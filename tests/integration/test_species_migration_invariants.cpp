#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/core/time_integration.hpp"
#include "cosmosim/gravity/tree_softening.hpp"

namespace {

cosmosim::core::HierarchicalTimeBinScheduler makeMigrationScheduler(
    const cosmosim::core::SimulationState& state) {
  std::uint8_t max_bin = 0;
  for (const std::uint8_t bin : state.particles.time_bin) {
    if (bin > max_bin) {
      max_bin = bin;
    }
  }
  cosmosim::core::HierarchicalTimeBinScheduler scheduler(max_bin);
  scheduler.reset(static_cast<std::uint32_t>(state.particles.size()), 0U, 0U);
  for (std::uint32_t i = 0; i < state.particles.size(); ++i) {
    scheduler.setElementBin(i, state.particles.time_bin[i], scheduler.currentTick());
  }
  return scheduler;
}


using cosmosim::core::HierarchicalTimeBinScheduler;
using cosmosim::core::ParticleMigrationCommit;
using cosmosim::core::ParticleMigrationRecord;
using cosmosim::core::ParticleReorderMode;
using cosmosim::core::ParticleSpecies;
using cosmosim::core::SimulationState;

constexpr std::uint32_t speciesTag(ParticleSpecies species) {
  return static_cast<std::uint32_t>(species);
}

void seedState(SimulationState& state) {
  state.resizeParticles(6);
  state.particle_sidecar.gravity_softening_comoving.resize(6);

  const std::array<std::uint64_t, 6> ids{1101, 1102, 1103, 1104, 1105, 1106};
  const std::array<std::uint32_t, 6> species{
      speciesTag(ParticleSpecies::kGas),
      speciesTag(ParticleSpecies::kDarkMatter),
      speciesTag(ParticleSpecies::kStar),
      speciesTag(ParticleSpecies::kGas),
      speciesTag(ParticleSpecies::kBlackHole),
      speciesTag(ParticleSpecies::kTracer)};

  for (std::size_t i = 0; i < state.particles.size(); ++i) {
    state.particle_sidecar.particle_id[i] = ids[i];
    state.particle_sidecar.species_tag[i] = species[i];
    state.particle_sidecar.owning_rank[i] = 0;
    state.particle_sidecar.sfc_key[i] = static_cast<std::uint64_t>(1000 - i * 10);
    state.particles.time_bin[i] = static_cast<std::uint8_t>(i % 3U);
    state.particles.mass_code[i] = 10.0 + static_cast<double>(i);
    state.particle_sidecar.gravity_softening_comoving[i] = 0.001 * static_cast<double>(i + 1);
  }

  state.star_particles.resize(1);
  state.star_particles.particle_index[0] = 2;
  state.star_particles.birth_mass_code[0] = 3.25;
  state.black_holes.resize(1);
  state.black_holes.particle_index[0] = 4;
  state.black_holes.host_cell_index[0] = 17;
  state.tracers.resize(1);
  state.tracers.particle_index[0] = 5;
  state.tracers.host_cell_index[0] = 9;

  state.species.count_by_species = {1, 2, 1, 1, 1};
  state.rebuildSpeciesIndex();
  assert(state.validateOwnershipInvariants());
}

std::uint32_t findParticleIndexById(const SimulationState& state, std::uint64_t id) {
  for (std::size_t i = 0; i < state.particles.size(); ++i) {
    if (state.particle_sidecar.particle_id[i] == id) {
      return static_cast<std::uint32_t>(i);
    }
  }
  assert(false && "particle id not found");
  return 0;
}

void test_species_migration_identity_invariants() {
  SimulationState state;
  seedState(state);

  const std::uint64_t migrated_id = 1103;
  const std::uint32_t outbound_index = findParticleIndexById(state, migrated_id);
  const std::array<std::uint32_t, 1> single_outbound{outbound_index};
  auto records = state.packParticleMigrationRecords(single_outbound, makeMigrationScheduler(state));
  assert(records.size() == 1);

  records[0].species_tag = speciesTag(ParticleSpecies::kGas);
  records[0].has_star_fields = false;
  records[0].has_gravity_softening_value = true;
  records[0].has_gravity_softening_override = true;
  records[0].gravity_softening_comoving = 0.0725;

  ParticleMigrationCommit commit;
  commit.world_rank = 0;
  commit.outbound_local_indices = {outbound_index};
  commit.inbound_records = records;
  state.commitParticleMigration(commit);

  const auto migrated_index = findParticleIndexById(state, migrated_id);
  assert(state.particle_sidecar.species_tag[migrated_index] == speciesTag(ParticleSpecies::kGas));
  assert((state.species.count_by_species == std::array<std::uint64_t, 5>{1, 3, 0, 1, 1}));
  assert(state.particles.time_bin[migrated_index] == records[0].time_bin);
  assert(state.particle_sidecar.gravity_softening_comoving[migrated_index] == 0.0725);
  assert(state.particle_sidecar.hasGravitySofteningOverride(migrated_index));
  assert(state.star_particles.size() == 0);

  const auto canonical = cosmosim::core::buildParticleReorderMap(state, ParticleReorderMode::kBySpecies);
  cosmosim::core::reorderParticles(state, canonical);
  assert(state.particle_sidecar.species_tag[findParticleIndexById(state, migrated_id)] == speciesTag(ParticleSpecies::kGas));

  state.resizeParticles(state.particles.size() + 2);
  for (std::size_t i = 0; i < state.particles.size(); ++i) {
    if (state.particle_sidecar.particle_id[i] == 0) {
      state.particle_sidecar.particle_id[i] = 2000 + static_cast<std::uint64_t>(i);
      state.particle_sidecar.species_tag[i] = speciesTag(ParticleSpecies::kGas);
      state.particle_sidecar.owning_rank[i] = 0;
      state.particle_sidecar.sfc_key[i] = 5000 + static_cast<std::uint64_t>(i);
      state.particles.time_bin[i] = 0;
      state.particles.mass_code[i] = 1.0;
      state.particle_sidecar.gravity_softening_comoving[i] = 0.003;
    }
  }
  state.species.count_by_species = {1, 5, 0, 1, 1};
  state.rebuildSpeciesIndex();

  assert(findParticleIndexById(state, migrated_id) < state.particles.size());
  assert(state.validateOwnershipInvariants());
}

void test_species_migration_sidecar_invariants() {
  SimulationState state;
  seedState(state);

  const std::uint32_t star_out = findParticleIndexById(state, 1103);
  const std::uint32_t bh_out = findParticleIndexById(state, 1105);

  const std::array<std::uint32_t, 2> outbound_indices{star_out, bh_out};
  auto outbound = state.packParticleMigrationRecords(outbound_indices, makeMigrationScheduler(state));
  assert(outbound.size() == 2);

  outbound[0].species_tag = speciesTag(ParticleSpecies::kTracer);
  outbound[0].has_star_fields = false;
  outbound[0].has_tracer_fields = true;
  outbound[0].tracer_fields.host_cell_index = 41;
  outbound[0].tracer_fields.parent_particle_id = 990001;
  outbound[0].tracer_fields.mass_fraction_of_host = 0.125;

  outbound[1].species_tag = speciesTag(ParticleSpecies::kGas);
  outbound[1].has_black_hole_fields = false;

  ParticleMigrationRecord inbound_star;
  inbound_star.particle_id = 7777;
  inbound_star.species_tag = speciesTag(ParticleSpecies::kStar);
  inbound_star.owning_rank = 0;
  inbound_star.sfc_key = 777;
  inbound_star.time_bin = 2;
  inbound_star.mass_code = 4.5;
  inbound_star.has_star_fields = true;
  inbound_star.star_fields.birth_mass_code = 8.0;
  inbound_star.has_gravity_softening_value = true;
  inbound_star.has_gravity_softening_override = true;
  inbound_star.gravity_softening_comoving = 0.025;

  ParticleMigrationCommit commit;
  commit.world_rank = 0;
  commit.outbound_local_indices = {star_out, bh_out};
  commit.inbound_records = std::vector<ParticleMigrationRecord>{outbound[0], outbound[1], inbound_star};
  state.commitParticleMigration(commit);

  assert(state.validateOwnershipInvariants());
  assert((state.species.count_by_species == std::array<std::uint64_t, 5>{1, 3, 1, 0, 2}));

  const auto tracer_idx = findParticleIndexById(state, 1103);
  assert(state.particle_sidecar.species_tag[tracer_idx] == speciesTag(ParticleSpecies::kTracer));
  bool found_tracer_sidecar = false;
  for (std::size_t row = 0; row < state.tracers.size(); ++row) {
    if (state.tracers.particle_index[row] == tracer_idx) {
      found_tracer_sidecar = true;
      assert(state.tracers.host_cell_index[row] == 41);
      assert(state.tracers.parent_particle_id[row] == 990001);
    }
  }
  assert(found_tracer_sidecar);

  const auto new_star_idx = findParticleIndexById(state, 7777);
  assert(state.particle_sidecar.species_tag[new_star_idx] == speciesTag(ParticleSpecies::kStar));
  bool found_star_sidecar = false;
  for (std::size_t row = 0; row < state.star_particles.size(); ++row) {
    if (state.star_particles.particle_index[row] == new_star_idx) {
      found_star_sidecar = true;
      assert(state.star_particles.birth_mass_code[row] == 8.0);
    }
  }
  assert(found_star_sidecar);
  assert(state.black_holes.size() == 0);
}


void test_species_migration_preserves_materialized_softening_without_promoting_default_to_override() {
  SimulationState state;
  seedState(state);

  const std::uint32_t outbound_index = findParticleIndexById(state, 1101);
  const std::array<std::uint32_t, 1> outbound_indices{outbound_index};
  auto records = state.packParticleMigrationRecords(outbound_indices, makeMigrationScheduler(state));
  assert(records.size() == 1);
  assert(records[0].has_gravity_softening_value);
  assert(!records[0].has_gravity_softening_override);

  records[0].species_tag = speciesTag(ParticleSpecies::kDarkMatter);

  ParticleMigrationCommit commit;
  commit.world_rank = 0;
  commit.outbound_local_indices = {outbound_index};
  commit.inbound_records = records;
  state.commitParticleMigration(commit);

  const std::uint32_t migrated_index = findParticleIndexById(state, 1101);
  assert(state.particle_sidecar.gravity_softening_comoving[migrated_index] == records[0].gravity_softening_comoving);
  assert(!state.particle_sidecar.hasGravitySofteningOverride(migrated_index));
  assert(state.validateOwnershipInvariants());
}

void test_species_migration_rejects_inbound_sidecar_mismatch() {
  SimulationState state;
  seedState(state);

  ParticleMigrationRecord invalid_star;
  invalid_star.particle_id = 8801;
  invalid_star.species_tag = speciesTag(ParticleSpecies::kStar);
  invalid_star.owning_rank = 0;
  invalid_star.mass_code = 1.0;
  invalid_star.time_bin = 0;
  invalid_star.has_star_fields = false;

  ParticleMigrationCommit commit;
  commit.world_rank = 0;
  commit.inbound_records = {invalid_star};

  bool threw = false;
  try {
    state.commitParticleMigration(commit);
  } catch (const std::invalid_argument&) {
    threw = true;
  }
  assert(threw);
  assert(state.validateOwnershipInvariants());

  ParticleMigrationRecord invalid_gas;
  invalid_gas.particle_id = 8802;
  invalid_gas.species_tag = speciesTag(ParticleSpecies::kGas);
  invalid_gas.owning_rank = 0;
  invalid_gas.mass_code = 1.0;
  invalid_gas.time_bin = 0;
  invalid_gas.has_tracer_fields = true;
  invalid_gas.tracer_fields.host_cell_index = 0;

  commit.inbound_records = {invalid_gas};
  threw = false;
  try {
    state.commitParticleMigration(commit);
  } catch (const std::invalid_argument&) {
    threw = true;
  }
  assert(threw);
  assert(state.validateOwnershipInvariants());

  ParticleMigrationRecord invalid_softening_override;
  invalid_softening_override.particle_id = 8803;
  invalid_softening_override.species_tag = speciesTag(ParticleSpecies::kDarkMatter);
  invalid_softening_override.owning_rank = 0;
  invalid_softening_override.mass_code = 1.0;
  invalid_softening_override.time_bin = 0;
  invalid_softening_override.has_gravity_softening_value = false;
  invalid_softening_override.has_gravity_softening_override = true;

  commit.inbound_records = {invalid_softening_override};
  threw = false;
  try {
    state.commitParticleMigration(commit);
  } catch (const std::invalid_argument&) {
    threw = true;
  }
  ParticleMigrationRecord missing_materialized_softening;
  missing_materialized_softening.particle_id = 8804;
  missing_materialized_softening.species_tag = speciesTag(ParticleSpecies::kDarkMatter);
  missing_materialized_softening.owning_rank = 0;
  missing_materialized_softening.mass_code = 1.0;
  missing_materialized_softening.time_bin = 0;
  missing_materialized_softening.has_gravity_softening_value = false;
  missing_materialized_softening.has_gravity_softening_override = false;

  commit.inbound_records = {missing_materialized_softening};
  threw = false;
  try {
    state.commitParticleMigration(commit);
  } catch (const std::invalid_argument&) {
    threw = true;
  }
  assert(threw);
  assert(state.validateOwnershipInvariants());
}

void test_species_migration_rejects_duplicate_final_particle_ids() {
  SimulationState state;
  seedState(state);

  ParticleMigrationRecord duplicate_kept_id;
  duplicate_kept_id.particle_id = 1101;
  duplicate_kept_id.species_tag = speciesTag(ParticleSpecies::kGas);
  duplicate_kept_id.owning_rank = 0;
  duplicate_kept_id.mass_code = 1.0;
  duplicate_kept_id.time_bin = 0;

  ParticleMigrationCommit commit;
  commit.world_rank = 0;
  commit.inbound_records = {duplicate_kept_id};

  bool threw = false;
  try {
    state.commitParticleMigration(commit);
  } catch (const std::invalid_argument&) {
    threw = true;
  }
  assert(threw);
  assert(state.validateUniqueParticleIds());
  assert(state.validateOwnershipInvariants());

  ParticleMigrationRecord duplicate_a;
  duplicate_a.particle_id = 9909;
  duplicate_a.species_tag = speciesTag(ParticleSpecies::kDarkMatter);
  duplicate_a.owning_rank = 0;
  duplicate_a.mass_code = 1.0;
  duplicate_a.time_bin = 0;
  ParticleMigrationRecord duplicate_b = duplicate_a;
  commit.inbound_records = {duplicate_a, duplicate_b};

  threw = false;
  try {
    state.commitParticleMigration(commit);
  } catch (const std::invalid_argument&) {
    threw = true;
  }
  assert(threw);
  assert(state.validateUniqueParticleIds());
  assert(state.validateOwnershipInvariants());
}

void test_species_migration_rejects_wrong_inbound_owner_rank() {
  SimulationState state;
  seedState(state);

  ParticleMigrationRecord wrong_owner;
  wrong_owner.particle_id = 8811;
  wrong_owner.species_tag = speciesTag(ParticleSpecies::kDarkMatter);
  wrong_owner.owning_rank = 1;
  wrong_owner.mass_code = 1.0;
  wrong_owner.time_bin = 0;

  ParticleMigrationCommit commit;
  commit.world_rank = 0;
  commit.inbound_records = {wrong_owner};

  bool threw = false;
  try {
    state.commitParticleMigration(commit);
  } catch (const std::invalid_argument&) {
    threw = true;
  }
  assert(threw);
  assert(state.validateOwnershipInvariants());
}

void test_species_migration_removes_stale_ghost_sidecar_rows_once() {
  SimulationState state;
  seedState(state);

  const std::uint32_t stale_star_index = findParticleIndexById(state, 1103);
  state.particle_sidecar.owning_rank[stale_star_index] = 1;
  const std::uint64_t generation_before = state.particleIndexGeneration();

  ParticleMigrationCommit commit;
  commit.world_rank = 0;
  commit.stale_local_ghost_indices = {stale_star_index};
  state.commitParticleMigration(commit);

  assert(state.particleIndexGeneration() == generation_before + 1);
  assert(state.star_particles.size() == 0);
  assert((state.species.count_by_species == std::array<std::uint64_t, 5>{1, 2, 0, 1, 1}));
  assert(state.validateOwnershipInvariants());
}

void test_species_migration_rejects_owned_or_duplicate_stale_ghost_indices() {
  SimulationState state;
  seedState(state);

  ParticleMigrationCommit commit;
  commit.world_rank = 0;
  commit.stale_local_ghost_indices = {findParticleIndexById(state, 1103)};

  bool threw = false;
  try {
    state.commitParticleMigration(commit);
  } catch (const std::invalid_argument&) {
    threw = true;
  }
  assert(threw);
  assert(state.validateOwnershipInvariants());

  const std::uint32_t remote_bh_index = findParticleIndexById(state, 1105);
  state.particle_sidecar.owning_rank[remote_bh_index] = 1;
  commit.stale_local_ghost_indices = {remote_bh_index, remote_bh_index};
  threw = false;
  try {
    state.commitParticleMigration(commit);
  } catch (const std::invalid_argument&) {
    threw = true;
  }
  assert(threw);
}

void test_species_migration_rejects_stale_sidecar_indices_before_commit() {
  SimulationState state;
  seedState(state);
  state.tracers.particle_index[0] = static_cast<std::uint32_t>(state.particles.size() + 8);

  ParticleMigrationCommit commit;
  commit.world_rank = 0;

  bool threw = false;
  try {
    state.commitParticleMigration(commit);
  } catch (const std::runtime_error&) {
    threw = true;
  }
  assert(threw);
}

void test_species_migration_softening_timestep_invariants() {
  SimulationState state;
  seedState(state);

  const std::uint64_t migrated_id = 1102;
  const std::uint32_t outbound_index = findParticleIndexById(state, migrated_id);
  const std::array<std::uint32_t, 1> outbound_indices{outbound_index};
  auto records = state.packParticleMigrationRecords(outbound_indices, makeMigrationScheduler(state));
  records[0].species_tag = speciesTag(ParticleSpecies::kStar);
  records[0].has_star_fields = true;
  records[0].star_fields.birth_mass_code = 12.0;
  records[0].has_gravity_softening_value = true;
  records[0].has_gravity_softening_override = true;
  records[0].gravity_softening_comoving = 0.077;

  ParticleMigrationCommit commit;
  commit.world_rank = 0;
  commit.outbound_local_indices = {outbound_index};
  commit.inbound_records = records;
  state.commitParticleMigration(commit);

  const std::uint32_t migrated_index = findParticleIndexById(state, migrated_id);
  assert(state.particles.time_bin[migrated_index] == records[0].time_bin);

  HierarchicalTimeBinScheduler scheduler(3);
  scheduler.reset(static_cast<std::uint32_t>(state.particles.size()), 0, 0);
  for (std::uint32_t i = 0; i < static_cast<std::uint32_t>(state.particles.size()); ++i) {
    scheduler.setElementBin(i, state.particles.time_bin[i], 0);
  }
  const auto active = scheduler.beginSubstep();
  assert(!active.empty());
  assert(std::find(active.begin(), active.end(), migrated_index) != active.end());

  cosmosim::gravity::TreeSofteningPolicy fallback;
  fallback.epsilon_comoving = 0.5;
  cosmosim::gravity::TreeSofteningView softening_view{
      .source_species_tag = std::span<const std::uint32_t>(state.particle_sidecar.species_tag.data(), state.particle_sidecar.species_tag.size()),
      .source_particle_epsilon_comoving = std::span<const double>(state.particle_sidecar.gravity_softening_comoving.data(),
                                                                  state.particle_sidecar.gravity_softening_comoving.size()),
      .source_particle_epsilon_override_mask = std::span<const std::uint8_t>(
          state.particle_sidecar.has_gravity_softening_override.data(),
          state.particle_sidecar.has_gravity_softening_override.size()),
      .target_particle_epsilon_comoving = std::span<const double>(),
      .species_policy = {.epsilon_comoving_by_species = {0.002, 0.004, 0.008, 0.016, 0.032}, .enabled = true},
  };
  const double resolved = cosmosim::gravity::resolveSourceSofteningEpsilon(migrated_index, fallback, softening_view);
  assert(resolved == 0.077);

  const auto reorder = cosmosim::core::buildParticleReorderMap(state, ParticleReorderMode::kBySfcKey);
  cosmosim::core::reorderParticles(state, reorder);
  const std::uint32_t reordered_index = findParticleIndexById(state, migrated_id);
  assert(state.particles.time_bin[reordered_index] == records[0].time_bin);
  assert(state.particle_sidecar.gravity_softening_comoving[reordered_index] == 0.077);
  assert(state.particle_sidecar.hasGravitySofteningOverride(reordered_index));
}


void test_particle_migration_carries_scheduler_authority_fields() {
  SimulationState state;
  seedState(state);

  HierarchicalTimeBinScheduler scheduler(4);
  scheduler.reset(static_cast<std::uint32_t>(state.particles.size()), 0, 5);
  for (std::uint32_t i = 0; i < static_cast<std::uint32_t>(state.particles.size()); ++i) {
    scheduler.setElementBin(i, static_cast<std::uint8_t>((i + 1U) % 4U), 5);
  }
  const std::uint32_t outbound_index = findParticleIndexById(state, 1105);
  scheduler.submitCandidateBin(outbound_index, 1, cosmosim::core::TimeStepCandidateSource::kUserClamp, "migration-authority-test");
  (void)scheduler.reconcileCandidateTransitions();

  const std::array<std::uint32_t, 1> outbound{outbound_index};
  auto records = state.packParticleMigrationRecords(outbound, scheduler);
  assert(records[0].has_scheduler_fields);
  assert(records[0].scheduler_fields.bin_index == records[0].time_bin);
  assert(records[0].scheduler_fields.pending_bin_index == 1U);

  const std::array<std::uint64_t, 1> destination_ids{records[0].particle_id};
  const auto rebuilt = cosmosim::core::rebuildSchedulerPersistentStateFromMigrationRecords(
      scheduler.currentTick(), scheduler.maxBin(), records, destination_ids);
  assert(rebuilt.bin_index.size() == 1);
  assert(rebuilt.bin_index[0] == records[0].scheduler_fields.bin_index);
  assert(rebuilt.next_activation_tick[0] == records[0].scheduler_fields.next_activation_tick);
  assert(rebuilt.pending_bin_index[0] == 1U);

  records[0].time_bin = static_cast<std::uint8_t>(records[0].time_bin + 1U);
  bool stale_mirror_threw = false;
  try {
    (void)cosmosim::core::rebuildSchedulerPersistentStateFromMigrationRecords(
        scheduler.currentTick(), scheduler.maxBin(), records, destination_ids);
  } catch (const std::invalid_argument&) {
    stale_mirror_threw = true;
  }
  assert(stale_mirror_threw);

  records[0].time_bin = records[0].scheduler_fields.bin_index;
  records[0].has_scheduler_fields = false;
  bool missing_authority_threw = false;
  try {
    (void)cosmosim::core::rebuildSchedulerPersistentStateFromMigrationRecords(
        scheduler.currentTick(), scheduler.maxBin(), records, destination_ids);
  } catch (const std::invalid_argument&) {
    missing_authority_threw = true;
  }
  assert(missing_authority_threw);
}

void test_gas_cell_migration_rebuilds_hydro_fields_by_particle_id() {
  SimulationState state;
  state.resizeParticles(4);
  state.resizeCells(3);
  state.particle_sidecar.gravity_softening_comoving.resize(4);

  const std::array<std::uint64_t, 4> ids{501, 900, 502, 503};
  const std::array<std::uint32_t, 4> species{
      speciesTag(ParticleSpecies::kGas),
      speciesTag(ParticleSpecies::kDarkMatter),
      speciesTag(ParticleSpecies::kGas),
      speciesTag(ParticleSpecies::kGas)};
  for (std::size_t i = 0; i < state.particles.size(); ++i) {
    state.particle_sidecar.particle_id[i] = ids[i];
    state.particle_sidecar.species_tag[i] = species[i];
    state.particle_sidecar.owning_rank[i] = 0;
    state.particle_sidecar.sfc_key[i] = static_cast<std::uint64_t>(100 + i);
    state.particles.mass_code[i] = 1.0 + static_cast<double>(i);
    state.particles.time_bin[i] = static_cast<std::uint8_t>(i);
    state.particle_sidecar.gravity_softening_comoving[i] = 0.01;
  }
  state.species.count_by_species = {1, 3, 0, 0, 0};
  state.rebuildSpeciesIndex();
  state.refreshGasCellIdentityFromParticleOrder();
  for (std::size_t cell = 0; cell < state.cells.size(); ++cell) {
    state.cells.mass_code[cell] = 10.0 + static_cast<double>(cell);
    state.cells.time_bin[cell] = static_cast<std::uint8_t>(cell + 1U);
    state.gas_cells.density_code[cell] = 1000.0 + static_cast<double>(cell);
    state.gas_cells.internal_energy_code[cell] = 2000.0 + static_cast<double>(cell);
  }
  cosmosim::core::requireParticleBoundGasCellContract(state, "migration gas fixture");

  struct GasFieldRecord {
    double mass_code = 0.0;
    std::uint8_t time_bin = 0;
    double density_code = 0.0;
    double internal_energy_code = 0.0;
  };
  std::unordered_map<std::uint64_t, GasFieldRecord> gas_fields_by_id;
  for (std::uint32_t cell = 0; cell < state.cells.size(); ++cell) {
    gas_fields_by_id.emplace(
        state.parentParticleIdForGasCellRow(cell).value(),
        GasFieldRecord{
            .mass_code = state.cells.mass_code[cell],
            .time_bin = state.cells.time_bin[cell],
            .density_code = state.gas_cells.density_code[cell],
            .internal_energy_code = state.gas_cells.internal_energy_code[cell],
        });
  }

  const std::uint32_t outbound_index = findParticleIndexById(state, 502);
  std::array<std::uint32_t, 1> outbound{outbound_index};
  auto records = state.packParticleMigrationRecords(outbound, makeMigrationScheduler(state));
  ParticleMigrationCommit commit;
  commit.world_rank = 0;
  commit.outbound_local_indices = {outbound_index};
  commit.inbound_records = records;
  state.commitParticleMigration(commit);

  (void)gas_fields_by_id;
  cosmosim::core::requireParticleBoundGasCellContract(state, "migration gas commit");
  assert(state.gasCellRowForParticleId(501) == 0);
  assert(state.gasCellRowForParticleId(503) == 1);
  assert(state.gasCellRowForParticleId(502) == 2);
  assert(state.gas_cells.density_code[state.gasCellRowForParticleId(501)] == 1000.0);
  assert(state.gas_cells.density_code[state.gasCellRowForParticleId(503)] == 1002.0);
  assert(state.gas_cells.density_code[state.gasCellRowForParticleId(502)] == 1001.0);
  assert(state.cells.time_bin[state.gasCellRowForParticleId(501)] == 1U);
  assert(state.cells.time_bin[state.gasCellRowForParticleId(503)] == 3U);
  assert(state.cells.time_bin[state.gasCellRowForParticleId(502)] == 2U);
  assert(state.gas_cells.internal_energy_code[state.gasCellRowForParticleId(502)] == 2001.0);
  assert(state.validateOwnershipInvariants());
}


void test_particle_indexed_module_sidecars_migrate_with_particles() {
  SimulationState state;
  state.resizeParticles(3);
  state.particle_sidecar.gravity_softening_comoving.resize(3);
  const std::array<std::uint64_t, 3> ids{7001, 7002, 7003};
  for (std::size_t i = 0; i < state.particles.size(); ++i) {
    state.particle_sidecar.particle_id[i] = ids[i];
    state.particle_sidecar.species_tag[i] = speciesTag(ParticleSpecies::kDarkMatter);
    state.particle_sidecar.owning_rank[i] = 0;
    state.particle_sidecar.sfc_key[i] = 100 + i;
    state.particles.mass_code[i] = 1.0 + static_cast<double>(i);
    state.particles.time_bin[i] = static_cast<std::uint8_t>(i);
    state.particle_sidecar.gravity_softening_comoving[i] = 0.01;
  }
  state.species.count_by_species = {3, 0, 0, 0, 0};
  state.rebuildSpeciesIndex();

  cosmosim::core::ModuleSidecarBlock indexed;
  indexed.module_name = "module_particle_payload";
  indexed.schema_version = 3;
  indexed.particle_indexed = true;
  indexed.row_stride_bytes = 2;
  indexed.particle_id_by_row = {7001, 7002, 7003};
  indexed.payload = {
      std::byte{0xA1}, std::byte{0xA2},
      std::byte{0xB1}, std::byte{0xB2},
      std::byte{0xC1}, std::byte{0xC2}};
  state.sidecars.upsert(indexed);

  cosmosim::core::ModuleSidecarBlock run_level;
  run_level.module_name = "run_level_payload";
  run_level.schema_version = 1;
  run_level.payload = {std::byte{0x44}, std::byte{0x55}};
  state.sidecars.upsert(run_level);

  const std::array<std::uint32_t, 1> outbound{1};
  auto records = state.packParticleMigrationRecords(outbound, makeMigrationScheduler(state));
  assert(records.size() == 1);
  assert(records[0].module_sidecar_payloads.size() == 1);
  assert(records[0].module_sidecar_payloads[0].module_name == "module_particle_payload");
  assert(records[0].module_sidecar_payloads[0].schema_version == 3U);
  assert(records[0].module_sidecar_payloads[0].payload[0] == std::byte{0xB1});
  assert(records[0].module_sidecar_payloads[0].payload[1] == std::byte{0xB2});

  ParticleMigrationCommit commit;
  commit.world_rank = 0;
  commit.outbound_local_indices = {1};
  commit.inbound_records = records;
  state.commitParticleMigration(commit);

  const auto* rebuilt = state.sidecars.find("module_particle_payload");
  assert(rebuilt != nullptr);
  assert(rebuilt->isParticleIndexed());
  assert(rebuilt->particle_id_by_row == std::vector<std::uint64_t>({7001, 7003, 7002}));
  assert(rebuilt->payload == std::vector<std::byte>({
      std::byte{0xA1}, std::byte{0xA2},
      std::byte{0xC1}, std::byte{0xC2},
      std::byte{0xB1}, std::byte{0xB2}}));

  const auto* preserved = state.sidecars.find("run_level_payload");
  assert(preserved != nullptr);
  assert(!preserved->isParticleIndexed());
  assert(preserved->payload == std::vector<std::byte>({std::byte{0x44}, std::byte{0x55}}));
}

}  // namespace

int main() {
  test_species_migration_identity_invariants();
  test_species_migration_sidecar_invariants();
  test_species_migration_softening_timestep_invariants();
  test_species_migration_preserves_materialized_softening_without_promoting_default_to_override();
  test_species_migration_rejects_inbound_sidecar_mismatch();
  test_species_migration_rejects_duplicate_final_particle_ids();
  test_species_migration_rejects_wrong_inbound_owner_rank();
  test_species_migration_removes_stale_ghost_sidecar_rows_once();
  test_species_migration_rejects_owned_or_duplicate_stale_ghost_indices();
  test_species_migration_rejects_stale_sidecar_indices_before_commit();
  test_particle_migration_carries_scheduler_authority_fields();
  test_gas_cell_migration_rebuilds_hydro_fields_by_particle_id();
  test_particle_indexed_module_sidecars_migrate_with_particles();
  return 0;
}
