#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <unordered_map>
#include <vector>

#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/core/time_integration.hpp"
#include "cosmosim/gravity/tree_softening.hpp"

namespace {

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
  auto records = state.packParticleMigrationRecords(single_outbound);
  assert(records.size() == 1);

  records[0].species_tag = speciesTag(ParticleSpecies::kGas);
  records[0].has_star_fields = false;
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
  auto outbound = state.packParticleMigrationRecords(outbound_indices);
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

void test_species_migration_softening_timestep_invariants() {
  SimulationState state;
  seedState(state);

  const std::uint64_t migrated_id = 1102;
  const std::uint32_t outbound_index = findParticleIndexById(state, migrated_id);
  const std::array<std::uint32_t, 1> outbound_indices{outbound_index};
  auto records = state.packParticleMigrationRecords(outbound_indices);
  records[0].species_tag = speciesTag(ParticleSpecies::kStar);
  records[0].has_star_fields = true;
  records[0].star_fields.birth_mass_code = 12.0;
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
}

}  // namespace

int main() {
  test_species_migration_identity_invariants();
  test_species_migration_sidecar_invariants();
  test_species_migration_softening_timestep_invariants();
  return 0;
}
