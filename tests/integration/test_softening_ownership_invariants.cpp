#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <unordered_map>
#include <vector>

#include "cosmosim/analysis/diagnostics.hpp"
#include "cosmosim/core/build_config.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/gravity/tree_softening.hpp"
#include "cosmosim/io/restart_checkpoint.hpp"

namespace {

using cosmosim::core::ParticleMigrationCommit;
using cosmosim::core::ParticleSpecies;
using cosmosim::core::SimulationState;

constexpr std::uint32_t speciesTag(ParticleSpecies species) {
  return static_cast<std::uint32_t>(species);
}

void rebuildSpeciesLedger(SimulationState& state) {
  state.species.count_by_species.fill(0);
  for (const auto tag : state.particle_sidecar.species_tag) {
    ++state.species.count_by_species[tag];
  }
  state.rebuildSpeciesIndex();
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

cosmosim::gravity::TreeSofteningSpeciesPolicy makeSpeciesPolicy() {
  cosmosim::gravity::TreeSofteningSpeciesPolicy policy;
  policy.enabled = true;
  policy.epsilon_comoving_by_species = {0.020, 0.030, 0.040, 0.050, 0.060};
  return policy;
}

SimulationState seedStateWithSelectiveOverrides() {
  SimulationState state;
  state.resizeParticles(6);

  const std::array<std::uint64_t, 6> ids{501, 502, 503, 504, 505, 506};
  const std::array<std::uint32_t, 6> species{
      speciesTag(ParticleSpecies::kDarkMatter),
      speciesTag(ParticleSpecies::kGas),
      speciesTag(ParticleSpecies::kStar),
      speciesTag(ParticleSpecies::kGas),
      speciesTag(ParticleSpecies::kBlackHole),
      speciesTag(ParticleSpecies::kTracer)};

  for (std::size_t i = 0; i < state.particles.size(); ++i) {
    state.particle_sidecar.particle_id[i] = ids[i];
    state.particle_sidecar.species_tag[i] = species[i];
    state.particle_sidecar.owning_rank[i] = 0;
    state.particle_sidecar.sfc_key[i] = 1000U - static_cast<std::uint64_t>(i * 100U);
  }

  state.particle_sidecar.gravity_softening_comoving.resize(state.particles.size(), 0.0);
  const auto species_policy = makeSpeciesPolicy();
  for (std::size_t i = 0; i < state.particles.size(); ++i) {
    const std::size_t sp = static_cast<std::size_t>(state.particle_sidecar.species_tag[i]);
    state.particle_sidecar.gravity_softening_comoving[i] = species_policy.epsilon_comoving_by_species[sp];
  }
  // Selective per-particle overrides over species defaults.
  state.particle_sidecar.setGravitySofteningOverride(1, 0.03125);  // gas override
  state.particle_sidecar.setGravitySofteningOverride(4, 0.07250);  // BH override

  state.star_particles.resize(1);
  state.star_particles.particle_index[0] = 2;
  state.black_holes.resize(1);
  state.black_holes.particle_index[0] = 4;
  state.tracers.resize(1);
  state.tracers.particle_index[0] = 5;

  rebuildSpeciesLedger(state);
  assert(state.validateOwnershipInvariants());
  return state;
}

void assertSofteningById(const SimulationState& state, std::uint64_t id, double expected) {
  const auto index = findParticleIndexById(state, id);
  assert(std::abs(state.particle_sidecar.gravity_softening_comoving[index] - expected) < 1.0e-15);
}

void test_softening_priority_invariants() {
  const auto species_policy = makeSpeciesPolicy();
  cosmosim::gravity::TreeSofteningPolicy global_policy;
  global_policy.epsilon_comoving = 0.125;

  const std::array<std::uint32_t, 3> source_species{
      speciesTag(ParticleSpecies::kGas),
      speciesTag(ParticleSpecies::kStar),
      99U};

  // Global-only fallback.
  {
    const cosmosim::gravity::TreeSofteningView view{
        .source_species_tag = source_species,
        .source_particle_epsilon_comoving = {},
        .target_particle_epsilon_comoving = {},
        .species_policy = {.enabled = false},
    };
    assert(std::abs(cosmosim::gravity::resolveSourceSofteningEpsilon(0, global_policy, view) - 0.125) < 1.0e-15);
  }

  // Species-default fallback.
  {
    const cosmosim::gravity::TreeSofteningView view{
        .source_species_tag = source_species,
        .source_particle_epsilon_comoving = {},
        .target_particle_epsilon_comoving = {},
        .species_policy = species_policy,
    };
    assert(std::abs(cosmosim::gravity::resolveSourceSofteningEpsilon(0, global_policy, view) - 0.030) < 1.0e-15);
    assert(std::abs(cosmosim::gravity::resolveSourceSofteningEpsilon(1, global_policy, view) - 0.040) < 1.0e-15);
    assert(std::abs(cosmosim::gravity::resolveSourceSofteningEpsilon(2, global_policy, view) - 0.125) < 1.0e-15);
  }

  // Per-particle override precedence.
  {
    const std::array<double, 3> per_particle_eps{0.333, 0.222, 0.111};
    const cosmosim::gravity::TreeSofteningView view{
        .source_species_tag = source_species,
        .source_particle_epsilon_comoving = per_particle_eps,
        .target_particle_epsilon_comoving = {},
        .species_policy = species_policy,
    };
    assert(std::abs(cosmosim::gravity::resolveSourceSofteningEpsilon(0, global_policy, view) - 0.333) < 1.0e-15);
    assert(std::abs(cosmosim::gravity::resolveSourceSofteningEpsilon(1, global_policy, view) - 0.222) < 1.0e-15);
    assert(std::abs(cosmosim::gravity::resolveSourceSofteningEpsilon(2, global_policy, view) - 0.111) < 1.0e-15);
  }

  // Materialized per-particle values are not authoritative overrides when an override mask is present.
  {
    const std::array<double, 3> materialized_eps{0.333, 0.222, 0.111};
    const std::array<std::uint8_t, 3> override_mask{1U, 0U, 0U};
    const cosmosim::gravity::TreeSofteningView view{
        .source_species_tag = source_species,
        .source_particle_epsilon_comoving = materialized_eps,
        .source_particle_epsilon_override_mask = override_mask,
        .target_particle_epsilon_comoving = {},
        .species_policy = species_policy,
    };
    assert(std::abs(cosmosim::gravity::resolveSourceSofteningEpsilon(0, global_policy, view) - 0.333) < 1.0e-15);
    assert(std::abs(cosmosim::gravity::resolveSourceSofteningEpsilon(1, global_policy, view) - 0.040) < 1.0e-15);
    assert(std::abs(cosmosim::gravity::resolveSourceSofteningEpsilon(2, global_policy, view) - 0.125) < 1.0e-15);
  }

  // Diagnostics are observers only.
  SimulationState state = seedStateWithSelectiveOverrides();
  const auto before = state.particle_sidecar.gravity_softening_comoving;

  cosmosim::core::SimulationConfig cfg;
  cfg.analysis.enable_diagnostics = true;
  cosmosim::analysis::DiagnosticsEngine engine(cfg);
  const auto health = engine.computeRunHealth(state);
  assert(health.gravity_softening_sidecar_size_ok);
  assert(state.particle_sidecar.gravity_softening_comoving == before);
}

void test_softening_override_resize_reorder_preservation() {
  SimulationState state = seedStateWithSelectiveOverrides();

  const double id502_eps = state.particle_sidecar.gravity_softening_comoving[findParticleIndexById(state, 502)];
  const double id505_eps = state.particle_sidecar.gravity_softening_comoving[findParticleIndexById(state, 505)];

  // Grow and confirm pre-existing rows preserved.
  const std::size_t original_count = state.particles.size();
  state.resizeParticles(original_count + 2);
  for (std::size_t i = original_count; i < state.particles.size(); ++i) {
    state.particle_sidecar.particle_id[i] = 800 + static_cast<std::uint64_t>(i);
    state.particle_sidecar.species_tag[i] = speciesTag(ParticleSpecies::kGas);
    state.particle_sidecar.owning_rank[i] = 0;
    state.particle_sidecar.sfc_key[i] = 5000U + static_cast<std::uint64_t>(i);
    state.particle_sidecar.gravity_softening_comoving[i] = 0.030;
  }
  rebuildSpeciesLedger(state);

  assertSofteningById(state, 502, id502_eps);
  assertSofteningById(state, 505, id505_eps);

  // Reorder and ensure override follows particle identity.
  const auto reorder = cosmosim::core::buildParticleReorderMap(state, cosmosim::core::ParticleReorderMode::kBySfcKey);
  cosmosim::core::reorderParticles(state, reorder);
  assertSofteningById(state, 502, id502_eps);
  assertSofteningById(state, 505, id505_eps);

  // Species migration preserves explicit overrides.
  const std::uint32_t outbound = findParticleIndexById(state, 505);
  auto packet = state.packParticleMigrationRecords(std::array<std::uint32_t, 1>{outbound});
  assert(packet.size() == 1);
  packet[0].species_tag = speciesTag(ParticleSpecies::kGas);
  packet[0].has_black_hole_fields = false;
  packet[0].has_gravity_softening_value = true;
  packet[0].has_gravity_softening_override = true;
  packet[0].gravity_softening_comoving = id505_eps;

  ParticleMigrationCommit commit;
  commit.world_rank = 0;
  commit.outbound_local_indices = {outbound};
  commit.inbound_records = packet;
  state.commitParticleMigration(commit);
  assertSofteningById(state, 505, id505_eps);

  // Active-set extraction lane mirrors authoritative sidecar data and follows reorder identity.
  const std::vector<std::uint32_t> active_indices{
      findParticleIndexById(state, 502),
      findParticleIndexById(state, 505)};
  std::vector<double> active_softening(active_indices.size(), 0.0);
  for (std::size_t i = 0; i < active_indices.size(); ++i) {
    active_softening[i] = state.particle_sidecar.gravity_softening_comoving[active_indices[i]];
  }
  assert(std::abs(active_softening[0] - id502_eps) < 1.0e-15);
  assert(std::abs(active_softening[1] - id505_eps) < 1.0e-15);

  // Shrink after reorder and ensure retained rows keep overrides.
  state.resizeParticles(6);
  assert(state.particle_sidecar.gravity_softening_comoving.size() == 6);
  assert(findParticleIndexById(state, 502) < state.particles.size());
  assertSofteningById(state, 502, id502_eps);
}

void test_softening_override_restart_roundtrip() {
#if COSMOSIM_ENABLE_HDF5
  SimulationState state = seedStateWithSelectiveOverrides();
  cosmosim::core::IntegratorState integrator_state;
  cosmosim::core::HierarchicalTimeBinScheduler scheduler(2);
  scheduler.reset(static_cast<std::uint32_t>(state.particles.size()), 0, 0);

  cosmosim::io::RestartWritePayload payload;
  payload.state = &state;
  payload.integrator_state = &integrator_state;
  payload.scheduler = &scheduler;
  payload.normalized_config_text = "schema_version = 1\n[mode]\nmode = zoom_in\n";
  payload.normalized_config_hash_hex = "a";

  const std::filesystem::path restart_path =
      std::filesystem::temp_directory_path() / "cosmosim_softening_ownership_roundtrip.hdf5";

  cosmosim::io::writeRestartCheckpointHdf5(restart_path, payload);
  const auto restored = cosmosim::io::readRestartCheckpointHdf5(restart_path);

  assert(restored.state.validateOwnershipInvariants());
  assert(restored.state.particle_sidecar.gravity_softening_comoving == state.particle_sidecar.gravity_softening_comoving);
  assert(restored.state.particle_sidecar.has_gravity_softening_override == state.particle_sidecar.has_gravity_softening_override);

  const auto species_policy = makeSpeciesPolicy();
  cosmosim::gravity::TreeSofteningPolicy global_policy;
  global_policy.epsilon_comoving = 0.125;
  const cosmosim::gravity::TreeSofteningView view{
      .source_species_tag = std::span<const std::uint32_t>(
          restored.state.particle_sidecar.species_tag.data(),
          restored.state.particle_sidecar.species_tag.size()),
      .source_particle_epsilon_comoving = std::span<const double>(
          restored.state.particle_sidecar.gravity_softening_comoving.data(),
          restored.state.particle_sidecar.gravity_softening_comoving.size()),
      .source_particle_epsilon_override_mask = std::span<const std::uint8_t>(
          restored.state.particle_sidecar.has_gravity_softening_override.data(),
          restored.state.particle_sidecar.has_gravity_softening_override.size()),
      .target_particle_epsilon_comoving = {},
      .species_policy = species_policy,
  };

  const auto id505 = findParticleIndexById(restored.state, 505);
  assert(std::abs(cosmosim::gravity::resolveSourceSofteningEpsilon(id505, global_policy, view) - 0.07250) < 1.0e-15);
#else
  // No-op in non-HDF5 builds; restart roundtrip is covered by HDF5-enabled test presets.
#endif
}

}  // namespace

int main() {
  test_softening_priority_invariants();
  test_softening_override_resize_reorder_preservation();
  test_softening_override_restart_roundtrip();
  return 0;
}
