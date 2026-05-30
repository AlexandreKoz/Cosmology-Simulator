#include <cassert>
#include <cstdint>
#include <filesystem>
#include <utility>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/core/time_integration.hpp"
#include "cosmosim/io/restart_checkpoint.hpp"
#include "restart_equivalence_harness.hpp"

namespace {

cosmosim::core::SimulationState makeToyDmState() {
  cosmosim::core::SimulationState state;
  state.resizeParticles(6);
  state.metadata.run_name = "restart_equivalence_harness";
  state.metadata.normalized_config_hash_hex = "stage8_restart_equivalence_harness";
  state.metadata.snapshot_stem = "snapshot";
  state.metadata.restart_stem = "restart";
  state.species.count_by_species[static_cast<std::size_t>(cosmosim::core::ParticleSpecies::kDarkMatter)] =
      state.particles.size();
  for (std::size_t pidx = 0; pidx < state.particles.size(); ++pidx) {
    state.particles.position_x_comoving[pidx] = 0.1 * static_cast<double>(pidx + 1U);
    state.particles.position_y_comoving[pidx] = 0.2 * static_cast<double>(pidx + 1U);
    state.particles.position_z_comoving[pidx] = 0.3 * static_cast<double>(pidx + 1U);
    state.particles.velocity_x_peculiar[pidx] = 0.01 * static_cast<double>(pidx + 1U);
    state.particles.velocity_y_peculiar[pidx] = -0.02 * static_cast<double>(pidx + 1U);
    state.particles.velocity_z_peculiar[pidx] = 0.03 * static_cast<double>(pidx + 1U);
    state.particles.mass_code[pidx] = 10.0 + static_cast<double>(pidx);
    state.particle_sidecar.particle_id[pidx] = 1000U + pidx;
    state.particle_sidecar.species_tag[pidx] = static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kDarkMatter);
    state.particle_sidecar.owning_rank[pidx] = 0;
    state.particle_sidecar.setGravitySofteningOverride(pidx, 0.01);
  }
  state.rebuildSpeciesIndex();
  return state;
}

cosmosim::core::HierarchicalTimeBinScheduler makeToyScheduler(std::uint32_t element_count) {
  cosmosim::core::HierarchicalTimeBinScheduler scheduler(2);
  scheduler.reset(element_count, 0, 0);
  for (std::uint32_t pidx = 0; pidx < element_count; ++pidx) {
    scheduler.setElementBin(pidx, static_cast<std::uint8_t>(pidx % 3U), scheduler.currentTick());
  }
  return scheduler;
}

cosmosim::core::IntegratorState makeToyIntegratorState() {
  cosmosim::core::IntegratorState integrator_state;
  integrator_state.current_time_code = 0.25;
  integrator_state.current_scale_factor = 0.5;
  integrator_state.current_redshift = 1.0;
  integrator_state.time_si_per_code = 1.0;
  integrator_state.current_boundary_kind = cosmosim::core::StepBoundaryKind::kGlobalSynchronizationPoint;
  integrator_state.last_completed_boundary_kind = cosmosim::core::StepBoundaryKind::kCheckpointPoint;
  integrator_state.last_completed_restart_safe = true;
  integrator_state.time_bins.hierarchical_enabled = true;
  integrator_state.time_bins.max_bin = 2;
  integrator_state.pm_refresh_enabled = true;
  integrator_state.pm_long_range_field_valid = true;
  integrator_state.pm_sync_state.reset(1);
  return integrator_state;
}

cosmosim::io::OutputCadencePersistentState makeOutputCadenceState() {
  cosmosim::io::OutputCadencePersistentState output_state;
  output_state.output_enabled = true;
  output_state.write_restarts = true;
  output_state.snapshot_due = false;
  output_state.checkpoint_due = false;
  output_state.snapshot_interval_steps = 3;
  output_state.next_snapshot_step_index = 3;
  output_state.snapshot_stem = "snapshot";
  output_state.restart_stem = "restart";
  return output_state;
}

cosmosim::io::StochasticPersistentState makeStochasticState() {
  cosmosim::io::StochasticPersistentState stochastic_state;
  stochastic_state.modules.push_back(cosmosim::io::StochasticModulePersistentState{
      .module_name = "star_formation",
      .schema_version = 1,
      .rng_policy = "stateless_splitmix64(seed,step_index,cell_index,rank_local_seed_offset)",
      .random_seed = 123456789ull,
      .rank_local_seed_offset = 0,
      .last_committed_step_index = 0,
      .deterministic_from_serialized_inputs = true,
  });
  stochastic_state.modules.push_back(cosmosim::io::StochasticModulePersistentState{
      .module_name = "stellar_feedback",
      .schema_version = 1,
      .rng_policy = "stateless_splitmix64(seed,step_index,star_index)",
      .random_seed = 42424242ull,
      .rank_local_seed_offset = 0,
      .last_committed_step_index = 0,
      .deterministic_from_serialized_inputs = true,
  });
  return stochastic_state;
}

void testRestartEquivalenceHarness() {
#if COSMOSIM_ENABLE_HDF5
  cosmosim::tests::RestartEquivalenceScenario scenario;
  scenario.initial_state = makeToyDmState();
  scenario.initial_integrator_state = makeToyIntegratorState();
  scenario.initial_scheduler = makeToyScheduler(static_cast<std::uint32_t>(scenario.initial_state.particles.size()));
  scenario.initial_output_cadence_state = makeOutputCadenceState();
  scenario.initial_stochastic_state = makeStochasticState();
  scenario.total_steps = 10;
  scenario.restart_step = 4;
  scenario.restart_path =
      std::filesystem::temp_directory_path() / "cosmosim_restart_equivalence_harness.hdf5";

  const cosmosim::tests::RestartEquivalenceResult result =
      cosmosim::tests::runRestartEquivalenceScenario(std::move(scenario));
  assert(result.direct_integrator_state.step_index == 10);
  assert(result.restarted_integrator_state.step_index == 10);
  assert(result.direct_scheduler_state.current_tick == result.restarted_scheduler_state.current_tick);
  std::filesystem::remove(std::filesystem::temp_directory_path() / "cosmosim_restart_equivalence_harness.hdf5");
#endif
}

}  // namespace

int main() {
  testRestartEquivalenceHarness();
  return 0;
}
