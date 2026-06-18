#pragma once

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <string>
#include <utility>

#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/core/time_integration.hpp"
#include "cosmosim/io/restart_checkpoint.hpp"
#include "restart_equivalence_harness.hpp"

namespace cosmosim::tests {

inline core::SimulationState makeStage8ParticleState(
    std::size_t particle_count,
    core::ParticleSpecies species,
    const std::string& run_name) {
  core::SimulationState state;
  state.resizeParticles(particle_count);
  state.metadata.run_name = run_name;
  state.metadata.normalized_config_hash_hex = run_name;
  state.metadata.snapshot_stem = "snapshot";
  state.metadata.restart_stem = "restart";
  state.species.count_by_species[static_cast<std::size_t>(species)] = particle_count;
  for (std::size_t pidx = 0; pidx < state.particles.size(); ++pidx) {
    const double weight = static_cast<double>(pidx + 1U);
    state.particles.position_x_comoving[pidx] = 0.05 * weight;
    state.particles.position_y_comoving[pidx] = 0.07 * weight;
    state.particles.position_z_comoving[pidx] = 0.11 * weight;
    state.particles.velocity_x_peculiar[pidx] = 0.003 * weight;
    state.particles.velocity_y_peculiar[pidx] = -0.002 * weight;
    state.particles.velocity_z_peculiar[pidx] = 0.001 * weight;
    state.particles.mass_code[pidx] = 1.0 + 0.25 * weight;
    state.particle_sidecar.particle_id[pidx] = 900000U + pidx;
    state.particle_sidecar.sfc_key[pidx] = 7000U + pidx;
    state.particle_sidecar.species_tag[pidx] = static_cast<std::uint32_t>(species);
    state.particle_sidecar.owning_rank[pidx] = 0;
    if ((pidx % 2U) == 0U) {
      state.particle_sidecar.setGravitySofteningOverride(pidx, 0.0025 * weight);
    }
  }
  state.rebuildSpeciesIndex();
  return state;
}

inline core::SimulationState makeStage8DmState(std::size_t particle_count, const std::string& run_name) {
  return makeStage8ParticleState(particle_count, core::ParticleSpecies::kDarkMatter, run_name);
}

inline core::SimulationState makeStage8HydroToyState(std::size_t cell_count, const std::string& run_name) {
  core::SimulationState state = makeStage8ParticleState(cell_count, core::ParticleSpecies::kGas, run_name);
  state.resizeCells(cell_count);
  state.resizePatches(1);
  state.rebuildSpeciesIndex();
  state.refreshGasCellIdentityFromParticleOrder();
  for (std::size_t cidx = 0; cidx < state.cells.size(); ++cidx) {
    const double weight = static_cast<double>(cidx + 1U);
    state.cells.center_x_comoving[cidx] = state.particles.position_x_comoving[cidx];
    state.cells.center_y_comoving[cidx] = state.particles.position_y_comoving[cidx];
    state.cells.center_z_comoving[cidx] = state.particles.position_z_comoving[cidx];
    state.cells.mass_code[cidx] = state.particles.mass_code[cidx];
    state.cells.patch_index[cidx] = 0;
    state.gas_cells.density_code[cidx] = 2.0 + 0.1 * weight;
    state.gas_cells.pressure_code[cidx] = 3.0 + 0.2 * weight;
    state.gas_cells.internal_energy_code[cidx] = 4.0 + 0.3 * weight;
    state.gas_cells.temperature_code[cidx] = 5.0 + 0.4 * weight;
    state.gas_cells.sound_speed_code[cidx] = 6.0 + 0.5 * weight;
  }
  state.patches.patch_id[0] = 100;
  state.patches.level[0] = 0;
  state.patches.first_cell[0] = 0;
  state.patches.cell_count[0] = static_cast<std::uint32_t>(cell_count);
  // Patch ownership is part of the GasCellIdentityMap restart contract.  The
  // patch descriptor is filled after the legacy particle-bound identity lanes,
  // so rebuild the identity records once PatchSoa carries its authoritative id.
  state.refreshGasCellIdentityMapFromSidecarLanes();
  state.rebuildSpeciesIndex();
  return state;
}

inline core::HierarchicalTimeBinScheduler makeStage8Scheduler(
    std::uint32_t element_count,
    std::uint8_t max_bin,
    std::uint64_t current_tick = 0) {
  core::HierarchicalTimeBinScheduler scheduler(max_bin);
  scheduler.reset(element_count, 0, current_tick);
  for (std::uint32_t pidx = 0; pidx < element_count; ++pidx) {
    scheduler.setElementBin(pidx, static_cast<std::uint8_t>(pidx % (static_cast<std::uint32_t>(max_bin) + 1U)), scheduler.currentTick());
  }
  return scheduler;
}

inline core::IntegratorState makeStage8IntegratorState(
    std::uint64_t pm_cadence_steps = 1,
    std::uint8_t max_bin = 2) {
  core::IntegratorState integrator_state;
  integrator_state.current_time_code = 0.125;
  integrator_state.current_scale_factor = 0.5;
  integrator_state.current_redshift = 1.0;
  integrator_state.time_si_per_code = 1.0;
  integrator_state.current_boundary_kind = core::StepBoundaryKind::kGlobalSynchronizationPoint;
  integrator_state.last_completed_boundary_kind = core::StepBoundaryKind::kCheckpointPoint;
  integrator_state.last_completed_restart_safe = true;
  integrator_state.time_bins.hierarchical_enabled = true;
  integrator_state.time_bins.max_bin = max_bin;
  integrator_state.pm_refresh_enabled = true;
  integrator_state.pm_long_range_field_valid = true;
  integrator_state.pm_sync_state.reset(pm_cadence_steps);
  return integrator_state;
}

inline io::OutputCadencePersistentState makeStage8OutputCadenceState(
    bool output_enabled,
    std::uint64_t interval_steps = 5) {
  io::OutputCadencePersistentState output_state;
  output_state.output_enabled = output_enabled;
  output_state.write_restarts = output_enabled;
  output_state.snapshot_due = false;
  output_state.checkpoint_due = false;
  output_state.snapshot_interval_steps = output_enabled ? interval_steps : 0U;
  output_state.next_snapshot_step_index = output_enabled ? interval_steps : 0U;
  output_state.snapshot_stem = "snapshot";
  output_state.restart_stem = "restart";
  return output_state;
}

inline io::StochasticPersistentState makeStage8StochasticState() {
  io::StochasticPersistentState stochastic_state;
  stochastic_state.modules.push_back(io::StochasticModulePersistentState{
      .module_name = "star_formation",
      .schema_version = 1,
      .rng_policy = "stateless_splitmix64(seed,step_index,cell_index,rank_local_seed_offset)",
      .random_seed = 123456789ull,
      .rank_local_seed_offset = 0,
      .last_committed_step_index = 0,
      .deterministic_from_serialized_inputs = true,
  });
  stochastic_state.modules.push_back(io::StochasticModulePersistentState{
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

inline RestartEquivalenceScenario makeStage8Scenario(
    core::SimulationState state,
    core::IntegratorState integrator_state,
    core::HierarchicalTimeBinScheduler scheduler,
    io::OutputCadencePersistentState output_state,
    std::filesystem::path restart_path,
    std::uint64_t total_steps,
    std::uint64_t restart_step) {
  RestartEquivalenceScenario scenario;
  scenario.initial_state = std::move(state);
  scenario.initial_integrator_state = integrator_state;
  scenario.initial_scheduler = std::move(scheduler);
  scenario.initial_output_cadence_state = std::move(output_state);
  scenario.initial_stochastic_state = makeStage8StochasticState();
  scenario.total_steps = total_steps;
  scenario.restart_step = restart_step;
  scenario.restart_path = std::move(restart_path);
  scenario.tolerances.position_abs = 1.0e-13;
  scenario.tolerances.velocity_abs = 1.0e-13;
  scenario.tolerances.scalar_abs = 1.0e-13;
  syncParticleTimeBinMirror(scenario.initial_state, scenario.initial_scheduler);
  return scenario;
}

inline std::filesystem::path stage8RestartPath(const std::string& stem) {
  return std::filesystem::temp_directory_path() / ("cosmosim_" + stem + ".hdf5");
}

}  // namespace cosmosim::tests
