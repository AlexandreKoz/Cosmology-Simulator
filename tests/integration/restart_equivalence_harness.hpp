#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <functional>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "cosmosim/core/provenance.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/core/time_integration.hpp"
#include "cosmosim/io/restart_checkpoint.hpp"

namespace cosmosim::tests {

struct RestartEquivalenceTolerances {
  double position_abs = 1.0e-14;
  double velocity_abs = 1.0e-14;
  double scalar_abs = 1.0e-14;
};

struct RestartEquivalenceOutputEventLog {
  std::vector<std::uint64_t> snapshot_steps;
  std::vector<std::uint64_t> checkpoint_steps;
};

struct RestartEquivalenceStepContext {
  core::SimulationState& state;
  core::IntegratorState& integrator_state;
  core::HierarchicalTimeBinScheduler& scheduler;
  io::OutputCadencePersistentState& output_state;
  io::StochasticPersistentState& stochastic_state;
  std::span<const std::uint32_t> active_particle_indices;
  double dt_code = 0.0;
};

using RestartEquivalenceStepKernel = std::function<void(RestartEquivalenceStepContext&)>;

struct RestartEquivalenceScenario {
  core::SimulationState initial_state;
  core::IntegratorState initial_integrator_state;
  core::HierarchicalTimeBinScheduler initial_scheduler{0};
  io::OutputCadencePersistentState initial_output_cadence_state;
  io::StochasticPersistentState initial_stochastic_state;
  std::uint64_t total_steps = 10;
  std::uint64_t restart_step = 4;
  std::filesystem::path restart_path;
  RestartEquivalenceTolerances tolerances;
  RestartEquivalenceStepKernel step_kernel;
};

struct RestartEquivalenceResult {
  core::SimulationState direct_state;
  core::SimulationState restarted_state;
  core::IntegratorState direct_integrator_state;
  core::IntegratorState restarted_integrator_state;
  core::TimeBinPersistentState direct_scheduler_state;
  core::TimeBinPersistentState restarted_scheduler_state;
  io::OutputCadencePersistentState direct_output_cadence_state;
  io::OutputCadencePersistentState restarted_output_cadence_state;
  io::StochasticPersistentState direct_stochastic_state;
  io::StochasticPersistentState restarted_stochastic_state;
  RestartEquivalenceOutputEventLog direct_output_events;
  RestartEquivalenceOutputEventLog restarted_output_events;
};

inline void failRestartEquivalence(const std::string& message) {
  throw std::runtime_error("restart equivalence mismatch: " + message);
}

inline void requireNear(double lhs, double rhs, double tolerance, const std::string& label) {
  if (std::abs(lhs - rhs) > tolerance) {
    std::ostringstream os;
    os << label << " lhs=" << lhs << " rhs=" << rhs << " tolerance=" << tolerance;
    failRestartEquivalence(os.str());
  }
}

template <class TValue>
inline void requireVectorEqual(const std::vector<TValue>& lhs, const std::vector<TValue>& rhs, const std::string& label) {
  if (lhs != rhs) {
    failRestartEquivalence(label);
  }
}

template <class TAlignedVector>
inline void requireAlignedExact(const TAlignedVector& lhs, const TAlignedVector& rhs, const std::string& label) {
  if (lhs.size() != rhs.size()) {
    failRestartEquivalence(label + " size");
  }
  for (std::size_t i = 0; i < lhs.size(); ++i) {
    if (lhs[i] != rhs[i]) {
      failRestartEquivalence(label + " element " + std::to_string(i));
    }
  }
}

template <class TAlignedVector>
inline void requireAlignedNear(
    const TAlignedVector& lhs,
    const TAlignedVector& rhs,
    double tolerance,
    const std::string& label) {
  if (lhs.size() != rhs.size()) {
    failRestartEquivalence(label + " size");
  }
  for (std::size_t i = 0; i < lhs.size(); ++i) {
    requireNear(lhs[i], rhs[i], tolerance, label + "[" + std::to_string(i) + "]");
  }
}

inline void syncParticleTimeBinMirror(
    core::SimulationState& state,
    const core::HierarchicalTimeBinScheduler& scheduler) {
  const core::TimeBinPersistentState scheduler_state = scheduler.exportPersistentState();
  if (scheduler_state.bin_index.size() != state.particles.size()) {
    throw std::runtime_error("restart equivalence harness requires one scheduler element per particle");
  }
  state.particles.time_bin.assign(scheduler_state.bin_index.begin(), scheduler_state.bin_index.end());
  if (state.cells.size() == scheduler_state.bin_index.size()) {
    state.cells.time_bin.assign(scheduler_state.bin_index.begin(), scheduler_state.bin_index.end());
  }
}

inline io::StochasticPersistentState committedStochasticStateForStep(
    io::StochasticPersistentState stochastic_state,
    std::uint64_t step_index) {
  for (io::StochasticModulePersistentState& module_state : stochastic_state.modules) {
    module_state.last_committed_step_index = step_index;
  }
  return stochastic_state;
}

inline void commitOutputCadenceForStep(
    io::OutputCadencePersistentState& output_state,
    std::uint64_t step_index) {
  output_state.last_completed_step_index = step_index;
  if (!output_state.output_enabled || output_state.snapshot_interval_steps == 0U) {
    output_state.snapshot_due = false;
    output_state.checkpoint_due = false;
    return;
  }
  if (step_index >= output_state.next_snapshot_step_index) {
    output_state.snapshot_due = true;
    output_state.checkpoint_due = output_state.write_restarts;
    while (step_index >= output_state.next_snapshot_step_index) {
      output_state.next_snapshot_step_index += output_state.snapshot_interval_steps;
    }
  } else {
    output_state.snapshot_due = false;
    output_state.checkpoint_due = false;
  }
}

inline void applyDefaultToyParticleStep(RestartEquivalenceStepContext& context) {
  for (const std::uint32_t particle_index : context.active_particle_indices) {
    const double weight = static_cast<double>(particle_index + 1U);
    context.state.particles.velocity_x_peculiar[particle_index] += 0.001 * weight;
    context.state.particles.velocity_y_peculiar[particle_index] -= 0.0005 * weight;
    context.state.particles.velocity_z_peculiar[particle_index] += 0.00025 * weight;
    context.state.particles.position_x_comoving[particle_index] +=
        context.dt_code * context.state.particles.velocity_x_peculiar[particle_index];
    context.state.particles.position_y_comoving[particle_index] +=
        context.dt_code * context.state.particles.velocity_y_peculiar[particle_index];
    context.state.particles.position_z_comoving[particle_index] +=
        context.dt_code * context.state.particles.velocity_z_peculiar[particle_index];
  }
}

inline void recordOutputEventsForStep(
    const io::OutputCadencePersistentState& output_state,
    std::uint64_t step_index,
    RestartEquivalenceOutputEventLog* output_events) {
  if (output_events == nullptr) {
    return;
  }
  if (output_state.snapshot_due) {
    output_events->snapshot_steps.push_back(step_index);
  }
  if (output_state.checkpoint_due) {
    output_events->checkpoint_steps.push_back(step_index);
  }
}

inline void advanceRuntimeStep(
    core::SimulationState& state,
    core::IntegratorState& integrator_state,
    core::HierarchicalTimeBinScheduler& scheduler,
    io::OutputCadencePersistentState& output_state,
    io::StochasticPersistentState& stochastic_state,
    const RestartEquivalenceStepKernel& step_kernel,
    RestartEquivalenceOutputEventLog* output_events) {
  auto active = scheduler.beginSubstep();
  const double dt_code = 0.03125;
  RestartEquivalenceStepContext context{
      .state = state,
      .integrator_state = integrator_state,
      .scheduler = scheduler,
      .output_state = output_state,
      .stochastic_state = stochastic_state,
      .active_particle_indices = active,
      .dt_code = dt_code,
  };
  if (step_kernel) {
    step_kernel(context);
  } else {
    applyDefaultToyParticleStep(context);
  }
  scheduler.endSubstep();
  ++integrator_state.step_index;
  integrator_state.current_time_code += dt_code;
  integrator_state.dt_time_code = dt_code;
  integrator_state.current_scale_factor += 0.0001;
  integrator_state.current_redshift = (1.0 / integrator_state.current_scale_factor) - 1.0;
  integrator_state.current_boundary_kind = core::StepBoundaryKind::kGlobalSynchronizationPoint;
  integrator_state.last_completed_boundary_kind = core::StepBoundaryKind::kCheckpointPoint;
  integrator_state.inside_kdk_step = false;
  integrator_state.last_completed_restart_safe = true;
  integrator_state.time_bins.hierarchical_enabled = true;
  integrator_state.time_bins.max_bin = scheduler.maxBin();
  integrator_state.time_bins.active_bin = 0;
  integrator_state.pm_refresh_enabled = true;
  integrator_state.pm_long_range_field_valid = true;
  const auto pm_event = integrator_state.pm_sync_state.registerKickOpportunity(
      integrator_state.step_index,
      integrator_state.current_scale_factor,
      true);
  if (pm_event.refresh_long_range_field) {
    integrator_state.pm_sync_state.commitRefresh(pm_event);
  }
  state.metadata.step_index = integrator_state.step_index;
  state.metadata.scale_factor = integrator_state.current_scale_factor;
  commitOutputCadenceForStep(output_state, integrator_state.step_index);
  recordOutputEventsForStep(output_state, integrator_state.step_index, output_events);
  stochastic_state = committedStochasticStateForStep(stochastic_state, integrator_state.step_index);
  syncParticleTimeBinMirror(state, scheduler);
}

inline io::RestartWritePayload makeRestartEquivalencePayload(
    const core::SimulationState& state,
    const core::IntegratorState& integrator_state,
    const core::HierarchicalTimeBinScheduler& scheduler,
    const io::OutputCadencePersistentState& output_state,
    const io::StochasticPersistentState& stochastic_state) {
  io::RestartWritePayload payload;
  payload.persistent_state.simulation_state = &state;
  payload.integrator_state = &integrator_state;
  payload.scheduler = &scheduler;
  payload.normalized_config_text = "stage8_restart_equivalence_harness=1\n";
  payload.normalized_config_hash_hex = core::stableConfigHashHex(payload.normalized_config_text);
  payload.provenance = core::makeProvenanceRecord(payload.normalized_config_hash_hex, "restart_equivalence", 0U);
  payload.provenance.normalized_config = payload.normalized_config_text;
  payload.provenance.gravity_treepm_update_cadence_steps = 1U;
  payload.distributed_gravity_state.decomposition_epoch = integrator_state.step_index;
  payload.distributed_gravity_state.world_size = 1;
  payload.distributed_gravity_state.pm_grid_nx = 8;
  payload.distributed_gravity_state.pm_grid_ny = 8;
  payload.distributed_gravity_state.pm_grid_nz = 8;
  payload.distributed_gravity_state.gravity_kick_opportunity =
      integrator_state.pm_sync_state.exportPersistentState().gravity_kick_opportunity;
  payload.distributed_gravity_state.pm_update_cadence_steps =
      integrator_state.pm_sync_state.exportPersistentState().cadence_steps;
  payload.distributed_gravity_state.long_range_field_version =
      integrator_state.pm_sync_state.exportPersistentState().field_version;
  payload.distributed_gravity_state.last_long_range_refresh_opportunity =
      integrator_state.pm_sync_state.exportPersistentState().last_refresh_opportunity;
  payload.distributed_gravity_state.long_range_field_built_step_index = integrator_state.step_index;
  payload.distributed_gravity_state.long_range_field_built_scale_factor = integrator_state.current_scale_factor;
  payload.distributed_gravity_state.owning_rank_by_item.assign(state.particles.size(), 0);
  payload.distributed_gravity_state.pm_slab_begin_x_by_rank = {0};
  payload.distributed_gravity_state.pm_slab_end_x_by_rank = {8};
  payload.output_cadence_state = output_state;
  payload.stochastic_state = stochastic_state;
  return payload;
}

inline void runSteps(
    core::SimulationState& state,
    core::IntegratorState& integrator_state,
    core::HierarchicalTimeBinScheduler& scheduler,
    io::OutputCadencePersistentState& output_state,
    io::StochasticPersistentState& stochastic_state,
    std::uint64_t step_count,
    const RestartEquivalenceStepKernel& step_kernel = {},
    RestartEquivalenceOutputEventLog* output_events = nullptr) {
  for (std::uint64_t step = 0; step < step_count; ++step) {
    advanceRuntimeStep(
        state, integrator_state, scheduler, output_state, stochastic_state, step_kernel, output_events);
  }
}

inline std::vector<io::StochasticModulePersistentState> sortedStochasticModules(
    io::StochasticPersistentState stochastic_state) {
  std::sort(
      stochastic_state.modules.begin(),
      stochastic_state.modules.end(),
      [](const io::StochasticModulePersistentState& lhs, const io::StochasticModulePersistentState& rhs) {
        return lhs.module_name < rhs.module_name;
      });
  return stochastic_state.modules;
}

inline void compareStochasticState(
    const io::StochasticPersistentState& lhs,
    const io::StochasticPersistentState& rhs) {
  const auto lhs_modules = sortedStochasticModules(lhs);
  const auto rhs_modules = sortedStochasticModules(rhs);
  if (lhs_modules.size() != rhs_modules.size()) {
    failRestartEquivalence("stochastic module count");
  }
  for (std::size_t i = 0; i < lhs_modules.size(); ++i) {
    if (lhs_modules[i].module_name != rhs_modules[i].module_name ||
        lhs_modules[i].schema_version != rhs_modules[i].schema_version ||
        lhs_modules[i].rng_policy != rhs_modules[i].rng_policy ||
        lhs_modules[i].random_seed != rhs_modules[i].random_seed ||
        lhs_modules[i].rank_local_seed_offset != rhs_modules[i].rank_local_seed_offset ||
        lhs_modules[i].last_committed_step_index != rhs_modules[i].last_committed_step_index ||
        lhs_modules[i].deterministic_from_serialized_inputs != rhs_modules[i].deterministic_from_serialized_inputs) {
      failRestartEquivalence("stochastic module state " + std::to_string(i));
    }
  }
}

inline void compareSchedulerState(
    const core::TimeBinPersistentState& lhs,
    const core::TimeBinPersistentState& rhs) {
  if (lhs.current_tick != rhs.current_tick || lhs.max_bin != rhs.max_bin) {
    failRestartEquivalence("scheduler scalar state");
  }
  requireVectorEqual(lhs.bin_index, rhs.bin_index, "scheduler bin_index");
  requireVectorEqual(lhs.next_activation_tick, rhs.next_activation_tick, "scheduler next_activation_tick");
  requireVectorEqual(lhs.active_flag, rhs.active_flag, "scheduler active_flag");
  requireVectorEqual(lhs.pending_bin_index, rhs.pending_bin_index, "scheduler pending_bin_index");
}

inline void compareOutputCadenceState(
    const io::OutputCadencePersistentState& lhs,
    const io::OutputCadencePersistentState& rhs) {
  if (lhs.output_enabled != rhs.output_enabled || lhs.write_restarts != rhs.write_restarts ||
      lhs.snapshot_due != rhs.snapshot_due || lhs.checkpoint_due != rhs.checkpoint_due ||
      lhs.last_completed_step_index != rhs.last_completed_step_index ||
      lhs.snapshot_interval_steps != rhs.snapshot_interval_steps ||
      lhs.next_snapshot_step_index != rhs.next_snapshot_step_index ||
      lhs.snapshot_stem != rhs.snapshot_stem || lhs.restart_stem != rhs.restart_stem) {
    failRestartEquivalence("output cadence state");
  }
}

inline void compareIntegratorState(
    const core::IntegratorState& lhs,
    const core::IntegratorState& rhs,
    const RestartEquivalenceTolerances& tolerances) {
  requireNear(lhs.current_time_code, rhs.current_time_code, tolerances.scalar_abs, "integrator current_time_code");
  requireNear(lhs.current_scale_factor, rhs.current_scale_factor, tolerances.scalar_abs, "integrator current_scale_factor");
  requireNear(lhs.current_redshift, rhs.current_redshift, tolerances.scalar_abs, "integrator current_redshift");
  if (lhs.step_index != rhs.step_index || lhs.current_boundary_kind != rhs.current_boundary_kind ||
      lhs.last_completed_boundary_kind != rhs.last_completed_boundary_kind || lhs.inside_kdk_step != rhs.inside_kdk_step ||
      lhs.last_completed_restart_safe != rhs.last_completed_restart_safe ||
      lhs.time_bins.hierarchical_enabled != rhs.time_bins.hierarchical_enabled ||
      lhs.time_bins.max_bin != rhs.time_bins.max_bin || lhs.pm_refresh_enabled != rhs.pm_refresh_enabled ||
      lhs.pm_long_range_field_valid != rhs.pm_long_range_field_valid) {
    failRestartEquivalence("integrator phase/scheduler metadata");
  }
  const auto lhs_pm = lhs.pm_sync_state.exportPersistentState();
  const auto rhs_pm = rhs.pm_sync_state.exportPersistentState();
  if (lhs_pm.cadence_steps != rhs_pm.cadence_steps ||
      lhs_pm.gravity_kick_opportunity != rhs_pm.gravity_kick_opportunity ||
      lhs_pm.last_refresh_opportunity != rhs_pm.last_refresh_opportunity ||
      lhs_pm.field_version != rhs_pm.field_version ||
      lhs_pm.last_refresh_step_index != rhs_pm.last_refresh_step_index ||
      lhs_pm.refresh_commit_pending != rhs_pm.refresh_commit_pending ||
      lhs_pm.pending_refresh_opportunity != rhs_pm.pending_refresh_opportunity ||
      lhs_pm.pending_refresh_field_version != rhs_pm.pending_refresh_field_version) {
    failRestartEquivalence("PM sync persistent state");
  }
  requireNear(lhs_pm.last_refresh_scale_factor, rhs_pm.last_refresh_scale_factor, tolerances.scalar_abs, "PM refresh scale factor");
}

inline void compareSimulationState(
    const core::SimulationState& lhs,
    const core::SimulationState& rhs,
    const RestartEquivalenceTolerances& tolerances) {
  if (lhs.particles.size() != rhs.particles.size()) {
    failRestartEquivalence("particle count");
  }
  requireAlignedNear(lhs.particles.position_x_comoving, rhs.particles.position_x_comoving, tolerances.position_abs, "pos_x");
  requireAlignedNear(lhs.particles.position_y_comoving, rhs.particles.position_y_comoving, tolerances.position_abs, "pos_y");
  requireAlignedNear(lhs.particles.position_z_comoving, rhs.particles.position_z_comoving, tolerances.position_abs, "pos_z");
  requireAlignedNear(lhs.particles.velocity_x_peculiar, rhs.particles.velocity_x_peculiar, tolerances.velocity_abs, "vel_x");
  requireAlignedNear(lhs.particles.velocity_y_peculiar, rhs.particles.velocity_y_peculiar, tolerances.velocity_abs, "vel_y");
  requireAlignedNear(lhs.particles.velocity_z_peculiar, rhs.particles.velocity_z_peculiar, tolerances.velocity_abs, "vel_z");
  requireAlignedNear(lhs.particles.mass_code, rhs.particles.mass_code, tolerances.scalar_abs, "mass");
  requireAlignedExact(lhs.particles.time_bin, rhs.particles.time_bin, "time_bin mirror");
  requireAlignedExact(lhs.particle_sidecar.particle_id, rhs.particle_sidecar.particle_id, "particle_id");
  requireAlignedExact(lhs.particle_sidecar.species_tag, rhs.particle_sidecar.species_tag, "species_tag");
  requireAlignedExact(lhs.particle_sidecar.owning_rank, rhs.particle_sidecar.owning_rank, "owning_rank");
  requireAlignedNear(lhs.particle_sidecar.gravity_softening_comoving, rhs.particle_sidecar.gravity_softening_comoving, tolerances.scalar_abs, "gravity_softening");
  requireAlignedExact(lhs.particle_sidecar.has_gravity_softening_override, rhs.particle_sidecar.has_gravity_softening_override, "has_gravity_softening_override");
  if (lhs.cells.size() != rhs.cells.size() || lhs.gas_cells.size() != rhs.gas_cells.size()) {
    failRestartEquivalence("gas cell count");
  }
  requireAlignedNear(lhs.cells.center_x_comoving, rhs.cells.center_x_comoving, tolerances.position_abs, "cell_center_x");
  requireAlignedNear(lhs.cells.center_y_comoving, rhs.cells.center_y_comoving, tolerances.position_abs, "cell_center_y");
  requireAlignedNear(lhs.cells.center_z_comoving, rhs.cells.center_z_comoving, tolerances.position_abs, "cell_center_z");
  requireAlignedNear(lhs.cells.mass_code, rhs.cells.mass_code, tolerances.scalar_abs, "cell_mass");
  requireAlignedExact(lhs.cells.time_bin, rhs.cells.time_bin, "cell_time_bin mirror");
  requireAlignedExact(lhs.cells.patch_index, rhs.cells.patch_index, "cell_patch_index");
  requireAlignedExact(lhs.patches.patch_id, rhs.patches.patch_id, "patch_id");
  requireAlignedExact(lhs.patches.level, rhs.patches.level, "patch_level");
  requireAlignedExact(lhs.patches.first_cell, rhs.patches.first_cell, "patch_first_cell");
  requireAlignedExact(lhs.patches.cell_count, rhs.patches.cell_count, "patch_cell_count");
  requireAlignedExact(lhs.patches.owning_rank, rhs.patches.owning_rank, "patch_owning_rank");
  requireAlignedExact(lhs.gas_cells.gas_cell_id, rhs.gas_cells.gas_cell_id, "gas_cell_id");
  requireAlignedExact(lhs.gas_cells.parent_particle_id, rhs.gas_cells.parent_particle_id, "gas_parent_particle_id");
  requireAlignedNear(lhs.gas_cells.density_code, rhs.gas_cells.density_code, tolerances.scalar_abs, "gas_density");
  requireAlignedNear(lhs.gas_cells.pressure_code, rhs.gas_cells.pressure_code, tolerances.scalar_abs, "gas_pressure");
  requireAlignedNear(lhs.gas_cells.internal_energy_code, rhs.gas_cells.internal_energy_code, tolerances.scalar_abs, "gas_internal_energy");
  requireAlignedNear(lhs.gas_cells.temperature_code, rhs.gas_cells.temperature_code, tolerances.scalar_abs, "gas_temperature");
  requireAlignedNear(lhs.gas_cells.sound_speed_code, rhs.gas_cells.sound_speed_code, tolerances.scalar_abs, "gas_sound_speed");
  if (lhs.star_particles.size() != rhs.star_particles.size()) {
    failRestartEquivalence("star particle sidecar count");
  }
  requireAlignedExact(lhs.star_particles.particle_index, rhs.star_particles.particle_index, "star_particle_index");
  requireAlignedNear(lhs.star_particles.formation_scale_factor, rhs.star_particles.formation_scale_factor, tolerances.scalar_abs, "star_formation_scale_factor");
  requireAlignedNear(lhs.star_particles.birth_mass_code, rhs.star_particles.birth_mass_code, tolerances.scalar_abs, "star_birth_mass");
  requireAlignedNear(lhs.star_particles.metallicity_mass_fraction, rhs.star_particles.metallicity_mass_fraction, tolerances.scalar_abs, "star_metallicity");
  requireAlignedNear(lhs.star_particles.stellar_age_years_last, rhs.star_particles.stellar_age_years_last, tolerances.scalar_abs, "star_age_last");
  requireAlignedNear(lhs.star_particles.stellar_returned_mass_cumulative_code, rhs.star_particles.stellar_returned_mass_cumulative_code, tolerances.scalar_abs, "star_returned_mass");
  requireAlignedNear(lhs.star_particles.stellar_returned_metals_cumulative_code, rhs.star_particles.stellar_returned_metals_cumulative_code, tolerances.scalar_abs, "star_returned_metals");
  requireAlignedNear(lhs.star_particles.stellar_feedback_energy_cumulative_erg, rhs.star_particles.stellar_feedback_energy_cumulative_erg, tolerances.scalar_abs, "star_feedback_energy");
  for (std::size_t channel = 0; channel < lhs.star_particles.stellar_returned_mass_channel_cumulative_code.size(); ++channel) {
    requireAlignedNear(lhs.star_particles.stellar_returned_mass_channel_cumulative_code[channel], rhs.star_particles.stellar_returned_mass_channel_cumulative_code[channel], tolerances.scalar_abs, "star_returned_mass_channel" + std::to_string(channel));
    requireAlignedNear(lhs.star_particles.stellar_returned_metals_channel_cumulative_code[channel], rhs.star_particles.stellar_returned_metals_channel_cumulative_code[channel], tolerances.scalar_abs, "star_returned_metals_channel" + std::to_string(channel));
    requireAlignedNear(lhs.star_particles.stellar_feedback_energy_channel_cumulative_erg[channel], rhs.star_particles.stellar_feedback_energy_channel_cumulative_erg[channel], tolerances.scalar_abs, "star_feedback_energy_channel" + std::to_string(channel));
  }
  if (lhs.metadata.step_index != rhs.metadata.step_index) {
    failRestartEquivalence("state metadata step_index");
  }
  requireNear(lhs.metadata.scale_factor, rhs.metadata.scale_factor, tolerances.scalar_abs, "state metadata scale_factor");
}

inline void compareOutputEventLog(
    const RestartEquivalenceOutputEventLog& lhs,
    const RestartEquivalenceOutputEventLog& rhs) {
  requireVectorEqual(lhs.snapshot_steps, rhs.snapshot_steps, "output snapshot event sequence");
  requireVectorEqual(lhs.checkpoint_steps, rhs.checkpoint_steps, "output checkpoint event sequence");
}

inline RestartEquivalenceResult runRestartEquivalenceScenario(RestartEquivalenceScenario scenario) {
  if (scenario.restart_step > scenario.total_steps) {
    throw std::invalid_argument("restart_step must be <= total_steps");
  }
  if (scenario.restart_path.empty()) {
    throw std::invalid_argument("restart_path must be set");
  }

  RestartEquivalenceResult result;
  result.direct_state = scenario.initial_state;
  result.direct_integrator_state = scenario.initial_integrator_state;
  core::HierarchicalTimeBinScheduler direct_scheduler = scenario.initial_scheduler;
  result.direct_output_cadence_state = scenario.initial_output_cadence_state;
  result.direct_stochastic_state = committedStochasticStateForStep(
      scenario.initial_stochastic_state,
      result.direct_integrator_state.step_index);
  syncParticleTimeBinMirror(result.direct_state, direct_scheduler);
  runSteps(
      result.direct_state,
      result.direct_integrator_state,
      direct_scheduler,
      result.direct_output_cadence_state,
      result.direct_stochastic_state,
      scenario.total_steps,
      scenario.step_kernel,
      &result.direct_output_events);
  result.direct_scheduler_state = direct_scheduler.exportPersistentState();

  result.restarted_state = scenario.initial_state;
  result.restarted_integrator_state = scenario.initial_integrator_state;
  core::HierarchicalTimeBinScheduler restarted_scheduler = scenario.initial_scheduler;
  result.restarted_output_cadence_state = scenario.initial_output_cadence_state;
  result.restarted_stochastic_state = committedStochasticStateForStep(
      scenario.initial_stochastic_state,
      result.restarted_integrator_state.step_index);
  syncParticleTimeBinMirror(result.restarted_state, restarted_scheduler);
  runSteps(
      result.restarted_state,
      result.restarted_integrator_state,
      restarted_scheduler,
      result.restarted_output_cadence_state,
      result.restarted_stochastic_state,
      scenario.restart_step,
      scenario.step_kernel,
      &result.restarted_output_events);

  const io::RestartWritePayload payload = makeRestartEquivalencePayload(
      result.restarted_state,
      result.restarted_integrator_state,
      restarted_scheduler,
      result.restarted_output_cadence_state,
      result.restarted_stochastic_state);
  io::writeRestartCheckpointHdf5(scenario.restart_path, payload);
  const io::RestartReadResult readback = io::readRestartCheckpointHdf5(scenario.restart_path);
  result.restarted_state = readback.state;
  result.restarted_integrator_state = readback.integrator_state;
  result.restarted_output_cadence_state = readback.output_cadence_state;
  result.restarted_stochastic_state = readback.stochastic_state;
  restarted_scheduler = core::HierarchicalTimeBinScheduler(readback.scheduler_state.max_bin);
  restarted_scheduler.importPersistentState(readback.scheduler_state);
  runSteps(
      result.restarted_state,
      result.restarted_integrator_state,
      restarted_scheduler,
      result.restarted_output_cadence_state,
      result.restarted_stochastic_state,
      scenario.total_steps - scenario.restart_step,
      scenario.step_kernel,
      &result.restarted_output_events);
  result.restarted_scheduler_state = restarted_scheduler.exportPersistentState();

  compareSimulationState(result.direct_state, result.restarted_state, scenario.tolerances);
  compareIntegratorState(result.direct_integrator_state, result.restarted_integrator_state, scenario.tolerances);
  compareSchedulerState(result.direct_scheduler_state, result.restarted_scheduler_state);
  compareOutputCadenceState(result.direct_output_cadence_state, result.restarted_output_cadence_state);
  compareOutputEventLog(result.direct_output_events, result.restarted_output_events);
  compareStochasticState(result.direct_stochastic_state, result.restarted_stochastic_state);
  return result;
}

}  // namespace cosmosim::tests
