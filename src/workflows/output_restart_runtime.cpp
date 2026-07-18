#include "cosmosim/workflows/output_restart_runtime.hpp"
#include "workflows/internal/runtime_stage_resource_access.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

#include "cosmosim/core/config.hpp"
#include "cosmosim/core/cuda_runtime.hpp"
#include "cosmosim/core/profiling.hpp"
#include "cosmosim/core/provenance.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/gravity/pm_solver.hpp"
#include "cosmosim/gravity/tree_pm_coupling.hpp"
#include "cosmosim/io/restart_checkpoint.hpp"
#include "cosmosim/io/snapshot_hdf5.hpp"
#include "cosmosim/parallel/distributed_memory.hpp"
#include "workflows/internal/output_verification.hpp"

namespace cosmosim::workflows::internal {
namespace {

[[nodiscard]] std::string formatRuntimeDouble(double value) {
  std::ostringstream stream;
  stream << std::scientific
         << std::setprecision(std::numeric_limits<double>::max_digits10)
         << value;
  return stream.str();
}

[[nodiscard]] std::string treePmAssignmentSchemeName(
    core::TreePmAssignmentScheme assignment_scheme) {
  switch (assignment_scheme) {
    case core::TreePmAssignmentScheme::kCic: return "cic";
    case core::TreePmAssignmentScheme::kTsc: return "tsc";
  }
  throw std::runtime_error("unhandled TreePM assignment scheme enum value");
}

[[nodiscard]] std::string treePmOpeningCriterionName(
    core::TreePmOpeningCriterion opening_criterion) {
  switch (opening_criterion) {
    case core::TreePmOpeningCriterion::kGeometric: return "geometric";
    case core::TreePmOpeningCriterion::kComDistance: return "com_distance";
    case core::TreePmOpeningCriterion::kRelativeForceError:
      return "relative_force_error";
  }
  throw std::runtime_error("unhandled TreePM opening criterion enum value");
}

[[nodiscard]] std::string pmDecompositionModeName(
    core::PmDecompositionMode mode) {
  switch (mode) {
    case core::PmDecompositionMode::kSlab: return "slab";
    case core::PmDecompositionMode::kPencil: return "pencil";
  }
  throw std::runtime_error("unhandled PM decomposition mode enum value");
}

[[nodiscard]] bool hasSpeciesSpecificSoftening(
    const core::SimulationConfig& config) {
  return config.numerics.gravity_softening_gas_kpc_comoving > 0.0 ||
      config.numerics.gravity_softening_dark_matter_kpc_comoving > 0.0 ||
      config.numerics.gravity_softening_star_kpc_comoving > 0.0 ||
      config.numerics.gravity_softening_black_hole_kpc_comoving > 0.0 ||
      config.numerics.gravity_softening_tracer_kpc_comoving > 0.0;
}

[[nodiscard]] std::string describeSofteningPolicy(
    const core::SimulationConfig& config) {
  return hasSpeciesSpecificSoftening(config)
      ? "comoving_species"
      : "comoving_fixed";
}

[[nodiscard]] io::StochasticPersistentState buildStochasticPersistentState(
    const core::SimulationConfig& config,
    const core::IntegratorState& integrator_state,
    std::uint32_t rank_local_seed_offset) {
  io::StochasticPersistentState stochastic_state;
  if (config.physics.enable_star_formation && config.physics.sf_stochastic_spawning) {
    stochastic_state.modules.push_back(io::StochasticModulePersistentState{
        .module_name = "star_formation",
        .schema_version = 1,
        .rng_policy = "stateless_splitmix64(seed,step_index,cell_index,rank_local_seed_offset)",
        .random_seed = config.physics.sf_random_seed,
        .rank_local_seed_offset = rank_local_seed_offset,
        .last_committed_step_index = integrator_state.step_index,
        .deterministic_from_serialized_inputs = true,
    });
  }
  if (config.physics.enable_feedback && config.physics.fb_variant == core::FeedbackVariant::kStochastic) {
    stochastic_state.modules.push_back(io::StochasticModulePersistentState{
        .module_name = "stellar_feedback",
        .schema_version = 1,
        .rng_policy = "stateless_splitmix64(seed,step_index,star_index)",
        .random_seed = config.physics.fb_random_seed,
        .rank_local_seed_offset = rank_local_seed_offset,
        .last_committed_step_index = integrator_state.step_index,
        .deterministic_from_serialized_inputs = true,
    });
  }
  return stochastic_state;
}

[[nodiscard]] bool stochasticStatesEquivalent(
    io::StochasticPersistentState lhs,
    io::StochasticPersistentState rhs) {
  const auto less_by_name = [](
      const io::StochasticModulePersistentState& a,
      const io::StochasticModulePersistentState& b) {
    return a.module_name < b.module_name;
  };
  std::sort(lhs.modules.begin(), lhs.modules.end(), less_by_name);
  std::sort(rhs.modules.begin(), rhs.modules.end(), less_by_name);
  if (lhs.modules.size() != rhs.modules.size()) {
    return false;
  }
  for (std::size_t i = 0; i < lhs.modules.size(); ++i) {
    const auto& a = lhs.modules[i];
    const auto& b = rhs.modules[i];
    if (a.module_name != b.module_name || a.schema_version != b.schema_version ||
        a.rng_policy != b.rng_policy || a.random_seed != b.random_seed ||
        a.rank_local_seed_offset != b.rank_local_seed_offset ||
        a.last_committed_step_index != b.last_committed_step_index ||
        a.deterministic_from_serialized_inputs != b.deterministic_from_serialized_inputs) {
      return false;
    }
  }
  return true;
}

[[nodiscard]] bool gravityForceCachesEquivalent(
    const io::GravityForceCachePersistentState& lhs,
    const io::GravityForceCachePersistentState& rhs) {
  return lhs.valid == rhs.valid &&
      lhs.particle_id == rhs.particle_id &&
      lhs.gas_cell_id == rhs.gas_cell_id &&
      lhs.particle_accel_x_comoving == rhs.particle_accel_x_comoving &&
      lhs.particle_accel_y_comoving == rhs.particle_accel_y_comoving &&
      lhs.particle_accel_z_comoving == rhs.particle_accel_z_comoving &&
      lhs.cell_accel_x_comoving == rhs.cell_accel_x_comoving &&
      lhs.cell_accel_y_comoving == rhs.cell_accel_y_comoving &&
      lhs.cell_accel_z_comoving == rhs.cell_accel_z_comoving;
}

[[nodiscard]] std::string formatIndexedRankedFileStem(
    std::string_view stem,
    std::uint64_t index,
    int world_size,
    int world_rank) {
  std::ostringstream out;
  out << stem << '_' << std::setw(3) << std::setfill('0') << index;
  if (world_size > 1) {
    out << "_rank" << std::setw(3) << std::setfill('0') << world_rank;
  }
  out << ".hdf5";
  return out.str();
}


[[nodiscard]] core::ProvenanceRecord makeGravityAwareProvenanceRecord(
    const core::FrozenConfig& frozen_config,
    const core::SimulationConfig& config) {
  core::ProvenanceRecord record = core::makeProvenanceRecord(
      frozen_config.provenance.config_hash_hex, "unknown");
  record.config_schema_name = "cosmosim_config";
  record.config_schema_version = std::to_string(frozen_config.config.schema_version);
  record.raw_input_config = frozen_config.raw_text;
  record.normalized_config = frozen_config.normalized_text;
  record.derived_runtime_state =
      core::serializeDerivedRuntimeConfig(core::deriveRuntimeConfig(frozen_config));
  const double dx = config.cosmology.box_size_x_mpc_comoving /
      static_cast<double>(config.numerics.treepm_pm_grid_nx);
  const double dy = config.cosmology.box_size_y_mpc_comoving /
      static_cast<double>(config.numerics.treepm_pm_grid_ny);
  const double dz = config.cosmology.box_size_z_mpc_comoving /
      static_cast<double>(config.numerics.treepm_pm_grid_nz);
  const double mesh_spacing_mpc_comoving = std::cbrt(dx * dy * dz);
  const gravity::TreePmSplitPolicy split_policy = gravity::makeTreePmSplitPolicyFromMeshSpacing(
      config.numerics.treepm_asmth_cells,
      config.numerics.treepm_rcut_cells,
      mesh_spacing_mpc_comoving);
  record.gravity_treepm_pm_grid = config.numerics.treepm_pm_grid_nx;
  record.gravity_treepm_pm_grid_nx = config.numerics.treepm_pm_grid_nx;
  record.gravity_treepm_pm_grid_ny = config.numerics.treepm_pm_grid_ny;
  record.gravity_treepm_pm_grid_nz = config.numerics.treepm_pm_grid_nz;
  record.gravity_treepm_assignment_scheme =
      treePmAssignmentSchemeName(config.numerics.treepm_assignment_scheme);
  record.gravity_treepm_window_deconvolution =
      config.numerics.treepm_enable_window_deconvolution;
  record.gravity_treepm_asmth_cells = config.numerics.treepm_asmth_cells;
  record.gravity_treepm_rcut_cells = config.numerics.treepm_rcut_cells;
  record.gravity_treepm_tree_opening_criterion =
      treePmOpeningCriterionName(config.numerics.treepm_tree_opening_criterion);
  record.gravity_treepm_tree_opening_theta = config.numerics.treepm_tree_opening_theta;
  record.gravity_treepm_tree_relative_force_tolerance =
      config.numerics.treepm_tree_relative_force_tolerance;
  record.gravity_treepm_tree_relative_force_acceleration_floor =
      config.numerics.treepm_tree_relative_force_acceleration_floor;
  record.gravity_treepm_mesh_spacing_mpc_comoving = mesh_spacing_mpc_comoving;
  record.gravity_treepm_mesh_spacing_x_mpc_comoving = dx;
  record.gravity_treepm_mesh_spacing_y_mpc_comoving = dy;
  record.gravity_treepm_mesh_spacing_z_mpc_comoving = dz;
  record.gravity_treepm_split_scale_mpc_comoving = split_policy.split_scale_comoving;
  record.gravity_treepm_cutoff_radius_mpc_comoving = split_policy.cutoff_radius_comoving;
  record.gravity_treepm_update_cadence_steps = config.numerics.treepm_update_cadence_steps;
  record.gravity_treepm_pm_decomposition_mode =
      pmDecompositionModeName(config.numerics.treepm_pm_decomposition_mode);
  record.gravity_treepm_tree_exchange_batch_bytes =
      config.numerics.treepm_tree_exchange_batch_bytes;
  record.gravity_softening_policy = describeSofteningPolicy(config);
  record.gravity_softening_kernel = "plummer";
  record.gravity_softening_epsilon_kpc_comoving = config.numerics.gravity_softening_kpc_comoving;
  record.gravity_pm_fft_backend = gravity::PmSolver::fftBackendName();
  switch (config.mode.zoom_long_range_strategy) {
    case core::ZoomLongRangeStrategy::kDisabled:
      record.zoom_long_range_strategy = "disabled";
      break;
    case core::ZoomLongRangeStrategy::kGlobalCoarsePlusFocusedHighResCorrection:
      record.zoom_long_range_strategy = "global_coarse_plus_focused_highres_correction";
      break;
  }
  record.zoom_region_center_x_mpc_comoving = config.mode.zoom_region_center_x_mpc_comoving;
  record.zoom_region_center_y_mpc_comoving = config.mode.zoom_region_center_y_mpc_comoving;
  record.zoom_region_center_z_mpc_comoving = config.mode.zoom_region_center_z_mpc_comoving;
  record.zoom_region_radius_mpc_comoving = config.mode.zoom_region_radius_mpc_comoving;
  record.zoom_focused_pm_grid = std::to_string(config.mode.zoom_focused_pm_grid_nx) + "x" +
      std::to_string(config.mode.zoom_focused_pm_grid_ny) + "x" +
      std::to_string(config.mode.zoom_focused_pm_grid_nz);
  record.zoom_contamination_radius_mpc_comoving = config.mode.zoom_contamination_radius_mpc_comoving;
  return record;
}

[[nodiscard]] std::string pmSlabSignature(const parallel::DistributedRestartState& distributed_state) {
  std::ostringstream stream;
  for (std::size_t rank = 0; rank < distributed_state.pm_slab_begin_x_by_rank.size(); ++rank) {
    if (rank > 0) {
      stream << ';';
    }
    stream << rank << ':' << distributed_state.pm_slab_begin_x_by_rank[rank] << '-' <<
        distributed_state.pm_slab_end_x_by_rank[rank];
  }
  return stream.str();
}

bool maybeWriteOutputs(
    const core::FrozenConfig& frozen_config,
    const core::SimulationConfig& config,
    const core::SimulationState& state,
    const core::IntegratorState& integrator_state,
    const core::HierarchicalTimeBinScheduler& scheduler,
    const core::HierarchicalTimeBinScheduler& gas_cell_scheduler,
    const workflows::GravityRestartStateProvider& gravity_state,
    ReferenceWorkflowReport& report,
    core::ProfilerSession& profiler,
    bool write_outputs_enabled,
    bool snapshot_due,
    bool checkpoint_due,
    std::uint64_t snapshot_interval_steps,
    std::uint64_t next_snapshot_step_index,
    double snapshot_interval_time_code,
    double next_snapshot_time_code,
    double restart_resume_dt_time_code) {
  if (!write_outputs_enabled || (!snapshot_due && !checkpoint_due)) {
    return false;
  }

#if !COSMOSIM_ENABLE_HDF5
  throw std::runtime_error(
      "runtime outputs requested, but this build lacks HDF5 support. Reconfigure with COSMOSIM_ENABLE_HDF5=ON.");
#else
  if (integrator_state.step_index == 0) {
    return false;
  }
  if (checkpoint_due) {
    core::assertCanWriteCheckpointAtBoundary(integrator_state, scheduler.currentTick());
    if (gas_cell_scheduler.currentTick() != scheduler.currentTick()) {
      throw std::runtime_error("output checkpoint requires particle and gas-cell schedulers on one integer tick");
    }
  }
  if (snapshot_due) {
    core::assertCanWriteSnapshotAtBoundary(integrator_state);
  }
  if ((snapshot_due || checkpoint_due) && !core::isOutputSafeBoundary(integrator_state.last_completed_boundary_kind)) {
    return false;
  }

  bool output_flushed = false;
  if (snapshot_due) {
    core::ProvenanceRecord snapshot_provenance =
        makeGravityAwareProvenanceRecord(frozen_config, config);
    snapshot_provenance.gravity_treepm_decomposition_epoch =
        gravity_state.decompositionEpoch();
    io::SnapshotWritePayload snapshot_payload;
    snapshot_payload.state = &state;
    snapshot_payload.config = &config;
    snapshot_payload.normalized_config_text = frozen_config.normalized_text;
    snapshot_payload.provenance = snapshot_provenance;
    report.snapshot_path = report.run_directory / formatIndexedRankedFileStem(config.output.output_stem, integrator_state.step_index, gravity_state.runtimeTopology().world_size, gravity_state.runtimeTopology().world_rank);
    io::writeGadgetArepoSnapshotHdf5(report.snapshot_path, snapshot_payload);
    report.snapshot_roundtrip_executed = true;
    const io::SnapshotReadResult snapshot_read = io::readGadgetArepoSnapshotHdf5(report.snapshot_path, config);
    const SnapshotRoundtripVerification snapshot_verification =
        verifySnapshotRoundtrip(
        snapshot_read,
        state,
        config,
        frozen_config.normalized_text,
        snapshot_provenance);
    report.snapshot_roundtrip_ok = snapshot_verification.ok;
    profiler.recordEvent(core::RuntimeEvent{
        .event_kind = "snapshot.write.complete",
        .severity = report.snapshot_roundtrip_ok ? core::RuntimeEventSeverity::kInfo : core::RuntimeEventSeverity::kWarning,
        .subsystem = "io.snapshot",
        .step_index = integrator_state.step_index,
        .simulation_time_code = integrator_state.current_time_code,
        .scale_factor = integrator_state.current_scale_factor,
        .message = report.snapshot_roundtrip_ok
            ? "snapshot output written and scientifically verified"
            : "snapshot output failed scientific roundtrip verification",
        .payload = {{"path", report.snapshot_path.string()},
                    {"verification", snapshot_verification.detail}},
    });
    if (!report.snapshot_roundtrip_ok) {
      throw std::runtime_error(
          "snapshot scientific roundtrip verification failed: " +
          snapshot_verification.detail);
    }
    output_flushed = true;
  }

  if (checkpoint_due) {
    core::IntegratorState restart_integrator_state = integrator_state;
    if (restart_resume_dt_time_code > 0.0) {
      // The completed step may have been shortened to land on an output
      // event.  Persist the pre-clip proposal as the next-step interval so a
      // resumed workflow follows the same timeline as an uninterrupted run.
      restart_integrator_state.dt_time_code = restart_resume_dt_time_code;
    }
    io::RestartWritePayload restart_payload;
    restart_payload.persistent_state.simulation_state = &state;
    restart_payload.integrator_state = &restart_integrator_state;
    restart_payload.scheduler = &scheduler;
    restart_payload.gas_cell_scheduler = &gas_cell_scheduler;
    const io::GravityForceCachePersistentState gravity_force_cache =
        gravity_state.exportRestartForceCache(state);
    restart_payload.gravity_force_cache = &gravity_force_cache;
    restart_payload.provenance =
        makeGravityAwareProvenanceRecord(frozen_config, config);
    restart_payload.normalized_config_text = frozen_config.normalized_text;
    restart_payload.normalized_config_hash_hex = frozen_config.provenance.config_hash_hex;
    restart_payload.distributed_gravity_state.schema_version = 2;
    restart_payload.distributed_gravity_state.decomposition_epoch = gravity_state.decompositionEpoch();
    restart_payload.distributed_gravity_state.world_size = gravity_state.runtimeTopology().world_size;
    restart_payload.distributed_gravity_state.pm_grid_nx = gravity_state.pmGridShape().nx;
    restart_payload.distributed_gravity_state.pm_grid_ny = gravity_state.pmGridShape().ny;
    restart_payload.distributed_gravity_state.pm_grid_nz = gravity_state.pmGridShape().nz;
    restart_payload.distributed_gravity_state.pm_decomposition_mode =
        pmDecompositionModeName(config.numerics.treepm_pm_decomposition_mode);
    restart_payload.distributed_gravity_state.gravity_kick_opportunity =
        integrator_state.pm_sync_state.gravityKickOpportunity();
    restart_payload.distributed_gravity_state.pm_update_cadence_steps =
        integrator_state.pm_sync_state.cadenceSteps();
    restart_payload.distributed_gravity_state.long_range_field_version =
        integrator_state.pm_sync_state.fieldVersion();
    restart_payload.distributed_gravity_state.last_long_range_refresh_opportunity =
        integrator_state.pm_sync_state.lastRefreshOpportunity();
    restart_payload.distributed_gravity_state.long_range_field_built_step_index =
        integrator_state.pm_sync_state.lastRefreshStepIndex();
    restart_payload.distributed_gravity_state.long_range_field_built_scale_factor =
        integrator_state.pm_sync_state.lastRefreshScaleFactor();
    restart_payload.distributed_gravity_state.long_range_restart_policy = "deterministic_rebuild";
    restart_payload.output_cadence_state.output_enabled = write_outputs_enabled;
    restart_payload.output_cadence_state.write_restarts = config.output.write_restarts;
    restart_payload.output_cadence_state.snapshot_due = false;
    restart_payload.output_cadence_state.checkpoint_due = false;
    restart_payload.output_cadence_state.last_completed_step_index = integrator_state.step_index;
    restart_payload.output_cadence_state.snapshot_interval_steps = snapshot_interval_steps;
    restart_payload.output_cadence_state.next_snapshot_step_index = next_snapshot_step_index;
    restart_payload.output_cadence_state.snapshot_interval_time_code =
        snapshot_interval_time_code;
    restart_payload.output_cadence_state.next_snapshot_time_code =
        next_snapshot_time_code;
    restart_payload.output_cadence_state.snapshot_stem = config.output.output_stem;
    restart_payload.output_cadence_state.restart_stem = config.output.restart_stem;
    restart_payload.stochastic_state = buildStochasticPersistentState(
        config,
        integrator_state,
        static_cast<std::uint32_t>(std::max(gravity_state.runtimeTopology().world_rank, 0)));
    restart_payload.distributed_gravity_state.owning_rank_by_item.reserve(state.particle_sidecar.owning_rank.size());
    for (const std::uint32_t owner : state.particle_sidecar.owning_rank) {
      restart_payload.distributed_gravity_state.owning_rank_by_item.push_back(static_cast<int>(owner));
    }
    const std::size_t world_size = static_cast<std::size_t>(gravity_state.runtimeTopology().world_size);
    restart_payload.distributed_gravity_state.pm_slab_begin_x_by_rank.resize(world_size, 0);
    restart_payload.distributed_gravity_state.pm_slab_end_x_by_rank.resize(world_size, 0);
    for (std::size_t rank = 0; rank < world_size; ++rank) {
      const parallel::PmSlabRange range = parallel::pmOwnedXRangeForRank(
          gravity_state.pmGridShape().nx,
          static_cast<int>(world_size),
          static_cast<int>(rank));
      restart_payload.distributed_gravity_state.pm_slab_begin_x_by_rank[rank] = range.begin_x;
      restart_payload.distributed_gravity_state.pm_slab_end_x_by_rank[rank] = range.end_x;
    }
    restart_payload.provenance.gravity_treepm_decomposition_epoch =
        restart_payload.distributed_gravity_state.decomposition_epoch;
    restart_payload.provenance.gravity_treepm_restart_world_size =
        restart_payload.distributed_gravity_state.world_size;
    restart_payload.provenance.gravity_treepm_restart_pm_grid =
        std::to_string(restart_payload.distributed_gravity_state.pm_grid_nx) + "x" +
        std::to_string(restart_payload.distributed_gravity_state.pm_grid_ny) + "x" +
        std::to_string(restart_payload.distributed_gravity_state.pm_grid_nz);
    restart_payload.provenance.gravity_treepm_restart_slab_signature =
        pmSlabSignature(restart_payload.distributed_gravity_state);
    restart_payload.provenance.gravity_treepm_restart_kick_opportunity =
        restart_payload.distributed_gravity_state.gravity_kick_opportunity;
    restart_payload.provenance.gravity_treepm_restart_field_version =
        restart_payload.distributed_gravity_state.long_range_field_version;
    restart_payload.provenance.gravity_treepm_long_range_restart_policy =
        restart_payload.distributed_gravity_state.long_range_restart_policy;

    report.restart_path = report.run_directory / formatIndexedRankedFileStem(config.output.restart_stem, integrator_state.step_index, gravity_state.runtimeTopology().world_size, gravity_state.runtimeTopology().world_rank);
    io::writeRestartCheckpointHdf5(report.restart_path, restart_payload);
    report.restart_roundtrip_executed = true;
    const io::RestartReadResult restart_read = io::readRestartCheckpointHdf5(report.restart_path);
    const auto compatibility = parallel::evaluateDistributedRestartCompatibility(
        restart_read.distributed_gravity_state,
        gravity_state.runtimeTopology());
    const bool restart_rank_qualified_name =
        gravity_state.runtimeTopology().world_size == 1 ||
        report.restart_path.filename().string().find("_rank") != std::string::npos;
    report.restart_roundtrip_ok = restartRuntimeStateExactlyEquivalent(
        restart_read.state, state) &&
        restart_read.integrator_state.pm_refresh_enabled == integrator_state.pm_refresh_enabled &&
        gravityForceCachesEquivalent(restart_read.gravity_force_cache, gravity_force_cache) &&
        restart_read.scheduler_state.current_tick == scheduler.currentTick() &&
        restart_read.distributed_gravity_state.decomposition_epoch ==
            restart_payload.distributed_gravity_state.decomposition_epoch &&
        restart_read.distributed_gravity_state.owning_rank_by_item.size() == state.particle_sidecar.owning_rank.size() &&
        restart_read.distributed_gravity_state.owning_rank_by_item == restart_payload.distributed_gravity_state.owning_rank_by_item &&
        restart_read.distributed_gravity_state.pm_slab_begin_x_by_rank == restart_payload.distributed_gravity_state.pm_slab_begin_x_by_rank &&
        restart_read.distributed_gravity_state.pm_slab_end_x_by_rank == restart_payload.distributed_gravity_state.pm_slab_end_x_by_rank &&
        restart_read.distributed_gravity_state.gravity_kick_opportunity == restart_payload.distributed_gravity_state.gravity_kick_opportunity &&
        restart_read.distributed_gravity_state.pm_update_cadence_steps == restart_payload.distributed_gravity_state.pm_update_cadence_steps &&
        restart_read.distributed_gravity_state.long_range_field_version == restart_payload.distributed_gravity_state.long_range_field_version &&
        restart_read.distributed_gravity_state.last_long_range_refresh_opportunity ==
            restart_payload.distributed_gravity_state.last_long_range_refresh_opportunity &&
        restart_read.distributed_gravity_state.long_range_field_built_step_index ==
            restart_payload.distributed_gravity_state.long_range_field_built_step_index &&
        std::abs(
            restart_read.distributed_gravity_state.long_range_field_built_scale_factor -
            restart_payload.distributed_gravity_state.long_range_field_built_scale_factor) <= 1.0e-12 &&
        restart_read.distributed_gravity_state.long_range_restart_policy ==
            restart_payload.distributed_gravity_state.long_range_restart_policy &&
        restart_read.output_cadence_state.output_enabled == restart_payload.output_cadence_state.output_enabled &&
        restart_read.output_cadence_state.write_restarts == restart_payload.output_cadence_state.write_restarts &&
        restart_read.output_cadence_state.next_snapshot_step_index ==
            restart_payload.output_cadence_state.next_snapshot_step_index &&
        restart_read.output_cadence_state.snapshot_interval_time_code ==
            restart_payload.output_cadence_state.snapshot_interval_time_code &&
        restart_read.output_cadence_state.next_snapshot_time_code ==
            restart_payload.output_cadence_state.next_snapshot_time_code &&
        restart_read.output_cadence_state.snapshot_stem == restart_payload.output_cadence_state.snapshot_stem &&
        restart_read.output_cadence_state.restart_stem == restart_payload.output_cadence_state.restart_stem &&
        stochasticStatesEquivalent(restart_read.stochastic_state, restart_payload.stochastic_state) &&
        restart_rank_qualified_name &&
        compatibility.compatible();
    profiler.recordEvent(core::RuntimeEvent{
        .event_kind = "restart.write.complete",
        .severity = report.restart_roundtrip_ok ? core::RuntimeEventSeverity::kInfo : core::RuntimeEventSeverity::kWarning,
        .subsystem = "io.restart",
        .step_index = integrator_state.step_index,
        .simulation_time_code = integrator_state.current_time_code,
        .scale_factor = integrator_state.current_scale_factor,
        .message = "restart checkpoint written and verified",
        .payload = {{"path", report.restart_path.string()},
                    {"restart_schema", io::restartSchema().name},
                    {"restart_schema_version", std::to_string(io::restartSchema().version)},
                    {"boundary_kind", std::string(core::stepBoundaryKindName(integrator_state.last_completed_boundary_kind))},
                    {"scheduler_current_tick", std::to_string(scheduler.currentTick())},
                    {"scheduler_max_bin", std::to_string(scheduler.maxBin())},
                    {"pm_cadence_steps", std::to_string(integrator_state.pm_sync_state.cadenceSteps())},
                    {"pm_gravity_kick_opportunity", std::to_string(integrator_state.pm_sync_state.gravityKickOpportunity())},
                    {"pm_field_version", std::to_string(integrator_state.pm_sync_state.fieldVersion())},
                    {"pm_long_range_field_valid", integrator_state.pm_long_range_field_valid ? "true" : "false"},
                    {"output_next_snapshot_step_index", std::to_string(restart_payload.output_cadence_state.next_snapshot_step_index)},
                    {"output_next_snapshot_time_code", formatRuntimeDouble(
                        restart_payload.output_cadence_state.next_snapshot_time_code)},
                    {"stochastic_module_count", std::to_string(restart_payload.stochastic_state.modules.size())}},
    });
    profiler.recordEvent(core::RuntimeEvent{
        .event_kind = "restart.read.complete",
        .severity = report.restart_roundtrip_ok ? core::RuntimeEventSeverity::kInfo : core::RuntimeEventSeverity::kWarning,
        .subsystem = "io.restart",
        .step_index = restart_read.integrator_state.step_index,
        .simulation_time_code = restart_read.integrator_state.current_time_code,
        .scale_factor = restart_read.integrator_state.current_scale_factor,
        .message = "restart checkpoint read and validated",
        .payload = {{"path", report.restart_path.string()},
                    {"restart_schema", restart_read.diagnostics.restart_schema_name},
                    {"restart_schema_version", std::to_string(restart_read.diagnostics.restart_schema_version)},
                    {"boundary_kind", restart_read.diagnostics.last_completed_boundary_kind},
                    {"restart_safe", restart_read.diagnostics.restart_safe ? "true" : "false"},
                    {"scheduler_current_tick", std::to_string(restart_read.diagnostics.scheduler_current_tick)},
                    {"scheduler_max_bin", std::to_string(restart_read.diagnostics.scheduler_max_bin)},
                    {"scheduler_element_count", std::to_string(restart_read.diagnostics.scheduler_element_count)},
                    {"scheduler_active_count", std::to_string(restart_read.diagnostics.scheduler_active_count)},
                    {"scheduler_pending_transition_count", std::to_string(restart_read.diagnostics.scheduler_pending_transition_count)},
                    {"pm_cadence_steps", std::to_string(restart_read.diagnostics.pm_cadence_steps)},
                    {"pm_gravity_kick_opportunity", std::to_string(restart_read.diagnostics.pm_gravity_kick_opportunity)},
                    {"pm_field_version", std::to_string(restart_read.diagnostics.pm_field_version)},
                    {"pm_long_range_field_valid", restart_read.diagnostics.pm_long_range_field_valid ? "true" : "false"},
                    {"output_next_snapshot_step_index", std::to_string(restart_read.diagnostics.output_next_snapshot_step_index)},
                    {"stochastic_module_count", std::to_string(restart_read.diagnostics.stochastic_module_count)},
                    {"payload_hash_hex", restart_read.payload_hash_hex}},
    });
    output_flushed = true;
  }
  return output_flushed;
#endif
}

[[nodiscard]] std::string joinRestartCompatibilityMessages(
    std::span<const std::string> mismatch_messages) {
  std::ostringstream stream;
  for (std::size_t i = 0; i < mismatch_messages.size(); ++i) {
    if (i != 0U) {
      stream << "; ";
    }
    stream << mismatch_messages[i];
  }
  return stream.str();
}

[[nodiscard]] parallel::LocalOwnershipIdentitySummary reduceLocalParticleIdentitySummary(
    const core::SimulationState& state,
    const parallel::MpiContext& mpi_context) {
  const parallel::LocalOwnershipIdentitySummary local =
      parallel::summarizeLocalOwnedParticleIds(state.particle_sidecar.particle_id);
  const std::uint64_t unique_rank_count =
      mpi_context.allreduceSumUint64(local.local_particle_ids_unique ? 1ULL : 0ULL);
  return parallel::LocalOwnershipIdentitySummary{
      .local_owned_count = mpi_context.allreduceSumUint64(local.local_owned_count),
      .local_particle_id_sum = mpi_context.allreduceSumUint64(local.local_particle_id_sum),
      .local_particle_id_square_sum = mpi_context.allreduceSumUint64(local.local_particle_id_square_sum),
      .local_particle_id_xor = mpi_context.allreduceXorUint64(local.local_particle_id_xor),
      .local_particle_ids_unique =
          unique_rank_count == static_cast<std::uint64_t>(mpi_context.worldSize()),
  };
}

void validateRestartResumeTopologyOrThrowImpl(
    const io::RestartReadResult& restart,
    const core::SimulationConfig& config,
    const core::FrozenConfig& frozen_config,
    const parallel::MpiContext& mpi_context) {
  const auto require_bin_zero_scheduler = [](
      const core::TimeBinPersistentState& scheduler_state,
      std::string_view label) {
    if (scheduler_state.max_bin != 0U) {
      throw std::runtime_error(
          "ReferenceWorkflow restart topology validation failed: " +
          std::string(label) + " max_bin must be zero for production KDK");
    }
    for (const std::uint8_t bin : scheduler_state.bin_index) {
      if (bin != 0U) {
        throw std::runtime_error(
            "ReferenceWorkflow restart topology validation failed: " +
            std::string(label) + " contains a nonzero committed time bin");
      }
    }
    for (const std::uint8_t pending : scheduler_state.pending_bin_index) {
      if (pending != 0U &&
          pending != core::HierarchicalTimeBinScheduler::k_unset_pending_bin) {
        throw std::runtime_error(
            "ReferenceWorkflow restart topology validation failed: " +
            std::string(label) + " contains a nonzero pending time bin");
      }
    }
  };
  require_bin_zero_scheduler(restart.scheduler_state, "particle scheduler");
  require_bin_zero_scheduler(restart.gas_cell_scheduler_state, "gas-cell scheduler");
  if (restart.normalized_config_hash_hex != frozen_config.provenance.config_hash_hex) {
    throw std::runtime_error(
        "ReferenceWorkflow restart topology validation failed: normalized config hash mismatch: expected=" +
        frozen_config.provenance.config_hash_hex + ", observed=" + restart.normalized_config_hash_hex);
  }
  const core::CudaRuntimeInfo cuda_runtime = core::queryCudaRuntime();
  const parallel::DistributedExecutionTopology runtime_topology =
      parallel::buildDistributedExecutionTopology(
          static_cast<std::size_t>(config.numerics.treepm_pm_grid_nx),
          static_cast<std::size_t>(config.numerics.treepm_pm_grid_ny),
          static_cast<std::size_t>(config.numerics.treepm_pm_grid_nz),
          mpi_context,
          config.parallel.mpi_ranks_expected,
          config.parallel.gpu_devices,
          cuda_runtime.runtime_available,
          cuda_runtime.visible_device_count,
          pmDecompositionModeName(config.numerics.treepm_pm_decomposition_mode));
  const parallel::DistributedRestartCompatibilityReport compatibility =
      parallel::evaluateDistributedRestartCompatibility(
          restart.distributed_gravity_state,
          runtime_topology);
  if (!compatibility.compatible()) {
    throw std::runtime_error(
        "ReferenceWorkflow restart topology validation failed for rank " +
        std::to_string(mpi_context.worldRank()) + "/" + std::to_string(mpi_context.worldSize()) +
        ": " + joinRestartCompatibilityMessages(compatibility.mismatch_messages));
  }
  if (restart.distributed_gravity_state.owning_rank_by_item.size() !=
      restart.state.particle_sidecar.owning_rank.size()) {
    throw std::runtime_error(
        "ReferenceWorkflow restart topology validation failed: /distributed_gravity owning_rank_by_item count " +
        std::to_string(restart.distributed_gravity_state.owning_rank_by_item.size()) +
        " does not match local particle rows " +
        std::to_string(restart.state.particle_sidecar.owning_rank.size()));
  }
  for (std::size_t row = 0; row < restart.state.particle_sidecar.owning_rank.size(); ++row) {
    const int observed_owner = static_cast<int>(restart.state.particle_sidecar.owning_rank[row]);
    const int restart_owner = restart.distributed_gravity_state.owning_rank_by_item[row];
    if (observed_owner != restart_owner) {
      throw std::runtime_error(
          "ReferenceWorkflow restart topology validation failed: local particle row " +
          std::to_string(row) + " owner mismatch: state=" + std::to_string(observed_owner) +
          ", /distributed_gravity=" + std::to_string(restart_owner));
    }
    if (restart_owner < 0 || restart_owner >= mpi_context.worldSize()) {
      throw std::runtime_error(
          "ReferenceWorkflow restart topology validation failed: local particle row " +
          std::to_string(row) + " owner rank is outside runtime world: owner=" +
          std::to_string(restart_owner) + ", world_size=" + std::to_string(mpi_context.worldSize()));
    }
  }
}

}  // namespace

void validateRestartResumeTopologyOrThrow(
    const io::RestartReadResult& restart,
    const core::SimulationConfig& config,
    const core::FrozenConfig& frozen_config,
    const parallel::MpiContext& mpi_context) {
  validateRestartResumeTopologyOrThrowImpl(
      restart, config, frozen_config, mpi_context);
}

PendingOutputBoundary initializeOutputCadence(
    const core::SimulationConfig& config,
    const ReferenceWorkflowOptions& options,
    const core::IntegratorState& integrator_state,
    const io::OutputCadencePersistentState* restored_output) {
  PendingOutputBoundary pending;
  if (options.write_outputs && config.output.snapshot_interval_steps > 0) {
    pending.snapshot_interval_steps =
        static_cast<std::uint64_t>(config.output.snapshot_interval_steps);
    if (restored_output != nullptr && restored_output->output_enabled) {
      if (restored_output->snapshot_interval_steps !=
          pending.snapshot_interval_steps) {
        throw std::runtime_error(
            "restart output step interval does not match the active normalized configuration");
      }
      pending.next_snapshot_step_index =
          restored_output->next_snapshot_step_index;
    } else {
      const std::uint64_t completed_intervals =
          integrator_state.step_index / pending.snapshot_interval_steps;
      if (completed_intervals == std::numeric_limits<std::uint64_t>::max() ||
          completed_intervals + 1U >
              std::numeric_limits<std::uint64_t>::max() /
                  pending.snapshot_interval_steps) {
        throw std::overflow_error(
            "output step cadence could not represent its next event");
      }
      pending.next_snapshot_step_index =
          (completed_intervals + 1U) * pending.snapshot_interval_steps;
    }
    if (pending.next_snapshot_step_index <= integrator_state.step_index) {
      throw std::runtime_error(
          "restart output step cadence does not provide a future ordered event");
    }
  } else if (restored_output != nullptr && options.write_outputs &&
             restored_output->output_enabled &&
             restored_output->snapshot_interval_steps > 0U) {
    throw std::runtime_error(
        "restart contains an active step output cadence that the current configuration disables");
  }

  if (options.write_outputs &&
      config.output.snapshot_interval_time_code > 0.0) {
    pending.snapshot_interval_time_code =
        config.output.snapshot_interval_time_code;
    if (restored_output != nullptr &&
        restored_output->snapshot_interval_time_code > 0.0) {
      if (!timelineTimesEqual(
              restored_output->snapshot_interval_time_code,
              config.output.snapshot_interval_time_code)) {
        throw std::runtime_error(
            "restart output code-time interval does not match the active normalized configuration");
      }
      pending.next_snapshot_time_code =
          restored_output->next_snapshot_time_code;
    } else {
      pending.next_snapshot_time_code = nextOrderedTimeEvent(
          config.numerics.t_code_begin,
          integrator_state.current_time_code,
          config.output.snapshot_interval_time_code);
    }
    if (!std::isfinite(pending.next_snapshot_time_code) ||
        pending.next_snapshot_time_code <= integrator_state.current_time_code) {
      throw std::runtime_error(
          "restart output code-time cadence does not provide a future ordered event");
    }
  } else if (restored_output != nullptr && options.write_outputs &&
             restored_output->snapshot_interval_time_code > 0.0) {
    throw std::runtime_error(
        "restart contains an active code-time output cadence that the current configuration disables");
  }
  return pending;
}

[[nodiscard]] bool timelineTimesEqual(double lhs, double rhs) {
  const double scale = std::max({1.0, std::abs(lhs), std::abs(rhs)});
  return std::abs(lhs - rhs) <= 8.0 * std::numeric_limits<double>::epsilon() * scale;
}

[[nodiscard]] double nextOrderedTimeEvent(
    double timeline_begin_time_code,
    double current_time_code,
    double interval_time_code) {
  if (!std::isfinite(timeline_begin_time_code) || !std::isfinite(current_time_code) ||
      !std::isfinite(interval_time_code) || interval_time_code <= 0.0) {
    throw std::invalid_argument("output code-time cadence requires finite times and a positive interval");
  }
  const double completed_intervals = std::max(
      0.0,
      std::floor((current_time_code - timeline_begin_time_code) / interval_time_code));
  double next_time_code =
      std::fma(completed_intervals + 1.0, interval_time_code, timeline_begin_time_code);
  if (next_time_code <= current_time_code || timelineTimesEqual(next_time_code, current_time_code)) {
    next_time_code = std::fma(completed_intervals + 2.0, interval_time_code, timeline_begin_time_code);
  }
  if (!std::isfinite(next_time_code) || next_time_code <= current_time_code) {
    throw std::overflow_error("output code-time cadence could not represent its next ordered event");
  }
  return next_time_code;
}

void latchOutputRequestForCompletedStep(
    const core::SimulationConfig& config,
    const ReferenceWorkflowOptions& options,
    std::uint64_t next_step_index,
    double next_time_code,
    PendingOutputBoundary& pending) {
  if (!options.write_outputs) {
    return;
  }
  const bool step_event_due = pending.snapshot_interval_steps > 0 &&
      next_step_index >= pending.next_snapshot_step_index;
  const bool time_event_due = pending.snapshot_interval_time_code > 0.0 &&
      timelineTimesEqual(next_time_code, pending.next_snapshot_time_code);
  if (!step_event_due && !time_event_due) {
    return;
  }
  pending.snapshot_due = true;
  pending.checkpoint_due = pending.checkpoint_due || config.output.write_restarts;
  pending.step_event_due = pending.step_event_due || step_event_due;
  pending.time_event_due = pending.time_event_due || time_event_due;
}

[[nodiscard]] core::StepBoundaryKind requestedBoundaryForPendingOutput(const PendingOutputBoundary& pending) {
  if (pending.checkpoint_due) {
    return core::StepBoundaryKind::kCheckpointPoint;
  }
  if (pending.snapshot_due) {
    return core::StepBoundaryKind::kSnapshotPoint;
  }
  return core::StepBoundaryKind::kGlobalSynchronizationPoint;
}

OutputRestartRuntime::OutputRestartRuntime(
    const core::FrozenConfig& frozen_config,
    const core::SimulationConfig& config,
    const core::HierarchicalTimeBinScheduler& scheduler,
    const core::HierarchicalTimeBinScheduler& gas_cell_scheduler,
    const GravityRestartStateProvider& gravity_state,
    ReferenceWorkflowReport& report,
    core::ProfilerSession& profiler,
    PendingOutputBoundary& pending_output,
    bool write_outputs_enabled)
    : m_frozen_config(frozen_config),
      m_config(config),
      m_scheduler(scheduler),
      m_gas_cell_scheduler(gas_cell_scheduler),
      m_gravity_state(gravity_state),
      m_report(report),
      m_profiler(profiler),
      m_pending_output(pending_output),
      m_write_outputs_enabled(write_outputs_enabled) {}

void OutputRestartRuntime::execute(OutputRestartStageView& view) {
  view.requireFresh();
  core::StepContext& context = internal::RuntimeStageAccess::outputRestartContext(
      view,
      {{RuntimeResourceKey::kParticlePosition, RuntimeResourceAccessMode::kRead},
       {RuntimeResourceKey::kParticleVelocity, RuntimeResourceAccessMode::kRead},
       {RuntimeResourceKey::kHydroConservedState, RuntimeResourceAccessMode::kRead},
       {RuntimeResourceKey::kSourceMutationState, RuntimeResourceAccessMode::kRead},
       {RuntimeResourceKey::kMigrationOwnership, RuntimeResourceAccessMode::kRead},
       {RuntimeResourceKey::kSchedulerTruth, RuntimeResourceAccessMode::kRead},
       {RuntimeResourceKey::kIntegratorTruth, RuntimeResourceAccessMode::kRead},
       {RuntimeResourceKey::kOutputRestartState, RuntimeResourceAccessMode::kReadWrite},
       {RuntimeResourceKey::kDiagnostics, RuntimeResourceAccessMode::kWrite}});
  if (context.stage != core::IntegrationStage::kOutputCheck) {
    throw std::logic_error(
        "output/restart runtime received an unregistered stage");
  }
  if (!m_pending_output.snapshot_due && !m_pending_output.checkpoint_due) {
    return;
  }
  if (!context.boundary.output_safe || !context.boundary.restart_safe) {
    if (m_pending_output.time_event_due) {
      throw std::runtime_error(
          "required code-time output event reached a boundary that is not output/restart safe");
    }
    if (m_pending_output.checkpoint_due) {
      core::assertCanWriteCheckpointAtBoundary(
          context.integrator_state, m_scheduler.currentTick());
    }
    return;
  }

  const double persisted_next_snapshot_time_code =
      m_pending_output.time_event_due
      ? nextOrderedTimeEvent(
            m_config.numerics.t_code_begin,
            context.integrator_state.current_time_code,
            m_pending_output.snapshot_interval_time_code)
      : m_pending_output.next_snapshot_time_code;
  std::uint64_t persisted_next_snapshot_step_index =
      m_pending_output.next_snapshot_step_index;
  if (m_pending_output.step_event_due) {
    do {
      if (persisted_next_snapshot_step_index >
          std::numeric_limits<std::uint64_t>::max() -
              m_pending_output.snapshot_interval_steps) {
        throw std::overflow_error(
            "output step cadence overflowed its ordered timeline");
      }
      persisted_next_snapshot_step_index +=
          m_pending_output.snapshot_interval_steps;
    } while (persisted_next_snapshot_step_index <=
             context.integrator_state.step_index);
  }

  const bool output_flushed = maybeWriteOutputs(
      m_frozen_config,
      m_config,
      context.state,
      context.integrator_state,
      m_scheduler,
      m_gas_cell_scheduler,
      m_gravity_state,
      m_report,
      m_profiler,
      m_write_outputs_enabled,
      m_pending_output.snapshot_due,
      m_pending_output.checkpoint_due,
      m_pending_output.snapshot_interval_steps,
      persisted_next_snapshot_step_index,
      m_pending_output.snapshot_interval_time_code,
      persisted_next_snapshot_time_code,
      m_pending_output.restart_resume_dt_time_code);
  if (output_flushed) {
    m_pending_output.snapshot_due = false;
    m_pending_output.checkpoint_due = false;
    m_pending_output.step_event_due = false;
    m_pending_output.time_event_due = false;
    m_pending_output.next_snapshot_step_index =
        persisted_next_snapshot_step_index;
    m_pending_output.next_snapshot_time_code =
        persisted_next_snapshot_time_code;
    m_pending_output.restart_resume_dt_time_code = 0.0;
  }
}

}  // namespace cosmosim::workflows::internal
