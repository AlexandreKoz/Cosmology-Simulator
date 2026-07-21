#include "cosmosim/workflows/reference_workflow.hpp"
#include "cosmosim/workflows/gravity_runtime.hpp"
#include "cosmosim/workflows/hydro_amr_runtime.hpp"
#include "cosmosim/workflows/runtime_capabilities.hpp"
#include "cosmosim/workflows/runtime_services.hpp"
#include "cosmosim/workflows/time_coordinator.hpp"
#include "workflows/internal/initial_condition_runtime.hpp"
#include "cosmosim/workflows/migration_balance_runtime.hpp"
#include "workflows/internal/output_restart_runtime.hpp"
#include "workflows/internal/reference_runtime_composition.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstddef>
#include <filesystem>
#include <iomanip>
#include <limits>
#include <numeric>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <cstring>
#include <exception>
#include <utility>
#include <vector>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/core/cosmology.hpp"
#include "cosmosim/core/cuda_runtime.hpp"
#include "cosmosim/core/profiling.hpp"
#include "cosmosim/core/simulation_mode.hpp"
#include "cosmosim/core/time_integration.hpp"
#include "cosmosim/core/units.hpp"
#include "cosmosim/gravity/tree_pm_coupling.hpp"
#include "cosmosim/io/io_contract.hpp"
#include "cosmosim/io/restart_checkpoint.hpp"
#include "cosmosim/io/snapshot_hdf5.hpp"

namespace cosmosim::workflows {
namespace {

[[nodiscard]] std::string formatRuntimeDouble(double value) {
  std::ostringstream stream;
  stream << std::scientific
         << std::setprecision(std::numeric_limits<double>::max_digits10)
         << value;
  return stream.str();
}


[[nodiscard]] std::string treePmAssignmentSchemeName(core::TreePmAssignmentScheme assignment_scheme) {
  switch (assignment_scheme) {
    case core::TreePmAssignmentScheme::kCic:
      return "cic";
    case core::TreePmAssignmentScheme::kTsc:
      return "tsc";
  }
  throw std::runtime_error("unhandled TreePm assignment scheme enum value");
}


[[nodiscard]] std::filesystem::path computeRunDirectory(
    const core::SimulationConfig& config,
    const std::filesystem::path* output_root_override) {
  const std::filesystem::path output_root =
      (output_root_override != nullptr && !output_root_override->empty())
      ? *output_root_override
      : std::filesystem::path(config.output.output_directory);
  return output_root / config.output.run_name;
}

[[nodiscard]] std::filesystem::path rankQualifiedRunDirectory(
    std::filesystem::path run_directory,
    int world_size,
    int world_rank) {
  if (world_size <= 1) {
    return run_directory;
  }
  std::ostringstream rank_suffix;
  rank_suffix << "_rank" << std::setw(3) << std::setfill('0') << world_rank;
  return run_directory.parent_path() / (run_directory.filename().string() + rank_suffix.str());
}

[[nodiscard]] std::filesystem::path resolveConfigRelativePath(
    const core::FrozenConfig& frozen_config,
    const std::filesystem::path& candidate) {
  if (candidate.is_absolute()) {
    return candidate;
  }

  const std::filesystem::path source_path(frozen_config.provenance.source_name);
  if (!source_path.empty() && source_path.has_parent_path()) {
    return source_path.parent_path() / candidate;
  }
  return candidate;
}

[[nodiscard]] bool repeatedCanonicalOrder(const std::vector<std::string>& observed) {
  const auto canonical = core::StageScheduler::kickDriftKickOrder();
  if (observed.empty() || observed.size() % canonical.size() != 0) {
    return false;
  }

  for (std::size_t chunk = 0; chunk < observed.size(); chunk += canonical.size()) {
    for (std::size_t i = 0; i < canonical.size(); ++i) {
      if (observed[chunk + i] != core::integrationStageName(canonical[i])) {
        return false;
      }
    }
  }
  return true;
}

void fnv1aMix(std::uint64_t& hash, const void* data, std::size_t size) {
  constexpr std::uint64_t k_fnv_offset_basis = 1469598103934665603ULL;
  constexpr std::uint64_t k_fnv_prime = 1099511628211ULL;
  if (hash == 0) {
    hash = k_fnv_offset_basis;
  }
  const auto* bytes = static_cast<const unsigned char*>(data);
  for (std::size_t i = 0; i < size; ++i) {
    hash ^= static_cast<std::uint64_t>(bytes[i]);
    hash *= k_fnv_prime;
  }
}

void fnv1aMix(std::uint64_t& hash, double value) { fnv1aMix(hash, &value, sizeof(value)); }

void fnv1aMix(std::uint64_t& hash, std::uint64_t value) { fnv1aMix(hash, &value, sizeof(value)); }

[[nodiscard]] std::uint64_t computeStateDigest(const core::SimulationState& state, const core::IntegratorState& integrator_state) {
  std::uint64_t hash = 0;
  fnv1aMix(hash, static_cast<std::uint64_t>(state.particles.size()));
  fnv1aMix(hash, static_cast<std::uint64_t>(state.cells.size()));
  for (const double value : state.particles.position_x_comoving) {
    fnv1aMix(hash, value);
  }
  for (const double value : state.particles.position_y_comoving) {
    fnv1aMix(hash, value);
  }
  for (const double value : state.particles.position_z_comoving) {
    fnv1aMix(hash, value);
  }
  for (const double value : state.particles.velocity_x_peculiar) {
    fnv1aMix(hash, value);
  }
  for (const double value : state.particles.velocity_y_peculiar) {
    fnv1aMix(hash, value);
  }
  for (const double value : state.particles.velocity_z_peculiar) {
    fnv1aMix(hash, value);
  }
  for (const double value : state.particles.mass_code) {
    fnv1aMix(hash, value);
  }
  for (const double value : state.gas_cells.density_code) {
    fnv1aMix(hash, value);
  }
  for (const double value : state.gas_cells.internal_energy_code) {
    fnv1aMix(hash, value);
  }
  fnv1aMix(hash, integrator_state.current_time_code);
  fnv1aMix(hash, integrator_state.current_scale_factor);
  fnv1aMix(hash, integrator_state.step_index);
  return hash;
}

[[nodiscard]] std::uint64_t computeParticleIdSum(const core::SimulationState& state) {
  std::uint64_t sum = 0;
  for (const std::uint64_t particle_id : state.particle_sidecar.particle_id) {
    sum += particle_id;
  }
  return sum;
}

[[nodiscard]] std::uint64_t computeParticleIdSquareSum(const core::SimulationState& state) {
  std::uint64_t sum = 0;
  for (const std::uint64_t particle_id : state.particle_sidecar.particle_id) {
    sum += particle_id * particle_id;
  }
  return sum;
}

[[nodiscard]] std::uint64_t computeParticleIdXor(const core::SimulationState& state) {
  std::uint64_t value = 0;
  for (const std::uint64_t particle_id : state.particle_sidecar.particle_id) {
    value ^= particle_id;
  }
  return value;
}


void ensureRunDirectory(const std::filesystem::path& run_directory) {
  std::filesystem::create_directories(run_directory);
}

void flushCommonArtifacts(
    const core::FrozenConfig& frozen_config,
    core::ProfilerSession& profiler,
    ReferenceWorkflowReport& report) {
  ensureRunDirectory(report.run_directory);

  report.normalized_config_snapshot_path = report.run_directory / "normalized_config.param.txt";
  if (!std::filesystem::exists(report.normalized_config_snapshot_path)) {
    core::writeNormalizedConfigSnapshot(frozen_config, report.run_directory);
  }
  report.normalized_config_snapshot_written = true;

  report.profiler_json_path = report.run_directory / "profile.json";
  report.profiler_csv_path = report.run_directory / "profile.csv";
  report.operational_report_json_path = report.run_directory / "operational_events.json";
  core::writeProfilerReportJson(profiler, report.profiler_json_path);
  core::writeProfilerReportCsv(profiler, report.profiler_csv_path);
  core::writeOperationalReportJson(
      profiler,
      report.operational_report_json_path,
      frozen_config.config.output.run_name,
      frozen_config.provenance.config_hash_hex);
}


}  // namespace

ReferenceWorkflowRunner::ReferenceWorkflowRunner(core::FrozenConfig frozen_config)
    : m_frozen_config(std::move(frozen_config)) {}

const core::FrozenConfig& ReferenceWorkflowRunner::frozenConfig() const noexcept {
  return m_frozen_config;
}

ReferenceWorkflowReport ReferenceWorkflowRunner::run(const ReferenceWorkflowOptions& options) const {
  return runImpl(nullptr, options);
}

ReferenceWorkflowReport ReferenceWorkflowRunner::run(
    const std::filesystem::path& output_root_override,
    const ReferenceWorkflowOptions& options) const {
  return runImpl(&output_root_override, options);
}


ReferenceWorkflowReport ReferenceWorkflowRunner::runImpl(
    const std::filesystem::path* output_root_override,
    const ReferenceWorkflowOptions& options) const {
  const core::SimulationConfig& config = m_frozen_config.config;
  parallel::MpiContext mpi_context;
  mpi_context.validateExpectedWorldSizeOrThrow(config.parallel.mpi_ranks_expected);

  ReferenceWorkflowReport report;
  report.world_size = mpi_context.worldSize();
  report.world_rank = mpi_context.worldRank();
  report.run_directory = rankQualifiedRunDirectory(
      computeRunDirectory(config, output_root_override),
      report.world_size,
      report.world_rank);
  report.config_compatible = true;
  report.schema_compatible =
      config.schema_version == 1 && io::gadgetArepoSchemaMap().schema_version >= 2 &&
      io::isRestartSchemaCompatible(io::restartSchema().version);

  core::ProfilerSession profiler(true);
  const RuntimeServices runtime_services{
      .mpi_context = mpi_context,
      .profiler = profiler,
      .deterministic_execution = true};
  const FailureCoordinator failure_coordinator(runtime_services);
  const internal::MigrationBalanceRuntime migration_balance(
      config, runtime_services);
  profiler.recordEvent(core::RuntimeEvent{
      .event_kind = "config.freeze",
      .severity = core::RuntimeEventSeverity::kInfo,
      .subsystem = "core.config",
      .step_index = options.step_index,
      .simulation_time_code = config.numerics.t_code_begin,
      .scale_factor = 1.0,
      .message = "frozen configuration accepted for runtime workflow",
      .payload = {{"schema_version", std::to_string(config.schema_version)},
                  {"config_hash_hex", m_frozen_config.provenance.config_hash_hex},
                  {"source_name", m_frozen_config.provenance.source_name}},
  });

  try {
    ensureRunDirectory(report.run_directory);
    const RuntimeCapabilityReport runtime_capabilities =
        buildRuntimeCapabilityReport(config);
    validateRequestedRuntimeCapabilities(config, runtime_capabilities);
    report.runtime_capability_report_path =
        report.run_directory / "runtime_capabilities.json";
    writeRuntimeCapabilityReportJson(
        runtime_capabilities, report.runtime_capability_report_path);
    core::writeNormalizedConfigSnapshot(m_frozen_config, report.run_directory);
    report.normalized_config_snapshot_path = report.run_directory / "normalized_config.param.txt";
    report.normalized_config_snapshot_written = true;

    if (!report.schema_compatible) {
      throw std::runtime_error("runtime workflow schema compatibility validation failed");
    }

#if !COSMOSIM_ENABLE_HDF5
    if (options.write_outputs) {
      throw std::runtime_error(
          "this build lacks HDF5 support, so the config-driven runtime cannot emit snapshots/restarts; reconfigure with COSMOSIM_ENABLE_HDF5=ON or use a no-output internal test path");
    }
#endif

    const core::ModePolicy mode_policy = core::buildModePolicy(config.mode);
    core::validateModePolicy(config, mode_policy);

    const internal::InitialConditionRuntime initial_conditions(
        m_frozen_config, runtime_services);
    internal::InitialConditionStartupResult startup =
        initial_conditions.materialize(options, report.run_directory);
    report.ic_manifest_path = startup.manifest_path;
    core::SimulationState state = std::move(startup.state);
    profiler.setMemoryReport(core::collectSimulationMemoryReport(state));
    const core::MemoryReport* startup_memory_report = profiler.memoryReport();
    if (startup_memory_report != nullptr) {
      profiler.recordEvent(core::RuntimeEvent{
          .event_kind = "memory.startup_snapshot",
          .severity = core::RuntimeEventSeverity::kInfo,
          .subsystem = "core.memory",
          .step_index = options.step_index,
          .simulation_time_code = config.numerics.t_code_begin,
          .scale_factor = 1.0,
          .message = "startup memory accounting snapshot from owned host buffers",
          .payload = {{"persistent_total_bytes", std::to_string(startup_memory_report->totals.persistent_total_bytes)},
                      {"transient_total_bytes", std::to_string(startup_memory_report->totals.transient_total_bytes)},
                      {"unknown_external_allocations", "true"}},
      });
    }
    initializeGravityState(state, config);
    const bool restoring_from_restart = startup.restoring_from_restart;
    if (restoring_from_restart) {
      std::exception_ptr local_restart_validation_failure;
      try {
        internal::validateRestartResumeTopologyOrThrow(
            *options.restart_state_override, config, m_frozen_config,
            mpi_context);
      } catch (...) {
        local_restart_validation_failure = std::current_exception();
      }
      failure_coordinator.rethrowCollectiveFailure(
          local_restart_validation_failure, "restart topology validation");
    } else if (!startup.already_partitioned) {
      migration_balance.initializeOwnership(state);
    }
    const parallel::LocalOwnershipIdentitySummary expected_global_identity =
        migration_balance.reduceIdentity(state);
    const std::vector<std::uint64_t> expected_global_particle_ids(
        state.particle_sidecar.particle_id.begin(), state.particle_sidecar.particle_id.end());

    report.local_particle_count = static_cast<std::uint64_t>(state.particles.size());
    report.local_cell_count = static_cast<std::uint64_t>(state.cells.size());
    report.local_particle_id_sum = computeParticleIdSum(state);
    const std::uint64_t local_particle_id_square_sum = computeParticleIdSquareSum(state);
    report.local_particle_id_xor = computeParticleIdXor(state);
    const auto local_identity = parallel::summarizeLocalOwnedParticleIds(state.particle_sidecar.particle_id);
    report.local_particle_ids_unique = local_identity.local_particle_ids_unique;
    report.global_particle_count = mpi_context.allreduceSumUint64(report.local_particle_count);
    report.global_cell_count = mpi_context.allreduceSumUint64(report.local_cell_count);
    report.global_particle_id_sum = mpi_context.allreduceSumUint64(report.local_particle_id_sum);
    const std::uint64_t global_particle_id_square_sum =
        mpi_context.allreduceSumUint64(local_particle_id_square_sum);
    report.global_particle_id_xor = mpi_context.allreduceXorUint64(report.local_particle_id_xor);
    const bool all_ranks_have_unique_local_ids =
        mpi_context.allreduceSumUint64(report.local_particle_ids_unique ? 1ULL : 0ULL) ==
        static_cast<std::uint64_t>(mpi_context.worldSize());
    const parallel::LocalOwnershipIdentitySummary reduced_identity{
        .local_owned_count = report.global_particle_count,
        .local_particle_id_sum = report.global_particle_id_sum,
        .local_particle_id_square_sum = global_particle_id_square_sum,
        .local_particle_id_xor = report.global_particle_id_xor,
        .local_particle_ids_unique = all_ranks_have_unique_local_ids,
    };
    report.global_particle_partition_identity_match = parallel::partitionIdentityMatchesGeneratedSet(
        reduced_identity,
        expected_global_identity.local_owned_count,
        expected_global_identity.local_particle_id_sum,
        expected_global_identity.local_particle_id_square_sum,
        expected_global_identity.local_particle_id_xor);
    if (!report.global_particle_partition_identity_match) {
      throw std::runtime_error(
          "distributed ownership invariant failed after initial decomposition: duplicate, missing, or extra authoritative particle IDs");
    }
    profiler.recordEvent(core::RuntimeEvent{
        .event_kind = "parallel.ownership.partition_summary",
        .severity = core::RuntimeEventSeverity::kInfo,
        .subsystem = "parallel.ownership",
        .step_index = options.step_index,
        .simulation_time_code = config.numerics.t_code_begin,
        .scale_factor = 1.0,
        .message = "initial ownership-partition summary recorded",
        .payload = {{"world_size", std::to_string(report.world_size)},
                    {"world_rank", std::to_string(report.world_rank)},
                    {"local_particle_count", std::to_string(report.local_particle_count)},
                    {"global_particle_count", std::to_string(report.global_particle_count)},
                    {"local_cell_count", std::to_string(report.local_cell_count)},
                    {"global_cell_count", std::to_string(report.global_cell_count)},
                    {"local_particle_id_sum", std::to_string(report.local_particle_id_sum)},
                    {"global_particle_id_sum", std::to_string(report.global_particle_id_sum)},
                    {"local_particle_id_xor", std::to_string(report.local_particle_id_xor)},
                    {"global_particle_id_xor", std::to_string(report.global_particle_id_xor)},
                    {"local_particle_ids_unique", report.local_particle_ids_unique ? "true" : "false"},
                    {"global_particle_partition_identity_match",
                     report.global_particle_partition_identity_match ? "true" : "false"}},
    });

    core::CosmologyBackgroundConfig background_config;
    background_config.hubble_param = config.cosmology.hubble_param;
    background_config.omega_matter = config.cosmology.omega_matter;
    background_config.omega_lambda = config.cosmology.omega_lambda;
    const std::optional<core::LambdaCdmBackground> background = mode_policy.cosmological_comoving_frame
        ? std::optional<core::LambdaCdmBackground>(core::LambdaCdmBackground(background_config))
        : std::nullopt;

    const core::UnitSystem runtime_units = core::makeUnitSystem(
        config.units.length_unit,
        config.units.mass_unit,
        config.units.velocity_unit);
    RungZeroTimeState time_state = initializeRungZeroTimeState(
        config,
        options,
        state,
        runtime_units,
        background.has_value() ? &background.value() : nullptr,
        restoring_from_restart ? options.restart_state_override : nullptr);
    core::IntegratorState& integrator_state = time_state.integratorState();

    const std::filesystem::path zoom_region_path =
        config.mode.zoom_region_file.empty()
        ? std::filesystem::path{}
        : resolveConfigRelativePath(m_frozen_config, std::filesystem::path(config.mode.zoom_region_file));
    internal::ReferenceRuntimeComposition runtime_composition =
        internal::buildReferenceRuntimeComposition(
            internal::ReferenceRuntimeCompositionInputs{
                .frozen_config = m_frozen_config,
                .config = config,
                .mode_policy = mode_policy,
                .units = runtime_units,
                .services = runtime_services,
                .particle_scheduler = time_state.particleScheduler(),
                .gas_cell_scheduler = time_state.gasCellScheduler(),
                .report = report,
                .profiler = profiler,
                .pending_output = time_state.pendingOutput(),
                .options = options,
                .zoom_region_path = zoom_region_path,
                .world_rank = static_cast<std::uint32_t>(
                    std::max(mpi_context.worldRank(), 0)),
            });
    GravityRuntime& gravity_callback = *runtime_composition.gravity;
    if (restoring_from_restart) {
      const io::RestartReadResult& restart = *options.restart_state_override;
      gravity_callback.restoreDecompositionEpoch(
          restart.distributed_gravity_state.decomposition_epoch);
      if (restart.gravity_force_cache.valid) {
        gravity_callback.importRestartForceCache(restart.gravity_force_cache, state);
      } else {
        // A pre-v20 checkpoint has no cache that can safely service the next
        // KDK pre-kick.  Do not pretend its persisted PM validity bit implies
        // a live force cache in this process; force the legal bootstrap path.
        integrator_state.pm_long_range_field_valid = false;
      }
    }
    report.treepm_pm_grid = gravity_callback.pmGridSize();
    report.treepm_pm_grid_nx = gravity_callback.pmGridShape().nx;
    report.treepm_pm_grid_ny = gravity_callback.pmGridShape().ny;
    report.treepm_pm_grid_nz = gravity_callback.pmGridShape().nz;
    report.treepm_pm_grid_shape = std::to_string(report.treepm_pm_grid_nx) + "x" +
        std::to_string(report.treepm_pm_grid_ny) + "x" + std::to_string(report.treepm_pm_grid_nz);
    report.treepm_update_cadence_steps = static_cast<int>(integrator_state.pm_sync_state.cadenceSteps());
    profiler.recordEvent(core::RuntimeEvent{
        .event_kind = "gravity.treepm_setup",
        .severity = core::RuntimeEventSeverity::kInfo,
        .subsystem = "gravity.treepm",
        .step_index = integrator_state.step_index,
        .simulation_time_code = integrator_state.current_time_code,
        .scale_factor = integrator_state.current_scale_factor,
        .message = "TreePM runtime configuration initialized",
        .payload = {
            {"pm_grid", std::to_string(config.numerics.treepm_pm_grid_nx) + "x" +
                    std::to_string(config.numerics.treepm_pm_grid_ny) + "x" +
                    std::to_string(config.numerics.treepm_pm_grid_nz)},
            {"pm_assignment_scheme", treePmAssignmentSchemeName(config.numerics.treepm_assignment_scheme)},
            {"pm_window_deconvolution", config.numerics.treepm_enable_window_deconvolution ? "true" : "false"},
            {"asmth_cells", formatRuntimeDouble(config.numerics.treepm_asmth_cells)},
            {"rcut_cells", formatRuntimeDouble(config.numerics.treepm_rcut_cells)},
            {"mesh_spacing_x_mpc_comoving", formatRuntimeDouble(config.cosmology.box_size_x_mpc_comoving / static_cast<double>(config.numerics.treepm_pm_grid_nx))},
            {"mesh_spacing_y_mpc_comoving", formatRuntimeDouble(config.cosmology.box_size_y_mpc_comoving / static_cast<double>(config.numerics.treepm_pm_grid_ny))},
            {"mesh_spacing_z_mpc_comoving", formatRuntimeDouble(config.cosmology.box_size_z_mpc_comoving / static_cast<double>(config.numerics.treepm_pm_grid_nz))},
            {"split_scale_mpc_comoving", formatRuntimeDouble(config.numerics.treepm_asmth_cells *
                std::cbrt((config.cosmology.box_size_x_mpc_comoving / static_cast<double>(config.numerics.treepm_pm_grid_nx)) *
                          (config.cosmology.box_size_y_mpc_comoving / static_cast<double>(config.numerics.treepm_pm_grid_ny)) *
                          (config.cosmology.box_size_z_mpc_comoving / static_cast<double>(config.numerics.treepm_pm_grid_nz))))},
            {"cutoff_radius_mpc_comoving", formatRuntimeDouble(config.numerics.treepm_rcut_cells *
                std::cbrt((config.cosmology.box_size_x_mpc_comoving / static_cast<double>(config.numerics.treepm_pm_grid_nx)) *
                          (config.cosmology.box_size_y_mpc_comoving / static_cast<double>(config.numerics.treepm_pm_grid_ny)) *
                          (config.cosmology.box_size_z_mpc_comoving / static_cast<double>(config.numerics.treepm_pm_grid_nz))))},
            {"pm_update_cadence_steps", std::to_string(config.numerics.treepm_update_cadence_steps)},
            {"gravitational_constant_code", formatRuntimeDouble(gravity_callback.gravitationalConstantCode())},
            {"zoom_long_range_strategy",
                config.mode.zoom_long_range_strategy == core::ZoomLongRangeStrategy::kDisabled
                    ? "disabled"
                    : "global_coarse_plus_focused_highres_correction"},
            {"zoom_region_center_x_mpc_comoving", formatRuntimeDouble(config.mode.zoom_region_center_x_mpc_comoving)},
            {"zoom_region_center_y_mpc_comoving", formatRuntimeDouble(config.mode.zoom_region_center_y_mpc_comoving)},
            {"zoom_region_center_z_mpc_comoving", formatRuntimeDouble(config.mode.zoom_region_center_z_mpc_comoving)},
            {"zoom_region_radius_mpc_comoving", formatRuntimeDouble(config.mode.zoom_region_radius_mpc_comoving)},
            {"zoom_focused_pm_grid", std::to_string(config.mode.zoom_focused_pm_grid_nx) + "x" +
                std::to_string(config.mode.zoom_focused_pm_grid_ny) + "x" +
                std::to_string(config.mode.zoom_focused_pm_grid_nz)},
            {"zoom_contamination_radius_mpc_comoving", formatRuntimeDouble(config.mode.zoom_contamination_radius_mpc_comoving)},
            {"softening_policy", describeGravitySofteningPolicy(config, state)},
            {"softening_kernel", "plummer"},
            {"softening_epsilon_kpc_comoving", formatRuntimeDouble(config.numerics.gravity_softening_kpc_comoving)},
            {"pm_fft_backend", gravity::PmSolver::fftBackendName()},
        },
    });
    {
      const std::array startup_reports{
          core::collectSimulationMemoryReport(state),
          gravity_callback.memoryReport()};
      profiler.setMemoryReport(core::mergeMemoryReports(startup_reports));
      const core::MemoryReport* runtime_memory_report = profiler.memoryReport();
      if (runtime_memory_report != nullptr) {
        profiler.recordEvent(core::RuntimeEvent{
            .event_kind = "memory.runtime_startup_snapshot",
            .severity = core::RuntimeEventSeverity::kInfo,
            .subsystem = "core.memory",
            .step_index = integrator_state.step_index,
            .simulation_time_code = integrator_state.current_time_code,
            .scale_factor = integrator_state.current_scale_factor,
            .message = "startup memory accounting snapshot including gravity solver workspaces",
            .payload = {{"persistent_total_bytes", std::to_string(runtime_memory_report->totals.persistent_total_bytes)},
                        {"transient_total_bytes", std::to_string(runtime_memory_report->totals.transient_total_bytes)},
                        {"unknown_total_bytes", std::to_string(runtime_memory_report->totals.unknown_total_bytes)}},
      });
    }
  }
    HydroAmrRuntime& hydro_callback = *runtime_composition.hydro_amr;
    TimeCoordinator time_coordinator(
        runtime_services,
        time_state,
        gravity_callback,
        hydro_callback,
        std::move(runtime_composition.execution_plan),
        migration_balance);
    time_coordinator.runRungZeroSegment(
        config,
        options,
        state,
        background.has_value() ? &background.value() : nullptr,
        expected_global_particle_ids,
        report,
        profiler,
        mode_policy,
        restoring_from_restart);

    report.final_state_digest = computeStateDigest(state, integrator_state);
    report.local_particle_count = static_cast<std::uint64_t>(state.particles.size());
    report.local_cell_count = static_cast<std::uint64_t>(state.cells.size());
    report.local_particle_id_sum = computeParticleIdSum(state);
    const std::uint64_t final_local_particle_id_square_sum = computeParticleIdSquareSum(state);
    report.local_particle_id_xor = computeParticleIdXor(state);
    const auto final_local_identity = parallel::summarizeLocalOwnedParticleIds(state.particle_sidecar.particle_id);
    report.local_particle_ids_unique = final_local_identity.local_particle_ids_unique;
    report.global_particle_count = mpi_context.allreduceSumUint64(report.local_particle_count);
    report.global_cell_count = mpi_context.allreduceSumUint64(report.local_cell_count);
    report.global_particle_id_sum = mpi_context.allreduceSumUint64(report.local_particle_id_sum);
    const std::uint64_t final_global_particle_id_square_sum =
        mpi_context.allreduceSumUint64(final_local_particle_id_square_sum);
    report.global_particle_id_xor = mpi_context.allreduceXorUint64(report.local_particle_id_xor);
    const bool final_all_ranks_have_unique_local_ids =
        mpi_context.allreduceSumUint64(report.local_particle_ids_unique ? 1ULL : 0ULL) ==
        static_cast<std::uint64_t>(mpi_context.worldSize());
    const parallel::LocalOwnershipIdentitySummary final_reduced_identity{
        .local_owned_count = report.global_particle_count,
        .local_particle_id_sum = report.global_particle_id_sum,
        .local_particle_id_square_sum = final_global_particle_id_square_sum,
        .local_particle_id_xor = report.global_particle_id_xor,
        .local_particle_ids_unique = final_all_ranks_have_unique_local_ids,
    };
    report.global_particle_partition_identity_match = parallel::partitionIdentityMatchesGeneratedSet(
        final_reduced_identity,
        expected_global_identity.local_owned_count,
        expected_global_identity.local_particle_id_sum,
        expected_global_identity.local_particle_id_square_sum,
        expected_global_identity.local_particle_id_xor);
    if (!report.global_particle_partition_identity_match) {
      throw std::runtime_error(
          "distributed ownership invariant failed at run completion: duplicate, missing, or extra authoritative particle IDs; "
          "actual_count=" + std::to_string(report.global_particle_count) +
          ", expected_count=" + std::to_string(expected_global_identity.local_owned_count) +
          ", actual_sum=" + std::to_string(report.global_particle_id_sum) +
          ", expected_sum=" + std::to_string(expected_global_identity.local_particle_id_sum) +
          ", actual_square_sum=" + std::to_string(final_global_particle_id_square_sum) +
          ", expected_square_sum=" + std::to_string(expected_global_identity.local_particle_id_square_sum) +
          ", actual_xor=" + std::to_string(report.global_particle_id_xor) +
          ", expected_xor=" + std::to_string(expected_global_identity.local_particle_id_xor));
    }
    report.treepm_long_range_refresh_count = gravity_callback.longRangeRefreshCount();
    report.treepm_long_range_reuse_count = gravity_callback.longRangeReuseCount();
    report.treepm_cadence_records.assign(
        gravity_callback.cadenceRecords().begin(),
        gravity_callback.cadenceRecords().end());
    report.canonical_stage_order = repeatedCanonicalOrder(report.stage_sequence);
    profiler.recordEvent(core::RuntimeEvent{
        .event_kind = "run.complete",
        .severity = core::RuntimeEventSeverity::kInfo,
        .subsystem = "workflows.reference",
        .step_index = integrator_state.step_index,
        .simulation_time_code = integrator_state.current_time_code,
        .scale_factor = integrator_state.current_scale_factor,
        .message = "runtime workflow completed",
        .payload = {{"completed_steps", std::to_string(report.completed_steps)},
                    {"final_time_code", formatRuntimeDouble(report.final_time_code)},
                    {"final_scale_factor", formatRuntimeDouble(report.final_scale_factor)},
                    {"run_directory", report.run_directory.string()},
                    {"treepm_update_cadence_steps", std::to_string(report.treepm_update_cadence_steps)},
                    {"treepm_long_range_refresh_count", std::to_string(report.treepm_long_range_refresh_count)},
                    {"treepm_long_range_reuse_count", std::to_string(report.treepm_long_range_reuse_count)},
                    {"local_particle_count", std::to_string(report.local_particle_count)},
                    {"global_particle_count", std::to_string(report.global_particle_count)},
                    {"local_particle_id_sum", std::to_string(report.local_particle_id_sum)},
                    {"global_particle_id_sum", std::to_string(report.global_particle_id_sum)},
                    {"local_particle_id_xor", std::to_string(report.local_particle_id_xor)},
                    {"global_particle_id_xor", std::to_string(report.global_particle_id_xor)},
                    {"local_particle_ids_unique", report.local_particle_ids_unique ? "true" : "false"},
                    {"global_particle_partition_identity_match",
                     report.global_particle_partition_identity_match ? "true" : "false"},
                    {"final_state_digest", std::to_string(report.final_state_digest)}},
    });
    flushCommonArtifacts(m_frozen_config, profiler, report);
    return report;
  } catch (const std::exception& ex) {
    profiler.recordEvent(core::RuntimeEvent{
        .event_kind = "run.failure",
        .severity = core::RuntimeEventSeverity::kFatal,
        .subsystem = "workflows.reference",
        .step_index = std::nullopt,
        .simulation_time_code = std::nullopt,
        .scale_factor = std::nullopt,
        .message = "runtime workflow failed",
        .payload = {{"error", ex.what()}, {"run_directory", report.run_directory.string()}},
    });
    try {
      flushCommonArtifacts(m_frozen_config, profiler, report);
    } catch (...) {
    }
    throw;
  }
}

}  // namespace cosmosim::workflows
