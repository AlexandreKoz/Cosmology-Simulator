#include "workflows/internal/initial_condition_runtime.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <utility>
#include <vector>

#include "workflows/internal/cartesian_gas_cell_layout.hpp"
#include "cosmosim/io/restart_checkpoint.hpp"
#include "cosmosim/parallel/distributed_memory.hpp"
#include "cosmosim/core/profiling.hpp"

namespace cosmosim::workflows::internal {
namespace {

constexpr std::size_t k_generated_particle_axis = 6;

[[nodiscard]] std::filesystem::path resolveConfigRelativePath(
    const core::FrozenConfig& frozen_config,
    const std::filesystem::path& candidate) {
  if (candidate.is_absolute()) {
    return candidate;
  }
  const std::filesystem::path source_path(
      frozen_config.provenance.source_name);
  if (!source_path.empty() && source_path.has_parent_path()) {
    return source_path.parent_path() / candidate;
  }
  return candidate;
}

[[nodiscard]] io::IcReadResult loadInitialConditions(
    const core::FrozenConfig& frozen_config,
    const RuntimeServices& services) {
  const core::SimulationConfig& config = frozen_config.config;
  if (config.mode.ic_convention ==
      core::InitialConditionConvention::kGenerated) {
    return io::convertGeneratedIsolatedIcToState(
        config, k_generated_particle_axis);
  }

  io::IcImportOptions import_options{
      .require_velocities = true,
      .require_particle_ids = true,
      .allow_mass_table_fallback = true,
      .chunk_particle_count = static_cast<std::size_t>(
          config.mode.ic_chunk_particle_count),
      .manifest = nullptr,
  };
  std::optional<io::IcManifest> loaded_manifest;
  std::filesystem::path ic_path;
  if (config.mode.ic_convention ==
      core::InitialConditionConvention::kManifestV1) {
    const std::filesystem::path manifest_path = resolveConfigRelativePath(
        frozen_config,
        std::filesystem::path(config.mode.ic_manifest_file));
    loaded_manifest = io::readIcManifestJson(manifest_path);
    for (std::filesystem::path& source_path : loaded_manifest->source_files) {
      if (source_path.is_relative()) {
        source_path = (manifest_path.parent_path() / source_path).lexically_normal();
      }
    }
    io::validateIcManifest(*loaded_manifest);
    import_options.manifest = &*loaded_manifest;
    ic_path = loaded_manifest->source_files.front();
  } else {
    ic_path = resolveConfigRelativePath(
        frozen_config, std::filesystem::path(config.mode.ic_file));
  }

  if (services.mpi_context.worldSize() > 1) {
    return io::readDistributedGadgetArepoHdf5Ic(
        ic_path, config, services.mpi_context, import_options);
  }
  return io::readGadgetArepoHdf5Ic(ic_path, config, import_options);
}

void materializeRootHydroPatchIfMissing(
    core::SimulationState& state,
    const core::SimulationConfig& config,
    int world_rank) {
  if (state.cells.size() == 0U || state.patches.size() != 0U) {
    return;
  }
  state.patches.resize(1U);
  state.patches.patch_id[0] = static_cast<std::uint64_t>(world_rank) + 1ULL;
  state.patches.level[0] = 0;
  state.patches.first_cell[0] = 0U;
  state.patches.cell_count[0] =
      static_cast<std::uint32_t>(state.cells.size());
  state.patches.owning_rank[0] = static_cast<std::uint32_t>(world_rank);
  if (state.cells.patch_index.size() != state.cells.size()) {
    state.cells.patch_index.resize(state.cells.size());
  }
  std::fill(state.cells.patch_index.begin(), state.cells.patch_index.end(), 0U);

  if (state.gas_cells.size() == state.cells.size()) {
    if (state.gas_cell_identity.empty()) {
      state.refreshGasCellIdentityMapFromSidecarLanes();
    } else {
      std::vector<core::GasCellIdentityRecord> records(
          state.gas_cell_identity.records().begin(),
          state.gas_cell_identity.records().end());
      for (core::GasCellIdentityRecord& record : records) {
        record.owning_patch_id = state.patches.patch_id[0];
      }
      state.replaceGasCellIdentityRecords(std::move(records));
    }
  }

  const CartesianGasCellLayoutBuildResult layout =
      buildCartesianGasCellRowLayout(state, config);
  if (layout.ok()) {
    state.patches.origin_x_comoving[0] =
        layout.layout.spec.origin_x_comoving;
    state.patches.origin_y_comoving[0] =
        layout.layout.spec.origin_y_comoving;
    state.patches.origin_z_comoving[0] =
        layout.layout.spec.origin_z_comoving;
    state.patches.extent_x_comoving[0] =
        layout.layout.spec.cell_width_x_comoving *
        static_cast<double>(layout.layout.spec.nx);
    state.patches.extent_y_comoving[0] =
        layout.layout.spec.cell_width_y_comoving *
        static_cast<double>(layout.layout.spec.ny);
    state.patches.extent_z_comoving[0] =
        layout.layout.spec.cell_width_z_comoving *
        static_cast<double>(layout.layout.spec.nz);
    state.patches.cell_dim_x[0] =
        static_cast<std::uint16_t>(layout.layout.spec.nx);
    state.patches.cell_dim_y[0] =
        static_cast<std::uint16_t>(layout.layout.spec.ny);
    state.patches.cell_dim_z[0] =
        static_cast<std::uint16_t>(layout.layout.spec.nz);
  }
  state.bumpCellIndexGeneration();
}

void finalizeStateMetadata(
    const core::FrozenConfig& frozen_config,
    core::SimulationState& state,
    int world_rank) {
  state.metadata.run_name = frozen_config.config.output.run_name;
  state.metadata.normalized_config_hash =
      frozen_config.provenance.config_hash;
  state.metadata.normalized_config_hash_hex =
      frozen_config.provenance.config_hash_hex;
  state.metadata.snapshot_stem = frozen_config.config.output.output_stem;
  state.metadata.restart_stem = frozen_config.config.output.restart_stem;
  state.metadata.scale_factor = state.metadata.scale_factor > 0.0
      ? state.metadata.scale_factor
      : std::max(1.0, frozen_config.config.numerics.t_code_begin);
  state.rebuildSpeciesIndex();
  materializeRootHydroPatchIfMissing(
      state, frozen_config.config, world_rank);
}

}  // namespace

InitialConditionRuntime::InitialConditionRuntime(
    const core::FrozenConfig& frozen_config,
    const RuntimeServices& services) noexcept
    : m_frozen_config(frozen_config), m_services(services) {}

InitialConditionStartupResult InitialConditionRuntime::materialize(
    const ReferenceWorkflowOptions& options,
    const std::filesystem::path& run_directory) const {
  if (options.initial_state_override != nullptr &&
      options.restart_state_override != nullptr) {
    throw std::invalid_argument(
        "ReferenceWorkflowOptions cannot provide both initial_state_override and restart_state_override");
  }
  if (options.restart_state_override != nullptr &&
      !options.initial_particle_scheduler_identity_records.empty()) {
    throw std::invalid_argument(
        "ReferenceWorkflowOptions cannot override scheduler identity records while resuming a restart payload");
  }

  io::IcReadResult ic_result;
  const bool restoring = options.restart_state_override != nullptr;
  if (restoring) {
    ic_result.state = options.restart_state_override->state;
    ic_result.report.defaulted_fields.push_back(
        "restart_state_override=checkpoint_payload");
  } else if (options.initial_state_override != nullptr) {
    ic_result.state = *options.initial_state_override;
    ic_result.report.defaulted_fields.push_back(
        "initial_state_override=caller_supplied");
  } else {
    ic_result = loadInitialConditions(m_frozen_config, m_services);
  }

  std::filesystem::path manifest_path;
  if (ic_result.report.manifest.has_value()) {
    manifest_path = run_directory / "ic_manifest.json";
    if (m_services.mpi_context.isRoot()) {
      io::writeIcManifestJson(*ic_result.report.manifest, manifest_path);
    }
  }
  finalizeStateMetadata(
      m_frozen_config,
      ic_result.state,
      m_services.mpi_context.worldRank());
  const io::IcImportCounters& counters = ic_result.report.counters;
  m_services.profiler.recordEvent(core::RuntimeEvent{
      .event_kind = "io.ic_ingestion.summary",
      .severity = core::RuntimeEventSeverity::kInfo,
      .subsystem = "io.initial_conditions",
      .message = "initial-condition ingestion completed with explicit provenance and bounded staging counters",
      .payload = {{"files_assigned", std::to_string(counters.files_assigned)},
                  {"chunks_assigned", std::to_string(counters.chunks_assigned)},
                  {"bytes_read", std::to_string(counters.bytes_read)},
                  {"records_converted", std::to_string(counters.records_converted)},
                  {"records_routed", std::to_string(counters.records_routed)},
                  {"bytes_sent", std::to_string(counters.bytes_sent)},
                  {"bytes_received", std::to_string(counters.bytes_received)},
                  {"peak_staging_bytes", std::to_string(counters.peak_staging_bytes)},
                  {"final_local_particle_count", std::to_string(counters.final_local_particle_count)},
                  {"already_partitioned", ic_result.report.already_partitioned ? "true" : "false"}},
  });
  const bool already_partitioned = ic_result.report.already_partitioned;
  return InitialConditionStartupResult{
      .state = std::move(ic_result.state),
      .import_report = std::move(ic_result.report),
      .manifest_path = std::move(manifest_path),
      .restoring_from_restart = restoring,
      .already_partitioned = already_partitioned,
  };
}

}  // namespace cosmosim::workflows::internal
