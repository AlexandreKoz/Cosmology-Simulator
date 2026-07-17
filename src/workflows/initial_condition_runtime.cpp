#include "workflows/internal/initial_condition_runtime.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <utility>
#include <vector>

#include "workflows/internal/cartesian_gas_cell_layout.hpp"
#include "cosmosim/io/restart_checkpoint.hpp"

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
    const core::FrozenConfig& frozen_config) {
  const core::SimulationConfig& config = frozen_config.config;
  if (config.mode.ic_file == "generated") {
    return io::convertGeneratedIsolatedIcToState(
        config, k_generated_particle_axis);
  }
  const std::filesystem::path ic_path = resolveConfigRelativePath(
      frozen_config, std::filesystem::path(config.mode.ic_file));
  return io::readGadgetArepoHdf5Ic(ic_path, config);
}

void materializeRootHydroPatchIfMissing(
    core::SimulationState& state,
    const core::SimulationConfig& config) {
  if (state.cells.size() == 0U || state.patches.size() != 0U) {
    return;
  }
  state.patches.resize(1U);
  state.patches.patch_id[0] = 1ULL;
  state.patches.level[0] = 0;
  state.patches.first_cell[0] = 0U;
  state.patches.cell_count[0] =
      static_cast<std::uint32_t>(state.cells.size());
  state.patches.owning_rank[0] = 0U;
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
    core::SimulationState& state) {
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
  materializeRootHydroPatchIfMissing(state, frozen_config.config);
}

}  // namespace

InitialConditionRuntime::InitialConditionRuntime(
    const core::FrozenConfig& frozen_config) noexcept
    : m_frozen_config(frozen_config) {}

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
    ic_result = loadInitialConditions(m_frozen_config);
  }

  std::filesystem::path manifest_path;
  if (ic_result.report.manifest.has_value()) {
    manifest_path = run_directory / "ic_manifest.json";
    io::writeIcManifestJson(*ic_result.report.manifest, manifest_path);
  }
  finalizeStateMetadata(m_frozen_config, ic_result.state);
  return InitialConditionStartupResult{
      .state = std::move(ic_result.state),
      .import_report = std::move(ic_result.report),
      .manifest_path = std::move(manifest_path),
      .restoring_from_restart = restoring,
  };
}

}  // namespace cosmosim::workflows::internal
