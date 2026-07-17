#include "cosmosim/workflows/migration_balance_runtime.hpp"
#include "workflows/internal/amr_migration_payload.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <iterator>
#include <limits>
#include <numeric>
#include <optional>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#if COSMOSIM_ENABLE_MPI
#include <mpi.h>
#endif

namespace cosmosim::workflows::internal {
namespace {

[[nodiscard]] std::optional<std::uint32_t> parentParticleRowForGasCellRow(
    const core::SimulationState& state,
    std::uint32_t cell_row,
    const std::unordered_map<std::uint64_t, std::uint32_t>& particle_row_by_id,
    std::string_view caller) {
  const auto* record = state.gas_cell_identity.findByLocalRow(cell_row);
  if (record == nullptr) {
    throw std::runtime_error(
        std::string(caller) +
        ": gas-cell identity map is missing a dense local row");
  }
  if (!record->parent_particle_id.has_value()) {
    return std::nullopt;
  }
  const auto parent_it = particle_row_by_id.find(*record->parent_particle_id);
  return parent_it == particle_row_by_id.end()
      ? std::nullopt
      : std::optional<std::uint32_t>(parent_it->second);
}

[[nodiscard]] const core::GasCellIdentityRecord&
gasCellIdentityRecordForLocalRow(
    const core::SimulationState& state,
    std::uint32_t cell_row,
    std::string_view caller) {
  const auto* record = state.gas_cell_identity.findByLocalRow(cell_row);
  if (record == nullptr) {
    throw std::runtime_error(
        std::string(caller) +
        ": gas-cell identity map is missing a dense local row");
  }
  return *record;
}

void syncTimeBinsFromScheduler(
    const core::HierarchicalTimeBinScheduler& scheduler,
    core::SimulationState& state) {
  core::syncTimeBinMirrorsFromScheduler(
      scheduler, state, core::TimeBinMirrorDomain::kParticles);
}

struct GasCellMigrationRecord {
  std::uint64_t particle_id = 0;
  std::uint64_t patch_id = 0;
  std::uint32_t patch_local_cell_offset = 0;
  double center_x_comoving = 0.0;
  double center_y_comoving = 0.0;
  double center_z_comoving = 0.0;
  double mass_code = 0.0;
  std::uint8_t time_bin = 0;
  double velocity_x_peculiar = 0.0;
  double velocity_y_peculiar = 0.0;
  double velocity_z_peculiar = 0.0;
  double density_code = 0.0;
  double pressure_code = 0.0;
  double internal_energy_code = 0.0;
  double temperature_code = 0.0;
  double sound_speed_code = 0.0;
};

struct LocalPatchMetadata {
  std::uint64_t patch_id = 0;
  std::int32_t level = 0;
  std::uint32_t original_cell_count = 0;
  std::uint64_t parent_patch_id = 0;
  std::uint64_t morton_key = 0;
  double origin_x_comoving = 0.0;
  double origin_y_comoving = 0.0;
  double origin_z_comoving = 0.0;
  double extent_x_comoving = 0.0;
  double extent_y_comoving = 0.0;
  double extent_z_comoving = 0.0;
  std::uint16_t cell_dim_x = 0;
  std::uint16_t cell_dim_y = 0;
  std::uint16_t cell_dim_z = 0;
  std::uint32_t owning_rank = 0;

  [[nodiscard]] bool hasExplicitCartesianGeometry() const noexcept {
    return cell_dim_x != 0U || cell_dim_y != 0U || cell_dim_z != 0U ||
        extent_x_comoving != 0.0 || extent_y_comoving != 0.0 || extent_z_comoving != 0.0;
  }
};

struct LocalGasCellMigrationState {
  std::unordered_map<std::uint64_t, GasCellMigrationRecord> gas_records_by_particle_id;
  std::unordered_map<std::uint64_t, LocalPatchMetadata> patch_metadata_by_patch_id;
};

static_assert(std::is_trivially_copyable_v<GasCellMigrationRecord>);
static_assert(std::is_trivially_copyable_v<LocalPatchMetadata>);



[[nodiscard]] LocalGasCellMigrationState collectLocalGasCellRecords(
    const core::SimulationState& state,
    std::span<const std::uint32_t> gas_particle_indices) {
  core::legacyRequireParticleBoundGasCellContract(state, "collectLocalGasCellRecords legacy import path");
  if (!state.patches.isConsistent()) {
    throw std::runtime_error("collectLocalGasCellRecords requires consistent retained AMR patch metadata");
  }

  LocalGasCellMigrationState retained{};
  retained.gas_records_by_particle_id.reserve(gas_particle_indices.size());
  retained.patch_metadata_by_patch_id.reserve(state.patches.size());
  const auto gas_globals = state.particle_species_index.globalIndices(core::ParticleSpecies::kGas);
  std::vector<std::uint8_t> keep_mask(state.particles.size(), 0U);
  for (const std::uint32_t particle_index : gas_particle_indices) {
    if (particle_index >= state.particles.size()) {
      throw std::out_of_range("gas particle index out of range while collecting gas migration records");
    }
    keep_mask[particle_index] = 1U;
  }

  for (std::size_t cell_index = 0; cell_index < gas_globals.size(); ++cell_index) {
    const std::uint32_t particle_index = gas_globals[cell_index];
    if (keep_mask[particle_index] == 0U) {
      continue;
    }
    const std::uint32_t old_patch_index = state.cells.patch_index[cell_index];
    if (old_patch_index >= state.patches.size()) {
      throw std::runtime_error("retained gas cell refers to missing AMR patch metadata during ownership compaction");
    }
    const std::uint32_t first_cell = state.patches.first_cell[old_patch_index];
    const std::uint32_t cell_count = state.patches.cell_count[old_patch_index];
    if (cell_index < first_cell || cell_index >= static_cast<std::size_t>(first_cell) + cell_count) {
      throw std::runtime_error("retained gas cell lies outside its AMR patch cell range during ownership compaction");
    }
    const std::uint64_t patch_id = state.patches.patch_id[old_patch_index];
    if (patch_id == 0ULL) {
      throw std::runtime_error("retained gas cell refers to an AMR patch with an invalid stable patch ID");
    }

    const LocalPatchMetadata metadata{
        .patch_id = patch_id,
        .level = state.patches.level[old_patch_index],
        .original_cell_count = cell_count,
        .parent_patch_id = state.patches.parent_patch_id[old_patch_index],
        .morton_key = state.patches.morton_key[old_patch_index],
        .origin_x_comoving = state.patches.origin_x_comoving[old_patch_index],
        .origin_y_comoving = state.patches.origin_y_comoving[old_patch_index],
        .origin_z_comoving = state.patches.origin_z_comoving[old_patch_index],
        .extent_x_comoving = state.patches.extent_x_comoving[old_patch_index],
        .extent_y_comoving = state.patches.extent_y_comoving[old_patch_index],
        .extent_z_comoving = state.patches.extent_z_comoving[old_patch_index],
        .cell_dim_x = state.patches.cell_dim_x[old_patch_index],
        .cell_dim_y = state.patches.cell_dim_y[old_patch_index],
        .cell_dim_z = state.patches.cell_dim_z[old_patch_index],
        .owning_rank = state.patches.owning_rank[old_patch_index],
    };
    const auto [metadata_it, metadata_inserted] =
        retained.patch_metadata_by_patch_id.emplace(patch_id, metadata);
    if (!metadata_inserted &&
        (metadata_it->second.level != metadata.level ||
         metadata_it->second.original_cell_count != metadata.original_cell_count ||
         metadata_it->second.parent_patch_id != metadata.parent_patch_id ||
         metadata_it->second.morton_key != metadata.morton_key ||
         metadata_it->second.owning_rank != metadata.owning_rank)) {
      throw std::runtime_error("duplicate stable AMR patch ID has conflicting metadata during ownership compaction");
    }

    GasCellMigrationRecord record;
    record.particle_id = state.particle_sidecar.particle_id[particle_index];
    record.patch_id = patch_id;
    record.patch_local_cell_offset = static_cast<std::uint32_t>(cell_index - first_cell);
    record.center_x_comoving = state.cells.center_x_comoving[cell_index];
    record.center_y_comoving = state.cells.center_y_comoving[cell_index];
    record.center_z_comoving = state.cells.center_z_comoving[cell_index];
    record.mass_code = state.cells.mass_code[cell_index];
    record.time_bin = state.cells.time_bin[cell_index];
    record.velocity_x_peculiar = state.gas_cells.velocity_x_peculiar[cell_index];
    record.velocity_y_peculiar = state.gas_cells.velocity_y_peculiar[cell_index];
    record.velocity_z_peculiar = state.gas_cells.velocity_z_peculiar[cell_index];
    record.density_code = state.gas_cells.density_code[cell_index];
    record.pressure_code = state.gas_cells.pressure_code[cell_index];
    record.internal_energy_code = state.gas_cells.internal_energy_code[cell_index];
    record.temperature_code = state.gas_cells.temperature_code[cell_index];
    record.sound_speed_code = state.gas_cells.sound_speed_code[cell_index];
    const auto [record_it, record_inserted] =
        retained.gas_records_by_particle_id.emplace(record.particle_id, record);
    if (!record_inserted) {
      throw std::runtime_error("duplicate retained gas particle ID during ownership compaction");
    }
  }
  return retained;
}

[[nodiscard]] std::vector<std::uint64_t> gasParticleIdByOldCellIndex(const core::SimulationState& state) {
  core::legacyRequireParticleBoundGasCellContract(state, "gasParticleIdByOldCellIndex legacy import path");
  std::vector<std::uint64_t> cell_particle_ids(state.cells.size(), 0ULL);
  for (std::size_t cell_index = 0; cell_index < state.cells.size(); ++cell_index) {
    cell_particle_ids[cell_index] =
        core::parentParticleIdForGasCellRow(state, static_cast<std::uint32_t>(cell_index)).value();
  }
  return cell_particle_ids;
}

void rebuildLocalGasStateFromParticleIds(
    core::SimulationState& state,
    const LocalGasCellMigrationState& retained,
    std::span<const std::uint64_t> old_cell_particle_id) {
  const auto gas_globals = state.particle_species_index.globalIndices(core::ParticleSpecies::kGas);
  if (gas_globals.size() > static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max())) {
    throw std::overflow_error("rebuildLocalGasStateFromParticleIds exceeds uint32_t cell-index capacity");
  }

  core::CellSoa rebuilt_cells;
  core::GasCellSidecar rebuilt_gas;
  core::PatchSoa rebuilt_patches;
  rebuilt_cells.resize(gas_globals.size());
  rebuilt_gas.resize(gas_globals.size());

  struct RebuiltPatchRange {
    const LocalPatchMetadata* metadata = nullptr;
    std::uint32_t first_cell = 0;
    std::uint32_t cell_count = 0;
    std::uint32_t previous_patch_local_offset = 0;
    bool has_previous_patch_local_offset = false;
  };
  std::vector<RebuiltPatchRange> rebuilt_patch_ranges;
  rebuilt_patch_ranges.reserve(retained.patch_metadata_by_patch_id.size());
  std::unordered_map<std::uint64_t, std::uint32_t> rebuilt_patch_index_by_patch_id;
  rebuilt_patch_index_by_patch_id.reserve(retained.patch_metadata_by_patch_id.size());

  std::unordered_map<std::uint64_t, std::uint32_t> new_cell_index_by_particle_id;
  new_cell_index_by_particle_id.reserve(gas_globals.size());
  for (std::size_t cell_index = 0; cell_index < gas_globals.size(); ++cell_index) {
    const std::uint64_t particle_id = state.particle_sidecar.particle_id[gas_globals[cell_index]];
    const auto record_it = retained.gas_records_by_particle_id.find(particle_id);
    if (record_it == retained.gas_records_by_particle_id.end()) {
      throw std::runtime_error("missing gas-cell migration record for local gas particle after ownership compaction");
    }
    const GasCellMigrationRecord& record = record_it->second;
    const auto metadata_it = retained.patch_metadata_by_patch_id.find(record.patch_id);
    if (metadata_it == retained.patch_metadata_by_patch_id.end()) {
      throw std::runtime_error("retained gas cell lost its stable AMR patch metadata during ownership compaction");
    }

    std::uint32_t rebuilt_patch_index = 0;
    const auto patch_index_it = rebuilt_patch_index_by_patch_id.find(record.patch_id);
    if (patch_index_it == rebuilt_patch_index_by_patch_id.end()) {
      if (rebuilt_patch_ranges.size() > static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max())) {
        throw std::overflow_error("rebuildLocalGasStateFromParticleIds exceeds uint32_t patch-index capacity");
      }
      rebuilt_patch_index = static_cast<std::uint32_t>(rebuilt_patch_ranges.size());
      rebuilt_patch_index_by_patch_id.emplace(record.patch_id, rebuilt_patch_index);
      rebuilt_patch_ranges.push_back(RebuiltPatchRange{
          .metadata = &metadata_it->second,
          .first_cell = static_cast<std::uint32_t>(cell_index),
      });
    } else {
      rebuilt_patch_index = patch_index_it->second;
      if (rebuilt_patch_index + 1U != rebuilt_patch_ranges.size()) {
        throw std::runtime_error(
            "retained gas cells from one AMR patch are non-contiguous after particle compaction; "
            "the legacy particle-bound gas-cell contract cannot preserve this topology");
      }
    }
    RebuiltPatchRange& patch_range = rebuilt_patch_ranges[rebuilt_patch_index];
    if (patch_range.has_previous_patch_local_offset &&
        record.patch_local_cell_offset <= patch_range.previous_patch_local_offset) {
      throw std::runtime_error("retained AMR patch-local cell ordering is not strictly increasing after compaction");
    }
    patch_range.previous_patch_local_offset = record.patch_local_cell_offset;
    patch_range.has_previous_patch_local_offset = true;
    ++patch_range.cell_count;

    rebuilt_gas.gas_cell_id[cell_index] = particle_id;
    rebuilt_gas.parent_particle_id[cell_index] = particle_id;
    rebuilt_cells.center_x_comoving[cell_index] = record.center_x_comoving;
    rebuilt_cells.center_y_comoving[cell_index] = record.center_y_comoving;
    rebuilt_cells.center_z_comoving[cell_index] = record.center_z_comoving;
    rebuilt_cells.mass_code[cell_index] = record.mass_code;
    rebuilt_cells.time_bin[cell_index] = record.time_bin;
    rebuilt_cells.patch_index[cell_index] = rebuilt_patch_index;
    rebuilt_gas.velocity_x_peculiar[cell_index] = record.velocity_x_peculiar;
    rebuilt_gas.velocity_y_peculiar[cell_index] = record.velocity_y_peculiar;
    rebuilt_gas.velocity_z_peculiar[cell_index] = record.velocity_z_peculiar;
    rebuilt_gas.density_code[cell_index] = record.density_code;
    rebuilt_gas.pressure_code[cell_index] = record.pressure_code;
    rebuilt_gas.internal_energy_code[cell_index] = record.internal_energy_code;
    rebuilt_gas.temperature_code[cell_index] = record.temperature_code;
    rebuilt_gas.sound_speed_code[cell_index] = record.sound_speed_code;
    const auto [index_it, inserted] =
        new_cell_index_by_particle_id.emplace(particle_id, static_cast<std::uint32_t>(cell_index));
    if (!inserted) {
      throw std::runtime_error("duplicate local gas particle ID after ownership compaction");
    }
  }

  rebuilt_patches.resize(rebuilt_patch_ranges.size());
  for (std::size_t patch_index = 0; patch_index < rebuilt_patch_ranges.size(); ++patch_index) {
    const RebuiltPatchRange& range = rebuilt_patch_ranges[patch_index];
    if (range.metadata == nullptr || range.metadata->patch_id == 0ULL || range.cell_count == 0U) {
      throw std::runtime_error("invalid retained AMR patch range while rebuilding local gas state");
    }
    if (range.metadata->hasExplicitCartesianGeometry() &&
        range.cell_count != range.metadata->original_cell_count) {
      throw std::runtime_error(
          "ownership compaction would split a geometry-bearing AMR patch; migrate complete patches or preserve "
          "explicit sparse-cell topology before rebuilding local gas state");
    }
    const LocalPatchMetadata& metadata = *range.metadata;
    rebuilt_patches.patch_id[patch_index] = metadata.patch_id;
    rebuilt_patches.level[patch_index] = metadata.level;
    rebuilt_patches.first_cell[patch_index] = range.first_cell;
    rebuilt_patches.cell_count[patch_index] = range.cell_count;
    rebuilt_patches.parent_patch_id[patch_index] = metadata.parent_patch_id;
    rebuilt_patches.morton_key[patch_index] = metadata.morton_key;
    rebuilt_patches.origin_x_comoving[patch_index] = metadata.origin_x_comoving;
    rebuilt_patches.origin_y_comoving[patch_index] = metadata.origin_y_comoving;
    rebuilt_patches.origin_z_comoving[patch_index] = metadata.origin_z_comoving;
    rebuilt_patches.extent_x_comoving[patch_index] = metadata.extent_x_comoving;
    rebuilt_patches.extent_y_comoving[patch_index] = metadata.extent_y_comoving;
    rebuilt_patches.extent_z_comoving[patch_index] = metadata.extent_z_comoving;
    rebuilt_patches.cell_dim_x[patch_index] = metadata.cell_dim_x;
    rebuilt_patches.cell_dim_y[patch_index] = metadata.cell_dim_y;
    rebuilt_patches.cell_dim_z[patch_index] = metadata.cell_dim_z;
    rebuilt_patches.owning_rank[patch_index] = metadata.owning_rank;
  }

  state.cells = std::move(rebuilt_cells);
  state.gas_cells = std::move(rebuilt_gas);
  state.patches = std::move(rebuilt_patches);
  state.bumpCellIndexGeneration();
  if (!state.patches.isConsistent()) {
    throw std::runtime_error("rebuildLocalGasStateFromParticleIds produced inconsistent AMR patch metadata");
  }
  state.refreshGasCellIdentityMapFromSidecarLanes();
  core::legacyRequireParticleBoundGasCellContract(state, "rebuildLocalGasStateFromParticleIds legacy import path");

  const auto remap_host_cell = [&](std::uint32_t old_cell_index) {
    if (old_cell_index >= old_cell_particle_id.size()) {
      throw std::runtime_error(
          "rebuildLocalGasStateFromParticleIds: sidecar host_cell_index does not refer to an old gas cell");
    }
    const std::uint64_t host_particle_id = old_cell_particle_id[old_cell_index];
    const auto found = new_cell_index_by_particle_id.find(host_particle_id);
    if (found == new_cell_index_by_particle_id.end()) {
      throw std::runtime_error(
          "rebuildLocalGasStateFromParticleIds: sidecar host gas cell was removed during ownership compaction");
    }
    return found->second;
  };

  for (std::size_t row = 0; row < state.black_holes.size(); ++row) {
    state.black_holes.host_cell_index[row] = remap_host_cell(state.black_holes.host_cell_index[row]);
  }
  for (std::size_t row = 0; row < state.tracers.size(); ++row) {
    state.tracers.host_cell_index[row] = remap_host_cell(state.tracers.host_cell_index[row]);
  }
}

void compactStateToCurrentOwner(
    core::SimulationState& state,
    int world_rank) {
  if (world_rank < 0) {
    throw std::invalid_argument("compactStateToCurrentOwner requires non-negative world rank");
  }
  std::vector<std::uint32_t> stale_ghost_indices;
  stale_ghost_indices.reserve(state.particles.size());
  for (std::size_t particle_index = 0; particle_index < state.particles.size(); ++particle_index) {
    const bool owned_here = state.particle_sidecar.owning_rank[particle_index] == static_cast<std::uint32_t>(world_rank);
    if (!owned_here) {
      stale_ghost_indices.push_back(static_cast<std::uint32_t>(particle_index));
    }
  }
  core::ParticleMigrationCommit commit;
  commit.world_rank = world_rank;
  commit.stale_local_ghost_indices = std::move(stale_ghost_indices);
  commit.preserve_gas_cell_state = true;
  state.commitParticleMigration(commit);

  if (state.patches.size() != 0 || state.cells.size() != 0) {
    if (!state.patches.isConsistent()) {
      throw std::runtime_error("compactStateToCurrentOwner requires consistent AMR patch metadata");
    }
    std::vector<std::uint32_t> stale_patch_indices;
    stale_patch_indices.reserve(state.patches.size());
    for (std::uint32_t patch_index = 0; patch_index < state.patches.size(); ++patch_index) {
      if (state.patches.owning_rank[patch_index] != static_cast<std::uint32_t>(world_rank)) {
        stale_patch_indices.push_back(patch_index);
      }
    }
    if (!stale_patch_indices.empty()) {
      core::AmrPatchMigrationCommit patch_commit;
      patch_commit.world_rank = world_rank;
      patch_commit.stale_local_ghost_patch_indices = std::move(stale_patch_indices);
      state.commitAmrPatchMigration(patch_commit);
    }
  }
}


[[nodiscard]] parallel::DecompositionConfig makeWorkflowDecompositionConfig(
    const core::SimulationConfig& config,
    int world_size) {
  parallel::DecompositionConfig decomposition_config;
  decomposition_config.world_size = world_size;
  decomposition_config.domain_x_min_comov = 0.0;
  decomposition_config.domain_x_max_comov = config.cosmology.box_size_x_mpc_comoving;
  decomposition_config.domain_y_min_comov = 0.0;
  decomposition_config.domain_y_max_comov = config.cosmology.box_size_y_mpc_comoving;
  decomposition_config.domain_z_min_comov = 0.0;
  decomposition_config.domain_z_max_comov = config.cosmology.box_size_z_mpc_comoving;
  decomposition_config.owned_particle_weight = 0.0;
  decomposition_config.active_target_weight = 0.0;
  decomposition_config.remote_tree_interaction_weight = 0.0;
  decomposition_config.work_weight = 0.0;
  decomposition_config.memory_weight = 0.0;
  decomposition_config.component_weights = parallel::DecompositionWeightCoefficients{
      .particle_count = config.parallel.decomposition_particle_count_weight,
      .gas_cell = config.parallel.decomposition_gas_cell_weight,
      .tree_interaction = config.parallel.decomposition_tree_interaction_weight,
      .pm_mesh = config.parallel.decomposition_pm_mesh_weight,
      .amr_patch = config.parallel.decomposition_amr_patch_weight,
      .active_fraction = config.parallel.decomposition_active_fraction_weight,
      .memory_pressure = config.parallel.decomposition_memory_pressure_weight,
      .gpu_occupancy = config.parallel.decomposition_gpu_occupancy_weight,
      .generic_work = config.parallel.decomposition_generic_work_weight,
  };
  return decomposition_config;
}

[[nodiscard]] parallel::DecompositionFeedbackCoefficients makeWorkflowFeedbackCoefficients(
    const core::SimulationConfig& config) {
  return parallel::DecompositionFeedbackCoefficients{
      .measured_tree_pair = config.parallel.decomposition_measured_tree_pair_weight,
      .measured_pm_cell = config.parallel.decomposition_measured_pm_cell_weight,
      .measured_amr_cell = config.parallel.decomposition_measured_amr_cell_weight,
      .measured_hydro_face = config.parallel.decomposition_measured_hydro_face_weight,
      .measured_wall_ms = config.parallel.decomposition_measured_wall_ms_weight,
  };
}

[[nodiscard]] std::uint64_t estimateParticleMemoryBytesForDecomposition(
    const core::SimulationState& state,
    std::uint32_t species_tag) {
  std::uint64_t bytes = sizeof(double) * 7U + sizeof(std::uint64_t) * 2U + sizeof(std::uint32_t) * 3U;
  if (!state.particle_sidecar.gravity_softening_comoving.empty()) {
    bytes += sizeof(double);
  }
  if (!state.particle_sidecar.has_gravity_softening_override.empty()) {
    bytes += sizeof(std::uint8_t);
  }
  if (species_tag == static_cast<std::uint32_t>(core::ParticleSpecies::kGas)) {
    bytes += sizeof(double) * 8U + sizeof(std::uint64_t) * 2U + sizeof(std::uint32_t) * 2U;
  } else if (species_tag == static_cast<std::uint32_t>(core::ParticleSpecies::kStar)) {
    bytes += sizeof(std::uint32_t) + sizeof(double) * 13U;
  } else if (species_tag == static_cast<std::uint32_t>(core::ParticleSpecies::kBlackHole)) {
    bytes += sizeof(std::uint32_t) * 2U + sizeof(double) * 8U;
  } else if (species_tag == static_cast<std::uint32_t>(core::ParticleSpecies::kTracer)) {
    bytes += sizeof(std::uint64_t) * 2U + sizeof(std::uint32_t) * 2U + sizeof(double) * 3U;
  }
  for (const core::ModuleSidecarBlock* block_ptr : state.sidecars.blocksSortedByName()) {
    const core::ModuleSidecarBlock& block = *block_ptr;
    if (!block.particle_indexed || block.row_stride_bytes == 0U) {
      continue;
    }
    const bool species_mask_requires_row = (block.required_species_mask & (1U << species_tag)) != 0U ||
        (block.requirement.kind == core::ModuleSidecarRequirementKind::kSpeciesMask &&
         (block.requirement.species_mask & (1U << species_tag)) != 0U);
    const bool predicate_may_require_row =
        (block.requirement.kind == core::ModuleSidecarRequirementKind::kGasDensityAtLeast &&
         species_tag == static_cast<std::uint32_t>(core::ParticleSpecies::kGas)) ||
        (block.requirement.kind == core::ModuleSidecarRequirementKind::kBlackHoleAccretionAtLeast &&
         species_tag == static_cast<std::uint32_t>(core::ParticleSpecies::kBlackHole)) ||
        block.requirement.kind == core::ModuleSidecarRequirementKind::kParticleFlagMask;
    if (species_mask_requires_row || predicate_may_require_row) {
      bytes += block.row_stride_bytes;
    }
  }
  return bytes;
}

[[nodiscard]] std::vector<parallel::DecompositionItem> buildRuntimeDecompositionItems(
    const core::SimulationState& state,
    const core::SimulationConfig& config,
    int world_rank,
    std::span<const std::uint32_t> active_particle_indices) {
  std::vector<std::uint8_t> active_mask(state.particles.size(), 0U);
  for (const std::uint32_t pidx : active_particle_indices) {
    if (pidx < active_mask.size()) {
      active_mask[pidx] = 1U;
    }
  }
  std::vector<std::uint32_t> patch_cell_count(state.patches.size(), 0U);
  for (std::size_t cell_index = 0; cell_index < state.cells.size(); ++cell_index) {
    const std::uint32_t patch_index = state.cells.patch_index[cell_index];
    if (patch_index < patch_cell_count.size()) {
      ++patch_cell_count[patch_index];
    }
  }

  std::vector<parallel::DecompositionItem> items;
  items.reserve(state.particles.size() + state.patches.size());
  for (std::size_t particle_index = 0; particle_index < state.particles.size(); ++particle_index) {
    const std::uint32_t species_tag = state.particle_sidecar.species_tag[particle_index];
    parallel::DecompositionItem item;
    item.entity_id = state.particle_sidecar.particle_id[particle_index];
    item.kind = parallel::DecompositionEntityKind::kParticle;
    item.current_owner_rank = state.particle_sidecar.owning_rank.empty()
        ? world_rank
        : static_cast<int>(state.particle_sidecar.owning_rank[particle_index]);
    item.x_comov = state.particles.position_x_comoving[particle_index];
    item.y_comov = state.particles.position_y_comoving[particle_index];
    item.z_comov = state.particles.position_z_comoving[particle_index];
    item.active_target_count_recent = active_mask[particle_index] != 0U ? 1U : 0U;
    item.remote_tree_interactions_recent = 1U;
    item.memory_bytes = estimateParticleMemoryBytesForDecomposition(state, species_tag);
    double amr_patch_cost = 0.0;
    if (species_tag == static_cast<std::uint32_t>(core::ParticleSpecies::kGas) && !state.cells.patch_index.empty()) {
      std::vector<std::uint32_t> seen_patch_indices;
      for (const std::uint32_t cell_row : state.gas_cell_identity.rowsForParentParticleId(item.entity_id)) {
        if (cell_row >= state.cells.patch_index.size()) {
          continue;
        }
        const std::uint32_t patch_index = state.cells.patch_index[cell_row];
        if (patch_index >= patch_cell_count.size() ||
            std::find(seen_patch_indices.begin(), seen_patch_indices.end(), patch_index) != seen_patch_indices.end()) {
          continue;
        }
        seen_patch_indices.push_back(patch_index);
        amr_patch_cost += static_cast<double>(patch_cell_count[patch_index]) *
            (1.0 + static_cast<double>(std::max<std::int32_t>(state.patches.level[patch_index], 0)));
      }
    }
    item.work_components = parallel::DecompositionWorkComponents{
        .particle_count_cost = 1.0,
        .gas_cell_cost = species_tag == static_cast<std::uint32_t>(core::ParticleSpecies::kGas) ? 1.0 : 0.0,
        .tree_interaction_cost = 1.0,
        .pm_mesh_cost = 1.0,
        .amr_patch_cost = amr_patch_cost,
        .active_fraction_cost = static_cast<double>(item.active_target_count_recent),
        .memory_pressure_cost = static_cast<double>(item.memory_bytes),
        .gpu_occupancy_cost = 0.0,
        .generic_work_cost = 1.0,
        .has_explicit_components = true,
    };
    items.push_back(item);
  }

  for (std::size_t patch_index = 0; patch_index < state.patches.size(); ++patch_index) {
    if (state.patches.cell_count[patch_index] == 0U) {
      continue;
    }
    const std::uint32_t first_cell = state.patches.first_cell[patch_index];
    if (first_cell >= state.cells.size()) {
      continue;
    }
    parallel::DecompositionItem item;
    item.entity_id = state.patches.patch_id[patch_index];
    item.kind = parallel::DecompositionEntityKind::kAmrPatch;
    item.current_owner_rank = state.patches.owning_rank.empty()
        ? world_rank
        : static_cast<int>(state.patches.owning_rank[patch_index]);
    item.x_comov = state.cells.center_x_comoving[first_cell];
    item.y_comov = state.cells.center_y_comoving[first_cell];
    item.z_comov = state.cells.center_z_comoving[first_cell];
    item.memory_bytes = static_cast<std::uint64_t>(state.patches.cell_count[patch_index]) *
        static_cast<std::uint64_t>(sizeof(double) * 8U + sizeof(std::uint32_t) * 2U);
    item.work_components = parallel::DecompositionWorkComponents{
        .amr_patch_cost = static_cast<double>(state.patches.cell_count[patch_index]) *
            (1.0 + static_cast<double>(std::max<std::int32_t>(state.patches.level[patch_index], 0))),
        .memory_pressure_cost = static_cast<double>(item.memory_bytes),
        .generic_work_cost = static_cast<double>(state.patches.cell_count[patch_index]),
        .has_explicit_components = true,
    };
    items.push_back(item);
  }
  parallel::applyRuntimeDecompositionFeedback(items, parallel::DecompositionRuntimeMeasurements{}, makeWorkflowFeedbackCoefficients(config));
  return items;
}


void syncTimeBinsFromScheduler(
    const core::HierarchicalTimeBinScheduler& scheduler,
    core::SimulationState& state);

namespace migration_wire {

void appendBytes(std::vector<std::uint8_t>& out, const void* data, std::size_t bytes) {
  const auto* first = static_cast<const std::uint8_t*>(data);
  out.insert(out.end(), first, first + bytes);
}

template <typename T>
void appendPod(std::vector<std::uint8_t>& out, const T& value) {
  static_assert(std::is_trivially_copyable_v<T>);
  appendBytes(out, &value, sizeof(T));
}

template <typename T>
T readPod(std::span<const std::uint8_t> bytes, std::size_t& offset, std::string_view label) {
  static_assert(std::is_trivially_copyable_v<T>);
  if (offset + sizeof(T) > bytes.size()) {
    throw std::runtime_error("particle migration wire packet truncated while reading " + std::string(label));
  }
  T value{};
  std::memcpy(&value, bytes.data() + offset, sizeof(T));
  offset += sizeof(T);
  return value;
}

void appendString(std::vector<std::uint8_t>& out, std::string_view value) {
  const std::uint64_t size = static_cast<std::uint64_t>(value.size());
  appendPod(out, size);
  appendBytes(out, value.data(), value.size());
}

std::string readString(std::span<const std::uint8_t> bytes, std::size_t& offset, std::string_view label) {
  const std::uint64_t size = readPod<std::uint64_t>(bytes, offset, label);
  if (size > static_cast<std::uint64_t>(bytes.size() - offset)) {
    throw std::runtime_error("particle migration wire packet truncated while reading string " + std::string(label));
  }
  std::string value(reinterpret_cast<const char*>(bytes.data() + offset), static_cast<std::size_t>(size));
  offset += static_cast<std::size_t>(size);
  return value;
}

void appendBytePayload(std::vector<std::uint8_t>& out, std::span<const std::byte> payload) {
  const std::uint64_t size = static_cast<std::uint64_t>(payload.size());
  appendPod(out, size);
  if (!payload.empty()) {
    appendBytes(out, payload.data(), payload.size());
  }
}

std::vector<std::byte> readBytePayload(std::span<const std::uint8_t> bytes, std::size_t& offset, std::string_view label) {
  const std::uint64_t size = readPod<std::uint64_t>(bytes, offset, label);
  if (size > static_cast<std::uint64_t>(bytes.size() - offset)) {
    throw std::runtime_error("particle migration wire packet truncated while reading payload " + std::string(label));
  }
  std::vector<std::byte> payload(static_cast<std::size_t>(size));
  if (!payload.empty()) {
    std::memcpy(payload.data(), bytes.data() + offset, payload.size());
  }
  offset += static_cast<std::size_t>(size);
  return payload;
}

void appendModulePayload(std::vector<std::uint8_t>& out, const core::ModuleParticleSidecarPayload& payload) {
  appendString(out, payload.module_name);
  appendPod(out, payload.schema_version);
  appendPod(out, payload.row_stride_bytes);
  appendPod(out, payload.required_species_mask);
  appendPod(out, static_cast<std::uint32_t>(payload.requirement.kind));
  appendPod(out, payload.requirement.species_mask);
  appendPod(out, payload.requirement.particle_flags_mask);
  appendPod(out, payload.requirement.threshold_code);
  appendBytePayload(out, std::span<const std::byte>(payload.payload.data(), payload.payload.size()));
}

core::ModuleParticleSidecarPayload readModulePayload(std::span<const std::uint8_t> bytes, std::size_t& offset) {
  core::ModuleParticleSidecarPayload payload;
  payload.module_name = readString(bytes, offset, "module_name");
  payload.schema_version = readPod<std::uint32_t>(bytes, offset, "module_schema_version");
  payload.row_stride_bytes = readPod<std::uint32_t>(bytes, offset, "module_row_stride_bytes");
  payload.required_species_mask = readPod<std::uint32_t>(bytes, offset, "module_required_species_mask");
  payload.requirement.kind = static_cast<core::ModuleSidecarRequirementKind>(
      readPod<std::uint32_t>(bytes, offset, "module_requirement_kind"));
  payload.requirement.species_mask = readPod<std::uint32_t>(bytes, offset, "module_requirement_species_mask");
  payload.requirement.particle_flags_mask = readPod<std::uint32_t>(bytes, offset, "module_requirement_particle_flags_mask");
  payload.requirement.threshold_code = readPod<double>(bytes, offset, "module_requirement_threshold_code");
  payload.payload = readBytePayload(bytes, offset, "module_payload");
  return payload;
}

void appendMigrationRecord(std::vector<std::uint8_t>& out, const core::ParticleMigrationRecord& record) {
  appendPod(out, record.particle_id);
  appendPod(out, record.sfc_key);
  appendPod(out, record.species_tag);
  appendPod(out, record.particle_flags);
  appendPod(out, record.owning_rank);
  appendPod(out, record.position_x_comoving);
  appendPod(out, record.position_y_comoving);
  appendPod(out, record.position_z_comoving);
  appendPod(out, record.velocity_x_peculiar);
  appendPod(out, record.velocity_y_peculiar);
  appendPod(out, record.velocity_z_peculiar);
  appendPod(out, record.mass_code);
  appendPod(out, record.time_bin);
  appendPod(out, static_cast<std::uint8_t>(record.has_scheduler_fields ? 1U : 0U));
  appendPod(out, record.scheduler_fields);
  appendPod(out, record.last_drift_time_code);
  appendPod(out, record.last_drift_scale_factor);
  appendPod(out, static_cast<std::uint8_t>(record.has_gravity_softening_value ? 1U : 0U));
  appendPod(out, static_cast<std::uint8_t>(record.has_gravity_softening_override ? 1U : 0U));
  appendPod(out, record.gravity_softening_comoving);
  appendPod(out, static_cast<std::uint8_t>(record.has_gas_cell_fields ? 1U : 0U));
  appendPod(out, record.gas_cell_fields);
  appendPod(out, static_cast<std::uint8_t>(record.has_star_fields ? 1U : 0U));
  appendPod(out, record.star_fields);
  appendPod(out, static_cast<std::uint8_t>(record.has_black_hole_fields ? 1U : 0U));
  appendPod(out, record.black_hole_fields);
  appendPod(out, static_cast<std::uint8_t>(record.has_tracer_fields ? 1U : 0U));
  appendPod(out, record.tracer_fields);
  const std::uint64_t module_count = static_cast<std::uint64_t>(record.module_sidecar_payloads.size());
  appendPod(out, module_count);
  for (const auto& payload : record.module_sidecar_payloads) {
    appendModulePayload(out, payload);
  }
}

core::ParticleMigrationRecord readMigrationRecord(std::span<const std::uint8_t> bytes, std::size_t& offset) {
  core::ParticleMigrationRecord record;
  record.particle_id = readPod<std::uint64_t>(bytes, offset, "particle_id");
  record.sfc_key = readPod<std::uint64_t>(bytes, offset, "sfc_key");
  record.species_tag = readPod<std::uint32_t>(bytes, offset, "species_tag");
  record.particle_flags = readPod<std::uint32_t>(bytes, offset, "particle_flags");
  record.owning_rank = readPod<std::uint32_t>(bytes, offset, "owning_rank");
  record.position_x_comoving = readPod<double>(bytes, offset, "position_x_comoving");
  record.position_y_comoving = readPod<double>(bytes, offset, "position_y_comoving");
  record.position_z_comoving = readPod<double>(bytes, offset, "position_z_comoving");
  record.velocity_x_peculiar = readPod<double>(bytes, offset, "velocity_x_peculiar");
  record.velocity_y_peculiar = readPod<double>(bytes, offset, "velocity_y_peculiar");
  record.velocity_z_peculiar = readPod<double>(bytes, offset, "velocity_z_peculiar");
  record.mass_code = readPod<double>(bytes, offset, "mass_code");
  record.time_bin = readPod<std::uint8_t>(bytes, offset, "time_bin");
  record.has_scheduler_fields = readPod<std::uint8_t>(bytes, offset, "has_scheduler_fields") != 0U;
  record.scheduler_fields = readPod<core::SchedulerMigrationFields>(bytes, offset, "scheduler_fields");
  record.last_drift_time_code = readPod<double>(bytes, offset, "last_drift_time_code");
  record.last_drift_scale_factor = readPod<double>(bytes, offset, "last_drift_scale_factor");
  record.has_gravity_softening_value = readPod<std::uint8_t>(bytes, offset, "has_gravity_softening_value") != 0U;
  record.has_gravity_softening_override = readPod<std::uint8_t>(bytes, offset, "has_gravity_softening_override") != 0U;
  record.gravity_softening_comoving = readPod<double>(bytes, offset, "gravity_softening_comoving");
  record.has_gas_cell_fields = readPod<std::uint8_t>(bytes, offset, "has_gas_cell_fields") != 0U;
  record.gas_cell_fields = readPod<core::GasCellMigrationFields>(bytes, offset, "gas_cell_fields");
  record.has_star_fields = readPod<std::uint8_t>(bytes, offset, "has_star_fields") != 0U;
  record.star_fields = readPod<core::StarParticleMigrationFields>(bytes, offset, "star_fields");
  record.has_black_hole_fields = readPod<std::uint8_t>(bytes, offset, "has_black_hole_fields") != 0U;
  record.black_hole_fields = readPod<core::BlackHoleParticleMigrationFields>(bytes, offset, "black_hole_fields");
  record.has_tracer_fields = readPod<std::uint8_t>(bytes, offset, "has_tracer_fields") != 0U;
  record.tracer_fields = readPod<core::TracerParticleMigrationFields>(bytes, offset, "tracer_fields");
  const std::uint64_t module_count = readPod<std::uint64_t>(bytes, offset, "module_count");
  if (module_count > 1'000'000ULL) {
    throw std::runtime_error("particle migration wire packet has unreasonable module payload count");
  }
  record.module_sidecar_payloads.reserve(static_cast<std::size_t>(module_count));
  for (std::uint64_t i = 0; i < module_count; ++i) {
    record.module_sidecar_payloads.push_back(readModulePayload(bytes, offset));
  }
  return record;
}

std::vector<std::uint8_t> serializeMigrationRecords(std::span<const core::ParticleMigrationRecord> records) {
  std::vector<std::uint8_t> bytes;
  appendPod(bytes, static_cast<std::uint64_t>(records.size()));
  for (const core::ParticleMigrationRecord& record : records) {
    appendMigrationRecord(bytes, record);
  }
  return bytes;
}

std::vector<core::ParticleMigrationRecord> deserializeMigrationRecords(std::span<const std::uint8_t> bytes) {
  std::size_t offset = 0;
  const std::uint64_t count = readPod<std::uint64_t>(bytes, offset, "migration_record_count");
  if (count > 100'000'000ULL) {
    throw std::runtime_error("particle migration wire packet has unreasonable record count");
  }
  std::vector<core::ParticleMigrationRecord> records;
  records.reserve(static_cast<std::size_t>(count));
  for (std::uint64_t i = 0; i < count; ++i) {
    records.push_back(readMigrationRecord(bytes, offset));
  }
  if (offset != bytes.size()) {
    throw std::runtime_error("particle migration wire packet has trailing bytes");
  }
  return records;
}

void appendGasCellMigrationRecord(std::vector<std::uint8_t>& out, const core::GasCellMigrationRecord& record) {
  appendPod(out, record.owning_rank);
  appendPod(out, record.fields);
}

core::GasCellMigrationRecord readGasCellMigrationRecord(std::span<const std::uint8_t> bytes, std::size_t& offset) {
  core::GasCellMigrationRecord record;
  record.owning_rank = readPod<std::uint32_t>(bytes, offset, "gas_cell_owning_rank");
  record.fields = readPod<core::GasCellMigrationFields>(bytes, offset, "gas_cell_fields");
  return record;
}

void appendAmrPatchMigrationRecord(std::vector<std::uint8_t>& out, const core::AmrPatchMigrationRecord& record) {
  appendPod(out, record.patch);
  appendPod(out, static_cast<std::uint64_t>(record.gas_cell_records.size()));
  for (const core::GasCellMigrationRecord& gas_record : record.gas_cell_records) {
    appendGasCellMigrationRecord(out, gas_record);
  }
  appendPod(out, static_cast<std::uint64_t>(record.gas_cell_scheduler_records.size()));
  for (const core::GasCellSchedulerMigrationRecord& scheduler_record : record.gas_cell_scheduler_records) {
    appendPod(out, scheduler_record);
  }
}

core::AmrPatchMigrationRecord readAmrPatchMigrationRecord(std::span<const std::uint8_t> bytes, std::size_t& offset) {
  core::AmrPatchMigrationRecord record;
  record.patch = readPod<core::AmrPatchMigrationFields>(bytes, offset, "amr_patch_fields");
  const std::uint64_t gas_count = readPod<std::uint64_t>(bytes, offset, "amr_patch_gas_cell_count");
  if (gas_count > 100'000'000ULL) {
    throw std::runtime_error("AMR patch migration wire packet has unreasonable gas-cell count");
  }
  record.gas_cell_records.reserve(static_cast<std::size_t>(gas_count));
  for (std::uint64_t i = 0; i < gas_count; ++i) {
    record.gas_cell_records.push_back(readGasCellMigrationRecord(bytes, offset));
  }
  const std::uint64_t scheduler_count = readPod<std::uint64_t>(bytes, offset, "amr_patch_scheduler_record_count");
  if (scheduler_count != gas_count) {
    throw std::runtime_error("AMR patch migration wire packet scheduler record count does not match gas-cell count");
  }
  record.gas_cell_scheduler_records.reserve(static_cast<std::size_t>(scheduler_count));
  for (std::uint64_t i = 0; i < scheduler_count; ++i) {
    record.gas_cell_scheduler_records.push_back(
        readPod<core::GasCellSchedulerMigrationRecord>(bytes, offset, "gas_cell_scheduler_record"));
  }
  return record;
}

std::vector<std::uint8_t> serializeAmrPatchMigrationRecords(
    std::span<const core::AmrPatchMigrationRecord> records) {
  std::vector<std::uint8_t> bytes;
  appendPod(bytes, static_cast<std::uint64_t>(records.size()));
  for (const core::AmrPatchMigrationRecord& record : records) {
    appendAmrPatchMigrationRecord(bytes, record);
  }
  return bytes;
}

std::vector<core::AmrPatchMigrationRecord> deserializeAmrPatchMigrationRecords(
    std::span<const std::uint8_t> bytes) {
  std::size_t offset = 0;
  const std::uint64_t count = readPod<std::uint64_t>(bytes, offset, "amr_patch_migration_record_count");
  if (count > 10'000'000ULL) {
    throw std::runtime_error("AMR patch migration wire packet has unreasonable record count");
  }
  std::vector<core::AmrPatchMigrationRecord> records;
  records.reserve(static_cast<std::size_t>(count));
  for (std::uint64_t i = 0; i < count; ++i) {
    records.push_back(readAmrPatchMigrationRecord(bytes, offset));
  }
  if (offset != bytes.size()) {
    throw std::runtime_error("AMR patch migration wire packet has trailing bytes");
  }
  return records;
}

}  // namespace migration_wire

std::vector<core::ParticleMigrationRecord> exchangeRuntimeParticleMigrationRecords(
    const parallel::MpiContext& mpi_context,
    const std::vector<std::vector<core::ParticleMigrationRecord>>& records_by_rank) {
  if (records_by_rank.size() != static_cast<std::size_t>(mpi_context.worldSize())) {
    throw std::invalid_argument("particle migration exchange records_by_rank size must match world size");
  }
  if (mpi_context.worldSize() == 1) {
    return records_by_rank.empty() ? std::vector<core::ParticleMigrationRecord>{} : records_by_rank[0];
  }
  if (!mpi_context.isEnabled()) {
    throw std::runtime_error("runtime particle migration execution requires MPI for world_size > 1");
  }
#if COSMOSIM_ENABLE_MPI
  std::vector<std::vector<std::uint8_t>> send_payloads(records_by_rank.size());
  std::vector<int> send_counts(records_by_rank.size(), 0);
  for (std::size_t rank = 0; rank < records_by_rank.size(); ++rank) {
    send_payloads[rank] = migration_wire::serializeMigrationRecords(records_by_rank[rank]);
    if (send_payloads[rank].size() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
      throw std::overflow_error("particle migration payload exceeds MPI int byte-count limit");
    }
    send_counts[rank] = static_cast<int>(send_payloads[rank].size());
  }
  std::vector<int> recv_counts(send_counts.size(), 0);
  MPI_Alltoall(send_counts.data(), 1, MPI_INT, recv_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);

  std::vector<int> send_offsets(send_counts.size(), 0);
  std::vector<int> recv_offsets(recv_counts.size(), 0);
  int total_send = 0;
  int total_recv = 0;
  for (std::size_t rank = 0; rank < send_counts.size(); ++rank) {
    send_offsets[rank] = total_send;
    recv_offsets[rank] = total_recv;
    total_send += send_counts[rank];
    total_recv += recv_counts[rank];
  }
  std::vector<std::uint8_t> send_bytes(static_cast<std::size_t>(total_send));
  for (std::size_t rank = 0; rank < send_payloads.size(); ++rank) {
    if (!send_payloads[rank].empty()) {
      std::memcpy(send_bytes.data() + send_offsets[rank], send_payloads[rank].data(), send_payloads[rank].size());
    }
  }
  std::vector<std::uint8_t> recv_bytes(static_cast<std::size_t>(total_recv));
  MPI_Alltoallv(
      send_bytes.data(),
      send_counts.data(),
      send_offsets.data(),
      MPI_BYTE,
      recv_bytes.data(),
      recv_counts.data(),
      recv_offsets.data(),
      MPI_BYTE,
      MPI_COMM_WORLD);

  std::vector<core::ParticleMigrationRecord> inbound;
  for (std::size_t rank = 0; rank < recv_counts.size(); ++rank) {
    const auto offset = static_cast<std::size_t>(recv_offsets[rank]);
    const auto count = static_cast<std::size_t>(recv_counts[rank]);
    auto decoded = migration_wire::deserializeMigrationRecords(
        std::span<const std::uint8_t>(recv_bytes.data() + offset, count));
    inbound.insert(inbound.end(), std::make_move_iterator(decoded.begin()), std::make_move_iterator(decoded.end()));
  }
  return inbound;
#else
  throw std::runtime_error("runtime particle migration execution requires an MPI-enabled build");
#endif
}

std::vector<core::AmrPatchMigrationRecord> exchangeRuntimeAmrPatchMigrationRecords(
    const parallel::MpiContext& mpi_context,
    const std::vector<std::vector<core::AmrPatchMigrationRecord>>& records_by_rank) {
  if (records_by_rank.size() != static_cast<std::size_t>(mpi_context.worldSize())) {
    throw std::invalid_argument("AMR patch migration exchange records_by_rank size must match world size");
  }
  if (mpi_context.worldSize() == 1) {
    return records_by_rank.empty() ? std::vector<core::AmrPatchMigrationRecord>{} : records_by_rank[0];
  }
  if (!mpi_context.isEnabled()) {
    throw std::runtime_error("runtime AMR patch migration execution requires MPI for world_size > 1");
  }
#if COSMOSIM_ENABLE_MPI
  std::vector<std::vector<std::uint8_t>> send_payloads(records_by_rank.size());
  std::vector<int> send_counts(records_by_rank.size(), 0);
  for (std::size_t rank = 0; rank < records_by_rank.size(); ++rank) {
    send_payloads[rank] = migration_wire::serializeAmrPatchMigrationRecords(records_by_rank[rank]);
    if (send_payloads[rank].size() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
      throw std::overflow_error("AMR patch migration payload exceeds MPI int byte-count limit");
    }
    send_counts[rank] = static_cast<int>(send_payloads[rank].size());
  }
  std::vector<int> recv_counts(send_counts.size(), 0);
  MPI_Alltoall(send_counts.data(), 1, MPI_INT, recv_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);

  std::vector<int> send_offsets(send_counts.size(), 0);
  std::vector<int> recv_offsets(recv_counts.size(), 0);
  int total_send = 0;
  int total_recv = 0;
  for (std::size_t rank = 0; rank < send_counts.size(); ++rank) {
    send_offsets[rank] = total_send;
    recv_offsets[rank] = total_recv;
    total_send += send_counts[rank];
    total_recv += recv_counts[rank];
  }
  std::vector<std::uint8_t> send_bytes(static_cast<std::size_t>(total_send));
  for (std::size_t rank = 0; rank < send_payloads.size(); ++rank) {
    if (!send_payloads[rank].empty()) {
      std::memcpy(send_bytes.data() + send_offsets[rank], send_payloads[rank].data(), send_payloads[rank].size());
    }
  }
  std::vector<std::uint8_t> recv_bytes(static_cast<std::size_t>(total_recv));
  MPI_Alltoallv(
      send_bytes.data(),
      send_counts.data(),
      send_offsets.data(),
      MPI_BYTE,
      recv_bytes.data(),
      recv_counts.data(),
      recv_offsets.data(),
      MPI_BYTE,
      MPI_COMM_WORLD);

  std::vector<core::AmrPatchMigrationRecord> inbound;
  for (std::size_t rank = 0; rank < recv_counts.size(); ++rank) {
    const auto offset = static_cast<std::size_t>(recv_offsets[rank]);
    const auto count = static_cast<std::size_t>(recv_counts[rank]);
    auto decoded = migration_wire::deserializeAmrPatchMigrationRecords(
        std::span<const std::uint8_t>(recv_bytes.data() + offset, count));
    inbound.insert(inbound.end(), std::make_move_iterator(decoded.begin()), std::make_move_iterator(decoded.end()));
  }
  return inbound;
#else
  throw std::runtime_error("runtime AMR patch migration execution requires an MPI-enabled build");
#endif
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

void requireGlobalOwnedParticlePartitionIdentity(
    const core::SimulationState& state,
    const parallel::MpiContext& mpi_context,
    std::span<const std::uint64_t> expected_global_particle_ids,
    std::string_view caller) {
  std::vector<std::uint64_t> local_owned_ids;
  local_owned_ids.reserve(state.particles.size());
  const std::uint32_t world_rank = static_cast<std::uint32_t>(mpi_context.worldRank());
  for (std::size_t particle_index = 0; particle_index < state.particles.size(); ++particle_index) {
    if (state.particle_sidecar.owning_rank[particle_index] == world_rank) {
      local_owned_ids.push_back(state.particle_sidecar.particle_id[particle_index]);
    }
  }
  const parallel::ExactOwnershipPartitionReport report =
      parallel::validateExactGlobalOwnershipPartition(mpi_context, local_owned_ids, expected_global_particle_ids);
  if (!report.valid()) {
    throw std::runtime_error(std::string(caller) +
        ": exact distributed authoritative ownership table has duplicate=" +
        std::to_string(report.duplicate_particle_ids.size()) +
        ", missing=" + std::to_string(report.missing_expected_particle_ids.size()) +
        ", extra=" + std::to_string(report.extra_particle_ids.size()) + " particle IDs");
  }
}

void recordRuntimeRebalanceDecision(
    core::ProfilerSession* profiler,
    const parallel::RuntimeRebalancePlan& rebalance,
    std::uint64_t step_index) {
  if (profiler == nullptr) {
    return;
  }
  const auto join_u64 = [](std::span<const std::uint64_t> values) {
    std::ostringstream stream;
    for (std::size_t i = 0; i < values.size(); ++i) {
      if (i != 0U) {
        stream << ',';
      }
      stream << values[i];
    }
    return stream.str();
  };
  profiler->recordEvent(core::RuntimeEvent{
      .event_kind = "parallel.decomposition.runtime_rebalance",
      .severity = rebalance.should_rebalance ? core::RuntimeEventSeverity::kWarning : core::RuntimeEventSeverity::kInfo,
      .subsystem = "parallel.domain_decomposition",
      .step_index = step_index,
      .message = rebalance.should_rebalance
          ? "runtime work-weighted decomposition produced migration intent"
          : "runtime work-weighted decomposition remained below rebalance threshold",
      .payload = {{"reason", rebalance.reason},
                  {"should_rebalance", rebalance.should_rebalance ? "true" : "false"},
                  {"particle_migration_count", std::to_string(rebalance.particle_migrations.size())},
                  {"amr_patch_ownership_update_count", std::to_string(rebalance.amr_patch_ownership_updates.size())},
                  {"migrated_load_fraction", std::to_string(rebalance.migrated_load_fraction)},
                  {"used_distributed_sfc_cuts", rebalance.used_distributed_sfc_cuts ? "true" : "false"},
                  {"exact_debug_audit_enabled", rebalance.exact_debug_audit_enabled ? "true" : "false"},
                  {"local_entities_considered", std::to_string(rebalance.local_entities_considered)},
                  {"global_entities_considered", std::to_string(rebalance.global_entities_considered)},
                  {"local_entities_moved", std::to_string(rebalance.local_entities_moved)},
                  {"global_entities_moved", std::to_string(rebalance.global_entities_moved)},
                  {"local_bytes_moved", std::to_string(rebalance.local_bytes_moved)},
                  {"global_bytes_moved", std::to_string(rebalance.global_bytes_moved)},
                  {"global_control_bytes", std::to_string(rebalance.global_control_bytes)},
                  {"peak_temporary_bytes", std::to_string(rebalance.peak_temporary_bytes)},
                  {"cut_displacement_fraction", std::to_string(rebalance.cut_displacement_fraction)},
                  {"sfc_cut_count", std::to_string(rebalance.sfc_cut_keys.size())},
                  {"sfc_cut_keys", join_u64(rebalance.sfc_cut_keys)},
                  {"sfc_cut_entity_ids", join_u64(rebalance.sfc_cut_entity_ids)},
                  {"current_weighted_imbalance_ratio", std::to_string(rebalance.current_metrics.weighted_imbalance_ratio)},
                  {"target_weighted_imbalance_ratio", std::to_string(rebalance.target_decomposition.metrics.weighted_imbalance_ratio)},
                  {"current_memory_imbalance_ratio", std::to_string(rebalance.current_metrics.memory_imbalance_ratio)},
                  {"target_memory_imbalance_ratio", std::to_string(rebalance.target_decomposition.metrics.memory_imbalance_ratio)}}});
}

void exchangeAndValidateAmrPatchPayloads(
    core::SimulationState& state,
    const parallel::MpiContext& mpi_context,
    int world_rank,
    std::uint64_t step_index,
    core::ProfilerSession* profiler) {
  const std::vector<parallel::AmrPatchPayloadRecord> local_records =
      buildMigrationAmrPatchPayloadRecords(state, world_rank);
  const std::vector<parallel::AmrPatchCellPayloadRecord> local_cell_records =
      buildMigrationAmrPatchCellPayloadRecords(state, world_rank);

  std::unordered_map<std::uint64_t, std::uint32_t> expected_cell_count_by_patch_id;
  std::unordered_map<std::uint64_t, std::uint32_t> observed_cell_count_by_patch_id;
  expected_cell_count_by_patch_id.reserve(local_records.size());
  observed_cell_count_by_patch_id.reserve(local_records.size());
  for (const parallel::AmrPatchPayloadRecord& record : local_records) {
    const auto [it, inserted] = expected_cell_count_by_patch_id.emplace(record.patch_id, record.cell_count);
    if (!inserted) {
      throw std::runtime_error("AMR owner-local validation detected duplicate authoritative patch payloads");
    }
    if (record.owner_rank != world_rank) {
      throw std::runtime_error("AMR owner-local validation found stale patch owner metadata");
    }
  }
  for (const parallel::AmrPatchCellPayloadRecord& record : local_cell_records) {
    const auto owner_it = expected_cell_count_by_patch_id.find(record.patch_id);
    if (owner_it == expected_cell_count_by_patch_id.end()) {
      throw std::runtime_error("AMR owner-local validation found a cell for an unknown local patch");
    }
    if (record.owner_rank != world_rank) {
      throw std::runtime_error("AMR owner-local validation found stale cell owner metadata");
    }
    if (record.local_cell_offset >= owner_it->second) {
      throw std::runtime_error("AMR owner-local validation found a cell offset beyond patch bounds");
    }
    ++observed_cell_count_by_patch_id[record.patch_id];
  }
  for (const auto& [patch_id, expected_count] : expected_cell_count_by_patch_id) {
    if (observed_cell_count_by_patch_id[patch_id] != expected_count) {
      throw std::runtime_error("AMR owner-local validation patch cell coverage mismatch for patch_id=" +
                               std::to_string(patch_id));
    }
  }

  parallel::DirectedAmrExchangeDiagnostics diagnostics;
  if (mpi_context.isEnabled() && mpi_context.worldSize() > 1) {
    const parallel::DirectedAmrPatchPayloadExchange directed = parallel::executeBlockingDirectedAmrPatchPayloadExchange(
        mpi_context,
        local_records,
        local_cell_records,
        step_index);
    diagnostics = directed.diagnostics;

    std::unordered_map<std::uint64_t, parallel::AmrPatchPayloadRecord> remote_patch_by_id;
    remote_patch_by_id.reserve(directed.patch_payloads_received.size());
    for (const parallel::AmrPatchPayloadRecord& record : directed.patch_payloads_received) {
      if (record.owner_rank == world_rank) {
        throw std::runtime_error("directed AMR validation received local authoritative patch as a remote ghost");
      }
      const auto [it, inserted] = remote_patch_by_id.emplace(record.patch_id, record);
      if (!inserted) {
        throw std::runtime_error("directed AMR validation received duplicate remote patch metadata");
      }
    }
    std::unordered_map<std::uint64_t, std::uint32_t> remote_cell_count_by_patch_id;
    remote_cell_count_by_patch_id.reserve(remote_patch_by_id.size());
    for (const parallel::AmrPatchCellPayloadRecord& record : directed.patch_cell_payloads_received) {
      const auto patch_it = remote_patch_by_id.find(record.patch_id);
      if (patch_it == remote_patch_by_id.end()) {
        throw std::runtime_error("directed AMR validation received remote cell payload without patch metadata");
      }
      if (record.owner_rank != patch_it->second.owner_rank) {
        throw std::runtime_error("directed AMR validation remote cell owner does not match patch owner");
      }
      if (record.local_cell_offset >= patch_it->second.cell_count) {
        throw std::runtime_error("directed AMR validation remote cell offset exceeds patch extent");
      }
      ++remote_cell_count_by_patch_id[record.patch_id];
    }
    for (const auto& [patch_id, record] : remote_patch_by_id) {
      if (remote_cell_count_by_patch_id[patch_id] != record.cell_count) {
        throw std::runtime_error("directed AMR validation remote patch cell coverage mismatch for patch_id=" +
                                 std::to_string(patch_id));
      }
    }
  }

  if (profiler != nullptr) {
    profiler->recordEvent(core::RuntimeEvent{
        .event_kind = "amr.directed_patch_payload_exchange",
        .severity = core::RuntimeEventSeverity::kInfo,
        .subsystem = "amr.patch_exchange",
        .step_index = step_index,
        .message = "validated owner-local AMR state and interface-scoped directed AMR payload coverage",
        .payload = {{"local_patch_payloads", std::to_string(local_records.size())},
                    {"local_patch_cell_payloads", std::to_string(local_cell_records.size())},
                    {"candidate_peer_count", std::to_string(diagnostics.candidate_peer_count)},
                    {"neighbor_peer_count", std::to_string(diagnostics.neighbor_peer_count)},
                    {"directed_patch_descriptor_records_sent", std::to_string(diagnostics.directed_patch_descriptor_records_sent)},
                    {"directed_patch_descriptor_records_received", std::to_string(diagnostics.directed_patch_descriptor_records_received)},
                    {"directed_patch_cell_records_sent", std::to_string(diagnostics.directed_patch_cell_records_sent)},
                    {"directed_patch_cell_records_received", std::to_string(diagnostics.directed_patch_cell_records_received)},
                    {"control_plane_bytes", std::to_string(diagnostics.control_plane_bytes)},
                    {"patch_descriptor_bytes", std::to_string(diagnostics.patch_descriptor_bytes)},
                    {"patch_cell_payload_bytes", std::to_string(diagnostics.patch_cell_payload_bytes)},
                    {"remote_patch_ghost_count", std::to_string(diagnostics.remote_patch_ghost_count)},
                    {"remote_interface_count", std::to_string(diagnostics.remote_interface_count)}}});
  }
}

[[nodiscard]] bool applyMeasuredRuntimeRebalancePlan(
    core::SimulationState& state,
    core::HierarchicalTimeBinScheduler& scheduler,
    core::HierarchicalTimeBinScheduler& gas_cell_scheduler,
    const core::SimulationConfig& config,
    const parallel::MpiContext& mpi_context,
    int world_rank,
    const parallel::DecompositionRuntimeMeasurements& measurements,
    std::span<const std::uint32_t> active_particle_indices,
    std::span<const std::uint64_t> expected_global_particle_ids,
    core::ProfilerSession* profiler,
    std::uint64_t step_index) {
  if (!config.parallel.decomposition_runtime_rebalance_enabled || mpi_context.worldSize() <= 1) {
    return false;
  }
  if (world_rank < 0 || world_rank >= mpi_context.worldSize()) {
    throw std::invalid_argument("runtime rebalance world_rank is outside MPI world");
  }
  auto local_items = buildRuntimeDecompositionItems(state, config, world_rank, active_particle_indices);
  parallel::applyRuntimeDecompositionFeedback(local_items, measurements, makeWorkflowFeedbackCoefficients(config));
  parallel::RuntimeRebalanceConfig rebalance_config{
      .world_size = mpi_context.worldSize(),
      .imbalance_trigger_ratio = config.parallel.decomposition_rebalance_imbalance_trigger,
      .memory_trigger_ratio = config.parallel.decomposition_rebalance_memory_trigger,
      .max_migrated_load_fraction = config.parallel.decomposition_rebalance_max_migrated_load_fraction,
      .allow_particle_migration = true,
      .allow_amr_patch_reassignment = true,
  };
  auto rebalance = parallel::buildDistributedRuntimeRebalancePlan(
      mpi_context,
      local_items,
      makeWorkflowDecompositionConfig(config, mpi_context.worldSize()),
      rebalance_config);
  rebalance.exact_debug_audit_enabled = config.parallel.decomposition_debug_exact_ownership_audit;
  if (!rebalance.should_rebalance) {
    exchangeAndValidateAmrPatchPayloads(state, mpi_context, world_rank, step_index, profiler);
    recordRuntimeRebalanceDecision(profiler, rebalance, step_index);
    return false;
  }

  const std::uint64_t particle_index_generation_before = state.particleIndexGeneration();

  std::unordered_map<std::uint64_t, std::uint32_t> local_index_by_particle_id;
  local_index_by_particle_id.reserve(state.particles.size());
  for (std::size_t particle_index = 0; particle_index < state.particles.size(); ++particle_index) {
    local_index_by_particle_id.emplace(
        state.particle_sidecar.particle_id[particle_index],
        static_cast<std::uint32_t>(particle_index));
  }

  std::unordered_map<std::uint32_t, int> outbound_target_by_local_index;
  outbound_target_by_local_index.reserve(rebalance.particle_migrations.size());
  std::unordered_map<std::uint32_t, int> outbound_patch_target_by_local_index;
  outbound_patch_target_by_local_index.reserve(rebalance.amr_patch_ownership_updates.size());
  const auto add_particle_migration = [&](std::uint32_t local_index, int new_owner_rank, std::string_view source_label) {
    if (local_index >= state.particles.size()) {
      throw std::out_of_range("runtime rebalance attempted to migrate a particle index outside SimulationState");
    }
    if (new_owner_rank < 0 || new_owner_rank >= mpi_context.worldSize()) {
      throw std::invalid_argument("runtime rebalance produced particle migration target outside MPI world");
    }
    if (new_owner_rank == world_rank) {
      return;
    }
    if (state.particle_sidecar.owning_rank[local_index] != static_cast<std::uint32_t>(world_rank)) {
      throw std::runtime_error("runtime rebalance attempted to migrate a non-authoritative local particle via " + std::string(source_label));
    }
    const auto [it, inserted] = outbound_target_by_local_index.emplace(local_index, new_owner_rank);
    if (!inserted && it->second != new_owner_rank) {
      throw std::runtime_error("runtime rebalance produced conflicting destinations for one particle");
    }
  };
  const auto add_patch_migration = [&](std::uint32_t local_patch_index, int new_owner_rank) {
    if (local_patch_index >= state.patches.size()) {
      throw std::out_of_range("runtime rebalance attempted to migrate an AMR patch index outside SimulationState");
    }
    if (new_owner_rank < 0 || new_owner_rank >= mpi_context.worldSize()) {
      throw std::invalid_argument("runtime rebalance produced AMR patch migration target outside MPI world");
    }
    if (new_owner_rank == world_rank) {
      return;
    }
    if (state.patches.owning_rank[local_patch_index] != static_cast<std::uint32_t>(world_rank)) {
      throw std::runtime_error("runtime rebalance attempted to migrate a non-authoritative local AMR patch");
    }
    const auto [it, inserted] = outbound_patch_target_by_local_index.emplace(local_patch_index, new_owner_rank);
    if (!inserted && it->second != new_owner_rank) {
      throw std::runtime_error("runtime rebalance produced conflicting destinations for one AMR patch");
    }
  };

  for (const auto& intent : rebalance.particle_migrations) {
    if (intent.old_owner_rank != world_rank || intent.new_owner_rank == world_rank) {
      continue;
    }
    const auto found = local_index_by_particle_id.find(intent.particle_id);
    if (found == local_index_by_particle_id.end()) {
      throw std::runtime_error("runtime rebalance migration intent references a particle ID missing on the owning rank");
    }
    add_particle_migration(found->second, intent.new_owner_rank, "particle_sfc_intent");
  }

  for (const auto& update : rebalance.amr_patch_ownership_updates) {
    if (update.new_owner_rank < 0 || update.new_owner_rank >= mpi_context.worldSize()) {
      throw std::runtime_error("runtime rebalance produced AMR patch owner outside MPI world");
    }
    for (std::size_t patch_index = 0; patch_index < state.patches.size(); ++patch_index) {
      if (state.patches.patch_id[patch_index] != update.patch_id) {
        continue;
      }
      if (state.patches.owning_rank[patch_index] == static_cast<std::uint32_t>(world_rank) &&
          update.new_owner_rank != world_rank) {
        add_patch_migration(static_cast<std::uint32_t>(patch_index), update.new_owner_rank);
        const std::uint32_t first_cell = state.patches.first_cell[patch_index];
        const std::uint32_t cell_count = state.patches.cell_count[patch_index];
        if (static_cast<std::uint64_t>(first_cell) + static_cast<std::uint64_t>(cell_count) > state.cells.size()) {
          throw std::runtime_error("runtime rebalance AMR patch cell range is outside CellSoa");
        }
        for (std::uint32_t cell_offset = 0; cell_offset < cell_count; ++cell_offset) {
          const std::uint32_t cell_index = first_cell + cell_offset;
          const auto gas_particle_index = parentParticleRowForGasCellRow(
              state,
              cell_index,
              local_index_by_particle_id,
              "runtime rebalance AMR patch ownership update");
          if (gas_particle_index.has_value() &&
              state.particle_sidecar.owning_rank[*gas_particle_index] == static_cast<std::uint32_t>(world_rank)) {
            add_particle_migration(*gas_particle_index, update.new_owner_rank, "amr_patch_ownership_update");
          }
        }
      }
    }
  }

  std::vector<std::uint32_t> outbound_local_indices;
  outbound_local_indices.reserve(outbound_target_by_local_index.size());
  for (const auto& [local_index, target_rank] : outbound_target_by_local_index) {
    (void)target_rank;
    outbound_local_indices.push_back(local_index);
  }
  std::sort(outbound_local_indices.begin(), outbound_local_indices.end());

  std::vector<std::uint32_t> outbound_patch_indices;
  outbound_patch_indices.reserve(outbound_patch_target_by_local_index.size());
  for (const auto& [patch_index, target_rank] : outbound_patch_target_by_local_index) {
    (void)target_rank;
    outbound_patch_indices.push_back(patch_index);
  }
  std::sort(outbound_patch_indices.begin(), outbound_patch_indices.end());

  std::vector<std::uint32_t> kept_cell_rows;
  kept_cell_rows.reserve(state.cells.size());
  for (std::uint32_t cell_row = 0; cell_row < state.cells.size(); ++cell_row) {
    const std::uint32_t patch_index = state.cells.patch_index[cell_row];
    if (patch_index >= state.patches.size()) {
      throw std::runtime_error("runtime rebalance found a gas cell with invalid patch_index");
    }
    if (!std::binary_search(outbound_patch_indices.begin(), outbound_patch_indices.end(), patch_index)) {
      kept_cell_rows.push_back(cell_row);
    }
  }
  std::vector<core::TimeBinSchedulerIdentityRecord> gas_scheduler_records =
      core::exportGasCellSchedulerIdentityRecords(gas_cell_scheduler, state, kept_cell_rows);

  std::vector<std::uint32_t> preserved_indices;
  preserved_indices.reserve(state.particles.size() - outbound_local_indices.size());
  for (std::uint32_t particle_index = 0; particle_index < state.particles.size(); ++particle_index) {
    if (!std::binary_search(outbound_local_indices.begin(), outbound_local_indices.end(), particle_index)) {
      preserved_indices.push_back(particle_index);
    }
  }
  std::vector<core::ParticleMigrationRecord> scheduler_records =
      state.packParticleMigrationRecords(preserved_indices, scheduler);

  std::vector<std::vector<core::ParticleMigrationRecord>> outbound_records_by_rank(
      static_cast<std::size_t>(mpi_context.worldSize()));
  for (const std::uint32_t local_index : outbound_local_indices) {
    const int target_rank = outbound_target_by_local_index.at(local_index);
    auto records = state.packParticleMigrationRecords(std::span<const std::uint32_t>(&local_index, 1), scheduler);
    if (records.size() != 1U) {
      throw std::runtime_error("runtime rebalance particle migration packing returned an unexpected record count");
    }
    records[0].owning_rank = static_cast<std::uint32_t>(target_rank);
    outbound_records_by_rank[static_cast<std::size_t>(target_rank)].push_back(std::move(records[0]));
  }

  std::vector<std::vector<core::AmrPatchMigrationRecord>> outbound_patch_records_by_rank(
      static_cast<std::size_t>(mpi_context.worldSize()));
  for (const std::uint32_t patch_index : outbound_patch_indices) {
    const int target_rank = outbound_patch_target_by_local_index.at(patch_index);
    auto records = state.packAmrPatchMigrationRecords(std::span<const std::uint32_t>(&patch_index, 1));
    if (records.size() != 1U) {
      throw std::runtime_error("runtime rebalance AMR patch migration packing returned an unexpected record count");
    }
    records[0].patch.owning_rank = static_cast<std::uint32_t>(target_rank);
    const std::uint32_t first_cell = state.patches.first_cell[patch_index];
    const std::uint32_t cell_count = state.patches.cell_count[patch_index];
    std::vector<std::uint32_t> patch_cell_rows;
    patch_cell_rows.reserve(cell_count);
    for (std::uint32_t offset = 0; offset < cell_count; ++offset) {
      patch_cell_rows.push_back(first_cell + offset);
    }
    const std::vector<core::TimeBinSchedulerIdentityRecord> patch_scheduler_records =
        core::exportGasCellSchedulerIdentityRecords(gas_cell_scheduler, state, patch_cell_rows);
    records[0].gas_cell_scheduler_records.reserve(patch_scheduler_records.size());
    for (const auto& scheduler_record : patch_scheduler_records) {
      records[0].gas_cell_scheduler_records.push_back(core::GasCellSchedulerMigrationRecord{
          .gas_cell_id = scheduler_record.element_id,
          .bin_index = scheduler_record.bin_index,
          .next_activation_tick = scheduler_record.next_activation_tick,
          .pending_bin_index = scheduler_record.pending_bin_index,
      });
    }
    for (core::GasCellMigrationRecord& gas_record : records[0].gas_cell_records) {
      gas_record.owning_rank = static_cast<std::uint32_t>(target_rank);
      gas_record.fields.patch_index = 0U;
      gas_record.fields.owning_patch_id = records[0].patch.patch_id;
    }
    outbound_patch_records_by_rank[static_cast<std::size_t>(target_rank)].push_back(std::move(records[0]));
  }

  std::vector<core::ParticleMigrationRecord> inbound_records =
      exchangeRuntimeParticleMigrationRecords(mpi_context, outbound_records_by_rank);
  std::vector<core::AmrPatchMigrationRecord> inbound_patch_records =
      exchangeRuntimeAmrPatchMigrationRecords(mpi_context, outbound_patch_records_by_rank);
  for (const core::ParticleMigrationRecord& record : inbound_records) {
    if (record.owning_rank != static_cast<std::uint32_t>(world_rank)) {
      throw std::runtime_error("runtime particle migration exchange delivered a record to the wrong destination rank");
    }
  }
  for (const core::AmrPatchMigrationRecord& record : inbound_patch_records) {
    if (record.patch.owning_rank != static_cast<std::uint32_t>(world_rank)) {
      throw std::runtime_error("runtime AMR patch migration exchange delivered a patch to the wrong destination rank");
    }
    if (record.gas_cell_scheduler_records.size() != record.gas_cell_records.size()) {
      throw std::runtime_error("runtime AMR patch migration payload is missing gas-cell scheduler identity records");
    }
    for (const core::GasCellSchedulerMigrationRecord& scheduler_record : record.gas_cell_scheduler_records) {
      gas_scheduler_records.push_back(core::TimeBinSchedulerIdentityRecord{
          .element_id = scheduler_record.gas_cell_id,
          .bin_index = scheduler_record.bin_index,
          .next_activation_tick = scheduler_record.next_activation_tick,
          .pending_bin_index = scheduler_record.pending_bin_index,
      });
    }
  }

  if (!outbound_patch_indices.empty() || !inbound_patch_records.empty()) {
    core::AmrPatchMigrationCommit patch_commit;
    patch_commit.world_rank = world_rank;
    patch_commit.outbound_local_patch_indices = outbound_patch_indices;
    patch_commit.inbound_records = inbound_patch_records;
    state.commitAmrPatchMigration(patch_commit);
    core::rebuildSchedulerFromGasCellIdentityRecords(gas_cell_scheduler, gas_scheduler_records, state);
    core::syncGasCellTimeBinMirrorsFromGasCellScheduler(gas_cell_scheduler, state);
  }

  if (!outbound_local_indices.empty() || !inbound_records.empty()) {
    scheduler_records.insert(scheduler_records.end(), inbound_records.begin(), inbound_records.end());
    core::ParticleMigrationCommit commit;
    commit.world_rank = world_rank;
    commit.outbound_local_indices = outbound_local_indices;
    commit.inbound_records = inbound_records;
    commit.preserve_gas_cell_state = true;
    state.commitParticleMigration(commit);

    std::vector<std::uint64_t> destination_particle_ids(state.particle_sidecar.particle_id.begin(),
                                                        state.particle_sidecar.particle_id.end());
    core::rebuildSchedulerFromParticleMigrationRecords(scheduler, scheduler_records, destination_particle_ids);
    syncTimeBinsFromScheduler(scheduler, state);
    if (config.parallel.decomposition_debug_exact_ownership_audit) {
      requireGlobalOwnedParticlePartitionIdentity(
          state, mpi_context, expected_global_particle_ids, "runtime rebalance particle migration commit");
    }
  }
  exchangeAndValidateAmrPatchPayloads(state, mpi_context, world_rank, step_index, profiler);

  recordRuntimeRebalanceDecision(profiler, rebalance, step_index);
  if (profiler != nullptr && (!outbound_local_indices.empty() || !inbound_records.empty())) {
    profiler->recordEvent(core::RuntimeEvent{
        .event_kind = "parallel.decomposition.runtime_migration_commit",
        .severity = core::RuntimeEventSeverity::kWarning,
        .subsystem = "parallel.domain_decomposition",
        .step_index = step_index,
        .message = "runtime rebalance committed authoritative particle migration at a safe scheduler boundary",
        .payload = {
            {"outbound_particle_count", std::to_string(outbound_local_indices.size())},
            {"inbound_particle_count", std::to_string(inbound_records.size())},
            {"current_weighted_imbalance_ratio", std::to_string(rebalance.current_metrics.weighted_imbalance_ratio)},
            {"target_weighted_imbalance_ratio", std::to_string(rebalance.target_decomposition.metrics.weighted_imbalance_ratio)},
            {"current_memory_imbalance_ratio", std::to_string(rebalance.current_metrics.memory_imbalance_ratio)},
            {"target_memory_imbalance_ratio", std::to_string(rebalance.target_decomposition.metrics.memory_imbalance_ratio)},
        },
    });
  }
  const bool local_particle_decomposition_changed =
      state.particleIndexGeneration() != particle_index_generation_before;
  const std::uint64_t changed_rank_count = mpi_context.allreduceSumUint64(
      local_particle_decomposition_changed ? 1ULL : 0ULL);
  return changed_rank_count > 0U;
}

void applyInitialGravityAwareDecomposition(
    core::SimulationState& state,
    const core::SimulationConfig& config,
    int world_size,
    int world_rank,
    core::ProfilerSession* profiler) {
  if (world_size <= 1) {
    return;
  }
  constexpr std::size_t k_density_grid = 16;
  const std::size_t grid_cells = k_density_grid * k_density_grid * k_density_grid;
  std::vector<std::uint32_t> occupancy(grid_cells, 0U);
  std::vector<std::uint32_t> active_occupancy(grid_cells, 0U);
  std::vector<std::uint32_t> gas_occupancy(grid_cells, 0U);
  std::vector<std::uint32_t> pm_x_occupancy(static_cast<std::size_t>(std::max(config.numerics.treepm_pm_grid_nx, 1)), 0U);
  const auto wrap = [](double x, double box) {
    if (box <= 0.0) {
      return x;
    }
    double wrapped = std::fmod(x, box);
    if (wrapped < 0.0) {
      wrapped += box;
    }
    return wrapped;
  };
  const auto density_cell_index = [&](double x, double y, double z) {
    const double box_x = config.cosmology.box_size_x_mpc_comoving;
    const double box_y = config.cosmology.box_size_y_mpc_comoving;
    const double box_z = config.cosmology.box_size_z_mpc_comoving;
    const std::size_t ix = (box_x > 0.0)
        ? std::min<std::size_t>(
              k_density_grid - 1U,
              static_cast<std::size_t>((wrap(x, box_x) / box_x) * static_cast<double>(k_density_grid)))
        : 0U;
    const std::size_t iy = (box_y > 0.0)
        ? std::min<std::size_t>(
              k_density_grid - 1U,
              static_cast<std::size_t>((wrap(y, box_y) / box_y) * static_cast<double>(k_density_grid)))
        : 0U;
    const std::size_t iz = (box_z > 0.0)
        ? std::min<std::size_t>(
              k_density_grid - 1U,
              static_cast<std::size_t>((wrap(z, box_z) / box_z) * static_cast<double>(k_density_grid)))
        : 0U;
    return (ix * k_density_grid + iy) * k_density_grid + iz;
  };
  const auto pm_x_index = [&](double x) {
    const double box_x = config.cosmology.box_size_x_mpc_comoving;
    const std::size_t nx = pm_x_occupancy.size();
    if (box_x <= 0.0 || nx == 0) {
      return std::size_t{0};
    }
    const double scaled = (wrap(x, box_x) / box_x) * static_cast<double>(nx);
    return std::min<std::size_t>(nx - 1U, static_cast<std::size_t>(scaled));
  };

  for (std::size_t particle_index = 0; particle_index < state.particles.size(); ++particle_index) {
    const std::size_t cell = density_cell_index(
        state.particles.position_x_comoving[particle_index],
        state.particles.position_y_comoving[particle_index],
        state.particles.position_z_comoving[particle_index]);
    ++occupancy[cell];
    if (!state.particles.time_bin.empty() && state.particles.time_bin[particle_index] == 0U) {
      ++active_occupancy[cell];
    }
    if (state.particle_sidecar.species_tag[particle_index] == static_cast<std::uint32_t>(core::ParticleSpecies::kGas)) {
      ++gas_occupancy[cell];
    }
    ++pm_x_occupancy[pm_x_index(state.particles.position_x_comoving[particle_index])];
  }

  std::vector<std::uint32_t> patch_cell_count(state.patches.size(), 0U);
  for (std::size_t cell_index = 0; cell_index < state.cells.size(); ++cell_index) {
    const std::uint32_t patch_index = state.cells.patch_index[cell_index];
    if (patch_index < patch_cell_count.size()) {
      ++patch_cell_count[patch_index];
    }
  }

  const auto particle_memory_bytes = [&](std::uint32_t species_tag) {
    std::uint64_t bytes = sizeof(double) * 7U + sizeof(std::uint64_t) * 2U + sizeof(std::uint32_t) * 3U;
    if (!state.particle_sidecar.gravity_softening_comoving.empty()) {
      bytes += sizeof(double);
    }
    if (!state.particle_sidecar.has_gravity_softening_override.empty()) {
      bytes += sizeof(std::uint8_t);
    }
    if (species_tag == static_cast<std::uint32_t>(core::ParticleSpecies::kGas)) {
      bytes += sizeof(double) * 8U + sizeof(std::uint64_t) * 2U + sizeof(std::uint32_t);
    } else if (species_tag == static_cast<std::uint32_t>(core::ParticleSpecies::kStar)) {
      bytes += sizeof(std::uint32_t) + sizeof(double) * 13U;
    } else if (species_tag == static_cast<std::uint32_t>(core::ParticleSpecies::kBlackHole)) {
      bytes += sizeof(std::uint32_t) * 2U + sizeof(double) * 8U;
    } else if (species_tag == static_cast<std::uint32_t>(core::ParticleSpecies::kTracer)) {
      bytes += sizeof(std::uint64_t) * 2U + sizeof(std::uint32_t) * 2U + sizeof(double) * 3U;
    }
    return bytes;
  };

  std::vector<parallel::DecompositionItem> items;
  items.reserve(state.particles.size() + state.patches.size());
  constexpr std::uint32_t k_invalid_patch_index = std::numeric_limits<std::uint32_t>::max();
  std::vector<std::uint32_t> patch_index_by_item;
  patch_index_by_item.reserve(state.particles.size() + state.patches.size());
  for (std::size_t particle_index = 0; particle_index < state.particles.size(); ++particle_index) {
    parallel::DecompositionItem item;
    item.entity_id = state.particle_sidecar.particle_id[particle_index];
    item.kind = parallel::DecompositionEntityKind::kParticle;
    item.current_owner_rank = static_cast<int>(state.particle_sidecar.owning_rank[particle_index]);
    item.x_comov = state.particles.position_x_comoving[particle_index];
    item.y_comov = state.particles.position_y_comoving[particle_index];
    item.z_comov = state.particles.position_z_comoving[particle_index];
    const std::size_t density_cell = density_cell_index(item.x_comov, item.y_comov, item.z_comov);
    const std::uint32_t local_density = occupancy[density_cell];
    const std::uint32_t local_active = active_occupancy[density_cell];
    const std::uint32_t local_gas = gas_occupancy[density_cell];
    const std::uint32_t pm_load = pm_x_occupancy[pm_x_index(item.x_comov)];
    const std::uint32_t species_tag = state.particle_sidecar.species_tag[particle_index];
    double amr_patch_cost = 0.0;
    if (species_tag == static_cast<std::uint32_t>(core::ParticleSpecies::kGas) && !state.cells.patch_index.empty()) {
      std::vector<std::uint32_t> seen_patch_indices;
      for (const std::uint32_t cell_row : state.gas_cell_identity.rowsForParentParticleId(item.entity_id)) {
        if (cell_row >= state.cells.patch_index.size()) {
          continue;
        }
        const std::uint32_t patch_index = state.cells.patch_index[cell_row];
        if (patch_index >= patch_cell_count.size() ||
            std::find(seen_patch_indices.begin(), seen_patch_indices.end(), patch_index) != seen_patch_indices.end()) {
          continue;
        }
        seen_patch_indices.push_back(patch_index);
        amr_patch_cost += static_cast<double>(patch_cell_count[patch_index]);
      }
    }
    const double local_density_d = static_cast<double>(std::max<std::uint32_t>(local_density, 1U));
    item.active_target_count_recent = local_active;
    item.remote_tree_interactions_recent = static_cast<std::uint64_t>(
        std::llround(local_density_d * std::log2(local_density_d + 1.0)));
    item.work_units = 1.0;
    item.memory_bytes = particle_memory_bytes(species_tag);
    item.work_components = parallel::DecompositionWorkComponents{
        .particle_count_cost = 1.0,
        .gas_cell_cost = (species_tag == static_cast<std::uint32_t>(core::ParticleSpecies::kGas))
            ? (1.0 + static_cast<double>(local_gas))
            : 0.0,
        .tree_interaction_cost = static_cast<double>(item.remote_tree_interactions_recent),
        .pm_mesh_cost = static_cast<double>(pm_load),
        .amr_patch_cost = amr_patch_cost,
        .active_fraction_cost = static_cast<double>(local_active),
        .memory_pressure_cost = static_cast<double>(item.memory_bytes),
        .gpu_occupancy_cost = 0.0,
        .generic_work_cost = 1.0 + std::sqrt(local_density_d),
        .has_explicit_components = true,
    };
    items.push_back(item);
    patch_index_by_item.push_back(k_invalid_patch_index);
  }

  if (state.patches.owning_rank.size() != state.patches.size()) {
    state.patches.owning_rank.assign(state.patches.size(), static_cast<std::uint32_t>(std::max(world_rank, 0)));
  }
  for (std::size_t patch_index = 0; patch_index < state.patches.size(); ++patch_index) {
    if (state.patches.cell_count[patch_index] == 0U) {
      continue;
    }
    const std::uint32_t first_cell = state.patches.first_cell[patch_index];
    if (first_cell >= state.cells.size()) {
      continue;
    }
    parallel::DecompositionItem patch_item;
    patch_item.entity_id = state.patches.patch_id[patch_index];
    patch_item.kind = parallel::DecompositionEntityKind::kAmrPatch;
    patch_item.current_owner_rank = static_cast<int>(state.patches.owning_rank[patch_index]);
    patch_item.x_comov = state.cells.center_x_comoving[first_cell];
    patch_item.y_comov = state.cells.center_y_comoving[first_cell];
    patch_item.z_comov = state.cells.center_z_comoving[first_cell];
    patch_item.memory_bytes = static_cast<std::uint64_t>(state.patches.cell_count[patch_index]) *
        static_cast<std::uint64_t>(sizeof(double) * 8U + sizeof(std::uint32_t));
    patch_item.work_components = parallel::DecompositionWorkComponents{
        .amr_patch_cost = static_cast<double>(state.patches.cell_count[patch_index]) *
            (1.0 + static_cast<double>(std::max(state.patches.level[patch_index], 0))),
        .memory_pressure_cost = static_cast<double>(patch_item.memory_bytes),
        .generic_work_cost = static_cast<double>(state.patches.cell_count[patch_index]),
        .has_explicit_components = true,
    };
    items.push_back(patch_item);
    patch_index_by_item.push_back(static_cast<std::uint32_t>(patch_index));
  }

  parallel::DecompositionConfig decomposition_config;
  decomposition_config.world_size = world_size;
  decomposition_config.domain_x_min_comov = 0.0;
  decomposition_config.domain_x_max_comov = config.cosmology.box_size_x_mpc_comoving;
  decomposition_config.domain_y_min_comov = 0.0;
  decomposition_config.domain_y_max_comov = config.cosmology.box_size_y_mpc_comoving;
  decomposition_config.domain_z_min_comov = 0.0;
  decomposition_config.domain_z_max_comov = config.cosmology.box_size_z_mpc_comoving;
  decomposition_config.owned_particle_weight = 0.0;
  decomposition_config.active_target_weight = 0.0;
  decomposition_config.remote_tree_interaction_weight = 0.0;
  decomposition_config.work_weight = 0.0;
  decomposition_config.memory_weight = 0.0;
  decomposition_config.component_weights = parallel::DecompositionWeightCoefficients{
      .particle_count = config.parallel.decomposition_particle_count_weight,
      .gas_cell = config.parallel.decomposition_gas_cell_weight,
      .tree_interaction = config.parallel.decomposition_tree_interaction_weight,
      .pm_mesh = config.parallel.decomposition_pm_mesh_weight,
      .amr_patch = config.parallel.decomposition_amr_patch_weight,
      .active_fraction = config.parallel.decomposition_active_fraction_weight,
      .memory_pressure = config.parallel.decomposition_memory_pressure_weight,
      .gpu_occupancy = config.parallel.decomposition_gpu_occupancy_weight,
      .generic_work = config.parallel.decomposition_generic_work_weight,
  };
  const auto plan = parallel::buildMortonSfcDecomposition(items, decomposition_config);
  for (std::size_t item_index = 0; item_index < state.particles.size(); ++item_index) {
    state.particle_sidecar.owning_rank[item_index] = static_cast<std::uint32_t>(plan.owning_rank_by_item[item_index]);
  }
  for (std::size_t item_index = 0; item_index < patch_index_by_item.size(); ++item_index) {
    const std::uint32_t patch_index = patch_index_by_item[item_index];
    if (patch_index != k_invalid_patch_index && patch_index < state.patches.owning_rank.size()) {
      state.patches.owning_rank[patch_index] = static_cast<std::uint32_t>(plan.owning_rank_by_item[item_index]);
    }
  }
  parallel::recordDistributedProfiling(profiler, plan.metrics, 0, 0);
  if (profiler != nullptr) {
    profiler->recordEvent(core::RuntimeEvent{
        .event_kind = "parallel.decomposition.work_weighted",
        .severity = core::RuntimeEventSeverity::kInfo,
        .subsystem = "parallel.domain_decomposition",
        .message = "initial domain decomposition used explicit work-weight components",
        .payload = {{"world_size", std::to_string(world_size)},
                    {"item_count", std::to_string(items.size())},
                    {"weighted_imbalance_ratio", std::to_string(plan.metrics.weighted_imbalance_ratio)},
                    {"memory_imbalance_ratio", std::to_string(plan.metrics.memory_imbalance_ratio)}}});
  }
  compactStateToCurrentOwner(state, world_rank);
}




}  // namespace

MigrationBalanceRuntime::MigrationBalanceRuntime(
    const core::SimulationConfig& config,
    const RuntimeServices& services) noexcept
    : m_config(config), m_services(services) {}

void MigrationBalanceRuntime::initializeOwnership(
    core::SimulationState& state) const {
  const std::size_t pm_grid_nx =
      static_cast<std::size_t>(m_config.numerics.treepm_pm_grid_nx);
  const int world_size = m_services.mpi_context.worldSize();
  const int world_rank = m_services.mpi_context.worldRank();
  if (pm_grid_nx == 0U || world_size <= 0) {
    throw std::invalid_argument(
        "initial ownership requires a nonzero PM x extent and world size");
  }
  const double box_size_x = m_config.cosmology.box_size_x_mpc_comoving;
  for (std::size_t index = 0; index < state.particles.size(); ++index) {
    double wrapped_x = state.particles.position_x_comoving[index];
    if (box_size_x > 0.0) {
      wrapped_x = std::fmod(wrapped_x, box_size_x);
      if (wrapped_x < 0.0) {
        wrapped_x += box_size_x;
      }
      if (wrapped_x >= box_size_x) {
        wrapped_x -= box_size_x;
      }
    }
    std::size_t global_x = 0U;
    if (box_size_x > 0.0) {
      const double scaled =
          (wrapped_x / box_size_x) * static_cast<double>(pm_grid_nx);
      global_x = std::min(
          pm_grid_nx - 1U, static_cast<std::size_t>(scaled));
    }
    state.particle_sidecar.owning_rank[index] =
        static_cast<std::uint32_t>(parallel::pmOwnerRankForGlobalX(
            pm_grid_nx, world_size, global_x));
  }
  applyInitialGravityAwareDecomposition(
      state,
      m_config,
      world_size,
      world_rank,
      &m_services.profiler);
}

parallel::LocalOwnershipIdentitySummary
MigrationBalanceRuntime::reduceIdentity(
    const core::SimulationState& state) const {
  return reduceLocalParticleIdentitySummary(state, m_services.mpi_context);
}

bool MigrationBalanceRuntime::rebalance(
    core::SimulationState& state,
    core::HierarchicalTimeBinScheduler& particle_scheduler,
    core::HierarchicalTimeBinScheduler& gas_cell_scheduler,
    const parallel::DecompositionRuntimeMeasurements& measurements,
    std::span<const std::uint32_t> active_particle_indices,
    std::span<const std::uint64_t> expected_global_particle_ids,
    std::uint64_t step_index) const {
  return applyMeasuredRuntimeRebalancePlan(
      state,
      particle_scheduler,
      gas_cell_scheduler,
      m_config,
      m_services.mpi_context,
      m_services.mpi_context.worldRank(),
      measurements,
      active_particle_indices,
      expected_global_particle_ids,
      &m_services.profiler,
      step_index);
}

}  // namespace cosmosim::workflows::internal
