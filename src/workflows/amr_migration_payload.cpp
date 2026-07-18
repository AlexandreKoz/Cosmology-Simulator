#include "workflows/internal/amr_migration_payload.hpp"
#include "workflows/internal/gas_cell_ownership.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace cosmosim::workflows::internal {

std::vector<parallel::AmrPatchPayloadRecord>
buildMigrationAmrPatchPayloadRecords(
    const core::SimulationState& state,
    int world_rank) {
  std::vector<parallel::AmrPatchPayloadRecord> records;
  if (world_rank < 0 || state.patches.size() == 0) {
    return records;
  }
  records.reserve(state.patches.size());
  for (std::size_t patch_index = 0; patch_index < state.patches.size();
       ++patch_index) {
    if (state.patches.owning_rank[patch_index] !=
        static_cast<std::uint32_t>(world_rank)) {
      continue;
    }
    const std::uint32_t first_cell = state.patches.first_cell[patch_index];
    const std::uint32_t cell_count = state.patches.cell_count[patch_index];
    if (cell_count == 0U) {
      continue;
    }
    if (static_cast<std::uint64_t>(first_cell) +
            static_cast<std::uint64_t>(cell_count) >
        state.cells.size()) {
      throw std::runtime_error(
          "AMR patch payload build found a patch range outside CellSoa");
    }
    parallel::AmrPatchPayloadRecord record;
    record.patch_id = state.patches.patch_id[patch_index];
    record.parent_patch_id = state.patches.parent_patch_id[patch_index];
    record.morton_key = state.patches.morton_key[patch_index];
    record.owner_rank = world_rank;
    record.level = static_cast<std::uint32_t>(
        std::max<std::int32_t>(state.patches.level[patch_index], 0));
    record.first_cell = first_cell;
    record.cell_count = cell_count;
    record.origin_x_comoving = state.patches.origin_x_comoving[patch_index];
    record.origin_y_comoving = state.patches.origin_y_comoving[patch_index];
    record.origin_z_comoving = state.patches.origin_z_comoving[patch_index];
    record.extent_x_comoving = state.patches.extent_x_comoving[patch_index];
    record.extent_y_comoving = state.patches.extent_y_comoving[patch_index];
    record.extent_z_comoving = state.patches.extent_z_comoving[patch_index];
    record.cell_dim_x = state.patches.cell_dim_x[patch_index];
    record.cell_dim_y = state.patches.cell_dim_y[patch_index];
    record.cell_dim_z = state.patches.cell_dim_z[patch_index];
    record.decomposition_epoch = state.gasCellIdentityGeneration();
    for (std::uint32_t offset = 0; offset < cell_count; ++offset) {
      const std::uint32_t cell_index = first_cell + offset;
      record.cell_mass_sum_code += state.cells.mass_code[cell_index];
      if (cell_index < state.gas_cells.internal_energy_code.size()) {
        record.gas_internal_energy_sum_code +=
            state.gas_cells.internal_energy_code[cell_index];
      }
    }
    parallel::validateAmrPatchPayloadRecord(record);
    records.push_back(record);
  }
  return records;
}

std::vector<parallel::AmrPatchCellPayloadRecord>
buildMigrationAmrPatchCellPayloadRecords(
    const core::SimulationState& state,
    int world_rank) {
  std::vector<parallel::AmrPatchCellPayloadRecord> records;
  if (world_rank < 0 || state.patches.size() == 0) {
    return records;
  }
  std::size_t total_owned_cells = 0;
  for (std::size_t patch_index = 0; patch_index < state.patches.size();
       ++patch_index) {
    if (state.patches.owning_rank[patch_index] ==
        static_cast<std::uint32_t>(world_rank)) {
      total_owned_cells += state.patches.cell_count[patch_index];
    }
  }
  records.reserve(total_owned_cells);
  for (std::size_t patch_index = 0; patch_index < state.patches.size();
       ++patch_index) {
    if (state.patches.owning_rank[patch_index] !=
        static_cast<std::uint32_t>(world_rank)) {
      continue;
    }
    const std::uint32_t first_cell = state.patches.first_cell[patch_index];
    const std::uint32_t cell_count = state.patches.cell_count[patch_index];
    if (static_cast<std::uint64_t>(first_cell) +
            static_cast<std::uint64_t>(cell_count) >
        state.cells.size()) {
      throw std::runtime_error(
          "AMR patch cell payload build found a patch range outside CellSoa");
    }
    for (std::uint32_t offset = 0; offset < cell_count; ++offset) {
      const std::uint32_t cell_index = first_cell + offset;
      parallel::AmrPatchCellPayloadRecord record;
      record.patch_id = state.patches.patch_id[patch_index];
      record.owner_rank = world_rank;
      record.local_cell_offset = offset;
      record.patch_index = static_cast<std::uint32_t>(patch_index);
      record.center_x_comoving = state.cells.center_x_comoving[cell_index];
      record.center_y_comoving = state.cells.center_y_comoving[cell_index];
      record.center_z_comoving = state.cells.center_z_comoving[cell_index];
      record.mass_code = state.cells.mass_code[cell_index];
      record.time_bin = state.cells.time_bin[cell_index];
      const core::GasCellIdentityRecord& identity =
          gasCellIdentityRecordForLocalRow(
              state, cell_index, "buildMigrationAmrPatchCellPayloadRecords");
      record.gas_cell_id = identity.gas_cell_id;
      record.parent_particle_id = identity.parent_particle_id.value_or(0U);
      record.velocity_x_peculiar =
          cell_index < state.gas_cells.velocity_x_peculiar.size()
              ? state.gas_cells.velocity_x_peculiar[cell_index]
              : 0.0;
      record.velocity_y_peculiar =
          cell_index < state.gas_cells.velocity_y_peculiar.size()
              ? state.gas_cells.velocity_y_peculiar[cell_index]
              : 0.0;
      record.velocity_z_peculiar =
          cell_index < state.gas_cells.velocity_z_peculiar.size()
              ? state.gas_cells.velocity_z_peculiar[cell_index]
              : 0.0;
      record.density_code = cell_index < state.gas_cells.density_code.size()
          ? state.gas_cells.density_code[cell_index]
          : 0.0;
      record.pressure_code = cell_index < state.gas_cells.pressure_code.size()
          ? state.gas_cells.pressure_code[cell_index]
          : 0.0;
      record.internal_energy_code =
          cell_index < state.gas_cells.internal_energy_code.size()
              ? state.gas_cells.internal_energy_code[cell_index]
              : 0.0;
      record.temperature_code =
          cell_index < state.gas_cells.temperature_code.size()
              ? state.gas_cells.temperature_code[cell_index]
              : 0.0;
      record.sound_speed_code =
          cell_index < state.gas_cells.sound_speed_code.size()
              ? state.gas_cells.sound_speed_code[cell_index]
              : 0.0;
      parallel::validateAmrPatchCellPayloadRecord(record);
      records.push_back(record);
    }
  }
  return records;
}

}  // namespace cosmosim::workflows::internal
