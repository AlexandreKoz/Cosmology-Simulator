#include "workflows/internal/gas_cell_ownership.hpp"

#include <algorithm>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace cosmosim::workflows::internal {

ParticleRowById buildParticleRowById(
    const core::SimulationState& state) {
  ParticleRowById row_by_id;
  row_by_id.reserve(state.particles.size());
  for (std::uint32_t row = 0; row < state.particles.size(); ++row) {
    const auto [_, inserted] = row_by_id.emplace(
        state.particle_sidecar.particle_id[row], row);
    if (!inserted) {
      throw std::runtime_error(
          "gas-cell ownership lookup encountered a duplicate local particle ID");
    }
  }
  return row_by_id;
}

std::optional<std::uint32_t> parentParticleRowForGasCellRow(
    const core::SimulationState& state,
    std::uint32_t cell_row,
    const ParticleRowById& particle_row_by_id,
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
  if (parent_it == particle_row_by_id.end()) {
    return std::nullopt;
  }
  return parent_it->second;
}

const core::GasCellIdentityRecord& gasCellIdentityRecordForLocalRow(
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

std::uint32_t gasCellOwnerRankForLocalRow(
    const core::SimulationState& state,
    std::uint32_t cell_row,
    const ParticleRowById& particle_row_by_id,
    std::string_view caller) {
  if (cell_row >= state.cells.size()) {
    throw std::out_of_range(
        std::string(caller) + ": gas-cell row is outside CellSoa");
  }
  const std::uint32_t patch_index = state.cells.patch_index[cell_row];
  if (patch_index < state.patches.size()) {
    if (patch_index >= state.patches.owning_rank.size()) {
      throw std::runtime_error(
          std::string(caller) +
          ": patch owning-rank lane is shorter than patch table");
    }
    return state.patches.owning_rank[patch_index];
  }

  const auto& identity =
      gasCellIdentityRecordForLocalRow(state, cell_row, caller);
  if (!identity.parent_particle_id.has_value()) {
    throw std::runtime_error(
        std::string(caller) +
        ": parentless gas cell has no valid patch owner");
  }
  const auto parent_it = particle_row_by_id.find(*identity.parent_particle_id);
  if (parent_it == particle_row_by_id.end()) {
    throw std::runtime_error(
        std::string(caller) +
        ": gas-cell parent_particle_id is not a local particle");
  }
  if (parent_it->second >= state.particle_sidecar.owning_rank.size()) {
    throw std::runtime_error(
        std::string(caller) +
        ": gas-cell parent row is outside owning-rank metadata");
  }
  return state.particle_sidecar.owning_rank[parent_it->second];
}

void synchronizeParentParticleCompatibilityMirrors(
    core::SimulationState& state,
    std::uint32_t world_rank,
    std::string_view caller) {
  state.requireGasCellIdentityMapCoversDenseRows(caller);
  const ParticleRowById particle_row_by_id = buildParticleRowById(state);

  struct ParentMirrorAccumulator {
    double mass_code = 0.0;
    double momentum_x_code = 0.0;
    double momentum_y_code = 0.0;
    double momentum_z_code = 0.0;
  };

  std::unordered_map<std::uint32_t, std::vector<std::uint32_t>> rows_by_parent;
  for (std::uint32_t cell_row = 0; cell_row < state.cells.size(); ++cell_row) {
    if (gasCellOwnerRankForLocalRow(
            state, cell_row, particle_row_by_id, caller) != world_rank) {
      continue;
    }
    const auto parent_row = parentParticleRowForGasCellRow(
        state, cell_row, particle_row_by_id, caller);
    if (parent_row.has_value()) {
      rows_by_parent[*parent_row].push_back(cell_row);
    }
  }

  for (auto& [parent_row, rows] : rows_by_parent) {
    std::sort(
        rows.begin(), rows.end(),
        [&](std::uint32_t lhs, std::uint32_t rhs) {
          return gasCellIdentityRecordForLocalRow(state, lhs, caller).gas_cell_id <
              gasCellIdentityRecordForLocalRow(state, rhs, caller).gas_cell_id;
        });
    ParentMirrorAccumulator aggregate;
    for (const std::uint32_t cell_row : rows) {
      const double mass_code = state.cells.mass_code[cell_row];
      aggregate.mass_code += mass_code;
      aggregate.momentum_x_code +=
          mass_code * state.gas_cells.velocity_x_peculiar[cell_row];
      aggregate.momentum_y_code +=
          mass_code * state.gas_cells.velocity_y_peculiar[cell_row];
      aggregate.momentum_z_code +=
          mass_code * state.gas_cells.velocity_z_peculiar[cell_row];
    }
    if (aggregate.mass_code > 0.0) {
      state.particles.mass_code[parent_row] = aggregate.mass_code;
      state.particles.velocity_x_peculiar[parent_row] =
          aggregate.momentum_x_code / aggregate.mass_code;
      state.particles.velocity_y_peculiar[parent_row] =
          aggregate.momentum_y_code / aggregate.mass_code;
      state.particles.velocity_z_peculiar[parent_row] =
          aggregate.momentum_z_code / aggregate.mass_code;
    }
  }
}

}  // namespace cosmosim::workflows::internal
