#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "cosmosim/amr/amr_framework.hpp"
#include "cosmosim/core/simulation_state.hpp"

namespace cosmosim::amr {

struct PatchLocalCellIndex {
  std::size_t i = 0;
  std::size_t j = 0;
  std::size_t k = 0;
  std::size_t linear_index = 0;
};

[[nodiscard]] inline std::size_t checkedPatchCellCount(
    const PatchDescriptor& patch,
    std::string_view context) {
  const std::size_t nx = patch.cell_dims[0];
  const std::size_t ny = patch.cell_dims[1];
  const std::size_t nz = patch.cell_dims[2];
  if (nx == 0U || ny == 0U || nz == 0U) {
    throw std::invalid_argument(std::string(context) + ": patch cell dimensions must be positive");
  }
  if (nx > std::numeric_limits<std::size_t>::max() / ny) {
    throw std::overflow_error(std::string(context) + ": patch cell dimensions overflow");
  }
  const std::size_t xy = nx * ny;
  if (xy > std::numeric_limits<std::size_t>::max() / nz) {
    throw std::overflow_error(std::string(context) + ": patch cell dimensions overflow");
  }
  return xy * nz;
}

[[nodiscard]] inline std::size_t linearPatchCellIndex(
    std::array<std::uint16_t, 3> dims,
    std::size_t i,
    std::size_t j,
    std::size_t k) {
  return i + static_cast<std::size_t>(dims[0]) * (j + static_cast<std::size_t>(dims[1]) * k);
}

[[nodiscard]] inline PatchLocalCellIndex patchLocalCellForPoint(
    const PatchDescriptor& patch,
    const std::array<double, 3>& point_comov,
    std::string_view context) {
  (void)checkedPatchCellCount(patch, context);
  std::array<std::size_t, 3> ijk{};
  for (std::size_t axis = 0; axis < 3; ++axis) {
    const double extent = patch.extent_comov[axis];
    const std::uint16_t dim = patch.cell_dims[axis];
    if (extent <= 0.0) {
      throw std::invalid_argument(std::string(context) + ": patch extents must be positive");
    }
    const double lower = patch.origin_comov[axis];
    const double upper = patch.origin_comov[axis] + extent;
    const double scale = std::max({1.0, std::abs(lower), std::abs(upper), std::abs(point_comov[axis])});
    if (point_comov[axis] < lower - 1.0e-10 * scale ||
        point_comov[axis] >= upper + 1.0e-10 * scale) {
      throw std::runtime_error(std::string(context) + ": gas-cell center lies outside explicit patch geometry");
    }
    const double cell_width = extent / static_cast<double>(dim);
    double coordinate = (point_comov[axis] - lower) / cell_width;
    coordinate = std::clamp(coordinate, 0.0, static_cast<double>(dim) - 1.0e-12);
    ijk[axis] = static_cast<std::size_t>(coordinate);
    if (ijk[axis] >= dim) {
      ijk[axis] = static_cast<std::size_t>(dim) - 1U;
    }
  }
  return PatchLocalCellIndex{
      .i = ijk[0],
      .j = ijk[1],
      .k = ijk[2],
      .linear_index = linearPatchCellIndex(patch.cell_dims, ijk[0], ijk[1], ijk[2])};
}

[[nodiscard]] inline PatchLocalCellIndex patchLocalCellForRow(
    const core::SimulationState& state,
    const PatchDescriptor& patch,
    std::uint32_t row,
    std::string_view context) {
  if (row >= state.cells.size()) {
    throw std::out_of_range(std::string(context) + ": gas-cell row is outside SimulationState cell storage");
  }
  const std::array<double, 3> point{
      state.cells.center_x_comoving.at(row),
      state.cells.center_y_comoving.at(row),
      state.cells.center_z_comoving.at(row)};
  return patchLocalCellForPoint(patch, point, context);
}

[[nodiscard]] inline std::vector<std::uint32_t> buildPatchLocalRowMap(
    const core::SimulationState& state,
    const PatchDescriptor& patch,
    std::span<const std::uint32_t> rows,
    std::string_view context) {
  const std::size_t expected_cells = checkedPatchCellCount(patch, context);
  if (rows.size() != expected_cells) {
    throw std::runtime_error(std::string(context) + ": patch gas-cell coverage does not match descriptor cell_dims");
  }
  std::vector<std::uint32_t> row_by_patch_cell(expected_cells, core::kInvalidGasCellRow);
  for (const std::uint32_t row : rows) {
    if (row >= state.cells.size() || row >= state.gas_cells.size()) {
      throw std::runtime_error(std::string(context) + ": patch row is outside SimulationState storage");
    }
    const auto* record = state.gas_cell_identity.findByLocalRow(row);
    if (record == nullptr || record->owning_patch_id != patch.patch_id ||
        record->gas_cell_id != state.gas_cells.gas_cell_id[row]) {
      throw std::runtime_error(std::string(context) + ": identity row does not belong to requested patch");
    }
    const std::size_t patch_cell = patchLocalCellForRow(state, patch, row, context).linear_index;
    if (patch_cell >= expected_cells || row_by_patch_cell[patch_cell] != core::kInvalidGasCellRow) {
      throw std::runtime_error(std::string(context) + ": duplicate or out-of-range patch-local gas-cell mapping");
    }
    row_by_patch_cell[patch_cell] = row;
  }
  for (const std::uint32_t row : row_by_patch_cell) {
    if (row == core::kInvalidGasCellRow) {
      throw std::runtime_error(std::string(context) + ": explicit patch geometry has a missing patch-local gas cell");
    }
  }
  return row_by_patch_cell;
}

[[nodiscard]] inline std::vector<std::uint32_t> rowsForPatchLocalCellOrder(
    const core::SimulationState& state,
    const PatchDescriptor& patch,
    std::string_view context) {
  const std::vector<std::uint32_t> rows = state.gas_cell_identity.rowsForPatch(patch.patch_id);
  return buildPatchLocalRowMap(state, patch, std::span<const std::uint32_t>(rows.data(), rows.size()), context);
}

}  // namespace cosmosim::amr
