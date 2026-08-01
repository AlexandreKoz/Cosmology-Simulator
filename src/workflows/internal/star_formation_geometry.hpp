#pragma once

#include <cstdint>
#include <limits>
#include <span>
#include <unordered_set>

#include "cosmosim/core/simulation_state.hpp"

namespace cosmosim::workflows::internal {

struct StarFormationPatchCellGeometry {
  bool valid = false;
  std::uint32_t patch_index = 0;
  std::uint32_t local_index = 0;
  std::uint32_t ix = 0;
  std::uint32_t iy = 0;
  std::uint32_t iz = 0;
  std::uint32_t nx = 0;
  std::uint32_t ny = 0;
  std::uint32_t nz = 0;
  std::uint32_t first_cell = 0;
  double dx_stored = 0.0;
  double dy_stored = 0.0;
  double dz_stored = 0.0;
};

[[nodiscard]] inline StarFormationPatchCellGeometry starFormationPatchCellGeometry(
    const core::SimulationState& state,
    std::uint32_t cell_index) {
  StarFormationPatchCellGeometry geometry;
  if (cell_index >= state.cells.size() || cell_index >= state.cells.patch_index.size()) {
    return geometry;
  }
  const std::uint32_t patch_index = state.cells.patch_index[cell_index];
  if (patch_index >= state.patches.size()) {
    return geometry;
  }
  const std::uint32_t first_cell = state.patches.first_cell[patch_index];
  const std::uint32_t cell_count = state.patches.cell_count[patch_index];
  if (cell_index < first_cell || cell_index >= first_cell + cell_count) {
    return geometry;
  }
  const std::uint32_t nx = state.patches.cell_dim_x[patch_index];
  const std::uint32_t ny = state.patches.cell_dim_y[patch_index];
  const std::uint32_t nz = state.patches.cell_dim_z[patch_index];
  const double extent_x = state.patches.extent_x_comoving[patch_index];
  const double extent_y = state.patches.extent_y_comoving[patch_index];
  const double extent_z = state.patches.extent_z_comoving[patch_index];
  if (nx == 0U || ny == 0U || nz == 0U || !(extent_x > 0.0) || !(extent_y > 0.0) ||
      !(extent_z > 0.0) || static_cast<std::uint64_t>(nx) * ny * nz != cell_count) {
    return geometry;
  }
  const std::uint32_t local_index = cell_index - first_cell;
  geometry.valid = true;
  geometry.patch_index = patch_index;
  geometry.local_index = local_index;
  geometry.ix = local_index % nx;
  geometry.iy = (local_index / nx) % ny;
  geometry.iz = local_index / (nx * ny);
  geometry.nx = nx;
  geometry.ny = ny;
  geometry.nz = nz;
  geometry.first_cell = first_cell;
  geometry.dx_stored = extent_x / static_cast<double>(nx);
  geometry.dy_stored = extent_y / static_cast<double>(ny);
  geometry.dz_stored = extent_z / static_cast<double>(nz);
  return geometry;
}

[[nodiscard]] inline std::uint32_t starFormationPatchRow(
    const StarFormationPatchCellGeometry& geometry,
    std::uint32_t ix,
    std::uint32_t iy,
    std::uint32_t iz) {
  return geometry.first_cell + ix + geometry.nx * (iy + geometry.ny * iz);
}

[[nodiscard]] inline double starFormationDerivativeAtCell(
    std::span<const double> field,
    const StarFormationPatchCellGeometry& geometry,
    int axis,
    double physical_spacing,
    std::uint32_t center_row) {
  if (!(physical_spacing > 0.0) || center_row >= field.size()) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  std::uint32_t coordinate = 0U;
  std::uint32_t dimension = 0U;
  if (axis == 0) {
    coordinate = geometry.ix;
    dimension = geometry.nx;
  } else if (axis == 1) {
    coordinate = geometry.iy;
    dimension = geometry.ny;
  } else {
    coordinate = geometry.iz;
    dimension = geometry.nz;
  }
  if (dimension <= 1U) {
    return 0.0;
  }

  const auto row_with_axis = [&](std::uint32_t value) {
    if (axis == 0) {
      return starFormationPatchRow(geometry, value, geometry.iy, geometry.iz);
    }
    if (axis == 1) {
      return starFormationPatchRow(geometry, geometry.ix, value, geometry.iz);
    }
    return starFormationPatchRow(geometry, geometry.ix, geometry.iy, value);
  };

  if (coordinate > 0U && coordinate + 1U < dimension) {
    const std::uint32_t left = row_with_axis(coordinate - 1U);
    const std::uint32_t right = row_with_axis(coordinate + 1U);
    if (left >= field.size() || right >= field.size()) {
      return std::numeric_limits<double>::quiet_NaN();
    }
    return (field[right] - field[left]) / (2.0 * physical_spacing);
  }
  if (coordinate == 0U) {
    const std::uint32_t right = row_with_axis(1U);
    return right < field.size() ? (field[right] - field[center_row]) / physical_spacing
                                : std::numeric_limits<double>::quiet_NaN();
  }
  const std::uint32_t left = row_with_axis(coordinate - 1U);
  return left < field.size() ? (field[center_row] - field[left]) / physical_spacing
                             : std::numeric_limits<double>::quiet_NaN();
}

[[nodiscard]] inline bool starFormationPatchIsLeaf(
    const core::SimulationState& state,
    std::uint32_t patch_index,
    const std::unordered_set<std::uint64_t>& patch_ids_with_children) {
  if (patch_index >= state.patches.size()) {
    return true;
  }
  return patch_ids_with_children.find(state.patches.patch_id[patch_index]) ==
      patch_ids_with_children.end();
}

}  // namespace cosmosim::workflows::internal
