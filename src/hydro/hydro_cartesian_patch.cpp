#include "cosmosim/hydro/hydro_cartesian_patch.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace cosmosim::hydro {
namespace {

void validateSpec(const HydroCartesianPatchSpec& spec) {
  if (spec.nx == 0 || spec.ny == 0 || spec.nz == 0) {
    throw std::invalid_argument("Cartesian hydro patch dimensions must be positive");
  }
  if (spec.cell_width_x_comoving <= 0.0 ||
      spec.cell_width_y_comoving <= 0.0 ||
      spec.cell_width_z_comoving <= 0.0) {
    throw std::invalid_argument("Cartesian hydro patch cell widths must be positive");
  }
}

[[nodiscard]] std::size_t productOrThrow(std::size_t nx, std::size_t ny, std::size_t nz) {
  if (nx > std::numeric_limits<std::size_t>::max() / ny) {
    throw std::overflow_error("Cartesian hydro patch dimensions overflow");
  }
  const std::size_t xy = nx * ny;
  if (xy > std::numeric_limits<std::size_t>::max() / nz) {
    throw std::overflow_error("Cartesian hydro patch dimensions overflow");
  }
  return xy * nz;
}

}  // namespace

std::size_t HydroPatchGeometry::cellCount() const noexcept {
  return nx * ny * nz;
}

std::size_t HydroPatchGeometry::totalCellStorageCount() const noexcept {
  return cellCount() + ghost_cells.size();
}

std::size_t HydroPatchGeometry::linearCellIndex(std::size_t i, std::size_t j, std::size_t k) const {
  if (i >= nx || j >= ny || k >= nz) {
    throw std::out_of_range("HydroPatchGeometry::linearCellIndex index out of range");
  }
  return (k * ny + j) * nx + i;
}

std::array<std::size_t, 3> HydroPatchGeometry::cellIjk(std::size_t row) const {
  if (row >= cellCount() || nx == 0 || ny == 0) {
    throw std::out_of_range("HydroPatchGeometry::cellIjk row out of range");
  }
  const std::size_t i = row % nx;
  const std::size_t j = (row / nx) % ny;
  const std::size_t k = row / (nx * ny);
  return {i, j, k};
}

std::size_t HydroPatchGeometry::neighborCell(
    std::size_t row,
    int di,
    int dj,
    int dk) const noexcept {
  if (row >= cellCount() || nx == 0 || ny == 0 || nz == 0) {
    return k_invalid_cell_index;
  }
  const std::size_t i = row % nx;
  const std::size_t j = (row / nx) % ny;
  const std::size_t k = row / (nx * ny);
  const int ni = static_cast<int>(i) + di;
  const int nj = static_cast<int>(j) + dj;
  const int nk = static_cast<int>(k) + dk;
  if (ni < 0 || nj < 0 || nk < 0 ||
      ni >= static_cast<int>(nx) ||
      nj >= static_cast<int>(ny) ||
      nk >= static_cast<int>(nz)) {
    return k_invalid_cell_index;
  }
  return (static_cast<std::size_t>(nk) * ny + static_cast<std::size_t>(nj)) * nx +
      static_cast<std::size_t>(ni);
}

HydroPatchGeometry makeCartesianPatchGeometry(const HydroCartesianPatchSpec& spec) {
  validateSpec(spec);
  const std::size_t cell_count = productOrThrow(spec.nx, spec.ny, spec.nz);

  HydroPatchGeometry geometry;
  geometry.nx = spec.nx;
  geometry.ny = spec.ny;
  geometry.nz = spec.nz;
  geometry.origin_x_comoving = spec.origin_x_comoving;
  geometry.origin_y_comoving = spec.origin_y_comoving;
  geometry.origin_z_comoving = spec.origin_z_comoving;
  geometry.cell_width_x_comoving = spec.cell_width_x_comoving;
  geometry.cell_width_y_comoving = spec.cell_width_y_comoving;
  geometry.cell_width_z_comoving = spec.cell_width_z_comoving;
  geometry.cell_volume_comoving =
      spec.cell_width_x_comoving * spec.cell_width_y_comoving * spec.cell_width_z_comoving;

  const std::size_t x_faces = (spec.nx > 1) ? (spec.nx - 1U) * spec.ny * spec.nz : 0U;
  const std::size_t y_faces = (spec.ny > 1) ? spec.nx * (spec.ny - 1U) * spec.nz : 0U;
  const std::size_t z_faces = (spec.nz > 1) ? spec.nx * spec.ny * (spec.nz - 1U) : 0U;
  geometry.faces.reserve(x_faces + y_faces + z_faces);

  for (std::size_t k = 0; k < spec.nz; ++k) {
    for (std::size_t j = 0; j < spec.ny; ++j) {
      for (std::size_t i = 0; i + 1U < spec.nx; ++i) {
        const std::size_t owner = (k * spec.ny + j) * spec.nx + i;
        const std::size_t neighbor = owner + 1U;
        geometry.faces.push_back(HydroFace{
            .owner_cell = owner,
            .neighbor_cell = neighbor,
            .owner_minus_cell = geometry.neighborCell(owner, -1, 0, 0),
            .neighbor_plus_cell = geometry.neighborCell(neighbor, 1, 0, 0),
            .area_comoving = spec.cell_width_y_comoving * spec.cell_width_z_comoving,
            .normal_x = 1.0,
            .normal_y = 0.0,
            .normal_z = 0.0,
            .axis = HydroFaceAxis::kX});
      }
    }
  }

  for (std::size_t k = 0; k < spec.nz; ++k) {
    for (std::size_t j = 0; j + 1U < spec.ny; ++j) {
      for (std::size_t i = 0; i < spec.nx; ++i) {
        const std::size_t owner = (k * spec.ny + j) * spec.nx + i;
        const std::size_t neighbor = owner + spec.nx;
        geometry.faces.push_back(HydroFace{
            .owner_cell = owner,
            .neighbor_cell = neighbor,
            .owner_minus_cell = geometry.neighborCell(owner, 0, -1, 0),
            .neighbor_plus_cell = geometry.neighborCell(neighbor, 0, 1, 0),
            .area_comoving = spec.cell_width_x_comoving * spec.cell_width_z_comoving,
            .normal_x = 0.0,
            .normal_y = 1.0,
            .normal_z = 0.0,
            .axis = HydroFaceAxis::kY});
      }
    }
  }

  for (std::size_t k = 0; k + 1U < spec.nz; ++k) {
    for (std::size_t j = 0; j < spec.ny; ++j) {
      for (std::size_t i = 0; i < spec.nx; ++i) {
        const std::size_t owner = (k * spec.ny + j) * spec.nx + i;
        const std::size_t neighbor = owner + spec.nx * spec.ny;
        geometry.faces.push_back(HydroFace{
            .owner_cell = owner,
            .neighbor_cell = neighbor,
            .owner_minus_cell = geometry.neighborCell(owner, 0, 0, -1),
            .neighbor_plus_cell = geometry.neighborCell(neighbor, 0, 0, 1),
            .area_comoving = spec.cell_width_x_comoving * spec.cell_width_y_comoving,
            .normal_x = 0.0,
            .normal_y = 0.0,
            .normal_z = 1.0,
            .axis = HydroFaceAxis::kZ});
      }
    }
  }

  if (geometry.cellCount() != cell_count) {
    throw std::logic_error("Cartesian hydro patch cell count mismatch");
  }
  return geometry;
}

std::array<std::size_t, 3> chooseNearCubicCartesianFactors(std::size_t cell_count) {
  if (cell_count == 0) {
    return {0, 0, 0};
  }

  std::array<std::size_t, 3> best{cell_count, 1U, 1U};
  double best_score = std::numeric_limits<double>::infinity();
  for (std::size_t nz = 1; nz * nz * nz <= cell_count; ++nz) {
    if (cell_count % nz != 0) {
      continue;
    }
    const std::size_t rem = cell_count / nz;
    for (std::size_t ny = nz; ny * ny <= rem; ++ny) {
      if (rem % ny != 0) {
        continue;
      }
      const std::size_t nx = rem / ny;
      if (nx < ny) {
        continue;
      }
      const double score =
          static_cast<double>(nx - ny) * static_cast<double>(nx - ny) +
          static_cast<double>(ny - nz) * static_cast<double>(ny - nz) +
          static_cast<double>(nx - nz) * static_cast<double>(nx - nz);
      if (score < best_score) {
        best = {nx, ny, nz};
        best_score = score;
      }
    }
  }
  return best;
}

}  // namespace cosmosim::hydro
