#include "cosmosim/amr/amr_hydro_geometry.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>

namespace cosmosim::amr {
namespace {

enum BoundarySlot : std::size_t {
  kLowerX = 0,
  kUpperX = 1,
  kLowerY = 2,
  kUpperY = 3,
  kLowerZ = 4,
  kUpperZ = 5,
};

[[nodiscard]] std::size_t checkedCellCount(const PatchDescriptor& patch) {
  const std::size_t nx = patch.cell_dims[0];
  const std::size_t ny = patch.cell_dims[1];
  const std::size_t nz = patch.cell_dims[2];
  if (nx == 0U || ny == 0U || nz == 0U) {
    throw std::invalid_argument("AMR hydro geometry requires positive patch cell dimensions");
  }
  if (nx > std::numeric_limits<std::size_t>::max() / ny) {
    throw std::overflow_error("AMR hydro patch dimensions overflow");
  }
  const std::size_t xy = nx * ny;
  if (xy > std::numeric_limits<std::size_t>::max() / nz) {
    throw std::overflow_error("AMR hydro patch dimensions overflow");
  }
  return xy * nz;
}

[[nodiscard]] double checkedCellWidth(double extent_comov, std::uint16_t cell_dim, const char* axis) {
  if (extent_comov <= 0.0) {
    throw std::invalid_argument(std::string("AMR hydro patch extent must be positive on ") + axis);
  }
  return extent_comov / static_cast<double>(cell_dim);
}

[[nodiscard]] hydro::HydroGhostMutationRights mutationRightsForBoundary(
    AmrHydroBoundaryClass boundary_class) {
  return boundary_class == AmrHydroBoundaryClass::kPhysical
      ? hydro::HydroGhostMutationRights::kWritablePhysicalBoundaryScratch
      : hydro::HydroGhostMutationRights::kReadOnlyImported;
}

[[nodiscard]] hydro::HydroBoundaryKind hydroBoundaryKindForBoundary(
    AmrHydroBoundaryClass boundary_class,
    hydro::HydroBoundaryKind physical_boundary_kind) {
  return boundary_class == AmrHydroBoundaryClass::kPhysical
      ? physical_boundary_kind
      : hydro::HydroBoundaryKind::kImportedMpi;
}

[[nodiscard]] AmrHydroGhostFillStatus fillStatusForBoundary(AmrHydroBoundaryClass boundary_class) {
  switch (boundary_class) {
    case AmrHydroBoundaryClass::kPhysical:
      return AmrHydroGhostFillStatus::kUnfilledPhysicalBoundary;
    case AmrHydroBoundaryClass::kSameLevel:
      return AmrHydroGhostFillStatus::kUnfilledSameLevel;
    case AmrHydroBoundaryClass::kCoarseFine:
      return AmrHydroGhostFillStatus::kUnfilledCoarseFine;
  }
  return AmrHydroGhostFillStatus::kUnfilledPhysicalBoundary;
}

void setNormal(hydro::HydroFace& face, hydro::HydroFaceAxis axis, hydro::HydroFaceSide side) {
  const double sign = side == hydro::HydroFaceSide::kLower ? -1.0 : 1.0;
  face.normal_x = axis == hydro::HydroFaceAxis::kX ? sign : 0.0;
  face.normal_y = axis == hydro::HydroFaceAxis::kY ? sign : 0.0;
  face.normal_z = axis == hydro::HydroFaceAxis::kZ ? sign : 0.0;
}

[[nodiscard]] double faceArea(const hydro::HydroPatchGeometry& geometry, hydro::HydroFaceAxis axis) {
  switch (axis) {
    case hydro::HydroFaceAxis::kX:
      return geometry.cell_width_y_comoving * geometry.cell_width_z_comoving;
    case hydro::HydroFaceAxis::kY:
      return geometry.cell_width_x_comoving * geometry.cell_width_z_comoving;
    case hydro::HydroFaceAxis::kZ:
      return geometry.cell_width_x_comoving * geometry.cell_width_y_comoving;
  }
  return 0.0;
}

[[nodiscard]] std::size_t boundaryOwner(
    const hydro::HydroPatchGeometry& geometry,
    hydro::HydroFaceAxis axis,
    hydro::HydroFaceSide side,
    std::size_t a,
    std::size_t b) {
  switch (axis) {
    case hydro::HydroFaceAxis::kX: {
      const std::size_t i = side == hydro::HydroFaceSide::kLower ? 0U : geometry.nx - 1U;
      return geometry.linearCellIndex(i, a, b);
    }
    case hydro::HydroFaceAxis::kY: {
      const std::size_t j = side == hydro::HydroFaceSide::kLower ? 0U : geometry.ny - 1U;
      return geometry.linearCellIndex(a, j, b);
    }
    case hydro::HydroFaceAxis::kZ: {
      const std::size_t k = side == hydro::HydroFaceSide::kLower ? 0U : geometry.nz - 1U;
      return geometry.linearCellIndex(a, b, k);
    }
  }
  return hydro::k_invalid_cell_index;
}

[[nodiscard]] BoundarySlot boundarySlot(hydro::HydroFaceAxis axis, hydro::HydroFaceSide side) {
  switch (axis) {
    case hydro::HydroFaceAxis::kX:
      return side == hydro::HydroFaceSide::kLower ? kLowerX : kUpperX;
    case hydro::HydroFaceAxis::kY:
      return side == hydro::HydroFaceSide::kLower ? kLowerY : kUpperY;
    case hydro::HydroFaceAxis::kZ:
      return side == hydro::HydroFaceSide::kLower ? kLowerZ : kUpperZ;
  }
  return kLowerX;
}

[[nodiscard]] std::uint64_t sourceGasCellIdForGhost(
    const AmrHydroPatchGeometry& patch_geometry,
    hydro::HydroFaceAxis axis,
    hydro::HydroFaceSide side,
    std::size_t a,
    std::size_t b,
    AmrHydroBoundaryClass boundary_class) {
  if (boundary_class != AmrHydroBoundaryClass::kPhysical) {
    return 0U;
  }

  const auto& geometry = patch_geometry.geometry;
  std::size_t source = hydro::k_invalid_cell_index;
  switch (axis) {
    case hydro::HydroFaceAxis::kX: {
      const std::size_t i = side == hydro::HydroFaceSide::kLower ? geometry.nx - 1U : 0U;
      source = geometry.linearCellIndex(i, a, b);
      break;
    }
    case hydro::HydroFaceAxis::kY: {
      const std::size_t j = side == hydro::HydroFaceSide::kLower ? geometry.ny - 1U : 0U;
      source = geometry.linearCellIndex(a, j, b);
      break;
    }
    case hydro::HydroFaceAxis::kZ: {
      const std::size_t k = side == hydro::HydroFaceSide::kLower ? geometry.nz - 1U : 0U;
      source = geometry.linearCellIndex(a, b, k);
      break;
    }
  }
  return patch_geometry.gas_cell_ids.at(source);
}

void appendBoundaryGhost(
    AmrHydroPatchGeometry& patch_geometry,
    const AmrHydroGeometryOptions& options,
    hydro::HydroFaceAxis axis,
    hydro::HydroFaceSide side,
    std::size_t a,
    std::size_t b) {
  auto& geometry = patch_geometry.geometry;
  const AmrHydroBoundaryClass boundary_class =
      options.boundary_classes[static_cast<std::size_t>(boundarySlot(axis, side))];
  const std::size_t ghost_slot = geometry.ghost_cells.size();
  const std::size_t ghost_cell = geometry.cellCount() + ghost_slot;
  const std::size_t owner = boundaryOwner(geometry, axis, side, a, b);
  const std::uint64_t owner_gas_cell_id = patch_geometry.gas_cell_ids.at(owner);
  const std::uint64_t source_gas_cell_id =
      sourceGasCellIdForGhost(patch_geometry, axis, side, a, b, boundary_class);

  geometry.ghost_cells.push_back(hydro::HydroGhostCell{
      .owner_real_cell = owner,
      .source_real_cell = owner,
      .ghost_cell = ghost_cell,
      .ghost_slot = ghost_slot,
      .boundary_kind = hydroBoundaryKindForBoundary(boundary_class, options.physical_boundary_kind),
      .axis = axis,
      .side = side,
      .mutation_rights = mutationRightsForBoundary(boundary_class)});

  hydro::HydroFace face{
      .owner_cell = owner,
      .neighbor_cell = ghost_cell,
      .owner_minus_cell = hydro::k_invalid_cell_index,
      .neighbor_plus_cell = hydro::k_invalid_cell_index,
      .ghost_cell_slot = ghost_slot,
      .area_comoving = faceArea(geometry, axis),
      .axis = axis};
  setNormal(face, axis, side);
  geometry.faces.push_back(face);

  patch_geometry.ghosts.push_back(AmrHydroGhostDescriptor{
      .patch_id = patch_geometry.patch.patch_id,
      .source_patch_id = 0U,
      .source_gas_cell_id = source_gas_cell_id,
      .owner_gas_cell_id = owner_gas_cell_id,
      .owner_local_cell_row = patch_geometry.local_cell_rows.at(owner),
      .ghost_slot = ghost_slot,
      .ghost_cell = ghost_cell,
      .boundary_class = boundary_class,
      .fill_status = fillStatusForBoundary(boundary_class)});
}

void appendAxisBoundaryGhosts(
    AmrHydroPatchGeometry& patch_geometry,
    const AmrHydroGeometryOptions& options,
    hydro::HydroFaceAxis axis) {
  const auto& geometry = patch_geometry.geometry;
  switch (axis) {
    case hydro::HydroFaceAxis::kX:
      for (std::size_t k = 0; k < geometry.nz; ++k) {
        for (std::size_t j = 0; j < geometry.ny; ++j) {
          appendBoundaryGhost(patch_geometry, options, axis, hydro::HydroFaceSide::kLower, j, k);
          appendBoundaryGhost(patch_geometry, options, axis, hydro::HydroFaceSide::kUpper, j, k);
        }
      }
      break;
    case hydro::HydroFaceAxis::kY:
      for (std::size_t k = 0; k < geometry.nz; ++k) {
        for (std::size_t i = 0; i < geometry.nx; ++i) {
          appendBoundaryGhost(patch_geometry, options, axis, hydro::HydroFaceSide::kLower, i, k);
          appendBoundaryGhost(patch_geometry, options, axis, hydro::HydroFaceSide::kUpper, i, k);
        }
      }
      break;
    case hydro::HydroFaceAxis::kZ:
      for (std::size_t j = 0; j < geometry.ny; ++j) {
        for (std::size_t i = 0; i < geometry.nx; ++i) {
          appendBoundaryGhost(patch_geometry, options, axis, hydro::HydroFaceSide::kLower, i, j);
          appendBoundaryGhost(patch_geometry, options, axis, hydro::HydroFaceSide::kUpper, i, j);
        }
      }
      break;
  }
}

[[nodiscard]] AmrHydroFaceDescriptor makeFaceDescriptor(
    const AmrHydroPatchGeometry& patch_geometry,
    std::size_t face_index) {
  const hydro::HydroFace& face = patch_geometry.geometry.faces.at(face_index);
  AmrHydroBoundaryClass boundary_class = AmrHydroBoundaryClass::kPhysical;
  std::uint64_t ghost_source_gas_cell_id = 0U;
  if (face.ghost_cell_slot != hydro::k_invalid_ghost_cell_slot) {
    const auto ghost_it = std::find_if(
        patch_geometry.ghosts.begin(),
        patch_geometry.ghosts.end(),
        [&face](const AmrHydroGhostDescriptor& ghost) {
          return ghost.ghost_slot == face.ghost_cell_slot;
        });
    if (ghost_it != patch_geometry.ghosts.end()) {
      boundary_class = ghost_it->boundary_class;
      ghost_source_gas_cell_id = ghost_it->source_gas_cell_id;
    }
  }

  return AmrHydroFaceDescriptor{
      .patch_id = patch_geometry.patch.patch_id,
      .face_id = (patch_geometry.patch.patch_id << 32U) | static_cast<std::uint64_t>(face_index),
      .owner_gas_cell_id = patch_geometry.gas_cell_ids.at(face.owner_cell),
      .neighbor_gas_cell_id = face.neighbor_cell < patch_geometry.gas_cell_ids.size()
          ? patch_geometry.gas_cell_ids[face.neighbor_cell]
          : 0U,
      .ghost_source_gas_cell_id = ghost_source_gas_cell_id,
      .owner_patch_cell = face.owner_cell,
      .neighbor_patch_cell = face.neighbor_cell,
      .ghost_slot = face.ghost_cell_slot,
      .area_comoving = face.area_comoving,
      .normal_x = face.normal_x,
      .normal_y = face.normal_y,
      .normal_z = face.normal_z,
      .axis = face.axis,
      .boundary_class = boundary_class};
}

[[nodiscard]] double positiveDensityForCell(const core::SimulationState& state, std::uint32_t row) {
  const double density = state.gas_cells.density_code.at(row);
  if (density > 0.0) {
    return density;
  }
  const double volume = 1.0;
  return std::max(state.cells.mass_code.at(row) / volume, 1.0e-14);
}

[[nodiscard]] std::size_t linearPatchCellIndex(
    std::array<std::uint16_t, 3> dims,
    std::size_t i,
    std::size_t j,
    std::size_t k) {
  return i + static_cast<std::size_t>(dims[0]) * (j + static_cast<std::size_t>(dims[1]) * k);
}

[[nodiscard]] std::size_t patchLocalCellForRow(
    const core::SimulationState& state,
    const PatchDescriptor& patch,
    std::uint32_t row) {
  const std::array<double, 3> point{
      state.cells.center_x_comoving.at(row),
      state.cells.center_y_comoving.at(row),
      state.cells.center_z_comoving.at(row)};
  std::array<std::size_t, 3> ijk{};
  for (std::size_t axis = 0; axis < 3; ++axis) {
    const double extent = patch.extent_comov[axis];
    const std::uint16_t dim = patch.cell_dims[axis];
    if (extent <= 0.0 || dim == 0U) {
      throw std::invalid_argument("buildAmrHydroPatchGeometry: invalid explicit patch geometry");
    }
    const double lower = patch.origin_comov[axis];
    const double upper = patch.origin_comov[axis] + extent;
    const double scale = std::max({1.0, std::abs(lower), std::abs(upper), std::abs(point[axis])});
    if (point[axis] < lower - 1.0e-10 * scale || point[axis] >= upper + 1.0e-10 * scale) {
      throw std::runtime_error("buildAmrHydroPatchGeometry: gas-cell center lies outside explicit patch geometry");
    }
    double coord = (point[axis] - lower) / (extent / static_cast<double>(dim));
    coord = std::clamp(coord, 0.0, static_cast<double>(dim) - 1.0e-12);
    ijk[axis] = static_cast<std::size_t>(coord);
    if (ijk[axis] >= dim) {
      ijk[axis] = dim - 1U;
    }
  }
  return linearPatchCellIndex(patch.cell_dims, ijk[0], ijk[1], ijk[2]);
}

}  // namespace

std::span<const std::uint64_t> AmrHydroPatchGeometry::gasCellIds() const noexcept {
  return gas_cell_ids;
}

std::vector<std::size_t> AmrHydroPatchGeometry::internalFaceIndices() const {
  std::vector<std::size_t> indices;
  indices.reserve(faces.size());
  for (std::size_t face_index = 0; face_index < faces.size(); ++face_index) {
    if (faces[face_index].ghost_slot == hydro::k_invalid_ghost_cell_slot) {
      indices.push_back(face_index);
    }
  }
  return indices;
}

AmrHydroPatchGeometry buildAmrHydroPatchGeometry(
    const core::SimulationState& state,
    const PatchDescriptor& patch,
    const AmrHydroGeometryOptions& options) {
  const std::size_t expected_cells = checkedCellCount(patch);
  state.requireGasCellIdentityMapCoversDenseRows("buildAmrHydroPatchGeometry");
  if (!state.gasCellIdentityMapMatchesSidecarLanes()) {
    throw std::runtime_error("buildAmrHydroPatchGeometry: gas-cell identity sidecar mirrors are stale");
  }

  std::vector<std::uint32_t> rows = state.gas_cell_identity.rowsForPatch(patch.patch_id);
  if (rows.size() != expected_cells) {
    throw std::runtime_error("buildAmrHydroPatchGeometry: patch gas-cell coverage does not match PatchDescriptor cell_dims");
  }
  std::vector<std::uint32_t> row_by_patch_cell(expected_cells, core::kInvalidGasCellRow);
  for (const std::uint32_t row : rows) {
    const std::size_t patch_cell = patchLocalCellForRow(state, patch, row);
    if (patch_cell >= expected_cells || row_by_patch_cell[patch_cell] != core::kInvalidGasCellRow) {
      throw std::runtime_error("buildAmrHydroPatchGeometry: duplicate or out-of-range patch-local gas-cell mapping");
    }
    row_by_patch_cell[patch_cell] = row;
  }
  for (const std::uint32_t row : row_by_patch_cell) {
    if (row == core::kInvalidGasCellRow) {
      throw std::runtime_error("buildAmrHydroPatchGeometry: explicit patch geometry has a missing patch-local gas cell");
    }
  }

  AmrHydroPatchGeometry result;
  result.patch = patch;
  result.source_gas_cell_identity_generation = state.gasCellIdentityGeneration();
  result.geometry.nx = patch.cell_dims[0];
  result.geometry.ny = patch.cell_dims[1];
  result.geometry.nz = patch.cell_dims[2];
  result.geometry.origin_x_comoving = patch.origin_comov[0];
  result.geometry.origin_y_comoving = patch.origin_comov[1];
  result.geometry.origin_z_comoving = patch.origin_comov[2];
  result.geometry.cell_width_x_comoving = checkedCellWidth(patch.extent_comov[0], patch.cell_dims[0], "x");
  result.geometry.cell_width_y_comoving = checkedCellWidth(patch.extent_comov[1], patch.cell_dims[1], "y");
  result.geometry.cell_width_z_comoving = checkedCellWidth(patch.extent_comov[2], patch.cell_dims[2], "z");
  result.geometry.cell_volume_comoving =
      result.geometry.cell_width_x_comoving *
      result.geometry.cell_width_y_comoving *
      result.geometry.cell_width_z_comoving;

  result.real_cells.reserve(row_by_patch_cell.size());
  result.gas_cell_ids.reserve(row_by_patch_cell.size());
  result.local_cell_rows.reserve(row_by_patch_cell.size());
  for (std::size_t patch_cell = 0; patch_cell < row_by_patch_cell.size(); ++patch_cell) {
    const std::uint32_t row = row_by_patch_cell[patch_cell];
    const auto* record = state.gas_cell_identity.findByLocalRow(row);
    if (record == nullptr || record->owning_patch_id != patch.patch_id) {
      throw std::runtime_error("buildAmrHydroPatchGeometry: identity row does not belong to requested patch");
    }
    if (row >= state.cells.size() || row >= state.gas_cells.size()) {
      throw std::runtime_error("buildAmrHydroPatchGeometry: identity row is outside gas-cell storage");
    }
    result.real_cells.push_back(AmrHydroCellDescriptor{
        .patch_id = patch.patch_id,
        .gas_cell_id = record->gas_cell_id,
        .local_cell_row = row,
        .patch_local_cell = patch_cell});
    result.gas_cell_ids.push_back(record->gas_cell_id);
    result.local_cell_rows.push_back(row);
  }

  const auto append_internal_face =
      [&result](std::size_t owner, std::size_t neighbor, hydro::HydroFaceAxis axis) {
        hydro::HydroFace face{
            .owner_cell = owner,
            .neighbor_cell = neighbor,
            .owner_minus_cell = hydro::k_invalid_cell_index,
            .neighbor_plus_cell = hydro::k_invalid_cell_index,
            .area_comoving = faceArea(result.geometry, axis),
            .axis = axis};
        switch (axis) {
          case hydro::HydroFaceAxis::kX:
            face.owner_minus_cell = result.geometry.neighborCell(owner, -1, 0, 0);
            face.neighbor_plus_cell = result.geometry.neighborCell(neighbor, 1, 0, 0);
            face.normal_x = 1.0;
            break;
          case hydro::HydroFaceAxis::kY:
            face.owner_minus_cell = result.geometry.neighborCell(owner, 0, -1, 0);
            face.neighbor_plus_cell = result.geometry.neighborCell(neighbor, 0, 1, 0);
            face.normal_y = 1.0;
            break;
          case hydro::HydroFaceAxis::kZ:
            face.owner_minus_cell = result.geometry.neighborCell(owner, 0, 0, -1);
            face.neighbor_plus_cell = result.geometry.neighborCell(neighbor, 0, 0, 1);
            face.normal_z = 1.0;
            break;
        }
        result.geometry.faces.push_back(face);
      };

  for (std::size_t k = 0; k < result.geometry.nz; ++k) {
    for (std::size_t j = 0; j < result.geometry.ny; ++j) {
      for (std::size_t i = 0; i + 1U < result.geometry.nx; ++i) {
        const std::size_t owner = result.geometry.linearCellIndex(i, j, k);
        append_internal_face(owner, owner + 1U, hydro::HydroFaceAxis::kX);
      }
    }
  }
  for (std::size_t k = 0; k < result.geometry.nz; ++k) {
    for (std::size_t j = 0; j + 1U < result.geometry.ny; ++j) {
      for (std::size_t i = 0; i < result.geometry.nx; ++i) {
        const std::size_t owner = result.geometry.linearCellIndex(i, j, k);
        append_internal_face(owner, owner + result.geometry.nx, hydro::HydroFaceAxis::kY);
      }
    }
  }
  for (std::size_t k = 0; k + 1U < result.geometry.nz; ++k) {
    for (std::size_t j = 0; j < result.geometry.ny; ++j) {
      for (std::size_t i = 0; i < result.geometry.nx; ++i) {
        const std::size_t owner = result.geometry.linearCellIndex(i, j, k);
        append_internal_face(owner, owner + result.geometry.nx * result.geometry.ny, hydro::HydroFaceAxis::kZ);
      }
    }
  }

  appendAxisBoundaryGhosts(result, options, hydro::HydroFaceAxis::kX);
  appendAxisBoundaryGhosts(result, options, hydro::HydroFaceAxis::kY);
  appendAxisBoundaryGhosts(result, options, hydro::HydroFaceAxis::kZ);

  result.faces.reserve(result.geometry.faces.size());
  for (std::size_t face_index = 0; face_index < result.geometry.faces.size(); ++face_index) {
    result.faces.push_back(makeFaceDescriptor(result, face_index));
  }
  return result;
}

hydro::HydroConservedStateSoa loadAmrHydroConservedState(
    const core::SimulationState& state,
    const AmrHydroPatchGeometry& patch_geometry,
    double adiabatic_index) {
  state.requireGasCellIdentityMapFresh(
      patch_geometry.source_gas_cell_identity_generation,
      "loadAmrHydroConservedState");
  hydro::HydroConservedStateSoa conserved(patch_geometry.geometry.totalCellStorageCount());
  for (std::size_t patch_cell = 0; patch_cell < patch_geometry.real_cells.size(); ++patch_cell) {
    const std::uint32_t row = patch_geometry.real_cells[patch_cell].local_cell_row;
    const hydro::HydroPrimitiveState primitive{
        .rho_comoving = positiveDensityForCell(state, row),
        .vel_x_peculiar = state.gas_cells.velocity_x_peculiar.at(row),
        .vel_y_peculiar = state.gas_cells.velocity_y_peculiar.at(row),
        .vel_z_peculiar = state.gas_cells.velocity_z_peculiar.at(row),
        .pressure_comoving = std::max(state.gas_cells.pressure_code.at(row), 1.0e-14)};
    conserved.storeCell(
        patch_cell,
        hydro::HydroCoreSolver::conservedFromPrimitive(primitive, adiabatic_index));
  }
  return conserved;
}

}  // namespace cosmosim::amr
