#include "cosmosim/amr/amr_ghost_fill.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "cosmosim/hydro/hydro_boundary_conditions.hpp"

namespace cosmosim::amr {
namespace {

constexpr double k_geometry_tol = 1.0e-10;

struct CellSelection {
  const AmrHydroGhostFillPatch* patch = nullptr;
  std::vector<std::size_t> cells;
};

[[nodiscard]] double cellWidth(
    const hydro::HydroPatchGeometry& geometry,
    hydro::HydroFaceAxis axis) {
  switch (axis) {
    case hydro::HydroFaceAxis::kX:
      return geometry.cell_width_x_comoving;
    case hydro::HydroFaceAxis::kY:
      return geometry.cell_width_y_comoving;
    case hydro::HydroFaceAxis::kZ:
      return geometry.cell_width_z_comoving;
  }
  return 0.0;
}

[[nodiscard]] int axisIndex(hydro::HydroFaceAxis axis) {
  switch (axis) {
    case hydro::HydroFaceAxis::kX:
      return 0;
    case hydro::HydroFaceAxis::kY:
      return 1;
    case hydro::HydroFaceAxis::kZ:
      return 2;
  }
  return 0;
}

[[nodiscard]] double sideSign(hydro::HydroFaceSide side) {
  return side == hydro::HydroFaceSide::kLower ? -1.0 : 1.0;
}

[[nodiscard]] std::array<double, 3> patchCellWidths(const PatchDescriptor& patch) {
  return {
      patch.extent_comov[0] / static_cast<double>(patch.cell_dims[0]),
      patch.extent_comov[1] / static_cast<double>(patch.cell_dims[1]),
      patch.extent_comov[2] / static_cast<double>(patch.cell_dims[2]),
  };
}

[[nodiscard]] std::array<double, 3> cellCenter(
    const PatchDescriptor& patch,
    std::size_t cell) {
  const std::size_t nx = patch.cell_dims[0];
  const std::size_t ny = patch.cell_dims[1];
  const std::size_t i = cell % nx;
  const std::size_t j = (cell / nx) % ny;
  const std::size_t k = cell / (nx * ny);
  const auto widths = patchCellWidths(patch);
  return {
      patch.origin_comov[0] + (static_cast<double>(i) + 0.5) * widths[0],
      patch.origin_comov[1] + (static_cast<double>(j) + 0.5) * widths[1],
      patch.origin_comov[2] + (static_cast<double>(k) + 0.5) * widths[2],
  };
}

[[nodiscard]] bool containsPoint(
    const PatchDescriptor& patch,
    const std::array<double, 3>& point) {
  for (std::size_t axis = 0; axis < 3; ++axis) {
    const double lower = patch.origin_comov[axis] - k_geometry_tol;
    const double upper = patch.origin_comov[axis] + patch.extent_comov[axis] + k_geometry_tol;
    if (point[axis] < lower || point[axis] >= upper) {
      return false;
    }
  }
  return true;
}

[[nodiscard]] std::size_t cellContainingPoint(
    const PatchDescriptor& patch,
    const std::array<double, 3>& point) {
  const auto widths = patchCellWidths(patch);
  std::array<std::size_t, 3> ijk{};
  for (std::size_t axis = 0; axis < 3; ++axis) {
    double coordinate = (point[axis] - patch.origin_comov[axis]) / widths[axis];
    coordinate = std::clamp(
        coordinate,
        0.0,
        static_cast<double>(patch.cell_dims[axis]) - k_geometry_tol);
    ijk[axis] = static_cast<std::size_t>(coordinate);
    if (ijk[axis] >= patch.cell_dims[axis]) {
      ijk[axis] = patch.cell_dims[axis] - 1U;
    }
  }
  return ijk[0] + static_cast<std::size_t>(patch.cell_dims[0]) *
                      (ijk[1] + static_cast<std::size_t>(patch.cell_dims[1]) * ijk[2]);
}

[[nodiscard]] bool centerInsideGhostVolume(
    const PatchDescriptor& target_patch,
    const hydro::HydroGhostCell& ghost,
    const std::array<double, 3>& center) {
  const auto owner_center = cellCenter(target_patch, ghost.owner_real_cell);
  const auto target_widths = patchCellWidths(target_patch);
  const int normal_axis = axisIndex(ghost.axis);
  for (int axis = 0; axis < 3; ++axis) {
    double mid = owner_center[axis];
    if (axis == normal_axis) {
      mid += sideSign(ghost.side) * target_widths[axis];
    }
    const double half_width = 0.5 * target_widths[axis] + k_geometry_tol;
    if (center[axis] < mid - half_width || center[axis] >= mid + half_width) {
      return false;
    }
  }
  return true;
}

[[nodiscard]] std::array<double, 3> ghostProbePoint(
    const AmrHydroPatchGeometry& patch_geometry,
    const hydro::HydroGhostCell& ghost) {
  std::array<double, 3> point = cellCenter(patch_geometry.patch, ghost.owner_real_cell);
  point[axisIndex(ghost.axis)] += sideSign(ghost.side) * cellWidth(patch_geometry.geometry, ghost.axis);
  return point;
}

void validatePatchView(const AmrHydroGhostFillPatch& patch) {
  if (patch.geometry == nullptr || patch.conserved == nullptr) {
    throw std::invalid_argument("fillAmrHydroGhostCells: patch view must provide geometry and conserved storage");
  }
  if (patch.conserved->size() < patch.geometry->geometry.totalCellStorageCount()) {
    throw std::invalid_argument("fillAmrHydroGhostCells: conserved storage is missing AMR ghost rows");
  }
}

[[nodiscard]] bool staleRemotePatch(const AmrHydroGhostFillPatch& patch) {
  return patch.residency == AmrGhostSourceResidency::kRemoteReadOnly &&
      patch.ghost_hydro_epoch != patch.expected_ghost_hydro_epoch;
}

[[nodiscard]] CellSelection selectSourceCells(
    const AmrHydroGhostFillPatch& target,
    const hydro::HydroGhostCell& ghost,
    AmrHydroBoundaryClass boundary_class,
    std::span<AmrHydroGhostFillPatch> patches) {
  const auto probe = ghostProbePoint(*target.geometry, ghost);
  CellSelection selection;

  for (const AmrHydroGhostFillPatch& candidate : patches) {
    if (candidate.geometry == nullptr || candidate.conserved == nullptr) {
      continue;
    }
    const PatchDescriptor& source_patch = candidate.geometry->patch;
    const PatchDescriptor& target_patch = target.geometry->patch;
    if (source_patch.patch_id == target_patch.patch_id) {
      continue;
    }
    if (!containsPoint(source_patch, probe)) {
      continue;
    }

    const bool same_level =
        source_patch.level == target_patch.level &&
        boundary_class == AmrHydroBoundaryClass::kSameLevel;
    const bool coarse_to_fine =
        source_patch.level < target_patch.level &&
        boundary_class == AmrHydroBoundaryClass::kCoarseFine;
    const bool fine_to_coarse =
        source_patch.level > target_patch.level &&
        boundary_class == AmrHydroBoundaryClass::kCoarseFine;
    if (!same_level && !coarse_to_fine && !fine_to_coarse) {
      continue;
    }

    selection.patch = &candidate;
    if (same_level || coarse_to_fine) {
      selection.cells.push_back(cellContainingPoint(source_patch, probe));
      return selection;
    }

    for (std::size_t cell = 0; cell < candidate.geometry->geometry.cellCount(); ++cell) {
      if (centerInsideGhostVolume(target_patch, ghost, cellCenter(source_patch, cell))) {
        selection.cells.push_back(cell);
      }
    }
    if (!selection.cells.empty()) {
      return selection;
    }
  }

  return selection;
}

[[nodiscard]] hydro::HydroConservedState averageSourceState(const CellSelection& selection) {
  hydro::HydroConservedState state;
  for (const std::size_t cell : selection.cells) {
    state += selection.patch->conserved->loadCell(cell);
  }
  const double inv_count = 1.0 / static_cast<double>(selection.cells.size());
  return inv_count * state;
}

void fillAmrGhost(
    AmrHydroGhostFillPatch& target,
    AmrHydroGhostDescriptor& descriptor,
    std::span<AmrHydroGhostFillPatch> patches,
    AmrHydroGhostFillDiagnostics& diagnostics) {
  const hydro::HydroGhostCell& ghost =
      target.geometry->geometry.ghost_cells.at(descriptor.ghost_slot);
  const CellSelection selection =
      selectSourceCells(target, ghost, descriptor.boundary_class, patches);
  if (selection.patch == nullptr || selection.cells.empty()) {
    descriptor.fill_status = AmrHydroGhostFillStatus::kMissingSource;
    ++diagnostics.missing_source_records;
    ++diagnostics.unresolved_ghosts;
    throw std::runtime_error(
        "fillAmrHydroGhostCells: unresolved AMR hydro ghost for patch " +
        std::to_string(target.geometry->patch.patch_id));
  }
  if (staleRemotePatch(*selection.patch)) {
    descriptor.fill_status = AmrHydroGhostFillStatus::kRejectedStaleRemote;
    ++diagnostics.stale_epoch_rejections;
    ++diagnostics.unresolved_ghosts;
    throw std::runtime_error("fillAmrHydroGhostCells: stale source AMR hydro ghost epoch");
  }

  target.conserved->storeCell(descriptor.ghost_cell, averageSourceState(selection));
  if (descriptor.boundary_class == AmrHydroBoundaryClass::kSameLevel) {
    descriptor.fill_status = AmrHydroGhostFillStatus::kFilledSameLevel;
    ++diagnostics.same_level_ghosts_filled;
  } else if (selection.patch->geometry->patch.level < target.geometry->patch.level) {
    descriptor.fill_status = AmrHydroGhostFillStatus::kFilledCoarseToFine;
    ++diagnostics.coarse_to_fine_ghosts_filled;
  } else {
    descriptor.fill_status = AmrHydroGhostFillStatus::kFilledFineToCoarse;
    ++diagnostics.fine_to_coarse_ghosts_filled;
  }
}

void markStaleRemoteTarget(
    AmrHydroGhostFillPatch& patch,
    AmrHydroGhostFillDiagnostics& diagnostics) {
  for (AmrHydroGhostDescriptor& ghost : patch.geometry->ghosts) {
    ghost.fill_status = AmrHydroGhostFillStatus::kRejectedStaleRemote;
    ++diagnostics.stale_epoch_rejections;
    ++diagnostics.unresolved_ghosts;
  }
}

}  // namespace

AmrHydroGhostFillDiagnostics fillAmrHydroGhostCells(
    std::span<AmrHydroGhostFillPatch> patches,
    double adiabatic_index) {
  AmrHydroGhostFillDiagnostics diagnostics;
  for (const AmrHydroGhostFillPatch& patch : patches) {
    validatePatchView(patch);
  }

  for (AmrHydroGhostFillPatch& patch : patches) {
    if (staleRemotePatch(patch)) {
      markStaleRemoteTarget(patch, diagnostics);
      continue;
    }
    hydro::fillHydroBoundaryGhostCells(*patch.conserved, patch.geometry->geometry, adiabatic_index);
    for (AmrHydroGhostDescriptor& ghost : patch.geometry->ghosts) {
      if (ghost.boundary_class == AmrHydroBoundaryClass::kPhysical) {
        ghost.fill_status = AmrHydroGhostFillStatus::kFilledPhysicalBoundary;
        ++diagnostics.physical_ghosts_filled;
        continue;
      }
      fillAmrGhost(patch, ghost, patches, diagnostics);
    }
  }
  return diagnostics;
}

}  // namespace cosmosim::amr
