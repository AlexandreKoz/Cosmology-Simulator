#include "cosmosim/amr/amr_ghost_fill.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "cosmosim/hydro/hydro_boundary_conditions.hpp"

namespace cosmosim::amr {
namespace {

constexpr double k_geometry_tol = 1.0e-10;
constexpr double k_time_tol = 1.0e-12;

struct CellSelection {
  const AmrHydroGhostFillPatch* patch = nullptr;
  const AmrHydroSparseGhostSource* sparse_remote = nullptr;
  std::vector<std::size_t> cells;
};

[[nodiscard]] const PatchDescriptor& selectionPatch(const CellSelection& selection) {
  if (selection.patch != nullptr) {
    return selection.patch->geometry->patch;
  }
  if (selection.sparse_remote != nullptr && selection.sparse_remote->patch != nullptr) {
    return *selection.sparse_remote->patch;
  }
  throw std::logic_error("fillAmrHydroGhostCells: empty source selection");
}

[[nodiscard]] const AmrHydroSparseRemoteCell* findSparseRemoteCell(
    const AmrHydroSparseGhostSource& source,
    std::size_t patch_local_cell) {
  const auto it = std::lower_bound(
      source.cells.begin(), source.cells.end(), patch_local_cell,
      [](const AmrHydroSparseRemoteCell& cell, std::size_t offset) {
        return static_cast<std::size_t>(cell.patch_local_cell) < offset;
      });
  if (it == source.cells.end() ||
      static_cast<std::size_t>(it->patch_local_cell) != patch_local_cell) {
    return nullptr;
  }
  return &*it;
}

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
  if (!patch.available_real_cells.empty() &&
      patch.available_real_cells.size() != patch.geometry->geometry.cellCount()) {
    throw std::invalid_argument("fillAmrHydroGhostCells: sparse source availability does not match real-cell count");
  }
  if (!std::isfinite(patch.target_state_time_code) || !std::isfinite(patch.ghost_fill_time_code) ||
      !std::isfinite(patch.source_current_state_time_code)) {
    throw std::invalid_argument("fillAmrHydroGhostCells: ghost-fill times must be finite");
  }
}

[[nodiscard]] bool staleRemotePatch(const AmrHydroGhostFillPatch& patch) {
  return patch.residency == AmrGhostSourceResidency::kRemoteReadOnly &&
      patch.ghost_hydro_epoch != patch.expected_ghost_hydro_epoch;
}

[[nodiscard]] bool timesMatch(double lhs, double rhs) {
  const double scale = std::max({1.0, std::abs(lhs), std::abs(rhs)});
  return std::abs(lhs - rhs) <= k_time_tol * scale;
}

[[nodiscard]] bool sourceCellAvailable(
    const AmrHydroGhostFillPatch& patch,
    std::size_t cell) {
  return patch.available_real_cells.empty() ||
      (cell < patch.available_real_cells.size() && patch.available_real_cells[cell] != 0U);
}

void requireSourceCellAvailable(
    const AmrHydroGhostFillPatch& patch,
    std::size_t cell) {
  if (!sourceCellAvailable(patch, cell)) {
    throw std::runtime_error(
        "fillAmrHydroGhostCells: sparse remote interface payload is missing a required source cell");
  }
}

[[nodiscard]] CellSelection selectSourceCells(
    const AmrHydroGhostFillPatch& target,
    const hydro::HydroGhostCell& ghost,
    AmrHydroBoundaryClass boundary_class,
    std::span<AmrHydroGhostFillPatch> patches,
    std::span<const AmrHydroSparseGhostSource> remote_sources) {
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
      const std::size_t source_cell = cellContainingPoint(source_patch, probe);
      requireSourceCellAvailable(candidate, source_cell);
      selection.cells.push_back(source_cell);
      return selection;
    }

    for (std::size_t cell = 0; cell < candidate.geometry->geometry.cellCount(); ++cell) {
      if (centerInsideGhostVolume(target_patch, ghost, cellCenter(source_patch, cell))) {
        requireSourceCellAvailable(candidate, cell);
        selection.cells.push_back(cell);
      }
    }
    if (!selection.cells.empty()) {
      return selection;
    }
  }

  for (const AmrHydroSparseGhostSource& candidate : remote_sources) {
    if (candidate.patch == nullptr) {
      continue;
    }
    const PatchDescriptor& source_patch = *candidate.patch;
    const PatchDescriptor& target_patch = target.geometry->patch;
    if (source_patch.patch_id == target_patch.patch_id ||
        !containsPoint(source_patch, probe)) {
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
    selection.sparse_remote = &candidate;
    if (same_level || coarse_to_fine) {
      const std::size_t source_cell = cellContainingPoint(source_patch, probe);
      if (findSparseRemoteCell(candidate, source_cell) == nullptr) {
        throw std::runtime_error(
            "fillAmrHydroGhostCells: sparse remote interface payload is missing a required source cell");
      }
      selection.cells.push_back(source_cell);
      return selection;
    }
    for (const AmrHydroSparseRemoteCell& cell : candidate.cells) {
      if (centerInsideGhostVolume(
              target_patch,
              ghost,
              cellCenter(source_patch, cell.patch_local_cell))) {
        selection.cells.push_back(cell.patch_local_cell);
      }
    }
    if (!selection.cells.empty()) {
      return selection;
    }
    selection.sparse_remote = nullptr;
  }

  return selection;
}

[[nodiscard]] hydro::HydroConservedState averageSourceState(const CellSelection& selection) {
  hydro::HydroConservedState state;
  for (const std::size_t cell : selection.cells) {
    if (selection.patch != nullptr) {
      state += selection.patch->conserved->loadCell(cell);
    } else {
      const AmrHydroSparseRemoteCell* remote_cell =
          findSparseRemoteCell(*selection.sparse_remote, cell);
      if (remote_cell == nullptr) {
        throw std::runtime_error(
            "fillAmrHydroGhostCells: sparse remote source disappeared during ghost averaging");
      }
      state += remote_cell->conserved;
    }
  }
  const double inv_count = 1.0 / static_cast<double>(selection.cells.size());
  return inv_count * state;
}

[[nodiscard]] bool selectionIsStale(const CellSelection& selection) {
  if (selection.patch != nullptr) {
    return staleRemotePatch(*selection.patch);
  }
  return selection.sparse_remote != nullptr &&
      selection.sparse_remote->ghost_hydro_epoch !=
          selection.sparse_remote->expected_ghost_hydro_epoch;
}

[[nodiscard]] double selectionSourceTime(const CellSelection& selection) {
  return selection.patch != nullptr
      ? selection.patch->source_current_state_time_code
      : selection.sparse_remote->source_current_state_time_code;
}

[[nodiscard]] bool selectionIsSparseRemote(const CellSelection& selection) noexcept {
  return selection.sparse_remote != nullptr;
}

[[nodiscard]] bool finiteConserved(const hydro::HydroConservedState& state) {
  return std::isfinite(state.mass_density_comoving) &&
      std::isfinite(state.momentum_density_x_comoving) &&
      std::isfinite(state.momentum_density_y_comoving) &&
      std::isfinite(state.momentum_density_z_comoving) &&
      std::isfinite(state.total_energy_density_comoving) &&
      std::isfinite(state.metal_mass_density_comoving);
}

[[nodiscard]] hydro::HydroConservedState temporalStateForCell(
    const AmrHydroGhostFillPatch& source,
    std::size_t source_cell,
    double fill_time_code,
    double adiabatic_index,
    bool& endpoint_fill,
    AmrHydroGhostFillDiagnostics& diagnostics) {
  if (source.temporal_boundary_history == nullptr) {
    ++diagnostics.temporal_history_missing_rejections;
    throw std::runtime_error("fillAmrHydroGhostCells: coarse-to-fine temporal history is required but unavailable");
  }
  const auto* record = source.temporal_boundary_history->findByPatchId(source.geometry->patch.patch_id);
  if (record == nullptr) {
    ++diagnostics.temporal_history_missing_rejections;
    throw std::runtime_error("fillAmrHydroGhostCells: no temporal history record for coarse source patch");
  }
  if (!record->end_state_valid || record->patch_level != source.geometry->patch.level) {
    ++diagnostics.temporal_history_invalid_rejections;
    throw std::runtime_error("fillAmrHydroGhostCells: coarse temporal history is incomplete or level-invalid");
  }
  if (record->patch_geometry_fingerprint != amrPatchGeometryFingerprint(source.geometry->patch)) {
    ++diagnostics.temporal_geometry_mismatch_rejections;
    throw std::runtime_error("fillAmrHydroGhostCells: coarse temporal history geometry fingerprint is stale");
  }
  if (record->gas_cell_identity_generation != source.geometry->source_gas_cell_identity_generation) {
    ++diagnostics.temporal_identity_mismatch_rejections;
    throw std::runtime_error("fillAmrHydroGhostCells: coarse temporal history gas identity generation is stale");
  }
  const double duration = record->interval_end_code - record->interval_start_code;
  if (!std::isfinite(duration) || duration <= 0.0 || !std::isfinite(fill_time_code)) {
    ++diagnostics.temporal_history_invalid_rejections;
    throw std::runtime_error("fillAmrHydroGhostCells: temporal history interval is invalid");
  }
  double alpha = (fill_time_code - record->interval_start_code) / duration;
  const double endpoint_tol = k_time_tol * std::max({1.0, std::abs(record->interval_start_code),
      std::abs(record->interval_end_code), std::abs(fill_time_code)});
  if (alpha < -endpoint_tol || alpha > 1.0 + endpoint_tol) {
    ++diagnostics.temporal_time_out_of_range_rejections;
    throw std::runtime_error("fillAmrHydroGhostCells: requested coarse temporal ghost lies outside validated interval");
  }
  if (alpha <= endpoint_tol) {
    alpha = 0.0;
    endpoint_fill = true;
  } else if (alpha >= 1.0 - endpoint_tol) {
    alpha = 1.0;
    endpoint_fill = true;
  }
  if (source_cell >= source.geometry->gas_cell_ids.size()) {
    ++diagnostics.temporal_history_invalid_rejections;
    throw std::runtime_error("fillAmrHydroGhostCells: temporal source cell is outside source geometry");
  }
  const std::uint64_t gas_cell_id = source.geometry->gas_cell_ids[source_cell];
  const auto cell_it = std::find_if(
      record->cells.begin(), record->cells.end(),
      [gas_cell_id, source_cell](const core::AmrTemporalBoundaryHistoryCellRecord& cell) {
        return cell.gas_cell_id == gas_cell_id && cell.patch_local_cell == source_cell;
      });
  if (cell_it == record->cells.end()) {
    ++diagnostics.temporal_identity_mismatch_rejections;
    throw std::runtime_error("fillAmrHydroGhostCells: temporal history cannot resolve coarse source by stable gas_cell_id");
  }
  const hydro::HydroConservedState start{
      .mass_density_comoving = cell_it->start_mass_density_comoving,
      .momentum_density_x_comoving = cell_it->start_momentum_density_x_comoving,
      .momentum_density_y_comoving = cell_it->start_momentum_density_y_comoving,
      .momentum_density_z_comoving = cell_it->start_momentum_density_z_comoving,
      .total_energy_density_comoving = cell_it->start_total_energy_density_comoving,
      .metal_mass_density_comoving = cell_it->start_metal_mass_density_comoving};
  const hydro::HydroConservedState end{
      .mass_density_comoving = cell_it->end_mass_density_comoving,
      .momentum_density_x_comoving = cell_it->end_momentum_density_x_comoving,
      .momentum_density_y_comoving = cell_it->end_momentum_density_y_comoving,
      .momentum_density_z_comoving = cell_it->end_momentum_density_z_comoving,
      .total_energy_density_comoving = cell_it->end_total_energy_density_comoving,
      .metal_mass_density_comoving = cell_it->end_metal_mass_density_comoving};
  const hydro::HydroConservedState result = (1.0 - alpha) * start + alpha * end;
  if (!finiteConserved(result) || result.mass_density_comoving <= 0.0) {
    ++diagnostics.temporal_history_invalid_rejections;
    throw std::runtime_error("fillAmrHydroGhostCells: temporally interpolated conserved state is non-finite or non-positive");
  }
  const hydro::HydroPrimitiveState primitive = hydro::HydroCoreSolver::primitiveFromConserved(result, adiabatic_index);
  if (!std::isfinite(primitive.rho_comoving) || !std::isfinite(primitive.pressure_comoving) ||
      !std::isfinite(primitive.vel_x_peculiar) || !std::isfinite(primitive.vel_y_peculiar) ||
      !std::isfinite(primitive.vel_z_peculiar) || primitive.rho_comoving <= 0.0 ||
      primitive.pressure_comoving <= 0.0) {
    ++diagnostics.temporal_history_invalid_rejections;
    throw std::runtime_error("fillAmrHydroGhostCells: primitive recovery of temporal conserved state is inadmissible");
  }
  return result;
}

[[noreturn]] void rejectTemporal(
    AmrHydroGhostDescriptor& descriptor,
    AmrHydroGhostFillStatus status,
    std::size_t& counter,
    AmrHydroGhostFillDiagnostics& diagnostics,
    std::string_view message) {
  descriptor.fill_status = status;
  ++counter;
  ++diagnostics.unresolved_ghosts;
  throw std::runtime_error(std::string("fillAmrHydroGhostCells: ") + std::string(message));
}

void fillAmrGhost(
    AmrHydroGhostFillPatch& target,
    AmrHydroGhostDescriptor& descriptor,
    std::span<AmrHydroGhostFillPatch> patches,
    std::span<const AmrHydroSparseGhostSource> remote_sources,
    double adiabatic_index,
    AmrHydroGhostFillDiagnostics& diagnostics) {
  const hydro::HydroGhostCell& ghost =
      target.geometry->geometry.ghost_cells.at(descriptor.ghost_slot);
  const CellSelection selection =
      selectSourceCells(target, ghost, descriptor.boundary_class, patches, remote_sources);
  if ((selection.patch == nullptr && selection.sparse_remote == nullptr) ||
      selection.cells.empty()) {
    descriptor.fill_status = AmrHydroGhostFillStatus::kMissingSource;
    ++diagnostics.missing_source_records;
    ++diagnostics.unresolved_ghosts;
    throw std::runtime_error(
        "fillAmrHydroGhostCells: unresolved AMR hydro ghost for patch " +
        std::to_string(target.geometry->patch.patch_id));
  }
  if (selectionIsStale(selection)) {
    descriptor.fill_status = AmrHydroGhostFillStatus::kRejectedStaleRemote;
    ++diagnostics.stale_epoch_rejections;
    ++diagnostics.unresolved_ghosts;
    throw std::runtime_error("fillAmrHydroGhostCells: stale source AMR hydro ghost epoch");
  }

  const PatchDescriptor& source_patch = selectionPatch(selection);
  const bool same_level = source_patch.level == target.geometry->patch.level;
  const bool coarse_to_fine = source_patch.level < target.geometry->patch.level;
  const bool fine_to_coarse = source_patch.level > target.geometry->patch.level;
  if (same_level) {
    if (!timesMatch(selectionSourceTime(selection), target.ghost_fill_time_code)) {
      rejectTemporal(
          descriptor,
          AmrHydroGhostFillStatus::kRejectedTemporalSameLevelMismatch,
          diagnostics.temporal_same_level_mismatch_rejections,
          diagnostics,
          "same-level source state time does not match requested ghost-fill time");
    }
    target.conserved->storeCell(descriptor.ghost_cell, averageSourceState(selection));
    descriptor.fill_status = AmrHydroGhostFillStatus::kFilledSameLevel;
    ++diagnostics.same_level_ghosts_filled;
    return;
  }
  if (fine_to_coarse) {
    if (!timesMatch(selectionSourceTime(selection), target.ghost_fill_time_code)) {
      rejectTemporal(
          descriptor,
          AmrHydroGhostFillStatus::kRejectedTemporalFineToCoarseMismatch,
          diagnostics.temporal_fine_to_coarse_misalignment_rejections,
          diagnostics,
          "fine-to-coarse ghost access is not at a synchronization time");
    }
    target.conserved->storeCell(descriptor.ghost_cell, averageSourceState(selection));
    descriptor.fill_status = AmrHydroGhostFillStatus::kFilledFineToCoarse;
    ++diagnostics.fine_to_coarse_ghosts_filled;
    return;
  }

  if (coarse_to_fine && target.enable_temporal_coarse_to_fine) {
    if (selectionIsSparseRemote(selection) ||
        (selection.patch->residency == AmrGhostSourceResidency::kRemoteReadOnly &&
         selection.patch->temporal_boundary_history == nullptr)) {
      rejectTemporal(
          descriptor,
          AmrHydroGhostFillStatus::kRejectedTemporalHistoryMissing,
          diagnostics.temporal_history_missing_rejections,
          diagnostics,
          "remote temporal interpolation requires an explicitly cached local temporal history");
    }
    hydro::HydroConservedState accumulated;
    bool endpoint_fill = false;
    const std::size_t missing_before = diagnostics.temporal_history_missing_rejections;
    const std::size_t invalid_before = diagnostics.temporal_history_invalid_rejections;
    const std::size_t out_of_range_before = diagnostics.temporal_time_out_of_range_rejections;
    const std::size_t geometry_before = diagnostics.temporal_geometry_mismatch_rejections;
    const std::size_t identity_before = diagnostics.temporal_identity_mismatch_rejections;
    try {
      for (const std::size_t cell : selection.cells) {
        bool cell_endpoint = false;
        accumulated += temporalStateForCell(
            *selection.patch,
            cell,
            target.ghost_fill_time_code,
            adiabatic_index,
            cell_endpoint,
            diagnostics);
        endpoint_fill = endpoint_fill || cell_endpoint;
      }
    } catch (const std::runtime_error&) {
      if (diagnostics.temporal_time_out_of_range_rejections > out_of_range_before) {
        descriptor.fill_status = AmrHydroGhostFillStatus::kRejectedTemporalTimeOutOfRange;
      } else if (diagnostics.temporal_history_missing_rejections > missing_before) {
        descriptor.fill_status = AmrHydroGhostFillStatus::kRejectedTemporalHistoryMissing;
      } else if (diagnostics.temporal_history_invalid_rejections > invalid_before ||
                 diagnostics.temporal_geometry_mismatch_rejections > geometry_before ||
                 diagnostics.temporal_identity_mismatch_rejections > identity_before) {
        descriptor.fill_status = AmrHydroGhostFillStatus::kRejectedTemporalHistoryInvalid;
      } else {
        descriptor.fill_status = AmrHydroGhostFillStatus::kRejectedTemporalHistoryInvalid;
      }
      ++diagnostics.unresolved_ghosts;
      throw;
    }
    target.conserved->storeCell(
        descriptor.ghost_cell,
        (1.0 / static_cast<double>(selection.cells.size())) * accumulated);
    descriptor.fill_status = AmrHydroGhostFillStatus::kFilledCoarseToFineTemporal;
    ++diagnostics.coarse_to_fine_ghosts_filled;
    ++diagnostics.temporal_coarse_to_fine_ghosts_filled;
    if (endpoint_fill) {
      ++diagnostics.temporal_endpoint_ghosts_filled;
    }
    return;
  }

  if (coarse_to_fine && !timesMatch(selectionSourceTime(selection), target.ghost_fill_time_code)) {
    rejectTemporal(
        descriptor,
        AmrHydroGhostFillStatus::kRejectedTemporalHistoryMissing,
        diagnostics.temporal_history_missing_rejections,
        diagnostics,
        "coarse-to-fine ghost fill is time-misaligned and lacks temporal interpolation");
  }
  target.conserved->storeCell(descriptor.ghost_cell, averageSourceState(selection));
  descriptor.fill_status = AmrHydroGhostFillStatus::kFilledCoarseToFine;
  ++diagnostics.coarse_to_fine_ghosts_filled;
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

[[nodiscard]] hydro::HydroConservedState conservedForStateRow(
    const core::SimulationState& state,
    std::uint32_t row,
    double adiabatic_index) {
  if (row >= state.cells.size() || row >= state.gas_cells.size()) {
    throw std::runtime_error("AMR temporal history source row is outside SimulationState");
  }
  const hydro::HydroPrimitiveState primitive{
      .rho_comoving = state.gas_cells.density_code[row],
      .vel_x_peculiar = state.gas_cells.velocity_x_peculiar[row],
      .vel_y_peculiar = state.gas_cells.velocity_y_peculiar[row],
      .vel_z_peculiar = state.gas_cells.velocity_z_peculiar[row],
      .pressure_comoving = state.gas_cells.pressure_code[row],
      .metallicity_mass_fraction = std::clamp(
          state.gas_cells.metal_mass_code[row] /
              std::max(state.cells.mass_code[row], 1.0e-30),
          0.0,
          1.0)};
  const auto conserved = hydro::HydroCoreSolver::conservedFromPrimitive(primitive, adiabatic_index);
  if (!finiteConserved(conserved) || conserved.mass_density_comoving <= 0.0) {
    throw std::runtime_error("AMR temporal history source state is non-finite or non-positive");
  }
  return conserved;
}

void assignStart(core::AmrTemporalBoundaryHistoryCellRecord& record, const hydro::HydroConservedState& state) {
  record.start_mass_density_comoving = state.mass_density_comoving;
  record.start_momentum_density_x_comoving = state.momentum_density_x_comoving;
  record.start_momentum_density_y_comoving = state.momentum_density_y_comoving;
  record.start_momentum_density_z_comoving = state.momentum_density_z_comoving;
  record.start_total_energy_density_comoving = state.total_energy_density_comoving;
  record.start_metal_mass_density_comoving = state.metal_mass_density_comoving;
}

void assignEnd(core::AmrTemporalBoundaryHistoryCellRecord& record, const hydro::HydroConservedState& state) {
  record.end_mass_density_comoving = state.mass_density_comoving;
  record.end_momentum_density_x_comoving = state.momentum_density_x_comoving;
  record.end_momentum_density_y_comoving = state.momentum_density_y_comoving;
  record.end_momentum_density_z_comoving = state.momentum_density_z_comoving;
  record.end_total_energy_density_comoving = state.total_energy_density_comoving;
  record.end_metal_mass_density_comoving = state.metal_mass_density_comoving;
}

}  // namespace

std::uint64_t amrPatchGeometryFingerprint(const PatchDescriptor& patch) {
  std::uint64_t hash = 1469598103934665603ULL;
  const auto mix = [&hash](std::uint64_t value) {
    hash ^= value;
    hash *= 1099511628211ULL;
  };
  mix(patch.patch_id);
  mix(patch.parent_patch_id);
  mix(static_cast<std::uint64_t>(patch.level));
  mix(patch.morton_key);
  for (const double value : patch.origin_comov) {
    mix(std::bit_cast<std::uint64_t>(value));
  }
  for (const double value : patch.extent_comov) {
    mix(std::bit_cast<std::uint64_t>(value));
  }
  for (const std::uint16_t value : patch.cell_dims) {
    mix(static_cast<std::uint64_t>(value));
  }
  return hash == 0U ? 1U : hash;
}

void captureAmrTemporalBoundaryHistoryStart(
    core::SimulationState& state,
    std::span<const PatchDescriptor> patches,
    double interval_start_code,
    double adiabatic_index) {
  if (!std::isfinite(interval_start_code)) {
    throw std::invalid_argument("captureAmrTemporalBoundaryHistoryStart requires finite interval_start_code");
  }
  if (!state.amr_temporal_boundary_history.empty()) {
    throw std::runtime_error("cannot begin a new AMR temporal boundary interval while an earlier interval is active");
  }
  std::vector<core::AmrTemporalBoundaryHistoryRecord> records;
  for (const PatchDescriptor& patch : patches) {
    const bool is_coarse_source = std::any_of(
        patches.begin(), patches.end(),
        [&patch](const PatchDescriptor& candidate) {
          return candidate.level > patch.level;
        });
    if (!is_coarse_source) {
      continue;
    }
    const AmrHydroPatchGeometry geometry = buildAmrHydroPatchGeometry(state, patch);
    core::AmrTemporalBoundaryHistoryRecord record;
    record.patch_id = patch.patch_id;
    record.patch_level = patch.level;
    record.patch_geometry_fingerprint = amrPatchGeometryFingerprint(patch);
    record.gas_cell_identity_generation = state.gasCellIdentityGeneration();
    record.interval_start_code = interval_start_code;
    record.interval_end_code = interval_start_code;
    record.end_state_valid = false;
    record.cells.reserve(geometry.real_cells.size());
    for (const AmrHydroCellDescriptor& descriptor : geometry.real_cells) {
      core::AmrTemporalBoundaryHistoryCellRecord cell;
      cell.gas_cell_id = descriptor.gas_cell_id;
      cell.patch_local_cell = descriptor.patch_local_cell;
      assignStart(cell, conservedForStateRow(state, descriptor.local_cell_row, adiabatic_index));
      record.cells.push_back(cell);
    }
    records.push_back(std::move(record));
  }
  state.amr_temporal_boundary_history.assign(std::move(records));
}

void captureAmrTemporalBoundaryHistoryEnd(
    core::SimulationState& state,
    std::span<const PatchDescriptor> patches,
    double interval_end_code,
    double adiabatic_index) {
  if (!std::isfinite(interval_end_code)) {
    throw std::invalid_argument("captureAmrTemporalBoundaryHistoryEnd requires finite interval_end_code");
  }
  for (core::AmrTemporalBoundaryHistoryRecord& record : state.amr_temporal_boundary_history.mutableRecords()) {
    const auto patch_it = std::find_if(
        patches.begin(), patches.end(),
        [&record](const PatchDescriptor& patch) { return patch.patch_id == record.patch_id; });
    if (patch_it == patches.end() || record.patch_geometry_fingerprint != amrPatchGeometryFingerprint(*patch_it)) {
      throw std::runtime_error("cannot complete AMR temporal history after patch geometry changed");
    }
    if (record.gas_cell_identity_generation != state.gasCellIdentityGeneration()) {
      throw std::runtime_error("cannot complete AMR temporal history after gas-cell identity generation changed");
    }
    if (interval_end_code <= record.interval_start_code) {
      throw std::invalid_argument("AMR temporal history interval must have positive duration");
    }
    for (core::AmrTemporalBoundaryHistoryCellRecord& cell : record.cells) {
      const auto row = state.rowForGasCellId(cell.gas_cell_id);
      if (!row.has_value()) {
        throw std::runtime_error("cannot complete AMR temporal history because coarse gas_cell_id is missing");
      }
      if (state.cells.patch_index[*row] >= state.patches.size() ||
          state.patches.patch_id[state.cells.patch_index[*row]] != record.patch_id) {
        throw std::runtime_error("cannot complete AMR temporal history because coarse gas cell changed patch ownership");
      }
      assignEnd(cell, conservedForStateRow(state, *row, adiabatic_index));
    }
    record.interval_end_code = interval_end_code;
    record.end_state_valid = true;
  }
}

void retireAmrTemporalBoundaryHistory(core::SimulationState& state) {
  state.amr_temporal_boundary_history.clear();
}

AmrHydroGhostFillDiagnostics fillAmrHydroGhostCells(
    std::span<AmrHydroGhostFillPatch> patches,
    std::span<const AmrHydroSparseGhostSource> remote_sources,
    double adiabatic_index) {
  AmrHydroGhostFillDiagnostics diagnostics;
  for (const AmrHydroGhostFillPatch& patch : patches) {
    validatePatchView(patch);
  }
  for (const AmrHydroSparseGhostSource& source : remote_sources) {
    if (source.patch == nullptr) {
      throw std::invalid_argument(
          "fillAmrHydroGhostCells: sparse remote source is missing its patch descriptor");
    }
    if (!std::isfinite(source.source_current_state_time_code)) {
      throw std::invalid_argument(
          "fillAmrHydroGhostCells: sparse remote source time must be finite");
    }
    const std::size_t patch_cell_count =
        static_cast<std::size_t>(source.patch->cell_dims[0]) *
        static_cast<std::size_t>(source.patch->cell_dims[1]) *
        static_cast<std::size_t>(source.patch->cell_dims[2]);
    std::size_t previous_offset = 0U;
    bool have_previous = false;
    for (const AmrHydroSparseRemoteCell& cell : source.cells) {
      const std::size_t offset = cell.patch_local_cell;
      if (offset >= patch_cell_count || cell.gas_cell_id == 0U) {
        throw std::invalid_argument(
            "fillAmrHydroGhostCells: sparse remote cell is outside its patch or has zero identity");
      }
      if (have_previous && offset <= previous_offset) {
        throw std::invalid_argument(
            "fillAmrHydroGhostCells: sparse remote cells must be strictly sorted by patch-local offset");
      }
      previous_offset = offset;
      have_previous = true;
    }
  }

  for (AmrHydroGhostFillPatch& patch : patches) {
    if (!patch.requires_ghost_fill) {
      continue;
    }
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
      fillAmrGhost(
          patch, ghost, patches, remote_sources, adiabatic_index, diagnostics);
    }
  }
  return diagnostics;
}

AmrHydroGhostFillDiagnostics fillAmrHydroGhostCells(
    std::span<AmrHydroGhostFillPatch> patches,
    double adiabatic_index) {
  return fillAmrHydroGhostCells(patches, {}, adiabatic_index);
}

}  // namespace cosmosim::amr
