#include "cosmosim/amr/amr_hydro_orchestrator.hpp"

#include "cosmosim/amr/amr_patch_indexing.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>

namespace cosmosim::amr {
namespace {

constexpr double k_geometry_tol = 1.0e-10;

[[nodiscard]] std::size_t product(std::array<std::uint16_t, 3> dims) {
  return static_cast<std::size_t>(dims[0]) * static_cast<std::size_t>(dims[1]) *
      static_cast<std::size_t>(dims[2]);
}

[[nodiscard]] std::size_t linearIndex(
    std::array<std::uint16_t, 3> dims,
    std::size_t i,
    std::size_t j,
    std::size_t k) {
  return i + static_cast<std::size_t>(dims[0]) * (j + static_cast<std::size_t>(dims[1]) * k);
}

[[nodiscard]] std::vector<double> sortedUnique(std::vector<double> values) {
  std::sort(values.begin(), values.end());
  values.erase(
      std::unique(
          values.begin(),
          values.end(),
          [](double lhs, double rhs) {
            const double scale = std::max({1.0, std::abs(lhs), std::abs(rhs)});
            return std::abs(lhs - rhs) <= k_geometry_tol * scale;
          }),
      values.end());
  return values;
}

[[nodiscard]] double spacingForUnique(const std::vector<double>& values, double fallback) {
  if (values.size() <= 1U) {
    return std::max(1.0e-12, fallback);
  }
  return std::max(1.0e-12, (values.back() - values.front()) / static_cast<double>(values.size() - 1U));
}

[[nodiscard]] std::array<std::uint16_t, 3> nearCubicFactors(std::size_t count) {
  if (count == 0U || count > std::numeric_limits<std::uint16_t>::max()) {
    throw std::invalid_argument("AMR production patch descriptor requires positive representable cell count");
  }
  std::array<std::uint16_t, 3> best{1, 1, static_cast<std::uint16_t>(count)};
  std::size_t best_surface = std::numeric_limits<std::size_t>::max();
  for (std::size_t nx = 1; nx <= count; ++nx) {
    if (count % nx != 0U) {
      continue;
    }
    const std::size_t yz = count / nx;
    for (std::size_t ny = 1; ny <= yz; ++ny) {
      if (yz % ny != 0U) {
        continue;
      }
      const std::size_t nz = yz / ny;
      const std::size_t surface = nx * ny + nx * nz + ny * nz;
      if (surface < best_surface && nx <= std::numeric_limits<std::uint16_t>::max() &&
          ny <= std::numeric_limits<std::uint16_t>::max() && nz <= std::numeric_limits<std::uint16_t>::max()) {
        best = {static_cast<std::uint16_t>(nx), static_cast<std::uint16_t>(ny), static_cast<std::uint16_t>(nz)};
        best_surface = surface;
      }
    }
  }
  return best;
}

[[nodiscard]] std::array<double, 3> cellWidths(const PatchDescriptor& patch) {
  return {
      patch.extent_comov[0] / static_cast<double>(patch.cell_dims[0]),
      patch.extent_comov[1] / static_cast<double>(patch.cell_dims[1]),
      patch.extent_comov[2] / static_cast<double>(patch.cell_dims[2]),
  };
}

[[nodiscard]] std::array<double, 3> cellCenter(const PatchDescriptor& patch, std::size_t cell) {
  const auto widths = cellWidths(patch);
  const std::size_t nx = patch.cell_dims[0];
  const std::size_t ny = patch.cell_dims[1];
  const std::size_t i = cell % nx;
  const std::size_t j = (cell / nx) % ny;
  const std::size_t k = cell / (nx * ny);
  return {
      patch.origin_comov[0] + (static_cast<double>(i) + 0.5) * widths[0],
      patch.origin_comov[1] + (static_cast<double>(j) + 0.5) * widths[1],
      patch.origin_comov[2] + (static_cast<double>(k) + 0.5) * widths[2],
  };
}

[[nodiscard]] bool containsPoint(const PatchDescriptor& patch, const std::array<double, 3>& point) {
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
  const auto widths = cellWidths(patch);
  std::array<std::size_t, 3> ijk{};
  for (std::size_t axis = 0; axis < 3; ++axis) {
    double coordinate = (point[axis] - patch.origin_comov[axis]) / widths[axis];
    coordinate = std::clamp(coordinate, 0.0, static_cast<double>(patch.cell_dims[axis]) - k_geometry_tol);
    ijk[axis] = static_cast<std::size_t>(coordinate);
    if (ijk[axis] >= patch.cell_dims[axis]) {
      ijk[axis] = patch.cell_dims[axis] - 1U;
    }
  }
  return linearIndex(patch.cell_dims, ijk[0], ijk[1], ijk[2]);
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

[[nodiscard]] hydro::HydroFaceSide oppositeSide(hydro::HydroFaceSide side) {
  return side == hydro::HydroFaceSide::kLower ? hydro::HydroFaceSide::kUpper : hydro::HydroFaceSide::kLower;
}

[[nodiscard]] std::array<double, 3> ghostProbePoint(
    const AmrHydroPatchGeometry& patch_geometry,
    const hydro::HydroGhostCell& ghost) {
  std::array<double, 3> point = cellCenter(patch_geometry.patch, ghost.owner_real_cell);
  const auto widths = cellWidths(patch_geometry.patch);
  point[axisIndex(ghost.axis)] += sideSign(ghost.side) * widths[axisIndex(ghost.axis)];
  return point;
}

[[nodiscard]] std::optional<PatchDescriptor> sourcePatchForGhost(
    const AmrHydroPatchGeometry& target,
    const hydro::HydroGhostCell& ghost,
    std::span<const PatchDescriptor> all_patches) {
  const std::array<double, 3> probe = ghostProbePoint(target, ghost);
  for (const PatchDescriptor& candidate : all_patches) {
    if (candidate.patch_id == target.patch.patch_id) {
      continue;
    }
    if (containsPoint(candidate, probe)) {
      return candidate;
    }
  }
  return std::nullopt;
}

[[nodiscard]] bool intervalsOverlap(double a0, double a1, double b0, double b1) {
  return std::max(a0, b0) < std::min(a1, b1) + k_geometry_tol;
}

[[nodiscard]] bool touchesOnSide(
    const PatchDescriptor& target,
    const PatchDescriptor& candidate,
    hydro::HydroFaceAxis axis,
    hydro::HydroFaceSide side) {
  const int normal_axis = axisIndex(axis);
  const double target_face = side == hydro::HydroFaceSide::kLower
      ? target.origin_comov[normal_axis]
      : target.origin_comov[normal_axis] + target.extent_comov[normal_axis];
  const double candidate_face = side == hydro::HydroFaceSide::kLower
      ? candidate.origin_comov[normal_axis] + candidate.extent_comov[normal_axis]
      : candidate.origin_comov[normal_axis];
  const double scale = std::max({1.0, std::abs(target_face), std::abs(candidate_face)});
  if (std::abs(target_face - candidate_face) > k_geometry_tol * scale) {
    return false;
  }
  for (int other_axis = 0; other_axis < 3; ++other_axis) {
    if (other_axis == normal_axis) {
      continue;
    }
    const double a0 = target.origin_comov[other_axis];
    const double a1 = a0 + target.extent_comov[other_axis];
    const double b0 = candidate.origin_comov[other_axis];
    const double b1 = b0 + candidate.extent_comov[other_axis];
    if (!intervalsOverlap(a0, a1, b0, b1)) {
      return false;
    }
  }
  return true;
}

[[nodiscard]] AmrHydroBoundaryClass boundaryClassForSide(
    const PatchDescriptor& target,
    std::span<const PatchDescriptor> patches,
    hydro::HydroFaceAxis axis,
    hydro::HydroFaceSide side) {
  for (const PatchDescriptor& candidate : patches) {
    if (candidate.patch_id == target.patch_id) {
      continue;
    }
    if (!touchesOnSide(target, candidate, axis, side)) {
      continue;
    }
    return candidate.level == target.level
        ? AmrHydroBoundaryClass::kSameLevel
        : AmrHydroBoundaryClass::kCoarseFine;
  }
  return AmrHydroBoundaryClass::kPhysical;
}

[[nodiscard]] std::uint64_t makeRegisterKey(
    std::uint64_t coarse_patch_id,
    std::size_t coarse_cell_index,
    hydro::HydroFaceAxis axis,
    hydro::HydroFaceSide orientation) {
  std::uint64_t key = 1469598103934665603ULL;
  const auto mix = [&key](std::uint64_t value) {
    key ^= value;
    key *= 1099511628211ULL;
  };
  mix(coarse_patch_id);
  mix(static_cast<std::uint64_t>(coarse_cell_index) + 1U);
  mix(static_cast<std::uint64_t>(axisIndex(axis)) + 11U);
  mix(orientation == hydro::HydroFaceSide::kLower ? 23U : 29U);
  return key == 0U ? 1U : key;
}

[[nodiscard]] double safeDensity(const core::SimulationState& state, std::uint32_t row, double floor) {
  return std::max(state.gas_cells.density_code[row], floor);
}

[[nodiscard]] hydro::HydroPrimitiveState primitiveForRow(
    const core::SimulationState& state,
    std::uint32_t row,
    double density_floor,
    double pressure_floor,
    double adiabatic_index) {
  const double rho = safeDensity(state, row, density_floor);
  double pressure = state.gas_cells.pressure_code[row];
  if (pressure <= pressure_floor) {
    const double internal_energy = std::max(state.gas_cells.internal_energy_code[row], pressure_floor);
    pressure = std::max((adiabatic_index - 1.0) * rho * internal_energy, pressure_floor);
  }
  return hydro::HydroPrimitiveState{
      .rho_comoving = rho,
      .vel_x_peculiar = state.gas_cells.velocity_x_peculiar[row],
      .vel_y_peculiar = state.gas_cells.velocity_y_peculiar[row],
      .vel_z_peculiar = state.gas_cells.velocity_z_peculiar[row],
      .pressure_comoving = pressure};
}

[[nodiscard]] ConservedState volumeIntegratedForRow(
    const core::SimulationState& state,
    std::uint32_t row,
    double volume,
    const ProductionAmrHydroOptions& options) {
  const hydro::HydroConservedState density = hydro::HydroCoreSolver::conservedFromPrimitive(
      primitiveForRow(state, row, options.density_floor, options.pressure_floor, options.adiabatic_index),
      options.adiabatic_index);
  return ConservedState{
      .mass_code = density.mass_density_comoving * volume,
      .momentum_x_code = density.momentum_density_x_comoving * volume,
      .momentum_y_code = density.momentum_density_y_comoving * volume,
      .momentum_z_code = density.momentum_density_z_comoving * volume,
      .total_energy_code = density.total_energy_density_comoving * volume};
}


void refreshTraceableGhostSourceIds(std::vector<AmrHydroPatchGeometry>& geometries) {
  std::unordered_map<std::uint64_t, std::size_t> patch_slot_by_id;
  patch_slot_by_id.reserve(geometries.size());
  for (std::size_t slot = 0; slot < geometries.size(); ++slot) {
    patch_slot_by_id.emplace(geometries[slot].patch.patch_id, slot);
  }
  for (AmrHydroPatchGeometry& target : geometries) {
    for (AmrHydroGhostDescriptor& ghost_descriptor : target.ghosts) {
      if (ghost_descriptor.source_patch_id == 0U) {
        continue;
      }
      const auto source_slot_it = patch_slot_by_id.find(ghost_descriptor.source_patch_id);
      if (source_slot_it == patch_slot_by_id.end()) {
        continue;
      }
      const AmrHydroPatchGeometry& source = geometries[source_slot_it->second];
      const hydro::HydroGhostCell& ghost = target.geometry.ghost_cells.at(ghost_descriptor.ghost_slot);
      const std::array<double, 3> probe = ghostProbePoint(target, ghost);
      if (!containsPoint(source.patch, probe)) {
        continue;
      }
      const std::size_t source_cell = cellContainingPoint(source.patch, probe);
      if (source_cell < source.gas_cell_ids.size()) {
        ghost_descriptor.source_gas_cell_id = source.gas_cell_ids[source_cell];
        for (std::size_t face_index = 0; face_index < target.geometry.faces.size(); ++face_index) {
          const hydro::HydroFace& face = target.geometry.faces[face_index];
          if (face.ghost_cell_slot == ghost_descriptor.ghost_slot &&
              face_index < target.geometry.flux_register_faces.size() &&
              target.geometry.flux_register_faces[face_index].role == hydro::HydroFluxRegisterFaceRole::kFine &&
              target.geometry.flux_register_faces[face_index].coarse_gas_cell_id == 0U) {
            target.geometry.flux_register_faces[face_index].coarse_gas_cell_id = ghost_descriptor.source_gas_cell_id;
          }
        }
      }
    }
  }
}

void writeVolumeIntegratedToRow(
    core::SimulationState& state,
    std::uint32_t row,
    const ConservedState& volume_state,
    double volume,
    const ProductionAmrHydroOptions& options) {
  if (volume <= 0.0) {
    throw std::invalid_argument("writeVolumeIntegratedToRow requires positive volume");
  }
  hydro::HydroConservedState density_state{
      .mass_density_comoving = volume_state.mass_code / volume,
      .momentum_density_x_comoving = volume_state.momentum_x_code / volume,
      .momentum_density_y_comoving = volume_state.momentum_y_code / volume,
      .momentum_density_z_comoving = volume_state.momentum_z_code / volume,
      .total_energy_density_comoving = volume_state.total_energy_code / volume};
  hydro::HydroPrimitiveState primitive =
      hydro::HydroCoreSolver::primitiveFromConserved(density_state, options.adiabatic_index);
  primitive.rho_comoving = std::max(primitive.rho_comoving, options.density_floor);
  primitive.pressure_comoving = std::max(primitive.pressure_comoving, options.pressure_floor);
  state.cells.mass_code[row] = primitive.rho_comoving * volume;
  state.gas_cells.density_code[row] = primitive.rho_comoving;
  state.gas_cells.pressure_code[row] = primitive.pressure_comoving;
  state.gas_cells.internal_energy_code[row] =
      primitive.pressure_comoving / ((options.adiabatic_index - 1.0) * primitive.rho_comoving);
  state.gas_cells.velocity_x_peculiar[row] = primitive.vel_x_peculiar;
  state.gas_cells.velocity_y_peculiar[row] = primitive.vel_y_peculiar;
  state.gas_cells.velocity_z_peculiar[row] = primitive.vel_z_peculiar;
  state.gas_cells.sound_speed_code[row] =
      std::sqrt(std::max(0.0, options.adiabatic_index * primitive.pressure_comoving / primitive.rho_comoving));
  state.gas_cells.temperature_code[row] = state.gas_cells.internal_energy_code[row];
}

[[nodiscard]] std::vector<std::uint32_t> patchLocalRowsForPatch(
    const core::SimulationState& state,
    const PatchDescriptor& patch,
    std::string_view context) {
  return rowsForPatchLocalCellOrder(state, patch, context);
}

[[nodiscard]] std::array<double, 5> totalsForRows(
    const core::SimulationState& state,
    std::span<const std::uint32_t> rows,
    double volume,
    const ProductionAmrHydroOptions& options) {
  std::array<double, 5> totals{};
  for (const std::uint32_t row : rows) {
    const ConservedState conserved = volumeIntegratedForRow(state, row, volume, options);
    totals[0] += conserved.mass_code;
    totals[1] += conserved.momentum_x_code;
    totals[2] += conserved.momentum_y_code;
    totals[3] += conserved.momentum_z_code;
    totals[4] += conserved.total_energy_code;
  }
  return totals;
}

void copyCellRow(core::SimulationState& state, std::uint32_t dst, const core::SimulationState& old, std::uint32_t src) {
  state.cells.center_x_comoving[dst] = old.cells.center_x_comoving[src];
  state.cells.center_y_comoving[dst] = old.cells.center_y_comoving[src];
  state.cells.center_z_comoving[dst] = old.cells.center_z_comoving[src];
  state.cells.mass_code[dst] = old.cells.mass_code[src];
  state.cells.time_bin[dst] = old.cells.time_bin[src];
  state.gas_cells.gas_cell_id[dst] = old.gas_cells.gas_cell_id[src];
  state.gas_cells.parent_particle_id[dst] = old.gas_cells.parent_particle_id[src];
  state.gas_cells.velocity_x_peculiar[dst] = old.gas_cells.velocity_x_peculiar[src];
  state.gas_cells.velocity_y_peculiar[dst] = old.gas_cells.velocity_y_peculiar[src];
  state.gas_cells.velocity_z_peculiar[dst] = old.gas_cells.velocity_z_peculiar[src];
  state.gas_cells.density_code[dst] = old.gas_cells.density_code[src];
  state.gas_cells.pressure_code[dst] = old.gas_cells.pressure_code[src];
  state.gas_cells.internal_energy_code[dst] = old.gas_cells.internal_energy_code[src];
  state.gas_cells.temperature_code[dst] = old.gas_cells.temperature_code[src];
  state.gas_cells.sound_speed_code[dst] = old.gas_cells.sound_speed_code[src];
}

void rebuildIdentityFromSidecars(core::SimulationState& state) {
  std::vector<core::GasCellIdentityRecord> records;
  records.reserve(state.cells.size());
  for (std::uint32_t row = 0; row < state.cells.size(); ++row) {
    if (state.gas_cells.gas_cell_id[row] == 0U) {
      throw std::runtime_error("production AMR hydro generated a zero gas_cell_id");
    }
    const std::uint32_t patch_index = state.cells.patch_index[row];
    if (patch_index >= state.patches.size()) {
      throw std::runtime_error("production AMR hydro generated an invalid patch index");
    }
    const std::uint64_t parent_id = state.gas_cells.parent_particle_id[row];
    records.push_back(core::GasCellIdentityRecord{
        .gas_cell_id = state.gas_cells.gas_cell_id[row],
        .parent_particle_id = parent_id == 0U ? std::optional<std::uint64_t>{} : std::optional<std::uint64_t>{parent_id},
        .owning_patch_id = state.patches.patch_id[patch_index],
        .local_cell_row = row});
  }
  state.replaceGasCellIdentityRecords(std::move(records));
}

[[nodiscard]] bool idRangeOverflows(std::uint64_t first, std::size_t count) {
  return count == 0U || first == 0U || first > std::numeric_limits<std::uint64_t>::max() - static_cast<std::uint64_t>(count - 1U);
}

void validatePatchIdRangeAvailable(
    const core::SimulationState& state,
    std::uint64_t first_patch_id,
    std::size_t count,
    std::uint64_t allowed_existing_patch_id = 0U) {
  if (idRangeOverflows(first_patch_id, count)) {
    throw std::invalid_argument("production AMR regrid patch-id range is zero or overflows");
  }
  std::unordered_set<std::uint64_t> requested;
  requested.reserve(count);
  for (std::size_t i = 0; i < count; ++i) {
    requested.insert(first_patch_id + static_cast<std::uint64_t>(i));
  }
  for (const std::uint64_t patch_id : state.patches.patch_id) {
    if (patch_id != allowed_existing_patch_id && requested.contains(patch_id)) {
      throw std::invalid_argument("production AMR regrid patch-id range collides with existing state");
    }
  }
}

void validateGasCellIdRangeAvailable(
    const core::SimulationState& state,
    std::uint64_t first_gas_cell_id,
    std::size_t count) {
  if (idRangeOverflows(first_gas_cell_id, count)) {
    throw std::invalid_argument("production AMR regrid gas-cell-id range is zero or overflows");
  }
  std::unordered_set<std::uint64_t> requested;
  requested.reserve(count);
  for (std::size_t i = 0; i < count; ++i) {
    requested.insert(first_gas_cell_id + static_cast<std::uint64_t>(i));
  }
  for (const std::uint64_t gas_cell_id : state.gas_cells.gas_cell_id) {
    if (requested.contains(gas_cell_id)) {
      throw std::invalid_argument("production AMR regrid gas-cell-id range collides with existing state");
    }
  }
  for (const core::GasCellIdentityRecord& record : state.gas_cell_identity.records()) {
    if (requested.contains(record.gas_cell_id)) {
      throw std::invalid_argument("production AMR regrid gas-cell-id range collides with identity map");
    }
  }
}

[[nodiscard]] std::uint64_t nextFreePatchId(const core::SimulationState& state, std::size_t count) {
  std::uint64_t candidate = 1U;
  for (const std::uint64_t patch_id : state.patches.patch_id) {
    candidate = std::max(candidate, patch_id + 1U);
  }
  validatePatchIdRangeAvailable(state, candidate, count);
  return candidate;
}

[[nodiscard]] std::uint64_t nextFreeGasCellId(const core::SimulationState& state, std::size_t count) {
  std::uint64_t candidate = 1U;
  for (const std::uint64_t gas_cell_id : state.gas_cells.gas_cell_id) {
    candidate = std::max(candidate, gas_cell_id + 1U);
  }
  validateGasCellIdRangeAvailable(state, candidate, count);
  return candidate;
}

[[nodiscard]] std::optional<std::uint32_t> patchIndexById(
    const core::SimulationState& state,
    std::uint64_t patch_id) {
  for (std::uint32_t index = 0; index < state.patches.size(); ++index) {
    if (state.patches.patch_id[index] == patch_id) {
      return index;
    }
  }
  return std::nullopt;
}

}  // namespace

bool patchStateRowHasExplicitGeometry(
    const core::SimulationState& state,
    std::size_t patch_index) {
  if (patch_index >= state.patches.size()) {
    return false;
  }
  return state.patches.extent_x_comoving[patch_index] > 0.0 &&
      state.patches.extent_y_comoving[patch_index] > 0.0 &&
      state.patches.extent_z_comoving[patch_index] > 0.0 &&
      state.patches.cell_dim_x[patch_index] != 0U &&
      state.patches.cell_dim_y[patch_index] != 0U &&
      state.patches.cell_dim_z[patch_index] != 0U;
}

PatchDescriptor patchDescriptorFromStateRow(
    const core::SimulationState& state,
    std::size_t patch_index) {
  if (patch_index >= state.patches.size()) {
    throw std::out_of_range("patchDescriptorFromStateRow: patch_index out of range");
  }
  if (!patchStateRowHasExplicitGeometry(state, patch_index)) {
    throw std::runtime_error("patchDescriptorFromStateRow: patch row lacks explicit restart-authoritative geometry");
  }
  PatchDescriptor descriptor;
  descriptor.patch_id = state.patches.patch_id[patch_index];
  descriptor.parent_patch_id = state.patches.parent_patch_id[patch_index];
  descriptor.level = static_cast<std::uint8_t>(std::max<std::int32_t>(state.patches.level[patch_index], 0));
  descriptor.morton_key = state.patches.morton_key[patch_index];
  descriptor.origin_comov = {
      state.patches.origin_x_comoving[patch_index],
      state.patches.origin_y_comoving[patch_index],
      state.patches.origin_z_comoving[patch_index]};
  descriptor.extent_comov = {
      state.patches.extent_x_comoving[patch_index],
      state.patches.extent_y_comoving[patch_index],
      state.patches.extent_z_comoving[patch_index]};
  descriptor.cell_dims = {
      state.patches.cell_dim_x[patch_index],
      state.patches.cell_dim_y[patch_index],
      state.patches.cell_dim_z[patch_index]};
  if (descriptor.patch_id == 0U || product(descriptor.cell_dims) != state.patches.cell_count[patch_index]) {
    throw std::runtime_error("patchDescriptorFromStateRow: invalid explicit patch descriptor lanes");
  }
  return descriptor;
}

void writePatchDescriptorToStateRow(
    core::SimulationState& state,
    std::size_t patch_index,
    const PatchDescriptor& patch) {
  if (patch_index >= state.patches.size()) {
    throw std::out_of_range("writePatchDescriptorToStateRow: patch_index out of range");
  }
  if (patch.patch_id == 0U || patch.extent_comov[0] <= 0.0 || patch.extent_comov[1] <= 0.0 ||
      patch.extent_comov[2] <= 0.0 || patch.cell_dims[0] == 0U || patch.cell_dims[1] == 0U ||
      patch.cell_dims[2] == 0U) {
    throw std::invalid_argument("writePatchDescriptorToStateRow: invalid explicit patch geometry");
  }
  state.patches.patch_id[patch_index] = patch.patch_id;
  state.patches.parent_patch_id[patch_index] = patch.parent_patch_id;
  state.patches.level[patch_index] = static_cast<std::int32_t>(patch.level);
  state.patches.morton_key[patch_index] = patch.morton_key;
  state.patches.origin_x_comoving[patch_index] = patch.origin_comov[0];
  state.patches.origin_y_comoving[patch_index] = patch.origin_comov[1];
  state.patches.origin_z_comoving[patch_index] = patch.origin_comov[2];
  state.patches.extent_x_comoving[patch_index] = patch.extent_comov[0];
  state.patches.extent_y_comoving[patch_index] = patch.extent_comov[1];
  state.patches.extent_z_comoving[patch_index] = patch.extent_comov[2];
  state.patches.cell_dim_x[patch_index] = patch.cell_dims[0];
  state.patches.cell_dim_y[patch_index] = patch.cell_dims[1];
  state.patches.cell_dim_z[patch_index] = patch.cell_dims[2];
}

bool hasProductionAmrHydroCoverage(const core::SimulationState& state) {
  if (state.cells.size() == 0U || state.patches.size() == 0U) {
    return false;
  }
  if (!state.gas_cell_identity.isConsistent() ||
      !state.gas_cell_identity.coversDenseLocalRows(state.cells.size()) ||
      !state.gasCellIdentityMapMatchesSidecarLanes()) {
    return false;
  }
  std::vector<std::uint8_t> covered(state.cells.size(), 0U);
  for (std::size_t patch = 0; patch < state.patches.size(); ++patch) {
    const std::uint32_t first = state.patches.first_cell[patch];
    const std::uint32_t count = state.patches.cell_count[patch];
    if (count == 0U) {
      continue;
    }
    if (!patchStateRowHasExplicitGeometry(state, patch)) {
      return false;
    }
    if (static_cast<std::size_t>(state.patches.cell_dim_x[patch]) *
            static_cast<std::size_t>(state.patches.cell_dim_y[patch]) *
            static_cast<std::size_t>(state.patches.cell_dim_z[patch]) != count) {
      return false;
    }
    if (first + count > state.cells.size()) {
      return false;
    }
    for (std::uint32_t offset = 0; offset < count; ++offset) {
      const std::uint32_t row = first + offset;
      if (covered[row] != 0U || state.cells.patch_index[row] != patch) {
        return false;
      }
      const auto* record = state.gas_cell_identity.findByLocalRow(row);
      if (record == nullptr || record->owning_patch_id != state.patches.patch_id[patch]) {
        return false;
      }
      covered[row] = 1U;
    }
  }
  return std::all_of(covered.begin(), covered.end(), [](std::uint8_t flag) { return flag != 0U; });
}

std::vector<PatchDescriptor> buildProductionAmrPatchDescriptors(const core::SimulationState& state) {
  if (!hasProductionAmrHydroCoverage(state)) {
    throw std::runtime_error(
        "buildProductionAmrPatchDescriptors requires complete SimulationState AMR patch coverage with explicit patch geometry");
  }
  std::vector<PatchDescriptor> descriptors;
  descriptors.reserve(state.patches.size());
  for (std::size_t patch_index = 0; patch_index < state.patches.size(); ++patch_index) {
    if (state.patches.cell_count[patch_index] == 0U) {
      continue;
    }
    descriptors.push_back(patchDescriptorFromStateRow(state, patch_index));
  }
  return descriptors;
}

void populateAmrHydroFluxRegisterFaces(
    AmrHydroPatchGeometry& patch_geometry,
    std::span<const PatchDescriptor> all_patches) {
  patch_geometry.geometry.flux_register_faces.assign(
      patch_geometry.geometry.faces.size(),
      hydro::HydroFluxRegisterFace{});
  for (std::size_t face_index = 0; face_index < patch_geometry.geometry.faces.size(); ++face_index) {
    const hydro::HydroFace& face = patch_geometry.geometry.faces[face_index];
    if (face.ghost_cell_slot == hydro::k_invalid_ghost_cell_slot) {
      continue;
    }
    const auto ghost_it = std::find_if(
        patch_geometry.ghosts.begin(),
        patch_geometry.ghosts.end(),
        [&face](const AmrHydroGhostDescriptor& ghost) { return ghost.ghost_slot == face.ghost_cell_slot; });
    if (ghost_it == patch_geometry.ghosts.end() ||
        ghost_it->boundary_class == AmrHydroBoundaryClass::kPhysical) {
      continue;
    }
    const hydro::HydroGhostCell& ghost = patch_geometry.geometry.ghost_cells.at(face.ghost_cell_slot);
    const std::optional<PatchDescriptor> source_patch = sourcePatchForGhost(patch_geometry, ghost, all_patches);
    if (!source_patch.has_value()) {
      continue;
    }
    AmrHydroGhostDescriptor& ghost_descriptor = *ghost_it;
    ghost_descriptor.source_patch_id = source_patch->patch_id;
    const std::array<double, 3> probe = ghostProbePoint(patch_geometry, ghost);
    [[maybe_unused]] const std::size_t source_cell = cellContainingPoint(*source_patch, probe);
    // Source gas-cell IDs are patched after all local AMR hydro geometries exist, because this
    // metadata belongs to the source patch's gas_cell_id vector, not the target patch.
    if (ghost_descriptor.boundary_class != AmrHydroBoundaryClass::kCoarseFine) {
      continue;
    }

    const bool target_is_coarse = patch_geometry.patch.level < source_patch->level;
    const PatchDescriptor& coarse_patch = target_is_coarse ? patch_geometry.patch : *source_patch;
    const std::size_t coarse_cell = target_is_coarse ? face.owner_cell : cellContainingPoint(coarse_patch, probe);
    const hydro::HydroFaceSide coarse_orientation = target_is_coarse ? ghost.side : oppositeSide(ghost.side);
    patch_geometry.geometry.flux_register_faces[face_index] = hydro::HydroFluxRegisterFace{
        .role = target_is_coarse ? hydro::HydroFluxRegisterFaceRole::kCoarse : hydro::HydroFluxRegisterFaceRole::kFine,
        .register_key = makeRegisterKey(coarse_patch.patch_id, coarse_cell, ghost.axis, coarse_orientation),
        .coarse_patch_id = coarse_patch.patch_id,
        .coarse_gas_cell_id = target_is_coarse ? patch_geometry.gas_cell_ids.at(coarse_cell) : ghost_descriptor.source_gas_cell_id,
        .coarse_cell_index = coarse_cell,
        .level = static_cast<int>(coarse_patch.level),
        .axis = ghost.axis,
        .orientation = coarse_orientation,
        .coarse_orientation_sign = target_is_coarse ? 1.0 : -1.0};
  }
}

void scatterAmrHydroConservedState(
    core::SimulationState& state,
    const AmrHydroPatchGeometry& patch_geometry,
    const hydro::HydroConservedStateSoa& conserved,
    double adiabatic_index) {
  state.requireGasCellIdentityMapFresh(
      patch_geometry.source_gas_cell_identity_generation,
      "scatterAmrHydroConservedState");
  if (conserved.size() < patch_geometry.geometry.cellCount()) {
    throw std::invalid_argument("scatterAmrHydroConservedState: conserved storage lacks real AMR cells");
  }
  std::unordered_map<std::uint64_t, std::uint32_t> parent_row_by_id;
  parent_row_by_id.reserve(state.particles.size());
  for (std::uint32_t row = 0; row < state.particles.size(); ++row) {
    parent_row_by_id.emplace(state.particle_sidecar.particle_id[row], row);
  }
  std::unordered_map<std::uint64_t, std::size_t> parent_use_count;
  for (const auto& record : state.gas_cell_identity.records()) {
    if (record.parent_particle_id.has_value()) {
      parent_use_count[*record.parent_particle_id] += 1U;
    }
  }

  for (std::size_t patch_cell = 0; patch_cell < patch_geometry.real_cells.size(); ++patch_cell) {
    const AmrHydroCellDescriptor& cell = patch_geometry.real_cells[patch_cell];
    const auto row = state.rowForGasCellId(cell.gas_cell_id);
    if (!row.has_value() || *row != cell.local_cell_row) {
      throw std::runtime_error("scatterAmrHydroConservedState: gas_cell_id resolved to a different row");
    }
    const hydro::HydroPrimitiveState primitive =
        hydro::HydroCoreSolver::primitiveFromConserved(conserved.loadCell(patch_cell), adiabatic_index);
    const double rho = std::max(primitive.rho_comoving, 1.0e-14);
    const double pressure = std::max(primitive.pressure_comoving, 1.0e-14);
    state.gas_cells.density_code[*row] = rho;
    state.gas_cells.pressure_code[*row] = pressure;
    state.gas_cells.internal_energy_code[*row] = pressure / ((adiabatic_index - 1.0) * rho);
    state.gas_cells.velocity_x_peculiar[*row] = primitive.vel_x_peculiar;
    state.gas_cells.velocity_y_peculiar[*row] = primitive.vel_y_peculiar;
    state.gas_cells.velocity_z_peculiar[*row] = primitive.vel_z_peculiar;
    state.gas_cells.sound_speed_code[*row] = std::sqrt(std::max(0.0, adiabatic_index * pressure / rho));
    state.gas_cells.temperature_code[*row] = state.gas_cells.internal_energy_code[*row];
    state.cells.mass_code[*row] = rho * patch_geometry.geometry.cell_volume_comoving;

    const auto parent_id = state.parentParticleIdForGasCellId(cell.gas_cell_id);
    if (parent_id.has_value() && parent_use_count[*parent_id] == 1U) {
      const auto parent_it = parent_row_by_id.find(*parent_id);
      if (parent_it != parent_row_by_id.end()) {
        const std::uint32_t parent_row = parent_it->second;
        state.particles.mass_code[parent_row] = state.cells.mass_code[*row];
        state.particles.velocity_x_peculiar[parent_row] = primitive.vel_x_peculiar;
        state.particles.velocity_y_peculiar[parent_row] = primitive.vel_y_peculiar;
        state.particles.velocity_z_peculiar[parent_row] = primitive.vel_z_peculiar;
      }
    }
  }
}

RefluxDiagnostics applyFluxRegistersToSimulationState(
    core::SimulationState& state,
    std::span<const FluxRegisterEntry> entries,
    std::span<const PatchDescriptor> all_patches,
    double adiabatic_index) {
  RefluxDiagnostics diagnostics;
  ProductionAmrHydroOptions options;
  options.adiabatic_index = adiabatic_index;
  for (const FluxRegisterEntry& entry : entries) {
    if (!entry.isComplete()) {
      ++diagnostics.skipped_incomplete_register_count;
      continue;
    }
    const double area_scale = std::max({1.0, std::abs(entry.coarse_area_comov), std::abs(entry.fine_area_comov)});
    if (entry.coarse_area_comov <= 0.0 || entry.fine_area_comov <= 0.0 ||
        std::abs(entry.coarse_area_comov - entry.fine_area_comov) > 1.0e-10 * area_scale) {
      ++diagnostics.skipped_area_mismatch_count;
      continue;
    }
    if (entry.coarse_gas_cell_id == 0U) {
      ++diagnostics.skipped_missing_target_count;
      continue;
    }
    const auto patch_it = std::find_if(
        all_patches.begin(),
        all_patches.end(),
        [&entry](const PatchDescriptor& patch) { return patch.patch_id == entry.coarse_patch_id; });
    if (patch_it == all_patches.end()) {
      ++diagnostics.skipped_missing_target_count;
      continue;
    }
    const auto row_opt = state.gas_cell_identity.rowForGasCellId(entry.coarse_gas_cell_id);
    if (!row_opt.has_value()) {
      ++diagnostics.skipped_missing_target_count;
      continue;
    }
    const std::uint32_t row = *row_opt;
    const auto* record = state.gas_cell_identity.findByGasCellId(entry.coarse_gas_cell_id);
    if (record == nullptr || record->owning_patch_id != entry.coarse_patch_id) {
      ++diagnostics.skipped_missing_target_count;
      continue;
    }
    const auto patch_index = patchIndexById(state, entry.coarse_patch_id);
    if (!patch_index.has_value() || row >= state.cells.patch_index.size() ||
        state.cells.patch_index[row] != *patch_index ||
        state.gas_cells.gas_cell_id[row] != entry.coarse_gas_cell_id ||
        entry.coarse_cell_index >= product(patch_it->cell_dims) ||
        patchLocalCellForRow(state, *patch_it, row, "applyFluxRegistersToSimulationState").linear_index != entry.coarse_cell_index) {
      throw std::runtime_error("applyFluxRegistersToSimulationState rejected a stale AMR reflux target mapping");
    }
    const double volume = patch_it->extent_comov[0] * patch_it->extent_comov[1] * patch_it->extent_comov[2] /
        static_cast<double>(std::max<std::size_t>(product(patch_it->cell_dims), 1U));
    ++diagnostics.complete_register_count;
    const hydro::HydroPrimitiveState old_primitive = primitiveForRow(
        state, row, 1.0e-14, 1.0e-14, adiabatic_index);
    hydro::HydroConservedState conserved =
        hydro::HydroCoreSolver::conservedFromPrimitive(old_primitive, adiabatic_index);
    ConservedState delta_flux = entry.fine_face_flux_code - entry.coarse_face_flux_code;
    delta_flux *= (entry.coarse_area_comov * entry.dt_code / volume);
    conserved.mass_density_comoving -= delta_flux.mass_code;
    conserved.momentum_density_x_comoving -= delta_flux.momentum_x_code;
    conserved.momentum_density_y_comoving -= delta_flux.momentum_y_code;
    conserved.momentum_density_z_comoving -= delta_flux.momentum_z_code;
    conserved.total_energy_density_comoving -= delta_flux.total_energy_code;
    const double old_internal_density = old_primitive.pressure_comoving / (adiabatic_index - 1.0);
    const hydro::HydroPrimitiveState primitive =
        hydro::HydroCoreSolver::primitiveFromConserved(conserved, adiabatic_index);
    state.gas_cells.density_code[row] = primitive.rho_comoving;
    state.gas_cells.pressure_code[row] = primitive.pressure_comoving;
    state.gas_cells.internal_energy_code[row] =
        primitive.pressure_comoving / ((adiabatic_index - 1.0) * std::max(primitive.rho_comoving, 1.0e-14));
    state.gas_cells.velocity_x_peculiar[row] = primitive.vel_x_peculiar;
    state.gas_cells.velocity_y_peculiar[row] = primitive.vel_y_peculiar;
    state.gas_cells.velocity_z_peculiar[row] = primitive.vel_z_peculiar;
    state.gas_cells.sound_speed_code[row] = std::sqrt(std::max(0.0, adiabatic_index * primitive.pressure_comoving /
        std::max(primitive.rho_comoving, 1.0e-14)));
    state.gas_cells.temperature_code[row] = state.gas_cells.internal_energy_code[row];
    state.cells.mass_code[row] = primitive.rho_comoving * volume;
    const double new_internal_density = primitive.pressure_comoving / (adiabatic_index - 1.0);

    diagnostics.corrected_cells += 1U;
    diagnostics.corrected_mass_code += std::abs(delta_flux.mass_code * volume);
    diagnostics.corrected_momentum_x_code += std::abs(delta_flux.momentum_x_code * volume);
    diagnostics.corrected_momentum_y_code += std::abs(delta_flux.momentum_y_code * volume);
    diagnostics.corrected_momentum_z_code += std::abs(delta_flux.momentum_z_code * volume);
    diagnostics.corrected_total_energy_code += std::abs(delta_flux.total_energy_code * volume);
    diagnostics.corrected_energy_code += std::abs(delta_flux.total_energy_code * volume);
    diagnostics.corrected_internal_energy_code += std::abs((new_internal_density - old_internal_density) * volume);
  }

  return diagnostics;
}

namespace {

[[nodiscard]] std::uint8_t axisToStorage(hydro::HydroFaceAxis axis) {
  return static_cast<std::uint8_t>(axisIndex(axis));
}

[[nodiscard]] std::uint8_t sideToStorage(hydro::HydroFaceSide side) {
  return side == hydro::HydroFaceSide::kLower ? 0U : 1U;
}

[[nodiscard]] hydro::HydroFaceAxis axisFromStorage(std::uint8_t axis) {
  switch (axis) {
    case 0U:
      return hydro::HydroFaceAxis::kX;
    case 1U:
      return hydro::HydroFaceAxis::kY;
    case 2U:
      return hydro::HydroFaceAxis::kZ;
    default:
      throw std::runtime_error("pending AMR flux register has invalid face axis metadata");
  }
}

[[nodiscard]] hydro::HydroFaceSide sideFromStorage(std::uint8_t side) {
  if (side > 1U) {
    throw std::runtime_error("pending AMR flux register has invalid face-side metadata");
  }
  return side == 0U ? hydro::HydroFaceSide::kLower : hydro::HydroFaceSide::kUpper;
}

void addFluxIntegral(
    double& mass,
    double& momentum_x,
    double& momentum_y,
    double& momentum_z,
    double& total_energy,
    const ConservedState& flux,
    double area_comov,
    double dt_code) {
  const double scale = area_comov * dt_code;
  mass += flux.mass_code * scale;
  momentum_x += flux.momentum_x_code * scale;
  momentum_y += flux.momentum_y_code * scale;
  momentum_z += flux.momentum_z_code * scale;
  total_energy += flux.total_energy_code * scale;
}

[[nodiscard]] ConservedState averageFluxFromIntegral(
    double mass,
    double momentum_x,
    double momentum_y,
    double momentum_z,
    double total_energy,
    double area_comov,
    double dt_code) {
  const double denom = area_comov * dt_code;
  if (denom <= 0.0) {
    return ConservedState{};
  }
  return ConservedState{
      .mass_code = mass / denom,
      .momentum_x_code = momentum_x / denom,
      .momentum_y_code = momentum_y / denom,
      .momentum_z_code = momentum_z / denom,
      .total_energy_code = total_energy / denom};
}

void validatePendingCompatible(
    const core::PendingFluxRegisterRecord& pending,
    const FluxRegisterEntry& entry) {
  if (pending.coarse_patch_id != entry.coarse_patch_id ||
      pending.coarse_gas_cell_id != entry.coarse_gas_cell_id ||
      pending.coarse_cell_index != entry.coarse_cell_index ||
      pending.level != entry.level ||
      pending.axis != axisToStorage(entry.axis) ||
      pending.orientation != sideToStorage(entry.orientation)) {
    throw std::runtime_error("pending AMR flux register received incompatible stable identity metadata for an existing key");
  }
}

[[nodiscard]] core::PendingFluxRegisterRecord makePendingRecordSeed(
    const core::SimulationState& state,
    const FluxRegisterEntry& entry,
    const ProductionAmrHydroOptions& options) {
  const double coarse_dt = options.reflux_coarse_dt_code > 0.0 ? options.reflux_coarse_dt_code : entry.dt_code;
  const double interval_start = options.reflux_interval_start_code;
  const double interval_end = options.reflux_interval_end_code > interval_start
      ? options.reflux_interval_end_code
      : interval_start + coarse_dt;
  const double expected_area = std::max({entry.face_area_comov, entry.coarse_area_comov, entry.fine_area_comov, 0.0});
  return core::PendingFluxRegisterRecord{
      .register_key = entry.register_key,
      .coarse_patch_id = entry.coarse_patch_id,
      .coarse_gas_cell_id = entry.coarse_gas_cell_id,
      .coarse_cell_index = entry.coarse_cell_index,
      .level = entry.level,
      .axis = axisToStorage(entry.axis),
      .orientation = sideToStorage(entry.orientation),
      .expected_area_comov = expected_area,
      .interval_start_code = interval_start,
      .interval_end_code = interval_end,
      .coarse_dt_code = coarse_dt,
      .expected_fine_substeps = std::max<std::uint32_t>(1U, options.expected_fine_substeps),
      .gas_cell_identity_generation = state.gasCellIdentityGeneration(),
      .patch_geometry_generation = state.cellIndexGeneration()};
}

[[nodiscard]] FluxRegisterEntry fluxEntryFromPendingRecord(
    const core::PendingFluxRegisterRecord& pending) {
  const double area = pending.expected_area_comov;
  const double dt = pending.coarse_dt_code > 0.0
      ? pending.coarse_dt_code
      : std::max(0.0, pending.interval_end_code - pending.interval_start_code);
  return FluxRegisterEntry{
      .register_key = pending.register_key,
      .coarse_patch_id = pending.coarse_patch_id,
      .coarse_gas_cell_id = pending.coarse_gas_cell_id,
      .coarse_cell_index = pending.coarse_cell_index,
      .level = pending.level,
      .axis = axisFromStorage(pending.axis),
      .orientation = sideFromStorage(pending.orientation),
      .coarse_face_flux_code = averageFluxFromIntegral(
          pending.coarse_mass_flux_integral_code,
          pending.coarse_momentum_x_flux_integral_code,
          pending.coarse_momentum_y_flux_integral_code,
          pending.coarse_momentum_z_flux_integral_code,
          pending.coarse_total_energy_flux_integral_code,
          area,
          dt),
      .fine_face_flux_code = averageFluxFromIntegral(
          pending.fine_mass_flux_integral_code,
          pending.fine_momentum_x_flux_integral_code,
          pending.fine_momentum_y_flux_integral_code,
          pending.fine_momentum_z_flux_integral_code,
          pending.fine_total_energy_flux_integral_code,
          area,
          dt),
      .face_area_comov = area,
      .coarse_area_comov = pending.coarse_area_accumulated_comov,
      .fine_area_comov = pending.fine_area_accumulated_comov,
      .dt_code = dt,
      .coarse_face_count = pending.coarse_face_count,
      .fine_face_count = pending.fine_face_count};
}

void accumulateDiagnostics(ProductionAmrHydroDiagnostics& lhs, const ProductionAmrHydroDiagnostics& rhs) {
  lhs.patch_count = std::max(lhs.patch_count, rhs.patch_count);
  lhs.advanced_patch_count += rhs.advanced_patch_count;
  lhs.active_cell_count += rhs.active_cell_count;
  lhs.active_face_count += rhs.active_face_count;
  lhs.flux_register_entry_count += rhs.flux_register_entry_count;
  lhs.pending_register_created_count += rhs.pending_register_created_count;
  lhs.pending_register_completed_count += rhs.pending_register_completed_count;
  lhs.pending_register_deferred_count += rhs.pending_register_deferred_count;
  lhs.pending_register_applied_count += rhs.pending_register_applied_count;
  lhs.pending_register_rejected_count += rhs.pending_register_rejected_count;
  lhs.ghost_fill.physical_ghosts_filled += rhs.ghost_fill.physical_ghosts_filled;
  lhs.ghost_fill.same_level_ghosts_filled += rhs.ghost_fill.same_level_ghosts_filled;
  lhs.ghost_fill.coarse_to_fine_ghosts_filled += rhs.ghost_fill.coarse_to_fine_ghosts_filled;
  lhs.ghost_fill.fine_to_coarse_ghosts_filled += rhs.ghost_fill.fine_to_coarse_ghosts_filled;
  lhs.ghost_fill.temporal_coarse_to_fine_ghosts_filled += rhs.ghost_fill.temporal_coarse_to_fine_ghosts_filled;
  lhs.ghost_fill.temporal_endpoint_ghosts_filled += rhs.ghost_fill.temporal_endpoint_ghosts_filled;
  lhs.ghost_fill.skipped_remote_ghosts += rhs.ghost_fill.skipped_remote_ghosts;
  lhs.ghost_fill.stale_epoch_rejections += rhs.ghost_fill.stale_epoch_rejections;
  lhs.ghost_fill.temporal_same_level_mismatch_rejections += rhs.ghost_fill.temporal_same_level_mismatch_rejections;
  lhs.ghost_fill.temporal_history_missing_rejections += rhs.ghost_fill.temporal_history_missing_rejections;
  lhs.ghost_fill.temporal_history_invalid_rejections += rhs.ghost_fill.temporal_history_invalid_rejections;
  lhs.ghost_fill.temporal_time_out_of_range_rejections += rhs.ghost_fill.temporal_time_out_of_range_rejections;
  lhs.ghost_fill.temporal_fine_to_coarse_misalignment_rejections += rhs.ghost_fill.temporal_fine_to_coarse_misalignment_rejections;
  lhs.ghost_fill.temporal_geometry_mismatch_rejections += rhs.ghost_fill.temporal_geometry_mismatch_rejections;
  lhs.ghost_fill.temporal_identity_mismatch_rejections += rhs.ghost_fill.temporal_identity_mismatch_rejections;
  lhs.ghost_fill.missing_source_records += rhs.ghost_fill.missing_source_records;
  lhs.ghost_fill.unresolved_ghosts += rhs.ghost_fill.unresolved_ghosts;
  lhs.reflux.complete_register_count += rhs.reflux.complete_register_count;
  lhs.reflux.skipped_incomplete_register_count += rhs.reflux.skipped_incomplete_register_count;
  lhs.reflux.skipped_area_mismatch_count += rhs.reflux.skipped_area_mismatch_count;
  lhs.reflux.skipped_missing_target_count += rhs.reflux.skipped_missing_target_count;
  lhs.reflux.corrected_cells += rhs.reflux.corrected_cells;
  lhs.reflux.corrected_mass_code += rhs.reflux.corrected_mass_code;
  lhs.reflux.corrected_momentum_x_code += rhs.reflux.corrected_momentum_x_code;
  lhs.reflux.corrected_momentum_y_code += rhs.reflux.corrected_momentum_y_code;
  lhs.reflux.corrected_momentum_z_code += rhs.reflux.corrected_momentum_z_code;
  lhs.reflux.corrected_total_energy_code += rhs.reflux.corrected_total_energy_code;
  lhs.reflux.corrected_energy_code += rhs.reflux.corrected_energy_code;
  lhs.reflux.corrected_internal_energy_code += rhs.reflux.corrected_internal_energy_code;
}

[[nodiscard]] std::vector<std::uint32_t> activeRowsForLevel(
    const core::SimulationState& state,
    int level,
    std::span<const std::uint32_t> requested_rows) {
  std::unordered_set<std::uint32_t> requested_lookup(requested_rows.begin(), requested_rows.end());
  const bool all_requested = requested_rows.empty();
  std::vector<std::uint32_t> rows;
  for (std::uint32_t row = 0; row < state.cells.size(); ++row) {
    if (!all_requested && !requested_lookup.contains(row)) {
      continue;
    }
    if (row >= state.cells.patch_index.size()) {
      continue;
    }
    const std::uint32_t patch_index = state.cells.patch_index[row];
    if (patch_index < state.patches.size() && state.patches.level[patch_index] == level) {
      rows.push_back(row);
    }
  }
  return rows;
}

}  // namespace

std::size_t mergeFluxRegistersIntoPendingStore(
    core::SimulationState& state,
    std::span<const FluxRegisterEntry> entries,
    const ProductionAmrHydroOptions& options) {
  std::size_t created_count = 0;
  for (const FluxRegisterEntry& entry : entries) {
    if (entry.register_key == 0U) {
      continue;
    }
    auto* pending = state.pending_flux_registers.findByRegisterKey(entry.register_key);
    if (pending == nullptr) {
      core::PendingFluxRegisterRecord seed = makePendingRecordSeed(state, entry, options);
      pending = &state.pending_flux_registers.upsertByRegisterKey(seed);
      ++created_count;
    } else {
      validatePendingCompatible(*pending, entry);
      pending->expected_fine_substeps = std::max(pending->expected_fine_substeps, options.expected_fine_substeps);
      pending->expected_area_comov = std::max({pending->expected_area_comov, entry.face_area_comov, entry.coarse_area_comov, entry.fine_area_comov});
    }

    if (entry.coarse_face_count > 0U && entry.coarse_area_comov > 0.0) {
      pending->coarse_face_count += entry.coarse_face_count;
      pending->coarse_area_accumulated_comov = std::max(pending->coarse_area_accumulated_comov, entry.coarse_area_comov);
      addFluxIntegral(
          pending->coarse_mass_flux_integral_code,
          pending->coarse_momentum_x_flux_integral_code,
          pending->coarse_momentum_y_flux_integral_code,
          pending->coarse_momentum_z_flux_integral_code,
          pending->coarse_total_energy_flux_integral_code,
          entry.coarse_face_flux_code,
          entry.coarse_area_comov,
          entry.dt_code);
    }
    if (entry.fine_face_count > 0U && entry.fine_area_comov > 0.0) {
      pending->fine_face_count += entry.fine_face_count;
      pending->fine_area_accumulated_comov = std::max(pending->fine_area_accumulated_comov, entry.fine_area_comov);
      addFluxIntegral(
          pending->fine_mass_flux_integral_code,
          pending->fine_momentum_x_flux_integral_code,
          pending->fine_momentum_y_flux_integral_code,
          pending->fine_momentum_z_flux_integral_code,
          pending->fine_total_energy_flux_integral_code,
          entry.fine_face_flux_code,
          entry.fine_area_comov,
          entry.dt_code);
      const std::uint32_t expected = std::max<std::uint32_t>(1U, pending->expected_fine_substeps);
      const std::uint32_t substep = expected == 1U ? 0U : std::min(options.fine_substep_index, expected - 1U);
      if (substep < 64U) {
        const std::uint64_t bit = 1ULL << substep;
        if ((pending->fine_substep_coverage_mask & bit) == 0U) {
          pending->fine_substep_coverage_mask |= bit;
          pending->completed_fine_substeps += 1U;
        }
      } else {
        pending->completed_fine_substeps += 1U;
      }
    }
  }
  return created_count;
}

RefluxDiagnostics applyCompletePendingFluxRegistersToSimulationState(
    core::SimulationState& state,
    std::span<const PatchDescriptor> all_patches,
    double adiabatic_index) {
  RefluxDiagnostics diagnostics;
  std::vector<std::uint64_t> applied_keys;
  for (const core::PendingFluxRegisterRecord& pending : state.pending_flux_registers.records()) {
    if (!pending.isComplete()) {
      ++diagnostics.skipped_incomplete_register_count;
      continue;
    }
    if (pending.gas_cell_identity_generation != state.gasCellIdentityGeneration()) {
      ++diagnostics.skipped_missing_target_count;
      continue;
    }
    const FluxRegisterEntry entry = fluxEntryFromPendingRecord(pending);
    RefluxDiagnostics one = applyFluxRegistersToSimulationState(state, std::span<const FluxRegisterEntry>(&entry, 1), all_patches, adiabatic_index);
    diagnostics.complete_register_count += one.complete_register_count;
    diagnostics.skipped_incomplete_register_count += one.skipped_incomplete_register_count;
    diagnostics.skipped_area_mismatch_count += one.skipped_area_mismatch_count;
    diagnostics.skipped_missing_target_count += one.skipped_missing_target_count;
    diagnostics.corrected_cells += one.corrected_cells;
    diagnostics.corrected_mass_code += one.corrected_mass_code;
    diagnostics.corrected_momentum_x_code += one.corrected_momentum_x_code;
    diagnostics.corrected_momentum_y_code += one.corrected_momentum_y_code;
    diagnostics.corrected_momentum_z_code += one.corrected_momentum_z_code;
    diagnostics.corrected_total_energy_code += one.corrected_total_energy_code;
    diagnostics.corrected_energy_code += one.corrected_energy_code;
    diagnostics.corrected_internal_energy_code += one.corrected_internal_energy_code;
    if (one.complete_register_count == 1U) {
      applied_keys.push_back(pending.register_key);
    }
  }
  state.pending_flux_registers.eraseCompletedByKey(std::span<const std::uint64_t>(applied_keys.data(), applied_keys.size()));
  return diagnostics;
}

namespace {

[[nodiscard]] ProductionAmrHydroDiagnostics advanceProductionAmrHydroSynchronizedImpl(
    core::SimulationState& state,
    std::span<const std::uint32_t> active_cell_rows,
    const hydro::HydroUpdateContext& update,
    const hydro::HydroSourceContext& global_source_context,
    const hydro::HydroCoreSolver& solver,
    const hydro::HydroRiemannSolver& riemann_solver,
    std::span<const hydro::HydroSourceTerm* const> source_terms,
    const ProductionAmrHydroOptions& options,
    const DistributedAmrHydroExchange* distributed_exchange) {
  if (!hasProductionAmrHydroCoverage(state)) {
    throw std::runtime_error("advanceProductionAmrHydro requires complete AMR patch coverage in SimulationState");
  }
  if (distributed_exchange != nullptr && options.persist_incomplete_flux_registers) {
    throw std::runtime_error(
        "distributed AMR pending flux-register accumulation is gated off; synchronized remote reflux must complete in-step");
  }
  ProductionAmrHydroDiagnostics diagnostics;
  const std::vector<PatchDescriptor> descriptors = buildProductionAmrPatchDescriptors(state);
  std::vector<PatchDescriptor> all_descriptors = descriptors;
  if (distributed_exchange != nullptr) {
    all_descriptors.reserve(descriptors.size() + distributed_exchange->remote_patches.size());
    for (const DistributedAmrRemotePatch& remote : distributed_exchange->remote_patches) {
      all_descriptors.push_back(remote.patch);
    }
  }
  diagnostics.patch_count = all_descriptors.size();

  std::unordered_set<std::uint32_t> active_lookup(active_cell_rows.begin(), active_cell_rows.end());
  const bool all_cells_active = active_cell_rows.empty();

  std::vector<AmrHydroPatchGeometry> geometries;
  std::vector<hydro::HydroConservedStateSoa> conserved_states;
  std::vector<AmrHydroPatchGeometry> remote_geometries;
  std::vector<hydro::HydroConservedStateSoa> remote_conserved_states;
  geometries.reserve(descriptors.size());
  conserved_states.reserve(descriptors.size());
  for (const PatchDescriptor& patch : descriptors) {
    AmrHydroGeometryOptions geometry_options;
    geometry_options.physical_boundary_kind = options.physical_boundary_kind;
    geometry_options.boundary_classes = {
        boundaryClassForSide(patch, all_descriptors, hydro::HydroFaceAxis::kX, hydro::HydroFaceSide::kLower),
        boundaryClassForSide(patch, all_descriptors, hydro::HydroFaceAxis::kX, hydro::HydroFaceSide::kUpper),
        boundaryClassForSide(patch, all_descriptors, hydro::HydroFaceAxis::kY, hydro::HydroFaceSide::kLower),
        boundaryClassForSide(patch, all_descriptors, hydro::HydroFaceAxis::kY, hydro::HydroFaceSide::kUpper),
        boundaryClassForSide(patch, all_descriptors, hydro::HydroFaceAxis::kZ, hydro::HydroFaceSide::kLower),
        boundaryClassForSide(patch, all_descriptors, hydro::HydroFaceAxis::kZ, hydro::HydroFaceSide::kUpper),
    };
    geometries.push_back(buildAmrHydroPatchGeometry(state, patch, geometry_options));
    populateAmrHydroFluxRegisterFaces(geometries.back(), all_descriptors);
    conserved_states.push_back(loadAmrHydroConservedState(state, geometries.back(), options.adiabatic_index));
  }
  if (distributed_exchange != nullptr) {
    remote_geometries.reserve(distributed_exchange->remote_patches.size());
    remote_conserved_states.reserve(distributed_exchange->remote_patches.size());
    for (const DistributedAmrRemotePatch& remote : distributed_exchange->remote_patches) {
      if (remote.owner_rank == distributed_exchange->local_rank) {
        throw std::runtime_error("distributed AMR remote patch source cannot be owned by the local rank");
      }
      if (remote.conserved_cells.size() != remote.gas_cell_ids.size()) {
        throw std::runtime_error("distributed AMR remote patch conserved payload does not match gas-cell identity payload");
      }
      AmrHydroGeometryOptions geometry_options;
      geometry_options.physical_boundary_kind = options.physical_boundary_kind;
      geometry_options.boundary_classes = {
          boundaryClassForSide(remote.patch, all_descriptors, hydro::HydroFaceAxis::kX, hydro::HydroFaceSide::kLower),
          boundaryClassForSide(remote.patch, all_descriptors, hydro::HydroFaceAxis::kX, hydro::HydroFaceSide::kUpper),
          boundaryClassForSide(remote.patch, all_descriptors, hydro::HydroFaceAxis::kY, hydro::HydroFaceSide::kLower),
          boundaryClassForSide(remote.patch, all_descriptors, hydro::HydroFaceAxis::kY, hydro::HydroFaceSide::kUpper),
          boundaryClassForSide(remote.patch, all_descriptors, hydro::HydroFaceAxis::kZ, hydro::HydroFaceSide::kLower),
          boundaryClassForSide(remote.patch, all_descriptors, hydro::HydroFaceAxis::kZ, hydro::HydroFaceSide::kUpper),
      };
      remote_geometries.push_back(buildRemoteAmrHydroPatchGeometry(
          remote.patch,
          remote.gas_cell_ids,
          remote.expected_ghost_hydro_epoch,
          geometry_options));
      populateAmrHydroFluxRegisterFaces(remote_geometries.back(), all_descriptors);
      hydro::HydroConservedStateSoa remote_soa(remote_geometries.back().geometry.totalCellStorageCount());
      for (std::size_t cell = 0; cell < remote.conserved_cells.size(); ++cell) {
        remote_soa.storeCell(cell, remote.conserved_cells[cell]);
      }
      remote_conserved_states.push_back(std::move(remote_soa));
    }
  }

  refreshTraceableGhostSourceIds(geometries);
  refreshTraceableGhostSourceIds(remote_geometries);
  std::vector<AmrHydroPatchGeometry> trace_geometries;
  if (!remote_geometries.empty()) {
    trace_geometries.reserve(geometries.size() + remote_geometries.size());
    for (const auto& geometry : geometries) {
      trace_geometries.push_back(geometry);
    }
    for (const auto& geometry : remote_geometries) {
      trace_geometries.push_back(geometry);
    }
    refreshTraceableGhostSourceIds(trace_geometries);
    for (std::size_t i = 0; i < geometries.size(); ++i) {
      geometries[i].ghosts = trace_geometries[i].ghosts;
      geometries[i].geometry.flux_register_faces = trace_geometries[i].geometry.flux_register_faces;
    }
    for (std::size_t i = 0; i < remote_geometries.size(); ++i) {
      remote_geometries[i].ghosts = trace_geometries[geometries.size() + i].ghosts;
      remote_geometries[i].geometry.flux_register_faces =
          trace_geometries[geometries.size() + i].geometry.flux_register_faces;
    }
  }

  std::vector<AmrHydroGhostFillPatch> ghost_views;
  ghost_views.reserve(geometries.size() + remote_geometries.size());
  for (std::size_t i = 0; i < geometries.size(); ++i) {
    const bool patch_requires_ghost_fill = all_cells_active || std::any_of(
        geometries[i].real_cells.begin(), geometries[i].real_cells.end(),
        [&active_lookup](const AmrHydroCellDescriptor& cell) {
          return active_lookup.contains(cell.local_cell_row);
        });
    ghost_views.push_back(AmrHydroGhostFillPatch{
        .geometry = &geometries[i],
        .conserved = &conserved_states[i],
        .residency = AmrGhostSourceResidency::kLocal,
        .ghost_hydro_epoch = state.gasCellIdentityGeneration(),
        .expected_ghost_hydro_epoch = state.gasCellIdentityGeneration(),
        .target_state_time_code = options.state_time_code,
        .ghost_fill_time_code = options.ghost_fill_time_code,
        .source_current_state_time_code = options.state_time_code,
        .temporal_boundary_history = options.enable_temporal_coarse_to_fine
            ? &state.amr_temporal_boundary_history
            : nullptr,
        .enable_temporal_coarse_to_fine = options.enable_temporal_coarse_to_fine,
        .requires_ghost_fill = patch_requires_ghost_fill});
  }
  if (distributed_exchange != nullptr) {
    for (std::size_t i = 0; i < remote_geometries.size(); ++i) {
      const DistributedAmrRemotePatch& remote = distributed_exchange->remote_patches[i];
      ghost_views.push_back(AmrHydroGhostFillPatch{
          .geometry = &remote_geometries[i],
          .conserved = &remote_conserved_states[i],
          .residency = AmrGhostSourceResidency::kRemoteReadOnly,
          .ghost_hydro_epoch = remote.ghost_hydro_epoch,
          .expected_ghost_hydro_epoch = remote.expected_ghost_hydro_epoch,
          .target_state_time_code = options.state_time_code,
          .ghost_fill_time_code = options.ghost_fill_time_code,
          .source_current_state_time_code = options.state_time_code,
          .temporal_boundary_history = nullptr,
          .enable_temporal_coarse_to_fine = false,
          .requires_ghost_fill = false});
    }
  }
  diagnostics.ghost_fill = fillAmrHydroGhostCells(ghost_views, options.adiabatic_index);

  FluxRegisterAccumulator flux_registers;
  for (std::size_t patch_index = 0; patch_index < geometries.size(); ++patch_index) {
    const AmrHydroPatchGeometry& patch_geometry = geometries[patch_index];
    std::vector<std::size_t> active_patch_cells;
    active_patch_cells.reserve(patch_geometry.real_cells.size());
    for (const AmrHydroCellDescriptor& cell : patch_geometry.real_cells) {
      if (all_cells_active || active_lookup.contains(cell.local_cell_row)) {
        active_patch_cells.push_back(cell.patch_local_cell);
      }
    }
    if (active_patch_cells.empty()) {
      continue;
    }
    std::unordered_set<std::size_t> patch_active_lookup(active_patch_cells.begin(), active_patch_cells.end());
    std::vector<std::size_t> active_faces;
    for (std::size_t face_index = 0; face_index < patch_geometry.geometry.faces.size(); ++face_index) {
      const hydro::HydroFace& face = patch_geometry.geometry.faces[face_index];
      if (patch_active_lookup.contains(face.owner_cell) ||
          (face.neighbor_cell < patch_geometry.geometry.cellCount() && patch_active_lookup.contains(face.neighbor_cell))) {
        active_faces.push_back(face_index);
      }
    }

    std::vector<double> gravity_x(conserved_states[patch_index].size(), 0.0);
    std::vector<double> gravity_y(conserved_states[patch_index].size(), 0.0);
    std::vector<double> gravity_z(conserved_states[patch_index].size(), 0.0);
    std::vector<double> hydrogen(conserved_states[patch_index].size(), 0.0);
    std::vector<double> metallicity(conserved_states[patch_index].size(), 0.0);
    std::vector<double> temperature(conserved_states[patch_index].size(), 0.0);
    for (std::size_t local_cell = 0; local_cell < patch_geometry.real_cells.size(); ++local_cell) {
      const std::uint32_t row = patch_geometry.real_cells[local_cell].local_cell_row;
      if (row < global_source_context.gravity_accel_x_peculiar.size()) {
        gravity_x[local_cell] = global_source_context.gravity_accel_x_peculiar[row];
      }
      if (row < global_source_context.gravity_accel_y_peculiar.size()) {
        gravity_y[local_cell] = global_source_context.gravity_accel_y_peculiar[row];
      }
      if (row < global_source_context.gravity_accel_z_peculiar.size()) {
        gravity_z[local_cell] = global_source_context.gravity_accel_z_peculiar[row];
      }
      if (row < global_source_context.hydrogen_number_density_cgs.size()) {
        hydrogen[local_cell] = global_source_context.hydrogen_number_density_cgs[row];
      }
      if (row < global_source_context.metallicity_mass_fraction.size()) {
        metallicity[local_cell] = global_source_context.metallicity_mass_fraction[row];
      }
      if (row < global_source_context.temperature_k.size()) {
        temperature[local_cell] = global_source_context.temperature_k[row];
      }
    }
    const hydro::HydroSourceContext patch_source_context{
        .update = update,
        .gravity_accel_x_peculiar = gravity_x,
        .gravity_accel_y_peculiar = gravity_y,
        .gravity_accel_z_peculiar = gravity_z,
        .hydrogen_number_density_cgs = hydrogen,
        .metallicity_mass_fraction = metallicity,
        .temperature_k = temperature,
        .redshift = global_source_context.redshift};

    const auto widths = cellWidths(patch_geometry.patch);
    hydro::MusclHancockReconstruction reconstruction(hydro::HydroReconstructionPolicy{
        .limiter = hydro::HydroSlopeLimiter::kMonotonizedCentral,
        .dt_over_dx_code = update.dt_code / std::max(widths[0], 1.0e-12),
        .dt_over_cell_width_code = {
            update.dt_code / std::max(widths[0], 1.0e-12),
            update.dt_code / std::max(widths[1], 1.0e-12),
            update.dt_code / std::max(widths[2], 1.0e-12)},
        .rho_floor = options.density_floor,
        .pressure_floor = options.pressure_floor,
        .enable_muscl_hancock_predictor = true,
        .adiabatic_index = options.adiabatic_index});
    hydro::HydroScratchBuffers scratch;
    hydro::HydroPrimitiveCacheSoa cache(conserved_states[patch_index].size());
    hydro::HydroProfileEvent profile;
    solver.advancePatchActiveSetWithScratch(
        conserved_states[patch_index],
        patch_geometry.geometry,
        hydro::HydroActiveSetView{.active_cells = active_patch_cells, .active_faces = active_faces},
        update,
        reconstruction,
        riemann_solver,
        source_terms,
        patch_source_context,
        scratch,
        &cache,
        &profile,
        &flux_registers);
    diagnostics.advanced_patch_count += 1U;
    diagnostics.active_cell_count += active_patch_cells.size();
    diagnostics.active_face_count += active_faces.size();
  }

  for (std::size_t patch_index = 0; patch_index < geometries.size(); ++patch_index) {
    scatterAmrHydroConservedState(state, geometries[patch_index], conserved_states[patch_index], options.adiabatic_index);
  }
  const std::vector<FluxRegisterEntry> entries = flux_registers.entries();
  diagnostics.flux_register_entry_count = entries.size();
  std::vector<FluxRegisterEntry> local_entries;
  std::vector<FluxRegisterEntry> remote_entries;
  if (distributed_exchange != nullptr) {
    local_entries.reserve(entries.size());
    remote_entries.reserve(entries.size());
    for (const FluxRegisterEntry& entry : entries) {
      if (entry.coarse_gas_cell_id != 0U && state.rowForGasCellId(entry.coarse_gas_cell_id).has_value()) {
        local_entries.push_back(entry);
      } else {
        remote_entries.push_back(entry);
      }
    }
    if (distributed_exchange->outbound_remote_flux_registers != nullptr) {
      distributed_exchange->outbound_remote_flux_registers->insert(
          distributed_exchange->outbound_remote_flux_registers->end(),
          remote_entries.begin(),
          remote_entries.end());
    } else if (!remote_entries.empty()) {
      throw std::runtime_error("distributed AMR generated remote reflux entries without an exchange sink");
    }
  } else {
    local_entries = entries;
  }
  if (options.persist_incomplete_flux_registers) {
    diagnostics.pending_register_created_count = mergeFluxRegistersIntoPendingStore(state, local_entries, options);
    std::size_t complete_pending = 0;
    for (const core::PendingFluxRegisterRecord& pending : state.pending_flux_registers.records()) {
      if (pending.isComplete()) {
        ++complete_pending;
      }
    }
    diagnostics.pending_register_completed_count = complete_pending;
    diagnostics.pending_register_deferred_count = state.pending_flux_registers.size() - complete_pending;
    diagnostics.reflux = applyCompletePendingFluxRegistersToSimulationState(state, descriptors, options.adiabatic_index);
    diagnostics.pending_register_applied_count = diagnostics.reflux.complete_register_count;
    diagnostics.pending_register_rejected_count = diagnostics.reflux.skipped_area_mismatch_count +
        diagnostics.reflux.skipped_missing_target_count;
  } else {
    diagnostics.reflux = applyFluxRegistersToSimulationState(state, local_entries, descriptors, options.adiabatic_index);
  }
  return diagnostics;
}

}  // namespace

ProductionAmrHydroDiagnostics advanceProductionAmrHydro(
    core::SimulationState& state,
    std::span<const std::uint32_t> active_cell_rows,
    const hydro::HydroUpdateContext& update,
    const hydro::HydroSourceContext& global_source_context,
    const hydro::HydroCoreSolver& solver,
    const hydro::HydroRiemannSolver& riemann_solver,
    std::span<const hydro::HydroSourceTerm* const> source_terms,
    const ProductionAmrHydroOptions& options) {
  if (options.sweep_mode == ProductionAmrHydroSweepMode::kLocalSubcycling) {
    return advanceProductionAmrHydroSubcycled(
        state,
        active_cell_rows,
        update,
        global_source_context,
        solver,
        riemann_solver,
        source_terms,
        options);
  }
  return advanceProductionAmrHydroSynchronizedImpl(
      state,
      active_cell_rows,
      update,
      global_source_context,
      solver,
      riemann_solver,
      source_terms,
      options,
      nullptr);
}

ProductionAmrHydroDiagnostics advanceDistributedProductionAmrHydro(
    core::SimulationState& state,
    std::span<const std::uint32_t> active_cell_rows,
    const hydro::HydroUpdateContext& update,
    const hydro::HydroSourceContext& global_source_context,
    const hydro::HydroCoreSolver& solver,
    const hydro::HydroRiemannSolver& riemann_solver,
    std::span<const hydro::HydroSourceTerm* const> source_terms,
    const ProductionAmrHydroOptions& options,
    const DistributedAmrHydroExchange& exchange) {
  if (options.sweep_mode == ProductionAmrHydroSweepMode::kLocalSubcycling ||
      options.enable_temporal_coarse_to_fine) {
    throw std::runtime_error(
        "distributed AMR hydro currently supports synchronized patch execution only; "
        "MPI temporal subcycling is not implemented");
  }
  return advanceProductionAmrHydroSynchronizedImpl(
      state,
      active_cell_rows,
      update,
      global_source_context,
      solver,
      riemann_solver,
      source_terms,
      options,
      &exchange);
}

ProductionAmrHydroDiagnostics advanceProductionAmrHydroSubcycled(
    core::SimulationState& state,
    std::span<const std::uint32_t> active_cell_rows,
    const hydro::HydroUpdateContext& coarse_update,
    const hydro::HydroSourceContext& global_source_context,
    const hydro::HydroCoreSolver& solver,
    const hydro::HydroRiemannSolver& riemann_solver,
    std::span<const hydro::HydroSourceTerm* const> source_terms,
    const ProductionAmrHydroOptions& options) {
  if (!hasProductionAmrHydroCoverage(state)) {
    throw std::runtime_error("advanceProductionAmrHydroSubcycled requires complete AMR patch coverage in SimulationState");
  }
  if (options.refinement_ratio < 2U) {
    throw std::invalid_argument("AMR hydro subcycling requires refinement_ratio >= 2");
  }
  if (!std::isfinite(options.state_time_code)) {
    throw std::invalid_argument("AMR hydro subcycling requires finite state_time_code");
  }
  const std::vector<PatchDescriptor> descriptors = buildProductionAmrPatchDescriptors(state);
  int min_level = std::numeric_limits<int>::max();
  int max_level = std::numeric_limits<int>::min();
  for (const PatchDescriptor& patch : descriptors) {
    const int patch_level = static_cast<int>(patch.level);
    min_level = std::min(min_level, patch_level);
    max_level = std::max(max_level, patch_level);
  }
  if (min_level == std::numeric_limits<int>::max()) {
    return ProductionAmrHydroDiagnostics{};
  }
  // The temporal-history representation below records one coarse interval and
  // is intentionally scoped to a single local coarse/fine pair.  Rejecting a
  // deeper hierarchy is safer than silently reusing an interval from the wrong
  // parent fine substep; recursive multi-level temporal histories remain work
  // for a later scheduler-integrated implementation.
  if (max_level - min_level > 1) {
    throw std::runtime_error(
        "local temporal AMR subcycling currently supports one coarse/fine level pair; "
        "a deeper hierarchy requires nested temporal boundary histories");
  }
  if (!state.amr_temporal_boundary_history.empty()) {
    throw std::runtime_error(
        "cannot start a new local AMR subcycle while restart-resumable temporal boundary history is active");
  }

  ProductionAmrHydroDiagnostics diagnostics;
  diagnostics.patch_count = descriptors.size();
  diagnostics.levels_advanced = static_cast<std::size_t>(max_level - min_level + 1);
  diagnostics.max_level_advanced = static_cast<std::size_t>(std::max(0, max_level));
  diagnostics.substeps_by_level.assign(static_cast<std::size_t>(max_level + 1), 0U);

  const std::vector<std::uint32_t> coarse_rows = activeRowsForLevel(state, min_level, active_cell_rows);
  const std::vector<std::uint32_t> fine_rows = activeRowsForLevel(state, max_level, active_cell_rows);
  const double coarse_start_code = options.state_time_code;
  const double coarse_end_code = coarse_start_code + coarse_update.dt_code;
  captureAmrTemporalBoundaryHistoryStart(state, descriptors, coarse_start_code, options.adiabatic_index);

  if (!coarse_rows.empty()) {
    ProductionAmrHydroOptions coarse_options = options;
    coarse_options.sweep_mode = ProductionAmrHydroSweepMode::kSynchronized;
    coarse_options.persist_incomplete_flux_registers = true;
    coarse_options.reflux_interval_start_code = coarse_start_code;
    coarse_options.reflux_interval_end_code = coarse_end_code;
    coarse_options.reflux_coarse_dt_code = coarse_update.dt_code;
    coarse_options.expected_fine_substeps = options.refinement_ratio;
    coarse_options.fine_substep_index = 0U;
    coarse_options.state_time_code = coarse_start_code;
    coarse_options.ghost_fill_time_code = coarse_start_code;
    coarse_options.enable_temporal_coarse_to_fine = false;
    const ProductionAmrHydroDiagnostics coarse_diag = advanceProductionAmrHydro(
        state, coarse_rows, coarse_update, global_source_context, solver, riemann_solver, source_terms, coarse_options);
    accumulateDiagnostics(diagnostics, coarse_diag);
    diagnostics.substeps_by_level[static_cast<std::size_t>(min_level)] = 1U;
    diagnostics.subcycled_level_step_count += 1U;
  }
  captureAmrTemporalBoundaryHistoryEnd(state, descriptors, coarse_end_code, options.adiabatic_index);

  const double fine_dt_code = coarse_update.dt_code / static_cast<double>(options.refinement_ratio);
  for (std::uint32_t substep = 0; substep < options.refinement_ratio; ++substep) {
    if (fine_rows.empty()) {
      continue;
    }
    const double fine_start_code = coarse_start_code + static_cast<double>(substep) * fine_dt_code;
    hydro::HydroUpdateContext fine_update = coarse_update;
    fine_update.dt_code = fine_dt_code;
    ProductionAmrHydroOptions fine_options = options;
    fine_options.sweep_mode = ProductionAmrHydroSweepMode::kSynchronized;
    fine_options.persist_incomplete_flux_registers = true;
    fine_options.reflux_interval_start_code = coarse_start_code;
    fine_options.reflux_interval_end_code = coarse_end_code;
    fine_options.reflux_coarse_dt_code = coarse_update.dt_code;
    fine_options.expected_fine_substeps = options.refinement_ratio;
    fine_options.fine_substep_index = substep;
    fine_options.state_time_code = fine_start_code;
    // Ghosts are consumed before the predictor inside HydroCoreSolver, hence
    // this is the physical beginning of the fine substep.
    fine_options.ghost_fill_time_code = fine_start_code;
    fine_options.enable_temporal_coarse_to_fine = true;
    const ProductionAmrHydroDiagnostics fine_diag = advanceProductionAmrHydro(
        state, fine_rows, fine_update, global_source_context, solver, riemann_solver, source_terms, fine_options);
    accumulateDiagnostics(diagnostics, fine_diag);
    diagnostics.substeps_by_level[static_cast<std::size_t>(max_level)] += 1U;
    diagnostics.subcycled_level_step_count += 1U;
  }

  const RefluxDiagnostics post_child_reflux = applyCompletePendingFluxRegistersToSimulationState(
      state, buildProductionAmrPatchDescriptors(state), options.adiabatic_index);
  diagnostics.reflux.complete_register_count += post_child_reflux.complete_register_count;
  diagnostics.reflux.skipped_incomplete_register_count += post_child_reflux.skipped_incomplete_register_count;
  diagnostics.reflux.skipped_area_mismatch_count += post_child_reflux.skipped_area_mismatch_count;
  diagnostics.reflux.skipped_missing_target_count += post_child_reflux.skipped_missing_target_count;
  diagnostics.reflux.corrected_cells += post_child_reflux.corrected_cells;
  diagnostics.reflux.corrected_mass_code += post_child_reflux.corrected_mass_code;
  diagnostics.reflux.corrected_momentum_x_code += post_child_reflux.corrected_momentum_x_code;
  diagnostics.reflux.corrected_momentum_y_code += post_child_reflux.corrected_momentum_y_code;
  diagnostics.reflux.corrected_momentum_z_code += post_child_reflux.corrected_momentum_z_code;
  diagnostics.reflux.corrected_total_energy_code += post_child_reflux.corrected_total_energy_code;
  diagnostics.reflux.corrected_energy_code += post_child_reflux.corrected_energy_code;
  diagnostics.reflux.corrected_internal_energy_code += post_child_reflux.corrected_internal_energy_code;
  diagnostics.pending_register_applied_count += post_child_reflux.complete_register_count;
  diagnostics.pending_register_rejected_count += post_child_reflux.skipped_area_mismatch_count +
      post_child_reflux.skipped_missing_target_count;

  std::size_t complete_pending = 0;
  for (const core::PendingFluxRegisterRecord& pending : state.pending_flux_registers.records()) {
    if (pending.isComplete()) {
      ++complete_pending;
    }
  }
  diagnostics.pending_register_completed_count = std::max(diagnostics.pending_register_completed_count, complete_pending);
  diagnostics.pending_register_deferred_count = state.pending_flux_registers.size() - complete_pending;
  if (state.pending_flux_registers.empty()) {
    retireAmrTemporalBoundaryHistory(state);
  }
  return diagnostics;
}

ProductionAmrRegridDiagnostics refineProductionPatchInSimulationState(
    core::SimulationState& state,
    const PatchDescriptor& parent_patch,
    std::uint64_t first_child_patch_id,
    std::uint64_t first_child_gas_cell_id,
    const ProductionAmrHydroOptions& options) {
  if (!state.amr_temporal_boundary_history.empty()) {
    throw std::runtime_error(
        "production AMR refine is prohibited while a restart-resumable temporal boundary history is active");
  }
  if (first_child_patch_id == 0U || first_child_gas_cell_id == 0U) {
    throw std::invalid_argument("production AMR refine requires nonzero child patch and gas-cell ID seeds");
  }
  if (!hasProductionAmrHydroCoverage(state)) {
    throw std::runtime_error("production AMR refine requires complete SimulationState patch coverage");
  }
  const auto parent_patch_index = patchIndexById(state, parent_patch.patch_id);
  if (!parent_patch_index.has_value()) {
    throw std::runtime_error("production AMR refine target patch is absent from SimulationState");
  }
  const std::vector<std::uint32_t> parent_rows = patchLocalRowsForPatch(state, parent_patch, "production AMR refine");
  if (parent_rows.size() != product(parent_patch.cell_dims)) {
    throw std::runtime_error("production AMR refine parent row coverage does not match descriptor");
  }
  validatePatchIdRangeAvailable(state, first_child_patch_id, 8U);
  validateGasCellIdRangeAvailable(state, first_child_gas_cell_id, parent_rows.size() * 8U);
  const double parent_volume = parent_patch.extent_comov[0] * parent_patch.extent_comov[1] * parent_patch.extent_comov[2] /
      static_cast<double>(product(parent_patch.cell_dims));
  const auto before = totalsForRows(state, parent_rows, parent_volume, options);

  const core::SimulationState old = state;
  const std::array<double, 3> child_extent = {
      parent_patch.extent_comov[0] * 0.5,
      parent_patch.extent_comov[1] * 0.5,
      parent_patch.extent_comov[2] * 0.5};
  const double child_volume = parent_volume / 8.0;
  const std::size_t new_cell_count = old.cells.size() - parent_rows.size() + parent_rows.size() * 8U;
  const std::size_t new_patch_count = old.patches.size() - 1U + 8U;
  state.resizeCells(new_cell_count);
  state.resizePatches(new_patch_count);

  std::uint32_t write_row = 0;
  std::uint32_t write_patch = 0;
  std::vector<std::uint32_t> new_child_rows;
  new_child_rows.reserve(parent_rows.size() * 8U);
  for (std::size_t old_patch_index = 0; old_patch_index < old.patches.size(); ++old_patch_index) {
    if (old.patches.patch_id[old_patch_index] != parent_patch.patch_id) {
      const std::uint32_t first = old.patches.first_cell[old_patch_index];
      const std::uint32_t count = old.patches.cell_count[old_patch_index];
      writePatchDescriptorToStateRow(state, write_patch, patchDescriptorFromStateRow(old, old_patch_index));
      state.patches.first_cell[write_patch] = write_row;
      state.patches.cell_count[write_patch] = count;
      state.patches.owning_rank[write_patch] = old.patches.owning_rank[old_patch_index];
      for (std::uint32_t offset = 0; offset < count; ++offset) {
        copyCellRow(state, write_row, old, first + offset);
        state.cells.patch_index[write_row] = write_patch;
        ++write_row;
      }
      ++write_patch;
      continue;
    }

    for (std::uint8_t octant = 0; octant < 8U; ++octant) {
      const std::array<double, 3> child_origin = {
          parent_patch.origin_comov[0] + (((octant & 1U) != 0U) ? child_extent[0] : 0.0),
          parent_patch.origin_comov[1] + (((octant & 2U) != 0U) ? child_extent[1] : 0.0),
          parent_patch.origin_comov[2] + (((octant & 4U) != 0U) ? child_extent[2] : 0.0)};
      const std::array<double, 3> child_widths = {
          child_extent[0] / static_cast<double>(parent_patch.cell_dims[0]),
          child_extent[1] / static_cast<double>(parent_patch.cell_dims[1]),
          child_extent[2] / static_cast<double>(parent_patch.cell_dims[2])};
      PatchDescriptor child_patch = parent_patch;
      child_patch.patch_id = first_child_patch_id + octant;
      child_patch.parent_patch_id = parent_patch.patch_id;
      child_patch.level = parent_patch.level + 1U;
      child_patch.morton_key = parent_patch.morton_key * 8U + static_cast<std::uint64_t>(octant) + 1U;
      child_patch.origin_comov = child_origin;
      child_patch.extent_comov = child_extent;
      writePatchDescriptorToStateRow(state, write_patch, child_patch);
      state.patches.first_cell[write_patch] = write_row;
      state.patches.cell_count[write_patch] = static_cast<std::uint32_t>(parent_rows.size());
      state.patches.owning_rank[write_patch] = old.patches.owning_rank[old_patch_index];
      for (std::size_t child_cell = 0; child_cell < parent_rows.size(); ++child_cell) {
        const std::size_t i = child_cell % parent_patch.cell_dims[0];
        const std::size_t j = (child_cell / parent_patch.cell_dims[0]) % parent_patch.cell_dims[1];
        const std::size_t k = child_cell / (static_cast<std::size_t>(parent_patch.cell_dims[0]) * parent_patch.cell_dims[1]);
        const std::array<double, 3> center = {
            child_origin[0] + (static_cast<double>(i) + 0.5) * child_widths[0],
            child_origin[1] + (static_cast<double>(j) + 0.5) * child_widths[1],
            child_origin[2] + (static_cast<double>(k) + 0.5) * child_widths[2]};
        const std::size_t parent_cell = cellContainingPoint(parent_patch, center);
        const std::uint32_t parent_row = parent_rows[parent_cell];
        const ConservedState parent_conserved = volumeIntegratedForRow(old, parent_row, parent_volume, options);
        const ConservedState child_conserved = parent_conserved * (1.0 / 8.0);
        state.cells.center_x_comoving[write_row] = center[0];
        state.cells.center_y_comoving[write_row] = center[1];
        state.cells.center_z_comoving[write_row] = center[2];
        state.cells.time_bin[write_row] = old.cells.time_bin[parent_row];
        state.cells.patch_index[write_row] = write_patch;
        state.gas_cells.gas_cell_id[write_row] = first_child_gas_cell_id++;
        state.gas_cells.parent_particle_id[write_row] = old.gas_cells.parent_particle_id[parent_row];
        writeVolumeIntegratedToRow(state, write_row, child_conserved, child_volume, options);
        new_child_rows.push_back(write_row);
        ++write_row;
      }
      ++write_patch;
    }
  }
  if (write_row != new_cell_count || write_patch != new_patch_count) {
    throw std::runtime_error("production AMR refine internal row accounting failed");
  }
  rebuildIdentityFromSidecars(state);
  state.bumpCellIndexGeneration();
  const auto after = totalsForRows(state, new_child_rows, child_volume, options);
  return ProductionAmrRegridDiagnostics{
      .refined_patch_count = 1,
      .created_gas_cell_count = new_child_rows.size(),
      .retired_gas_cell_count = parent_rows.size(),
      .conserved_mass_before = before[0],
      .conserved_mass_after = after[0],
      .conserved_momentum_x_before = before[1],
      .conserved_momentum_x_after = after[1],
      .conserved_momentum_y_before = before[2],
      .conserved_momentum_y_after = after[2],
      .conserved_momentum_z_before = before[3],
      .conserved_momentum_z_after = after[3],
      .conserved_total_energy_before = before[4],
      .conserved_total_energy_after = after[4]};
}

ProductionAmrRegridDiagnostics refineProductionPatchInSimulationState(
    core::SimulationState& state,
    const PatchDescriptor& parent_patch,
    const ProductionAmrHydroOptions& options) {
  const std::size_t child_cell_count = product(parent_patch.cell_dims) * 8U;
  return refineProductionPatchInSimulationState(
      state, parent_patch, nextFreePatchId(state, 8U), nextFreeGasCellId(state, child_cell_count), options);
}

ProductionAmrRegridDiagnostics derefineProductionPatchInSimulationState(
    core::SimulationState& state,
    const PatchDescriptor& parent_patch,
    std::uint64_t replacement_gas_cell_id,
    const ProductionAmrHydroOptions& options) {
  if (!state.amr_temporal_boundary_history.empty()) {
    throw std::runtime_error(
        "production AMR derefine is prohibited while a restart-resumable temporal boundary history is active");
  }
  if (replacement_gas_cell_id == 0U) {
    throw std::invalid_argument("production AMR derefine requires a nonzero replacement gas_cell_id seed");
  }
  if (!hasProductionAmrHydroCoverage(state)) {
    throw std::runtime_error("production AMR derefine requires complete SimulationState patch coverage");
  }
  const std::vector<PatchDescriptor> descriptors = buildProductionAmrPatchDescriptors(state);
  std::vector<PatchDescriptor> children;
  for (const PatchDescriptor& patch : descriptors) {
    const bool parent_link_ok = patch.parent_patch_id == parent_patch.patch_id || patch.parent_patch_id == 0U;
    const bool same_parent_extent = parent_link_ok && patch.level == parent_patch.level + 1U &&
        patch.origin_comov[0] >= parent_patch.origin_comov[0] - k_geometry_tol &&
        patch.origin_comov[1] >= parent_patch.origin_comov[1] - k_geometry_tol &&
        patch.origin_comov[2] >= parent_patch.origin_comov[2] - k_geometry_tol &&
        patch.origin_comov[0] + patch.extent_comov[0] <= parent_patch.origin_comov[0] + parent_patch.extent_comov[0] + k_geometry_tol &&
        patch.origin_comov[1] + patch.extent_comov[1] <= parent_patch.origin_comov[1] + parent_patch.extent_comov[1] + k_geometry_tol &&
        patch.origin_comov[2] + patch.extent_comov[2] <= parent_patch.origin_comov[2] + parent_patch.extent_comov[2] + k_geometry_tol;
    if (same_parent_extent) {
      children.push_back(patch);
    }
  }
  const std::array<double, 3> expected_child_extent{
      parent_patch.extent_comov[0] * 0.5,
      parent_patch.extent_comov[1] * 0.5,
      parent_patch.extent_comov[2] * 0.5};
  std::array<bool, 8> octant_seen{};
  for (const PatchDescriptor& child : children) {
    if (child.level != parent_patch.level + 1U || child.cell_dims != parent_patch.cell_dims) {
      throw std::runtime_error("production AMR derefine child level/cell dimensions do not match exact octant contract");
    }
    std::uint8_t octant = 0U;
    for (std::size_t axis = 0; axis < 3; ++axis) {
      const double scale = std::max({1.0, std::abs(child.extent_comov[axis]), std::abs(expected_child_extent[axis])});
      if (std::abs(child.extent_comov[axis] - expected_child_extent[axis]) > k_geometry_tol * scale) {
        throw std::runtime_error("production AMR derefine child extent does not match exact octant contract");
      }
      const double lower = parent_patch.origin_comov[axis];
      const double upper = parent_patch.origin_comov[axis] + expected_child_extent[axis];
      const double child_origin = child.origin_comov[axis];
      const double lower_scale = std::max({1.0, std::abs(child_origin), std::abs(lower)});
      const double upper_scale = std::max({1.0, std::abs(child_origin), std::abs(upper)});
      if (std::abs(child_origin - lower) <= k_geometry_tol * lower_scale) {
        continue;
      }
      if (std::abs(child_origin - upper) <= k_geometry_tol * upper_scale) {
        octant |= static_cast<std::uint8_t>(1U << axis);
        continue;
      }
      throw std::runtime_error("production AMR derefine child origin does not match an exact parent octant");
    }
    if (octant_seen[octant]) {
      throw std::runtime_error("production AMR derefine detected duplicate child octant");
    }
    octant_seen[octant] = true;
  }
  if (children.size() != 8U ||
      !std::all_of(octant_seen.begin(), octant_seen.end(), [](bool seen) { return seen; })) {
    throw std::runtime_error("production AMR derefine requires exactly eight child patches covering the parent descriptor");
  }
  const double child_volume = children.front().extent_comov[0] * children.front().extent_comov[1] * children.front().extent_comov[2] /
      static_cast<double>(product(children.front().cell_dims));
  std::vector<std::uint32_t> child_rows;
  for (const PatchDescriptor& child : children) {
    const std::vector<std::uint32_t> rows = patchLocalRowsForPatch(state, child, "production AMR derefine child coverage");
    child_rows.insert(child_rows.end(), rows.begin(), rows.end());
  }
  const auto before = totalsForRows(state, child_rows, child_volume, options);
  validatePatchIdRangeAvailable(state, parent_patch.patch_id, 1U);
  validateGasCellIdRangeAvailable(state, replacement_gas_cell_id, product(parent_patch.cell_dims));

  const core::SimulationState old = state;
  const std::size_t parent_cell_count = product(parent_patch.cell_dims);
  const double parent_volume = parent_patch.extent_comov[0] * parent_patch.extent_comov[1] * parent_patch.extent_comov[2] /
      static_cast<double>(parent_cell_count);
  std::vector<ConservedState> restricted(parent_cell_count);
  std::vector<std::uint64_t> restricted_parent_particle_id(parent_cell_count, std::numeric_limits<std::uint64_t>::max());
  std::vector<std::uint8_t> restricted_time_bin(parent_cell_count, std::numeric_limits<std::uint8_t>::max());
  for (const PatchDescriptor& child : children) {
    const std::vector<std::uint32_t> rows = patchLocalRowsForPatch(old, child, "production AMR derefine child restriction");
    for (std::size_t child_cell = 0; child_cell < rows.size(); ++child_cell) {
      const std::array<double, 3> center = cellCenter(child, child_cell);
      const std::size_t parent_cell = cellContainingPoint(parent_patch, center);
      const std::uint32_t child_row = rows[child_cell];
      restricted[parent_cell] += volumeIntegratedForRow(old, child_row, child_volume, options);
      const std::uint64_t child_parent = old.gas_cells.parent_particle_id[child_row];
      std::uint64_t& parent_slot = restricted_parent_particle_id[parent_cell];
      if (parent_slot == std::numeric_limits<std::uint64_t>::max()) {
        parent_slot = child_parent;
      } else if (parent_slot != child_parent) {
        parent_slot = 0U;
      }
      restricted_time_bin[parent_cell] = std::min(restricted_time_bin[parent_cell], old.cells.time_bin[child_row]);
    }
  }

  std::unordered_set<std::uint64_t> child_patch_ids;
  for (const PatchDescriptor& child : children) {
    child_patch_ids.insert(child.patch_id);
  }
  const std::size_t kept_patch_count = old.patches.size() - child_patch_ids.size() + 1U;
  const std::size_t kept_cell_count = old.cells.size() - child_rows.size() + parent_cell_count;
  state.resizeCells(kept_cell_count);
  state.resizePatches(kept_patch_count);
  std::uint32_t write_row = 0;
  std::uint32_t write_patch = 0;
  std::vector<std::uint32_t> parent_rows;
  parent_rows.reserve(parent_cell_count);
  bool inserted_parent = false;
  for (std::size_t old_patch_index = 0; old_patch_index < old.patches.size(); ++old_patch_index) {
    const std::uint64_t old_patch_id = old.patches.patch_id[old_patch_index];
    if (child_patch_ids.contains(old_patch_id)) {
      if (!inserted_parent) {
        writePatchDescriptorToStateRow(state, write_patch, parent_patch);
        state.patches.first_cell[write_patch] = write_row;
        state.patches.cell_count[write_patch] = static_cast<std::uint32_t>(parent_cell_count);
        state.patches.owning_rank[write_patch] = old.patches.owning_rank[old_patch_index];
        for (std::size_t cell = 0; cell < parent_cell_count; ++cell) {
          const std::array<double, 3> center = cellCenter(parent_patch, cell);
          state.cells.center_x_comoving[write_row] = center[0];
          state.cells.center_y_comoving[write_row] = center[1];
          state.cells.center_z_comoving[write_row] = center[2];
          state.cells.time_bin[write_row] = restricted_time_bin[cell] == std::numeric_limits<std::uint8_t>::max()
              ? 0U
              : restricted_time_bin[cell];
          state.cells.patch_index[write_row] = write_patch;
          state.gas_cells.gas_cell_id[write_row] = replacement_gas_cell_id++;
          state.gas_cells.parent_particle_id[write_row] =
              restricted_parent_particle_id[cell] == std::numeric_limits<std::uint64_t>::max() ? 0U : restricted_parent_particle_id[cell];
          writeVolumeIntegratedToRow(state, write_row, restricted[cell], parent_volume, options);
          parent_rows.push_back(write_row);
          ++write_row;
        }
        ++write_patch;
        inserted_parent = true;
      }
      continue;
    }
    const std::uint32_t first = old.patches.first_cell[old_patch_index];
    const std::uint32_t count = old.patches.cell_count[old_patch_index];
    writePatchDescriptorToStateRow(state, write_patch, patchDescriptorFromStateRow(old, old_patch_index));
    state.patches.first_cell[write_patch] = write_row;
    state.patches.cell_count[write_patch] = count;
    state.patches.owning_rank[write_patch] = old.patches.owning_rank[old_patch_index];
    for (std::uint32_t offset = 0; offset < count; ++offset) {
      copyCellRow(state, write_row, old, first + offset);
      state.cells.patch_index[write_row] = write_patch;
      ++write_row;
    }
    ++write_patch;
  }
  if (!inserted_parent || write_row != kept_cell_count || write_patch != kept_patch_count) {
    throw std::runtime_error("production AMR derefine internal row accounting failed");
  }
  rebuildIdentityFromSidecars(state);
  state.bumpCellIndexGeneration();
  const auto after = totalsForRows(state, parent_rows, parent_volume, options);
  return ProductionAmrRegridDiagnostics{
      .derefined_patch_count = 1,
      .created_gas_cell_count = parent_rows.size(),
      .retired_gas_cell_count = child_rows.size(),
      .conserved_mass_before = before[0],
      .conserved_mass_after = after[0],
      .conserved_momentum_x_before = before[1],
      .conserved_momentum_x_after = after[1],
      .conserved_momentum_y_before = before[2],
      .conserved_momentum_y_after = after[2],
      .conserved_momentum_z_before = before[3],
      .conserved_momentum_z_after = after[3],
      .conserved_total_energy_before = before[4],
      .conserved_total_energy_after = after[4]};
}

ProductionAmrRegridDiagnostics derefineProductionPatchInSimulationState(
    core::SimulationState& state,
    const PatchDescriptor& parent_patch,
    const ProductionAmrHydroOptions& options) {
  return derefineProductionPatchInSimulationState(
      state, parent_patch, nextFreeGasCellId(state, product(parent_patch.cell_dims)), options);
}

}  // namespace cosmosim::amr
