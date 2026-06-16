#include "cosmosim/amr/amr_hydro_orchestrator.hpp"

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

[[nodiscard]] std::vector<std::uint32_t> sortedRowsForPatch(
    const core::SimulationState& state,
    std::uint64_t patch_id) {
  std::vector<std::uint32_t> rows = state.gas_cell_identity.rowsForPatch(patch_id);
  std::sort(rows.begin(), rows.end());
  return rows;
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
  state.gas_cell_identity.assign(std::move(records));
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
    throw std::runtime_error("buildProductionAmrPatchDescriptors requires complete SimulationState AMR patch coverage");
  }
  std::vector<PatchDescriptor> descriptors;
  descriptors.reserve(state.patches.size());
  for (std::size_t patch_index = 0; patch_index < state.patches.size(); ++patch_index) {
    const std::uint32_t first = state.patches.first_cell[patch_index];
    const std::uint32_t count = state.patches.cell_count[patch_index];
    if (count == 0U) {
      continue;
    }
    std::vector<double> xs;
    std::vector<double> ys;
    std::vector<double> zs;
    xs.reserve(count);
    ys.reserve(count);
    zs.reserve(count);
    for (std::uint32_t offset = 0; offset < count; ++offset) {
      const std::uint32_t row = first + offset;
      xs.push_back(state.cells.center_x_comoving[row]);
      ys.push_back(state.cells.center_y_comoving[row]);
      zs.push_back(state.cells.center_z_comoving[row]);
    }
    xs = sortedUnique(std::move(xs));
    ys = sortedUnique(std::move(ys));
    zs = sortedUnique(std::move(zs));
    std::array<std::uint16_t, 3> dims{};
    if (xs.size() * ys.size() * zs.size() == count && !xs.empty() && !ys.empty() && !zs.empty()) {
      dims = {static_cast<std::uint16_t>(xs.size()), static_cast<std::uint16_t>(ys.size()), static_cast<std::uint16_t>(zs.size())};
    } else {
      dims = nearCubicFactors(count);
    }
    const double dx = spacingForUnique(xs, 1.0 / static_cast<double>(std::max<std::uint16_t>(dims[0], 1U)));
    const double dy = spacingForUnique(ys, 1.0 / static_cast<double>(std::max<std::uint16_t>(dims[1], 1U)));
    const double dz = spacingForUnique(zs, 1.0 / static_cast<double>(std::max<std::uint16_t>(dims[2], 1U)));
    PatchDescriptor descriptor;
    descriptor.patch_id = state.patches.patch_id[patch_index];
    descriptor.parent_patch_id = 0U;
    descriptor.level = static_cast<std::uint8_t>(std::max<std::int32_t>(state.patches.level[patch_index], 0));
    descriptor.morton_key = descriptor.patch_id;
    descriptor.cell_dims = dims;
    descriptor.origin_comov = {
        xs.empty() ? 0.0 : xs.front() - 0.5 * dx,
        ys.empty() ? 0.0 : ys.front() - 0.5 * dy,
        zs.empty() ? 0.0 : zs.front() - 0.5 * dz,
    };
    descriptor.extent_comov = {
        dx * static_cast<double>(dims[0]),
        dy * static_cast<double>(dims[1]),
        dz * static_cast<double>(dims[2]),
    };
    if (product(descriptor.cell_dims) != count) {
      throw std::runtime_error("buildProductionAmrPatchDescriptors could not derive patch-local cell dimensions");
    }
    descriptors.push_back(descriptor);
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
    const auto patch_it = std::find_if(
        all_patches.begin(),
        all_patches.end(),
        [&entry](const PatchDescriptor& patch) { return patch.patch_id == entry.coarse_patch_id; });
    if (patch_it == all_patches.end()) {
      continue;
    }
    const std::vector<std::uint32_t> rows = sortedRowsForPatch(state, entry.coarse_patch_id);
    if (entry.coarse_cell_index >= rows.size()) {
      continue;
    }
    const std::uint32_t row = rows[entry.coarse_cell_index];
    const double volume = patch_it->extent_comov[0] * patch_it->extent_comov[1] * patch_it->extent_comov[2] /
        static_cast<double>(std::max<std::size_t>(product(patch_it->cell_dims), 1U));
    const hydro::HydroPrimitiveState old_primitive = primitiveForRow(
        state, row, 1.0e-14, 1.0e-14, adiabatic_index);
    hydro::HydroConservedState conserved =
        hydro::HydroCoreSolver::conservedFromPrimitive(old_primitive, adiabatic_index);
    ConservedState delta_flux = entry.fine_face_flux_code - entry.coarse_face_flux_code;
    delta_flux *= (entry.face_area_comov * entry.dt_code / volume);
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
    diagnostics.corrected_energy_code += std::abs(delta_flux.total_energy_code * volume);
    diagnostics.corrected_internal_energy_code += std::abs((new_internal_density - old_internal_density) * volume);
  }
  return diagnostics;
}

ProductionAmrHydroDiagnostics advanceProductionAmrHydro(
    core::SimulationState& state,
    std::span<const std::uint32_t> active_cell_rows,
    const hydro::HydroUpdateContext& update,
    const hydro::HydroSourceContext& global_source_context,
    const hydro::HydroCoreSolver& solver,
    const hydro::HydroRiemannSolver& riemann_solver,
    std::span<const hydro::HydroSourceTerm* const> source_terms,
    const ProductionAmrHydroOptions& options) {
  if (!hasProductionAmrHydroCoverage(state)) {
    throw std::runtime_error("advanceProductionAmrHydro requires complete AMR patch coverage in SimulationState");
  }
  ProductionAmrHydroDiagnostics diagnostics;
  const std::vector<PatchDescriptor> descriptors = buildProductionAmrPatchDescriptors(state);
  diagnostics.patch_count = descriptors.size();

  std::unordered_set<std::uint32_t> active_lookup(active_cell_rows.begin(), active_cell_rows.end());
  const bool all_cells_active = active_cell_rows.empty();

  std::vector<AmrHydroPatchGeometry> geometries;
  std::vector<hydro::HydroConservedStateSoa> conserved_states;
  geometries.reserve(descriptors.size());
  conserved_states.reserve(descriptors.size());
  for (const PatchDescriptor& patch : descriptors) {
    AmrHydroGeometryOptions geometry_options;
    geometry_options.physical_boundary_kind = options.physical_boundary_kind;
    geometry_options.boundary_classes = {
        boundaryClassForSide(patch, descriptors, hydro::HydroFaceAxis::kX, hydro::HydroFaceSide::kLower),
        boundaryClassForSide(patch, descriptors, hydro::HydroFaceAxis::kX, hydro::HydroFaceSide::kUpper),
        boundaryClassForSide(patch, descriptors, hydro::HydroFaceAxis::kY, hydro::HydroFaceSide::kLower),
        boundaryClassForSide(patch, descriptors, hydro::HydroFaceAxis::kY, hydro::HydroFaceSide::kUpper),
        boundaryClassForSide(patch, descriptors, hydro::HydroFaceAxis::kZ, hydro::HydroFaceSide::kLower),
        boundaryClassForSide(patch, descriptors, hydro::HydroFaceAxis::kZ, hydro::HydroFaceSide::kUpper),
    };
    geometries.push_back(buildAmrHydroPatchGeometry(state, patch, geometry_options));
    populateAmrHydroFluxRegisterFaces(geometries.back(), descriptors);
    conserved_states.push_back(loadAmrHydroConservedState(state, geometries.back(), options.adiabatic_index));
  }

  refreshTraceableGhostSourceIds(geometries);

  std::vector<AmrHydroGhostFillPatch> ghost_views;
  ghost_views.reserve(geometries.size());
  for (std::size_t i = 0; i < geometries.size(); ++i) {
    ghost_views.push_back(AmrHydroGhostFillPatch{
        .geometry = &geometries[i],
        .conserved = &conserved_states[i],
        .residency = AmrGhostSourceResidency::kLocal,
        .ghost_hydro_epoch = state.gasCellIdentityGeneration(),
        .expected_ghost_hydro_epoch = state.gasCellIdentityGeneration()});
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
  diagnostics.reflux = applyFluxRegistersToSimulationState(state, entries, descriptors, options.adiabatic_index);
  return diagnostics;
}

ProductionAmrRegridDiagnostics refineProductionPatchInSimulationState(
    core::SimulationState& state,
    const PatchDescriptor& parent_patch,
    std::uint64_t first_child_patch_id,
    std::uint64_t first_child_gas_cell_id,
    const ProductionAmrHydroOptions& options) {
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
  const std::vector<std::uint32_t> parent_rows = sortedRowsForPatch(state, parent_patch.patch_id);
  if (parent_rows.size() != product(parent_patch.cell_dims)) {
    throw std::runtime_error("production AMR refine parent row coverage does not match descriptor");
  }
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
      state.patches.patch_id[write_patch] = old.patches.patch_id[old_patch_index];
      state.patches.level[write_patch] = old.patches.level[old_patch_index];
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
      state.patches.patch_id[write_patch] = first_child_patch_id + octant;
      state.patches.level[write_patch] = static_cast<std::int32_t>(parent_patch.level) + 1;
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

ProductionAmrRegridDiagnostics derefineProductionPatchInSimulationState(
    core::SimulationState& state,
    const PatchDescriptor& parent_patch,
    std::uint64_t replacement_gas_cell_id,
    const ProductionAmrHydroOptions& options) {
  if (replacement_gas_cell_id == 0U) {
    throw std::invalid_argument("production AMR derefine requires a nonzero replacement gas_cell_id seed");
  }
  if (!hasProductionAmrHydroCoverage(state)) {
    throw std::runtime_error("production AMR derefine requires complete SimulationState patch coverage");
  }
  const std::vector<PatchDescriptor> descriptors = buildProductionAmrPatchDescriptors(state);
  std::vector<PatchDescriptor> children;
  for (const PatchDescriptor& patch : descriptors) {
    const bool same_parent_extent = patch.level == parent_patch.level + 1U &&
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
  if (children.size() != 8U) {
    throw std::runtime_error("production AMR derefine requires exactly eight child patches covering the parent descriptor");
  }
  const double child_volume = children.front().extent_comov[0] * children.front().extent_comov[1] * children.front().extent_comov[2] /
      static_cast<double>(product(children.front().cell_dims));
  std::vector<std::uint32_t> child_rows;
  for (const PatchDescriptor& child : children) {
    const std::vector<std::uint32_t> rows = sortedRowsForPatch(state, child.patch_id);
    child_rows.insert(child_rows.end(), rows.begin(), rows.end());
  }
  const auto before = totalsForRows(state, child_rows, child_volume, options);

  const core::SimulationState old = state;
  const std::size_t parent_cell_count = product(parent_patch.cell_dims);
  const double parent_volume = parent_patch.extent_comov[0] * parent_patch.extent_comov[1] * parent_patch.extent_comov[2] /
      static_cast<double>(parent_cell_count);
  std::vector<ConservedState> restricted(parent_cell_count);
  for (const PatchDescriptor& child : children) {
    const std::vector<std::uint32_t> rows = sortedRowsForPatch(old, child.patch_id);
    for (std::size_t child_cell = 0; child_cell < rows.size(); ++child_cell) {
      const std::array<double, 3> center = cellCenter(child, child_cell);
      const std::size_t parent_cell = cellContainingPoint(parent_patch, center);
      restricted[parent_cell] += volumeIntegratedForRow(old, rows[child_cell], child_volume, options);
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
        state.patches.patch_id[write_patch] = parent_patch.patch_id;
        state.patches.level[write_patch] = parent_patch.level;
        state.patches.first_cell[write_patch] = write_row;
        state.patches.cell_count[write_patch] = static_cast<std::uint32_t>(parent_cell_count);
        state.patches.owning_rank[write_patch] = old.patches.owning_rank[old_patch_index];
        for (std::size_t cell = 0; cell < parent_cell_count; ++cell) {
          const std::array<double, 3> center = cellCenter(parent_patch, cell);
          state.cells.center_x_comoving[write_row] = center[0];
          state.cells.center_y_comoving[write_row] = center[1];
          state.cells.center_z_comoving[write_row] = center[2];
          state.cells.time_bin[write_row] = 0U;
          state.cells.patch_index[write_row] = write_patch;
          state.gas_cells.gas_cell_id[write_row] = replacement_gas_cell_id++;
          state.gas_cells.parent_particle_id[write_row] = 0U;
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
    state.patches.patch_id[write_patch] = old.patches.patch_id[old_patch_index];
    state.patches.level[write_patch] = old.patches.level[old_patch_index];
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

}  // namespace cosmosim::amr
