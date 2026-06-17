#include <cassert>
#include <cmath>
#include <cstdint>
#include <optional>
#include <vector>

#include "cosmosim/amr/amr_ghost_fill.hpp"

namespace {

constexpr double k_gamma = 1.4;
constexpr double k_tol = 1.0e-12;

[[nodiscard]] std::size_t linearCellIndex(
    const cosmosim::amr::PatchDescriptor& patch,
    std::size_t i,
    std::size_t j,
    std::size_t k) {
  return i + static_cast<std::size_t>(patch.cell_dims[0]) *
                 (j + static_cast<std::size_t>(patch.cell_dims[1]) * k);
}

void assertStateNear(
    const cosmosim::hydro::HydroConservedState& actual,
    const cosmosim::hydro::HydroConservedState& expected) {
  assert(std::abs(actual.mass_density_comoving - expected.mass_density_comoving) < k_tol);
  assert(std::abs(actual.momentum_density_x_comoving - expected.momentum_density_x_comoving) < k_tol);
  assert(std::abs(actual.momentum_density_y_comoving - expected.momentum_density_y_comoving) < k_tol);
  assert(std::abs(actual.momentum_density_z_comoving - expected.momentum_density_z_comoving) < k_tol);
  assert(std::abs(actual.total_energy_density_comoving - expected.total_energy_density_comoving) < k_tol);
}

cosmosim::core::SimulationState makeState(
    const std::vector<cosmosim::amr::PatchDescriptor>& patches) {
  cosmosim::core::SimulationState state;
  std::size_t cell_count = 0;
  for (const auto& patch : patches) {
    cell_count += static_cast<std::size_t>(patch.cell_dims[0]) * patch.cell_dims[1] * patch.cell_dims[2];
  }
  state.resizeCells(cell_count);
  state.resizePatches(patches.size());

  std::vector<cosmosim::core::GasCellIdentityRecord> records;
  records.reserve(cell_count);
  std::uint32_t row = 0;
  std::uint64_t gas_cell_id = 5001;
  for (std::size_t patch_row = 0; patch_row < patches.size(); ++patch_row) {
    const auto& patch = patches[patch_row];
    const std::size_t patch_cell_count =
        static_cast<std::size_t>(patch.cell_dims[0]) * patch.cell_dims[1] * patch.cell_dims[2];
    state.patches.patch_id[patch_row] = patch.patch_id;
    state.patches.parent_patch_id[patch_row] = patch.parent_patch_id;
    state.patches.level[patch_row] = patch.level;
    state.patches.morton_key[patch_row] = patch.morton_key == 0U ? patch.patch_id : patch.morton_key;
    state.patches.origin_x_comoving[patch_row] = patch.origin_comov[0];
    state.patches.origin_y_comoving[patch_row] = patch.origin_comov[1];
    state.patches.origin_z_comoving[patch_row] = patch.origin_comov[2];
    state.patches.extent_x_comoving[patch_row] = patch.extent_comov[0];
    state.patches.extent_y_comoving[patch_row] = patch.extent_comov[1];
    state.patches.extent_z_comoving[patch_row] = patch.extent_comov[2];
    state.patches.cell_dim_x[patch_row] = patch.cell_dims[0];
    state.patches.cell_dim_y[patch_row] = patch.cell_dims[1];
    state.patches.cell_dim_z[patch_row] = patch.cell_dims[2];
    state.patches.first_cell[patch_row] = row;
    state.patches.cell_count[patch_row] = static_cast<std::uint32_t>(patch_cell_count);
    state.patches.owning_rank[patch_row] = 0;
    for (std::size_t patch_cell = 0; patch_cell < patch_cell_count; ++patch_cell) {
      const double value = 100.0 * static_cast<double>(patch_row + 1U) + static_cast<double>(patch_cell);
      const std::size_t i = patch_cell % patch.cell_dims[0];
      const std::size_t j = (patch_cell / patch.cell_dims[0]) % patch.cell_dims[1];
      const std::size_t k = patch_cell / (static_cast<std::size_t>(patch.cell_dims[0]) * patch.cell_dims[1]);
      state.cells.patch_index[row] = static_cast<std::uint32_t>(patch_row);
      state.cells.center_x_comoving[row] = patch.origin_comov[0] + (static_cast<double>(i) + 0.5) *
          patch.extent_comov[0] / static_cast<double>(patch.cell_dims[0]);
      state.cells.center_y_comoving[row] = patch.origin_comov[1] + (static_cast<double>(j) + 0.5) *
          patch.extent_comov[1] / static_cast<double>(patch.cell_dims[1]);
      state.cells.center_z_comoving[row] = patch.origin_comov[2] + (static_cast<double>(k) + 0.5) *
          patch.extent_comov[2] / static_cast<double>(patch.cell_dims[2]);
      state.cells.mass_code[row] = value;
      state.gas_cells.gas_cell_id[row] = gas_cell_id;
      state.gas_cells.parent_particle_id[row] = 0;
      state.gas_cells.density_code[row] = value;
      state.gas_cells.pressure_code[row] = 1.0;
      state.gas_cells.velocity_x_peculiar[row] = 0.0;
      state.gas_cells.velocity_y_peculiar[row] = 0.0;
      state.gas_cells.velocity_z_peculiar[row] = 0.0;
      records.push_back(cosmosim::core::GasCellIdentityRecord{
          .gas_cell_id = gas_cell_id,
          .parent_particle_id = std::nullopt,
          .owning_patch_id = patch.patch_id,
          .local_cell_row = row});
      ++row;
      ++gas_cell_id;
    }
  }
  state.gas_cell_identity.assign(std::move(records));
  return state;
}

cosmosim::amr::AmrHydroGeometryOptions sameLevelOptions(
    cosmosim::hydro::HydroFaceAxis axis,
    cosmosim::hydro::HydroFaceSide side) {
  cosmosim::amr::AmrHydroGeometryOptions options;
  options.physical_boundary_kind = cosmosim::hydro::HydroBoundaryKind::kOpen;
  const std::size_t slot =
      axis == cosmosim::hydro::HydroFaceAxis::kX
          ? (side == cosmosim::hydro::HydroFaceSide::kLower ? 0U : 1U)
          : axis == cosmosim::hydro::HydroFaceAxis::kY
              ? (side == cosmosim::hydro::HydroFaceSide::kLower ? 2U : 3U)
              : (side == cosmosim::hydro::HydroFaceSide::kLower ? 4U : 5U);
  options.boundary_classes[slot] = cosmosim::amr::AmrHydroBoundaryClass::kSameLevel;
  return options;
}

[[nodiscard]] const cosmosim::amr::AmrHydroGhostDescriptor& findGhost(
    const cosmosim::amr::AmrHydroPatchGeometry& geometry,
    cosmosim::hydro::HydroFaceAxis axis,
    cosmosim::hydro::HydroFaceSide side,
    std::size_t owner_cell) {
  for (const auto& descriptor : geometry.ghosts) {
    const auto& ghost = geometry.geometry.ghost_cells[descriptor.ghost_slot];
    if (ghost.axis == axis && ghost.side == side && ghost.owner_real_cell == owner_cell) {
      return descriptor;
    }
  }
  assert(false);
  return geometry.ghosts.front();
}

void testSameLevelGhostFillForAxis(cosmosim::hydro::HydroFaceAxis axis) {
  cosmosim::amr::PatchDescriptor left;
  left.patch_id = 11;
  left.level = 0;
  left.cell_dims = {2, 2, 2};
  left.extent_comov = {1.0, 1.0, 1.0};
  cosmosim::amr::PatchDescriptor right = left;
  right.patch_id = 12;
  if (axis == cosmosim::hydro::HydroFaceAxis::kX) {
    right.origin_comov[0] = 1.0;
  } else if (axis == cosmosim::hydro::HydroFaceAxis::kY) {
    right.origin_comov[1] = 1.0;
  } else {
    right.origin_comov[2] = 1.0;
  }

  const auto state = makeState({left, right});
  auto left_geometry = cosmosim::amr::buildAmrHydroPatchGeometry(
      state, left, sameLevelOptions(axis, cosmosim::hydro::HydroFaceSide::kUpper));
  auto right_geometry = cosmosim::amr::buildAmrHydroPatchGeometry(
      state, right, sameLevelOptions(axis, cosmosim::hydro::HydroFaceSide::kLower));

  auto left_conserved = cosmosim::amr::loadAmrHydroConservedState(state, left_geometry, k_gamma);
  auto right_conserved = cosmosim::amr::loadAmrHydroConservedState(state, right_geometry, k_gamma);
  const auto before_neighbor = right_conserved.loadCell(0);
  std::vector<cosmosim::amr::AmrHydroGhostFillPatch> patches{
      {.geometry = &left_geometry, .conserved = &left_conserved},
      {.geometry = &right_geometry, .conserved = &right_conserved}};

  const auto diagnostics = cosmosim::amr::fillAmrHydroGhostCells(patches, k_gamma);
  assert(diagnostics.same_level_ghosts_filled == 8U);

  const std::size_t left_owner =
      axis == cosmosim::hydro::HydroFaceAxis::kX
          ? linearCellIndex(left, 1, 0, 0)
          : axis == cosmosim::hydro::HydroFaceAxis::kY
              ? linearCellIndex(left, 0, 1, 0)
              : linearCellIndex(left, 0, 0, 1);
  const std::size_t right_source =
      axis == cosmosim::hydro::HydroFaceAxis::kX
          ? linearCellIndex(right, 0, 0, 0)
          : axis == cosmosim::hydro::HydroFaceAxis::kY
              ? linearCellIndex(right, 0, 0, 0)
              : linearCellIndex(right, 0, 0, 0);
  const auto& ghost = findGhost(
      left_geometry,
      axis,
      cosmosim::hydro::HydroFaceSide::kUpper,
      left_owner);
  assert(ghost.fill_status == cosmosim::amr::AmrHydroGhostFillStatus::kFilledSameLevel);
  assertStateNear(left_conserved.loadCell(ghost.ghost_cell), right_conserved.loadCell(right_source));
  assertStateNear(right_conserved.loadCell(0), before_neighbor);
}

void testSameLevelGhostFillInXyz() {
  testSameLevelGhostFillForAxis(cosmosim::hydro::HydroFaceAxis::kX);
  testSameLevelGhostFillForAxis(cosmosim::hydro::HydroFaceAxis::kY);
  testSameLevelGhostFillForAxis(cosmosim::hydro::HydroFaceAxis::kZ);
}

void testRejectsStaleRemoteGhostEpoch() {
  cosmosim::amr::PatchDescriptor patch;
  patch.patch_id = 21;
  patch.cell_dims = {2, 2, 2};
  const auto state = makeState({patch});
  auto geometry = cosmosim::amr::buildAmrHydroPatchGeometry(state, patch);
  auto conserved = cosmosim::amr::loadAmrHydroConservedState(state, geometry, k_gamma);
  std::vector<cosmosim::amr::AmrHydroGhostFillPatch> patches{{
      .geometry = &geometry,
      .conserved = &conserved,
      .residency = cosmosim::amr::AmrGhostSourceResidency::kRemoteReadOnly,
      .ghost_hydro_epoch = 4,
      .expected_ghost_hydro_epoch = 5,
  }};

  const auto diagnostics = cosmosim::amr::fillAmrHydroGhostCells(patches, k_gamma);
  assert(diagnostics.stale_epoch_rejections == geometry.ghosts.size());
  assert(diagnostics.unresolved_ghosts == geometry.ghosts.size());
  for (const auto& ghost : geometry.ghosts) {
    assert(ghost.fill_status == cosmosim::amr::AmrHydroGhostFillStatus::kRejectedStaleRemote);
  }
}

}  // namespace

int main() {
  testSameLevelGhostFillInXyz();
  testRejectsStaleRemoteGhostEpoch();
  return 0;
}
