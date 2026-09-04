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
  assert(std::abs(actual.metal_mass_density_comoving - expected.metal_mass_density_comoving) < k_tol);
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
      // Deliberately non-zero/non-uniform passive metal state so local/sparse
      // equivalence also covers distributed hydro metal transport semantics.
      state.gas_cells.metal_mass_code[row] =
          value * (0.01 + 0.001 * static_cast<double>((patch_cell % 7U) + 1U));
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

cosmosim::amr::AmrHydroGeometryOptions coarseFineOptions(
    cosmosim::hydro::HydroFaceSide side) {
  cosmosim::amr::AmrHydroGeometryOptions options;
  options.physical_boundary_kind = cosmosim::hydro::HydroBoundaryKind::kOpen;
  options.boundary_classes[
      side == cosmosim::hydro::HydroFaceSide::kLower ? 0U : 1U] =
      cosmosim::amr::AmrHydroBoundaryClass::kCoarseFine;
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

void testSparseSameLevelRemoteGhostMatchesDenseReference() {
  cosmosim::amr::PatchDescriptor left;
  left.patch_id = 31;
  left.level = 0;
  left.cell_dims = {2, 2, 2};
  left.extent_comov = {1.0, 1.0, 1.0};
  cosmosim::amr::PatchDescriptor right = left;
  right.patch_id = 32;
  right.origin_comov[0] = 1.0;

  // Dense two-patch reference using the established local/local path.
  const auto dense_state = makeState({left, right});
  auto dense_left_geometry = cosmosim::amr::buildAmrHydroPatchGeometry(
      dense_state, left, sameLevelOptions(
          cosmosim::hydro::HydroFaceAxis::kX,
          cosmosim::hydro::HydroFaceSide::kUpper));
  auto dense_right_geometry = cosmosim::amr::buildAmrHydroPatchGeometry(
      dense_state, right, sameLevelOptions(
          cosmosim::hydro::HydroFaceAxis::kX,
          cosmosim::hydro::HydroFaceSide::kLower));
  auto dense_left = cosmosim::amr::loadAmrHydroConservedState(
      dense_state, dense_left_geometry, k_gamma);
  auto dense_right = cosmosim::amr::loadAmrHydroConservedState(
      dense_state, dense_right_geometry, k_gamma);
  std::vector<cosmosim::amr::AmrHydroGhostFillPatch> dense_patches{
      {.geometry = &dense_left_geometry, .conserved = &dense_left},
      {.geometry = &dense_right_geometry, .conserved = &dense_right}};
  const auto dense_diagnostics = cosmosim::amr::fillAmrHydroGhostCells(
      dense_patches, k_gamma);
  assert(dense_diagnostics.same_level_ghosts_filled == 8U);

  // Sparse remote source contains only the X-lower boundary plane of the
  // right patch. No remote interior placeholders or remote solver geometry
  // exist in this path.
  const auto local_state = makeState({left});
  auto sparse_left_geometry = cosmosim::amr::buildAmrHydroPatchGeometry(
      local_state, left, sameLevelOptions(
          cosmosim::hydro::HydroFaceAxis::kX,
          cosmosim::hydro::HydroFaceSide::kUpper));
  auto sparse_left = cosmosim::amr::loadAmrHydroConservedState(
      local_state, sparse_left_geometry, k_gamma);
  std::vector<cosmosim::amr::AmrHydroSparseRemoteCell> sparse_cells;
  sparse_cells.reserve(4U);
  for (std::uint32_t k = 0; k < 2U; ++k) {
    for (std::uint32_t j = 0; j < 2U; ++j) {
      const std::uint32_t offset = static_cast<std::uint32_t>(
          linearCellIndex(right, 0U, j, k));
      sparse_cells.push_back(cosmosim::amr::AmrHydroSparseRemoteCell{
          .patch_local_cell = offset,
          .gas_cell_id = 9000U + offset,
          .conserved = dense_right.loadCell(offset)});
    }
  }
  const std::size_t sparse_capacity = sparse_cells.capacity();
  const cosmosim::amr::AmrHydroSparseGhostSource remote{
      .patch = &right,
      .owner_rank = 1,
      .ghost_hydro_epoch = 3U,
      .expected_ghost_hydro_epoch = 3U,
      .source_current_state_time_code = 0.0,
      .cells = sparse_cells};
  std::vector<cosmosim::amr::AmrHydroGhostFillPatch> sparse_patches{
      {.geometry = &sparse_left_geometry, .conserved = &sparse_left}};
  const auto sparse_diagnostics = cosmosim::amr::fillAmrHydroGhostCells(
      sparse_patches,
      std::span<const cosmosim::amr::AmrHydroSparseGhostSource>(&remote, 1U),
      k_gamma);
  assert(sparse_diagnostics.same_level_ghosts_filled == 4U);
  assert(sparse_cells.capacity() == sparse_capacity);

  for (std::size_t owner = 0; owner < 8U; ++owner) {
    const std::size_t i = owner % 2U;
    if (i != 1U) {
      continue;
    }
    const auto& dense_ghost = findGhost(
        dense_left_geometry,
        cosmosim::hydro::HydroFaceAxis::kX,
        cosmosim::hydro::HydroFaceSide::kUpper,
        owner);
    const auto& sparse_ghost = findGhost(
        sparse_left_geometry,
        cosmosim::hydro::HydroFaceAxis::kX,
        cosmosim::hydro::HydroFaceSide::kUpper,
        owner);
    assertStateNear(
        sparse_left.loadCell(sparse_ghost.ghost_cell),
        dense_left.loadCell(dense_ghost.ghost_cell));
  }
}

void testSparseFineToCoarseMatchesDenseReferenceAndRequiresCompleteStencil() {
  cosmosim::amr::PatchDescriptor coarse;
  coarse.patch_id = 41U;
  coarse.level = 0U;
  coarse.origin_comov = {0.0, 0.0, 0.0};
  coarse.extent_comov = {1.0, 1.0, 1.0};
  coarse.cell_dims = {2U, 2U, 2U};

  cosmosim::amr::PatchDescriptor fine;
  fine.patch_id = 42U;
  fine.parent_patch_id = coarse.patch_id;
  fine.level = 1U;
  fine.origin_comov = {1.0, 0.0, 0.0};
  fine.extent_comov = {0.5, 1.0, 1.0};
  fine.cell_dims = {2U, 4U, 4U};

  // Dense accepted reference. Both fine X layers deliberately carry distinct
  // values through makeState(), so averaging only the first layer cannot pass.
  const auto dense_state = makeState({coarse, fine});
  auto dense_coarse_geometry = cosmosim::amr::buildAmrHydroPatchGeometry(
      dense_state, coarse, coarseFineOptions(cosmosim::hydro::HydroFaceSide::kUpper));
  auto dense_fine_geometry = cosmosim::amr::buildAmrHydroPatchGeometry(
      dense_state, fine, coarseFineOptions(cosmosim::hydro::HydroFaceSide::kLower));
  auto dense_coarse = cosmosim::amr::loadAmrHydroConservedState(
      dense_state, dense_coarse_geometry, k_gamma);
  auto dense_fine = cosmosim::amr::loadAmrHydroConservedState(
      dense_state, dense_fine_geometry, k_gamma);
  std::vector<cosmosim::amr::AmrHydroGhostFillPatch> dense_patches{
      {.geometry = &dense_coarse_geometry, .conserved = &dense_coarse},
      {.geometry = &dense_fine_geometry, .conserved = &dense_fine}};
  const auto dense_diagnostics = cosmosim::amr::fillAmrHydroGhostCells(
      dense_patches, k_gamma);
  assert(dense_diagnostics.fine_to_coarse_ghosts_filled == 4U);

  const auto local_state = makeState({coarse});
  auto sparse_coarse_geometry = cosmosim::amr::buildAmrHydroPatchGeometry(
      local_state, coarse, coarseFineOptions(cosmosim::hydro::HydroFaceSide::kUpper));
  auto sparse_coarse = cosmosim::amr::loadAmrHydroConservedState(
      local_state, sparse_coarse_geometry, k_gamma);

  std::vector<cosmosim::amr::AmrHydroSparseRemoteCell> complete_cells;
  complete_cells.reserve(32U);
  for (std::uint32_t offset = 0U; offset < 32U; ++offset) {
    complete_cells.push_back(cosmosim::amr::AmrHydroSparseRemoteCell{
        .patch_local_cell = offset,
        .gas_cell_id = 12000U + offset,
        .conserved = dense_fine.loadCell(offset)});
  }
  const cosmosim::amr::AmrHydroSparseGhostSource complete_remote{
      .patch = &fine,
      .owner_rank = 1,
      .ghost_hydro_epoch = 9U,
      .expected_ghost_hydro_epoch = 9U,
      .source_current_state_time_code = 0.0,
      .cells = complete_cells};
  std::vector<cosmosim::amr::AmrHydroGhostFillPatch> sparse_patches{
      {.geometry = &sparse_coarse_geometry, .conserved = &sparse_coarse}};
  const auto sparse_diagnostics = cosmosim::amr::fillAmrHydroGhostCells(
      sparse_patches,
      std::span<const cosmosim::amr::AmrHydroSparseGhostSource>(&complete_remote, 1U),
      k_gamma);
  assert(sparse_diagnostics.fine_to_coarse_ghosts_filled == 4U);

  for (std::size_t owner = 0; owner < 8U; ++owner) {
    if ((owner % 2U) != 1U) {
      continue;
    }
    const auto& dense_ghost = findGhost(
        dense_coarse_geometry,
        cosmosim::hydro::HydroFaceAxis::kX,
        cosmosim::hydro::HydroFaceSide::kUpper,
        owner);
    const auto& sparse_ghost = findGhost(
        sparse_coarse_geometry,
        cosmosim::hydro::HydroFaceAxis::kX,
        cosmosim::hydro::HydroFaceSide::kUpper,
        owner);
    assertStateNear(
        sparse_coarse.loadCell(sparse_ghost.ghost_cell),
        dense_coarse.loadCell(dense_ghost.ghost_cell));
  }

  // Remove one cell from the second normal source layer. The sparse path must
  // reject incomplete restriction rather than averaging a partial stencil.
  std::vector<cosmosim::amr::AmrHydroSparseRemoteCell> incomplete_cells = complete_cells;
  const std::size_t missing_offset = linearCellIndex(fine, 1U, 0U, 0U);
  incomplete_cells.erase(
      incomplete_cells.begin() + static_cast<std::ptrdiff_t>(missing_offset));
  const cosmosim::amr::AmrHydroSparseGhostSource incomplete_remote{
      .patch = &fine,
      .owner_rank = 1,
      .ghost_hydro_epoch = 9U,
      .expected_ghost_hydro_epoch = 9U,
      .source_current_state_time_code = 0.0,
      .cells = incomplete_cells};
  auto rejected_coarse_geometry = cosmosim::amr::buildAmrHydroPatchGeometry(
      local_state, coarse, coarseFineOptions(cosmosim::hydro::HydroFaceSide::kUpper));
  auto rejected_coarse = cosmosim::amr::loadAmrHydroConservedState(
      local_state, rejected_coarse_geometry, k_gamma);
  std::vector<cosmosim::amr::AmrHydroGhostFillPatch> rejected_patches{
      {.geometry = &rejected_coarse_geometry, .conserved = &rejected_coarse}};
  bool threw = false;
  try {
    (void)cosmosim::amr::fillAmrHydroGhostCells(
        rejected_patches,
        std::span<const cosmosim::amr::AmrHydroSparseGhostSource>(&incomplete_remote, 1U),
        k_gamma);
  } catch (const std::runtime_error&) {
    threw = true;
  }
  assert(threw);
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
  testSparseSameLevelRemoteGhostMatchesDenseReference();
  testSparseFineToCoarseMatchesDenseReferenceAndRequiresCompleteStencil();
  testRejectsStaleRemoteGhostEpoch();
  return 0;
}
