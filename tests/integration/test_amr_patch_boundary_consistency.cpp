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

void assertNear(double actual, double expected) {
  assert(std::abs(actual - expected) < k_tol);
}

void assertStateNear(
    const cosmosim::hydro::HydroConservedState& actual,
    const cosmosim::hydro::HydroConservedState& expected) {
  assertNear(actual.mass_density_comoving, expected.mass_density_comoving);
  assertNear(actual.momentum_density_x_comoving, expected.momentum_density_x_comoving);
  assertNear(actual.momentum_density_y_comoving, expected.momentum_density_y_comoving);
  assertNear(actual.momentum_density_z_comoving, expected.momentum_density_z_comoving);
  assertNear(actual.total_energy_density_comoving, expected.total_energy_density_comoving);
}

cosmosim::core::SimulationState makeBoundaryState(
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
  std::uint64_t gas_cell_id = 7001;
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
      const double density = 100.0 * static_cast<double>(patch_row + 1U) + static_cast<double>(patch_cell);
      const double vx = 0.01 * static_cast<double>(patch_cell + 1U);
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
      state.cells.mass_code[row] = density;
      state.gas_cells.gas_cell_id[row] = gas_cell_id;
      state.gas_cells.parent_particle_id[row] = 0;
      state.gas_cells.density_code[row] = density;
      state.gas_cells.pressure_code[row] = 2.0 + 0.1 * static_cast<double>(patch_cell);
      state.gas_cells.velocity_x_peculiar[row] = vx;
      state.gas_cells.velocity_y_peculiar[row] = 0.5 * vx;
      state.gas_cells.velocity_z_peculiar[row] = -0.25 * vx;
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

[[nodiscard]] cosmosim::amr::AmrHydroGeometryOptions coarseFineOptions(
    cosmosim::hydro::HydroFaceSide side,
    cosmosim::hydro::HydroBoundaryKind physical_kind = cosmosim::hydro::HydroBoundaryKind::kOpen) {
  cosmosim::amr::AmrHydroGeometryOptions options;
  options.physical_boundary_kind = physical_kind;
  options.boundary_classes[side == cosmosim::hydro::HydroFaceSide::kLower ? 0U : 1U] =
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

void testCoarseFineBoundaryFill() {
  cosmosim::amr::PatchDescriptor coarse;
  coarse.patch_id = 31;
  coarse.level = 0;
  coarse.origin_comov = {0.0, 0.0, 0.0};
  coarse.extent_comov = {1.0, 1.0, 1.0};
  coarse.cell_dims = {2, 2, 2};

  cosmosim::amr::PatchDescriptor fine;
  fine.patch_id = 32;
  fine.level = 1;
  fine.origin_comov = {1.0, 0.0, 0.0};
  fine.extent_comov = {0.5, 1.0, 1.0};
  fine.cell_dims = {2, 4, 4};

  const auto state = makeBoundaryState({coarse, fine});
  auto coarse_geometry = cosmosim::amr::buildAmrHydroPatchGeometry(
      state, coarse, coarseFineOptions(cosmosim::hydro::HydroFaceSide::kUpper));
  auto fine_geometry = cosmosim::amr::buildAmrHydroPatchGeometry(
      state, fine, coarseFineOptions(cosmosim::hydro::HydroFaceSide::kLower));
  auto coarse_conserved = cosmosim::amr::loadAmrHydroConservedState(state, coarse_geometry, k_gamma);
  auto fine_conserved = cosmosim::amr::loadAmrHydroConservedState(state, fine_geometry, k_gamma);

  const auto coarse_before = coarse_conserved.loadCell(linearCellIndex(coarse, 1, 0, 0));
  const auto fine_before = fine_conserved.loadCell(linearCellIndex(fine, 0, 0, 0));
  std::vector<cosmosim::amr::AmrHydroGhostFillPatch> patches{
      {.geometry = &coarse_geometry, .conserved = &coarse_conserved},
      {.geometry = &fine_geometry, .conserved = &fine_conserved}};
  const auto diagnostics = cosmosim::amr::fillAmrHydroGhostCells(patches, k_gamma);

  assert(diagnostics.coarse_to_fine_ghosts_filled == 16U);
  assert(diagnostics.fine_to_coarse_ghosts_filled == 4U);

  const auto& fine_lower_ghost = findGhost(
      fine_geometry,
      cosmosim::hydro::HydroFaceAxis::kX,
      cosmosim::hydro::HydroFaceSide::kLower,
      linearCellIndex(fine, 0, 0, 0));
  assertStateNear(fine_conserved.loadCell(fine_lower_ghost.ghost_cell), coarse_before);

  const auto& coarse_upper_ghost = findGhost(
      coarse_geometry,
      cosmosim::hydro::HydroFaceAxis::kX,
      cosmosim::hydro::HydroFaceSide::kUpper,
      linearCellIndex(coarse, 1, 0, 0));
  cosmosim::hydro::HydroConservedState fine_average;
  for (std::size_t k = 0; k < 2; ++k) {
    for (std::size_t j = 0; j < 2; ++j) {
      for (std::size_t i = 0; i < 2; ++i) {
        fine_average += fine_conserved.loadCell(linearCellIndex(fine, i, j, k));
      }
    }
  }
  fine_average = 0.125 * fine_average;
  assertStateNear(coarse_conserved.loadCell(coarse_upper_ghost.ghost_cell), fine_average);

  assertStateNear(coarse_conserved.loadCell(linearCellIndex(coarse, 1, 0, 0)), coarse_before);
  assertStateNear(fine_conserved.loadCell(linearCellIndex(fine, 0, 0, 0)), fine_before);
}

void testPhysicalReflectiveBoundaryStillUsesH1Rules() {
  cosmosim::amr::PatchDescriptor patch;
  patch.patch_id = 41;
  patch.cell_dims = {2, 2, 2};
  const auto state = makeBoundaryState({patch});
  cosmosim::amr::AmrHydroGeometryOptions options;
  options.physical_boundary_kind = cosmosim::hydro::HydroBoundaryKind::kReflective;
  auto geometry = cosmosim::amr::buildAmrHydroPatchGeometry(state, patch, options);
  auto conserved = cosmosim::amr::loadAmrHydroConservedState(state, geometry, k_gamma);

  std::vector<cosmosim::amr::AmrHydroGhostFillPatch> patches{{.geometry = &geometry, .conserved = &conserved}};
  const auto source = cosmosim::hydro::HydroCoreSolver::primitiveFromConserved(conserved.loadCell(0), k_gamma);
  (void)cosmosim::amr::fillAmrHydroGhostCells(patches, k_gamma);
  const auto& ghost = findGhost(
      geometry,
      cosmosim::hydro::HydroFaceAxis::kX,
      cosmosim::hydro::HydroFaceSide::kLower,
      0);
  const auto primitive =
      cosmosim::hydro::HydroCoreSolver::primitiveFromConserved(conserved.loadCell(ghost.ghost_cell), k_gamma);

  assertNear(primitive.rho_comoving, source.rho_comoving);
  assertNear(primitive.vel_x_peculiar, -source.vel_x_peculiar);
  assertNear(primitive.vel_y_peculiar, source.vel_y_peculiar);
  assertNear(primitive.vel_z_peculiar, source.vel_z_peculiar);
  assertNear(primitive.pressure_comoving, source.pressure_comoving);
}

}  // namespace

int main() {
  testCoarseFineBoundaryFill();
  testPhysicalReflectiveBoundaryStillUsesH1Rules();
  return 0;
}
