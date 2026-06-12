#include <cassert>
#include <cmath>
#include <vector>

#include "cosmosim/hydro/hydro_boundary_conditions.hpp"
#include "cosmosim/hydro/hydro_cartesian_patch.hpp"
#include "cosmosim/hydro/hydro_riemann.hpp"
#include "cosmosim/hydro/hydro_reconstruction.hpp"

namespace {

constexpr double k_gamma = 1.4;
constexpr double k_tol = 1.0e-12;

[[nodiscard]] cosmosim::hydro::HydroPatchGeometry makeGeometry(
    cosmosim::hydro::HydroBoundaryKind boundary_kind) {
  auto geometry = cosmosim::hydro::makeCartesianPatchGeometry(cosmosim::hydro::HydroCartesianPatchSpec{
      .nx = 2,
      .ny = 2,
      .nz = 2,
      .origin_x_comoving = 0.0,
      .origin_y_comoving = 0.0,
      .origin_z_comoving = 0.0,
      .cell_width_x_comoving = 1.0,
      .cell_width_y_comoving = 2.0,
      .cell_width_z_comoving = 3.0});
  cosmosim::hydro::appendCartesianBoundaryGhostFaces(geometry, boundary_kind);
  return geometry;
}

[[nodiscard]] cosmosim::hydro::HydroPrimitiveState primitiveForCell(std::size_t cell_index) {
  const double row = static_cast<double>(cell_index);
  return cosmosim::hydro::HydroPrimitiveState{
      .rho_comoving = 1.0 + row,
      .vel_x_peculiar = 10.0 + row,
      .vel_y_peculiar = 20.0 + row,
      .vel_z_peculiar = 30.0 + row,
      .pressure_comoving = 2.0 + 0.5 * row};
}

void fillRealCells(
    cosmosim::hydro::HydroConservedStateSoa& conserved,
    const cosmosim::hydro::HydroPatchGeometry& geometry) {
  for (std::size_t cell = 0; cell < geometry.cellCount(); ++cell) {
    conserved.storeCell(
        cell,
        cosmosim::hydro::HydroCoreSolver::conservedFromPrimitive(primitiveForCell(cell), k_gamma));
  }
}

void assertClose(double lhs, double rhs) {
  assert(std::abs(lhs - rhs) < k_tol);
}

void assertPrimitiveMatches(
    const cosmosim::hydro::HydroPrimitiveState& actual,
    cosmosim::hydro::HydroPrimitiveState expected) {
  assertClose(actual.rho_comoving, expected.rho_comoving);
  assertClose(actual.vel_x_peculiar, expected.vel_x_peculiar);
  assertClose(actual.vel_y_peculiar, expected.vel_y_peculiar);
  assertClose(actual.vel_z_peculiar, expected.vel_z_peculiar);
  assertClose(actual.pressure_comoving, expected.pressure_comoving);
}

void reflectExpected(
    cosmosim::hydro::HydroPrimitiveState& expected,
    cosmosim::hydro::HydroFaceAxis axis) {
  switch (axis) {
    case cosmosim::hydro::HydroFaceAxis::kX:
      expected.vel_x_peculiar = -expected.vel_x_peculiar;
      break;
    case cosmosim::hydro::HydroFaceAxis::kY:
      expected.vel_y_peculiar = -expected.vel_y_peculiar;
      break;
    case cosmosim::hydro::HydroFaceAxis::kZ:
      expected.vel_z_peculiar = -expected.vel_z_peculiar;
      break;
  }
}

void assertBoundaryFill(cosmosim::hydro::HydroBoundaryKind boundary_kind) {
  const auto geometry = makeGeometry(boundary_kind);
  assert(geometry.ghost_cells.size() == 24U);
  cosmosim::hydro::HydroConservedStateSoa conserved(geometry.totalCellStorageCount());
  fillRealCells(conserved, geometry);

  const auto owner_before = conserved.loadCell(0);
  cosmosim::hydro::fillHydroBoundaryGhostCells(conserved, geometry, k_gamma);
  const auto owner_after = conserved.loadCell(0);
  assertClose(owner_before.mass_density_comoving, owner_after.mass_density_comoving);
  assertClose(owner_before.momentum_density_x_comoving, owner_after.momentum_density_x_comoving);
  assertClose(owner_before.momentum_density_y_comoving, owner_after.momentum_density_y_comoving);
  assertClose(owner_before.momentum_density_z_comoving, owner_after.momentum_density_z_comoving);
  assertClose(owner_before.total_energy_density_comoving, owner_after.total_energy_density_comoving);

  bool saw_x = false;
  bool saw_y = false;
  bool saw_z = false;
  for (const cosmosim::hydro::HydroGhostCell& ghost : geometry.ghost_cells) {
    assert(ghost.boundary_kind == boundary_kind);
    assert(ghost.owner_real_cell < geometry.cellCount());
    assert(ghost.source_real_cell < geometry.cellCount());
    assert(ghost.ghost_cell >= geometry.cellCount());
    assert(ghost.mutation_rights ==
        cosmosim::hydro::HydroGhostMutationRights::kWritablePhysicalBoundaryScratch);

    auto expected = primitiveForCell(ghost.source_real_cell);
    if (boundary_kind == cosmosim::hydro::HydroBoundaryKind::kReflective) {
      reflectExpected(expected, ghost.axis);
    }
    const auto actual = cosmosim::hydro::HydroCoreSolver::primitiveFromConserved(
        conserved.loadCell(ghost.ghost_cell),
        k_gamma);
    assertPrimitiveMatches(actual, expected);

    saw_x = saw_x || ghost.axis == cosmosim::hydro::HydroFaceAxis::kX;
    saw_y = saw_y || ghost.axis == cosmosim::hydro::HydroFaceAxis::kY;
    saw_z = saw_z || ghost.axis == cosmosim::hydro::HydroFaceAxis::kZ;
  }
  assert(saw_x && saw_y && saw_z);
}

void testPeriodicOpenReflectiveGhostRows() {
  assertBoundaryFill(cosmosim::hydro::HydroBoundaryKind::kPeriodic);
  assertBoundaryFill(cosmosim::hydro::HydroBoundaryKind::kOpen);
  assertBoundaryFill(cosmosim::hydro::HydroBoundaryKind::kReflective);
}

void testReconstructionConsumesGhostState() {
  auto geometry = makeGeometry(cosmosim::hydro::HydroBoundaryKind::kReflective);
  cosmosim::hydro::HydroConservedStateSoa conserved(geometry.totalCellStorageCount());
  fillRealCells(conserved, geometry);
  cosmosim::hydro::fillHydroBoundaryGhostCells(conserved, geometry, k_gamma);

  const cosmosim::hydro::HydroFace* boundary_face = nullptr;
  for (const cosmosim::hydro::HydroFace& face : geometry.faces) {
    if (face.ghost_cell_slot != cosmosim::hydro::k_invalid_ghost_cell_slot &&
        geometry.ghost_cells[face.ghost_cell_slot].axis == cosmosim::hydro::HydroFaceAxis::kY) {
      boundary_face = &face;
      break;
    }
  }
  assert(boundary_face != nullptr);

  cosmosim::hydro::PiecewiseConstantReconstruction reconstruction;
  cosmosim::hydro::HydroPrimitiveState left;
  cosmosim::hydro::HydroPrimitiveState right;
  reconstruction.reconstructFace(conserved, *boundary_face, left, right, k_gamma);

  const auto& ghost = geometry.ghost_cells[boundary_face->ghost_cell_slot];
  auto expected_right = primitiveForCell(ghost.source_real_cell);
  reflectExpected(expected_right, ghost.axis);
  assertPrimitiveMatches(left, primitiveForCell(ghost.owner_real_cell));
  assertPrimitiveMatches(right, expected_right);
}

void testImportedGhostMetadataIsReadOnly() {
  auto geometry = makeGeometry(cosmosim::hydro::HydroBoundaryKind::kOpen);
  geometry.ghost_cells.front().boundary_kind = cosmosim::hydro::HydroBoundaryKind::kImportedMpi;
  geometry.ghost_cells.front().mutation_rights =
      cosmosim::hydro::HydroGhostMutationRights::kReadOnlyImported;

  cosmosim::hydro::HydroConservedStateSoa conserved(geometry.totalCellStorageCount());
  fillRealCells(conserved, geometry);
  const auto sentinel = cosmosim::hydro::HydroCoreSolver::conservedFromPrimitive(
      {.rho_comoving = 99.0, .vel_x_peculiar = 1.0, .pressure_comoving = 7.0},
      k_gamma);
  conserved.storeCell(geometry.ghost_cells.front().ghost_cell, sentinel);

  cosmosim::hydro::fillHydroBoundaryGhostCells(conserved, geometry, k_gamma);
  const auto actual = conserved.loadCell(geometry.ghost_cells.front().ghost_cell);
  assertClose(actual.mass_density_comoving, sentinel.mass_density_comoving);
  assertClose(actual.total_energy_density_comoving, sentinel.total_energy_density_comoving);
}

[[nodiscard]] std::vector<std::size_t> indexRange(std::size_t count) {
  std::vector<std::size_t> indices(count);
  for (std::size_t i = 0; i < count; ++i) {
    indices[i] = i;
  }
  return indices;
}

void testPeriodicGhostFluxUpdatesWrappedRealCell() {
  auto geometry = makeGeometry(cosmosim::hydro::HydroBoundaryKind::kPeriodic);
  cosmosim::hydro::HydroConservedStateSoa conserved(geometry.totalCellStorageCount());
  fillRealCells(conserved, geometry);
  cosmosim::hydro::fillHydroBoundaryGhostCells(conserved, geometry, k_gamma);

  const std::vector<std::size_t> real_cells = indexRange(geometry.cellCount());
  const std::vector<std::size_t> active_faces = indexRange(geometry.faces.size());
  const cosmosim::hydro::HydroConservationTotals initial =
      cosmosim::hydro::HydroCoreSolver::conservationTotals(conserved, geometry, real_cells);

  cosmosim::hydro::HydroUpdateContext update{.dt_code = 2.5e-4, .scale_factor = 1.0, .hubble_rate_code = 0.0};
  cosmosim::hydro::HydroSourceContext source_context{.update = update};
  cosmosim::hydro::HydroCoreSolver solver(k_gamma);
  cosmosim::hydro::MusclHancockReconstruction reconstruction(cosmosim::hydro::HydroReconstructionPolicy{
      .limiter = cosmosim::hydro::HydroSlopeLimiter::kMonotonizedCentral,
      .dt_over_dx_code = update.dt_code / geometry.cell_width_x_comoving,
      .dt_over_cell_width_code = {
          update.dt_code / geometry.cell_width_x_comoving,
          update.dt_code / geometry.cell_width_y_comoving,
          update.dt_code / geometry.cell_width_z_comoving},
      .rho_floor = 1.0e-12,
      .pressure_floor = 1.0e-12,
      .enable_muscl_hancock_predictor = true,
      .adiabatic_index = k_gamma});
  cosmosim::hydro::HllcRiemannSolver riemann_solver;
  cosmosim::hydro::HydroScratchBuffers scratch;
  cosmosim::hydro::HydroPrimitiveCacheSoa primitive_cache(conserved.size());

  for (std::size_t step = 0; step < 4U; ++step) {
    cosmosim::hydro::fillHydroBoundaryGhostCells(conserved, geometry, k_gamma);
    solver.advancePatchActiveSetWithScratch(
        conserved,
        geometry,
        cosmosim::hydro::HydroActiveSetView{.active_cells = real_cells, .active_faces = active_faces},
        update,
        reconstruction,
        riemann_solver,
        {},
        source_context,
        scratch,
        &primitive_cache,
        nullptr);
  }

  const cosmosim::hydro::HydroConservationTotals final =
      cosmosim::hydro::HydroCoreSolver::conservationTotals(conserved, geometry, real_cells);
  assert(std::abs(final.mass - initial.mass) < 1.0e-10);
  assert(std::abs(final.momentum_x - initial.momentum_x) < 1.0e-10);
  assert(std::abs(final.momentum_y - initial.momentum_y) < 1.0e-10);
  assert(std::abs(final.momentum_z - initial.momentum_z) < 1.0e-10);
  assert(std::abs(final.total_energy - initial.total_energy) < 1.0e-10);
}

}  // namespace

int main() {
  testPeriodicOpenReflectiveGhostRows();
  testReconstructionConsumesGhostState();
  testImportedGhostMetadataIsReadOnly();
  testPeriodicGhostFluxUpdatesWrappedRealCell();
  return 0;
}
