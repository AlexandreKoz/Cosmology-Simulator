#include <cassert>
#include <cmath>
#include <vector>

#include "cosmosim/hydro/hydro_cartesian_patch.hpp"
#include "cosmosim/hydro/hydro_core_solver.hpp"

namespace {

constexpr double k_tol = 5.0e-10;

[[nodiscard]] cosmosim::hydro::HydroPatchGeometry makeGeometry() {
  cosmosim::hydro::HydroPatchGeometry geometry = cosmosim::hydro::makeCartesianPatchGeometry(cosmosim::hydro::HydroCartesianPatchSpec{
      .nx = 4,
      .ny = 3,
      .nz = 2,
      .origin_x_comoving = 0.0,
      .origin_y_comoving = 0.0,
      .origin_z_comoving = 0.0,
      .cell_width_x_comoving = 0.25,
      .cell_width_y_comoving = 1.0 / 3.0,
      .cell_width_z_comoving = 0.5});
  for (std::size_t k = 0; k < geometry.nz; ++k) {
    for (std::size_t j = 0; j < geometry.ny; ++j) {
      geometry.faces.push_back(cosmosim::hydro::HydroFace{
          .owner_cell = geometry.linearCellIndex(geometry.nx - 1U, j, k),
          .neighbor_cell = geometry.linearCellIndex(0U, j, k),
          .owner_minus_cell = geometry.linearCellIndex(geometry.nx - 2U, j, k),
          .neighbor_plus_cell = geometry.linearCellIndex(1U, j, k),
          .area_comoving = geometry.cell_width_y_comoving * geometry.cell_width_z_comoving,
          .normal_x = 1.0,
          .axis = cosmosim::hydro::HydroFaceAxis::kX});
    }
  }
  for (std::size_t k = 0; k < geometry.nz; ++k) {
    for (std::size_t i = 0; i < geometry.nx; ++i) {
      geometry.faces.push_back(cosmosim::hydro::HydroFace{
          .owner_cell = geometry.linearCellIndex(i, geometry.ny - 1U, k),
          .neighbor_cell = geometry.linearCellIndex(i, 0U, k),
          .owner_minus_cell = geometry.linearCellIndex(i, geometry.ny - 2U, k),
          .neighbor_plus_cell = geometry.linearCellIndex(i, 1U, k),
          .area_comoving = geometry.cell_width_x_comoving * geometry.cell_width_z_comoving,
          .normal_y = 1.0,
          .axis = cosmosim::hydro::HydroFaceAxis::kY});
    }
  }
  for (std::size_t j = 0; j < geometry.ny; ++j) {
    for (std::size_t i = 0; i < geometry.nx; ++i) {
      geometry.faces.push_back(cosmosim::hydro::HydroFace{
          .owner_cell = geometry.linearCellIndex(i, j, geometry.nz - 1U),
          .neighbor_cell = geometry.linearCellIndex(i, j, 0U),
          .owner_minus_cell = geometry.linearCellIndex(i, j, geometry.nz - 2U),
          .neighbor_plus_cell = geometry.linearCellIndex(i, j, 1U),
          .area_comoving = geometry.cell_width_x_comoving * geometry.cell_width_y_comoving,
          .normal_z = 1.0,
          .axis = cosmosim::hydro::HydroFaceAxis::kZ});
    }
  }
  return geometry;
}

void fillClosedBoxState(
    cosmosim::hydro::HydroConservedStateSoa& conserved,
    const cosmosim::hydro::HydroPatchGeometry& geometry,
    double gamma) {
  for (std::size_t cell = 0; cell < conserved.size(); ++cell) {
    cosmosim::hydro::HydroPrimitiveState primitive;
    primitive.rho_comoving = 1.25;
    primitive.vel_x_peculiar = 0.03;
    primitive.vel_y_peculiar = -0.02;
    primitive.vel_z_peculiar = 0.01;
    primitive.pressure_comoving = 1.1;
    conserved.storeCell(cell, cosmosim::hydro::HydroCoreSolver::conservedFromPrimitive(primitive, gamma));
  }
}

[[nodiscard]] std::vector<std::size_t> allCells(std::size_t count) {
  std::vector<std::size_t> cells(count);
  for (std::size_t i = 0; i < count; ++i) {
    cells[i] = i;
  }
  return cells;
}

void assertClose(double lhs, double rhs, double tolerance) {
  assert(std::abs(lhs - rhs) <= tolerance);
}

void testClosedBoxNoSourceConservation() {
  constexpr double gamma = 1.4;
  constexpr std::size_t step_count = 20;
  const cosmosim::hydro::HydroPatchGeometry geometry = makeGeometry();
  cosmosim::hydro::HydroConservedStateSoa conserved(geometry.cellCount());
  fillClosedBoxState(conserved, geometry, gamma);

  const std::vector<std::size_t> cells = allCells(conserved.size());
  const cosmosim::hydro::HydroConservationTotals initial =
      cosmosim::hydro::HydroCoreSolver::conservationTotals(conserved, geometry, cells);

  cosmosim::hydro::HydroUpdateContext update{.dt_code = 5.0e-4, .scale_factor = 1.0, .hubble_rate_code = 0.0};
  cosmosim::hydro::HydroSourceContext source_context{.update = update};
  cosmosim::hydro::HydroCoreSolver solver(gamma);
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
      .adiabatic_index = gamma});
  cosmosim::hydro::HllcRiemannSolver riemann_solver;
  cosmosim::hydro::HydroProfileEvent last_profile;

  for (std::size_t step = 0; step < step_count; ++step) {
    last_profile = {};
    solver.advancePatch(conserved, geometry, update, reconstruction, riemann_solver, {}, source_context, &last_profile);
    assert(last_profile.conservation.internal_energy_floor_count == 0U);
    assert(std::abs(last_profile.conservation.source_delta.mass) < k_tol);
    assert(std::abs(last_profile.conservation.source_delta.total_energy) < k_tol);
    assert(std::abs(last_profile.conservation.flux_delta.mass) < k_tol);
    assert(std::abs(last_profile.conservation.flux_delta.momentum_x) < k_tol);
    assert(std::abs(last_profile.conservation.flux_delta.momentum_y) < k_tol);
    assert(std::abs(last_profile.conservation.flux_delta.momentum_z) < k_tol);
    assert(std::abs(last_profile.conservation.flux_delta.total_energy) < k_tol);
    assert(std::abs(last_profile.conservation.residual.mass) < k_tol);
    assert(std::abs(last_profile.conservation.residual.total_energy) < k_tol);
    assert(std::abs(last_profile.conservation.residual.internal_energy) < k_tol);
  }

  const cosmosim::hydro::HydroConservationTotals final =
      cosmosim::hydro::HydroCoreSolver::conservationTotals(conserved, geometry, cells);
  assertClose(final.mass, initial.mass, k_tol);
  assertClose(final.momentum_x, initial.momentum_x, k_tol);
  assertClose(final.momentum_y, initial.momentum_y, k_tol);
  assertClose(final.momentum_z, initial.momentum_z, k_tol);
  assertClose(final.total_energy, initial.total_energy, k_tol);
  assertClose(final.internal_energy, initial.internal_energy, k_tol);
}

}  // namespace

int main() {
  testClosedBoxNoSourceConservation();
  return 0;
}
