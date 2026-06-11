#include <cassert>
#include <cmath>
#include <vector>

#include "cosmosim/hydro/hydro_cartesian_patch.hpp"
#include "cosmosim/hydro/hydro_core_solver.hpp"

namespace {

constexpr double k_tol = 1.0e-9;

cosmosim::hydro::HydroPatchGeometry makeCartesianGeometry(
    std::size_t nx,
    std::size_t ny,
    std::size_t nz) {
  return cosmosim::hydro::makeCartesianPatchGeometry(cosmosim::hydro::HydroCartesianPatchSpec{
      .nx = nx,
      .ny = ny,
      .nz = nz,
      .origin_x_comoving = 0.0,
      .origin_y_comoving = 0.0,
      .origin_z_comoving = 0.0,
      .cell_width_x_comoving = 1.0 / static_cast<double>(nx),
      .cell_width_y_comoving = 1.0 / static_cast<double>(ny),
      .cell_width_z_comoving = 1.0 / static_cast<double>(nz)});
}

void fillSodLikeInitialState(
    cosmosim::hydro::HydroConservedStateSoa& conserved,
    const cosmosim::hydro::HydroPatchGeometry& geometry,
    double gamma) {
  for (std::size_t i = 0; i < conserved.size(); ++i) {
    const auto ijk = geometry.cellIjk(i);
    cosmosim::hydro::HydroPrimitiveState primitive;
    if (ijk[0] < geometry.nx / 2U) {
      primitive.rho_comoving = 1.0;
      primitive.pressure_comoving = 1.0;
    } else {
      primitive.rho_comoving = 0.125;
      primitive.pressure_comoving = 0.1;
    }

    conserved.storeCell(i, cosmosim::hydro::HydroCoreSolver::conservedFromPrimitive(primitive, gamma));
  }
}

void testSodLikePeriodicConservationAndRegression() {
  constexpr std::size_t k_nx = 4;
  constexpr std::size_t k_ny = 4;
  constexpr std::size_t k_nz = 4;
  constexpr std::size_t k_cell_count = k_nx * k_ny * k_nz;
  constexpr std::size_t k_step_count = 30;
  constexpr double k_gamma = 1.4;

  const cosmosim::hydro::HydroPatchGeometry geometry = makeCartesianGeometry(k_nx, k_ny, k_nz);
  cosmosim::hydro::HydroConservedStateSoa conserved(k_cell_count);
  fillSodLikeInitialState(conserved, geometry, k_gamma);

  cosmosim::hydro::HydroUpdateContext update;
  update.dt_code = 1.0e-3;
  update.scale_factor = 1.0;

  cosmosim::hydro::HydroSourceContext source_context;
  source_context.update = update;

  cosmosim::hydro::HydroCoreSolver solver(k_gamma);
  cosmosim::hydro::MusclHancockReconstruction reconstruction(cosmosim::hydro::HydroReconstructionPolicy{
      .limiter = cosmosim::hydro::HydroSlopeLimiter::kMonotonizedCentral,
      .dt_over_dx_code = update.dt_code / geometry.cell_width_x_comoving,
      .dt_over_cell_width_code = {
          update.dt_code / geometry.cell_width_x_comoving,
          update.dt_code / geometry.cell_width_y_comoving,
          update.dt_code / geometry.cell_width_z_comoving},
      .rho_floor = 1.0e-10,
      .pressure_floor = 1.0e-10,
      .enable_muscl_hancock_predictor = true});
  cosmosim::hydro::HllcRiemannSolver riemann_solver;

  const double initial_mass = [&]() {
    double sum = 0.0;
    for (double rho : conserved.massDensityComoving()) {
      sum += rho;
    }
    return sum;
  }();

  for (std::size_t step = 0; step < k_step_count; ++step) {
    solver.advancePatch(conserved, geometry, update, reconstruction, riemann_solver, {}, source_context, nullptr);
  }

  double final_mass = 0.0;
  for (double rho : conserved.massDensityComoving()) {
    final_mass += rho;
  }

  assert(std::abs(final_mass - initial_mass) < k_tol);
  assert(conserved.massDensityComoving()[geometry.linearCellIndex(1U, 1U, 1U)] > 0.2);
  assert(conserved.massDensityComoving()[geometry.linearCellIndex(2U, 1U, 1U)] < 0.9);
}

void testContactAdvectionStability() {
  constexpr double k_gamma = 1.4;
  constexpr std::size_t k_nx = 4;
  constexpr std::size_t k_ny = 4;
  constexpr std::size_t k_nz = 4;
  constexpr std::size_t k_cells = k_nx * k_ny * k_nz;
  const auto geometry = makeCartesianGeometry(k_nx, k_ny, k_nz);
  cosmosim::hydro::HydroConservedStateSoa conserved(k_cells);

  for (std::size_t i = 0; i < k_cells; ++i) {
    const auto ijk = geometry.cellIjk(i);
    cosmosim::hydro::HydroPrimitiveState primitive;
    primitive.rho_comoving = (ijk[0] < geometry.nx / 2U) ? 2.0 : 1.0;
    primitive.vel_x_peculiar = 0.2;
    primitive.pressure_comoving = 1.0;
    conserved.storeCell(i, cosmosim::hydro::HydroCoreSolver::conservedFromPrimitive(primitive, k_gamma));
  }

  cosmosim::hydro::HydroUpdateContext update{.dt_code = 2.0e-3, .scale_factor = 1.0, .hubble_rate_code = 0.0};
  cosmosim::hydro::HydroSourceContext source_context{.update = update};

  cosmosim::hydro::HydroCoreSolver solver(k_gamma);
  cosmosim::hydro::MusclHancockReconstruction reconstruction(cosmosim::hydro::HydroReconstructionPolicy{
      .limiter = cosmosim::hydro::HydroSlopeLimiter::kVanLeer,
      .dt_over_dx_code = update.dt_code / geometry.cell_width_x_comoving,
      .dt_over_cell_width_code = {
          update.dt_code / geometry.cell_width_x_comoving,
          update.dt_code / geometry.cell_width_y_comoving,
          update.dt_code / geometry.cell_width_z_comoving},
      .rho_floor = 1.0e-10,
      .pressure_floor = 1.0e-10,
      .enable_muscl_hancock_predictor = true});
  cosmosim::hydro::HllcRiemannSolver riemann_solver;

  for (std::size_t step = 0; step < 20; ++step) {
    solver.advancePatch(conserved, geometry, update, reconstruction, riemann_solver, {}, source_context, nullptr);
  }

  for (double rho : conserved.massDensityComoving()) {
    assert(rho > 0.0);
  }
}

}  // namespace

int main() {
  testSodLikePeriodicConservationAndRegression();
  testContactAdvectionStability();
  return 0;
}
