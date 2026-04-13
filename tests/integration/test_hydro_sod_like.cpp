#include <cassert>
#include <cmath>
#include <vector>

#include "cosmosim/hydro/hydro_core_solver.hpp"

namespace {

constexpr double k_tol = 1.0e-9;

cosmosim::hydro::HydroPatchGeometry makePeriodic1dGeometry(std::size_t cell_count) {
  cosmosim::hydro::HydroPatchGeometry geometry;
  geometry.cell_volume_comoving = 1.0;
  geometry.faces.reserve(cell_count);

  for (std::size_t i = 0; i < cell_count; ++i) {
    geometry.faces.push_back(cosmosim::hydro::HydroFace{
        .owner_cell = i,
        .neighbor_cell = (i + 1U) % cell_count,
        .area_comoving = 1.0,
        .normal_x = 1.0,
        .normal_y = 0.0,
        .normal_z = 0.0});
  }

  return geometry;
}

void fillSodLikeInitialState(cosmosim::hydro::HydroConservedStateSoa& conserved, double gamma) {
  for (std::size_t i = 0; i < conserved.size(); ++i) {
    cosmosim::hydro::HydroPrimitiveState primitive;
    if (i < conserved.size() / 2U) {
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
  constexpr std::size_t k_cell_count = 64;
  constexpr std::size_t k_step_count = 30;
  constexpr double k_gamma = 1.4;

  cosmosim::hydro::HydroConservedStateSoa conserved(k_cell_count);
  fillSodLikeInitialState(conserved, k_gamma);

  const cosmosim::hydro::HydroPatchGeometry geometry = makePeriodic1dGeometry(k_cell_count);

  cosmosim::hydro::HydroUpdateContext update;
  update.dt_code = 1.0e-3;
  update.scale_factor = 1.0;

  cosmosim::hydro::HydroSourceContext source_context;
  source_context.update = update;

  cosmosim::hydro::HydroCoreSolver solver(k_gamma);
  cosmosim::hydro::MusclHancockReconstruction reconstruction(cosmosim::hydro::HydroReconstructionPolicy{
      .limiter = cosmosim::hydro::HydroSlopeLimiter::kMonotonizedCentral,
      .dt_over_dx_code = update.dt_code,
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
  assert(conserved.massDensityComoving()[31U] > 0.2);
  assert(conserved.massDensityComoving()[32U] < 0.9);
}

void testContactAdvectionStability() {
  constexpr double k_gamma = 1.4;
  constexpr std::size_t k_cells = 64;
  cosmosim::hydro::HydroConservedStateSoa conserved(k_cells);

  for (std::size_t i = 0; i < k_cells; ++i) {
    cosmosim::hydro::HydroPrimitiveState primitive;
    primitive.rho_comoving = (i < k_cells / 2U) ? 2.0 : 1.0;
    primitive.vel_x_peculiar = 0.2;
    primitive.pressure_comoving = 1.0;
    conserved.storeCell(i, cosmosim::hydro::HydroCoreSolver::conservedFromPrimitive(primitive, k_gamma));
  }

  const auto geometry = makePeriodic1dGeometry(k_cells);
  cosmosim::hydro::HydroUpdateContext update{.dt_code = 2.0e-3, .scale_factor = 1.0, .hubble_rate_code = 0.0};
  cosmosim::hydro::HydroSourceContext source_context{.update = update};

  cosmosim::hydro::HydroCoreSolver solver(k_gamma);
  cosmosim::hydro::MusclHancockReconstruction reconstruction(cosmosim::hydro::HydroReconstructionPolicy{
      .limiter = cosmosim::hydro::HydroSlopeLimiter::kVanLeer,
      .dt_over_dx_code = update.dt_code,
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
