#include <cassert>
#include <cmath>

#include "cosmosim/hydro/hydro_core_solver.hpp"
#include "cosmosim/hydro/hydro_riemann.hpp"

namespace {

constexpr double k_tol = 1.0e-10;

void testHlleSymmetricPressureFlux() {
  const cosmosim::hydro::HydroPrimitiveState left{.rho_comoving = 1.0, .pressure_comoving = 2.0};
  const cosmosim::hydro::HydroPrimitiveState right = left;
  const cosmosim::hydro::HydroFace face{.owner_cell = 0, .neighbor_cell = 1, .area_comoving = 1.0, .normal_x = 1.0};

  cosmosim::hydro::HlleRiemannSolver solver;
  const auto flux = solver.computeFlux(left, right, face, 1.4);
  assert(std::abs(flux.mass_density_comoving) < k_tol);
  assert(std::abs(flux.momentum_density_x_comoving - 2.0) < 1.0e-8);
}

void testHllcSymmetricContactRegression() {
  const cosmosim::hydro::HydroPrimitiveState left{.rho_comoving = 1.0, .vel_x_peculiar = 0.0, .pressure_comoving = 1.0};
  const cosmosim::hydro::HydroPrimitiveState right = left;
  const cosmosim::hydro::HydroFace face{.owner_cell = 0, .neighbor_cell = 1, .area_comoving = 1.0, .normal_x = 1.0};

  cosmosim::hydro::HllcRiemannSolver solver;
  const auto flux = solver.computeFlux(left, right, face, 1.4);
  assert(std::abs(flux.mass_density_comoving) < k_tol);
  assert(std::abs(flux.momentum_density_x_comoving - 1.0) < 1.0e-8);
  assert(solver.fallbackCount() == 0U);
}

}  // namespace

int main() {
  testHlleSymmetricPressureFlux();
  testHllcSymmetricContactRegression();
  return 0;
}
