#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>

#include "cosmosim/cosmosim.hpp"

int main() {
  const cosmosim::gravity::PmGridShape shape{16, 16, 16};
  cosmosim::gravity::PmGridStorage grid_a1(shape);
  cosmosim::gravity::PmGridStorage grid_a05(shape);
  cosmosim::gravity::PmSolver solver(shape);

  std::vector<double> x{0.25, 0.75};
  std::vector<double> y{0.5, 0.5};
  std::vector<double> z{0.5, 0.5};
  std::vector<double> m{1.0, -1.0};

  cosmosim::gravity::PmSolveOptions opts;
  opts.box_size_x_mpc_comoving = 1.0;
  opts.box_size_y_mpc_comoving = 1.0;
  opts.box_size_z_mpc_comoving = 1.0;
  opts.assignment_scheme = cosmosim::gravity::PmAssignmentScheme::kCic;
  opts.boundary_condition = cosmosim::gravity::PmBoundaryCondition::kPeriodic;

  auto solve_and_norm = [&](cosmosim::gravity::PmGridStorage& grid, double scale_factor) {
    opts.scale_factor = scale_factor;
    grid.clear();
    solver.assignDensity(grid, x, y, z, m, opts);
    solver.solvePoissonPeriodic(grid, opts);
    double accum = 0.0;
    for (double v : grid.force_x()) accum += v * v;
    for (double v : grid.force_y()) accum += v * v;
    for (double v : grid.force_z()) accum += v * v;
    return std::sqrt(accum);
  };

  const double norm_a1 = solve_and_norm(grid_a1, 1.0);
  const double norm_a05 = solve_and_norm(grid_a05, 0.5);
  assert(norm_a1 > 0.0);
  const double ratio = norm_a05 / norm_a1;
  assert(std::abs(ratio - 0.25) < 1.0e-6);
  return 0;
}
