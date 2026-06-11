#include <array>
#include <cassert>
#include <cmath>
#include <utility>

#include "cosmosim/hydro/hydro_cartesian_patch.hpp"
#include "cosmosim/hydro/hydro_core_solver.hpp"

namespace {

constexpr std::size_t k_n = 4;
constexpr std::size_t k_steps = 12;
constexpr double k_gamma = 1.4;
constexpr double k_tol = 2.0e-12;

enum class Axis {
  kX,
  kY,
  kZ,
};

cosmosim::hydro::HydroPatchGeometry makeGeometry() {
  return cosmosim::hydro::makeCartesianPatchGeometry(cosmosim::hydro::HydroCartesianPatchSpec{
      .nx = k_n,
      .ny = k_n,
      .nz = k_n,
      .origin_x_comoving = 0.0,
      .origin_y_comoving = 0.0,
      .origin_z_comoving = 0.0,
      .cell_width_x_comoving = 1.0 / static_cast<double>(k_n),
      .cell_width_y_comoving = 1.0 / static_cast<double>(k_n),
      .cell_width_z_comoving = 1.0 / static_cast<double>(k_n)});
}

std::size_t coordinateForAxis(const std::array<std::size_t, 3>& ijk, Axis axis) {
  switch (axis) {
    case Axis::kX:
      return ijk[0];
    case Axis::kY:
      return ijk[1];
    case Axis::kZ:
      return ijk[2];
  }
  return ijk[0];
}

cosmosim::hydro::HydroConservedStateSoa runAxis(Axis axis) {
  const cosmosim::hydro::HydroPatchGeometry geometry = makeGeometry();
  cosmosim::hydro::HydroConservedStateSoa conserved(geometry.cellCount());
  for (std::size_t row = 0; row < conserved.size(); ++row) {
    const std::array<std::size_t, 3> ijk = geometry.cellIjk(row);
    const bool left_state = coordinateForAxis(ijk, axis) < k_n / 2U;
    cosmosim::hydro::HydroPrimitiveState primitive;
    primitive.rho_comoving = left_state ? 1.0 : 0.125;
    primitive.pressure_comoving = left_state ? 1.0 : 0.1;
    conserved.storeCell(row, cosmosim::hydro::HydroCoreSolver::conservedFromPrimitive(primitive, k_gamma));
  }

  const cosmosim::hydro::HydroUpdateContext update{.dt_code = 7.5e-4, .scale_factor = 1.0, .hubble_rate_code = 0.0};
  const cosmosim::hydro::HydroSourceContext source_context{.update = update};
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
      .enable_muscl_hancock_predictor = true,
      .adiabatic_index = k_gamma});
  cosmosim::hydro::HllcRiemannSolver riemann_solver;

  for (std::size_t step = 0; step < k_steps; ++step) {
    solver.advancePatch(conserved, geometry, update, reconstruction, riemann_solver, {}, source_context, nullptr);
  }
  return conserved;
}

cosmosim::hydro::HydroPrimitiveState loadPermuted(
    const cosmosim::hydro::HydroConservedStateSoa& conserved,
    const cosmosim::hydro::HydroPatchGeometry& geometry,
    Axis axis,
    std::size_t i,
    std::size_t j,
    std::size_t k) {
  std::size_t row = geometry.linearCellIndex(i, j, k);
  if (axis == Axis::kY) {
    row = geometry.linearCellIndex(j, i, k);
  } else if (axis == Axis::kZ) {
    row = geometry.linearCellIndex(k, j, i);
  }
  cosmosim::hydro::HydroPrimitiveState primitive =
      cosmosim::hydro::HydroCoreSolver::primitiveFromConserved(conserved.loadCell(row), k_gamma);
  if (axis == Axis::kY) {
    std::swap(primitive.vel_x_peculiar, primitive.vel_y_peculiar);
  } else if (axis == Axis::kZ) {
    std::swap(primitive.vel_x_peculiar, primitive.vel_z_peculiar);
  }
  return primitive;
}

void assertEquivalentProfiles(
    const cosmosim::hydro::HydroConservedStateSoa& x_state,
    const cosmosim::hydro::HydroConservedStateSoa& other_state,
    Axis other_axis) {
  const cosmosim::hydro::HydroPatchGeometry geometry = makeGeometry();
  for (std::size_t k = 0; k < k_n; ++k) {
    for (std::size_t j = 0; j < k_n; ++j) {
      for (std::size_t i = 0; i < k_n; ++i) {
        const cosmosim::hydro::HydroPrimitiveState x =
            loadPermuted(x_state, geometry, Axis::kX, i, j, k);
        const cosmosim::hydro::HydroPrimitiveState other =
            loadPermuted(other_state, geometry, other_axis, i, j, k);
        assert(std::abs(x.rho_comoving - other.rho_comoving) < k_tol);
        assert(std::abs(x.vel_x_peculiar - other.vel_x_peculiar) < k_tol);
        assert(std::abs(x.vel_y_peculiar - other.vel_y_peculiar) < k_tol);
        assert(std::abs(x.vel_z_peculiar - other.vel_z_peculiar) < k_tol);
        assert(std::abs(x.pressure_comoving - other.pressure_comoving) < k_tol);
      }
    }
  }
}

void testAxisAlignedShockSymmetry() {
  const cosmosim::hydro::HydroConservedStateSoa x_state = runAxis(Axis::kX);
  const cosmosim::hydro::HydroConservedStateSoa y_state = runAxis(Axis::kY);
  const cosmosim::hydro::HydroConservedStateSoa z_state = runAxis(Axis::kZ);

  assertEquivalentProfiles(x_state, y_state, Axis::kY);
  assertEquivalentProfiles(x_state, z_state, Axis::kZ);
}

}  // namespace

int main() {
  testAxisAlignedShockSymmetry();
  return 0;
}
