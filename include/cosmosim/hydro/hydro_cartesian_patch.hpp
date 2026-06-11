#pragma once

#include <array>
#include <cstddef>

#include "cosmosim/hydro/hydro_core_solver.hpp"

namespace cosmosim::hydro {

struct HydroCartesianPatchSpec {
  std::size_t nx = 0;
  std::size_t ny = 0;
  std::size_t nz = 0;
  double origin_x_comoving = 0.0;
  double origin_y_comoving = 0.0;
  double origin_z_comoving = 0.0;
  double cell_width_x_comoving = 0.0;
  double cell_width_y_comoving = 0.0;
  double cell_width_z_comoving = 0.0;
};

[[nodiscard]] HydroPatchGeometry makeCartesianPatchGeometry(const HydroCartesianPatchSpec& spec);

// Return exact factors for cell_count, ordered nx >= ny >= nz and biased toward a cube.
[[nodiscard]] std::array<std::size_t, 3> chooseNearCubicCartesianFactors(std::size_t cell_count);

}  // namespace cosmosim::hydro
