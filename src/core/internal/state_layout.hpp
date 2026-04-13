#pragma once

#include <cstdint>
#include <vector>

namespace cosmosim::core::internal {

// Hot particle data: structure-of-arrays for solver-facing vectorization.
struct ParticleHotSoa {
  std::vector<double> position_x_comoving;
  std::vector<double> position_y_comoving;
  std::vector<double> position_z_comoving;
  std::vector<double> velocity_x_peculiar;
  std::vector<double> velocity_y_peculiar;
  std::vector<double> velocity_z_peculiar;
  std::vector<double> mass_code;
};

// Cold sidecar data: metadata and provenance that should not pollute hot loops.
struct ParticleColdSidecar {
  std::vector<std::uint64_t> particle_id;
  std::vector<std::uint32_t> species;
  std::vector<std::uint32_t> flags;
};

}  // namespace cosmosim::core::internal
