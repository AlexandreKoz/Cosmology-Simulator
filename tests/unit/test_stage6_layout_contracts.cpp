#include <cassert>
#include <cstdint>
#include <type_traits>

#include "cosmosim/core/simulation_state.hpp"

namespace {

template <typename T>
concept HasCanonicalOwners = requires(T s) {
  s.particles.position_x_comoving;
  s.particles.position_y_comoving;
  s.particles.position_z_comoving;
  s.particles.velocity_x_peculiar;
  s.particles.velocity_y_peculiar;
  s.particles.velocity_z_peculiar;
  s.particles.mass_code;
  s.particles.time_bin;
  s.particle_sidecar.particle_id;
  s.particle_sidecar.species_tag;
  s.particle_sidecar.owning_rank;
  s.particle_sidecar.gravity_softening_comoving;
  s.particle_sidecar.has_gravity_softening_override;
};

template <typename T>
concept HasPersistentAccelerationLane = requires(T p) {
  p.acceleration_x_comoving;
};

}  // namespace

int main() {
  static_assert(HasCanonicalOwners<cosmosim::core::SimulationState>);
  static_assert(!HasPersistentAccelerationLane<cosmosim::core::ParticleSoa>);
  static_assert(std::is_same_v<decltype(cosmosim::core::SimulationState{}.particles), cosmosim::core::ParticleSoa>);

  cosmosim::core::SimulationState state;
  state.resizeParticles(8);
  for (std::size_t i = 0; i < state.particles.size(); ++i) {
    state.particle_sidecar.particle_id[i] = 9000 + i;
    state.particle_sidecar.species_tag[i] = static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kDarkMatter);
    state.particle_sidecar.owning_rank[i] = 0;
    state.particles.mass_code[i] = 1.0 + static_cast<double>(i);
  }

  const auto before_sidecar = state.particle_sidecar.particle_id;
  cosmosim::core::TransientStepWorkspace workspace;
  const std::uint32_t active_idx[] = {1, 3, 5};
  auto gravity_view = cosmosim::core::buildGravityParticleKernelView(state, active_idx, workspace);
  gravity_view.mass_code[2] = 42.0;
  cosmosim::core::scatterGravityParticleKernelView(gravity_view, state);

  assert(state.particles.mass_code[5] == 42.0);
  assert(state.particle_sidecar.particle_id == before_sidecar);
  return 0;
}
