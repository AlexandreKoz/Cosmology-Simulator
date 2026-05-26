#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>

#include "cosmosim/core/simulation_state.hpp"

int main() {
  cosmosim::core::SimulationState state;
  state.resizeParticles(16);

  for (std::size_t i = 0; i < state.particles.size(); ++i) {
    state.particle_sidecar.particle_id[i] = 10000 + i;
    state.particle_sidecar.species_tag[i] = static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kDarkMatter);
    state.particle_sidecar.owning_rank[i] = 0;
    state.particles.position_x_comoving[i] = static_cast<double>(i) * 0.1;
    state.particles.position_y_comoving[i] = static_cast<double>(i) * 0.2;
    state.particles.position_z_comoving[i] = static_cast<double>(i) * 0.3;
    state.particles.mass_code[i] = 1.0 + 0.01 * static_cast<double>(i);
  }

  cosmosim::core::TransientStepWorkspace workspace;
  std::array<std::uint32_t, 6> active{0, 2, 4, 7, 10, 12};
  auto compact = cosmosim::core::buildGravityParticleKernelView(state, active, workspace);

  double subset_sum = 0.0;
  for (std::size_t i = 0; i < compact.size(); ++i) {
    subset_sum += compact.position_x_comoving[i] + compact.position_y_comoving[i] + compact.position_z_comoving[i] + compact.mass_code[i];
    compact.mass_code[i] += 1.0;
  }

  double full_sum_on_subset = 0.0;
  for (const std::uint32_t idx : active) {
    full_sum_on_subset += state.particles.position_x_comoving[idx] + state.particles.position_y_comoving[idx] +
                          state.particles.position_z_comoving[idx] + state.particles.mass_code[idx];
  }

  assert(std::abs(subset_sum - full_sum_on_subset) < 1.0e-14);

  cosmosim::core::scatterGravityParticleKernelView(compact, state);
  for (const std::uint32_t idx : active) {
    assert(std::abs(state.particles.mass_code[idx] - (2.0 + 0.01 * static_cast<double>(idx))) < 1.0e-14);
  }

  return 0;
}
