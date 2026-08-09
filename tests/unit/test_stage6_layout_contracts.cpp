#include <cassert>
#include <cstdint>
#include <cstddef>
#include <stdexcept>
#include <type_traits>
#include <utility>

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

template <typename T>
concept HasPersistentGasReconstructionGradientLane = requires(T gas) {
  gas.recon_gradient_x;
  gas.recon_gradient_y;
  gas.recon_gradient_z;
};

}  // namespace

int main() {
  static_assert(HasCanonicalOwners<cosmosim::core::SimulationState>);
  static_assert(cosmosim::core::k_is_canonical_particle_state_owner_v<cosmosim::core::SimulationState>);
  static_assert(!cosmosim::core::k_is_canonical_particle_state_owner_v<cosmosim::core::ParticleSoaStorage>);
  static_assert(!cosmosim::core::ParticleSoaStorage::k_owns_persistent_particle_truth);
  static_assert(!cosmosim::core::ParticleSoaStorage::k_is_restart_serializable);
  static_assert(!HasPersistentAccelerationLane<cosmosim::core::ParticleSoa>);
  static_assert(!HasPersistentGasReconstructionGradientLane<cosmosim::core::GasCellSidecar>);
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

  auto direct_gravity = cosmosim::core::buildGravityParticleKernelViewAllParticlesDirect(state, workspace);
  assert(direct_gravity.position_x_comoving.data() == state.particles.position_x_comoving.data());
  assert(direct_gravity.velocity_x_peculiar.data() == state.particles.velocity_x_peculiar.data());
  assert(direct_gravity.mass_code.data() == state.particles.mass_code.data());
  direct_gravity.mass_code[0] = 77.0;
  assert(state.particles.mass_code[0] == 77.0);

  const std::uint32_t duplicate_idx[] = {2, 2};
  auto duplicate_view = cosmosim::core::buildGravityParticleKernelView(state, duplicate_idx, workspace);
  duplicate_view.mass_code[0] = 123.0;
  duplicate_view.mass_code[1] = 456.0;
  const double mass_before_duplicate_scatter = state.particles.mass_code[2];
  bool duplicate_scatter_threw = false;
  try {
    cosmosim::core::scatterGravityParticleKernelView(duplicate_view, state);
  } catch (const std::invalid_argument&) {
    duplicate_scatter_threw = true;
  }
  assert(duplicate_scatter_threw);
  assert(state.particles.mass_code[2] == mass_before_duplicate_scatter);

  // Regression: growth of the monotonic scratch arena must never invalidate
  // allocations returned earlier in the same scratch epoch.
  cosmosim::core::MonotonicScratchAllocator scratch(8);
  auto* early = scratch.allocateArray<std::uint64_t>(1);
  assert(early != nullptr);
  *early = UINT64_C(0x123456789abcdef0);
  auto* aligned = scratch.allocateBytes(257, 64);
  assert(aligned != nullptr);
  assert(reinterpret_cast<std::uintptr_t>(aligned) % 64U == 0U);
  (void)scratch.allocateBytes(4096, alignof(std::max_align_t));
  assert(*early == UINT64_C(0x123456789abcdef0));
  const auto capacity_before_reset = scratch.capacityBytes();
  scratch.reset();
  assert(scratch.capacityBytes() == capacity_before_reset);
  auto* reused = scratch.allocateArray<std::uint64_t>(1);
  assert(reused != nullptr);

  bool bad_alignment_threw = false;
  try {
    (void)scratch.allocateBytes(16, 3);
  } catch (const std::invalid_argument&) {
    bad_alignment_threw = true;
  }
  assert(bad_alignment_threw);
  return 0;
}
