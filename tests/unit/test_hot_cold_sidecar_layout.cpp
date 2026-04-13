#include <array>
#include <cassert>
#include <cstdint>
#include <stdexcept>

#include "cosmosim/core/simulation_state.hpp"

int main() {
  cosmosim::core::SimulationState state;
  state.resizeParticles(6);
  state.species.count_by_species = {2, 0, 2, 1, 1};

  const std::array<std::uint32_t, 6> species{0, 2, 4, 2, 0, 3};
  const std::array<std::uint8_t, 6> time_bins{2, 0, 1, 0, 2, 1};
  const std::array<std::uint64_t, 6> sfc{30, 10, 60, 20, 40, 50};
  for (std::size_t i = 0; i < 6; ++i) {
    state.particle_sidecar.particle_id[i] = 100 + i;
    state.particle_sidecar.sfc_key[i] = sfc[i];
    state.particle_sidecar.species_tag[i] = species[i];
    state.particle_sidecar.owning_rank[i] = 0;
    state.particles.time_bin[i] = time_bins[i];
    state.particles.position_x_comoving[i] = static_cast<double>(i);
    state.particles.position_y_comoving[i] = static_cast<double>(i) + 10.0;
    state.particles.position_z_comoving[i] = static_cast<double>(i) + 20.0;
    state.particles.velocity_x_peculiar[i] = 0.1;
    state.particles.velocity_y_peculiar[i] = 0.2;
    state.particles.velocity_z_peculiar[i] = 0.3;
    state.particles.mass_code[i] = 1.0;
  }

  state.star_particles.resize(2);
  state.star_particles.particle_index = {1, 3};
  state.black_holes.resize(1);
  state.black_holes.particle_index[0] = 5;
  state.tracers.resize(1);
  state.tracers.particle_index[0] = 2;

  state.rebuildSpeciesIndex();

  const auto reorder_map = cosmosim::core::buildParticleReorderMap(
      state,
      cosmosim::core::ParticleReorderMode::kByTimeBin);
  assert(reorder_map.isConsistent(state.particles.size()));

  const auto old_star_first_id = state.particle_sidecar.particle_id[state.star_particles.particle_index[0]];
  cosmosim::core::reorderParticles(state, reorder_map);

  cosmosim::core::debugAssertNoStaleParticleIndices(state);
  assert(state.validateOwnershipInvariants());

  const auto new_star_first_id = state.particle_sidecar.particle_id[state.star_particles.particle_index[0]];
  assert(new_star_first_id == old_star_first_id);

  for (std::size_t i = 1; i < state.particles.size(); ++i) {
    assert(state.particles.time_bin[i - 1] <= state.particles.time_bin[i]);
  }

  cosmosim::core::TransientStepWorkspace workspace;
  const std::array<std::uint32_t, 3> active_particles{0, 2, 4};
  auto gravity_view = cosmosim::core::buildGravityParticleKernelView(state, active_particles, workspace);
  assert(gravity_view.size() == 3);
  gravity_view.velocity_x_peculiar[1] += 2.0;
  cosmosim::core::scatterGravityParticleKernelView(gravity_view, state);
  assert(state.particles.velocity_x_peculiar[active_particles[1]] == gravity_view.velocity_x_peculiar[1]);

  state.resizeCells(3);
  for (std::size_t i = 0; i < 3; ++i) {
    state.cells.center_x_comoving[i] = static_cast<double>(i);
    state.cells.center_y_comoving[i] = static_cast<double>(i);
    state.cells.center_z_comoving[i] = static_cast<double>(i);
    state.cells.mass_code[i] = 1.0;
    state.cells.patch_index[i] = 0;
    state.gas_cells.density_code[i] = 10.0 + static_cast<double>(i);
    state.gas_cells.pressure_code[i] = 20.0 + static_cast<double>(i);
  }

  const std::array<std::uint32_t, 2> active_cells{0, 1};
  auto hydro_view = cosmosim::core::buildHydroCellKernelView(state, active_cells, workspace);
  hydro_view.density_code[0] = 99.0;
  cosmosim::core::scatterHydroCellKernelView(hydro_view, state);
  assert(state.gas_cells.density_code[0] == 99.0);

  bool threw = false;
  try {
    state.star_particles.particle_index[0] = 99;
    cosmosim::core::debugAssertNoStaleParticleIndices(state);
  } catch (const std::runtime_error&) {
    threw = true;
  }
  assert(threw);

  return 0;
}
