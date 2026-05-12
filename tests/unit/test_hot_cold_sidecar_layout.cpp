#include <array>
#include <cassert>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/core/time_integration.hpp"

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

  bool mirror_reorder_threw = false;
  try {
    (void)cosmosim::core::buildParticleReorderMap(
        state,
        cosmosim::core::ParticleReorderMode::kByTimeBin);
  } catch (const std::invalid_argument&) {
    mirror_reorder_threw = true;
  }
  assert(mirror_reorder_threw);

  cosmosim::core::HierarchicalTimeBinScheduler scheduler(3);
  scheduler.reset(static_cast<std::uint32_t>(state.particles.size()), 0, 0);
  for (std::uint32_t i = 0; i < state.particles.size(); ++i) {
    scheduler.setElementBin(i, time_bins[i], 0);
  }
  cosmosim::core::syncTimeBinMirrorsFromScheduler(scheduler, state);
  const auto reorder_map = cosmosim::core::buildParticleReorderMapByScheduler(state, scheduler);
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

  const std::size_t gravity_capacity_before_clear = workspace.particle_position_x_comoving.capacity();
  const std::size_t gravity_index_capacity_before_clear = workspace.gravity_particle_index.capacity();
  workspace.clear();
  assert(workspace.particle_position_x_comoving.capacity() >= gravity_capacity_before_clear);
  assert(workspace.gravity_particle_index.capacity() >= gravity_index_capacity_before_clear);
  const std::array<std::uint32_t, 2> repeated_particles{1, 3};
  auto repeated_gravity_view = cosmosim::core::buildGravityParticleKernelView(state, repeated_particles, workspace);
  assert(repeated_gravity_view.source_particle_index_generation == state.particleIndexGeneration());
  assert(workspace.particle_position_x_comoving.capacity() >= gravity_capacity_before_clear);

  bool stale_particle_threw = false;
  auto stale_after_resize_view = repeated_gravity_view;
  state.resizeParticles(state.particles.size() + 1U);
  try {
    cosmosim::core::scatterGravityParticleKernelView(stale_after_resize_view, state);
  } catch (const std::runtime_error&) {
    stale_particle_threw = true;
  }
  assert(stale_particle_threw);
  state.resizeParticles(6);

  cosmosim::core::SimulationState hydro_state;
  hydro_state.resizeParticles(3);
  hydro_state.resizeCells(3);
  hydro_state.species.count_by_species = {0, 3, 0, 0, 0};
  for (std::size_t i = 0; i < 3; ++i) {
    hydro_state.particle_sidecar.particle_id[i] = 900 + i;
    hydro_state.particle_sidecar.species_tag[i] =
        static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kGas);
    hydro_state.particle_sidecar.owning_rank[i] = 0;
    hydro_state.cells.center_x_comoving[i] = static_cast<double>(i);
    hydro_state.cells.center_y_comoving[i] = static_cast<double>(i);
    hydro_state.cells.center_z_comoving[i] = static_cast<double>(i);
    hydro_state.cells.mass_code[i] = 1.0;
    hydro_state.cells.patch_index[i] = 0;
    hydro_state.gas_cells.density_code[i] = 10.0 + static_cast<double>(i);
    hydro_state.gas_cells.pressure_code[i] = 20.0 + static_cast<double>(i);
  }
  hydro_state.rebuildSpeciesIndex();
  assert(hydro_state.validateOwnershipInvariants());

  const std::array<std::uint32_t, 2> active_cells{0, 1};
  auto hydro_view = cosmosim::core::buildHydroCellKernelView(hydro_state, active_cells, workspace);
  hydro_view.density_code[0] = 99.0;
  cosmosim::core::scatterHydroCellKernelView(hydro_view, hydro_state);
  assert(hydro_state.gas_cells.density_code[0] == 99.0);

  const std::size_t hydro_capacity_before_clear = workspace.hydro_cell_density_code.capacity();
  workspace.clear();
  assert(workspace.hydro_cell_density_code.capacity() >= hydro_capacity_before_clear);
  auto stale_hydro_view = cosmosim::core::buildHydroCellKernelView(hydro_state, active_cells, workspace);
  hydro_state.resizeCells(hydro_state.cells.size() + 1U);
  bool stale_cell_threw = false;
  try {
    cosmosim::core::scatterHydroCellKernelView(stale_hydro_view, hydro_state);
  } catch (const std::runtime_error&) {
    stale_cell_threw = true;
  }
  assert(stale_cell_threw);

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
