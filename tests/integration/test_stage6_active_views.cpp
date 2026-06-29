#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <optional>

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

  cosmosim::core::SimulationState gas_state;
  gas_state.resizeCells(4);
  for (std::size_t row = 0; row < gas_state.cells.size(); ++row) {
    gas_state.cells.center_x_comoving[row] = static_cast<double>(row);
    gas_state.cells.mass_code[row] = 10.0 + static_cast<double>(row);
    gas_state.gas_cells.density_code[row] = 20.0 + static_cast<double>(row);
    gas_state.gas_cells.pressure_code[row] = 30.0 + static_cast<double>(row);
  }
  gas_state.replaceGasCellIdentityRecords({
      {.gas_cell_id = 401, .parent_particle_id = std::nullopt, .owning_patch_id = 0, .local_cell_row = 0},
      {.gas_cell_id = 404, .parent_particle_id = std::nullopt, .owning_patch_id = 0, .local_cell_row = 1},
      {.gas_cell_id = 402, .parent_particle_id = std::nullopt, .owning_patch_id = 0, .local_cell_row = 2},
      {.gas_cell_id = 403, .parent_particle_id = std::nullopt, .owning_patch_id = 0, .local_cell_row = 3},
  });

  cosmosim::core::TransientStepWorkspace hydro_workspace;
  std::array<std::uint32_t, 2> active_cells{1, 3};
  auto hydro_compact = cosmosim::core::buildHydroCellKernelView(gas_state, active_cells, hydro_workspace);
  assert(hydro_compact.gas_cell_id[0] == 404);
  assert(hydro_compact.gas_cell_id[1] == 403);
  hydro_compact.cell_index[0] = 0;
  hydro_compact.cell_index[1] = 0;
  hydro_compact.density_code[0] = 88.0;
  hydro_compact.pressure_code[1] = 99.0;
  cosmosim::core::scatterHydroCellKernelView(hydro_compact, gas_state);
  assert(gas_state.gas_cells.density_code[0] == 20.0);
  assert(gas_state.gas_cells.density_code[1] == 88.0);
  assert(gas_state.gas_cells.pressure_code[3] == 99.0);

  return 0;
}
