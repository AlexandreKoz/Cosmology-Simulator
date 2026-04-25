#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#include "cosmosim/core/simulation_state.hpp"

namespace {

cosmosim::core::SimulationState makeGasContractState() {
  cosmosim::core::SimulationState state;
  state.resizeParticles(4);
  state.resizeCells(3);
  state.species.count_by_species = {1, 3, 0, 0, 0};

  state.particle_sidecar.particle_id = {100, 200, 101, 102};
  state.particle_sidecar.species_tag = {
      static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kGas),
      static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kDarkMatter),
      static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kGas),
      static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kGas)};
  state.particle_sidecar.sfc_key = {30, 0, 20, 10};

  for (std::size_t cell = 0; cell < state.cells.size(); ++cell) {
    state.cells.center_x_comoving[cell] = 0.1 * static_cast<double>(cell + 1);
    state.cells.mass_code[cell] = 10.0 + static_cast<double>(cell);
    state.gas_cells.density_code[cell] = 100.0 + static_cast<double>(cell);
    state.gas_cells.pressure_code[cell] = 200.0 + static_cast<double>(cell);
    state.gas_cells.internal_energy_code[cell] = 300.0 + static_cast<double>(cell);
    state.gas_cells.recon_gradient_x[cell] = 0.5 + static_cast<double>(cell);
  }

  state.rebuildSpeciesIndex();
  return state;
}

std::unordered_map<std::uint64_t, double> densityByGasParticleId(const cosmosim::core::SimulationState& state) {
  std::unordered_map<std::uint64_t, double> by_id;
  const auto gas_globals = state.particle_species_index.globalIndices(cosmosim::core::ParticleSpecies::kGas);
  assert(gas_globals.size() == state.cells.size());
  for (std::size_t cell = 0; cell < state.cells.size(); ++cell) {
    by_id.emplace(state.particle_sidecar.particle_id[gas_globals[cell]], state.gas_cells.density_code[cell]);
  }
  return by_id;
}

void test_gas_cell_identity_invariants() {
  cosmosim::core::SimulationState state = makeGasContractState();
  cosmosim::core::debugAssertGasCellIdentityContract(state);

  const auto baseline_density = densityByGasParticleId(state);
  assert(baseline_density.at(100) == 100.0);
  assert(baseline_density.at(101) == 101.0);
  assert(baseline_density.at(102) == 102.0);

  cosmosim::core::TransientStepWorkspace workspace;
  const std::array<std::uint32_t, 2> active_cells{0, 2};
  auto hydro_view = cosmosim::core::buildHydroCellKernelView(state, active_cells, workspace);
  hydro_view.density_code[0] = 555.0;
  hydro_view.pressure_code[1] = 777.0;
  cosmosim::core::scatterHydroCellKernelView(hydro_view, state);

  assert(state.gas_cells.density_code[0] == 555.0);
  assert(state.gas_cells.pressure_code[2] == 777.0);
  // Reconstruction scratch lanes are persistent sidecar fields and are not part of HydroCellKernelView.
  assert(state.gas_cells.recon_gradient_x[0] == 0.5);
  assert(state.gas_cells.internal_energy_code[2] == 302.0);
}

void test_gas_cell_reorder_resize_invariants() {
  cosmosim::core::SimulationState state = makeGasContractState();
  const auto baseline_density = densityByGasParticleId(state);

  // Allowed reorder: gas relative order unchanged.
  const auto by_species = cosmosim::core::buildParticleReorderMap(state, cosmosim::core::ParticleReorderMode::kBySpecies);
  cosmosim::core::reorderParticles(state, by_species);
  assert(densityByGasParticleId(state) == baseline_density);

  // Forbidden reorder: gas relative order changes without gas-cell rebuild.
  const auto by_sfc = cosmosim::core::buildParticleReorderMap(state, cosmosim::core::ParticleReorderMode::kBySfcKey);
  bool reorder_threw = false;
  try {
    cosmosim::core::reorderParticles(state, by_sfc);
  } catch (const std::invalid_argument&) {
    reorder_threw = true;
  }
  assert(reorder_threw);

  cosmosim::core::TransientStepWorkspace workspace;
  const std::array<std::uint32_t, 1> active_cell{0};
  auto view = cosmosim::core::buildHydroCellKernelView(state, active_cell, workspace);
  state.resizeCells(2);

  bool stale_scatter_threw = false;
  try {
    cosmosim::core::scatterHydroCellKernelView(view, state);
  } catch (const std::runtime_error&) {
    stale_scatter_threw = true;
  }
  assert(stale_scatter_threw);

  bool contract_threw = false;
  try {
    cosmosim::core::debugAssertGasCellIdentityContract(state);
  } catch (const std::runtime_error&) {
    contract_threw = true;
  }
  assert(contract_threw);

  // Active extraction remains available as a generic cell view; contract checks are explicit.
  (void)cosmosim::core::buildHydroCellKernelView(state, active_cell, workspace);
}

}  // namespace

int main() {
  test_gas_cell_identity_invariants();
  test_gas_cell_reorder_resize_invariants();
  return 0;
}
