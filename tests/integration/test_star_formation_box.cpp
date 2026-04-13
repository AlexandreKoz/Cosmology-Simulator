#include <cassert>
#include <vector>

#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/physics/star_formation.hpp"

int main() {
  cosmosim::core::SimulationState state;
  state.resizeCells(8);
  for (std::size_t i = 0; i < state.cells.size(); ++i) {
    state.cells.center_x_comoving[i] = static_cast<double>(i);
    state.cells.mass_code[i] = 2.0;
    state.gas_cells.density_code[i] = 20.0;
    state.gas_cells.temperature_code[i] = 8.0e3;
  }

  cosmosim::physics::StarFormationConfig config;
  config.epsilon_ff = 0.02;
  config.stochastic_spawning = false;
  cosmosim::physics::StarFormationModel model(config);

  std::vector<std::uint32_t> active(state.cells.size());
  for (std::size_t i = 0; i < active.size(); ++i) {
    active[i] = static_cast<std::uint32_t>(i);
  }
  std::vector<double> div_v(state.cells.size(), -0.5);
  std::vector<double> metallicity(state.cells.size(), 0.01);

  const double gas_mass_before = 16.0;
  auto report = model.apply(state, active, div_v, metallicity, 3.0e8, 0.9, 42, 0);
  assert(report.counters.eligible_cells == state.cells.size());
  assert(report.counters.spawned_particles == state.cells.size());

  double gas_mass_after = 0.0;
  for (const double m : state.cells.mass_code) {
    gas_mass_after += m;
  }

  double star_mass = 0.0;
  for (const double m : state.particles.mass_code) {
    star_mass += m;
  }

  assert(star_mass > 0.0);
  assert(gas_mass_after < gas_mass_before);
  assert(gas_mass_after + star_mass <= gas_mass_before + 1.0e-9);
  assert(state.validateOwnershipInvariants());
  return 0;
}
