#include <cassert>
#include <cmath>
#include <vector>

#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/physics/stellar_feedback.hpp"

namespace {

void seedSingleStarCase(cosmosim::core::SimulationState& state) {
  state.resizeCells(4);
  for (std::size_t i = 0; i < state.cells.size(); ++i) {
    state.cells.center_x_comoving[i] = static_cast<double>(i);
    state.cells.center_y_comoving[i] = 0.0;
    state.cells.center_z_comoving[i] = 0.0;
    state.cells.mass_code[i] = 10.0;
    state.gas_cells.density_code[i] = 10.0;
    state.gas_cells.internal_energy_code[i] = 0.0;
  }

  state.resizeParticles(1);
  state.star_particles.resize(1);
  state.particles.position_x_comoving[0] = 0.1;
  state.particles.position_y_comoving[0] = 0.0;
  state.particles.position_z_comoving[0] = 0.0;
  state.star_particles.particle_index[0] = 0;
  state.star_particles.birth_mass_code[0] = 1.0;
}

void testBudgetConservationNoVariant() {
  cosmosim::core::SimulationState state;
  seedSingleStarCase(state);

  cosmosim::physics::StellarFeedbackConfig config;
  config.variant = cosmosim::physics::StellarFeedbackVariant::kNone;
  config.mode = cosmosim::physics::StellarFeedbackMode::kThermalKineticMomentum;
  config.epsilon_thermal = 1.0;
  config.epsilon_kinetic = 0.0;
  config.epsilon_momentum = 0.0;
  config.neighbor_count = 2;

  cosmosim::physics::StellarFeedbackModel model(config);
  cosmosim::physics::StellarFeedbackModuleState module_state;

  const std::vector<std::uint32_t> active = {0};
  const std::vector<double> returned_mass = {0.5};
  const std::vector<double> returned_metals = {0.05};

  const auto report = model.apply(state, module_state, active, returned_mass, returned_metals, 1.0e-3);
  assert(report.counters.feedback_stars == 1);
  assert(std::abs(report.counters.deposited_mass_code - 0.5) < 1.0e-12);
  assert(std::abs(report.counters.deposited_metals_code - 0.05) < 1.0e-12);
  assert(report.counters.unresolved_mass_code < 1.0e-14);
  assert(report.counters.unresolved_thermal_energy_erg < 1.0e-14);
}

void testModeSelectionMomentumOnly() {
  cosmosim::physics::StellarFeedbackConfig config;
  config.mode = cosmosim::physics::StellarFeedbackMode::kMomentum;
  cosmosim::physics::StellarFeedbackModel model(config);

  const auto budget = model.computeBudget(1.0, 0.1, 0.01);
  assert(budget.thermal_energy_erg == 0.0);
  assert(budget.kinetic_energy_erg == 0.0);
  assert(budget.momentum_budget_code > 0.0);
}

}  // namespace

int main() {
  testBudgetConservationNoVariant();
  testModeSelectionMomentumOnly();
  return 0;
}
