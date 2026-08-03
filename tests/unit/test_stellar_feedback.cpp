#include <cassert>
#include <cmath>
#include <vector>

#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/physics/stellar_feedback.hpp"

namespace {

void seedSingleStarCase(cosmosim::core::SimulationState& state, std::size_t cell_count = 4U) {
  state.resizeCells(cell_count);
  for (std::size_t i = 0; i < state.cells.size(); ++i) {
    state.cells.center_x_comoving[i] = static_cast<double>(i);
    state.cells.center_y_comoving[i] = 0.0;
    state.cells.center_z_comoving[i] = 0.0;
    state.cells.mass_code[i] = 10.0;
    state.gas_cells.gas_cell_id[i] = 100U + i;
    state.gas_cells.density_code[i] = 5.0;  // volume = 2
    state.gas_cells.internal_energy_code[i] = 2.0;
    state.gas_cells.metal_mass_code[i] = 0.1;
  }

  state.resizeParticles(1);
  state.star_particles.resize(1);
  state.particles.position_x_comoving[0] = 0.1;
  state.particles.position_y_comoving[0] = 0.0;
  state.particles.position_z_comoving[0] = 0.0;
  state.star_particles.particle_index[0] = 0;
  state.star_particles.birth_mass_code[0] = 1.0;
}

cosmosim::physics::StellarFeedbackConfig baseConfig() {
  cosmosim::physics::StellarFeedbackConfig config;
  config.variant = cosmosim::physics::StellarFeedbackVariant::kNone;
  config.mode = cosmosim::physics::StellarFeedbackMode::kThermal;
  config.epsilon_thermal = 1.0;
  config.epsilon_kinetic = 0.0;
  config.epsilon_momentum = 0.0;
  config.neighbor_count = 1;
  return config;
}

void testMassDensityEnergyAndMetalConservation() {
  cosmosim::core::SimulationState state;
  seedSingleStarCase(state);
  auto config = baseConfig();
  config.total_energy_code_per_erg = 2.0;
  cosmosim::physics::StellarFeedbackModel model(config);
  cosmosim::physics::StellarFeedbackModuleState module_state;

  const std::vector<std::uint32_t> active = {0};
  const std::vector<double> returned_mass = {0.5};
  const std::vector<double> returned_metals = {0.05};
  const std::vector<double> energy = {4.0};
  const auto report = model.apply(
      state, module_state, active, returned_mass, returned_metals, 1.0, energy);

  assert(report.counters.feedback_stars == 1);
  assert(std::abs(report.counters.deposited_mass_code - 0.5) < 1.0e-12);
  assert(std::abs(report.counters.deposited_metals_code - 0.05) < 1.0e-12);
  assert(std::abs(state.cells.mass_code[0] - 10.5) < 1.0e-12);
  assert(std::abs(state.gas_cells.density_code[0] - 5.25) < 1.0e-12);
  assert(std::abs(state.gas_cells.metal_mass_code[0] - 0.15) < 1.0e-12);
  const double expected_u = (10.0 * 2.0 + 4.0 * 2.0) / 10.5;
  assert(std::abs(state.gas_cells.internal_energy_code[0] - expected_u) < 1.0e-12);
  assert(state.star_particles.enrichment_carry_mass_code[0] == 0.0);
  assert(state.star_particles.stellar_deposited_mass_cumulative_code[0] == 0.5);
}

void testNoNeighborBudgetIsDurablyCarried() {
  cosmosim::core::SimulationState state;
  seedSingleStarCase(state, 0U);
  cosmosim::physics::StellarFeedbackModel model(baseConfig());
  cosmosim::physics::StellarFeedbackModuleState module_state;
  const std::vector<std::uint32_t> active = {0};
  const std::vector<double> returned_mass = {0.25};
  const std::vector<double> returned_metals = {0.02};
  const std::vector<double> energy = {3.0};
  const auto report = model.apply(
      state, module_state, active, returned_mass, returned_metals, 1.0, energy);
  assert(report.counters.deposited_mass_code == 0.0);
  assert(std::abs(state.star_particles.enrichment_carry_mass_code[0] - 0.25) < 1.0e-14);
  assert(std::abs(state.star_particles.enrichment_carry_metals_code[0] - 0.02) < 1.0e-14);
  assert(std::abs(state.star_particles.enrichment_carry_feedback_energy_erg[0] - 3.0) < 1.0e-14);
}

void testModeSelectionMomentumOnly() {
  auto config = baseConfig();
  config.mode = cosmosim::physics::StellarFeedbackMode::kMomentum;
  config.epsilon_thermal = 0.0;
  config.epsilon_momentum = 1.0;
  cosmosim::physics::StellarFeedbackModel model(config);
  const auto budget = model.computeBudget(1.0, 0.1, 0.01);
  assert(budget.thermal_energy_erg == 0.0);
  assert(budget.kinetic_energy_erg == 0.0);
  assert(budget.momentum_budget_code > 0.0);
}

}  // namespace

int main() {
  testMassDensityEnergyAndMetalConservation();
  testNoNeighborBudgetIsDurablyCarried();
  testModeSelectionMomentumOnly();
  return 0;
}
