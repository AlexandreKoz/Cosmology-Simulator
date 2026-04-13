#include <cassert>
#include <cmath>
#include <vector>

#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/physics/stellar_feedback.hpp"

int main() {
  cosmosim::core::SimulationState state;
  state.resizeCells(8);
  for (std::size_t i = 0; i < state.cells.size(); ++i) {
    state.cells.center_x_comoving[i] = static_cast<double>(i % 2);
    state.cells.center_y_comoving[i] = static_cast<double>((i / 2) % 2);
    state.cells.center_z_comoving[i] = static_cast<double>((i / 4) % 2);
    state.cells.mass_code[i] = 5.0;
    state.gas_cells.density_code[i] = 5.0;
    state.gas_cells.internal_energy_code[i] = 1.0;
  }

  state.resizeParticles(2);
  state.star_particles.resize(2);
  for (std::size_t star_index = 0; star_index < 2; ++star_index) {
    state.particles.position_x_comoving[star_index] = 0.5;
    state.particles.position_y_comoving[star_index] = 0.5;
    state.particles.position_z_comoving[star_index] = 0.5;
    state.star_particles.particle_index[star_index] = static_cast<std::uint32_t>(star_index);
    state.star_particles.birth_mass_code[star_index] = 1.0;
  }

  cosmosim::physics::StellarFeedbackConfig config;
  config.mode = cosmosim::physics::StellarFeedbackMode::kThermalKineticMomentum;
  config.variant = cosmosim::physics::StellarFeedbackVariant::kStochastic;
  config.stochastic_event_probability = 1.0;
  config.neighbor_count = 4;
  cosmosim::physics::StellarFeedbackModel model(config);
  cosmosim::physics::StellarFeedbackModuleState module_state;

  const std::vector<std::uint32_t> active = {0, 1};
  const std::vector<double> returned_mass = {0.3, 0.2};
  const std::vector<double> returned_metals = {0.03, 0.02};

  const auto report = model.apply(state, module_state, active, returned_mass, returned_metals, 1.0e-3);
  assert(report.counters.feedback_stars == 2);
  assert(std::abs(report.counters.deposited_mass_code - 0.5) < 1.0e-12);
  assert(std::abs(report.counters.deposited_metals_code - 0.05) < 1.0e-12);
  assert(report.counters.deposited_thermal_energy_erg > 0.0);
  assert(report.counters.deposited_kinetic_energy_erg > 0.0);
  assert(report.counters.unresolved_momentum_code > 0.0);

  return 0;
}
