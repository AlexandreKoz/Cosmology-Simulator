#include <array>
#include <cassert>

#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/physics/black_hole_agn.hpp"

int main() {
  cosmosim::core::SimulationState state;
  state.resizeCells(1);
  state.cells.center_x_comoving[0] = 0.0;
  state.cells.center_y_comoving[0] = 0.0;
  state.cells.center_z_comoving[0] = 0.0;
  state.cells.time_bin[0] = 0;
  state.gas_cells.density_code[0] = 50.0;
  state.gas_cells.sound_speed_code[0] = 10.0;
  state.gas_cells.internal_energy_code[0] = 1.0;

  cosmosim::physics::BlackHoleAgnConfig config;
  config.enabled = true;
  config.seed_halo_mass_threshold_code = 1.0;
  config.seed_mass_code = 10.0;
  config.feedback_coupling_efficiency = 0.7;
  cosmosim::physics::BlackHoleAgnModel model(config);

  std::array<cosmosim::physics::BlackHoleSeedCandidate, 1> seeds{{{0, 10.0, 0}}};
  auto report = model.apply(state, seeds, 1.0, 0);
  assert(report.counters.seeded_bh == 1);

  double previous_mass = state.black_holes.subgrid_mass_code[0];
  double previous_feedback = state.black_holes.cumulative_feedback_energy_code[0];

  for (std::uint64_t step = 1; step <= 4; ++step) {
    report = model.apply(state, std::array<cosmosim::physics::BlackHoleSeedCandidate, 0>{}, 0.5, step);
    assert(report.counters.active_bh == 1);
    assert(state.black_holes.subgrid_mass_code[0] >= previous_mass);
    assert(state.black_holes.cumulative_feedback_energy_code[0] >= previous_feedback);
    previous_mass = state.black_holes.subgrid_mass_code[0];
    previous_feedback = state.black_holes.cumulative_feedback_energy_code[0];
  }

  assert(state.black_holes.duty_cycle_total_time_code[0] >= 2.0);
  assert(state.gas_cells.internal_energy_code[0] > 1.0);
  assert(state.validateOwnershipInvariants());
  return 0;
}
