#include <array>
#include <cassert>
#include <cmath>
#include <string>

#include "cosmosim/core/time_integration.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/physics/star_formation.hpp"

namespace {

void testThresholdEligibility() {
  cosmosim::physics::StarFormationConfig config;
  config.density_threshold_code = 5.0;
  config.temperature_threshold_k = 1.0e4;
  config.min_converging_flow_rate_code = -0.2;
  cosmosim::physics::StarFormationModel model(config);

  cosmosim::physics::StarFormationCellInput cell;
  cell.gas_mass_code = 1.0;
  cell.gas_density_code = 5.2;
  cell.gas_temperature_k = 9.0e3;
  cell.velocity_divergence_code = -0.3;
  assert(model.isEligible(cell));

  cell.gas_temperature_k = 2.0e4;
  assert(!model.isEligible(cell));
}

void testSchmidtKennicuttRate() {
  cosmosim::physics::StarFormationConfig config;
  config.epsilon_ff = 0.02;
  cosmosim::physics::StarFormationModel model(config);

  const double rho = 12.0;
  const double t_ff = model.freeFallTimeCode(rho);
  const double expected_rate = config.epsilon_ff * rho / t_ff;
  assert(std::abs(model.sfrDensityRateCode(rho) - expected_rate) / expected_rate < 1.0e-12);
}

void testConservationAndMetadata() {
  cosmosim::core::SimulationState state;
  state.resizeCells(1);
  state.cells.center_x_comoving[0] = 1.0;
  state.cells.center_y_comoving[0] = 2.0;
  state.cells.center_z_comoving[0] = 3.0;
  state.cells.mass_code[0] = 10.0;
  state.gas_cells.density_code[0] = 50.0;
  state.gas_cells.temperature_code[0] = 1.0e3;

  cosmosim::physics::StarFormationConfig config;
  config.stochastic_spawning = false;
  config.epsilon_ff = 0.05;
  cosmosim::physics::StarFormationModel model(config);

  const std::array<std::uint32_t, 1> active{0};
  const std::array<double, 1> div_v{-1.0};
  const std::array<double, 1> metallicity{0.02};

  const auto report = model.apply(state, active, div_v, metallicity, 1.0e9, 0.5, 7, 0);
  assert(report.counters.spawn_events == 1);
  assert(report.counters.spawned_mass_code > 0.0);

  const double gas_mass_after = state.cells.mass_code[0];
  const double star_mass_after = state.particles.mass_code[0];
  assert(std::abs((gas_mass_after + star_mass_after) - 10.0) < 1.0e-10);
  assert(state.star_particles.size() == 1);
  assert(state.star_particles.formation_scale_factor[0] == 0.5);

  const auto* sidecar = state.sidecars.find("star_formation");
  assert(sidecar != nullptr);
  const std::string payload(reinterpret_cast<const char*>(sidecar->payload.data()), sidecar->payload.size());
  assert(payload.find("spawn_events=1") != std::string::npos);
}

void testTimeIntegrationCallbackHook() {
  cosmosim::core::SimulationState state;
  state.resizeCells(1);
  state.cells.mass_code[0] = 2.0;
  state.gas_cells.density_code[0] = 20.0;
  state.gas_cells.temperature_code[0] = 5.0e3;

  cosmosim::physics::StarFormationConfig config;
  config.stochastic_spawning = false;
  config.epsilon_ff = 0.1;
  cosmosim::physics::StarFormationCallback callback{cosmosim::physics::StarFormationModel(config)};
  callback.setVelocityDivergenceCode(std::array<double, 1>{-0.5});
  callback.setMetallicityMassFraction(std::array<double, 1>{0.02});

  cosmosim::core::IntegratorState integrator_state;
  integrator_state.dt_time_code = 1.0e8;
  const std::array<std::uint32_t, 1> active_cells{0};
  cosmosim::core::ActiveSetDescriptor active_set{
      .cell_indices = active_cells,
      .cells_are_subset = true,
  };
  cosmosim::core::StepContext context{
      .state = state,
      .integrator_state = integrator_state,
      .active_set = active_set,
      .stage = cosmosim::core::IntegrationStage::kSourceTerms,
  };

  callback.onStage(context);
  assert(callback.lastStepReport().counters.spawn_events == 1);
  assert(state.particles.size() == 1);
}

}  // namespace

int main() {
  testThresholdEligibility();
  testSchmidtKennicuttRate();
  testConservationAndMetadata();
  testTimeIntegrationCallbackHook();
  return 0;
}
