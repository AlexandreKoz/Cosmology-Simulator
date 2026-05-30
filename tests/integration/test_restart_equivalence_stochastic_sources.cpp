#include <cassert>
#include <cstdint>
#include <filesystem>
#include <utility>
#include <vector>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/physics/star_formation.hpp"
#include "restart_equivalence_harness.hpp"
#include "restart_equivalence_scenarios.hpp"

namespace {

void applyProductionStochasticStarFormation(cosmosim::tests::RestartEquivalenceStepContext& context) {
  cosmosim::physics::StarFormationConfig config;
  config.enabled = true;
  config.density_threshold_code = 10.0;
  config.temperature_threshold_k = 1.0e4;
  config.min_converging_flow_rate_code = 0.0;
  config.epsilon_ff = 1.0;
  config.min_star_particle_mass_code = 0.2;
  config.stochastic_spawning = true;
  config.random_seed = 123456789ull;
  config.newton_g_code = 1.0;
  cosmosim::physics::StarFormationModel model(config);

  std::vector<std::uint32_t> active_cells;
  active_cells.reserve(context.active_particle_indices.size());
  for (const std::uint32_t element_index : context.active_particle_indices) {
    if (element_index < context.state.cells.size()) {
      active_cells.push_back(element_index);
    }
  }
  std::vector<double> velocity_divergence(context.state.cells.size(), -1.0);
  std::vector<double> metallicity(context.state.cells.size(), 0.02);
  const std::size_t particle_count_before = context.state.particles.size();
  const auto report = model.apply(
      context.state,
      active_cells,
      velocity_divergence,
      metallicity,
      context.dt_code,
      context.integrator_state.current_scale_factor,
      context.integrator_state.step_index + 1U,
      0U);
  const std::size_t particle_count_after = context.state.particles.size();
  if (particle_count_after > particle_count_before) {
    context.scheduler.appendElements(
        static_cast<std::uint32_t>(particle_count_after - particle_count_before),
        0,
        context.scheduler.currentTick() + 1U);
  }
  (void)report;
}

}  // namespace

int main() {
#if COSMOSIM_ENABLE_HDF5
  const auto restart_path = cosmosim::tests::stage8RestartPath("restart_equivalence_stochastic_sources");
  auto state = cosmosim::tests::makeStage8HydroToyState(8, "restart_equivalence_stochastic_sources");
  for (std::size_t cidx = 0; cidx < state.cells.size(); ++cidx) {
    state.cells.mass_code[cidx] = 2.0;
    state.particles.mass_code[cidx] = 2.0;
    state.gas_cells.density_code[cidx] = 1000.0 + static_cast<double>(cidx);
    state.gas_cells.pressure_code[cidx] = 1.0;
    state.gas_cells.temperature_code[cidx] = 100.0;
    state.gas_cells.internal_energy_code[cidx] = 1.5;
    state.gas_cells.sound_speed_code[cidx] = 0.1;
  }
  auto scheduler = cosmosim::tests::makeStage8Scheduler(static_cast<std::uint32_t>(state.particles.size()), 1);
  auto integrator_state = cosmosim::tests::makeStage8IntegratorState(1, 1);
  auto output_state = cosmosim::tests::makeStage8OutputCadenceState(false);
  auto scenario = cosmosim::tests::makeStage8Scenario(
      std::move(state), integrator_state, std::move(scheduler), std::move(output_state), restart_path, 6, 3);
  scenario.step_kernel = applyProductionStochasticStarFormation;
  scenario.tolerances.position_abs = 1.0e-12;
  scenario.tolerances.velocity_abs = 1.0e-12;
  scenario.tolerances.scalar_abs = 1.0e-12;
  const auto result = cosmosim::tests::runRestartEquivalenceScenario(std::move(scenario));
  assert(result.direct_state.star_particles.size() > 0U);
  assert(result.direct_state.star_particles.size() == result.restarted_state.star_particles.size());
  std::filesystem::remove(restart_path);
#endif
  return 0;
}
