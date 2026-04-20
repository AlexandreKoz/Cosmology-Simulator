#include <cassert>
#include <cmath>
#include <cstdint>
#include <utility>
#include <vector>

#include "cosmosim/core/time_integration.hpp"
#include "cosmosim/physics/star_formation.hpp"

namespace {

class GravityKickMock final : public cosmosim::core::IntegrationCallback {
 public:
  std::string_view callbackName() const override { return "gravity_kick_mock"; }

  void onStage(cosmosim::core::StepContext& context) override {
    if (context.stage != cosmosim::core::IntegrationStage::kGravityKickPre &&
        context.stage != cosmosim::core::IntegrationStage::kGravityKickPost) {
      return;
    }

    auto& state = context.state;
    const double kick = 0.25 * context.integrator_state.dt_time_code;
    for (std::size_t i = 0; i < state.particles.size(); ++i) {
      state.particles.velocity_x_peculiar[i] += kick;
    }
  }
};

class ActiveSubsetKickMock final : public cosmosim::core::IntegrationCallback {
 public:
  std::string_view callbackName() const override { return "active_subset_kick_mock"; }

  void onStage(cosmosim::core::StepContext& context) override {
    if (context.stage != cosmosim::core::IntegrationStage::kGravityKickPre &&
        context.stage != cosmosim::core::IntegrationStage::kGravityKickPost) {
      return;
    }
    ++kick_stage_invocations;
    kicked_particles_per_stage.push_back(context.active_set.particle_indices.size());

    const double kick = 0.5 * context.integrator_state.dt_time_code;
    for (const std::uint32_t particle_index : context.active_set.particle_indices) {
      context.state.particles.velocity_x_peculiar[particle_index] += kick;
    }
  }

  int kick_stage_invocations = 0;
  std::vector<std::size_t> kicked_particles_per_stage;
};

void runNoPhysicsLoop() {
  cosmosim::core::SimulationState state;
  cosmosim::core::IntegratorState integrator_state;
  integrator_state.dt_time_code = 0.1;
  integrator_state.current_scale_factor = 1.0;

  cosmosim::core::StepOrchestrator orchestrator;
  for (int step = 0; step < 5; ++step) {
    orchestrator.executeSingleStep(state, integrator_state, {}, nullptr, nullptr);
  }

  assert(std::abs(integrator_state.current_time_code - 0.5) < 1.0e-12);
  assert(integrator_state.step_index == 5U);
  assert(std::abs(integrator_state.current_scale_factor - 1.0) < 1.0e-12);
}

void runGravityOnlyLoop() {
  cosmosim::core::SimulationState state;
  state.resizeParticles(4);

  for (std::size_t i = 0; i < state.particles.size(); ++i) {
    state.particles.velocity_x_peculiar[i] = 0.0;
  }

  cosmosim::core::IntegratorState integrator_state;
  integrator_state.dt_time_code = 5.0e16;
  integrator_state.current_scale_factor = 1.0;

  cosmosim::core::CosmologyBackgroundConfig cfg;
  cfg.hubble_param = 0.7;
  cfg.omega_matter = 0.3;
  cfg.omega_lambda = 0.7;
  cosmosim::core::LambdaCdmBackground background(cfg);

  GravityKickMock gravity_kick;
  cosmosim::core::StepOrchestrator orchestrator;
  orchestrator.registerCallback(gravity_kick);

  for (int step = 0; step < 4; ++step) {
    orchestrator.executeSingleStep(state, integrator_state, {}, &background, nullptr);
  }

  // Two kick stages per step with 0.25*dt each => net velocity increment per step is 0.5*dt.
  const double expected_velocity = 4.0 * 0.5 * integrator_state.dt_time_code;
  for (std::size_t i = 0; i < state.particles.size(); ++i) {
    assert(std::abs(state.particles.velocity_x_peculiar[i] - expected_velocity) < 1.0e-12);
  }

  assert(integrator_state.current_time_code > 0.0);
  assert(integrator_state.current_scale_factor > 1.0);
}

void runStarFormationSourceTermLoop() {
  cosmosim::core::SimulationState state;
  state.resizeCells(4);
  for (std::size_t i = 0; i < state.cells.size(); ++i) {
    state.cells.center_x_comoving[i] = static_cast<double>(i);
    state.cells.mass_code[i] = 4.0;
    state.gas_cells.density_code[i] = 30.0;
    state.gas_cells.temperature_code[i] = 5.0e3;
  }

  cosmosim::core::IntegratorState integrator_state;
  integrator_state.dt_time_code = 5.0e8;
  integrator_state.current_scale_factor = 1.0;

  cosmosim::physics::StarFormationConfig sf_config;
  sf_config.stochastic_spawning = false;
  sf_config.epsilon_ff = 0.03;
  cosmosim::physics::StarFormationModel sf_model(sf_config);
  cosmosim::physics::StarFormationCallback sf_callback(std::move(sf_model));

  std::vector<std::uint32_t> active_cells(state.cells.size());
  std::vector<double> velocity_divergence(state.cells.size(), -0.5);
  std::vector<double> metallicity(state.cells.size(), 0.01);
  for (std::size_t i = 0; i < active_cells.size(); ++i) {
    active_cells[i] = static_cast<std::uint32_t>(i);
  }
  sf_callback.setVelocityDivergenceCode(velocity_divergence);
  sf_callback.setMetallicityMassFraction(metallicity);

  cosmosim::core::StepOrchestrator orchestrator;
  orchestrator.registerCallback(sf_callback);

  double gas_mass_before = 0.0;
  for (const double mass : state.cells.mass_code) {
    gas_mass_before += mass;
  }

  cosmosim::core::ActiveSetDescriptor active_set{
      .cell_indices = active_cells,
      .cells_are_subset = true,
  };
  for (int step = 0; step < 3; ++step) {
    orchestrator.executeSingleStep(state, integrator_state, active_set, nullptr, nullptr);
  }

  double gas_mass_after = 0.0;
  for (const double mass : state.cells.mass_code) {
    gas_mass_after += mass;
  }

  double star_mass = 0.0;
  for (const double mass : state.particles.mass_code) {
    star_mass += mass;
  }

  assert(sf_callback.lastStepReport().counters.scanned_cells == state.cells.size());
  assert(gas_mass_after < gas_mass_before);
  assert(star_mass > 0.0);
  assert(std::abs((gas_mass_after + star_mass) - gas_mass_before) < 1.0e-9);
  assert(integrator_state.step_index == 3U);
}

void runActiveSubsetKickContractLoop() {
  cosmosim::core::SimulationState state;
  state.resizeParticles(6);
  for (std::size_t i = 0; i < state.particles.size(); ++i) {
    state.particles.velocity_x_peculiar[i] = 0.0;
  }

  cosmosim::core::IntegratorState integrator_state;
  integrator_state.dt_time_code = 0.2;
  integrator_state.current_scale_factor = 1.0;

  std::vector<std::uint32_t> active_particles = {0, 2, 5};
  cosmosim::core::ActiveSetDescriptor active_set{
      .particle_indices = active_particles,
      .particles_are_subset = true,
  };

  ActiveSubsetKickMock subset_kick;
  cosmosim::core::StepOrchestrator orchestrator;
  orchestrator.registerCallback(subset_kick);
  orchestrator.executeSingleStep(state, integrator_state, active_set, nullptr, nullptr);

  assert(subset_kick.kick_stage_invocations == 2);
  assert(subset_kick.kicked_particles_per_stage.size() == 2);
  assert(subset_kick.kicked_particles_per_stage[0] == active_particles.size());
  assert(subset_kick.kicked_particles_per_stage[1] == active_particles.size());

  // Two kick stages, each applies 0.5*dt to active particles only.
  const double expected_active_velocity = integrator_state.dt_time_code;
  for (std::size_t i = 0; i < state.particles.size(); ++i) {
    const bool is_active = (i == 0 || i == 2 || i == 5);
    if (is_active) {
      assert(std::abs(state.particles.velocity_x_peculiar[i] - expected_active_velocity) < 1.0e-12);
    } else {
      assert(std::abs(state.particles.velocity_x_peculiar[i]) < 1.0e-12);
    }
  }
}

}  // namespace

int main() {
  runNoPhysicsLoop();
  runGravityOnlyLoop();
  runStarFormationSourceTermLoop();
  runActiveSubsetKickContractLoop();
  return 0;
}
