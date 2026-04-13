#include <array>
#include <cassert>
#include <vector>

#include "cosmosim/core/config.hpp"
#include "cosmosim/core/simulation_mode.hpp"

namespace {

void testModePolicyDefaultsAndValidation() {
  cosmosim::core::SimulationConfig config;
  config.mode.mode = cosmosim::core::SimulationMode::kCosmoCube;
  config.mode.hydro_boundary = cosmosim::core::ModeHydroBoundary::kAuto;
  config.mode.gravity_boundary = cosmosim::core::ModeGravityBoundary::kAuto;

  const auto policy = cosmosim::core::buildModePolicy(config.mode);
  assert(policy.hydro_boundary == cosmosim::core::BoundaryCondition::kPeriodic);
  assert(policy.gravity_boundary == cosmosim::core::GravityBoundaryModel::kPeriodicPoisson);
  cosmosim::core::validateModePolicy(config, policy);
}

void testInvalidModeOverrideFailsValidation() {
  cosmosim::core::SimulationConfig config;
  config.mode.mode = cosmosim::core::SimulationMode::kZoomIn;
  config.mode.hydro_boundary = cosmosim::core::ModeHydroBoundary::kOpen;

  bool threw = false;
  try {
    const auto policy = cosmosim::core::buildModePolicy(config.mode);
    cosmosim::core::validateModePolicy(config, policy);
  } catch (const cosmosim::core::ConfigError&) {
    threw = true;
  }
  assert(threw);
}

void testHydroBoundaryFill() {
  const std::array<double, 4> density = {1.0, 2.0, 3.0, 4.0};
  const std::array<double, 4> velocity = {10.0, 20.0, 30.0, 40.0};
  const std::array<double, 4> pressure = {100.0, 200.0, 300.0, 400.0};

  std::vector<double> left_density(2, 0.0);
  std::vector<double> right_density(2, 0.0);
  std::vector<double> left_velocity(2, 0.0);
  std::vector<double> right_velocity(2, 0.0);
  std::vector<double> left_pressure(2, 0.0);
  std::vector<double> right_pressure(2, 0.0);

  cosmosim::core::fillHydroGhostCells1d(
      {.density_code = density, .velocity_normal_code = velocity, .pressure_code = pressure},
      {.left_density_code = left_density,
       .right_density_code = right_density,
       .left_velocity_normal_code = left_velocity,
       .right_velocity_normal_code = right_velocity,
       .left_pressure_code = left_pressure,
       .right_pressure_code = right_pressure},
      cosmosim::core::BoundaryCondition::kReflective);

  assert(left_density[0] == 4.0);
  assert(left_velocity[0] == -40.0);
  assert(right_density[0] == 4.0);
  assert(right_velocity[1] == -30.0);
}

void testGravityBoundaryFillPeriodicAndIsolated() {
  const std::array<double, 4> potential = {5.0, 6.0, 7.0, 8.0};
  std::vector<double> left(2, 0.0);
  std::vector<double> right(2, 0.0);

  cosmosim::core::fillGravityPotentialGhostCells1d(
      potential,
      left,
      right,
      cosmosim::core::GravityBoundaryModel::kPeriodicPoisson,
      -99.0);
  assert(left[0] == 7.0);
  assert(right[0] == 5.0);

  cosmosim::core::fillGravityPotentialGhostCells1d(
      potential,
      left,
      right,
      cosmosim::core::GravityBoundaryModel::kIsolatedMonopoleDirichlet,
      -99.0);
  assert(left[0] == -99.0);
  assert(right[1] == -99.0);
}

}  // namespace

int main() {
  testModePolicyDefaultsAndValidation();
  testInvalidModeOverrideFailsValidation();
  testHydroBoundaryFill();
  testGravityBoundaryFillPeriodicAndIsolated();
  return 0;
}
