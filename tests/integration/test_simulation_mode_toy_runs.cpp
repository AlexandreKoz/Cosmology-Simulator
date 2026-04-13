#include <array>
#include <stdexcept>
#include <string>
#include <vector>

#include "cosmosim/core/config.hpp"
#include "cosmosim/core/simulation_mode.hpp"

namespace {

void requireOrThrow(bool condition, const std::string& message) {
  if (!condition) {
    throw std::runtime_error(message);
  }
}

void testPeriodicToyRunPolicyAndGhosts() {
  cosmosim::core::SimulationConfig config;
  config.mode.mode = cosmosim::core::SimulationMode::kCosmoCube;
  config.mode.hydro_boundary = cosmosim::core::ModeHydroBoundary::kPeriodic;
  config.mode.gravity_boundary = cosmosim::core::ModeGravityBoundary::kPeriodic;

  const auto policy = cosmosim::core::buildModePolicy(config.mode);
  cosmosim::core::validateModePolicy(config, policy);

  requireOrThrow(policy.hydro_boundary == cosmosim::core::BoundaryCondition::kPeriodic,
      "cosmo_cube must map to periodic hydro boundary");
  requireOrThrow(policy.gravity_boundary == cosmosim::core::GravityBoundaryModel::kPeriodicPoisson,
      "cosmo_cube must map to periodic gravity boundary");

  const std::array<double, 4> potential = {11.0, 12.0, 13.0, 14.0};
  std::vector<double> left(2, 0.0);
  std::vector<double> right(2, 0.0);
  cosmosim::core::fillGravityPotentialGhostCells1d(
      potential,
      left,
      right,
      policy.gravity_boundary,
      -1.0);

  requireOrThrow(left[0] == 13.0 && left[1] == 14.0,
      "periodic gravity left ghosts must wrap from right interior edge");
  requireOrThrow(right[0] == 11.0 && right[1] == 12.0,
      "periodic gravity right ghosts must wrap from left interior edge");
}

void testIsolatedToyRunPolicyAndGhosts() {
  cosmosim::core::SimulationConfig config;
  config.mode.mode = cosmosim::core::SimulationMode::kIsolatedGalaxy;
  config.mode.hydro_boundary = cosmosim::core::ModeHydroBoundary::kOpen;
  config.mode.gravity_boundary = cosmosim::core::ModeGravityBoundary::kIsolatedMonopole;

  const auto policy = cosmosim::core::buildModePolicy(config.mode);
  cosmosim::core::validateModePolicy(config, policy);

  requireOrThrow(policy.hydro_boundary == cosmosim::core::BoundaryCondition::kOpen,
      "isolated_galaxy must map to open hydro boundary");
  requireOrThrow(policy.gravity_boundary == cosmosim::core::GravityBoundaryModel::kIsolatedMonopoleDirichlet,
      "isolated_galaxy must map to isolated gravity boundary");

  const std::array<double, 3> density = {1.0, 2.0, 3.0};
  const std::array<double, 3> velocity = {4.0, 5.0, 6.0};
  const std::array<double, 3> pressure = {7.0, 8.0, 9.0};

  std::vector<double> left_density(1, 0.0);
  std::vector<double> right_density(1, 0.0);
  std::vector<double> left_velocity(1, 0.0);
  std::vector<double> right_velocity(1, 0.0);
  std::vector<double> left_pressure(1, 0.0);
  std::vector<double> right_pressure(1, 0.0);

  cosmosim::core::fillHydroGhostCells1d(
      {.density_code = density, .velocity_normal_code = velocity, .pressure_code = pressure},
      {.left_density_code = left_density,
       .right_density_code = right_density,
       .left_velocity_normal_code = left_velocity,
       .right_velocity_normal_code = right_velocity,
       .left_pressure_code = left_pressure,
       .right_pressure_code = right_pressure},
      policy.hydro_boundary);

  requireOrThrow(left_density[0] == 1.0 && right_density[0] == 3.0,
      "open hydro ghosts must clamp to edge cells");
  requireOrThrow(left_velocity[0] == 4.0 && right_pressure[0] == 9.0,
      "open hydro ghosts must preserve edge primitive values");

  std::vector<double> left_potential(1, 0.0);
  std::vector<double> right_potential(1, 0.0);
  cosmosim::core::fillGravityPotentialGhostCells1d(
      density,
      left_potential,
      right_potential,
      policy.gravity_boundary,
      -42.0);
  requireOrThrow(left_potential[0] == -42.0 && right_potential[0] == -42.0,
      "isolated gravity ghosts must use fixed reference potential");
}

void testModeValidationRejectsBoundaryMismatch() {
  cosmosim::core::SimulationConfig config;
  config.mode.mode = cosmosim::core::SimulationMode::kCosmoCube;
  config.mode.hydro_boundary = cosmosim::core::ModeHydroBoundary::kPeriodic;
  config.mode.gravity_boundary = cosmosim::core::ModeGravityBoundary::kIsolatedMonopole;

  const auto policy = cosmosim::core::buildModePolicy(config.mode);
  bool threw = false;
  try {
    cosmosim::core::validateModePolicy(config, policy);
  } catch (const cosmosim::core::ConfigError&) {
    threw = true;
  }
  requireOrThrow(threw,
      "validateModePolicy must reject cosmological mode with isolated gravity boundary override");
}

}  // namespace

int main() {
  testPeriodicToyRunPolicyAndGhosts();
  testIsolatedToyRunPolicyAndGhosts();
  testModeValidationRejectsBoundaryMismatch();
  return 0;
}
