#pragma once

#include <cstddef>
#include <span>
#include <string>

#include "cosmosim/core/config.hpp"

namespace cosmosim::core {

enum class BoundaryCondition {
  kPeriodic,
  kOpen,
  kReflective,
};

enum class GravityBoundaryModel {
  kPeriodicPoisson,
  kIsolatedMonopoleDirichlet,
};

struct ModePolicy {
  SimulationMode simulation_mode = SimulationMode::kZoomIn;
  BoundaryCondition hydro_boundary = BoundaryCondition::kPeriodic;
  GravityBoundaryModel gravity_boundary = GravityBoundaryModel::kPeriodicPoisson;
  bool cosmological_comoving_frame = true;
  bool zoom_region_expected = false;
};

struct Hydro1dStateView {
  std::span<const double> density_code;
  std::span<const double> velocity_normal_code;
  std::span<const double> pressure_code;
};

struct Hydro1dGhostBufferView {
  std::span<double> left_density_code;
  std::span<double> right_density_code;
  std::span<double> left_velocity_normal_code;
  std::span<double> right_velocity_normal_code;
  std::span<double> left_pressure_code;
  std::span<double> right_pressure_code;
};

void fillHydroGhostCells1d(
    const Hydro1dStateView& interior,
    const Hydro1dGhostBufferView& ghosts,
    BoundaryCondition boundary_condition);

void fillGravityPotentialGhostCells1d(
    std::span<const double> interior_potential,
    std::span<double> left_ghost_potential,
    std::span<double> right_ghost_potential,
    GravityBoundaryModel boundary_model,
    double isolated_reference_potential = 0.0);

[[nodiscard]] ModePolicy buildModePolicy(const ModeConfig& mode_config);
void validateModePolicy(const SimulationConfig& config, const ModePolicy& policy);

[[nodiscard]] std::string boundaryConditionToString(BoundaryCondition boundary_condition);
[[nodiscard]] std::string gravityBoundaryModelToString(GravityBoundaryModel boundary_model);

}  // namespace cosmosim::core
