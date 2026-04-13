#include "cosmosim/core/simulation_mode.hpp"

#include <algorithm>
#include <stdexcept>

namespace cosmosim::core {
namespace {

[[nodiscard]] std::size_t requireSameSize(
    std::span<const double> density,
    std::span<const double> velocity,
    std::span<const double> pressure) {
  if (density.size() != velocity.size() || density.size() != pressure.size()) {
    throw std::invalid_argument("hydro state spans must have matching sizes");
  }
  return density.size();
}

void validateGhostBuffers(const Hydro1dGhostBufferView& ghosts) {
  const std::size_t left_count = ghosts.left_density_code.size();
  const std::size_t right_count = ghosts.right_density_code.size();
  if (left_count == 0 || right_count == 0) {
    throw std::invalid_argument("hydro ghost buffers must provide at least one left/right cell");
  }

  if (ghosts.left_velocity_normal_code.size() != left_count ||
      ghosts.left_pressure_code.size() != left_count ||
      ghosts.right_velocity_normal_code.size() != right_count ||
      ghosts.right_pressure_code.size() != right_count) {
    throw std::invalid_argument("hydro ghost buffers must have consistent component sizes");
  }
}

[[nodiscard]] bool isIsolatedMode(SimulationMode mode) {
  return mode == SimulationMode::kIsolatedGalaxy || mode == SimulationMode::kIsolatedCluster;
}

}  // namespace

void fillHydroGhostCells1d(
    const Hydro1dStateView& interior,
    const Hydro1dGhostBufferView& ghosts,
    BoundaryCondition boundary_condition) {
  const std::size_t interior_count = requireSameSize(
      interior.density_code,
      interior.velocity_normal_code,
      interior.pressure_code);
  if (interior_count == 0) {
    throw std::invalid_argument("hydro interior state must not be empty");
  }
  validateGhostBuffers(ghosts);

  const std::size_t left_ghost_count = ghosts.left_density_code.size();
  const std::size_t right_ghost_count = ghosts.right_density_code.size();

  for (std::size_t i = 0; i < left_ghost_count; ++i) {
    const std::size_t src = std::min(i, interior_count - 1U);
    const std::size_t mirrored_src = std::min(interior_count - 1U - src, interior_count - 1U);

    switch (boundary_condition) {
      case BoundaryCondition::kPeriodic: {
        const std::size_t periodic_src = (interior_count - left_ghost_count + i) % interior_count;
        ghosts.left_density_code[i] = interior.density_code[periodic_src];
        ghosts.left_velocity_normal_code[i] = interior.velocity_normal_code[periodic_src];
        ghosts.left_pressure_code[i] = interior.pressure_code[periodic_src];
        break;
      }
      case BoundaryCondition::kOpen: {
        ghosts.left_density_code[i] = interior.density_code.front();
        ghosts.left_velocity_normal_code[i] = interior.velocity_normal_code.front();
        ghosts.left_pressure_code[i] = interior.pressure_code.front();
        break;
      }
      case BoundaryCondition::kReflective: {
        ghosts.left_density_code[i] = interior.density_code[mirrored_src];
        ghosts.left_velocity_normal_code[i] = -interior.velocity_normal_code[mirrored_src];
        ghosts.left_pressure_code[i] = interior.pressure_code[mirrored_src];
        break;
      }
    }
  }

  for (std::size_t i = 0; i < right_ghost_count; ++i) {
    const std::size_t src = std::min(i, interior_count - 1U);
    const std::size_t mirrored_src = interior_count - 1U - src;

    switch (boundary_condition) {
      case BoundaryCondition::kPeriodic: {
        const std::size_t periodic_src = i % interior_count;
        ghosts.right_density_code[i] = interior.density_code[periodic_src];
        ghosts.right_velocity_normal_code[i] = interior.velocity_normal_code[periodic_src];
        ghosts.right_pressure_code[i] = interior.pressure_code[periodic_src];
        break;
      }
      case BoundaryCondition::kOpen: {
        ghosts.right_density_code[i] = interior.density_code.back();
        ghosts.right_velocity_normal_code[i] = interior.velocity_normal_code.back();
        ghosts.right_pressure_code[i] = interior.pressure_code.back();
        break;
      }
      case BoundaryCondition::kReflective: {
        ghosts.right_density_code[i] = interior.density_code[mirrored_src];
        ghosts.right_velocity_normal_code[i] = -interior.velocity_normal_code[mirrored_src];
        ghosts.right_pressure_code[i] = interior.pressure_code[mirrored_src];
        break;
      }
    }
  }
}

void fillGravityPotentialGhostCells1d(
    std::span<const double> interior_potential,
    std::span<double> left_ghost_potential,
    std::span<double> right_ghost_potential,
    GravityBoundaryModel boundary_model,
    double isolated_reference_potential) {
  if (interior_potential.empty()) {
    throw std::invalid_argument("gravity interior potential must not be empty");
  }
  if (left_ghost_potential.empty() || right_ghost_potential.empty()) {
    throw std::invalid_argument("gravity ghost buffers must provide at least one left/right cell");
  }

  for (std::size_t i = 0; i < left_ghost_potential.size(); ++i) {
    if (boundary_model == GravityBoundaryModel::kPeriodicPoisson) {
      const std::size_t src = (interior_potential.size() - left_ghost_potential.size() + i) %
          interior_potential.size();
      left_ghost_potential[i] = interior_potential[src];
    } else {
      left_ghost_potential[i] = isolated_reference_potential;
    }
  }

  for (std::size_t i = 0; i < right_ghost_potential.size(); ++i) {
    if (boundary_model == GravityBoundaryModel::kPeriodicPoisson) {
      const std::size_t src = i % interior_potential.size();
      right_ghost_potential[i] = interior_potential[src];
    } else {
      right_ghost_potential[i] = isolated_reference_potential;
    }
  }
}

ModePolicy buildModePolicy(const ModeConfig& mode_config) {
  ModePolicy policy;
  policy.simulation_mode = mode_config.mode;
  policy.zoom_region_expected = mode_config.zoom_high_res_region;

  switch (mode_config.mode) {
    case SimulationMode::kCosmoCube:
      policy.hydro_boundary = BoundaryCondition::kPeriodic;
      policy.gravity_boundary = GravityBoundaryModel::kPeriodicPoisson;
      policy.cosmological_comoving_frame = true;
      break;
    case SimulationMode::kZoomIn:
      policy.hydro_boundary = BoundaryCondition::kPeriodic;
      policy.gravity_boundary = GravityBoundaryModel::kPeriodicPoisson;
      policy.cosmological_comoving_frame = true;
      break;
    case SimulationMode::kIsolatedGalaxy:
      policy.hydro_boundary = BoundaryCondition::kOpen;
      policy.gravity_boundary = GravityBoundaryModel::kIsolatedMonopoleDirichlet;
      policy.cosmological_comoving_frame = false;
      break;
    case SimulationMode::kIsolatedCluster:
      policy.hydro_boundary = BoundaryCondition::kReflective;
      policy.gravity_boundary = GravityBoundaryModel::kIsolatedMonopoleDirichlet;
      policy.cosmological_comoving_frame = false;
      break;
  }

  switch (mode_config.hydro_boundary) {
    case ModeHydroBoundary::kPeriodic:
      policy.hydro_boundary = BoundaryCondition::kPeriodic;
      break;
    case ModeHydroBoundary::kOpen:
      policy.hydro_boundary = BoundaryCondition::kOpen;
      break;
    case ModeHydroBoundary::kReflective:
      policy.hydro_boundary = BoundaryCondition::kReflective;
      break;
    case ModeHydroBoundary::kAuto:
      break;
  }

  switch (mode_config.gravity_boundary) {
    case ModeGravityBoundary::kPeriodic:
      policy.gravity_boundary = GravityBoundaryModel::kPeriodicPoisson;
      break;
    case ModeGravityBoundary::kIsolatedMonopole:
      policy.gravity_boundary = GravityBoundaryModel::kIsolatedMonopoleDirichlet;
      break;
    case ModeGravityBoundary::kAuto:
      break;
  }

  return policy;
}

void validateModePolicy(const SimulationConfig& config, const ModePolicy& policy) {
  if ((config.mode.mode == SimulationMode::kCosmoCube || config.mode.mode == SimulationMode::kZoomIn) &&
      policy.gravity_boundary != GravityBoundaryModel::kPeriodicPoisson) {
    throw ConfigError("cosmological modes require mode.gravity_boundary=periodic");
  }

  if ((config.mode.mode == SimulationMode::kCosmoCube || config.mode.mode == SimulationMode::kZoomIn) &&
      policy.hydro_boundary != BoundaryCondition::kPeriodic) {
    throw ConfigError("cosmological modes require mode.hydro_boundary=periodic");
  }

  if (isIsolatedMode(config.mode.mode) && policy.gravity_boundary == GravityBoundaryModel::kPeriodicPoisson) {
    throw ConfigError("isolated modes require non-periodic gravity boundary treatment");
  }

  if (config.mode.mode == SimulationMode::kZoomIn && config.mode.zoom_high_res_region &&
      config.mode.zoom_region_file.empty()) {
    throw ConfigError("mode.zoom_region_file is required when mode.zoom_high_res_region is true");
  }
}

std::string boundaryConditionToString(BoundaryCondition boundary_condition) {
  switch (boundary_condition) {
    case BoundaryCondition::kPeriodic:
      return "periodic";
    case BoundaryCondition::kOpen:
      return "open";
    case BoundaryCondition::kReflective:
      return "reflective";
  }
  return "unknown";
}

std::string gravityBoundaryModelToString(GravityBoundaryModel boundary_model) {
  switch (boundary_model) {
    case GravityBoundaryModel::kPeriodicPoisson:
      return "periodic";
    case GravityBoundaryModel::kIsolatedMonopoleDirichlet:
      return "isolated_monopole";
  }
  return "unknown";
}

}  // namespace cosmosim::core
