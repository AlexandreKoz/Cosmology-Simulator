#pragma once

#include <algorithm>
#include <cmath>

#include "cosmosim/gravity/tree_gravity.hpp"
#include "cosmosim/gravity/tree_softening.hpp"

namespace cosmosim::gravity::internal {


[[nodiscard]] inline bool targetInsideNodeAabbFromCenterDelta(
    double delta_x,
    double delta_y,
    double delta_z,
    double half_size) noexcept {
  return std::abs(delta_x) <= half_size &&
      std::abs(delta_y) <= half_size &&
      std::abs(delta_z) <= half_size;
}

struct TreeNodeAcceptanceInput {
  bool is_leaf = false;
  bool target_inside_node = false;
  double half_size = 0.0;
  double com_center_offset = 0.0;
  double node_mass_code = 0.0;
  double r2 = 0.0;
  bool previous_acceleration_available = false;
  double previous_acceleration_magnitude_code = 0.0;
  double target_softening_comoving = 0.0;
  double node_softening_min_comoving = 0.0;
  double node_softening_max_comoving = 0.0;
};

[[nodiscard]] inline bool acceptNodeByComDistanceMac(
    double theta,
    double half_size,
    double com_center_offset,
    double r2) noexcept {
  const double r = std::sqrt(r2 + 1.0e-30);
  return ((2.0 * half_size + com_center_offset) / r) < theta;
}

[[nodiscard]] inline bool acceptNodeByMac(
    bool is_leaf,
    bool target_inside_node,
    double half_size,
    double com_center_offset,
    double node_mass_code,
    double r2,
    bool previous_acceleration_available,
    double previous_acceleration_magnitude_code,
    const TreeGravityOptions& options) noexcept {
  if (is_leaf) return true;
  if (target_inside_node) return false;
  const double width = 2.0 * half_size;
  if (options.opening_criterion == TreeOpeningCriterion::kBarnesHutGeometric) {
    return (width / std::sqrt(r2 + 1.0e-30)) < options.opening_theta;
  }
  if (options.opening_criterion == TreeOpeningCriterion::kBarnesHutComDistance) {
    return acceptNodeByComDistanceMac(options.opening_theta, half_size, com_center_offset, r2);
  }
  if (!previous_acceleration_available) {
    return acceptNodeByComDistanceMac(options.opening_theta, half_size, com_center_offset, r2);
  }
  const double acceleration_scale_code = std::max(
      std::abs(previous_acceleration_magnitude_code),
      options.relative_force_acceleration_floor_code);
  const double estimated_error_scale =
      options.gravitational_constant_code * node_mass_code * width * width;
  const double allowed_error_scale =
      options.relative_force_tolerance * acceleration_scale_code * r2 * r2;
  return estimated_error_scale <= allowed_error_scale;
}

[[nodiscard]] inline bool passesSofteningEnvelopeGuard(
    bool is_leaf,
    double half_size,
    double r,
    double target_softening_comoving,
    double node_softening_min_comoving,
    double node_softening_max_comoving) noexcept {
  if (is_leaf) return true;
  if (node_softening_max_comoving - node_softening_min_comoving <= 1.0e-12) return true;
  const double pair_softening_max = combineSofteningPairEpsilonUnchecked(
      node_softening_max_comoving, target_softening_comoving);
  return r > (2.0 * half_size + 2.0 * pair_softening_max);
}


[[nodiscard]] inline bool acceptNodeByCommonTreePolicy(
    const TreeNodeAcceptanceInput& input,
    const TreeGravityOptions& options) noexcept {
  const double r = std::sqrt(input.r2 + 1.0e-30);
  return acceptNodeByMac(
             input.is_leaf,
             input.target_inside_node,
             input.half_size,
             input.com_center_offset,
             input.node_mass_code,
             input.r2,
             input.previous_acceleration_available,
             input.previous_acceleration_magnitude_code,
             options) &&
      passesSofteningEnvelopeGuard(
          input.is_leaf,
          input.half_size,
          r,
          input.target_softening_comoving,
          input.node_softening_min_comoving,
          input.node_softening_max_comoving);
}


}  // namespace cosmosim::gravity::internal
