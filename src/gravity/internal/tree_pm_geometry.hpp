#pragma once

#include <algorithm>
#include <cmath>
#include <limits>

namespace cosmosim::gravity::internal {

// Certify once per target/axis from actual source or node bounds. Tree centers
// and COMs may live in an unwrapped frame; canonical particles do not certify
// their deltas. No normalization or new coordinate lanes are needed.
struct TreePmPeriodicMetric {
  double length = 0.0;
  bool bounded = false;
  [[nodiscard]] double operator()(double delta) const noexcept {
    if (!(length > 0.0)) return delta;
    const double half = 0.5 * length;
    // Division can round an adjacent representable ratio onto a half tie.
    // Retain the established nearbyint result in this narrow boundary band.
    if (!bounded || std::abs(std::abs(delta) - half) <=
        4.0 * std::numeric_limits<double>::epsilon() * length)
      return delta - length * std::nearbyint(delta / length);
    return delta > half ? delta - length : (delta < -half ? delta + length : delta);
  }
};

[[nodiscard]] inline TreePmPeriodicMetric treePmPeriodicMetric(
    double minimum, double maximum, double target, double length,
    bool round_to_nearest = true) noexcept {
  return {length, round_to_nearest && std::isfinite(minimum) && std::isfinite(maximum) &&
      minimum - target >= -length && maximum - target <= length};
}

// Deltas have already used each axis's general nearest-image operation.
// The maximum remains the established conservative unwrapped-cube envelope.
struct TreePmAabbDistances {
  double minimum_squared;
  double absolute_x, absolute_y, absolute_z;
  TreePmAabbDistances(double dx, double dy, double dz, double h) noexcept {
    const double x = std::abs(dx), y = std::abs(dy), z = std::abs(dz);
    const double ex = std::max(0.0, x-h), ey = std::max(0.0, y-h), ez = std::max(0.0, z-h);
    minimum_squared = ex*ex + ey*ey + ez*ez;
    absolute_x=x;absolute_y=y;absolute_z=z;
  }
  [[nodiscard]] double maximumSquared(double h) const noexcept {
    const double mx = absolute_x+h, my = absolute_y+h, mz = absolute_z+h;
    return mx*mx + my*my + mz*mz;
  }
};

[[nodiscard]] inline bool treePmSquaredDistanceWithinCutoff(
    double distance2, double cutoff, double cutoff2) noexcept {
  // Preserve the old sqrt comparison at rounding-sensitive boundaries and
  // extreme scales; ordinary rejections/containment need no square root.
  if (!std::isnormal(cutoff2) || !std::isfinite(distance2) ||
      std::abs(distance2-cutoff2) <= 8.0*std::numeric_limits<double>::epsilon()*cutoff2) {
    return std::sqrt(distance2) <= cutoff;
  }
  return distance2 <= cutoff2;
}

}  // namespace cosmosim::gravity::internal
