#pragma once

#include <algorithm>
#include <cmath>
#include <limits>

namespace cosmosim::gravity::internal {

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
