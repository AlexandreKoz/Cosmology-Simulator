#pragma once

#include <cmath>
#include <stdexcept>

namespace cosmosim::core {

[[nodiscard]] inline double canonicalPeriodicPosition(double value, double box_length) {
  if (!std::isfinite(value) || !std::isfinite(box_length) || box_length <= 0.0) {
    throw std::invalid_argument("periodic canonicalization requires finite value and positive box length");
  }
  double wrapped = std::fmod(value, box_length);
  if (wrapped < 0.0) wrapped += box_length;
  if (!(wrapped < box_length)) wrapped = 0.0;
  return wrapped;
}

[[nodiscard]] inline double minimumImageDelta(double delta, double box_length) {
  if (!std::isfinite(delta) || !std::isfinite(box_length) || box_length <= 0.0) {
    throw std::invalid_argument("minimum-image displacement requires finite delta and positive box length");
  }
  return delta - box_length * std::nearbyint(delta / box_length);
}

[[nodiscard]] inline bool containsCanonicalPeriodicPosition(
    double value, double box_length, double tolerance = 0.0) noexcept {
  return std::isfinite(value) && std::isfinite(box_length) && box_length > 0.0 &&
      std::isfinite(tolerance) && tolerance >= 0.0 &&
      value >= -tolerance && value < box_length + tolerance;
}

}  // namespace cosmosim::core
