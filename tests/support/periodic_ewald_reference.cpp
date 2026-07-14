#include "periodic_ewald_reference.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

namespace cosmosim::test_support {
namespace {

constexpr double k_pi = 3.141592653589793238462643383279502884;
constexpr double k_two_over_sqrt_pi = 1.128379167095512573896158903121545172;

class CompensatedSum {
 public:
  void add(double value) {
    if (!std::isfinite(value)) {
      throw std::overflow_error("periodic Ewald reference produced a non-finite contribution");
    }

    const double next = m_sum + value;
    if (!std::isfinite(next)) {
      throw std::overflow_error("periodic Ewald reference accumulation overflowed");
    }
    if (std::abs(m_sum) >= std::abs(value)) {
      m_correction += (m_sum - next) + value;
    } else {
      m_correction += (value - next) + m_sum;
    }
    if (!std::isfinite(m_correction)) {
      throw std::overflow_error("periodic Ewald reference compensation overflowed");
    }
    m_sum = next;
  }

  [[nodiscard]] double value() const {
    const double result = m_sum + m_correction;
    if (!std::isfinite(result)) {
      throw std::overflow_error("periodic Ewald reference produced a non-finite acceleration");
    }
    return result;
  }

 private:
  double m_sum = 0.0;
  double m_correction = 0.0;
};

struct VectorAccumulator {
  CompensatedSum x;
  CompensatedSum y;
  CompensatedSum z;

  void add(double value_x, double value_y, double value_z) {
    x.add(value_x);
    y.add(value_y);
    z.add(value_z);
  }

  [[nodiscard]] PeriodicEwaldVector3 value() const {
    return {.x = x.value(), .y = y.value(), .z = z.value()};
  }
};

[[nodiscard]] bool finiteVector(const PeriodicEwaldVector3& value) {
  return std::isfinite(value.x) && std::isfinite(value.y) && std::isfinite(value.z);
}

void requireFiniteVector(const PeriodicEwaldVector3& value, const std::string& label, std::size_t index) {
  if (!finiteVector(value)) {
    std::ostringstream message;
    message << "periodic Ewald " << label << '[' << index << "] position must be finite";
    throw std::invalid_argument(message.str());
  }
}

void validateLimits(const std::array<int, 3>& limits, const char* label) {
  for (std::size_t axis = 0; axis < limits.size(); ++axis) {
    if (limits[axis] < 0) {
      std::ostringstream message;
      message << "periodic Ewald " << label << '[' << axis << "] must be non-negative";
      throw std::invalid_argument(message.str());
    }
  }
}

[[nodiscard]] double canonicalCoordinate(double coordinate, double length) {
  double wrapped = std::fmod(coordinate, length);
  if (wrapped < 0.0) {
    wrapped += length;
  }
  if (wrapped == length || wrapped == 0.0) {
    return 0.0;
  }
  return wrapped;
}

[[nodiscard]] PeriodicEwaldVector3 canonicalPosition(
    const PeriodicEwaldVector3& position,
    const PeriodicEwaldBox& box) {
  return {
      .x = canonicalCoordinate(position.x, box.length_x),
      .y = canonicalCoordinate(position.y, box.length_y),
      .z = canonicalCoordinate(position.z, box.length_z),
  };
}

[[nodiscard]] double minimumImageDelta(double source, double target, double length) {
  const double delta = std::remainder(source - target, length);
  return delta == 0.0 ? 0.0 : delta;
}

[[nodiscard]] PeriodicEwaldVector3 minimumImageDisplacement(
    const PeriodicEwaldVector3& source,
    const PeriodicEwaldVector3& target,
    const PeriodicEwaldBox& box) {
  return {
      .x = minimumImageDelta(source.x, target.x, box.length_x),
      .y = minimumImageDelta(source.y, target.y, box.length_y),
      .z = minimumImageDelta(source.z, target.z, box.length_z),
  };
}

[[nodiscard]] bool isZeroModuloBox(
    const PeriodicEwaldVector3& displacement,
    const PeriodicEwaldBox& box) {
  constexpr double k_match_epsilon = 128.0 * std::numeric_limits<double>::epsilon();
  const double tolerance_x = k_match_epsilon * std::max(1.0, box.length_x);
  const double tolerance_y = k_match_epsilon * std::max(1.0, box.length_y);
  const double tolerance_z = k_match_epsilon * std::max(1.0, box.length_z);
  return std::abs(displacement.x) <= tolerance_x &&
      std::abs(displacement.y) <= tolerance_y &&
      std::abs(displacement.z) <= tolerance_z;
}

void validateInputs(
    std::span<const PeriodicEwaldSource> sources,
    std::span<const PeriodicEwaldTarget> targets,
    const PeriodicEwaldBox& box,
    const PeriodicEwaldOptions& options) {
  if (!std::isfinite(box.length_x) || box.length_x <= 0.0 ||
      !std::isfinite(box.length_y) || box.length_y <= 0.0 ||
      !std::isfinite(box.length_z) || box.length_z <= 0.0) {
    throw std::invalid_argument("periodic Ewald box lengths must be finite and positive");
  }

  const double volume = box.length_x * box.length_y * box.length_z;
  if (!std::isfinite(volume) || volume <= 0.0) {
    throw std::invalid_argument("periodic Ewald box volume must be finite and positive");
  }
  if (!std::isfinite(options.gravitational_constant) || options.gravitational_constant <= 0.0) {
    throw std::invalid_argument("periodic Ewald gravitational constant must be finite and positive");
  }
  if (!std::isfinite(options.alpha_inverse_length) || options.alpha_inverse_length <= 0.0) {
    throw std::invalid_argument("periodic Ewald alpha must be finite and positive");
  }
  const double alpha_squared = options.alpha_inverse_length * options.alpha_inverse_length;
  if (!std::isfinite(alpha_squared)) {
    throw std::invalid_argument("periodic Ewald alpha squared must be finite");
  }
  validateLimits(options.real_image_limits, "real_image_limits");
  validateLimits(options.reciprocal_mode_limits, "reciprocal_mode_limits");

  for (std::size_t source_index = 0; source_index < sources.size(); ++source_index) {
    requireFiniteVector(sources[source_index].position, "source", source_index);
    if (!std::isfinite(sources[source_index].mass) || sources[source_index].mass < 0.0) {
      std::ostringstream message;
      message << "periodic Ewald source[" << source_index << "] mass must be finite and non-negative";
      throw std::invalid_argument(message.str());
    }
  }

  for (std::size_t target_index = 0; target_index < targets.size(); ++target_index) {
    requireFiniteVector(targets[target_index].position, "target", target_index);
    const std::size_t self_source_index = targets[target_index].self_source_index;
    if (self_source_index != k_no_periodic_ewald_self_source && self_source_index >= sources.size()) {
      std::ostringstream message;
      message << "periodic Ewald target[" << target_index << "] self_source_index is out of range";
      throw std::invalid_argument(message.str());
    }
  }
}

[[nodiscard]] std::int64_t negativeLimit(int limit) {
  return -static_cast<std::int64_t>(limit);
}

}  // namespace

std::vector<PeriodicEwaldVector3> periodicEwaldAccelerations(
    std::span<const PeriodicEwaldSource> sources,
    std::span<const PeriodicEwaldTarget> targets,
    const PeriodicEwaldBox& box,
    const PeriodicEwaldOptions& options) {
  validateInputs(sources, targets, box, options);

  std::vector<PeriodicEwaldVector3> canonical_sources;
  canonical_sources.reserve(sources.size());
  for (const auto& source : sources) {
    canonical_sources.push_back(canonicalPosition(source.position, box));
  }
  std::vector<PeriodicEwaldVector3> canonical_targets;
  canonical_targets.reserve(targets.size());
  for (const auto& target : targets) {
    canonical_targets.push_back(canonicalPosition(target.position, box));
  }

  const double volume = box.length_x * box.length_y * box.length_z;
  const double alpha = options.alpha_inverse_length;
  const double alpha_squared = alpha * alpha;
  const double reciprocal_prefactor = 4.0 * k_pi * options.gravitational_constant / volume;
  const std::array<double, 3> reciprocal_fundamental{
      2.0 * k_pi / box.length_x,
      2.0 * k_pi / box.length_y,
      2.0 * k_pi / box.length_z,
  };
  if (!std::isfinite(reciprocal_prefactor) ||
      !std::isfinite(reciprocal_fundamental[0]) ||
      !std::isfinite(reciprocal_fundamental[1]) ||
      !std::isfinite(reciprocal_fundamental[2])) {
    throw std::invalid_argument("periodic Ewald reciprocal normalization must be finite");
  }

  std::vector<PeriodicEwaldVector3> accelerations;
  accelerations.reserve(targets.size());
  for (std::size_t target_index = 0; target_index < targets.size(); ++target_index) {
    const auto& target = targets[target_index];
    VectorAccumulator acceleration;

    for (std::size_t source_index = 0; source_index < sources.size(); ++source_index) {
      const auto& source = sources[source_index];
      PeriodicEwaldVector3 displacement = minimumImageDisplacement(
          canonical_sources[source_index], canonical_targets[target_index], box);
      const bool is_self_source = target.self_source_index == source_index;
      if (is_self_source) {
        if (!isZeroModuloBox(displacement, box)) {
          std::ostringstream message;
          message << "periodic Ewald target[" << target_index
                  << "] self source is not at the same periodic position";
          throw std::invalid_argument(message.str());
        }
        displacement = {};
      }
      if (source.mass == 0.0) {
        continue;
      }

      const std::int64_t real_limit_x = options.real_image_limits[0];
      const std::int64_t real_limit_y = options.real_image_limits[1];
      const std::int64_t real_limit_z = options.real_image_limits[2];
      for (std::int64_t image_x = negativeLimit(options.real_image_limits[0]);
           image_x <= real_limit_x;
           ++image_x) {
        for (std::int64_t image_y = negativeLimit(options.real_image_limits[1]);
             image_y <= real_limit_y;
             ++image_y) {
          for (std::int64_t image_z = negativeLimit(options.real_image_limits[2]);
               image_z <= real_limit_z;
               ++image_z) {
            if (is_self_source && image_x == 0 && image_y == 0 && image_z == 0) {
              continue;
            }

            const double separation_x = displacement.x + static_cast<double>(image_x) * box.length_x;
            const double separation_y = displacement.y + static_cast<double>(image_y) * box.length_y;
            const double separation_z = displacement.z + static_cast<double>(image_z) * box.length_z;
            const double radius_squared = separation_x * separation_x +
                separation_y * separation_y + separation_z * separation_z;
            if (radius_squared == 0.0) {
              std::ostringstream message;
              message << "periodic Ewald unsoftened singularity between target[" << target_index
                      << "] and non-self source[" << source_index << ']';
              throw std::invalid_argument(message.str());
            }
            if (!std::isfinite(radius_squared)) {
              throw std::overflow_error("periodic Ewald real-space separation overflowed");
            }

            const double radius = std::sqrt(radius_squared);
            const double alpha_radius = alpha * radius;
            const double screening = std::erfc(alpha_radius) +
                k_two_over_sqrt_pi * alpha_radius * std::exp(-alpha_radius * alpha_radius);
            const double factor = options.gravitational_constant * source.mass * screening /
                (radius_squared * radius);
            acceleration.add(
                factor * separation_x,
                factor * separation_y,
                factor * separation_z);
          }
        }
      }

      // The reciprocal self term is identically zero because sin(k . 0) = 0.
      // Skipping it explicitly avoids manufacturing roundoff self force.
      if (is_self_source) {
        continue;
      }

      const std::int64_t reciprocal_limit_x = options.reciprocal_mode_limits[0];
      const std::int64_t reciprocal_limit_y = options.reciprocal_mode_limits[1];
      const std::int64_t reciprocal_limit_z = options.reciprocal_mode_limits[2];
      for (std::int64_t mode_x = negativeLimit(options.reciprocal_mode_limits[0]);
           mode_x <= reciprocal_limit_x;
           ++mode_x) {
        const double wave_x = static_cast<double>(mode_x) * reciprocal_fundamental[0];
        for (std::int64_t mode_y = negativeLimit(options.reciprocal_mode_limits[1]);
             mode_y <= reciprocal_limit_y;
             ++mode_y) {
          const double wave_y = static_cast<double>(mode_y) * reciprocal_fundamental[1];
          for (std::int64_t mode_z = negativeLimit(options.reciprocal_mode_limits[2]);
               mode_z <= reciprocal_limit_z;
               ++mode_z) {
            if (mode_x == 0 && mode_y == 0 && mode_z == 0) {
              continue;
            }

            const double wave_z = static_cast<double>(mode_z) * reciprocal_fundamental[2];
            const double wave_number_squared =
                wave_x * wave_x + wave_y * wave_y + wave_z * wave_z;
            if (!std::isfinite(wave_number_squared) || wave_number_squared <= 0.0) {
              throw std::overflow_error("periodic Ewald reciprocal wave number overflowed");
            }
            const double damping = std::exp(-wave_number_squared / (4.0 * alpha_squared));
            const double phase = wave_x * displacement.x +
                wave_y * displacement.y + wave_z * displacement.z;
            const double factor = reciprocal_prefactor * source.mass * damping *
                std::sin(phase) / wave_number_squared;
            acceleration.add(factor * wave_x, factor * wave_y, factor * wave_z);
          }
        }
      }
    }

    accelerations.push_back(acceleration.value());
  }
  return accelerations;
}

std::vector<PeriodicEwaldVector3> periodicEwaldAccelerationsAtSources(
    std::span<const PeriodicEwaldSource> sources,
    const PeriodicEwaldBox& box,
    const PeriodicEwaldOptions& options) {
  std::vector<PeriodicEwaldTarget> targets;
  targets.reserve(sources.size());
  for (std::size_t source_index = 0; source_index < sources.size(); ++source_index) {
    targets.push_back({.position = sources[source_index].position, .self_source_index = source_index});
  }
  return periodicEwaldAccelerations(sources, targets, box, options);
}

}  // namespace cosmosim::test_support
