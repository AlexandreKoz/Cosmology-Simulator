#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "periodic_ewald_reference.hpp"

namespace {

using cosmosim::test_support::PeriodicEwaldBox;
using cosmosim::test_support::PeriodicEwaldOptions;
using cosmosim::test_support::PeriodicEwaldSource;
using cosmosim::test_support::PeriodicEwaldTarget;
using cosmosim::test_support::PeriodicEwaldVector3;
using cosmosim::test_support::periodicEwaldAccelerations;
using cosmosim::test_support::periodicEwaldAccelerationsAtSources;

void requireOrThrow(bool condition, const std::string& message) {
  if (!condition) {
    throw std::runtime_error(message);
  }
}

template <typename Callable>
void requireInvalidArgument(Callable&& callable, const std::string& label) {
  bool rejected = false;
  try {
    callable();
  } catch (const std::invalid_argument&) {
    rejected = true;
  }
  requireOrThrow(rejected, label + " must throw std::invalid_argument");
}

[[nodiscard]] double squaredNorm(const PeriodicEwaldVector3& value) {
  return value.x * value.x + value.y * value.y + value.z * value.z;
}

[[nodiscard]] double norm(const PeriodicEwaldVector3& value) {
  return std::sqrt(squaredNorm(value));
}

[[nodiscard]] PeriodicEwaldVector3 difference(
    const PeriodicEwaldVector3& lhs,
    const PeriodicEwaldVector3& rhs) {
  return {.x = lhs.x - rhs.x, .y = lhs.y - rhs.y, .z = lhs.z - rhs.z};
}

[[nodiscard]] double relativeL2(
    const std::vector<PeriodicEwaldVector3>& lhs,
    const std::vector<PeriodicEwaldVector3>& rhs) {
  requireOrThrow(lhs.size() == rhs.size(), "relativeL2 input sizes must match");
  double difference_squared = 0.0;
  double reference_squared = 0.0;
  for (std::size_t index = 0; index < lhs.size(); ++index) {
    difference_squared += squaredNorm(difference(lhs[index], rhs[index]));
    reference_squared += squaredNorm(rhs[index]);
  }
  return std::sqrt(difference_squared / std::max(reference_squared, 1.0e-300));
}

[[nodiscard]] std::string convergenceMessage(
    double loose_to_medium,
    double medium_to_strict,
    double strict_to_stricter,
    double alpha_independence) {
  std::ostringstream message;
  message << std::setprecision(17)
          << "periodic Ewald convergence metrics: loose_to_medium=" << loose_to_medium
          << ", medium_to_strict=" << medium_to_strict
          << ", strict_to_stricter=" << strict_to_stricter
          << ", strict_alpha_independence=" << alpha_independence;
  return message.str();
}

[[nodiscard]] PeriodicEwaldOptions strictOptions() {
  return PeriodicEwaldOptions{
      .gravitational_constant = 0.73,
      .alpha_inverse_length = 2.3,
      .real_image_limits = {2, 2, 2},
      .reciprocal_mode_limits = {7, 9, 12},
  };
}

[[nodiscard]] std::vector<PeriodicEwaldSource> convergenceFixture() {
  return {
      {.position = {.x = 0.08, .y = 0.17, .z = 0.29}, .mass = 0.7},
      {.position = {.x = 0.91, .y = 1.12, .z = 1.51}, .mass = 1.1},
      {.position = {.x = 0.43, .y = 0.76, .z = 0.64}, .mass = 0.9},
      {.position = {.x = 0.67, .y = 0.31, .z = 1.07}, .mass = 1.4},
      {.position = {.x = 0.22, .y = 1.01, .z = 0.93}, .mass = 0.6},
  };
}

void testIncreasingStrictnessConvergenceAndAlphaIndependence() {
  const PeriodicEwaldBox box{.length_x = 1.0, .length_y = 1.3, .length_z = 1.7};
  const std::vector<PeriodicEwaldSource> sources = convergenceFixture();

  PeriodicEwaldOptions loose = strictOptions();
  loose.real_image_limits = {0, 0, 0};
  loose.reciprocal_mode_limits = {1, 1, 1};

  PeriodicEwaldOptions medium = strictOptions();
  medium.real_image_limits = {1, 1, 1};
  medium.reciprocal_mode_limits = {4, 5, 7};

  const PeriodicEwaldOptions strict = strictOptions();

  PeriodicEwaldOptions stricter = strict;
  stricter.real_image_limits = {3, 3, 3};
  stricter.reciprocal_mode_limits = {10, 13, 17};

  PeriodicEwaldOptions alternate_alpha = stricter;
  alternate_alpha.alpha_inverse_length = 3.1;
  alternate_alpha.reciprocal_mode_limits = {12, 16, 21};

  const auto loose_acceleration = periodicEwaldAccelerationsAtSources(sources, box, loose);
  const auto medium_acceleration = periodicEwaldAccelerationsAtSources(sources, box, medium);
  const auto strict_acceleration = periodicEwaldAccelerationsAtSources(sources, box, strict);
  const auto stricter_acceleration = periodicEwaldAccelerationsAtSources(sources, box, stricter);
  const auto alternate_alpha_acceleration = periodicEwaldAccelerationsAtSources(sources, box, alternate_alpha);

  const double loose_to_medium = relativeL2(loose_acceleration, medium_acceleration);
  const double medium_to_strict = relativeL2(medium_acceleration, strict_acceleration);
  const double strict_to_stricter = relativeL2(strict_acceleration, stricter_acceleration);
  const double alpha_independence = relativeL2(stricter_acceleration, alternate_alpha_acceleration);
  const std::string metrics = convergenceMessage(
      loose_to_medium, medium_to_strict, strict_to_stricter, alpha_independence);
  std::cout << metrics << '\n';

  requireOrThrow(loose_to_medium > 1.0e-6, "loose Ewald configuration was not discriminating; " + metrics);
  requireOrThrow(
      medium_to_strict < 0.05 * loose_to_medium,
      "Ewald error did not contract from loose to medium/strict; " + metrics);
  requireOrThrow(
      strict_to_stricter < 0.05 * medium_to_strict,
      "Ewald error did not contract from strict to stricter; " + metrics);
  requireOrThrow(strict_to_stricter < 5.0e-12, "strict Ewald reference did not converge; " + metrics);
  requireOrThrow(
      alpha_independence < 5.0e-11,
      "converged Ewald force depends on splitting alpha; " + metrics);
}

void testIntegerBoxTranslationInvariance() {
  const PeriodicEwaldBox box{.length_x = 1.0, .length_y = 1.3, .length_z = 1.7};
  const std::vector<PeriodicEwaldSource> sources = convergenceFixture();
  std::vector<PeriodicEwaldSource> translated = sources;
  const int shifts[][3] = {
      {2, -3, 1},
      {-1, 4, -2},
      {5, 1, 3},
      {-4, -2, 2},
      {3, 5, -4},
  };
  for (std::size_t index = 0; index < translated.size(); ++index) {
    translated[index].position.x += static_cast<double>(shifts[index][0]) * box.length_x;
    translated[index].position.y += static_cast<double>(shifts[index][1]) * box.length_y;
    translated[index].position.z += static_cast<double>(shifts[index][2]) * box.length_z;
  }

  const PeriodicEwaldOptions options = strictOptions();
  const auto reference = periodicEwaldAccelerationsAtSources(sources, box, options);
  const auto shifted = periodicEwaldAccelerationsAtSources(translated, box, options);
  const double relative_error = relativeL2(shifted, reference);
  requireOrThrow(
      relative_error < 2.0e-13,
      "periodic Ewald acceleration changed under independent integer-box translations: " +
          std::to_string(relative_error));
}

void testZeroSelfForceAndSymmetricCancellation() {
  const PeriodicEwaldBox box{.length_x = 2.0, .length_y = 3.0, .length_z = 5.0};
  PeriodicEwaldOptions options{
      .gravitational_constant = 1.25,
      .alpha_inverse_length = 1.7,
      .real_image_limits = {3, 3, 3},
      .reciprocal_mode_limits = {9, 13, 20},
  };

  const std::vector<PeriodicEwaldSource> single_source{
      {.position = {.x = 0.37, .y = 1.21, .z = 3.44}, .mass = 2.5},
  };
  const auto self_acceleration = periodicEwaldAccelerationsAtSources(single_source, box, options);
  requireOrThrow(self_acceleration.size() == 1, "single-source Ewald output size mismatch");
  requireOrThrow(norm(self_acceleration[0]) < 2.0e-14, "periodic Ewald self force is nonzero");

  const std::vector<PeriodicEwaldTarget> translated_self_target{
      {.position = {.x = 2.37, .y = -1.79, .z = 13.44}, .self_source_index = 0},
  };
  const auto translated_self = periodicEwaldAccelerations(
      single_source, translated_self_target, box, options);
  requireOrThrow(norm(translated_self[0]) < 2.0e-14, "translated periodic self force is nonzero");

  const PeriodicEwaldVector3 center{.x = 1.0, .y = 1.5, .z = 2.5};
  const std::vector<PeriodicEwaldSource> symmetric_sources{
      {.position = {.x = 0.75, .y = 1.125, .z = 1.875}, .mass = 1.4},
      {.position = {.x = 1.25, .y = 1.875, .z = 3.125}, .mass = 1.4},
  };
  const std::vector<PeriodicEwaldTarget> center_target{{.position = center}};
  const auto center_acceleration = periodicEwaldAccelerations(
      symmetric_sources, center_target, box, options);
  requireOrThrow(
      norm(center_acceleration[0]) < 5.0e-14,
      "equal sources symmetric about a target did not cancel");
}

void testMassWeightedNetForceCancellation() {
  const PeriodicEwaldBox box{.length_x = 1.0, .length_y = 1.3, .length_z = 1.7};
  const std::vector<PeriodicEwaldSource> sources = convergenceFixture();
  const auto accelerations = periodicEwaldAccelerationsAtSources(sources, box, strictOptions());

  PeriodicEwaldVector3 net_force{};
  double force_scale = 0.0;
  for (std::size_t index = 0; index < sources.size(); ++index) {
    net_force.x += sources[index].mass * accelerations[index].x;
    net_force.y += sources[index].mass * accelerations[index].y;
    net_force.z += sources[index].mass * accelerations[index].z;
    force_scale += sources[index].mass * norm(accelerations[index]);
  }
  requireOrThrow(
      norm(net_force) <= 5.0e-14 * std::max(force_scale, 1.0),
      "periodic Ewald mass-weighted net force does not cancel");
}

void testGravitationalConstantNormalizationAndAttractionSign() {
  const PeriodicEwaldBox box{.length_x = 4.0, .length_y = 5.0, .length_z = 6.0};
  const std::vector<PeriodicEwaldSource> source{
      {.position = {.x = 1.2, .y = 2.0, .z = 3.0}, .mass = 1.7},
  };
  const std::vector<PeriodicEwaldTarget> target{
      {.position = {.x = 1.0, .y = 2.0, .z = 3.0}},
  };
  PeriodicEwaldOptions unit_g{
      .gravitational_constant = 1.0,
      .alpha_inverse_length = 1.1,
      .real_image_limits = {2, 2, 2},
      .reciprocal_mode_limits = {9, 11, 13},
  };
  PeriodicEwaldOptions scaled_g = unit_g;
  scaled_g.gravitational_constant = 2.75;

  const auto unit_acceleration = periodicEwaldAccelerations(source, target, box, unit_g);
  const auto scaled_acceleration = periodicEwaldAccelerations(source, target, box, scaled_g);
  const PeriodicEwaldVector3 expected_scaled{
      .x = 2.75 * unit_acceleration[0].x,
      .y = 2.75 * unit_acceleration[0].y,
      .z = 2.75 * unit_acceleration[0].z,
  };
  const double scaling_error = norm(difference(scaled_acceleration[0], expected_scaled)) /
      std::max(norm(expected_scaled), 1.0e-300);
  requireOrThrow(unit_acceleration[0].x > 0.0, "periodic Ewald force has the repulsive sign");
  requireOrThrow(scaling_error < 5.0e-15, "periodic Ewald force does not scale linearly with G");
}

void testInputValidationAndUnsoftenedEnvelope() {
  const PeriodicEwaldBox box{.length_x = 1.0, .length_y = 1.3, .length_z = 1.7};
  const std::vector<PeriodicEwaldSource> source{
      {.position = {.x = 0.2, .y = 0.3, .z = 0.4}, .mass = 1.0},
  };
  const std::vector<PeriodicEwaldTarget> probe{
      {.position = {.x = 0.6, .y = 0.7, .z = 0.8}},
  };
  const PeriodicEwaldOptions valid = strictOptions();

  PeriodicEwaldBox invalid_box = box;
  invalid_box.length_y = 0.0;
  requireInvalidArgument(
      [&]() { (void)periodicEwaldAccelerations(source, probe, invalid_box, valid); },
      "zero box length");

  PeriodicEwaldOptions invalid_alpha = valid;
  invalid_alpha.alpha_inverse_length = std::numeric_limits<double>::quiet_NaN();
  requireInvalidArgument(
      [&]() { (void)periodicEwaldAccelerations(source, probe, box, invalid_alpha); },
      "non-finite alpha");

  PeriodicEwaldOptions invalid_truncation = valid;
  invalid_truncation.real_image_limits[1] = -1;
  requireInvalidArgument(
      [&]() { (void)periodicEwaldAccelerations(source, probe, box, invalid_truncation); },
      "negative real-image truncation");

  std::vector<PeriodicEwaldSource> invalid_mass = source;
  invalid_mass[0].mass = -1.0;
  requireInvalidArgument(
      [&]() { (void)periodicEwaldAccelerations(invalid_mass, probe, box, valid); },
      "negative source mass");

  std::vector<PeriodicEwaldSource> invalid_position = source;
  invalid_position[0].position.z = std::numeric_limits<double>::infinity();
  requireInvalidArgument(
      [&]() { (void)periodicEwaldAccelerations(invalid_position, probe, box, valid); },
      "non-finite source position");

  std::vector<PeriodicEwaldTarget> invalid_self = probe;
  invalid_self[0].self_source_index = 1;
  requireInvalidArgument(
      [&]() { (void)periodicEwaldAccelerations(source, invalid_self, box, valid); },
      "out-of-range self source");

  invalid_self[0].self_source_index = 0;
  requireInvalidArgument(
      [&]() { (void)periodicEwaldAccelerations(source, invalid_self, box, valid); },
      "mislocated self source");

  const std::vector<PeriodicEwaldTarget> coincident_probe{{.position = source[0].position}};
  requireInvalidArgument(
      [&]() { (void)periodicEwaldAccelerations(source, coincident_probe, box, valid); },
      "coincident unsoftened non-self pair");

  const std::vector<PeriodicEwaldSource> empty_sources;
  const auto empty_acceleration = periodicEwaldAccelerations(empty_sources, probe, box, valid);
  requireOrThrow(
      empty_acceleration.size() == 1 && squaredNorm(empty_acceleration[0]) == 0.0,
      "empty source set must produce exact zero acceleration");
}

}  // namespace

int main() {
  testIncreasingStrictnessConvergenceAndAlphaIndependence();
  testIntegerBoxTranslationInvariance();
  testZeroSelfForceAndSymmetricCancellation();
  testMassWeightedNetForceCancellation();
  testGravitationalConstantNormalizationAndAttractionSign();
  testInputValidationAndUnsoftenedEnvelope();
  return 0;
}
