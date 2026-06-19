#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace cosmosim::tests::amr_hydro_campaign {

// Code-neutral line or radial profile.  Consumers must compare profiles only
// after metadata compatibility has been checked; AMR rows are never assumed to
// share topology or ordering.
struct ProfilePoint {
  double coordinate_code = 0.0;
  double rho_code = 0.0;
  double velocity_code = 0.0;
  double pressure_code = 0.0;
};

struct ComparisonMetadata {
  std::string schema_version;
  std::string case_name;
  std::string geometry;  // "line" (Sod) or "radial" (Sedov).
  std::string source_code_name;
  std::string source_code_version;
  double gamma = 0.0;
  double final_time_code = 0.0;
  std::string units;
  std::string boundary_conditions;
  bool gravity_enabled = false;
  bool cooling_enabled = false;
  bool cosmological_expansion_enabled = false;
};

struct ErrorNorms {
  double l1_density = 0.0;
  double l1_velocity = 0.0;
  double l1_pressure = 0.0;
  double l2_density = 0.0;
};

inline void require(bool condition, std::string_view message) {
  if (!condition) {
    throw std::runtime_error(std::string(message));
  }
}

inline void validateProfile(const std::vector<ProfilePoint>& points, std::string_view label) {
  require(points.size() >= 2U, std::string(label) + ": profile needs at least two samples");
  for (std::size_t i = 0; i < points.size(); ++i) {
    const auto& point = points[i];
    require(std::isfinite(point.coordinate_code) && std::isfinite(point.rho_code) &&
            std::isfinite(point.velocity_code) && std::isfinite(point.pressure_code),
            std::string(label) + ": non-finite profile value");
    require(point.rho_code > 0.0 && point.pressure_code > 0.0,
            std::string(label) + ": non-positive density or pressure");
    if (i > 0U) {
      require(points[i - 1U].coordinate_code < point.coordinate_code,
              std::string(label) + ": coordinates must be strictly increasing");
    }
  }
}

inline void validateComparable(const ComparisonMetadata& lhs, const ComparisonMetadata& rhs) {
  require(lhs.schema_version == "cosmosim_amr_profile_v1" && rhs.schema_version == "cosmosim_amr_profile_v1",
          "comparison metadata schema mismatch");
  require(lhs.case_name == rhs.case_name, "comparison case mismatch");
  require(lhs.geometry == rhs.geometry, "comparison geometry mismatch");
  require(lhs.units == rhs.units, "comparison units mismatch; explicit conversion is required");
  require(std::abs(lhs.gamma - rhs.gamma) <= 1.0e-14, "comparison EOS gamma mismatch");
  require(std::abs(lhs.final_time_code - rhs.final_time_code) <= 1.0e-12,
          "comparison final time mismatch");
  require(lhs.boundary_conditions == rhs.boundary_conditions, "comparison boundary-condition mismatch");
  require(lhs.gravity_enabled == rhs.gravity_enabled && lhs.cooling_enabled == rhs.cooling_enabled &&
              lhs.cosmological_expansion_enabled == rhs.cosmological_expansion_enabled,
          "comparison physics switch mismatch");
}

[[nodiscard]] inline ProfilePoint interpolateLinear(
    const std::vector<ProfilePoint>& profile,
    double coordinate_code) {
  validateProfile(profile, "interpolateLinear input");
  require(coordinate_code >= profile.front().coordinate_code && coordinate_code <= profile.back().coordinate_code,
          "interpolation coordinate outside profile support");
  auto upper = std::lower_bound(
      profile.begin(), profile.end(), coordinate_code,
      [](const ProfilePoint& point, double coordinate) { return point.coordinate_code < coordinate; });
  if (upper == profile.begin() || upper == profile.end()) {
    return upper == profile.end() ? profile.back() : profile.front();
  }
  const auto lower = upper - 1;
  const double alpha = (coordinate_code - lower->coordinate_code) /
      (upper->coordinate_code - lower->coordinate_code);
  return ProfilePoint{
      .coordinate_code = coordinate_code,
      .rho_code = (1.0 - alpha) * lower->rho_code + alpha * upper->rho_code,
      .velocity_code = (1.0 - alpha) * lower->velocity_code + alpha * upper->velocity_code,
      .pressure_code = (1.0 - alpha) * lower->pressure_code + alpha * upper->pressure_code};
}

[[nodiscard]] inline ErrorNorms compareOnCommonGrid(
    const std::vector<ProfilePoint>& numerical,
    const std::vector<ProfilePoint>& reference,
    const std::vector<double>& coordinates) {
  validateProfile(numerical, "numerical");
  validateProfile(reference, "reference");
  require(!coordinates.empty(), "comparison requires a common coordinate grid");
  ErrorNorms result;
  for (const double coordinate : coordinates) {
    const ProfilePoint n = interpolateLinear(numerical, coordinate);
    const ProfilePoint r = interpolateLinear(reference, coordinate);
    result.l1_density += std::abs(n.rho_code - r.rho_code);
    result.l1_velocity += std::abs(n.velocity_code - r.velocity_code);
    result.l1_pressure += std::abs(n.pressure_code - r.pressure_code);
    result.l2_density += (n.rho_code - r.rho_code) * (n.rho_code - r.rho_code);
  }
  const double inv_count = 1.0 / static_cast<double>(coordinates.size());
  result.l1_density *= inv_count;
  result.l1_velocity *= inv_count;
  result.l1_pressure *= inv_count;
  result.l2_density = std::sqrt(result.l2_density * inv_count);
  return result;
}

[[nodiscard]] inline double observedOrder(double coarse_error, double fine_error, double refinement_ratio) {
  require(std::isfinite(coarse_error) && std::isfinite(fine_error) && coarse_error > 0.0 && fine_error > 0.0,
          "observedOrder requires positive finite errors");
  require(std::isfinite(refinement_ratio) && refinement_ratio > 1.0,
          "observedOrder requires refinement ratio greater than one");
  return std::log(coarse_error / fine_error) / std::log(refinement_ratio);
}

[[nodiscard]] inline double strongestGradientCoordinate(const std::vector<ProfilePoint>& profile) {
  validateProfile(profile, "shock locator profile");
  double strongest = -1.0;
  double coordinate = profile.front().coordinate_code;
  for (std::size_t i = 1; i < profile.size(); ++i) {
    const double dx = profile[i].coordinate_code - profile[i - 1U].coordinate_code;
    const double gradient = std::abs((profile[i].rho_code - profile[i - 1U].rho_code) / dx);
    if (gradient > strongest) {
      strongest = gradient;
      coordinate = 0.5 * (profile[i].coordinate_code + profile[i - 1U].coordinate_code);
    }
  }
  return coordinate;
}

}  // namespace cosmosim::tests::amr_hydro_campaign
