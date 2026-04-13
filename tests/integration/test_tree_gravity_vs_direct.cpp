#include <algorithm>
#include <cmath>
#include <cstdint>
#include <span>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "cosmosim/gravity/tree_gravity.hpp"

namespace {

void requireOrThrow(bool condition, const std::string& message) {
  if (!condition) {
    throw std::runtime_error(message);
  }
}

void directSumAcceleration(
    std::span<const double> pos_x,
    std::span<const double> pos_y,
    std::span<const double> pos_z,
    std::span<const double> mass,
    std::span<const std::uint32_t> active,
    const cosmosim::gravity::TreeGravityOptions& options,
    std::span<double> ax,
    std::span<double> ay,
    std::span<double> az) {
  for (std::size_t i = 0; i < active.size(); ++i) {
    const std::uint32_t target = active[i];
    const double tx = pos_x[target];
    const double ty = pos_y[target];
    const double tz = pos_z[target];
    double acc_x = 0.0;
    double acc_y = 0.0;
    double acc_z = 0.0;
    for (std::size_t j = 0; j < mass.size(); ++j) {
      if (j == target) {
        continue;
      }
      const double dx = pos_x[j] - tx;
      const double dy = pos_y[j] - ty;
      const double dz = pos_z[j] - tz;
      const double r2 = dx * dx + dy * dy + dz * dz;
      const double factor = options.gravitational_constant_code * mass[j] *
          cosmosim::gravity::softenedInvR3(r2, options.softening);
      acc_x += factor * dx;
      acc_y += factor * dy;
      acc_z += factor * dz;
    }
    ax[i] = acc_x;
    ay[i] = acc_y;
    az[i] = acc_z;
  }
}

[[nodiscard]] double solveAndMeasureMaxRelativeError(
    std::span<const double> pos_x,
    std::span<const double> pos_y,
    std::span<const double> pos_z,
    std::span<const double> mass,
    std::span<const std::uint32_t> active,
    const cosmosim::gravity::TreeGravityOptions& options,
    cosmosim::gravity::TreeGravityProfile* profile,
    double* out_mean_relative_error) {
  cosmosim::gravity::TreeGravitySolver solver;
  solver.build(pos_x, pos_y, pos_z, mass, options, profile);

  std::vector<double> tree_ax(active.size(), 0.0);
  std::vector<double> tree_ay(active.size(), 0.0);
  std::vector<double> tree_az(active.size(), 0.0);
  solver.evaluateActiveSet(pos_x, pos_y, pos_z, mass, active, tree_ax, tree_ay, tree_az, options, profile);

  std::vector<double> direct_ax(active.size(), 0.0);
  std::vector<double> direct_ay(active.size(), 0.0);
  std::vector<double> direct_az(active.size(), 0.0);
  directSumAcceleration(pos_x, pos_y, pos_z, mass, active, options, direct_ax, direct_ay, direct_az);

  double max_relative_error = 0.0;
  double sum_relative_error = 0.0;
  for (std::size_t i = 0; i < active.size(); ++i) {
    const double direct_norm = std::sqrt(
        direct_ax[i] * direct_ax[i] + direct_ay[i] * direct_ay[i] + direct_az[i] * direct_az[i]);
    const double diff_norm = std::sqrt(
        (tree_ax[i] - direct_ax[i]) * (tree_ax[i] - direct_ax[i]) +
        (tree_ay[i] - direct_ay[i]) * (tree_ay[i] - direct_ay[i]) +
        (tree_az[i] - direct_az[i]) * (tree_az[i] - direct_az[i]));
    const double denom = std::max(direct_norm, 1.0e-12);
    const double relative_error = diff_norm / denom;
    max_relative_error = std::max(max_relative_error, relative_error);
    sum_relative_error += relative_error;
  }

  if (out_mean_relative_error != nullptr) {
    *out_mean_relative_error = sum_relative_error / static_cast<double>(active.size());
  }
  return max_relative_error;
}

}  // namespace

int main() {
  constexpr std::size_t particle_count = 96;
  std::vector<double> pos_x(particle_count, 0.0);
  std::vector<double> pos_y(particle_count, 0.0);
  std::vector<double> pos_z(particle_count, 0.0);
  std::vector<double> mass(particle_count, 0.0);

  for (std::size_t i = 0; i < particle_count; ++i) {
    pos_x[i] = std::fmod(static_cast<double>((29U * i + 7U) % 997U) * 0.013, 1.0);
    pos_y[i] = std::fmod(static_cast<double>((47U * i + 11U) % 991U) * 0.017, 1.0);
    pos_z[i] = std::fmod(static_cast<double>((71U * i + 3U) % 983U) * 0.019, 1.0);
    mass[i] = 0.5 + static_cast<double>((13U * i) % 9U) * 0.05;
  }

  std::vector<std::uint32_t> active(16);
  for (std::size_t i = 0; i < active.size(); ++i) {
    active[i] = static_cast<std::uint32_t>(i * 3U);
  }

  cosmosim::gravity::TreeGravityOptions tight_options;
  tight_options.opening_theta = 0.4;
  tight_options.gravitational_constant_code = 1.0;
  tight_options.max_leaf_size = 4;
  tight_options.softening.epsilon_comoving = 1.0e-3;

  cosmosim::gravity::TreeGravityProfile tight_profile;
  double tight_mean_relative_error = 0.0;
  const double tight_max_relative_error = solveAndMeasureMaxRelativeError(
      pos_x,
      pos_y,
      pos_z,
      mass,
      active,
      tight_options,
      &tight_profile,
      &tight_mean_relative_error);

  cosmosim::gravity::TreeGravityOptions loose_options = tight_options;
  loose_options.opening_theta = 0.8;

  cosmosim::gravity::TreeGravityProfile loose_profile;
  double loose_mean_relative_error = 0.0;
  const double loose_max_relative_error = solveAndMeasureMaxRelativeError(
      pos_x,
      pos_y,
      pos_z,
      mass,
      active,
      loose_options,
      &loose_profile,
      &loose_mean_relative_error);

  std::ostringstream diag;
  diag << "Tree-vs-direct validation failed. "
       << "Assumption: this test targets isolated/non-periodic tree traversal only (no minimum-image wrapping). "
       << "tight(theta=" << tight_options.opening_theta << "): max_rel=" << tight_max_relative_error
       << ", mean_rel=" << tight_mean_relative_error
       << "; loose(theta=" << loose_options.opening_theta << "): max_rel=" << loose_max_relative_error
       << ", mean_rel=" << loose_mean_relative_error
       << "; tight.accepted_nodes=" << tight_profile.accepted_nodes
       << ", tight.visited_nodes=" << tight_profile.visited_nodes;

  requireOrThrow(tight_max_relative_error < 0.03, diag.str());
  requireOrThrow(tight_mean_relative_error < 0.015, diag.str());
  requireOrThrow(loose_max_relative_error >= tight_max_relative_error, diag.str());
  requireOrThrow(tight_profile.build_ms >= 0.0, "Tree profile.build_ms must be non-negative");
  requireOrThrow(tight_profile.multipole_ms >= 0.0, "Tree profile.multipole_ms must be non-negative");
  requireOrThrow(tight_profile.traversal_ms >= 0.0, "Tree profile.traversal_ms must be non-negative");
  requireOrThrow(tight_profile.accepted_nodes > 0, "Tree profile.accepted_nodes must be positive");

  return 0;
}
