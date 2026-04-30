#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <numeric>
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

void directSumAccelerationWithSofteningView(
    std::span<const double> pos_x,
    std::span<const double> pos_y,
    std::span<const double> pos_z,
    std::span<const double> mass,
    std::span<const std::uint32_t> active,
    const cosmosim::gravity::TreeGravityOptions& options,
    const cosmosim::gravity::TreeSofteningView& softening_view,
    std::span<double> ax,
    std::span<double> ay,
    std::span<double> az) {
  for (std::size_t i = 0; i < active.size(); ++i) {
    const std::uint32_t target = active[i];
    const double tx = pos_x[target];
    const double ty = pos_y[target];
    const double tz = pos_z[target];
    const double target_eps = cosmosim::gravity::resolveTargetSofteningEpsilon(i, options.softening, softening_view);
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
      const double source_eps = cosmosim::gravity::resolveSourceSofteningEpsilon(j, options.softening, softening_view);
      const double pair_eps = cosmosim::gravity::combineSofteningPairEpsilon(source_eps, target_eps);
      const double factor = options.gravitational_constant_code * mass[j] * cosmosim::gravity::softenedInvR3(r2, pair_eps);
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

struct Distribution {
  std::vector<double> pos_x;
  std::vector<double> pos_y;
  std::vector<double> pos_z;
  std::vector<double> mass;
  std::vector<std::uint32_t> active;
};

[[nodiscard]] Distribution makeQuasiUniformDistribution() {
  constexpr std::size_t particle_count = 96;
  Distribution distribution;
  distribution.pos_x.resize(particle_count, 0.0);
  distribution.pos_y.resize(particle_count, 0.0);
  distribution.pos_z.resize(particle_count, 0.0);
  distribution.mass.resize(particle_count, 0.0);
  for (std::size_t i = 0; i < particle_count; ++i) {
    distribution.pos_x[i] = std::fmod(static_cast<double>((29U * i + 7U) % 997U) * 0.013, 1.0);
    distribution.pos_y[i] = std::fmod(static_cast<double>((47U * i + 11U) % 991U) * 0.017, 1.0);
    distribution.pos_z[i] = std::fmod(static_cast<double>((71U * i + 3U) % 983U) * 0.019, 1.0);
    distribution.mass[i] = 0.5 + static_cast<double>((13U * i) % 9U) * 0.05;
  }
  distribution.active.resize(16);
  for (std::size_t i = 0; i < distribution.active.size(); ++i) {
    distribution.active[i] = static_cast<std::uint32_t>(i * 3U);
  }
  return distribution;
}

[[nodiscard]] Distribution makeClusteredDistribution() {
  constexpr std::size_t particle_count = 96;
  Distribution distribution;
  distribution.pos_x.resize(particle_count, 0.0);
  distribution.pos_y.resize(particle_count, 0.0);
  distribution.pos_z.resize(particle_count, 0.0);
  distribution.mass.resize(particle_count, 0.0);
  for (std::size_t i = 0; i < particle_count; ++i) {
    const bool first_cluster = i < particle_count / 2U;
    const double phase = static_cast<double>(i % 48U);
    const double center = first_cluster ? 0.28 : 0.74;
    distribution.pos_x[i] = center + 0.045 * std::sin(0.31 * phase + (first_cluster ? 0.3 : 0.9));
    distribution.pos_y[i] = center + 0.038 * std::cos(0.27 * phase + (first_cluster ? 0.2 : 1.1));
    distribution.pos_z[i] = center + 0.041 * std::sin(0.19 * phase + (first_cluster ? 0.7 : 1.4));
    distribution.mass[i] = 0.4 + 0.03 * static_cast<double>((7U * i + 3U) % 11U);
  }
  distribution.active.resize(24);
  for (std::size_t i = 0; i < distribution.active.size(); ++i) {
    distribution.active[i] = static_cast<std::uint32_t>((5U * i + 1U) % particle_count);
  }
  return distribution;
}

void testMixedSpeciesSofteningAgainstDirect() {
  Distribution dist = makeQuasiUniformDistribution();
  dist.active.resize(dist.pos_x.size());
  for (std::size_t i = 0; i < dist.active.size(); ++i) {
    dist.active[i] = static_cast<std::uint32_t>(i);
  }
  std::vector<std::uint32_t> species_tag(dist.pos_x.size(), 0U);
  for (std::size_t i = 0; i < species_tag.size(); ++i) {
    species_tag[i] = static_cast<std::uint32_t>(i % 3U);
  }
  std::vector<double> source_override(species_tag.size(), 0.0);
  std::vector<std::uint8_t> source_override_mask(species_tag.size(), 0U);
  source_override[5] = 0.03;
  source_override_mask[5] = 1U;
  source_override[14] = 0.05;
  source_override_mask[14] = 1U;
  std::vector<double> target_override(dist.active.size(), 0.0);
  std::vector<std::uint8_t> target_override_mask(dist.active.size(), 0U);
  for (std::size_t i = 0; i < target_override.size(); ++i) {
    if ((i % 7U) == 0U) {
      target_override[i] = 0.04;
      target_override_mask[i] = 1U;
    }
  }

  cosmosim::gravity::TreeGravityOptions options;
  options.opening_theta = 0.35;
  options.multipole_order = cosmosim::gravity::TreeMultipoleOrder::kQuadrupole;
  options.gravitational_constant_code = 1.0;
  options.max_leaf_size = 1;
  options.softening.epsilon_comoving = 1.0e-3;
  cosmosim::gravity::TreeSofteningView softening_view{
      .source_species_tag = species_tag,
      .source_particle_epsilon_comoving = source_override,
      .source_particle_epsilon_override_mask = source_override_mask,
      .target_particle_epsilon_comoving = target_override,
      .target_particle_epsilon_override_mask = target_override_mask,
      .species_policy = {.epsilon_comoving_by_species = {0.001, 0.01, 0.02, 0.0, 0.0}, .enabled = true},
  };

  cosmosim::gravity::TreeGravitySolver solver;
  solver.build(dist.pos_x, dist.pos_y, dist.pos_z, dist.mass, options, nullptr, softening_view);
  std::vector<double> tree_ax(dist.active.size(), 0.0);
  std::vector<double> tree_ay(dist.active.size(), 0.0);
  std::vector<double> tree_az(dist.active.size(), 0.0);
  solver.evaluateActiveSet(
      dist.pos_x,
      dist.pos_y,
      dist.pos_z,
      dist.mass,
      dist.active,
      tree_ax,
      tree_ay,
      tree_az,
      options,
      nullptr,
      softening_view);

  std::vector<double> ref_ax(dist.active.size(), 0.0);
  std::vector<double> ref_ay(dist.active.size(), 0.0);
  std::vector<double> ref_az(dist.active.size(), 0.0);
  directSumAccelerationWithSofteningView(
      dist.pos_x, dist.pos_y, dist.pos_z, dist.mass, dist.active, options, softening_view, ref_ax, ref_ay, ref_az);

  double diff2 = 0.0;
  double ref2 = 0.0;
  for (std::size_t i = 0; i < tree_ax.size(); ++i) {
    const double dx = tree_ax[i] - ref_ax[i];
    const double dy = tree_ay[i] - ref_ay[i];
    const double dz = tree_az[i] - ref_az[i];
    diff2 += dx * dx + dy * dy + dz * dz;
    ref2 += ref_ax[i] * ref_ax[i] + ref_ay[i] * ref_ay[i] + ref_az[i] * ref_az[i];
  }
  const double rel_l2 = std::sqrt(diff2 / std::max(ref2, 1.0e-30));
  requireOrThrow(rel_l2 < 0.02, "mixed-species/per-particle softening tree vs direct mismatch");
}

}  // namespace

int main() {
  testMixedSpeciesSofteningAgainstDirect();
  const std::vector<Distribution> distributions = {makeQuasiUniformDistribution(), makeClusteredDistribution()};
  for (std::size_t dist_i = 0; dist_i < distributions.size(); ++dist_i) {
    const Distribution& dist = distributions[dist_i];
    cosmosim::gravity::TreeGravityOptions options;
    options.opening_theta = 0.6;
    options.opening_criterion = cosmosim::gravity::TreeOpeningCriterion::kBarnesHutComDistance;
    options.multipole_order = cosmosim::gravity::TreeMultipoleOrder::kQuadrupole;
    options.gravitational_constant_code = 1.0;
    options.max_leaf_size = 4;
    options.softening.epsilon_comoving = 1.0e-3;

    const std::array<double, 3> theta_schedule = {0.35, 0.6, 0.9};
    std::array<double, 3> max_error{};
    std::array<double, 3> mean_error{};
    std::array<std::uint64_t, 3> visited_nodes{};
    std::array<std::uint64_t, 3> pp_interactions{};

    for (std::size_t i = 0; i < theta_schedule.size(); ++i) {
      options.opening_theta = theta_schedule[i];
      cosmosim::gravity::TreeGravityProfile profile;
      max_error[i] = solveAndMeasureMaxRelativeError(
          dist.pos_x,
          dist.pos_y,
          dist.pos_z,
          dist.mass,
          dist.active,
          options,
          &profile,
          &mean_error[i]);
      visited_nodes[i] = profile.visited_nodes;
      pp_interactions[i] = profile.particle_particle_interactions;
      requireOrThrow(profile.build_ms >= 0.0, "Tree profile.build_ms must be non-negative");
      requireOrThrow(profile.multipole_ms >= 0.0, "Tree profile.multipole_ms must be non-negative");
      requireOrThrow(profile.traversal_ms >= 0.0, "Tree profile.traversal_ms must be non-negative");
      requireOrThrow(profile.accepted_nodes > 0, "Tree profile.accepted_nodes must be positive");
    }

    std::ostringstream diag;
    diag << "Tree-vs-direct validation failed for distribution=" << dist_i
         << " (isolated/non-periodic tree traversal only): "
         << "theta=[" << theta_schedule[0] << "," << theta_schedule[1] << "," << theta_schedule[2] << "]"
         << " max_error=[" << max_error[0] << "," << max_error[1] << "," << max_error[2] << "]"
         << " mean_error=[" << mean_error[0] << "," << mean_error[1] << "," << mean_error[2] << "]"
         << " visited_nodes=[" << visited_nodes[0] << "," << visited_nodes[1] << "," << visited_nodes[2] << "]"
         << " pp_interactions=[" << pp_interactions[0] << "," << pp_interactions[1] << "," << pp_interactions[2] << "]";

    requireOrThrow(max_error[0] <= max_error[1] + 1.0e-6, diag.str());
    requireOrThrow(max_error[1] <= max_error[2] + 1.0e-6, diag.str());
    requireOrThrow(mean_error[0] <= mean_error[1] + 1.0e-6, diag.str());
    requireOrThrow(mean_error[1] <= mean_error[2] + 1.0e-6, diag.str());
    requireOrThrow(visited_nodes[0] >= visited_nodes[1], diag.str());
    requireOrThrow(visited_nodes[1] >= visited_nodes[2], diag.str());
    requireOrThrow(pp_interactions[0] >= pp_interactions[1], diag.str());
    requireOrThrow(pp_interactions[1] >= pp_interactions[2], diag.str());
    requireOrThrow(max_error[0] < 0.02, diag.str());
    requireOrThrow(mean_error[0] < 0.01, diag.str());

    cosmosim::gravity::TreeGravityOptions geometric_options = options;
    geometric_options.opening_theta = 0.6;
    geometric_options.opening_criterion = cosmosim::gravity::TreeOpeningCriterion::kBarnesHutGeometric;
    cosmosim::gravity::TreeGravityProfile geometric_profile;
    double geometric_mean_error = 0.0;
    const double geometric_max_error = solveAndMeasureMaxRelativeError(
        dist.pos_x,
        dist.pos_y,
        dist.pos_z,
        dist.mass,
        dist.active,
        geometric_options,
        &geometric_profile,
        &geometric_mean_error);
    requireOrThrow(max_error[1] <= 0.95 * geometric_max_error + 1.0e-6, diag.str());
    requireOrThrow(mean_error[1] <= 0.95 * geometric_mean_error + 1.0e-6, diag.str());
    requireOrThrow(visited_nodes[1] >= geometric_profile.visited_nodes, diag.str());
  }

  return 0;
}
