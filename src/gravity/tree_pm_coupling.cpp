#include "cosmosim/gravity/tree_pm_coupling.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace cosmosim::gravity {
namespace {

[[nodiscard]] double minimumImageDelta(double delta, double box_size_comoving) {
  if (box_size_comoving <= 0.0) {
    return delta;
  }
  return delta - box_size_comoving * std::nearbyint(delta / box_size_comoving);
}

[[nodiscard]] bool forceAccumulatorShapeValid(const TreePmForceAccumulatorView& accumulator) {
  return accumulator.active_particle_index.size() == accumulator.accel_x_comoving.size() &&
      accumulator.active_particle_index.size() == accumulator.accel_y_comoving.size() &&
      accumulator.active_particle_index.size() == accumulator.accel_z_comoving.size();
}

void resizeCompactSidecars(std::vector<double>& first, std::vector<double>& second, std::vector<double>& third, std::size_t size) {
  first.assign(size, 0.0);
  second.assign(size, 0.0);
  third.assign(size, 0.0);
}

[[nodiscard]] double minimumDistanceToNodeAabb(
    double px,
    double py,
    double pz,
    double cx,
    double cy,
    double cz,
    double half_size_comoving,
    double box_size_comoving) {
  const double dx_abs = std::abs(minimumImageDelta(cx - px, box_size_comoving));
  const double dy_abs = std::abs(minimumImageDelta(cy - py, box_size_comoving));
  const double dz_abs = std::abs(minimumImageDelta(cz - pz, box_size_comoving));
  const double ex = std::max(0.0, dx_abs - half_size_comoving);
  const double ey = std::max(0.0, dy_abs - half_size_comoving);
  const double ez = std::max(0.0, dz_abs - half_size_comoving);
  return std::sqrt(ex * ex + ey * ey + ez * ez);
}

[[nodiscard]] double maximumDistanceToNodeAabb(
    double px,
    double py,
    double pz,
    double cx,
    double cy,
    double cz,
    double half_size_comoving,
    double box_size_comoving) {
  const double dx_abs = std::abs(minimumImageDelta(cx - px, box_size_comoving));
  const double dy_abs = std::abs(minimumImageDelta(cy - py, box_size_comoving));
  const double dz_abs = std::abs(minimumImageDelta(cz - pz, box_size_comoving));
  const double max_x = dx_abs + half_size_comoving;
  const double max_y = dy_abs + half_size_comoving;
  const double max_z = dz_abs + half_size_comoving;
  return std::sqrt(max_x * max_x + max_y * max_y + max_z * max_z);
}

}  // namespace

void TreePmForceAccumulatorView::reset() const {
  std::fill(accel_x_comoving.begin(), accel_x_comoving.end(), 0.0);
  std::fill(accel_y_comoving.begin(), accel_y_comoving.end(), 0.0);
  std::fill(accel_z_comoving.begin(), accel_z_comoving.end(), 0.0);
}

void TreePmForceAccumulatorView::addToActiveSlot(
    std::size_t active_slot,
    double ax_comoving,
    double ay_comoving,
    double az_comoving) const {
  accel_x_comoving[active_slot] += ax_comoving;
  accel_y_comoving[active_slot] += ay_comoving;
  accel_z_comoving[active_slot] += az_comoving;
}

TreePmCoordinator::TreePmCoordinator(PmGridShape pm_shape)
    : m_shape(pm_shape), m_grid(pm_shape), m_pm_solver(pm_shape), m_tree_solver() {}

void TreePmCoordinator::solveActiveSet(
    std::span<const double> pos_x_comoving,
    std::span<const double> pos_y_comoving,
    std::span<const double> pos_z_comoving,
    std::span<const double> mass_code,
    const TreePmForceAccumulatorView& accumulator,
    const TreePmOptions& options,
    TreePmProfileEvent* profile,
    TreePmDiagnostics* diagnostics) {
  if (!forceAccumulatorShapeValid(accumulator)) {
    throw std::invalid_argument("TreePM force accumulator spans must have matching active-set extent");
  }
  if (pos_x_comoving.size() != pos_y_comoving.size() || pos_x_comoving.size() != pos_z_comoving.size() ||
      pos_x_comoving.size() != mass_code.size()) {
    throw std::invalid_argument("TreePM requires position and mass spans with equal extent");
  }
  validateTreePmSplitPolicy(options.split_policy);

  const auto start = std::chrono::steady_clock::now();
  accumulator.reset();

  // PM owns long-range force via explicit Gaussian Fourier filter.
  PmSolveOptions pm_options = options.pm_options;
  pm_options.tree_pm_split_scale_comoving = options.split_policy.split_scale_comoving;
  m_pm_solver.assignDensity(m_grid, pos_x_comoving, pos_y_comoving, pos_z_comoving, mass_code, pm_options,
      profile != nullptr ? &profile->pm_profile : nullptr);
  m_pm_solver.solvePoissonPeriodic(m_grid, pm_options, profile != nullptr ? &profile->pm_profile : nullptr);

  const std::size_t active_count = accumulator.active_particle_index.size();
  resizeCompactSidecars(m_active_pos_x_comoving, m_active_pos_y_comoving, m_active_pos_z_comoving, active_count);
  resizeCompactSidecars(m_active_pm_ax_comoving, m_active_pm_ay_comoving, m_active_pm_az_comoving, active_count);

  for (std::size_t i = 0; i < accumulator.active_particle_index.size(); ++i) {
    const std::uint32_t particle_index = accumulator.active_particle_index[i];
    if (particle_index >= pos_x_comoving.size()) {
      throw std::out_of_range("TreePM active index exceeds particle count");
    }
    m_active_pos_x_comoving[i] = pos_x_comoving[particle_index];
    m_active_pos_y_comoving[i] = pos_y_comoving[particle_index];
    m_active_pos_z_comoving[i] = pos_z_comoving[particle_index];
  }

  m_pm_solver.interpolateForces(
      m_grid,
      m_active_pos_x_comoving,
      m_active_pos_y_comoving,
      m_active_pos_z_comoving,
      m_active_pm_ax_comoving,
      m_active_pm_ay_comoving,
      m_active_pm_az_comoving,
      pm_options,
      profile != nullptr ? &profile->pm_profile : nullptr);

  for (std::size_t i = 0; i < accumulator.active_particle_index.size(); ++i) {
    accumulator.addToActiveSlot(i, m_active_pm_ax_comoving[i], m_active_pm_ay_comoving[i], m_active_pm_az_comoving[i]);
  }

  // Tree owns short-range residual with the complementary real-space kernel.
  const auto tree_start = std::chrono::steady_clock::now();
  m_tree_solver.build(pos_x_comoving, pos_y_comoving, pos_z_comoving, mass_code, options.tree_options,
      profile != nullptr ? &profile->tree_profile : nullptr);
  evaluateShortRangeResidual(
      pos_x_comoving,
      pos_y_comoving,
      pos_z_comoving,
      mass_code,
      accumulator,
      options,
      profile != nullptr ? &profile->tree_profile : nullptr);
  const auto tree_stop = std::chrono::steady_clock::now();

  if (diagnostics != nullptr) {
    *diagnostics = computeTreePmDiagnostics(options.split_policy);
    diagnostics->residual_pruned_nodes = m_last_residual_stats.pruned_nodes;
    diagnostics->residual_pair_skips_cutoff = m_last_residual_stats.pair_skips_cutoff;
    diagnostics->residual_pair_evaluations = m_last_residual_stats.pair_evaluations;
  }

  if (profile != nullptr) {
    profile->tree_short_range_ms += std::chrono::duration<double, std::milli>(tree_stop - tree_start).count();
    profile->coupling_overhead_ms += std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start).count();
  }
}

void TreePmCoordinator::evaluateShortRangeResidual(
    std::span<const double> pos_x_comoving,
    std::span<const double> pos_y_comoving,
    std::span<const double> pos_z_comoving,
    std::span<const double> mass_code,
    const TreePmForceAccumulatorView& accumulator,
    const TreePmOptions& options,
    TreeGravityProfile* tree_profile) {
  std::uint64_t visited_nodes = 0;
  std::uint64_t accepted_nodes = 0;
  std::uint64_t pp_interactions = 0;
  std::uint64_t cutoff_pruned_nodes = 0;
  std::uint64_t cutoff_pair_skips = 0;

  const auto traversal_start = std::chrono::steady_clock::now();
  const TreeNodeSoa& nodes = m_tree_solver.nodes();
  const TreeMortonOrdering& ordering = m_tree_solver.ordering();
  const double box_size_comoving = options.pm_options.box_size_mpc_comoving;
  const double cutoff_radius_comoving = options.split_policy.cutoff_radius_comoving;
  const double cutoff_radius2_comoving = cutoff_radius_comoving * cutoff_radius_comoving;
  std::vector<std::uint32_t> stack;
  stack.reserve(256);

  for (std::size_t active_i = 0; active_i < accumulator.active_particle_index.size(); ++active_i) {
    const std::uint32_t particle_index = accumulator.active_particle_index[active_i];
    const double px = pos_x_comoving[particle_index];
    const double py = pos_y_comoving[particle_index];
    const double pz = pos_z_comoving[particle_index];

    double ax = 0.0;
    double ay = 0.0;
    double az = 0.0;

    stack.clear();
    stack.push_back(0U);

    while (!stack.empty()) {
      const std::uint32_t node_index = stack.back();
      stack.pop_back();
      ++visited_nodes;

      const double half_size = nodes.half_size_comoving[node_index];
      const double min_node_distance = minimumDistanceToNodeAabb(
          px,
          py,
          pz,
          nodes.center_x_comoving[node_index],
          nodes.center_y_comoving[node_index],
          nodes.center_z_comoving[node_index],
          half_size,
          box_size_comoving);
      if (min_node_distance > cutoff_radius_comoving) {
        ++cutoff_pruned_nodes;
        continue;
      }

      const double dx = minimumImageDelta(nodes.com_x_comoving[node_index] - px, box_size_comoving);
      const double dy = minimumImageDelta(nodes.com_y_comoving[node_index] - py, box_size_comoving);
      const double dz = minimumImageDelta(nodes.com_z_comoving[node_index] - pz, box_size_comoving);
      const double r2 = dx * dx + dy * dy + dz * dz;
      const double r = std::sqrt(std::max(r2, 1.0e-30));

      const double l_over_r = (2.0 * half_size) / r;
      const bool is_leaf = nodes.child_count[node_index] == 0;
      const bool geometric_accept = is_leaf || (l_over_r < options.tree_options.opening_theta);
      const bool node_within_cutoff = is_leaf || (maximumDistanceToNodeAabb(
          px,
          py,
          pz,
          nodes.center_x_comoving[node_index],
          nodes.center_y_comoving[node_index],
          nodes.center_z_comoving[node_index],
          half_size,
          box_size_comoving) <= cutoff_radius_comoving);
      const bool accept = geometric_accept && node_within_cutoff;

      if (accept) {
        ++accepted_nodes;
        if (is_leaf) {
          const std::uint32_t begin = nodes.particle_begin[node_index];
          const std::uint32_t end = begin + nodes.particle_count[node_index];
          for (std::uint32_t sorted_i = begin; sorted_i < end; ++sorted_i) {
            const std::uint32_t source_index = ordering.sorted_particle_index[sorted_i];
            if (source_index == particle_index) {
              continue;
            }
            const double sx = minimumImageDelta(pos_x_comoving[source_index] - px, box_size_comoving);
            const double sy = minimumImageDelta(pos_y_comoving[source_index] - py, box_size_comoving);
            const double sz = minimumImageDelta(pos_z_comoving[source_index] - pz, box_size_comoving);
            const double sr2 = sx * sx + sy * sy + sz * sz;
            if (sr2 > cutoff_radius2_comoving) {
              ++cutoff_pair_skips;
              continue;
            }
            const double sr = std::sqrt(std::max(sr2, 1.0e-30));
            const double split_factor = treePmGaussianShortRangeForceFactor(sr, options.split_policy.split_scale_comoving);
            // Contract: short-range residual is the softened tree force multiplied by the
            // Gaussian real-space residual factor so that tree+PM composes to the unsplit
            // softened force before explicit r_cut truncation.
            const double softened_factor =
                softenedInvR3(sr2, options.tree_options.softening) * split_factor * options.tree_options.gravitational_constant_code;
            ax += softened_factor * mass_code[source_index] * sx;
            ay += softened_factor * mass_code[source_index] * sy;
            az += softened_factor * mass_code[source_index] * sz;
            ++pp_interactions;
          }
        } else {
          const double split_factor = treePmGaussianShortRangeForceFactor(r, options.split_policy.split_scale_comoving);
          // Same softened-residual contract as the leaf pair path, applied to accepted nodes.
          const double softened_factor =
              softenedInvR3(r2, options.tree_options.softening) * split_factor * options.tree_options.gravitational_constant_code;
          ax += softened_factor * nodes.mass_code[node_index] * dx;
          ay += softened_factor * nodes.mass_code[node_index] * dy;
          az += softened_factor * nodes.mass_code[node_index] * dz;
        }
      } else {
        const std::size_t child_offset = static_cast<std::size_t>(node_index) * 8U;
        for (std::uint8_t octant = 0; octant < 8U; ++octant) {
          const std::uint32_t child = nodes.child_index[child_offset + octant];
          if (child != std::numeric_limits<std::uint32_t>::max()) {
            stack.push_back(child);
          }
        }
      }
    }

    accumulator.addToActiveSlot(active_i, ax, ay, az);
  }

  if (tree_profile != nullptr) {
    tree_profile->visited_nodes += visited_nodes;
    tree_profile->accepted_nodes += accepted_nodes;
    tree_profile->particle_particle_interactions += pp_interactions;
    tree_profile->traversal_ms +=
        std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - traversal_start).count();
    if (!accumulator.active_particle_index.empty()) {
      tree_profile->average_interactions_per_target = static_cast<double>(tree_profile->particle_particle_interactions) /
          static_cast<double>(accumulator.active_particle_index.size());
    }
  }
  m_last_residual_stats.pruned_nodes = cutoff_pruned_nodes;
  m_last_residual_stats.pair_skips_cutoff = cutoff_pair_skips;
  m_last_residual_stats.pair_evaluations = pp_interactions;
}

TreePmDiagnostics computeTreePmDiagnostics(const TreePmSplitPolicy& split_policy) {
  validateTreePmSplitPolicy(split_policy);

  TreePmDiagnostics diagnostics;
  diagnostics.mesh_spacing_comoving = split_policy.mesh_spacing_comoving;
  diagnostics.asmth_cells = split_policy.asmth_cells;
  diagnostics.rcut_cells = split_policy.rcut_cells;
  diagnostics.split_scale_comoving = split_policy.split_scale_comoving;
  diagnostics.cutoff_radius_comoving = split_policy.cutoff_radius_comoving;
  diagnostics.short_range_factor_at_split =
      treePmGaussianShortRangeForceFactor(split_policy.split_scale_comoving, split_policy.split_scale_comoving);
  diagnostics.long_range_factor_at_split =
      treePmGaussianLongRangeForceFactor(split_policy.split_scale_comoving, split_policy.split_scale_comoving);
  diagnostics.short_range_factor_at_cutoff =
      treePmGaussianShortRangeForceFactor(split_policy.cutoff_radius_comoving, split_policy.split_scale_comoving);
  diagnostics.long_range_factor_at_cutoff =
      treePmGaussianLongRangeForceFactor(split_policy.cutoff_radius_comoving, split_policy.split_scale_comoving);
  diagnostics.composition_error_at_split = std::abs(
      diagnostics.short_range_factor_at_split + diagnostics.long_range_factor_at_split - 1.0);

  const double radii[] = {
      0.25 * split_policy.split_scale_comoving,
      0.5 * split_policy.split_scale_comoving,
      split_policy.split_scale_comoving,
      2.0 * split_policy.split_scale_comoving,
      4.0 * split_policy.split_scale_comoving,
  };
  for (const double radius_comoving : radii) {
    const double composed = treePmGaussianShortRangeForceFactor(radius_comoving, split_policy.split_scale_comoving) +
        treePmGaussianLongRangeForceFactor(radius_comoving, split_policy.split_scale_comoving);
    diagnostics.max_relative_composition_error = std::max(
        diagnostics.max_relative_composition_error,
        std::abs(composed - 1.0) / std::max(1.0, std::abs(composed)));
  }

  return diagnostics;
}

}  // namespace cosmosim::gravity
