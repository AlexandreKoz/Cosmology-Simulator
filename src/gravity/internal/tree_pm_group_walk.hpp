#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <vector>

#include "cosmosim/gravity/tree_gravity.hpp"
#include "cosmosim/gravity/tree_pm_coupling.hpp"

namespace cosmosim::gravity::internal {

// Group walk for the TreePM short-range residual (adaptive relative_monopole
// policy). A group is the set of active targets inside one tree leaf. One
// tree walk per group builds a shared interaction list that is valid for every
// target in the group:
//
// * node pruning uses the minimum box-box distance between the node AABB and
//   the group AABB, so no node any target could reach inside r_cut is lost;
// * a node is accepted only if G M l^2 <= alpha a_min r_min^4 and
//   l < theta_max r_min, where r_min is the distance from the node COM to the
//   group AABB (a lower bound on every target's COM distance) and a_min the
//   smallest previous-acceleration magnitude in the group. Every target
//   therefore satisfies the per-target criterion with margin;
// * a node whose AABB touches the group AABB is opened (no target may lie
//   inside an accepted node), and under kContainNode a node must lie inside
//   r_cut of every group point.
//
// Accepted leaves are listed one entry per leaf with the leaf cell geometry
// (minimum-image center offset from the group center plus half width). Each
// target then culls leaves whose cell lies entirely beyond r_cut by the exact
// per-axis periodic distance max(0, |minimg(c - p)| - h), which bounds every
// source-target separation in that leaf from below. Surviving leaves that are
// contiguous in tree order are coalesced into one range evaluated with a
// branch-free pair loop suited to compiler vectorization. Culling only skips
// pairs the kernel would mask (r2 > r_cut^2); forces differ from the unculled
// list at roundoff level through summation order only.

struct TreePmGroupBox {
  std::array<double, 3> center{};
  std::array<double, 3> half{};
};

using TreePmGroupInteractionList = TreePmGroupWalkWorkerScratch;

inline void clearTreePmGroupInteractionList(TreePmGroupInteractionList& list) noexcept {
  list.multipole_nodes.clear();
  list.range_begin.clear();
  list.range_end.clear();
  list.leaf_offset_x.clear();
  list.leaf_offset_y.clear();
  list.leaf_offset_z.clear();
  list.leaf_half_size.clear();
}

// Leaf culling threshold. A leaf is rejected for a target only when its cell's
// squared distance exceeds r_cut^2 by this relative margin, so rounding in the
// cell-offset arithmetic can never drop a pair the pair kernel would keep
// (r2 <= r_cut^2). Pairs inside the margin are still masked exactly by the kernel.
inline constexpr double k_tree_pm_group_leaf_cull_slack = 1.0 + 1.0e-9;

struct TreePmGroupWalkParameters {
  double cutoff_radius2 = 0.0;
  double relative_tolerance_over_g = 0.0;  // alpha / G folded: test M l^2 <= this * a_min r^4
  double maximum_opening_angle2 = 0.0;
  bool contain_cutoff = true;
};

// Minimum-image delta for node/group centers. Tree node coordinates may live
// in an unwrapped frame, so the general nearest-image form is used here (one
// call per visited node, not per pair).
[[nodiscard]] inline double treePmGroupNearestImage(double delta, double length) noexcept {
  if (!(length > 0.0)) return delta;
  return delta - length * std::nearbyint(delta / length);
}

template <class Stack>
inline void buildTreePmGroupInteractionList(
    const TreeNodeSoa& nodes, const TreePmGroupBox& group, double group_min_acceleration,
    const std::array<double, 3>& box_lengths, const TreePmGroupWalkParameters& p,
    Stack& stack, TreePmGroupInteractionList& out, TreePmTraversalCounters& counters) {
  clearTreePmGroupInteractionList(out);
  if (nodes.size() == 0U) return;
  stack.clear();
  stack.push_back(0U);
  while (!stack.empty()) {
    const TreeLocalIndex node = stack.back();
    stack.pop_back();
    ++counters.visited_nodes;
    const double h = nodes.half_size_comoving[node];
    const std::array<double, 3> dc{
        treePmGroupNearestImage(nodes.center_x_comoving[node] - group.center[0], box_lengths[0]),
        treePmGroupNearestImage(nodes.center_y_comoving[node] - group.center[1], box_lengths[1]),
        treePmGroupNearestImage(nodes.center_z_comoving[node] - group.center[2], box_lengths[2])};
    double gap2 = 0.0, far2 = 0.0;
    for (std::size_t axis = 0; axis < 3U; ++axis) {
      const double a = std::abs(dc[axis]);
      const double gap = std::max(0.0, a - h - group.half[axis]);
      const double far = a + h + group.half[axis];
      gap2 += gap * gap;
      far2 += far * far;
    }
    if (gap2 > p.cutoff_radius2) {
      ++counters.cutoff_pruned_nodes;
      continue;
    }
    if (nodes.child_count[node] == 0U) {
      ++counters.accepted_leaves;
      ++counters.accepted_nodes;
      const TreeLocalIndex begin = nodes.particle_begin[node];
      const TreeLocalIndex end = begin + static_cast<TreeLocalIndex>(nodes.particle_count[node]);
      // One entry per leaf keeps the geometry needed for per-target culling;
      // touching ranges are coalesced per target after culling instead.
      out.range_begin.push_back(begin);
      out.range_end.push_back(end);
      out.leaf_offset_x.push_back(dc[0]);
      out.leaf_offset_y.push_back(dc[1]);
      out.leaf_offset_z.push_back(dc[2]);
      out.leaf_half_size.push_back(h);
      continue;
    }
    bool accept = gap2 > 0.0 && (!p.contain_cutoff || far2 <= p.cutoff_radius2);
    if (accept) {
      const std::array<double, 3> dm{
          treePmGroupNearestImage(nodes.com_x_comoving[node] - group.center[0], box_lengths[0]),
          treePmGroupNearestImage(nodes.com_y_comoving[node] - group.center[1], box_lengths[1]),
          treePmGroupNearestImage(nodes.com_z_comoving[node] - group.center[2], box_lengths[2])};
      double r_min2 = 0.0;
      for (std::size_t axis = 0; axis < 3U; ++axis) {
        const double gap = std::max(0.0, std::abs(dm[axis]) - group.half[axis]);
        r_min2 += gap * gap;
      }
      const double width2 = 4.0 * h * h;
      accept = r_min2 > 0.0 && width2 < p.maximum_opening_angle2 * r_min2 &&
          nodes.mass_code[node] * width2 <=
              p.relative_tolerance_over_g * group_min_acceleration * r_min2 * r_min2;
      if (!accept) ++counters.relative_mac_rejections;
    } else {
      ++counters.near_node_rejections;
    }
    if (accept) {
      ++counters.accepted_internal_multipoles;
      ++counters.accepted_nodes;
      out.multipole_nodes.push_back(node);
    } else {
      ++counters.opened_nodes;
      const std::size_t child_offset = static_cast<std::size_t>(node) * 8U;
      for (std::size_t octant = 8U; octant > 0U; --octant) {
        const TreeLocalIndex child = nodes.child_index[child_offset + octant - 1U];
        if (child != kInvalidTreeLocalIndex) stack.push_back(child);
      }
    }
  }
}

struct TreePmGroupPairKernel {
  const double* x = nullptr;
  const double* y = nullptr;
  const double* z = nullptr;
  const double* mass = nullptr;
  std::array<double, 3> length{};
  double epsilon2 = 0.0;
  double cutoff_radius2 = 0.0;
  double inverse_2a_squared = 0.0;   // t = r^2 / (2 r_s)^2
  double inverse_8a_cubed = 0.0;     // long-range coefficient scale
  double table_inverse_step = 0.0;
  double table_t_limit = 0.0;        // strictly below the last table knot
  const double* table_value = nullptr;
  const double* table_slope = nullptr;
  double gravitational_constant = 0.0;
};

// Accumulates the softened-minus-Gaussian residual of tree-order sources
// [begin,end) on one target. Requires bounded periodic deltas (every wrapped
// source and target lie in [0,L)), so one conditional shift per axis is the
// nearest image; pairs near half a box lie beyond r_cut < L/2 and are masked.
// Self and coincident pairs contribute exactly zero. nearest_pair_r2 is
// min-combined with the squared separation of every evaluated pair; the
// group path cannot tell an exactly coincident distinct source from the
// target itself, so r2 == 0 is excluded there (the per-target walk keeps it).
inline void accumulateTreePmGroupPairs(
    const TreePmGroupPairKernel& k, double px, double py, double pz,
    std::size_t begin, std::size_t end,
    double& ax, double& ay, double& az, std::uint64_t& pairs,
    double& nearest_pair_r2) noexcept {
  const double lx = k.length[0], ly = k.length[1], lz = k.length[2];
  const double hx = 0.5 * lx, hy = 0.5 * ly, hz = 0.5 * lz;
  double sx = 0.0, sy = 0.0, sz = 0.0;
  double nearest = nearest_pair_r2;
  std::uint64_t counted = 0U;
#if defined(_OPENMP)
#pragma omp simd reduction(+:sx,sy,sz,counted) reduction(min:nearest)
#endif
  for (std::size_t j = begin; j < end; ++j) {
    double dx = k.x[j] - px, dy = k.y[j] - py, dz = k.z[j] - pz;
    dx = dx > hx ? dx - lx : (dx < -hx ? dx + lx : dx);
    dy = dy > hy ? dy - ly : (dy < -hy ? dy + ly : dy);
    dz = dz > hz ? dz - lz : (dz < -hz ? dz + lz : dz);
    const double r2 = dx * dx + dy * dy + dz * dz;
    const bool inside = r2 > 0.0 && r2 <= k.cutoff_radius2;
    const double t = std::min(r2 * k.inverse_2a_squared, k.table_t_limit);
    const double cell = t * k.table_inverse_step;
    const auto i = static_cast<std::size_t>(cell);
    const double u = cell - static_cast<double>(i);
    const double u2 = u * u, u3 = u2 * u;
    const double coefficient_long =
        ((2.0 * u3 - 3.0 * u2 + 1.0) * k.table_value[i] + (u3 - 2.0 * u2 + u) * k.table_slope[i] +
         (-2.0 * u3 + 3.0 * u2) * k.table_value[i + 1U] + (u3 - u2) * k.table_slope[i + 1U]) *
        k.inverse_8a_cubed;
    const double d = r2 + k.epsilon2;
    const double softened = 1.0 / (d * std::sqrt(d));
    const double factor = inside ? (softened - coefficient_long) * k.gravitational_constant * k.mass[j] : 0.0;
    sx += factor * dx;
    sy += factor * dy;
    sz += factor * dz;
    counted += inside ? 1U : 0U;
    nearest = std::min(nearest, inside ? r2 : nearest);
  }
  nearest_pair_r2 = nearest;
  ax += sx;
  ay += sy;
  az += sz;
  pairs += counted;
}

struct TreePmGroupLeafCullCounts {
  std::uint64_t candidate_leaves = 0;
  std::uint64_t culled_leaves = 0;
  std::uint64_t candidate_particles = 0;
};

// Evaluates the list's source leaves on one target with per-leaf culling.
// target_offset is the target's minimum-image offset from the group center
// used to build the list (|component| <= L/2), so one conditional shift of
// (leaf_offset - target_offset) is the minimum image of (leaf center - target).
// No square root is taken; the comparison is on squared distances.
inline void accumulateTreePmGroupCulledLeaves(
    const TreePmGroupPairKernel& k, const TreePmGroupInteractionList& list,
    const std::array<double, 3>& target_offset, double px, double py, double pz,
    double& ax, double& ay, double& az, std::uint64_t& pairs, double& nearest_pair_r2,
    TreePmGroupLeafCullCounts& counts) noexcept {
  const double lx = k.length[0], ly = k.length[1], lz = k.length[2];
  const double hx = 0.5 * lx, hy = 0.5 * ly, hz = 0.5 * lz;
  const double cull_radius2 = k.cutoff_radius2 * k_tree_pm_group_leaf_cull_slack;
  const auto wrap = [](double d, double length, double half) noexcept {
    return d > half ? d - length : (d < -half ? d + length : d);
  };
  const std::size_t leaf_count = list.range_begin.size();
  std::size_t run_begin = 0U, run_end = 0U;
  bool run_open = false;
  std::uint64_t culled = 0U, candidates = 0U;
  for (std::size_t leaf = 0; leaf < leaf_count; ++leaf) {
    const double h = list.leaf_half_size[leaf];
    const double gx = std::max(0.0, std::abs(wrap(list.leaf_offset_x[leaf] - target_offset[0], lx, hx)) - h);
    const double gy = std::max(0.0, std::abs(wrap(list.leaf_offset_y[leaf] - target_offset[1], ly, hy)) - h);
    const double gz = std::max(0.0, std::abs(wrap(list.leaf_offset_z[leaf] - target_offset[2], lz, hz)) - h);
    if (gx * gx + gy * gy + gz * gz > cull_radius2) {
      ++culled;
      continue;
    }
    const std::size_t begin = list.range_begin[leaf];
    const std::size_t end = list.range_end[leaf];
    candidates += end - begin;
    if (run_open && run_end == begin) {
      run_end = end;
      continue;
    }
    if (run_open) {
      accumulateTreePmGroupPairs(k, px, py, pz, run_begin, run_end, ax, ay, az, pairs, nearest_pair_r2);
    }
    run_begin = begin;
    run_end = end;
    run_open = true;
  }
  if (run_open) {
    accumulateTreePmGroupPairs(k, px, py, pz, run_begin, run_end, ax, ay, az, pairs, nearest_pair_r2);
  }
  counts.candidate_leaves += leaf_count;
  counts.culled_leaves += culled;
  counts.candidate_particles += candidates;
}

}  // namespace cosmosim::gravity::internal
