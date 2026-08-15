#include "cosmosim/gravity/tree_gravity.hpp"
#include "cosmosim/core/build_config.hpp"

#include <algorithm>
#include <bit>
#include <array>
#include <chrono>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// shared hot tree interaction invariants
#include "internal/tree_interaction_common.hpp"

namespace cosmosim::gravity {
namespace {

struct OctantSpan {
  TreeLocalCount begin = 0;
  TreeLocalCount end = 0;
};

[[nodiscard]] std::uint64_t sourceFingerprint(
    std::span<const double> pos_x_comoving,
    std::span<const double> pos_y_comoving,
    std::span<const double> pos_z_comoving,
    std::span<const double> mass_code) noexcept {
  // FNV-1a over the exact source bits. This is intentionally only the legacy
  // safety path; production workflow callers should provide a non-zero source
  // generation and avoid an O(N) validation scan on every traversal.
  std::uint64_t hash = 1469598103934665603ULL;
  const auto mix = [&hash](std::uint64_t value) {
    hash ^= value;
    hash *= 1099511628211ULL;
  };
  mix(static_cast<std::uint64_t>(pos_x_comoving.size()));
  for (std::size_t i = 0; i < pos_x_comoving.size(); ++i) {
    mix(std::bit_cast<std::uint64_t>(pos_x_comoving[i]));
    mix(std::bit_cast<std::uint64_t>(pos_y_comoving[i]));
    mix(std::bit_cast<std::uint64_t>(pos_z_comoving[i]));
    mix(std::bit_cast<std::uint64_t>(mass_code[i]));
  }
  return hash;
}

[[nodiscard]] std::uint8_t octantForParticle(
    double x_comoving,
    double y_comoving,
    double z_comoving,
    double center_x_comoving,
    double center_y_comoving,
    double center_z_comoving) {
  const std::uint8_t ox = x_comoving >= center_x_comoving ? 1U : 0U;
  const std::uint8_t oy = y_comoving >= center_y_comoving ? 1U : 0U;
  const std::uint8_t oz = z_comoving >= center_z_comoving ? 1U : 0U;
  return static_cast<std::uint8_t>((ox << 2U) | (oy << 1U) | oz);
}

void validateInputSpans(
    std::span<const double> pos_x_comoving,
    std::span<const double> pos_y_comoving,
    std::span<const double> pos_z_comoving,
    std::span<const double> mass_code) {
  if (pos_x_comoving.size() != pos_y_comoving.size() || pos_x_comoving.size() != pos_z_comoving.size() ||
      pos_x_comoving.size() != mass_code.size()) {
    throw std::invalid_argument("Tree gravity requires equal particle span lengths");
  }
  for (std::size_t i = 0; i < mass_code.size(); ++i) {
    if (!std::isfinite(pos_x_comoving[i]) || !std::isfinite(pos_y_comoving[i]) ||
        !std::isfinite(pos_z_comoving[i]) || !std::isfinite(mass_code[i]) ||
        mass_code[i] < 0.0) {
      throw std::invalid_argument(
          "Tree gravity requires finite coordinates and finite non-negative masses");
    }
  }
}

void validateOptions(const TreeGravityOptions& options) {
  if (!std::isfinite(options.opening_theta) || options.opening_theta <= 0.0 ||
      !std::isfinite(options.relative_force_tolerance) || options.relative_force_tolerance <= 0.0 ||
      !std::isfinite(options.relative_force_acceleration_floor_code) ||
      options.relative_force_acceleration_floor_code <= 0.0 ||
      !std::isfinite(options.gravitational_constant_code) || options.gravitational_constant_code <= 0.0 ||
      !std::isfinite(options.softening.epsilon_comoving) || options.softening.epsilon_comoving < 0.0 ||
      options.softening.kernel != TreeSofteningKernel::kPlummer ||
      options.max_leaf_size == 0) {
    throw std::invalid_argument("Invalid tree gravity options");
  }

  switch (options.opening_criterion) {
    case TreeOpeningCriterion::kBarnesHutGeometric:
    case TreeOpeningCriterion::kBarnesHutComDistance:
    case TreeOpeningCriterion::kRelativeForceError:
      break;
    default:
      throw std::invalid_argument("Invalid tree gravity opening criterion");
  }
  switch (options.multipole_order) {
    case TreeMultipoleOrder::kMonopole:
    case TreeMultipoleOrder::kQuadrupole:
      break;
    default:
      throw std::invalid_argument("Invalid tree gravity multipole order");
  }
}

[[nodiscard]] std::array<double, 3> monopolePlusQuadrupoleAccel(
    const TreeNodeSoa& nodes,
    TreeLocalIndex node_index,
    double dx,
    double dy,
    double dz,
    double target_softening_comoving,
    const TreeGravityOptions& options) {
  const double r2 = dx * dx + dy * dy + dz * dz;
  const double pair_epsilon =
      combineSofteningPairEpsilonUnchecked(nodes.softening_max_comoving[node_index], target_softening_comoving);
  const double eps2 = pair_epsilon * pair_epsilon;
  const double denom = std::max(r2 + eps2, 1.0e-30);

  const double inv_r3 = 1.0 / (denom * std::sqrt(denom));

  const double gm = options.gravitational_constant_code * nodes.mass_code[node_index];
  double ax = gm * inv_r3 * dx;
  double ay = gm * inv_r3 * dy;
  double az = gm * inv_r3 * dz;

  if (options.multipole_order != TreeMultipoleOrder::kQuadrupole) {
    return {ax, ay, az};
  }

  const double qxx = nodes.quad_xx[node_index];
  const double qxy = nodes.quad_xy[node_index];
  const double qxz = nodes.quad_xz[node_index];
  const double qyy = nodes.quad_yy[node_index];
  const double qyz = nodes.quad_yz[node_index];
  const double qzz = nodes.quad_zz[node_index];


  // Expand g_i(d) = d_i (|d|^2 + eps^2)^(-3/2) through second
  // order about the node COM.  Q alone is sufficient only for a harmonic
  // unsoftened kernel; the trace lane supplies the isotropic second moment
  // needed by the softened kernel.
  const double moment_trace = nodes.second_moment_trace[node_index];
  const double mxx = (qxx + moment_trace) / 3.0;
  const double mxy = qxy / 3.0;
  const double mxz = qxz / 3.0;
  const double myy = (qyy + moment_trace) / 3.0;
  const double myz = qyz / 3.0;
  const double mzz = (qzz + moment_trace) / 3.0;
  const double mdx = mxx * dx + mxy * dy + mxz * dz;
  const double mdy = mxy * dx + myy * dy + myz * dz;
  const double mdz = mxz * dx + myz * dy + mzz * dz;
  const double dmd = dx * mdx + dy * mdy + dz * mdz;
  const double r = std::sqrt(std::max(r2, 1.0e-30));
  const double radial_first = -3.0 * r / (denom * denom * std::sqrt(denom));
  const double radial_second =
      -3.0 / (denom * denom * std::sqrt(denom)) +
      15.0 * r2 / (denom * denom * denom * std::sqrt(denom));
  const double first_over_r = radial_first / r;
  const double contracted_radial = radial_second / r2 - radial_first / (r2 * r);
  const double correction_scale = options.gravitational_constant_code;
  ax += correction_scale *
      (mdx * first_over_r + 0.5 * moment_trace * dx * first_over_r +
       0.5 * dx * dmd * contracted_radial);
  ay += correction_scale *
      (mdy * first_over_r + 0.5 * moment_trace * dy * first_over_r +
       0.5 * dy * dmd * contracted_radial);
  az += correction_scale *
      (mdz * first_over_r + 0.5 * moment_trace * dz * first_over_r +
       0.5 * dz * dmd * contracted_radial);
  return {ax, ay, az};
}

void pushChildrenNearFirst(
    const TreeNodeSoa& nodes,
    TreeLocalIndex node_index,
    double px,
    double py,
    double pz,
    std::vector<TreeLocalIndex>& stack) {
  std::array<std::pair<double, TreeLocalIndex>, 8> child_dist2{};
  std::size_t count = 0;
  const std::size_t child_offset = static_cast<std::size_t>(node_index) * 8U;
  for (std::uint8_t octant = 0; octant < 8U; ++octant) {
    const TreeLocalIndex child = nodes.child_index[child_offset + octant];
    if (child == kInvalidTreeLocalIndex) {
      continue;
    }
    const double dx = nodes.center_x_comoving[child] - px;
    const double dy = nodes.center_y_comoving[child] - py;
    const double dz = nodes.center_z_comoving[child] - pz;
    child_dist2[count++] = {dx * dx + dy * dy + dz * dz, child};
  }

  std::sort(child_dist2.begin(), child_dist2.begin() + static_cast<std::ptrdiff_t>(count),
      [](const auto& lhs, const auto& rhs) { return lhs.first < rhs.first; });
  for (std::size_t i = count; i > 0; --i) {
    stack.push_back(child_dist2[i - 1].second);
  }
}

}  // namespace

std::size_t TreeNodeSoa::size() const {
  return center_x_comoving.size();
}

void TreeNodeSoa::clear() {
  center_x_comoving.clear();
  center_y_comoving.clear();
  center_z_comoving.clear();
  half_size_comoving.clear();
  mass_code.clear();
  com_x_comoving.clear();
  com_y_comoving.clear();
  com_z_comoving.clear();
  quad_xx.clear();
  quad_xy.clear();
  quad_xz.clear();
  quad_yy.clear();
  quad_yz.clear();
  quad_zz.clear();
  second_moment_trace.clear();
  softening_min_comoving.clear();
  softening_max_comoving.clear();
  child_base.clear();
  child_count.clear();
  child_index.clear();
  particle_begin.clear();
  particle_count.clear();
}

void TreeNodeSoa::reserve(std::size_t count, bool include_quadrupoles) {
  center_x_comoving.reserve(count);
  center_y_comoving.reserve(count);
  center_z_comoving.reserve(count);
  half_size_comoving.reserve(count);
  mass_code.reserve(count);
  com_x_comoving.reserve(count);
  com_y_comoving.reserve(count);
  com_z_comoving.reserve(count);
  if (include_quadrupoles) {
    quad_xx.reserve(count);
    quad_xy.reserve(count);
    quad_xz.reserve(count);
    quad_yy.reserve(count);
    quad_yz.reserve(count);
    quad_zz.reserve(count);
    second_moment_trace.reserve(count);
  } else {
    // Quadrupoles are cold state for monopole builds. Release a previous
    // quadrupole build's retained capacity rather than charging every future
    // monopole traversal for lanes it cannot read.
    std::vector<double>{}.swap(quad_xx);
    std::vector<double>{}.swap(quad_xy);
    std::vector<double>{}.swap(quad_xz);
    std::vector<double>{}.swap(quad_yy);
    std::vector<double>{}.swap(quad_yz);
    std::vector<double>{}.swap(quad_zz);
    std::vector<double>{}.swap(second_moment_trace);
  }
  softening_min_comoving.reserve(count);
  softening_max_comoving.reserve(count);
  child_base.reserve(count);
  child_count.reserve(count);
  child_index.reserve(count * 8U);
  particle_begin.reserve(count);
  particle_count.reserve(count);
}


void TreeNodeSoa::appendMemoryReport(core::MemoryReportBuilder& builder) const {
  const auto add = [&builder](std::string label, const auto& container) {
    const std::uint64_t bytes = core::ownedCapacityBytesForContainer(container);
    builder.addEntry(core::MemoryEntry{.subsystem = core::MemorySubsystem::kTree,
                                       .lifetime = core::MemoryLifetime::kTransient,
                                       .label = std::move(label),
                                       .current_size_bytes = core::currentSizeBytesForContainer(container),
                                       .owned_capacity_bytes = bytes,
                                       .high_water_bytes = 0U,
                                       .estimated_next_step_bytes = 0U,
                                       .uncertainty_note = "historical lane high-water is not retained after optional lane release; solver reports node-capacity high-water separately"});
  };
  add("tree.nodes.center_x_comoving", center_x_comoving);
  add("tree.nodes.center_y_comoving", center_y_comoving);
  add("tree.nodes.center_z_comoving", center_z_comoving);
  add("tree.nodes.half_size_comoving", half_size_comoving);
  add("tree.nodes.mass_code", mass_code);
  add("tree.nodes.com_x_comoving", com_x_comoving);
  add("tree.nodes.com_y_comoving", com_y_comoving);
  add("tree.nodes.com_z_comoving", com_z_comoving);
  add("tree.nodes.quad_xx", quad_xx);
  add("tree.nodes.quad_xy", quad_xy);
  add("tree.nodes.quad_xz", quad_xz);
  add("tree.nodes.quad_yy", quad_yy);
  add("tree.nodes.quad_yz", quad_yz);
  add("tree.nodes.quad_zz", quad_zz);
  add("tree.nodes.second_moment_trace", second_moment_trace);
  add("tree.nodes.softening_min_comoving", softening_min_comoving);
  add("tree.nodes.softening_max_comoving", softening_max_comoving);
  add("tree.nodes.child_base", child_base);
  add("tree.nodes.child_count", child_count);
  add("tree.nodes.child_index", child_index);
  add("tree.nodes.particle_begin", particle_begin);
  add("tree.nodes.particle_count", particle_count);
}

void TreeGravitySolver::build(
    const TreeGravitySourceView& source_view,
    const TreeGravityOptions& options,
    TreeGravityProfile* profile,
    const TreeSofteningView& softening_view) {
  build(
      source_view.pos_x_comoving,
      source_view.pos_y_comoving,
      source_view.pos_z_comoving,
      source_view.mass_code,
      options,
      profile,
      softening_view,
      source_view.source_generation);
}

void TreeGravitySolver::build(
    std::span<const double> pos_x_comoving,
    std::span<const double> pos_y_comoving,
    std::span<const double> pos_z_comoving,
    std::span<const double> mass_code,
    const TreeGravityOptions& options,
    TreeGravityProfile* profile,
    const TreeSofteningView& softening_view,
    GravitySourceGeneration source_generation) {
  const auto build_start = std::chrono::steady_clock::now();
  validateInputSpans(pos_x_comoving, pos_y_comoving, pos_z_comoving, mass_code);
  validateOptions(options);
  if (pos_x_comoving.size() > static_cast<std::size_t>(kInvalidTreeLocalIndex)) {
    throw std::overflow_error("Tree gravity source count exceeds the 32-bit tree-index contract");
  }

  m_build_source_count = pos_x_comoving.size();
  m_build_source_generation = source_generation;
  m_build_source_fingerprint = !source_generation.valid()
      ? sourceFingerprint(pos_x_comoving, pos_y_comoving, pos_z_comoving, mass_code)
      : 0U;
  m_build_multipole_order = options.multipole_order;
  m_build_max_leaf_size = options.max_leaf_size;
  m_build_softening = options.softening;

  m_nodes.clear();
  m_source_softening_epsilon_comoving.clear();
  m_source_softening_epsilon_comoving.resize(pos_x_comoving.size(), options.softening.epsilon_comoving);
  if (!softening_view.source_particle_epsilon_comoving.empty() &&
      softening_view.source_particle_epsilon_comoving.size() != pos_x_comoving.size()) {
    throw std::invalid_argument("Tree gravity source softening sidecar size must match source particle count");
  }
  if (!softening_view.source_particle_epsilon_override_mask.empty() &&
      softening_view.source_particle_epsilon_override_mask.size() != pos_x_comoving.size()) {
    throw std::invalid_argument("Tree gravity source softening override mask size must match source particle count");
  }
  if (!softening_view.source_species_tag.empty() && softening_view.source_species_tag.size() != pos_x_comoving.size()) {
    throw std::invalid_argument("Tree gravity source species sidecar size must match source particle count");
  }
  for (std::size_t i = 0; i < pos_x_comoving.size(); ++i) {
    m_source_softening_epsilon_comoving[i] = resolveSourceSofteningEpsilon(i, options.softening, softening_view);
  }
  m_partition_scratch.assign(pos_x_comoving.size(), 0U);
  const auto ordering_start = std::chrono::steady_clock::now();
  m_ordering = buildMortonOrdering(pos_x_comoving, pos_y_comoving, pos_z_comoving);
  const auto ordering_stop = std::chrono::steady_clock::now();
  if (pos_x_comoving.empty()) {
    m_tree_build_generation = nextGravityIdentity(m_tree_build_generation, "Tree build generation overflow");
    if (profile != nullptr) {
      *profile = {};
      profile->build_count = 1U;
      profile->morton_ordering_ms =
          std::chrono::duration<double, std::milli>(ordering_stop - ordering_start).count();
      profile->build_ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - build_start).count();
    }
    return;
  }

  const TreeBounds bounds = computeTreeBounds(pos_x_comoving, pos_y_comoving, pos_z_comoving);
  const double max_extent = std::max(bounds.maxExtentComoving(), 1.0e-10);
  const double center_x_comoving = 0.5 * (bounds.min_x_comoving + bounds.max_x_comoving);
  const double center_y_comoving = 0.5 * (bounds.min_y_comoving + bounds.max_y_comoving);
  const double center_z_comoving = 0.5 * (bounds.min_z_comoving + bounds.max_z_comoving);

  // A full octree with non-trivial leaf occupancy needs O(N/leaf_size) nodes,
  // not 2*N nodes. Seed from the expected leaf count and the previous build's
  // actual capacity high-water; std::vector growth remains the correctness
  // fallback for adversarial clustered geometry.
  const std::size_t leaf_estimate =
      (pos_x_comoving.size() + options.max_leaf_size - 1U) / options.max_leaf_size;
  const std::size_t structural_estimate = leaf_estimate > ((std::numeric_limits<std::size_t>::max() - 1U) / 2U)
      ? pos_x_comoving.size()
      : 1U + 2U * leaf_estimate;
  std::size_t reserve_nodes = std::max<std::size_t>(1U, structural_estimate);
  if (m_node_capacity_high_water > 0U) {
    const std::size_t headroom = m_node_capacity_high_water / 8U + 8U;
    const std::size_t prior_with_headroom = m_node_capacity_high_water >
            std::numeric_limits<std::size_t>::max() - headroom
        ? m_node_capacity_high_water
        : m_node_capacity_high_water + headroom;
    reserve_nodes = std::max(reserve_nodes, prior_with_headroom);
  }
  m_nodes.reserve(reserve_nodes, options.multipole_order == TreeMultipoleOrder::kQuadrupole);
  const auto topology_start = std::chrono::steady_clock::now();
  const TreeLocalIndex root_index = buildNodeRecursive(
      pos_x_comoving,
      pos_y_comoving,
      pos_z_comoving,
      mass_code,
      0,
      checkedTreeLocalIndex(m_ordering.sorted_particle_index.size(), "tree ordering exceeds local index policy"),
      center_x_comoving,
      center_y_comoving,
      center_z_comoving,
      0.5 * max_extent * (1.0 + 1.0e-8),
      options);
  (void)root_index;
  const auto topology_stop = std::chrono::steady_clock::now();
  m_node_capacity_high_water = std::max(m_node_capacity_high_water, m_nodes.center_x_comoving.capacity());
  m_tree_build_generation = nextGravityIdentity(m_tree_build_generation, "Tree build generation overflow");

  const auto multipole_start = std::chrono::steady_clock::now();
  accumulateMultipoles(pos_x_comoving, pos_y_comoving, pos_z_comoving, mass_code, 0, options.multipole_order);
  const auto multipole_stop = std::chrono::steady_clock::now();

  if (profile != nullptr) {
    profile->build_count = 1U;
    profile->multipole_refresh_count = 1U;
    profile->build_ms = std::chrono::duration<double, std::milli>(multipole_start - build_start).count();
    profile->multipole_ms = std::chrono::duration<double, std::milli>(multipole_stop - multipole_start).count();
    profile->morton_ordering_ms =
        std::chrono::duration<double, std::milli>(ordering_stop - ordering_start).count();
    profile->topology_build_ms =
        std::chrono::duration<double, std::milli>(topology_stop - topology_start).count();
    profile->estimated_node_reserve = static_cast<std::uint64_t>(reserve_nodes);
    profile->actual_node_count = static_cast<std::uint64_t>(m_nodes.size());
    profile->node_capacity_high_water = static_cast<std::uint64_t>(m_node_capacity_high_water);
  }
}

void TreeGravitySolver::evaluateActiveSet(
    const TreeGravitySourceView& source_view,
    const TreeGravityTargetView& target_view,
    const TreeGravityOptions& options,
    TreeGravityProfile* profile,
    const TreeSofteningView& softening_view) const {
  evaluateActiveSet(
      source_view.pos_x_comoving,
      source_view.pos_y_comoving,
      source_view.pos_z_comoving,
      source_view.mass_code,
      target_view.active_particle_index,
      target_view.accel_x_comoving,
      target_view.accel_y_comoving,
      target_view.accel_z_comoving,
      options,
      profile,
      softening_view,
      target_view.previous_acceleration_magnitude_code,
      source_view.source_generation);
}

void TreeGravitySolver::evaluateActiveSet(
    std::span<const double> pos_x_comoving,
    std::span<const double> pos_y_comoving,
    std::span<const double> pos_z_comoving,
    std::span<const double> mass_code,
    std::span<const TreeLocalIndex> active_particle_index,
    std::span<double> accel_x_comoving,
    std::span<double> accel_y_comoving,
    std::span<double> accel_z_comoving,
    const TreeGravityOptions& options,
    TreeGravityProfile* profile,
    const TreeSofteningView& softening_view,
    std::span<const double> previous_acceleration_magnitude_code,
    GravitySourceGeneration source_generation) const {
  validateInputSpans(pos_x_comoving, pos_y_comoving, pos_z_comoving, mass_code);
  validateOptions(options);
  if (!built()) {
    throw std::runtime_error("Tree must be built before traversal");
  }
  if (pos_x_comoving.size() != m_build_source_count) {
    throw std::invalid_argument(
        "Tree traversal source count differs from the source state used to build the tree");
  }
  if (m_build_source_generation.valid()) {
    if (!source_generation.valid() || source_generation != m_build_source_generation) {
      throw std::invalid_argument(
          "Tree traversal source generation differs from the source state used to build the tree");
    }
  } else {
    if (sourceFingerprint(pos_x_comoving, pos_y_comoving, pos_z_comoving, mass_code) !=
        m_build_source_fingerprint) {
      throw std::invalid_argument(
          "Tree traversal source content differs from the legacy source state used to build the tree");
    }
    // Legacy callers without a generation token pay an O(N) identity check.
    // Include resolved source softening so a same-sized sidecar/policy mutation
    // cannot silently traverse a tree whose cached pair softenings are stale.
    for (std::size_t source_index = 0; source_index < pos_x_comoving.size(); ++source_index) {
      if (resolveSourceSofteningEpsilon(source_index, options.softening, softening_view) !=
          m_source_softening_epsilon_comoving[source_index]) {
        throw std::invalid_argument(
            "Tree traversal source softening differs from the legacy source state used to build the tree");
      }
    }
  }
  if (options.multipole_order != m_build_multipole_order ||
      options.max_leaf_size != m_build_max_leaf_size ||
      options.softening.kernel != m_build_softening.kernel ||
      options.softening.epsilon_comoving != m_build_softening.epsilon_comoving) {
    throw std::invalid_argument(
        "Tree traversal options changed a build-time topology/multipole/softening contract; rebuild the tree");
  }
  if (active_particle_index.size() != accel_x_comoving.size() || active_particle_index.size() != accel_y_comoving.size() ||
      active_particle_index.size() != accel_z_comoving.size()) {
    throw std::invalid_argument("Active-set and acceleration spans must match");
  }
  if (!previous_acceleration_magnitude_code.empty() &&
      previous_acceleration_magnitude_code.size() != active_particle_index.size()) {
    throw std::invalid_argument("Previous acceleration magnitude span must match the active-set size");
  }
  if (!softening_view.target_particle_epsilon_comoving.empty() &&
      softening_view.target_particle_epsilon_comoving.size() != active_particle_index.size()) {
    throw std::invalid_argument("Tree gravity target softening sidecar size must match active-set size");
  }
  if (!softening_view.target_particle_epsilon_override_mask.empty() &&
      softening_view.target_particle_epsilon_override_mask.size() != active_particle_index.size()) {
    throw std::invalid_argument("Tree gravity target softening override mask size must match active-set size");
  }
  if (!softening_view.target_species_tag.empty() &&
      softening_view.target_species_tag.size() != active_particle_index.size()) {
    throw std::invalid_argument("Tree gravity target species sidecar size must match active-set size");
  }

  // Validate all failure-capable per-target inputs before entering the OpenMP
  // traversal. The hot tree walk therefore remains exception-free and each
  // active target owns a disjoint acceleration output slot.
  std::vector<double> target_softening_by_active(active_particle_index.size(), 0.0);
  for (std::size_t active_i = 0; active_i < active_particle_index.size(); ++active_i) {
    const TreeLocalIndex particle_index = active_particle_index[active_i];
    if (particle_index >= pos_x_comoving.size()) {
      throw std::out_of_range("Active particle index exceeds particle count");
    }
    target_softening_by_active[active_i] =
        resolveTargetSofteningEpsilon(active_i, particle_index, options.softening, softening_view);
  }

  const auto traversal_start = std::chrono::steady_clock::now();
  std::uint64_t accepted_nodes = 0;
  std::uint64_t opened_nodes = 0;
  std::uint64_t visited_nodes = 0;
  std::uint64_t pp_interactions = 0;

#if COSMOSIM_HAVE_OPENMP
#pragma omp parallel reduction(+ : accepted_nodes, opened_nodes, visited_nodes, pp_interactions)
#endif
  {
    std::vector<TreeLocalIndex> stack;
    stack.reserve(256);

#if COSMOSIM_HAVE_OPENMP
#pragma omp for schedule(dynamic, 8)
#endif
    for (std::ptrdiff_t active_slot = 0;
         active_slot < static_cast<std::ptrdiff_t>(active_particle_index.size());
         ++active_slot) {
      const std::size_t active_i = static_cast<std::size_t>(active_slot);
      const TreeLocalIndex particle_index = active_particle_index[active_i];
      const double px = pos_x_comoving[particle_index];
      const double py = pos_y_comoving[particle_index];
      const double pz = pos_z_comoving[particle_index];
      const double target_softening_comoving = target_softening_by_active[active_i];
      const bool previous_acceleration_supplied = !previous_acceleration_magnitude_code.empty();
      const double previous_acceleration_code = previous_acceleration_supplied
          ? previous_acceleration_magnitude_code[active_i]
          : 0.0;
      const bool previous_acceleration_available =
          previous_acceleration_supplied && std::isfinite(previous_acceleration_code);

      double ax = 0.0;
      double ay = 0.0;
      double az = 0.0;

      stack.clear();
      stack.push_back(0U);

      while (!stack.empty()) {
        const TreeLocalIndex node_index = stack.back();
        stack.pop_back();
        ++visited_nodes;

        const double dx = m_nodes.com_x_comoving[node_index] - px;
        const double dy = m_nodes.com_y_comoving[node_index] - py;
        const double dz = m_nodes.com_z_comoving[node_index] - pz;
        const double r2 = dx * dx + dy * dy + dz * dz;
        const double half_size = m_nodes.half_size_comoving[node_index];
        const double center_dx = m_nodes.center_x_comoving[node_index] - m_nodes.com_x_comoving[node_index];
        const double center_dy = m_nodes.center_y_comoving[node_index] - m_nodes.com_y_comoving[node_index];
        const double center_dz = m_nodes.center_z_comoving[node_index] - m_nodes.com_z_comoving[node_index];
        const double com_offset = std::sqrt(center_dx * center_dx + center_dy * center_dy + center_dz * center_dz);
        const bool is_leaf = m_nodes.child_count[node_index] == 0;
        const bool target_inside_node = !is_leaf &&
            std::abs(px - m_nodes.center_x_comoving[node_index]) <= half_size &&
            std::abs(py - m_nodes.center_y_comoving[node_index]) <= half_size &&
            std::abs(pz - m_nodes.center_z_comoving[node_index]) <= half_size;
        const double r = std::sqrt(r2 + 1.0e-30);
        const bool mac_accept = internal::acceptNodeByMac(
            is_leaf,
            target_inside_node,
            half_size,
            com_offset,
            m_nodes.mass_code[node_index],
            r2,
            previous_acceleration_available,
            previous_acceleration_code,
            options);
        const bool softening_accept = internal::passesSofteningEnvelopeGuard(
            is_leaf,
            half_size,
            r,
            target_softening_comoving,
            m_nodes.softening_min_comoving[node_index],
            m_nodes.softening_max_comoving[node_index]);
        const bool accept = mac_accept && softening_accept;

        if (accept) {
          ++accepted_nodes;
          if (is_leaf) {
            const TreeLocalCount begin = m_nodes.particle_begin[node_index];
            const TreeLocalCount leaf_end = begin + m_nodes.particle_count[node_index];
            for (TreeLocalIndex sorted_i = begin; sorted_i < leaf_end; ++sorted_i) {
              const TreeLocalIndex source_index = m_ordering.sorted_particle_index[sorted_i];
              if (source_index == particle_index) {
                continue;
              }
              const double sx = pos_x_comoving[source_index] - px;
              const double sy = pos_y_comoving[source_index] - py;
              const double sz = pos_z_comoving[source_index] - pz;
              const double sr2 = sx * sx + sy * sy + sz * sz;
              const double pair_epsilon = combineSofteningPairEpsilonUnchecked(
                  m_source_softening_epsilon_comoving[source_index], target_softening_comoving);
              const double factor = options.gravitational_constant_code * mass_code[source_index] *
                  softenedInvR3Unchecked(sr2, pair_epsilon);
              ax += factor * sx;
              ay += factor * sy;
              az += factor * sz;
              ++pp_interactions;
            }
          } else {
            const auto contrib = monopolePlusQuadrupoleAccel(
                m_nodes, node_index, dx, dy, dz, target_softening_comoving, options);
            ax += contrib[0];
            ay += contrib[1];
            az += contrib[2];
          }
        } else {
          ++opened_nodes;
          pushChildrenNearFirst(m_nodes, node_index, px, py, pz, stack);
        }
      }

      accel_x_comoving[active_i] = ax;
      accel_y_comoving[active_i] = ay;
      accel_z_comoving[active_i] = az;
    }
  }

  if (profile != nullptr) {
    profile->visited_nodes += visited_nodes;
    profile->accepted_nodes += accepted_nodes;
    profile->opened_nodes += opened_nodes;
    profile->particle_particle_interactions += pp_interactions;
    profile->traversed_target_count += static_cast<std::uint64_t>(active_particle_index.size());
    profile->traversal_ms +=
        std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - traversal_start).count();
    if (profile->traversed_target_count > 0U) {
      profile->average_interactions_per_target =
          static_cast<double>(profile->particle_particle_interactions) /
          static_cast<double>(profile->traversed_target_count);
    }
  }
}

const TreeNodeSoa& TreeGravitySolver::nodes() const {
  return m_nodes;
}

const TreeMortonOrdering& TreeGravitySolver::ordering() const {
  return m_ordering;
}

TreeBuildGeneration TreeGravitySolver::treeBuildGeneration() const noexcept {
  return m_tree_build_generation;
}

void TreeGravitySolver::appendMemoryReport(core::MemoryReportBuilder& builder) const {
  m_nodes.appendMemoryReport(builder);
  const auto add = [&builder](core::MemorySubsystem subsystem, std::string label, const auto& container) {
    const std::uint64_t capacity_bytes = core::ownedCapacityBytesForContainer(container);
    builder.addEntry(core::MemoryEntry{.subsystem = subsystem,
                                       .lifetime = core::MemoryLifetime::kTransient,
                                       .label = std::move(label),
                                       .current_size_bytes = core::currentSizeBytesForContainer(container),
                                       .owned_capacity_bytes = capacity_bytes,
                                       .high_water_bytes = capacity_bytes,
                                       .estimated_next_step_bytes = 0U,
                                       .uncertainty_note = "retained vector capacity is the observed allocation high-water; next-step requirement is not predicted"});
  };
  add(core::MemorySubsystem::kTree, "tree.ordering.sorted_particle_index", m_ordering.sorted_particle_index);
  add(core::MemorySubsystem::kTree, "tree.ordering.morton_key", m_ordering.morton_key);
  add(core::MemorySubsystem::kTree, "tree.source_softening", m_source_softening_epsilon_comoving);
  add(core::MemorySubsystem::kScratch, "tree.partition_scratch", m_partition_scratch);
}

bool TreeGravitySolver::built() const {
  return m_nodes.size() > 0;
}

TreeLocalIndex TreeGravitySolver::buildNodeRecursive(
    std::span<const double> pos_x_comoving,
    std::span<const double> pos_y_comoving,
    std::span<const double> pos_z_comoving,
    std::span<const double> mass_code,
    TreeLocalIndex begin,
    TreeLocalIndex end,
    double center_x_comoving,
    double center_y_comoving,
    double center_z_comoving,
    double half_size_comoving,
    const TreeGravityOptions& options) {
  (void)mass_code;
  if (m_nodes.size() >= static_cast<std::size_t>(kInvalidTreeLocalIndex)) {
    throw std::overflow_error("Tree node count exceeds the 32-bit node-index contract");
  }
  const TreeLocalIndex node_index = checkedTreeLocalIndex(m_nodes.size(), "tree node count exceeds local index policy");
  m_nodes.center_x_comoving.push_back(center_x_comoving);
  m_nodes.center_y_comoving.push_back(center_y_comoving);
  m_nodes.center_z_comoving.push_back(center_z_comoving);
  m_nodes.half_size_comoving.push_back(half_size_comoving);
  m_nodes.mass_code.push_back(0.0);
  m_nodes.com_x_comoving.push_back(center_x_comoving);
  m_nodes.com_y_comoving.push_back(center_y_comoving);
  m_nodes.com_z_comoving.push_back(center_z_comoving);
  if (options.multipole_order == TreeMultipoleOrder::kQuadrupole) {
    m_nodes.quad_xx.push_back(0.0);
    m_nodes.quad_xy.push_back(0.0);
    m_nodes.quad_xz.push_back(0.0);
    m_nodes.quad_yy.push_back(0.0);
    m_nodes.quad_yz.push_back(0.0);
    m_nodes.quad_zz.push_back(0.0);
    m_nodes.second_moment_trace.push_back(0.0);
  }
  m_nodes.softening_min_comoving.push_back(std::numeric_limits<double>::infinity());
  m_nodes.softening_max_comoving.push_back(0.0);
  m_nodes.child_base.push_back(0U);
  m_nodes.child_count.push_back(0U);
  for (std::uint8_t i = 0; i < 8U; ++i) {
    m_nodes.child_index.push_back(kInvalidTreeLocalIndex);
  }
  m_nodes.particle_begin.push_back(begin);
  m_nodes.particle_count.push_back(end - begin);

  if ((end - begin) <= options.max_leaf_size || half_size_comoving <= 1.0e-12) {
    return node_index;
  }

  std::array<TreeLocalCount, 8> octant_count{};
  for (TreeLocalCount i = begin; i < end; ++i) {
    const TreeLocalIndex particle = m_ordering.sorted_particle_index[i];
    const std::uint8_t octant = octantForParticle(
        pos_x_comoving[particle], pos_y_comoving[particle], pos_z_comoving[particle], center_x_comoving, center_y_comoving,
        center_z_comoving);
    ++octant_count[octant];
  }

  std::array<TreeLocalCount, 9> octant_offsets{};
  octant_offsets[0] = begin;
  for (std::size_t i = 0; i < 8; ++i) {
    octant_offsets[i + 1] = octant_offsets[i] + octant_count[i];
  }

  if (m_partition_scratch.size() < m_ordering.sorted_particle_index.size()) {
    throw std::logic_error("Tree partition workspace is smaller than the source ordering");
  }
  std::array<TreeLocalCount, 8> cursor{};
  for (std::size_t i = 0; i < 8; ++i) {
    cursor[i] = octant_offsets[i];
  }

  for (TreeLocalCount i = begin; i < end; ++i) {
    const TreeLocalIndex particle = m_ordering.sorted_particle_index[i];
    const std::uint8_t octant = octantForParticle(
        pos_x_comoving[particle], pos_y_comoving[particle], pos_z_comoving[particle], center_x_comoving, center_y_comoving,
        center_z_comoving);
    m_partition_scratch[cursor[octant]++] = particle;
  }

  std::copy(
      m_partition_scratch.begin() + begin,
      m_partition_scratch.begin() + end,
      m_ordering.sorted_particle_index.begin() + begin);

  const TreeLocalIndex child_base = checkedTreeLocalIndex(m_nodes.size(), "tree child base exceeds local index policy");
  std::uint8_t non_empty_children = 0;
  const std::size_t child_slot_offset = static_cast<std::size_t>(node_index) * 8U;

  for (std::uint8_t octant = 0; octant < 8U; ++octant) {
    const TreeLocalIndex child_begin = octant_offsets[octant];
    const TreeLocalIndex child_end = octant_offsets[octant + 1];
    if (child_begin == child_end) {
      continue;
    }

    const double child_half = 0.5 * half_size_comoving;
    const double child_center_x = center_x_comoving + ((octant & 4U) != 0U ? child_half : -child_half);
    const double child_center_y = center_y_comoving + ((octant & 2U) != 0U ? child_half : -child_half);
    const double child_center_z = center_z_comoving + ((octant & 1U) != 0U ? child_half : -child_half);

    const TreeLocalIndex built_child_index = buildNodeRecursive(
        pos_x_comoving,
        pos_y_comoving,
        pos_z_comoving,
        mass_code,
        child_begin,
        child_end,
        child_center_x,
        child_center_y,
        child_center_z,
        child_half,
        options);
    m_nodes.child_index[child_slot_offset + octant] = built_child_index;
    ++non_empty_children;
  }

  m_nodes.child_base[node_index] = child_base;
  m_nodes.child_count[node_index] = non_empty_children;
  return node_index;
}

void TreeGravitySolver::accumulateMultipoles(
    std::span<const double> pos_x_comoving,
    std::span<const double> pos_y_comoving,
    std::span<const double> pos_z_comoving,
    std::span<const double> mass_code,
    TreeLocalIndex node_index,
    TreeMultipoleOrder multipole_order) {
  if (m_nodes.child_count[node_index] == 0U) {
    const TreeLocalCount begin = m_nodes.particle_begin[node_index];
    const TreeLocalCount end = begin + m_nodes.particle_count[node_index];
    double total_mass = 0.0;
    double wx = 0.0;
    double wy = 0.0;
    double wz = 0.0;
    for (TreeLocalCount i = begin; i < end; ++i) {
      const TreeLocalIndex particle = m_ordering.sorted_particle_index[i];
      const double m = mass_code[particle];
      m_nodes.softening_min_comoving[node_index] =
          std::min(m_nodes.softening_min_comoving[node_index], m_source_softening_epsilon_comoving[particle]);
      m_nodes.softening_max_comoving[node_index] =
          std::max(m_nodes.softening_max_comoving[node_index], m_source_softening_epsilon_comoving[particle]);
      total_mass += m;
      wx += m * pos_x_comoving[particle];
      wy += m * pos_y_comoving[particle];
      wz += m * pos_z_comoving[particle];
    }
    m_nodes.mass_code[node_index] = total_mass;
    if (total_mass > 0.0) {
      m_nodes.com_x_comoving[node_index] = wx / total_mass;
      m_nodes.com_y_comoving[node_index] = wy / total_mass;
      m_nodes.com_z_comoving[node_index] = wz / total_mass;
    }
    if (multipole_order == TreeMultipoleOrder::kQuadrupole && total_mass > 0.0) {
      const double cx = m_nodes.com_x_comoving[node_index];
      const double cy = m_nodes.com_y_comoving[node_index];
      const double cz = m_nodes.com_z_comoving[node_index];
      double i_xx = 0.0;
      double i_xy = 0.0;
      double i_xz = 0.0;
      double i_yy = 0.0;
      double i_yz = 0.0;
      double i_zz = 0.0;
      for (TreeLocalCount i = begin; i < end; ++i) {
        const TreeLocalIndex particle = m_ordering.sorted_particle_index[i];
        const double m = mass_code[particle];
        const double rx = pos_x_comoving[particle] - cx;
        const double ry = pos_y_comoving[particle] - cy;
        const double rz = pos_z_comoving[particle] - cz;
        i_xx += m * rx * rx;
        i_xy += m * rx * ry;
        i_xz += m * rx * rz;
        i_yy += m * ry * ry;
        i_yz += m * ry * rz;
        i_zz += m * rz * rz;
      }
      const double trace = i_xx + i_yy + i_zz;
      m_nodes.second_moment_trace[node_index] = trace;
      m_nodes.quad_xx[node_index] = 3.0 * i_xx - trace;
      m_nodes.quad_xy[node_index] = 3.0 * i_xy;
      m_nodes.quad_xz[node_index] = 3.0 * i_xz;
      m_nodes.quad_yy[node_index] = 3.0 * i_yy - trace;
      m_nodes.quad_yz[node_index] = 3.0 * i_yz;
      m_nodes.quad_zz[node_index] = 3.0 * i_zz - trace;
    }
    return;
  }

  double total_mass = 0.0;
  double wx = 0.0;
  double wy = 0.0;
  double wz = 0.0;

  const std::size_t child_offset = static_cast<std::size_t>(node_index) * 8U;
  for (std::uint8_t octant = 0; octant < 8U; ++octant) {
    const TreeLocalIndex child = m_nodes.child_index[child_offset + octant];
    if (child == kInvalidTreeLocalIndex) {
      continue;
    }
    accumulateMultipoles(pos_x_comoving, pos_y_comoving, pos_z_comoving, mass_code, child, multipole_order);
    m_nodes.softening_min_comoving[node_index] =
        std::min(m_nodes.softening_min_comoving[node_index], m_nodes.softening_min_comoving[child]);
    m_nodes.softening_max_comoving[node_index] =
        std::max(m_nodes.softening_max_comoving[node_index], m_nodes.softening_max_comoving[child]);
    total_mass += m_nodes.mass_code[child];
    wx += m_nodes.mass_code[child] * m_nodes.com_x_comoving[child];
    wy += m_nodes.mass_code[child] * m_nodes.com_y_comoving[child];
    wz += m_nodes.mass_code[child] * m_nodes.com_z_comoving[child];
  }

  m_nodes.mass_code[node_index] = total_mass;
  if (total_mass > 0.0) {
    m_nodes.com_x_comoving[node_index] = wx / total_mass;
    m_nodes.com_y_comoving[node_index] = wy / total_mass;
    m_nodes.com_z_comoving[node_index] = wz / total_mass;
  }

  if (multipole_order != TreeMultipoleOrder::kQuadrupole || total_mass <= 0.0) {
    return;
  }

  const double cx = m_nodes.com_x_comoving[node_index];
  const double cy = m_nodes.com_y_comoving[node_index];
  const double cz = m_nodes.com_z_comoving[node_index];

  double qxx = 0.0;
  double qxy = 0.0;
  double qxz = 0.0;
  double qyy = 0.0;
  double qyz = 0.0;
  double qzz = 0.0;
  double second_moment_trace = 0.0;

  for (std::uint8_t octant = 0; octant < 8U; ++octant) {
    const TreeLocalIndex child = m_nodes.child_index[child_offset + octant];
    if (child == kInvalidTreeLocalIndex) {
      continue;
    }
    const double dx = m_nodes.com_x_comoving[child] - cx;
    const double dy = m_nodes.com_y_comoving[child] - cy;
    const double dz = m_nodes.com_z_comoving[child] - cz;
    const double m = m_nodes.mass_code[child];
    const double d2 = dx * dx + dy * dy + dz * dz;

    qxx += m_nodes.quad_xx[child] + m * (3.0 * dx * dx - d2);
    qxy += m_nodes.quad_xy[child] + m * (3.0 * dx * dy);
    qxz += m_nodes.quad_xz[child] + m * (3.0 * dx * dz);
    qyy += m_nodes.quad_yy[child] + m * (3.0 * dy * dy - d2);
    qyz += m_nodes.quad_yz[child] + m * (3.0 * dy * dz);
    qzz += m_nodes.quad_zz[child] + m * (3.0 * dz * dz - d2);
    second_moment_trace += m_nodes.second_moment_trace[child] + m * d2;
  }

  m_nodes.quad_xx[node_index] = qxx;
  m_nodes.quad_xy[node_index] = qxy;
  m_nodes.quad_xz[node_index] = qxz;
  m_nodes.quad_yy[node_index] = qyy;
  m_nodes.quad_yz[node_index] = qyz;
  m_nodes.quad_zz[node_index] = qzz;
  m_nodes.second_moment_trace[node_index] = second_moment_trace;
}

}  // namespace cosmosim::gravity
