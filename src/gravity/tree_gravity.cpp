#include "cosmosim/gravity/tree_gravity.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace cosmosim::gravity {
namespace {

struct OctantSpan {
  std::uint32_t begin = 0;
  std::uint32_t end = 0;
};

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
  child_base.clear();
  child_count.clear();
  child_index.clear();
  particle_begin.clear();
  particle_count.clear();
}

void TreeNodeSoa::reserve(std::size_t count) {
  center_x_comoving.reserve(count);
  center_y_comoving.reserve(count);
  center_z_comoving.reserve(count);
  half_size_comoving.reserve(count);
  mass_code.reserve(count);
  com_x_comoving.reserve(count);
  com_y_comoving.reserve(count);
  com_z_comoving.reserve(count);
  child_base.reserve(count);
  child_count.reserve(count);
  child_index.reserve(count * 8U);
  particle_begin.reserve(count);
  particle_count.reserve(count);
}

void TreeGravitySolver::build(
    std::span<const double> pos_x_comoving,
    std::span<const double> pos_y_comoving,
    std::span<const double> pos_z_comoving,
    std::span<const double> mass_code,
    const TreeGravityOptions& options,
    TreeGravityProfile* profile) {
  const auto build_start = std::chrono::steady_clock::now();
  validateInputSpans(pos_x_comoving, pos_y_comoving, pos_z_comoving, mass_code);
  if (options.opening_theta <= 0.0 || options.gravitational_constant_code <= 0.0 || options.max_leaf_size == 0) {
    throw std::invalid_argument("Invalid tree gravity options");
  }

  m_nodes.clear();
  m_ordering = buildMortonOrdering(pos_x_comoving, pos_y_comoving, pos_z_comoving);
  if (pos_x_comoving.empty()) {
    if (profile != nullptr) {
      *profile = {};
      profile->build_ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - build_start).count();
    }
    return;
  }

  const TreeBounds bounds = computeTreeBounds(pos_x_comoving, pos_y_comoving, pos_z_comoving);
  const double max_extent = std::max(bounds.maxExtentComoving(), 1.0e-10);
  const double center_x_comoving = 0.5 * (bounds.min_x_comoving + bounds.max_x_comoving);
  const double center_y_comoving = 0.5 * (bounds.min_y_comoving + bounds.max_y_comoving);
  const double center_z_comoving = 0.5 * (bounds.min_z_comoving + bounds.max_z_comoving);

  m_nodes.reserve(pos_x_comoving.size() * 2U);
  const std::uint32_t root_index = buildNodeRecursive(
      pos_x_comoving,
      pos_y_comoving,
      pos_z_comoving,
      mass_code,
      0,
      static_cast<std::uint32_t>(m_ordering.sorted_particle_index.size()),
      center_x_comoving,
      center_y_comoving,
      center_z_comoving,
      0.5 * max_extent * (1.0 + 1.0e-8),
      options);
  (void)root_index;

  const auto multipole_start = std::chrono::steady_clock::now();
  accumulateMultipoles(pos_x_comoving, pos_y_comoving, pos_z_comoving, mass_code, 0);
  const auto multipole_stop = std::chrono::steady_clock::now();

  if (profile != nullptr) {
    profile->build_ms = std::chrono::duration<double, std::milli>(multipole_start - build_start).count();
    profile->multipole_ms = std::chrono::duration<double, std::milli>(multipole_stop - multipole_start).count();
  }
}

void TreeGravitySolver::evaluateActiveSet(
    std::span<const double> pos_x_comoving,
    std::span<const double> pos_y_comoving,
    std::span<const double> pos_z_comoving,
    std::span<const double> mass_code,
    std::span<const std::uint32_t> active_particle_index,
    std::span<double> accel_x_comoving,
    std::span<double> accel_y_comoving,
    std::span<double> accel_z_comoving,
    const TreeGravityOptions& options,
    TreeGravityProfile* profile) const {
  validateInputSpans(pos_x_comoving, pos_y_comoving, pos_z_comoving, mass_code);
  if (!built()) {
    throw std::runtime_error("Tree must be built before traversal");
  }
  if (active_particle_index.size() != accel_x_comoving.size() || active_particle_index.size() != accel_y_comoving.size() ||
      active_particle_index.size() != accel_z_comoving.size()) {
    throw std::invalid_argument("Active-set and acceleration spans must match");
  }

  const auto traversal_start = std::chrono::steady_clock::now();
  std::uint64_t accepted_nodes = 0;
  std::uint64_t visited_nodes = 0;
  std::uint64_t pp_interactions = 0;

  for (std::size_t active_i = 0; active_i < active_particle_index.size(); ++active_i) {
    const std::uint32_t particle_index = active_particle_index[active_i];
    if (particle_index >= pos_x_comoving.size()) {
      throw std::out_of_range("Active particle index exceeds particle count");
    }

    const double px = pos_x_comoving[particle_index];
    const double py = pos_y_comoving[particle_index];
    const double pz = pos_z_comoving[particle_index];

    double ax = 0.0;
    double ay = 0.0;
    double az = 0.0;

    std::vector<std::uint32_t> stack;
    stack.push_back(0U);

    while (!stack.empty()) {
      const std::uint32_t node_index = stack.back();
      stack.pop_back();
      ++visited_nodes;

      const double dx = m_nodes.com_x_comoving[node_index] - px;
      const double dy = m_nodes.com_y_comoving[node_index] - py;
      const double dz = m_nodes.com_z_comoving[node_index] - pz;
      const double r2 = dx * dx + dy * dy + dz * dz;
      const double half_size = m_nodes.half_size_comoving[node_index];
      const double l_over_r = (2.0 * half_size) / std::sqrt(r2 + 1.0e-30);
      const bool is_leaf = m_nodes.child_count[node_index] == 0;
      const bool accept = is_leaf || (l_over_r < options.opening_theta);

      if (accept) {
        ++accepted_nodes;
        if (is_leaf) {
          const std::uint32_t begin = m_nodes.particle_begin[node_index];
          const std::uint32_t end = begin + m_nodes.particle_count[node_index];
          for (std::uint32_t sorted_i = begin; sorted_i < end; ++sorted_i) {
            const std::uint32_t source_index = m_ordering.sorted_particle_index[sorted_i];
            if (source_index == particle_index) {
              continue;
            }
            const double sx = pos_x_comoving[source_index] - px;
            const double sy = pos_y_comoving[source_index] - py;
            const double sz = pos_z_comoving[source_index] - pz;
            const double sr2 = sx * sx + sy * sy + sz * sz;
            const double factor = options.gravitational_constant_code * mass_code[source_index] *
                softenedInvR3(sr2, options.softening);
            ax += factor * sx;
            ay += factor * sy;
            az += factor * sz;
            ++pp_interactions;
          }
        } else {
          const double factor = options.gravitational_constant_code * m_nodes.mass_code[node_index] *
              softenedInvR3(r2, options.softening);
          ax += factor * dx;
          ay += factor * dy;
          az += factor * dz;
        }
      } else {
        const std::size_t child_offset = static_cast<std::size_t>(node_index) * 8U;
        for (std::uint8_t octant = 0; octant < 8U; ++octant) {
          const std::uint32_t child = m_nodes.child_index[child_offset + octant];
          if (child != std::numeric_limits<std::uint32_t>::max()) {
            stack.push_back(child);
          }
        }
      }
    }

    accel_x_comoving[active_i] = ax;
    accel_y_comoving[active_i] = ay;
    accel_z_comoving[active_i] = az;
  }

  if (profile != nullptr) {
    profile->visited_nodes += visited_nodes;
    profile->accepted_nodes += accepted_nodes;
    profile->particle_particle_interactions += pp_interactions;
    profile->traversal_ms +=
        std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - traversal_start).count();
    if (!active_particle_index.empty()) {
      profile->average_interactions_per_target =
          static_cast<double>(profile->particle_particle_interactions) / static_cast<double>(active_particle_index.size());
    }
  }
}

const TreeNodeSoa& TreeGravitySolver::nodes() const {
  return m_nodes;
}

const TreeMortonOrdering& TreeGravitySolver::ordering() const {
  return m_ordering;
}

bool TreeGravitySolver::built() const {
  return m_nodes.size() > 0;
}

std::uint32_t TreeGravitySolver::buildNodeRecursive(
    std::span<const double> pos_x_comoving,
    std::span<const double> pos_y_comoving,
    std::span<const double> pos_z_comoving,
    std::span<const double> mass_code,
    std::uint32_t begin,
    std::uint32_t end,
    double center_x_comoving,
    double center_y_comoving,
    double center_z_comoving,
    double half_size_comoving,
    const TreeGravityOptions& options) {
  (void)mass_code;
  const std::uint32_t node_index = static_cast<std::uint32_t>(m_nodes.size());
  m_nodes.center_x_comoving.push_back(center_x_comoving);
  m_nodes.center_y_comoving.push_back(center_y_comoving);
  m_nodes.center_z_comoving.push_back(center_z_comoving);
  m_nodes.half_size_comoving.push_back(half_size_comoving);
  m_nodes.mass_code.push_back(0.0);
  m_nodes.com_x_comoving.push_back(center_x_comoving);
  m_nodes.com_y_comoving.push_back(center_y_comoving);
  m_nodes.com_z_comoving.push_back(center_z_comoving);
  m_nodes.child_base.push_back(0U);
  m_nodes.child_count.push_back(0U);
  for (std::uint8_t i = 0; i < 8U; ++i) {
    m_nodes.child_index.push_back(std::numeric_limits<std::uint32_t>::max());
  }
  m_nodes.particle_begin.push_back(begin);
  m_nodes.particle_count.push_back(end - begin);

  if ((end - begin) <= options.max_leaf_size || half_size_comoving <= 1.0e-12) {
    return node_index;
  }

  std::array<std::uint32_t, 8> octant_count{};
  for (std::uint32_t i = begin; i < end; ++i) {
    const std::uint32_t particle = m_ordering.sorted_particle_index[i];
    const std::uint8_t octant = octantForParticle(
        pos_x_comoving[particle], pos_y_comoving[particle], pos_z_comoving[particle], center_x_comoving, center_y_comoving,
        center_z_comoving);
    ++octant_count[octant];
  }

  std::array<std::uint32_t, 9> octant_offsets{};
  octant_offsets[0] = begin;
  for (std::size_t i = 0; i < 8; ++i) {
    octant_offsets[i + 1] = octant_offsets[i] + octant_count[i];
  }

  std::vector<std::uint32_t> scratch(end - begin, 0U);
  std::array<std::uint32_t, 8> cursor{};
  for (std::size_t i = 0; i < 8; ++i) {
    cursor[i] = octant_offsets[i] - begin;
  }

  for (std::uint32_t i = begin; i < end; ++i) {
    const std::uint32_t particle = m_ordering.sorted_particle_index[i];
    const std::uint8_t octant = octantForParticle(
        pos_x_comoving[particle], pos_y_comoving[particle], pos_z_comoving[particle], center_x_comoving, center_y_comoving,
        center_z_comoving);
    scratch[cursor[octant]++] = particle;
  }

  std::copy(scratch.begin(), scratch.end(), m_ordering.sorted_particle_index.begin() + begin);

  const std::uint32_t child_base = static_cast<std::uint32_t>(m_nodes.size());
  std::uint8_t non_empty_children = 0;
  const std::size_t child_slot_offset = static_cast<std::size_t>(node_index) * 8U;

  for (std::uint8_t octant = 0; octant < 8U; ++octant) {
    const std::uint32_t child_begin = octant_offsets[octant];
    const std::uint32_t child_end = octant_offsets[octant + 1];
    if (child_begin == child_end) {
      continue;
    }

    const double child_half = 0.5 * half_size_comoving;
    const double child_center_x = center_x_comoving + ((octant & 4U) != 0U ? child_half : -child_half);
    const double child_center_y = center_y_comoving + ((octant & 2U) != 0U ? child_half : -child_half);
    const double child_center_z = center_z_comoving + ((octant & 1U) != 0U ? child_half : -child_half);

    const std::uint32_t built_child_index = buildNodeRecursive(
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
    std::uint32_t node_index) {
  if (m_nodes.child_count[node_index] == 0U) {
    const std::uint32_t begin = m_nodes.particle_begin[node_index];
    const std::uint32_t end = begin + m_nodes.particle_count[node_index];
    double total_mass = 0.0;
    double wx = 0.0;
    double wy = 0.0;
    double wz = 0.0;
    for (std::uint32_t i = begin; i < end; ++i) {
      const std::uint32_t particle = m_ordering.sorted_particle_index[i];
      const double m = mass_code[particle];
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
    return;
  }

  double total_mass = 0.0;
  double wx = 0.0;
  double wy = 0.0;
  double wz = 0.0;

  const std::size_t child_offset = static_cast<std::size_t>(node_index) * 8U;
  for (std::uint8_t octant = 0; octant < 8U; ++octant) {
    const std::uint32_t child = m_nodes.child_index[child_offset + octant];
    if (child == std::numeric_limits<std::uint32_t>::max()) {
      continue;
    }
    accumulateMultipoles(pos_x_comoving, pos_y_comoving, pos_z_comoving, mass_code, child);
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
}

}  // namespace cosmosim::gravity
