#include "cosmosim/gravity/tree_ordering.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace cosmosim::gravity {
namespace {

[[nodiscard]] std::uint64_t expandBits21(std::uint32_t value) {
  std::uint64_t x = value & 0x1FFFFFU;
  x = (x | (x << 32U)) & 0x1F00000000FFFFULL;
  x = (x | (x << 16U)) & 0x1F0000FF0000FFULL;
  x = (x | (x << 8U)) & 0x100F00F00F00F00FULL;
  x = (x | (x << 4U)) & 0x10C30C30C30C30C3ULL;
  x = (x | (x << 2U)) & 0x1249249249249249ULL;
  return x;
}

[[nodiscard]] std::uint64_t morton3D(std::uint32_t x, std::uint32_t y, std::uint32_t z) {
  return (expandBits21(x) << 2U) | (expandBits21(y) << 1U) | expandBits21(z);
}

}  // namespace

double TreeBounds::maxExtentComoving() const {
  return std::max({max_x_comoving - min_x_comoving, max_y_comoving - min_y_comoving, max_z_comoving - min_z_comoving});
}

TreeBounds computeTreeBounds(
    std::span<const double> pos_x_comoving,
    std::span<const double> pos_y_comoving,
    std::span<const double> pos_z_comoving) {
  if (pos_x_comoving.size() != pos_y_comoving.size() || pos_x_comoving.size() != pos_z_comoving.size()) {
    throw std::invalid_argument("Tree bounds requires equal position spans");
  }
  if (pos_x_comoving.empty()) {
    return {};
  }

  TreeBounds bounds;
  bounds.min_x_comoving = std::numeric_limits<double>::infinity();
  bounds.min_y_comoving = std::numeric_limits<double>::infinity();
  bounds.min_z_comoving = std::numeric_limits<double>::infinity();
  bounds.max_x_comoving = -std::numeric_limits<double>::infinity();
  bounds.max_y_comoving = -std::numeric_limits<double>::infinity();
  bounds.max_z_comoving = -std::numeric_limits<double>::infinity();

  for (std::size_t i = 0; i < pos_x_comoving.size(); ++i) {
    bounds.min_x_comoving = std::min(bounds.min_x_comoving, pos_x_comoving[i]);
    bounds.min_y_comoving = std::min(bounds.min_y_comoving, pos_y_comoving[i]);
    bounds.min_z_comoving = std::min(bounds.min_z_comoving, pos_z_comoving[i]);
    bounds.max_x_comoving = std::max(bounds.max_x_comoving, pos_x_comoving[i]);
    bounds.max_y_comoving = std::max(bounds.max_y_comoving, pos_y_comoving[i]);
    bounds.max_z_comoving = std::max(bounds.max_z_comoving, pos_z_comoving[i]);
  }
  return bounds;
}

TreeMortonOrdering buildMortonOrdering(
    std::span<const double> pos_x_comoving,
    std::span<const double> pos_y_comoving,
    std::span<const double> pos_z_comoving) {
  const TreeBounds bounds = computeTreeBounds(pos_x_comoving, pos_y_comoving, pos_z_comoving);

  TreeMortonOrdering ordering;
  ordering.sorted_particle_index.resize(pos_x_comoving.size());
  ordering.morton_key.resize(pos_x_comoving.size());

  if (pos_x_comoving.empty()) {
    return ordering;
  }

  const double extent = std::max(bounds.maxExtentComoving(), 1.0e-12);
  constexpr double k_grid = static_cast<double>((1U << 21U) - 1U);

  for (std::size_t i = 0; i < pos_x_comoving.size(); ++i) {
    const double nx = std::clamp((pos_x_comoving[i] - bounds.min_x_comoving) / extent, 0.0, 1.0);
    const double ny = std::clamp((pos_y_comoving[i] - bounds.min_y_comoving) / extent, 0.0, 1.0);
    const double nz = std::clamp((pos_z_comoving[i] - bounds.min_z_comoving) / extent, 0.0, 1.0);
    const std::uint32_t qx = static_cast<std::uint32_t>(std::llround(nx * k_grid));
    const std::uint32_t qy = static_cast<std::uint32_t>(std::llround(ny * k_grid));
    const std::uint32_t qz = static_cast<std::uint32_t>(std::llround(nz * k_grid));
    ordering.sorted_particle_index[i] = static_cast<std::uint32_t>(i);
    ordering.morton_key[i] = morton3D(qx, qy, qz);
  }

  std::vector<std::size_t> order(ordering.sorted_particle_index.size());
  for (std::size_t i = 0; i < order.size(); ++i) {
    order[i] = i;
  }
  std::stable_sort(order.begin(), order.end(), [&ordering](std::size_t lhs, std::size_t rhs) {
    return ordering.morton_key[lhs] < ordering.morton_key[rhs];
  });

  std::vector<std::uint32_t> sorted_particle(order.size());
  std::vector<std::uint64_t> sorted_key(order.size());
  for (std::size_t i = 0; i < order.size(); ++i) {
    sorted_particle[i] = ordering.sorted_particle_index[order[i]];
    sorted_key[i] = ordering.morton_key[order[i]];
  }
  ordering.sorted_particle_index = std::move(sorted_particle);
  ordering.morton_key = std::move(sorted_key);
  return ordering;
}

}  // namespace cosmosim::gravity
