#include "cosmosim/gravity/tree_ordering.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace cosmosim::gravity {
namespace {

void radixSortMortonOrdering(TreeMortonOrdering& ordering) {
  const std::size_t count = ordering.morton_key.size();
  if (count < 2U) {
    return;
  }

  std::vector<std::uint64_t> scratch_key(count);
  std::vector<std::uint32_t> scratch_index(count);
  bool source_is_primary = true;
  for (unsigned shift = 0U; shift < 64U; shift += 8U) {
    std::array<std::size_t, 256> bucket_count{};
    const auto& source_key = source_is_primary ? ordering.morton_key : scratch_key;
    const auto& source_index = source_is_primary ? ordering.sorted_particle_index : scratch_index;
    auto& destination_key = source_is_primary ? scratch_key : ordering.morton_key;
    auto& destination_index = source_is_primary ? scratch_index : ordering.sorted_particle_index;

    for (const std::uint64_t key : source_key) {
      ++bucket_count[static_cast<std::size_t>((key >> shift) & 0xFFU)];
    }
    std::array<std::size_t, 256> bucket_offset{};
    for (std::size_t bucket = 1U; bucket < bucket_offset.size(); ++bucket) {
      bucket_offset[bucket] = bucket_offset[bucket - 1U] + bucket_count[bucket - 1U];
    }
    for (std::size_t i = 0U; i < count; ++i) {
      const std::size_t bucket = static_cast<std::size_t>((source_key[i] >> shift) & 0xFFU);
      const std::size_t destination = bucket_offset[bucket]++;
      destination_key[destination] = source_key[i];
      destination_index[destination] = source_index[i];
    }
    source_is_primary = !source_is_primary;
  }
  // Eight byte passes return the result to the primary arrays. Keep this
  // assertion-like guard explicit in case the radix width is changed later.
  if (!source_is_primary) {
    ordering.morton_key = std::move(scratch_key);
    ordering.sorted_particle_index = std::move(scratch_index);
  }
}

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
    if (!std::isfinite(pos_x_comoving[i]) || !std::isfinite(pos_y_comoving[i]) ||
        !std::isfinite(pos_z_comoving[i])) {
      throw std::invalid_argument(
          "Tree bounds requires finite particle coordinates");
    }
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
  if (pos_x_comoving.size() >
      static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max())) {
    throw std::overflow_error(
        "Morton ordering exceeds the 32-bit particle-index contract");
  }
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
    ordering.sorted_particle_index[i] = checkedTreeLocalIndex(i, "Tree Morton particle index exceeds local-index policy");
    ordering.morton_key[i] = morton3D(qx, qy, qz);
  }

  // Morton keys are fixed-width integers, so a stable byte-radix pass avoids
  // the former O(N log N) comparison sort and its extra size_t permutation.
  radixSortMortonOrdering(ordering);
  return ordering;
}

}  // namespace cosmosim::gravity
