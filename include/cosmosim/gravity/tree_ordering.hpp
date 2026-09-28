#pragma once

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include "cosmosim/gravity/tree_index.hpp"

namespace cosmosim::gravity {

// Morton ordering helper for locality-friendly particle reindexing.
struct TreeMortonOrdering {
  std::vector<TreeLocalIndex> sorted_particle_index;
  std::vector<std::uint64_t> morton_key;
};

struct TreeBounds {
  double min_x_comoving = 0.0;
  double min_y_comoving = 0.0;
  double min_z_comoving = 0.0;
  double max_x_comoving = 0.0;
  double max_y_comoving = 0.0;
  double max_z_comoving = 0.0;

  [[nodiscard]] double maxExtentComoving() const;
};

[[nodiscard]] TreeBounds computeTreeBounds(
    std::span<const double> pos_x_comoving,
    std::span<const double> pos_y_comoving,
    std::span<const double> pos_z_comoving);

// Build the deterministic 21-bit Morton ordering into retained caller-owned
// storage. The scratch spans must already cover the source population; this
// routine never allocates a second population-scale ordering internally.
void buildMortonOrderingInPlace(
    std::span<const double> pos_x_comoving,
    std::span<const double> pos_y_comoving,
    std::span<const double> pos_z_comoving,
    TreeMortonOrdering& ordering,
    std::span<std::uint64_t> scratch_key,
    std::span<TreeLocalIndex> scratch_index);

}  // namespace cosmosim::gravity
