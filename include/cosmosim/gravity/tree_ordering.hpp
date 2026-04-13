#pragma once

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace cosmosim::gravity {

// Morton ordering helper for locality-friendly particle reindexing.
struct TreeMortonOrdering {
  std::vector<std::uint32_t> sorted_particle_index;
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

[[nodiscard]] TreeMortonOrdering buildMortonOrdering(
    std::span<const double> pos_x_comoving,
    std::span<const double> pos_y_comoving,
    std::span<const double> pos_z_comoving);

}  // namespace cosmosim::gravity
