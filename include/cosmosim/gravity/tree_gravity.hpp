#pragma once

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include "cosmosim/gravity/tree_ordering.hpp"
#include "cosmosim/gravity/tree_softening.hpp"

namespace cosmosim::gravity {

enum class TreeOpeningCriterion {
  kBarnesHutGeometric,
  kBarnesHutComDistance,
};

enum class TreeMultipoleOrder {
  kMonopole = 0,
  kQuadrupole = 2,
};

struct TreeGravityOptions {
  double opening_theta = 0.6;
  TreeOpeningCriterion opening_criterion = TreeOpeningCriterion::kBarnesHutComDistance;
  TreeMultipoleOrder multipole_order = TreeMultipoleOrder::kQuadrupole;
  double gravitational_constant_code = 1.0;
  std::size_t max_leaf_size = 16;
  TreeSofteningPolicy softening{};
};

struct TreeGravityProfile {
  double build_ms = 0.0;
  double multipole_ms = 0.0;
  double traversal_ms = 0.0;
  std::uint64_t visited_nodes = 0;
  std::uint64_t accepted_nodes = 0;
  std::uint64_t particle_particle_interactions = 0;
  double average_interactions_per_target = 0.0;
};

// Compact SoA tree node representation with fixed 8-way child fanout sidecars.
struct TreeNodeSoa {
  std::vector<double> center_x_comoving;
  std::vector<double> center_y_comoving;
  std::vector<double> center_z_comoving;
  std::vector<double> half_size_comoving;
  std::vector<double> mass_code;
  std::vector<double> com_x_comoving;
  std::vector<double> com_y_comoving;
  std::vector<double> com_z_comoving;
  // Symmetric, traceless quadrupole tensor components Q_ij around each node COM.
  std::vector<double> quad_xx;
  std::vector<double> quad_xy;
  std::vector<double> quad_xz;
  std::vector<double> quad_yy;
  std::vector<double> quad_yz;
  std::vector<double> quad_zz;
  std::vector<double> softening_max_comoving;
  std::vector<std::uint32_t> child_base;
  std::vector<std::uint8_t> child_count;
  std::vector<std::uint32_t> child_index;
  std::vector<std::uint32_t> particle_begin;
  std::vector<std::uint32_t> particle_count;

  [[nodiscard]] std::size_t size() const;
  void clear();
  void reserve(std::size_t count);
};

class TreeGravitySolver {
 public:
  TreeGravitySolver() = default;

  void build(
      std::span<const double> pos_x_comoving,
      std::span<const double> pos_y_comoving,
      std::span<const double> pos_z_comoving,
      std::span<const double> mass_code,
      const TreeGravityOptions& options,
      TreeGravityProfile* profile = nullptr,
      const TreeSofteningView& softening_view = {});

  void evaluateActiveSet(
      std::span<const double> pos_x_comoving,
      std::span<const double> pos_y_comoving,
      std::span<const double> pos_z_comoving,
      std::span<const double> mass_code,
      std::span<const std::uint32_t> active_particle_index,
      std::span<double> accel_x_comoving,
      std::span<double> accel_y_comoving,
      std::span<double> accel_z_comoving,
      const TreeGravityOptions& options,
      TreeGravityProfile* profile = nullptr,
      const TreeSofteningView& softening_view = {}) const;

  [[nodiscard]] const TreeNodeSoa& nodes() const;
  [[nodiscard]] const TreeMortonOrdering& ordering() const;

 private:
  [[nodiscard]] bool built() const;
  [[nodiscard]] std::uint32_t buildNodeRecursive(
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
      const TreeGravityOptions& options);
  void accumulateMultipoles(
      std::span<const double> pos_x_comoving,
      std::span<const double> pos_y_comoving,
      std::span<const double> pos_z_comoving,
      std::span<const double> mass_code,
      std::uint32_t node_index,
      TreeMultipoleOrder multipole_order);

  TreeNodeSoa m_nodes;
  TreeMortonOrdering m_ordering;
  std::vector<double> m_source_softening_epsilon_comoving;
};

}  // namespace cosmosim::gravity
