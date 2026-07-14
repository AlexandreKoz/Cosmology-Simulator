#pragma once

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include "cosmosim/core/memory_accounting.hpp"
#include "cosmosim/gravity/tree_ordering.hpp"
#include "cosmosim/gravity/tree_softening.hpp"

namespace cosmosim::gravity {

enum class TreeOpeningCriterion {
  kBarnesHutGeometric,
  kBarnesHutComDistance,
  kRelativeForceError,
};

enum class TreeMultipoleOrder {
  kMonopole = 0,
  kQuadrupole = 2,
};

struct TreeGravityOptions {
  // Keep the direct solver API aligned with the typed production-config
  // default certified by the periodic Ewald validation suite.
  double opening_theta = 0.7;
  TreeOpeningCriterion opening_criterion = TreeOpeningCriterion::kBarnesHutComDistance;
  TreeMultipoleOrder multipole_order = TreeMultipoleOrder::kQuadrupole;
  double gravitational_constant_code = 1.0;
  std::size_t max_leaf_size = 16;
  TreeSofteningPolicy softening{};
  // Relative force-error MAC parameters. For a node of mass M and width l,
  // accept when G M l^2 <= alpha max(|a_prev|, a_floor) r^4.
  double relative_force_tolerance = 0.005;
  double relative_force_acceleration_floor_code = 1.0e-30;
};

struct TreeGravityProfile {
  double build_ms = 0.0;
  double multipole_ms = 0.0;
  double traversal_ms = 0.0;
  std::uint64_t build_count = 0;
  std::uint64_t multipole_refresh_count = 0;
  std::uint64_t visited_nodes = 0;
  std::uint64_t accepted_nodes = 0;
  std::uint64_t opened_nodes = 0;
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
  // Trace of the raw second central moment I_ij = sum(m * dx_i * dx_j).
  // The trace is redundant for an unsoftened 1/r harmonic kernel, but is
  // required to evaluate a true second-order expansion of softened and
  // TreePM-screened kernels.
  std::vector<double> second_moment_trace;
  std::vector<double> softening_min_comoving;
  std::vector<double> softening_max_comoving;
  std::vector<std::uint32_t> child_base;
  std::vector<std::uint8_t> child_count;
  std::vector<std::uint32_t> child_index;
  std::vector<std::uint32_t> particle_begin;
  std::vector<std::uint32_t> particle_count;

  [[nodiscard]] std::size_t size() const;
  void clear();
  void reserve(std::size_t count);
  void appendMemoryReport(core::MemoryReportBuilder& builder) const;
};

class TreeGravitySolver {
 public:
  struct TreeGravitySourceView {
    std::span<const double> pos_x_comoving;
    std::span<const double> pos_y_comoving;
    std::span<const double> pos_z_comoving;
    std::span<const double> mass_code;
  };

  struct TreeGravityTargetView {
    std::span<const std::uint32_t> active_particle_index;
    std::span<double> accel_x_comoving;
    std::span<double> accel_y_comoving;
    std::span<double> accel_z_comoving;
    // Optional magnitude of the previous total acceleration for each active
    // target. An empty span, or a non-finite entry, selects the deterministic
    // COM-distance fallback for that target.
    std::span<const double> previous_acceleration_magnitude_code{};
  };

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
      const TreeGravitySourceView& source_view,
      const TreeGravityTargetView& target_view,
      const TreeGravityOptions& options,
      TreeGravityProfile* profile = nullptr,
      const TreeSofteningView& softening_view = {}) const;

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
      const TreeSofteningView& softening_view = {},
      std::span<const double> previous_acceleration_magnitude_code = {}) const;

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
