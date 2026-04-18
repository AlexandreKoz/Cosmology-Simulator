#pragma once

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include "cosmosim/gravity/pm_solver.hpp"
#include "cosmosim/gravity/tree_gravity.hpp"
#include "cosmosim/gravity/tree_pm_split_kernel.hpp"

namespace cosmosim::gravity {

// Shared compact active-set view for force accumulation ownership.
struct TreePmForceAccumulatorView {
  std::span<const std::uint32_t> active_particle_index;
  std::span<double> accel_x_comoving;
  std::span<double> accel_y_comoving;
  std::span<double> accel_z_comoving;

  void reset() const;
  void addToActiveSlot(std::size_t active_slot, double ax_comoving, double ay_comoving, double az_comoving) const;
};

struct TreePmOptions {
  PmSolveOptions pm_options{};
  TreeGravityOptions tree_options{};
  TreePmSplitPolicy split_policy{};
};

struct TreePmDiagnostics {
  double mesh_spacing_comoving = 0.0;
  double asmth_cells = 0.0;
  double rcut_cells = 0.0;
  double split_scale_comoving = 0.0;
  double cutoff_radius_comoving = 0.0;
  double short_range_factor_at_split = 0.0;
  double long_range_factor_at_split = 0.0;
  double short_range_factor_at_cutoff = 0.0;
  double long_range_factor_at_cutoff = 0.0;
  double composition_error_at_split = 0.0;
  double max_relative_composition_error = 0.0;
  std::uint64_t residual_pruned_nodes = 0;
  std::uint64_t residual_pair_skips_cutoff = 0;
  std::uint64_t residual_pair_evaluations = 0;
};

struct TreePmProfileEvent {
  PmProfileEvent pm_profile{};
  TreeGravityProfile tree_profile{};
  double tree_short_range_ms = 0.0;
  double coupling_overhead_ms = 0.0;
};

// Thin coordinator that makes TreePM ownership explicit and auditable.
class TreePmCoordinator {
 public:
  explicit TreePmCoordinator(PmGridShape pm_shape);

  void solveActiveSetWithPmCadence(
      std::span<const double> pos_x_comoving,
      std::span<const double> pos_y_comoving,
      std::span<const double> pos_z_comoving,
      std::span<const double> mass_code,
      const TreePmForceAccumulatorView& accumulator,
      const TreePmOptions& options,
      bool refresh_long_range_field,
      TreePmProfileEvent* profile = nullptr,
      TreePmDiagnostics* diagnostics = nullptr);

  void solveActiveSet(
      std::span<const double> pos_x_comoving,
      std::span<const double> pos_y_comoving,
      std::span<const double> pos_z_comoving,
      std::span<const double> mass_code,
      const TreePmForceAccumulatorView& accumulator,
      const TreePmOptions& options,
      TreePmProfileEvent* profile = nullptr,
      TreePmDiagnostics* diagnostics = nullptr);

 private:
  void evaluateShortRangeResidual(
      std::span<const double> pos_x_comoving,
      std::span<const double> pos_y_comoving,
      std::span<const double> pos_z_comoving,
      std::span<const double> mass_code,
      const TreePmForceAccumulatorView& accumulator,
      const TreePmOptions& options,
      TreeGravityProfile* tree_profile);

  struct ResidualTraversalStats {
    std::uint64_t pruned_nodes = 0;
    std::uint64_t pair_skips_cutoff = 0;
    std::uint64_t pair_evaluations = 0;
  };

  PmGridShape m_shape;
  PmGridStorage m_grid;
  PmSolver m_pm_solver;
  TreeGravitySolver m_tree_solver;

  // Reused compact-sidecar scratch arrays to avoid per-step heap churn.
  std::vector<double> m_active_pos_x_comoving;
  std::vector<double> m_active_pos_y_comoving;
  std::vector<double> m_active_pos_z_comoving;
  std::vector<double> m_active_pm_ax_comoving;
  std::vector<double> m_active_pm_ay_comoving;
  std::vector<double> m_active_pm_az_comoving;
  ResidualTraversalStats m_last_residual_stats;
  bool m_has_cached_long_range_field = false;
};

[[nodiscard]] TreePmDiagnostics computeTreePmDiagnostics(const TreePmSplitPolicy& split_policy);

}  // namespace cosmosim::gravity
