#pragma once

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include "cosmosim/core/memory_accounting.hpp"
#include "cosmosim/gravity/pm_solver.hpp"
#include "cosmosim/gravity/tree_gravity.hpp"
#include "cosmosim/gravity/tree_pm_split_kernel.hpp"

namespace cosmosim::gravity {

// Shared compact active-set view for force accumulation ownership.
struct TreePmForceAccumulatorView {
  // By default each active index identifies a local source particle and its
  // coordinates are read from the source arrays. When all three explicit
  // target-position spans are present, they are authoritative instead;
  // UINT32_MAX then denotes a target with no local source/self identity.
  std::span<const std::uint32_t> active_particle_index;
  std::span<double> accel_x_comoving;
  std::span<double> accel_y_comoving;
  std::span<double> accel_z_comoving;
  // Optional previous total-acceleration magnitude used by the relative
  // force-error MAC. Missing or non-finite entries select the deterministic
  // COM-distance fallback for that target.
  std::span<const double> previous_acceleration_magnitude_code{};
  std::span<const double> target_pos_x_comoving{};
  std::span<const double> target_pos_y_comoving{};
  std::span<const double> target_pos_z_comoving{};

  void reset() const;
  void addToActiveSlot(std::size_t active_slot, double ax_comoving, double ay_comoving, double az_comoving) const;
};

struct TreePmOptions {
  PmSolveOptions pm_options{};
  TreeGravityOptions tree_options{};
  TreePmSplitPolicy split_policy{};
  bool enable_zoom_long_range_correction = false;
  PmGridShape zoom_focused_pm_shape{};
  std::span<const std::uint8_t> source_is_high_res;
  std::span<const std::uint8_t> active_is_high_res;
  double zoom_region_center_x_comoving = 0.0;
  double zoom_region_center_y_comoving = 0.0;
  double zoom_region_center_z_comoving = 0.0;
  double zoom_region_radius_comoving = 0.0;
  double zoom_contamination_radius_comoving = 0.0;
  // Runtime identity carried by every distributed short-range request and
  // response. The workflow owns these epochs; the coordinator owns only its
  // per-instance exchange sequence.
  std::uint64_t decomposition_epoch = 0;
  std::uint64_t force_epoch = 0;
  std::uint64_t tree_exchange_batch_bytes = 4ULL * 1024ULL * 1024ULL;
  std::uint64_t zoom_high_res_allgather_limit_bytes = 256ULL * 1024ULL * 1024ULL;
};

struct TreePmDiagnostics {
  std::uint64_t local_source_count = 0;
  std::uint64_t local_active_target_count = 0;
  std::uint64_t empty_source_rank_count = 0;
  std::uint64_t empty_target_rank_count = 0;
  std::uint64_t local_tree_node_count = 0;
  std::uint64_t remote_hierarchy_packet_count = 0;
  std::uint64_t communicating_peer_count = 0;
  std::uint64_t pm_solve_count = 0;
  std::uint64_t pm_reuse_count = 0;
  std::uint64_t pm_halo_value_count = 0;
  std::uint64_t pm_local_nx = 0;
  std::uint64_t pm_local_ny = 0;
  std::uint64_t pm_local_nz = 0;
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
  std::uint64_t residual_remote_request_packets = 0;
  std::uint64_t residual_remote_response_packets = 0;
  std::uint64_t residual_remote_request_bytes = 0;
  std::uint64_t residual_remote_response_bytes = 0;
  std::uint64_t residual_remote_request_batches = 0;
  std::uint64_t residual_remote_peer_participations = 0;
  std::uint64_t residual_remote_targets_with_requests = 0;
  std::uint64_t residual_remote_targets_without_requests = 0;
  std::uint64_t residual_remote_pairs_pruned_by_bounds = 0;
  std::uint64_t residual_remote_request_packets_max_peer = 0;
  std::uint64_t residual_remote_response_packets_max_peer = 0;
  double residual_remote_request_packet_imbalance_ratio = 0.0;
  std::uint64_t zoom_high_res_source_count = 0;
  std::uint64_t zoom_low_res_source_count = 0;
  std::uint64_t zoom_low_res_contamination_count = 0;
  std::uint64_t zoom_high_res_allgather_bytes = 0;
  std::uint64_t zoom_high_res_allgather_limit_bytes = 0;
  double zoom_low_res_contamination_mass_code = 0.0;
  double force_l2_pm_global = 0.0;
  double force_l2_pm_zoom_correction = 0.0;
  double force_l2_tree_short_range = 0.0;
  double force_l2_tree_short_range_local = 0.0;
  double force_l2_tree_short_range_remote = 0.0;
  double force_l2_total = 0.0;
  // Periodic TreePM builds its transient tree in one contiguous, seam-safe
  // unwrapped frame per axis. These root diagnostics make that geometry
  // contract directly testable without exposing mutable tree storage.
  double tree_root_half_size_comoving = 0.0;
  double tree_root_com_x_comoving = 0.0;
  double tree_root_com_y_comoving = 0.0;
  double tree_root_com_z_comoving = 0.0;
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
  TreePmCoordinator(PmGridShape pm_shape, parallel::PmSlabLayout pm_layout);
  TreePmCoordinator(PmGridShape pm_shape, parallel::PmSlabLayout pm_layout, parallel::MpiContext mpi_context);

  [[nodiscard]] const parallel::PmSlabLayout& slabLayout() const noexcept;
  [[nodiscard]] bool ownsFullPmDomain() const noexcept;
  [[nodiscard]] const parallel::PmSlabHaloExchangeResult& lastPmSlabHaloExchange() const noexcept;
  [[nodiscard]] core::MemoryReport memoryReport() const;

  void solveActiveSetWithPmCadence(
      std::span<const double> pos_x_comoving,
      std::span<const double> pos_y_comoving,
      std::span<const double> pos_z_comoving,
      std::span<const double> mass_code,
      const TreePmForceAccumulatorView& accumulator,
      const TreePmOptions& options,
      bool refresh_long_range_field,
      TreePmProfileEvent* profile = nullptr,
      TreePmDiagnostics* diagnostics = nullptr,
      const TreeSofteningView& softening_view = {});

  void solveActiveSet(
      std::span<const double> pos_x_comoving,
      std::span<const double> pos_y_comoving,
      std::span<const double> pos_z_comoving,
      std::span<const double> mass_code,
      const TreePmForceAccumulatorView& accumulator,
      const TreePmOptions& options,
      TreePmProfileEvent* profile = nullptr,
      TreePmDiagnostics* diagnostics = nullptr,
      const TreeSofteningView& softening_view = {});

 private:
  void evaluateShortRangeResidual(
      std::span<const double> pos_x_comoving,
      std::span<const double> pos_y_comoving,
      std::span<const double> pos_z_comoving,
      std::span<const double> mass_code,
      const TreePmForceAccumulatorView& accumulator,
      const TreePmOptions& options,
      const TreeSofteningView& softening_view,
      bool rank_local_serial_mode,
      TreeGravityProfile* tree_profile);

  struct ResidualTraversalStats {
    std::uint64_t pruned_nodes = 0;
    std::uint64_t pair_skips_cutoff = 0;
    std::uint64_t pair_evaluations = 0;
    std::uint64_t remote_request_packets = 0;
    std::uint64_t remote_response_packets = 0;
    std::uint64_t remote_hierarchy_packets = 0;
    std::uint64_t communicating_peer_count = 0;
    std::uint64_t remote_request_bytes = 0;
    std::uint64_t remote_response_bytes = 0;
    std::uint64_t remote_request_batches = 0;
    std::uint64_t remote_peer_participations = 0;
    std::uint64_t remote_targets_with_requests = 0;
    std::uint64_t remote_targets_without_requests = 0;
    std::uint64_t remote_pairs_pruned_by_bounds = 0;
    std::uint64_t remote_request_packets_max_peer = 0;
    std::uint64_t remote_response_packets_max_peer = 0;
    double remote_request_packet_imbalance_ratio = 0.0;
    double local_short_range_sum_sq = 0.0;
    double remote_short_range_sum_sq = 0.0;
  };

  PmGridShape m_shape;
  parallel::MpiContext m_mpi_context;
  PmGridStorage m_grid;
  PmSolver m_pm_solver;
  TreeGravitySolver m_tree_solver;

  // Periodic tree-build coordinates are transient derived state. PM assignment
  // and particle truth continue to use the caller-owned wrapped coordinates.
  std::vector<double> m_tree_source_x_comoving;
  std::vector<double> m_tree_source_y_comoving;
  std::vector<double> m_tree_source_z_comoving;

  // Reused compact-sidecar scratch arrays to avoid per-step heap churn.
  std::vector<double> m_active_pos_x_comoving;
  std::vector<double> m_active_pos_y_comoving;
  std::vector<double> m_active_pos_z_comoving;
  std::vector<double> m_active_pm_ax_comoving;
  std::vector<double> m_active_pm_ay_comoving;
  std::vector<double> m_active_pm_az_comoving;
  std::vector<double> m_active_zoom_corr_ax_comoving;
  std::vector<double> m_active_zoom_corr_ay_comoving;
  std::vector<double> m_active_zoom_corr_az_comoving;
  struct TreeExchangeWorkspace {
    int world_size = 1;
    std::vector<int> send_counts;
    std::vector<int> recv_counts;
    std::vector<int> send_displs;
    std::vector<int> recv_displs;
    std::vector<int> response_send_counts;
    std::vector<int> response_recv_counts;
    std::vector<int> response_send_displs;
    std::vector<int> response_recv_displs;
    std::vector<std::uint8_t> send_payload;
    std::vector<std::uint8_t> recv_payload;
    std::vector<std::uint8_t> response_send_payload;
    std::vector<std::uint8_t> response_recv_payload;
    std::vector<double> remote_batch_ax;
    std::vector<double> remote_batch_ay;
    std::vector<double> remote_batch_az;
    std::vector<std::uint32_t> expected_response_count;
    std::vector<std::uint32_t> received_response_count;
  } m_tree_exchange_workspace;
  ResidualTraversalStats m_last_residual_stats;
  parallel::PmSlabHaloExchangeResult m_last_pm_slab_halo_exchange{};
  std::uint64_t m_pm_halo_exchange_sequence = 0;
  std::uint64_t m_tree_exchange_sequence = 0;
  struct LongRangeFieldValidity {
    bool valid = false;
    std::uint64_t decomposition_epoch = 0;
    std::uint64_t force_epoch = 0;
    double scale_factor = 0.0;
    double gravitational_constant_code = 0.0;
    double split_scale_comoving = 0.0;
    double box_size_x_comoving = 0.0;
    double box_size_y_comoving = 0.0;
    double box_size_z_comoving = 0.0;
    PmAssignmentScheme assignment_scheme = PmAssignmentScheme::kCic;
    PmBoundaryCondition boundary_condition = PmBoundaryCondition::kPeriodic;
    core::PmDecompositionMode decomposition_mode = core::PmDecompositionMode::kSlab;
    bool window_deconvolution = false;
  } m_long_range_field_validity;
};

[[nodiscard]] TreePmDiagnostics computeTreePmDiagnostics(const TreePmSplitPolicy& split_policy);

}  // namespace cosmosim::gravity
