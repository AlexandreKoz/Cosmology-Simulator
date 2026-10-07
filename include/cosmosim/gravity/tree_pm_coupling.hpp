#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <span>
#include <vector>

#include "cosmosim/core/memory_accounting.hpp"
#include "cosmosim/gravity/gravity_communication_arena.hpp"
#include "cosmosim/gravity/gravity_state_identity.hpp"
#include "cosmosim/gravity/pm_solver.hpp"
#include "cosmosim/gravity/tree_gravity.hpp"
#include "cosmosim/gravity/tree_pm_split_kernel.hpp"

namespace cosmosim::gravity {

// Short-range exchange wire record widths. The coordinator's host packet
// structs are layout-compatible with these widths on supported ABIs; that
// identity is compile-time enforced by static_assert next to the packet
// definitions in tree_pm_coupling.cpp, so a padding change on a supported ABI
// fails the build instead of silently letting wire bytes stand in for host
// object bytes. Gravity memory estimation reuses the same constants so
// protocol and budget arithmetic cannot drift apart.
inline constexpr std::size_t kTreePmShortRangeRequestWireBytes = 96U;
inline constexpr std::size_t kTreePmShortRangeResponseWireBytes = 80U;
inline constexpr std::size_t kTreePmResidualBlockSize = 64U;
// One bundle per worker; no pair/node hot-path atomics or shared cache lines.
struct alignas(64) TreePmTraversalCounters {
  std::uint64_t visited_nodes = 0;
  std::uint64_t accepted_nodes = 0;  // compatibility: leaves + internal nodes
  std::uint64_t opened_nodes = 0;
  std::uint64_t direct_pair_evaluations = 0;
  std::uint64_t cutoff_pruned_nodes = 0;
  std::uint64_t cutoff_skipped_pairs = 0;
  std::uint64_t remote_pairs_pruned_by_bounds = 0;
  std::uint64_t accepted_internal_multipoles = 0;
  std::uint64_t accepted_leaves = 0;
  std::uint64_t selected_mac_rejections = 0;
  std::uint64_t relative_mac_rejections = 0;
  std::uint64_t maximum_angle_rejections = 0;
  std::uint64_t strict_envelope_rejections = 0;
  std::uint64_t softening_rejections = 0;
  std::uint64_t near_node_rejections = 0;
  std::uint64_t cutoff_containment_rejections = 0;
  std::uint64_t geometric_history_fallbacks = 0;
  std::uint64_t targets = 0;
  // Summed logical-block work time; no per-target clock reads.
  double elapsed_work_ms = 0.0;
  double block_work_ms_max = 0.0;
};
inline constexpr std::size_t kTreePmResidualCounterBytes = sizeof(TreePmTraversalCounters);
// Optional worker scratch, separate from the aligned hot counter bundle.
inline constexpr std::size_t kTreePmResidualSpatialCounterBytes =
    2U * sizeof(std::array<std::uint64_t, parallel::k_spatial_work_bin_count>);

// Why TreePM declined the installed authoritative top-domain geometry and
// used the conservative local-tree root packet instead.
enum class TreePmDomainGeometryFallbackReason : std::uint32_t {
  kNone = 0,
  kNoGeometryInstalled = 1,
  kStaleSourceGeneration = 2,
  kDecompositionEpochMismatch = 3,
  kSourceCoverageFailure = 4,
  kGeometryPreparationFailure = 5,
  kCollectivePeerFallback = 6,
};

[[nodiscard]] inline const char* treePmDomainGeometryFallbackReasonName(
    TreePmDomainGeometryFallbackReason reason) noexcept {
  switch (reason) {
    case TreePmDomainGeometryFallbackReason::kNone:
      return "none";
    case TreePmDomainGeometryFallbackReason::kNoGeometryInstalled:
      return "no_geometry_installed";
    case TreePmDomainGeometryFallbackReason::kStaleSourceGeneration:
      return "stale_source_generation";
    case TreePmDomainGeometryFallbackReason::kDecompositionEpochMismatch:
      return "decomposition_epoch_mismatch";
    case TreePmDomainGeometryFallbackReason::kSourceCoverageFailure:
      return "source_coverage_failure";
    case TreePmDomainGeometryFallbackReason::kGeometryPreparationFailure:
      return "geometry_preparation_failure";
    case TreePmDomainGeometryFallbackReason::kCollectivePeerFallback:
      return "collective_peer_fallback";
  }
  return "unknown";
}

// Shared compact active-set view for force accumulation ownership.
struct TreePmForceAccumulatorView {
  // By default each active index identifies a local source particle and its
  // coordinates are read from the source arrays. When all three explicit
  // target-position spans are present, they are authoritative instead;
  // UINT32_MAX then denotes a target with no local source/self identity.
  std::span<const TreeLocalIndex> active_particle_index;
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
  // Optional compact PM-only output, captured before the residual is added.
  // All three spans must cover the active set. Empty for ordinary KDK.
  std::span<double> long_range_accel_x_comoving{};
  std::span<double> long_range_accel_y_comoving{};
  std::span<double> long_range_accel_z_comoving{};
  // Optional source-indexed Tree-only lanes. Each active source row is reset
  // and accumulated directly, avoiding subtraction of nearly equal PM/total
  // values. Requires indexed local source targets, not independent positions.
  std::span<double> short_range_accel_x_comoving{};
  std::span<double> short_range_accel_y_comoving{};
  std::span<double> short_range_accel_z_comoving{};

  void reset() const;
  void addToActiveSlot(std::size_t active_slot, double ax_comoving, double ay_comoving, double az_comoving) const;
  void addShortRangeToActiveSlot(std::size_t active_slot, double ax_comoving, double ay_comoving, double az_comoving) const;
};

enum class TreePmAcceptancePolicy : std::uint8_t {
  kStrictReference,
  kAdaptiveRelative,
};

struct TreePmOptions {
  PmSolveOptions pm_options{};
  TreeGravityOptions tree_options{};
  TreePmSplitPolicy split_policy{};
  TreePmAcceptancePolicy acceptance_policy = TreePmAcceptancePolicy::kStrictReference;
  double adaptive_maximum_opening_angle = 0.25;
  bool spatial_work_history_enabled = false;
  // Exact snapshot reuse and certified motion refit are independent policies.
  bool identical_source_tree_reuse_enabled = false;
  bool topology_refit_enabled = false;
  // Nonzero certifies unchanged dense source-row identity. Rank-local token,
  // not a collective generation; unknown/reordered views reject motion refit.
  std::uint64_t source_layout_generation = 0U;
  // Integrator split operator: evaluate ONLY the complementary short-range
  // force, without interpolating or reusing a stale PM field. Long-range kicks
  // have separate synchronized endpoint authority in hierarchical KDK.
  bool short_range_only = false;
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
  DecompositionEpoch decomposition_epoch{};
  // Authoritative compact top-domain geometry owned by src/parallel. Empty is
  // an explicit lower-capability fallback to the derived local tree-root
  // envelope used by older/reduced workflow paths.
  std::span<const parallel::TopDomainLeaf> authoritative_domain_leaves;
  // Monotonic identity of the physical source snapshot represented by a PM
  // field. Reuse is legal only when this generation still matches. The
  // workflow owns this token; it must advance when source positions/masses or
  // source membership change.
  GravitySourceGeneration source_generation{};
  // Source generation the installed authoritative top-domain leaves were
  // refit against. TreePM uses authoritative geometry only when this equals
  // source_generation; otherwise it falls back conservatively and reports
  // TreePmDomainGeometryFallbackReason::kStaleSourceGeneration.
  GravitySourceGeneration authoritative_geometry_source_generation{};
  PmFieldVersion pm_field_version{};
  ForceEvaluationEpoch force_epoch{};
  std::uint64_t tree_exchange_batch_bytes = 4ULL * 1024ULL * 1024ULL;
  std::uint64_t zoom_high_res_allgather_limit_bytes = 256ULL * 1024ULL * 1024ULL;
};

struct TreePmDiagnostics {
  TreePmTraversalCounters local_traversal{};
  TreePmTraversalCounters incoming_traversal{};
  std::array<double, parallel::k_spatial_work_bin_count> spatial_work_per_target{};
  std::uint64_t spatial_work_history_solves = 0;
  // Nonempty worker-region summaries for the current force only. Work time
  // is summed logical-block time, not an additive force-phase wall timer.
  std::uint64_t worker_region_count = 0;
  double worker_targets_min = 0.0;
  double worker_targets_max = 0.0;
  double worker_visits_max = 0.0;
  double worker_pairs_max = 0.0;
  double worker_multipoles_max = 0.0;
  double worker_work_ms_min = 0.0;
  double worker_work_ms_max = 0.0;
  double worker_work_ms_sum = 0.0;
  std::uint64_t local_source_count = 0;
  std::uint64_t local_active_target_count = 0;
  std::uint64_t empty_source_rank_count = 0;
  std::uint64_t empty_target_rank_count = 0;
  std::uint64_t local_tree_node_count = 0;
  std::uint64_t remote_hierarchy_packet_count = 0;
  std::uint64_t communicating_peer_count = 0;
  std::uint64_t top_level_domain_leaf_count = 0;
  std::uint64_t authoritative_domain_leaf_count = 0;
  // Authoritative top-domain geometry lifecycle for this force solve.
  std::uint64_t domain_geometry_source_generation = 0;
  std::uint64_t current_gravity_source_generation = 0;
  std::uint64_t domain_geometry_fresh = 0;
  std::uint64_t domain_geometry_fallback_used = 0;
  std::uint64_t domain_geometry_fallback_reason = 0;
  std::uint64_t domain_geometry_uncovered_source_count = 0;
  std::uint64_t domain_hierarchy_node_count = 0;
  std::uint64_t domain_cache_hit_count = 0;
  std::uint64_t domain_cache_miss_count = 0;
  std::uint64_t graph_cache_hit_count = 0;
  std::uint64_t graph_cache_miss_count = 0;
  std::uint64_t let_candidate_peer_count = 0;
  std::uint64_t let_exported_target_count = 0;
  std::uint64_t let_imported_target_count = 0;
  std::uint64_t let_wire_bytes_sent = 0;
  std::uint64_t let_wire_bytes_received = 0;
  // Compatibility alias: four-wire-buffer high-water only.
  std::uint64_t let_high_water_bytes = 0;
  std::uint64_t let_wire_buffer_high_water_bytes = 0;
  std::uint64_t let_known_workspace_high_water_bytes = 0;
  std::uint64_t communication_arena_capacity_bytes = 0;
  std::uint64_t communication_arena_logical_high_water_bytes = 0;
  double exported_targets_per_requested_target = 0.0;
  double let_discovery_ms = 0.0;
  double let_graph_setup_ms = 0.0;
  double let_communication_ms = 0.0;
  double let_overlap_local_work_ms = 0.0;
  double let_communication_wait_ms = 0.0;
  double let_overlap_efficiency = 0.0;
  // Incoming short-range remote-phase timer split. The compute timer
  // surrounds only validated incoming target force evaluation against this
  // rank's local tree; request decode/validation and response encode/pack are
  // separate timers, as are consensus and response exchange. The compatibility
  // alias let_remote_traversal_ms carries incoming compute only and never
  // decode+validation+hash+tree+encode as "traversal."
  double let_remote_traversal_ms = 0.0;
  double incoming_request_decode_validation_ms = 0.0;
  double incoming_remote_target_compute_ms = 0.0;
  double incoming_response_encode_pack_ms = 0.0;
  // Response count/displacement layout and response payload buffer sizing.
  double protocol_validation_ms = 0.0;
  double protocol_consensus_ms = 0.0;
  double response_exchange_ms = 0.0;
  std::uint64_t pm_solve_count = 0;
  std::uint64_t pm_reuse_count = 0;
  // Logical left+right halo values per scalar force component (not XYZ total).
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
  // Exact pair-evaluation truth: residual_pair_evaluations =
  // local_pair_evaluations + incoming_remote_pair_evaluations. Local counts
  // targets owned by this rank; incoming counts remote targets evaluated
  // against this rank's tree. tree_profile.particle_particle_interactions
  // remains the combined residual total for downstream profile consumers.
  std::uint64_t residual_pair_evaluations = 0;
  std::uint64_t local_pair_evaluations = 0;
  std::uint64_t incoming_remote_pair_evaluations = 0;
  std::uint64_t residual_remote_request_packets = 0;
  std::uint64_t residual_remote_response_packets = 0;
  std::uint64_t residual_remote_request_bytes = 0;
  std::uint64_t residual_remote_response_bytes = 0;
  std::uint64_t residual_remote_request_batches = 0;
  std::uint64_t residual_remote_peer_participations = 0;
  std::uint64_t residual_remote_targets_with_requests = 0;
  std::uint64_t residual_remote_targets_without_requests = 0;
  std::uint64_t residual_remote_pairs_pruned_by_bounds = 0;
  // OpenMP residual execution provenance for this solve. Scratch high-water is
  // the retained contiguous worker-stack capacity, not a per-target estimate.
  std::uint64_t openmp_compiled = 0;
  std::uint64_t openmp_configured_workers = 0;
  std::uint64_t openmp_observed_workers = 0;
  std::uint64_t residual_local_target_count = 0;
  std::uint64_t residual_incoming_target_count = 0;
  std::uint64_t residual_worker_scratch_high_water_bytes = 0;
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
  double source_preprocess_ms = 0.0;
  double let_discovery_ms = 0.0;
  double let_graph_setup_ms = 0.0;
  double let_communication_ms = 0.0;
  double let_overlap_local_work_ms = 0.0;
  double let_communication_wait_ms = 0.0;
  double let_overlap_efficiency = 0.0;
  double remote_traversal_ms = 0.0;
  double incoming_request_decode_validation_ms = 0.0;
  double incoming_remote_target_compute_ms = 0.0;
  double incoming_response_encode_pack_ms = 0.0;
  double protocol_validation_ms = 0.0;
  double protocol_consensus_ms = 0.0;
  double response_exchange_ms = 0.0;
};

// Thin coordinator that makes TreePM ownership explicit and auditable.
class TreePmCoordinator {
 public:
  explicit TreePmCoordinator(PmGridShape pm_shape);
  TreePmCoordinator(PmGridShape pm_shape, parallel::PmSlabLayout pm_layout);
  TreePmCoordinator(
      PmGridShape pm_shape,
      parallel::PmSlabLayout pm_layout,
      parallel::MpiContext mpi_context,
      core::MemoryGovernor* memory_governor = nullptr);
  ~TreePmCoordinator();

  [[nodiscard]] const parallel::PmSlabLayout& slabLayout() const noexcept;
  [[nodiscard]] bool ownsFullPmDomain() const noexcept;
  [[nodiscard]] const parallel::PmSlabHaloExchangeResult& lastPmSlabHaloExchange() const noexcept;
  [[nodiscard]] core::MemoryReport memoryReport() const;
  // Explicitly release cached MPI-owned communicators while the MPI session is active.
  // Normal application/benchmark lifecycle must call this before MPI_Finalize; the
  // destructor retains only a noexcept emergency fallback.
  void shutdownMpiResources();

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
  // Per-target-family residual traversal bundle. The coordinator maintains
  // one instance for locally owned targets and one for incoming remote
  // targets so pair evaluations can be reported without double counting.
  using ResidualTraversalCounters = TreePmTraversalCounters;
   static_assert(sizeof(ResidualTraversalCounters) == kTreePmResidualCounterBytes,
                 "TreePM residual counter storage contract changed");

  // Prepared once before traversal only when spatial feedback is enabled.
  // Owned targets accumulate here across batches; incoming targets never do.
  struct alignas(64) ResidualSpatialWorkCounters {
    std::array<std::uint64_t, parallel::k_spatial_work_bin_count> spatial_work{};
    std::array<std::uint64_t, parallel::k_spatial_work_bin_count> spatial_targets{};
  };
  static_assert(sizeof(ResidualSpatialWorkCounters) == kTreePmResidualSpatialCounterBytes,
                "TreePM spatial worker scratch storage contract changed");

  void evaluateShortRangeResidual(
      std::span<const double> pos_x_comoving,
      std::span<const double> pos_y_comoving,
      std::span<const double> pos_z_comoving,
      std::span<const double> mass_code,
      const TreePmForceAccumulatorView& accumulator,
      const TreePmOptions& options,
       const TreeSofteningView& softening_view,
       bool distributed_payload_communication_required,
       TreeGravityProfile* tree_profile);


  struct ResidualTraversalStats {
    TreePmDiagnostics traversal_summary{};
    std::uint64_t pruned_nodes = 0;
    std::uint64_t pair_skips_cutoff = 0;
    // Exact identity: pair_evaluations =
    // local_pair_evaluations + incoming_remote_pair_evaluations.
    std::uint64_t pair_evaluations = 0;
    std::uint64_t local_pair_evaluations = 0;
    std::uint64_t incoming_remote_pair_evaluations = 0;
    std::uint64_t remote_request_packets = 0;
    std::uint64_t remote_response_packets = 0;
    std::uint64_t remote_hierarchy_packets = 0;
    std::uint64_t communicating_peer_count = 0;
    std::uint64_t top_level_domain_leaf_count = 0;
     std::uint64_t authoritative_domain_leaf_count = 0;
     std::uint64_t domain_geometry_source_generation = 0;
     std::uint64_t current_gravity_source_generation = 0;
     std::uint64_t domain_geometry_fresh = 0;
     std::uint64_t domain_geometry_fallback_used = 0;
     std::uint64_t domain_geometry_fallback_reason = 0;
     std::uint64_t domain_geometry_uncovered_source_count = 0;
     std::uint64_t domain_hierarchy_node_count = 0;

    std::uint64_t domain_cache_hit_count = 0;
    std::uint64_t domain_cache_miss_count = 0;
    std::uint64_t graph_cache_hit_count = 0;
    std::uint64_t graph_cache_miss_count = 0;
    std::uint64_t let_candidate_peer_count = 0;
    std::uint64_t let_exported_target_count = 0;
    std::uint64_t let_imported_target_count = 0;
    std::uint64_t let_wire_bytes_sent = 0;
    std::uint64_t let_wire_bytes_received = 0;
    // Compatibility alias: four-wire-buffer high-water only.
    std::uint64_t let_high_water_bytes = 0;
    std::uint64_t let_wire_buffer_high_water_bytes = 0;
    std::uint64_t let_known_workspace_high_water_bytes = 0;
    double let_discovery_ms = 0.0;
    double let_graph_setup_ms = 0.0;
    double let_communication_ms = 0.0;
    double let_overlap_local_work_ms = 0.0;
    double let_communication_wait_ms = 0.0;
    double let_overlap_efficiency = 0.0;
    double let_remote_traversal_ms = 0.0;
    double incoming_request_decode_validation_ms = 0.0;
    double incoming_remote_target_compute_ms = 0.0;
    double incoming_response_encode_pack_ms = 0.0;
    double protocol_validation_ms = 0.0;
    double protocol_consensus_ms = 0.0;
    double response_exchange_ms = 0.0;
    std::uint64_t remote_request_bytes = 0;
    std::uint64_t remote_response_bytes = 0;
    std::uint64_t remote_request_batches = 0;
    std::uint64_t remote_peer_participations = 0;
    std::uint64_t remote_targets_with_requests = 0;
    std::uint64_t remote_targets_without_requests = 0;
    std::uint64_t remote_pairs_pruned_by_bounds = 0;
    std::uint64_t incoming_remote_target_evaluations = 0;
    std::uint64_t openmp_observed_workers = 0;
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
  // Bounded derived planner feedback; survives ownership changes because bins
  // address fixed coordinates. Never force truth or restart continuation state.
  std::array<double, parallel::k_spatial_work_bin_count> m_spatial_work_history{};
  std::array<double, parallel::k_spatial_work_bin_count> m_spatial_target_history{};
  std::uint64_t m_spatial_work_history_solves = 0;
  DecompositionEpoch m_tree_decomposition_epoch{};
  PmBoundaryCondition m_tree_boundary = PmBoundaryCondition::kPeriodic;
  std::array<double, 3U> m_tree_box_lengths{};
  std::array<double, 3U> m_tree_unwrap_anchor{};
  std::uint64_t m_tree_source_layout_generation = 0U;

  // Periodic tree-build coordinates are transient derived state. PM assignment
  // and particle truth continue to use the caller-owned wrapped coordinates.
  std::vector<double> m_tree_source_x_comoving;
  std::vector<double> m_tree_source_y_comoving;
  std::vector<double> m_tree_source_z_comoving;

  // Zoom correction is optional and owns compact lanes only while the selected
  // profile enables it. Ordinary source-index targets alias authoritative
  // source coordinates and PM writes directly into the active force buffer.
  std::vector<double> m_active_zoom_corr_ax_comoving;
  std::vector<double> m_active_zoom_corr_ay_comoving;
  std::vector<double> m_active_zoom_corr_az_comoving;
  std::array<std::uint64_t, 3> m_zoom_corr_high_water_bytes{};
  // One physical rank-local backing allocation reused by the mutually
  // exclusive PM density, halo, PM interpolation, and short-range Tree
  // communication phases. Logical protocol storage is phase-local and
  // arena-backed; persistent topology/cache state remains outside it.
  GravityCommunicationArena m_communication_arena;
  std::uint64_t m_tree_exchange_logical_high_water_bytes = 0U;

  // Contiguous OpenMP residual DFS slots: worker_count * (1 + 7 * max_depth)
  // TreeLocalIndex entries. Retained between residual evaluations so the
  // MemoryGovernor sees one stable high-water instead of per-call growth.
  std::vector<TreeLocalIndex> m_worker_stack_storage;
  std::uint64_t m_worker_stack_high_water_bytes = 0;
  std::vector<ResidualTraversalCounters> m_worker_counter_storage;
  std::uint64_t m_worker_counter_high_water_bytes = 0;
  std::vector<ResidualSpatialWorkCounters> m_worker_spatial_counter_storage;
  std::uint64_t m_worker_spatial_counter_high_water_bytes = 0;
  std::vector<double> m_block_sum_sq_storage;
  std::uint64_t m_block_sum_sq_high_water_bytes = 0;
  ResidualTraversalStats m_last_residual_stats;
  parallel::PmSlabHaloExchangeResult m_last_pm_slab_halo_exchange{};
  std::uint64_t m_pm_halo_exchange_sequence = 0;
  std::uint64_t m_tree_exchange_sequence = 0;
  struct LetDomainCache {
    bool valid = false;
    DecompositionEpoch decomposition_epoch{};
    int world_size = 1;
    std::uint64_t local_geometry_fingerprint = 0U;
    std::uint64_t geometry_fingerprint = 0U;
    bool authoritative_geometry = false;
    std::vector<parallel::TreePseudoParticlePacket> top_level_domain_leaves;
  } m_let_domain_cache;
  struct LetDomainHierarchyCacheOpaque;
  std::unique_ptr<LetDomainHierarchyCacheOpaque> m_let_domain_hierarchy_cache;
  struct SparsePeerGraphCacheOpaque;
  std::unique_ptr<SparsePeerGraphCacheOpaque> m_sparse_peer_graph_cache;
  struct LongRangeFieldValidity {
    bool valid = false;
    DecompositionEpoch decomposition_epoch{};
    GravitySourceGeneration source_generation{};
    PmFieldVersion pm_field_version{};
    ForceEvaluationEpoch last_force_epoch{};
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
