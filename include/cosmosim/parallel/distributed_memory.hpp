#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <iosfwd>
#include <exception>
#include <limits>
#include <stdexcept>
#include <span>
#include <string>
#include <unordered_map>
#include <vector>

#include "cosmosim/core/checked_arithmetic.hpp"
#include "cosmosim/core/profiling.hpp"
#include "cosmosim/parallel/distributed_mesh.hpp"

namespace cosmosim::core {

class SimulationState;

}  // namespace cosmosim::core

namespace cosmosim::parallel {

struct DecompositionConfig;

enum class DecompositionEntityKind : std::uint8_t {
  kParticle = 0,
  kHydroCell = 1,
  kAmrPatch = 2,
  kPmMeshCell = 3,
};

struct DecompositionWorkComponents {
  double particle_count_cost = 0.0;
  double gas_cell_cost = 0.0;
  double tree_interaction_cost = 0.0;
  double pm_mesh_cost = 0.0;
  double amr_patch_cost = 0.0;
  double active_fraction_cost = 0.0;
  double memory_pressure_cost = 0.0;
  double transient_memory_cost = 0.0;
  double source_event_cost = 0.0;
  double communication_cost = 0.0;
  double gpu_occupancy_cost = 0.0;
  double generic_work_cost = 0.0;
  bool has_explicit_components = false;

  [[nodiscard]] double rawTotal() const noexcept;
};

struct DecompositionWeightCoefficients {
  double particle_count = 1.0;
  double gas_cell = 1.0;
  double tree_interaction = 1.0;
  double pm_mesh = 1.0;
  double amr_patch = 1.0;
  double active_fraction = 1.0;
  double memory_pressure = 1.0;
  double gpu_occupancy = 0.0;
  double generic_work = 1.0;
};

struct DecompositionItem {
  std::uint64_t entity_id = 0;
  DecompositionEntityKind kind = DecompositionEntityKind::kParticle;
  int current_owner_rank = -1;
  double x_comov = 0.0;
  double y_comov = 0.0;
  double z_comov = 0.0;
  // Optional conservative spatial footprint of this decomposition unit. Point
  // particles use a degenerate box. AMR patches may span many authoritative
  // gas-cell centres and therefore publish the complete owned footprint. The
  // decomposition layer, not gravity, owns this geometry contract.
  bool has_spatial_bounds = false;
  double min_x_comov = 0.0;
  double max_x_comov = 0.0;
  double min_y_comov = 0.0;
  double max_y_comov = 0.0;
  double min_z_comov = 0.0;
  double max_z_comov = 0.0;
  std::uint64_t active_target_count_recent = 0;
  std::uint64_t remote_tree_interactions_recent = 0;
  double work_units = 1.0;
  std::uint64_t memory_bytes = 0;
  DecompositionWorkComponents work_components{};
};

// Compact production planner record. Replaces the former dual N-sized
// DecompositionItem + LocalKeyedItem populations during cut planning, cut
// samples, migration-intent emission, and streaming metric recompute.
// Geometry and full work_components remain in the source DecompositionItem;
// this record carries only SFC order, ownership, load, and index fields.
// local_index has two documented contracts: source-view streaming records
// (makeSourceRecords) carry the canonical kind-specific source row (particle
// rows and included patch rows, unique within their kind), while the rich
// compatibility adapter and the startup record stream carry a combined
// ordinal over their local record population. The shared compact planner
// core (buildMortonSfcDecompositionFromCompact) accepts only the combined
// contract and enforces {local_index} as a permutation of [0, record count).
struct CompactRuntimeDecompositionRecord {
  std::uint64_t entity_id = 0;
  std::uint64_t sfc_key = 0;
  std::uint64_t memory_bytes = 0;
  double weighted_load = 0.0;
  std::size_t local_index = 0;
  std::uint64_t active_target_count_recent = 0;
  std::uint64_t remote_tree_interactions_recent = 0;
  int current_owner_rank = -1;
  DecompositionEntityKind kind = DecompositionEntityKind::kParticle;
};
static_assert(sizeof(CompactRuntimeDecompositionRecord) <= 64U,
              "CompactRuntimeDecompositionRecord must fit in 64 bytes");

struct RuntimeDecompositionSourceView {
  int world_rank = 0;
  std::size_t particle_count = 0U;
  std::size_t patch_count = 0U;
  std::span<const std::uint64_t> particle_ids{};
  std::span<const double> particle_x_comoving{};
  std::span<const double> particle_y_comoving{};
  std::span<const double> particle_z_comoving{};
  std::span<const std::uint32_t> particle_species_tag{};
  std::span<const std::uint32_t> particle_owning_rank{};
  std::span<const std::uint8_t> active_particle_mask{};
  std::span<const std::uint64_t> patch_ids{};
  std::span<const std::int32_t> patch_levels{};
  std::span<const std::uint32_t> patch_owning_rank{};
  std::span<const std::uint32_t> patch_first_cells{};
  std::span<const std::uint32_t> patch_cell_counts{};
  std::span<const std::uint32_t> compact_patch_indices{};
  std::span<const std::uint16_t> patch_cell_dim_x{};
  std::span<const std::uint16_t> patch_cell_dim_y{};
  std::span<const std::uint16_t> patch_cell_dim_z{};
  std::span<const double> cell_x_comoving{};
  std::span<const double> cell_y_comoving{};
  std::span<const double> cell_z_comoving{};
  std::span<const std::uint32_t> cell_patch_indices{};
  std::span<const std::uint32_t> gas_patch_list_offsets{};
  std::span<const std::uint32_t> gas_patch_indices{};
  std::array<std::uint64_t, 5U> particle_memory_bytes_by_species{};
  std::uint32_t gas_species_tag = 1U;
  std::uint32_t star_species_tag = 2U;
  std::uint32_t black_hole_species_tag = 3U;
  std::uint64_t gas_transient_memory_bytes_per_cell = 0U;
  std::uint64_t source_scratch_bytes = 0U;

  [[nodiscard]] std::size_t localEntityCount() const noexcept {
    return particle_count + compact_patch_indices.size();
  }
};

// One authoritative gas-incidence existence condition shared by
// RuntimeDecompositionSourceStorage construction and every MemoryGovernor
// admission that charges gas-incidence structures. Gas-empty/DMO states have
// no gas identity rows or no cells, so no offsets, indices, or construction
// scratch can exist for them.
[[nodiscard]] bool runtimeDecompositionHasGasIncidenceSource(
    const core::SimulationState& state) noexcept;

// Authoritative pre-construction byte model for the structures
// RuntimeDecompositionSourceStorage retains or temporarily requires while
// building a RuntimeDecompositionSourceView. Formulas follow the storage
// constructor's actual allocation semantics (capacity-based where the storage
// reports capacity): DMO/gas-empty states charge zero gas-incidence bytes
// through runtimeDecompositionHasGasIncidenceSource, the compact patch index
// is reserved once per patch row, and no removed legacy vector (such as a
// derived patch cell-count array) is charged. The total is a conservative
// envelope: every component is an upper bound on the corresponding live
// allocation, so estimate.total_bytes >= bytes the storage actually retains
// or transiently needs for those components.
struct RuntimeDecompositionSourceMemoryEstimate {
  // Zero when active_particle_indices is empty (no mask is allocated).
  std::uint64_t active_mask_bytes = 0U;
  std::uint64_t compact_patch_index_bytes = 0U;
  // Zero unless runtimeDecompositionHasGasIncidenceSource(state).
  std::uint64_t gas_incidence_offset_bytes = 0U;
  std::uint64_t gas_incidence_index_bytes = 0U;
  std::uint64_t gas_incidence_construction_bytes = 0U;
  // Fixed species memory table owned inline by the storage.
  std::uint64_t other_source_owned_bytes = 0U;
  std::uint64_t total_bytes = 0U;
};

// Shared MemoryGovernor admission estimate for both the runtime rebalance
// plan reservation and the top-domain seed reservation. Checked arithmetic;
// throws std::overflow_error on impossible populations.
[[nodiscard]] RuntimeDecompositionSourceMemoryEstimate
estimateRuntimeDecompositionSourceStorage(
    const core::SimulationState& state,
    std::span<const std::uint32_t> active_particle_indices);

// Authoritative MemoryGovernor admission estimate for the compact startup
// gravity-aware planner's transient live set. This answers the planner
// question ("can the planner safely allocate its temporary working set
// now?") and is additional to, never a replacement for, the resulting-cut
// feasibility checks (max_rank_memory_bytes / rank_transient_reserve_bytes).
// The total models the simultaneous live peak: compact records and the
// startup occupancy/patch-mapping scratch coexist with the planner's
// MemoryGroup population, owner vector, sorted-index vector, plan metric and
// range lanes, and the bounded per-entity gas-incidence lookup temporaries.
// MemoryGroup is modeled as {size_t begin, end; uint64_t memory_bytes;
// double weighted_load}; a static_assert inside the planner TU enforces that
// layout against the actual cut-algorithm struct. entity_upper_bound is
// particles + patches (records.reserve upper bound, which also bounds the
// one-group-per-record worst case: groups.reserve(entity_upper_bound)).
// density_grid_cell_count is the startup density grid cell count and
// pm_x_bins the startup PM-x occupancy lane count, both supplied by the
// startup caller so the grid policy stays owned by the startup workflow.
// Checked arithmetic; throws std::overflow_error on impossible populations.
struct CompactStartupPlannerMemoryEstimate {
  std::uint64_t compact_record_bytes = 0U;
  std::uint64_t memory_group_bytes = 0U;
  std::uint64_t owning_rank_bytes = 0U;
  std::uint64_t sorted_index_bytes = 0U;
  std::uint64_t occupancy_bytes = 0U;
  std::uint64_t patch_mapping_bytes = 0U;
  std::uint64_t other_known_scratch_bytes = 0U;
  std::uint64_t total_bytes = 0U;
};

[[nodiscard]] CompactStartupPlannerMemoryEstimate
estimateCompactStartupPlannerTransientBytes(
    std::size_t entity_upper_bound,
    std::size_t patch_count,
    std::size_t cell_count,
    std::size_t world_size,
    std::size_t pm_x_bins,
    std::size_t density_grid_cell_count);

struct CompactTopDomainSeedRecord {
  std::uint64_t sfc_key = 0;
  std::uint64_t entity_id = 0;
  std::size_t local_index = 0;
  DecompositionEntityKind kind = DecompositionEntityKind::kParticle;
};
static_assert(sizeof(CompactTopDomainSeedRecord) <= 32U,
              "CompactTopDomainSeedRecord must fit in 32 bytes");

[[nodiscard]] std::uint64_t sfcKeyForCompactRuntimeRecord(
    const CompactRuntimeDecompositionRecord& record) noexcept;

// Canonical 10-bit-quantized Morton key for one coordinate triple under the
// configured decomposition domain. This is the single SFC key rule shared by
// the rich item planner, compact records, source-view records, and startup
// initial placement.
[[nodiscard]] std::uint64_t sfcKeyForPosition(
    double x_comov,
    double y_comov,
    double z_comov,
    const DecompositionConfig& config);

// Streaming construction: one pass over items, no second geometry/component
// population. weighted_load preserves the exact weightedLoad() contract.
[[nodiscard]] std::vector<CompactRuntimeDecompositionRecord>
makeCompactRuntimeDecompositionRecords(
    std::span<const DecompositionItem> items,
    const DecompositionConfig& config);

struct DecompositionConfig {
  int world_size = 1;
  double domain_x_min_comov = 0.0;
  double domain_x_max_comov = 1.0;
  double domain_y_min_comov = 0.0;
  double domain_y_max_comov = 1.0;
  double domain_z_min_comov = 0.0;
  double domain_z_max_comov = 1.0;
  double owned_particle_weight = 1.0;
  double active_target_weight = 0.0;
  double remote_tree_interaction_weight = 0.0;
  double work_weight = 0.0;
  double memory_weight = 0.0;
  DecompositionWeightCoefficients component_weights{};
  // Hard rank-local envelope for decomposition-owned persistent state plus the
  // declared transient reserve. Zero disables this guard for compatibility.
  std::uint64_t max_rank_memory_bytes = 0;
  std::uint64_t rank_transient_reserve_bytes = 0;
  bool prefer_component_work_model = true;
};

// Single weighted-load authority for entities that publish explicit work
// components (the startup gravity-aware path and any caller holding resolved
// components without a full DecompositionItem). Applies the same
// component-vs-legacy selection, fallback chain, and final default as the
// rich item contract; components are clamped to nonnegative values first.
[[nodiscard]] double weightedLoadFromExplicitComponents(
    const DecompositionWorkComponents& components,
    DecompositionEntityKind kind,
    std::uint64_t active_target_count_recent,
    std::uint64_t remote_tree_interactions_recent,
    double work_units,
    std::uint64_t memory_bytes,
    const DecompositionConfig& config);

struct RankRange {
  std::size_t begin_sorted = 0;
  std::size_t end_sorted = 0;
};

struct LoadBalanceMetrics {
  std::vector<double> weighted_load_by_rank;
  std::vector<std::uint64_t> memory_bytes_by_rank;
  std::vector<std::uint64_t> owned_particles_by_rank;
  std::vector<std::uint64_t> active_targets_by_rank;
  std::vector<std::uint64_t> remote_tree_interactions_by_rank;
  std::vector<double> particle_count_cost_by_rank;
  std::vector<double> gas_cell_cost_by_rank;
  std::vector<double> tree_interaction_cost_by_rank;
  std::vector<double> pm_mesh_cost_by_rank;
  std::vector<double> amr_patch_cost_by_rank;
  std::vector<double> active_fraction_cost_by_rank;
  std::vector<double> memory_pressure_cost_by_rank;
  std::vector<double> transient_memory_cost_by_rank;
  std::vector<double> source_event_cost_by_rank;
  std::vector<double> communication_cost_by_rank;
  std::vector<double> gpu_occupancy_cost_by_rank;
  std::vector<double> generic_work_cost_by_rank;
  double mean_weighted_load = 0.0;
  double max_weighted_load = 0.0;
  double weighted_imbalance_ratio = 0.0;
  std::uint64_t total_memory_bytes = 0;
  std::uint64_t max_memory_bytes = 0;
  double memory_imbalance_ratio = 0.0;
  std::vector<std::uint64_t> peak_memory_bytes_by_rank;
  std::uint64_t max_peak_memory_bytes = 0;
  double peak_memory_imbalance_ratio = 0.0;
};

// Accumulates one entity's clamped work components into the rank's component
// metric lanes with sign +1 (accumulate) or -1 (withdraw). Shared by the
// compact planner, the rich planner adapter, current-ownership metrics, and
// the startup component pass.
void addWorkComponentsToMetrics(
    LoadBalanceMetrics& metrics,
    std::size_t rank,
    const DecompositionWorkComponents& components,
    double sign);

struct DecompositionPlan {
  std::vector<int> owning_rank_by_item;
  std::vector<std::size_t> sorted_indices;
  std::vector<RankRange> ranges_by_rank;
  LoadBalanceMetrics metrics;
};

// Compact authoritative top-domain leaf exported by the decomposition layer
// for locality-aware consumers such as TreePM. A rank may own multiple leaves.
// Bounds describe authoritative decomposition units rather than a gravity-tree
// root, so topology validity is tied to the decomposition epoch/geometry, not
// the frequency with which a local gravity tree is rebuilt.
struct TopDomainLeaf {
  std::uint64_t domain_leaf_id = 0;
  int owner_rank = -1;
  std::uint64_t decomposition_epoch = 0;
  std::uint64_t sfc_key_begin = 0;
  std::uint64_t sfc_key_end = 0;
  double min_x_comov = 0.0;
  double max_x_comov = 0.0;
  double min_y_comov = 0.0;
  double max_y_comov = 0.0;
  double min_z_comov = 0.0;
  double max_z_comov = 0.0;
  double work_weight = 0.0;
  std::uint64_t entity_count = 0;
  bool periodic_geometry = false;
};

[[nodiscard]] std::vector<TopDomainLeaf> buildAuthoritativeTopDomainLeaves(
    std::span<const DecompositionItem> local_items,
    const DecompositionConfig& config,
    int owner_rank,
    std::uint64_t decomposition_epoch,
    std::size_t max_leaves_per_rank = 8U);

// Compact seed builder for the canonical production source view: one
// minimal-record pass, one in-place total-order sort, bounded leaf grouping,
// and streaming kind + local_index geometry resolution.
[[nodiscard]] std::vector<TopDomainLeaf> buildAuthoritativeTopDomainLeavesFromCompact(
    std::span<const DecompositionItem> local_items,
    const DecompositionConfig& config,
    int owner_rank,
    std::uint64_t decomposition_epoch,
    std::size_t max_leaves_per_rank = 8U);

[[nodiscard]] std::vector<TopDomainLeaf> buildAuthoritativeTopDomainLeavesFromSource(
    const RuntimeDecompositionSourceView& source,
    const DecompositionConfig& config,
    int owner_rank,
    std::uint64_t decomposition_epoch,
    std::size_t max_leaves_per_rank = 8U);

// Observability for a compact leaf refit against the current source snapshot.
struct TopDomainGeometryRefitDiagnostics {
  std::uint64_t refreshed_leaf_count = 0;
  std::uint64_t out_of_seed_range_source_count = 0;
  std::uint64_t empty_seed_leaf_count = 0;
  std::uint64_t source_count = 0;
};

// O(N_local) refresh of authoritative top-domain leaf bounds without
// reconstructing DecompositionItems. Seed leaves supply ownership, epoch, and
// SFC intervals; a single scan of the provided source coordinates expands
// each interval's AABB and assigns out-of-interval sources to the nearest
// seed leaf (still covered; ownership unchanged). Empty leaves are omitted
// from the published geometry so no non-finite bounds are exported.
[[nodiscard]] std::vector<TopDomainLeaf> refitAuthoritativeTopDomainLeaves(
    std::span<const TopDomainLeaf> seed_leaves,
    std::span<const double> pos_x_comoving,
    std::span<const double> pos_y_comoving,
    std::span<const double> pos_z_comoving,
    const DecompositionConfig& config,
    int owner_rank,
    std::uint64_t decomposition_epoch,
    TopDomainGeometryRefitDiagnostics* diagnostics = nullptr);

[[nodiscard]] std::uint64_t topDomainGeometryFingerprint(
    std::span<const TopDomainLeaf> leaves) noexcept;

[[nodiscard]] DecompositionPlan buildMortonSfcDecomposition(
    std::span<const DecompositionItem> items,
    const DecompositionConfig& config);

// Compact SFC decomposition planner over ≤64-byte records. One in-place
// deterministic total-order std::sort on (sfc_key, entity_id, local_index),
// indivisible (key, ID) memory-group cuts, hard rank-memory feasibility, and
// streaming metric recompute — the single cut algorithm shared by the rich
// reference adapter and startup initial placement. Precondition: record
// {local_index} must form a permutation of [0, records.size()): every index
// in range, every index exactly once. The planner enforces this contract
// (range, uniqueness, and full coverage) before ownership write-back, so
// owning_rank_by_item and sorted_indices are expressed in that local_index
// space with no default-owner-zero masking of missing indices.
// components_by_local_index, when non-empty, must cover [0, records.size()) and
// supplies per-rank work-component metrics during the sorted pass; startup
// callers that fold components into weighted_load at record creation may leave
// it empty and accumulate component diagnostics afterwards via
// addWorkComponentsToMetrics.
[[nodiscard]] DecompositionPlan buildMortonSfcDecompositionFromCompact(
    std::span<CompactRuntimeDecompositionRecord> records,
    const DecompositionConfig& config,
    std::span<const DecompositionWorkComponents> components_by_local_index = {});

struct DecompositionRuntimeMeasurements {
  // Measured feedback from the previous solver window. These are aggregate
  // counters/timings that are distributed across decomposition units according
  // to their existing spatial/proxy contribution. A zero/empty frame leaves
  // proxy-only costs unchanged.
  std::uint64_t tree_pair_evaluations_recent = 0;
  std::uint64_t tree_remote_request_bytes_recent = 0;
  std::uint64_t pm_mesh_cells_touched_recent = 0;
  std::uint64_t pm_fft_transpose_bytes_recent = 0;
  std::uint64_t amr_patch_cells_updated_recent = 0;
  std::uint64_t hydro_face_fluxes_recent = 0;
  std::uint64_t ghost_exchange_bytes_recent = 0;
  double tree_wall_ms_recent = 0.0;
  double pm_wall_ms_recent = 0.0;
  double amr_wall_ms_recent = 0.0;
  double hydro_wall_ms_recent = 0.0;
  double gpu_kernel_ms_recent = 0.0;
  double accelerator_occupancy_fraction_recent = 0.0;
  bool has_measurements = false;
};

struct DecompositionFeedbackCoefficients {
  double measured_tree_pair = 1.0;
  double measured_pm_cell = 1.0;
  double measured_amr_cell = 1.0;
  double measured_hydro_face = 1.0;
  double measured_wall_ms = 1.0;
};

void applyRuntimeDecompositionFeedback(
    std::span<DecompositionItem> items,
    const DecompositionRuntimeMeasurements& measurements,
    const DecompositionFeedbackCoefficients& coefficients);

struct RuntimeRebalanceConfig {
  int world_size = 1;
  double imbalance_trigger_ratio = 1.25;
  double memory_trigger_ratio = 1.50;
  double max_migrated_load_fraction = 0.25;
  bool allow_particle_migration = true;
  bool allow_amr_patch_reassignment = true;
};

struct ParticleMigrationIntent {
  std::uint64_t particle_id = 0;
  std::size_t item_index = 0;
  int old_owner_rank = 0;
  int new_owner_rank = 0;
  double work_units = 0.0;
};

struct AmrPatchOwnershipUpdate {
  std::uint64_t patch_id = 0;
  int old_owner_rank = 0;
  int new_owner_rank = 0;
};

class MpiContext;

struct RuntimeRebalancePlan {
  bool should_rebalance = false;
  std::string reason;
  LoadBalanceMetrics current_metrics{};
  DecompositionPlan target_decomposition{};
  std::vector<ParticleMigrationIntent> particle_migrations;
  std::vector<AmrPatchOwnershipUpdate> amr_patch_ownership_updates;
  double migrated_load = 0.0;
  double migrated_load_fraction = 0.0;
  std::vector<std::uint64_t> sfc_cut_keys;
  std::vector<std::uint64_t> sfc_cut_entity_ids;
  std::uint64_t local_entities_considered = 0;
  std::uint64_t global_entities_considered = 0;
  std::uint64_t local_entities_moved = 0;
  std::uint64_t global_entities_moved = 0;
  std::uint64_t local_bytes_moved = 0;
  std::uint64_t global_bytes_moved = 0;
  std::uint64_t local_control_bytes = 0;
  std::uint64_t global_control_bytes = 0;
  std::uint64_t peak_temporary_bytes = 0;
  double cut_displacement_fraction = 0.0;
  bool used_distributed_sfc_cuts = false;
  bool exact_debug_audit_enabled = false;
  // Compact-planner provenance (zero on the legacy reference path).
  bool used_compact_planner = false;
  std::uint64_t planner_record_bytes = 0;
  std::uint64_t planner_sample_bytes = 0;
  std::uint64_t planner_prefix_bytes = 0;
   std::uint64_t planner_migration_intent_bytes = 0;
   std::uint64_t planner_other_known_scratch_bytes = 0;
   std::uint64_t planner_local_entity_count = 0;
   std::uint64_t planner_peak_temporary_bytes = 0;
   std::uint64_t planner_local_peak_temporary_bytes = 0;
   double planner_bytes_per_entity = 0.0;

};

[[nodiscard]] LoadBalanceMetrics computeCurrentOwnershipLoadBalanceMetrics(
    std::span<const DecompositionItem> items,
    const DecompositionConfig& config);

[[nodiscard]] RuntimeRebalancePlan buildRuntimeRebalancePlan(
    std::span<const DecompositionItem> items,
    const DecompositionConfig& decomposition_config,
    const RuntimeRebalanceConfig& rebalance_config);

[[nodiscard]] RuntimeRebalancePlan buildDistributedRuntimeRebalancePlan(
    const MpiContext& mpi_context,
    std::span<const DecompositionItem> local_items,
    const DecompositionConfig& decomposition_config,
    const RuntimeRebalanceConfig& rebalance_config);

// Production large-N planner: compact SFC records, single in-place total-order
// sort (sfc_key, entity_id, local_index), bounded cut samples, migration
// intents only for owner changes, streaming metric recompute. Leaves
// target_decomposition.owning_rank_by_item / sorted_indices empty (not used
// by the production migration path). Legacy rich builders remain available
// only as reference/debug/test compatibility paths; they are not used by
// runtime rebalance, top-domain seed construction, or startup initial
// placement, all three of which are compact.
[[nodiscard]] RuntimeRebalancePlan buildCompactDistributedRuntimeRebalancePlan(
    const MpiContext& mpi_context,
    std::span<const DecompositionItem> local_items,
    const DecompositionConfig& decomposition_config,
    const RuntimeRebalanceConfig& rebalance_config);

[[nodiscard]] RuntimeRebalancePlan buildCompactDistributedRuntimeRebalancePlan(
    const MpiContext& mpi_context,
    const RuntimeDecompositionSourceView& source,
    const DecompositionConfig& decomposition_config,
    const RuntimeRebalanceConfig& rebalance_config,
    const DecompositionRuntimeMeasurements& measurements,
    const DecompositionFeedbackCoefficients& feedback_coefficients);

[[nodiscard]] std::vector<DecompositionItem> gatherDecompositionItemsAcrossRanks(
    const MpiContext& mpi_context,
    std::span<const DecompositionItem> local_items);

// Classic MPI collectives expose signed-int counts/displacements. Logical
// transfers are planned in size_t and split into bounded round-local layouts so
// no individual MPI call depends on a population-scale aggregate fitting in
// int. The public planner is intentionally small so overflow/boundary behavior
// can be unit-tested without allocating large payloads.
inline constexpr std::size_t k_default_mpi_transport_round_bytes =
    16U * 1024U * 1024U;

struct BoundedMpiRoundLayout {
  std::vector<int> counts;
  std::vector<int> displacements;
  std::vector<std::size_t> logical_offsets;
  std::size_t round_count = 0U;
};

struct BoundedMpiTransferPlan {
  std::vector<std::size_t> logical_counts;
  std::vector<std::size_t> logical_displacements;
  std::vector<BoundedMpiRoundLayout> rounds;
  std::size_t logical_total_count = 0U;
  std::size_t per_peer_count_limit = 0U;
};

[[nodiscard]] BoundedMpiTransferPlan planBoundedMpiTransferRounds(
    std::span<const std::size_t> logical_counts,
    std::size_t mpi_count_limit = static_cast<std::size_t>(std::numeric_limits<int>::max()),
    std::size_t round_count_limit = k_default_mpi_transport_round_bytes);

// Production uses the conservative internal round bound above. Test-enabled
// builds may override it through a test-only environment seam so MPI
// integration tests can force multiple rounds using tiny payloads.
[[nodiscard]] std::size_t mpiTransportRoundLimitBytes();

[[nodiscard]] std::vector<std::vector<std::uint8_t>> exchangeBoundedAlltoallBytes(
    const MpiContext& mpi_context,
    const std::vector<std::vector<std::uint8_t>>& send_payloads);


enum class LocalIndexResidency : std::uint8_t {
  kOwned = 0,
  kGhost = 1,
};

enum class ExchangeObjectKind : std::uint8_t {
  kLocalParticle = 0,
  kImportedGhostParticle = 1,
  kTreePseudoParticle = 2,
  kPmMeshCell = 3,
  kHydroGhostCell = 4,
  kAmrPatchMetadata = 5,
};

class MpiContext;
struct PmSlabLayout;

struct OwnershipDescriptor {
  ExchangeObjectKind kind = ExchangeObjectKind::kLocalParticle;
  std::uint64_t object_id = 0;
  int owner_rank = 0;
  int local_rank = 0;
  std::uint64_t decomposition_epoch = 0;
  bool is_authoritative = true;
  bool is_mutable = true;
};

void validateOwnershipDescriptor(const OwnershipDescriptor& descriptor);

struct GhostLayerEpoch {
  std::uint64_t decomposition_epoch = 0;
  std::uint64_t ghost_sync_epoch = 0;
  std::uint64_t particle_index_generation = 0;

  [[nodiscard]] bool matches(const GhostLayerEpoch& expected) const noexcept;
};

struct LocalGhostDescriptor {
  LocalIndexResidency residency = LocalIndexResidency::kOwned;
  int owning_rank = 0;
  std::uint64_t particle_id = 0;
  GhostLayerEpoch epoch{};
};

enum class GhostTransferRole : std::uint8_t {
  kOutboundSend = 0,
  kInboundReceive = 1,
};

enum class GhostTransferIntent : std::uint8_t {
  kGhostRefreshRequest = 0,
  kGhostRefreshReceiveStaging = 1,
  kOwnershipMigrationSend = 2,
  kOwnershipMigrationReceiveStaging = 3,
};

struct GhostTransferDescriptor {
  GhostTransferRole role = GhostTransferRole::kOutboundSend;
  GhostTransferIntent intent = GhostTransferIntent::kGhostRefreshRequest;
  int peer_rank = 0;
  std::size_t neighbor_slot = 0;
  LocalIndexResidency expected_post_transfer_residency = LocalIndexResidency::kGhost;
  std::vector<std::uint32_t> local_indices;
};

struct GhostExchangePlan {
  std::vector<int> neighbor_ranks;
  std::vector<std::vector<std::uint32_t>> send_local_indices_by_neighbor;
  std::vector<std::vector<std::uint32_t>> recv_local_indices_by_neighbor;
  std::vector<GhostTransferDescriptor> outbound_transfers;
  std::vector<GhostTransferDescriptor> inbound_transfers;
  std::uint64_t send_bytes = 0;
  std::uint64_t recv_bytes = 0;
  GhostLayerEpoch epoch{};
  // Monotonic phase discriminator for future nonblocking/overlap layers.
  // Blocking code derives it from GhostLayerEpoch::ghost_sync_epoch so rank-local
  // neighbor ordering can never influence message matching.
  std::uint64_t exchange_sequence = 0;
  bool uses_blocking_exchange = true;
  bool nonblocking_overlap_enabled = false;
};

[[nodiscard]] GhostExchangePlan buildGhostExchangePlan(
    int world_rank,
    std::span<const LocalGhostDescriptor> local_ghost_descriptors,
    std::size_t bytes_per_ghost);

[[nodiscard]] GhostExchangePlan buildGhostExchangePlan(
    int world_rank,
    std::span<const int> ghost_owner_rank_by_local_index,
    std::size_t bytes_per_ghost);

[[nodiscard]] GhostExchangePlan buildExplicitGhostExchangePlan(
    int world_rank,
    std::span<const int> neighbor_ranks,
    std::span<const std::vector<std::uint32_t>> send_local_indices_by_neighbor,
    std::span<const std::vector<std::uint32_t>> recv_local_indices_by_neighbor,
    std::size_t bytes_per_ghost,
    const GhostLayerEpoch& epoch,
    bool enable_nonblocking_overlap = false);

void validateGhostExchangePlan(const GhostExchangePlan& plan);
void validateGhostTransferAgainstResidency(
    const GhostTransferDescriptor& descriptor,
    std::span<const LocalGhostDescriptor> local_ghost_descriptors,
    int world_rank);
void validateBlockingGhostExchangeContracts(
    const GhostExchangePlan& plan,
    std::span<const LocalGhostDescriptor> local_ghost_descriptors,
    int world_rank,
    const GhostLayerEpoch& expected_epoch);

struct ReductionAgreement {
  double deterministic_baseline_sum = 0.0;
  double measured_sum = 0.0;
  double absolute_error = 0.0;
  double relative_error = 0.0;
};

enum class ReductionAgreementMode : std::uint8_t {
  kAbsoluteOnly = 0,
  kRelativeOnly = 1,
  kAbsoluteAndRelative = 2,
  kAbsoluteOrRelative = 3,
};

struct ReductionAgreementPolicy {
  ReductionAgreementMode mode = ReductionAgreementMode::kAbsoluteOrRelative;
  double absolute_tolerance = 0.0;
  double relative_tolerance = 0.0;
};


struct LocalOwnershipIdentitySummary {
  std::uint64_t local_owned_count = 0;
  std::uint64_t local_particle_id_sum = 0;
  std::uint64_t local_particle_id_square_sum = 0;
  std::uint64_t local_particle_id_xor = 0;
  bool local_particle_ids_unique = true;
};

struct ExactOwnershipPartitionReport {
  std::uint64_t global_owned_count = 0;
  bool local_particle_ids_unique = true;
  bool globally_unique = true;
  bool matches_expected_ids = true;
  std::vector<std::uint64_t> duplicate_particle_ids;
  std::vector<std::uint64_t> missing_expected_particle_ids;
  std::vector<std::uint64_t> extra_particle_ids;

  [[nodiscard]] bool valid() const noexcept {
    return local_particle_ids_unique && globally_unique && matches_expected_ids &&
        duplicate_particle_ids.empty() && missing_expected_particle_ids.empty() && extra_particle_ids.empty();
  }
};

[[nodiscard]] LocalOwnershipIdentitySummary summarizeLocalOwnedParticleIds(
    std::span<const std::uint64_t> local_particle_ids);

[[nodiscard]] ExactOwnershipPartitionReport validateExactGlobalOwnershipPartition(
    const MpiContext& mpi_context,
    std::span<const std::uint64_t> local_owned_particle_ids,
    std::span<const std::uint64_t> expected_local_reference_particle_ids);

[[nodiscard]] bool partitionIdentityMatchesGeneratedSet(
    const LocalOwnershipIdentitySummary& reduced_global_summary,
    std::uint64_t expected_global_count,
    std::uint64_t expected_particle_id_sum,
    std::uint64_t expected_particle_id_square_sum,
    std::uint64_t expected_particle_id_xor);

[[nodiscard]] bool partitionIdentityMatchesGeneratedSet(
    const LocalOwnershipIdentitySummary& reduced_global_summary,
    std::uint64_t expected_global_count,
    std::uint64_t expected_particle_id_sum,
    std::uint64_t expected_particle_id_xor);

[[nodiscard]] double deterministicRankOrderedSum(std::span<const double> per_rank_values);
[[nodiscard]] ReductionAgreement compareReductionAgreement(
    std::span<const double> per_rank_values,
    double measured_sum);
[[nodiscard]] bool satisfiesReductionAgreement(
    const ReductionAgreement& agreement,
    const ReductionAgreementPolicy& policy);

struct RankConfigDigest {
  int world_rank = 0;
  std::uint64_t normalized_config_hash = 0;
  int mpi_ranks_expected = 1;
  bool deterministic_reduction = true;
};

enum class RankConfigMismatchProperty : std::uint8_t {
  kNormalizedConfigHash = 0,
  kMpiRanksExpected = 1,
  kDeterministicReduction = 2,
};

struct RankConfigMismatch {
  RankConfigMismatchProperty property = RankConfigMismatchProperty::kNormalizedConfigHash;
  int baseline_rank = 0;
  int rank = 0;
  std::string baseline_value;
  std::string rank_value;
};

struct RankConfigConsensus {
  bool normalized_config_hash_match = true;
  bool mpi_ranks_expected_match = true;
  bool deterministic_reduction_match = true;
  std::vector<int> mismatched_ranks;
  std::vector<RankConfigMismatch> mismatches;

  [[nodiscard]] bool allConsistent() const noexcept;
};

[[nodiscard]] RankConfigConsensus evaluateRankConfigConsensus(
    std::span<const RankConfigDigest> digests);

struct GhostExchangeBufferSoA {
  GhostLayerEpoch epoch{};
  std::vector<std::uint64_t> entity_id;
  // Gravity/PM/tree-facing lanes. Empty optional lanes are encoded as zero only
  // by legacy callers; production solver refreshes should populate them.
  std::vector<double> position_x_comoving;
  std::vector<double> position_y_comoving;
  std::vector<double> position_z_comoving;
  std::vector<double> mass_code;
  // Hydro-facing boundary state lanes.
  std::vector<double> density_code;
  std::vector<double> velocity_x_code;
  std::vector<double> velocity_y_code;
  std::vector<double> velocity_z_code;
  std::vector<double> pressure_code;
  std::vector<double> internal_energy_code;

  [[nodiscard]] bool isConsistent() const noexcept;
  [[nodiscard]] std::size_t size() const noexcept;
  [[nodiscard]] bool hasGravityPayload() const noexcept;
  [[nodiscard]] bool hasHydroPayload() const noexcept;
};

struct ReadOnlyGhostExchangeView {
  GhostLayerEpoch epoch{};
  std::span<const std::uint64_t> entity_id;
  std::span<const double> position_x_comoving;
  std::span<const double> position_y_comoving;
  std::span<const double> position_z_comoving;
  std::span<const double> mass_code;
  std::span<const double> density_code;
  std::span<const double> velocity_x_code;
  std::span<const double> velocity_y_code;
  std::span<const double> velocity_z_code;
  std::span<const double> pressure_code;
  std::span<const double> internal_energy_code;

  [[nodiscard]] std::size_t size() const noexcept;
  [[nodiscard]] bool isConsistent() const noexcept;
  [[nodiscard]] bool isFresh(const GhostLayerEpoch& expected_epoch) const noexcept;
};

[[nodiscard]] ReadOnlyGhostExchangeView makeReadOnlyGhostExchangeView(
    const GhostExchangeBufferSoA& storage);
void requireFreshGhostExchangeView(
    const ReadOnlyGhostExchangeView& view,
    const GhostLayerEpoch& expected_epoch);

class GhostExchangeBuffer {
 public:
  void clear();
  [[nodiscard]] std::size_t byteSize() const noexcept;

  void packFrom(
      const GhostExchangeBufferSoA& source,
      std::span<const std::uint32_t> local_indices);

  void packFrom(
      const GhostTransferDescriptor& descriptor,
      const GhostExchangeBufferSoA& source,
      std::span<const std::uint32_t> local_indices);

  void unpackAppendTo(GhostExchangeBufferSoA& destination) const;

  void unpackAppendTo(
      const GhostTransferDescriptor& descriptor,
      GhostExchangeBufferSoA& destination) const;

  [[nodiscard]] std::span<const std::uint8_t> encodedBytes() const noexcept;
  void replaceEncodedBytes(std::vector<std::uint8_t> bytes);

 private:
  std::vector<std::uint8_t> m_bytes;
};

[[nodiscard]] std::size_t ghostRefreshPayloadRecordBytes() noexcept;
void validateGhostRefreshPayloadDescriptor(const GhostTransferDescriptor& descriptor);
[[nodiscard]] int ghostExchangePairStableTag(int tag_base, int local_rank, int peer_rank);
[[nodiscard]] int ghostExchangeSequencedTag(
    int tag_base,
    int local_rank,
    int peer_rank,
    std::uint64_t exchange_sequence);

struct DistributedRestartState {
  std::uint32_t schema_version = 2;
  std::uint64_t decomposition_epoch = 0;
  int world_size = 1;
  std::size_t pm_grid_nx = 0;
  std::size_t pm_grid_ny = 0;
  std::size_t pm_grid_nz = 0;
  std::string pm_decomposition_mode = "slab";
  std::uint64_t gravity_kick_opportunity = 0;
  std::uint64_t pm_update_cadence_steps = 1;
  std::uint64_t long_range_field_version = 0;
  std::uint64_t last_long_range_refresh_opportunity = 0;
  std::uint64_t long_range_field_built_step_index = 0;
  double long_range_field_built_scale_factor = 1.0;
  std::string long_range_restart_policy = "deterministic_rebuild";
  std::vector<int> owning_rank_by_item;
  std::vector<std::size_t> pm_slab_begin_x_by_rank;
  std::vector<std::size_t> pm_slab_end_x_by_rank;

  // The stream path preserves the historical text representation without a
  // second population-sized string. serialize() performs a checked counting
  // pass and one exact-size destination allocation.
  void serializeTo(std::ostream& stream) const;
  [[nodiscard]] std::size_t serializedSizeBytes() const;
  [[nodiscard]] std::string serialize() const;
  [[nodiscard]] static DistributedRestartState deserialize(const std::string& encoded);
};

struct DistributedRestartCompatibilityReport {
  bool supported_schema_match = true;
  bool world_size_match = true;
  bool pm_grid_shape_match = true;
  bool pm_decomposition_mode_match = true;
  bool pm_slab_table_shape_match = true;
  bool pm_local_slab_match = true;
  bool pm_cadence_steps_match = true;
  bool gravity_kick_state_match = true;
  bool long_range_field_state_match = true;
  std::vector<std::string> mismatch_messages;

  [[nodiscard]] bool compatible() const noexcept {
    return supported_schema_match && world_size_match && pm_grid_shape_match && pm_decomposition_mode_match &&
        pm_slab_table_shape_match && pm_local_slab_match && pm_cadence_steps_match && gravity_kick_state_match &&
        long_range_field_state_match;
  }
};

struct DistributedExecutionTopology;

[[nodiscard]] DistributedRestartCompatibilityReport evaluateDistributedRestartCompatibility(
    const DistributedRestartState& restart_state,
    const DistributedExecutionTopology& runtime_topology);

class MpiContext {
 public:
  MpiContext();
  MpiContext(bool is_enabled, int world_size, int world_rank);

  [[nodiscard]] bool isEnabled() const noexcept;
  [[nodiscard]] bool isRoot() const noexcept;
  [[nodiscard]] int worldSize() const noexcept;
  [[nodiscard]] int worldRank() const noexcept;
  [[nodiscard]] int localRank() const noexcept;
  [[nodiscard]] int localSize() const noexcept;
  void validateExpectedWorldSizeOrThrow(int expected_world_size) const;

  [[nodiscard]] double allreduceSumDouble(double local_value) const;
  void allreduceSumDoublesInPlace(std::span<double> values) const;
  [[nodiscard]] double allreduceMinDouble(double local_value) const;
  [[nodiscard]] std::uint64_t allreduceSumUint64(std::uint64_t local_value) const;
  void allreduceSumUint64sInPlace(std::span<std::uint64_t> values) const;
  [[nodiscard]] std::uint64_t exclusiveScanSumUint64(std::uint64_t local_value) const;
  [[nodiscard]] std::uint64_t allreduceMaxUint64(std::uint64_t local_value) const;
  [[nodiscard]] std::uint64_t allreduceMinUint64(std::uint64_t local_value) const;
  [[nodiscard]] std::uint64_t allreduceXorUint64(std::uint64_t local_value) const;

  // Collective preflight gate used after any rank-local preparation that can
  // throw. All ranks must call this for the same phase. If any participant
  // failed, every participant throws before entering the following payload
  // communication; the first failing rank and its bounded diagnostic are
  // propagated deterministically.
  void rethrowCollectivePreparationFailure(
      const std::exception_ptr& local_failure,
      std::string_view phase_name) const;

  // Variable-size byte collectives used by correctness-first distributed
  // metadata/catalog assembly. Only the root receives gathered payload bytes;
  // broadcast returns the root payload on every rank.
  [[nodiscard]] std::vector<std::uint8_t> gatherBytesToRoot(
      std::span<const std::uint8_t> local_bytes, int root_rank = 0) const;
  [[nodiscard]] std::vector<std::uint8_t> broadcastBytesFromRoot(
      std::span<const std::uint8_t> root_bytes, int root_rank = 0) const;

  // Bounded byte all-gather used by distributed metadata paths. Rank order and
  // byte order within each rank are preserved exactly.
  [[nodiscard]] std::vector<std::uint8_t> allgatherBytesBounded(
      std::span<const std::uint8_t> local_bytes) const;

 private:
  bool m_is_enabled = false;
  int m_world_size = 1;
  int m_world_rank = 0;
  int m_local_rank = 0;
  int m_local_size = 1;
};

struct BlockingGhostExchangeResult {
  GhostExchangeBufferSoA received_ghosts;
  std::uint64_t sent_bytes = 0;
  std::uint64_t received_bytes = 0;
};

struct BlockingGhostRefreshExchange {
  GhostExchangePlan plan;
  BlockingGhostExchangeResult result;
};

struct GhostCacheLifecycle {
  GhostLayerEpoch epoch{};
  bool valid = false;
  std::uint64_t refresh_count = 0;
  std::uint64_t invalidation_count = 0;
};

void invalidateGhostCache(
    GhostCacheLifecycle& lifecycle,
    const GhostLayerEpoch& next_epoch);
void markGhostCacheCommitted(
    GhostCacheLifecycle& lifecycle,
    const GhostLayerEpoch& committed_epoch);
void requireValidGhostCache(
    const GhostCacheLifecycle& lifecycle,
    const GhostLayerEpoch& expected_epoch,
    std::string_view caller);

struct GhostRefreshCommitReport {
  std::size_t updated_ghost_slots = 0;
  std::uint64_t committed_payload_bytes = 0;
};

[[nodiscard]] GhostRefreshCommitReport commitBlockingGhostRefreshResult(
    GhostExchangeBufferSoA& ghost_storage,
    std::span<const LocalGhostDescriptor> local_ghost_descriptors,
    const GhostExchangePlan& plan,
    const BlockingGhostExchangeResult& result,
    const GhostLayerEpoch& expected_epoch);

[[nodiscard]] BlockingGhostExchangeResult executeBlockingGhostRefreshExchange(
    const MpiContext& mpi_context,
    const GhostExchangePlan& plan,
    std::span<const LocalGhostDescriptor> local_ghost_descriptors,
    const GhostExchangeBufferSoA& authoritative_local_state,
    const GhostLayerEpoch& expected_epoch);

// Correctness-first high-level ghost refresh. Local ghost descriptors declare
// which remote particle IDs this rank needs. The blocking path first exchanges
// those demands, derives outbound send rows by matching peer requests against
// authoritative owned local descriptors, then uses the same payload validation
// and commit contract as explicit plans.
[[nodiscard]] BlockingGhostRefreshExchange executeBlockingGhostRefreshExchangeFromDescriptors(
    const MpiContext& mpi_context,
    std::span<const LocalGhostDescriptor> local_ghost_descriptors,
    const GhostExchangeBufferSoA& authoritative_local_state,
    const GhostLayerEpoch& expected_epoch);

struct RankDeviceAssignment {
  int requested_device_count = 0;
  int visible_device_count = 0;
  int active_device_count = 0;
  int assigned_device_index = -1;
  bool uses_cuda = false;

  [[nodiscard]] bool isValid() const noexcept;
};

struct DistributedExecutionTopology {
  int world_size = 1;
  int world_rank = 0;
  int local_rank = 0;
  bool mpi_enabled = false;
  std::string pm_decomposition_mode = "slab";
  PmSlabLayout pm_slab{};
  RankDeviceAssignment device_assignment{};

  [[nodiscard]] bool isDistributed() const noexcept { return world_size > 1; }
  [[nodiscard]] bool usesCuda() const noexcept { return device_assignment.uses_cuda; }
};

[[nodiscard]] RankDeviceAssignment selectRankDeviceAssignment(
    int local_rank,
    int configured_gpu_devices,
    bool cuda_runtime_available,
    int visible_device_count);

[[nodiscard]] DistributedExecutionTopology buildDistributedExecutionTopology(
    std::size_t global_nx,
    std::size_t global_ny,
    std::size_t global_nz,
    const MpiContext& mpi_context,
    int mpi_ranks_expected,
    int configured_gpu_devices,
    bool cuda_runtime_available,
    int visible_device_count,
    std::string pm_decomposition_mode = "slab");

}  // namespace cosmosim::parallel
