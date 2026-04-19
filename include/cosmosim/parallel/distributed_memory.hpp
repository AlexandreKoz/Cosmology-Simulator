#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <span>
#include <string>
#include <vector>

#include "cosmosim/core/profiling.hpp"

namespace cosmosim::parallel {

enum class DecompositionEntityKind : std::uint8_t {
  kParticle = 0,
  kAmrPatch = 1,
};

struct DecompositionItem {
  std::uint64_t entity_id = 0;
  DecompositionEntityKind kind = DecompositionEntityKind::kParticle;
  double x_comov = 0.0;
  double y_comov = 0.0;
  double z_comov = 0.0;
  std::uint64_t active_target_count_recent = 0;
  std::uint64_t remote_tree_interactions_recent = 0;
  double work_units = 1.0;
  std::uint64_t memory_bytes = 0;
};

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
};

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
  double mean_weighted_load = 0.0;
  double max_weighted_load = 0.0;
  double weighted_imbalance_ratio = 0.0;
  std::uint64_t total_memory_bytes = 0;
  std::uint64_t max_memory_bytes = 0;
  double memory_imbalance_ratio = 0.0;
};

struct DecompositionPlan {
  std::vector<int> owning_rank_by_item;
  std::vector<std::size_t> sorted_indices;
  std::vector<RankRange> ranges_by_rank;
  LoadBalanceMetrics metrics;
};

[[nodiscard]] DecompositionPlan buildMortonSfcDecomposition(
    std::span<const DecompositionItem> items,
    const DecompositionConfig& config);

enum class LocalIndexResidency : std::uint8_t {
  kOwned = 0,
  kGhost = 1,
};

struct LocalGhostDescriptor {
  LocalIndexResidency residency = LocalIndexResidency::kOwned;
  int owning_rank = 0;
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
};

[[nodiscard]] GhostExchangePlan buildGhostExchangePlan(
    int world_rank,
    std::span<const LocalGhostDescriptor> local_ghost_descriptors,
    std::size_t bytes_per_ghost);

[[nodiscard]] GhostExchangePlan buildGhostExchangePlan(
    int world_rank,
    std::span<const int> ghost_owner_rank_by_local_index,
    std::size_t bytes_per_ghost);

void validateGhostExchangePlan(const GhostExchangePlan& plan);

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
  std::vector<std::uint64_t> entity_id;
  std::vector<double> density_code;
  std::vector<double> velocity_x_code;
  std::vector<double> pressure_code;

  [[nodiscard]] bool isConsistent() const noexcept;
  [[nodiscard]] std::size_t size() const noexcept;
};

class GhostExchangeBuffer {
 public:
  void clear();
  [[nodiscard]] std::size_t byteSize() const noexcept;

  void packFrom(
      const GhostExchangeBufferSoA& source,
      std::span<const std::uint32_t> local_indices);

  void unpackAppendTo(GhostExchangeBufferSoA& destination) const;

 private:
  std::vector<std::uint8_t> m_bytes;
};

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

  [[nodiscard]] std::string serialize() const;
  [[nodiscard]] static DistributedRestartState deserialize(const std::string& encoded);
};

struct DistributedRestartCompatibilityReport {
  bool world_size_match = true;
  bool pm_grid_shape_match = true;
  bool pm_decomposition_mode_match = true;
  bool pm_local_slab_match = true;
  std::vector<std::string> mismatch_messages;

  [[nodiscard]] bool compatible() const noexcept {
    return world_size_match && pm_grid_shape_match && pm_decomposition_mode_match && pm_local_slab_match;
  }
};

struct DistributedExecutionTopology;

[[nodiscard]] DistributedRestartCompatibilityReport evaluateDistributedRestartCompatibility(
    const DistributedRestartState& restart_state,
    const DistributedExecutionTopology& runtime_topology);

struct PmSlabRange {
  std::size_t begin_x = 0;
  std::size_t end_x = 0;

  [[nodiscard]] std::size_t extentX() const noexcept {
    return (end_x >= begin_x) ? (end_x - begin_x) : 0;
  }
  [[nodiscard]] bool contains(std::size_t global_x) const noexcept {
    return global_x >= begin_x && global_x < end_x;
  }
};

struct PmSlabLayout {
  std::size_t global_nx = 0;
  std::size_t global_ny = 0;
  std::size_t global_nz = 0;
  int world_size = 1;
  int world_rank = 0;
  PmSlabRange owned_x{};

  [[nodiscard]] std::size_t local_nx() const noexcept {
    return owned_x.extentX();
  }
  [[nodiscard]] std::size_t localCellCount() const noexcept {
    return local_nx() * global_ny * global_nz;
  }
  [[nodiscard]] bool isValid() const noexcept;
  [[nodiscard]] bool ownsGlobalX(std::size_t global_x) const noexcept {
    return owned_x.contains(global_x);
  }
  [[nodiscard]] bool ownsGlobalCell(std::size_t global_x, std::size_t global_y, std::size_t global_z) const noexcept {
    if (global_y >= global_ny || global_z >= global_nz) {
      return false;
    }
    return ownsGlobalX(global_x);
  }
  [[nodiscard]] std::size_t localXFromGlobal(std::size_t global_x) const;
  [[nodiscard]] std::size_t globalXFromLocal(std::size_t local_x) const;
  [[nodiscard]] std::size_t localLinearIndex(std::size_t global_x, std::size_t global_y, std::size_t global_z) const;
  [[nodiscard]] bool ownsFullDomain() const noexcept {
    return owned_x.begin_x == 0 && owned_x.end_x == global_nx;
  }
};

[[nodiscard]] inline PmSlabRange pmOwnedXRangeForRank(std::size_t global_nx, int world_size, int rank) {
  if (global_nx == 0) {
    throw std::invalid_argument("global_nx must be positive");
  }
  if (world_size <= 0) {
    throw std::invalid_argument("world_size must be positive");
  }
  if (rank < 0 || rank >= world_size) {
    throw std::invalid_argument("rank must be within [0, world_size)");
  }
  const std::size_t world = static_cast<std::size_t>(world_size);
  const std::size_t rank_u = static_cast<std::size_t>(rank);
  const std::size_t base = global_nx / world;
  const std::size_t remainder = global_nx % world;
  const std::size_t begin = rank_u * base + std::min(rank_u, remainder);
  const std::size_t extent = base + (rank_u < remainder ? 1U : 0U);
  return {.begin_x = begin, .end_x = begin + extent};
}

[[nodiscard]] inline int pmOwnerRankForGlobalX(std::size_t global_nx, int world_size, std::size_t global_x) {
  if (global_x >= global_nx) {
    throw std::out_of_range("global_x must be within [0, global_nx)");
  }
  for (int rank = 0; rank < world_size; ++rank) {
    const PmSlabRange owned = pmOwnedXRangeForRank(global_nx, world_size, rank);
    if (owned.contains(global_x)) {
      return rank;
    }
  }
  throw std::logic_error("no owner rank found for global_x");
}

[[nodiscard]] inline int pmOwnerRankForGlobalCell(
    std::size_t global_nx,
    std::size_t global_ny,
    std::size_t global_nz,
    int world_size,
    std::size_t global_x,
    std::size_t global_y,
    std::size_t global_z) {
  if (global_y >= global_ny || global_z >= global_nz) {
    throw std::out_of_range("global cell y/z index out of range");
  }
  return pmOwnerRankForGlobalX(global_nx, world_size, global_x);
}
[[nodiscard]] inline PmSlabLayout makePmSlabLayout(
    std::size_t global_nx,
    std::size_t global_ny,
    std::size_t global_nz,
    int world_size,
    int world_rank) {
  PmSlabLayout layout{
      .global_nx = global_nx,
      .global_ny = global_ny,
      .global_nz = global_nz,
      .world_size = world_size,
      .world_rank = world_rank,
      .owned_x = pmOwnedXRangeForRank(global_nx, world_size, world_rank),
  };
  if (!layout.isValid()) {
    throw std::invalid_argument("constructed PM slab layout is invalid");
  }
  return layout;
}

inline bool PmSlabLayout::isValid() const noexcept {
  if (global_nx == 0 || global_ny == 0 || global_nz == 0) {
    return false;
  }
  if (world_size <= 0 || world_rank < 0 || world_rank >= world_size) {
    return false;
  }
  if (owned_x.begin_x > owned_x.end_x || owned_x.end_x > global_nx) {
    return false;
  }
  const PmSlabRange expected = pmOwnedXRangeForRank(global_nx, world_size, world_rank);
  return expected.begin_x == owned_x.begin_x && expected.end_x == owned_x.end_x;
}

inline std::size_t PmSlabLayout::localXFromGlobal(std::size_t global_x) const {
  if (!ownsGlobalX(global_x)) {
    throw std::out_of_range("global x index is not owned by this PM slab");
  }
  return global_x - owned_x.begin_x;
}

inline std::size_t PmSlabLayout::globalXFromLocal(std::size_t local_x) const {
  if (local_x >= local_nx()) {
    throw std::out_of_range("local PM slab x index out of range");
  }
  return owned_x.begin_x + local_x;
}

inline std::size_t PmSlabLayout::localLinearIndex(std::size_t global_x, std::size_t global_y, std::size_t global_z) const {
  if (!ownsGlobalCell(global_x, global_y, global_z)) {
    throw std::out_of_range("global PM cell is not owned by this slab");
  }
  const std::size_t local_x = localXFromGlobal(global_x);
  return (local_x * global_ny + global_y) * global_nz + global_z;
}

void recordDistributedProfiling(
    core::ProfilerSession* profiler,
    const LoadBalanceMetrics& metrics,
    std::uint64_t ghost_exchange_send_bytes,
    std::uint64_t ghost_exchange_recv_bytes);

class MpiContext {
 public:
  MpiContext();
  MpiContext(bool is_enabled, int world_size, int world_rank);

  [[nodiscard]] bool isEnabled() const noexcept;
  [[nodiscard]] bool isRoot() const noexcept;
  [[nodiscard]] int worldSize() const noexcept;
  [[nodiscard]] int worldRank() const noexcept;
  void validateExpectedWorldSizeOrThrow(int expected_world_size) const;

  [[nodiscard]] double allreduceSumDouble(double local_value) const;
  [[nodiscard]] std::uint64_t allreduceSumUint64(std::uint64_t local_value) const;

 private:
  bool m_is_enabled = false;
  int m_world_size = 1;
  int m_world_rank = 0;
};

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
  bool mpi_enabled = false;
  PmSlabLayout pm_slab{};
  RankDeviceAssignment device_assignment{};

  [[nodiscard]] bool isDistributed() const noexcept { return world_size > 1; }
  [[nodiscard]] bool usesCuda() const noexcept { return device_assignment.uses_cuda; }
};

[[nodiscard]] RankDeviceAssignment selectRankDeviceAssignment(
    int world_rank,
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
    int visible_device_count);

}  // namespace cosmosim::parallel
