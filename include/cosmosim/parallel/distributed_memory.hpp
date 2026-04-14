#pragma once

#include <cstddef>
#include <cstdint>
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
  double work_weight = 1.0;
  double memory_weight = 0.0;
};

struct RankRange {
  std::size_t begin_sorted = 0;
  std::size_t end_sorted = 0;
};

struct LoadBalanceMetrics {
  std::vector<double> weighted_load_by_rank;
  std::vector<std::uint64_t> memory_bytes_by_rank;
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

struct GhostExchangePlan {
  std::vector<int> neighbor_ranks;
  std::vector<std::vector<std::uint32_t>> send_local_indices_by_neighbor;
  std::vector<std::vector<std::uint32_t>> recv_local_indices_by_neighbor;
  std::uint64_t send_bytes = 0;
  std::uint64_t recv_bytes = 0;
};

enum class LocalIndexResidency : std::uint8_t {
  kOwned = 0,
  kGhost = 1,
};

struct LocalGhostDescriptor {
  LocalIndexResidency residency = LocalIndexResidency::kOwned;
  int owning_rank = 0;
};

[[nodiscard]] GhostExchangePlan buildGhostExchangePlan(
    int world_rank,
    std::span<const LocalGhostDescriptor> local_ghost_descriptors,
    std::size_t bytes_per_ghost);

[[nodiscard]] GhostExchangePlan buildGhostExchangePlan(
    int world_rank,
    std::span<const int> ghost_owner_rank_by_local_index,
    std::size_t bytes_per_ghost);

struct ReductionAgreement {
  double deterministic_sum = 0.0;
  double reference_sum = 0.0;
  double absolute_error = 0.0;
  double relative_error = 0.0;
};

[[nodiscard]] double deterministicRankOrderedSum(std::span<const double> per_rank_values);
[[nodiscard]] ReductionAgreement compareReductionAgreement(
    std::span<const double> per_rank_values,
    double measured_sum);

struct RankConfigDigest {
  int world_rank = 0;
  std::uint64_t normalized_config_hash = 0;
  int mpi_ranks_expected = 1;
  bool deterministic_reduction = true;
};

struct RankConfigConsensus {
  bool normalized_config_hash_match = true;
  bool mpi_ranks_expected_match = true;
  bool deterministic_reduction_match = true;
  std::vector<int> mismatched_ranks;

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
  std::uint32_t schema_version = 1;
  std::uint64_t decomposition_epoch = 0;
  int world_size = 1;
  std::vector<int> owning_rank_by_item;

  [[nodiscard]] std::string serialize() const;
  [[nodiscard]] static DistributedRestartState deserialize(const std::string& encoded);
};

void recordDistributedProfiling(
    core::ProfilerSession* profiler,
    const LoadBalanceMetrics& metrics,
    std::uint64_t ghost_exchange_send_bytes,
    std::uint64_t ghost_exchange_recv_bytes);

class MpiContext {
 public:
  MpiContext();

  [[nodiscard]] bool isEnabled() const noexcept;
  [[nodiscard]] int worldSize() const noexcept;
  [[nodiscard]] int worldRank() const noexcept;

  [[nodiscard]] double allreduceSumDouble(double local_value) const;
  [[nodiscard]] std::uint64_t allreduceSumUint64(std::uint64_t local_value) const;

 private:
  bool m_is_enabled = false;
  int m_world_size = 1;
  int m_world_rank = 0;
};

}  // namespace cosmosim::parallel
