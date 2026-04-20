#include "cosmosim/parallel/distributed_memory.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

#include "cosmosim/core/build_config.hpp"

#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
#include <mpi.h>
#endif

namespace cosmosim::parallel {
namespace {

[[nodiscard]] double clampUnit(double value) {
  if (value <= 0.0) {
    return 0.0;
  }
  if (value >= 1.0) {
    return std::nextafter(1.0, 0.0);
  }
  return value;
}

[[nodiscard]] std::uint32_t quantize10bit(double coordinate, double min_coord, double max_coord) {
  const double extent = max_coord - min_coord;
  if (!(extent > 0.0)) {
    throw std::invalid_argument("decomposition domain extents must be positive");
  }
  const double normalized = clampUnit((coordinate - min_coord) / extent);
  constexpr double k_scale = 1024.0;
  const double scaled = std::floor(normalized * k_scale);
  const auto q = static_cast<std::uint32_t>(scaled);
  return std::min<std::uint32_t>(q, 1023U);
}

[[nodiscard]] std::uint64_t expandBits3d(std::uint32_t x) {
  std::uint64_t value = x & 0x3ffU;
  value = (value | (value << 16)) & 0x30000ffU;
  value = (value | (value << 8)) & 0x300f00fU;
  value = (value | (value << 4)) & 0x30c30c3U;
  value = (value | (value << 2)) & 0x9249249U;
  return value;
}

[[nodiscard]] std::uint64_t mortonKey3d(std::uint32_t x, std::uint32_t y, std::uint32_t z) {
  return (expandBits3d(z) << 2U) | (expandBits3d(y) << 1U) | expandBits3d(x);
}

[[nodiscard]] double weightedLoad(const DecompositionItem& item, const DecompositionConfig& config) {
  const double owned_particle_term =
      (item.kind == DecompositionEntityKind::kParticle) ? 1.0 : 0.0;
  const double active_target_term = static_cast<double>(item.active_target_count_recent);
  const double remote_tree_term = static_cast<double>(item.remote_tree_interactions_recent);
  const double work_term = std::max(0.0, item.work_units);
  const double memory_term = static_cast<double>(item.memory_bytes);
  return config.owned_particle_weight * owned_particle_term +
      config.active_target_weight * active_target_term +
      config.remote_tree_interaction_weight * remote_tree_term +
      config.work_weight * work_term + config.memory_weight * memory_term;
}

template <typename T>
void appendPod(std::vector<std::uint8_t>& bytes, const T& value) {
  const std::size_t old_size = bytes.size();
  bytes.resize(old_size + sizeof(T));
  std::memcpy(bytes.data() + old_size, &value, sizeof(T));
}

template <typename T>
[[nodiscard]] T readPod(const std::vector<std::uint8_t>& bytes, std::size_t* offset) {
  if (*offset + sizeof(T) > bytes.size()) {
    throw std::runtime_error("ghost buffer decode overflow");
  }
  T value{};
  std::memcpy(&value, bytes.data() + *offset, sizeof(T));
  *offset += sizeof(T);
  return value;
}

[[nodiscard]] std::vector<std::string> splitLines(const std::string& text) {
  std::vector<std::string> lines;
  std::stringstream stream(text);
  std::string line;
  while (std::getline(stream, line)) {
    lines.push_back(line);
  }
  return lines;
}

[[nodiscard]] double absoluteValue(double value) {
  return (value < 0.0) ? -value : value;
}

[[nodiscard]] double stableRelativeError(double measured, double reference, double absolute_error) {
  const double denom = std::max(absoluteValue(reference), std::numeric_limits<double>::min());
  return absolute_error / denom;
}

[[nodiscard]] std::string rankConfigValueString(std::uint64_t value) {
  return std::to_string(value);
}

[[nodiscard]] std::string rankConfigValueString(int value) {
  return std::to_string(value);
}

[[nodiscard]] std::string rankConfigValueString(bool value) {
  return value ? "true" : "false";
}

void appendRankConfigMismatch(
    RankConfigConsensus* consensus,
    RankConfigMismatchProperty property,
    int baseline_rank,
    int rank,
    std::string baseline_value,
    std::string rank_value) {
  consensus->mismatches.push_back(RankConfigMismatch{
      .property = property,
      .baseline_rank = baseline_rank,
      .rank = rank,
      .baseline_value = std::move(baseline_value),
      .rank_value = std::move(rank_value),
  });
}

[[nodiscard]] constexpr std::size_t ghostExchangeRecordBytes() {
  return sizeof(std::uint64_t) + sizeof(double) * 3U;
}

void validateTransferDescriptor(
    const GhostTransferDescriptor& descriptor,
    GhostTransferRole expected_role,
    int expected_peer_rank,
    std::size_t expected_neighbor_slot,
    std::span<const std::uint32_t> expected_indices) {
  if (descriptor.role != expected_role) {
    throw std::invalid_argument("ghost transfer descriptor role does not match container");
  }
  if (descriptor.peer_rank != expected_peer_rank) {
    throw std::invalid_argument("ghost transfer descriptor peer_rank does not match neighbor slot");
  }
  if (descriptor.neighbor_slot != expected_neighbor_slot) {
    throw std::invalid_argument("ghost transfer descriptor neighbor_slot mismatch");
  }
  if (descriptor.local_indices.size() != expected_indices.size() ||
      !std::equal(descriptor.local_indices.begin(), descriptor.local_indices.end(), expected_indices.begin())) {
    throw std::invalid_argument("ghost transfer descriptor indices drift from canonical plan indices");
  }
  if (expected_role == GhostTransferRole::kOutboundSend &&
      descriptor.intent != GhostTransferIntent::kGhostRefreshRequest &&
      descriptor.intent != GhostTransferIntent::kOwnershipMigrationSend) {
    throw std::invalid_argument("outbound transfer intent must be ghost refresh request or migration send");
  }
  if (expected_role == GhostTransferRole::kInboundReceive &&
      descriptor.intent != GhostTransferIntent::kGhostRefreshReceiveStaging &&
      descriptor.intent != GhostTransferIntent::kOwnershipMigrationReceiveStaging) {
    throw std::invalid_argument("inbound transfer intent must be receive-staging intent");
  }
  if (descriptor.local_indices.empty() &&
      !(expected_role == GhostTransferRole::kOutboundSend &&
        descriptor.intent == GhostTransferIntent::kGhostRefreshRequest)) {
    throw std::invalid_argument("ghost transfer descriptor local_indices must be non-empty");
  }
  if (descriptor.intent == GhostTransferIntent::kGhostRefreshRequest ||
      descriptor.intent == GhostTransferIntent::kGhostRefreshReceiveStaging) {
    if (descriptor.expected_post_transfer_residency != LocalIndexResidency::kGhost) {
      throw std::invalid_argument("ghost refresh transfers must keep ghost post-transfer residency");
    }
  }
}

}  // namespace

DecompositionPlan buildMortonSfcDecomposition(std::span<const DecompositionItem> items, const DecompositionConfig& config) {
  if (config.world_size <= 0) {
    throw std::invalid_argument("world_size must be positive");
  }
  if (config.owned_particle_weight < 0.0 || config.active_target_weight < 0.0 ||
      config.remote_tree_interaction_weight < 0.0 || config.work_weight < 0.0 || config.memory_weight < 0.0) {
    throw std::invalid_argument("decomposition weights must be non-negative");
  }

  struct KeyedItem {
    std::size_t index = 0;
    std::uint64_t morton_key = 0;
    double weighted_load = 0.0;
  };

  std::vector<KeyedItem> keyed(items.size());
  for (std::size_t i = 0; i < items.size(); ++i) {
    const std::uint32_t qx = quantize10bit(items[i].x_comov, config.domain_x_min_comov, config.domain_x_max_comov);
    const std::uint32_t qy = quantize10bit(items[i].y_comov, config.domain_y_min_comov, config.domain_y_max_comov);
    const std::uint32_t qz = quantize10bit(items[i].z_comov, config.domain_z_min_comov, config.domain_z_max_comov);
    keyed[i] = KeyedItem{
        .index = i,
        .morton_key = mortonKey3d(qx, qy, qz),
        .weighted_load = weightedLoad(items[i], config),
    };
  }

  std::stable_sort(keyed.begin(), keyed.end(), [&](const KeyedItem& a, const KeyedItem& b) {
    if (a.morton_key != b.morton_key) {
      return a.morton_key < b.morton_key;
    }
    return items[a.index].entity_id < items[b.index].entity_id;
  });

  DecompositionPlan plan;
  plan.owning_rank_by_item.assign(items.size(), 0);
  plan.sorted_indices.resize(items.size());
  plan.ranges_by_rank.assign(static_cast<std::size_t>(config.world_size), RankRange{});
  plan.metrics.weighted_load_by_rank.assign(static_cast<std::size_t>(config.world_size), 0.0);
  plan.metrics.memory_bytes_by_rank.assign(static_cast<std::size_t>(config.world_size), 0ULL);
  plan.metrics.owned_particles_by_rank.assign(static_cast<std::size_t>(config.world_size), 0ULL);
  plan.metrics.active_targets_by_rank.assign(static_cast<std::size_t>(config.world_size), 0ULL);
  plan.metrics.remote_tree_interactions_by_rank.assign(static_cast<std::size_t>(config.world_size), 0ULL);

  const double total_load = std::accumulate(
      keyed.begin(), keyed.end(), 0.0, [](double acc, const KeyedItem& entry) { return acc + entry.weighted_load; });
  const double target_per_rank = (config.world_size > 0) ? total_load / static_cast<double>(config.world_size) : 0.0;

  int current_rank = 0;
  double current_prefix_load = 0.0;
  std::size_t rank_begin = 0;

  for (std::size_t sorted_pos = 0; sorted_pos < keyed.size(); ++sorted_pos) {
    const std::size_t original_index = keyed[sorted_pos].index;
    plan.sorted_indices[sorted_pos] = original_index;
    plan.owning_rank_by_item[original_index] = current_rank;
    plan.metrics.weighted_load_by_rank[static_cast<std::size_t>(current_rank)] += keyed[sorted_pos].weighted_load;
    plan.metrics.memory_bytes_by_rank[static_cast<std::size_t>(current_rank)] += items[original_index].memory_bytes;
    if (items[original_index].kind == DecompositionEntityKind::kParticle) {
      ++plan.metrics.owned_particles_by_rank[static_cast<std::size_t>(current_rank)];
    }
    plan.metrics.active_targets_by_rank[static_cast<std::size_t>(current_rank)] +=
        items[original_index].active_target_count_recent;
    plan.metrics.remote_tree_interactions_by_rank[static_cast<std::size_t>(current_rank)] +=
        items[original_index].remote_tree_interactions_recent;
    const double item_load = keyed[sorted_pos].weighted_load;
    current_prefix_load += item_load;

    if (current_rank + 1 >= config.world_size) {
      continue;
    }
    if (sorted_pos + 1 >= keyed.size()) {
      continue;
    }

    const std::size_t items_remaining = keyed.size() - (sorted_pos + 1U);
    const std::size_t ranks_remaining = static_cast<std::size_t>(config.world_size - (current_rank + 1));
    if (items_remaining < ranks_remaining) {
      continue;
    }

    const double next_target_prefix = target_per_rank * static_cast<double>(current_rank + 1);
    const bool crossed_target = current_prefix_load >= next_target_prefix;
    const bool rank_has_multiple_items = sorted_pos > rank_begin;
    if (!crossed_target || !rank_has_multiple_items) {
      continue;
    }

    const double before_distance = std::abs((current_prefix_load - item_load) - next_target_prefix);
    const double after_distance = std::abs(current_prefix_load - next_target_prefix);
    const bool cut_before_current = before_distance < after_distance;
    if (cut_before_current) {
      plan.ranges_by_rank[static_cast<std::size_t>(current_rank)] = RankRange{
          .begin_sorted = rank_begin,
          .end_sorted = sorted_pos};
      ++current_rank;
      rank_begin = sorted_pos;
      plan.owning_rank_by_item[original_index] = current_rank;
      plan.metrics.weighted_load_by_rank[static_cast<std::size_t>(current_rank)] += item_load;
      plan.metrics.weighted_load_by_rank[static_cast<std::size_t>(current_rank - 1)] -= item_load;
      plan.metrics.memory_bytes_by_rank[static_cast<std::size_t>(current_rank)] += items[original_index].memory_bytes;
      plan.metrics.memory_bytes_by_rank[static_cast<std::size_t>(current_rank - 1)] -= items[original_index].memory_bytes;
      if (items[original_index].kind == DecompositionEntityKind::kParticle) {
        ++plan.metrics.owned_particles_by_rank[static_cast<std::size_t>(current_rank)];
        --plan.metrics.owned_particles_by_rank[static_cast<std::size_t>(current_rank - 1)];
      }
      plan.metrics.active_targets_by_rank[static_cast<std::size_t>(current_rank)] +=
          items[original_index].active_target_count_recent;
      plan.metrics.active_targets_by_rank[static_cast<std::size_t>(current_rank - 1)] -=
          items[original_index].active_target_count_recent;
      plan.metrics.remote_tree_interactions_by_rank[static_cast<std::size_t>(current_rank)] +=
          items[original_index].remote_tree_interactions_recent;
      plan.metrics.remote_tree_interactions_by_rank[static_cast<std::size_t>(current_rank - 1)] -=
          items[original_index].remote_tree_interactions_recent;
      continue;
    }

    plan.ranges_by_rank[static_cast<std::size_t>(current_rank)] = RankRange{
        .begin_sorted = rank_begin,
        .end_sorted = sorted_pos + 1U};
    ++current_rank;
    rank_begin = sorted_pos + 1U;
  }

  plan.ranges_by_rank[static_cast<std::size_t>(current_rank)] = RankRange{.begin_sorted = rank_begin, .end_sorted = keyed.size()};
  for (int rank = current_rank + 1; rank < config.world_size; ++rank) {
    plan.ranges_by_rank[static_cast<std::size_t>(rank)] = RankRange{.begin_sorted = keyed.size(), .end_sorted = keyed.size()};
  }

  const auto max_load_it = std::max_element(plan.metrics.weighted_load_by_rank.begin(), plan.metrics.weighted_load_by_rank.end());
  plan.metrics.max_weighted_load =
      (max_load_it == plan.metrics.weighted_load_by_rank.end()) ? 0.0 : *max_load_it;
  plan.metrics.mean_weighted_load =
      plan.metrics.weighted_load_by_rank.empty()
          ? 0.0
          : (std::accumulate(plan.metrics.weighted_load_by_rank.begin(), plan.metrics.weighted_load_by_rank.end(), 0.0) /
             static_cast<double>(plan.metrics.weighted_load_by_rank.size()));
  plan.metrics.weighted_imbalance_ratio =
      (plan.metrics.mean_weighted_load > 0.0) ? (plan.metrics.max_weighted_load / plan.metrics.mean_weighted_load) : 0.0;

  plan.metrics.total_memory_bytes = std::accumulate(
      plan.metrics.memory_bytes_by_rank.begin(), plan.metrics.memory_bytes_by_rank.end(), 0ULL);
  const auto max_mem_it = std::max_element(plan.metrics.memory_bytes_by_rank.begin(), plan.metrics.memory_bytes_by_rank.end());
  plan.metrics.max_memory_bytes = (max_mem_it == plan.metrics.memory_bytes_by_rank.end()) ? 0ULL : *max_mem_it;
  const double mean_memory = plan.metrics.memory_bytes_by_rank.empty()
                                 ? 0.0
                                 : (static_cast<double>(plan.metrics.total_memory_bytes) /
                                    static_cast<double>(plan.metrics.memory_bytes_by_rank.size()));
  plan.metrics.memory_imbalance_ratio = (mean_memory > 0.0) ? (static_cast<double>(plan.metrics.max_memory_bytes) / mean_memory) : 0.0;

  return plan;
}

GhostExchangePlan buildGhostExchangePlan(
    int world_rank,
    std::span<const LocalGhostDescriptor> local_ghost_descriptors,
    std::size_t bytes_per_ghost) {
  if (bytes_per_ghost == 0) {
    throw std::invalid_argument("bytes_per_ghost must be positive");
  }
  GhostExchangePlan plan;
  std::vector<int> owners;
  owners.reserve(local_ghost_descriptors.size());

  for (const LocalGhostDescriptor descriptor : local_ghost_descriptors) {
    if (descriptor.owning_rank < 0) {
      throw std::invalid_argument("ghost owner rank must be non-negative");
    }
    if (descriptor.residency == LocalIndexResidency::kOwned) {
      if (descriptor.owning_rank != world_rank) {
        throw std::invalid_argument("owned local index must have world_rank ownership");
      }
      continue;
    }
    if (descriptor.owning_rank == world_rank) {
      throw std::invalid_argument("ghost local index cannot be owned by world_rank");
    }
    owners.push_back(descriptor.owning_rank);
  }

  std::sort(owners.begin(), owners.end());
  owners.erase(std::unique(owners.begin(), owners.end()), owners.end());

  plan.neighbor_ranks = owners;
  plan.send_local_indices_by_neighbor.assign(owners.size(), {});
  plan.recv_local_indices_by_neighbor.assign(owners.size(), {});
  plan.outbound_transfers.assign(owners.size(), {});
  plan.inbound_transfers.assign(owners.size(), {});

  for (std::uint32_t local_index = 0; local_index < local_ghost_descriptors.size(); ++local_index) {
    const LocalGhostDescriptor descriptor = local_ghost_descriptors[local_index];
    if (descriptor.residency == LocalIndexResidency::kOwned) {
      continue;
    }
    const auto it = std::lower_bound(owners.begin(), owners.end(), descriptor.owning_rank);
    if (it == owners.end() || *it != descriptor.owning_rank) {
      throw std::logic_error("owner rank map mismatch");
    }
    const std::size_t neighbor_slot = static_cast<std::size_t>(std::distance(owners.begin(), it));
    plan.recv_local_indices_by_neighbor[neighbor_slot].push_back(local_index);
  }

  for (std::size_t i = 0; i < owners.size(); ++i) {
    // Descriptor-only planning can identify local ghost import slots, but not the peer-owned
    // source rows that must be exported back. Keep outbound payload indices empty rather than
    // pretending the local ghost rows themselves are valid send sources.
    plan.send_local_indices_by_neighbor[i].clear();
    plan.outbound_transfers[i] = GhostTransferDescriptor{
        .role = GhostTransferRole::kOutboundSend,
        .intent = GhostTransferIntent::kGhostRefreshRequest,
        .peer_rank = owners[i],
        .neighbor_slot = i,
        .expected_post_transfer_residency = LocalIndexResidency::kGhost,
        .local_indices = plan.send_local_indices_by_neighbor[i],
    };
    plan.inbound_transfers[i] = GhostTransferDescriptor{
        .role = GhostTransferRole::kInboundReceive,
        .intent = GhostTransferIntent::kGhostRefreshReceiveStaging,
        .peer_rank = owners[i],
        .neighbor_slot = i,
        .expected_post_transfer_residency = LocalIndexResidency::kGhost,
        .local_indices = plan.recv_local_indices_by_neighbor[i],
    };
    plan.recv_bytes +=
        static_cast<std::uint64_t>(plan.recv_local_indices_by_neighbor[i].size()) * bytes_per_ghost;
    plan.send_bytes +=
        static_cast<std::uint64_t>(plan.send_local_indices_by_neighbor[i].size()) * bytes_per_ghost;
  }

  validateGhostExchangePlan(plan);
  return plan;
}

GhostExchangePlan buildGhostExchangePlan(
    int world_rank,
    std::span<const int> ghost_owner_rank_by_local_index,
    std::size_t bytes_per_ghost) {
  std::vector<LocalGhostDescriptor> descriptors;
  descriptors.reserve(ghost_owner_rank_by_local_index.size());
  for (const int owner_rank : ghost_owner_rank_by_local_index) {
    descriptors.push_back(LocalGhostDescriptor{
        .residency = (owner_rank == world_rank) ? LocalIndexResidency::kOwned : LocalIndexResidency::kGhost,
        .owning_rank = owner_rank,
    });
  }
  return buildGhostExchangePlan(world_rank, descriptors, bytes_per_ghost);
}

void validateGhostExchangePlan(const GhostExchangePlan& plan) {
  const std::size_t neighbor_count = plan.neighbor_ranks.size();
  if (plan.send_local_indices_by_neighbor.size() != neighbor_count ||
      plan.recv_local_indices_by_neighbor.size() != neighbor_count ||
      plan.outbound_transfers.size() != neighbor_count ||
      plan.inbound_transfers.size() != neighbor_count) {
    throw std::invalid_argument("ghost exchange plan containers must have matching neighbor counts");
  }
  for (std::size_t i = 0; i < neighbor_count; ++i) {
    if (i > 0 && plan.neighbor_ranks[i - 1] >= plan.neighbor_ranks[i]) {
      throw std::invalid_argument("ghost exchange plan neighbor_ranks must be strictly increasing");
    }
    validateTransferDescriptor(
        plan.outbound_transfers[i],
        GhostTransferRole::kOutboundSend,
        plan.neighbor_ranks[i],
        i,
        plan.send_local_indices_by_neighbor[i]);
    validateTransferDescriptor(
        plan.inbound_transfers[i],
        GhostTransferRole::kInboundReceive,
        plan.neighbor_ranks[i],
        i,
        plan.recv_local_indices_by_neighbor[i]);
  }
}

double deterministicRankOrderedSum(std::span<const double> per_rank_values) {
  long double sum = 0.0;
  for (const double value : per_rank_values) {
    sum += static_cast<long double>(value);
  }
  return static_cast<double>(sum);
}

ReductionAgreement compareReductionAgreement(
    std::span<const double> per_rank_values,
    double measured_sum) {
  const double deterministic_baseline_sum = deterministicRankOrderedSum(per_rank_values);
  const double absolute_error = absoluteValue(measured_sum - deterministic_baseline_sum);
  return ReductionAgreement{
      .deterministic_baseline_sum = deterministic_baseline_sum,
      .measured_sum = measured_sum,
      .absolute_error = absolute_error,
      .relative_error = stableRelativeError(measured_sum, deterministic_baseline_sum, absolute_error),
  };
}

bool satisfiesReductionAgreement(
    const ReductionAgreement& agreement,
    const ReductionAgreementPolicy& policy) {
  if (policy.absolute_tolerance < 0.0 || policy.relative_tolerance < 0.0) {
    throw std::invalid_argument("reduction agreement tolerances must be non-negative");
  }
  const bool absolute_ok = agreement.absolute_error <= policy.absolute_tolerance;
  const bool relative_ok = agreement.relative_error <= policy.relative_tolerance;
  switch (policy.mode) {
    case ReductionAgreementMode::kAbsoluteOnly:
      return absolute_ok;
    case ReductionAgreementMode::kRelativeOnly:
      return relative_ok;
    case ReductionAgreementMode::kAbsoluteAndRelative:
      return absolute_ok && relative_ok;
    case ReductionAgreementMode::kAbsoluteOrRelative:
      return absolute_ok || relative_ok;
  }
  throw std::invalid_argument("unknown reduction agreement mode");
}

bool RankConfigConsensus::allConsistent() const noexcept {
  return normalized_config_hash_match && mpi_ranks_expected_match && deterministic_reduction_match;
}

RankConfigConsensus evaluateRankConfigConsensus(std::span<const RankConfigDigest> digests) {
  RankConfigConsensus consensus;
  if (digests.empty()) {
    return consensus;
  }

  const RankConfigDigest baseline = digests.front();
  for (const RankConfigDigest& digest : digests) {
    bool rank_matches = true;
    if (digest.normalized_config_hash != baseline.normalized_config_hash) {
      consensus.normalized_config_hash_match = false;
      rank_matches = false;
      appendRankConfigMismatch(
          &consensus,
          RankConfigMismatchProperty::kNormalizedConfigHash,
          baseline.world_rank,
          digest.world_rank,
          rankConfigValueString(baseline.normalized_config_hash),
          rankConfigValueString(digest.normalized_config_hash));
    }
    if (digest.mpi_ranks_expected != baseline.mpi_ranks_expected) {
      consensus.mpi_ranks_expected_match = false;
      rank_matches = false;
      appendRankConfigMismatch(
          &consensus,
          RankConfigMismatchProperty::kMpiRanksExpected,
          baseline.world_rank,
          digest.world_rank,
          rankConfigValueString(baseline.mpi_ranks_expected),
          rankConfigValueString(digest.mpi_ranks_expected));
    }
    if (digest.deterministic_reduction != baseline.deterministic_reduction) {
      consensus.deterministic_reduction_match = false;
      rank_matches = false;
      appendRankConfigMismatch(
          &consensus,
          RankConfigMismatchProperty::kDeterministicReduction,
          baseline.world_rank,
          digest.world_rank,
          rankConfigValueString(baseline.deterministic_reduction),
          rankConfigValueString(digest.deterministic_reduction));
    }
    if (!rank_matches) {
      consensus.mismatched_ranks.push_back(digest.world_rank);
    }
  }
  return consensus;
}

bool GhostExchangeBufferSoA::isConsistent() const noexcept {
  return entity_id.size() == density_code.size() && entity_id.size() == velocity_x_code.size() &&
         entity_id.size() == pressure_code.size();
}

std::size_t GhostExchangeBufferSoA::size() const noexcept { return entity_id.size(); }

void GhostExchangeBuffer::clear() { m_bytes.clear(); }

std::size_t GhostExchangeBuffer::byteSize() const noexcept { return m_bytes.size(); }

void GhostExchangeBuffer::packFrom(const GhostExchangeBufferSoA& source, std::span<const std::uint32_t> local_indices) {
  if (!source.isConsistent()) {
    throw std::invalid_argument("ghost source SoA fields must have matching sizes");
  }

  m_bytes.clear();
  appendPod<std::uint64_t>(m_bytes, static_cast<std::uint64_t>(local_indices.size()));

  for (const std::uint32_t index : local_indices) {
    if (index >= source.size()) {
      throw std::out_of_range("ghost pack local index out of range");
    }

    appendPod<std::uint64_t>(m_bytes, source.entity_id[index]);
    appendPod<double>(m_bytes, source.density_code[index]);
    appendPod<double>(m_bytes, source.velocity_x_code[index]);
    appendPod<double>(m_bytes, source.pressure_code[index]);
  }
}

void GhostExchangeBuffer::unpackAppendTo(GhostExchangeBufferSoA& destination) const {
  if (!destination.isConsistent()) {
    throw std::invalid_argument("ghost destination SoA fields must have matching sizes");
  }
  if (m_bytes.size() < sizeof(std::uint64_t)) {
    throw std::runtime_error("ghost buffer is too small");
  }

  std::size_t offset = 0;
  const std::uint64_t count = readPod<std::uint64_t>(m_bytes, &offset);
  const std::uint64_t expected_payload_bytes =
      static_cast<std::uint64_t>(sizeof(std::uint64_t)) + count * static_cast<std::uint64_t>(ghostExchangeRecordBytes());
  if (expected_payload_bytes != static_cast<std::uint64_t>(m_bytes.size())) {
    throw std::runtime_error("ghost buffer payload shape does not match encoded count");
  }

  destination.entity_id.reserve(destination.entity_id.size() + static_cast<std::size_t>(count));
  destination.density_code.reserve(destination.density_code.size() + static_cast<std::size_t>(count));
  destination.velocity_x_code.reserve(destination.velocity_x_code.size() + static_cast<std::size_t>(count));
  destination.pressure_code.reserve(destination.pressure_code.size() + static_cast<std::size_t>(count));

  for (std::uint64_t i = 0; i < count; ++i) {
    destination.entity_id.push_back(readPod<std::uint64_t>(m_bytes, &offset));
    destination.density_code.push_back(readPod<double>(m_bytes, &offset));
    destination.velocity_x_code.push_back(readPod<double>(m_bytes, &offset));
    destination.pressure_code.push_back(readPod<double>(m_bytes, &offset));
  }

  if (offset != m_bytes.size()) {
    throw std::runtime_error("ghost buffer decode found trailing bytes");
  }
}

std::string DistributedRestartState::serialize() const {
  std::ostringstream stream;
  stream << "schema_version=" << schema_version << '\n';
  stream << "decomposition_epoch=" << decomposition_epoch << '\n';
  stream << "world_size=" << world_size << '\n';
  stream << "pm_grid_nx=" << pm_grid_nx << '\n';
  stream << "pm_grid_ny=" << pm_grid_ny << '\n';
  stream << "pm_grid_nz=" << pm_grid_nz << '\n';
  stream << "pm_decomposition_mode=" << pm_decomposition_mode << '\n';
  stream << "gravity_kick_opportunity=" << gravity_kick_opportunity << '\n';
  stream << "pm_update_cadence_steps=" << pm_update_cadence_steps << '\n';
  stream << "long_range_field_version=" << long_range_field_version << '\n';
  stream << "last_long_range_refresh_opportunity=" << last_long_range_refresh_opportunity << '\n';
  stream << "long_range_field_built_step_index=" << long_range_field_built_step_index << '\n';
  stream << "long_range_field_built_scale_factor=" << long_range_field_built_scale_factor << '\n';
  stream << "long_range_restart_policy=" << long_range_restart_policy << '\n';
  stream << "item_count=" << owning_rank_by_item.size() << '\n';
  for (std::size_t i = 0; i < owning_rank_by_item.size(); ++i) {
    stream << "rank[" << i << "]=" << owning_rank_by_item[i] << '\n';
  }
  stream << "pm_slab_rank_count=" << pm_slab_begin_x_by_rank.size() << '\n';
  for (std::size_t rank = 0; rank < pm_slab_begin_x_by_rank.size(); ++rank) {
    stream << "pm_slab_begin_x[" << rank << "]=" << pm_slab_begin_x_by_rank[rank] << '\n';
    stream << "pm_slab_end_x[" << rank << "]=" << pm_slab_end_x_by_rank[rank] << '\n';
  }
  return stream.str();
}

DistributedRestartState DistributedRestartState::deserialize(const std::string& encoded) {
  DistributedRestartState state;
  const std::vector<std::string> lines = splitLines(encoded);
  std::size_t expected_item_count = 0;
  std::vector<bool> seen_rank_entry;
  std::size_t expected_slab_rank_count = 0;
  std::vector<bool> seen_slab_begin;
  std::vector<bool> seen_slab_end;

  for (const std::string& line : lines) {
    if (line.empty()) {
      continue;
    }

    const std::size_t eq = line.find('=');
    if (eq == std::string::npos) {
      throw std::invalid_argument("invalid restart encoding line");
    }
    const std::string key = line.substr(0, eq);
    const std::string value = line.substr(eq + 1);

    if (key == "schema_version") {
      state.schema_version = static_cast<std::uint32_t>(std::stoul(value));
    } else if (key == "decomposition_epoch") {
      state.decomposition_epoch = std::stoull(value);
    } else if (key == "world_size") {
      state.world_size = std::stoi(value);
    } else if (key == "pm_grid_nx") {
      state.pm_grid_nx = static_cast<std::size_t>(std::stoull(value));
    } else if (key == "pm_grid_ny") {
      state.pm_grid_ny = static_cast<std::size_t>(std::stoull(value));
    } else if (key == "pm_grid_nz") {
      state.pm_grid_nz = static_cast<std::size_t>(std::stoull(value));
    } else if (key == "pm_decomposition_mode") {
      state.pm_decomposition_mode = value;
    } else if (key == "gravity_kick_opportunity") {
      state.gravity_kick_opportunity = std::stoull(value);
    } else if (key == "pm_update_cadence_steps") {
      state.pm_update_cadence_steps = std::stoull(value);
    } else if (key == "long_range_field_version") {
      state.long_range_field_version = std::stoull(value);
    } else if (key == "last_long_range_refresh_opportunity") {
      state.last_long_range_refresh_opportunity = std::stoull(value);
    } else if (key == "long_range_field_built_step_index") {
      state.long_range_field_built_step_index = std::stoull(value);
    } else if (key == "long_range_field_built_scale_factor") {
      state.long_range_field_built_scale_factor = std::stod(value);
    } else if (key == "long_range_restart_policy") {
      state.long_range_restart_policy = value;
    } else if (key == "item_count") {
      expected_item_count = static_cast<std::size_t>(std::stoull(value));
      state.owning_rank_by_item.assign(expected_item_count, 0);
      seen_rank_entry.assign(expected_item_count, false);
    } else if (key == "pm_slab_rank_count") {
      expected_slab_rank_count = static_cast<std::size_t>(std::stoull(value));
      state.pm_slab_begin_x_by_rank.assign(expected_slab_rank_count, 0);
      state.pm_slab_end_x_by_rank.assign(expected_slab_rank_count, 0);
      seen_slab_begin.assign(expected_slab_rank_count, false);
      seen_slab_end.assign(expected_slab_rank_count, false);
    } else if (key.rfind("rank[", 0) == 0) {
      const std::size_t open = key.find('[');
      const std::size_t close = key.find(']');
      if (open == std::string::npos || close == std::string::npos || close <= open + 1) {
        throw std::invalid_argument("invalid rank entry in restart encoding");
      }
      const std::size_t index = static_cast<std::size_t>(std::stoull(key.substr(open + 1, close - open - 1)));
      if (index >= state.owning_rank_by_item.size()) {
        throw std::out_of_range("restart rank index out of bounds");
      }
      if (seen_rank_entry[index]) {
        throw std::invalid_argument("duplicate restart rank entry");
      }
      state.owning_rank_by_item[index] = std::stoi(value);
      seen_rank_entry[index] = true;
    } else if (key.rfind("pm_slab_begin_x[", 0) == 0 || key.rfind("pm_slab_end_x[", 0) == 0) {
      const bool is_begin = key.rfind("pm_slab_begin_x[", 0) == 0;
      const std::size_t open = key.find('[');
      const std::size_t close = key.find(']');
      if (open == std::string::npos || close == std::string::npos || close <= open + 1) {
        throw std::invalid_argument("invalid PM slab entry in restart encoding");
      }
      const std::size_t rank_index = static_cast<std::size_t>(std::stoull(key.substr(open + 1, close - open - 1)));
      if (rank_index >= expected_slab_rank_count) {
        throw std::out_of_range("restart PM slab rank index out of bounds");
      }
      if (is_begin) {
        if (seen_slab_begin[rank_index]) {
          throw std::invalid_argument("duplicate PM slab begin entry");
        }
        state.pm_slab_begin_x_by_rank[rank_index] = static_cast<std::size_t>(std::stoull(value));
        seen_slab_begin[rank_index] = true;
      } else {
        if (seen_slab_end[rank_index]) {
          throw std::invalid_argument("duplicate PM slab end entry");
        }
        state.pm_slab_end_x_by_rank[rank_index] = static_cast<std::size_t>(std::stoull(value));
        seen_slab_end[rank_index] = true;
      }
    }
  }

  if (state.owning_rank_by_item.size() != expected_item_count) {
    throw std::runtime_error("restart decode item count mismatch");
  }
  if (state.world_size <= 0) {
    throw std::invalid_argument("restart world_size must be positive");
  }
  if (state.pm_update_cadence_steps == 0) {
    throw std::invalid_argument("restart PM cadence must be >= 1");
  }
  if (state.last_long_range_refresh_opportunity > state.gravity_kick_opportunity) {
    throw std::invalid_argument("restart cadence state is inconsistent: last refresh opportunity exceeds kick opportunity");
  }
  if (state.long_range_field_version == 0 && state.last_long_range_refresh_opportunity != 0) {
    throw std::invalid_argument("restart cadence state is inconsistent: non-zero refresh opportunity with zero field version");
  }
  if (state.long_range_restart_policy != "deterministic_rebuild") {
    throw std::invalid_argument("restart long-range policy is unsupported: " + state.long_range_restart_policy);
  }
  if (state.schema_version >= 2) {
    if (state.pm_grid_nx == 0 || state.pm_grid_ny == 0 || state.pm_grid_nz == 0) {
      throw std::invalid_argument("restart PM grid dimensions must be > 0 for schema_version >= 2");
    }
    if (state.pm_decomposition_mode.empty()) {
      throw std::invalid_argument("restart PM decomposition mode must be non-empty");
    }
    if (expected_slab_rank_count != static_cast<std::size_t>(state.world_size)) {
      throw std::invalid_argument("restart PM slab rank count must match world_size");
    }
    for (std::size_t rank = 0; rank < expected_slab_rank_count; ++rank) {
      if (!seen_slab_begin[rank] || !seen_slab_end[rank]) {
        throw std::runtime_error("restart decode missing PM slab ownership entry");
      }
      if (state.pm_slab_end_x_by_rank[rank] < state.pm_slab_begin_x_by_rank[rank]) {
        throw std::invalid_argument("restart PM slab end_x must be >= begin_x");
      }
    }
  }
  for (bool seen : seen_rank_entry) {
    if (!seen) {
      throw std::runtime_error("restart decode missing ownership entry");
    }
  }
  return state;
}

DistributedRestartCompatibilityReport evaluateDistributedRestartCompatibility(
    const DistributedRestartState& restart_state,
    const DistributedExecutionTopology& runtime_topology) {
  DistributedRestartCompatibilityReport report;
  if (restart_state.world_size != runtime_topology.world_size) {
    report.world_size_match = false;
    report.mismatch_messages.push_back(
        "world_size mismatch: restart=" + std::to_string(restart_state.world_size) +
        ", runtime=" + std::to_string(runtime_topology.world_size));
  }
  if (restart_state.pm_grid_nx != runtime_topology.pm_slab.global_nx ||
      restart_state.pm_grid_ny != runtime_topology.pm_slab.global_ny ||
      restart_state.pm_grid_nz != runtime_topology.pm_slab.global_nz) {
    report.pm_grid_shape_match = false;
    report.mismatch_messages.push_back(
        "PM grid mismatch: restart=(" + std::to_string(restart_state.pm_grid_nx) + "," +
        std::to_string(restart_state.pm_grid_ny) + "," + std::to_string(restart_state.pm_grid_nz) +
        "), runtime=(" + std::to_string(runtime_topology.pm_slab.global_nx) + "," +
        std::to_string(runtime_topology.pm_slab.global_ny) + "," +
        std::to_string(runtime_topology.pm_slab.global_nz) + ")");
  }
  if (restart_state.pm_decomposition_mode != runtime_topology.pm_decomposition_mode) {
    report.pm_decomposition_mode_match = false;
    report.mismatch_messages.push_back(
        "PM decomposition mode mismatch: restart=" + restart_state.pm_decomposition_mode +
        ", runtime=" + runtime_topology.pm_decomposition_mode);
  }
  if (restart_state.pm_update_cadence_steps == 0) {
    report.pm_cadence_steps_match = false;
    report.mismatch_messages.push_back(
        "PM cadence mismatch: restart cadence must be >= 1, got 0");
  }
  if (restart_state.last_long_range_refresh_opportunity > restart_state.gravity_kick_opportunity) {
    report.gravity_kick_state_match = false;
    report.mismatch_messages.push_back(
        "gravity kick mismatch: last refresh opportunity exceeds current kick opportunity");
  }
  if ((restart_state.long_range_field_version == 0) !=
      (restart_state.last_long_range_refresh_opportunity == 0)) {
    report.long_range_field_state_match = false;
    report.mismatch_messages.push_back(
        "long-range field mismatch: field version and refresh opportunity are inconsistent");
  }
  if (runtime_topology.world_rank < 0 ||
      runtime_topology.world_rank >= static_cast<int>(restart_state.pm_slab_begin_x_by_rank.size())) {
    report.pm_local_slab_match = false;
    report.mismatch_messages.push_back("runtime world_rank is outside restart PM slab ownership table");
    return report;
  }
  const std::size_t rank = static_cast<std::size_t>(runtime_topology.world_rank);
  const std::size_t begin_restart = restart_state.pm_slab_begin_x_by_rank[rank];
  const std::size_t end_restart = restart_state.pm_slab_end_x_by_rank[rank];
  if (begin_restart != runtime_topology.pm_slab.owned_x.begin_x ||
      end_restart != runtime_topology.pm_slab.owned_x.end_x) {
    report.pm_local_slab_match = false;
    report.mismatch_messages.push_back(
        "PM slab mismatch for rank " + std::to_string(runtime_topology.world_rank) +
        ": restart=[" + std::to_string(begin_restart) + "," + std::to_string(end_restart) +
        "), runtime=[" + std::to_string(runtime_topology.pm_slab.owned_x.begin_x) + "," +
        std::to_string(runtime_topology.pm_slab.owned_x.end_x) + ")");
  }
  return report;
}

void recordDistributedProfiling(
    core::ProfilerSession* profiler,
    const LoadBalanceMetrics& metrics,
    std::uint64_t ghost_exchange_send_bytes,
    std::uint64_t ghost_exchange_recv_bytes) {
  if (profiler == nullptr) {
    return;
  }

  profiler->counters().setCount("parallel.ghost_exchange_send_bytes", ghost_exchange_send_bytes);
  profiler->counters().setCount("parallel.ghost_exchange_recv_bytes", ghost_exchange_recv_bytes);
  profiler->counters().setCount("parallel.total_memory_bytes", metrics.total_memory_bytes);

  const std::uint64_t imbalance_ppm = (metrics.weighted_imbalance_ratio <= 0.0)
                                          ? 0ULL
                                          : static_cast<std::uint64_t>(std::llround(metrics.weighted_imbalance_ratio * 1.0e6));
  profiler->counters().setCount("parallel.weighted_imbalance_ratio_ppm", imbalance_ppm);

  const std::uint64_t bytes_moved = ghost_exchange_send_bytes + ghost_exchange_recv_bytes;
  profiler->addBytesMoved(bytes_moved);
}

MpiContext::MpiContext() {
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  int initialized = 0;
  MPI_Initialized(&initialized);
  if (initialized != 0) {
    m_is_enabled = true;
    MPI_Comm_size(MPI_COMM_WORLD, &m_world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &m_world_rank);
  }
#endif
}

MpiContext::MpiContext(bool is_enabled, int world_size, int world_rank)
    : m_is_enabled(is_enabled), m_world_size(world_size), m_world_rank(world_rank) {
  if (world_size <= 0) {
    throw std::invalid_argument("MpiContext world_size must be positive");
  }
  if (world_rank < 0 || world_rank >= world_size) {
    throw std::invalid_argument("MpiContext world_rank must be within [0, world_size)");
  }
}

bool MpiContext::isEnabled() const noexcept { return m_is_enabled; }

bool MpiContext::isRoot() const noexcept { return m_world_rank == 0; }

int MpiContext::worldSize() const noexcept { return m_world_size; }

int MpiContext::worldRank() const noexcept { return m_world_rank; }

void MpiContext::validateExpectedWorldSizeOrThrow(int expected_world_size) const {
  if (expected_world_size <= 0) {
    throw std::invalid_argument("expected_world_size must be positive");
  }
  if (expected_world_size != m_world_size) {
    throw std::runtime_error(
        "parallel.mpi_ranks_expected does not match runtime world size: expected=" +
        std::to_string(expected_world_size) + ", runtime=" + std::to_string(m_world_size));
  }
}

double MpiContext::allreduceSumDouble(double local_value) const {
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  if (m_is_enabled) {
    double global = 0.0;
    MPI_Allreduce(&local_value, &global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return global;
  }
#endif
  return local_value;
}

std::uint64_t MpiContext::allreduceSumUint64(std::uint64_t local_value) const {
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  if (m_is_enabled) {
    std::uint64_t global = 0;
    MPI_Allreduce(&local_value, &global, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
    return global;
  }
#endif
  return local_value;
}

std::uint64_t MpiContext::allreduceXorUint64(std::uint64_t local_value) const {
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  if (m_is_enabled) {
    std::uint64_t global = 0;
    MPI_Allreduce(&local_value, &global, 1, MPI_UINT64_T, MPI_BXOR, MPI_COMM_WORLD);
    return global;
  }
#endif
  return local_value;
}



bool RankDeviceAssignment::isValid() const noexcept {
  if (requested_device_count < 0 || visible_device_count < 0 || active_device_count < 0) {
    return false;
  }
  if (!uses_cuda) {
    return assigned_device_index == -1;
  }
  return active_device_count > 0 && assigned_device_index >= 0 && assigned_device_index < active_device_count &&
      visible_device_count >= active_device_count;
}

RankDeviceAssignment selectRankDeviceAssignment(
    int world_rank,
    int configured_gpu_devices,
    bool cuda_runtime_available,
    int visible_device_count) {
  if (world_rank < 0) {
    throw std::invalid_argument("world_rank must be non-negative");
  }
  if (configured_gpu_devices < 0) {
    throw std::invalid_argument("configured_gpu_devices must be >= 0");
  }
  if (visible_device_count < 0) {
    throw std::invalid_argument("visible_device_count must be >= 0");
  }

  RankDeviceAssignment assignment;
  assignment.requested_device_count = configured_gpu_devices;
  assignment.visible_device_count = visible_device_count;

  if (configured_gpu_devices == 0) {
    return assignment;
  }
  if (!cuda_runtime_available || visible_device_count == 0) {
    throw std::runtime_error(
        "parallel.gpu_devices requested CUDA PM execution, but no CUDA runtime devices are available");
  }
  if (configured_gpu_devices > visible_device_count) {
    throw std::runtime_error(
        "parallel.gpu_devices exceeds visible CUDA devices: requested=" + std::to_string(configured_gpu_devices) +
        ", visible=" + std::to_string(visible_device_count));
  }

  assignment.uses_cuda = true;
  assignment.active_device_count = configured_gpu_devices;
  assignment.assigned_device_index = world_rank % configured_gpu_devices;
  return assignment;
}

DistributedExecutionTopology buildDistributedExecutionTopology(
    std::size_t global_nx,
    std::size_t global_ny,
    std::size_t global_nz,
    const MpiContext& mpi_context,
    int mpi_ranks_expected,
    int configured_gpu_devices,
    bool cuda_runtime_available,
    int visible_device_count,
    std::string pm_decomposition_mode) {
  mpi_context.validateExpectedWorldSizeOrThrow(mpi_ranks_expected);

  DistributedExecutionTopology topology;
  topology.world_size = mpi_context.worldSize();
  topology.world_rank = mpi_context.worldRank();
  topology.mpi_enabled = mpi_context.isEnabled();
  topology.pm_decomposition_mode = std::move(pm_decomposition_mode);
  topology.pm_slab = makePmSlabLayout(global_nx, global_ny, global_nz, mpi_context.worldSize(), mpi_context.worldRank());
  topology.device_assignment =
      selectRankDeviceAssignment(mpi_context.worldRank(), configured_gpu_devices, cuda_runtime_available, visible_device_count);
  if (!topology.device_assignment.isValid()) {
    throw std::runtime_error("constructed distributed execution topology is invalid");
  }
  return topology;
}


}  // namespace cosmosim::parallel
