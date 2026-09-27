#include "cosmosim/parallel/distributed_memory.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <iomanip>
#include <limits>
#include <numeric>
#include <optional>
#include <sstream>
#include <streambuf>
#include <unordered_map>
#include <unordered_set>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/core/memory_governor.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "parallel/internal/distributed_memory_support.hpp"

#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
#include <mpi.h>
#endif

namespace cosmosim::parallel {
namespace {

using internal::injectMpiTestFault;

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
BoundedMpiTransferPlan planBoundedMpiTransferRounds(
    std::span<const std::size_t> logical_counts,
    std::size_t mpi_count_limit,
    std::size_t round_count_limit) {
  if (mpi_count_limit == 0U ||
      mpi_count_limit > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    throw std::invalid_argument(
        "bounded MPI planner count limit must be within classic MPI int range");
  }
  if (round_count_limit == 0U || round_count_limit > mpi_count_limit) {
    throw std::invalid_argument(
        "bounded MPI planner round limit must be positive and no larger than the MPI count limit");
  }

  BoundedMpiTransferPlan plan;
  plan.logical_counts.assign(logical_counts.begin(), logical_counts.end());
  plan.logical_displacements.resize(logical_counts.size(), 0U);
  if (logical_counts.empty()) {
    return plan;
  }
  if (logical_counts.size() > round_count_limit) {
    throw std::invalid_argument(
        "bounded MPI planner round limit is too small to address every participant safely");
  }

  std::size_t logical_total = 0U;
  for (std::size_t peer = 0; peer < logical_counts.size(); ++peer) {
    plan.logical_displacements[peer] = logical_total;
    logical_total = core::checkedSizeAdd(
        logical_total, logical_counts[peer], "bounded MPI logical prefix");
  }
  plan.logical_total_count = logical_total;

  const std::size_t fair_round_share = round_count_limit / logical_counts.size();
  plan.per_peer_count_limit = std::min(mpi_count_limit, fair_round_share);
  if (plan.per_peer_count_limit == 0U) {
    throw std::invalid_argument("bounded MPI planner computed a zero per-peer round limit");
  }

  std::size_t round_count = 0U;
  for (const std::size_t count : logical_counts) {
    const std::size_t peer_rounds = count / plan.per_peer_count_limit +
        (count % plan.per_peer_count_limit == 0U ? 0U : 1U);
    round_count = std::max(round_count, peer_rounds);
  }
  plan.rounds.reserve(round_count);
  std::vector<std::size_t> consumed(logical_counts.size(), 0U);
  for (std::size_t round_index = 0; round_index < round_count; ++round_index) {
    BoundedMpiRoundLayout round;
    round.counts.resize(logical_counts.size(), 0);
    round.displacements.resize(logical_counts.size(), 0);
    round.logical_offsets = consumed;

    std::size_t round_total = 0U;
    for (std::size_t peer = 0; peer < logical_counts.size(); ++peer) {
      const std::size_t remaining = logical_counts[peer] - consumed[peer];
      const std::size_t count = std::min(plan.per_peer_count_limit, remaining);
      round.displacements[peer] = core::checkedIntegralNarrow<int>(
          round_total, "bounded MPI round displacement");
      round.counts[peer] = core::checkedIntegralNarrow<int>(
          count, "bounded MPI round count");
      round_total = core::checkedSizeAdd(
          round_total, count, "bounded MPI round aggregate");
      if (round_total > round_count_limit || round_total > mpi_count_limit) {
        throw std::overflow_error("bounded MPI round aggregate exceeds representability limit");
      }
      consumed[peer] = core::checkedSizeAdd(
          consumed[peer], count, "bounded MPI logical coverage");
    }
    round.round_count = round_total;
    plan.rounds.push_back(std::move(round));
  }
  if (consumed != plan.logical_counts) {
    throw std::logic_error("bounded MPI planner did not cover the complete logical payload");
  }
  return plan;
}

std::size_t mpiTransportRoundLimitBytes() {
#if COSMOSIM_ENABLE_TESTS
  const char* raw = std::getenv("COSMOSIM_MPI_TEST_TRANSPORT_LIMIT_BYTES");
  if (raw != nullptr && *raw != '\0') {
    try {
      const unsigned long long parsed = std::stoull(raw);
      if (parsed > 0ULL && parsed <= static_cast<unsigned long long>(std::numeric_limits<int>::max())) {
        return static_cast<std::size_t>(parsed);
      }
    } catch (...) {
      // Invalid test-only overrides fall back to the production-safe default.
    }
  }
#endif
  return k_default_mpi_transport_round_bytes;
}

std::vector<std::vector<std::uint8_t>> exchangeBoundedAlltoallBytes(
    const MpiContext& mpi_context,
    const std::vector<std::vector<std::uint8_t>>& send_payloads) {
  const int world_size = mpi_context.worldSize();
  if (world_size <= 0) {
    throw std::invalid_argument("bounded byte all-to-all requires positive world size");
  }
  if (!mpi_context.isEnabled()) {
    if (world_size != 1 || send_payloads.size() != 1U) {
      throw std::runtime_error(
          "bounded byte all-to-all requires MPI when world size exceeds one");
    }
    return {send_payloads.front()};
  }
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  std::vector<std::uint64_t> send_counts64;
  std::vector<std::uint64_t> recv_counts64;
  std::exception_ptr local_failure;
  try {
    if (send_payloads.size() != static_cast<std::size_t>(world_size)) {
      throw std::invalid_argument(
          "bounded byte all-to-all payload rank extent must match MPI world size");
    }
    send_counts64.resize(static_cast<std::size_t>(world_size), 0U);
    recv_counts64.resize(static_cast<std::size_t>(world_size), 0U);
    for (std::size_t rank = 0; rank < send_payloads.size(); ++rank) {
      send_counts64[rank] = core::checkedIntegralNarrow<std::uint64_t>(
          send_payloads[rank].size(), "bounded byte all-to-all send count");
    }
  } catch (...) {
    local_failure = std::current_exception();
  }
  mpi_context.rethrowCollectivePreparationFailure(
      local_failure, "bounded byte all-to-all pre-count preparation");

  if (MPI_Alltoall(
          send_counts64.data(), 1, MPI_UINT64_T,
          recv_counts64.data(), 1, MPI_UINT64_T,
          MPI_COMM_WORLD) != MPI_SUCCESS) {
    throw std::runtime_error("bounded byte all-to-all count exchange failed");
  }

  BoundedMpiTransferPlan send_plan;
  BoundedMpiTransferPlan recv_plan;
  std::vector<std::vector<std::uint8_t>> recv_payloads;
  std::vector<std::uint8_t> send_round_buffer;
  std::vector<std::uint8_t> recv_round_buffer;
  BoundedMpiRoundLayout zero_round;
  local_failure = nullptr;
  try {
    std::vector<std::size_t> send_counts(send_counts64.size(), 0U);
    std::vector<std::size_t> recv_counts(recv_counts64.size(), 0U);
    for (std::size_t rank = 0; rank < send_counts.size(); ++rank) {
      send_counts[rank] = core::checkedIntegralNarrow<std::size_t>(
          send_counts64[rank], "bounded byte all-to-all logical send count");
      recv_counts[rank] = core::checkedIntegralNarrow<std::size_t>(
          recv_counts64[rank], "bounded byte all-to-all logical receive count");
    }
    const std::size_t round_limit = mpiTransportRoundLimitBytes();
    send_plan = planBoundedMpiTransferRounds(
        send_counts,
        static_cast<std::size_t>(std::numeric_limits<int>::max()),
        round_limit);
    recv_plan = planBoundedMpiTransferRounds(
        recv_counts,
        static_cast<std::size_t>(std::numeric_limits<int>::max()),
        round_limit);

    recv_payloads.resize(static_cast<std::size_t>(world_size));
    for (std::size_t rank = 0; rank < recv_payloads.size(); ++rank) {
      recv_payloads[rank].resize(recv_counts[rank]);
    }
    std::size_t maximum_send_round = 0U;
    for (const auto& round : send_plan.rounds) {
      maximum_send_round = std::max(maximum_send_round, round.round_count);
    }
    std::size_t maximum_recv_round = 0U;
    for (const auto& round : recv_plan.rounds) {
      maximum_recv_round = std::max(maximum_recv_round, round.round_count);
    }
    send_round_buffer.resize(maximum_send_round);
    recv_round_buffer.resize(maximum_recv_round);
    zero_round.counts.resize(static_cast<std::size_t>(world_size), 0);
    zero_round.displacements.resize(static_cast<std::size_t>(world_size), 0);
    zero_round.logical_offsets.resize(static_cast<std::size_t>(world_size), 0U);
    injectMpiTestFault(mpi_context, "alltoall_post_count");
  } catch (...) {
    local_failure = std::current_exception();
  }
  mpi_context.rethrowCollectivePreparationFailure(
      local_failure, "bounded byte all-to-all payload preparation");

  const std::uint64_t local_round_count = static_cast<std::uint64_t>(
      std::max(send_plan.rounds.size(), recv_plan.rounds.size()));
  const std::uint64_t global_round_count =
      mpi_context.allreduceMaxUint64(local_round_count);
  for (std::uint64_t round_index = 0U; round_index < global_round_count;
       ++round_index) {
    const BoundedMpiRoundLayout& send_round =
        round_index < send_plan.rounds.size()
            ? send_plan.rounds[static_cast<std::size_t>(round_index)]
            : zero_round;
    const BoundedMpiRoundLayout& recv_round =
        round_index < recv_plan.rounds.size()
            ? recv_plan.rounds[static_cast<std::size_t>(round_index)]
            : zero_round;

    for (std::size_t peer = 0; peer < send_round.counts.size(); ++peer) {
      const std::size_t count = static_cast<std::size_t>(send_round.counts[peer]);
      if (count == 0U) {
        continue;
      }
      std::memcpy(
          send_round_buffer.data() +
              static_cast<std::size_t>(send_round.displacements[peer]),
          send_payloads[peer].data() + send_round.logical_offsets[peer],
          count);
    }
    if (MPI_Alltoallv(
            send_round_buffer.empty() ? nullptr : send_round_buffer.data(),
            send_round.counts.data(), send_round.displacements.data(), MPI_BYTE,
            recv_round_buffer.empty() ? nullptr : recv_round_buffer.data(),
            recv_round.counts.data(), recv_round.displacements.data(), MPI_BYTE,
            MPI_COMM_WORLD) != MPI_SUCCESS) {
      throw std::runtime_error("bounded byte all-to-all payload exchange failed");
    }
    for (std::size_t peer = 0; peer < recv_round.counts.size(); ++peer) {
      const std::size_t count = static_cast<std::size_t>(recv_round.counts[peer]);
      if (count == 0U) {
        continue;
      }
      std::memcpy(
          recv_payloads[peer].data() + recv_round.logical_offsets[peer],
          recv_round_buffer.data() +
              static_cast<std::size_t>(recv_round.displacements[peer]),
          count);
    }
  }
  return recv_payloads;
#else
  throw std::runtime_error(
      "bounded byte all-to-all requires an MPI-enabled build");
#endif
}

std::vector<DecompositionItem> gatherDecompositionItemsAcrossRanks(
    const MpiContext& mpi_context,
    std::span<const DecompositionItem> local_items) {
  static_assert(std::is_trivially_copyable_v<DecompositionItem>);
  if (mpi_context.worldSize() == 1) {
    return std::vector<DecompositionItem>(local_items.begin(), local_items.end());
  }
  if (!mpi_context.isEnabled()) {
    throw std::runtime_error("multi-rank decomposition item gather requires MPI to be enabled");
  }
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  std::size_t local_byte_count = 0U;
  std::exception_ptr local_preparation_failure;
  try {
    local_byte_count = core::checkedSizeMultiply(
        local_items.size(), sizeof(DecompositionItem),
        "decomposition item gather local byte count");
  } catch (...) {
    local_preparation_failure = std::current_exception();
  }
  mpi_context.rethrowCollectivePreparationFailure(
      local_preparation_failure,
      "decomposition item gather local byte framing");

  const auto local_bytes = std::span<const std::uint8_t>(
      reinterpret_cast<const std::uint8_t*>(local_items.data()),
      local_byte_count);
  const std::vector<std::uint8_t> recv_bytes =
      mpi_context.allgatherBytesBounded(local_bytes);

  std::vector<DecompositionItem> gathered;
  std::exception_ptr local_reassembly_failure;
  try {
    if (recv_bytes.size() % sizeof(DecompositionItem) != 0U) {
      throw std::runtime_error(
          "decomposition item gather returned a non-record-aligned byte count");
    }
    gathered.resize(recv_bytes.size() / sizeof(DecompositionItem));
    if (!recv_bytes.empty()) {
      std::memcpy(gathered.data(), recv_bytes.data(), recv_bytes.size());
    }
  } catch (...) {
    local_reassembly_failure = std::current_exception();
  }
  mpi_context.rethrowCollectivePreparationFailure(
      local_reassembly_failure,
      "decomposition item gather receive reassembly");
  return gathered;
#else
  throw std::runtime_error("multi-rank decomposition item gather requires an MPI build");
#endif
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

  bool has_epoch = false;
  GhostLayerEpoch common_epoch{};
  for (const LocalGhostDescriptor descriptor : local_ghost_descriptors) {
    if (!has_epoch) {
      common_epoch = descriptor.epoch;
      has_epoch = true;
    } else if (!descriptor.epoch.matches(common_epoch)) {
      throw std::invalid_argument("ghost descriptors in one exchange plan must share a common epoch");
    }
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

  plan.epoch = common_epoch;
  plan.exchange_sequence = common_epoch.ghost_sync_epoch;
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

GhostExchangePlan buildExplicitGhostExchangePlan(
    int world_rank,
    std::span<const int> neighbor_ranks,
    std::span<const std::vector<std::uint32_t>> send_local_indices_by_neighbor,
    std::span<const std::vector<std::uint32_t>> recv_local_indices_by_neighbor,
    std::size_t bytes_per_ghost,
    const GhostLayerEpoch& epoch,
    bool enable_nonblocking_overlap) {
  if (world_rank < 0) {
    throw std::invalid_argument("world_rank must be non-negative");
  }
  if (bytes_per_ghost == 0) {
    throw std::invalid_argument("bytes_per_ghost must be positive");
  }
  if (neighbor_ranks.size() != send_local_indices_by_neighbor.size() ||
      neighbor_ranks.size() != recv_local_indices_by_neighbor.size()) {
    throw std::invalid_argument("explicit ghost exchange plan container sizes must match");
  }

  GhostExchangePlan plan;
  plan.neighbor_ranks.assign(neighbor_ranks.begin(), neighbor_ranks.end());
  plan.send_local_indices_by_neighbor.assign(
      send_local_indices_by_neighbor.begin(), send_local_indices_by_neighbor.end());
  plan.recv_local_indices_by_neighbor.assign(
      recv_local_indices_by_neighbor.begin(), recv_local_indices_by_neighbor.end());
  plan.outbound_transfers.resize(neighbor_ranks.size());
  plan.inbound_transfers.resize(neighbor_ranks.size());
  plan.epoch = epoch;
  plan.exchange_sequence = epoch.ghost_sync_epoch;
  plan.uses_blocking_exchange = true;
  plan.nonblocking_overlap_enabled = enable_nonblocking_overlap;

  for (std::size_t i = 0; i < neighbor_ranks.size(); ++i) {
    if (neighbor_ranks[i] < 0 || neighbor_ranks[i] == world_rank) {
      throw std::invalid_argument("explicit ghost exchange neighbor rank must be a remote non-negative rank");
    }
    plan.outbound_transfers[i] = GhostTransferDescriptor{
        .role = GhostTransferRole::kOutboundSend,
        .intent = GhostTransferIntent::kGhostRefreshRequest,
        .peer_rank = neighbor_ranks[i],
        .neighbor_slot = i,
        .expected_post_transfer_residency = LocalIndexResidency::kGhost,
        .local_indices = plan.send_local_indices_by_neighbor[i],
    };
    plan.inbound_transfers[i] = GhostTransferDescriptor{
        .role = GhostTransferRole::kInboundReceive,
        .intent = GhostTransferIntent::kGhostRefreshReceiveStaging,
        .peer_rank = neighbor_ranks[i],
        .neighbor_slot = i,
        .expected_post_transfer_residency = LocalIndexResidency::kGhost,
        .local_indices = plan.recv_local_indices_by_neighbor[i],
    };
    plan.send_bytes += static_cast<std::uint64_t>(plan.send_local_indices_by_neighbor[i].size()) * bytes_per_ghost;
    plan.recv_bytes += static_cast<std::uint64_t>(plan.recv_local_indices_by_neighbor[i].size()) * bytes_per_ghost;
  }

  validateGhostExchangePlan(plan);
  return plan;
}

void validateGhostExchangePlan(const GhostExchangePlan& plan) {
  if (!plan.uses_blocking_exchange && !plan.nonblocking_overlap_enabled) {
    throw std::invalid_argument("ghost exchange plan must expose either the default blocking path or an explicit overlap path");
  }
  if (plan.nonblocking_overlap_enabled && !plan.uses_blocking_exchange) {
    throw std::invalid_argument("nonblocking ghost exchange overlap must share the blocking ownership contract");
  }
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

void validateGhostTransferAgainstResidency(
    const GhostTransferDescriptor& descriptor,
    std::span<const LocalGhostDescriptor> local_ghost_descriptors,
    int world_rank) {
  if (world_rank < 0) {
    throw std::invalid_argument("world_rank must be non-negative");
  }
  for (const std::uint32_t local_index : descriptor.local_indices) {
    if (local_index >= local_ghost_descriptors.size()) {
      throw std::out_of_range("ghost transfer descriptor local index out of residency table range");
    }
    const LocalGhostDescriptor local = local_ghost_descriptors[local_index];
    if (descriptor.role == GhostTransferRole::kOutboundSend) {
      if (local.residency != LocalIndexResidency::kOwned || local.owning_rank != world_rank) {
        throw std::invalid_argument("outbound ghost or migration payload must be packed from authoritative local state");
      }
    } else {
      if (descriptor.intent == GhostTransferIntent::kGhostRefreshReceiveStaging) {
        if (local.residency != LocalIndexResidency::kGhost || local.owning_rank == world_rank) {
          throw std::invalid_argument("ghost refresh receive staging must unpack into remote-owned ghost slots");
        }
      } else if (descriptor.intent == GhostTransferIntent::kOwnershipMigrationReceiveStaging) {
        if (descriptor.expected_post_transfer_residency != LocalIndexResidency::kOwned) {
          throw std::invalid_argument("ownership migration receive staging must produce owned local state");
        }
      }
    }
  }
}

void validateBlockingGhostExchangeContracts(
    const GhostExchangePlan& plan,
    std::span<const LocalGhostDescriptor> local_ghost_descriptors,
    int world_rank,
    const GhostLayerEpoch& expected_epoch) {
  validateGhostExchangePlan(plan);
  if (!plan.uses_blocking_exchange) {
    throw std::invalid_argument("default ghost exchange path must be blocking and correctness-first");
  }
  if (!plan.epoch.matches(expected_epoch)) {
    throw std::invalid_argument("ghost exchange plan epoch is stale for the current decomposition/sync generation");
  }
  for (const GhostTransferDescriptor& descriptor : plan.outbound_transfers) {
    validateGhostTransferAgainstResidency(descriptor, local_ghost_descriptors, world_rank);
  }
  for (const GhostTransferDescriptor& descriptor : plan.inbound_transfers) {
    validateGhostTransferAgainstResidency(descriptor, local_ghost_descriptors, world_rank);
  }
}



}  // namespace cosmosim::parallel
