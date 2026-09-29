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

using internal::commitOptionalGhostLane;
using internal::ghostExchangeRecordBytes;
using internal::injectMpiTestFault;
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
using internal::queryActiveMpiWorld;
#endif

}  // namespace
namespace {

template <typename Record>
[[nodiscard]] std::vector<Record> allgatherTrivialRecordsBounded(
    const MpiContext& mpi_context,
    std::span<const Record> local_records,
    std::string_view phase) {
  static_assert(std::is_trivially_copyable_v<Record>);

  std::size_t local_byte_count = 0U;
  std::exception_ptr local_preparation_failure;
  try {
    local_byte_count = core::checkedSizeMultiply(
        local_records.size(), sizeof(Record),
        std::string(phase) + " local byte count");
  } catch (...) {
    local_preparation_failure = std::current_exception();
  }
  mpi_context.rethrowCollectivePreparationFailure(
      local_preparation_failure,
      std::string(phase) + " local wire preparation");

  const auto local_wire = std::span<const std::uint8_t>(
      reinterpret_cast<const std::uint8_t*>(local_records.data()),
      local_byte_count);
  std::vector<std::uint8_t> gathered_wire =
      mpi_context.allgatherBytesBounded(local_wire);

  std::vector<Record> result;
  std::exception_ptr local_reassembly_failure;
  try {
    if (gathered_wire.size() % sizeof(Record) != 0U) {
      throw std::runtime_error(
          std::string(phase) + " returned partial record bytes");
    }
    result.resize(gathered_wire.size() / sizeof(Record));
    if (!gathered_wire.empty()) {
      std::memcpy(result.data(), gathered_wire.data(), gathered_wire.size());
    }
  } catch (...) {
    local_reassembly_failure = std::current_exception();
  }
  mpi_context.rethrowCollectivePreparationFailure(
      local_reassembly_failure,
      std::string(phase) + " receive reassembly");
  return result;
}

}  // namespace

std::vector<HydroConservativeFluxCorrectionRecord> executeBlockingHydroConservativeFluxCorrectionExchange(
    const MpiContext& mpi_context,
    std::span<const HydroConservativeFluxCorrectionRecord> local_records,
    std::uint64_t exchange_sequence) {
  (void)exchange_sequence;
  std::exception_ptr local_validation_failure;
  try {
    for (const HydroConservativeFluxCorrectionRecord& record : local_records) {
      validateHydroConservativeFluxCorrectionRecord(record);
      if (record.source_rank != mpi_context.worldRank()) {
        throw std::invalid_argument(
            "hydro conservative flux correction source rank does not match MPI context");
      }
    }
  } catch (...) {
    local_validation_failure = std::current_exception();
  }
  if (mpi_context.isEnabled()) {
    mpi_context.rethrowCollectivePreparationFailure(
        local_validation_failure,
        "hydro conservative flux correction local validation");
  } else if (local_validation_failure != nullptr) {
    std::rethrow_exception(local_validation_failure);
  }
  if (!mpi_context.isEnabled()) {
    return std::vector<HydroConservativeFluxCorrectionRecord>(
        local_records.begin(), local_records.end());
  }
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  const int world_size = mpi_context.worldSize();
  std::vector<HydroConservativeFluxCorrectionRecord> result =
      allgatherTrivialRecordsBounded(
          mpi_context, local_records,
          "hydro conservative flux correction exchange");

  std::exception_ptr local_result_failure;
  try {
    for (const HydroConservativeFluxCorrectionRecord& record : result) {
      validateHydroConservativeFluxCorrectionRecord(record);
      if (record.source_rank < 0 || record.source_rank >= world_size ||
          record.owner_rank < 0 || record.owner_rank >= world_size) {
        throw std::runtime_error(
            "hydro conservative flux correction exchange returned invalid rank metadata");
      }
    }
  } catch (...) {
    local_result_failure = std::current_exception();
  }
  mpi_context.rethrowCollectivePreparationFailure(
      local_result_failure,
      "hydro conservative flux correction receive validation");
  return result;
#else
  throw std::runtime_error(
      "hydro conservative flux correction exchange requires MPI support when MPI context is enabled");
#endif
}

std::vector<HydroGhostCellRequest> executeBlockingHydroGhostCellRequestExchange(
    const MpiContext& mpi_context,
    std::span<const HydroGhostCellRequest> local_requests,
    std::uint64_t exchange_sequence) {
  (void)exchange_sequence;
  std::exception_ptr local_validation_failure;
  try {
    for (const HydroGhostCellRequest& request : local_requests) {
      validateHydroGhostCellRequest(request);
      if (request.descriptor.consumer_rank != mpi_context.worldRank()) {
        throw std::invalid_argument(
            "hydro ghost cell request consumer rank does not match MPI context");
      }
    }
  } catch (...) {
    local_validation_failure = std::current_exception();
  }
  if (mpi_context.isEnabled()) {
    mpi_context.rethrowCollectivePreparationFailure(
        local_validation_failure,
        "hydro ghost cell request local validation");
  } else if (local_validation_failure != nullptr) {
    std::rethrow_exception(local_validation_failure);
  }
  if (!mpi_context.isEnabled()) {
    return std::vector<HydroGhostCellRequest>(
        local_requests.begin(), local_requests.end());
  }
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  const int world_size = mpi_context.worldSize();
  std::vector<HydroGhostCellRequest> result =
      allgatherTrivialRecordsBounded(
          mpi_context, local_requests,
          "hydro ghost cell request exchange");

  std::exception_ptr local_result_failure;
  try {
    for (const HydroGhostCellRequest& request : result) {
      validateHydroGhostCellRequest(request);
      if (request.descriptor.owner_rank < 0 ||
          request.descriptor.owner_rank >= world_size ||
          request.descriptor.consumer_rank < 0 ||
          request.descriptor.consumer_rank >= world_size) {
        throw std::runtime_error(
            "hydro ghost cell request exchange returned invalid rank metadata");
      }
    }
  } catch (...) {
    local_result_failure = std::current_exception();
  }
  mpi_context.rethrowCollectivePreparationFailure(
      local_result_failure,
      "hydro ghost cell request receive validation");
  return result;
#else
  throw std::runtime_error(
      "hydro ghost cell request exchange requires MPI support when MPI context is enabled");
#endif
}

std::vector<HydroGhostCellPayloadRecord> executeBlockingHydroGhostCellPayloadExchange(
    const MpiContext& mpi_context,
    std::span<const HydroGhostCellPayloadRecord> local_records,
    std::uint64_t exchange_sequence) {
  (void)exchange_sequence;
  std::exception_ptr local_validation_failure;
  try {
    for (const HydroGhostCellPayloadRecord& record : local_records) {
      validateHydroGhostCellPayloadRecord(record);
      if (record.descriptor.owner_rank != mpi_context.worldRank()) {
        throw std::invalid_argument(
            "hydro ghost cell payload owner rank does not match MPI context");
      }
    }
  } catch (...) {
    local_validation_failure = std::current_exception();
  }
  if (mpi_context.isEnabled()) {
    mpi_context.rethrowCollectivePreparationFailure(
        local_validation_failure,
        "hydro ghost cell payload local validation");
  } else if (local_validation_failure != nullptr) {
    std::rethrow_exception(local_validation_failure);
  }
  if (!mpi_context.isEnabled()) {
    return std::vector<HydroGhostCellPayloadRecord>(
        local_records.begin(), local_records.end());
  }
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  const int world_size = mpi_context.worldSize();
  std::vector<HydroGhostCellPayloadRecord> result =
      allgatherTrivialRecordsBounded(
          mpi_context, local_records,
          "hydro ghost cell payload exchange");

  std::exception_ptr local_result_failure;
  try {
    for (const HydroGhostCellPayloadRecord& record : result) {
      validateHydroGhostCellPayloadRecord(record);
      if (record.descriptor.owner_rank < 0 ||
          record.descriptor.owner_rank >= world_size ||
          record.descriptor.consumer_rank < 0 ||
          record.descriptor.consumer_rank >= world_size) {
        throw std::runtime_error(
            "hydro ghost cell payload exchange returned invalid rank metadata");
      }
    }
  } catch (...) {
    local_result_failure = std::current_exception();
  }
  mpi_context.rethrowCollectivePreparationFailure(
      local_result_failure,
      "hydro ghost cell payload receive validation");
  return result;
#else
  throw std::runtime_error(
      "hydro ghost cell payload exchange requires MPI support when MPI context is enabled");
#endif
}

PmSlabHaloExchangeResult executeBlockingPmSlabHaloExchange(
    const MpiContext& mpi_context,
    const PmSlabLayout& layout,
    std::span<const double> local_scalar_field,
    std::size_t halo_depth_x,
    bool periodic_x,
    std::uint64_t exchange_sequence) {
#if !defined(COSMOSIM_ENABLE_MPI) || !COSMOSIM_ENABLE_MPI
  (void)exchange_sequence;
#endif
  PmSlabHaloExchangeResult result;
  std::vector<double> send_left;
  std::vector<double> send_right;
  std::size_t halo_value_count = 0U;
  std::uint64_t payload_bytes = 0U;
  int left_peer = -1;
  int right_peer = -1;
  bool no_exchange = false;
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  int communicator_world_size = 1;
  int communicator_world_rank = 0;
  const bool communicator_mpi_active =
      queryActiveMpiWorld(communicator_world_size, communicator_world_rank);
#endif

  std::exception_ptr local_preparation_failure;
  try {
    if (!layout.isValid()) {
      throw std::invalid_argument("PM slab halo exchange requires a valid slab layout");
    }
    if (layout.world_size != mpi_context.worldSize() ||
        layout.world_rank != mpi_context.worldRank()) {
      throw std::invalid_argument(
          "PM slab halo exchange layout world metadata must match MPI context");
    }
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
    if (!communicator_mpi_active &&
        (layout.world_size > 1 || mpi_context.isEnabled())) {
      throw std::invalid_argument(
          "PM slab halo exchange requires an active MPI_COMM_WORLD for an enabled or distributed context");
    }
    if (communicator_mpi_active &&
        (layout.world_size != communicator_world_size ||
         layout.world_rank != communicator_world_rank)) {
      throw std::invalid_argument(
          "PM slab halo exchange layout world metadata must match MPI_COMM_WORLD");
    }
#endif
    if (layout.global_ny >
        std::numeric_limits<std::size_t>::max() / layout.global_nz) {
      throw std::overflow_error("PM slab halo exchange plane size overflows size_t");
    }
    const std::size_t plane_size = layout.global_ny * layout.global_nz;
    if (layout.local_nx() >
        std::numeric_limits<std::size_t>::max() / plane_size) {
      throw std::overflow_error("PM slab halo exchange local field size overflows size_t");
    }
    const std::size_t expected_local_values = layout.local_nx() * plane_size;
    if (local_scalar_field.size() != expected_local_values) {
      throw std::invalid_argument(
          "PM slab halo exchange field size does not match local slab cell count");
    }
    if (halo_depth_x == 0 || layout.world_size == 1 || layout.local_nx() == 0) {
      no_exchange = true;
    } else {
      if (!mpi_context.isEnabled()) {
        throw std::runtime_error(
            "PM slab halo exchange requires MPI for distributed layouts");
      }
      // A single-neighbor payload cannot span a second slab. Use the smallest
      // non-empty slab extent so every communicating rank posts matching counts.
      std::size_t minimum_nonempty_slab_nx =
          std::numeric_limits<std::size_t>::max();
      for (int rank = 0; rank < layout.world_size; ++rank) {
        const PmSlabRange owned =
            pmOwnedXRangeForRank(layout.global_nx, layout.world_size, rank);
        if (owned.extentX() > 0U) {
          minimum_nonempty_slab_nx =
              std::min(minimum_nonempty_slab_nx, owned.extentX());
        }
      }
      if (minimum_nonempty_slab_nx ==
          std::numeric_limits<std::size_t>::max()) {
        throw std::logic_error("PM slab halo exchange layout has no non-empty owner");
      }
      const std::size_t depth =
          std::min(halo_depth_x, minimum_nonempty_slab_nx);
      if (depth > std::numeric_limits<std::size_t>::max() / plane_size) {
        throw std::overflow_error("PM slab halo exchange payload size overflows size_t");
      }
      halo_value_count = depth * plane_size;
      if (halo_value_count >
          static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw std::overflow_error(
            "PM slab halo exchange payload count exceeds MPI int limit");
      }
      if (halo_value_count >
          std::numeric_limits<std::uint64_t>::max() / sizeof(double)) {
        throw std::overflow_error(
            "PM slab halo exchange byte diagnostics overflow uint64_t");
      }
      payload_bytes =
          static_cast<std::uint64_t>(halo_value_count) * sizeof(double);
      result.halo_depth_x = depth;
      if (layout.owned_x.begin_x > 0) {
        left_peer = pmOwnerRankForGlobalX(
            layout.global_nx,
            layout.world_size,
            layout.owned_x.begin_x - 1U);
      } else if (periodic_x) {
        left_peer = pmOwnerRankForGlobalX(
            layout.global_nx,
            layout.world_size,
            layout.global_nx - 1U);
      }
      if (layout.owned_x.end_x < layout.global_nx) {
        right_peer = pmOwnerRankForGlobalX(
            layout.global_nx,
            layout.world_size,
            layout.owned_x.end_x);
      } else if (periodic_x) {
        right_peer = pmOwnerRankForGlobalX(
            layout.global_nx, layout.world_size, 0U);
      }
      result.left_peer_rank = left_peer;
      result.right_peer_rank = right_peer;
      result.left_halo.assign(halo_value_count, 0.0);
      result.right_halo.assign(halo_value_count, 0.0);
      send_left.assign(halo_value_count, 0.0);
      send_right.assign(halo_value_count, 0.0);
      const std::span<const double> left_source =
          local_scalar_field.first(halo_value_count);
      const std::span<const double> right_source =
          local_scalar_field.last(halo_value_count);
      std::copy(left_source.begin(), left_source.end(), send_left.begin());
      std::copy(right_source.begin(), right_source.end(), send_right.begin());

      const int local_rank = mpi_context.worldRank();
      const auto is_remote_peer = [&](int peer) {
        return peer >= 0 && peer != local_rank;
      };
      const std::uint64_t remote_side_count =
          static_cast<std::uint64_t>(is_remote_peer(left_peer)) +
          static_cast<std::uint64_t>(is_remote_peer(right_peer));
      if (remote_side_count > 0 &&
          payload_bytes >
              std::numeric_limits<std::uint64_t>::max() / remote_side_count) {
        throw std::overflow_error(
            "PM slab halo exchange aggregate byte diagnostics overflow uint64_t");
      }
      result.sent_bytes = payload_bytes * remote_side_count;
      result.received_bytes = payload_bytes * remote_side_count;
    }
  } catch (...) {
    local_preparation_failure = std::current_exception();
  }

#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  if (communicator_mpi_active && communicator_world_size > 1) {
    const std::uint64_t local_failure_vote =
        local_preparation_failure ? 1U : 0U;
    std::uint64_t failure_count = 0U;
    MPI_Allreduce(
        &local_failure_vote,
        &failure_count,
        1,
        MPI_UINT64_T,
        MPI_SUM,
        MPI_COMM_WORLD);
    if (failure_count != 0U) {
      if (local_preparation_failure) {
        std::rethrow_exception(local_preparation_failure);
      }
      throw std::runtime_error(
          "PM slab halo exchange peer rejected protocol preparation");
    }
    const std::array<std::uint64_t, 7> local_protocol_identity{
        static_cast<std::uint64_t>(layout.global_nx),
        static_cast<std::uint64_t>(layout.global_ny),
        static_cast<std::uint64_t>(layout.global_nz),
        static_cast<std::uint64_t>(halo_depth_x),
        periodic_x ? 1U : 0U,
        exchange_sequence,
        static_cast<std::uint64_t>(communicator_world_size),
    };
    std::array<std::uint64_t, 7> minimum_protocol_identity{};
    std::array<std::uint64_t, 7> maximum_protocol_identity{};
    MPI_Allreduce(
        local_protocol_identity.data(),
        minimum_protocol_identity.data(),
        static_cast<int>(local_protocol_identity.size()),
        MPI_UINT64_T,
        MPI_MIN,
        MPI_COMM_WORLD);
    MPI_Allreduce(
        local_protocol_identity.data(),
        maximum_protocol_identity.data(),
        static_cast<int>(local_protocol_identity.size()),
        MPI_UINT64_T,
        MPI_MAX,
        MPI_COMM_WORLD);
    if (minimum_protocol_identity != maximum_protocol_identity) {
      throw std::runtime_error(
          "PM slab halo exchange ranks disagree on global shape, halo depth, "
          "boundary mode, or exchange sequence");
    }
  }
#endif
  if (local_preparation_failure) {
    std::rethrow_exception(local_preparation_failure);
  }
  if (no_exchange) {
    return result;
  }

#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  constexpr int k_pm_halo_tag_base = 8810;
  constexpr int k_send_left_side = 0;
  constexpr int k_send_right_side = 1;
  const auto edge_index = [&](int peer) {
    const int local = mpi_context.worldRank();
    return (std::abs(local - peer) == 1) ? std::min(local, peer) : (layout.world_size - 1);
  };
  const auto side_tag = [&](int peer, int side) {
    return k_pm_halo_tag_base + edge_index(peer) * 2 + side;
  };

  const int local_rank = mpi_context.worldRank();
  const auto is_remote_peer = [&](int peer) {
    return peer >= 0 && peer != local_rank;
  };
  if (left_peer == local_rank) {
    std::copy(send_right.begin(), send_right.end(), result.left_halo.begin());
  }
  if (right_peer == local_rank) {
    std::copy(send_left.begin(), send_left.end(), result.right_halo.begin());
  }

  std::array<MPI_Request, 4> requests{};
  int request_count = 0;
  const int mpi_value_count = static_cast<int>(halo_value_count);
  const auto post_receive = [&](int peer, std::vector<double>& receive, int sender_side) {
    if (!is_remote_peer(peer)) {
      return;
    }
    MPI_Irecv(
        receive.data(),
        mpi_value_count,
        MPI_DOUBLE,
        peer,
        ghostExchangeSequencedTag(side_tag(peer, sender_side), local_rank, peer, exchange_sequence),
        MPI_COMM_WORLD,
        &requests[static_cast<std::size_t>(request_count++)]);
  };
  const auto post_send = [&](int peer, const std::vector<double>& send, int sender_side) {
    if (!is_remote_peer(peer)) {
      return;
    }
    MPI_Isend(
        const_cast<double*>(send.data()),
        mpi_value_count,
        MPI_DOUBLE,
        peer,
        ghostExchangeSequencedTag(side_tag(peer, sender_side), local_rank, peer, exchange_sequence),
        MPI_COMM_WORLD,
        &requests[static_cast<std::size_t>(request_count++)]);
  };

  // Post both receives before either send. Distinct side tags keep two-sided
  // traffic unambiguous when the periodic left and right owner are the same rank.
  post_receive(left_peer, result.left_halo, k_send_right_side);
  post_receive(right_peer, result.right_halo, k_send_left_side);
  post_send(left_peer, send_left, k_send_left_side);
  post_send(right_peer, send_right, k_send_right_side);
  if (request_count > 0) {
    MPI_Waitall(request_count, requests.data(), MPI_STATUSES_IGNORE);
  }
  return result;
#else
  throw std::runtime_error("PM slab halo exchange requires MPI support when MPI context is enabled");
#endif
}

GhostRefreshCommitReport commitBlockingGhostRefreshResult(
    GhostExchangeBufferSoA& ghost_storage,
    std::span<const LocalGhostDescriptor> local_ghost_descriptors,
    const GhostExchangePlan& plan,
    const BlockingGhostExchangeResult& result,
    const GhostLayerEpoch& expected_epoch) {
  validateGhostExchangePlan(plan);
  if (!plan.epoch.matches(expected_epoch)) {
    throw std::invalid_argument("commitBlockingGhostRefreshResult: ghost exchange plan epoch is stale");
  }
  if (!ghost_storage.isConsistent() || !result.received_ghosts.isConsistent()) {
    throw std::invalid_argument("commitBlockingGhostRefreshResult: ghost storage and result payloads must be component-consistent");
  }
  if (ghost_storage.size() < local_ghost_descriptors.size()) {
    throw std::invalid_argument("commitBlockingGhostRefreshResult: ghost storage must expose one slot per local descriptor");
  }
  if (!result.received_ghosts.epoch.matches(expected_epoch)) {
    throw std::invalid_argument("commitBlockingGhostRefreshResult: received ghost payload epoch is stale");
  }

  std::size_t expected_count = 0;
  for (const auto& indices : plan.recv_local_indices_by_neighbor) {
    expected_count += indices.size();
  }
  if (result.received_ghosts.size() != expected_count) {
    throw std::invalid_argument("commitBlockingGhostRefreshResult: received payload count does not match plan receive slots");
  }

  GhostRefreshCommitReport report;
  std::size_t result_row = 0;
  for (std::size_t slot = 0; slot < plan.recv_local_indices_by_neighbor.size(); ++slot) {
    for (const std::uint32_t local_index : plan.recv_local_indices_by_neighbor[slot]) {
      if (local_index >= local_ghost_descriptors.size() || local_index >= ghost_storage.size()) {
        throw std::out_of_range("commitBlockingGhostRefreshResult: receive slot index out of range");
      }
      const LocalGhostDescriptor descriptor = local_ghost_descriptors[local_index];
      if (descriptor.residency != LocalIndexResidency::kGhost || descriptor.owning_rank != plan.neighbor_ranks[slot]) {
        throw std::invalid_argument("commitBlockingGhostRefreshResult: receive slot is not a ghost owned by the exchange peer");
      }
      if (!descriptor.epoch.matches(expected_epoch)) {
        throw std::invalid_argument("commitBlockingGhostRefreshResult: local ghost descriptor is stale");
      }
      if (result.received_ghosts.entity_id[result_row] != descriptor.particle_id) {
        throw std::invalid_argument("commitBlockingGhostRefreshResult: received entity_id does not match ghost slot particle_id");
      }
      ghost_storage.entity_id[local_index] = result.received_ghosts.entity_id[result_row];
      const std::size_t storage_size = ghost_storage.size();
      commitOptionalGhostLane(&ghost_storage.position_x_comoving, result.received_ghosts.position_x_comoving, storage_size, local_index, result_row);
      commitOptionalGhostLane(&ghost_storage.position_y_comoving, result.received_ghosts.position_y_comoving, storage_size, local_index, result_row);
      commitOptionalGhostLane(&ghost_storage.position_z_comoving, result.received_ghosts.position_z_comoving, storage_size, local_index, result_row);
      commitOptionalGhostLane(&ghost_storage.mass_code, result.received_ghosts.mass_code, storage_size, local_index, result_row);
      commitOptionalGhostLane(&ghost_storage.density_code, result.received_ghosts.density_code, storage_size, local_index, result_row);
      commitOptionalGhostLane(&ghost_storage.velocity_x_code, result.received_ghosts.velocity_x_code, storage_size, local_index, result_row);
      commitOptionalGhostLane(&ghost_storage.velocity_y_code, result.received_ghosts.velocity_y_code, storage_size, local_index, result_row);
      commitOptionalGhostLane(&ghost_storage.velocity_z_code, result.received_ghosts.velocity_z_code, storage_size, local_index, result_row);
      commitOptionalGhostLane(&ghost_storage.pressure_code, result.received_ghosts.pressure_code, storage_size, local_index, result_row);
      commitOptionalGhostLane(&ghost_storage.internal_energy_code, result.received_ghosts.internal_energy_code, storage_size, local_index, result_row);
      ++result_row;
      ++report.updated_ghost_slots;
    }
  }
  ghost_storage.epoch = expected_epoch;
  report.committed_payload_bytes = static_cast<std::uint64_t>(report.updated_ghost_slots) *
      static_cast<std::uint64_t>(ghostExchangeRecordBytes());
  return report;
}

void invalidateGhostCache(GhostCacheLifecycle& lifecycle, const GhostLayerEpoch& next_epoch) {
  lifecycle.epoch = next_epoch;
  lifecycle.valid = false;
  ++lifecycle.invalidation_count;
}

void markGhostCacheCommitted(GhostCacheLifecycle& lifecycle, const GhostLayerEpoch& committed_epoch) {
  lifecycle.epoch = committed_epoch;
  lifecycle.valid = true;
  ++lifecycle.refresh_count;
}

void requireValidGhostCache(
    const GhostCacheLifecycle& lifecycle,
    const GhostLayerEpoch& expected_epoch,
    std::string_view caller) {
  if (!lifecycle.valid || !lifecycle.epoch.matches(expected_epoch)) {
    throw std::runtime_error(std::string(caller) + ": stale or invalid ghost cache used by solver");
  }
}


BlockingGhostRefreshExchange executeBlockingGhostRefreshExchangeFromDescriptors(
    const MpiContext& mpi_context,
    std::span<const LocalGhostDescriptor> local_ghost_descriptors,
    const GhostExchangeBufferSoA& authoritative_local_state,
    const GhostLayerEpoch& expected_epoch) {
  return executeBlockingGhostRefreshExchangeFromDescriptors(
      mpi_context,
      local_ghost_descriptors,
      makeReadOnlyGhostExchangeView(authoritative_local_state),
      expected_epoch);
}

BlockingGhostRefreshExchange executeBlockingGhostRefreshExchangeFromDescriptors(
    const MpiContext& mpi_context,
    std::span<const LocalGhostDescriptor> local_ghost_descriptors,
    const ReadOnlyGhostExchangeView& authoritative_local_state,
    const GhostLayerEpoch& expected_epoch) {
  const int world_rank = mpi_context.worldRank();
  const int world_size = mpi_context.worldSize();
  std::unordered_map<std::uint64_t, std::uint32_t> owned_index_by_particle_id;
  std::vector<std::vector<std::uint64_t>> requested_particle_ids_by_rank;
  std::vector<std::vector<std::uint32_t>> recv_indices_by_rank;
  std::exception_ptr local_preparation_failure;
  try {
    if (!authoritative_local_state.isConsistent()) {
      throw std::invalid_argument(
          "executeBlockingGhostRefreshExchangeFromDescriptors: authoritative payload state is inconsistent");
    }
    if (!authoritative_local_state.epoch.matches(expected_epoch)) {
      throw std::invalid_argument(
          "executeBlockingGhostRefreshExchangeFromDescriptors: authoritative payload epoch is stale");
    }
    if (authoritative_local_state.size() < local_ghost_descriptors.size()) {
      throw std::invalid_argument(
          "executeBlockingGhostRefreshExchangeFromDescriptors: payload state must expose one row per local descriptor");
    }
    if (world_rank < 0 || world_size <= 0 || world_rank >= world_size) {
      throw std::invalid_argument(
          "executeBlockingGhostRefreshExchangeFromDescriptors: invalid MPI context");
    }

    owned_index_by_particle_id.reserve(local_ghost_descriptors.size());
    requested_particle_ids_by_rank.resize(static_cast<std::size_t>(world_size));
    recv_indices_by_rank.resize(static_cast<std::size_t>(world_size));
    for (std::uint32_t local_index = 0; local_index < local_ghost_descriptors.size(); ++local_index) {
      const LocalGhostDescriptor descriptor = local_ghost_descriptors[local_index];
      if (!descriptor.epoch.matches(expected_epoch)) {
        throw std::invalid_argument(
            "executeBlockingGhostRefreshExchangeFromDescriptors: stale local ghost descriptor epoch");
      }
      if (descriptor.owning_rank < 0 || descriptor.owning_rank >= world_size) {
        throw std::invalid_argument(
            "executeBlockingGhostRefreshExchangeFromDescriptors: descriptor owning_rank outside MPI world");
      }
      if (authoritative_local_state.entity_id[local_index] != descriptor.particle_id) {
        throw std::invalid_argument(
            "executeBlockingGhostRefreshExchangeFromDescriptors: payload entity_id does not match descriptor particle_id");
      }
      if (descriptor.residency == LocalIndexResidency::kOwned) {
        if (descriptor.owning_rank != world_rank) {
          throw std::invalid_argument(
              "executeBlockingGhostRefreshExchangeFromDescriptors: owned descriptor has nonlocal owner");
        }
        auto [_, inserted] = owned_index_by_particle_id.emplace(descriptor.particle_id, local_index);
        if (!inserted) {
          throw std::invalid_argument(
              "executeBlockingGhostRefreshExchangeFromDescriptors: duplicate owned particle_id in descriptor table");
        }
      } else {
        if (descriptor.owning_rank == world_rank) {
          throw std::invalid_argument(
              "executeBlockingGhostRefreshExchangeFromDescriptors: ghost descriptor is owned by local rank");
        }
        requested_particle_ids_by_rank[static_cast<std::size_t>(descriptor.owning_rank)].push_back(
            descriptor.particle_id);
        recv_indices_by_rank[static_cast<std::size_t>(descriptor.owning_rank)].push_back(local_index);
      }
    }
  } catch (...) {
    local_preparation_failure = std::current_exception();
  }
  mpi_context.rethrowCollectivePreparationFailure(
      local_preparation_failure,
      "ghost refresh descriptor/request preparation");
  const bool has_local_ghost_demands = std::any_of(
      recv_indices_by_rank.begin(), recv_indices_by_rank.end(), [](const auto& rows) { return !rows.empty(); });
  if (!mpi_context.isEnabled()) {
    if (has_local_ghost_demands) {
      throw std::runtime_error(
          "executeBlockingGhostRefreshExchangeFromDescriptors: non-empty ghost demand requires MPI");
    }
    return BlockingGhostRefreshExchange{
        .plan = buildExplicitGhostExchangePlan(
            world_rank,
            std::span<const int>{},
            std::span<const std::vector<std::uint32_t>>{},
            std::span<const std::vector<std::uint32_t>>{},
            ghostRefreshPayloadRecordBytes(),
            expected_epoch),
        .result = BlockingGhostExchangeResult{.received_ghosts = GhostExchangeBufferSoA{
            .epoch = expected_epoch,
            .entity_id = {},
            .position_x_comoving = {},
            .position_y_comoving = {},
            .position_z_comoving = {},
            .mass_code = {},
            .density_code = {},
            .velocity_x_code = {},
            .velocity_y_code = {},
            .velocity_z_code = {},
            .pressure_code = {},
            .internal_energy_code = {},
        }},
    };
  }

#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  constexpr int k_request_count_tag_base = 4810;
  constexpr int k_request_ready_tag_base = 5310;
  constexpr int k_request_payload_tag_base = 5810;
  constexpr int k_request_mapping_ready_tag_base = 6310;
  std::vector<std::vector<std::uint32_t>> send_indices_by_rank;
  local_preparation_failure = nullptr;
  try {
    send_indices_by_rank.resize(static_cast<std::size_t>(world_size));
  } catch (...) {
    local_preparation_failure = std::current_exception();
  }
  mpi_context.rethrowCollectivePreparationFailure(
      local_preparation_failure,
      "ghost refresh request exchange storage preparation");

  std::exception_ptr first_protocol_failure;
  const auto record_peer_failure = [&](int peer_rank, std::string_view phase) {
    if (first_protocol_failure != nullptr) {
      return;
    }
    try {
      throw std::runtime_error(
          "ghost refresh peer " + std::to_string(peer_rank) +
          " rejected " + std::string(phase));
    } catch (...) {
      first_protocol_failure = std::current_exception();
    }
  };

  for (int peer_rank = 0; peer_rank < world_size; ++peer_rank) {
    if (peer_rank == world_rank) {
      continue;
    }
    const auto& request_ids =
        requested_particle_ids_by_rank[static_cast<std::size_t>(peer_rank)];
    const std::uint64_t send_count = core::checkedIntegralNarrow<std::uint64_t>(
        request_ids.size(), "ghost refresh request send count");
    std::uint64_t recv_count = 0U;
    if (MPI_Sendrecv(
            const_cast<std::uint64_t*>(&send_count), 1, MPI_UINT64_T,
            peer_rank,
            ghostExchangeSequencedTag(
                k_request_count_tag_base, world_rank, peer_rank,
                expected_epoch.ghost_sync_epoch),
            &recv_count, 1, MPI_UINT64_T,
            peer_rank,
            ghostExchangeSequencedTag(
                k_request_count_tag_base, world_rank, peer_rank,
                expected_epoch.ghost_sync_epoch),
            MPI_COMM_WORLD, MPI_STATUS_IGNORE) != MPI_SUCCESS) {
      throw std::runtime_error(
          "executeBlockingGhostRefreshExchangeFromDescriptors: request-count Sendrecv failed");
    }

    std::vector<std::uint64_t> received_request_ids;
    int request_send_count = 0;
    int request_recv_count = 0;
    std::exception_ptr peer_local_failure = first_protocol_failure;
    if (peer_local_failure == nullptr) {
      try {
        request_send_count = core::checkedIntegralNarrow<int>(
            send_count, "ghost refresh request send count");
        const std::size_t checked_recv_count = core::checkedIntegralNarrow<std::size_t>(
            recv_count, "ghost refresh request receive count");
        request_recv_count = core::checkedIntegralNarrow<int>(
            checked_recv_count, "ghost refresh request receive MPI count");
        received_request_ids.resize(checked_recv_count);
        injectMpiTestFault(mpi_context, "ghost_request_post_metadata");
      } catch (...) {
        peer_local_failure = std::current_exception();
        if (first_protocol_failure == nullptr) {
          first_protocol_failure = peer_local_failure;
        }
      }
    }

    int local_ready = peer_local_failure == nullptr ? 1 : 0;
    int peer_ready = 0;
    if (MPI_Sendrecv(
            &local_ready, 1, MPI_INT, peer_rank,
            ghostExchangeSequencedTag(
                k_request_ready_tag_base, world_rank, peer_rank,
                expected_epoch.ghost_sync_epoch),
            &peer_ready, 1, MPI_INT, peer_rank,
            ghostExchangeSequencedTag(
                k_request_ready_tag_base, world_rank, peer_rank,
                expected_epoch.ghost_sync_epoch),
            MPI_COMM_WORLD, MPI_STATUS_IGNORE) != MPI_SUCCESS) {
      throw std::runtime_error(
          "executeBlockingGhostRefreshExchangeFromDescriptors: request-readiness Sendrecv failed");
    }
    if (local_ready == 0 || peer_ready == 0) {
      if (peer_ready == 0) {
        record_peer_failure(peer_rank, "request payload preparation");
      }
      continue;
    }

    if (MPI_Sendrecv(
            request_ids.empty() ? nullptr : const_cast<std::uint64_t*>(request_ids.data()),
            request_send_count,
            MPI_UINT64_T, peer_rank,
            ghostExchangeSequencedTag(
                k_request_payload_tag_base, world_rank, peer_rank,
                expected_epoch.ghost_sync_epoch),
            received_request_ids.empty() ? nullptr : received_request_ids.data(),
            request_recv_count,
            MPI_UINT64_T, peer_rank,
            ghostExchangeSequencedTag(
                k_request_payload_tag_base, world_rank, peer_rank,
                expected_epoch.ghost_sync_epoch),
            MPI_COMM_WORLD, MPI_STATUS_IGNORE) != MPI_SUCCESS) {
      throw std::runtime_error(
          "executeBlockingGhostRefreshExchangeFromDescriptors: request-payload Sendrecv failed");
    }

    peer_local_failure = nullptr;
    try {
      auto& send_rows = send_indices_by_rank[static_cast<std::size_t>(peer_rank)];
      send_rows.reserve(received_request_ids.size());
      for (const std::uint64_t particle_id : received_request_ids) {
        const auto it = owned_index_by_particle_id.find(particle_id);
        if (it == owned_index_by_particle_id.end()) {
          throw std::runtime_error(
              "executeBlockingGhostRefreshExchangeFromDescriptors: peer requested a particle_id not owned by this rank");
        }
        send_rows.push_back(it->second);
      }
      injectMpiTestFault(mpi_context, "ghost_request_post_payload");
    } catch (...) {
      peer_local_failure = std::current_exception();
      if (first_protocol_failure == nullptr) {
        first_protocol_failure = peer_local_failure;
      }
    }

    local_ready = peer_local_failure == nullptr ? 1 : 0;
    peer_ready = 0;
    if (MPI_Sendrecv(
            &local_ready, 1, MPI_INT, peer_rank,
            ghostExchangeSequencedTag(
                k_request_mapping_ready_tag_base, world_rank, peer_rank,
                expected_epoch.ghost_sync_epoch),
            &peer_ready, 1, MPI_INT, peer_rank,
            ghostExchangeSequencedTag(
                k_request_mapping_ready_tag_base, world_rank, peer_rank,
                expected_epoch.ghost_sync_epoch),
            MPI_COMM_WORLD, MPI_STATUS_IGNORE) != MPI_SUCCESS) {
      throw std::runtime_error(
          "executeBlockingGhostRefreshExchangeFromDescriptors: request-mapping readiness Sendrecv failed");
    }
    if (peer_ready == 0) {
      record_peer_failure(peer_rank, "request mapping validation");
    }
  }

  mpi_context.rethrowCollectivePreparationFailure(
      first_protocol_failure,
      "ghost refresh request discovery");
  std::vector<int> neighbor_ranks;
  std::vector<std::vector<std::uint32_t>> send_indices_by_neighbor;
  std::vector<std::vector<std::uint32_t>> recv_indices_by_neighbor;
  BlockingGhostRefreshExchange exchange;
  local_preparation_failure = nullptr;
  try {
    for (int peer_rank = 0; peer_rank < world_size; ++peer_rank) {
      if (peer_rank == world_rank) {
        continue;
      }
      const auto& send_rows = send_indices_by_rank[static_cast<std::size_t>(peer_rank)];
      const auto& recv_rows = recv_indices_by_rank[static_cast<std::size_t>(peer_rank)];
      if (send_rows.empty() && recv_rows.empty()) {
        continue;
      }
      neighbor_ranks.push_back(peer_rank);
      send_indices_by_neighbor.push_back(send_rows);
      recv_indices_by_neighbor.push_back(recv_rows);
    }
    exchange.plan = buildExplicitGhostExchangePlan(
        world_rank,
        neighbor_ranks,
        send_indices_by_neighbor,
        recv_indices_by_neighbor,
        ghostRefreshPayloadRecordBytes(),
        expected_epoch);
    exchange.plan.exchange_sequence = expected_epoch.ghost_sync_epoch;
  } catch (...) {
    local_preparation_failure = std::current_exception();
  }
  mpi_context.rethrowCollectivePreparationFailure(
      local_preparation_failure,
      "ghost refresh payload-plan preparation");
  exchange.result = executeBlockingGhostRefreshExchange(
      mpi_context,
      exchange.plan,
      local_ghost_descriptors,
      authoritative_local_state,
      expected_epoch);
  return exchange;
#else
  throw std::runtime_error(
      "executeBlockingGhostRefreshExchangeFromDescriptors: MPI support is not compiled in");
#endif
}

BlockingGhostExchangeResult executeBlockingGhostRefreshExchange(
    const MpiContext& mpi_context,
    const GhostExchangePlan& plan,
    std::span<const LocalGhostDescriptor> local_ghost_descriptors,
    const GhostExchangeBufferSoA& authoritative_local_state,
    const GhostLayerEpoch& expected_epoch) {
  return executeBlockingGhostRefreshExchange(
      mpi_context,
      plan,
      local_ghost_descriptors,
      makeReadOnlyGhostExchangeView(authoritative_local_state),
      expected_epoch);
}

BlockingGhostExchangeResult executeBlockingGhostRefreshExchange(
    const MpiContext& mpi_context,
    const GhostExchangePlan& plan,
    std::span<const LocalGhostDescriptor> local_ghost_descriptors,
    const ReadOnlyGhostExchangeView& authoritative_local_state,
    const GhostLayerEpoch& expected_epoch) {
  BlockingGhostExchangeResult result;
  result.received_ghosts.epoch = expected_epoch;

  std::vector<GhostExchangeBuffer> send_buffers;
  std::exception_ptr local_preparation_failure;
  try {
    validateBlockingGhostExchangeContracts(
        plan, local_ghost_descriptors, mpi_context.worldRank(), expected_epoch);
    if (!authoritative_local_state.isConsistent()) {
      throw std::invalid_argument(
          "executeBlockingGhostRefreshExchange: authoritative local ghost payload state is inconsistent");
    }
    if (!authoritative_local_state.epoch.matches(expected_epoch)) {
      throw std::invalid_argument(
          "executeBlockingGhostRefreshExchange: authoritative local payload epoch is stale");
    }
    if (!plan.neighbor_ranks.empty() && !mpi_context.isEnabled()) {
      throw std::runtime_error(
          "executeBlockingGhostRefreshExchange: non-empty ghost exchange requires MPI; serial path must have no neighbors");
    }

    send_buffers.resize(plan.neighbor_ranks.size());
    std::size_t total_receive_records = 0U;
    for (std::size_t slot = 0; slot < plan.neighbor_ranks.size(); ++slot) {
      for (const std::uint32_t local_index : plan.send_local_indices_by_neighbor[slot]) {
        if (local_index >= authoritative_local_state.size() ||
            local_index >= local_ghost_descriptors.size()) {
          throw std::out_of_range(
              "executeBlockingGhostRefreshExchange: send descriptor index is outside local payload state");
        }
        if (authoritative_local_state.entity_id[local_index] !=
            local_ghost_descriptors[local_index].particle_id) {
          throw std::invalid_argument(
              "executeBlockingGhostRefreshExchange: send payload entity_id does not match owned descriptor particle_id");
        }
      }
      send_buffers[slot].packFrom(
          plan.outbound_transfers[slot],
          authoritative_local_state,
          plan.send_local_indices_by_neighbor[slot]);
      total_receive_records = core::checkedSizeAdd(
          total_receive_records,
          plan.recv_local_indices_by_neighbor[slot].size(),
          "ghost refresh total receive record count");
    }

    // Reserve only lanes that are semantically present in the authoritative
    // payload. Optional-lane absence is part of the wire contract: in
    // particular, generic DMO particle ghosts must not allocate or advertise
    // gas/hydro lanes merely because an exchange occurs.
    result.received_ghosts.entity_id.reserve(total_receive_records);
    const auto reserve_if_source_present = [total_receive_records](
                                                std::span<const double> source,
                                                std::vector<double>* destination) {
      if (!source.empty()) {
        destination->reserve(total_receive_records);
      }
    };
    reserve_if_source_present(authoritative_local_state.position_x_comoving, &result.received_ghosts.position_x_comoving);
    reserve_if_source_present(authoritative_local_state.position_y_comoving, &result.received_ghosts.position_y_comoving);
    reserve_if_source_present(authoritative_local_state.position_z_comoving, &result.received_ghosts.position_z_comoving);
    reserve_if_source_present(authoritative_local_state.mass_code, &result.received_ghosts.mass_code);
    reserve_if_source_present(authoritative_local_state.density_code, &result.received_ghosts.density_code);
    reserve_if_source_present(authoritative_local_state.velocity_x_code, &result.received_ghosts.velocity_x_code);
    reserve_if_source_present(authoritative_local_state.velocity_y_code, &result.received_ghosts.velocity_y_code);
    reserve_if_source_present(authoritative_local_state.velocity_z_code, &result.received_ghosts.velocity_z_code);
    reserve_if_source_present(authoritative_local_state.pressure_code, &result.received_ghosts.pressure_code);
    reserve_if_source_present(authoritative_local_state.internal_energy_code, &result.received_ghosts.internal_energy_code);
  } catch (...) {
    local_preparation_failure = std::current_exception();
  }

  if (mpi_context.isEnabled()) {
    mpi_context.rethrowCollectivePreparationFailure(
        local_preparation_failure,
        "ghost refresh payload pre-communication preparation");
  } else if (local_preparation_failure != nullptr) {
    std::rethrow_exception(local_preparation_failure);
  }
  // A rank with no local peer payload still remains a participant in the
  // MPI-world protocol: MPI-enabled ranks must reach the final distributed
  // failure agreement below even when this peer loop is empty. The serial
  // no-neighbor path can return locally because it has no world collective.
  if (!mpi_context.isEnabled() && plan.neighbor_ranks.empty()) {
    return result;
  }

#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  constexpr int k_size_tag_base = 6810;
  constexpr int k_ready_tag_base = 7310;
  constexpr int k_payload_tag_base = 7810;
  constexpr int k_decode_ready_tag_base = 8310;
  std::exception_ptr first_protocol_failure;
  const auto record_peer_failure = [&](int peer_rank, std::string_view phase) {
    if (first_protocol_failure != nullptr) {
      return;
    }
    try {
      throw std::runtime_error(
          "ghost payload peer " + std::to_string(peer_rank) +
          " rejected " + std::string(phase));
    } catch (...) {
      first_protocol_failure = std::current_exception();
    }
  };

  for (std::size_t slot = 0; slot < plan.neighbor_ranks.size(); ++slot) {
    const auto send_bytes = send_buffers[slot].encodedBytes();
    const std::uint64_t send_size = core::checkedIntegralNarrow<std::uint64_t>(
        send_bytes.size(), "ghost refresh payload send size");
    std::uint64_t recv_size = 0U;
    const int peer_rank = plan.neighbor_ranks[slot];
    if (MPI_Sendrecv(
            const_cast<std::uint64_t*>(&send_size), 1, MPI_UINT64_T,
            peer_rank,
            ghostExchangeSequencedTag(
                k_size_tag_base, mpi_context.worldRank(), peer_rank,
                plan.exchange_sequence),
            &recv_size, 1, MPI_UINT64_T,
            peer_rank,
            ghostExchangeSequencedTag(
                k_size_tag_base, mpi_context.worldRank(), peer_rank,
                plan.exchange_sequence),
            MPI_COMM_WORLD, MPI_STATUS_IGNORE) != MPI_SUCCESS) {
      throw std::runtime_error(
          "executeBlockingGhostRefreshExchange: payload-size Sendrecv failed");
    }

    std::vector<std::uint8_t> recv_bytes;
    int payload_send_count = 0;
    int payload_recv_count = 0;
    std::exception_ptr peer_local_failure = first_protocol_failure;
    if (peer_local_failure == nullptr) {
      try {
        payload_send_count = core::checkedIntegralNarrow<int>(
            send_size, "ghost refresh payload send MPI count");
        const std::size_t checked_recv_size = core::checkedIntegralNarrow<std::size_t>(
            recv_size, "ghost refresh payload receive size");
        payload_recv_count = core::checkedIntegralNarrow<int>(
            checked_recv_size, "ghost refresh payload receive MPI count");
        recv_bytes.resize(checked_recv_size);
        injectMpiTestFault(mpi_context, "ghost_payload_post_metadata");
      } catch (...) {
        peer_local_failure = std::current_exception();
        if (first_protocol_failure == nullptr) {
          first_protocol_failure = peer_local_failure;
        }
      }
    }

    int local_ready = peer_local_failure == nullptr ? 1 : 0;
    int peer_ready = 0;
    if (MPI_Sendrecv(
            &local_ready, 1, MPI_INT, peer_rank,
            ghostExchangeSequencedTag(
                k_ready_tag_base, mpi_context.worldRank(), peer_rank,
                plan.exchange_sequence),
            &peer_ready, 1, MPI_INT, peer_rank,
            ghostExchangeSequencedTag(
                k_ready_tag_base, mpi_context.worldRank(), peer_rank,
                plan.exchange_sequence),
            MPI_COMM_WORLD, MPI_STATUS_IGNORE) != MPI_SUCCESS) {
      throw std::runtime_error(
          "executeBlockingGhostRefreshExchange: payload-readiness Sendrecv failed");
    }
    if (local_ready == 0 || peer_ready == 0) {
      if (peer_ready == 0) {
        record_peer_failure(peer_rank, "payload preparation");
      }
      continue;
    }

    if (MPI_Sendrecv(
            send_bytes.empty() ? nullptr : const_cast<std::uint8_t*>(send_bytes.data()),
            payload_send_count,
            MPI_BYTE, peer_rank,
            ghostExchangeSequencedTag(
                k_payload_tag_base, mpi_context.worldRank(), peer_rank,
                plan.exchange_sequence),
            recv_bytes.empty() ? nullptr : recv_bytes.data(),
            payload_recv_count,
            MPI_BYTE, peer_rank,
            ghostExchangeSequencedTag(
                k_payload_tag_base, mpi_context.worldRank(), peer_rank,
                plan.exchange_sequence),
            MPI_COMM_WORLD, MPI_STATUS_IGNORE) != MPI_SUCCESS) {
      throw std::runtime_error(
          "executeBlockingGhostRefreshExchange: ghost payload Sendrecv failed");
    }

    peer_local_failure = nullptr;
    try {
      GhostExchangeBuffer recv_buffer;
      recv_buffer.replaceEncodedBytes(std::move(recv_bytes));
      const std::size_t old_received_count = result.received_ghosts.size();
      recv_buffer.unpackAppendTo(
          plan.inbound_transfers[slot], result.received_ghosts);
      for (std::size_t i = 0;
           i < plan.recv_local_indices_by_neighbor[slot].size(); ++i) {
        const std::uint32_t local_index =
            plan.recv_local_indices_by_neighbor[slot][i];
        if (local_index >= local_ghost_descriptors.size()) {
          throw std::out_of_range(
              "executeBlockingGhostRefreshExchange: receive descriptor index is outside residency table");
        }
        if (result.received_ghosts.entity_id[old_received_count + i] !=
            local_ghost_descriptors[local_index].particle_id) {
          throw std::invalid_argument(
              "executeBlockingGhostRefreshExchange: received ghost entity_id does not match receive descriptor particle_id");
        }
      }
      injectMpiTestFault(mpi_context, "ghost_payload_post_payload");
      result.sent_bytes = core::checkedSizeAdd(
          core::checkedIntegralNarrow<std::size_t>(
              result.sent_bytes, "ghost refresh sent byte accumulator"),
          core::checkedIntegralNarrow<std::size_t>(
              send_size, "ghost refresh sent byte count"),
          "ghost refresh sent byte aggregate");
      result.received_bytes = core::checkedSizeAdd(
          core::checkedIntegralNarrow<std::size_t>(
              result.received_bytes, "ghost refresh received byte accumulator"),
          core::checkedIntegralNarrow<std::size_t>(
              recv_size, "ghost refresh received byte count"),
          "ghost refresh received byte aggregate");
    } catch (...) {
      peer_local_failure = std::current_exception();
      if (first_protocol_failure == nullptr) {
        first_protocol_failure = peer_local_failure;
      }
    }

    local_ready = peer_local_failure == nullptr ? 1 : 0;
    peer_ready = 0;
    if (MPI_Sendrecv(
            &local_ready, 1, MPI_INT, peer_rank,
            ghostExchangeSequencedTag(
                k_decode_ready_tag_base, mpi_context.worldRank(), peer_rank,
                plan.exchange_sequence),
            &peer_ready, 1, MPI_INT, peer_rank,
            ghostExchangeSequencedTag(
                k_decode_ready_tag_base, mpi_context.worldRank(), peer_rank,
                plan.exchange_sequence),
            MPI_COMM_WORLD, MPI_STATUS_IGNORE) != MPI_SUCCESS) {
      throw std::runtime_error(
          "executeBlockingGhostRefreshExchange: payload-decode readiness Sendrecv failed");
    }
    if (peer_ready == 0) {
      record_peer_failure(peer_rank, "payload decode/validation");
    }
  }

  mpi_context.rethrowCollectivePreparationFailure(
      first_protocol_failure,
      "ghost refresh sequential peer payload exchange");
  return result;
#else
  throw std::runtime_error(
      "executeBlockingGhostRefreshExchange: MPI support is not compiled in");
#endif
}


PmSlabHaloExchangeStatus executeBlockingPmSlabHaloExchangeInto(
    const MpiContext& mpi_context,
    const PmSlabLayout& layout,
    std::span<const double> local_scalar_field,
    std::size_t halo_depth_x,
    bool periodic_x,
    std::span<double> left_halo_out,
    std::span<double> right_halo_out,
    std::pmr::memory_resource* communication_resource,
    std::uint64_t exchange_sequence) {
#if !defined(COSMOSIM_ENABLE_MPI) || !COSMOSIM_ENABLE_MPI
  (void)exchange_sequence;
#endif
  if (communication_resource == nullptr) {
    throw std::invalid_argument(
        "PM slab halo arena exchange requires a non-null communication resource");
  }
  PmSlabHaloExchangeStatus result;
  std::pmr::vector<double> send_left(communication_resource);
  std::pmr::vector<double> send_right(communication_resource);
  std::size_t halo_value_count = 0U;
  std::uint64_t payload_bytes = 0U;
  int left_peer = -1;
  int right_peer = -1;
  bool no_exchange = false;
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  int communicator_world_size = 1;
  int communicator_world_rank = 0;
  const bool communicator_mpi_active =
      queryActiveMpiWorld(communicator_world_size, communicator_world_rank);
#endif

  std::exception_ptr local_preparation_failure;
  try {
    if (!layout.isValid()) {
      throw std::invalid_argument("PM slab halo exchange requires a valid slab layout");
    }
    if (layout.world_size != mpi_context.worldSize() ||
        layout.world_rank != mpi_context.worldRank()) {
      throw std::invalid_argument(
          "PM slab halo exchange layout world metadata must match MPI context");
    }
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
    if (!communicator_mpi_active &&
        (layout.world_size > 1 || mpi_context.isEnabled())) {
      throw std::invalid_argument(
          "PM slab halo exchange requires an active MPI_COMM_WORLD for an enabled or distributed context");
    }
    if (communicator_mpi_active &&
        (layout.world_size != communicator_world_size ||
         layout.world_rank != communicator_world_rank)) {
      throw std::invalid_argument(
          "PM slab halo exchange layout world metadata must match MPI_COMM_WORLD");
    }
#endif
    if (layout.global_ny >
        std::numeric_limits<std::size_t>::max() / layout.global_nz) {
      throw std::overflow_error("PM slab halo exchange plane size overflows size_t");
    }
    const std::size_t plane_size = layout.global_ny * layout.global_nz;
    if (layout.local_nx() >
        std::numeric_limits<std::size_t>::max() / plane_size) {
      throw std::overflow_error("PM slab halo exchange local field size overflows size_t");
    }
    const std::size_t expected_local_values = layout.local_nx() * plane_size;
    if (local_scalar_field.size() != expected_local_values) {
      throw std::invalid_argument(
          "PM slab halo exchange field size does not match local slab cell count");
    }
    if (halo_depth_x == 0 || layout.world_size == 1 || layout.local_nx() == 0) {
      no_exchange = true;
    } else {
      if (!mpi_context.isEnabled()) {
        throw std::runtime_error(
            "PM slab halo exchange requires MPI for distributed layouts");
      }
      std::size_t minimum_nonempty_slab_nx =
          std::numeric_limits<std::size_t>::max();
      for (int rank = 0; rank < layout.world_size; ++rank) {
        const PmSlabRange owned =
            pmOwnedXRangeForRank(layout.global_nx, layout.world_size, rank);
        if (owned.extentX() > 0U) {
          minimum_nonempty_slab_nx =
              std::min(minimum_nonempty_slab_nx, owned.extentX());
        }
      }
      if (minimum_nonempty_slab_nx ==
          std::numeric_limits<std::size_t>::max()) {
        throw std::logic_error("PM slab halo exchange layout has no non-empty owner");
      }
      const std::size_t depth =
          std::min(halo_depth_x, minimum_nonempty_slab_nx);
      if (depth > std::numeric_limits<std::size_t>::max() / plane_size) {
        throw std::overflow_error("PM slab halo exchange payload size overflows size_t");
      }
      halo_value_count = depth * plane_size;
      if (halo_value_count >
          static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw std::overflow_error(
            "PM slab halo exchange payload count exceeds MPI int limit");
      }
      if (left_halo_out.size() < halo_value_count ||
          right_halo_out.size() < halo_value_count) {
        throw std::invalid_argument(
            "PM slab halo exchange output cache spans are smaller than the halo payload");
      }
      if (halo_value_count >
          std::numeric_limits<std::uint64_t>::max() / sizeof(double)) {
        throw std::overflow_error(
            "PM slab halo exchange byte diagnostics overflow uint64_t");
      }
      payload_bytes =
          static_cast<std::uint64_t>(halo_value_count) * sizeof(double);
      result.halo_depth_x = depth;
      if (layout.owned_x.begin_x > 0) {
        left_peer = pmOwnerRankForGlobalX(
            layout.global_nx, layout.world_size, layout.owned_x.begin_x - 1U);
      } else if (periodic_x) {
        left_peer = pmOwnerRankForGlobalX(
            layout.global_nx, layout.world_size, layout.global_nx - 1U);
      }
      if (layout.owned_x.end_x < layout.global_nx) {
        right_peer = pmOwnerRankForGlobalX(
            layout.global_nx, layout.world_size, layout.owned_x.end_x);
      } else if (periodic_x) {
        right_peer = pmOwnerRankForGlobalX(
            layout.global_nx, layout.world_size, 0U);
      }
      result.left_peer_rank = left_peer;
      result.right_peer_rank = right_peer;
      send_left.resize(halo_value_count);
      send_right.resize(halo_value_count);
      const std::span<const double> left_source =
          local_scalar_field.first(halo_value_count);
      const std::span<const double> right_source =
          local_scalar_field.last(halo_value_count);
      std::copy(left_source.begin(), left_source.end(), send_left.begin());
      std::copy(right_source.begin(), right_source.end(), send_right.begin());

      const int local_rank = mpi_context.worldRank();
      const auto is_remote_peer = [&](int peer) {
        return peer >= 0 && peer != local_rank;
      };
      const std::uint64_t remote_side_count =
          static_cast<std::uint64_t>(is_remote_peer(left_peer)) +
          static_cast<std::uint64_t>(is_remote_peer(right_peer));
      if (remote_side_count > 0 &&
          payload_bytes >
              std::numeric_limits<std::uint64_t>::max() / remote_side_count) {
        throw std::overflow_error(
            "PM slab halo exchange aggregate byte diagnostics overflow uint64_t");
      }
      result.sent_bytes = payload_bytes * remote_side_count;
      result.received_bytes = payload_bytes * remote_side_count;
    }
  } catch (...) {
    local_preparation_failure = std::current_exception();
  }

#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  if (communicator_mpi_active && communicator_world_size > 1) {
    const std::uint64_t local_failure_vote =
        local_preparation_failure ? 1U : 0U;
    std::uint64_t failure_count = 0U;
    MPI_Allreduce(
        &local_failure_vote, &failure_count, 1, MPI_UINT64_T, MPI_SUM,
        MPI_COMM_WORLD);
    if (failure_count != 0U) {
      if (local_preparation_failure) {
        std::rethrow_exception(local_preparation_failure);
      }
      throw std::runtime_error(
          "PM slab halo exchange peer rejected protocol preparation");
    }
    const std::array<std::uint64_t, 7> local_protocol_identity{
        static_cast<std::uint64_t>(layout.global_nx),
        static_cast<std::uint64_t>(layout.global_ny),
        static_cast<std::uint64_t>(layout.global_nz),
        static_cast<std::uint64_t>(halo_depth_x), periodic_x ? 1U : 0U,
        exchange_sequence, static_cast<std::uint64_t>(communicator_world_size),
    };
    std::array<std::uint64_t, 7> minimum_protocol_identity{};
    std::array<std::uint64_t, 7> maximum_protocol_identity{};
    MPI_Allreduce(
        local_protocol_identity.data(), minimum_protocol_identity.data(),
        static_cast<int>(local_protocol_identity.size()), MPI_UINT64_T,
        MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(
        local_protocol_identity.data(), maximum_protocol_identity.data(),
        static_cast<int>(local_protocol_identity.size()), MPI_UINT64_T,
        MPI_MAX, MPI_COMM_WORLD);
    if (minimum_protocol_identity != maximum_protocol_identity) {
      throw std::runtime_error(
          "PM slab halo exchange ranks disagree on global shape, halo depth, boundary mode, or exchange sequence");
    }
  }
#endif
  if (local_preparation_failure) {
    std::rethrow_exception(local_preparation_failure);
  }
  if (no_exchange) {
    return result;
  }

#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  constexpr int k_pm_halo_tag_base = 8810;
  constexpr int k_send_left_side = 0;
  constexpr int k_send_right_side = 1;
  const auto edge_index = [&](int peer) {
    const int local = mpi_context.worldRank();
    return (std::abs(local - peer) == 1) ? std::min(local, peer)
                                         : (layout.world_size - 1);
  };
  const auto side_tag = [&](int peer, int side) {
    return k_pm_halo_tag_base + edge_index(peer) * 2 + side;
  };

  const int local_rank = mpi_context.worldRank();
  const auto is_remote_peer = [&](int peer) {
    return peer >= 0 && peer != local_rank;
  };
  auto left_receive = left_halo_out.first(halo_value_count);
  auto right_receive = right_halo_out.first(halo_value_count);
  if (left_peer == local_rank) {
    std::copy(send_right.begin(), send_right.end(), left_receive.begin());
  }
  if (right_peer == local_rank) {
    std::copy(send_left.begin(), send_left.end(), right_receive.begin());
  }

  std::array<MPI_Request, 4> requests{};
  int request_count = 0;
  const int mpi_value_count = static_cast<int>(halo_value_count);
  const auto post_receive = [&](int peer, std::span<double> receive,
                                int sender_side) {
    if (!is_remote_peer(peer)) {
      return;
    }
    MPI_Irecv(
        receive.data(), mpi_value_count, MPI_DOUBLE, peer,
        ghostExchangeSequencedTag(
            side_tag(peer, sender_side), local_rank, peer, exchange_sequence),
        MPI_COMM_WORLD, &requests[static_cast<std::size_t>(request_count++)]);
  };
  const auto post_send = [&](int peer, std::span<const double> send,
                             int sender_side) {
    if (!is_remote_peer(peer)) {
      return;
    }
    MPI_Isend(
        const_cast<double*>(send.data()), mpi_value_count, MPI_DOUBLE, peer,
        ghostExchangeSequencedTag(
            side_tag(peer, sender_side), local_rank, peer, exchange_sequence),
        MPI_COMM_WORLD, &requests[static_cast<std::size_t>(request_count++)]);
  };

  post_receive(left_peer, left_receive, k_send_right_side);
  post_receive(right_peer, right_receive, k_send_left_side);
  post_send(left_peer, send_left, k_send_left_side);
  post_send(right_peer, send_right, k_send_right_side);
  if (request_count > 0) {
    MPI_Waitall(request_count, requests.data(), MPI_STATUSES_IGNORE);
  }
  return result;
#else
  throw std::runtime_error(
      "PM slab halo exchange requires MPI support when MPI context is enabled");
#endif
}

}  // namespace cosmosim::parallel
