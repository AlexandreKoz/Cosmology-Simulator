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
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
using internal::queryActiveMpiWorld;
using internal::queryNodeLocalTopology;
#endif

}  // namespace
MpiContext::MpiContext() {
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  m_is_enabled = queryActiveMpiWorld(m_world_size, m_world_rank);
  if (m_is_enabled) {
    queryNodeLocalTopology(m_world_rank, m_local_rank, m_local_size);
  }
#endif
}

MpiContext::MpiContext(bool is_enabled, int world_size, int world_rank)
    : m_is_enabled(is_enabled), m_world_size(world_size), m_world_rank(world_rank), m_local_rank(world_rank),
      m_local_size(is_enabled ? world_size : 1) {
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

int MpiContext::localRank() const noexcept { return m_local_rank; }

int MpiContext::localSize() const noexcept { return m_local_size; }

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

void MpiContext::allreduceSumDoublesInPlace(std::span<double> values) const {
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  if (m_is_enabled && !values.empty()) {
    if (values.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
      throw std::overflow_error("MPI double-vector reduction exceeds int count range");
    }
    if (MPI_Allreduce(
            MPI_IN_PLACE, values.data(), static_cast<int>(values.size()),
            MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD) != MPI_SUCCESS) {
      throw std::runtime_error("MPI_Allreduce failed for distributed double diagnostics");
    }
  }
#else
  (void)values;
#endif
}

double MpiContext::allreduceMinDouble(double local_value) const {
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  if (m_is_enabled) {
    double global = 0.0;
    MPI_Allreduce(&local_value, &global, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
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

void MpiContext::allreduceSumUint64sInPlace(std::span<std::uint64_t> values) const {
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  if (m_is_enabled && !values.empty()) {
    if (values.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
      throw std::overflow_error("MPI uint64-vector reduction exceeds int count range");
    }
    if (MPI_Allreduce(
            MPI_IN_PLACE, values.data(), static_cast<int>(values.size()),
            MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD) != MPI_SUCCESS) {
      throw std::runtime_error("MPI_Allreduce failed for distributed uint64 diagnostics");
    }
  }
#else
  (void)values;
#endif
}

std::uint64_t MpiContext::exclusiveScanSumUint64(std::uint64_t local_value) const {
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  if (m_is_enabled) {
    std::uint64_t prefix = 0U;
    const int rc = MPI_Exscan(
        &local_value, &prefix, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
    if (rc != MPI_SUCCESS) {
      throw std::runtime_error(
          "MPI_Exscan failed while computing distributed snapshot row offsets");
    }
    return m_world_rank == 0 ? 0U : prefix;
  }
#endif
  return 0U;
}

std::uint64_t MpiContext::allreduceMaxUint64(std::uint64_t local_value) const {
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  if (m_is_enabled) {
    std::uint64_t global = 0;
    MPI_Allreduce(&local_value, &global, 1, MPI_UINT64_T, MPI_MAX, MPI_COMM_WORLD);
    return global;
  }
#endif
  return local_value;
}

std::uint64_t MpiContext::allreduceMinUint64(std::uint64_t local_value) const {
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  if (m_is_enabled) {
    std::uint64_t global = 0;
    MPI_Allreduce(&local_value, &global, 1, MPI_UINT64_T, MPI_MIN, MPI_COMM_WORLD);
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


void MpiContext::rethrowCollectivePreparationFailure(
    const std::exception_ptr& local_failure,
    std::string_view phase_name) const {
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  if (m_is_enabled) {
    const int local_failed = local_failure != nullptr ? 1 : 0;
    int failed_count = 0;
    if (MPI_Allreduce(
            &local_failed, &failed_count, 1, MPI_INT, MPI_SUM,
            MPI_COMM_WORLD) != MPI_SUCCESS) {
      throw std::runtime_error(
          "MPI readiness Allreduce failed during collective preparation phase");
    }
    if (failed_count == 0) {
      return;
    }

    const int local_candidate = local_failure != nullptr ? m_world_rank : m_world_size;
    int failure_rank = m_world_size;
    if (MPI_Allreduce(
            &local_candidate, &failure_rank, 1, MPI_INT, MPI_MIN,
            MPI_COMM_WORLD) != MPI_SUCCESS) {
      throw std::runtime_error(
          "MPI failure-rank Allreduce failed during collective preparation phase");
    }

    constexpr std::size_t k_maximum_message_bytes = 2047U;
    std::array<char, k_maximum_message_bytes + 1U> message_buffer{};
    std::uint32_t message_length = 0U;
    if (m_world_rank == failure_rank) {
      const char* message = "unknown non-standard exception";
      try {
        std::rethrow_exception(local_failure);
      } catch (const std::exception& error) {
        message = error.what();
      } catch (...) {
      }
      const std::size_t raw_length = std::char_traits<char>::length(message);
      message_length = static_cast<std::uint32_t>(
          std::min(raw_length, k_maximum_message_bytes));
      std::copy_n(message, message_length, message_buffer.data());
    }
    if (MPI_Bcast(
            &message_length, 1, MPI_UINT32_T, failure_rank,
            MPI_COMM_WORLD) != MPI_SUCCESS ||
        MPI_Bcast(
            message_buffer.data(), static_cast<int>(message_buffer.size()), MPI_CHAR,
            failure_rank, MPI_COMM_WORLD) != MPI_SUCCESS) {
      throw std::runtime_error(
          "MPI diagnostic broadcast failed during collective preparation phase");
    }
    throw std::runtime_error(
        "collective preparation phase '" + std::string(phase_name) +
        "' failed on rank " + std::to_string(failure_rank) + ": " +
        std::string(message_buffer.data(),
                    message_buffer.data() + message_length));
  }
#endif
  if (local_failure != nullptr) {
    std::rethrow_exception(local_failure);
  }
}

std::vector<std::uint8_t> MpiContext::gatherBytesToRoot(
    std::span<const std::uint8_t> local_bytes,
    int root_rank) const {
  if (root_rank < 0 || root_rank >= m_world_size) {
    throw std::invalid_argument("MpiContext::gatherBytesToRoot invalid root rank");
  }
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  if (m_is_enabled) {
    std::uint64_t local_count64 = 0U;
    std::vector<std::uint64_t> counts64;
    std::exception_ptr local_failure;
    try {
      local_count64 = core::checkedIntegralNarrow<std::uint64_t>(
          local_bytes.size(), "gatherBytesToRoot local byte count");
      counts64.resize(static_cast<std::size_t>(m_world_size), 0U);
    } catch (...) {
      local_failure = std::current_exception();
    }
    rethrowCollectivePreparationFailure(
        local_failure, "gatherBytesToRoot count-buffer preparation");

    if (MPI_Allgather(
            &local_count64, 1, MPI_UINT64_T,
            counts64.data(), 1, MPI_UINT64_T,
            MPI_COMM_WORLD) != MPI_SUCCESS) {
      throw std::runtime_error("gatherBytesToRoot count Allgather failed");
    }

    BoundedMpiTransferPlan plan;
    std::vector<std::uint8_t> gathered;
    std::vector<std::uint8_t> round_receive_buffer;
    local_failure = nullptr;
    try {
      std::vector<std::size_t> logical_counts(counts64.size(), 0U);
      for (std::size_t rank = 0; rank < counts64.size(); ++rank) {
        logical_counts[rank] = core::checkedIntegralNarrow<std::size_t>(
            counts64[rank], "gatherBytesToRoot received byte count");
      }
      plan = planBoundedMpiTransferRounds(
          logical_counts,
          static_cast<std::size_t>(std::numeric_limits<int>::max()),
          mpiTransportRoundLimitBytes());
      if (m_world_rank == root_rank) {
        gathered.resize(plan.logical_total_count);
        std::size_t maximum_round_count = 0U;
        for (const auto& round : plan.rounds) {
          maximum_round_count = std::max(maximum_round_count, round.round_count);
        }
        round_receive_buffer.resize(maximum_round_count);
      }
      injectMpiTestFault(*this, "gather_post_count");
    } catch (...) {
      local_failure = std::current_exception();
    }
    rethrowCollectivePreparationFailure(
        local_failure, "gatherBytesToRoot payload preparation");

    for (const BoundedMpiRoundLayout& round : plan.rounds) {
      const std::size_t local_rank = static_cast<std::size_t>(m_world_rank);
      const int send_count = round.counts[local_rank];
      const std::size_t local_offset = round.logical_offsets[local_rank];
      const std::uint8_t* send_pointer = send_count == 0
          ? nullptr
          : local_bytes.data() + local_offset;
      const int status = MPI_Gatherv(
          const_cast<std::uint8_t*>(send_pointer), send_count, MPI_BYTE,
          m_world_rank == root_rank && !round_receive_buffer.empty()
              ? round_receive_buffer.data()
              : nullptr,
          m_world_rank == root_rank ? round.counts.data() : nullptr,
          m_world_rank == root_rank ? round.displacements.data() : nullptr,
          MPI_BYTE, root_rank, MPI_COMM_WORLD);
      if (status != MPI_SUCCESS) {
        throw std::runtime_error("gatherBytesToRoot bounded payload Gatherv failed");
      }
      if (m_world_rank == root_rank) {
        for (std::size_t rank = 0; rank < round.counts.size(); ++rank) {
          const std::size_t count = static_cast<std::size_t>(round.counts[rank]);
          if (count == 0U) {
            continue;
          }
          // The planner checked the logical prefix and every consumed peer offset
          // before any payload round, so this addition cannot overflow here.
          const std::size_t destination_offset =
              plan.logical_displacements[rank] + round.logical_offsets[rank];
          std::memcpy(
              gathered.data() + destination_offset,
              round_receive_buffer.data() +
                  static_cast<std::size_t>(round.displacements[rank]),
              count);
        }
      }
    }
    return gathered;
  }
#endif
  if (root_rank != 0) {
    throw std::invalid_argument("serial MpiContext only supports root rank 0");
  }
  return {local_bytes.begin(), local_bytes.end()};
}

std::vector<std::uint8_t> MpiContext::broadcastBytesFromRoot(
    std::span<const std::uint8_t> root_bytes,
    int root_rank) const {
  if (root_rank < 0 || root_rank >= m_world_size) {
    throw std::invalid_argument("MpiContext::broadcastBytesFromRoot invalid root rank");
  }
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  if (m_is_enabled) {
    std::uint64_t byte_count64 = 0U;
    std::exception_ptr local_failure;
    try {
      if (m_world_rank == root_rank) {
        byte_count64 = core::checkedIntegralNarrow<std::uint64_t>(
            root_bytes.size(), "broadcastBytesFromRoot root byte count");
      }
    } catch (...) {
      local_failure = std::current_exception();
    }
    rethrowCollectivePreparationFailure(
        local_failure, "broadcastBytesFromRoot size preparation");
    if (MPI_Bcast(
            &byte_count64, 1, MPI_UINT64_T, root_rank,
            MPI_COMM_WORLD) != MPI_SUCCESS) {
      throw std::runtime_error("broadcastBytesFromRoot size Bcast failed");
    }

    std::vector<std::uint8_t> bytes;
    local_failure = nullptr;
    try {
      const std::size_t byte_count = core::checkedIntegralNarrow<std::size_t>(
          byte_count64, "broadcastBytesFromRoot receive byte count");
      bytes.resize(byte_count);
      if (m_world_rank == root_rank && !root_bytes.empty()) {
        std::copy(root_bytes.begin(), root_bytes.end(), bytes.begin());
      }
      injectMpiTestFault(*this, "broadcast_post_count");
    } catch (...) {
      local_failure = std::current_exception();
    }
    rethrowCollectivePreparationFailure(
        local_failure, "broadcastBytesFromRoot payload preparation");

    const std::size_t round_limit = std::min(
        mpiTransportRoundLimitBytes(),
        static_cast<std::size_t>(std::numeric_limits<int>::max()));
    for (std::size_t offset = 0U; offset < bytes.size();) {
      const std::size_t chunk_size = std::min(round_limit, bytes.size() - offset);
      const int chunk_count = core::checkedIntegralNarrow<int>(
          chunk_size, "broadcastBytesFromRoot bounded chunk count");
      if (MPI_Bcast(
              bytes.data() + offset, chunk_count, MPI_BYTE,
              root_rank, MPI_COMM_WORLD) != MPI_SUCCESS) {
        throw std::runtime_error("broadcastBytesFromRoot bounded payload Bcast failed");
      }
      // chunk_size is bounded by bytes.size() - offset.
      offset += chunk_size;
    }
    return bytes;
  }
#endif
  if (root_rank != 0) {
    throw std::invalid_argument("serial MpiContext only supports root rank 0");
  }
  return {root_bytes.begin(), root_bytes.end()};
}

std::vector<std::uint8_t> MpiContext::allgatherBytesBounded(
    std::span<const std::uint8_t> local_bytes) const {
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  if (m_is_enabled) {
    std::uint64_t local_count64 = 0U;
    std::vector<std::uint64_t> counts64;
    std::exception_ptr local_failure;
    try {
      local_count64 = core::checkedIntegralNarrow<std::uint64_t>(
          local_bytes.size(), "allgatherBytesBounded local byte count");
      counts64.resize(static_cast<std::size_t>(m_world_size), 0U);
    } catch (...) {
      local_failure = std::current_exception();
    }
    rethrowCollectivePreparationFailure(
        local_failure, "allgatherBytesBounded count-buffer preparation");
    if (MPI_Allgather(
            &local_count64, 1, MPI_UINT64_T,
            counts64.data(), 1, MPI_UINT64_T,
            MPI_COMM_WORLD) != MPI_SUCCESS) {
      throw std::runtime_error("allgatherBytesBounded count Allgather failed");
    }

    BoundedMpiTransferPlan plan;
    std::vector<std::uint8_t> gathered;
    std::vector<std::uint8_t> round_receive_buffer;
    local_failure = nullptr;
    try {
      std::vector<std::size_t> logical_counts(counts64.size(), 0U);
      for (std::size_t rank = 0; rank < counts64.size(); ++rank) {
        logical_counts[rank] = core::checkedIntegralNarrow<std::size_t>(
            counts64[rank], "allgatherBytesBounded received byte count");
      }
      plan = planBoundedMpiTransferRounds(
          logical_counts,
          static_cast<std::size_t>(std::numeric_limits<int>::max()),
          mpiTransportRoundLimitBytes());
      gathered.resize(plan.logical_total_count);
      std::size_t maximum_round_count = 0U;
      for (const auto& round : plan.rounds) {
        maximum_round_count = std::max(maximum_round_count, round.round_count);
      }
      round_receive_buffer.resize(maximum_round_count);
      injectMpiTestFault(*this, "allgather_post_count");
    } catch (...) {
      local_failure = std::current_exception();
    }
    rethrowCollectivePreparationFailure(
        local_failure, "allgatherBytesBounded payload preparation");

    for (const BoundedMpiRoundLayout& round : plan.rounds) {
      const std::size_t local_rank = static_cast<std::size_t>(m_world_rank);
      const int send_count = round.counts[local_rank];
      const std::uint8_t* send_pointer = send_count == 0
          ? nullptr
          : local_bytes.data() + round.logical_offsets[local_rank];
      if (MPI_Allgatherv(
              const_cast<std::uint8_t*>(send_pointer), send_count, MPI_BYTE,
              round_receive_buffer.empty() ? nullptr : round_receive_buffer.data(),
              round.counts.data(), round.displacements.data(), MPI_BYTE,
              MPI_COMM_WORLD) != MPI_SUCCESS) {
        throw std::runtime_error("allgatherBytesBounded bounded payload Allgatherv failed");
      }
      for (std::size_t rank = 0; rank < round.counts.size(); ++rank) {
        const std::size_t count = static_cast<std::size_t>(round.counts[rank]);
        if (count == 0U) {
          continue;
        }
        // The planner checked the logical prefix and every consumed peer offset
        // before any payload round, so this addition cannot overflow here.
        const std::size_t destination_offset =
            plan.logical_displacements[rank] + round.logical_offsets[rank];
        std::memcpy(
            gathered.data() + destination_offset,
            round_receive_buffer.data() +
                static_cast<std::size_t>(round.displacements[rank]),
            count);
      }
    }
    return gathered;
  }
#endif
  return {local_bytes.begin(), local_bytes.end()};
}


}  // namespace cosmosim::parallel
