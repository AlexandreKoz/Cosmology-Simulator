#include "io/internal/ic_failure_protocol.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <stdexcept>

#include "io/internal/ic_record_codec.hpp"
#include "io/internal/ic_mpi_collectives.hpp"

#if COSMOSIM_ENABLE_MPI
#include <mpi.h>
#endif

namespace cosmosim::io::failure_protocol_internal {

#if COSMOSIM_ENABLE_HDF5 && COSMOSIM_ENABLE_MPI
using mpi_collective_internal::mpiAllreduce;
using mpi_collective_internal::mpiAlltoall;
using mpi_collective_internal::mpiAlltoallv;
using mpi_collective_internal::mpiBcast;

namespace {

thread_local std::uint64_t* g_collective_phase_counter = nullptr;

struct AlltoallLayout {
  std::vector<int> send_counts;
  std::vector<int> receive_counts;
  std::vector<int> send_displacements;
  std::vector<int> receive_displacements;
  std::uint64_t send_total = 0U;
  std::uint64_t receive_total = 0U;
};

}  // namespace

std::string collectiveFailureMessage(
    const parallel::MpiContext& mpi_context,
    const std::exception_ptr& local_failure,
    std::string_view phase) {
  const int local_failed = local_failure ? 1 : 0;
  int failed_count = 0;
  mpiAllreduce(
      &local_failed, &failed_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if (failed_count == 0) {
    return {};
  }
  int candidate = local_failure ? mpi_context.worldRank()
                                : mpi_context.worldSize();
  int failure_rank = mpi_context.worldSize();
  mpiAllreduce(
      &candidate, &failure_rank, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

  static constexpr std::size_t kMaximumMessageBytes = 4095U;
  std::array<char, kMaximumMessageBytes + 1U> buffer{};
  std::uint32_t length = 0U;
  if (mpi_context.worldRank() == failure_rank) {
    const char* message = "unknown non-standard exception";
    try {
      std::rethrow_exception(local_failure);
    } catch (const std::exception& error) {
      message = error.what();
    } catch (...) {
    }
    const std::size_t message_length = std::char_traits<char>::length(message);
    length = static_cast<std::uint32_t>(
        std::min(message_length, kMaximumMessageBytes));
    std::copy_n(message, length, buffer.data());
  }
  mpiBcast(&length, 1, MPI_UINT32_T, failure_rank, MPI_COMM_WORLD);
  mpiBcast(
      buffer.data(), static_cast<int>(buffer.size()), MPI_CHAR, failure_rank,
      MPI_COMM_WORLD);
  return std::string(phase) + " failed on rank " +
      std::to_string(failure_rank) + ": " +
      std::string(buffer.data(), buffer.data() + length);
}

void incrementCollectivePhaseCounter() noexcept {
  if (g_collective_phase_counter == nullptr) {
    return;
  }
  if (*g_collective_phase_counter <
      std::numeric_limits<std::uint64_t>::max()) {
    ++(*g_collective_phase_counter);
  }
}

CollectivePhaseCounterScope::CollectivePhaseCounterScope(
    std::uint64_t& counter) noexcept
    : m_previous(g_collective_phase_counter) {
  g_collective_phase_counter = &counter;
}

CollectivePhaseCounterScope::~CollectivePhaseCounterScope() {
  g_collective_phase_counter = m_previous;
}

void injectIcTestFault(
    const parallel::MpiContext& mpi_context,
    std::string_view phase) {
#if COSMOSIM_ENABLE_TESTS
  const char* raw = std::getenv("COSMOSIM_IC_TEST_FAULT");
  if (raw == nullptr || *raw == '\0') {
    return;
  }
  const std::string specification(raw);
  const std::size_t separator = specification.rfind(':');
  if (separator == std::string::npos) {
    return;
  }
  const std::string configured_phase = specification.substr(0, separator);
  int configured_rank = -1;
  try {
    configured_rank = std::stoi(specification.substr(separator + 1U));
  } catch (...) {
    return;
  }
  if (configured_phase == phase &&
      configured_rank == mpi_context.worldRank()) {
    throw std::runtime_error(
        "test-only injected IC failure at phase " + std::string(phase));
  }
#else
  static_cast<void>(mpi_context);
  static_cast<void>(phase);
#endif
}

void mutateIcTestRoute(
    const parallel::MpiContext& mpi_context,
    std::vector<std::vector<std::uint8_t>>& per_rank) {
#if COSMOSIM_ENABLE_TESTS
  static bool mutation_applied = false;
  if (mutation_applied) {
    return;
  }
  const char* raw = std::getenv("COSMOSIM_IC_TEST_ROUTE_MUTATION");
  if (raw == nullptr || *raw == '\0') {
    return;
  }
  const std::string specification(raw);
  const std::size_t separator = specification.rfind(':');
  if (separator == std::string::npos) {
    return;
  }
  int configured_rank = -1;
  try {
    configured_rank = std::stoi(specification.substr(separator + 1U));
  } catch (...) {
    return;
  }
  if (configured_rank != mpi_context.worldRank()) {
    return;
  }
  const std::string operation = specification.substr(0U, separator);
  auto bucket = std::find_if(
      per_rank.begin(), per_rank.end(), [](const auto& candidate) {
        return !candidate.empty();
      });
  if (bucket == per_rank.end()) {
    return;
  }
  const auto [record_begin, record_size] = internal::lastIcWireRecordSpan(*bucket);
  if (record_size == 0U) return;
  if (operation == "drop") {
    bucket->resize(record_begin);
  } else if (operation == "duplicate") {
    const std::vector<std::uint8_t> copy(
        bucket->begin() + static_cast<std::ptrdiff_t>(record_begin), bucket->end());
    bucket->insert(bucket->end(), copy.begin(), copy.end());
  } else {
    return;
  }
  mutation_applied = true;
#else
  static_cast<void>(mpi_context);
  static_cast<void>(per_rank);
#endif
}

std::string broadcastRootString(
    const parallel::MpiContext& mpi_context,
    std::string root_value) {
  std::uint64_t length = mpi_context.isRoot() ? root_value.size() : 0U;
  mpiBcast(&length, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
  runCollectivePhaseVoid(
      mpi_context, "root-string receive allocation", [&]() {
        if (length > std::numeric_limits<std::size_t>::max()) {
          throw std::overflow_error("root string exceeds local size_t");
        }
        root_value.resize(static_cast<std::size_t>(length));
      });
  static constexpr std::uint64_t kBroadcastChunk = 64U * 1024U * 1024U;
  for (std::uint64_t offset = 0U; offset < length;
       offset += kBroadcastChunk) {
    const int chunk = static_cast<int>(
        std::min(kBroadcastChunk, length - offset));
    mpiBcast(
        root_value.data() + static_cast<std::size_t>(offset), chunk, MPI_CHAR,
        0, MPI_COMM_WORLD);
  }
  return root_value;
}

std::vector<std::uint8_t> alltoallBytes(
    const parallel::MpiContext& mpi_context,
    const std::vector<std::vector<std::uint8_t>>& per_rank,
    std::uint64_t& sent_bytes,
    std::uint64_t& received_bytes,
    std::uint64_t* exchange_peak_bytes) {
  const int world_size = mpi_context.worldSize();
  AlltoallLayout layout = runCollectivePhase<AlltoallLayout>(
      mpi_context, "IC all-to-all send-layout preparation", [&]() {
        if (per_rank.size() != static_cast<std::size_t>(world_size)) {
          throw std::invalid_argument(
              "IC all-to-all requires one send bucket per rank");
        }
        AlltoallLayout prepared;
        prepared.send_counts.assign(static_cast<std::size_t>(world_size), 0);
        prepared.receive_counts.assign(
            static_cast<std::size_t>(world_size), 0);
        prepared.send_displacements.assign(
            static_cast<std::size_t>(world_size), 0);
        prepared.receive_displacements.assign(
            static_cast<std::size_t>(world_size), 0);
        for (int rank = 0; rank < world_size; ++rank) {
          const std::size_t size =
              per_rank[static_cast<std::size_t>(rank)].size();
          if (size > static_cast<std::size_t>(
                         std::numeric_limits<int>::max()) ||
              prepared.send_total >
                  static_cast<std::uint64_t>(
                      std::numeric_limits<int>::max()) - size) {
            throw std::overflow_error(
                "bounded IC send payload exceeds MPI int count/displacement");
          }
          prepared.send_displacements[static_cast<std::size_t>(rank)] =
              static_cast<int>(prepared.send_total);
          prepared.send_counts[static_cast<std::size_t>(rank)] =
              static_cast<int>(size);
          prepared.send_total += size;
        }
        return prepared;
      });

  mpiAlltoall(
      layout.send_counts.data(), 1, MPI_INT, layout.receive_counts.data(), 1,
      MPI_INT, MPI_COMM_WORLD);

  struct ExchangeBuffers {
    std::vector<std::uint8_t> send;
    std::vector<std::uint8_t> receive;
  };
  ExchangeBuffers buffers = runCollectivePhase<ExchangeBuffers>(
      mpi_context, "IC all-to-all receive-layout and buffer preparation", [&]() {
        for (int rank = 0; rank < world_size; ++rank) {
          const int count =
              layout.receive_counts[static_cast<std::size_t>(rank)];
          if (count < 0 ||
              layout.receive_total >
                  static_cast<std::uint64_t>(
                      std::numeric_limits<int>::max()) -
                      static_cast<std::uint64_t>(count)) {
            throw std::overflow_error(
                "bounded IC receive payload exceeds MPI int count/displacement");
          }
          layout.receive_displacements[static_cast<std::size_t>(rank)] =
              static_cast<int>(layout.receive_total);
          layout.receive_total += static_cast<std::uint64_t>(count);
        }
        ExchangeBuffers prepared;
        prepared.send.resize(static_cast<std::size_t>(layout.send_total));
        prepared.receive.resize(
            static_cast<std::size_t>(layout.receive_total));
        for (int rank = 0; rank < world_size; ++rank) {
          const auto& bucket = per_rank[static_cast<std::size_t>(rank)];
          std::copy(
              bucket.begin(), bucket.end(),
              prepared.send.begin() +
                  layout.send_displacements[static_cast<std::size_t>(rank)]);
        }
        return prepared;
      });

  mpiAlltoallv(
      buffers.send.data(), layout.send_counts.data(),
      layout.send_displacements.data(), MPI_BYTE, buffers.receive.data(),
      layout.receive_counts.data(), layout.receive_displacements.data(),
      MPI_BYTE, MPI_COMM_WORLD);
  sent_bytes += layout.send_total;
  received_bytes += layout.receive_total;
  if (exchange_peak_bytes != nullptr) {
    const std::uint64_t layout_bytes =
        static_cast<std::uint64_t>(layout.send_counts.capacity() +
                                   layout.receive_counts.capacity() +
                                   layout.send_displacements.capacity() +
                                   layout.receive_displacements.capacity()) *
        sizeof(int);
    *exchange_peak_bytes = layout_bytes +
        static_cast<std::uint64_t>(buffers.send.capacity()) +
        static_cast<std::uint64_t>(buffers.receive.capacity());
  }
  return std::move(buffers.receive);
}

#endif

}  // namespace cosmosim::io::failure_protocol_internal
