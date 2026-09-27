#pragma once

#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/parallel/distributed_memory.hpp"

#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
#include <mpi.h>
#endif

namespace cosmosim::parallel::internal {

#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
[[nodiscard]] inline bool queryActiveMpiWorld(int& world_size, int& world_rank) noexcept {
  world_size = 1;
  world_rank = 0;
  int initialized = 0;
  MPI_Initialized(&initialized);
  if (initialized == 0) {
    return false;
  }
  int finalized = 0;
  MPI_Finalized(&finalized);
  if (finalized != 0) {
    return false;
  }
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  return true;
}

inline void queryNodeLocalTopology(int fallback_world_rank, int& local_rank, int& local_size) {
  MPI_Comm local_comm = MPI_COMM_NULL;
  const int split_result = MPI_Comm_split_type(
      MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &local_comm);
  if (split_result != MPI_SUCCESS || local_comm == MPI_COMM_NULL) {
    throw std::runtime_error(
        "MPI_Comm_split_type(MPI_COMM_TYPE_SHARED) failed while determining node-local topology");
  }
  local_rank = fallback_world_rank;
  local_size = 1;
  const int rank_result = MPI_Comm_rank(local_comm, &local_rank);
  const int size_result = MPI_Comm_size(local_comm, &local_size);
  const int free_result = MPI_Comm_free(&local_comm);
  if (rank_result != MPI_SUCCESS || size_result != MPI_SUCCESS || free_result != MPI_SUCCESS) {
    throw std::runtime_error("MPI node-local communicator query failed");
  }
}
#endif

inline void injectMpiTestFault(const MpiContext& mpi_context, std::string_view phase) {
#if COSMOSIM_ENABLE_TESTS
  const char* raw = std::getenv("COSMOSIM_MPI_TEST_FAULT");
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
  if (configured_rank == mpi_context.worldRank() &&
      specification.substr(0U, separator) == phase) {
    throw std::runtime_error(
        "test-only injected MPI preparation failure at phase " + std::string(phase));
  }
#else
  static_cast<void>(mpi_context);
  static_cast<void>(phase);
#endif
}

[[nodiscard]] constexpr std::size_t ghostExchangeRecordBytes() {
  return sizeof(std::uint64_t) + sizeof(double) * 10U;
}

[[nodiscard]] inline bool laneIsPresentOrEmpty(
    std::size_t size, std::size_t expected) noexcept {
  return size == 0U || size == expected;
}

[[nodiscard]] inline double optionalLaneValue(
    const std::vector<double>& lane, std::size_t index) {
  return lane.empty() ? 0.0 : lane[index];
}

inline void commitOptionalGhostLane(
    std::vector<double>* destination,
    const std::vector<double>& source,
    std::size_t destination_size,
    std::size_t destination_index,
    std::size_t source_index) {
  if (source.empty()) {
    return;
  }
  if (source_index >= source.size()) {
    throw std::invalid_argument("ghost source optional lane does not cover committed row");
  }
  if (destination->empty()) {
    destination->assign(destination_size, 0.0);
  } else if (destination->size() != destination_size) {
    throw std::invalid_argument("ghost optional lane size does not match storage size");
  }
  (*destination)[destination_index] = source[source_index];
}

}  // namespace cosmosim::parallel::internal
