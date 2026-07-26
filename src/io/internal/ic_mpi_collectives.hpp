#pragma once

#include <cstdint>

#include "cosmosim/core/build_config.hpp"

#if COSMOSIM_ENABLE_MPI
#include <mpi.h>
#endif

namespace cosmosim::io::mpi_collective_internal {

#if COSMOSIM_ENABLE_HDF5 && COSMOSIM_ENABLE_MPI

struct MpiCollectiveCallCounts {
  std::uint64_t total = 0U;
  std::uint64_t routing_total = 0U;
  std::uint64_t allreduce = 0U;
  std::uint64_t bcast = 0U;
  std::uint64_t gather = 0U;
  std::uint64_t gatherv = 0U;
  std::uint64_t alltoall = 0U;
  std::uint64_t alltoallv = 0U;
};

class MpiCollectiveCounterScope {
 public:
  explicit MpiCollectiveCounterScope(MpiCollectiveCallCounts& counts) noexcept;
  MpiCollectiveCounterScope(const MpiCollectiveCounterScope&) = delete;
  MpiCollectiveCounterScope& operator=(const MpiCollectiveCounterScope&) = delete;
  ~MpiCollectiveCounterScope();

 private:
  MpiCollectiveCallCounts* m_previous = nullptr;
};

class RoutingMpiCollectiveScope {
 public:
  RoutingMpiCollectiveScope() noexcept;
  RoutingMpiCollectiveScope(const RoutingMpiCollectiveScope&) = delete;
  RoutingMpiCollectiveScope& operator=(const RoutingMpiCollectiveScope&) = delete;
  ~RoutingMpiCollectiveScope();
};

int mpiAllreduce(
    const void* send_buffer,
    void* receive_buffer,
    int count,
    MPI_Datatype datatype,
    MPI_Op operation,
    MPI_Comm communicator) noexcept;

int mpiBcast(
    void* buffer,
    int count,
    MPI_Datatype datatype,
    int root,
    MPI_Comm communicator) noexcept;

int mpiGather(
    const void* send_buffer,
    int send_count,
    MPI_Datatype send_datatype,
    void* receive_buffer,
    int receive_count,
    MPI_Datatype receive_datatype,
    int root,
    MPI_Comm communicator) noexcept;

int mpiGatherv(
    const void* send_buffer,
    int send_count,
    MPI_Datatype send_datatype,
    void* receive_buffer,
    const int* receive_counts,
    const int* displacements,
    MPI_Datatype receive_datatype,
    int root,
    MPI_Comm communicator) noexcept;

int mpiAlltoall(
    const void* send_buffer,
    int send_count,
    MPI_Datatype send_datatype,
    void* receive_buffer,
    int receive_count,
    MPI_Datatype receive_datatype,
    MPI_Comm communicator) noexcept;

int mpiAlltoallv(
    const void* send_buffer,
    const int* send_counts,
    const int* send_displacements,
    MPI_Datatype send_datatype,
    void* receive_buffer,
    const int* receive_counts,
    const int* receive_displacements,
    MPI_Datatype receive_datatype,
    MPI_Comm communicator) noexcept;

#endif

}  // namespace cosmosim::io::mpi_collective_internal
