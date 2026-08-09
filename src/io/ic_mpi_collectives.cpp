#include "io/internal/ic_mpi_collectives.hpp"

#include <cstdio>
#include <exception>
#include <limits>

namespace cosmosim::io::mpi_collective_internal {

#if COSMOSIM_ENABLE_HDF5 && COSMOSIM_ENABLE_MPI
namespace {

thread_local MpiCollectiveCallCounts* g_counts = nullptr;
thread_local std::uint32_t g_routing_scope_depth = 0U;

enum class CollectiveKind {
  kAllreduce,
  kBcast,
  kGather,
  kGatherv,
  kAlltoall,
  kAlltoallv,
};

void saturatingIncrement(std::uint64_t& value) noexcept {
  if (value < std::numeric_limits<std::uint64_t>::max()) {
    ++value;
  }
}

[[noreturn]] void abortForMpiFailure(
    MPI_Comm communicator,
    int status,
    const char* operation) noexcept {
  char error_string[MPI_MAX_ERROR_STRING]{};
  int error_length = 0;
  const int describe_status =
      MPI_Error_string(status, error_string, &error_length);
  if (describe_status == MPI_SUCCESS && error_length > 0) {
    std::fprintf(
        stderr, "CosmoSim fatal MPI I/O error in %s: %.*s (status=%d)\n",
        operation, error_length, error_string, status);
  } else {
    std::fprintf(
        stderr, "CosmoSim fatal MPI I/O error in %s (status=%d)\n",
        operation, status);
  }
  std::fflush(stderr);
  MPI_Abort(communicator, status);
  std::terminate();
}

int requireMpiSuccess(
    MPI_Comm communicator,
    int status,
    const char* operation) noexcept {
  if (status != MPI_SUCCESS) {
    abortForMpiFailure(communicator, status, operation);
  }
  return status;
}

void recordCall(CollectiveKind kind) noexcept {
  if (g_counts == nullptr) {
    return;
  }
  saturatingIncrement(g_counts->total);
  if (g_routing_scope_depth > 0U) {
    saturatingIncrement(g_counts->routing_total);
  }
  switch (kind) {
    case CollectiveKind::kAllreduce:
      saturatingIncrement(g_counts->allreduce);
      break;
    case CollectiveKind::kBcast:
      saturatingIncrement(g_counts->bcast);
      break;
    case CollectiveKind::kGather:
      saturatingIncrement(g_counts->gather);
      break;
    case CollectiveKind::kGatherv:
      saturatingIncrement(g_counts->gatherv);
      break;
    case CollectiveKind::kAlltoall:
      saturatingIncrement(g_counts->alltoall);
      break;
    case CollectiveKind::kAlltoallv:
      saturatingIncrement(g_counts->alltoallv);
      break;
  }
}

}  // namespace

void configureMpiErrorsReturn(MPI_Comm communicator) noexcept {
  const int status = MPI_Comm_set_errhandler(communicator, MPI_ERRORS_RETURN);
  requireMpiSuccess(communicator, status, "MPI_Comm_set_errhandler");
}

MpiCollectiveCounterScope::MpiCollectiveCounterScope(
    MpiCollectiveCallCounts& counts) noexcept
    : m_previous(g_counts) {
  g_counts = &counts;
}

MpiCollectiveCounterScope::~MpiCollectiveCounterScope() {
  g_counts = m_previous;
}

RoutingMpiCollectiveScope::RoutingMpiCollectiveScope() noexcept {
  if (g_routing_scope_depth < std::numeric_limits<std::uint32_t>::max()) {
    ++g_routing_scope_depth;
  }
}

RoutingMpiCollectiveScope::~RoutingMpiCollectiveScope() {
  if (g_routing_scope_depth > 0U) {
    --g_routing_scope_depth;
  }
}

int mpiAllreduce(
    const void* send_buffer,
    void* receive_buffer,
    int count,
    MPI_Datatype datatype,
    MPI_Op operation,
    MPI_Comm communicator) noexcept {
  recordCall(CollectiveKind::kAllreduce);
  return requireMpiSuccess(
      communicator,
      MPI_Allreduce(
          send_buffer, receive_buffer, count, datatype, operation, communicator),
      "MPI_Allreduce");
}

int mpiBcast(
    void* buffer,
    int count,
    MPI_Datatype datatype,
    int root,
    MPI_Comm communicator) noexcept {
  recordCall(CollectiveKind::kBcast);
  return requireMpiSuccess(
      communicator, MPI_Bcast(buffer, count, datatype, root, communicator),
      "MPI_Bcast");
}

int mpiGather(
    const void* send_buffer,
    int send_count,
    MPI_Datatype send_datatype,
    void* receive_buffer,
    int receive_count,
    MPI_Datatype receive_datatype,
    int root,
    MPI_Comm communicator) noexcept {
  recordCall(CollectiveKind::kGather);
  return requireMpiSuccess(
      communicator,
      MPI_Gather(
          send_buffer, send_count, send_datatype, receive_buffer, receive_count,
          receive_datatype, root, communicator),
      "MPI_Gather");
}

int mpiGatherv(
    const void* send_buffer,
    int send_count,
    MPI_Datatype send_datatype,
    void* receive_buffer,
    const int* receive_counts,
    const int* displacements,
    MPI_Datatype receive_datatype,
    int root,
    MPI_Comm communicator) noexcept {
  recordCall(CollectiveKind::kGatherv);
  return requireMpiSuccess(
      communicator,
      MPI_Gatherv(
          send_buffer, send_count, send_datatype, receive_buffer, receive_counts,
          displacements, receive_datatype, root, communicator),
      "MPI_Gatherv");
}

int mpiAlltoall(
    const void* send_buffer,
    int send_count,
    MPI_Datatype send_datatype,
    void* receive_buffer,
    int receive_count,
    MPI_Datatype receive_datatype,
    MPI_Comm communicator) noexcept {
  recordCall(CollectiveKind::kAlltoall);
  return requireMpiSuccess(
      communicator,
      MPI_Alltoall(
          send_buffer, send_count, send_datatype, receive_buffer, receive_count,
          receive_datatype, communicator),
      "MPI_Alltoall");
}

int mpiAlltoallv(
    const void* send_buffer,
    const int* send_counts,
    const int* send_displacements,
    MPI_Datatype send_datatype,
    void* receive_buffer,
    const int* receive_counts,
    const int* receive_displacements,
    MPI_Datatype receive_datatype,
    MPI_Comm communicator) noexcept {
  recordCall(CollectiveKind::kAlltoallv);
  return requireMpiSuccess(
      communicator,
      MPI_Alltoallv(
          send_buffer, send_counts, send_displacements, send_datatype,
          receive_buffer, receive_counts, receive_displacements,
          receive_datatype, communicator),
      "MPI_Alltoallv");
}

#endif

}  // namespace cosmosim::io::mpi_collective_internal
