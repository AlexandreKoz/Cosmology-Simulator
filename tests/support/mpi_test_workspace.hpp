#pragma once

#include <cstdint>
#include <string_view>

#include <mpi.h>

#include "test_temp_workspace.hpp"

namespace cosmosim::test_support {

// Creates one collision-resistant filesystem workspace shared by every rank in
// `communicator`. Rank zero creates the run token and owns eventual cleanup.
// Callers must ensure all ranks are finished with files in the workspace (for
// example with MPI_Barrier) before the owning workspace is destroyed.
[[nodiscard]] inline TestTempWorkspace createMpiSharedWorkspace(
    std::string_view test_stem,
    MPI_Comm communicator = MPI_COMM_WORLD,
    int root_rank = 0) {
  int world_rank = 0;
  MPI_Comm_rank(communicator, &world_rank);
  std::uint64_t token = world_rank == root_rank
      ? TestTempWorkspace::uniqueRunToken()
      : 0U;
  MPI_Bcast(&token, 1, MPI_UINT64_T, root_rank, communicator);
  return TestTempWorkspace::createMpiShared(
      test_stem, token, world_rank == root_rank);
}

}  // namespace cosmosim::test_support
