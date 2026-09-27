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
    int local_rank,
    int configured_gpu_devices,
    bool cuda_runtime_available,
    int visible_device_count) {
  if (local_rank < 0) {
    throw std::invalid_argument("local_rank must be non-negative");
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
  assignment.assigned_device_index = local_rank % configured_gpu_devices;
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
  topology.local_rank = mpi_context.localRank();
  topology.mpi_enabled = mpi_context.isEnabled();
  topology.pm_decomposition_mode = std::move(pm_decomposition_mode);
  topology.pm_slab = makePmSlabLayout(global_nx, global_ny, global_nz, mpi_context.worldSize(), mpi_context.worldRank());
  topology.device_assignment =
      selectRankDeviceAssignment(mpi_context.localRank(), configured_gpu_devices, cuda_runtime_available, visible_device_count);
  if (!topology.device_assignment.isValid()) {
    throw std::runtime_error("constructed distributed execution topology is invalid");
  }
  return topology;
}



}  // namespace cosmosim::parallel
