#include "cosmosim/workflows/runtime_services.hpp"

#include <cstdint>
#include <stdexcept>
#include <string>

#include "cosmosim/parallel/distributed_memory.hpp"

namespace cosmosim::workflows {

FailureCoordinator::FailureCoordinator(const RuntimeServices& services) noexcept
    : m_services(services) {}

std::uint64_t FailureCoordinator::failedRankCount(bool local_failed) const {
  return m_services.mpi_context.allreduceSumUint64(local_failed ? 1ULL : 0ULL);
}

void FailureCoordinator::rethrowCollectiveFailure(
    const std::exception_ptr& local_failure,
    std::string_view phase_name) const {
  const std::uint64_t failed_rank_count =
      failedRankCount(local_failure != nullptr);
  if (failed_rank_count == 0U) {
    return;
  }
  if (local_failure != nullptr) {
    std::rethrow_exception(local_failure);
  }
  throw std::runtime_error(
      "collective phase '" + std::string(phase_name) + "' aborted because " +
      std::to_string(failed_rank_count) + " peer rank(s) reported failure");
}

}  // namespace cosmosim::workflows
