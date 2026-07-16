#include "cosmosim/workflows/runtime_services.hpp"

#include <cassert>
#include <exception>
#include <stdexcept>
#include <string>

#include "cosmosim/core/profiling.hpp"
#include "cosmosim/parallel/distributed_memory.hpp"

int main() {
  const cosmosim::parallel::MpiContext mpi_context(false, 1, 0);
  cosmosim::core::ProfilerSession profiler(true);
  const cosmosim::workflows::RuntimeServices services{
      .mpi_context = mpi_context,
      .profiler = profiler,
      .deterministic_execution = true};
  const cosmosim::workflows::FailureCoordinator failures(services);

  assert(failures.failedRankCount(false) == 0U);
  failures.rethrowCollectiveFailure(nullptr, "healthy_phase");

  std::exception_ptr failure;
  try {
    throw std::invalid_argument("rank-local validation failure");
  } catch (...) {
    failure = std::current_exception();
  }
  bool rethrew_local_cause = false;
  try {
    failures.rethrowCollectiveFailure(failure, "failed_phase");
  } catch (const std::invalid_argument& error) {
    rethrew_local_cause =
        std::string(error.what()) == "rank-local validation failure";
  }
  assert(rethrew_local_cause);
  return 0;
}
