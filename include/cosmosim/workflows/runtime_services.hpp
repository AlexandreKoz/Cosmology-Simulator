#pragma once

#include <cstdint>
#include <exception>
#include <string_view>

namespace cosmosim::core {
class ProfilerSession;
}

namespace cosmosim::parallel {
class MpiContext;
}

namespace cosmosim::workflows {

// Process-lifetime dependencies are created by the composition root and
// borrowed by runtime components. Services must never create replacement MPI
// contexts or profiler authorities internally.
struct RuntimeServices {
  const parallel::MpiContext& mpi_context;
  core::ProfilerSession& profiler;
  bool deterministic_execution = true;
};

// Collective-safe phase gate. Every rank must call the gate for a given phase;
// no rank may enter the following collective phase after a peer reports a
// local failure.
class FailureCoordinator {
 public:
  explicit FailureCoordinator(const RuntimeServices& services) noexcept;

  [[nodiscard]] std::uint64_t failedRankCount(bool local_failed) const;
  void rethrowCollectiveFailure(
      const std::exception_ptr& local_failure,
      std::string_view phase_name) const;

 private:
  const RuntimeServices& m_services;
};

}  // namespace cosmosim::workflows
