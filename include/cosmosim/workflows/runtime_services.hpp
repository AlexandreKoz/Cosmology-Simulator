#pragma once

#include <cstdint>
#include <exception>
#include <string_view>

namespace cosmosim::core {
class MemoryGovernor;
class ProfilerSession;
struct MemoryReport;
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
  // Optional only for standalone/unit composition. The production reference
  // workflow always supplies the one process-level governor authority.
  core::MemoryGovernor* memory_governor = nullptr;
  bool deterministic_execution = true;
};

// Collective-safe phase gate. Every rank must call the gate for a given phase;
// no rank may enter the following collective phase after a peer reports a
// local failure.
// Reduce the current memory snapshot across ranks at an explicit workflow
// boundary. Optional OS observations are reported only when every rank exposes
// that metric; rank maximum remains the safety-relevant value.
void attachDistributedMemoryTelemetry(
    core::MemoryReport& report,
    const RuntimeServices& services);

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
