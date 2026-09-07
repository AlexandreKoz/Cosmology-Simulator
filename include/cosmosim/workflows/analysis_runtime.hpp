#pragma once

#include <memory>
#include <string>
#include <vector>

#include "cosmosim/core/config.hpp"
#include "cosmosim/core/time_integration.hpp"
#include "cosmosim/workflows/runtime_resources.hpp"

namespace cosmosim::workflows {

struct RuntimeServices;

class AnalysisRuntime {
 public:
  virtual ~AnalysisRuntime() = default;

  virtual void audit(AnalysisStageView& view) = 0;
  virtual void executeDiagnostics(AnalysisStageView& view) = 0;
  // Required-only owner preflight. Optional products have separate, deferrable
  // physical leases and must never turn this callback into a hard stage gate.
  [[nodiscard]] virtual std::uint64_t estimateRequiredIncrementalBytes(
      const core::SimulationState& state, std::uint64_t completed_step) const = 0;
  // Best-effort optional cadence is coalesced, not historical replay.
  // Checkpoint metadata reports pending work; the pending queue is deliberately
  // not restart truth. A resumed run starts a new optional cadence epoch.
  [[nodiscard]] virtual std::string optionalCadenceProvenance() const = 0;
  virtual void finalizePending(
      std::uint64_t completed_step, double time_code, double scale_factor) = 0;

};

[[nodiscard]] std::unique_ptr<AnalysisRuntime> makeAnalysisRuntime(
    const core::SimulationConfig& config,
    std::vector<std::string>& stage_sequence,
    const RuntimeServices& services);

}  // namespace cosmosim::workflows
