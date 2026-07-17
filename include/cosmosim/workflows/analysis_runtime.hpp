#pragma once

#include <memory>
#include <string>
#include <vector>

#include "cosmosim/core/config.hpp"
#include "cosmosim/core/time_integration.hpp"
#include "cosmosim/workflows/runtime_resources.hpp"

namespace cosmosim::workflows {

class AnalysisRuntime {
 public:
  virtual ~AnalysisRuntime() = default;

  virtual void audit(AnalysisStageView& view) = 0;
  virtual void executeDiagnostics(AnalysisStageView& view) = 0;

 protected:
  [[nodiscard]] static core::StepContext& stageContext(
      AnalysisStageView& view) {
    return view.ownerContext();
  }
};

[[nodiscard]] std::unique_ptr<AnalysisRuntime> makeAnalysisRuntime(
    const core::SimulationConfig& config,
    std::vector<std::string>& stage_sequence);

}  // namespace cosmosim::workflows
