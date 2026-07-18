#include "cosmosim/workflows/analysis_runtime.hpp"

class ExternalAnalysis final : public cosmosim::workflows::AnalysisRuntime {
 public:
  void audit(cosmosim::workflows::AnalysisStageView& view) override {
    auto& context = stageContext(view);
    context.state.particles.mass_code[0] = 0.0;
  }

  void executeDiagnostics(
      cosmosim::workflows::AnalysisStageView&) override {}
};
