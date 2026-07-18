#include "cosmosim/workflows/runtime_resources.hpp"

void escape(cosmosim::workflows::AnalysisStageView& view) {
  (void)view.ownerContext();
}
