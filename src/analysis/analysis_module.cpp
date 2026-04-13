#include "cosmosim/analysis/analysis_module.hpp"

#include "analysis/internal/analysis_boundary.hpp"

namespace cosmosim::analysis {

std::string_view AnalysisModule::name() {
  return "analysis";
}

std::string_view AnalysisModule::ownershipBoundary() {
  return internal::k_analysis_boundary;
}

}  // namespace cosmosim::analysis
