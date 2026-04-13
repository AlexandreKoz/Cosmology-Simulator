#pragma once

#include <string_view>

namespace cosmosim::analysis {

class AnalysisModule {
 public:
  static std::string_view name();
  static std::string_view ownershipBoundary();
};

}  // namespace cosmosim::analysis
