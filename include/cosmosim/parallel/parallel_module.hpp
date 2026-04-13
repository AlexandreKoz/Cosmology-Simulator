#pragma once

#include <string_view>

namespace cosmosim::parallel {

class ParallelModule {
 public:
  static std::string_view name();
  static std::string_view ownershipBoundary();
};

}  // namespace cosmosim::parallel
