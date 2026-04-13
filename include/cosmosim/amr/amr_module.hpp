#pragma once

#include <string_view>

namespace cosmosim::amr {

class AmrModule {
 public:
  static std::string_view name();
  static std::string_view ownershipBoundary();
};

}  // namespace cosmosim::amr
