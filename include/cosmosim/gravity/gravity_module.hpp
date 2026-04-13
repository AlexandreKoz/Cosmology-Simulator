#pragma once

#include <string_view>

namespace cosmosim::gravity {

class GravityModule {
 public:
  static std::string_view name();
  static std::string_view ownershipBoundary();
};

}  // namespace cosmosim::gravity
