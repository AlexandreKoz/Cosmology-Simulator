#pragma once

#include <string_view>

namespace cosmosim::hydro {

class HydroModule {
 public:
  static std::string_view name();
  static std::string_view ownershipBoundary();
};

}  // namespace cosmosim::hydro
