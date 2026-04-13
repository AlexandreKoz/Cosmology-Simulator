#pragma once

#include <string_view>

namespace cosmosim::physics {

class PhysicsModule {
 public:
  static std::string_view name();
  static std::string_view ownershipBoundary();
};

}  // namespace cosmosim::physics
