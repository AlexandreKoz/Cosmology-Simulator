#pragma once

#include <string_view>

namespace cosmosim::io {

class IoModule {
 public:
  static std::string_view name();
  static std::string_view ownershipBoundary();
};

}  // namespace cosmosim::io
