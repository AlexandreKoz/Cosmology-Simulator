#pragma once

#include <string_view>

namespace cosmosim::utils {

class UtilsModule {
 public:
  static std::string_view name();
  static std::string_view ownershipBoundary();
};

}  // namespace cosmosim::utils
