#include "cosmosim/utils/utils_module.hpp"

#include "utils/internal/utils_boundary.hpp"

namespace cosmosim::utils {

std::string_view UtilsModule::name() {
  return "utils";
}

std::string_view UtilsModule::ownershipBoundary() {
  return internal::k_utils_boundary;
}

}  // namespace cosmosim::utils
