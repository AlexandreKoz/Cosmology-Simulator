#include "cosmosim/amr/amr_module.hpp"

#include "amr/internal/amr_boundary.hpp"

namespace cosmosim::amr {

std::string_view AmrModule::name() {
  return "amr";
}

std::string_view AmrModule::ownershipBoundary() {
  return internal::k_amr_boundary;
}

}  // namespace cosmosim::amr
