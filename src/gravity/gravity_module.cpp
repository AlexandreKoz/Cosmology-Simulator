#include "cosmosim/gravity/gravity_module.hpp"

#include "gravity/internal/gravity_boundary.hpp"

namespace cosmosim::gravity {

std::string_view GravityModule::name() {
  return "gravity";
}

std::string_view GravityModule::ownershipBoundary() {
  return internal::k_gravity_boundary;
}

}  // namespace cosmosim::gravity
