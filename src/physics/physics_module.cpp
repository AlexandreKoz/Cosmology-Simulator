#include "cosmosim/physics/physics_module.hpp"

#include "physics/internal/physics_boundary.hpp"

namespace cosmosim::physics {

std::string_view PhysicsModule::name() {
  return "physics";
}

std::string_view PhysicsModule::ownershipBoundary() {
  return internal::k_physics_boundary;
}

}  // namespace cosmosim::physics
