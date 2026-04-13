#include "cosmosim/hydro/hydro_module.hpp"

#include "hydro/internal/hydro_boundary.hpp"

namespace cosmosim::hydro {

std::string_view HydroModule::name() {
  return "hydro";
}

std::string_view HydroModule::ownershipBoundary() {
  return internal::k_hydro_boundary;
}

}  // namespace cosmosim::hydro
