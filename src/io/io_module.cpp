#include "cosmosim/io/io_module.hpp"

#include "io/internal/io_boundary.hpp"

namespace cosmosim::io {

std::string_view IoModule::name() {
  return "io";
}

std::string_view IoModule::ownershipBoundary() {
  return internal::k_io_boundary;
}

}  // namespace cosmosim::io
