#include "cosmosim/parallel/parallel_module.hpp"

#include "parallel/internal/parallel_boundary.hpp"

namespace cosmosim::parallel {

std::string_view ParallelModule::name() {
  return "parallel";
}

std::string_view ParallelModule::ownershipBoundary() {
  return internal::k_parallel_boundary;
}

}  // namespace cosmosim::parallel
