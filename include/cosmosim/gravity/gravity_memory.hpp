#pragma once

#include <cstddef>
#include <cstdint>

#include "cosmosim/core/config.hpp"
#include "cosmosim/core/memory_accounting.hpp"
#include "cosmosim/gravity/pm_solver.hpp"
#include "cosmosim/gravity/tree_gravity.hpp"

namespace cosmosim::gravity {

struct GravityMemoryEstimateInput {
  std::uint64_t local_source_count = 0U;
  std::uint64_t local_target_count = 0U;
  std::size_t tree_leaf_size = 16U;
  TreeMultipoleOrder multipole_order = TreeMultipoleOrder::kQuadrupole;
  PmGridShape pm_shape{};
  PmAssignmentScheme assignment_scheme = PmAssignmentScheme::kCic;
  core::PmDecompositionMode decomposition_mode = core::PmDecompositionMode::kSlab;
  std::uint32_t mpi_rank_count = 1U;
  bool zoom_enabled = false;
  PmGridShape zoom_pm_shape{};
  bool periodic_tree_coordinates = true;
  bool indexed_target_coordinates = true;
  bool cuda_resident = false;
  std::uint64_t tree_exchange_batch_bytes = 4ULL * 1024ULL * 1024ULL;
};

struct GravityMemoryEstimate {
  core::MemoryReport report;
  std::uint64_t known_peak_bytes = 0U;
  std::uint64_t external_backend_unknown_bytes = 0U;
  std::uint64_t estimated_tree_nodes = 0U;
};

[[nodiscard]] GravityMemoryEstimate estimateGravityMemory(
    const GravityMemoryEstimateInput& input);

void enforceGravityMemoryBudget(
    const GravityMemoryEstimate& estimate,
    std::uint64_t budget_bytes);

}  // namespace cosmosim::gravity
