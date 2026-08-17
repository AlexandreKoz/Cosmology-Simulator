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
  std::uint64_t local_particle_count = 0U;
  std::uint64_t local_cell_count = 0U;
  std::size_t tree_leaf_size = 16U;
  TreeMultipoleOrder multipole_order = TreeMultipoleOrder::kQuadrupole;
  PmGridShape pm_shape{};
  PmAssignmentScheme assignment_scheme = PmAssignmentScheme::kCic;
  core::PmDecompositionMode decomposition_mode = core::PmDecompositionMode::kSlab;
  std::uint32_t mpi_rank_count = 1U;
  int mpi_world_rank = 0;
  bool zoom_enabled = false;
  PmGridShape zoom_pm_shape{};
  bool periodic_tree_coordinates = true;
  bool indexed_target_coordinates = true;
  bool cuda_resident = false;
  std::uint64_t tree_exchange_batch_bytes = 4ULL * 1024ULL * 1024ULL;
  std::uint64_t pm_exchange_batch_bytes = 16ULL * 1024ULL * 1024ULL;
  std::uint64_t backend_unknown_reserve_bytes = 0U;
  double safety_margin_fraction = 0.0;
};

struct GravityMemoryEstimate {
  core::MemoryReport report;
  std::uint64_t known_peak_bytes = 0U;
  std::uint64_t external_backend_unknown_bytes = 0U;
  std::uint64_t budget_required_bytes = 0U;
  std::uint64_t estimated_tree_nodes = 0U;
  std::uint64_t pm_plan_owned_bytes = 0U;
};

struct DmoProcessMemoryPolicy {
  std::uint32_t mpi_rank_count = 1U;
  std::uint64_t scheduler_owned_bytes = 0U;
  std::uint64_t output_restart_overlap_bytes = 0U;
  std::uint64_t mpi_external_reserve_bytes = 0U;
  std::uint64_t fftw_external_reserve_bytes = 0U;
  std::uint64_t hdf5_external_reserve_bytes = 0U;
  std::uint64_t allocator_external_reserve_bytes = 0U;
  double safety_margin_fraction = 0.0;
};

struct DmoProcessMemoryEstimate {
  core::MemoryReport report;
  std::uint64_t known_owned_peak_bytes = 0U;
  std::uint64_t external_unknown_reserve_bytes = 0U;
  std::uint64_t modeled_subtotal_bytes = 0U;
  std::uint64_t safety_margin_bytes = 0U;
  std::uint64_t budget_required_bytes = 0U;
  std::uint64_t aggregate_required_bytes = 0U;
};

[[nodiscard]] GravityMemoryEstimate estimateGravityMemory(
    const GravityMemoryEstimateInput& input);

[[nodiscard]] DmoProcessMemoryEstimate estimateDmoProcessMemory(
    const core::MemoryReport& canonical_runtime_report,
    const GravityMemoryEstimate& gravity_estimate,
    const DmoProcessMemoryPolicy& policy);

void enforceGravityMemoryBudget(
    const GravityMemoryEstimate& estimate,
    std::uint64_t budget_bytes);

void enforceDmoProcessMemoryBudget(
    const DmoProcessMemoryEstimate& estimate,
    std::uint64_t budget_bytes);

}  // namespace cosmosim::gravity
