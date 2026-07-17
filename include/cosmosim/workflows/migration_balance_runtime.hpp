#pragma once

#include <cstdint>
#include <span>

#include "cosmosim/core/config.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/core/time_integration.hpp"
#include "cosmosim/parallel/distributed_memory.hpp"
#include "cosmosim/workflows/runtime_services.hpp"

namespace cosmosim::workflows::internal {

// Sole workflow owner for initial placement, runtime rebalance, migration
// transport/commit, compaction, scheduler remap, and ownership-generation
// invalidation. It borrows the one process RuntimeServices bundle.
class MigrationBalanceRuntime {
 public:
  MigrationBalanceRuntime(
      const core::SimulationConfig& config,
      const RuntimeServices& services) noexcept;

  void initializeOwnership(core::SimulationState& state) const;

  [[nodiscard]] parallel::LocalOwnershipIdentitySummary reduceIdentity(
      const core::SimulationState& state) const;

  [[nodiscard]] bool rebalance(
      core::SimulationState& state,
      core::HierarchicalTimeBinScheduler& particle_scheduler,
      core::HierarchicalTimeBinScheduler& gas_cell_scheduler,
      const parallel::DecompositionRuntimeMeasurements& measurements,
      std::span<const std::uint32_t> active_particle_indices,
      std::span<const std::uint64_t> expected_global_particle_ids,
      std::uint64_t step_index) const;

 private:
  const core::SimulationConfig& m_config;
  const RuntimeServices& m_services;
};

}  // namespace cosmosim::workflows::internal
