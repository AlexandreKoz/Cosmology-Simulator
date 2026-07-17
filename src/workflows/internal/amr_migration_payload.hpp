#pragma once

#include <vector>

#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/parallel/distributed_memory.hpp"

namespace cosmosim::workflows::internal {

[[nodiscard]] std::vector<parallel::AmrPatchPayloadRecord>
buildMigrationAmrPatchPayloadRecords(
    const core::SimulationState& state,
    int world_rank);
[[nodiscard]] std::vector<parallel::AmrPatchCellPayloadRecord>
buildMigrationAmrPatchCellPayloadRecords(
    const core::SimulationState& state,
    int world_rank);

}  // namespace cosmosim::workflows::internal
