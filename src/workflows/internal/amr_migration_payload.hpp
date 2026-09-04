#pragma once

#include <span>
#include <vector>

#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/parallel/distributed_memory.hpp"

namespace cosmosim::workflows::internal {

[[nodiscard]] std::vector<parallel::AmrPatchPayloadRecord>
buildMigrationAmrPatchPayloadRecords(
    const core::SimulationState& state,
    int world_rank);
void fillMigrationAmrPatchBoundaryCellPayloadChunk(
    const core::SimulationState& state,
    int world_rank,
    std::span<const parallel::AmrPatchBoundaryCellRequest> requests,
    std::uint64_t first_record,
    std::size_t max_records,
    std::vector<parallel::AmrPatchCellPayloadRecord>& output);

}  // namespace cosmosim::workflows::internal
