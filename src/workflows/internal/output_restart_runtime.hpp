#pragma once

#include "cosmosim/workflows/output_restart_runtime.hpp"

namespace cosmosim::core {
class ProfilerSession;
struct FrozenConfig;
struct SimulationConfig;
struct SimulationState;
}

namespace cosmosim::io {
struct GravityForceCachePersistentState;
struct OutputCadencePersistentState;
struct RestartReadResult;
}

namespace cosmosim::parallel {
struct DistributedExecutionTopology;
class MpiContext;
}

namespace cosmosim::workflows::internal {

void validateRestartResumeTopologyOrThrow(
    const io::RestartReadResult& restart,
    const core::SimulationConfig& config,
    const core::FrozenConfig& frozen_config,
    const parallel::MpiContext& mpi_context);

}  // namespace cosmosim::workflows::internal
