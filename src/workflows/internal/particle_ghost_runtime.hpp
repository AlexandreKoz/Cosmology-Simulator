#pragma once

#include <cstddef>
#include <cstdint>
#include <string_view>

#include "cosmosim/core/time_integration.hpp"
#include "cosmosim/parallel/distributed_memory.hpp"

namespace cosmosim::workflows::internal {

struct SolverGhostRefreshReport {
  std::uint64_t sent_bytes = 0;
  std::uint64_t received_bytes = 0;
  std::size_t committed_slots = 0;
};

[[nodiscard]] SolverGhostRefreshReport refreshParticleGhostsForSolver(
    core::StepContext& context,
    const parallel::MpiContext& mpi_context,
    std::string_view subsystem_name,
    parallel::GhostCacheLifecycle* lifecycle = nullptr);

}  // namespace cosmosim::workflows::internal
