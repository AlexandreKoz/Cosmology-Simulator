#pragma once

#include <memory>

namespace cosmosim::core { class MemoryGovernor; }
namespace cosmosim::parallel { class MpiContext; }
namespace cosmosim::physics { class ParticleIdPrecommit; }

namespace cosmosim::workflows::internal {

// Internal test/owner construction surface. The production SourceRuntime owns
// the same implementation and remains the sole source-ID coordination authority.
[[nodiscard]] std::unique_ptr<physics::ParticleIdPrecommit> makeDistributedParticleIdRegistry(
    const parallel::MpiContext& mpi_context, core::MemoryGovernor* governor);

}  // namespace cosmosim::workflows::internal
