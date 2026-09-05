#pragma once

#include <cstdint>
#include <memory>

#include "cosmosim/core/config.hpp"
#include "cosmosim/core/simulation_mode.hpp"
#include "cosmosim/core/time_integration.hpp"
#include "cosmosim/core/units.hpp"
#include "cosmosim/workflows/runtime_resources.hpp"
#include "cosmosim/parallel/distributed_memory.hpp"

namespace cosmosim::physics {
class EffectiveMultiphaseEosTable;
}  // namespace cosmosim::physics

namespace cosmosim::core {
struct MemoryReport;
}

namespace cosmosim::workflows {

struct RuntimeServices;

// Owns construction and deterministic registration order for the existing
// rung-zero source callbacks. No source model or update law is changed here.
class SourceRuntime {
 public:
  virtual ~SourceRuntime() = default;

  virtual void execute(SourceMutationStageView& view) = 0;
  [[nodiscard]] virtual core::MemoryReport memoryReport() const = 0;

};

[[nodiscard]] std::unique_ptr<SourceRuntime> makeSourceRuntime(
    const core::SimulationConfig& config,
    const core::ModePolicy& mode_policy,
    const core::UnitSystem& units,
    std::uint32_t world_rank,
    const parallel::MpiContext& mpi_context,
    std::shared_ptr<const physics::EffectiveMultiphaseEosTable> effective_eos_table = nullptr,
    const RuntimeServices* runtime_services = nullptr);

}  // namespace cosmosim::workflows
