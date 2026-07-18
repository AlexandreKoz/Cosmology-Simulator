#pragma once

#include <cstdint>
#include <memory>

#include "cosmosim/core/config.hpp"
#include "cosmosim/core/time_integration.hpp"
#include "cosmosim/core/units.hpp"
#include "cosmosim/workflows/runtime_resources.hpp"

namespace cosmosim::workflows {

// Owns construction and deterministic registration order for the existing
// rung-zero source callbacks. No source model or update law is changed here.
class SourceRuntime {
 public:
  virtual ~SourceRuntime() = default;

  virtual void execute(SourceMutationStageView& view) = 0;

};

[[nodiscard]] std::unique_ptr<SourceRuntime> makeSourceRuntime(
    const core::SimulationConfig& config,
    const core::UnitSystem& units,
    std::uint32_t world_rank);

}  // namespace cosmosim::workflows
