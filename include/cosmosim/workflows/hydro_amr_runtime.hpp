#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>

#include "cosmosim/core/time_integration.hpp"
#include "cosmosim/hydro/hydro_core_solver.hpp"
#include "cosmosim/workflows/gravity_runtime.hpp"

namespace cosmosim::workflows {

struct RuntimeServices;

// Owner interface for the currently supported rung-zero hydro/AMR path.
// Communication counters are exposed as values so callers cannot reach the
// owner's ghost, patch, conserved-state, or geometry workspaces.
class HydroAmrRuntime {
 public:
  virtual ~HydroAmrRuntime() = default;

  virtual void execute(HydroAmrStageView& view) = 0;

  [[nodiscard]] virtual const hydro::HydroProfileEvent&
  lastHydroProfile() const noexcept = 0;
  [[nodiscard]] virtual const core::HydroCflDiagnostics&
  lastHydroCflDiagnostics() const noexcept = 0;
  [[nodiscard]] virtual std::uint64_t
  ghostExchangeBytesRecent() const noexcept = 0;
  [[nodiscard]] virtual std::size_t
  remoteImportedGhostCount() const noexcept = 0;
  [[nodiscard]] virtual std::size_t
  remoteInterfaceFaceCount() const noexcept = 0;
  [[nodiscard]] virtual std::size_t
  remoteStaleInvalidPayloadCount() const noexcept = 0;

 protected:
  [[nodiscard]] static core::StepContext& stageContext(
      HydroAmrStageView& view) {
    return view.ownerContext();
  }
};

[[nodiscard]] std::unique_ptr<HydroAmrRuntime> makeHydroAmrRuntime(
    const core::SimulationConfig& config,
    const core::ModePolicy& mode_policy,
    const GravityAccelerationProvider& gravity_acceleration,
    const RuntimeServices& services);

}  // namespace cosmosim::workflows
