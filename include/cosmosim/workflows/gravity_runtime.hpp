#pragma once

#include <cstdint>
#include <filesystem>
#include <memory>
#include <span>

#include "cosmosim/core/profiling.hpp"
#include "cosmosim/core/simulation_mode.hpp"
#include "cosmosim/core/time_integration.hpp"
#include "cosmosim/gravity/tree_pm_coupling.hpp"
#include "cosmosim/io/restart_checkpoint.hpp"
#include "cosmosim/parallel/distributed_memory.hpp"
#include "cosmosim/workflows/reference_workflow.hpp"
#include "cosmosim/workflows/runtime_resources.hpp"

namespace cosmosim::workflows {

struct RuntimeServices;

// Gravity-owned startup normalization and diagnostic policy description.
// These are called by the composition lifecycle before the first gravity task.
void initializeGravityState(
    core::SimulationState& state,
    const core::SimulationConfig& config);
[[nodiscard]] std::string describeGravitySofteningPolicy(
    const core::SimulationConfig& config,
    const core::SimulationState& state);

// Narrow read surface consumed by hydro/source orchestration. Dense-row
// acceleration mirrors are invalidated by particle migration through
// commitParticleDecompositionChange().
class GravityAccelerationProvider {
 public:
  virtual ~GravityAccelerationProvider() = default;

  [[nodiscard]] virtual std::span<const double> cellAccelX() const noexcept = 0;
  [[nodiscard]] virtual std::span<const double> cellAccelY() const noexcept = 0;
  [[nodiscard]] virtual std::span<const double> cellAccelZ() const noexcept = 0;
  [[nodiscard]] virtual std::span<const double> particleAccelX() const noexcept = 0;
  [[nodiscard]] virtual std::span<const double> particleAccelY() const noexcept = 0;
  [[nodiscard]] virtual std::span<const double> particleAccelZ() const noexcept = 0;
};

// Persistence-facing gravity surface. Output/restart depends on this contract,
// never on the concrete solver implementation.
class GravityRestartStateProvider {
 public:
  virtual ~GravityRestartStateProvider() = default;

  [[nodiscard]] virtual std::uint64_t decompositionEpoch() const noexcept = 0;
  [[nodiscard]] virtual const gravity::PmGridShape& pmGridShape() const noexcept = 0;
  [[nodiscard]] virtual const parallel::DistributedExecutionTopology&
  runtimeTopology() const noexcept = 0;
  [[nodiscard]] virtual io::GravityForceCachePersistentState
  exportRestartForceCache(const core::SimulationState& state) const = 0;
};

// Owner interface for rung-zero TreePM/PM cadence and force-cache lifecycle.
// Stage dispatch enters only through the typed GravityStageView contract.
class GravityRuntime : public GravityAccelerationProvider,
                       public GravityRestartStateProvider {
 public:
  ~GravityRuntime() override = default;

  virtual void execute(GravityStageView& view) = 0;

  [[nodiscard]] virtual std::size_t pmGridSize() const noexcept = 0;
  [[nodiscard]] virtual std::uint64_t longRangeRefreshCount() const noexcept = 0;
  [[nodiscard]] virtual std::uint64_t longRangeReuseCount() const noexcept = 0;
  [[nodiscard]] virtual int pmCadenceSteps() const noexcept = 0;
  [[nodiscard]] virtual double gravitationalConstantCode() const noexcept = 0;
  [[nodiscard]] virtual std::span<const ReferenceWorkflowReport::TreePmCadenceRecord>
  cadenceRecords() const noexcept = 0;
  [[nodiscard]] virtual core::MemoryReport memoryReport() const = 0;
  [[nodiscard]] virtual const parallel::DecompositionRuntimeMeasurements&
  lastRuntimeDecompositionMeasurements() const noexcept = 0;

  virtual void restoreDecompositionEpoch(std::uint64_t decomposition_epoch) = 0;
  virtual void commitParticleDecompositionChange() = 0;
  virtual void importRestartForceCache(
      const io::GravityForceCachePersistentState& cache,
      const core::SimulationState& state) = 0;

};

[[nodiscard]] std::unique_ptr<GravityRuntime> makeGravityRuntime(
    const core::SimulationConfig& config,
    const core::ModePolicy& mode_policy,
    const RuntimeServices& services,
    std::filesystem::path zoom_region_path = {});

}  // namespace cosmosim::workflows
