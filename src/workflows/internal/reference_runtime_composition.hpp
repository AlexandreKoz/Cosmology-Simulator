#pragma once

#include <cstdint>
#include <filesystem>
#include <memory>

#include "cosmosim/workflows/runtime_module_registry.hpp"

namespace cosmosim::core {
class HierarchicalTimeBinScheduler;
class ProfilerSession;
struct FrozenConfig;
struct ModePolicy;
struct SimulationConfig;
struct UnitSystem;
}

namespace cosmosim::workflows {
class GravityRuntime;
class HydroAmrRuntime;
struct ReferenceWorkflowOptions;
struct ReferenceWorkflowReport;
struct RuntimeServices;

namespace internal {
struct PendingOutputBoundary;

struct ReferenceRuntimeCompositionInputs {
  const core::FrozenConfig& frozen_config;
  const core::SimulationConfig& config;
  const core::ModePolicy& mode_policy;
  const core::UnitSystem& units;
  const RuntimeServices& services;
  const core::HierarchicalTimeBinScheduler& particle_scheduler;
  const core::HierarchicalTimeBinScheduler& gas_cell_scheduler;
  ReferenceWorkflowReport& report;
  core::ProfilerSession& profiler;
  PendingOutputBoundary& pending_output;
  const ReferenceWorkflowOptions& options;
  std::filesystem::path zoom_region_path;
  std::uint32_t world_rank = 0;
};

// The execution plan owns every stage-contributing service. The exposed
// gravity/hydro handles are narrow lifecycle/diagnostic surfaces needed by the
// time coordinator outside stage execution.
struct ReferenceRuntimeComposition {
  std::shared_ptr<GravityRuntime> gravity;
  std::shared_ptr<HydroAmrRuntime> hydro_amr;
  RuntimeExecutionPlan execution_plan;
};

[[nodiscard]] ReferenceRuntimeComposition buildReferenceRuntimeComposition(
    ReferenceRuntimeCompositionInputs inputs);

}  // namespace internal
}  // namespace cosmosim::workflows
