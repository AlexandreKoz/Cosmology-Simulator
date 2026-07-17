#pragma once

#include <cstdint>
#include <span>

#include "cosmosim/core/time_integration.hpp"
#include "cosmosim/workflows/gravity_runtime.hpp"
#include "cosmosim/workflows/reference_workflow.hpp"

namespace cosmosim::core {
class ProfilerSession;
struct FrozenConfig;
struct SimulationConfig;
}

namespace cosmosim::io {
struct OutputCadencePersistentState;
}

namespace cosmosim::workflows::internal {

struct PendingOutputBoundary {
  bool snapshot_due = false;
  bool checkpoint_due = false;
  bool step_event_due = false;
  bool time_event_due = false;
  std::uint64_t snapshot_interval_steps = 0;
  std::uint64_t next_snapshot_step_index = 0;
  double snapshot_interval_time_code = 0.0;
  double next_snapshot_time_code = 0.0;
  double restart_resume_dt_time_code = 0.0;
};

[[nodiscard]] PendingOutputBoundary initializeOutputCadence(
    const core::SimulationConfig& config,
    const ReferenceWorkflowOptions& options,
    const core::IntegratorState& integrator_state,
    const io::OutputCadencePersistentState* restored_output);

[[nodiscard]] bool timelineTimesEqual(double lhs, double rhs);
[[nodiscard]] double nextOrderedTimeEvent(
    double timeline_begin_time_code,
    double current_time_code,
    double interval_time_code);
void latchOutputRequestForCompletedStep(
    const core::SimulationConfig& config,
    const ReferenceWorkflowOptions& options,
    std::uint64_t next_step_index,
    double next_time_code,
    PendingOutputBoundary& pending);
[[nodiscard]] core::StepBoundaryKind requestedBoundaryForPendingOutput(
    const PendingOutputBoundary& pending);

// Persistence owner for cadence, schema-preserving snapshot/checkpoint writes,
// field-aware readback verification, and restart state assembly.
class OutputRestartRuntime final {
 public:
  OutputRestartRuntime(
      const core::FrozenConfig& frozen_config,
      const core::SimulationConfig& config,
      const core::HierarchicalTimeBinScheduler& scheduler,
      const core::HierarchicalTimeBinScheduler& gas_cell_scheduler,
      const workflows::GravityRestartStateProvider& gravity_state,
      ReferenceWorkflowReport& report,
      core::ProfilerSession& profiler,
      PendingOutputBoundary& pending_output,
      bool write_outputs_enabled);

  void execute(OutputRestartStageView& view);

 private:
  const core::FrozenConfig& m_frozen_config;
  const core::SimulationConfig& m_config;
  const core::HierarchicalTimeBinScheduler& m_scheduler;
  const core::HierarchicalTimeBinScheduler& m_gas_cell_scheduler;
  const workflows::GravityRestartStateProvider& m_gravity_state;
  ReferenceWorkflowReport& m_report;
  core::ProfilerSession& m_profiler;
  PendingOutputBoundary& m_pending_output;
  bool m_write_outputs_enabled = false;
};

}  // namespace cosmosim::workflows::internal
