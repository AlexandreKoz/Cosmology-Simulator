// 20261025 optimization campaign: authored, NOT RUN — deferred to maintainer.
#include "../support/test_require.hpp"

#include <cstdint>
#include <stdexcept>
#include <string>

#include "cosmosim/core/config.hpp"
#include "cosmosim/core/time_integration.hpp"
#include "cosmosim/workflows/output_restart_runtime.hpp"
#include "cosmosim/workflows/reference_workflow.hpp"

namespace {

using cosmosim::workflows::internal::PendingOutputBoundary;
using cosmosim::workflows::internal::initializeOutputCadence;
using cosmosim::workflows::internal::latchOutputRequestForCompletedStep;
using cosmosim::workflows::internal::nextRestartStepIndexAfter;

cosmosim::core::FrozenConfig loadOutput(const std::string& output_body) {
  return cosmosim::core::loadFrozenConfigFromString(
      "[mode]\nmode = zoom_in\n[output]\n" + output_body, "output_cadence_test");
}

// Legacy policy: every snapshot event also latches a checkpoint.
void testLegacyCoupledCheckpointCadence() {
  const auto frozen = loadOutput("snapshot_interval_steps = 4\nwrite_restarts = true\n");
  cosmosim::workflows::ReferenceWorkflowOptions options;
  cosmosim::core::IntegratorState integrator;
  integrator.step_index = 0U;
  PendingOutputBoundary pending = initializeOutputCadence(frozen.config, options, integrator, nullptr);
  CHUI_REQUIRE(pending.restart_interval_steps == 0U);
  for (std::uint64_t next = 1U; next < 4U; ++next) {
    latchOutputRequestForCompletedStep(frozen.config, options, next, 0.0, pending);
    CHUI_REQUIRE(!pending.snapshot_due && !pending.checkpoint_due);
  }
  latchOutputRequestForCompletedStep(frozen.config, options, 4U, 0.0, pending);
  CHUI_REQUIRE(pending.snapshot_due && pending.checkpoint_due);
}

// Independent policy: snapshots and checkpoints latch on their own cadences.
void testIndependentCheckpointCadence() {
  const auto frozen = loadOutput(
      "snapshot_interval_steps = 4\nwrite_restarts = true\nrestart_interval_steps = 6\n");
  cosmosim::workflows::ReferenceWorkflowOptions options;
  cosmosim::core::IntegratorState integrator;
  integrator.step_index = 0U;
  PendingOutputBoundary pending = initializeOutputCadence(frozen.config, options, integrator, nullptr);
  CHUI_REQUIRE(pending.restart_interval_steps == 6U);
  CHUI_REQUIRE(pending.next_restart_step_index == 6U);

  latchOutputRequestForCompletedStep(frozen.config, options, 4U, 0.0, pending);
  CHUI_REQUIRE(pending.snapshot_due);
  CHUI_REQUIRE(!pending.checkpoint_due);  // a science snapshot no longer implies a checkpoint
  pending = initializeOutputCadence(frozen.config, options, integrator, nullptr);

  latchOutputRequestForCompletedStep(frozen.config, options, 6U, 0.0, pending);
  CHUI_REQUIRE(!pending.snapshot_due);  // a checkpoint does not require a science snapshot
  CHUI_REQUIRE(pending.checkpoint_due && pending.restart_step_event_due);
}

// Resumed runs derive the next checkpoint from the persisted step index.
void testResumedCheckpointCadenceIsDerived() {
  CHUI_REQUIRE(nextRestartStepIndexAfter(0U, 6U) == 6U);
  CHUI_REQUIRE(nextRestartStepIndexAfter(5U, 6U) == 6U);
  CHUI_REQUIRE(nextRestartStepIndexAfter(6U, 6U) == 12U);
  CHUI_REQUIRE(nextRestartStepIndexAfter(13U, 6U) == 18U);
  const auto frozen = loadOutput("restart_interval_steps = 6\n");
  cosmosim::workflows::ReferenceWorkflowOptions options;
  cosmosim::core::IntegratorState integrator;
  integrator.step_index = 12U;
  const PendingOutputBoundary pending =
      initializeOutputCadence(frozen.config, options, integrator, nullptr);
  CHUI_REQUIRE(pending.next_restart_step_index == 18U);
  bool rejected = false;
  try {
    (void)nextRestartStepIndexAfter(1U, 0U);
  } catch (const std::invalid_argument&) {
    rejected = true;
  }
  CHUI_REQUIRE(rejected);
}

}  // namespace

int main() {
  testLegacyCoupledCheckpointCadence();
  testIndependentCheckpointCadence();
  testResumedCheckpointCadenceIsDerived();
  return 0;
}
