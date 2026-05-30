#include <cassert>
#include <filesystem>
#include <utility>

#include "cosmosim/core/build_config.hpp"
#include "restart_equivalence_harness.hpp"
#include "restart_equivalence_scenarios.hpp"

int main() {
#if COSMOSIM_ENABLE_HDF5
  const auto restart_path = cosmosim::tests::stage8RestartPath("restart_equivalence_multirate_bins");
  auto state = cosmosim::tests::makeStage8DmState(18, "restart_equivalence_multirate_bins");
  auto scheduler = cosmosim::tests::makeStage8Scheduler(static_cast<std::uint32_t>(state.particles.size()), 4);
  scheduler.submitCandidateBin(0, 4, cosmosim::core::TimeStepCandidateSource::kUserClamp);
  scheduler.submitCandidateBin(7, 3, cosmosim::core::TimeStepCandidateSource::kUserClamp);
  auto integrator_state = cosmosim::tests::makeStage8IntegratorState(2, 4);
  auto output_state = cosmosim::tests::makeStage8OutputCadenceState(false);
  auto scenario = cosmosim::tests::makeStage8Scenario(
      std::move(state), integrator_state, std::move(scheduler), std::move(output_state), restart_path, 100, 40);
  const auto result = cosmosim::tests::runRestartEquivalenceScenario(std::move(scenario));
  assert(result.direct_scheduler_state.bin_index == result.restarted_scheduler_state.bin_index);
  assert(result.direct_scheduler_state.next_activation_tick == result.restarted_scheduler_state.next_activation_tick);
  assert(result.direct_scheduler_state.pending_bin_index == result.restarted_scheduler_state.pending_bin_index);
  std::filesystem::remove(restart_path);
#endif
  return 0;
}
