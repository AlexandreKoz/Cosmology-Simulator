#include <cassert>
#include <filesystem>
#include <utility>

#include "cosmosim/core/build_config.hpp"
#include "restart_equivalence_harness.hpp"
#include "restart_equivalence_scenarios.hpp"

int main() {
#if COSMOSIM_ENABLE_HDF5
  const auto restart_path = cosmosim::tests::stage8RestartPath("restart_equivalence_dm_only");
  auto state = cosmosim::tests::makeStage8DmState(12, "restart_equivalence_dm_only");
  auto scheduler = cosmosim::tests::makeStage8Scheduler(static_cast<std::uint32_t>(state.particles.size()), 2);
  auto integrator_state = cosmosim::tests::makeStage8IntegratorState(1, 2);
  auto output_state = cosmosim::tests::makeStage8OutputCadenceState(false);
  auto scenario = cosmosim::tests::makeStage8Scenario(
      std::move(state), integrator_state, std::move(scheduler), std::move(output_state), restart_path, 100, 40);
  const auto result = cosmosim::tests::runRestartEquivalenceScenario(std::move(scenario));
  assert(result.direct_integrator_state.step_index == 100);
  assert(result.restarted_integrator_state.step_index == 100);
  assert(result.direct_state.particles.size() == result.restarted_state.particles.size());
  std::filesystem::remove(restart_path);
#endif
  return 0;
}
