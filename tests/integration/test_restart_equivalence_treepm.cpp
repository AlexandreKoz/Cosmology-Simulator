#include <cassert>
#include <filesystem>
#include <utility>

#include "cosmosim/core/build_config.hpp"
#include "restart_equivalence_harness.hpp"
#include "restart_equivalence_scenarios.hpp"

int main() {
#if COSMOSIM_ENABLE_HDF5
  const auto restart_path = cosmosim::tests::stage8RestartPath("restart_equivalence_treepm");
  auto state = cosmosim::tests::makeStage8DmState(16, "restart_equivalence_treepm");
  auto scheduler = cosmosim::tests::makeStage8Scheduler(static_cast<std::uint32_t>(state.particles.size()), 3);
  auto integrator_state = cosmosim::tests::makeStage8IntegratorState(3, 3);
  auto output_state = cosmosim::tests::makeStage8OutputCadenceState(false);
  auto scenario = cosmosim::tests::makeStage8Scenario(
      std::move(state), integrator_state, std::move(scheduler), std::move(output_state), restart_path, 100, 40);
  const auto result = cosmosim::tests::runRestartEquivalenceScenario(std::move(scenario));
  const auto pm_state = result.restarted_integrator_state.pm_sync_state.exportPersistentState();
  assert(pm_state.cadence_steps == 3);
  assert(pm_state.field_version == result.direct_integrator_state.pm_sync_state.exportPersistentState().field_version);
  assert(result.restarted_integrator_state.pm_long_range_field_valid);
  std::filesystem::remove(restart_path);
#endif
  return 0;
}
