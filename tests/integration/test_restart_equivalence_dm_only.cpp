#include <cassert>
#include <filesystem>
#include <utility>

#include "cosmosim/core/build_config.hpp"
#include "restart_equivalence_harness.hpp"
#include "restart_equivalence_scenarios.hpp"
#include "workflows/internal/output_verification.hpp"

int main() {
#if COSMOSIM_ENABLE_HDF5
  const auto restart_path = cosmosim::tests::stage8RestartPath("restart_equivalence_dm_only");
  auto state = cosmosim::tests::makeStage8DmState(12, "restart_equivalence_dm_only");
  auto scheduler = cosmosim::tests::makeStage8Scheduler(static_cast<std::uint32_t>(state.particles.size()), 2);
  auto integrator_state = cosmosim::tests::makeStage8IntegratorState(1, 2);
  auto output_state = cosmosim::tests::makeStage8OutputCadenceState(false);
  auto scenario = cosmosim::tests::makeStage8Scenario(
      std::move(state), integrator_state, std::move(scheduler), std::move(output_state), restart_path, 100, 40);
  auto result = cosmosim::tests::runRestartEquivalenceScenario(std::move(scenario));
  assert(result.direct_integrator_state.step_index == 100);
  assert(result.restarted_integrator_state.step_index == 100);
  assert(result.direct_state.particles.size() == result.restarted_state.particles.size());
  assert(cosmosim::workflows::internal::restartRuntimeStateExactlyEquivalent(
      result.restarted_state, result.direct_state));
  assert(!result.restarted_state.particle_sidecar.last_drift_time_code.empty());
  assert(!result.restarted_state.particle_sidecar.last_drift_scale_factor.empty());

  const double saved_drift_time =
      result.restarted_state.particle_sidecar.last_drift_time_code.front();
  result.restarted_state.particle_sidecar.last_drift_time_code.front() =
      saved_drift_time + 0.125;
  assert(!cosmosim::workflows::internal::restartRuntimeStateExactlyEquivalent(
      result.restarted_state, result.direct_state));
  result.restarted_state.particle_sidecar.last_drift_time_code.front() = saved_drift_time;

  const double saved_drift_scale_factor =
      result.restarted_state.particle_sidecar.last_drift_scale_factor.front();
  result.restarted_state.particle_sidecar.last_drift_scale_factor.front() =
      saved_drift_scale_factor + 0.03125;
  assert(!cosmosim::workflows::internal::restartRuntimeStateExactlyEquivalent(
      result.restarted_state, result.direct_state));
  result.restarted_state.particle_sidecar.last_drift_scale_factor.front() =
      saved_drift_scale_factor;
  assert(cosmosim::workflows::internal::restartRuntimeStateExactlyEquivalent(
      result.restarted_state, result.direct_state));
  std::filesystem::remove(restart_path);
#endif
  return 0;
}
