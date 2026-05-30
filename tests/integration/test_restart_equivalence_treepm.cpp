#include <cassert>
#include <cmath>
#include <filesystem>
#include <utility>
#include <vector>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/gravity/tree_pm_coupling.hpp"
#include "restart_equivalence_harness.hpp"
#include "restart_equivalence_scenarios.hpp"

namespace {

[[nodiscard]] constexpr cosmosim::gravity::PmGridShape pmShapeForAvailableBackend(
    std::size_t fftw_n,
    std::size_t fallback_n) noexcept {
#if COSMOSIM_ENABLE_FFTW
  return {fftw_n, fftw_n, fftw_n};
#else
  return {fallback_n, fallback_n, fallback_n};
#endif
}

[[nodiscard]] double wrapUnitBox(double x) {
  x = std::fmod(x, 1.0);
  return x < 0.0 ? x + 1.0 : x;
}

void applyProductionTreePmKick(cosmosim::tests::RestartEquivalenceStepContext& context) {
  std::vector<double> acc_x(context.active_particle_indices.size(), 0.0);
  std::vector<double> acc_y(context.active_particle_indices.size(), 0.0);
  std::vector<double> acc_z(context.active_particle_indices.size(), 0.0);

  cosmosim::gravity::TreePmForceAccumulatorView accumulator{
      context.active_particle_indices,
      acc_x,
      acc_y,
      acc_z,
  };

  cosmosim::gravity::TreePmOptions options;
  options.pm_options.box_size_mpc_comoving = 1.0;
  options.pm_options.scale_factor = context.integrator_state.current_scale_factor;
  options.pm_options.gravitational_constant_code = 1.0;
  options.pm_options.enable_window_deconvolution = true;
  options.tree_options.opening_theta = 0.55;
  options.tree_options.max_leaf_size = 4;
  options.tree_options.gravitational_constant_code = 1.0;
  options.tree_options.softening.epsilon_comoving = 0.015;
  const auto pm_shape = pmShapeForAvailableBackend(16, 8);
  const double mesh_spacing = options.pm_options.box_size_mpc_comoving / static_cast<double>(pm_shape.nx);
  options.split_policy = cosmosim::gravity::makeTreePmSplitPolicyFromMeshSpacing(1.25, 4.5, mesh_spacing);

  cosmosim::gravity::TreePmCoordinator coordinator(pm_shape);
  cosmosim::gravity::TreePmDiagnostics diagnostics;
  coordinator.solveActiveSet(
      context.state.particles.position_x_comoving,
      context.state.particles.position_y_comoving,
      context.state.particles.position_z_comoving,
      context.state.particles.mass_code,
      accumulator,
      options,
      nullptr,
      &diagnostics);
  assert(std::isfinite(diagnostics.force_l2_total));

  for (std::size_t active_slot = 0; active_slot < context.active_particle_indices.size(); ++active_slot) {
    const std::uint32_t particle_index = context.active_particle_indices[active_slot];
    context.state.particles.velocity_x_peculiar[particle_index] += context.dt_code * acc_x[active_slot];
    context.state.particles.velocity_y_peculiar[particle_index] += context.dt_code * acc_y[active_slot];
    context.state.particles.velocity_z_peculiar[particle_index] += context.dt_code * acc_z[active_slot];
    context.state.particles.position_x_comoving[particle_index] = wrapUnitBox(
        context.state.particles.position_x_comoving[particle_index] +
        context.dt_code * context.state.particles.velocity_x_peculiar[particle_index]);
    context.state.particles.position_y_comoving[particle_index] = wrapUnitBox(
        context.state.particles.position_y_comoving[particle_index] +
        context.dt_code * context.state.particles.velocity_y_peculiar[particle_index]);
    context.state.particles.position_z_comoving[particle_index] = wrapUnitBox(
        context.state.particles.position_z_comoving[particle_index] +
        context.dt_code * context.state.particles.velocity_z_peculiar[particle_index]);
  }
}

}  // namespace

int main() {
#if COSMOSIM_ENABLE_HDF5
  const auto restart_path = cosmosim::tests::stage8RestartPath("restart_equivalence_treepm");
  auto state = cosmosim::tests::makeStage8DmState(16, "restart_equivalence_treepm");
  for (std::size_t pidx = 0; pidx < state.particles.size(); ++pidx) {
    state.particles.position_x_comoving[pidx] = std::fmod(0.071 * static_cast<double>(pidx + 1U), 1.0);
    state.particles.position_y_comoving[pidx] = std::fmod(0.113 * static_cast<double>(pidx + 2U), 1.0);
    state.particles.position_z_comoving[pidx] = std::fmod(0.173 * static_cast<double>(pidx + 3U), 1.0);
    state.particles.mass_code[pidx] = 0.75 + 0.03 * static_cast<double>(pidx % 5U);
  }
  auto scheduler = cosmosim::tests::makeStage8Scheduler(static_cast<std::uint32_t>(state.particles.size()), 3);
  auto integrator_state = cosmosim::tests::makeStage8IntegratorState(3, 3);
  auto output_state = cosmosim::tests::makeStage8OutputCadenceState(false);
  auto scenario = cosmosim::tests::makeStage8Scenario(
      std::move(state), integrator_state, std::move(scheduler), std::move(output_state), restart_path, 40, 16);
  scenario.step_kernel = applyProductionTreePmKick;
  scenario.tolerances.position_abs = 1.0e-12;
  scenario.tolerances.velocity_abs = 1.0e-12;
  scenario.tolerances.scalar_abs = 1.0e-12;
  const auto result = cosmosim::tests::runRestartEquivalenceScenario(std::move(scenario));
  const auto pm_state = result.restarted_integrator_state.pm_sync_state.exportPersistentState();
  assert(pm_state.cadence_steps == 3);
  assert(pm_state.field_version == result.direct_integrator_state.pm_sync_state.exportPersistentState().field_version);
  assert(result.restarted_integrator_state.pm_long_range_field_valid);
  std::filesystem::remove(restart_path);
#endif
  return 0;
}
