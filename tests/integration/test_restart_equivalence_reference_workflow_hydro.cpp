#include <array>
#include <cassert>
#include <filesystem>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/io/restart_checkpoint.hpp"
#include "reference_workflow_hydro_test_fixture.hpp"

int main() {
#if !COSMOSIM_ENABLE_HDF5
  // Established guarded convention: this target remains CTest-registerable in
  // CPU builds, but full workflow restart equivalence is executed only with HDF5.
  return 0;
#else
  namespace fixture = cosmosim::test::workflow_hydro_fixture;
  constexpr std::array<std::uint32_t, fixture::k_cell_count> canonical{{0, 1, 2, 3, 4, 5, 6, 7}};
  constexpr std::array<std::uint32_t, fixture::k_cell_count> shuffled{{5, 0, 7, 2, 4, 1, 6, 3}};
  const auto initial_state = fixture::makeState(canonical);
  const auto frozen = fixture::makeFrozenConfig("h1_reference_workflow_hydro_restart");
  cosmosim::workflows::ReferenceWorkflowRunner runner(frozen);
  const std::filesystem::path root =
      std::filesystem::temp_directory_path() / "cosmosim_h1_reference_workflow_restart";

  const auto direct_report = runner.run(
      root / "direct",
      cosmosim::workflows::ReferenceWorkflowOptions{
          .write_outputs = true,
          .initial_state_override = &initial_state,
          .max_steps_override = 2});
  assert(direct_report.completed_steps == 2U);
  assert(direct_report.restart_roundtrip_ok);
  const auto direct_restart = cosmosim::io::readRestartCheckpointHdf5(direct_report.restart_path);

  const auto split_first_report = runner.run(
      root / "split_first",
      cosmosim::workflows::ReferenceWorkflowOptions{
          .write_outputs = true,
          .initial_state_override = &initial_state,
          .max_steps_override = 1});
  assert(split_first_report.completed_steps == 1U);
  assert(split_first_report.restart_roundtrip_ok);
  auto checkpoint = cosmosim::io::readRestartCheckpointHdf5(split_first_report.restart_path);
  const std::uint64_t identity_generation_before_reorder = checkpoint.state.gasCellIdentityGeneration();
  fixture::reorderGasRowsAndScheduler(checkpoint, shuffled);
  assert(checkpoint.state.gasCellIdentityGeneration() > identity_generation_before_reorder);
  assert(checkpoint.state.rowForGasCellId(70001U).value() != 0U);

  const auto resumed_report = runner.run(
      root / "resumed",
      cosmosim::workflows::ReferenceWorkflowOptions{
          .write_outputs = true,
          .restart_state_override = &checkpoint,
          .max_steps_override = 1});
  assert(resumed_report.completed_steps == 1U);
  assert(resumed_report.restart_roundtrip_ok);
  const auto resumed_restart = cosmosim::io::readRestartCheckpointHdf5(resumed_report.restart_path);

  fixture::assertEquivalentByStableId(direct_restart.state, resumed_restart.state);
  fixture::assertEquivalentGasSchedulerById(direct_restart, resumed_restart);
  fixture::assertConservationEquivalent(
      fixture::globalConservation(direct_restart.state),
      fixture::globalConservation(resumed_restart.state));
  fixture::assertEquivalentCfl(
      direct_report.final_hydro_cfl_diagnostics,
      resumed_report.final_hydro_cfl_diagnostics);

  assert(direct_restart.integrator_state.step_index == resumed_restart.integrator_state.step_index);
  fixture::assertNear(
      direct_restart.integrator_state.current_time_code,
      resumed_restart.integrator_state.current_time_code);
  fixture::assertNear(
      direct_restart.integrator_state.current_scale_factor,
      resumed_restart.integrator_state.current_scale_factor);
  assert(direct_restart.scheduler_state.current_tick == resumed_restart.scheduler_state.current_tick);
  assert(direct_restart.scheduler_state.bin_index == resumed_restart.scheduler_state.bin_index);
  assert(direct_restart.gas_cell_scheduler_state.current_tick ==
         resumed_restart.gas_cell_scheduler_state.current_tick);
  // The resumed branch deliberately rebuilt the dense-row map at the safe
  // checkpoint boundary.  Generation must record that topology/identity-map
  // event even though every physical quantity is equivalent by gas_cell_id.
  assert(resumed_restart.state.gasCellIdentityGeneration() ==
         direct_restart.state.gasCellIdentityGeneration() + 1U);
  direct_restart.state.requireGasCellIdentityMapCoversDenseRows(
      "direct workflow restart-equivalence final state");
  resumed_restart.state.requireGasCellIdentityMapCoversDenseRows(
      "reordered workflow restart-equivalence final state");
  assert(direct_restart.state.patches.patch_id == resumed_restart.state.patches.patch_id);
  assert(direct_restart.state.patches.cell_dim_x == resumed_restart.state.patches.cell_dim_x);
  assert(direct_restart.output_cadence_state.last_completed_step_index ==
         resumed_restart.output_cadence_state.last_completed_step_index);
  assert(direct_restart.output_cadence_state.next_snapshot_step_index ==
         resumed_restart.output_cadence_state.next_snapshot_step_index);

  fixture::cleanup(root);
  return 0;
#endif
}
