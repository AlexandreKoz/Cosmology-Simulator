#include <array>
#include <cassert>
#include <filesystem>
#include <stdexcept>
#include <string>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/io/restart_checkpoint.hpp"
#include "reference_workflow_hydro_test_fixture.hpp"

int main() {
#if !COSMOSIM_ENABLE_HDF5
  // HDF5 is required here solely to read final workflow state without adding a
  // production-only state-export seam. The workflow target still compiles in
  // CPU/no-HDF5 configurations; HDF5 acceptance is exercised in hdf5-debug.
  return 0;
#else
  namespace fixture = cosmosim::test::workflow_hydro_fixture;
  constexpr std::array<std::uint32_t, fixture::k_cell_count> canonical{{0, 1, 2, 3, 4, 5, 6, 7}};
  constexpr std::array<std::uint32_t, fixture::k_cell_count> shuffled{{5, 0, 7, 2, 4, 1, 6, 3}};

  const auto canonical_state = fixture::makeState(canonical);
  const auto shuffled_state = fixture::makeState(shuffled);
  const auto canonical_frozen = fixture::makeFrozenConfig("h1_workflow_row_order_canonical");
  const auto shuffled_frozen = fixture::makeFrozenConfig("h1_workflow_row_order_shuffled");
  const std::filesystem::path output_root =
      std::filesystem::temp_directory_path() / "cosmosim_h1_workflow_row_order";

  cosmosim::workflows::ReferenceWorkflowRunner canonical_runner(canonical_frozen);
  cosmosim::workflows::ReferenceWorkflowRunner shuffled_runner(shuffled_frozen);
  const auto canonical_report = canonical_runner.run(
      output_root,
      cosmosim::workflows::ReferenceWorkflowOptions{
          .write_outputs = true,
          .initial_state_override = &canonical_state,
          .max_steps_override = 2});
  const auto shuffled_report = shuffled_runner.run(
      output_root,
      cosmosim::workflows::ReferenceWorkflowOptions{
          .write_outputs = true,
          .initial_state_override = &shuffled_state,
          .max_steps_override = 2});

  assert(canonical_report.completed_steps == 2U);
  assert(shuffled_report.completed_steps == 2U);
  assert(canonical_report.restart_roundtrip_ok);
  assert(shuffled_report.restart_roundtrip_ok);
  fixture::assertEquivalentCfl(
      canonical_report.final_hydro_cfl_diagnostics,
      shuffled_report.final_hydro_cfl_diagnostics);

  const auto canonical_restart = cosmosim::io::readRestartCheckpointHdf5(canonical_report.restart_path);
  const auto shuffled_restart = cosmosim::io::readRestartCheckpointHdf5(shuffled_report.restart_path);
  fixture::assertEquivalentByStableId(canonical_restart.state, shuffled_restart.state);
  fixture::assertConservationEquivalent(
      fixture::globalConservation(canonical_restart.state),
      fixture::globalConservation(shuffled_restart.state));
  fixture::assertEquivalentGasSchedulerById(canonical_restart, shuffled_restart);

  // The dense storage rows differ by construction. Any equivalence therefore
  // proves gather/solve/scatter and periodic ghost ownership use physical layout
  // plus stable ID rather than canonical dense order.
  assert(canonical_restart.state.rowForGasCellId(70001U).value() !=
         shuffled_restart.state.rowForGasCellId(70001U).value());

  // Negative geometry coverage: an incomplete Cartesian product must fail
  // loudly in the actual workflow rather than entering a near-cubic fallback.
  auto invalid_state = fixture::makeState(canonical);
  invalid_state.cells.center_x_comoving[7] = invalid_state.cells.center_x_comoving[6];
  cosmosim::workflows::ReferenceWorkflowRunner invalid_runner(
      fixture::makeFrozenConfig("h1_workflow_invalid_cartesian_layout"));
  bool rejected = false;
  try {
    (void)invalid_runner.run(
        output_root,
        cosmosim::workflows::ReferenceWorkflowOptions{
            .write_outputs = false,
            .initial_state_override = &invalid_state,
            .max_steps_override = 1});
  } catch (const std::runtime_error& error) {
    const std::string message = error.what();
    rejected = message.find("Cartesian") != std::string::npos ||
        message.find("cartesian") != std::string::npos;
  }
  assert(rejected);

  fixture::cleanup(output_root);
  return 0;
#endif
}
