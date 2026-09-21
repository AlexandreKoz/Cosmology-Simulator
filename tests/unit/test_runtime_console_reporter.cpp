#include <cstdlib>
#include <iostream>
#include <optional>
#include <sstream>
#include <string>

#include "cosmosim/workflows/runtime_console_reporter.hpp"

namespace {

void require(bool condition, const char* message) {
  if (!condition) {
    std::cerr << "test_runtime_console_reporter failed: " << message << '\n';
    std::exit(1);
  }
}

std::size_t countOccurrences(const std::string& text, const std::string& needle) {
  std::size_t count = 0U;
  std::size_t pos = 0U;
  while ((pos = text.find(needle, pos)) != std::string::npos) {
    ++count;
    pos += needle.size();
  }
  return count;
}

}  // namespace

int main() {
  using cosmosim::workflows::RuntimeConsoleOptions;
  using cosmosim::workflows::RuntimeConsoleReporter;
  using cosmosim::workflows::RuntimeConsoleStartupStatus;
  using cosmosim::workflows::RuntimeConsoleStepStatus;

  require(
      RuntimeConsoleReporter::shouldEmitStepForCadence(1U, std::nullopt, 10U, 0.0, 30.0),
      "first productive step should emit when cadence is enabled");
  require(
      !RuntimeConsoleReporter::shouldEmitStepForCadence(1U, std::nullopt, 0U, 0.0, 0.0),
      "both zero cadence triggers should disable productive-step output");
  require(
      !RuntimeConsoleReporter::shouldEmitStepForCadence(5U, 1U, 10U, 1.0, 30.0),
      "step cadence should not fire early");
  require(
      RuntimeConsoleReporter::shouldEmitStepForCadence(11U, 1U, 10U, 1.0, 30.0),
      "step cadence should fire at the configured interval");
  require(
      RuntimeConsoleReporter::shouldEmitStepForCadence(2U, 1U, 10U, 31.0, 30.0),
      "wall-time cadence should fire after maximum silence");
  require(
      !RuntimeConsoleReporter::shouldEmitStepForCadence(1U, 1U, 1U, 31.0, 30.0),
      "the same productive boundary must not emit twice");

  std::ostringstream info;
  std::ostringstream diag;
  RuntimeConsoleOptions options;
  options.enabled = true;
  options.status_every_steps = 1U;
  options.status_seconds = 30.0;
  options.config_path = "configs/test config.param.txt";
  RuntimeConsoleReporter reporter(options, 0, info, diag);

  RuntimeConsoleStartupStatus noncosmo;
  noncosmo.run_name = "isolated";
  noncosmo.run_directory = "outputs/isolated";
  noncosmo.rank_directory = "outputs/isolated_rank000";
  noncosmo.simulation_mode = "isolated_galaxy";
  noncosmo.global_particle_count = 4U;
  noncosmo.t_code = 0.0;
  noncosmo.endpoint_t_code = 1.0;
  reporter.emitStartup(noncosmo);
  require(info.str().find("[CHUI][START]") != std::string::npos,
          "startup record should be emitted on rank zero");
  require(info.str().find(" z=") == std::string::npos,
          "non-cosmological startup must not fabricate redshift");
  require(info.str().find("run_directory=\"outputs/isolated\"") != std::string::npos,
          "startup should report the logical run directory");
  require(info.str().find("rank_directory=\"outputs/isolated_rank000\"") != std::string::npos,
          "startup should distinguish the rank-local artifact directory");

  RuntimeConsoleStepStatus step;
  step.step_index = 1U;
  step.t_code = 0.1;
  step.dt_time_code = 0.1;
  step.a_scale = 0.5;
  step.redshift = 1.0;
  step.active_particle_count = 4U;
  step.total_particle_count = 4U;
  step.wall_step_seconds = 0.01;
  step.pm_activity = "refresh";
  reporter.emitStep(step);
  require(info.str().find("[CHUI][STEP]") != std::string::npos,
          "productive-step record should be emitted");
  require(info.str().find(" a=5.00000000e-01 z=1.00000000e+00") != std::string::npos,
          "cosmological step should print finite a/z");
  require(countOccurrences(info.str(), "[CHUI][STEP]") == 1U,
          "one boundary must emit at most one step record");

  reporter.emitSnapshotCommitted(
      1U, "outputs/isolated/snapshots/snap_001.0.hdf5",
      "outputs/isolated/snapshots/snap_001.complete", 2U);
  require(info.str().find("[CHUI][SNAPSHOT] step=1 state=committed members=2") !=
              std::string::npos,
          "snapshot record should expose logical-set member count");
  require(info.str().find("manifest=\"outputs/isolated/snapshots/snap_001.complete\"") !=
              std::string::npos,
          "snapshot record should expose the transactional completion manifest");
  require(info.str().find("rank0_member=\"outputs/isolated/snapshots/snap_001.0.hdf5\"") !=
              std::string::npos,
          "multi-rank snapshot record should label rank-zero member explicitly");

  reporter.emitDone(
      1U, 0.1, std::nullopt, std::nullopt,
      "outputs/isolated", "outputs/isolated_rank000",
      "outputs/isolated_rank000/normalized_config.param.txt",
      "outputs/isolated_rank000/operational_events.json");
  require(info.str().find("[CHUI][DONE]") != std::string::npos,
          "completion record should be emitted");
  require(info.str().find("run_directory=\"outputs/isolated\"") != std::string::npos,
          "completion should retain logical run-directory truth");

  std::ostringstream quiet_info;
  std::ostringstream quiet_diag;
  RuntimeConsoleOptions quiet_options;
  quiet_options.enabled = true;
  quiet_options.quiet = true;
  RuntimeConsoleReporter quiet_reporter(quiet_options, 0, quiet_info, quiet_diag);
  quiet_reporter.emitRuntimePhase("test");
  quiet_reporter.emitWarning("test", "warning survives quiet mode");
  require(quiet_info.str().empty(), "quiet mode should suppress routine information");
  require(quiet_diag.str().find("[CHUI][WARN]") != std::string::npos,
          "quiet mode must preserve warnings");

  std::ostringstream rank1_info;
  std::ostringstream rank1_diag;
  RuntimeConsoleReporter rank1_reporter(options, 1, rank1_info, rank1_diag);
  rank1_reporter.emitRuntimePhase("test");
  rank1_reporter.emitWarning("test", "root-only reporter warning");
  require(rank1_info.str().empty() && rank1_diag.str().empty(),
          "normal console reporter must be rank-zero owned");

  return 0;
}
