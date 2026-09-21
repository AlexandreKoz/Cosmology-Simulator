#include "cosmosim/workflows/runtime_console_reporter.hpp"

#include <cmath>
#include <iomanip>
#include <limits>
#include <ostream>
#include <sstream>
#include <utility>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/core/memory_accounting.hpp"
#include "cosmosim/core/memory_governor.hpp"
#include "cosmosim/core/profiling.hpp"
#include "cosmosim/core/version.hpp"

namespace cosmosim::workflows {
namespace {

[[nodiscard]] std::string formatDouble(double value) {
  std::ostringstream out;
  out << std::scientific << std::setprecision(8) << value;
  return out.str();
}

[[nodiscard]] std::string formatSeconds(double value) {
  std::ostringstream out;
  out << std::fixed << std::setprecision(3) << value;
  return out.str();
}

[[nodiscard]] std::string optionalDouble(
    const std::optional<double>& value,
    std::string_view unavailable = "n/a") {
  if (!value.has_value() || !std::isfinite(*value)) {
    return std::string(unavailable);
  }
  return formatDouble(*value);
}

[[nodiscard]] std::string optionalBytes(
    const std::optional<std::uint64_t>& value) {
  if (!value.has_value()) {
    return "n/a";
  }
  if (*value == std::numeric_limits<std::uint64_t>::max()) {
    return "unlimited";
  }
  return std::to_string(*value);
}

[[nodiscard]] std::optional<std::uint64_t> memoryHeadroom(
    const core::MemoryReport& report) {
  if (!report.governor_snapshot.has_value()) {
    return std::nullopt;
  }
  return report.governor_snapshot->headroom_bytes;
}

[[nodiscard]] std::string memoryPressure(const core::MemoryReport& report) {
  if (!report.governor_snapshot.has_value()) {
    return "unavailable";
  }
  return std::string(core::memoryPressureLabel(report.governor_snapshot->pressure));
}

}  // namespace

RuntimeConsoleReporter::RuntimeConsoleReporter(
    RuntimeConsoleOptions options,
    int world_rank,
    std::ostream& info_stream,
    std::ostream& diagnostic_stream)
    : m_options(std::move(options)),
      m_world_rank(world_rank),
      m_info_stream(&info_stream),
      m_diagnostic_stream(&diagnostic_stream),
      m_start_time(Clock::now()),
      m_last_output_time(m_start_time) {}

bool RuntimeConsoleReporter::enabled() const noexcept {
  return m_options.enabled;
}

bool RuntimeConsoleReporter::informationEnabled() const noexcept {
  return enabled() && isRoot() && !m_options.quiet;
}

bool RuntimeConsoleReporter::isRoot() const noexcept {
  return m_world_rank == 0;
}

const RuntimeConsoleOptions& RuntimeConsoleReporter::options() const noexcept {
  return m_options;
}

bool RuntimeConsoleReporter::shouldEmitStepForCadence(
    std::uint64_t step_index,
    std::optional<std::uint64_t> last_emitted_step,
    std::uint64_t status_every_steps,
    double seconds_since_last_output,
    double status_seconds) noexcept {
  const bool step_trigger_enabled = status_every_steps > 0U;
  const bool time_trigger_enabled = std::isfinite(status_seconds) && status_seconds > 0.0;
  if (!step_trigger_enabled && !time_trigger_enabled) {
    return false;
  }
  if (!last_emitted_step.has_value()) {
    return true;
  }
  if (step_index <= *last_emitted_step) {
    return false;
  }
  const bool step_due = step_trigger_enabled &&
      step_index - *last_emitted_step >= status_every_steps;
  const bool time_due = time_trigger_enabled &&
      std::isfinite(seconds_since_last_output) &&
      seconds_since_last_output >= status_seconds;
  return step_due || time_due;
}

void RuntimeConsoleReporter::emitRuntimePhase(
    std::string_view phase,
    std::string_view detail) {
  if (!informationEnabled()) {
    return;
  }
  std::ostringstream payload;
  payload << "phase=" << phase;
  if (!detail.empty()) {
    payload << " detail=" << detail;
  }
  emitInfoLine("RUNTIME", payload.str());
}

void RuntimeConsoleReporter::emitInitialConditions(
    std::uint64_t global_particle_count,
    std::uint64_t global_cell_count,
    bool restoring_from_restart,
    const std::filesystem::path& manifest_path) {
  if (!informationEnabled()) {
    return;
  }
  std::ostringstream payload;
  payload << "status=ready"
          << " particles=" << global_particle_count
          << " cells=" << global_cell_count
          << " source=" << (restoring_from_restart ? "restart" : "initial_conditions");
  if (!manifest_path.empty()) {
    payload << " manifest=" << std::quoted(manifest_path.string());
  }
  emitInfoLine("IC", payload.str());
}

void RuntimeConsoleReporter::emitStartup(
    const RuntimeConsoleStartupStatus& status) {
  if (!informationEnabled()) {
    return;
  }
  std::ostringstream payload;
  payload << "project=CHUI"
          << " version=" << core::versionString()
          << " config=" << std::quoted(m_options.config_path)
          << " run_name=" << std::quoted(status.run_name)
          << " run_directory=" << std::quoted(status.run_directory.string());
  if (!status.rank_directory.empty() && status.rank_directory != status.run_directory) {
    payload << " rank_directory=" << std::quoted(status.rank_directory.string());
  }
  payload << " mode=" << status.simulation_mode
          << " mpi_ranks=" << status.mpi_world_size
          << " openmp_compiled=" << (status.openmp_compiled ? "true" : "false")
          << " openmp_threads=" << status.openmp_threads
          << " build_mpi=" << (COSMOSIM_ENABLE_MPI ? "true" : "false")
          << " build_hdf5=" << (COSMOSIM_ENABLE_HDF5 ? "true" : "false")
          << " build_fftw=" << (COSMOSIM_ENABLE_FFTW ? "true" : "false")
          << " build_fftw_mpi="
          << ((COSMOSIM_ENABLE_MPI && COSMOSIM_ENABLE_FFTW) ? "true" : "false")
          << " particles=" << status.global_particle_count
          << " cells=" << status.global_cell_count
          << " t_code=" << formatDouble(status.t_code)
          << " dt_time_code=" << formatDouble(status.dt_time_code);
  if (status.a_scale.has_value()) {
    payload << " a=" << optionalDouble(status.a_scale)
            << " z=" << optionalDouble(status.redshift);
  }
  payload << " endpoint_t_code=" << formatDouble(status.endpoint_t_code);
  if (status.endpoint_a_scale.has_value()) {
    payload << " endpoint_a=" << optionalDouble(status.endpoint_a_scale);
  }
  if (status.pm_grid_nx > 0U && status.pm_grid_ny > 0U && status.pm_grid_nz > 0U) {
    payload << " pm_grid=" << status.pm_grid_nx << 'x'
            << status.pm_grid_ny << 'x' << status.pm_grid_nz
            << " pm_cadence_steps=" << status.pm_update_cadence_steps;
  }
  payload << " snapshot_interval_steps=" << status.snapshot_interval_steps
          << " snapshot_interval_time_code=" << formatDouble(status.snapshot_interval_time_code)
          << " restarts=" << (status.write_restarts ? "enabled" : "disabled")
          << " restart_state=" << (status.restoring_from_restart ? "resumed" : "fresh")
          << " memory_pressure="
          << (status.memory_pressure.empty() ? "unavailable" : status.memory_pressure)
          << " memory_headroom_bytes=" << optionalBytes(status.memory_headroom_bytes);
  emitInfoLine("START", payload.str());
}

void RuntimeConsoleReporter::emitPmSetup(
    std::size_t pm_grid_nx,
    std::size_t pm_grid_ny,
    std::size_t pm_grid_nz,
    int update_cadence_steps,
    std::string_view backend) {
  if (!informationEnabled()) {
    return;
  }
  std::ostringstream payload;
  payload << "status=ready"
          << " grid=" << pm_grid_nx << 'x' << pm_grid_ny << 'x' << pm_grid_nz
          << " cadence_steps=" << update_cadence_steps
          << " backend=" << backend;
  emitInfoLine("PM", payload.str());
}

void RuntimeConsoleReporter::emitMemory(
    const core::MemoryReport& report,
    std::string_view phase) {
  if (!informationEnabled()) {
    return;
  }
  std::ostringstream payload;
  payload << "phase=" << phase
          << " persistent_bytes=" << report.totals.persistent_total_bytes
          << " transient_bytes=" << report.totals.transient_total_bytes
          << " rank_max_owned_bytes=" << report.distributed.rank_max_owned_bytes
          << " pressure=" << memoryPressure(report)
          << " headroom_bytes=" << optionalBytes(memoryHeadroom(report));
  if (report.process_memory_reconciliation.has_value() &&
      report.process_memory_reconciliation->observed_rss_bytes.has_value()) {
    payload << " rss_bytes="
            << *report.process_memory_reconciliation->observed_rss_bytes;
  }
  emitInfoLine("MEMORY", payload.str());
}

void RuntimeConsoleReporter::emitStep(const RuntimeConsoleStepStatus& status) {
  if (!informationEnabled()) {
    return;
  }
  const double silence_seconds = secondsSinceLastOutput();
  if (!shouldEmitStepForCadence(
          status.step_index,
          m_last_emitted_step,
          m_options.status_every_steps,
          silence_seconds,
          m_options.status_seconds)) {
    return;
  }

  std::ostringstream payload;
  payload << "step=" << status.step_index
          << " t_code=" << formatDouble(status.t_code)
          << " dt_time_code=" << formatDouble(status.dt_time_code);
  if (status.a_scale.has_value()) {
    payload << " a=" << optionalDouble(status.a_scale)
            << " z=" << optionalDouble(status.redshift);
  }
  payload << " active_particles=" << status.active_particle_count << '/'
          << status.total_particle_count;
  if (status.total_cell_count > 0U || status.active_cell_count > 0U) {
    payload << " active_cells=" << status.active_cell_count << '/'
            << status.total_cell_count;
  }
  payload << " wall_step_s=" << formatSeconds(status.wall_step_seconds)
          << " wall_total_s=" << formatSeconds(secondsSince(m_start_time))
          << " pm=" << (status.pm_activity.empty() ? "none" : status.pm_activity)
          << " output=" << (status.output_activity.empty() ? "none" : status.output_activity)
          << " memory_pressure="
          << (status.memory_pressure.empty() ? "unavailable" : status.memory_pressure)
          << " memory_headroom_bytes=" << optionalBytes(status.memory_headroom_bytes);
  emitInfoLine("STEP", payload.str());
  m_last_emitted_step = status.step_index;
}

void RuntimeConsoleReporter::emitDecomposition(
    std::uint64_t step_index,
    std::uint64_t decomposition_epoch) {
  if (!informationEnabled()) {
    return;
  }
  std::ostringstream payload;
  payload << "step=" << step_index
          << " status=committed"
          << " epoch=" << decomposition_epoch;
  emitInfoLine("DECOMP", payload.str());
}

void RuntimeConsoleReporter::emitSnapshotCommitted(
    std::uint64_t step_index,
    const std::filesystem::path& member_path,
    const std::filesystem::path& set_path,
    std::uint32_t member_count) {
  if (!informationEnabled()) {
    return;
  }
  std::ostringstream payload;
  payload << "step=" << step_index
          << " state=committed"
          << " members=" << member_count
          << " manifest=" << std::quoted(set_path.string());
  if (member_count <= 1U) {
    payload << " file=" << std::quoted(member_path.string());
  } else {
    payload << " rank0_member=" << std::quoted(member_path.string());
  }
  emitInfoLine("SNAPSHOT", payload.str());
}

void RuntimeConsoleReporter::emitRestartCommitted(
    std::uint64_t step_index,
    const std::filesystem::path& restart_path) {
  if (!informationEnabled()) {
    return;
  }
  std::ostringstream payload;
  payload << "step=" << step_index
          << " status=verified"
          << " path=" << std::quoted(restart_path.string());
  emitInfoLine("RESTART", payload.str());
}

void RuntimeConsoleReporter::emitWarning(
    std::string_view subsystem,
    std::string_view message) {
  if (!enabled() || !isRoot()) {
    return;
  }
  std::ostringstream payload;
  payload << "subsystem=" << subsystem << " message=" << std::quoted(std::string(message));
  emitDiagnosticLine("WARN", payload.str());
}

void RuntimeConsoleReporter::emitDone(
    std::uint64_t completed_steps,
    double final_time_code,
    std::optional<double> final_scale_factor,
    std::optional<double> final_redshift,
    const std::filesystem::path& run_directory,
    const std::filesystem::path& rank_directory,
    const std::filesystem::path& normalized_config_path,
    const std::filesystem::path& operational_report_path) {
  if (!informationEnabled()) {
    return;
  }
  std::ostringstream payload;
  payload << "completed_steps=" << completed_steps
          << " final_t_code=" << formatDouble(final_time_code);
  if (final_scale_factor.has_value()) {
    payload << " final_a=" << optionalDouble(final_scale_factor)
            << " final_z=" << optionalDouble(final_redshift);
  }
  payload << " wall_total_s=" << formatSeconds(secondsSince(m_start_time))
          << " run_directory=" << std::quoted(run_directory.string());
  if (!rank_directory.empty() && rank_directory != run_directory) {
    payload << " rank_directory=" << std::quoted(rank_directory.string());
  }
  payload << " normalized_config=" << std::quoted(normalized_config_path.string())
          << " operational_report=" << std::quoted(operational_report_path.string());
  emitInfoLine("DONE", payload.str());
}

void RuntimeConsoleReporter::emitInfoLine(
    std::string_view tag,
    std::string_view payload) {
  if (!informationEnabled() || m_info_stream == nullptr) {
    return;
  }
  *m_info_stream << "[CHUI][" << tag << "] " << payload << '\n' << std::flush;
  m_last_output_time = Clock::now();
}

void RuntimeConsoleReporter::emitDiagnosticLine(
    std::string_view tag,
    std::string_view payload) {
  if (!enabled() || !isRoot() || m_diagnostic_stream == nullptr) {
    return;
  }
  *m_diagnostic_stream << "[CHUI][" << tag << "] " << payload << '\n' << std::flush;
  m_last_output_time = Clock::now();
}

double RuntimeConsoleReporter::secondsSince(Clock::time_point start) const noexcept {
  return std::chrono::duration<double>(Clock::now() - start).count();
}

double RuntimeConsoleReporter::secondsSinceLastOutput() const noexcept {
  return secondsSince(m_last_output_time);
}

}  // namespace cosmosim::workflows
