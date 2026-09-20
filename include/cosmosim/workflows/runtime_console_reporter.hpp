#pragma once

#include <chrono>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <iosfwd>
#include <optional>
#include <string>
#include <string_view>

namespace cosmosim::core {
struct MemoryReport;
}

namespace cosmosim::workflows {

// Operational presentation controls only. These values are deliberately kept
// outside the persistent simulation configuration and never participate in
// scientific/restart authority.
struct RuntimeConsoleOptions {
  bool enabled = false;
  bool quiet = false;
  std::uint64_t status_every_steps = 10U;
  double status_seconds = 30.0;
  std::string config_path;
};

struct RuntimeConsoleStartupStatus {
  std::string run_name;
  std::filesystem::path run_directory;
  std::string simulation_mode;
  int mpi_world_size = 1;
  bool openmp_compiled = false;
  int openmp_threads = 1;
  std::uint64_t global_particle_count = 0U;
  std::uint64_t global_cell_count = 0U;
  double t_code = 0.0;
  double dt_time_code = 0.0;
  std::optional<double> a_scale;
  std::optional<double> redshift;
  double endpoint_t_code = 0.0;
  std::optional<double> endpoint_a_scale;
  std::size_t pm_grid_nx = 0U;
  std::size_t pm_grid_ny = 0U;
  std::size_t pm_grid_nz = 0U;
  int pm_update_cadence_steps = 0;
  int snapshot_interval_steps = 0;
  double snapshot_interval_time_code = 0.0;
  bool write_restarts = false;
  bool restoring_from_restart = false;
  std::optional<std::uint64_t> memory_headroom_bytes;
  std::string memory_pressure;
};

struct RuntimeConsoleStepStatus {
  std::uint64_t step_index = 0U;
  double t_code = 0.0;
  double dt_time_code = 0.0;
  std::optional<double> a_scale;
  std::optional<double> redshift;
  std::uint64_t active_particle_count = 0U;
  std::uint64_t total_particle_count = 0U;
  std::uint64_t active_cell_count = 0U;
  std::uint64_t total_cell_count = 0U;
  double wall_step_seconds = 0.0;
  std::string pm_activity;
  std::string output_activity;
  std::optional<std::uint64_t> memory_headroom_bytes;
  std::string memory_pressure;
};

class RuntimeConsoleReporter final {
 public:
  using Clock = std::chrono::steady_clock;

  RuntimeConsoleReporter(
      RuntimeConsoleOptions options,
      int world_rank,
      std::ostream& info_stream,
      std::ostream& diagnostic_stream);

  [[nodiscard]] bool enabled() const noexcept;
  [[nodiscard]] bool informationEnabled() const noexcept;
  [[nodiscard]] bool isRoot() const noexcept;
  [[nodiscard]] const RuntimeConsoleOptions& options() const noexcept;

  // Pure cadence predicate used by unit tests and the live reporter. A first
  // productive record is emitted when at least one cadence trigger is enabled.
  [[nodiscard]] static bool shouldEmitStepForCadence(
      std::uint64_t step_index,
      std::optional<std::uint64_t> last_emitted_step,
      std::uint64_t status_every_steps,
      double seconds_since_last_output,
      double status_seconds) noexcept;

  void emitRuntimePhase(std::string_view phase, std::string_view detail = {});
  void emitInitialConditions(
      std::uint64_t global_particle_count,
      std::uint64_t global_cell_count,
      bool restoring_from_restart,
      const std::filesystem::path& manifest_path);
  void emitStartup(const RuntimeConsoleStartupStatus& status);
  void emitPmSetup(
      std::size_t pm_grid_nx,
      std::size_t pm_grid_ny,
      std::size_t pm_grid_nz,
      int update_cadence_steps,
      std::string_view backend);
  void emitMemory(const core::MemoryReport& report, std::string_view phase);
  void emitStep(const RuntimeConsoleStepStatus& status);
  void emitDecomposition(std::uint64_t step_index, std::uint64_t decomposition_epoch);
  void emitSnapshotCommitted(
      std::uint64_t step_index,
      const std::filesystem::path& member_path,
      const std::filesystem::path& set_path);
  void emitRestartCommitted(
      std::uint64_t step_index,
      const std::filesystem::path& restart_path);
  void emitWarning(std::string_view subsystem, std::string_view message);
  void emitDone(
      std::uint64_t completed_steps,
      double final_time_code,
      std::optional<double> final_scale_factor,
      std::optional<double> final_redshift,
      const std::filesystem::path& run_directory,
      const std::filesystem::path& normalized_config_path,
      const std::filesystem::path& operational_report_path);

 private:
  void emitInfoLine(std::string_view tag, std::string_view payload);
  void emitDiagnosticLine(std::string_view tag, std::string_view payload);
  [[nodiscard]] double secondsSince(Clock::time_point start) const noexcept;
  [[nodiscard]] double secondsSinceLastOutput() const noexcept;

  RuntimeConsoleOptions m_options;
  int m_world_rank = 0;
  std::ostream* m_info_stream = nullptr;
  std::ostream* m_diagnostic_stream = nullptr;
  Clock::time_point m_start_time;
  Clock::time_point m_last_output_time;
  std::optional<std::uint64_t> m_last_emitted_step;
};

}  // namespace cosmosim::workflows
