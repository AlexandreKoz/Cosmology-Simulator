#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

#include "cosmosim/core/config.hpp"
#include "cosmosim/core/time_integration.hpp"

namespace cosmosim::analysis {

enum class DiagnosticClass : std::uint8_t {
  kRunHealth = 0,
  kScienceLight = 1,
  kScienceHeavy = 2,
};

struct PowerSpectrumBin {
  double k_center_code = 0.0;
  double power_code_volume = 0.0;
  std::uint64_t mode_count = 0;
};

struct StarFormationHistoryBin {
  double scale_factor_center = 0.0;
  double formed_mass_code = 0.0;
};

struct AngularMomentumBudget {
  std::array<double, 3> total_l_code{};
  std::array<double, 3> gas_l_code{};
  std::array<double, 3> star_l_code{};
  std::array<double, 3> dark_matter_l_code{};
  std::array<double, 3> black_hole_l_code{};
};

struct RunHealthCounters {
  std::uint64_t particle_count = 0;
  std::uint64_t cell_count = 0;
  std::uint64_t star_count = 0;
  bool ownership_invariants_ok = false;
  bool unique_particle_ids_ok = false;
  std::uint64_t non_finite_particles = 0;
  std::uint64_t non_finite_cells = 0;
};

struct DiagnosticsBundle {
  std::uint64_t step_index = 0;
  double scale_factor = 1.0;
  DiagnosticClass diagnostic_class = DiagnosticClass::kRunHealth;
  RunHealthCounters health;
  std::vector<PowerSpectrumBin> power_spectrum;
  std::vector<StarFormationHistoryBin> star_formation_history;
  AngularMomentumBudget angular_momentum;
  std::vector<double> xy_slice_density_code;
  std::vector<double> xy_projection_density_code;
  std::size_t quicklook_grid_n = 0;
  std::string quicklook_projection_csv_path;
};

struct DiagnosticsTiming {
  double cumulative_run_health_ms = 0.0;
  double cumulative_light_ms = 0.0;
  double cumulative_heavy_ms = 0.0;
  std::uint64_t run_health_calls = 0;
  std::uint64_t light_calls = 0;
  std::uint64_t heavy_calls = 0;
};

class DiagnosticsEngine {
 public:
  explicit DiagnosticsEngine(core::SimulationConfig config);

  [[nodiscard]] RunHealthCounters computeRunHealth(const core::SimulationState& state) const;

  [[nodiscard]] std::vector<PowerSpectrumBin> computePowerSpectrum(
      const core::SimulationState& state,
      std::size_t mesh_n,
      std::size_t bin_count) const;

  [[nodiscard]] std::vector<StarFormationHistoryBin> computeStarFormationHistory(
      const core::SimulationState& state,
      std::size_t bin_count) const;

  [[nodiscard]] AngularMomentumBudget computeAngularMomentumBudget(
      const core::SimulationState& state) const;

  [[nodiscard]] std::vector<double> computeGasXySliceDensity(
      const core::SimulationState& state,
      std::size_t grid_n) const;

  [[nodiscard]] std::vector<double> computeGasXyProjectionDensity(
      const core::SimulationState& state,
      std::size_t grid_n) const;

  [[nodiscard]] DiagnosticsBundle generateBundle(
      const core::SimulationState& state,
      std::uint64_t step_index,
      double scale_factor,
      DiagnosticClass diagnostic_class) const;

  void writeBundle(const DiagnosticsBundle& bundle) const;
  void enforceRetentionPolicy() const;

 private:
  [[nodiscard]] std::filesystem::path diagnosticsOutputDirectory() const;
  [[nodiscard]] std::filesystem::path bundlePath(const DiagnosticsBundle& bundle) const;
  [[nodiscard]] std::filesystem::path quicklookPath(const DiagnosticsBundle& bundle) const;

  core::SimulationConfig m_config;
};

class DiagnosticsCallback final : public core::IntegrationCallback {
 public:
  explicit DiagnosticsCallback(core::SimulationConfig config);

  [[nodiscard]] std::string_view callbackName() const override;
  void onStage(core::StepContext& context) override;

  [[nodiscard]] const DiagnosticsTiming& timing() const noexcept;

 private:
  void runDiagnostics(core::StepContext& context, DiagnosticClass diagnostic_class);

  DiagnosticsEngine m_engine;
  core::SimulationConfig m_config;
  DiagnosticsTiming m_timing;
};

}  // namespace cosmosim::analysis
