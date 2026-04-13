#pragma once

#include <cstddef>
#include <cstdint>
#include <span>
#include <string>
#include <vector>

#include "cosmosim/core/config.hpp"
#include "cosmosim/core/simulation_state.hpp"

namespace cosmosim::physics {

enum class StellarFeedbackMode : std::uint8_t {
  kThermal = 0,
  kKinetic = 1,
  kMomentum = 2,
  kThermalKineticMomentum = 3,
};

enum class StellarFeedbackVariant : std::uint8_t {
  kNone = 0,
  kDelayedCooling = 1,
  kStochastic = 2,
};

struct StellarFeedbackConfig {
  bool enabled = true;
  StellarFeedbackMode mode = StellarFeedbackMode::kThermalKineticMomentum;
  StellarFeedbackVariant variant = StellarFeedbackVariant::kNone;
  bool use_returned_mass_budget = true;
  double epsilon_thermal = 0.6;
  double epsilon_kinetic = 0.3;
  double epsilon_momentum = 0.1;
  double sn_energy_erg_per_mass_code = 1.0e49;
  double momentum_code_per_mass_code = 3.0e3;
  std::uint32_t neighbor_count = 8;
  double delayed_cooling_time_code = 0.0;
  double stochastic_event_probability = 0.25;
  std::uint64_t random_seed = 42424242ull;
  std::uint32_t metadata_schema_version = 1;
};

struct StellarFeedbackBudget {
  double source_mass_code = 0.0;
  double returned_mass_code = 0.0;
  double returned_metals_code = 0.0;
  double total_energy_erg = 0.0;
  double thermal_energy_erg = 0.0;
  double kinetic_energy_erg = 0.0;
  double momentum_budget_code = 0.0;
};

struct StellarFeedbackTarget {
  std::uint32_t cell_index = 0;
  double weight = 0.0;
  double radial_dx_comoving = 0.0;
  double radial_dy_comoving = 0.0;
  double radial_dz_comoving = 0.0;
};

struct StellarFeedbackStarReport {
  std::uint32_t star_index = 0;
  std::uint32_t particle_index = 0;
  StellarFeedbackBudget budget;
  std::size_t target_count = 0;
  bool stochastic_event_fired = true;
  bool delayed_cooling_applied = false;
  double deposited_mass_code = 0.0;
  double deposited_metals_code = 0.0;
  double deposited_thermal_energy_erg = 0.0;
  double deposited_kinetic_energy_erg = 0.0;
  double deposited_momentum_code = 0.0;
  double unresolved_mass_code = 0.0;
  double unresolved_metals_code = 0.0;
  double unresolved_thermal_energy_erg = 0.0;
  double unresolved_kinetic_energy_erg = 0.0;
  double unresolved_momentum_code = 0.0;
};

struct StellarFeedbackStepCounters {
  std::uint64_t scanned_stars = 0;
  std::uint64_t feedback_stars = 0;
  std::uint64_t target_cells_visited = 0;
  double source_mass_code = 0.0;
  double deposited_mass_code = 0.0;
  double deposited_metals_code = 0.0;
  double deposited_thermal_energy_erg = 0.0;
  double deposited_kinetic_energy_erg = 0.0;
  double deposited_momentum_code = 0.0;
  double unresolved_mass_code = 0.0;
  double unresolved_metals_code = 0.0;
  double unresolved_thermal_energy_erg = 0.0;
  double unresolved_kinetic_energy_erg = 0.0;
  double unresolved_momentum_code = 0.0;
};

struct StellarFeedbackStepReport {
  StellarFeedbackStepCounters counters;
  std::vector<StellarFeedbackStarReport> star_reports;
};

struct StellarFeedbackModuleState {
  // Sidecar-style bookkeeping kept outside core hot particle layout.
  std::vector<double> last_returned_mass_cumulative_code;
  std::vector<double> carry_mass_code;
  std::vector<double> carry_metals_code;
  std::vector<double> carry_thermal_energy_erg;
  std::vector<double> carry_kinetic_energy_erg;
  std::vector<double> carry_momentum_code;

  void ensureStarCapacity(std::size_t star_count);
};

class StellarFeedbackModel {
 public:
  explicit StellarFeedbackModel(StellarFeedbackConfig config);

  [[nodiscard]] const StellarFeedbackConfig& config() const noexcept;
  [[nodiscard]] StellarFeedbackBudget computeBudget(
      double source_mass_code,
      double returned_mass_code,
      double returned_metals_code) const;

  [[nodiscard]] std::vector<StellarFeedbackTarget> selectTargets(
      const core::SimulationState& state,
      std::uint32_t particle_index) const;

  [[nodiscard]] StellarFeedbackStepReport apply(
      core::SimulationState& state,
      StellarFeedbackModuleState& module_state,
      std::span<const std::uint32_t> active_star_indices,
      std::span<const double> returned_mass_delta_code,
      std::span<const double> returned_metals_delta_code,
      double dt_code) const;

  [[nodiscard]] core::ModuleSidecarBlock buildMetadataSidecar(const StellarFeedbackStepReport& report) const;

 private:
  [[nodiscard]] static std::string modeToString(StellarFeedbackMode mode);
  [[nodiscard]] static std::string variantToString(StellarFeedbackVariant variant);
  [[nodiscard]] bool stochasticEventFires(std::uint32_t star_index, std::uint64_t step_seed) const;

  StellarFeedbackConfig m_config;
};

[[nodiscard]] StellarFeedbackConfig makeStellarFeedbackConfig(const core::PhysicsConfig& physics_config);

}  // namespace cosmosim::physics
