#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <string>
#include <vector>

#include "cosmosim/core/config.hpp"
#include "cosmosim/core/simulation_state.hpp"

namespace cosmosim::physics {

enum class StellarYieldChannel : std::uint8_t {
  kWinds = 0,
  kCoreCollapseSn = 1,
  kTypeIaSn = 2,
};

constexpr std::size_t k_stellar_yield_channel_count = 3;

struct StellarEvolutionCumulative {
  double return_fraction_total = 0.0;
  double metal_yield_fraction_total = 0.0;
  double energy_erg_per_initial_mass_code = 0.0;
  std::array<double, k_stellar_yield_channel_count> return_fraction_channel{};
  std::array<double, k_stellar_yield_channel_count> metal_yield_fraction_channel{};
  std::array<double, k_stellar_yield_channel_count> energy_erg_per_initial_mass_channel{};
};

struct StellarEvolutionIntervalBudget {
  double returned_mass_code = 0.0;
  double returned_metals_code = 0.0;
  double feedback_energy_erg = 0.0;
  double remnant_change_code = 0.0;
  std::array<double, k_stellar_yield_channel_count> returned_mass_channel_code{};
  std::array<double, k_stellar_yield_channel_count> returned_metals_channel_code{};
  std::array<double, k_stellar_yield_channel_count> feedback_energy_channel_erg{};
};

struct StellarEvolutionTable {
  std::string table_id = "builtin_reference";
  std::string table_version = "v1";
  std::string source_path = "builtin";
  std::vector<double> age_yr;
  std::vector<double> return_fraction_total;
  std::vector<double> metal_yield_fraction_total;
  std::vector<double> energy_erg_per_initial_mass_code;
  std::array<std::vector<double>, k_stellar_yield_channel_count> return_fraction_channel;
  std::array<std::vector<double>, k_stellar_yield_channel_count> metal_yield_fraction_channel;
  std::array<std::vector<double>, k_stellar_yield_channel_count> energy_erg_per_initial_mass_channel;

  [[nodiscard]] bool isConsistent() const noexcept;
  [[nodiscard]] StellarEvolutionCumulative evaluateAtAgeYears(double age_yr_value) const;
  [[nodiscard]] StellarEvolutionIntervalBudget integrateInterval(
      double age_begin_yr,
      double age_end_yr,
      double initial_mass_code) const;

  [[nodiscard]] static StellarEvolutionTable loadFromTextFile(
      const std::string& path,
      const std::string& source_tag = "stellar_evolution_table");
  [[nodiscard]] static StellarEvolutionTable makeBuiltinReference();
};

struct StellarEvolutionConfig {
  bool enabled = true;
  std::string table_path;
  std::uint32_t metadata_schema_version = 1;
  double hubble_time_years = 1.44e10;
};

struct StellarEvolutionStarBudget {
  std::uint32_t star_index = 0;
  std::uint32_t particle_index = 0;
  double star_age_begin_years = 0.0;
  double star_age_end_years = 0.0;
  double mass_old_code = 0.0;
  double mass_new_code = 0.0;
  StellarEvolutionIntervalBudget interval;
};

struct StellarEvolutionStepCounters {
  std::uint64_t scanned_stars = 0;
  std::uint64_t evolved_stars = 0;
  double returned_mass_code = 0.0;
  double returned_metals_code = 0.0;
  double feedback_energy_erg = 0.0;
};

struct StellarEvolutionStepReport {
  StellarEvolutionStepCounters counters;
  std::vector<StellarEvolutionStarBudget> budgets;
};

class StellarEvolutionBookkeeper {
 public:
  StellarEvolutionBookkeeper(StellarEvolutionConfig config, StellarEvolutionTable table);

  [[nodiscard]] const StellarEvolutionConfig& config() const noexcept;
  [[nodiscard]] const StellarEvolutionTable& table() const noexcept;

  [[nodiscard]] double evaluateStarAgeYears(double formation_scale_factor, double current_scale_factor) const;

  [[nodiscard]] StellarEvolutionStepReport apply(
      core::SimulationState& state,
      std::span<const std::uint32_t> active_star_indices,
      double current_scale_factor,
      double dt_code) const;

  [[nodiscard]] core::ModuleSidecarBlock buildMetadataSidecar(const StellarEvolutionStepReport& report) const;

 private:
  StellarEvolutionConfig m_config;
  StellarEvolutionTable m_table;
};

[[nodiscard]] StellarEvolutionConfig makeStellarEvolutionConfig(const core::PhysicsConfig& physics_config);
[[nodiscard]] StellarEvolutionTable loadStellarEvolutionTable(const core::PhysicsConfig& physics_config);

}  // namespace cosmosim::physics
