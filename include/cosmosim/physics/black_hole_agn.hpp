#pragma once

#include <cstdint>
#include <span>
#include <string_view>
#include <vector>

#include "cosmosim/core/config.hpp"
#include "cosmosim/core/constants.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/core/time_integration.hpp"

namespace cosmosim::physics {

struct BlackHoleAgnConfig {
  bool enabled = false;
  double seed_halo_mass_threshold_code = 1.0e3;
  double seed_mass_code = 1.0;
  std::uint32_t seed_max_per_cell = 1;
  double alpha_bondi = 1.0;
  bool use_eddington_cap = true;
  double epsilon_r = 0.1;
  double epsilon_f = 0.05;
  double feedback_coupling_efficiency = 1.0;
  double duty_cycle_active_edd_ratio_threshold = 0.01;
  double proton_mass_si = 1.67262192369e-27;
  double thomson_cross_section_si = 6.6524587321e-29;
  double newton_g_si = core::constants::k_newton_g_si;
  double speed_of_light_si = core::constants::k_speed_of_light_si;
  std::uint32_t metadata_schema_version = 1;
};

struct BlackHoleSeedCandidate {
  std::uint32_t cell_index = 0;
  double host_halo_mass_code = 0.0;
  std::uint32_t owning_rank = 0;
};

struct BlackHoleRates {
  double mdot_bondi_code = 0.0;
  double mdot_edd_code = 0.0;
  double mdot_acc_code = 0.0;
  double eddington_ratio = 0.0;
  double feedback_power_code = 0.0;
};

struct BlackHoleAgnCounters {
  std::uint64_t scanned_bh = 0;
  std::uint64_t active_bh = 0;
  std::uint64_t seed_candidates = 0;
  std::uint64_t seeded_bh = 0;
  double integrated_accreted_mass_code = 0.0;
  double integrated_feedback_energy_code = 0.0;
  double deposited_feedback_energy_code = 0.0;
  double integrated_duty_cycle_active_time_code = 0.0;
  double integrated_duty_cycle_total_time_code = 0.0;
};

struct BlackHoleAgnStepReport {
  BlackHoleAgnCounters counters;
  std::vector<std::uint32_t> seeded_cell_indices;
};

class BlackHoleAgnModel {
 public:
  explicit BlackHoleAgnModel(BlackHoleAgnConfig config);

  [[nodiscard]] const BlackHoleAgnConfig& config() const noexcept;
  [[nodiscard]] bool isSeedEligible(
      const core::SimulationState& state,
      const BlackHoleSeedCandidate& candidate) const;
  [[nodiscard]] BlackHoleRates computeAccretionRates(
      double bh_mass_code,
      double gas_density_code,
      double sound_speed_code,
      double relative_velocity_code) const;

  [[nodiscard]] BlackHoleAgnStepReport apply(
      core::SimulationState& state,
      std::span<const BlackHoleSeedCandidate> seed_candidates,
      double dt_code,
      std::uint64_t step_index) const;

  [[nodiscard]] core::ModuleSidecarBlock buildMetadataSidecar(const BlackHoleAgnCounters& counters) const;

 private:
  BlackHoleAgnConfig m_config;
};

[[nodiscard]] BlackHoleAgnConfig makeBlackHoleAgnConfig(const core::PhysicsConfig& physics_config);

class BlackHoleAgnCallback final : public core::IntegrationCallback {
 public:
  explicit BlackHoleAgnCallback(BlackHoleAgnModel model);

  [[nodiscard]] std::string_view callbackName() const override;
  void onStage(core::StepContext& context) override;

  void setSeedCandidates(std::span<const BlackHoleSeedCandidate> seed_candidates);
  [[nodiscard]] const BlackHoleAgnStepReport& lastStepReport() const noexcept;

 private:
  BlackHoleAgnModel m_model;
  std::vector<BlackHoleSeedCandidate> m_seed_candidates;
  BlackHoleAgnStepReport m_last_step_report;
};

}  // namespace cosmosim::physics
