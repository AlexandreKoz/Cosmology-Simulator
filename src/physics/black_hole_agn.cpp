#include "cosmosim/physics/black_hole_agn.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

#include "cosmosim/core/constants.hpp"

namespace cosmosim::physics {
namespace {

constexpr double k_mass_floor = 1.0e-20;
constexpr double k_density_floor = 1.0e-20;
constexpr double k_speed_floor = 1.0e-20;
constexpr double k_time_floor = 1.0e-30;

[[nodiscard]] std::uint64_t nextParticleId(const core::SimulationState& state) {
  std::uint64_t max_id = 0;
  for (const std::uint64_t particle_id : state.particle_sidecar.particle_id) {
    max_id = std::max(max_id, particle_id);
  }
  return max_id + 1;
}

[[nodiscard]] std::uint32_t countBhInCell(const core::SimulationState& state, std::uint32_t cell_index) {
  std::uint32_t count = 0;
  for (std::size_t bh_index = 0; bh_index < state.black_holes.size(); ++bh_index) {
    if (state.black_holes.host_cell_index[bh_index] == cell_index) {
      ++count;
    }
  }
  return count;
}

}  // namespace

BlackHoleAgnModel::BlackHoleAgnModel(BlackHoleAgnConfig config) : m_config(std::move(config)) {
  if (m_config.seed_halo_mass_threshold_code <= 0.0 || m_config.seed_mass_code <= 0.0) {
    throw std::invalid_argument("BlackHoleAgnModel: seed thresholds and masses must be > 0");
  }
  if (m_config.alpha_bondi <= 0.0) {
    throw std::invalid_argument("BlackHoleAgnModel: alpha_bondi must be > 0");
  }
  if (m_config.epsilon_r <= 0.0 || m_config.epsilon_r > 1.0 || m_config.epsilon_f < 0.0 ||
      m_config.epsilon_f > 1.0 || m_config.feedback_coupling_efficiency < 0.0 ||
      m_config.feedback_coupling_efficiency > 1.0) {
    throw std::invalid_argument("BlackHoleAgnModel: efficiencies must be in physically conservative bounds");
  }
  if (m_config.proton_mass_si <= 0.0 || m_config.thomson_cross_section_si <= 0.0 ||
      m_config.newton_g_si <= 0.0 || m_config.speed_of_light_si <= 0.0) {
    throw std::invalid_argument("BlackHoleAgnModel: constants must be > 0");
  }
}

const BlackHoleAgnConfig& BlackHoleAgnModel::config() const noexcept { return m_config; }

bool BlackHoleAgnModel::isSeedEligible(
    const core::SimulationState& state,
    const BlackHoleSeedCandidate& candidate) const {
  if (!m_config.enabled || candidate.cell_index >= state.cells.size()) {
    return false;
  }
  if (candidate.host_halo_mass_code < m_config.seed_halo_mass_threshold_code) {
    return false;
  }
  return countBhInCell(state, candidate.cell_index) < m_config.seed_max_per_cell;
}

BlackHoleRates BlackHoleAgnModel::computeAccretionRates(
    double bh_mass_code,
    double gas_density_code,
    double sound_speed_code,
    double relative_velocity_code) const {
  BlackHoleRates rates;
  const double mass = std::max(bh_mass_code, k_mass_floor);
  const double rho = std::max(gas_density_code, k_density_floor);
  const double cs2 = sound_speed_code * sound_speed_code;
  const double v2 = relative_velocity_code * relative_velocity_code;
  const double denom = std::pow(std::max(cs2 + v2, k_speed_floor * k_speed_floor), 1.5);

  rates.mdot_bondi_code =
      m_config.alpha_bondi * 4.0 * core::constants::k_pi * m_config.newton_g_si * m_config.newton_g_si *
      mass * mass * rho / denom;

  rates.mdot_edd_code =
      4.0 * core::constants::k_pi * m_config.newton_g_si * mass * m_config.proton_mass_si /
      (m_config.epsilon_r * m_config.thomson_cross_section_si * m_config.speed_of_light_si);

  rates.mdot_acc_code = m_config.use_eddington_cap ? std::min(rates.mdot_bondi_code, rates.mdot_edd_code)
                                                    : rates.mdot_bondi_code;
  rates.eddington_ratio = rates.mdot_acc_code / std::max(rates.mdot_edd_code, k_time_floor);
  rates.feedback_power_code = m_config.epsilon_f * m_config.epsilon_r * rates.mdot_acc_code *
                              m_config.speed_of_light_si * m_config.speed_of_light_si;
  return rates;
}

BlackHoleAgnStepReport BlackHoleAgnModel::apply(
    core::SimulationState& state,
    std::span<const BlackHoleSeedCandidate> seed_candidates,
    double dt_code,
    std::uint64_t /*step_index*/) const {
  BlackHoleAgnStepReport report;
  if (!m_config.enabled || dt_code <= 0.0) {
    state.sidecars.upsert(buildMetadataSidecar(report.counters));
    return report;
  }

  for (std::size_t bh_local_index = 0; bh_local_index < state.black_holes.size(); ++bh_local_index) {
    ++report.counters.scanned_bh;
    const std::uint32_t particle_index = state.black_holes.particle_index[bh_local_index];
    const std::uint32_t host_cell_index = state.black_holes.host_cell_index[bh_local_index];
    if (particle_index >= state.particles.size() || host_cell_index >= state.gas_cells.size()) {
      continue;
    }

    ++report.counters.active_bh;
    const BlackHoleRates rates = computeAccretionRates(
        state.black_holes.subgrid_mass_code[bh_local_index],
        state.gas_cells.density_code[host_cell_index],
        state.gas_cells.sound_speed_code[host_cell_index],
        0.0);

    const double delta_mass = std::max(rates.mdot_acc_code * dt_code, 0.0);
    const double delta_feedback_energy = std::max(rates.feedback_power_code * dt_code, 0.0);
    const double delta_feedback_deposited = delta_feedback_energy * m_config.feedback_coupling_efficiency;

    state.black_holes.accretion_rate_code[bh_local_index] = rates.mdot_acc_code;
    state.black_holes.feedback_energy_code[bh_local_index] = delta_feedback_energy;
    state.black_holes.eddington_ratio[bh_local_index] = rates.eddington_ratio;
    state.black_holes.subgrid_mass_code[bh_local_index] += delta_mass;
    state.black_holes.cumulative_accreted_mass_code[bh_local_index] += delta_mass;
    state.black_holes.cumulative_feedback_energy_code[bh_local_index] += delta_feedback_energy;
    state.black_holes.duty_cycle_total_time_code[bh_local_index] += dt_code;

    if (rates.eddington_ratio >= m_config.duty_cycle_active_edd_ratio_threshold) {
      state.black_holes.duty_cycle_active_time_code[bh_local_index] += dt_code;
      report.counters.integrated_duty_cycle_active_time_code += dt_code;
    }
    report.counters.integrated_duty_cycle_total_time_code += dt_code;

    // We deposit feedback energy to host-cell internal energy in the same frame as internal_energy_code.
    if (host_cell_index < state.gas_cells.internal_energy_code.size()) {
      state.gas_cells.internal_energy_code[host_cell_index] += delta_feedback_deposited;
    }

    report.counters.integrated_accreted_mass_code += delta_mass;
    report.counters.integrated_feedback_energy_code += delta_feedback_energy;
    report.counters.deposited_feedback_energy_code += delta_feedback_deposited;

    // Keep gravity-hot particle mass synchronized with BH subgrid mass for restart consistency.
    state.particles.mass_code[particle_index] = state.black_holes.subgrid_mass_code[bh_local_index];
  }

  std::uint64_t next_particle_id = nextParticleId(state);
  for (const BlackHoleSeedCandidate& candidate : seed_candidates) {
    ++report.counters.seed_candidates;
    if (!isSeedEligible(state, candidate)) {
      continue;
    }

    const std::size_t particle_index = state.particles.size();
    state.resizeParticles(particle_index + 1);
    state.particles.position_x_comoving[particle_index] = state.cells.center_x_comoving[candidate.cell_index];
    state.particles.position_y_comoving[particle_index] = state.cells.center_y_comoving[candidate.cell_index];
    state.particles.position_z_comoving[particle_index] = state.cells.center_z_comoving[candidate.cell_index];
    state.particles.velocity_x_peculiar[particle_index] = 0.0;
    state.particles.velocity_y_peculiar[particle_index] = 0.0;
    state.particles.velocity_z_peculiar[particle_index] = 0.0;
    state.particles.mass_code[particle_index] = m_config.seed_mass_code;
    state.particles.time_bin[particle_index] = state.cells.time_bin[candidate.cell_index];

    state.particle_sidecar.particle_id[particle_index] = next_particle_id++;
    state.particle_sidecar.sfc_key[particle_index] = state.particle_sidecar.particle_id[particle_index];
    state.particle_sidecar.species_tag[particle_index] =
        static_cast<std::uint32_t>(core::ParticleSpecies::kBlackHole);
    state.particle_sidecar.particle_flags[particle_index] = 0;
    state.particle_sidecar.owning_rank[particle_index] = candidate.owning_rank;

    const std::size_t bh_local_index = state.black_holes.size();
    state.black_holes.resize(bh_local_index + 1);
    state.black_holes.particle_index[bh_local_index] = static_cast<std::uint32_t>(particle_index);
    state.black_holes.host_cell_index[bh_local_index] = candidate.cell_index;
    state.black_holes.subgrid_mass_code[bh_local_index] = m_config.seed_mass_code;
    state.black_holes.accretion_rate_code[bh_local_index] = 0.0;
    state.black_holes.feedback_energy_code[bh_local_index] = 0.0;
    state.black_holes.eddington_ratio[bh_local_index] = 0.0;
    state.black_holes.cumulative_accreted_mass_code[bh_local_index] = 0.0;
    state.black_holes.cumulative_feedback_energy_code[bh_local_index] = 0.0;
    state.black_holes.duty_cycle_active_time_code[bh_local_index] = 0.0;
    state.black_holes.duty_cycle_total_time_code[bh_local_index] = 0.0;

    ++report.counters.seeded_bh;
    report.seeded_cell_indices.push_back(candidate.cell_index);
  }

  state.species.count_by_species[static_cast<std::size_t>(core::ParticleSpecies::kBlackHole)] +=
      report.counters.seeded_bh;
  state.rebuildSpeciesIndex();
  state.sidecars.upsert(buildMetadataSidecar(report.counters));
  return report;
}

core::ModuleSidecarBlock BlackHoleAgnModel::buildMetadataSidecar(const BlackHoleAgnCounters& counters) const {
  std::ostringstream stream;
  stream << "module=black_hole_agn\n";
  stream << "schema_version=" << m_config.metadata_schema_version << "\n";
  stream << "seed_halo_mass_threshold_code=" << m_config.seed_halo_mass_threshold_code << "\n";
  stream << "seed_mass_code=" << m_config.seed_mass_code << "\n";
  stream << "seed_max_per_cell=" << m_config.seed_max_per_cell << "\n";
  stream << "alpha_bondi=" << m_config.alpha_bondi << "\n";
  stream << "use_eddington_cap=" << (m_config.use_eddington_cap ? "true" : "false") << "\n";
  stream << "epsilon_r=" << m_config.epsilon_r << "\n";
  stream << "epsilon_f=" << m_config.epsilon_f << "\n";
  stream << "feedback_coupling_efficiency=" << m_config.feedback_coupling_efficiency << "\n";
  stream << "duty_cycle_active_edd_ratio_threshold=" << m_config.duty_cycle_active_edd_ratio_threshold << "\n";
  stream << "scanned_bh=" << counters.scanned_bh << "\n";
  stream << "active_bh=" << counters.active_bh << "\n";
  stream << "seed_candidates=" << counters.seed_candidates << "\n";
  stream << "seeded_bh=" << counters.seeded_bh << "\n";
  stream << "integrated_accreted_mass_code=" << counters.integrated_accreted_mass_code << "\n";
  stream << "integrated_feedback_energy_code=" << counters.integrated_feedback_energy_code << "\n";
  stream << "deposited_feedback_energy_code=" << counters.deposited_feedback_energy_code << "\n";
  stream << "integrated_duty_cycle_active_time_code=" << counters.integrated_duty_cycle_active_time_code << "\n";
  stream << "integrated_duty_cycle_total_time_code=" << counters.integrated_duty_cycle_total_time_code << "\n";

  const std::string text = stream.str();
  core::ModuleSidecarBlock block;
  block.module_name = "black_hole_agn";
  block.schema_version = m_config.metadata_schema_version;
  block.payload.resize(text.size());
  for (std::size_t i = 0; i < text.size(); ++i) {
    block.payload[i] = static_cast<std::byte>(text[i]);
  }
  return block;
}

BlackHoleAgnConfig makeBlackHoleAgnConfig(const core::PhysicsConfig& physics_config) {
  BlackHoleAgnConfig config;
  config.enabled = physics_config.enable_black_hole_agn;
  config.seed_halo_mass_threshold_code = physics_config.bh_seed_halo_mass_threshold_code;
  config.seed_mass_code = physics_config.bh_seed_mass_code;
  config.seed_max_per_cell = physics_config.bh_seed_max_per_cell;
  config.alpha_bondi = physics_config.bh_alpha_bondi;
  config.use_eddington_cap = physics_config.bh_use_eddington_cap;
  config.epsilon_r = physics_config.bh_epsilon_r;
  config.epsilon_f = physics_config.bh_epsilon_f;
  config.feedback_coupling_efficiency = physics_config.bh_feedback_coupling_efficiency;
  config.duty_cycle_active_edd_ratio_threshold = physics_config.bh_duty_cycle_active_edd_ratio_threshold;
  config.proton_mass_si = physics_config.bh_proton_mass_si;
  config.thomson_cross_section_si = physics_config.bh_thomson_cross_section_si;
  config.newton_g_si = physics_config.bh_newton_g_si;
  config.speed_of_light_si = physics_config.bh_speed_of_light_si;
  return config;
}

BlackHoleAgnCallback::BlackHoleAgnCallback(BlackHoleAgnModel model) : m_model(std::move(model)) {}

std::string_view BlackHoleAgnCallback::callbackName() const { return "black_hole_agn_callback"; }

void BlackHoleAgnCallback::onStage(core::StepContext& context) {
  if (context.stage != core::IntegrationStage::kSourceTerms) {
    return;
  }
  m_last_step_report =
      m_model.apply(context.state, m_seed_candidates, context.integrator_state.dt_time_code, context.integrator_state.step_index);
}

void BlackHoleAgnCallback::setSeedCandidates(std::span<const BlackHoleSeedCandidate> seed_candidates) {
  m_seed_candidates.assign(seed_candidates.begin(), seed_candidates.end());
}

const BlackHoleAgnStepReport& BlackHoleAgnCallback::lastStepReport() const noexcept {
  return m_last_step_report;
}

}  // namespace cosmosim::physics
