#include "cosmosim/physics/black_hole_agn.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

#include "cosmosim/core/constants.hpp"
#include "cosmosim/physics/star_formation.hpp"

namespace cosmosim::physics {
namespace {

constexpr double k_mass_floor = 1.0e-20;
constexpr double k_density_floor = 1.0e-20;
constexpr double k_speed_floor = 1.0e-20;
constexpr double k_time_floor = 1.0e-30;

[[nodiscard]] std::uint64_t blackHoleSeedBirthKey(
    std::uint64_t gas_cell_id,
    std::uint64_t step_index) noexcept {
  // SplitMix64-style mixing gives a stable decomposition-independent key from
  // authoritative gas identity and the global integration step.
  std::uint64_t value = gas_cell_id ^ (step_index + 0x9e3779b97f4a7c15ULL);
  value = (value ^ (value >> 30U)) * 0xbf58476d1ce4e5b9ULL;
  value = (value ^ (value >> 27U)) * 0x94d049bb133111ebULL;
  value ^= value >> 31U;
  return value == 0U ? 0x6a09e667f3bcc909ULL : value;
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
      m_config.newton_g_si <= 0.0 || m_config.speed_of_light_si <= 0.0 ||
      m_config.proton_mass_code <= 0.0 || m_config.thomson_cross_section_code <= 0.0 ||
      m_config.newton_g_code <= 0.0 || m_config.speed_of_light_code <= 0.0) {
    throw std::invalid_argument("BlackHoleAgnModel: SI and code-unit constants must be > 0");
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
  const double gas_mass = state.cells.mass_code[candidate.cell_index];
  if (!std::isfinite(gas_mass) || gas_mass <= m_config.seed_mass_code + k_mass_floor) {
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
      m_config.alpha_bondi * 4.0 * core::constants::k_pi * m_config.newton_g_code * m_config.newton_g_code *
      mass * mass * rho / denom;

  rates.mdot_edd_code =
      4.0 * core::constants::k_pi * m_config.newton_g_code * mass * m_config.proton_mass_code /
      (m_config.epsilon_r * m_config.thomson_cross_section_code * m_config.speed_of_light_code);

  rates.mdot_acc_code = m_config.use_eddington_cap ? std::min(rates.mdot_bondi_code, rates.mdot_edd_code)
                                                    : rates.mdot_bondi_code;
  rates.eddington_ratio = rates.mdot_acc_code / std::max(rates.mdot_edd_code, k_time_floor);
  rates.feedback_power_code = m_config.epsilon_f * m_config.epsilon_r * rates.mdot_acc_code *
                              m_config.speed_of_light_code * m_config.speed_of_light_code;
  return rates;
}

BlackHoleAgnCounters BlackHoleAgnModel::applyAccretionFromView(
    BlackHoleAgnAccretionView view,
    double dt_code,
    double scale_factor,
    bool density_is_comoving) const {
  BlackHoleAgnCounters counters;
  if (!m_config.enabled || dt_code <= 0.0) {
    return counters;
  }
  if (!std::isfinite(dt_code) || !std::isfinite(scale_factor) || scale_factor <= 0.0) {
    throw std::invalid_argument(
        "BlackHoleAgnModel: accretion requires finite positive dt and scale factor");
  }

  const auto require_same_size = [](std::size_t expected, std::size_t actual, const char* lane) {
    if (actual != expected) {
      throw std::invalid_argument(std::string("BlackHoleAgnModel: inconsistent accretion view lane ") + lane);
    }
  };
  const std::size_t bh_count = view.particle_index.size();
  require_same_size(bh_count, view.host_cell_index.size(), "host_cell_index");
  require_same_size(bh_count, view.subgrid_mass_code.size(), "subgrid_mass_code");
  require_same_size(bh_count, view.accretion_rate_code.size(), "accretion_rate_code");
  require_same_size(bh_count, view.feedback_energy_code.size(), "feedback_energy_code");
  require_same_size(bh_count, view.eddington_ratio.size(), "eddington_ratio");
  require_same_size(bh_count, view.cumulative_accreted_mass_code.size(), "cumulative_accreted_mass_code");
  require_same_size(bh_count, view.cumulative_feedback_energy_code.size(), "cumulative_feedback_energy_code");
  require_same_size(bh_count, view.duty_cycle_active_time_code.size(), "duty_cycle_active_time_code");
  require_same_size(bh_count, view.duty_cycle_total_time_code.size(), "duty_cycle_total_time_code");

  const std::size_t gas_count = view.gas_mass_code.size();
  require_same_size(gas_count, view.gas_density_code.size(), "gas_density_code");
  require_same_size(gas_count, view.gas_metal_mass_code.size(), "gas_metal_mass_code");
  require_same_size(gas_count, view.gas_sound_speed_code.size(), "gas_sound_speed_code");
  require_same_size(gas_count, view.gas_velocity_x_peculiar.size(), "gas_velocity_x_peculiar");
  require_same_size(gas_count, view.gas_velocity_y_peculiar.size(), "gas_velocity_y_peculiar");
  require_same_size(gas_count, view.gas_velocity_z_peculiar.size(), "gas_velocity_z_peculiar");
  require_same_size(gas_count, view.gas_internal_energy_code.size(), "gas_internal_energy_code");

  const std::size_t particle_count = view.particle_mass_code.size();
  require_same_size(particle_count, view.particle_velocity_x_peculiar.size(), "particle_velocity_x_peculiar");
  require_same_size(particle_count, view.particle_velocity_y_peculiar.size(), "particle_velocity_y_peculiar");
  require_same_size(particle_count, view.particle_velocity_z_peculiar.size(), "particle_velocity_z_peculiar");

  const double density_scale = density_is_comoving
      ? 1.0 / (scale_factor * scale_factor * scale_factor)
      : 1.0;

  for (const std::uint32_t bh_local_index : view.active_black_hole_indices) {
    ++counters.scanned_bh;
    if (bh_local_index >= bh_count) {
      throw std::out_of_range("BlackHoleAgnModel: active BH index is out of range");
    }
    const std::uint32_t particle_index = view.particle_index[bh_local_index];
    const std::uint32_t host_cell_index = view.host_cell_index[bh_local_index];
    if (particle_index >= particle_count || host_cell_index >= gas_count) {
      throw std::out_of_range("BlackHoleAgnModel: BH particle/host-cell index is out of range");
    }

    const double gas_mass_before = view.gas_mass_code[host_cell_index];
    const double gas_density_stored = view.gas_density_code[host_cell_index];
    const double gas_metal_before = view.gas_metal_mass_code[host_cell_index];
    const double sound_speed = view.gas_sound_speed_code[host_cell_index];
    const double bh_mass_before = view.subgrid_mass_code[bh_local_index];
    if (!std::isfinite(gas_mass_before) || gas_mass_before < 0.0 ||
        !std::isfinite(gas_density_stored) || gas_density_stored < 0.0 ||
        !std::isfinite(gas_metal_before) || gas_metal_before < 0.0 ||
        !std::isfinite(sound_speed) || sound_speed < 0.0 ||
        !std::isfinite(bh_mass_before) || bh_mass_before < 0.0) {
      throw std::invalid_argument("BlackHoleAgnModel: non-finite or negative physical accretion state");
    }

    ++counters.active_bh;
    const double dv_x = view.particle_velocity_x_peculiar[particle_index] -
        view.gas_velocity_x_peculiar[host_cell_index];
    const double dv_y = view.particle_velocity_y_peculiar[particle_index] -
        view.gas_velocity_y_peculiar[host_cell_index];
    const double dv_z = view.particle_velocity_z_peculiar[particle_index] -
        view.gas_velocity_z_peculiar[host_cell_index];
    const double relative_velocity = std::sqrt(dv_x * dv_x + dv_y * dv_y + dv_z * dv_z);
    const double gas_density_physical = gas_density_stored * density_scale;
    const BlackHoleRates requested_rates = computeAccretionRates(
        bh_mass_before, gas_density_physical, sound_speed, relative_velocity);

    const double requested_mass = std::max(requested_rates.mdot_acc_code * dt_code, 0.0);
    const double available_mass = std::max(gas_mass_before - k_mass_floor, 0.0);
    const double delta_mass = std::min(requested_mass, available_mass);
    const double actual_mdot = delta_mass / dt_code;
    const double actual_eddington_ratio =
        actual_mdot / std::max(requested_rates.mdot_edd_code, k_time_floor);
    const double actual_feedback_power = m_config.epsilon_f * m_config.epsilon_r * actual_mdot *
        m_config.speed_of_light_code * m_config.speed_of_light_code;
    const double delta_feedback_energy = std::max(actual_feedback_power * dt_code, 0.0);
    const double requested_deposit = delta_feedback_energy * m_config.feedback_coupling_efficiency;

    const double gas_mass_after = gas_mass_before - delta_mass;
    const double retained_fraction = gas_mass_before > 0.0 ? gas_mass_after / gas_mass_before : 1.0;
    const double transferred_metal_mass = gas_metal_before * (1.0 - retained_fraction);
    view.gas_mass_code[host_cell_index] = gas_mass_after;
    view.gas_density_code[host_cell_index] = gas_density_stored * retained_fraction;
    view.gas_metal_mass_code[host_cell_index] = std::max(gas_metal_before - transferred_metal_mass, 0.0);

    // The tracked baseline uses a Newtonian mass ledger: mdot_acc is the gas
    // rest mass transferred into the BH dynamical/subgrid mass. epsilon_r is
    // an unresolved radiative-energy yield, not an implicit deletion of mass.
    const double bh_mass_after = bh_mass_before + delta_mass;
    if (delta_mass > 0.0 && bh_mass_after > 0.0) {
      view.particle_velocity_x_peculiar[particle_index] =
          (bh_mass_before * view.particle_velocity_x_peculiar[particle_index] +
           delta_mass * view.gas_velocity_x_peculiar[host_cell_index]) / bh_mass_after;
      view.particle_velocity_y_peculiar[particle_index] =
          (bh_mass_before * view.particle_velocity_y_peculiar[particle_index] +
           delta_mass * view.gas_velocity_y_peculiar[host_cell_index]) / bh_mass_after;
      view.particle_velocity_z_peculiar[particle_index] =
          (bh_mass_before * view.particle_velocity_z_peculiar[particle_index] +
           delta_mass * view.gas_velocity_z_peculiar[host_cell_index]) / bh_mass_after;
    }
    view.subgrid_mass_code[bh_local_index] = bh_mass_after;
    view.particle_mass_code[particle_index] = bh_mass_after;

    // internal_energy_code is specific internal energy. Deposit total feedback
    // energy by dividing by the remaining host-gas mass instead of mixing units.
    double deposited_feedback_energy = 0.0;
    if (requested_deposit > 0.0 && gas_mass_after > k_mass_floor) {
      view.gas_internal_energy_code[host_cell_index] += requested_deposit / gas_mass_after;
      deposited_feedback_energy = requested_deposit;
    }

    view.accretion_rate_code[bh_local_index] = actual_mdot;
    view.feedback_energy_code[bh_local_index] = delta_feedback_energy;
    view.eddington_ratio[bh_local_index] = actual_eddington_ratio;
    view.cumulative_accreted_mass_code[bh_local_index] += delta_mass;
    view.cumulative_feedback_energy_code[bh_local_index] += delta_feedback_energy;
    view.duty_cycle_total_time_code[bh_local_index] += dt_code;

    if (actual_eddington_ratio >= m_config.duty_cycle_active_edd_ratio_threshold) {
      view.duty_cycle_active_time_code[bh_local_index] += dt_code;
      counters.integrated_duty_cycle_active_time_code += dt_code;
    }
    counters.integrated_duty_cycle_total_time_code += dt_code;
    counters.integrated_accreted_mass_code += delta_mass;
    counters.gas_mass_removed_code += delta_mass;
    counters.metal_mass_transferred_code += transferred_metal_mass;
    counters.integrated_feedback_energy_code += delta_feedback_energy;
    counters.deposited_feedback_energy_code += deposited_feedback_energy;
  }
  return counters;
}

BlackHoleAgnStepReport BlackHoleAgnModel::apply(
    core::SimulationState& state,
    std::span<const BlackHoleSeedCandidate> seed_candidates,
    double dt_code,
    double scale_factor,
    bool density_is_comoving,
    std::uint64_t step_index,
    ParticleIdPrecommit* id_precommit) const {
  BlackHoleAgnStepReport report;
  if (!m_config.enabled || dt_code <= 0.0) {
    state.sidecars.upsert(buildMetadataSidecar(report.counters));
    return report;
  }

  std::vector<std::uint32_t> active_black_hole_indices(state.black_holes.size());
  for (std::size_t i = 0; i < active_black_hole_indices.size(); ++i) {
    if (i > static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max())) {
      throw std::overflow_error("BlackHoleAgnModel: local BH count exceeds uint32 index capacity");
    }
    active_black_hole_indices[i] = static_cast<std::uint32_t>(i);
  }
  BlackHoleAgnAccretionView accretion_view{
      .active_black_hole_indices = active_black_hole_indices,
      .particle_index = state.black_holes.particle_index,
      .host_cell_index = state.black_holes.host_cell_index,
      .subgrid_mass_code = state.black_holes.subgrid_mass_code,
      .accretion_rate_code = state.black_holes.accretion_rate_code,
      .feedback_energy_code = state.black_holes.feedback_energy_code,
      .eddington_ratio = state.black_holes.eddington_ratio,
      .cumulative_accreted_mass_code = state.black_holes.cumulative_accreted_mass_code,
      .cumulative_feedback_energy_code = state.black_holes.cumulative_feedback_energy_code,
      .duty_cycle_active_time_code = state.black_holes.duty_cycle_active_time_code,
      .duty_cycle_total_time_code = state.black_holes.duty_cycle_total_time_code,
      .gas_mass_code = state.cells.mass_code,
      .gas_density_code = state.gas_cells.density_code,
      .gas_metal_mass_code = state.gas_cells.metal_mass_code,
      .gas_sound_speed_code = state.gas_cells.sound_speed_code,
      .gas_velocity_x_peculiar = state.gas_cells.velocity_x_peculiar,
      .gas_velocity_y_peculiar = state.gas_cells.velocity_y_peculiar,
      .gas_velocity_z_peculiar = state.gas_cells.velocity_z_peculiar,
      .gas_internal_energy_code = state.gas_cells.internal_energy_code,
      .particle_mass_code = state.particles.mass_code,
      .particle_velocity_x_peculiar = state.particles.velocity_x_peculiar,
      .particle_velocity_y_peculiar = state.particles.velocity_y_peculiar,
      .particle_velocity_z_peculiar = state.particles.velocity_z_peculiar,
  };
  const BlackHoleAgnCounters accretion_counters =
      applyAccretionFromView(accretion_view, dt_code, scale_factor, density_is_comoving);
  report.counters.scanned_bh += accretion_counters.scanned_bh;
  report.counters.active_bh += accretion_counters.active_bh;
  report.counters.integrated_accreted_mass_code += accretion_counters.integrated_accreted_mass_code;
  report.counters.gas_mass_removed_code += accretion_counters.gas_mass_removed_code;
  report.counters.metal_mass_transferred_code += accretion_counters.metal_mass_transferred_code;
  report.counters.integrated_feedback_energy_code += accretion_counters.integrated_feedback_energy_code;
  report.counters.deposited_feedback_energy_code += accretion_counters.deposited_feedback_energy_code;
  report.counters.integrated_duty_cycle_active_time_code += accretion_counters.integrated_duty_cycle_active_time_code;
  report.counters.integrated_duty_cycle_total_time_code += accretion_counters.integrated_duty_cycle_total_time_code;

  std::vector<const BlackHoleSeedCandidate*> accepted_candidates;
  std::vector<std::uint64_t> birth_keys;
  accepted_candidates.reserve(seed_candidates.size());
  birth_keys.reserve(seed_candidates.size());
  for (const BlackHoleSeedCandidate& candidate : seed_candidates) {
    ++report.counters.seed_candidates;
    if (!isSeedEligible(state, candidate)) {
      continue;
    }
    const std::uint64_t gas_cell_id = state.gas_cells.gas_cell_id[candidate.cell_index];
    if (gas_cell_id == 0U) {
      throw std::runtime_error("BlackHoleAgnModel: BH seeding requires stable nonzero gas_cell_id");
    }
    accepted_candidates.push_back(&candidate);
    birth_keys.push_back(blackHoleSeedBirthKey(gas_cell_id, step_index));
  }

  LocalParticleIdRegistry local_registry;
  ParticleIdPrecommit& registry = id_precommit != nullptr
      ? *id_precommit
      : static_cast<ParticleIdPrecommit&>(local_registry);
  const std::vector<std::uint64_t> seeded_ids = registry.precommit(state, birth_keys);
  if (seeded_ids.size() != accepted_candidates.size()) {
    throw std::runtime_error("BlackHoleAgnModel: particle-ID precommit returned the wrong seed count");
  }

  for (std::size_t seed_index = 0; seed_index < accepted_candidates.size(); ++seed_index) {
    const BlackHoleSeedCandidate& candidate = *accepted_candidates[seed_index];
    const std::size_t cell_index = candidate.cell_index;
    const double gas_mass_before = state.cells.mass_code[cell_index];
    if (gas_mass_before <= m_config.seed_mass_code + k_mass_floor) {
      throw std::runtime_error("BlackHoleAgnModel: host gas mass changed below seed transaction requirement");
    }
    const double retained_fraction = (gas_mass_before - m_config.seed_mass_code) / gas_mass_before;
    const double metal_before = state.gas_cells.metal_mass_code[cell_index];
    const double swallowed_metal = metal_before * (1.0 - retained_fraction);
    state.cells.mass_code[cell_index] -= m_config.seed_mass_code;
    state.gas_cells.density_code[cell_index] *= retained_fraction;
    state.gas_cells.metal_mass_code[cell_index] = std::max(metal_before - swallowed_metal, 0.0);

    const std::size_t particle_index = state.particles.size();
    if (particle_index > static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max())) {
      throw std::overflow_error("BlackHoleAgnModel: local particle count exceeds BH sidecar index capacity");
    }
    state.resizeParticles(particle_index + 1);
    state.particles.position_x_comoving[particle_index] = state.cells.center_x_comoving[cell_index];
    state.particles.position_y_comoving[particle_index] = state.cells.center_y_comoving[cell_index];
    state.particles.position_z_comoving[particle_index] = state.cells.center_z_comoving[cell_index];
    state.particles.velocity_x_peculiar[particle_index] = state.gas_cells.velocity_x_peculiar[cell_index];
    state.particles.velocity_y_peculiar[particle_index] = state.gas_cells.velocity_y_peculiar[cell_index];
    state.particles.velocity_z_peculiar[particle_index] = state.gas_cells.velocity_z_peculiar[cell_index];
    state.particles.mass_code[particle_index] = m_config.seed_mass_code;
    state.particle_sidecar.particle_id[particle_index] = seeded_ids[seed_index];
    state.particle_sidecar.sfc_key[particle_index] = seeded_ids[seed_index];
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
    report.counters.seeded_mass_code += m_config.seed_mass_code;
    report.counters.gas_mass_removed_code += m_config.seed_mass_code;
    report.counters.metal_mass_transferred_code += swallowed_metal;
    report.seeded_cell_indices.push_back(candidate.cell_index);
    report.seeded_particle_ids.push_back(seeded_ids[seed_index]);
  }

  state.species.count_by_species[static_cast<std::size_t>(core::ParticleSpecies::kBlackHole)] +=
      report.counters.seeded_bh;
  state.rebuildSpeciesIndex();
  state.sidecars.upsert(buildMetadataSidecar(report.counters));
  return report;
}

BlackHoleAgnStepReport BlackHoleAgnModel::apply(
    core::SimulationState& state,
    std::span<const BlackHoleSeedCandidate> seed_candidates,
    double dt_code,
    std::uint64_t step_index) const {
  return apply(state, seed_candidates, dt_code, 1.0, false, step_index, nullptr);
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
  stream << "gas_mass_removed_code=" << counters.gas_mass_removed_code << "\n";
  stream << "metal_mass_transferred_code=" << counters.metal_mass_transferred_code << "\n";
  stream << "seeded_mass_code=" << counters.seeded_mass_code << "\n";
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

std::span<const core::IntegrationStage> BlackHoleAgnCallback::integrationStages() const {
  static constexpr std::array stages{core::IntegrationStage::kSourceTerms};
  return stages;
}

std::span<const core::StageContract> BlackHoleAgnCallback::stageContracts() const {
  static constexpr std::array contracts{core::StageContract{
      .stage = core::IntegrationStage::kSourceTerms,
      .required_inputs = core::StageDataDomain::kParticles | core::StageDataDomain::kGasCells,
      .mutated_state = core::StageDataDomain::kParticles | core::StageDataDomain::kGasCells | core::StageDataDomain::kDiagnostics,
      .produced_outputs = core::StageDataDomain::kParticles | core::StageDataDomain::kGasCells | core::StageDataDomain::kDiagnostics,
      .allowed_side_effects = core::StageDataDomain::kDiagnostics,
      .sync_requirements = core::StageSyncRequirement::kLocalOnly,
      .active_set_family = core::StageActiveSetFamily::kActiveParticles,
      .restart_safety = core::StageSafety::kUnsafe,
      .output_safety = core::StageSafety::kUnsafe,
      .owner = core::StageSubsystem::kSources,
  }};
  return contracts;
}

void BlackHoleAgnCallback::onStage(core::StepContext& context) {
  if (context.stage != core::IntegrationStage::kSourceTerms) {
    throw std::logic_error("black-hole AGN handler received an unregistered stage");
  }
  const double evaluation_scale_factor = context.timeline_step.scale_factor_end;
  m_last_step_report = m_model.apply(
      context.state,
      m_seed_candidates,
      context.integrator_state.dt_time_code,
      evaluation_scale_factor,
      context.cosmology_background != nullptr,
      context.integrator_state.step_index,
      nullptr);
}

void BlackHoleAgnCallback::setSeedCandidates(std::span<const BlackHoleSeedCandidate> seed_candidates) {
  m_seed_candidates.assign(seed_candidates.begin(), seed_candidates.end());
}

const BlackHoleAgnStepReport& BlackHoleAgnCallback::lastStepReport() const noexcept {
  return m_last_step_report;
}

}  // namespace cosmosim::physics
