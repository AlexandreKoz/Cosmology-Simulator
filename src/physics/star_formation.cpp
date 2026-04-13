#include "cosmosim/physics/star_formation.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

#include "cosmosim/core/constants.hpp"

namespace cosmosim::physics {
namespace {

constexpr double k_density_floor = 1.0e-20;
constexpr double k_mass_floor = 1.0e-20;
constexpr double k_u01_norm = 1.0 / static_cast<double>(std::numeric_limits<std::uint64_t>::max());

[[nodiscard]] std::uint64_t splitmix64(std::uint64_t x) {
  x += 0x9e3779b97f4a7c15ull;
  x = (x ^ (x >> 30u)) * 0xbf58476d1ce4e5b9ull;
  x = (x ^ (x >> 27u)) * 0x94d049bb133111ebull;
  return x ^ (x >> 31u);
}

[[nodiscard]] double uniform01(
    std::uint64_t seed,
    std::uint64_t step_index,
    std::uint32_t cell_index,
    std::uint32_t rank_local_seed_offset) {
  std::uint64_t mixed = seed;
  mixed ^= splitmix64(step_index + 0x9e3779b97f4a7c15ull);
  mixed ^= splitmix64(static_cast<std::uint64_t>(cell_index) + 0xbf58476d1ce4e5b9ull);
  mixed ^= splitmix64(static_cast<std::uint64_t>(rank_local_seed_offset) + 0x94d049bb133111ebull);
  const std::uint64_t bits = splitmix64(mixed);
  return static_cast<double>(bits) * k_u01_norm;
}

[[nodiscard]] std::uint64_t nextParticleId(const core::SimulationState& state) {
  std::uint64_t max_id = 0;
  for (const std::uint64_t particle_id : state.particle_sidecar.particle_id) {
    max_id = std::max(max_id, particle_id);
  }
  return max_id + 1;
}

}  // namespace

StarFormationModel::StarFormationModel(StarFormationConfig config) : m_config(std::move(config)) {
  if (m_config.density_threshold_code <= 0.0) {
    throw std::invalid_argument("StarFormationModel: density_threshold_code must be > 0");
  }
  if (m_config.temperature_threshold_k <= 0.0) {
    throw std::invalid_argument("StarFormationModel: temperature_threshold_k must be > 0");
  }
  if (m_config.epsilon_ff < 0.0 || m_config.epsilon_ff > 1.0) {
    throw std::invalid_argument("StarFormationModel: epsilon_ff must be in [0, 1]");
  }
  if (m_config.min_star_particle_mass_code <= 0.0) {
    throw std::invalid_argument("StarFormationModel: min_star_particle_mass_code must be > 0");
  }
}

const StarFormationConfig& StarFormationModel::config() const noexcept {
  return m_config;
}

bool StarFormationModel::isEligible(const StarFormationCellInput& cell) const {
  if (!m_config.enabled) {
    return false;
  }
  if (cell.gas_mass_code <= k_mass_floor || cell.gas_density_code < m_config.density_threshold_code) {
    return false;
  }
  if (cell.gas_temperature_k > m_config.temperature_threshold_k) {
    return false;
  }
  if (cell.velocity_divergence_code > m_config.min_converging_flow_rate_code) {
    return false;
  }
  return true;
}

double StarFormationModel::freeFallTimeCode(double gas_density_code) const {
  const double safe_density = std::max(gas_density_code, k_density_floor);
  return std::sqrt(3.0 * core::constants::k_pi / (32.0 * core::constants::k_newton_g_si * safe_density));
}

double StarFormationModel::sfrDensityRateCode(double gas_density_code) const {
  if (gas_density_code <= 0.0) {
    return 0.0;
  }
  const double t_ff = freeFallTimeCode(gas_density_code);
  return m_config.epsilon_ff * gas_density_code / std::max(t_ff, 1.0e-30);
}

double StarFormationModel::expectedSpawnMassCode(const StarFormationCellInput& cell, double dt_code) const {
  if (dt_code <= 0.0 || !isEligible(cell)) {
    return 0.0;
  }
  const double t_ff = freeFallTimeCode(cell.gas_density_code);
  const double expected = m_config.epsilon_ff * cell.gas_mass_code * dt_code / std::max(t_ff, 1.0e-30);
  return std::clamp(expected, 0.0, cell.gas_mass_code);
}

StarFormationCellOutcome StarFormationModel::sampleCellOutcome(
    const StarFormationCellInput& cell,
    double dt_code,
    std::uint64_t step_index,
    std::uint32_t rank_local_seed_offset) const {
  StarFormationCellOutcome outcome;
  outcome.eligible = isEligible(cell);
  if (!outcome.eligible || dt_code <= 0.0) {
    return outcome;
  }

  outcome.free_fall_time_code = freeFallTimeCode(cell.gas_density_code);
  outcome.sfr_density_rate_code = sfrDensityRateCode(cell.gas_density_code);
  outcome.expected_spawn_mass_code = expectedSpawnMassCode(cell, dt_code);
  if (outcome.expected_spawn_mass_code <= k_mass_floor) {
    return outcome;
  }

  if (!m_config.stochastic_spawning) {
    outcome.spawned_mass_code = outcome.expected_spawn_mass_code;
    outcome.spawned_particle_count = outcome.spawned_mass_code > 0.0 ? 1u : 0u;
    return outcome;
  }

  const double lambda = outcome.expected_spawn_mass_code / m_config.min_star_particle_mass_code;
  std::uint32_t spawn_count = static_cast<std::uint32_t>(std::floor(lambda));
  const double remainder = std::clamp(lambda - static_cast<double>(spawn_count), 0.0, 1.0);
  outcome.random_u01 = uniform01(m_config.random_seed, step_index, cell.cell_index, rank_local_seed_offset);
  if (outcome.random_u01 < remainder) {
    ++spawn_count;
  }
  outcome.spawned_particle_count = spawn_count;

  const double spawned_mass = static_cast<double>(spawn_count) * m_config.min_star_particle_mass_code;
  outcome.spawned_mass_code = std::min(spawned_mass, cell.gas_mass_code);
  return outcome;
}

StarFormationStepReport StarFormationModel::apply(
    core::SimulationState& state,
    std::span<const std::uint32_t> active_cell_indices,
    std::span<const double> velocity_divergence_code,
    std::span<const double> metallicity_mass_fraction,
    double dt_code,
    double scale_factor,
    std::uint64_t step_index,
    std::uint32_t rank_local_seed_offset) const {
  StarFormationStepReport report;
  if (!m_config.enabled || dt_code <= 0.0) {
    state.sidecars.upsert(buildMetadataSidecar(report.counters));
    return report;
  }

  std::uint64_t next_particle_id = nextParticleId(state);
  for (const std::uint32_t cell_index : active_cell_indices) {
    ++report.counters.scanned_cells;
    if (cell_index >= state.cells.size() || cell_index >= state.gas_cells.size()) {
      continue;
    }

    StarFormationCellInput cell;
    cell.cell_index = cell_index;
    cell.gas_mass_code = state.cells.mass_code[cell_index];
    cell.gas_density_code = state.gas_cells.density_code[cell_index];
    cell.gas_temperature_k = state.gas_cells.temperature_code[cell_index];
    cell.velocity_divergence_code =
        (cell_index < velocity_divergence_code.size()) ? velocity_divergence_code[cell_index] : 0.0;
    cell.metallicity_mass_fraction =
        (cell_index < metallicity_mass_fraction.size()) ? metallicity_mass_fraction[cell_index] : 0.0;

    const StarFormationCellOutcome outcome = sampleCellOutcome(cell, dt_code, step_index, rank_local_seed_offset);
    if (!outcome.eligible) {
      continue;
    }
    ++report.counters.eligible_cells;
    report.counters.expected_spawn_mass_code += outcome.expected_spawn_mass_code;
    if (outcome.spawned_mass_code <= k_mass_floor) {
      continue;
    }

    const double gas_mass_before = state.cells.mass_code[cell_index];
    const double gas_density_before = state.gas_cells.density_code[cell_index];
    const double transfer_mass = std::min(outcome.spawned_mass_code, gas_mass_before);
    const double transfer_fraction = transfer_mass / std::max(gas_mass_before, k_mass_floor);

    state.cells.mass_code[cell_index] = std::max(gas_mass_before - transfer_mass, 0.0);
    state.gas_cells.density_code[cell_index] = std::max(gas_density_before * (1.0 - transfer_fraction), 0.0);

    const std::size_t particle_index = state.particles.size();
    state.resizeParticles(particle_index + 1);
    state.particles.position_x_comoving[particle_index] = state.cells.center_x_comoving[cell_index];
    state.particles.position_y_comoving[particle_index] = state.cells.center_y_comoving[cell_index];
    state.particles.position_z_comoving[particle_index] = state.cells.center_z_comoving[cell_index];
    state.particles.velocity_x_peculiar[particle_index] = 0.0;
    state.particles.velocity_y_peculiar[particle_index] = 0.0;
    state.particles.velocity_z_peculiar[particle_index] = 0.0;
    state.particles.mass_code[particle_index] = transfer_mass;
    state.particles.time_bin[particle_index] = state.cells.time_bin[cell_index];

    state.particle_sidecar.particle_id[particle_index] = next_particle_id++;
    state.particle_sidecar.sfc_key[particle_index] = state.particle_sidecar.particle_id[particle_index];
    state.particle_sidecar.species_tag[particle_index] =
        static_cast<std::uint32_t>(core::ParticleSpecies::kStar);
    state.particle_sidecar.particle_flags[particle_index] = 0u;
    state.particle_sidecar.owning_rank[particle_index] = rank_local_seed_offset;

    const std::size_t star_index = state.star_particles.size();
    state.star_particles.resize(star_index + 1);
    state.star_particles.particle_index[star_index] = static_cast<std::uint32_t>(particle_index);
    state.star_particles.formation_scale_factor[star_index] = scale_factor;
    state.star_particles.birth_mass_code[star_index] = transfer_mass;
    state.star_particles.metallicity_mass_fraction[star_index] = cell.metallicity_mass_fraction;
    state.star_particles.stellar_age_years_last[star_index] = 0.0;
    state.star_particles.stellar_returned_mass_cumulative_code[star_index] = 0.0;
    state.star_particles.stellar_returned_metals_cumulative_code[star_index] = 0.0;
    state.star_particles.stellar_feedback_energy_cumulative_erg[star_index] = 0.0;
    for (std::size_t channel = 0; channel < state.star_particles.stellar_returned_mass_channel_cumulative_code.size();
         ++channel) {
      state.star_particles.stellar_returned_mass_channel_cumulative_code[channel][star_index] = 0.0;
      state.star_particles.stellar_returned_metals_channel_cumulative_code[channel][star_index] = 0.0;
      state.star_particles.stellar_feedback_energy_channel_cumulative_erg[channel][star_index] = 0.0;
    }

    ++report.counters.spawn_events;
    report.counters.spawned_particles += 1;
    report.counters.spawned_mass_code += transfer_mass;
    report.spawned_from_cells.push_back(cell_index);
  }

  state.species.count_by_species[static_cast<std::size_t>(core::ParticleSpecies::kStar)] +=
      report.counters.spawned_particles;
  state.rebuildSpeciesIndex();
  state.sidecars.upsert(buildMetadataSidecar(report.counters));
  return report;
}

core::ModuleSidecarBlock StarFormationModel::buildMetadataSidecar(const StarFormationCounters& counters) const {
  std::ostringstream stream;
  stream << "module=star_formation\n";
  stream << "schema_version=" << m_config.metadata_schema_version << "\n";
  stream << "density_threshold_code=" << m_config.density_threshold_code << "\n";
  stream << "temperature_threshold_k=" << m_config.temperature_threshold_k << "\n";
  stream << "epsilon_ff=" << m_config.epsilon_ff << "\n";
  stream << "min_star_particle_mass_code=" << m_config.min_star_particle_mass_code << "\n";
  stream << "stochastic_spawning=" << (m_config.stochastic_spawning ? "true" : "false") << "\n";
  stream << "random_seed=" << m_config.random_seed << "\n";
  stream << "scanned_cells=" << counters.scanned_cells << "\n";
  stream << "eligible_cells=" << counters.eligible_cells << "\n";
  stream << "spawn_events=" << counters.spawn_events << "\n";
  stream << "spawned_particles=" << counters.spawned_particles << "\n";
  stream << "expected_spawn_mass_code=" << counters.expected_spawn_mass_code << "\n";
  stream << "spawned_mass_code=" << counters.spawned_mass_code << "\n";

  const std::string text = stream.str();
  core::ModuleSidecarBlock block;
  block.module_name = "star_formation";
  block.schema_version = m_config.metadata_schema_version;
  block.payload.resize(text.size());
  for (std::size_t i = 0; i < text.size(); ++i) {
    block.payload[i] = static_cast<std::byte>(text[i]);
  }
  return block;
}

StarFormationConfig makeStarFormationConfig(const core::PhysicsConfig& physics_config) {
  StarFormationConfig config;
  config.enabled = physics_config.enable_star_formation;
  config.density_threshold_code = physics_config.sf_density_threshold_code;
  config.temperature_threshold_k = physics_config.sf_temperature_threshold_k;
  config.min_converging_flow_rate_code = physics_config.sf_min_converging_flow_rate_code;
  config.epsilon_ff = physics_config.sf_epsilon_ff;
  config.min_star_particle_mass_code = physics_config.sf_min_star_particle_mass_code;
  config.stochastic_spawning = physics_config.sf_stochastic_spawning;
  config.random_seed = physics_config.sf_random_seed;
  return config;
}

StarFormationCallback::StarFormationCallback(StarFormationModel model, std::uint32_t rank_local_seed_offset)
    : m_model(std::move(model)), m_rank_local_seed_offset(rank_local_seed_offset) {}

std::string_view StarFormationCallback::callbackName() const { return "star_formation_callback"; }

void StarFormationCallback::onStage(core::StepContext& context) {
  if (context.stage != core::IntegrationStage::kSourceTerms) {
    return;
  }

  const std::size_t cell_count = context.state.cells.size();
  if (cell_count == 0) {
    m_last_step_report = {};
    return;
  }
  ensureFieldSizes(cell_count);

  std::span<const std::uint32_t> active_cells = context.active_set.cell_indices;
  if (!context.active_set.cells_are_subset && active_cells.empty()) {
    m_full_cell_indices.resize(cell_count);
    std::iota(m_full_cell_indices.begin(), m_full_cell_indices.end(), 0U);
    active_cells = m_full_cell_indices;
  }

  m_last_step_report = m_model.apply(
      context.state,
      active_cells,
      m_velocity_divergence_code,
      m_metallicity_mass_fraction,
      context.integrator_state.dt_time_code,
      context.integrator_state.current_scale_factor,
      context.integrator_state.step_index,
      m_rank_local_seed_offset);
}

void StarFormationCallback::setVelocityDivergenceCode(std::span<const double> velocity_divergence_code) {
  m_velocity_divergence_code.assign(velocity_divergence_code.begin(), velocity_divergence_code.end());
}

void StarFormationCallback::setMetallicityMassFraction(std::span<const double> metallicity_mass_fraction) {
  m_metallicity_mass_fraction.assign(metallicity_mass_fraction.begin(), metallicity_mass_fraction.end());
}

void StarFormationCallback::setRankLocalSeedOffset(std::uint32_t rank_local_seed_offset) {
  m_rank_local_seed_offset = rank_local_seed_offset;
}

const StarFormationStepReport& StarFormationCallback::lastStepReport() const noexcept {
  return m_last_step_report;
}

void StarFormationCallback::ensureFieldSizes(std::size_t cell_count) {
  if (m_velocity_divergence_code.size() < cell_count) {
    m_velocity_divergence_code.resize(cell_count, 0.0);
  }
  if (m_metallicity_mass_fraction.size() < cell_count) {
    m_metallicity_mass_fraction.resize(cell_count, 0.0);
  }
}

}  // namespace cosmosim::physics
