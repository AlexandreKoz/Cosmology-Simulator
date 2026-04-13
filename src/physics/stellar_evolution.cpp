#include "cosmosim/physics/stellar_evolution.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace cosmosim::physics {
namespace {

constexpr double k_age_floor_yr = 0.0;
constexpr double k_mass_floor = 1.0e-20;

[[nodiscard]] double linearInterpolate(double x0, double x1, double y0, double y1, double x) {
  if (x1 <= x0) {
    return y0;
  }
  const double t = std::clamp((x - x0) / (x1 - x0), 0.0, 1.0);
  return y0 + t * (y1 - y0);
}

[[nodiscard]] std::array<double, k_stellar_yield_channel_count> differenceArray(
    const std::array<double, k_stellar_yield_channel_count>& end,
    const std::array<double, k_stellar_yield_channel_count>& begin) {
  std::array<double, k_stellar_yield_channel_count> result{};
  for (std::size_t channel = 0; channel < k_stellar_yield_channel_count; ++channel) {
    result[channel] = std::max(end[channel] - begin[channel], 0.0);
  }
  return result;
}

[[nodiscard]] std::string trim(const std::string& input) {
  const auto begin = input.find_first_not_of(" \t\r\n");
  if (begin == std::string::npos) {
    return {};
  }
  const auto end = input.find_last_not_of(" \t\r\n");
  return input.substr(begin, end - begin + 1);
}

void normalizeChannelSums(
    double total,
    std::array<double, k_stellar_yield_channel_count>& channel_values) {
  const double channel_sum = channel_values[0] + channel_values[1] + channel_values[2];
  if (channel_sum <= 0.0 || total <= 0.0) {
    return;
  }
  const double scale = total / channel_sum;
  for (double& value : channel_values) {
    value *= scale;
  }
}

}  // namespace

bool StellarEvolutionTable::isConsistent() const noexcept {
  const std::size_t row_count = age_yr.size();
  if (row_count < 2 || return_fraction_total.size() != row_count ||
      metal_yield_fraction_total.size() != row_count || energy_erg_per_initial_mass_code.size() != row_count) {
    return false;
  }

  for (std::size_t i = 1; i < row_count; ++i) {
    if (age_yr[i] < age_yr[i - 1]) {
      return false;
    }
  }

  for (std::size_t channel = 0; channel < k_stellar_yield_channel_count; ++channel) {
    if (return_fraction_channel[channel].size() != row_count ||
        metal_yield_fraction_channel[channel].size() != row_count ||
        energy_erg_per_initial_mass_channel[channel].size() != row_count) {
      return false;
    }
  }
  return true;
}

StellarEvolutionCumulative StellarEvolutionTable::evaluateAtAgeYears(double age_yr_value) const {
  if (!isConsistent()) {
    throw std::runtime_error("StellarEvolutionTable: inconsistent table state");
  }

  const double clamped_age = std::clamp(age_yr_value, age_yr.front(), age_yr.back());
  const auto upper = std::lower_bound(age_yr.begin(), age_yr.end(), clamped_age);
  const std::size_t upper_index = static_cast<std::size_t>(std::distance(age_yr.begin(), upper));
  const std::size_t lower_index = (upper_index == 0) ? 0 : (upper_index - 1);
  const std::size_t right_index = std::min(upper_index, age_yr.size() - 1);

  StellarEvolutionCumulative result;
  result.return_fraction_total = linearInterpolate(
      age_yr[lower_index], age_yr[right_index], return_fraction_total[lower_index],
      return_fraction_total[right_index], clamped_age);
  result.metal_yield_fraction_total = linearInterpolate(
      age_yr[lower_index], age_yr[right_index], metal_yield_fraction_total[lower_index],
      metal_yield_fraction_total[right_index], clamped_age);
  result.energy_erg_per_initial_mass_code = linearInterpolate(
      age_yr[lower_index], age_yr[right_index], energy_erg_per_initial_mass_code[lower_index],
      energy_erg_per_initial_mass_code[right_index], clamped_age);

  for (std::size_t channel = 0; channel < k_stellar_yield_channel_count; ++channel) {
    result.return_fraction_channel[channel] = linearInterpolate(
        age_yr[lower_index], age_yr[right_index], return_fraction_channel[channel][lower_index],
        return_fraction_channel[channel][right_index], clamped_age);
    result.metal_yield_fraction_channel[channel] = linearInterpolate(
        age_yr[lower_index], age_yr[right_index], metal_yield_fraction_channel[channel][lower_index],
        metal_yield_fraction_channel[channel][right_index], clamped_age);
    result.energy_erg_per_initial_mass_channel[channel] = linearInterpolate(
        age_yr[lower_index], age_yr[right_index], energy_erg_per_initial_mass_channel[channel][lower_index],
        energy_erg_per_initial_mass_channel[channel][right_index], clamped_age);
  }

  return result;
}

StellarEvolutionIntervalBudget StellarEvolutionTable::integrateInterval(
    double age_begin_yr,
    double age_end_yr,
    double initial_mass_code) const {
  const double clamped_begin = std::max(age_begin_yr, k_age_floor_yr);
  const double clamped_end = std::max(age_end_yr, clamped_begin);

  const StellarEvolutionCumulative begin_state = evaluateAtAgeYears(clamped_begin);
  const StellarEvolutionCumulative end_state = evaluateAtAgeYears(clamped_end);

  StellarEvolutionIntervalBudget interval;
  interval.returned_mass_code =
      std::max(end_state.return_fraction_total - begin_state.return_fraction_total, 0.0) * initial_mass_code;
  interval.returned_metals_code =
      std::max(end_state.metal_yield_fraction_total - begin_state.metal_yield_fraction_total, 0.0) *
      initial_mass_code;
  interval.feedback_energy_erg =
      std::max(end_state.energy_erg_per_initial_mass_code - begin_state.energy_erg_per_initial_mass_code, 0.0) *
      initial_mass_code;

  auto return_fraction_channel_delta =
      differenceArray(end_state.return_fraction_channel, begin_state.return_fraction_channel);
  auto metal_fraction_channel_delta =
      differenceArray(end_state.metal_yield_fraction_channel, begin_state.metal_yield_fraction_channel);
  auto energy_channel_delta =
      differenceArray(end_state.energy_erg_per_initial_mass_channel, begin_state.energy_erg_per_initial_mass_channel);
  normalizeChannelSums(interval.returned_mass_code / std::max(initial_mass_code, k_mass_floor), return_fraction_channel_delta);
  normalizeChannelSums(interval.returned_metals_code / std::max(initial_mass_code, k_mass_floor), metal_fraction_channel_delta);
  normalizeChannelSums(interval.feedback_energy_erg / std::max(initial_mass_code, k_mass_floor), energy_channel_delta);

  for (std::size_t channel = 0; channel < k_stellar_yield_channel_count; ++channel) {
    interval.returned_mass_channel_code[channel] = return_fraction_channel_delta[channel] * initial_mass_code;
    interval.returned_metals_channel_code[channel] = metal_fraction_channel_delta[channel] * initial_mass_code;
    interval.feedback_energy_channel_erg[channel] = energy_channel_delta[channel] * initial_mass_code;
  }

  return interval;
}

StellarEvolutionTable StellarEvolutionTable::loadFromTextFile(
    const std::string& path,
    const std::string& source_tag) {
  std::ifstream input(path);
  if (!input) {
    throw std::runtime_error("StellarEvolutionTable: failed to open table path '" + path + "'");
  }

  StellarEvolutionTable table;
  table.source_path = path;

  std::string line;
  while (std::getline(input, line)) {
    const std::string normalized = trim(line);
    if (normalized.empty()) {
      continue;
    }

    if (normalized[0] == '#') {
      const std::string metadata = trim(normalized.substr(1));
      const std::size_t equal_pos = metadata.find('=');
      if (equal_pos != std::string::npos) {
        const std::string key = trim(metadata.substr(0, equal_pos));
        const std::string value = trim(metadata.substr(equal_pos + 1));
        if (key == "table_id") {
          table.table_id = value;
        } else if (key == "table_version") {
          table.table_version = value;
        }
      }
      continue;
    }

    std::istringstream row_stream(normalized);
    double age_yr_row = 0.0;
    double return_total_row = 0.0;
    double metal_total_row = 0.0;
    double energy_total_row = 0.0;
    std::array<double, k_stellar_yield_channel_count> return_channel_row{};
    std::array<double, k_stellar_yield_channel_count> metal_channel_row{};
    std::array<double, k_stellar_yield_channel_count> energy_channel_row{};

    row_stream >> age_yr_row >> return_total_row >> metal_total_row >> energy_total_row >> return_channel_row[0] >>
        return_channel_row[1] >> return_channel_row[2] >> metal_channel_row[0] >> metal_channel_row[1] >>
        metal_channel_row[2] >> energy_channel_row[0] >> energy_channel_row[1] >> energy_channel_row[2];
    row_stream >> std::ws;
    if (!row_stream || !row_stream.eof()) {
      throw std::runtime_error(
          "StellarEvolutionTable: malformed row in '" + path + "' for source tag '" + source_tag + "'");
    }

    table.age_yr.push_back(age_yr_row);
    table.return_fraction_total.push_back(return_total_row);
    table.metal_yield_fraction_total.push_back(metal_total_row);
    table.energy_erg_per_initial_mass_code.push_back(energy_total_row);
    for (std::size_t channel = 0; channel < k_stellar_yield_channel_count; ++channel) {
      table.return_fraction_channel[channel].push_back(return_channel_row[channel]);
      table.metal_yield_fraction_channel[channel].push_back(metal_channel_row[channel]);
      table.energy_erg_per_initial_mass_channel[channel].push_back(energy_channel_row[channel]);
    }
  }

  if (!table.isConsistent()) {
    throw std::runtime_error("StellarEvolutionTable: loaded table is inconsistent for path '" + path + "'");
  }

  return table;
}

StellarEvolutionTable StellarEvolutionTable::makeBuiltinReference() {
  StellarEvolutionTable table;
  table.table_id = "builtin_reference";
  table.table_version = "v1_2026_04";
  table.source_path = "builtin";

  const std::array<double, 7> age_yr = {0.0, 3.0e6, 1.0e7, 1.0e8, 1.0e9, 5.0e9, 1.4e10};
  const std::array<double, 7> return_total = {0.0, 0.02, 0.10, 0.18, 0.28, 0.34, 0.40};
  const std::array<double, 7> metal_total = {0.0, 0.001, 0.004, 0.008, 0.013, 0.017, 0.021};
  const std::array<double, 7> energy_total = {0.0, 1.0e46, 6.0e48, 1.0e49, 1.6e49, 1.9e49, 2.1e49};

  const std::array<std::array<double, 7>, k_stellar_yield_channel_count> return_channel = {
      std::array<double, 7>{0.0, 0.015, 0.055, 0.110, 0.170, 0.205, 0.240},
      std::array<double, 7>{0.0, 0.005, 0.045, 0.060, 0.070, 0.070, 0.070},
      std::array<double, 7>{0.0, 0.000, 0.000, 0.010, 0.040, 0.065, 0.090}};

  const std::array<std::array<double, 7>, k_stellar_yield_channel_count> metal_channel = {
      std::array<double, 7>{0.0, 0.0003, 0.0013, 0.0033, 0.0060, 0.0080, 0.0100},
      std::array<double, 7>{0.0, 0.0007, 0.0027, 0.0032, 0.0030, 0.0028, 0.0025},
      std::array<double, 7>{0.0, 0.0000, 0.0000, 0.0015, 0.0040, 0.0062, 0.0085}};

  const std::array<std::array<double, 7>, k_stellar_yield_channel_count> energy_channel = {
      std::array<double, 7>{0.0, 1.0e46, 3.5e47, 7.0e47, 1.2e48, 1.6e48, 2.0e48},
      std::array<double, 7>{0.0, 0.0, 5.65e48, 8.3e48, 9.3e48, 9.5e48, 9.6e48},
      std::array<double, 7>{0.0, 0.0, 0.0, 1.0e48, 5.0e48, 7.9e48, 9.5e48}};

  table.age_yr.assign(age_yr.begin(), age_yr.end());
  table.return_fraction_total.assign(return_total.begin(), return_total.end());
  table.metal_yield_fraction_total.assign(metal_total.begin(), metal_total.end());
  table.energy_erg_per_initial_mass_code.assign(energy_total.begin(), energy_total.end());
  for (std::size_t channel = 0; channel < k_stellar_yield_channel_count; ++channel) {
    table.return_fraction_channel[channel].assign(return_channel[channel].begin(), return_channel[channel].end());
    table.metal_yield_fraction_channel[channel].assign(metal_channel[channel].begin(), metal_channel[channel].end());
    table.energy_erg_per_initial_mass_channel[channel].assign(
        energy_channel[channel].begin(), energy_channel[channel].end());
  }
  return table;
}

StellarEvolutionBookkeeper::StellarEvolutionBookkeeper(StellarEvolutionConfig config, StellarEvolutionTable table)
    : m_config(std::move(config)), m_table(std::move(table)) {
  if (!m_table.isConsistent()) {
    throw std::invalid_argument("StellarEvolutionBookkeeper: stellar evolution table is inconsistent");
  }
  if (m_config.hubble_time_years <= 0.0) {
    throw std::invalid_argument("StellarEvolutionBookkeeper: hubble_time_years must be > 0");
  }
}

const StellarEvolutionConfig& StellarEvolutionBookkeeper::config() const noexcept {
  return m_config;
}

const StellarEvolutionTable& StellarEvolutionBookkeeper::table() const noexcept {
  return m_table;
}

double StellarEvolutionBookkeeper::evaluateStarAgeYears(
    double formation_scale_factor,
    double current_scale_factor) const {
  const double safe_birth_a = std::clamp(formation_scale_factor, 1.0e-6, 1.0);
  const double safe_current_a = std::clamp(current_scale_factor, safe_birth_a, 1.0);
  const double age_proxy = std::log(safe_current_a / safe_birth_a) * m_config.hubble_time_years;
  return std::max(age_proxy, 0.0);
}

StellarEvolutionStepReport StellarEvolutionBookkeeper::apply(
    core::SimulationState& state,
    std::span<const std::uint32_t> active_star_indices,
    double current_scale_factor,
    double dt_code) const {
  StellarEvolutionStepReport report;
  if (!m_config.enabled || dt_code <= 0.0 || state.star_particles.size() == 0) {
    state.sidecars.upsert(buildMetadataSidecar(report));
    return report;
  }

  const double dt_years = std::max(dt_code * m_config.hubble_time_years, 0.0);
  for (const std::uint32_t star_index : active_star_indices) {
    ++report.counters.scanned_stars;
    if (star_index >= state.star_particles.size()) {
      continue;
    }

    const std::uint32_t particle_index = state.star_particles.particle_index[star_index];
    if (particle_index >= state.particles.size()) {
      continue;
    }

    const double birth_mass = state.star_particles.birth_mass_code[star_index];
    if (birth_mass <= k_mass_floor) {
      continue;
    }

    const double age_begin_years = evaluateStarAgeYears(
        state.star_particles.formation_scale_factor[star_index], current_scale_factor);
    const double age_end_years = age_begin_years + dt_years;
    StellarEvolutionIntervalBudget interval = m_table.integrateInterval(age_begin_years, age_end_years, birth_mass);

    const double mass_old = state.particles.mass_code[particle_index];
    const double returned_mass = std::clamp(interval.returned_mass_code, 0.0, mass_old);
    const double mass_new = std::max(mass_old - returned_mass, 0.0);
    interval.returned_mass_code = returned_mass;
    interval.remnant_change_code = mass_old - mass_new - returned_mass;

    state.particles.mass_code[particle_index] = mass_new;
    state.star_particles.stellar_age_years_last[star_index] = age_end_years;
    state.star_particles.stellar_returned_mass_cumulative_code[star_index] += returned_mass;
    state.star_particles.stellar_returned_metals_cumulative_code[star_index] += interval.returned_metals_code;
    state.star_particles.stellar_feedback_energy_cumulative_erg[star_index] += interval.feedback_energy_erg;
    for (std::size_t channel = 0; channel < k_stellar_yield_channel_count; ++channel) {
      state.star_particles.stellar_returned_mass_channel_cumulative_code[channel][star_index] +=
          interval.returned_mass_channel_code[channel];
      state.star_particles.stellar_returned_metals_channel_cumulative_code[channel][star_index] +=
          interval.returned_metals_channel_code[channel];
      state.star_particles.stellar_feedback_energy_channel_cumulative_erg[channel][star_index] +=
          interval.feedback_energy_channel_erg[channel];
    }

    ++report.counters.evolved_stars;
    report.counters.returned_mass_code += returned_mass;
    report.counters.returned_metals_code += interval.returned_metals_code;
    report.counters.feedback_energy_erg += interval.feedback_energy_erg;

    report.budgets.push_back(StellarEvolutionStarBudget{
        .star_index = star_index,
        .particle_index = particle_index,
        .star_age_begin_years = age_begin_years,
        .star_age_end_years = age_end_years,
        .mass_old_code = mass_old,
        .mass_new_code = mass_new,
        .interval = interval,
    });
  }

  state.sidecars.upsert(buildMetadataSidecar(report));
  return report;
}

core::ModuleSidecarBlock StellarEvolutionBookkeeper::buildMetadataSidecar(
    const StellarEvolutionStepReport& report) const {
  std::ostringstream stream;
  stream << "module=stellar_evolution\n";
  stream << "schema_version=" << m_config.metadata_schema_version << "\n";
  stream << "table_id=" << m_table.table_id << "\n";
  stream << "table_version=" << m_table.table_version << "\n";
  stream << "table_source_path=" << m_table.source_path << "\n";
  stream << "scanned_stars=" << report.counters.scanned_stars << "\n";
  stream << "evolved_stars=" << report.counters.evolved_stars << "\n";
  stream << "returned_mass_code=" << report.counters.returned_mass_code << "\n";
  stream << "returned_metals_code=" << report.counters.returned_metals_code << "\n";
  stream << "feedback_energy_erg=" << report.counters.feedback_energy_erg << "\n";

  const std::string text = stream.str();
  core::ModuleSidecarBlock block;
  block.module_name = "stellar_evolution";
  block.schema_version = m_config.metadata_schema_version;
  block.payload.resize(text.size());
  for (std::size_t i = 0; i < text.size(); ++i) {
    block.payload[i] = static_cast<std::byte>(text[i]);
  }
  return block;
}

StellarEvolutionConfig makeStellarEvolutionConfig(const core::PhysicsConfig& physics_config) {
  StellarEvolutionConfig config;
  config.enabled = physics_config.enable_stellar_evolution;
  config.table_path = physics_config.stellar_evolution_table_path;
  config.hubble_time_years = physics_config.stellar_evolution_hubble_time_years;
  return config;
}

StellarEvolutionTable loadStellarEvolutionTable(const core::PhysicsConfig& physics_config) {
  if (physics_config.stellar_evolution_table_path.empty()) {
    return StellarEvolutionTable::makeBuiltinReference();
  }
  return StellarEvolutionTable::loadFromTextFile(
      physics_config.stellar_evolution_table_path,
      "physics.stellar_evolution_table_path");
}

}  // namespace cosmosim::physics
