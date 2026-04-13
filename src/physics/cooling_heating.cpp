#include "cosmosim/physics/cooling_heating.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace cosmosim::physics {
namespace {

constexpr double k_small = 1.0e-30;

[[nodiscard]] double clampTemperature(double temperature_k) {
  return std::max(temperature_k, 1.0);
}

[[nodiscard]] double safeLog10(double value) {
  return std::log10(std::max(value, k_small));
}

[[nodiscard]] double linearInterpolate(double x0, double x1, double y0, double y1, double x) {
  if (std::abs(x1 - x0) < 1.0e-16) {
    return y0;
  }
  const double t = (x - x0) / (x1 - x0);
  return (1.0 - t) * y0 + t * y1;
}

[[nodiscard]] double uvHeatingScale(core::UvBackgroundModel uv_background_model, double redshift) {
  switch (uv_background_model) {
    case core::UvBackgroundModel::kNone:
      return 0.0;
    case core::UvBackgroundModel::kHm12: {
      const double z_peak = 2.0;
      const double width = 2.5;
      const double dz = (redshift - z_peak) / width;
      return std::exp(-dz * dz);
    }
    case core::UvBackgroundModel::kFg20:
      return 1.0 / (1.0 + std::pow((1.0 + redshift) / 4.0, 2.0));
  }
  throw std::invalid_argument("CoolingRateProvider: unhandled UvBackgroundModel enum value");
}

[[nodiscard]] double shieldingFactor(
    core::SelfShieldingModel self_shielding_model,
    double hydrogen_number_density_cgs) {
  switch (self_shielding_model) {
    case core::SelfShieldingModel::kNone:
      return 1.0;
    case core::SelfShieldingModel::kRahmati13Like: {
      const double threshold = 1.0e-2;
      return 1.0 / (1.0 + std::pow(hydrogen_number_density_cgs / threshold, 1.6));
    }
  }
  throw std::invalid_argument("CoolingRateProvider: unhandled SelfShieldingModel enum value");
}

[[nodiscard]] double effectiveHydrogenNumberDensity(
    double mass_density_comoving,
    double default_hydrogen_number_density_cgs) {
  if (default_hydrogen_number_density_cgs > 0.0) {
    return default_hydrogen_number_density_cgs;
  }
  return std::max(mass_density_comoving, 1.0e-10);
}

}  // namespace

double CoolingRateResult::netCoolingRateErgCm3S() const {
  return cooling_rate_erg_cm3_s - heating_rate_erg_cm3_s;
}

MetalLineCoolingTable MetalLineCoolingTable::loadFromTextFile(
    const std::string& path,
    const std::string& table_label) {
  std::ifstream input(path);
  if (!input) {
    throw std::runtime_error("MetalLineCoolingTable: failed to open file '" + path + "'");
  }

  MetalLineCoolingTable table;
  table.m_provenance.source_path = path;
  table.m_provenance.source_format = "text_two_column";
  table.m_provenance.table_label = table_label;

  std::string line;
  while (std::getline(input, line)) {
    if (line.empty() || line.front() == '#') {
      continue;
    }
    std::istringstream stream(line);
    double log10_temperature_k = 0.0;
    double log10_cooling_rate_erg_cm3_s = 0.0;
    stream >> log10_temperature_k >> log10_cooling_rate_erg_cm3_s;
    if (stream.fail()) {
      throw std::runtime_error("MetalLineCoolingTable: malformed row in '" + path + "'");
    }
    table.m_log10_temperature_k.push_back(log10_temperature_k);
    table.m_log10_cooling_rate_erg_cm3_s.push_back(log10_cooling_rate_erg_cm3_s);
  }

  if (table.m_log10_temperature_k.empty()) {
    throw std::runtime_error("MetalLineCoolingTable: no data rows in '" + path + "'");
  }

  for (std::size_t i = 1; i < table.m_log10_temperature_k.size(); ++i) {
    if (table.m_log10_temperature_k[i] <= table.m_log10_temperature_k[i - 1]) {
      throw std::runtime_error("MetalLineCoolingTable: log10(T) must be strictly increasing");
    }
  }

  table.m_provenance.row_count = table.m_log10_temperature_k.size();
  return table;
}

bool MetalLineCoolingTable::empty() const noexcept {
  return m_log10_temperature_k.empty();
}

double MetalLineCoolingTable::lookupCoolingRateErgCm3S(double temperature_k) const {
  if (m_log10_temperature_k.empty()) {
    return 0.0;
  }

  const double target = safeLog10(clampTemperature(temperature_k));
  const auto upper = std::upper_bound(m_log10_temperature_k.begin(), m_log10_temperature_k.end(), target);
  if (upper == m_log10_temperature_k.begin()) {
    return std::pow(10.0, m_log10_cooling_rate_erg_cm3_s.front());
  }
  if (upper == m_log10_temperature_k.end()) {
    return std::pow(10.0, m_log10_cooling_rate_erg_cm3_s.back());
  }

  const std::size_t right_index = static_cast<std::size_t>(std::distance(m_log10_temperature_k.begin(), upper));
  const std::size_t left_index = right_index - 1;
  const double interpolated_log10_rate = linearInterpolate(
      m_log10_temperature_k[left_index],
      m_log10_temperature_k[right_index],
      m_log10_cooling_rate_erg_cm3_s[left_index],
      m_log10_cooling_rate_erg_cm3_s[right_index],
      target);
  return std::pow(10.0, interpolated_log10_rate);
}

const CoolingTableProvenance& MetalLineCoolingTable::provenance() const noexcept {
  return m_provenance;
}

CoolingRateProvider::CoolingRateProvider(CoolingModelConfig config) : m_config(std::move(config)) {
  if (m_config.temperature_floor_k <= 0.0) {
    throw std::invalid_argument("CoolingRateProvider: temperature_floor_k must be > 0");
  }
  if (m_config.max_fractional_energy_change_per_substep <= 0.0 ||
      m_config.max_fractional_energy_change_per_substep >= 1.0) {
    throw std::invalid_argument(
        "CoolingRateProvider: max_fractional_energy_change_per_substep must be in (0, 1)");
  }
  if (m_config.max_subcycles <= 0) {
    throw std::invalid_argument("CoolingRateProvider: max_subcycles must be > 0");
  }
  if (m_config.temperature_per_internal_energy_k <= 0.0) {
    throw std::invalid_argument("CoolingRateProvider: temperature_per_internal_energy_k must be > 0");
  }
}

void CoolingRateProvider::setMetalLineTable(MetalLineCoolingTable table) {
  m_metal_line_table = std::move(table);
}

bool CoolingRateProvider::hasMetalLineTable() const noexcept {
  return m_metal_line_table.has_value();
}

CoolingRateResult CoolingRateProvider::lookupRates(const CoolingRateQuery& query) const {
  const double temperature_k = std::max(query.temperature_k, m_config.temperature_floor_k);
  const double n_h = std::max(query.hydrogen_number_density_cgs, 1.0e-12);

  // Primordial baseline: compact analytic fit with explicit UV/shielding controls.
  const double primordial_cooling_coeff = 1.0e-22 * std::sqrt(temperature_k / 1.0e4);
  const double primordial_cooling = primordial_cooling_coeff * n_h * n_h;

  const double uv_scale = uvHeatingScale(m_config.uv_background_model, query.redshift);
  const double shielding = shieldingFactor(m_config.self_shielding_model, n_h);
  const double primordial_heating = 2.0e-24 * uv_scale * shielding * n_h;

  double metal_line_cooling = 0.0;
  if (m_config.enable_metal_line_cooling && m_metal_line_table.has_value()) {
    metal_line_cooling = m_metal_line_table->lookupCoolingRateErgCm3S(temperature_k) *
        std::max(query.metallicity_mass_fraction, 0.0);
  }

  CoolingRateResult result;
  result.cooling_rate_erg_cm3_s = primordial_cooling + metal_line_cooling;
  result.heating_rate_erg_cm3_s = primordial_heating;
  return result;
}

const CoolingModelConfig& CoolingRateProvider::config() const noexcept {
  return m_config;
}

std::optional<CoolingTableProvenance> CoolingRateProvider::metalLineTableProvenance() const {
  if (!m_metal_line_table.has_value()) {
    return std::nullopt;
  }
  return m_metal_line_table->provenance();
}

CoolingSourceIntegrator::CoolingSourceIntegrator(double energy_floor_code)
    : m_energy_floor_code(energy_floor_code) {
  if (m_energy_floor_code <= 0.0) {
    throw std::invalid_argument("CoolingSourceIntegrator: energy_floor_code must be > 0");
  }
}

CoolingIntegrationResult CoolingSourceIntegrator::integrateSpecificInternalEnergy(
    double specific_internal_energy_code,
    double mass_density_comoving,
    double dt_code,
    const CoolingRateQuery& query,
    const CoolingRateProvider& rate_provider) const {
  if (dt_code <= 0.0) {
    throw std::invalid_argument("CoolingSourceIntegrator: dt_code must be > 0");
  }

  CoolingIntegrationResult result;
  double energy = std::max(specific_internal_energy_code, m_energy_floor_code);
  double t_remaining = dt_code;
  const CoolingModelConfig& config = rate_provider.config();

  for (int subcycle = 0; subcycle < config.max_subcycles && t_remaining > 0.0; ++subcycle) {
    CoolingRateQuery local_query = query;
    local_query.temperature_k = std::max(
        energy * config.temperature_per_internal_energy_k,
        config.temperature_floor_k);

    const CoolingRateResult rates = rate_provider.lookupRates(local_query);
    const double volumetric_net_rate = rates.heating_rate_erg_cm3_s - rates.cooling_rate_erg_cm3_s;
    const double de_dt = volumetric_net_rate / std::max(mass_density_comoving, 1.0e-16);

    const double cooling_time = std::max(energy, m_energy_floor_code) / std::max(std::abs(de_dt), 1.0e-30);
    const double stable_dt = std::max(
        config.max_fractional_energy_change_per_substep * cooling_time,
        1.0e-14 * dt_code);
    const double dt_sub = std::min(t_remaining, stable_dt);

    energy = std::max(energy + de_dt * dt_sub, m_energy_floor_code);
    t_remaining -= dt_sub;
    result.diagnostics.subcycles_used = subcycle + 1;
  }

  result.diagnostics.hit_subcycle_limit = (t_remaining > 1.0e-14 * dt_code);

  CoolingRateQuery final_query = query;
  final_query.temperature_k = std::max(
      energy * config.temperature_per_internal_energy_k,
      config.temperature_floor_k);
  const CoolingRateResult final_rates = rate_provider.lookupRates(final_query);
  const double final_de_dt = (final_rates.heating_rate_erg_cm3_s - final_rates.cooling_rate_erg_cm3_s) /
      std::max(mass_density_comoving, 1.0e-16);
  result.diagnostics.suggested_next_dt_code = config.max_fractional_energy_change_per_substep *
      std::max(energy, m_energy_floor_code) / std::max(std::abs(final_de_dt), 1.0e-30);
  result.specific_internal_energy_code = energy;
  return result;
}

CoolingHeatingSource::CoolingHeatingSource(CoolingRateProvider provider, CoolingSourceIntegrator integrator)
    : m_rate_provider(std::move(provider)), m_integrator(std::move(integrator)) {}

hydro::HydroConservedState CoolingHeatingSource::sourceForCell(
    std::size_t cell_index,
    const hydro::HydroConservedState& conserved,
    const hydro::HydroPrimitiveState& primitive,
    const hydro::HydroSourceContext& context) const {
  const double kinetic_density = 0.5 *
      (conserved.momentum_density_x_comoving * conserved.momentum_density_x_comoving +
       conserved.momentum_density_y_comoving * conserved.momentum_density_y_comoving +
       conserved.momentum_density_z_comoving * conserved.momentum_density_z_comoving) /
      std::max(conserved.mass_density_comoving, 1.0e-16);
  const double internal_density = std::max(conserved.total_energy_density_comoving - kinetic_density, 0.0);
  const double specific_internal_energy = internal_density / std::max(conserved.mass_density_comoving, 1.0e-16);

  const double hydrogen_density =
      (cell_index < context.hydrogen_number_density_cgs.size()) ? context.hydrogen_number_density_cgs[cell_index]
                                                                 : effectiveHydrogenNumberDensity(primitive.rho_comoving, 0.0);
  const double metallicity =
      (cell_index < context.metallicity_mass_fraction.size()) ? context.metallicity_mass_fraction[cell_index] : 0.0;
  const double temperature_k =
      (cell_index < context.temperature_k.size()) ? context.temperature_k[cell_index] : 1.0e4;

  const CoolingRateQuery query{
      .temperature_k = temperature_k,
      .hydrogen_number_density_cgs = hydrogen_density,
      .metallicity_mass_fraction = metallicity,
      .redshift = context.redshift};

  const CoolingIntegrationResult integrated = m_integrator.integrateSpecificInternalEnergy(
      specific_internal_energy,
      conserved.mass_density_comoving,
      context.update.dt_code,
      query,
      m_rate_provider);

  hydro::HydroConservedState source{};
  source.total_energy_density_comoving =
      (integrated.specific_internal_energy_code - specific_internal_energy) *
      conserved.mass_density_comoving / std::max(context.update.dt_code, 1.0e-30);
  return source;
}

}  // namespace cosmosim::physics
