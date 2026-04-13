#include "cosmosim/core/units.hpp"

#include <algorithm>
#include <cctype>
#include <stdexcept>

#include "cosmosim/core/constants.hpp"

namespace cosmosim::core {
namespace {

[[nodiscard]] std::string toLower(std::string value) {
  std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
    return static_cast<char>(std::tolower(c));
  });
  return value;
}

[[nodiscard]] double lengthSiPerCodeFromName(const std::string& name) {
  const std::string unit = toLower(name);
  if (unit == "m") {
    return 1.0;
  }
  if (unit == "cm") {
    return constants::k_centimeter_si;
  }
  if (unit == "kpc") {
    return constants::k_kiloparsec_si;
  }
  if (unit == "mpc") {
    return constants::k_megaparsec_si;
  }
  throw std::invalid_argument("unsupported length unit: " + name);
}

[[nodiscard]] double massSiPerCodeFromName(const std::string& name) {
  const std::string unit = toLower(name);
  if (unit == "kg") {
    return 1.0;
  }
  if (unit == "g") {
    return constants::k_gram_si;
  }
  if (unit == "msun") {
    return constants::k_solar_mass_si;
  }
  throw std::invalid_argument("unsupported mass unit: " + name);
}

[[nodiscard]] double velocitySiPerCodeFromName(const std::string& name) {
  const std::string unit = toLower(name);
  if (unit == "m_s") {
    return 1.0;
  }
  if (unit == "cm_s") {
    return constants::k_centimeter_si;
  }
  if (unit == "km_s") {
    return constants::k_kilometer_si;
  }
  throw std::invalid_argument("unsupported velocity unit: " + name);
}

}  // namespace

double UnitSystem::lengthCodeToSi(double length_code) const { return length_code * length_si_per_code; }

double UnitSystem::lengthSiToCode(double length_si) const { return length_si / length_si_per_code; }

double UnitSystem::massCodeToSi(double mass_code) const { return mass_code * mass_si_per_code; }

double UnitSystem::massSiToCode(double mass_si) const { return mass_si / mass_si_per_code; }

double UnitSystem::velocityCodeToSi(double velocity_code) const {
  return velocity_code * velocity_si_per_code;
}

double UnitSystem::velocitySiToCode(double velocity_si) const { return velocity_si / velocity_si_per_code; }

double UnitSystem::densityCodeToSi(double density_code) const {
  const double density_si_per_code = mass_si_per_code / (length_si_per_code * length_si_per_code * length_si_per_code);
  return density_code * density_si_per_code;
}

double UnitSystem::densitySiToCode(double density_si) const {
  const double density_si_per_code = mass_si_per_code / (length_si_per_code * length_si_per_code * length_si_per_code);
  return density_si / density_si_per_code;
}

double UnitSystem::densityCodeToCgs(double density_code) const {
  const double density_si = densityCodeToSi(density_code);
  return density_si * 1.0e-3;
}

UnitSystem makeUnitSystem(
    const std::string& length_unit,
    const std::string& mass_unit,
    const std::string& velocity_unit) {
  UnitSystem units;
  units.length_unit = toLower(length_unit);
  units.mass_unit = toLower(mass_unit);
  units.velocity_unit = toLower(velocity_unit);

  units.length_si_per_code = lengthSiPerCodeFromName(units.length_unit);
  units.mass_si_per_code = massSiPerCodeFromName(units.mass_unit);
  units.velocity_si_per_code = velocitySiPerCodeFromName(units.velocity_unit);
  return units;
}

double comovingToPhysicalLength(double x_comoving, double scale_factor) {
  return scale_factor * x_comoving;
}

double physicalToComovingLength(double r_phys, double scale_factor) { return r_phys / scale_factor; }

double peculiarVelocityFromComovingRate(double dx_comoving_dt, double scale_factor) {
  return scale_factor * dx_comoving_dt;
}

double comovingRateFromPeculiarVelocity(double v_pec, double scale_factor) {
  return v_pec / scale_factor;
}

}  // namespace cosmosim::core
