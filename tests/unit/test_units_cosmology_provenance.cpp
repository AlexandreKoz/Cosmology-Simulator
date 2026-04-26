#include <cassert>
#include <cmath>
#include <string>

#include "cosmosim/core/constants.hpp"
#include "cosmosim/core/config.hpp"
#include "cosmosim/core/cosmology.hpp"
#include "cosmosim/core/provenance.hpp"
#include "cosmosim/core/units.hpp"

namespace {

constexpr double k_tolerance = 1.0e-12;

void testHubbleAndCriticalDensityAtAOne() {
  cosmosim::core::CosmologyBackgroundConfig cfg;
  cfg.hubble_param = 0.7;
  cfg.omega_matter = 0.3;
  cfg.omega_lambda = 0.7;
  cfg.omega_radiation = 0.0;
  cfg.omega_curvature = 0.0;

  cosmosim::core::LambdaCdmBackground background(cfg);

  const double expected_h0 = 0.7 * cosmosim::core::constants::k_hubble_100_km_s_mpc_si;
  const double h_a1 = background.hubbleSi(1.0);
  assert(std::abs((h_a1 - expected_h0) / expected_h0) < k_tolerance);

  const double rho_crit = background.criticalDensitySi(1.0);
  const double expected_rho =
      (3.0 * expected_h0 * expected_h0) /
      (8.0 * cosmosim::core::constants::k_pi * cosmosim::core::constants::k_newton_g_si);
  assert(std::abs((rho_crit - expected_rho) / expected_rho) < k_tolerance);
}

void testScaleFactorDependence() {
  cosmosim::core::CosmologyBackgroundConfig cfg;
  cfg.hubble_param = 0.67;
  cfg.omega_matter = 0.31;
  cfg.omega_lambda = 0.68;
  cfg.omega_radiation = 0.01;
  cfg.omega_curvature = 0.0;

  cosmosim::core::LambdaCdmBackground background(cfg);
  const double a = 0.5;
  const double e = std::sqrt(
      cfg.omega_radiation * std::pow(a, -4.0) + cfg.omega_matter * std::pow(a, -3.0) +
      cfg.omega_curvature * std::pow(a, -2.0) + cfg.omega_lambda);
  assert(std::abs(background.eFactor(a) - e) < 1.0e-12);
}

void testComovingPhysicalConversions() {
  const double a = 0.25;
  const double x_comoving = 16.0;
  const double r_phys = cosmosim::core::comovingToPhysicalLength(x_comoving, a);
  assert(std::abs(r_phys - 4.0) < k_tolerance);
  assert(std::abs(cosmosim::core::physicalToComovingLength(r_phys, a) - x_comoving) < k_tolerance);

  const double dx_dt = 20.0;
  const double v_pec = cosmosim::core::peculiarVelocityFromComovingRate(dx_dt, a);
  assert(std::abs(v_pec - 5.0) < k_tolerance);
  assert(std::abs(cosmosim::core::comovingRateFromPeculiarVelocity(v_pec, a) - dx_dt) <
         k_tolerance);
}

void testUnitsConversions() {
  const auto units = cosmosim::core::makeUnitSystem("mpc", "msun", "km_s");

  const double length_si = units.lengthCodeToSi(1.0);
  assert(std::abs(length_si - cosmosim::core::constants::k_megaparsec_si) < 1.0e-3);

  const double rho_code = 2.5;
  const double rho_si = units.densityCodeToSi(rho_code);
  const double rho_roundtrip = units.densitySiToCode(rho_si);
  assert(std::abs(rho_roundtrip - rho_code) < k_tolerance);

  const double rho_cgs = units.densityCodeToCgs(1.0);
  assert(rho_cgs > 0.0);
}

void testStableHash() {
  const std::string normalized = "a=1\nb=2\n";
  const std::string hash_a = cosmosim::core::stableConfigHashHex(normalized);
  const std::string hash_b = cosmosim::core::stableConfigHashHex(normalized);
  assert(hash_a == hash_b);
}

void testDerivedConstantsConsistencyFromNormalizedConfig() {
  const std::string text = R"(
[mode]
mode = cosmo_cube
[units]
length_unit = mpc
mass_unit = msun
velocity_unit = km_s
[cosmology]
hubble_param = 0.7
omega_matter = 0.3
omega_lambda = 0.7
)";
  const auto frozen = cosmosim::core::loadFrozenConfigFromString(text, "derived_constants");
  const auto units_from_frozen = cosmosim::core::makeUnitSystem(
      frozen.config.units.length_unit,
      frozen.config.units.mass_unit,
      frozen.config.units.velocity_unit);
  const auto units_from_roundtrip = cosmosim::core::makeUnitSystem("mpc", "msun", "km_s");
  assert(std::abs(units_from_frozen.length_si_per_code - units_from_roundtrip.length_si_per_code) < k_tolerance);

  cosmosim::core::CosmologyBackgroundConfig bg_cfg;
  bg_cfg.hubble_param = frozen.config.cosmology.hubble_param;
  bg_cfg.omega_matter = frozen.config.cosmology.omega_matter;
  bg_cfg.omega_lambda = frozen.config.cosmology.omega_lambda;
  bg_cfg.omega_radiation = 0.0;
  bg_cfg.omega_curvature = 0.0;
  const cosmosim::core::LambdaCdmBackground background(bg_cfg);
  const double expected_h0 = frozen.config.cosmology.hubble_param * cosmosim::core::constants::k_hubble_100_km_s_mpc_si;
  assert(std::abs(background.hubble0Si() - expected_h0) / expected_h0 < k_tolerance);
}

}  // namespace

int main() {
  testHubbleAndCriticalDensityAtAOne();
  testScaleFactorDependence();
  testComovingPhysicalConversions();
  testUnitsConversions();
  testStableHash();
  testDerivedConstantsConsistencyFromNormalizedConfig();
  return 0;
}
