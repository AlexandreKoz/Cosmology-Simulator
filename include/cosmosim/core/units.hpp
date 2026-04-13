#pragma once

#include <string>

namespace cosmosim::core {

// Canonical code-unit mapping used by low-level physics helpers.
// The *_si_per_code members define one code unit expressed in SI.
struct UnitSystem {
  std::string length_unit;
  std::string mass_unit;
  std::string velocity_unit;

  double length_si_per_code = 1.0;
  double mass_si_per_code = 1.0;
  double velocity_si_per_code = 1.0;

  [[nodiscard]] double lengthCodeToSi(double length_code) const;
  [[nodiscard]] double lengthSiToCode(double length_si) const;

  [[nodiscard]] double massCodeToSi(double mass_code) const;
  [[nodiscard]] double massSiToCode(double mass_si) const;

  [[nodiscard]] double velocityCodeToSi(double velocity_code) const;
  [[nodiscard]] double velocitySiToCode(double velocity_si) const;

  [[nodiscard]] double densityCodeToSi(double density_code) const;
  [[nodiscard]] double densitySiToCode(double density_si) const;

  [[nodiscard]] double densityCodeToCgs(double density_code) const;
};

[[nodiscard]] UnitSystem makeUnitSystem(
    const std::string& length_unit,
    const std::string& mass_unit,
    const std::string& velocity_unit);

// Frame-conversion helpers:
// r_phys = a * x_comoving
[[nodiscard]] double comovingToPhysicalLength(double x_comoving, double scale_factor);
[[nodiscard]] double physicalToComovingLength(double r_phys, double scale_factor);
// Peculiar velocity conversion:
// v_pec = a * dx_comoving/dt
[[nodiscard]] double peculiarVelocityFromComovingRate(double dx_comoving_dt, double scale_factor);
[[nodiscard]] double comovingRateFromPeculiarVelocity(double v_pec, double scale_factor);

}  // namespace cosmosim::core
