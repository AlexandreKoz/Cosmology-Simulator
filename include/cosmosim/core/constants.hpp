#pragma once

namespace cosmosim::core::constants {

// Mathematical constant pi.
constexpr double k_pi = 3.141592653589793238462643383279502884;
// Newtonian gravitational constant in SI units [m^3 kg^-1 s^-2].
constexpr double k_newton_g_si = 6.67430e-11;                // m^3 kg^-1 s^-2
// Speed of light in vacuum [m s^-1].
constexpr double k_speed_of_light_si = 2.99792458e8;         // m s^-1
// Parsec in SI units [m].
constexpr double k_parsec_si = 3.0856775814913673e16;        // m
// Kiloparsec in SI units [m].
constexpr double k_kiloparsec_si = 1.0e3 * k_parsec_si;      // m
// Megaparsec in SI units [m].
constexpr double k_megaparsec_si = 1.0e6 * k_parsec_si;      // m
// Solar mass in SI units [kg].
constexpr double k_solar_mass_si = 1.98847e30;               // kg
// Kilometer in SI units [m].
constexpr double k_kilometer_si = 1.0e3;                     // m
// Centimeter in SI units [m].
constexpr double k_centimeter_si = 1.0e-2;                   // m
// Gram in SI units [kg].
constexpr double k_gram_si = 1.0e-3;                         // kg
// Seconds per megayear [s].
constexpr double k_seconds_per_megayear = 3.15576e13;        // s
// Hubble scale conversion 100 km/s/Mpc expressed in SI [s^-1].
constexpr double k_hubble_100_km_s_mpc_si =
    (100.0 * k_kilometer_si) / k_megaparsec_si;              // s^-1

}  // namespace cosmosim::core::constants
