#pragma once

namespace cosmosim::core {

// LambdaCDM background parameters.
// hubble_param is the dimensionless little-h: H0 = h * 100 km/s/Mpc.
struct CosmologyBackgroundConfig {
  double hubble_param = 0.674;  // h = H0 / (100 km/s/Mpc)
  double omega_matter = 0.315;
  double omega_lambda = 0.685;
  double omega_radiation = 0.0;
  double omega_curvature = 0.0;
};

class LambdaCdmBackground {
 public:
  explicit LambdaCdmBackground(CosmologyBackgroundConfig config = {});

  // Hubble constant at z=0 in SI units [s^-1].
  [[nodiscard]] const CosmologyBackgroundConfig& config() const;
  [[nodiscard]] double hubble0Si() const;
  // E(a) = H(a) / H0.
  [[nodiscard]] double eFactor(double scale_factor) const;
  // H(a) = H0 * sqrt(omega_r a^-4 + omega_m a^-3 + omega_k a^-2 + omega_lambda).
  [[nodiscard]] double hubbleSi(double scale_factor) const;
  // rho_crit(a) = 3 H(a)^2 / (8 pi G), SI units [kg m^-3].
  [[nodiscard]] double criticalDensitySi(double scale_factor) const;
  // Exact FLRW elapsed proper time: integral da / (a H(a)), SI seconds.
  [[nodiscard]] double cosmicTimeIntervalSi(
      double scale_factor_begin,
      double scale_factor_end) const;
  // Proper cosmic time since a configurable positive lower bound.
  [[nodiscard]] double cosmicTimeSinceScaleFactorSi(
      double scale_factor,
      double lower_scale_factor = 1.0e-8) const;

 private:
  CosmologyBackgroundConfig m_config;
};

}  // namespace cosmosim::core
