#include "cosmosim/core/cosmology.hpp"

#include <cmath>
#include <stdexcept>

#include "cosmosim/core/constants.hpp"

namespace cosmosim::core {

LambdaCdmBackground::LambdaCdmBackground(CosmologyBackgroundConfig config) : m_config(config) {
  if (m_config.hubble_param <= 0.0) {
    throw std::invalid_argument("hubble_param must be positive");
  }
}

const CosmologyBackgroundConfig& LambdaCdmBackground::config() const { return m_config; }

double LambdaCdmBackground::hubble0Si() const {
  return m_config.hubble_param * constants::k_hubble_100_km_s_mpc_si;
}

double LambdaCdmBackground::eFactor(double scale_factor) const {
  if (scale_factor <= 0.0) {
    throw std::invalid_argument("scale_factor must be positive");
  }
  const double a2 = scale_factor * scale_factor;
  const double a3 = a2 * scale_factor;
  const double a4 = a2 * a2;

  // E(a)^2 = omega_r a^-4 + omega_m a^-3 + omega_k a^-2 + omega_lambda.
  const double density_sum = m_config.omega_radiation / a4 + m_config.omega_matter / a3 +
                             m_config.omega_curvature / a2 + m_config.omega_lambda;
  if (density_sum <= 0.0) {
    throw std::invalid_argument("density sum for H(a) must be positive");
  }
  return std::sqrt(density_sum);
}

double LambdaCdmBackground::hubbleSi(double scale_factor) const {
  return hubble0Si() * eFactor(scale_factor);
}

double LambdaCdmBackground::criticalDensitySi(double scale_factor) const {
  const double h = hubbleSi(scale_factor);
  // rho_crit(a) = 3 H(a)^2 / (8 pi G).
  return (3.0 * h * h) / (8.0 * constants::k_pi * constants::k_newton_g_si);
}

}  // namespace cosmosim::core
