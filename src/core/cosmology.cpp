#include "cosmosim/core/cosmology.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <stdexcept>

#include "cosmosim/core/constants.hpp"

namespace cosmosim::core {

namespace {

[[nodiscard]] double adaptiveSimpsonRecursive(
    const std::function<double(double)>& integrand,
    double left,
    double right,
    double f_left,
    double f_mid,
    double f_right,
    double whole,
    double tolerance,
    int depth_remaining) {
  const double mid = 0.5 * (left + right);
  const double left_mid = 0.5 * (left + mid);
  const double right_mid = 0.5 * (mid + right);
  const double f_left_mid = integrand(left_mid);
  const double f_right_mid = integrand(right_mid);
  const double left_part = (mid - left) * (f_left + 4.0 * f_left_mid + f_mid) / 6.0;
  const double right_part = (right - mid) * (f_mid + 4.0 * f_right_mid + f_right) / 6.0;
  const double refined = left_part + right_part;
  const double error = refined - whole;
  if (!std::isfinite(refined) || !std::isfinite(error)) {
    throw std::runtime_error("adaptive cosmology integration produced a non-finite estimate");
  }
  if (std::abs(error) <= 15.0 * tolerance) {
    return refined + error / 15.0;
  }
  if (depth_remaining <= 0) {
    throw std::runtime_error("adaptive cosmology integration failed to converge within recursion limit");
  }
  return adaptiveSimpsonRecursive(
             integrand, left, mid, f_left, f_left_mid, f_mid, left_part,
             0.5 * tolerance, depth_remaining - 1) +
         adaptiveSimpsonRecursive(
             integrand, mid, right, f_mid, f_right_mid, f_right, right_part,
             0.5 * tolerance, depth_remaining - 1);
}

[[nodiscard]] double adaptiveSimpson(
    const std::function<double(double)>& integrand,
    double left,
    double right,
    double tolerance) {
  const double mid = 0.5 * (left + right);
  const double f_left = integrand(left);
  const double f_mid = integrand(mid);
  const double f_right = integrand(right);
  const double whole = (right - left) * (f_left + 4.0 * f_mid + f_right) / 6.0;
  return adaptiveSimpsonRecursive(
      integrand, left, right, f_left, f_mid, f_right, whole, tolerance, 24);
}

}  // namespace


LambdaCdmBackground::LambdaCdmBackground(CosmologyBackgroundConfig config) : m_config(config) {
  if (!std::isfinite(m_config.hubble_param) || !(m_config.hubble_param > 0.0)) {
    throw std::invalid_argument("hubble_param must be finite and positive");
  }
  if (!std::isfinite(m_config.omega_matter) || !std::isfinite(m_config.omega_lambda) ||
      !std::isfinite(m_config.omega_radiation) || !std::isfinite(m_config.omega_curvature)) {
    throw std::invalid_argument("cosmology density parameters must be finite");
  }
  if (m_config.omega_matter < 0.0 || m_config.omega_lambda < 0.0 || m_config.omega_radiation < 0.0) {
    throw std::invalid_argument("matter, lambda, and radiation density parameters must be non-negative");
  }
}

const CosmologyBackgroundConfig& LambdaCdmBackground::config() const { return m_config; }

double LambdaCdmBackground::hubble0Si() const {
  return m_config.hubble_param * constants::k_hubble_100_km_s_mpc_si;
}

double LambdaCdmBackground::eFactor(double scale_factor) const {
  if (!std::isfinite(scale_factor) || !(scale_factor > 0.0)) {
    throw std::invalid_argument("scale_factor must be finite and positive");
  }
  const double a2 = scale_factor * scale_factor;
  const double a3 = a2 * scale_factor;
  const double a4 = a2 * a2;

  // E(a)^2 = omega_r a^-4 + omega_m a^-3 + omega_k a^-2 + omega_lambda.
  const double density_sum = m_config.omega_radiation / a4 + m_config.omega_matter / a3 +
                             m_config.omega_curvature / a2 + m_config.omega_lambda;
  if (!std::isfinite(density_sum) || !(density_sum > 0.0)) {
    throw std::invalid_argument("density sum for H(a) must be finite and positive");
  }
  return std::sqrt(density_sum);
}

double LambdaCdmBackground::hubbleSi(double scale_factor) const {
  const double hubble = hubble0Si() * eFactor(scale_factor);
  if (!std::isfinite(hubble) || !(hubble > 0.0)) {
    throw std::runtime_error("H(a) is non-finite or non-positive");
  }
  return hubble;
}

double LambdaCdmBackground::criticalDensitySi(double scale_factor) const {
  const double h = hubbleSi(scale_factor);
  // rho_crit(a) = 3 H(a)^2 / (8 pi G).
  return (3.0 * h * h) / (8.0 * constants::k_pi * constants::k_newton_g_si);
}


double LambdaCdmBackground::cosmicTimeIntervalSi(
    double scale_factor_begin,
    double scale_factor_end) const {
  if (!std::isfinite(scale_factor_begin) || !std::isfinite(scale_factor_end) ||
      scale_factor_begin <= 0.0 || scale_factor_end <= 0.0) {
    throw std::invalid_argument("scale factors must be finite and positive");
  }
  if (scale_factor_begin == scale_factor_end) {
    return 0.0;
  }
  const double sign = scale_factor_end > scale_factor_begin ? 1.0 : -1.0;
  const double left = std::min(scale_factor_begin, scale_factor_end);
  const double right = std::max(scale_factor_begin, scale_factor_end);
  const auto integrand = [this](double scale_factor) {
    return 1.0 / (scale_factor * hubbleSi(scale_factor));
  };
  // Relative target is tightened below double-precision needs of the scheduler.
  const double estimate = (right - left) * integrand(0.5 * (left + right));
  const double tolerance = std::max(std::abs(estimate) * 1.0e-11, 1.0e-6);
  return sign * adaptiveSimpson(integrand, left, right, tolerance);
}

double LambdaCdmBackground::cosmicTimeSinceScaleFactorSi(
    double scale_factor,
    double lower_scale_factor) const {
  if (scale_factor < lower_scale_factor) {
    throw std::invalid_argument("scale_factor must be >= lower_scale_factor");
  }
  return cosmicTimeIntervalSi(lower_scale_factor, scale_factor);
}

}  // namespace cosmosim::core
