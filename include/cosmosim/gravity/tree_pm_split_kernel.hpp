#pragma once

#include <array>
#include <cmath>
#include <stdexcept>

#include "cosmosim/gravity/tree_softening.hpp"

namespace cosmosim::gravity {

// Explicitly documented TreePM Gaussian split metadata and utility API.
enum class TreePmSplitKernel {
  kGaussianErfc,
};

struct TreePmSplitPolicy {
  TreePmSplitKernel kernel = TreePmSplitKernel::kGaussianErfc;
  double mesh_spacing_comoving = 0.0;
  double asmth_cells = 0.0;
  double rcut_cells = 0.0;
  double split_scale_comoving = 0.0;
  double cutoff_radius_comoving = 0.0;
};

inline void validateTreePmSplitPolicy(const TreePmSplitPolicy& policy) {
  if (!std::isfinite(policy.mesh_spacing_comoving) || policy.mesh_spacing_comoving <= 0.0) {
    throw std::invalid_argument("TreePM mesh_spacing_comoving must be > 0");
  }
  if (!std::isfinite(policy.asmth_cells) || policy.asmth_cells <= 0.0) {
    throw std::invalid_argument("TreePM asmth_cells must be > 0");
  }
  if (!std::isfinite(policy.rcut_cells) || policy.rcut_cells <= 0.0) {
    throw std::invalid_argument("TreePM rcut_cells must be > 0");
  }
  if (!std::isfinite(policy.split_scale_comoving) || policy.split_scale_comoving <= 0.0) {
    throw std::invalid_argument("TreePM split_scale_comoving must be > 0");
  }
  if (!std::isfinite(policy.cutoff_radius_comoving) || policy.cutoff_radius_comoving <= 0.0) {
    throw std::invalid_argument("TreePM cutoff_radius_comoving must be > 0");
  }
  if (policy.kernel != TreePmSplitKernel::kGaussianErfc) {
    throw std::invalid_argument("Unsupported TreePM split kernel");
  }
}

[[nodiscard]] inline TreePmSplitPolicy makeTreePmSplitPolicyFromMeshSpacing(
    double asmth_cells,
    double rcut_cells,
    double mesh_spacing_comoving,
    TreePmSplitKernel kernel = TreePmSplitKernel::kGaussianErfc) {
  TreePmSplitPolicy policy;
  policy.kernel = kernel;
  policy.mesh_spacing_comoving = mesh_spacing_comoving;
  policy.asmth_cells = asmth_cells;
  policy.rcut_cells = rcut_cells;
  policy.split_scale_comoving = asmth_cells * mesh_spacing_comoving;
  policy.cutoff_radius_comoving = rcut_cells * mesh_spacing_comoving;
  validateTreePmSplitPolicy(policy);
  return policy;
}

[[nodiscard]] inline double treePmGaussianShortRangeForceFactorUnchecked(
    double distance_comoving,
    double split_scale_comoving) noexcept {
  if (distance_comoving <= 0.0) {
    return 1.0;
  }
  const double q = (0.5 * distance_comoving) / split_scale_comoving;
  constexpr double inv_sqrt_pi = 0.564189583547756286948079451560772586;
  if (std::isinf(q)) return 0.0;
  return std::erfc(q) + 2.0 * inv_sqrt_pi * q * std::exp(-q * q);
}

[[nodiscard]] inline double treePmGaussianShortRangeForceFactor(
    double distance_comoving,
    double split_scale_comoving) {
  if (!std::isfinite(distance_comoving) || distance_comoving < 0.0) {
    throw std::invalid_argument("TreePM distance_comoving must be finite and non-negative");
  }
  if (!std::isfinite(split_scale_comoving) || split_scale_comoving <= 0.0) {
    throw std::invalid_argument("TreePM split_scale_comoving must be finite and positive");
  }
  return treePmGaussianShortRangeForceFactorUnchecked(distance_comoving, split_scale_comoving);
}

[[nodiscard]] inline double treePmGaussianLongRangeForceFactorUnchecked(
    double distance_comoving,
    double split_scale_comoving) noexcept {
  if (distance_comoving <= 0.0) return 0.0;
  const double q = (0.5 * distance_comoving) / split_scale_comoving;
  if (std::isinf(q)) return 1.0;
  constexpr double k_inv_sqrt_pi = 0.564189583547756286948079451560772586;
  if (q < 0.125) {
    // L(q)/q^3 = 4/sqrt(pi) sum_n (-q^2)^n/[n! (2n+3)].
    // Through n=7: omitted term < 2e-20 on this interval. Do not form 1-S.
    const double t = q * q;
    const double h = 4.0 * k_inv_sqrt_pi * (1.0/3.0 + t * (-1.0/5.0 +
        t * (1.0/14.0 + t * (-1.0/54.0 + t * (1.0/264.0 +
        t * (-1.0/1560.0 + t * (1.0/10800.0 - t/85680.0)))))));
    return q * q * q * h;
  }
  return std::erf(q) - 2.0 * k_inv_sqrt_pi * q * std::exp(-q * q);
}

// Dimensionless H(t)=L(sqrt(t))/t^(3/2), t=(r/(2a))^2. The same
// polynomial also supplies a stable coefficient when L itself underflows.
[[nodiscard]] inline double treePmGaussianCoefficientSeries(double t) noexcept {
  constexpr double k_four_inv_sqrt_pi = 2.25675833419102514779231780624309034;
  return k_four_inv_sqrt_pi * (1.0/3.0 + t * (-1.0/5.0 + t * (1.0/14.0 +
      t * (-1.0/54.0 + t * (1.0/264.0 + t * (-1.0/1560.0 +
      t * (1.0/10800.0 - t/85680.0)))))));
}

[[nodiscard]] inline double treePmGaussianLongRangeInvR3Unchecked(
    double distance_comoving, double split_scale_comoving) noexcept {
  const double q = (0.5 * distance_comoving) / split_scale_comoving;
  if (q < 0.125) {
    // Sequential divisions avoid premature a^3 overflow/underflow.
    return ((0.125 * treePmGaussianCoefficientSeries(q*q) /
        split_scale_comoving) / split_scale_comoving) / split_scale_comoving;
  }
  return ((treePmGaussianLongRangeForceFactorUnchecked(distance_comoving,
      split_scale_comoving) / distance_comoving) / distance_comoving) / distance_comoving;
}

// Experimental direct-pair coefficient only. One immutable bounded table per enabled coordinator,
// reserved before allocation and initialized before worker launch.
// Cubic Hermite on t in [0,64], 4096 cells. Qualification target: absolute
// dimensionless H error <= 5e-11 (NOT a relative residual-force guarantee).
// Outside the table use the analytic law. Split changes require no invalidation:
// t and the 1/(8a^3) normalization are supplied by each call.
class TreePmGaussianCoefficientTable {
 public:
  static constexpr std::size_t k_cells = 4096U;
  static constexpr double k_step = 64.0 / static_cast<double>(k_cells);
  [[nodiscard]] double longRangeInvR3(double r2, double a) const noexcept {
    const double t = 0.25 * (r2 / a) / a;
    if (!(t >= 0.0) || t >= 64.0) {
      return treePmGaussianLongRangeInvR3Unchecked(std::sqrt(r2), a);
    }
    return ((0.125*dimensionlessCoefficient(t)/a)/a)/a;
  }
  [[nodiscard]] double dimensionlessCoefficient(double t) const noexcept {
    if (t < 0.015625) return treePmGaussianCoefficientSeries(t);
    if (!(t < 64.0)) {
      const double q = std::sqrt(t);
      return treePmGaussianLongRangeInvR3Unchecked(q, 0.5);
    }
    const double cell = t / k_step;
    const auto i = static_cast<std::size_t>(cell);
    const double u = cell - static_cast<double>(i);
    const double u2 = u*u;
    const double u3 = u2*u;
    return (2*u3 - 3*u2 + 1)*m_value[i] + (u3 - 2*u2 + u)*m_slope[i] +
        (-2*u3 + 3*u2)*m_value[i+1U] + (u3 - u2)*m_slope[i+1U];
  }
  TreePmGaussianCoefficientTable() noexcept {
    constexpr double k_two_inv_sqrt_pi = 1.12837916709551257389615890312154517;
    for (std::size_t i = 0U; i <= k_cells; ++i) {
      const double t = static_cast<double>(i)*k_step;
      const double q = std::sqrt(t);
      m_value[i] = i == 0U ? treePmGaussianCoefficientSeries(0.0) :
          treePmGaussianLongRangeForceFactorUnchecked(q, 0.5)/(t*q);
      const double derivative = i == 0U ? -4.0*0.564189583547756286948079451560772586/5.0 :
          (k_two_inv_sqrt_pi*std::exp(-t) - 1.5*m_value[i])/t;
      m_slope[i] = k_step*derivative;
    }
  }
 private:
  std::array<double, k_cells+1U> m_value{};
  std::array<double, k_cells+1U> m_slope{};
};

[[nodiscard]] inline double treePmGaussianLongRangeForceFactor(
    double distance_comoving,
    double split_scale_comoving) {
  if (!std::isfinite(distance_comoving) || distance_comoving < 0.0) {
    throw std::invalid_argument("TreePM distance_comoving must be finite and non-negative");
  }
  if (!std::isfinite(split_scale_comoving) || split_scale_comoving <= 0.0) {
    throw std::invalid_argument("TreePM split_scale_comoving must be finite and positive");
  }
  return treePmGaussianLongRangeForceFactorUnchecked(distance_comoving, split_scale_comoving);
}

// PM carries the unsoftened Gaussian long-range field.  Therefore a finite-
// softening TreePM solve must use a real-space residual equal to
//
//   K_short = K_softened - K_long,newtonian
//
// rather than K_softened multiplied by the Newtonian short-range factor.
// This definition composes Tree + PM to the requested Plummer-softened force
// before the explicit short-range cutoff is applied.
[[nodiscard]] inline double treePmSoftenedShortRangeInvR3Unchecked(
    double squared_distance,
    double epsilon_comoving,
    double split_scale_comoving) noexcept {
  if (squared_distance <= 0.0) {
    return softenedInvR3Unchecked(squared_distance, epsilon_comoving);
  }
  const double distance = std::sqrt(squared_distance);
  return softenedInvR3Unchecked(squared_distance, epsilon_comoving) -
      treePmGaussianLongRangeInvR3Unchecked(distance, split_scale_comoving);
}

[[nodiscard]] inline double treePmSoftenedShortRangeInvR3(
    double squared_distance,
    double epsilon_comoving,
    double split_scale_comoving) {
  if (!std::isfinite(squared_distance) || squared_distance < 0.0) {
    throw std::invalid_argument("TreePM squared_distance must be finite and non-negative");
  }
  validatedSofteningEpsilon(epsilon_comoving, "TreePM pair softening");
  if (!std::isfinite(split_scale_comoving) || split_scale_comoving <= 0.0) {
    throw std::invalid_argument("TreePM split_scale_comoving must be finite and positive");
  }
  return treePmSoftenedShortRangeInvR3Unchecked(
      squared_distance, epsilon_comoving, split_scale_comoving);
}

[[nodiscard]] inline double treePmGaussianFourierLongRangeFilterUnchecked(
    double wave_number_comoving,
    double split_scale_comoving) noexcept {
  const double product = wave_number_comoving * split_scale_comoving;
  return std::exp(-(product * product));
}

[[nodiscard]] inline double treePmGaussianFourierLongRangeFilter(
    double wave_number_comoving,
    double split_scale_comoving) {
  if (!std::isfinite(wave_number_comoving) || wave_number_comoving < 0.0) {
    throw std::invalid_argument("TreePM wave_number_comoving must be finite and non-negative");
  }
  if (!std::isfinite(split_scale_comoving) || split_scale_comoving <= 0.0) {
    throw std::invalid_argument("TreePM split_scale_comoving must be finite and positive");
  }
  return treePmGaussianFourierLongRangeFilterUnchecked(wave_number_comoving, split_scale_comoving);
}

}  // namespace cosmosim::gravity
