#pragma once

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <limits>
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
  // Hot-path preconditions: finite r2 >= 0, finite a > 0, validated by the
  // coordinator. Like the Unchecked analytic helpers, this does not throw.
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
    // Preserve the established zero-softening/coincident convention (zero).
    // For finite positive Plummer softening the radial coefficient is the
    // continuous softened-minus-Gaussian limit, even though its vector is zero.
    if (epsilon_comoving == 0.0) return softenedInvR3Unchecked(squared_distance, epsilon_comoving);
    // Scale before subtracting: two individually overflowing coefficients
    // must not produce inf-inf/NaN. The dimensionless ratio is always <=1.
    constexpr double k_inverse_six_sqrt_pi = 0.094031597257959381158013241926795431;
    if (epsilon_comoving <= split_scale_comoving) {
      const double ratio = epsilon_comoving / split_scale_comoving;
      const double coefficient = 1.0 - k_inverse_six_sqrt_pi * ratio * ratio * ratio;
      return ((coefficient / epsilon_comoving) / epsilon_comoving) / epsilon_comoving;
    }
    const double ratio = split_scale_comoving / epsilon_comoving;
    const double coefficient = ratio * ratio * ratio - k_inverse_six_sqrt_pi;
    return ((coefficient / split_scale_comoving) / split_scale_comoving) / split_scale_comoving;
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

// Derivatives of f(x), x=r^2, for acceleration d*f(x). The contractions are
// f'(r)/r=2 df/dx and (f''(r)-f'(r)/r)/r^2=4 d^2f/dx^2. One exponential and
// one erf supply all three coefficients; the small-t series avoids cancellation.
struct TreePmResidualRadialCoefficients {
  double value;
  double first_over_r;
  double contracted_second;
};

[[nodiscard]] inline TreePmResidualRadialCoefficients treePmResidualRadialCoefficients(
    double x, double epsilon, double a) noexcept {
  const double d = x + epsilon*epsilon;
  const double softened = 1.0/(d*std::sqrt(d));
  const double t = (0.25*x/a)/a;
  constexpr double c = 2.25675833419102514779231780624309034;
  double h, dh, ddh;
  if (t < 0.015625) {
    // Simultaneous Horner derivatives of the established n<=7 series.
    constexpr std::array<double, 8> coefficient{
        1.0/3.0, -1.0/5.0, 1.0/14.0, -1.0/54.0,
        1.0/264.0, -1.0/1560.0, 1.0/10800.0, -1.0/85680.0};
    h = coefficient.back(); dh = 0.0; ddh = 0.0;
    for (std::size_t i = 7U; i > 0U; --i) {
      ddh = ddh*t + 2.0*dh;
      dh = dh*t + h;
      h = h*t + coefficient[i-1U];
    }
    h *= c; dh *= c; ddh *= c;
  } else {
    const double q = std::sqrt(t);
    const double exponential = std::exp(-t);
    h = (std::erf(q) - 0.5*c*q*exponential)/(t*q);
    dh = (0.5*c*exponential - 1.5*h)/t;
    ddh = (-0.5*c*exponential - 2.5*dh)/t;
  }
  const double scale = ((0.125/a)/a)/a;
  const double dt_dx = (0.25/a)/a;
  return {softened - scale*h,
      -3.0*softened/d - 2.0*scale*dt_dx*dh,
      15.0*softened/(d*d) - 4.0*scale*dt_dx*dt_dx*ddh};
}

// Optional homogeneous-pair accelerator. Geometric bins in s=r^2+eps^2,
// indexed by exponent/mantissa bits, resolve the Plummer core without log or
// sqrt in eligible pair evaluations. Each cell stores cubic Hermite data and
// an absolute fourth-derivative error envelope. Ill-conditioned cancellation,
// zero/underflow softening and unsupported dynamic ranges use the analytic law.
// Fixed storage, independent of particle/node/thread count. No per-pair allocation.
class TreePmResidualCoefficientTable {
 public:
  static constexpr std::size_t k_cells = 16384U;
  static constexpr unsigned k_shift = 43U; // 512 mantissa bins per exponent
  static constexpr double k_relative_error_limit = 5.0e-11;

  [[nodiscard]] static bool eligible(double epsilon, double a, double cutoff) noexcept {
    const double e2 = epsilon*epsilon;
    const double end = cutoff*cutoff + e2;
    if (!std::numeric_limits<double>::is_iec559 || sizeof(double) != sizeof(std::uint64_t) ||
        !std::isnormal(e2) || !std::isnormal(end) || !std::isnormal(a*a) ||
        !std::isnormal(((0.125/a)/a)/a) || end <= e2) return false;
    const auto begin_code = std::bit_cast<std::uint64_t>(e2) >> k_shift;
    const auto end_code = std::bit_cast<std::uint64_t>(end) >> k_shift;
    return end_code >= begin_code && end_code - begin_code + 1U <= k_cells;
  }
  TreePmResidualCoefficientTable(double epsilon, double a, double cutoff) noexcept
      : m_epsilon(epsilon), m_a(a), m_cutoff(cutoff), m_epsilon2(epsilon*epsilon) {
    if (!eligible(epsilon, a, cutoff)) return;
    // Start at the first knot >= eps^2; the tiny omitted core interval falls back.
    const auto bits = std::bit_cast<std::uint64_t>(m_epsilon2);
    m_begin = (bits >> k_shift) + ((bits & ((1ULL << k_shift)-1U)) != 0U);
    const auto end = std::bit_cast<std::uint64_t>(cutoff*cutoff+m_epsilon2) >> k_shift;
    m_count = end >= m_begin ? static_cast<std::size_t>(end-m_begin+1U) : 0U;
    for (std::size_t i = 0; i < m_count; ++i) {
      const double left = knot(m_begin+i), right = knot(m_begin+i+1U);
      const double x = left-m_epsilon2, h = right-left;
      const auto l = treePmResidualRadialCoefficients(x, epsilon, a);
      const auto r = treePmResidualRadialCoefficients(right-m_epsilon2, epsilon, a);
      auto& cell = m_cell[i];
      cell = {l.value, r.value, 0.5*h*l.first_over_r, 0.5*h*r.first_over_r, 0.0};
      // Hermite error <= h^4 max|f''''(x)|/384. Plummer derivative is
      // (945/16)s^(-11/2). Gaussian integral gives a global bound and,
      // for x>0, the tighter unscreened Newtonian derivative bound.
      const double hs = h/left;
      const double plummer = (945.0/16.0)/384.0/(left*std::sqrt(left))*hs*hs*hs*hs;
      const double ha = (0.25*h/a)/a;
      double gaussian = ((0.5*0.564189583547756286948079451560772586/a)/a)/a /
          (11.0*384.0)*ha*ha*ha*ha;
      if (x > 0.0) {
        const double hx = h/x;
        const double newton = (945.0/16.0)/384.0/(x*std::sqrt(x))*hx*hx*hx*hx;
        gaussian = std::min(gaussian, newton);
      }
      cell.error = 2.0*(plummer+gaussian) + 128.0*std::numeric_limits<double>::epsilon()*
          (std::abs(l.value)+std::abs(r.value)+1.0/(left*std::sqrt(left)));
    }
  }
  [[nodiscard]] bool matches(double epsilon, double a, double cutoff) const noexcept {
    return epsilon == m_epsilon && a == m_a && cutoff == m_cutoff;
  }
  // False means the caller must use its analytic reference, not another table.
  [[nodiscard]] bool evaluate(double x, double& value) const noexcept {
    if (!(x > 0.0) || !std::isfinite(x)) return false;
    const double s = x+m_epsilon2;
    const auto code = std::bit_cast<std::uint64_t>(s) >> k_shift;
    if (code < m_begin || code-m_begin >= m_count) return false;
    const auto& cell = m_cell[static_cast<std::size_t>(code-m_begin)];
    const double left = knot(code), h = knot(code+1U)-left;
    const double u = (s-left)/h, u2=u*u, u3=u2*u;
    value = (2*u3-3*u2+1)*cell.left + (u3-2*u2+u)*cell.left_slope +
        (-2*u3+3*u2)*cell.right + (u3-u2)*cell.right_slope;
    return std::isfinite(value) && std::isfinite(cell.error) &&
        cell.error <= k_relative_error_limit*std::abs(value);
  }
 private:
  struct Cell { double left, right, left_slope, right_slope, error; };
  [[nodiscard]] static double knot(std::uint64_t code) noexcept {
    return std::bit_cast<double>(code << k_shift);
  }
  double m_epsilon, m_a, m_cutoff, m_epsilon2;
  std::uint64_t m_begin = 0U;
  std::size_t m_count = 0U;
  std::array<Cell, k_cells> m_cell{};
};

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
