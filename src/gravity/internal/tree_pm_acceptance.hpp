#pragma once

#include <algorithm>
#include <cmath>
#include <limits>

#include "cosmosim/gravity/tree_pm_coupling.hpp"
#include "tree_interaction_common.hpp"

namespace cosmosim::gravity::internal {

// Accuracy proxy for CHUI's actual residual F_i = G M x_i f(r), where
// f = (r^2+eps^2)^(-3/2) + (S(r)-1)/r^3. The quadrupole path retains the
// complete second moment (including its trace); this proxy deliberately
// charges the second-order scale even there. It is NOT a certified bound on
// the quadrupole remainder. Qualification must calibrate the requested alpha.
// E = G M max(l^2/r^4, rho^2 [3 |f'| + |r f''-f'|]/2),
// rho = sqrt(3) h + |COM-center|. Accept E <= alpha |A_previous|.
// A is the unscaled TOTAL Tree+PM kernel, not A/a^2 or A/a^3.
[[nodiscard]] inline double residualAccuracyProxy(
    const TreeNodeAcceptanceInput& in, double split_scale, double g_code) noexcept {
  if (!(in.r2 > 0.0) || !std::isfinite(in.r2)) {
    return std::numeric_limits<double>::infinity();
  }
  const double r = std::sqrt(in.r2);
  const double eps = combineSofteningPairEpsilonUnchecked(
      in.node_softening_max_comoving, in.target_softening_comoving);
  const double d = in.r2 + eps * eps;
  const double inv_d5 = 1.0 / (d * d * std::sqrt(d));
  const double inv_r3 = 1.0 / (in.r2 * r);
  const double inv_r4 = 1.0 / (in.r2 * in.r2);
  const double s = treePmGaussianShortRangeForceFactorUnchecked(r, split_scale);
  const double q = r / (2.0 * split_scale);
  const double exponential = std::exp(-q * q);
  constexpr double k_inv_sqrt_pi = 0.564189583547756286948079451560772586;
  const double scale3 = split_scale * split_scale * split_scale;
  const double s_first = -0.5 * k_inv_sqrt_pi * in.r2 * exponential / scale3;
  const double s_second = k_inv_sqrt_pi * exponential *
      (-r / scale3 + 0.25 * r * in.r2 / (scale3 * split_scale * split_scale));
  const double first = -3.0 * r * inv_d5 + s_first * inv_r3 -
      3.0 * (s - 1.0) * inv_r4;
  const double second = -3.0 * inv_d5 + 15.0 * in.r2 * inv_d5 / d +
      s_second * inv_r3 - 6.0 * s_first * inv_r4 +
      12.0 * (s - 1.0) * inv_r4 / r;
  const double width = 2.0 * in.half_size;
  const double rho = std::sqrt(3.0) * in.half_size + in.com_center_offset;
  const double proxy = g_code * in.node_mass_code * std::max(
      width * width * inv_r4,
      0.5 * rho * rho * (3.0 * std::abs(first) + std::abs(r * second - first)));
  return std::isfinite(proxy) && proxy >= 0.0 ? proxy : std::numeric_limits<double>::infinity();
}

[[nodiscard]] inline bool acceptTreePmNode(
    const TreeNodeAcceptanceInput& in, bool within_cutoff,
    const TreePmOptions& options, TreePmTraversalCounters& counters) noexcept {
  if (in.is_leaf) return true;
  const double r = std::sqrt(std::max(in.r2, 1.0e-30));
  const bool softening_ok = passesSofteningEnvelopeGuard(
      false, in.half_size, r, in.target_softening_comoving,
      in.node_softening_min_comoving, in.node_softening_max_comoving);
  bool mac_ok = false;
  bool envelope_ok = false;
  if (options.acceptance_policy == TreePmAcceptancePolicy::kStrictReference) {
    mac_ok = acceptNodeByMac(false, in.target_inside_node, in.half_size,
        in.com_center_offset, in.node_mass_code, in.r2,
        in.previous_acceleration_available, in.previous_acceleration_magnitude_code,
        options.tree_options);
    envelope_ok = options.tree_options.multipole_order == TreeMultipoleOrder::kQuadrupole &&
        (2.0 * in.half_size / r) < 0.08;
    if (!envelope_ok) ++counters.strict_envelope_rejections;
    if (!mac_ok) ++counters.selected_mac_rejections;
  } else {
    // No fabricated acceleration floor. Tiny/zero, missing, negative or
    // nonfinite history uses geometric fallback, deterministic per target.
    const bool history_ok = in.previous_acceleration_available &&
        std::isfinite(in.previous_acceleration_magnitude_code) &&
        in.previous_acceleration_magnitude_code >
            options.tree_options.relative_force_acceleration_floor_code;
    if (history_ok) {
      const double estimate = residualAccuracyProxy(in,
          options.split_policy.split_scale_comoving,
          options.tree_options.gravitational_constant_code);
      const double allowance = options.tree_options.relative_force_tolerance *
          in.previous_acceleration_magnitude_code;
      mac_ok = std::isfinite(estimate) && std::isfinite(allowance) && estimate <= allowance;
      if (!mac_ok) ++counters.relative_mac_rejections;
    } else {
      ++counters.geometric_history_fallbacks;
      mac_ok = acceptNodeByComDistanceMac(options.tree_options.opening_theta,
          in.half_size, in.com_center_offset, in.r2);
      if (!mac_ok) ++counters.selected_mac_rejections;
    }
    envelope_ok = (2.0 * in.half_size / r) < options.adaptive_maximum_opening_angle;
    if (!envelope_ok) ++counters.maximum_angle_rejections;
  }
  if (!softening_ok) ++counters.softening_rejections;
  if (in.target_inside_node) ++counters.near_node_rejections;
  if (!within_cutoff) ++counters.cutoff_containment_rejections;
  // Rejection counters overlap: each failed guard is counted. Opened nodes
  // remain the unique number of descent decisions.
  return mac_ok && envelope_ok && softening_ok && within_cutoff &&
      !in.target_inside_node && std::isfinite(in.r2);
}

}  // namespace cosmosim::gravity::internal
