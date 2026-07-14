#pragma once

#include <array>
#include <cstddef>
#include <limits>
#include <span>
#include <vector>

namespace cosmosim::test_support {

struct PeriodicEwaldVector3 {
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;
};

struct PeriodicEwaldBox {
  double length_x = 0.0;
  double length_y = 0.0;
  double length_z = 0.0;
};

struct PeriodicEwaldSource {
  PeriodicEwaldVector3 position{};
  double mass = 0.0;
};

inline constexpr std::size_t k_no_periodic_ewald_self_source =
    std::numeric_limits<std::size_t>::max();

struct PeriodicEwaldTarget {
  PeriodicEwaldVector3 position{};
  std::size_t self_source_index = k_no_periodic_ewald_self_source;
};

struct PeriodicEwaldOptions {
  // G in the caller's units. The returned particle kernel is scale-free in
  // comoving coordinates; cosmological kick factors belong to the integrator.
  double gravitational_constant = 1.0;

  // Ewald alpha, with units inverse to the position/box-length unit.
  double alpha_inverse_length = 1.0;

  // Inclusive image ranges [-limit_i, +limit_i] on real-space {x, y, z}.
  std::array<int, 3> real_image_limits{1, 1, 1};

  // Inclusive integer-mode ranges [-limit_i, +limit_i] on reciprocal
  // {x, y, z}. The (0, 0, 0) mode is always omitted.
  std::array<int, 3> reciprocal_mode_limits{4, 4, 4};
};

// Test-only, unsoftened, double-precision periodic acceleration reference.
// It does not call or share kernels with the production PM/TreePM solvers.
//
// With d_ij = x_j - x_i reduced to a minimum-image representative,
// R_n = d_ij + (n_x L_x, n_y L_y, n_z L_z), and
// k_h = 2*pi*(h_x/L_x, h_y/L_y, h_z/L_z), this evaluates
//
//   a_i(real) = G sum_j m_j sum_n' R_n/|R_n|^3
//       * [erfc(alpha |R_n|)
//          + 2 alpha |R_n|/sqrt(pi) exp(-alpha^2 |R_n|^2)]
//
//   a_i(recip) = 4*pi*G/V sum_j m_j sum_{h != 0}
//       exp(-|k_h|^2/(4 alpha^2)) k_h/|k_h|^2 sin(k_h . d_ij).
//
// The omitted reciprocal zero mode is the homogeneous compensating-background
// convention for the periodic Poisson equation sourced by rho - mean(rho).
// In the corresponding periodic Green function its Ewald term is
// -pi/(alpha^2 V); the gravitational potential therefore receives the constant
// +pi G M/(alpha^2 V). Its gradient is zero, so it has no acceleration term
// here (and the potential gauge remains outside this acceleration-only API).
//
// If self_source_index names source j, the singular n=(0,0,0) real term is
// omitted and the analytically zero reciprocal self force is skipped. Periodic
// self images remain in the symmetric real sum. Distinct, coincident unsoftened
// source-target pairs are rejected rather than silently regularized.
// The direct O(N_target N_source (N_image + N_mode)) loops and compensated
// sums are deliberate: this is an auditable small-N reference, not a runtime
// solver or benchmark implementation.
[[nodiscard]] std::vector<PeriodicEwaldVector3> periodicEwaldAccelerations(
    std::span<const PeriodicEwaldSource> sources,
    std::span<const PeriodicEwaldTarget> targets,
    const PeriodicEwaldBox& box,
    const PeriodicEwaldOptions& options);

// Convenience entry point for acceleration on every source. Target i is marked
// as the self image of source i, so the central self interaction is excluded.
[[nodiscard]] std::vector<PeriodicEwaldVector3> periodicEwaldAccelerationsAtSources(
    std::span<const PeriodicEwaldSource> sources,
    const PeriodicEwaldBox& box,
    const PeriodicEwaldOptions& options);

}  // namespace cosmosim::test_support
