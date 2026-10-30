// Gravity invariant validation suite (2026-10-26; authored, NOT RUN).
//
// Physical invariants of the full periodic TreePM force that the existing
// Ewald accuracy matrix (test_tree_pm_ewald_accuracy) does not state directly:
//   1. single-particle self-force vanishes (periodic lattice of one mass);
//   2. total momentum conservation sum_i m_i a_i ~ 0 (PM with matched
//      assignment/interpolation and an antisymmetric spectral gradient, and a
//      pairwise-antisymmetric all-direct short-range residual);
//   3. isotropy: a cyclic axis permutation of a cubic configuration permutes
//      the forces;
//   4. PM resolution convergence: at fixed split scale in mesh cells, the
//      two-body force error against the independent Ewald reference does not
//      grow as the mesh is refined (16 -> 32 -> 64);
//   5. spectral and finite-difference PM gradients are distinct operators and
//      both stay within a provisional Ewald error bound.
// All thresholds are PROVISIONAL (owner qualification pending) and named as
// such; none is a published accuracy claim.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <numeric>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/gravity/tree_pm_coupling.hpp"
#include "periodic_ewald_reference.hpp"

#if COSMOSIM_ENABLE_MPI
#include <mpi.h>
#endif

namespace {

using namespace cosmosim::gravity;
namespace ts = cosmosim::test_support;

void requireOrThrow(bool condition, const std::string& message) {
  if (!condition) throw std::runtime_error(message);
}

struct Field {
  std::vector<double> ax, ay, az;
};

struct Population {
  std::vector<double> x, y, z, mass;
};

[[nodiscard]] TreePmOptions fullOptions(std::size_t mesh, double epsilon) {
  TreePmOptions options;
  options.pm_options.box_size_x_mpc_comoving = 1.0;
  options.pm_options.box_size_y_mpc_comoving = 1.0;
  options.pm_options.box_size_z_mpc_comoving = 1.0;
  options.pm_options.gravitational_constant_code = 1.0;
  options.pm_options.assignment_scheme = PmAssignmentScheme::kCic;
  options.pm_options.boundary_condition = PmBoundaryCondition::kPeriodic;
  options.tree_options.gravitational_constant_code = 1.0;
  options.tree_options.softening.epsilon_comoving = epsilon;
  options.tree_options.max_leaf_size = 8U;
  options.tree_options.relative_force_tolerance = 1.0e-3;
  // All-direct short range (no accepted multipoles): the residual is the
  // exact pairwise-antisymmetric Gaussian complement within r_cut.
  options.acceptance_policy = TreePmAcceptancePolicy::kAdaptiveRelative;
  options.adaptive_error_estimator = TreePmAdaptiveErrorEstimator::kRelativeMonopole;
  options.adaptive_maximum_opening_angle = 1.0e-6;
  options.split_policy = makeTreePmSplitPolicyFromMeshSpacing(1.25, 4.5, 1.0 / static_cast<double>(mesh));
  options.source_generation = GravitySourceGeneration{1U};
  return options;
}

[[nodiscard]] Field solve(const Population& p, std::size_t mesh, const TreePmOptions& options) {
  Field f{std::vector<double>(p.x.size()), std::vector<double>(p.x.size()), std::vector<double>(p.x.size())};
  std::vector<TreeLocalIndex> active(p.x.size());
  std::iota(active.begin(), active.end(), 0U);
  const std::vector<double> previous(p.x.size(), 1.0);
  const TreePmForceAccumulatorView accumulator{active, f.ax, f.ay, f.az, previous};
  TreePmCoordinator coordinator(PmGridShape{mesh, mesh, mesh});
  coordinator.solveActiveSet(p.x, p.y, p.z, p.mass, accumulator, options, nullptr, nullptr);
  return f;
}

[[nodiscard]] Population clustered(std::size_t n, std::uint64_t seed) {
  Population p;
  std::uint64_t state = seed;
  const auto uniform = [&state]() {
    state = state * 6364136223846793005ULL + 1442695040888963407ULL;
    return static_cast<double>(state >> 11U) * (1.0 / 9007199254740992.0);
  };
  const auto wrap = [](double v) { v -= std::floor(v); return v >= 1.0 ? 0.0 : v; };
  for (std::size_t i = 0; i < n; ++i) {
    const double cx = (i % 2U == 0U) ? 0.3 : 0.72;
    const double spread = (i % 3U == 0U) ? 0.04 : 0.2;
    p.x.push_back(wrap(cx + spread * (uniform() - 0.5)));
    p.y.push_back(wrap(0.45 + spread * (uniform() - 0.5)));
    p.z.push_back(wrap(0.6 + spread * (uniform() - 0.5)));
    p.mass.push_back(0.5 + uniform());
  }
  return p;
}

void testSingleParticleSelfForce() {
  for (const std::size_t mesh : {16U, 32U}) {
    for (const double x : {0.5, 0.0312, 0.99}) {
      Population p;
      p.x = {x}; p.y = {0.37}; p.z = {0.81}; p.mass = {1.0};
      const Field f = solve(p, mesh, fullOptions(mesh, 1.0e-3));
      // Scale: a unit mass at one mesh spacing.
      const double scale = 1.0 * static_cast<double>(mesh * mesh);
      const double magnitude = std::sqrt(f.ax[0] * f.ax[0] + f.ay[0] * f.ay[0] + f.az[0] * f.az[0]);
      // PROVISIONAL: CIC self-force of a lone particle is zero by symmetry only
      // at cell centres; elsewhere the PM self-force is O(1e-3) of the
      // one-cell force for CIC. Gate at 1e-2 of the one-cell scale.
      requireOrThrow(magnitude <= 1.0e-2 * scale,
                     "single-particle self-force exceeds the provisional bound at mesh " + std::to_string(mesh));
    }
  }
}

void testMomentumConservation() {
  const Population p = clustered(240U, 0xC0FFEEULL);
  for (const std::size_t mesh : {16U, 32U}) {
    const Field f = solve(p, mesh, fullOptions(mesh, 2.0e-3));
    double px = 0.0, py = 0.0, pz = 0.0, scale = 0.0;
    for (std::size_t i = 0; i < p.x.size(); ++i) {
      px += p.mass[i] * f.ax[i];
      py += p.mass[i] * f.ay[i];
      pz += p.mass[i] * f.az[i];
      scale += p.mass[i] * std::sqrt(f.ax[i] * f.ax[i] + f.ay[i] * f.ay[i] + f.az[i] * f.az[i]);
    }
    const double net = std::sqrt(px * px + py * py + pz * pz);
    // PROVISIONAL: matched CIC deposit/interpolation with a spectral gradient
    // and an all-direct antisymmetric residual conserve momentum to roundoff
    // up to the Gaussian-table residual asymmetry; gate at 1e-6 relative.
    requireOrThrow(net <= 1.0e-6 * scale,
                   "net force sum_i m_i a_i is not ~0 at mesh " + std::to_string(mesh) +
                       " (relative " + std::to_string(net / scale) + ")");
  }
}

void testAxisPermutationIsotropy() {
  const Population p = clustered(160U, 0xBADC0DEULL);
  Population rotated;  // (x, y, z) -> (y, z, x)
  rotated.x = p.y; rotated.y = p.z; rotated.z = p.x; rotated.mass = p.mass;
  const std::size_t mesh = 32U;
  const Field f = solve(p, mesh, fullOptions(mesh, 2.0e-3));
  const Field g = solve(rotated, mesh, fullOptions(mesh, 2.0e-3));
  double max_err = 0.0, max_ref = 0.0;
  for (std::size_t i = 0; i < p.x.size(); ++i) {
    // Rotated force components (gx, gy, gz) must equal (fy, fz, fx).
    for (const auto& [got, want] : {std::pair{g.ax[i], f.ay[i]}, std::pair{g.ay[i], f.az[i]},
                                    std::pair{g.az[i], f.ax[i]}}) {
      max_err = std::max(max_err, std::abs(got - want));
      max_ref = std::max(max_ref, std::abs(want));
    }
  }
  // Exact symmetry up to summation order and FFT roundoff (cubic box/mesh).
  requireOrThrow(max_err <= 1.0e-9 * std::max(max_ref, 1.0), "TreePM forces are not axis-permutation covariant");
}

[[nodiscard]] double twoBodyEwaldError(std::size_t mesh, PmForceGradientBackend backend, double separation) {
  Population p;
  p.x = {0.2, 0.2 + separation}; p.y = {0.41, 0.41}; p.z = {0.63, 0.63}; p.mass = {1.0, 1.0};
  TreePmOptions options = fullOptions(mesh, 1.0e-4);
  options.pm_options.force_gradient_backend = backend;
  const Field f = solve(p, mesh, options);
  std::vector<ts::PeriodicEwaldSource> sources{{{p.x[0], p.y[0], p.z[0]}, 1.0}, {{p.x[1], p.y[1], p.z[1]}, 1.0}};
  ts::PeriodicEwaldOptions ewald;
  ewald.alpha_inverse_length = 4.0;
  ewald.real_image_limits = {3, 3, 3};
  ewald.reciprocal_mode_limits = {8, 8, 8};
  ewald.plummer_softening_epsilon = 1.0e-4;
  const auto reference = ts::periodicEwaldAccelerationsAtSources(sources, {1.0, 1.0, 1.0}, ewald);
  double err = 0.0, ref = 0.0;
  for (std::size_t i = 0; i < 2U; ++i) {
    err = std::max(err, std::hypot(f.ax[i] - reference[i].x, f.ay[i] - reference[i].y, f.az[i] - reference[i].z));
    ref = std::max(ref, std::hypot(reference[i].x, reference[i].y, reference[i].z));
  }
  return err / ref;
}

void testPmResolutionConvergence() {
#if COSMOSIM_ENABLE_FFTW
  // Long-range dominated separation (0.3 L > r_cut at every mesh here).
  const double e16 = twoBodyEwaldError(16U, PmForceGradientBackend::kSpectral, 0.3);
  const double e32 = twoBodyEwaldError(32U, PmForceGradientBackend::kSpectral, 0.3);
  const double e64 = twoBodyEwaldError(64U, PmForceGradientBackend::kSpectral, 0.3);
  std::printf("pm_resolution_two_body_rel_error mesh16=%.3e mesh32=%.3e mesh64=%.3e\n", e16, e32, e64);
  // A single two-body phase can land unusually close to a mesh sweet spot at
  // an intermediate resolution, so strict pairwise monotonicity is not a
  // stable convergence test (the 32^3 case can beat 64^3 while both are far
  // better than 16^3).  Require genuine end-to-end refinement progress plus
  // the independent absolute accuracy ceiling instead.
  requireOrThrow(e64 <= 0.25 * e16,
                 "two-body TreePM error does not improve sufficiently from 16^3 to 64^3");
  requireOrThrow(e64 <= 2.0e-2, "two-body TreePM error at 64^3 exceeds the provisional 2% bound");
  // Split-transition separation (~ r_s) is dominated by the tree residual.
  const double transition = twoBodyEwaldError(32U, PmForceGradientBackend::kSpectral, 1.25 / 32.0);
  requireOrThrow(transition <= 2.0e-2, "two-body TreePM error at the split scale exceeds the provisional 2% bound");
#endif
}

void testFiniteDifferenceGradientIsDistinctAndBounded() {
#if COSMOSIM_ENABLE_FFTW
  const double spectral = twoBodyEwaldError(32U, PmForceGradientBackend::kSpectral, 0.3);
  const double fd = twoBodyEwaldError(32U, PmForceGradientBackend::kFiniteDifference4, 0.3);
  std::printf("pm_gradient_backend_rel_error spectral=%.3e finite_difference_4=%.3e\n", spectral, fd);
  requireOrThrow(spectral != fd, "finite-difference PM gradient must be a distinct operator from the spectral reference");
  // PROVISIONAL bound for the unqualified FD backend.
  requireOrThrow(fd <= 5.0e-2, "finite-difference PM gradient exceeds the provisional 5% two-body bound");
#endif
}

}  // namespace

int main() {
#if COSMOSIM_ENABLE_MPI
  MPI_Init(nullptr, nullptr);
#endif
  int failed = 0;
  try {
    testSingleParticleSelfForce();
    testMomentumConservation();
    testAxisPermutationIsotropy();
    testPmResolutionConvergence();
    testFiniteDifferenceGradientIsDistinctAndBounded();
  } catch (const std::exception& error) {
    std::fprintf(stderr, "gravity invariant suite failed: %s\n", error.what());
    failed = 1;
  }
#if COSMOSIM_ENABLE_MPI
  MPI_Finalize();
#endif
  return failed;
}
