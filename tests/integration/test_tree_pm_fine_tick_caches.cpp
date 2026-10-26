// Fine-tick cache regression (2026-10-26), single rank.
//
// A retained TreePmCoordinator reused across "fine ticks" (refit-enabled,
// topology-generation caches, certified uniform-mass reuse) must reproduce
// an eager reference — a fresh coordinator per tick with refit and reuse
// disabled — for active fractions from all particles down to one, small
// motions (certified refit), large motions and periodic seam crossings
// (refit rejected → rebuild), and dense-row reorders (layout change →
// rebuild). The configuration opens every node, so both sides evaluate the
// identical direct residual and differ only by summation order (1e-12 rel).
// Also: uniform-mass scalar path versus explicit mass lane, exact mass
// scaling, and OpenMP thread-count determinism of per-target forces.

#include <algorithm>
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
#include "cosmosim/parallel/distributed_memory.hpp"

#if defined(_OPENMP)
#include <omp.h>
#endif
#if COSMOSIM_ENABLE_MPI
#include <mpi.h>
#endif

namespace {

using namespace cosmosim::gravity;

void requireOrThrow(bool condition, const std::string& message) {
  if (!condition) throw std::runtime_error(message);
}

struct Sources {
  std::vector<double> x, y, z, mass;
};

[[nodiscard]] Sources makeSources(std::size_t n, std::uint64_t seed) {
  Sources s;
  std::uint64_t state = seed;
  const auto uniform = [&state]() {
    state = state * 6364136223846793005ULL + 1442695040888963407ULL;
    return static_cast<double>(state >> 11U) * (1.0 / 9007199254740992.0);
  };
  const auto wrap = [](double v) { v -= std::floor(v); return v >= 1.0 ? 0.0 : v; };
  for (std::size_t i = 0; i < n; ++i) {
    const double cx = (i % 3U == 0U) ? 0.002 : (i % 3U == 1U ? 0.5 : 0.75);
    const double spread = (i % 2U == 0U) ? 0.03 : 0.25;
    s.x.push_back(wrap(cx + spread * (uniform() - 0.5)));
    s.y.push_back(wrap(0.5 + spread * (uniform() - 0.5)));
    s.z.push_back(wrap(0.25 + spread * (uniform() - 0.5)));
  }
  s.mass.assign(n, 1.0);
  return s;
}

[[nodiscard]] TreePmOptions baseOptions() {
  TreePmOptions options;
  options.pm_options.box_size_x_mpc_comoving = 1.0;
  options.pm_options.box_size_y_mpc_comoving = 1.0;
  options.pm_options.box_size_z_mpc_comoving = 1.0;
  options.tree_options.softening.epsilon_comoving = 2.0e-3;
  options.tree_options.max_leaf_size = 8U;
  options.tree_options.relative_force_tolerance = 1.0e-3;
  options.split_policy = makeTreePmSplitPolicyFromMeshSpacing(1.25, 6.25, 0.03125);
  options.short_range_only = true;
  options.acceptance_policy = TreePmAcceptancePolicy::kAdaptiveRelative;
  options.adaptive_error_estimator = TreePmAdaptiveErrorEstimator::kRelativeMonopole;
  options.adaptive_maximum_opening_angle = 1.0e-6;
  options.gaussian_pair_lookup_enabled = true;
  return options;
}

struct Forces {
  std::vector<double> ax, ay, az;
  TreePmDiagnostics diagnostics{};
  TreePmProfileEvent profile{};
};

[[nodiscard]] Forces solve(TreePmCoordinator& coordinator, const Sources& s, std::span<const TreeLocalIndex> active,
                           const TreePmOptions& options) {
  Forces f;
  f.ax.assign(active.size(), 0.0);
  f.ay.assign(active.size(), 0.0);
  f.az.assign(active.size(), 0.0);
  const std::vector<double> previous(active.size(), 1.0);
  const TreePmForceAccumulatorView accumulator{active, f.ax, f.ay, f.az, previous};
  coordinator.solveActiveSetWithPmCadence(s.x, s.y, s.z, s.mass, accumulator, options, false, &f.profile,
                                          &f.diagnostics);
  return f;
}

[[nodiscard]] Forces eagerReference(const Sources& s, std::span<const TreeLocalIndex> active, TreePmOptions options) {
  options.topology_refit_enabled = false;
  options.identical_source_tree_reuse_enabled = false;
  options.source_mass_generation = 0U;
  TreePmCoordinator fresh(PmGridShape{32, 32, 32});
  return solve(fresh, s, active, options);
}

void requireClose(const Forces& got, const Forces& want, double tolerance, const std::string& label) {
  requireOrThrow(got.ax.size() == want.ax.size(), label + ": extent mismatch");
  for (std::size_t i = 0; i < got.ax.size(); ++i) {
    for (const auto& [g, w] : {std::pair{got.ax[i], want.ax[i]}, std::pair{got.ay[i], want.ay[i]},
                               std::pair{got.az[i], want.az[i]}}) {
      requireOrThrow(std::isfinite(g) && std::abs(g - w) <= tolerance * std::max(1.0, std::abs(w)),
                     label + ": force mismatch at target " + std::to_string(i));
    }
  }
}

[[nodiscard]] std::vector<TreeLocalIndex> activeSubset(std::size_t n, std::size_t stride) {
  std::vector<TreeLocalIndex> active;
  for (std::size_t i = 0; i < n; i += stride) active.push_back(static_cast<TreeLocalIndex>(i));
  return active;
}

void testRetainedCachesMatchEagerReferenceAcrossTicks() {
  constexpr std::size_t n = 600U;
  Sources s = makeSources(n, 0x5DEECE66DULL);
  TreePmOptions options = baseOptions();
  options.topology_refit_enabled = true;
  options.source_layout_generation = 7U;
  options.source_mass_generation = 11U;
  TreePmCoordinator retained(PmGridShape{32, 32, 32});
  std::uint64_t generation = 1U;
  std::uint64_t refit_ticks = 0U, cache_reuses = 0U;
  const std::vector<std::size_t> strides{1U, 4U, 37U, n};  // 100%, 25%, ~3%, one target
  for (std::size_t tick = 0; tick < 12U; ++tick) {
    // Small coherent motion keeps every source inside its original leaf cell
    // most ticks (certified refit); the canonical periodic wrap is kept.
    for (std::size_t i = 0; i < n; ++i) {
      s.x[i] += 1.0e-6 * std::sin(0.37 * static_cast<double>(i + tick));
      s.y[i] += 1.0e-6 * std::cos(0.11 * static_cast<double>(i * 3U + tick));
      s.x[i] -= std::floor(s.x[i]);
      s.y[i] -= std::floor(s.y[i]);
    }
    options.source_generation = GravitySourceGeneration{++generation};
    const auto active = activeSubset(n, strides[tick % strides.size()]);
    const Forces cached = solve(retained, s, active, options);
    const Forces eager = eagerReference(s, active, options);
    requireClose(cached, eager, 1.0e-12, "retained fine-tick caches, tick " + std::to_string(tick));
    refit_ticks += cached.profile.tree_profile.refit_accepted;
    cache_reuses += cached.diagnostics.topology_cache_reuses;
    if (tick > 0U) {
      requireOrThrow(cached.diagnostics.mass_uniformity_scans == 0U && cached.diagnostics.mass_uniformity_reuses == 1U,
                     "certified mass generation did not reuse the uniformity result");
    }
    requireOrThrow(cached.diagnostics.uniform_source_mass_lane_released == 1U,
                   "uniform population did not release the tree-order mass lane");
  }
  requireOrThrow(refit_ticks > 0U, "no fine tick used a certified refit");
  requireOrThrow(cache_reuses > 0U, "topology caches were never reused across a certified refit");
}

void testRefitRejectionSeamCrossingAndReorderRebuild() {
  constexpr std::size_t n = 400U;
  Sources s = makeSources(n, 0x9E3779B97F4A7C15ULL);
  TreePmOptions options = baseOptions();
  options.topology_refit_enabled = true;
  options.source_layout_generation = 3U;
  options.source_generation = GravitySourceGeneration{1U};
  TreePmCoordinator retained(PmGridShape{32, 32, 32});
  const auto all = activeSubset(n, 1U);
  (void)solve(retained, s, all, options);
  // Large motion and seam crossing: sources near x = 0 move across the seam.
  for (std::size_t i = 0; i < n; ++i) {
    s.x[i] -= 0.01;
    s.x[i] -= std::floor(s.x[i]);
  }
  options.source_generation = GravitySourceGeneration{2U};
  const Forces moved = solve(retained, s, all, options);
  requireOrThrow(moved.profile.tree_profile.build_count == 1U && moved.diagnostics.topology_cache_rebuilds > 0U,
                 "leaf escape across the seam did not rebuild the tree and topology caches");
  requireClose(moved, eagerReference(s, all, options), 1.0e-12, "rebuild after refit rejection");
  // Dense-row reorder: permute rows and advance the layout generation.
  std::vector<std::size_t> perm(n);
  std::iota(perm.begin(), perm.end(), 0U);
  std::reverse(perm.begin(), perm.end());
  Sources reordered;
  for (const std::size_t p : perm) {
    reordered.x.push_back(s.x[p]); reordered.y.push_back(s.y[p]); reordered.z.push_back(s.z[p]);
    reordered.mass.push_back(s.mass[p]);
  }
  options.source_layout_generation = 4U;
  options.source_generation = GravitySourceGeneration{3U};
  const Forces after_reorder = solve(retained, reordered, all, options);
  requireOrThrow(after_reorder.profile.tree_profile.build_count == 1U &&
                     after_reorder.profile.tree_profile.rebuild_reason == 8U,
                 "dense-row reorder did not invalidate the topology");
  for (std::size_t i = 0; i < n; ++i) {
    // Row i of the reordered view is row perm[i] of the original.
    const double want = moved.ax[perm[i]];
    requireOrThrow(std::abs(after_reorder.ax[i] - want) <= 1.0e-12 * std::max(1.0, std::abs(want)),
                   "reordered rows did not carry their own forces");
  }
}

void testUniformScalarPathMatchesExplicitMassLane() {
  constexpr std::size_t n = 300U;
  Sources uniform = makeSources(n, 0x2545F4914F6CDD1DULL);
  // Explicit-lane twin: identical except one extra heavy source far outside
  // every target's cutoff; the active set excludes it.
  Sources explicit_lane = uniform;
  explicit_lane.x.push_back(0.5); explicit_lane.y.push_back(0.0); explicit_lane.z.push_back(0.75);
  explicit_lane.mass.push_back(3.0);
  Sources uniform_with_far = uniform;
  uniform_with_far.x.push_back(0.5); uniform_with_far.y.push_back(0.0); uniform_with_far.z.push_back(0.75);
  uniform_with_far.mass.push_back(1.0);
  TreePmOptions options = baseOptions();
  options.source_generation = GravitySourceGeneration{1U};
  std::vector<TreeLocalIndex> active;
  for (std::size_t i = 0; i < n; ++i) {
    const double dy = std::min(uniform.y[i], 1.0 - uniform.y[i]);
    if (dy > 0.3) active.push_back(static_cast<TreeLocalIndex>(i));  // far from the extra source
  }
  requireOrThrow(!active.empty(), "fixture has no target far from the extra source");
  TreePmCoordinator a(PmGridShape{32, 32, 32});
  TreePmCoordinator b(PmGridShape{32, 32, 32});
  const Forces scalar = solve(a, uniform_with_far, active, options);
  const Forces lane = solve(b, explicit_lane, active, options);
  requireOrThrow(scalar.diagnostics.uniform_source_mass_lane_released == 1U &&
                     lane.diagnostics.uniform_source_mass_lane_released == 0U,
                 "fixture did not exercise both mass representations");
  requireClose(scalar, lane, 0.0, "uniform scalar mass versus explicit mass lane");
  // Exact scaling: doubling a uniform mass doubles every pair term exactly.
  Sources doubled = uniform;
  std::fill(doubled.mass.begin(), doubled.mass.end(), 2.0);
  TreePmCoordinator c(PmGridShape{32, 32, 32});
  TreePmCoordinator d(PmGridShape{32, 32, 32});
  const Forces one = solve(c, uniform, active, options);
  const Forces two = solve(d, doubled, active, options);
  for (std::size_t i = 0; i < active.size(); ++i) {
    requireOrThrow(two.ax[i] == 2.0 * one.ax[i] && two.ay[i] == 2.0 * one.ay[i] && two.az[i] == 2.0 * one.az[i],
                   "uniform mass scaling is not exact");
  }
}

void testOpenMpThreadCountDeterminism() {
#if defined(_OPENMP)
  constexpr std::size_t n = 500U;
  const Sources s = makeSources(n, 0xA0761D6478BD642FULL);
  const auto active = activeSubset(n, 2U);
  for (const bool group_walk : {false, true}) {
    TreePmOptions options = baseOptions();
    options.group_walk_enabled = group_walk;
    options.group_walk_target_group_size = group_walk ? 32U : 0U;
    options.source_generation = GravitySourceGeneration{1U};
    const int saved = omp_get_max_threads();
    omp_set_num_threads(1);
    TreePmCoordinator serial(PmGridShape{32, 32, 32});
    const Forces one = solve(serial, s, active, options);
    omp_set_num_threads(std::max(4, saved));
    TreePmCoordinator threaded(PmGridShape{32, 32, 32});
    const Forces many = solve(threaded, s, active, options);
    omp_set_num_threads(saved);
    requireClose(many, one, 0.0, std::string("OpenMP thread-count determinism, group_walk=") +
                                     (group_walk ? "on" : "off"));
  }
#endif
}

}  // namespace

int main() {
#if COSMOSIM_ENABLE_MPI
  MPI_Init(nullptr, nullptr);
#endif
  int failed = 0;
  try {
#if COSMOSIM_ENABLE_MPI
    // Rank-local serial coordinators; every rank runs the same fixture.
#endif
    testRetainedCachesMatchEagerReferenceAcrossTicks();
    testRefitRejectionSeamCrossingAndReorderRebuild();
    testUniformScalarPathMatchesExplicitMassLane();
    testOpenMpThreadCountDeterminism();
  } catch (const std::exception& error) {
    std::fprintf(stderr, "fine-tick cache regression failed: %s\n", error.what());
    failed = 1;
  }
#if COSMOSIM_ENABLE_MPI
  MPI_Finalize();
#endif
  return failed;
}
