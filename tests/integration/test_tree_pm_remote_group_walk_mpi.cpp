// Distributed remote (incoming LET) group-walk regression (2026-10-26).
//
// Owner-executed on 2 and 4 ranks. For each configuration the same
// distributed decomposition is solved twice — per-target remote walk
// (group_walk_enabled = false, the reference) and grouped remote walk — and
// every owned target's short-range acceleration is compared component-wise.
// The configuration opens every node (adaptive_maximum_opening_angle ~ 0) so
// both walks evaluate the identical in-cutoff direct pair set; summed incoming
// pair counts must therefore agree exactly, which also pins the coincident-
// pair convention of the remote group kernel. The test refuses to pass on a
// single rank and requires the grouped remote path to have run.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/gravity/tree_pm_coupling.hpp"
#include "cosmosim/parallel/distributed_memory.hpp"

#if COSMOSIM_ENABLE_MPI
#include <mpi.h>
#endif

namespace {

#if COSMOSIM_ENABLE_MPI

using namespace cosmosim::gravity;

void requireOrThrow(bool condition, const std::string& message) {
  if (!condition) throw std::runtime_error(message);
}

[[nodiscard]] std::uint64_t globalSum(std::uint64_t value) {
  std::uint64_t total = 0U;
  MPI_Allreduce(&value, &total, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
  return total;
}

[[nodiscard]] bool allRanks(bool value) {
  int local = value ? 1 : 0, global = 0;
  MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  return global == 1;
}

enum class Ownership { kSlabByX, kAsymmetric, kRoundRobin };

struct Case {
  std::string name;
  std::size_t group_size = 0U;      // 0 = legacy one-leaf groups
  std::size_t max_leaf_size = 8U;
  double epsilon = 2.0e-3;
  Ownership ownership = Ownership::kSlabByX;
  bool partial_history = false;     // mixed eligible / history-fallback groups
};

struct Population {
  std::vector<double> x, y, z, mass, previous;
};

// Clusters straddling the periodic seams (x=0/1, corner), a dense interior
// core, a sparse uniform background and one exactly coincident pair.
[[nodiscard]] Population makePopulation(bool partial_history) {
  Population p;
  std::uint64_t state = 0x243F6A8885A308D3ULL;
  const auto uniform = [&state]() {
    state = state * 6364136223846793005ULL + 1442695040888963407ULL;
    return static_cast<double>(state >> 11U) * (1.0 / 9007199254740992.0);
  };
  const auto wrap = [](double v) { v -= std::floor(v); return v >= 1.0 ? 0.0 : v; };
  const std::array<std::array<double, 4>, 5> clusters{{
      {0.005, 0.5, 0.5, 0.02},   // across the x seam
      {0.995, 0.995, 0.005, 0.02},  // corner seam
      {0.5, 0.5, 0.5, 0.015},    // dense core
      {0.26, 0.74, 0.21, 0.05},
      {0.75, 0.25, 0.8, 0.03}}};
  for (const auto& c : clusters) {
    for (int i = 0; i < 72; ++i) {
      p.x.push_back(wrap(c[0] + c[3] * (uniform() - 0.5)));
      p.y.push_back(wrap(c[1] + c[3] * (uniform() - 0.5)));
      p.z.push_back(wrap(c[2] + c[3] * (uniform() - 0.5)));
    }
  }
  for (int i = 0; i < 48; ++i) {  // sparse background
    p.x.push_back(uniform()); p.y.push_back(uniform()); p.z.push_back(uniform());
  }
  // Exactly coincident distinct pair (will usually land on different ranks
  // under round-robin ownership).
  p.x.push_back(0.4321); p.y.push_back(0.1234); p.z.push_back(0.5678);
  p.x.push_back(0.4321); p.y.push_back(0.1234); p.z.push_back(0.5678);
  p.mass.assign(p.x.size(), 1.0);
  p.previous.assign(p.x.size(), 50.0);
  if (partial_history) {
    for (std::size_t i = 0; i < p.previous.size(); i += 5U) {
      p.previous[i] = std::numeric_limits<double>::quiet_NaN();
    }
  }
  return p;
}

[[nodiscard]] int ownerOf(std::size_t i, double x, Ownership ownership, int world_size) {
  switch (ownership) {
    case Ownership::kSlabByX:
      return std::min(world_size - 1, static_cast<int>(x * static_cast<double>(world_size)));
    case Ownership::kAsymmetric:
      // ~80% on rank 0, the rest spread over the others.
      return (i % 5U != 0U || world_size == 1) ? 0 : 1 + static_cast<int>((i / 5U) % static_cast<std::size_t>(world_size - 1));
    case Ownership::kRoundRobin:
      return static_cast<int>(i % static_cast<std::size_t>(world_size));
  }
  return 0;
}

[[nodiscard]] TreePmOptions makeOptions(const Case& c) {
  TreePmOptions options;
  options.pm_options.box_size_x_mpc_comoving = 1.0;
  options.pm_options.box_size_y_mpc_comoving = 1.0;
  options.pm_options.box_size_z_mpc_comoving = 1.0;
  options.tree_options.softening.epsilon_comoving = c.epsilon;
  options.tree_options.max_leaf_size = c.max_leaf_size;
  options.tree_options.relative_force_tolerance = 1.0e-3;
  options.split_policy = makeTreePmSplitPolicyFromMeshSpacing(1.25, 6.25, 0.03125);
  options.short_range_only = true;
  options.acceptance_policy = TreePmAcceptancePolicy::kAdaptiveRelative;
  options.adaptive_error_estimator = TreePmAdaptiveErrorEstimator::kRelativeMonopole;
  options.adaptive_maximum_opening_angle = 1.0e-6;  // open every node: all-direct pair sets
  options.gaussian_pair_lookup_enabled = true;
  options.source_generation = GravitySourceGeneration{1U};
  options.source_layout_generation = 1U;
  options.group_walk_target_group_size = c.group_size;
  return options;
}

struct Solve {
  std::vector<double> ax, ay, az;
  TreePmDiagnostics diagnostics{};
};

[[nodiscard]] Solve solveDistributed(const Population& local, const TreePmOptions& options) {
  int world_size = 1, world_rank = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  const PmGridShape shape{32, 32, 32};
  const auto layout = cosmosim::parallel::makePmSlabLayout(shape.nx, shape.ny, shape.nz, world_size, world_rank);
  TreePmCoordinator coordinator(shape, layout);
  Solve s;
  s.ax.assign(local.x.size(), 0.0);
  s.ay.assign(local.x.size(), 0.0);
  s.az.assign(local.x.size(), 0.0);
  std::vector<TreeLocalIndex> active(local.x.size());
  for (std::size_t i = 0; i < active.size(); ++i) active[i] = static_cast<TreeLocalIndex>(i);
  const TreePmForceAccumulatorView accumulator{active, s.ax, s.ay, s.az, local.previous};
  coordinator.solveActiveSetWithPmCadence(local.x, local.y, local.z, local.mass, accumulator, options,
                                          /*refresh_long_range_field=*/false, nullptr, &s.diagnostics);
  return s;
}

void runCase(const Case& c, std::uint64_t& global_grouped_remote, std::uint64_t& global_fallback_remote) {
  int world_size = 1, world_rank = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  const Population all = makePopulation(c.partial_history);
  Population local;
  for (std::size_t i = 0; i < all.x.size(); ++i) {
    if (ownerOf(i, all.x[i], c.ownership, world_size) != world_rank) continue;
    local.x.push_back(all.x[i]); local.y.push_back(all.y[i]); local.z.push_back(all.z[i]);
    local.mass.push_back(all.mass[i]); local.previous.push_back(all.previous[i]);
  }
  TreePmOptions options = makeOptions(c);
  options.group_walk_enabled = false;
  const Solve reference = solveDistributed(local, options);
  options.group_walk_enabled = true;
  const Solve grouped = solveDistributed(local, options);

  double max_err = 0.0;
  std::size_t worst = 0U;
  for (std::size_t i = 0; i < local.x.size(); ++i) {
    for (const auto& [got, want] : {std::pair{grouped.ax[i], reference.ax[i]}, std::pair{grouped.ay[i], reference.ay[i]},
                                    std::pair{grouped.az[i], reference.az[i]}}) {
      const double err = std::abs(got - want) / std::max(1.0, std::abs(want));
      if (!std::isfinite(got) || err > max_err) { max_err = std::isfinite(got) ? err : 1.0e300; worst = i; }
    }
  }
  const std::string label = c.name + " (ranks=" + std::to_string(world_size) + ")";
  requireOrThrow(allRanks(max_err <= 1.0e-9),
                 label + ": grouped remote walk differs from per-target walk (rank " + std::to_string(world_rank) +
                     " target " + std::to_string(worst) + " rel=" + std::to_string(max_err) + ")");

  const auto& g_in = grouped.diagnostics.incoming_traversal;
  const auto& r_in = reference.diagnostics.incoming_traversal;
  // Identical in-cutoff pair set (all-direct): exact incoming pair agreement,
  // including the coincident distinct pair counted on the remote side.
  requireOrThrow(globalSum(g_in.direct_pair_evaluations) == globalSum(r_in.direct_pair_evaluations),
                 label + ": incoming direct pair count differs between grouped and per-target walks");
  requireOrThrow(globalSum(grouped.diagnostics.local_traversal.direct_pair_evaluations) ==
                     globalSum(reference.diagnostics.local_traversal.direct_pair_evaluations),
                 label + ": owned direct pair count differs between grouped and per-target walks");
  // Exactly-once accounting of incoming targets.
  requireOrThrow(globalSum(g_in.targets) == globalSum(r_in.targets),
                 label + ": incoming target evaluation count changed");
  requireOrThrow(allRanks(g_in.group_walk_grouped_targets + g_in.group_walk_fallback_targets <= g_in.targets),
                 label + ": incoming grouped + fallback targets exceed evaluated targets");
  requireOrThrow(allRanks(r_in.group_walk_groups == 0U), label + ": reference run used the group walk");
  if (c.group_size != 0U) {
    requireOrThrow(allRanks(g_in.group_walk_target_group_size_max <= c.group_size),
                   label + ": remote group exceeded the configured target group size");
  }
  global_grouped_remote += globalSum(g_in.group_walk_grouped_targets);
  global_fallback_remote += globalSum(g_in.group_walk_fallback_history);
}

void testRemoteGroupWalkMatchesPerTargetWalk() {
  int world_size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  requireOrThrow(world_size >= 2, "remote group-walk regression requires at least two MPI ranks");
  const std::vector<Case> cases{
      {"slab ownership, legacy groups", 0U, 8U, 2.0e-3, Ownership::kSlabByX, false},
      {"slab ownership, group 16", 16U, 8U, 2.0e-3, Ownership::kSlabByX, false},
      {"slab ownership, group 64, leaf 4", 64U, 4U, 2.0e-3, Ownership::kSlabByX, false},
      {"asymmetric ownership, group 32", 32U, 8U, 2.0e-3, Ownership::kAsymmetric, false},
      {"round-robin ownership, group 32, wide softening", 32U, 8U, 6.0e-3, Ownership::kRoundRobin, false},
      {"round-robin ownership, mixed history fallback", 32U, 8U, 2.0e-3, Ownership::kRoundRobin, true},
      {"slab ownership, leaf 16, group 16", 16U, 16U, 1.0e-3, Ownership::kSlabByX, true},
  };
  std::uint64_t grouped_remote = 0U, fallback_remote = 0U;
  for (const Case& c : cases) runCase(c, grouped_remote, fallback_remote);
  requireOrThrow(grouped_remote > 0U, "the grouped remote (incoming LET) path never ran");
  requireOrThrow(fallback_remote > 0U, "mixed-history cases never exercised the remote per-target fallback");
}

#endif

}  // namespace

int main() {
#if COSMOSIM_ENABLE_MPI
  MPI_Init(nullptr, nullptr);
  int failed = 0;
  try {
    testRemoteGroupWalkMatchesPerTargetWalk();
  } catch (const std::exception& error) {
    std::fprintf(stderr, "remote group-walk regression failed: %s\n", error.what());
    failed = 1;
  }
  int any_failed = 0;
  MPI_Allreduce(&failed, &any_failed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Finalize();
  return any_failed == 0 ? 0 : 1;
#else
  std::fprintf(stderr, "remote group-walk regression requires COSMOSIM_ENABLE_MPI\n");
  return 1;
#endif
}
