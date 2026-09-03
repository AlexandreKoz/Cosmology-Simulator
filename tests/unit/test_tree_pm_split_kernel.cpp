#include <cassert>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>

#include "cosmosim/gravity/tree_pm_coupling.hpp"
#include "gravity/internal/tree_pm_transport_planner.hpp"

namespace {

void testSplitKernelComplementarity() {
  const cosmosim::gravity::TreePmSplitPolicy split_policy =
      cosmosim::gravity::makeTreePmSplitPolicyFromMeshSpacing(1.6, 5.0, 0.125);

  const double radii[] = {0.01, 0.05, 0.2, 0.4, 0.8};
  for (const double radius : radii) {
    const double short_factor =
        cosmosim::gravity::treePmGaussianShortRangeForceFactor(radius, split_policy.split_scale_comoving);
    const double long_factor =
        cosmosim::gravity::treePmGaussianLongRangeForceFactor(radius, split_policy.split_scale_comoving);
    assert(short_factor >= 0.0);
    assert(short_factor <= 1.0 + 1.0e-12);
    assert(long_factor >= -1.0e-12);
    assert(long_factor <= 1.0);
    assert(std::abs(short_factor + long_factor - 1.0) < 1.0e-12);
  }
}


void testFiniteSofteningResidualComposesWithPmLongRange() {
  constexpr double split_scale = 0.2;
  const double radii[] = {0.02, 0.05, 0.1, 0.2, 0.5, 1.0};
  const double epsilon_over_split[] = {0.0, 0.008, 0.016637952, 0.025, 0.05, 0.10, 0.20};

  for (const double epsilon_ratio : epsilon_over_split) {
    const double epsilon = epsilon_ratio * split_scale;
    for (const double radius : radii) {
      const double r2 = radius * radius;
      const double newton_inv_r3 = 1.0 / (r2 * radius);
      const double pm_long =
          cosmosim::gravity::treePmGaussianLongRangeForceFactor(radius, split_scale) *
          newton_inv_r3;
      const double tree_residual =
          cosmosim::gravity::treePmSoftenedShortRangeInvR3(r2, epsilon, split_scale);
      const double expected = cosmosim::gravity::softenedInvR3(r2, epsilon);
      const double scale = std::max(std::abs(expected), 1.0e-30);
      assert(std::abs((tree_residual + pm_long) - expected) / scale < 5.0e-13);
    }
  }
}

void testMeshCellDerivedSplitSemantics() {
  const double mesh_spacing = 0.025;
  const double asmth_cells = 1.25;
  const double rcut_cells = 4.5;
  const cosmosim::gravity::TreePmSplitPolicy split_policy =
      cosmosim::gravity::makeTreePmSplitPolicyFromMeshSpacing(asmth_cells, rcut_cells, mesh_spacing);

  assert(std::abs(split_policy.mesh_spacing_comoving - mesh_spacing) < 1.0e-15);
  assert(std::abs(split_policy.asmth_cells - asmth_cells) < 1.0e-15);
  assert(std::abs(split_policy.rcut_cells - rcut_cells) < 1.0e-15);
  assert(std::abs(split_policy.split_scale_comoving - asmth_cells * mesh_spacing) < 1.0e-15);
  assert(std::abs(split_policy.cutoff_radius_comoving - rcut_cells * mesh_spacing) < 1.0e-15);
}

void testDiagnosticsContinuityAtSplitScale() {
  const cosmosim::gravity::TreePmSplitPolicy split_policy =
      cosmosim::gravity::makeTreePmSplitPolicyFromMeshSpacing(2.0, 6.0, 0.0625);

  const cosmosim::gravity::TreePmDiagnostics diagnostics = cosmosim::gravity::computeTreePmDiagnostics(split_policy);
  assert(std::abs(diagnostics.mesh_spacing_comoving - split_policy.mesh_spacing_comoving) < 1.0e-15);
  assert(std::abs(diagnostics.asmth_cells - split_policy.asmth_cells) < 1.0e-15);
  assert(std::abs(diagnostics.rcut_cells - split_policy.rcut_cells) < 1.0e-15);
  assert(std::abs(diagnostics.split_scale_comoving - split_policy.split_scale_comoving) < 1.0e-15);
  assert(std::abs(diagnostics.cutoff_radius_comoving - split_policy.cutoff_radius_comoving) < 1.0e-15);
  assert(diagnostics.short_range_factor_at_split > 0.0);
  assert(diagnostics.long_range_factor_at_split > 0.0);
  assert(diagnostics.short_range_factor_at_cutoff >= 0.0);
  assert(diagnostics.short_range_factor_at_cutoff < 0.5);
  assert(diagnostics.long_range_factor_at_cutoff > 0.5);
  assert(diagnostics.composition_error_at_split < 1.0e-12);
  assert(diagnostics.max_relative_composition_error < 1.0e-12);
}

void testSparseTreePmAggregateRoundPlanner() {
  using cosmosim::gravity::internal::planSparseTreePmRound;
  using cosmosim::gravity::internal::sparseTreePmPhysicalRoundCount;

  const auto zero_peer_plan = planSparseTreePmRound(4096U, 0U, 96U, 80U, 1024U);
  assert(zero_peer_plan.targets_per_peer_per_round == 4096U);
  assert(zero_peer_plan.aggregate_bytes_per_target == 0U);

  const auto many_peer_plan = planSparseTreePmRound(4096U, 8U, 96U, 80U, 4096U);
  assert(many_peer_plan.aggregate_bytes_per_target == 8U * 96U);
  assert(many_peer_plan.targets_per_peer_per_round == 5U);
  assert(many_peer_plan.targets_per_peer_per_round *
             many_peer_plan.aggregate_bytes_per_target <= 4096U);

  const std::uint64_t logical_targets =
      static_cast<std::uint64_t>(std::numeric_limits<int>::max()) + 1000000ULL;
  const std::uint64_t physical_rounds = sparseTreePmPhysicalRoundCount(
      logical_targets, many_peer_plan.targets_per_peer_per_round);
  assert(physical_rounds > 1U);
  assert(physical_rounds ==
         1U + (logical_targets - 1U) /
             static_cast<std::uint64_t>(many_peer_plan.targets_per_peer_per_round));

  bool threw = false;
  try {
    (void)planSparseTreePmRound(1U, 64U, 96U, 80U, 4096U);
  } catch (const std::overflow_error&) {
    threw = true;
  }
  assert(threw);
}

}  // namespace

int main() {
  testSplitKernelComplementarity();
  testFiniteSofteningResidualComposesWithPmLongRange();
  testMeshCellDerivedSplitSemantics();
  testDiagnosticsContinuityAtSplitScale();
  testSparseTreePmAggregateRoundPlanner();
  return 0;
}
