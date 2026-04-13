#include <cassert>
#include <cmath>

#include "cosmosim/gravity/tree_pm_coupling.hpp"

namespace {

void testSplitKernelComplementarity() {
  cosmosim::gravity::TreePmSplitPolicy split_policy;
  split_policy.split_scale_comoving = 0.2;

  const double radii[] = {0.01, 0.05, 0.2, 0.4, 1.0};
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

void testDiagnosticsContinuityAtSplitScale() {
  cosmosim::gravity::TreePmSplitPolicy split_policy;
  split_policy.split_scale_comoving = 0.125;

  const cosmosim::gravity::TreePmDiagnostics diagnostics = cosmosim::gravity::computeTreePmDiagnostics(split_policy);
  assert(std::abs(diagnostics.split_scale_comoving - split_policy.split_scale_comoving) < 1.0e-15);
  assert(diagnostics.short_range_factor_at_split > 0.0);
  assert(diagnostics.long_range_factor_at_split > 0.0);
  assert(diagnostics.composition_error_at_split < 1.0e-12);
  assert(diagnostics.max_relative_composition_error < 1.0e-12);
}

}  // namespace

int main() {
  testSplitKernelComplementarity();
  testDiagnosticsContinuityAtSplitScale();
  return 0;
}
