#include <cassert>
#include <cmath>

#include "cosmosim/gravity/tree_pm_coupling.hpp"

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

}  // namespace

int main() {
  testSplitKernelComplementarity();
  testMeshCellDerivedSplitSemantics();
  testDiagnosticsContinuityAtSplitScale();
  return 0;
}
