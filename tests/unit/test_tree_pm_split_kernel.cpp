#include "../support/test_require.hpp"
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>

#include "cosmosim/gravity/tree_pm_coupling.hpp"
#include "gravity/internal/tree_pm_transport_planner.hpp"
#include "gravity/internal/tree_pm_geometry.hpp"
#include "gravity/internal/tree_pm_acceptance.hpp"

namespace {

// Independent long-double integral series; production uses a fixed polynomial
// or erf, whereas this reference sums the convergent series through tolerance.
long double referenceH(long double t) {
  const long double four_inv_sqrt_pi = 4.0L/std::sqrt(std::acos(-1.0L));
  if (t < 1.0L) {
    long double sum = 0.0L, term = 1.0L;
    for (unsigned n = 0; n < 80U; ++n) {
      const long double add = term/static_cast<long double>(2U*n+3U);
      sum += add;
      if (std::abs(add) < 1e-30L) break;
      term *= -t/static_cast<long double>(n+1U);
    }
    return four_inv_sqrt_pi*sum;
  }
  const long double q = std::sqrt(t);
  return (std::erf(q) - 0.5L*four_inv_sqrt_pi*q*std::exp(-t))/(t*q);
}

void testStableZeroAndTinyRadius() {
  using namespace cosmosim::gravity;
  for (const double a : {1e-90, 1e-10, 0.2, 1.0, 1e10, 1e90}) {
    CHUI_REQUIRE(treePmGaussianShortRangeForceFactor(0.0,a) == 1.0);
    CHUI_REQUIRE(treePmGaussianLongRangeForceFactor(0.0,a) == 0.0);
    const long double scale = static_cast<long double>(a)*a*a;
    const long double limit = 1.0L/(6.0L*std::sqrt(std::acos(-1.0L)));
    CHUI_REQUIRE(std::abs(treePmGaussianLongRangeInvR3Unchecked(0.0,a)*scale-limit) < 3e-16L);
    for (const double ratio : {1e-16, 1e-10, 1e-5, 0.249999999, 0.25, 0.250000001, 1.0, 5.0, 16.0}) {
      const double r = a*ratio;
      const long double q = static_cast<long double>(r)/(2.0L*a);
      const long double expected = referenceH(q*q)/8.0L;
      const double coefficient = treePmGaussianLongRangeInvR3Unchecked(r,a);
      CHUI_REQUIRE(std::isfinite(coefficient));
      CHUI_REQUIRE(std::abs(coefficient*scale-expected) < 3e-14L*expected);
    }
    for (const double eps_ratio : {0.001, 0.01, 0.05, 0.2}) {
      const double eps = a*eps_ratio;
      const double coincident = treePmSoftenedShortRangeInvR3(0.0,eps,a);
      const long double expected_zero = std::pow(static_cast<long double>(eps),-3.0L)-limit/scale;
      CHUI_REQUIRE(std::abs((coincident-expected_zero)*scale) < 1e-6L);
      CHUI_REQUIRE(treePmSoftenedShortRangeInvR3(0.0,0.0,a) == softenedInvR3(0.0,0.0));
      const double r = a*1e-12;
      const double residual = treePmSoftenedShortRangeInvR3(r*r,eps,a);
      CHUI_REQUIRE(std::abs((residual-coincident)*scale) < 1e-6L);
      const long double softened = std::pow(static_cast<long double>(r)*r +
          static_cast<long double>(eps)*eps, -1.5L);
      CHUI_REQUIRE(std::abs((residual-softened)*scale + limit) < 1e-6L);
    }
  }
}

void testCoincidentResidualLookupConsistency() {
  using namespace cosmosim::gravity;
  const TreePmGaussianCoefficientTable first, second;
  for (const double a : {1e-90, 1e-5, 0.25, 1.0, 1e90}) {
    const long double scale = static_cast<long double>(a)*a*a;
    const long double c = 1.0L/(6.0L*std::sqrt(std::acos(-1.0L)));
    for (const double ratio : {0.05, 0.5, 1.0, 3.0}) {
      const double eps = a*ratio;
      const double zero = treePmSoftenedShortRangeInvR3(0.0, eps, a);
      const long double expected = std::pow(static_cast<long double>(eps)/a, -3.0L) - c;
      CHUI_REQUIRE(std::abs(zero*scale-expected) < 8e-12L);
      const double lookup = softenedInvR3(0.0, eps)-first.longRangeInvR3(0.0, a);
      CHUI_REQUIRE(std::abs((zero-lookup)*scale) < 8e-12L);
      for (const double r_ratio : {1e-16, 1e-12}) {
        const double r = a*r_ratio;
        CHUI_REQUIRE(std::abs((treePmSoftenedShortRangeInvR3(r*r, eps, a)-zero)*scale) < 8e-12L);
      }
    }
    for (const double t : {0.0, 0.015625, 0.5, 63.99, 64.0, 100.0}) {
      const double r2 = 4*a*a*t;
      const double before = first.longRangeInvR3(r2, a);
      CHUI_REQUIRE(before == second.longRangeInvR3(r2, a));
      (void)first.longRangeInvR3(4*a*a, a);
      CHUI_REQUIRE(before == first.longRangeInvR3(r2, a));
    }
  }
}

void testGaussianLookupAcrossInterval() {
  using namespace cosmosim::gravity;
  const TreePmGaussianCoefficientTable table;
  for (const double a : {1e-50, 1e-3, 0.2, 1.0, 1e50}) {
    const long double scale = static_cast<long double>(a)*a*a;
    // Exercise cell interiors, endpoints, near-zero analytic branch and the
    // last cell. Absolute budget is on H, independent of a and epsilon.
    for (unsigned i = 0; i <= 131072U; ++i) {
      const double t = 64.0*static_cast<double>(i)/131072.0;
      const double r2 = 4.0*a*a*t;
      const double lr = table.longRangeInvR3(r2,a);
      const long double actual_t = 0.25L*static_cast<long double>(r2)/a/a;
      const long double expected_h = referenceH(actual_t);
      CHUI_REQUIRE(std::abs(8.0L*lr*scale - expected_h) <= 5e-11L);
      if (i % 127U == 0U && r2 > 0.0) {
        for (const double eps_ratio : {0.0, 0.008, 0.05, 0.2}) {
          const double eps = a*eps_ratio;
          const double lookup_residual = softenedInvR3(r2,eps)-lr;
          const double analytic = treePmSoftenedShortRangeInvR3(r2,eps,a);
          CHUI_REQUIRE(std::abs((lookup_residual-analytic)*scale) < 6.3e-12L);
        }
      }
    }
    for (const double ratio : {std::nextafter(16.0,0.0), 16.0,
         std::nextafter(16.0,17.0), 32.0}) {
      const double r = ratio*a;
      CHUI_REQUIRE(std::abs((table.longRangeInvR3(r*r,a) -
          treePmGaussianLongRangeInvR3Unchecked(r,a))*scale) < 6.3e-12L);
    }
  }
}

void testInvalidKernelArguments() {
  using namespace cosmosim::gravity;
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double inf = std::numeric_limits<double>::infinity();
  // The true zero-radius coefficient is outside double range here; preserve
  // its sign rather than producing NaN through separately overflowing terms.
  const double positive_overflow = treePmSoftenedShortRangeInvR3(0.0, 1e-110, 1e-110);
  const double negative_overflow = treePmSoftenedShortRangeInvR3(0.0, 1e-108, 1e-110);
  CHUI_REQUIRE(std::isinf(positive_overflow) && positive_overflow > 0.0);
  CHUI_REQUIRE(std::isinf(negative_overflow) && negative_overflow < 0.0);
  CHUI_REQUIRE(treePmSoftenedShortRangeInvR3(0.0, 1e110, 1e110) == 0.0);
  for (const double r : {-1.0,nan,inf}) {
    bool rejected = false;
    try { (void)treePmGaussianShortRangeForceFactor(r,1.0); }
    catch (const std::invalid_argument&) { rejected = true; }
    CHUI_REQUIRE(rejected);
    rejected = false;
    try { (void)treePmGaussianLongRangeForceFactor(r,1.0); }
    catch (const std::invalid_argument&) { rejected = true; }
    CHUI_REQUIRE(rejected);
  }
  for (const double a : {0.0,-1.0,nan,inf}) {
    bool rejected = false;
    try { (void)treePmSoftenedShortRangeInvR3(1.0,0.01,a); }
    catch (const std::invalid_argument&) { rejected = true; }
    CHUI_REQUIRE(rejected);
  }
  for (const double invalid : {-1.0,nan,inf}) {
    bool rejected = false;
    try { (void)treePmSoftenedShortRangeInvR3(invalid,0.01,1.0); }
    catch (const std::invalid_argument&) { rejected = true; }
    CHUI_REQUIRE(rejected);
    rejected = false;
    try { (void)treePmSoftenedShortRangeInvR3(1.0,invalid,1.0); }
    catch (const std::invalid_argument&) { rejected = true; }
    CHUI_REQUIRE(rejected);
  }
  for (const double a : {0.0, -1.0, nan, inf}) {
    bool rejected = false;
    try { (void)treePmSoftenedShortRangeInvR3(0.0, 0.1, a); }
    catch (const std::invalid_argument&) { rejected = true; }
    CHUI_REQUIRE(rejected);
  }
  for (const double eps : {-1.0, nan, inf}) {
    bool rejected = false;
    try { (void)treePmSoftenedShortRangeInvR3(0.0, eps, 1.0); }
    catch (const std::invalid_argument&) { rejected = true; }
    CHUI_REQUIRE(rejected);
  }
  auto unsupported = makeTreePmSplitPolicyFromMeshSpacing(1.25, 6.25, 0.01);
  unsupported.kernel = static_cast<TreePmSplitKernel>(99);
  bool rejected = false;
  try { validateTreePmSplitPolicy(unsupported); }
  catch (const std::invalid_argument& error) {
    rejected = std::string(error.what()).find("Unsupported TreePM split kernel") != std::string::npos;
  }
  CHUI_REQUIRE(rejected);
  // Finite valid inputs with enormous dimensionless q must avoid inf*0.
  CHUI_REQUIRE(treePmGaussianShortRangeForceFactor(1e300,1e-300) == 0.0);
  CHUI_REQUIRE(treePmGaussianLongRangeForceFactor(1e300,1e-300) == 1.0);
}

void testSquaredGeometryAndStagedMac() {
  using namespace cosmosim::gravity;
  using namespace cosmosim::gravity::internal;
  for (const double cutoff : {1e-150,0.25,1.0,1e150}) {
    const double c2 = cutoff*cutoff;
    for (const double d2 : {std::nextafter(c2,0.0),c2,std::nextafter(c2,
         std::numeric_limits<double>::infinity()),c2*0.5,c2*2.0}) {
      CHUI_REQUIRE(treePmSquaredDistanceWithinCutoff(d2,cutoff,c2) == (std::sqrt(d2) <= cutoff));
    }
  }
  const TreePmAabbDistances box(0.3,-0.4,0.5,0.1);
  CHUI_REQUIRE(std::abs(box.minimum_squared-(0.04+0.09+0.16)) < 1e-15);
  CHUI_REQUIRE(std::abs(box.maximumSquared(0.1)-(0.16+0.25+0.36)) < 1e-15);
  TreePmOptions full;
  full.split_policy = makeTreePmSplitPolicyFromMeshSpacing(1.25,6.25,0.01);
  full.tree_options.multipole_order = TreeMultipoleOrder::kQuadrupole;
  full.acceptance_policy = TreePmAcceptancePolicy::kAdaptiveRelative;
  auto fast = full;
  fast.full_mac_diagnostics = false;
  for (unsigned i = 1U; i <= 10000U; ++i) {
    TreeNodeAcceptanceInput in{.target_inside_node = i%7U == 0U,
        .half_size = 0.001*static_cast<double>(i%17U+1U), .com_center_offset = 0.0002,
        .node_mass_code = 1.0, .r2 = 0.001*static_cast<double>(i%101U+1U),
        .previous_acceleration_available = i%3U != 0U,
        .previous_acceleration_magnitude_code = 0.01*static_cast<double>(i),
        .target_softening_comoving = 0.001, .node_softening_min_comoving = 0.0001,
        .node_softening_max_comoving = 0.002};
    TreePmTraversalCounters a,b;
    const bool inside_cutoff = i%5U != 0U;
    const bool accepted = acceptTreePmNode(in,inside_cutoff,full,a);
    CHUI_REQUIRE(accepted == acceptTreePmNode(in,inside_cutoff,fast,b));
    CHUI_REQUIRE(a.skipped_mac_evaluations == 0U);
    CHUI_REQUIRE(a.maximum_angle_rejections == b.maximum_angle_rejections);
    CHUI_REQUIRE(a.softening_rejections == b.softening_rejections);
    CHUI_REQUIRE(a.near_node_rejections == b.near_node_rejections);
    if (b.skipped_mac_evaluations) CHUI_REQUIRE(!accepted);
    else CHUI_REQUIRE(a.relative_mac_rejections == b.relative_mac_rejections);
  }
}

void testInvalidAdaptiveHistoryUsesGeometricFallback() {
  using namespace cosmosim::gravity;
  using namespace cosmosim::gravity::internal;
  TreePmOptions options;
  options.split_policy = makeTreePmSplitPolicyFromMeshSpacing(1.25, 6.25, 0.01);
  options.acceptance_policy = TreePmAcceptancePolicy::kAdaptiveRelative;
  TreeNodeAcceptanceInput input{.half_size=0.001, .com_center_offset=0.0,
      .node_mass_code=1.0, .r2=0.01, .previous_acceleration_available=true,
      .target_softening_comoving=0.0, .node_softening_min_comoving=0.0,
      .node_softening_max_comoving=0.0};
  for (const double history : {0.0, -1.0, std::numeric_limits<double>::quiet_NaN(),
                               std::numeric_limits<double>::infinity()}) {
    input.previous_acceleration_magnitude_code = history;
    for (const bool available : {false, true}) {
      input.previous_acceleration_available = available;
      TreePmTraversalCounters counters;
      const bool accepted = acceptTreePmNode(input, true, options, counters);
      CHUI_REQUIRE(accepted == acceptNodeByComDistanceMac(options.tree_options.opening_theta,
          input.half_size, input.com_center_offset, input.r2));
      CHUI_REQUIRE(counters.geometric_history_fallbacks == 1U);
      CHUI_REQUIRE(counters.relative_mac_rejections == 0U);
      input.target_inside_node = true;
      TreePmTraversalCounters inside;
      CHUI_REQUIRE(!acceptTreePmNode(input, true, options, inside));
      input.target_inside_node = false;
    }
  }
}

void testSplitKernelComplementarity() {
  const cosmosim::gravity::TreePmSplitPolicy split_policy =
      cosmosim::gravity::makeTreePmSplitPolicyFromMeshSpacing(1.6, 5.0, 0.125);

  const double radii[] = {0.01, 0.05, 0.2, 0.4, 0.8};
  for (const double radius : radii) {
    const double short_factor =
        cosmosim::gravity::treePmGaussianShortRangeForceFactor(radius, split_policy.split_scale_comoving);
    const double long_factor =
        cosmosim::gravity::treePmGaussianLongRangeForceFactor(radius, split_policy.split_scale_comoving);
    CHUI_REQUIRE(short_factor >= 0.0);
    CHUI_REQUIRE(short_factor <= 1.0 + 1.0e-12);
    CHUI_REQUIRE(long_factor >= -1.0e-12);
    CHUI_REQUIRE(long_factor <= 1.0);
    CHUI_REQUIRE(std::abs(short_factor + long_factor - 1.0) < 1.0e-12);
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
      CHUI_REQUIRE(std::abs((tree_residual + pm_long) - expected) / scale < 5.0e-13);
    }
  }
}

void testMeshCellDerivedSplitSemantics() {
  const double mesh_spacing = 0.025;
  const double asmth_cells = 1.25;
  const double rcut_cells = 4.5;
  const cosmosim::gravity::TreePmSplitPolicy split_policy =
      cosmosim::gravity::makeTreePmSplitPolicyFromMeshSpacing(asmth_cells, rcut_cells, mesh_spacing);

  CHUI_REQUIRE(std::abs(split_policy.mesh_spacing_comoving - mesh_spacing) < 1.0e-15);
  CHUI_REQUIRE(std::abs(split_policy.asmth_cells - asmth_cells) < 1.0e-15);
  CHUI_REQUIRE(std::abs(split_policy.rcut_cells - rcut_cells) < 1.0e-15);
  CHUI_REQUIRE(std::abs(split_policy.split_scale_comoving - asmth_cells * mesh_spacing) < 1.0e-15);
  CHUI_REQUIRE(std::abs(split_policy.cutoff_radius_comoving - rcut_cells * mesh_spacing) < 1.0e-15);
}

void testDiagnosticsContinuityAtSplitScale() {
  const cosmosim::gravity::TreePmSplitPolicy split_policy =
      cosmosim::gravity::makeTreePmSplitPolicyFromMeshSpacing(2.0, 6.0, 0.0625);

  const cosmosim::gravity::TreePmDiagnostics diagnostics = cosmosim::gravity::computeTreePmDiagnostics(split_policy);
  CHUI_REQUIRE(std::abs(diagnostics.mesh_spacing_comoving - split_policy.mesh_spacing_comoving) < 1.0e-15);
  CHUI_REQUIRE(std::abs(diagnostics.asmth_cells - split_policy.asmth_cells) < 1.0e-15);
  CHUI_REQUIRE(std::abs(diagnostics.rcut_cells - split_policy.rcut_cells) < 1.0e-15);
  CHUI_REQUIRE(std::abs(diagnostics.split_scale_comoving - split_policy.split_scale_comoving) < 1.0e-15);
  CHUI_REQUIRE(std::abs(diagnostics.cutoff_radius_comoving - split_policy.cutoff_radius_comoving) < 1.0e-15);
  CHUI_REQUIRE(diagnostics.short_range_factor_at_split > 0.0);
  CHUI_REQUIRE(diagnostics.long_range_factor_at_split > 0.0);
  CHUI_REQUIRE(diagnostics.short_range_factor_at_cutoff >= 0.0);
  CHUI_REQUIRE(diagnostics.short_range_factor_at_cutoff < 0.5);
  CHUI_REQUIRE(diagnostics.long_range_factor_at_cutoff > 0.5);
  CHUI_REQUIRE(diagnostics.composition_error_at_split < 1.0e-12);
  CHUI_REQUIRE(diagnostics.max_relative_composition_error < 1.0e-12);
}

void testSparseTreePmAggregateRoundPlanner() {
  using cosmosim::gravity::internal::planSparseTreePmRound;
  using cosmosim::gravity::internal::sparseTreePmPhysicalRoundCount;

  const auto zero_peer_plan = planSparseTreePmRound(4096U, 0U, 96U, 80U, 1024U);
  CHUI_REQUIRE(zero_peer_plan.targets_per_peer_per_round == 4096U);
  CHUI_REQUIRE(zero_peer_plan.aggregate_bytes_per_target == 0U);

  const auto many_peer_plan = planSparseTreePmRound(4096U, 8U, 96U, 80U, 4096U);
  CHUI_REQUIRE(many_peer_plan.aggregate_bytes_per_target == 8U * 96U);
  CHUI_REQUIRE(many_peer_plan.targets_per_peer_per_round == 5U);
  CHUI_REQUIRE(many_peer_plan.targets_per_peer_per_round *
             many_peer_plan.aggregate_bytes_per_target <= 4096U);

  const std::uint64_t logical_targets =
      static_cast<std::uint64_t>(std::numeric_limits<int>::max()) + 1000000ULL;
  const std::uint64_t physical_rounds = sparseTreePmPhysicalRoundCount(
      logical_targets, many_peer_plan.targets_per_peer_per_round);
  CHUI_REQUIRE(physical_rounds > 1U);
  CHUI_REQUIRE(physical_rounds ==
         1U + (logical_targets - 1U) /
             static_cast<std::uint64_t>(many_peer_plan.targets_per_peer_per_round));

  bool threw = false;
  try {
    (void)planSparseTreePmRound(1U, 64U, 96U, 80U, 4096U);
  } catch (const std::overflow_error&) {
    threw = true;
  }
  CHUI_REQUIRE(threw);
}

}  // namespace

int main() {
  testStableZeroAndTinyRadius();
  testCoincidentResidualLookupConsistency();
  testGaussianLookupAcrossInterval();
  testInvalidKernelArguments();
  testSquaredGeometryAndStagedMac();
  testInvalidAdaptiveHistoryUsesGeometricFallback();
  testSplitKernelComplementarity();
  testFiniteSofteningResidualComposesWithPmLongRange();
  testMeshCellDerivedSplitSemantics();
  testDiagnosticsContinuityAtSplitScale();
  testSparseTreePmAggregateRoundPlanner();
  return 0;
}
