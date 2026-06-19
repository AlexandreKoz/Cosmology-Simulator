#include <cassert>
#include <cmath>
#include <vector>

#include "amr_hydro_validation_campaign.hpp"

namespace {

using namespace cosmosim::tests::amr_hydro_campaign;

std::vector<ProfilePoint> makeReference() {
  return {{0.0, 1.0, 0.0, 1.0}, {0.25, 0.9, 0.1, 0.8}, {0.5, 0.5, 0.4, 0.4},
          {0.75, 0.2, 0.2, 0.15}, {1.0, 0.125, 0.0, 0.1}};
}

std::vector<ProfilePoint> makeResolutionProfile(double amplitude) {
  auto points = makeReference();
  for (std::size_t i = 0; i < points.size(); ++i) {
    const double signed_error = amplitude * (static_cast<double>(i) - 2.0);
    points[i].rho_code = std::max(1.0e-6, points[i].rho_code + signed_error);
    points[i].pressure_code = std::max(1.0e-6, points[i].pressure_code + 0.5 * signed_error);
  }
  return points;
}

void testConvergenceMetricContract() {
  // This is a fast deterministic contract for the campaign metric machinery,
  // not a claim that a physical Sod/Sedov study has been executed.
  const auto reference = makeReference();
  const std::vector<double> common_grid{0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0};
  const auto coarse = compareOnCommonGrid(makeResolutionProfile(0.04), reference, common_grid);
  const auto medium = compareOnCommonGrid(makeResolutionProfile(0.02), reference, common_grid);
  const auto fine = compareOnCommonGrid(makeResolutionProfile(0.01), reference, common_grid);
  assert(coarse.l1_density > medium.l1_density && medium.l1_density > fine.l1_density);
  assert(observedOrder(coarse.l1_density, medium.l1_density, 2.0) > 0.9);
  assert(observedOrder(medium.l1_density, fine.l1_density, 2.0) > 0.9);
  assert(std::abs(strongestGradientCoordinate(reference) - 0.375) < 1.0e-12);
}

}  // namespace

int main() {
  testConvergenceMetricContract();
  return 0;
}
