#include <cassert>
#include <stdexcept>
#include <vector>

#include "amr_hydro_validation_campaign.hpp"

namespace {

using namespace cosmosim::tests::amr_hydro_campaign;

ComparisonMetadata metadata(const char* code) {
  return ComparisonMetadata{
      .schema_version = "cosmosim_amr_profile_v1",
      .case_name = "sod_1d_temporal_amr",
      .geometry = "line",
      .source_code_name = code,
      .source_code_version = "test-fixture",
      .gamma = 1.4,
      .final_time_code = 0.2,
      .units = "code_units_v1",
      .boundary_conditions = "outflow",
      .gravity_enabled = false,
      .cooling_enabled = false,
      .cosmological_expansion_enabled = false};
}

void testSyntheticContractFixture() {
  const auto chui = metadata("CHUI");
  const auto external = metadata("synthetic_contract_fixture");
  validateComparable(chui, external);
  const std::vector<ProfilePoint> lhs{{0.0, 1.0, 0.0, 1.0}, {0.5, 0.5, 0.3, 0.4}, {1.0, 0.125, 0.0, 0.1}};
  const std::vector<ProfilePoint> rhs{{0.0, 1.0, 0.0, 1.0}, {0.5, 0.52, 0.29, 0.41}, {1.0, 0.125, 0.0, 0.1}};
  const auto norms = compareOnCommonGrid(lhs, rhs, {0.0, 0.25, 0.5, 0.75, 1.0});
  assert(norms.l1_density > 0.0);
  assert(norms.l2_density > 0.0);
}

void testRejectsIncompatiblePhysics() {
  auto lhs = metadata("CHUI");
  auto rhs = metadata("external");
  rhs.gamma = 5.0 / 3.0;
  bool threw = false;
  try {
    validateComparable(lhs, rhs);
  } catch (const std::runtime_error&) {
    threw = true;
  }
  assert(threw);
}

}  // namespace

int main() {
  testSyntheticContractFixture();
  testRejectsIncompatiblePhysics();
  return 0;
}
