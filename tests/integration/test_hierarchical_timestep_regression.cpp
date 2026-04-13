#include <cassert>

#include "cosmosim/core/time_integration.hpp"

namespace {

void testMappingRegressionReference() {
  const cosmosim::core::TimeStepLimits limits{
      .min_dt_time_code = 0.0625,
      .max_dt_time_code = 1.0,
      .max_bin = 4,
  };

  const auto coarse = cosmosim::core::mapDtToTimeBin(0.999, limits);
  const auto medium = cosmosim::core::mapDtToTimeBin(0.2, limits);
  const auto fine = cosmosim::core::mapDtToTimeBin(0.0625, limits);

  assert(coarse.bin_index == 3);
  assert(medium.bin_index == 1);
  assert(fine.bin_index == 0);

  assert(cosmosim::core::binIndexToDt(coarse.bin_index, limits) == 0.5);
  assert(cosmosim::core::binIndexToDt(medium.bin_index, limits) == 0.125);
}

}  // namespace

int main() {
  testMappingRegressionReference();
  return 0;
}
