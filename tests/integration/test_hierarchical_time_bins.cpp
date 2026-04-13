#include <array>
#include <cassert>
#include <cstdint>
#include <vector>

#include "cosmosim/core/time_integration.hpp"

namespace {

void testDeterministicActivationOrder() {
  cosmosim::core::HierarchicalTimeBinScheduler scheduler(2);
  scheduler.reset(6, 2, 0);

  scheduler.setElementBin(0, 0, 0);
  scheduler.setElementBin(1, 0, 0);
  scheduler.setElementBin(2, 1, 0);
  scheduler.setElementBin(3, 1, 0);
  scheduler.setElementBin(4, 2, 0);
  scheduler.setElementBin(5, 2, 0);

  const std::array<std::vector<std::uint32_t>, 8> expected = {
      std::vector<std::uint32_t>{0, 1, 2, 3, 4, 5},
      std::vector<std::uint32_t>{0, 1},
      std::vector<std::uint32_t>{0, 1, 2, 3},
      std::vector<std::uint32_t>{0, 1},
      std::vector<std::uint32_t>{0, 1, 2, 3, 4, 5},
      std::vector<std::uint32_t>{0, 1},
      std::vector<std::uint32_t>{0, 1, 2, 3},
      std::vector<std::uint32_t>{0, 1},
  };

  for (std::size_t step = 0; step < expected.size(); ++step) {
    const auto active = scheduler.beginSubstep();
    assert(active.size() == expected[step].size());
    for (std::size_t i = 0; i < active.size(); ++i) {
      assert(active[i] == expected[step][i]);
    }
    scheduler.endSubstep();
  }

  const auto& diagnostics = scheduler.diagnostics();
  assert(diagnostics.occupancy_by_bin[0] == 2);
  assert(diagnostics.occupancy_by_bin[1] == 2);
  assert(diagnostics.occupancy_by_bin[2] == 2);
}

}  // namespace

int main() {
  testDeterministicActivationOrder();
  return 0;
}
