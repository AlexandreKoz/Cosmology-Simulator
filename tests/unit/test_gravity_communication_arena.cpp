// 20261025 optimization campaign: authored, NOT RUN — deferred to maintainer.
#include "../support/test_require.hpp"

#include <cstdint>
#include <memory_resource>
#include <stdexcept>
#include <vector>

#include "cosmosim/core/memory_governor.hpp"
#include "cosmosim/gravity/gravity_communication_arena.hpp"

namespace {

using cosmosim::gravity::GravityCommunicationArena;
using Phase = GravityCommunicationArena::Phase;

// Growth is legal only at an idle boundary, keeps governor accounting exact,
// and never shrinks implicitly.
void testIdleBoundaryGrowthIsGoverned() {
  cosmosim::core::MemoryGovernor governor({.hard_limit_bytes = 1U << 20U});
  {
    GravityCommunicationArena arena(&governor);
    arena.configure(4096U);
    CHUI_REQUIRE(arena.capacityBytes() == 4096U && arena.growthEvents() == 0U);
    CHUI_REQUIRE(governor.snapshot().committed_bytes == 4096U);
    {
      auto lease = arena.begin(Phase::kTreeExchange);
      std::pmr::vector<std::uint8_t> bytes(lease.resource());
      bytes.resize(2048U);
      lease.recordLogicalHighWater(2048U);
      bool grew_while_leased = false;
      try {
        arena.configure(8192U);
        grew_while_leased = true;
      } catch (const std::logic_error&) {
      }
      CHUI_REQUIRE(!grew_while_leased);
    }
    arena.configure(2048U);  // smaller request: no-op, no shrink
    CHUI_REQUIRE(arena.capacityBytes() == 4096U);
    arena.configure(16384U);
    CHUI_REQUIRE(arena.capacityBytes() == 16384U && arena.growthEvents() == 1U);
    CHUI_REQUIRE(governor.snapshot().committed_bytes == 16384U);
    CHUI_REQUIRE(arena.maxLogicalHighWater() == 2048U);
    CHUI_REQUIRE(arena.unusedCapacityBytes() == 16384U - 2048U);
    {
      auto lease = arena.begin(Phase::kPmDensity);
      std::pmr::vector<std::uint8_t> bytes(lease.resource());
      bytes.resize(12000U);  // fits only the grown storage
    }
  }
  CHUI_REQUIRE(governor.snapshot().committed_bytes == 0U);
}

// A refused growth leaves the previous capacity intact and usable.
void testRefusedGrowthKeepsPreviousCapacity() {
  cosmosim::core::MemoryGovernor governor({.hard_limit_bytes = 8192U});
  GravityCommunicationArena arena(&governor);
  arena.configure(4096U);
  bool refused = false;
  try {
    arena.configure(8192U);  // old 4096 + new 8192 exceeds the 8192 hard limit
  } catch (const std::exception&) {
    refused = true;
  }
  CHUI_REQUIRE(refused);
  CHUI_REQUIRE(arena.capacityBytes() == 4096U && arena.growthEvents() == 0U);
  CHUI_REQUIRE(governor.snapshot().committed_bytes == 4096U);
  auto lease = arena.begin(Phase::kPmHalo);
  std::pmr::vector<std::uint8_t> bytes(lease.resource());
  bytes.resize(1024U);
}

}  // namespace

int main() {
  testIdleBoundaryGrowthIsGoverned();
  testRefusedGrowthKeepsPreviousCapacity();
  return 0;
}
