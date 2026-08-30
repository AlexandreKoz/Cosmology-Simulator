#include <cassert>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <new>
#include <stdexcept>

#include "cosmosim/core/memory_governor.hpp"
#include "cosmosim/core/simulation_state.hpp"

namespace {

using cosmosim::core::MemoryClass;
using cosmosim::core::MemoryGovernor;
using cosmosim::core::MemoryGovernorPolicy;
using cosmosim::core::MonotonicScratchAllocator;

void testInitialCapacityRejectedBeforeCommit() {
  MemoryGovernor governor(MemoryGovernorPolicy{.hard_limit_bytes = 512U});
  bool rejected = false;
  try {
    MonotonicScratchAllocator scratch(&governor, 1024U);
    (void)scratch;
  } catch (const std::runtime_error&) {
    rejected = true;
  }
  assert(rejected);
  const auto snapshot = governor.snapshot();
  assert(snapshot.reserved_bytes == 0U);
  assert(snapshot.committed_bytes == 0U);
  assert(snapshot.rejection_count == 1U);
}

void testSuccessfulGrowthCommitsPhysicalCapacityAndResetRetainsIt() {
  MemoryGovernor governor(MemoryGovernorPolicy{.hard_limit_bytes = 8192U});
  {
    MonotonicScratchAllocator scratch(&governor, 1024U);
    assert(scratch.governed());
    assert(scratch.capacityBytes() == 1024U);
    assert(governor.snapshot().committed_bytes == 1024U);

    auto* first = scratch.allocateBytes(128U, 64U);
    assert(first != nullptr);
    assert(reinterpret_cast<std::uintptr_t>(first) % 64U == 0U);
    assert(scratch.usedBytes() >= 128U);
    const std::size_t logical_high_water = scratch.logicalHighWaterBytes();
    const std::size_t retained_capacity = scratch.capacityBytes();
    const std::uint64_t committed_before_reset = governor.snapshot().committed_bytes;

    scratch.reset();
    assert(scratch.usedBytes() == 0U);
    assert(scratch.logicalHighWaterBytes() == logical_high_water);
    assert(scratch.capacityBytes() == retained_capacity);
    assert(governor.snapshot().committed_bytes == committed_before_reset);

    auto* reused = scratch.allocateBytes(128U, 64U);
    assert(reused == first);
    assert(governor.snapshot().committed_bytes == committed_before_reset);
  }
  assert(governor.snapshot().committed_bytes == 0U);
}

void testMultipleBlocksAccountAndKeepPointersStable() {
  MemoryGovernor governor(MemoryGovernorPolicy{.hard_limit_bytes = 16384U});
  MonotonicScratchAllocator scratch(&governor, 1024U);
  auto* early = scratch.allocateArray<std::uint64_t>(100U);
  assert(early != nullptr);
  early[0] = UINT64_C(0x123456789abcdef0);
  (void)scratch.allocateBytes(512U, alignof(std::max_align_t));
  // First 1024-byte block cannot satisfy both requests, so growth follows the
  // allocator's geometric 2048-byte next-block policy.
  assert(scratch.capacityBytes() == 3072U);
  assert(governor.snapshot().committed_bytes == 3072U);
  assert(*early == UINT64_C(0x123456789abcdef0));
  assert(scratch.capacityHighWaterBytes() == 3072U);
}

void testBackingAllocationFailureReturnsPendingReservation() {
  MemoryGovernor governor;
  MonotonicScratchAllocator scratch(&governor);
  bool allocation_failed = false;
  try {
    (void)scratch.allocateBytes(
        std::numeric_limits<std::size_t>::max(), 1U);
  } catch (const std::bad_array_new_length&) {
    allocation_failed = true;
  } catch (const std::bad_alloc&) {
    allocation_failed = true;
  } catch (const std::length_error&) {
    // Some standard libraries reject impossible array extents before operator
    // new. This is still an allocation-path failure after reservation.
    allocation_failed = true;
  }
  assert(allocation_failed);
  const auto snapshot = governor.snapshot();
  assert(snapshot.reserved_bytes == 0U);
  assert(snapshot.committed_bytes == 0U);
}

void testAlignmentValidationDoesNotChangeGovernorState() {
  MemoryGovernor governor(MemoryGovernorPolicy{.hard_limit_bytes = 4096U});
  MonotonicScratchAllocator scratch(&governor);
  bool bad_alignment = false;
  try {
    (void)scratch.allocateBytes(16U, 3U);
  } catch (const std::invalid_argument&) {
    bad_alignment = true;
  }
  assert(bad_alignment);
  assert(governor.snapshot().committed_bytes == 0U);
  assert(governor.snapshot().reserved_bytes == 0U);
}

void testStandaloneCompatibilityRemainsUngoverned() {
  MonotonicScratchAllocator scratch(64U);
  assert(!scratch.governed());
  assert(scratch.capacityBytes() >= 64U);
  auto* value = scratch.allocateArray<double>(4U);
  assert(value != nullptr);
}

}  // namespace

int main() {
  testInitialCapacityRejectedBeforeCommit();
  testSuccessfulGrowthCommitsPhysicalCapacityAndResetRetainsIt();
  testMultipleBlocksAccountAndKeepPointersStable();
  testBackingAllocationFailureReturnsPendingReservation();
  testAlignmentValidationDoesNotChangeGovernorState();
  testStandaloneCompatibilityRemainsUngoverned();
  return 0;
}
