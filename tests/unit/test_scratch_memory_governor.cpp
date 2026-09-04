#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <new>
#include <stdexcept>

#include "cosmosim/core/memory_governor.hpp"
#include "cosmosim/core/simulation_state.hpp"

namespace {
bool g_fail_next_array_allocation = false;
}

void* operator new[](std::size_t size) {
  if (g_fail_next_array_allocation) {
    g_fail_next_array_allocation = false;
    throw std::bad_alloc();
  }
  if (void* storage = std::malloc(size); storage != nullptr) {
    return storage;
  }
  throw std::bad_alloc();
}

void operator delete[](void* storage) noexcept {
  std::free(storage);
}

void operator delete[](void* storage, std::size_t) noexcept {
  std::free(storage);
}

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
  g_fail_next_array_allocation = true;
  try {
    (void)scratch.allocateBytes(2048U, alignof(std::max_align_t));
  } catch (const std::bad_alloc&) {
    allocation_failed = true;
  }
  assert(allocation_failed);
  assert(!g_fail_next_array_allocation);
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

void testPreparedGravityIndexLaneUsesGovernedScratch() {
  MemoryGovernor governor(MemoryGovernorPolicy{.hard_limit_bytes = 8192U});
  cosmosim::core::SimulationState state;
  state.resizeParticles(8U);
  cosmosim::core::TransientStepWorkspace workspace(&governor);

  workspace.prepareGravityParticleIndexScratch(state.particles.size());
  assert(!workspace.gravity_particle_index_scratch.empty());
  assert(governor.snapshot().committed_bytes == workspace.scratch.capacityBytes());
  assert(workspace.gravity_particle_index.empty());

  auto view = cosmosim::core::buildGravityParticleKernelViewAllParticlesDirect(
      state, workspace);
  assert(view.particle_index.data() == workspace.gravity_particle_index_scratch.data());
  assert(workspace.gravity_particle_index.empty());
  for (std::size_t i = 0; i < view.particle_index.size(); ++i) {
    assert(view.particle_index[i] == i);
  }

  const std::uint64_t committed_before_reset = governor.snapshot().committed_bytes;
  workspace.clear();
  assert(workspace.gravity_particle_index_scratch.empty());
  assert(governor.snapshot().committed_bytes == committed_before_reset);
  workspace.prepareGravityParticleIndexScratch(state.particles.size());
  assert(governor.snapshot().committed_bytes == committed_before_reset);
}

}  // namespace

int main() {
  testInitialCapacityRejectedBeforeCommit();
  testSuccessfulGrowthCommitsPhysicalCapacityAndResetRetainsIt();
  testMultipleBlocksAccountAndKeepPointersStable();
  testBackingAllocationFailureReturnsPendingReservation();
  testAlignmentValidationDoesNotChangeGovernorState();
  testStandaloneCompatibilityRemainsUngoverned();
  testPreparedGravityIndexLaneUsesGovernedScratch();
  return 0;
}
