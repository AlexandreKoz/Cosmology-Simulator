#include <cassert>
#include <cstdint>
#include <limits>
#include <stdexcept>

#include "cosmosim/core/memory_governor.hpp"
#include "cosmosim/io/restart_checkpoint.hpp"

namespace {

void testCandidateBytesZero() {
  cosmosim::io::RestartReadCandidateDimensions dims;
  assert(cosmosim::io::restartReadCandidateStagingBytes(dims) == 0U);
}

void testCandidateBytesScalesWithCounts() {
  cosmosim::io::RestartReadCandidateDimensions small{
      .particle_count = 10U, .cell_count = 8U, .patch_count = 1U};
  cosmosim::io::RestartReadCandidateDimensions large{
      .particle_count = 20U, .cell_count = 8U, .patch_count = 1U};
  assert(
      cosmosim::io::restartReadCandidateStagingBytes(large) >
      cosmosim::io::restartReadCandidateStagingBytes(small));
}

void testOverflowRejectionBeforeAllocation() {
  cosmosim::io::RestartReadCandidateDimensions huge{
      .particle_count = std::numeric_limits<std::uint64_t>::max() / 2U};
  bool rejected = false;
  try {
    (void)cosmosim::io::restartReadCandidateStagingBytes(huge);
  } catch (const std::overflow_error&) {
    rejected = true;
  }
  assert(rejected);
}

void testTightHeadroomRejectionAndRetry() {
  cosmosim::io::RestartReadCandidateDimensions dims{
      .particle_count = 16U, .cell_count = 8U, .patch_count = 2U,
      .pending_flux_count = 4U, .temporal_cell_count = 2U};
  const std::uint64_t need = cosmosim::io::restartReadCandidateStagingBytes(dims);
  assert(need > 0U);
  cosmosim::core::MemoryGovernor tight({.hard_limit_bytes = need - 1U});
  tight.setBaselineOwnedBytes(0U);
  bool rejected = false;
  try {
    auto r = tight.reserve(
        cosmosim::core::MemoryClass::kPhaseResident, need, "restart.readback.candidate");
    (void)r;
  } catch (const cosmosim::core::MemoryAdmissionError&) {
    rejected = true;
  }
  assert(rejected);
  cosmosim::core::MemoryGovernor ample({.hard_limit_bytes = need + 1024U});
  ample.setBaselineOwnedBytes(0U);
  auto reservation = ample.reserve(
      cosmosim::core::MemoryClass::kPhaseResident, need, "restart.readback.candidate");
  reservation.commit();
  assert(ample.snapshot().committed_bytes == need);
}

}  // namespace

int main() {
  testCandidateBytesZero();
  testCandidateBytesScalesWithCounts();
  testOverflowRejectionBeforeAllocation();
  testTightHeadroomRejectionAndRetry();
  return 0;
}
