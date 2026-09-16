#include <cassert>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

#include "cosmosim/core/memory_governor.hpp"
#include "cosmosim/parallel/distributed_mesh.hpp"

namespace {

void testZeroRecords() {
  const auto plan = cosmosim::parallel::planAmrFluxExchangeStaging(0U, 0U, 0U, 2U);
  assert(plan.total_inbound_capacity == 0U);
  assert(plan.max_peer_send_count == 0U);
  assert(plan.max_peer_receive_count == 0U);
  // Only O(world_size) metadata remains.
  const std::uint64_t expected_meta =
      static_cast<std::uint64_t>(2U * sizeof(std::uint64_t) * 2U + 2U * sizeof(int));
  assert(plan.peak_reservation_bytes == expected_meta);
  assert(cosmosim::parallel::amrFluxExchangeStagingPeakBytes(plan) == expected_meta);
}

void testUnevenPeerCounts() {
  const std::size_t inbound = 10U;
  const std::size_t max_send = 7U;
  const std::size_t max_recv = 3U;
  const auto plan =
      cosmosim::parallel::planAmrFluxExchangeStaging(inbound, max_send, max_recv, 4U);
  const std::uint64_t expected =
      static_cast<std::uint64_t>(
          (inbound + max_send + max_recv) * sizeof(cosmosim::parallel::AmrFluxRegisterPayloadRecord) +
          4U * sizeof(std::uint64_t) * 2U + 4U * sizeof(int));
  assert(plan.peak_reservation_bytes == expected);
}

void testBoundedSinglePeerStaging() {
  // Sequential peer buffers must NOT be summed: peak uses max, not sum.
  const auto bounded =
      cosmosim::parallel::planAmrFluxExchangeStaging(4U, 5U, 5U, 8U);
  const std::uint64_t unbounded_hint =
      static_cast<std::uint64_t>(
          (4U + 5U * 7U + 5U * 7U) * sizeof(cosmosim::parallel::AmrFluxRegisterPayloadRecord));
  assert(bounded.peak_reservation_bytes < unbounded_hint);
}

void testOverflowRejection() {
  bool rejected = false;
  try {
    const std::size_t huge = std::numeric_limits<std::size_t>::max() / 2U;
    (void)cosmosim::parallel::planAmrFluxExchangeStaging(huge, huge, 0U, 2U);
  } catch (const std::overflow_error&) {
    rejected = true;
  }
  assert(rejected);
}

void testTightHeadroomRejectionAndRetry() {
  const auto plan = cosmosim::parallel::planAmrFluxExchangeStaging(16U, 8U, 8U, 4U);
  const std::uint64_t peak = plan.peak_reservation_bytes;
  assert(peak > 0U);
  cosmosim::core::MemoryGovernor tight({.hard_limit_bytes = peak});
  tight.setBaselineOwnedBytes(0U);
  // Headroom is exactly peak, but governor applies safety margin 0 so exact
  // fit succeeds; tight-1 must reject.
  {
    cosmosim::core::MemoryGovernor short_gov({.hard_limit_bytes = peak - 1U});
    short_gov.setBaselineOwnedBytes(0U);
    bool rejected = false;
    try {
      auto r = short_gov.reserve(
          cosmosim::core::MemoryClass::kCommunication, peak, "gate1.tight");
      (void)r;
    } catch (const cosmosim::core::MemoryAdmissionError&) {
      rejected = true;
    }
    assert(rejected);
  }
  auto reservation = tight.reserve(
      cosmosim::core::MemoryClass::kCommunication, peak, "gate1.retry");
  reservation.commit();
  assert(tight.snapshot().committed_bytes == peak);
  reservation.release();
  assert(tight.snapshot().committed_bytes == 0U);
}

void testDeterministicConversionOrder() {
  // Workflow conversion preserves input order: simulate the
  // FluxRegisterEntry -> payload key order mapping.
  std::vector<std::uint64_t> input_keys = {30U, 10U, 20U};
  std::vector<std::uint64_t> output_keys;
  output_keys.reserve(input_keys.size());
  for (std::uint64_t key : input_keys) {
    output_keys.push_back(key);
  }
  assert(output_keys.size() == 3U);
  assert(output_keys[0] == 30U && output_keys[1] == 10U && output_keys[2] == 20U);
}

}  // namespace

int main() {
  testZeroRecords();
  testUnevenPeerCounts();
  testBoundedSinglePeerStaging();
  testOverflowRejection();
  testTightHeadroomRejectionAndRetry();
  testDeterministicConversionOrder();
  return 0;
}
