#include <atomic>
#include <cassert>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include "cosmosim/core/memory_governor.hpp"

namespace {

using cosmosim::core::MemoryClass;
using cosmosim::core::MemoryGovernor;
using cosmosim::core::MemoryGovernorPolicy;
using cosmosim::core::MemoryPressure;

void testUnlimitedBudgetAccounting() {
  MemoryGovernor governor;
  governor.setBaselineOwnedBytes(123U);
  auto reservation = governor.reserve(MemoryClass::kScratchArena, 4096U, "unit.unlimited");
  auto pending = governor.snapshot();
  assert(pending.hard_limit_bytes == 0U);
  assert(pending.baseline_owned_bytes == 123U);
  assert(pending.reserved_bytes == 4096U);
  assert(pending.committed_bytes == 0U);
  assert(pending.accounted_bytes == 4219U);
  assert(pending.headroom_bytes == std::numeric_limits<std::uint64_t>::max());
  assert(pending.pressure == MemoryPressure::kGreen);

  reservation.commit();
  const auto committed = governor.snapshot();
  assert(committed.reserved_bytes == 0U);
  assert(committed.committed_bytes == 4096U);
  assert(committed.pressure == MemoryPressure::kGreen);
  reservation.release();
  assert(governor.snapshot().committed_bytes == 0U);
}

void testAdmissionBoundariesAndZeroByteReservation() {
  MemoryGovernor governor(MemoryGovernorPolicy{.hard_limit_bytes = 1024U});
  auto zero = governor.reserve(MemoryClass::kDiagnostic, 0U, "unit.zero");
  assert(zero.pending());
  assert(governor.snapshot().reserved_bytes == 0U);
  zero.commit();
  assert(zero.committed());
  zero.release();

  auto exact = governor.reserve(MemoryClass::kPhaseResident, 1024U, "unit.exact");
  assert(governor.snapshot().reserved_bytes == 1024U);
  exact.commit();
  assert(governor.snapshot().committed_bytes == 1024U);

  bool rejected = false;
  try {
    (void)governor.reserve(MemoryClass::kCommunication, 1U, "unit.one_over");
  } catch (const std::runtime_error& error) {
    rejected = true;
    const std::string message = error.what();
    assert(message.find("owner=unit.one_over") != std::string::npos);
    assert(message.find("class=communication") != std::string::npos);
    assert(message.find("requested_bytes=1") != std::string::npos);
    assert(message.find("hard_limit_bytes=1024") != std::string::npos);
    assert(message.find("pressure=trip") != std::string::npos);
  }
  assert(rejected);
  assert(governor.snapshot().rejection_count == 1U);
  exact.release();
}

void testSafetyMarginMatchesPreflightStyleAdmission() {
  MemoryGovernor governor(MemoryGovernorPolicy{
      .hard_limit_bytes = 1250U,
      .external_runtime_reserve_bytes = 100U,
      .planned_overlap_reserve_bytes = 100U,
      .safety_margin_fraction = 0.25,
  });
  // Raw 1000 -> ceil(1000*1.25) = 1250, exactly admissible.
  governor.setBaselineOwnedBytes(800U);
  const auto at_limit = governor.snapshot();
  assert(at_limit.accounted_bytes == 1000U);
  assert(at_limit.safety_adjusted_accounted_bytes == 1250U);
  assert(at_limit.headroom_bytes == 0U);

  bool rejected = false;
  try {
    (void)governor.reserve(MemoryClass::kScratchArena, 1U, "unit.margin_over");
  } catch (const std::runtime_error&) {
    rejected = true;
  }
  assert(rejected);
}

void testPressureBoundaries() {
  MemoryGovernor governor(MemoryGovernorPolicy{.hard_limit_bytes = 1000U});
  governor.setBaselineOwnedBytes(849U);
  assert(governor.snapshot().pressure == MemoryPressure::kGreen);
  governor.setBaselineOwnedBytes(850U);
  assert(governor.snapshot().pressure == MemoryPressure::kAmber);
  governor.setBaselineOwnedBytes(949U);
  assert(governor.snapshot().pressure == MemoryPressure::kAmber);
  governor.setBaselineOwnedBytes(950U);
  assert(governor.snapshot().pressure == MemoryPressure::kRed);
  governor.setBaselineOwnedBytes(1000U);
  assert(governor.snapshot().pressure == MemoryPressure::kRed);
  governor.setBaselineOwnedBytes(1001U);
  assert(governor.snapshot().pressure == MemoryPressure::kTrip);
  assert(governor.snapshot().headroom_bytes == 0U);
}

void testTinyFiniteBudgetPressureIsTotal() {
  MemoryGovernor one_byte(MemoryGovernorPolicy{.hard_limit_bytes = 1U});
  assert(one_byte.snapshot().pressure == MemoryPressure::kGreen);
  one_byte.setBaselineOwnedBytes(1U);
  assert(one_byte.snapshot().pressure == MemoryPressure::kRed);
  one_byte.setBaselineOwnedBytes(2U);
  assert(one_byte.snapshot().pressure == MemoryPressure::kTrip);

  MemoryGovernor two_bytes(MemoryGovernorPolicy{.hard_limit_bytes = 2U});
  assert(two_bytes.snapshot().pressure == MemoryPressure::kGreen);
  two_bytes.setBaselineOwnedBytes(1U);
  assert(two_bytes.snapshot().pressure == MemoryPressure::kAmber);
  two_bytes.setBaselineOwnedBytes(2U);
  assert(two_bytes.snapshot().pressure == MemoryPressure::kRed);

  MemoryGovernor three_bytes(MemoryGovernorPolicy{.hard_limit_bytes = 3U});
  assert(three_bytes.snapshot().pressure == MemoryPressure::kGreen);
  three_bytes.setBaselineOwnedBytes(1U);
  assert(three_bytes.snapshot().pressure == MemoryPressure::kAmber);
  three_bytes.setBaselineOwnedBytes(2U);
  assert(three_bytes.snapshot().pressure == MemoryPressure::kRed);
}

void testPendingDestructionAndCommitLifecycle() {
  MemoryGovernor governor(MemoryGovernorPolicy{.hard_limit_bytes = 4096U});
  {
    auto pending = governor.reserve(MemoryClass::kScratchArena, 256U, "unit.pending_destroy");
    assert(pending.pending());
    assert(governor.snapshot().reserved_bytes == 256U);
  }
  assert(governor.snapshot().reserved_bytes == 0U);

  {
    auto committed = governor.reserve(MemoryClass::kPersistentCache, 512U, "unit.committed_destroy");
    committed.commit();
    assert(governor.snapshot().committed_bytes == 512U);
  }
  assert(governor.snapshot().committed_bytes == 0U);
  const auto snapshot = governor.snapshot();
  assert(snapshot.peak_reserved_bytes >= 512U);
  assert(snapshot.peak_committed_bytes >= 512U);
  assert(snapshot.peak_accounted_bytes >= 512U);
}

void testMoveConstructionMoveAssignmentAndDoubleCommit() {
  MemoryGovernor governor(MemoryGovernorPolicy{.hard_limit_bytes = 4096U});
  auto first = governor.reserve(MemoryClass::kScratchArena, 300U, "unit.move_first");
  auto moved(std::move(first));
  assert(!first.valid());
  assert(moved.pending());
  moved.commit();

  bool double_commit_rejected = false;
  try {
    moved.commit();
  } catch (const std::logic_error&) {
    double_commit_rejected = true;
  }
  assert(double_commit_rejected);

  auto second = governor.reserve(MemoryClass::kCommunication, 200U, "unit.move_second");
  assert(governor.snapshot().reserved_bytes == 200U);
  second = std::move(moved);
  // Move assignment releases second's pending 200 bytes before taking the
  // committed 300-byte reservation.
  const auto after_move = governor.snapshot();
  assert(after_move.reserved_bytes == 0U);
  assert(after_move.committed_bytes == 300U);
  assert(!moved.valid());
  second.release();
  second.release();
  assert(governor.snapshot().committed_bytes == 0U);
}

void testMultipleReservationsClassesAndReleaseOrder() {
  MemoryGovernor governor(MemoryGovernorPolicy{.hard_limit_bytes = 1000U});
  auto scratch = governor.reserve(MemoryClass::kScratchArena, 300U, "unit.scratch");
  auto comm = governor.reserve(MemoryClass::kCommunication, 200U, "unit.comm");
  auto diag = governor.reserve(MemoryClass::kDiagnostic, 100U, "unit.diag");
  auto pending = governor.snapshot();
  assert(pending.reserved_bytes == 600U);
  assert(pending.reserved_by_class[cosmosim::core::memoryClassIndex(MemoryClass::kScratchArena)] == 300U);
  assert(pending.reserved_by_class[cosmosim::core::memoryClassIndex(MemoryClass::kCommunication)] == 200U);

  comm.commit();
  scratch.commit();
  diag.release();
  auto mixed = governor.snapshot();
  assert(mixed.reserved_bytes == 0U);
  assert(mixed.committed_bytes == 500U);
  assert(mixed.committed_by_class[cosmosim::core::memoryClassIndex(MemoryClass::kScratchArena)] == 300U);
  assert(mixed.committed_by_class[cosmosim::core::memoryClassIndex(MemoryClass::kCommunication)] == 200U);

  scratch.release();
  assert(governor.snapshot().committed_bytes == 200U);
  comm.release();
  assert(governor.snapshot().committed_bytes == 0U);
}

void testCheckedArithmeticAndOversizedRequest() {
  bool checked_add_threw = false;
  try {
    (void)cosmosim::core::checkedMemoryBytesAdd(
        std::numeric_limits<std::uint64_t>::max(), 1U, "unit.checked_add");
  } catch (const std::overflow_error&) {
    checked_add_threw = true;
  }
  assert(checked_add_threw);

  MemoryGovernor governor;
  governor.setBaselineOwnedBytes(1U);
  bool oversized_rejected = false;
  try {
    (void)governor.reserve(
        MemoryClass::kScratchArena,
        std::numeric_limits<std::uint64_t>::max(),
        "unit.oversized");
  } catch (const std::overflow_error& error) {
    oversized_rejected = true;
    assert(std::string(error.what()).find("uint64_accounting_overflow") != std::string::npos);
  }
  assert(oversized_rejected);
  assert(governor.snapshot().rejection_count == 1U);
  assert(governor.snapshot().reserved_bytes == 0U);

  MemoryGovernor baseline_governor(MemoryGovernorPolicy{
      .external_runtime_reserve_bytes = 2U});
  baseline_governor.setBaselineOwnedBytes(10U);
  bool baseline_overflow_rejected = false;
  try {
    baseline_governor.setBaselineOwnedBytes(
        std::numeric_limits<std::uint64_t>::max());
  } catch (const std::overflow_error&) {
    baseline_overflow_rejected = true;
  }
  assert(baseline_overflow_rejected);
  assert(baseline_governor.snapshot().baseline_owned_bytes == 10U);
}

void testConcurrentControlPlaneReservations() {
  MemoryGovernor governor(MemoryGovernorPolicy{.hard_limit_bytes = 1U << 20U});
  constexpr int k_threads = 8;
  constexpr int k_iterations = 200;
  std::atomic<int> completed{0};
  std::vector<std::thread> workers;
  workers.reserve(k_threads);
  for (int thread_index = 0; thread_index < k_threads; ++thread_index) {
    workers.emplace_back([&]() {
      for (int i = 0; i < k_iterations; ++i) {
        auto reservation = governor.reserve(
            MemoryClass::kPhaseResident, 64U, "unit.concurrent");
        reservation.commit();
      }
      completed.fetch_add(1, std::memory_order_relaxed);
    });
  }
  for (auto& worker : workers) {
    worker.join();
  }
  assert(completed.load(std::memory_order_relaxed) == k_threads);
  const auto snapshot = governor.snapshot();
  assert(snapshot.reserved_bytes == 0U);
  assert(snapshot.committed_bytes == 0U);
  assert(snapshot.rejection_count == 0U);
  assert(snapshot.peak_committed_bytes > 0U);
}

}  // namespace

int main() {
  testUnlimitedBudgetAccounting();
  testAdmissionBoundariesAndZeroByteReservation();
  testSafetyMarginMatchesPreflightStyleAdmission();
  testPressureBoundaries();
  testTinyFiniteBudgetPressureIsTotal();
  testPendingDestructionAndCommitLifecycle();
  testMoveConstructionMoveAssignmentAndDoubleCommit();
  testMultipleReservationsClassesAndReleaseOrder();
  testCheckedArithmeticAndOversizedRequest();
  testConcurrentControlPlaneReservations();
  return 0;
}
