#pragma once

#include <cstdint>
#include <limits>
#include <stdexcept>

namespace cosmosim::workflows::internal {

// Coalesced optional science cadence. This is scheduling metadata, not a
// snapshot of historical physical state. No heap allocation or per-element
// state is retained. All counters are exact and checked.
struct OptionalDiagnosticCadence {
  std::uint64_t first_due_step = 0U;
  std::uint64_t latest_due_step = 0U;
  std::uint64_t missed_count = 0U;
  bool pending = false;

  void recordDue(std::uint64_t step) {
    if (step == 0U || (pending && step <= latest_due_step)) {
      throw std::logic_error("optional diagnostic due epochs must increase");
    }
    if (missed_count == std::numeric_limits<std::uint64_t>::max()) {
      throw std::overflow_error("optional diagnostic missed count overflows uint64");
    }
    if (!pending) { first_due_step = step; pending = true; }
    latest_due_step = step;
    ++missed_count;
  }

  void clear() noexcept {
    first_due_step = 0U;
    latest_due_step = 0U;
    missed_count = 0U;
    pending = false;
  }
};

}  // namespace cosmosim::workflows::internal
