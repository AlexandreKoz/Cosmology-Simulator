#include <cassert>
#include <span>
#include <vector>

#include "cosmosim/core/memory_accounting.hpp"
#include "cosmosim/core/simulation_state.hpp"

namespace {

void testCapacityBasedAccountingUsesCapacityNotSize() {
  std::vector<std::uint64_t> values;
  values.reserve(64);
  values.resize(3);
  const std::uint64_t bytes = cosmosim::core::ownedCapacityBytesForContainer(values);
  assert(bytes == static_cast<std::uint64_t>(64 * sizeof(std::uint64_t)));
}

void testSpanViewReportsNoOwnedBytes() {
  std::vector<double> values(16, 1.0);
  const std::span<const double> view(values.data(), 5);
  const std::uint64_t referenced = cosmosim::core::referencedBytesForSpan(view);
  assert(referenced == static_cast<std::uint64_t>(5 * sizeof(double)));
}

void testAllCategoriesPresentEvenIfZero() {
  cosmosim::core::SimulationState state;
  const cosmosim::core::MemoryReport report = cosmosim::core::collectSimulationMemoryReport(state, nullptr);
  bool seen[static_cast<std::size_t>(cosmosim::core::MemorySubsystem::kCount)]{};
  for (const auto& entry : report.entries) {
    seen[cosmosim::core::memorySubsystemIndex(entry.subsystem)] = true;
  }
  for (bool present : seen) {
    assert(present);
  }
}

}  // namespace

int main() {
  testCapacityBasedAccountingUsesCapacityNotSize();
  testSpanViewReportsNoOwnedBytes();
  testAllCategoriesPresentEvenIfZero();
  return 0;
}
