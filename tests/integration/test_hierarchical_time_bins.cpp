#include <array>
#include <cassert>
#include <cstdint>
#include <stdexcept>
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


void testDiagnosticsDoNotAlterSchedulerTruth() {
  cosmosim::core::HierarchicalTimeBinScheduler scheduler(3);
  scheduler.reset(5, 2, 4);
  scheduler.setElementBin(0, 0, scheduler.currentTick());
  scheduler.setElementBin(1, 1, scheduler.currentTick());
  scheduler.setElementBin(2, 2, scheduler.currentTick());
  scheduler.setElementBin(3, 3, scheduler.currentTick());
  scheduler.setElementBin(4, 1, scheduler.currentTick());

  const auto before = scheduler.exportPersistentState();
  const auto active = scheduler.beginSubstep();
  (void)active;
  const auto& diagnostics = scheduler.diagnostics();
  assert(diagnostics.occupancy_by_bin.size() == 4);
  assert(diagnostics.active_count_by_bin.size() == 4);
  const auto after_diagnostics = scheduler.exportPersistentState();
  assert(after_diagnostics.current_tick == before.current_tick);
  assert(after_diagnostics.max_bin == before.max_bin);
  assert(after_diagnostics.bin_index == before.bin_index);
  assert(after_diagnostics.next_activation_tick == before.next_activation_tick);
  assert(after_diagnostics.pending_bin_index == before.pending_bin_index);

  const std::uint32_t active_sum = diagnostics.active_elements;
  (void)active_sum;
  const auto after_read = scheduler.exportPersistentState();
  assert(after_read.bin_index == before.bin_index);
  assert(after_read.next_activation_tick == before.next_activation_tick);
  assert(after_read.pending_bin_index == before.pending_bin_index);
}

void testImportRejectsSchedulerSchemaInconsistency() {
  cosmosim::core::TimeBinPersistentState persistent;
  persistent.current_tick = 1;
  persistent.max_bin = 2;
  persistent.bin_index = {0, 1};
  persistent.next_activation_tick = {1};
  persistent.active_flag = {0, 0};
  persistent.pending_bin_index = {
      cosmosim::core::HierarchicalTimeBinScheduler::k_unset_pending_bin,
      cosmosim::core::HierarchicalTimeBinScheduler::k_unset_pending_bin};

  cosmosim::core::HierarchicalTimeBinScheduler scheduler(2);
  bool threw = false;
  try {
    scheduler.importPersistentState(persistent);
  } catch (const std::invalid_argument&) {
    threw = true;
  }
  assert(threw);
}

}  // namespace

int main() {
  testDeterministicActivationOrder();
  testDiagnosticsDoNotAlterSchedulerTruth();
  testImportRejectsSchedulerSchemaInconsistency();
  return 0;
}
