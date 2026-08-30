#include <cassert>
#include <cstdint>
#include <limits>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include "cosmosim/core/device_buffer.hpp"
#include "cosmosim/core/memory_accounting.hpp"
#include "cosmosim/core/memory_governor.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/core/time_scheduler.hpp"

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


void testRuntimeAccountingCoversAllPersistentLanesAndWorkspaceCapacity() {
  cosmosim::core::SimulationState state;
  state.resizeParticles(4);
  state.resizeCells(3);
  state.star_particles.resize(1);
  state.black_holes.resize(1);
  state.tracers.resize(1);
  state.rebuildSpeciesIndex();

  cosmosim::core::TransientStepWorkspace workspace;
  workspace.particle_position_x_comoving.reserve(9);
  workspace.hydro_cell_index.reserve(7);
  workspace.hydro_recon_gradient_x.reserve(11);
  static_cast<void>(workspace.scratch.allocateBytes(128, alignof(double)));

  const cosmosim::core::MemoryReport report = cosmosim::core::collectSimulationMemoryReport(state, &workspace);
  assert(report.totals.persistent_by_subsystem[cosmosim::core::memorySubsystemIndex(cosmosim::core::MemorySubsystem::kParticles)] > 0);
  assert(report.totals.persistent_by_subsystem[cosmosim::core::memorySubsystemIndex(cosmosim::core::MemorySubsystem::kGasHydro)] > 0);
  assert(report.totals.persistent_by_subsystem[cosmosim::core::memorySubsystemIndex(cosmosim::core::MemorySubsystem::kSidecars)] > 0);
  assert(report.totals.transient_by_subsystem[cosmosim::core::memorySubsystemIndex(cosmosim::core::MemorySubsystem::kActiveSets)] > 0);
  assert(report.totals.transient_by_subsystem[cosmosim::core::memorySubsystemIndex(cosmosim::core::MemorySubsystem::kScratch)] >= 128);

  bool saw_velocity_z = false;
  bool saw_sound_speed = false;
  bool saw_hydro_gradient_scratch = false;
  for (const auto& entry : report.entries) {
    if (entry.label == "particles.velocity_z_peculiar") { saw_velocity_z = true; }
    if (entry.label == "gas_cells.sound_speed_code") { saw_sound_speed = true; }
    if (entry.label == "workspace.hydro_recon_gradient_x") { saw_hydro_gradient_scratch = true; }
  }
  assert(saw_velocity_z);
  assert(saw_sound_speed);
  assert(saw_hydro_gradient_scratch);
}

void testDeviceBufferPreservesHistoricalHighWaterAcrossShrink() {
  cosmosim::core::DeviceBufferDouble buffer;
  buffer.resize(64U);
  assert(buffer.size() == 64U);
  assert(buffer.highWaterSize() == 64U);
  assert(buffer.highWaterBytes() == 64U * sizeof(double));

  buffer.resize(8U);
  assert(buffer.size() == 8U);
  assert(buffer.sizeBytes() == 8U * sizeof(double));
  assert(buffer.highWaterSize() == 64U);
  assert(buffer.highWaterBytes() == 64U * sizeof(double));

  cosmosim::core::DeviceBufferDouble moved(std::move(buffer));
  assert(moved.size() == 8U);
  assert(moved.highWaterSize() == 64U);
}

void testPreRunEstimateReportsRequiredSubsystemsAndUncertainty() {
  cosmosim::core::MemoryBudgetEstimateInput input;
  input.particle_capacity = 100;
  input.gas_cell_capacity = 80;
  input.star_capacity = 10;
  input.black_hole_capacity = 2;
  input.tracer_capacity = 5;
  input.active_particle_capacity = 25;
  input.active_cell_capacity = 20;
  input.tree_node_capacity = 160;
  input.pm_grid_cells = 64;
  input.mpi_exchange_particle_capacity = 12;
  input.output_buffer_bytes = 4096;

  const cosmosim::core::MemoryReport estimate = cosmosim::core::estimatePreRunMemoryBudget(input);
  assert(estimate.totals.persistent_total_bytes > 0);
  assert(estimate.totals.transient_total_bytes > 0);
  bool saw_pm = false;
  bool saw_output = false;
  bool saw_estimate_note = false;
  for (const auto& entry : estimate.entries) {
    if (entry.subsystem == cosmosim::core::MemorySubsystem::kPmMesh && entry.estimate_only) { saw_pm = true; }
    if (entry.subsystem == cosmosim::core::MemorySubsystem::kOutputBuffers && entry.owned_capacity_bytes == 4096) { saw_output = true; }
    if (entry.estimate_only && !entry.uncertainty_note.empty()) { saw_estimate_note = true; }
  }
  assert(saw_pm);
  assert(saw_output);
  assert(saw_estimate_note);
}

void testMemorySubsystemAndClassCoexistAndSerialize() {
  cosmosim::core::MemoryReportBuilder builder;
  builder.addEntry(cosmosim::core::MemoryEntry{
      .subsystem = cosmosim::core::MemorySubsystem::kMpiBuffers,
      .lifetime = cosmosim::core::MemoryLifetime::kTransient,
      .memory_class = cosmosim::core::MemoryClass::kCommunication,
      .label = "unit.communication",
      .current_size_bytes = 32U,
      .owned_capacity_bytes = 64U,
      .high_water_bytes = 64U,
  });
  cosmosim::core::MemoryReport report = std::move(builder).finish();
  const auto& entry = report.entries.front();
  assert(entry.subsystem == cosmosim::core::MemorySubsystem::kMpiBuffers);
  assert(entry.memory_class.has_value());
  assert(*entry.memory_class == cosmosim::core::MemoryClass::kCommunication);
  const std::string formatted = cosmosim::core::formatMemoryReportHumanReadable(report);
  assert(formatted.find("subsystem=mpi_buffers") != std::string::npos);
  assert(formatted.find("class=communication") != std::string::npos);
}

void testGovernorReconciliationExcludesGovernedScratchCommitment() {
  cosmosim::core::MemoryGovernor governor(
      cosmosim::core::MemoryGovernorPolicy{.hard_limit_bytes = 1U << 20U});
  cosmosim::core::SimulationState state;
  state.resizeParticles(2U);
  cosmosim::core::TransientStepWorkspace workspace(&governor);
  static_cast<void>(workspace.scratch.allocateBytes(128U, alignof(double)));

  cosmosim::core::MemoryReport report =
      cosmosim::core::collectSimulationMemoryReport(state, &workspace);
  bool saw_governed_scratch = false;
  for (const auto& entry : report.entries) {
    if (entry.label == "workspace.scratch_arena") {
      saw_governed_scratch = true;
      assert(entry.governed_commitment);
      assert(entry.memory_class == cosmosim::core::MemoryClass::kScratchArena);
      assert(entry.owned_capacity_bytes == workspace.scratch.capacityBytes());
      assert(entry.current_size_bytes == workspace.scratch.usedBytes());
    }
  }
  assert(saw_governed_scratch);

  const std::uint64_t baseline =
      cosmosim::core::memoryReportBaselineOwnedBytes(report);
  governor.setBaselineOwnedBytes(baseline);
  cosmosim::core::attachMemoryGovernorSnapshot(report, governor);
  assert(report.governor_snapshot.has_value());
  const auto snapshot = *report.governor_snapshot;
  assert(snapshot.baseline_owned_bytes == baseline);
  assert(snapshot.committed_bytes == workspace.scratch.capacityBytes());
  assert(snapshot.accounted_bytes ==
         baseline + static_cast<std::uint64_t>(workspace.scratch.capacityBytes()));

  const std::string formatted = cosmosim::core::formatMemoryReportHumanReadable(report);
  assert(formatted.find("governor hard_limit_bytes=") != std::string::npos);
  assert(formatted.find("governed_commitment=true") != std::string::npos);
}

void testUnlimitedGovernorReporting() {
  cosmosim::core::MemoryGovernor governor;
  cosmosim::core::MemoryReport report;
  cosmosim::core::attachMemoryGovernorSnapshot(report, governor);
  assert(report.governor_snapshot->hard_limit_bytes == 0U);
  assert(report.governor_snapshot->pressure == cosmosim::core::MemoryPressure::kGreen);
  assert(report.governor_snapshot->headroom_bytes ==
         std::numeric_limits<std::uint64_t>::max());
}

void testSchedulerOwnershipParticipatesInBaselineReport() {
  cosmosim::core::HierarchicalTimeBinScheduler particle_scheduler(2U);
  cosmosim::core::HierarchicalTimeBinScheduler gas_scheduler(2U);
  particle_scheduler.reset(16U, 0U);
  gas_scheduler.reset(7U, 0U);

  const cosmosim::core::MemoryReport report =
      cosmosim::core::collectSchedulerMemoryReport(
          particle_scheduler, gas_scheduler);
  bool saw_particle = false;
  bool saw_gas = false;
  for (const auto& entry : report.entries) {
    if (entry.label == "scheduler.particle_owned_state") {
      saw_particle = true;
      assert(entry.owned_capacity_bytes == particle_scheduler.ownedCapacityBytes());
      assert(entry.memory_class == cosmosim::core::MemoryClass::kCanonicalPersistent);
    }
    if (entry.label == "scheduler.gas_cell_owned_state") {
      saw_gas = true;
      assert(entry.owned_capacity_bytes == gas_scheduler.ownedCapacityBytes());
      assert(entry.memory_class == cosmosim::core::MemoryClass::kCanonicalPersistent);
    }
  }
  assert(saw_particle && saw_gas);
  assert(cosmosim::core::memoryReportBaselineOwnedBytes(report) ==
         particle_scheduler.ownedCapacityBytes() + gas_scheduler.ownedCapacityBytes());
}

}  // namespace

int main() {
  testCapacityBasedAccountingUsesCapacityNotSize();
  testSpanViewReportsNoOwnedBytes();
  testAllCategoriesPresentEvenIfZero();
  testRuntimeAccountingCoversAllPersistentLanesAndWorkspaceCapacity();
  testDeviceBufferPreservesHistoricalHighWaterAcrossShrink();
  testPreRunEstimateReportsRequiredSubsystemsAndUncertainty();
  testMemorySubsystemAndClassCoexistAndSerialize();
  testGovernorReconciliationExcludesGovernedScratchCommitment();
  testUnlimitedGovernorReporting();
  testSchedulerOwnershipParticipatesInBaselineReport();
  return 0;
}
