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

}  // namespace

int main() {
  testCapacityBasedAccountingUsesCapacityNotSize();
  testSpanViewReportsNoOwnedBytes();
  testAllCategoriesPresentEvenIfZero();
  testRuntimeAccountingCoversAllPersistentLanesAndWorkspaceCapacity();
  testPreRunEstimateReportsRequiredSubsystemsAndUncertainty();
  return 0;
}
