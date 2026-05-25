#include "cosmosim/core/memory_accounting.hpp"

#include <sstream>

#include "cosmosim/core/simulation_state.hpp"

namespace cosmosim::core {
namespace {

void addTotals(MemoryReport& report, const MemoryEntry& entry) {
  const std::size_t index = memorySubsystemIndex(entry.subsystem);
  if (entry.lifetime == MemoryLifetime::kPersistent) {
    report.totals.persistent_by_subsystem[index] += entry.owned_capacity_bytes;
    report.totals.persistent_total_bytes += entry.owned_capacity_bytes;
  } else if (entry.lifetime == MemoryLifetime::kTransient) {
    report.totals.transient_by_subsystem[index] += entry.owned_capacity_bytes;
    report.totals.transient_total_bytes += entry.owned_capacity_bytes;
  }
}

template <typename T>
void addOwned(MemoryReportBuilder& builder, MemorySubsystem subsystem, MemoryLifetime lifetime, std::string_view label, const T& container) {
  builder.addEntry(MemoryEntry{.subsystem = subsystem,
                               .lifetime = lifetime,
                               .label = std::string(label),
                               .owned_capacity_bytes = ownedCapacityBytesForContainer(container)});
}

}  // namespace

void MemoryReportBuilder::addEntry(MemoryEntry entry) {
  addTotals(m_report, entry);
  m_report.entries.push_back(std::move(entry));
}

MemoryReport MemoryReportBuilder::finish() && {
  for (std::size_t i = 0; i < static_cast<std::size_t>(MemorySubsystem::kCount); ++i) {
    if (m_report.totals.persistent_by_subsystem[i] == 0 && m_report.totals.transient_by_subsystem[i] == 0) {
      m_report.entries.push_back(MemoryEntry{.subsystem = static_cast<MemorySubsystem>(i),
                                             .lifetime = MemoryLifetime::kUnknown,
                                             .label = "category_present",
                                             .owned_capacity_bytes = 0});
    }
  }
  return std::move(m_report);
}

std::size_t memorySubsystemIndex(MemorySubsystem subsystem) noexcept { return static_cast<std::size_t>(subsystem); }

std::string_view memorySubsystemLabel(MemorySubsystem subsystem) noexcept {
  switch (subsystem) {
    case MemorySubsystem::kParticles:
      return "particles";
    case MemorySubsystem::kGasHydro:
      return "gas_hydro";
    case MemorySubsystem::kTree:
      return "tree";
    case MemorySubsystem::kPmMesh:
      return "pm_mesh";
    case MemorySubsystem::kSidecars:
      return "sidecars";
    case MemorySubsystem::kActiveSets:
      return "active_sets";
    case MemorySubsystem::kMpiBuffers:
      return "mpi_buffers";
    case MemorySubsystem::kScratch:
      return "scratch";
    case MemorySubsystem::kOutputBuffers:
      return "output_buffers";
    case MemorySubsystem::kCount:
      return "invalid";
  }
  return "invalid";
}
std::string_view memoryLifetimeLabel(MemoryLifetime lifetime) noexcept {
  switch (lifetime) {
    case MemoryLifetime::kPersistent:
      return "persistent";
    case MemoryLifetime::kTransient:
      return "transient";
    case MemoryLifetime::kUnknown:
      return "unknown";
  }
  return "unknown";
}

MemoryReport collectSimulationMemoryReport(const SimulationState& state, const TransientStepWorkspace* workspace) {
  MemoryReportBuilder builder;
  addOwned(builder, MemorySubsystem::kParticles, MemoryLifetime::kPersistent, "particles.position_x_comoving", state.particles.position_x_comoving);
  addOwned(builder, MemorySubsystem::kParticles, MemoryLifetime::kPersistent, "particles.position_y_comoving", state.particles.position_y_comoving);
  addOwned(builder, MemorySubsystem::kParticles, MemoryLifetime::kPersistent, "particles.position_z_comoving", state.particles.position_z_comoving);
  addOwned(builder, MemorySubsystem::kParticles, MemoryLifetime::kPersistent, "particles.velocity_x_peculiar", state.particles.velocity_x_peculiar);
  addOwned(builder, MemorySubsystem::kParticles, MemoryLifetime::kPersistent, "particles.velocity_y_peculiar", state.particles.velocity_y_peculiar);
  addOwned(builder, MemorySubsystem::kParticles, MemoryLifetime::kPersistent, "particles.velocity_z_peculiar", state.particles.velocity_z_peculiar);
  addOwned(builder, MemorySubsystem::kParticles, MemoryLifetime::kPersistent, "particles.mass_code", state.particles.mass_code);

  addOwned(builder, MemorySubsystem::kGasHydro, MemoryLifetime::kPersistent, "cells.center_x_comoving", state.cells.center_x_comoving);
  addOwned(builder, MemorySubsystem::kGasHydro, MemoryLifetime::kPersistent, "cells.mass_code", state.cells.mass_code);
  addOwned(builder, MemorySubsystem::kGasHydro, MemoryLifetime::kPersistent, "gas_cells.density_code", state.gas_cells.density_code);

  addOwned(builder, MemorySubsystem::kSidecars, MemoryLifetime::kPersistent, "particle_sidecar.particle_id", state.particle_sidecar.particle_id);
  addOwned(builder, MemorySubsystem::kSidecars, MemoryLifetime::kPersistent, "particle_sidecar.owning_rank", state.particle_sidecar.owning_rank);
  addOwned(builder, MemorySubsystem::kSidecars, MemoryLifetime::kPersistent, "star_particles.particle_index", state.star_particles.particle_index);
  addOwned(builder, MemorySubsystem::kSidecars, MemoryLifetime::kPersistent, "black_holes.particle_index", state.black_holes.particle_index);
  addOwned(builder, MemorySubsystem::kSidecars, MemoryLifetime::kPersistent, "tracers.particle_index", state.tracers.particle_index);

  builder.addEntry(MemoryEntry{.subsystem = MemorySubsystem::kTree,
                               .lifetime = MemoryLifetime::kUnknown,
                               .label = "tree.external_allocations",
                               .estimate_only = true,
                               .uncertainty_note = "Tree node storage is subsystem-internal and not yet wired into core report hooks."});
  builder.addEntry(MemoryEntry{.subsystem = MemorySubsystem::kPmMesh,
                               .lifetime = MemoryLifetime::kUnknown,
                               .label = "pm_mesh.external_allocations",
                               .estimate_only = true,
                               .uncertainty_note = "PM mesh/FFTW/GPU allocation accounting not yet exposed to central memory report."});
  builder.addEntry(MemoryEntry{.subsystem = MemorySubsystem::kMpiBuffers,
                               .lifetime = MemoryLifetime::kUnknown,
                               .label = "mpi_buffers.external_allocations",
                               .estimate_only = true,
                               .uncertainty_note = "MPI transfer buffers vary by backend and are not fully observable from core state."});
  builder.addEntry(MemoryEntry{.subsystem = MemorySubsystem::kOutputBuffers,
                               .lifetime = MemoryLifetime::kUnknown,
                               .label = "output_buffers.external_allocations",
                               .estimate_only = true,
                               .uncertainty_note = "HDF5/IO staging buffers are intentionally not treated as restart truth and are currently reported as unknown."});

  if (workspace != nullptr) {
    addOwned(builder, MemorySubsystem::kActiveSets, MemoryLifetime::kTransient, "workspace.particle_id", workspace->particle_id);
    addOwned(builder, MemorySubsystem::kActiveSets, MemoryLifetime::kTransient, "workspace.gravity_particle_index", workspace->gravity_particle_index);
    addOwned(builder, MemorySubsystem::kActiveSets, MemoryLifetime::kTransient, "workspace.cell_patch_index", workspace->cell_patch_index);
    builder.addEntry(MemoryEntry{.subsystem = MemorySubsystem::kScratch,
                                 .lifetime = MemoryLifetime::kTransient,
                                 .label = "workspace.scratch",
                                 .owned_capacity_bytes = static_cast<std::uint64_t>(workspace->scratch.capacityBytes())});
  }

  MemoryReport report = std::move(builder).finish();
  report.notes.push_back("Owned bytes use container capacity*sizeof(T); spans/views are non-owning and report zero owned bytes.");
  report.notes.push_back("Unknown/estimate entries intentionally avoid claiming allocator, MPI, GPU, FFTW, or HDF5 internal precision.");
  return report;
}

std::string formatMemoryReportHumanReadable(const MemoryReport& report) {
  std::ostringstream out;
  out << "memory_report persistent_total_bytes=" << report.totals.persistent_total_bytes
      << " transient_total_bytes=" << report.totals.transient_total_bytes << "\n";
  for (const MemoryEntry& entry : report.entries) {
    out << " - subsystem=" << memorySubsystemLabel(entry.subsystem) << " lifetime=" << memoryLifetimeLabel(entry.lifetime)
        << " label=" << entry.label << " owned_capacity_bytes=" << entry.owned_capacity_bytes;
    if (entry.referenced_bytes > 0) {
      out << " referenced_bytes=" << entry.referenced_bytes;
    }
    if (entry.estimate_only) {
      out << " estimate_only=true";
    }
    if (!entry.uncertainty_note.empty()) {
      out << " note=\"" << entry.uncertainty_note << "\"";
    }
    out << "\n";
  }
  return out.str();
}

}  // namespace cosmosim::core
