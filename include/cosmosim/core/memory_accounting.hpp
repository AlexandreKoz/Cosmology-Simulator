#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "cosmosim/core/memory_governor.hpp"
#include "cosmosim/core/process_memory.hpp"

namespace cosmosim::core {

class SimulationState;
struct TransientStepWorkspace;
class HierarchicalTimeBinScheduler;

enum class MemorySubsystem : std::uint8_t {
  kParticles = 0,
  kGasHydro = 1,
  kTree = 2,
  kPmMesh = 3,
  kSidecars = 4,
  kActiveSets = 5,
  kMpiBuffers = 6,
  kScratch = 7,
  kOutputBuffers = 8,
  kCount = 9,
};

enum class MemoryLifetime : std::uint8_t { kPersistent = 0, kTransient = 1, kUnknown = 2 };

struct MemoryEntry {
  MemorySubsystem subsystem = MemorySubsystem::kParticles;
  MemoryLifetime lifetime = MemoryLifetime::kUnknown;
  // Orthogonal to subsystem ownership: subsystem answers who owns a range;
  // memory_class answers what residency/lifetime purpose it serves.
  std::optional<MemoryClass> memory_class;
  std::string label;
  // Bytes currently occupied by live elements. Capacity is reported
  // separately because reusable HPC workspaces intentionally retain storage
  // between steps.
  std::uint64_t current_size_bytes = 0;
  std::uint64_t owned_capacity_bytes = 0;
  std::uint64_t referenced_bytes = 0;
  std::uint64_t high_water_bytes = 0;
  // Best estimate for the next force/integration boundary. A zero value means
  // that the owning component does not have a meaningful prediction yet.
  std::uint64_t estimated_next_step_bytes = 0;
  bool estimate_only = false;
  // True only when owned_capacity_bytes describes a byte range already
  // represented by MemoryGovernor.committed_bytes. Baseline reconciliation
  // excludes these entries so governed physical memory is counted once.
  bool governed_commitment = false;
  std::string uncertainty_note;
};

struct MemoryTotals {
  std::array<std::uint64_t, static_cast<std::size_t>(MemorySubsystem::kCount)> persistent_by_subsystem{};
  std::array<std::uint64_t, static_cast<std::size_t>(MemorySubsystem::kCount)> transient_by_subsystem{};
  std::array<std::uint64_t, static_cast<std::size_t>(MemorySubsystem::kCount)> unknown_by_subsystem{};
  std::uint64_t persistent_total_bytes = 0;
  std::uint64_t transient_total_bytes = 0;
  std::uint64_t unknown_total_bytes = 0;
};

struct DistributedMemorySummary {
  bool valid = false;
  int rank_count = 1;
  std::uint64_t local_owned_bytes = 0;
  std::uint64_t global_sum_owned_bytes = 0;
  std::uint64_t rank_max_owned_bytes = 0;
  double max_to_mean_imbalance_ratio = 1.0;
};

struct ProcessMemoryReconciliation {
  // Estimated bytes that are currently resident by CHUI contract: retained
  // baseline ownership, committed governed allocations, and the configured
  // opaque external-runtime allowance.  Future policy headroom and merely
  // reserved-but-uncommitted allocations are intentionally excluded.
  std::uint64_t current_declared_residency_bytes = 0U;
  // Backward-compatible report key.  This now aliases
  // current_declared_residency_bytes rather than governor policy demand.
  std::uint64_t known_accounted_bytes = 0U;
  std::optional<std::uint64_t> observed_rss_bytes;
  std::optional<std::uint64_t> observed_peak_rss_bytes;
  std::optional<std::uint64_t> observed_pss_bytes;
  std::optional<std::uint64_t> unexplained_resident_bytes;
  std::optional<double> observed_to_known_ratio;
};

struct DistributedMemoryMetricSummary {
  bool valid = false;
  std::uint64_t local_bytes = 0U;
  std::uint64_t global_sum_bytes = 0U;
  std::uint64_t rank_max_bytes = 0U;
  double rank_mean_bytes = 0.0;
  double max_to_mean_imbalance_ratio = 0.0;
};

struct DistributedProcessMemorySummary {
  int rank_count = 1;
  DistributedMemoryMetricSummary known_accounted{};
  DistributedMemoryMetricSummary current_rss{};
  DistributedMemoryMetricSummary peak_rss{};
  DistributedMemoryMetricSummary communication_high_water{};
};

struct MemoryReport {
  std::vector<MemoryEntry> entries;
  MemoryTotals totals;
  std::vector<std::string> notes;
  DistributedMemorySummary distributed{};
  std::optional<MemoryGovernorSnapshot> governor_snapshot;
  std::optional<ProcessMemoryObservation> process_memory_observation;
  std::optional<ProcessMemoryReconciliation> process_memory_reconciliation;
  DistributedProcessMemorySummary distributed_process_memory{};
};

struct MemoryBudgetEstimateInput {
  std::uint64_t particle_capacity = 0;
  std::uint64_t gas_cell_capacity = 0;
  std::uint64_t star_capacity = 0;
  std::uint64_t black_hole_capacity = 0;
  std::uint64_t tracer_capacity = 0;
  std::uint64_t active_particle_capacity = 0;
  std::uint64_t active_cell_capacity = 0;
  std::uint64_t tree_node_capacity = 0;
  std::uint64_t pm_grid_cells = 0;
  std::uint64_t mpi_exchange_particle_capacity = 0;
  std::uint64_t output_buffer_bytes = 0;
};

class MemoryReportBuilder {
 public:
  void addEntry(MemoryEntry entry);
  [[nodiscard]] MemoryReport finish() &&;

 private:
  MemoryReport m_report;
};

[[nodiscard]] std::size_t memorySubsystemIndex(MemorySubsystem subsystem) noexcept;
[[nodiscard]] std::string_view memorySubsystemLabel(MemorySubsystem subsystem) noexcept;
[[nodiscard]] std::string_view memoryLifetimeLabel(MemoryLifetime lifetime) noexcept;

template <typename T>
[[nodiscard]] std::uint64_t ownedCapacityBytesForContainer(const T& container) {
  const std::uint64_t count = static_cast<std::uint64_t>(container.capacity());
  constexpr std::uint64_t k_element_bytes = static_cast<std::uint64_t>(sizeof(typename T::value_type));
  if (count > std::numeric_limits<std::uint64_t>::max() / k_element_bytes) {
    throw std::overflow_error("memory-accounting container capacity byte count overflow");
  }
  return count * k_element_bytes;
}

template <typename T>
[[nodiscard]] std::uint64_t currentSizeBytesForContainer(const T& container) {
  const std::uint64_t count = static_cast<std::uint64_t>(container.size());
  constexpr std::uint64_t k_element_bytes = static_cast<std::uint64_t>(sizeof(typename T::value_type));
  if (count > std::numeric_limits<std::uint64_t>::max() / k_element_bytes) {
    throw std::overflow_error("memory-accounting container size byte count overflow");
  }
  return count * k_element_bytes;
}

template <typename T>
[[nodiscard]] std::uint64_t referencedBytesForSpan(const T& view_span) {
  const std::uint64_t count = static_cast<std::uint64_t>(view_span.size());
  constexpr std::uint64_t k_element_bytes = static_cast<std::uint64_t>(sizeof(typename T::value_type));
  if (count > std::numeric_limits<std::uint64_t>::max() / k_element_bytes) {
    throw std::overflow_error("memory-accounting span byte count overflow");
  }
  return count * k_element_bytes;
}

[[nodiscard]] MemoryReport collectSimulationMemoryReport(
    const SimulationState& state,
    const TransientStepWorkspace* workspace = nullptr);

// Authoritative scheduler storage is process-persistent runtime truth and must
// appear in the same ownership tree used to seed/reconcile MemoryGovernor.
[[nodiscard]] MemoryReport collectSchedulerMemoryReport(
    const HierarchicalTimeBinScheduler& particle_scheduler,
    const HierarchicalTimeBinScheduler& gas_cell_scheduler);

[[nodiscard]] MemoryReport mergeMemoryReports(std::span<const MemoryReport> reports);

// Sum actual CHUI-owned capacity that is not already represented by a
// governor commitment. Estimate-only and unknown/external entries are excluded.
[[nodiscard]] std::uint64_t memoryReportBaselineOwnedBytes(const MemoryReport& report);

void attachMemoryGovernorSnapshot(
    MemoryReport& report,
    const MemoryGovernor& governor);

// Attach a real process-residency sample and reconcile it against the explicit
// CHUI accounting model. The governor snapshot, when present, is the preferred
// known-demand authority because it already includes baseline ownership,
// active reservations/commitments, and configured opaque-runtime reserves.
void attachProcessMemoryObservation(
    MemoryReport& report,
    ProcessMemoryObservation observation);

[[nodiscard]] std::uint64_t memoryReportCommunicationHighWaterBytes(
    const MemoryReport& report);

[[nodiscard]] MemoryReport estimatePreRunMemoryBudget(const MemoryBudgetEstimateInput& input);

[[nodiscard]] std::string formatMemoryReportHumanReadable(const MemoryReport& report);

}  // namespace cosmosim::core
