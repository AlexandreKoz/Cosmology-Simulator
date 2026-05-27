#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

namespace cosmosim::core {

class SimulationState;
struct TransientStepWorkspace;

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
  std::string label;
  std::uint64_t owned_capacity_bytes = 0;
  std::uint64_t referenced_bytes = 0;
  std::uint64_t high_water_bytes = 0;
  bool estimate_only = false;
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

struct MemoryReport {
  std::vector<MemoryEntry> entries;
  MemoryTotals totals;
  std::vector<std::string> notes;
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
[[nodiscard]] constexpr std::uint64_t ownedCapacityBytesForContainer(const T& container) noexcept {
  return static_cast<std::uint64_t>(container.capacity()) * static_cast<std::uint64_t>(sizeof(typename T::value_type));
}

template <typename T>
[[nodiscard]] constexpr std::uint64_t referencedBytesForSpan(const T& view_span) noexcept {
  return static_cast<std::uint64_t>(view_span.size()) * static_cast<std::uint64_t>(sizeof(typename T::value_type));
}

[[nodiscard]] MemoryReport collectSimulationMemoryReport(
    const SimulationState& state,
    const TransientStepWorkspace* workspace = nullptr);

[[nodiscard]] MemoryReport estimatePreRunMemoryBudget(const MemoryBudgetEstimateInput& input);

[[nodiscard]] std::string formatMemoryReportHumanReadable(const MemoryReport& report);

}  // namespace cosmosim::core
