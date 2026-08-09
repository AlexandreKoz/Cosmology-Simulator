#pragma once

#include <atomic>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <memory>
#include <mutex>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/core/memory_accounting.hpp"

namespace cosmosim::core {

struct AllocatorStatsSnapshot {
  std::uint64_t alloc_calls = 0;
  std::uint64_t free_calls = 0;
  std::uint64_t bytes_allocated = 0;
  std::uint64_t bytes_freed = 0;
  std::uint64_t peak_live_bytes = 0;
  std::uint64_t live_bytes = 0;
};

class AllocatorStats {
 public:
  void recordAllocate(std::uint64_t bytes);
  void recordFree(std::uint64_t bytes);
  [[nodiscard]] AllocatorStatsSnapshot snapshot() const noexcept;
  void reset() noexcept;

 private:
  std::atomic<std::uint64_t> m_alloc_calls{0};
  std::atomic<std::uint64_t> m_free_calls{0};
  std::atomic<std::uint64_t> m_bytes_allocated{0};
  std::atomic<std::uint64_t> m_bytes_freed{0};
  std::atomic<std::uint64_t> m_peak_live_bytes{0};
  std::atomic<std::uint64_t> m_live_bytes{0};
};

class CounterRegistry {
 public:
  CounterRegistry();

  // addCount is lock-free after a thread's first registration: each worker
  // writes only to its thread-local shard. setCount/reset are synchronization-
  // boundary operations and must not race with addCount for the same registry.
  void addCount(std::string_view key, std::uint64_t delta = 1);
  void setCount(std::string_view key, std::uint64_t value);
  [[nodiscard]] std::uint64_t count(std::string_view key) const;
  [[nodiscard]] std::unordered_map<std::string, std::uint64_t> snapshotEntries() const;
  void reset();

 private:
  struct ThreadShard {
    std::unordered_map<std::string, std::uint64_t> counts;
  };

  [[nodiscard]] ThreadShard& localShard();

  std::uint64_t m_registry_id = 0;
  mutable std::mutex m_registration_mutex;
  std::vector<std::unique_ptr<ThreadShard>> m_shards;
  std::unordered_map<std::string, std::uint64_t> m_base_counts;
};

struct ProfileNode {
  std::string name;
  std::uint64_t call_count = 0;
  double inclusive_ms = 0.0;
  double exclusive_ms = 0.0;
  std::uint64_t bytes_moved = 0;
  std::unordered_map<std::string, std::size_t> child_lookup;
  std::vector<std::size_t> children;
};

enum class RuntimeEventSeverity : std::uint8_t {
  kInfo = 0,
  kWarning = 1,
  kError = 2,
  kFatal = 3,
};

struct RuntimeEvent {
  std::string event_kind;
  RuntimeEventSeverity severity = RuntimeEventSeverity::kInfo;
  std::string subsystem;
  std::optional<std::uint64_t> step_index;
  std::optional<double> simulation_time_code;
  std::optional<double> scale_factor;
  std::string message;
  std::unordered_map<std::string, std::string> payload;
};

class ProfilerSession {
 public:
  using Clock = std::chrono::steady_clock;

  explicit ProfilerSession(bool enabled = false);

  void setEnabled(bool enabled) noexcept;
  [[nodiscard]] bool enabled() const noexcept;

  void beginScope(std::string_view phase_name);
  void endScope();
  void addBytesMoved(std::uint64_t bytes);

  CounterRegistry& counters() noexcept;
  const CounterRegistry& counters() const noexcept;

  AllocatorStats& allocatorStats() noexcept;
  const AllocatorStats& allocatorStats() const noexcept;

  // Snapshot accessors merge per-thread shards and therefore require a quiescent
  // profiling boundary (for example, after an OpenMP region has joined).
  [[nodiscard]] const std::vector<ProfileNode>& nodes() const;
  [[nodiscard]] std::size_t rootNodeIndex() const noexcept;
  void recordEvent(RuntimeEvent event);
  [[nodiscard]] const std::vector<RuntimeEvent>& events() const;
  void setMemoryReport(MemoryReport report);
  [[nodiscard]] const MemoryReport* memoryReport() const noexcept;

  void reset();

 private:
  struct ActiveScope {
    std::size_t node_index = 0;
    Clock::time_point start;
    double child_elapsed_ms = 0.0;
  };

  struct SequencedEvent {
    std::uint64_t sequence = 0;
    RuntimeEvent event;
  };

  struct ThreadShard {
    std::vector<ProfileNode> nodes;
    std::vector<ActiveScope> scope_stack;
    std::vector<SequencedEvent> events;

    ThreadShard();
  };

  [[nodiscard]] ThreadShard& localShard();
  [[nodiscard]] std::size_t findOrCreateChild(
      ThreadShard& shard,
      std::size_t parent_index,
      std::string_view phase_name);
  void rebuildMergedNodes() const;
  void rebuildMergedEvents() const;

  std::uint64_t m_session_id = 0;
  std::atomic<bool> m_enabled{false};
  mutable std::mutex m_registration_mutex;
  std::vector<std::unique_ptr<ThreadShard>> m_shards;
  CounterRegistry m_counters;
  AllocatorStats m_allocator_stats;
  std::atomic<std::uint64_t> m_next_event_sequence{0};
  mutable std::vector<ProfileNode> m_merged_nodes;
  mutable std::vector<RuntimeEvent> m_merged_events;
  std::optional<MemoryReport> m_memory_report;

};

class ScopedProfile {
 public:
  ScopedProfile(ProfilerSession* session, std::string_view phase_name);
  ~ScopedProfile() noexcept;

  ScopedProfile(const ScopedProfile&) = delete;
  ScopedProfile& operator=(const ScopedProfile&) = delete;

 private:
  ProfilerSession* m_session = nullptr;
  bool m_active = false;
};

void writeProfilerReportJson(const ProfilerSession& session, const std::filesystem::path& output_path);
void writeProfilerReportCsv(const ProfilerSession& session, const std::filesystem::path& output_path);
void writeOperationalReportJson(
    const ProfilerSession& session,
    const std::filesystem::path& output_path,
    std::string_view run_label,
    std::string_view provenance_config_hash_hex);

}  // namespace cosmosim::core

#if defined(COSMOSIM_ENABLE_PROFILING) && COSMOSIM_ENABLE_PROFILING
#define COSMOSIM_DETAIL_CONCAT_IMPL(lhs, rhs) lhs##rhs
#define COSMOSIM_DETAIL_CONCAT(lhs, rhs) COSMOSIM_DETAIL_CONCAT_IMPL(lhs, rhs)
#if defined(__COUNTER__)
#define COSMOSIM_PROFILE_SCOPE(session_ptr, phase_name) \
  ::cosmosim::core::ScopedProfile COSMOSIM_DETAIL_CONCAT(cosmosim_profile_scope_, __COUNTER__)((session_ptr), (phase_name))
#else
#define COSMOSIM_PROFILE_SCOPE(session_ptr, phase_name) \
  ::cosmosim::core::ScopedProfile COSMOSIM_DETAIL_CONCAT(cosmosim_profile_scope_, __LINE__)((session_ptr), (phase_name))
#endif
#define COSMOSIM_PROFILE_BYTES(session_ptr, bytes) \
  do {                                            \
    if ((session_ptr) != nullptr) {               \
      (session_ptr)->addBytesMoved((bytes));      \
    }                                             \
  } while (false)
#else
#define COSMOSIM_PROFILE_SCOPE(session_ptr, phase_name) ((void)0)
#define COSMOSIM_PROFILE_BYTES(session_ptr, bytes) ((void)0)
#endif
