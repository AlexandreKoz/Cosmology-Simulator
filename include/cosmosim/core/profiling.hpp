#pragma once

#include <chrono>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

#include "cosmosim/core/build_config.hpp"

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
  AllocatorStatsSnapshot m_stats;
};

class CounterRegistry {
 public:
  void addCount(std::string_view key, std::uint64_t delta = 1);
  void setCount(std::string_view key, std::uint64_t value);
  [[nodiscard]] std::uint64_t count(std::string_view key) const;
  [[nodiscard]] const std::unordered_map<std::string, std::uint64_t>& entries() const noexcept;
  void reset();

 private:
  std::unordered_map<std::string, std::uint64_t> m_counts;
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

  [[nodiscard]] const std::vector<ProfileNode>& nodes() const noexcept;
  [[nodiscard]] std::size_t rootNodeIndex() const noexcept;
  void recordEvent(RuntimeEvent event);
  [[nodiscard]] const std::vector<RuntimeEvent>& events() const noexcept;

  void reset();

 private:
  struct ActiveScope {
    std::size_t node_index = 0;
    Clock::time_point start;
    double child_elapsed_ms = 0.0;
  };

  bool m_enabled = false;
  std::vector<ProfileNode> m_nodes;
  std::vector<ActiveScope> m_scope_stack;
  CounterRegistry m_counters;
  AllocatorStats m_allocator_stats;
  std::vector<RuntimeEvent> m_events;

  [[nodiscard]] std::size_t findOrCreateChild(std::size_t parent_index, std::string_view phase_name);
};

class ScopedProfile {
 public:
  ScopedProfile(ProfilerSession* session, std::string_view phase_name);
  ~ScopedProfile();

  ScopedProfile(const ScopedProfile&) = delete;
  ScopedProfile& operator=(const ScopedProfile&) = delete;

 private:
  ProfilerSession* m_session = nullptr;
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
#define COSMOSIM_PROFILE_SCOPE(session_ptr, phase_name) \
  ::cosmosim::core::ScopedProfile COSMOSIM_profile_scope_##__LINE__((session_ptr), (phase_name))
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
