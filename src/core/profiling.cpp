#include "cosmosim/core/profiling.hpp"

#include <algorithm>
#include <atomic>
#include <fstream>
#include <iomanip>
#include <optional>
#include <sstream>
#include <stdexcept>

namespace cosmosim::core {
namespace {

[[nodiscard]] std::string escapeJson(const std::string& input) {
  std::string out;
  out.reserve(input.size());
  for (const char c : input) {
    switch (c) {
      case '"':
        out += "\\\"";
        break;
      case '\\':
        out += "\\\\";
        break;
      case '\n':
        out += "\\n";
        break;
      default:
        out.push_back(c);
        break;
    }
  }
  return out;
}

[[nodiscard]] const char* runtimeEventSeverityLabel(RuntimeEventSeverity severity) {
  switch (severity) {
    case RuntimeEventSeverity::kInfo:
      return "info";
    case RuntimeEventSeverity::kWarning:
      return "warning";
    case RuntimeEventSeverity::kError:
      return "error";
    case RuntimeEventSeverity::kFatal:
      return "fatal";
  }
  return "unknown";
}

void writeOptionalU64(std::ostream& out, const std::optional<std::uint64_t>& value) {
  if (value.has_value()) {
    out << *value;
    return;
  }
  out << "null";
}

void writeOptionalDouble(std::ostream& out, const std::optional<double>& value) {
  if (value.has_value()) {
    out << std::fixed << std::setprecision(6) << *value;
    return;
  }
  out << "null";
}

void writeEventJson(std::ostream& out, const RuntimeEvent& event, int indent_level) {
  const std::string indent(static_cast<std::size_t>(indent_level), ' ');
  const std::string field_indent(static_cast<std::size_t>(indent_level + 2), ' ');
  out << indent << "{\n";
  out << field_indent << "\"event_kind\": \"" << escapeJson(event.event_kind) << "\",\n";
  out << field_indent << "\"severity\": \"" << runtimeEventSeverityLabel(event.severity) << "\",\n";
  out << field_indent << "\"subsystem\": \"" << escapeJson(event.subsystem) << "\",\n";
  out << field_indent << "\"step_index\": ";
  writeOptionalU64(out, event.step_index);
  out << ",\n";
  out << field_indent << "\"simulation_time_code\": ";
  writeOptionalDouble(out, event.simulation_time_code);
  out << ",\n";
  out << field_indent << "\"scale_factor\": ";
  writeOptionalDouble(out, event.scale_factor);
  out << ",\n";
  out << field_indent << "\"message\": \"" << escapeJson(event.message) << "\",\n";
  out << field_indent << "\"payload\": {";
  bool first_payload = true;
  for (const auto& [key, value] : event.payload) {
    out << (first_payload ? "\n" : ",\n");
    out << field_indent << "  \"" << escapeJson(key) << "\": \"" << escapeJson(value) << "\"";
    first_payload = false;
  }
  if (!first_payload) {
    out << "\n" << field_indent;
  }
  out << "}\n";
  out << indent << "}";
}

void writeMemoryReportJson(std::ostream& out, const MemoryReport& report, int indent_level) {
  const std::string indent(static_cast<std::size_t>(indent_level), ' ');
  const std::string field_indent(static_cast<std::size_t>(indent_level + 2), ' ');
  out << indent << "{\n";
  out << field_indent << "\"persistent_total_bytes\": " << report.totals.persistent_total_bytes << ",\n";
  out << field_indent << "\"transient_total_bytes\": " << report.totals.transient_total_bytes << ",\n";
  out << field_indent << "\"unknown_external_allocations\": true\n";
  out << indent << "}";
}

void writeNodeJson(std::ostream& out, const std::vector<ProfileNode>& nodes, std::size_t node_index, int indent_level) {
  const ProfileNode& node = nodes.at(node_index);
  const std::string indent(static_cast<std::size_t>(indent_level), ' ');
  const std::string child_indent(static_cast<std::size_t>(indent_level + 2), ' ');

  out << indent << "{\n";
  out << child_indent << "\"name\": \"" << escapeJson(node.name) << "\",\n";
  out << child_indent << "\"call_count\": " << node.call_count << ",\n";
  out << child_indent << "\"inclusive_ms\": " << std::fixed << std::setprecision(6) << node.inclusive_ms << ",\n";
  out << child_indent << "\"exclusive_ms\": " << std::fixed << std::setprecision(6) << node.exclusive_ms << ",\n";
  out << child_indent << "\"bytes_moved\": " << node.bytes_moved << ",\n";
  out << child_indent << "\"children\": [";
  if (!node.children.empty()) {
    out << "\n";
    for (std::size_t i = 0; i < node.children.size(); ++i) {
      writeNodeJson(out, nodes, node.children[i], indent_level + 4);
      if (i + 1 < node.children.size()) {
        out << ",";
      }
      out << "\n";
    }
    out << child_indent;
  }
  out << "]\n";
  out << indent << "}";
}

}  // namespace

namespace {
std::atomic<std::uint64_t> g_next_counter_registry_id{1};
std::atomic<std::uint64_t> g_next_profiler_session_id{1};
}

void AllocatorStats::recordAllocate(std::uint64_t bytes) {
  m_alloc_calls.fetch_add(1, std::memory_order_relaxed);
  m_bytes_allocated.fetch_add(bytes, std::memory_order_relaxed);
  const std::uint64_t live = m_live_bytes.fetch_add(bytes, std::memory_order_relaxed) + bytes;
  std::uint64_t peak = m_peak_live_bytes.load(std::memory_order_relaxed);
  while (peak < live &&
         !m_peak_live_bytes.compare_exchange_weak(
             peak, live, std::memory_order_relaxed, std::memory_order_relaxed)) {
  }
}

void AllocatorStats::recordFree(std::uint64_t bytes) {
  m_free_calls.fetch_add(1, std::memory_order_relaxed);
  m_bytes_freed.fetch_add(bytes, std::memory_order_relaxed);
  std::uint64_t live = m_live_bytes.load(std::memory_order_relaxed);
  while (true) {
    const std::uint64_t replacement = bytes > live ? 0 : live - bytes;
    if (m_live_bytes.compare_exchange_weak(
            live, replacement, std::memory_order_relaxed, std::memory_order_relaxed)) {
      break;
    }
  }
}

AllocatorStatsSnapshot AllocatorStats::snapshot() const noexcept {
  return AllocatorStatsSnapshot{
      .alloc_calls = m_alloc_calls.load(std::memory_order_relaxed),
      .free_calls = m_free_calls.load(std::memory_order_relaxed),
      .bytes_allocated = m_bytes_allocated.load(std::memory_order_relaxed),
      .bytes_freed = m_bytes_freed.load(std::memory_order_relaxed),
      .peak_live_bytes = m_peak_live_bytes.load(std::memory_order_relaxed),
      .live_bytes = m_live_bytes.load(std::memory_order_relaxed),
  };
}

void AllocatorStats::reset() noexcept {
  m_alloc_calls.store(0, std::memory_order_relaxed);
  m_free_calls.store(0, std::memory_order_relaxed);
  m_bytes_allocated.store(0, std::memory_order_relaxed);
  m_bytes_freed.store(0, std::memory_order_relaxed);
  m_peak_live_bytes.store(0, std::memory_order_relaxed);
  m_live_bytes.store(0, std::memory_order_relaxed);
}

CounterRegistry::CounterRegistry()
    : m_registry_id(g_next_counter_registry_id.fetch_add(1, std::memory_order_relaxed)) {}

CounterRegistry::ThreadShard& CounterRegistry::localShard() {
  struct TlsBinding {
    std::uint64_t registry_id;
    ThreadShard* shard;
  };
  thread_local std::vector<TlsBinding> bindings;
  for (const TlsBinding& binding : bindings) {
    if (binding.registry_id == m_registry_id) {
      return *binding.shard;
    }
  }
  auto shard = std::make_unique<ThreadShard>();
  ThreadShard* shard_ptr = shard.get();
  {
    std::lock_guard lock(m_registration_mutex);
    m_shards.push_back(std::move(shard));
  }
  bindings.push_back(TlsBinding{m_registry_id, shard_ptr});
  return *shard_ptr;
}

void CounterRegistry::addCount(std::string_view key, std::uint64_t delta) {
  localShard().counts[std::string(key)] += delta;
}

void CounterRegistry::setCount(std::string_view key, std::uint64_t value) {
  std::lock_guard lock(m_registration_mutex);
  const std::string owned_key(key);
  m_base_counts[owned_key] = value;
  for (auto& shard : m_shards) {
    shard->counts.erase(owned_key);
  }
}

std::uint64_t CounterRegistry::count(std::string_view key) const {
  const auto snapshot = snapshotEntries();
  const auto it = snapshot.find(std::string(key));
  return it == snapshot.end() ? 0 : it->second;
}

std::unordered_map<std::string, std::uint64_t> CounterRegistry::snapshotEntries() const {
  // Called at reporting/synchronization boundaries after worker threads join.
  std::lock_guard lock(m_registration_mutex);
  auto merged = m_base_counts;
  for (const auto& shard : m_shards) {
    for (const auto& [key, value] : shard->counts) {
      merged[key] += value;
    }
  }
  return merged;
}

void CounterRegistry::reset() {
  std::lock_guard lock(m_registration_mutex);
  m_base_counts.clear();
  for (auto& shard : m_shards) {
    shard->counts.clear();
  }
}

ProfilerSession::ThreadShard::ThreadShard() {
  nodes.push_back(ProfileNode{.name = "root", .call_count = 0});
}

ProfilerSession::ProfilerSession(bool enabled)
    : m_session_id(g_next_profiler_session_id.fetch_add(1, std::memory_order_relaxed)),
      m_enabled(enabled) {
  m_merged_nodes.push_back(ProfileNode{.name = "root", .call_count = 0});
}

ProfilerSession::ThreadShard& ProfilerSession::localShard() {
  struct TlsBinding {
    std::uint64_t session_id;
    ThreadShard* shard;
  };
  thread_local std::vector<TlsBinding> bindings;
  for (const TlsBinding& binding : bindings) {
    if (binding.session_id == m_session_id) {
      return *binding.shard;
    }
  }
  auto shard = std::make_unique<ThreadShard>();
  ThreadShard* shard_ptr = shard.get();
  {
    std::lock_guard lock(m_registration_mutex);
    m_shards.push_back(std::move(shard));
  }
  bindings.push_back(TlsBinding{m_session_id, shard_ptr});
  return *shard_ptr;
}

void ProfilerSession::setEnabled(bool enabled) noexcept {
  m_enabled.store(enabled, std::memory_order_relaxed);
}

bool ProfilerSession::enabled() const noexcept {
  return m_enabled.load(std::memory_order_relaxed);
}

void ProfilerSession::beginScope(std::string_view phase_name) {
  if (!enabled()) {
    return;
  }
  ThreadShard& shard = localShard();
  const std::size_t parent_index =
      shard.scope_stack.empty() ? rootNodeIndex() : shard.scope_stack.back().node_index;
  const std::size_t node_index = findOrCreateChild(shard, parent_index, phase_name);
  ++shard.nodes[node_index].call_count;
  shard.scope_stack.push_back(ActiveScope{
      .node_index = node_index,
      .start = Clock::now(),
      .child_elapsed_ms = 0.0,
  });
}

void ProfilerSession::endScope() {
  if (!enabled()) {
    return;
  }
  ThreadShard& shard = localShard();
  if (shard.scope_stack.empty()) {
    throw std::logic_error("ProfilerSession::endScope called with empty thread-local scope stack");
  }
  ActiveScope scope = shard.scope_stack.back();
  shard.scope_stack.pop_back();
  const double elapsed_ms = std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(
      Clock::now() - scope.start).count();
  ProfileNode& node = shard.nodes[scope.node_index];
  node.inclusive_ms += elapsed_ms;
  node.exclusive_ms += std::max(0.0, elapsed_ms - scope.child_elapsed_ms);
  if (!shard.scope_stack.empty()) {
    shard.scope_stack.back().child_elapsed_ms += elapsed_ms;
  }
}

void ProfilerSession::addBytesMoved(std::uint64_t bytes) {
  if (!enabled()) {
    return;
  }
  ThreadShard& shard = localShard();
  const std::size_t node_index =
      shard.scope_stack.empty() ? rootNodeIndex() : shard.scope_stack.back().node_index;
  shard.nodes[node_index].bytes_moved += bytes;
}

CounterRegistry& ProfilerSession::counters() noexcept { return m_counters; }
const CounterRegistry& ProfilerSession::counters() const noexcept { return m_counters; }
AllocatorStats& ProfilerSession::allocatorStats() noexcept { return m_allocator_stats; }
const AllocatorStats& ProfilerSession::allocatorStats() const noexcept { return m_allocator_stats; }

std::size_t ProfilerSession::findOrCreateChild(
    ThreadShard& shard,
    std::size_t parent_index,
    std::string_view phase_name) {
  const std::string phase_key(phase_name);
  const auto it = shard.nodes[parent_index].child_lookup.find(phase_key);
  if (it != shard.nodes[parent_index].child_lookup.end()) {
    return it->second;
  }
  const std::size_t child_index = shard.nodes.size();
  shard.nodes.push_back(ProfileNode{.name = phase_key, .call_count = 0});
  shard.nodes[parent_index].children.push_back(child_index);
  shard.nodes[parent_index].child_lookup.emplace(phase_key, child_index);
  return child_index;
}

void ProfilerSession::rebuildMergedNodes() const {
  std::lock_guard lock(m_registration_mutex);
  m_merged_nodes.clear();
  m_merged_nodes.push_back(ProfileNode{.name = "root", .call_count = 0});

  const auto merge_node = [&](auto&& self,
                              const std::vector<ProfileNode>& source,
                              std::size_t source_index,
                              std::size_t target_index) -> void {
    const ProfileNode& source_node = source[source_index];
    ProfileNode& target = m_merged_nodes[target_index];
    target.call_count += source_node.call_count;
    target.inclusive_ms += source_node.inclusive_ms;
    target.exclusive_ms += source_node.exclusive_ms;
    target.bytes_moved += source_node.bytes_moved;
    for (const std::size_t source_child : source_node.children) {
      const std::string& name = source[source_child].name;
      std::size_t target_child = 0;
      const auto existing = m_merged_nodes[target_index].child_lookup.find(name);
      if (existing == m_merged_nodes[target_index].child_lookup.end()) {
        target_child = m_merged_nodes.size();
        m_merged_nodes.push_back(ProfileNode{.name = name, .call_count = 0});
        m_merged_nodes[target_index].children.push_back(target_child);
        m_merged_nodes[target_index].child_lookup.emplace(name, target_child);
      } else {
        target_child = existing->second;
      }
      self(self, source, source_child, target_child);
    }
  };

  for (const auto& shard : m_shards) {
    merge_node(merge_node, shard->nodes, rootNodeIndex(), rootNodeIndex());
  }
}

const std::vector<ProfileNode>& ProfilerSession::nodes() const {
  rebuildMergedNodes();
  return m_merged_nodes;
}

std::size_t ProfilerSession::rootNodeIndex() const noexcept { return 0; }

void ProfilerSession::recordEvent(RuntimeEvent event) {
  ThreadShard& shard = localShard();
  shard.events.push_back(SequencedEvent{
      .sequence = m_next_event_sequence.fetch_add(1, std::memory_order_relaxed),
      .event = std::move(event),
  });
}

void ProfilerSession::rebuildMergedEvents() const {
  std::lock_guard lock(m_registration_mutex);
  std::vector<const SequencedEvent*> ordered;
  for (const auto& shard : m_shards) {
    for (const SequencedEvent& event : shard->events) {
      ordered.push_back(&event);
    }
  }
  std::sort(ordered.begin(), ordered.end(), [](const SequencedEvent* lhs, const SequencedEvent* rhs) {
    return lhs->sequence < rhs->sequence;
  });
  m_merged_events.clear();
  m_merged_events.reserve(ordered.size());
  for (const SequencedEvent* event : ordered) {
    m_merged_events.push_back(event->event);
  }
}

const std::vector<RuntimeEvent>& ProfilerSession::events() const {
  rebuildMergedEvents();
  return m_merged_events;
}

void ProfilerSession::setMemoryReport(MemoryReport report) {
  m_memory_report = std::move(report);
}

const MemoryReport* ProfilerSession::memoryReport() const noexcept {
  return m_memory_report.has_value() ? &(*m_memory_report) : nullptr;
}

void ProfilerSession::reset() {
  std::lock_guard lock(m_registration_mutex);
  for (auto& shard : m_shards) {
    shard->nodes.clear();
    shard->nodes.push_back(ProfileNode{.name = "root", .call_count = 0});
    shard->scope_stack.clear();
    shard->events.clear();
  }
  m_counters.reset();
  m_allocator_stats.reset();
  m_next_event_sequence.store(0, std::memory_order_relaxed);
  m_merged_nodes.clear();
  m_merged_nodes.push_back(ProfileNode{.name = "root", .call_count = 0});
  m_merged_events.clear();
  m_memory_report.reset();
}

ScopedProfile::ScopedProfile(ProfilerSession* session, std::string_view phase_name)
    : m_session(session) {
  if (m_session != nullptr) {
    m_session->beginScope(phase_name);
  }
}

ScopedProfile::~ScopedProfile() {
  if (m_session != nullptr) {
    m_session->endScope();
  }
}

void writeProfilerReportJson(const ProfilerSession& session, const std::filesystem::path& output_path) {
  std::ofstream out(output_path);
  if (!out) {
    throw std::runtime_error("failed to open profiler JSON output path");
  }

  out << "{\n";
  out << "  \"schema_version\": 1,\n";
  out << "  \"enabled\": " << (session.enabled() ? "true" : "false") << ",\n";

  const AllocatorStatsSnapshot allocator = session.allocatorStats().snapshot();
  out << "  \"allocator\": {\n";
  out << "    \"alloc_calls\": " << allocator.alloc_calls << ",\n";
  out << "    \"free_calls\": " << allocator.free_calls << ",\n";
  out << "    \"bytes_allocated\": " << allocator.bytes_allocated << ",\n";
  out << "    \"bytes_freed\": " << allocator.bytes_freed << ",\n";
  out << "    \"peak_live_bytes\": " << allocator.peak_live_bytes << ",\n";
  out << "    \"live_bytes\": " << allocator.live_bytes << "\n";
  out << "  },\n";
  if (const MemoryReport* report = session.memoryReport(); report != nullptr) {
    out << "  \"memory_report\": ";
    writeMemoryReportJson(out, *report, 2);
    out << ",\n";
  }

  out << "  \"counters\": {";
  bool first = true;
  for (const auto& [key, value] : session.counters().snapshotEntries()) {
    out << (first ? "\n" : ",\n");
    out << "    \"" << escapeJson(key) << "\": " << value;
    first = false;
  }
  if (!first) {
    out << "\n  ";
  }
  out << "},\n";

  out << "  \"phases\": ";
  writeNodeJson(out, session.nodes(), session.rootNodeIndex(), 2);
  out << "\n}\n";
}

void writeProfilerReportCsv(const ProfilerSession& session, const std::filesystem::path& output_path) {
  std::ofstream out(output_path);
  if (!out) {
    throw std::runtime_error("failed to open profiler CSV output path");
  }

  out << "path,call_count,inclusive_ms,exclusive_ms,bytes_moved\n";

  std::vector<std::pair<std::size_t, std::string>> stack;
  stack.emplace_back(session.rootNodeIndex(), "root");

  while (!stack.empty()) {
    const auto [node_index, path] = stack.back();
    stack.pop_back();

    const ProfileNode& node = session.nodes().at(node_index);
    out << path << ',' << node.call_count << ',' << std::fixed << std::setprecision(6) << node.inclusive_ms << ','
        << node.exclusive_ms << ',' << node.bytes_moved << '\n';

    for (auto it = node.children.rbegin(); it != node.children.rend(); ++it) {
      const ProfileNode& child = session.nodes().at(*it);
      stack.emplace_back(*it, path + "/" + child.name);
    }
  }
}

void writeOperationalReportJson(
    const ProfilerSession& session,
    const std::filesystem::path& output_path,
    std::string_view run_label,
    std::string_view provenance_config_hash_hex) {
  std::ofstream out(output_path);
  if (!out) {
    throw std::runtime_error("failed to open operational diagnostics JSON output path");
  }

  std::uint64_t warning_count = 0;
  std::uint64_t error_count = 0;
  std::uint64_t fatal_count = 0;
  for (const RuntimeEvent& event : session.events()) {
    if (event.severity == RuntimeEventSeverity::kWarning) {
      ++warning_count;
    } else if (event.severity == RuntimeEventSeverity::kError) {
      ++error_count;
    } else if (event.severity == RuntimeEventSeverity::kFatal) {
      ++fatal_count;
    }
  }

  out << "{\n";
  out << "  \"schema_version\": 1,\n";
  out << "  \"run_label\": \"" << escapeJson(std::string(run_label)) << "\",\n";
  out << "  \"provenance_config_hash_hex\": \"" << escapeJson(std::string(provenance_config_hash_hex)) << "\",\n";
  out << "  \"summary\": {\n";
  out << "    \"event_count\": " << session.events().size() << ",\n";
  out << "    \"warning_count\": " << warning_count << ",\n";
  out << "    \"error_count\": " << error_count << ",\n";
  out << "    \"fatal_count\": " << fatal_count << ",\n";
  out << "    \"status\": \"" << (fatal_count > 0 || error_count > 0 ? "error" : "ok") << "\"\n";
  out << "  },\n";
  if (const MemoryReport* report = session.memoryReport(); report != nullptr) {
    out << "  \"memory_report\": ";
    writeMemoryReportJson(out, *report, 2);
    out << ",\n";
  }
  out << "  \"events\": [";
  if (!session.events().empty()) {
    out << "\n";
    for (std::size_t i = 0; i < session.events().size(); ++i) {
      writeEventJson(out, session.events()[i], 4);
      if (i + 1 < session.events().size()) {
        out << ",";
      }
      out << "\n";
    }
    out << "  ";
  }
  out << "]\n";
  out << "}\n";
}

}  // namespace cosmosim::core
