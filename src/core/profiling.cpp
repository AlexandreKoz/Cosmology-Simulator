#include "cosmosim/core/profiling.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <system_error>
#include <utility>

namespace cosmosim::core {
namespace {

[[nodiscard]] ProfileNode makeProfileNode(std::string name) {
  ProfileNode node;
  node.name = std::move(name);
  return node;
}

[[nodiscard]] std::string escapeJson(const std::string& input) {
  std::string out;
  out.reserve(input.size());
  constexpr char k_hex[] = "0123456789abcdef";
  for (const char raw_c : input) {
    const auto c = static_cast<unsigned char>(raw_c);
    switch (c) {
      case '"':
        out += "\\\"";
        break;
      case '\\':
        out += "\\\\";
        break;
      case '\b':
        out += "\\b";
        break;
      case '\f':
        out += "\\f";
        break;
      case '\n':
        out += "\\n";
        break;
      case '\r':
        out += "\\r";
        break;
      case '\t':
        out += "\\t";
        break;
      default:
        if (c < 0x20U) {
          out += "\\u00";
          out.push_back(k_hex[(c >> 4U) & 0x0FU]);
          out.push_back(k_hex[c & 0x0FU]);
        } else {
          out.push_back(static_cast<char>(c));
        }
        break;
    }
  }
  return out;
}

[[nodiscard]] std::string escapeCsvField(std::string_view input) {
  if (input.find_first_of(",\"\r\n") == std::string_view::npos) {
    return std::string(input);
  }
  std::string out;
  out.reserve(input.size() + 2U);
  out.push_back('"');
  for (const char c : input) {
    if (c == '"') {
      out += "\"\"";
    } else {
      out.push_back(c);
    }
  }
  out.push_back('"');
  return out;
}

template <typename Writer>
void writeAtomically(const std::filesystem::path& output_path, Writer&& writer) {
  std::filesystem::path part_path = output_path;
  part_path += ".part";
  {
    std::ofstream out(part_path, std::ios::out | std::ios::trunc);
    if (!out) {
      throw std::runtime_error("failed to open temporary profiler output path: " + part_path.string());
    }
    writer(out);
    out.flush();
    if (!out) {
      throw std::runtime_error("failed while writing profiler output: " + part_path.string());
    }
    out.close();
    if (!out) {
      throw std::runtime_error("failed while closing profiler output: " + part_path.string());
    }
  }

  std::error_code ec;
  std::filesystem::rename(part_path, output_path, ec);
  if (!ec) {
    return;
  }

  // Never delete a previously published report merely to make rename succeed.
  // POSIX rename replaces atomically; platforms that cannot atomically replace
  // an existing destination fail closed and preserve the old complete file.
  std::error_code remove_ec;
  std::filesystem::remove(part_path, remove_ec);
  throw std::runtime_error(
      "failed to atomically finalize profiler output '" + output_path.string() + "': " + ec.message());
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
  if (value.has_value() && std::isfinite(*value)) {
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
  std::vector<std::pair<std::string, std::string>> ordered_payload(event.payload.begin(), event.payload.end());
  std::sort(ordered_payload.begin(), ordered_payload.end(), [](const auto& lhs, const auto& rhs) {
    return lhs.first < rhs.first;
  });
  for (const auto& [key, value] : ordered_payload) {
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
  const auto write_optional_u64 = [&out](const std::optional<std::uint64_t>& value) {
    if (value.has_value()) {
      out << *value;
    } else {
      out << "null";
    }
  };
  const auto write_distributed_metric = [&out](
      const DistributedMemoryMetricSummary& metric) {
    out << "{\"valid\": " << (metric.valid ? "true" : "false")
        << ", \"local_bytes\": " << metric.local_bytes
        << ", \"global_sum_bytes\": " << metric.global_sum_bytes
        << ", \"rank_max_bytes\": " << metric.rank_max_bytes
        << ", \"rank_mean_bytes\": " << metric.rank_mean_bytes
        << ", \"max_to_mean_imbalance_ratio\": "
        << metric.max_to_mean_imbalance_ratio << "}";
  };
  out << indent << "{\n";
  out << field_indent << "\"persistent_total_bytes\": " << report.totals.persistent_total_bytes << ",\n";
  out << field_indent << "\"transient_total_bytes\": " << report.totals.transient_total_bytes << ",\n";
  out << field_indent << "\"unknown_external_allocations\": true,\n";
  out << field_indent << "\"distributed\": {\"valid\": "
      << (report.distributed.valid ? "true" : "false")
      << ", \"rank_count\": " << report.distributed.rank_count
      << ", \"local_owned_bytes\": " << report.distributed.local_owned_bytes
      << ", \"global_sum_owned_bytes\": " << report.distributed.global_sum_owned_bytes
      << ", \"rank_max_owned_bytes\": " << report.distributed.rank_max_owned_bytes
      << ", \"max_to_mean_imbalance_ratio\": "
      << report.distributed.max_to_mean_imbalance_ratio << "}";
  if (report.governor_snapshot.has_value()) {
    const MemoryGovernorSnapshot& governor = *report.governor_snapshot;
    out << ",\n";
    out << field_indent << "\"governor\": {"
        << "\"hard_limit_bytes\": " << governor.hard_limit_bytes
        << ", \"baseline_owned_bytes\": " << governor.baseline_owned_bytes
        << ", \"external_runtime_reserve_bytes\": " << governor.external_runtime_reserve_bytes
        << ", \"planned_overlap_reserve_bytes\": " << governor.planned_overlap_reserve_bytes
        << ", \"committed_bytes\": " << governor.committed_bytes
        << ", \"reserved_bytes\": " << governor.reserved_bytes
        << ", \"accounted_bytes\": " << governor.accounted_bytes
        << ", \"safety_adjusted_accounted_bytes\": "
        << governor.safety_adjusted_accounted_bytes
        << ", \"headroom_bytes\": " << governor.headroom_bytes
        << ", \"pressure\": \"" << memoryPressureLabel(governor.pressure) << "\""
        << ", \"peak_committed_bytes\": " << governor.peak_committed_bytes
        << ", \"peak_reserved_bytes\": " << governor.peak_reserved_bytes
        << ", \"peak_accounted_bytes\": " << governor.peak_accounted_bytes
        << ", \"reservation_rejection_count\": " << governor.rejection_count
        << "}";
  }
  if (report.process_memory_reconciliation.has_value()) {
    const ProcessMemoryReconciliation& process =
        *report.process_memory_reconciliation;
    out << ",\n" << field_indent << "\"process_memory\": {"
        << "\"current_declared_residency_bytes\": "
        << process.current_declared_residency_bytes
        << ", \"known_accounted_bytes\": " << process.known_accounted_bytes
        << ", \"observed_rss_bytes\": ";
    write_optional_u64(process.observed_rss_bytes);
    out << ", \"observed_peak_rss_bytes\": ";
    write_optional_u64(process.observed_peak_rss_bytes);
    out << ", \"observed_pss_bytes\": ";
    write_optional_u64(process.observed_pss_bytes);
    out << ", \"unexplained_resident_bytes\": ";
    write_optional_u64(process.unexplained_resident_bytes);
    out << ", \"observed_to_known_ratio\": ";
    if (process.observed_to_known_ratio.has_value()) {
      out << *process.observed_to_known_ratio;
    } else {
      out << "null";
    }
    out << "}";
  }
  const DistributedProcessMemorySummary& distributed_process =
      report.distributed_process_memory;
  if (distributed_process.known_accounted.valid ||
      distributed_process.current_rss.valid ||
      distributed_process.peak_rss.valid ||
      distributed_process.communication_high_water.valid) {
    out << ",\n" << field_indent << "\"distributed_process_memory\": {"
        << "\"rank_count\": " << distributed_process.rank_count
        << ", \"known_accounted\": ";
    write_distributed_metric(distributed_process.known_accounted);
    out << ", \"current_rss\": ";
    write_distributed_metric(distributed_process.current_rss);
    out << ", \"peak_rss\": ";
    write_distributed_metric(distributed_process.peak_rss);
    out << ", \"communication_high_water\": ";
    write_distributed_metric(distributed_process.communication_high_water);
    out << "}";
  }
  out << "\n";
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
  nodes.push_back(makeProfileNode("root"));
}

ProfilerSession::ProfilerSession(bool enabled)
    : m_session_id(g_next_profiler_session_id.fetch_add(1, std::memory_order_relaxed)),
      m_enabled(enabled) {
  m_merged_nodes.push_back(makeProfileNode("root"));
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
  shard.nodes.push_back(makeProfileNode(phase_key));
  shard.nodes[parent_index].children.push_back(child_index);
  shard.nodes[parent_index].child_lookup.emplace(phase_key, child_index);
  return child_index;
}

void ProfilerSession::rebuildMergedNodes() const {
  std::lock_guard lock(m_registration_mutex);
  m_merged_nodes.clear();
  m_merged_nodes.push_back(makeProfileNode("root"));

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
        m_merged_nodes.push_back(makeProfileNode(name));
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
  for (ProfileNode& node : m_merged_nodes) {
    std::sort(node.children.begin(), node.children.end(), [this](std::size_t lhs, std::size_t rhs) {
      return m_merged_nodes[lhs].name < m_merged_nodes[rhs].name;
    });
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
  for (const auto& shard : m_shards) {
    if (!shard->scope_stack.empty()) {
      throw std::logic_error("ProfilerSession::reset called while profiling scopes are active");
    }
  }
  for (auto& shard : m_shards) {
    shard->nodes.clear();
    shard->nodes.push_back(makeProfileNode("root"));
    shard->scope_stack.clear();
    shard->events.clear();
  }
  m_counters.reset();
  m_allocator_stats.reset();
  m_next_event_sequence.store(0, std::memory_order_relaxed);
  m_merged_nodes.clear();
  m_merged_nodes.push_back(makeProfileNode("root"));
  m_merged_events.clear();
  m_memory_report.reset();
}

ScopedProfile::ScopedProfile(ProfilerSession* session, std::string_view phase_name)
    : m_session(session), m_active(session != nullptr && session->enabled()) {
  if (m_active) {
    m_session->beginScope(phase_name);
  }
}

ScopedProfile::~ScopedProfile() noexcept {
  if (m_active) {
    try {
      m_session->endScope();
    } catch (...) {
      // Profiling is diagnostic infrastructure. Never turn exception unwinding
      // into std::terminate because a profiling scope was misused.
    }
  }
}

void writeProfilerReportJson(const ProfilerSession& session, const std::filesystem::path& output_path) {
  writeAtomically(output_path, [&](std::ostream& out) {
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

    const auto counter_snapshot = session.counters().snapshotEntries();
    std::vector<std::pair<std::string, std::uint64_t>> counters(
        counter_snapshot.begin(), counter_snapshot.end());
    std::sort(counters.begin(), counters.end(), [](const auto& lhs, const auto& rhs) {
      return lhs.first < rhs.first;
    });
    out << "  \"counters\": {";
    bool first = true;
    for (const auto& [key, value] : counters) {
      out << (first ? "\n" : ",\n");
      out << "    \"" << escapeJson(key) << "\": " << value;
      first = false;
    }
    if (!first) {
      out << "\n  ";
    }
    out << "},\n";

    const std::vector<ProfileNode>& nodes = session.nodes();
    out << "  \"phases\": ";
    writeNodeJson(out, nodes, session.rootNodeIndex(), 2);
    out << "\n}\n";
  });
}

void writeProfilerReportCsv(const ProfilerSession& session, const std::filesystem::path& output_path) {
  writeAtomically(output_path, [&](std::ostream& out) {
    out << "path,call_count,inclusive_ms,exclusive_ms,bytes_moved\n";

    const std::vector<ProfileNode>& nodes = session.nodes();
    std::vector<std::pair<std::size_t, std::string>> stack;
    stack.emplace_back(session.rootNodeIndex(), "root");

    while (!stack.empty()) {
      const auto [node_index, path] = stack.back();
      stack.pop_back();

      const ProfileNode& node = nodes.at(node_index);
      out << escapeCsvField(path) << ',' << node.call_count << ',' << std::fixed << std::setprecision(6)
          << node.inclusive_ms << ',' << node.exclusive_ms << ',' << node.bytes_moved << '\n';

      for (auto it = node.children.rbegin(); it != node.children.rend(); ++it) {
        const ProfileNode& child = nodes.at(*it);
        stack.emplace_back(*it, path + "/" + child.name);
      }
    }
  });
}

void writeOperationalReportJson(
    const ProfilerSession& session,
    const std::filesystem::path& output_path,
    std::string_view run_label,
    std::string_view provenance_config_hash_hex) {
  writeAtomically(output_path, [&](std::ostream& out) {
    const std::vector<RuntimeEvent>& events = session.events();
    std::uint64_t warning_count = 0;
    std::uint64_t error_count = 0;
    std::uint64_t fatal_count = 0;
    for (const RuntimeEvent& event : events) {
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
    out << "    \"event_count\": " << events.size() << ",\n";
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
    if (!events.empty()) {
      out << "\n";
      for (std::size_t i = 0; i < events.size(); ++i) {
        writeEventJson(out, events[i], 4);
        if (i + 1 < events.size()) {
          out << ',';
        }
        out << "\n";
      }
      out << "  ";
    }
    out << "]\n";
    out << "}\n";
  });
}


}  // namespace cosmosim::core
