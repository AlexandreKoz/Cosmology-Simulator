#include "cosmosim/core/profiling.hpp"

#include <algorithm>
#include <fstream>
#include <iomanip>
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

void AllocatorStats::recordAllocate(std::uint64_t bytes) {
  ++m_stats.alloc_calls;
  m_stats.bytes_allocated += bytes;
  m_stats.live_bytes += bytes;
  m_stats.peak_live_bytes = std::max(m_stats.peak_live_bytes, m_stats.live_bytes);
}

void AllocatorStats::recordFree(std::uint64_t bytes) {
  ++m_stats.free_calls;
  m_stats.bytes_freed += bytes;
  m_stats.live_bytes = (bytes > m_stats.live_bytes) ? 0 : (m_stats.live_bytes - bytes);
}

AllocatorStatsSnapshot AllocatorStats::snapshot() const noexcept { return m_stats; }

void AllocatorStats::reset() noexcept { m_stats = {}; }

void CounterRegistry::addCount(std::string_view key, std::uint64_t delta) {
  m_counts[std::string(key)] += delta;
}

void CounterRegistry::setCount(std::string_view key, std::uint64_t value) {
  m_counts[std::string(key)] = value;
}

std::uint64_t CounterRegistry::count(std::string_view key) const {
  const auto it = m_counts.find(std::string(key));
  return (it == m_counts.end()) ? 0 : it->second;
}

const std::unordered_map<std::string, std::uint64_t>& CounterRegistry::entries() const noexcept { return m_counts; }

void CounterRegistry::reset() { m_counts.clear(); }

ProfilerSession::ProfilerSession(bool enabled) : m_enabled(enabled) {
  m_nodes.push_back(ProfileNode{.name = "root", .call_count = 0});
}

void ProfilerSession::setEnabled(bool enabled) noexcept { m_enabled = enabled; }

bool ProfilerSession::enabled() const noexcept { return m_enabled; }

void ProfilerSession::beginScope(std::string_view phase_name) {
  if (!m_enabled) {
    return;
  }

  const std::size_t parent_index = m_scope_stack.empty() ? rootNodeIndex() : m_scope_stack.back().node_index;
  const std::size_t node_index = findOrCreateChild(parent_index, phase_name);

  ProfileNode& node = m_nodes[node_index];
  ++node.call_count;
  m_scope_stack.push_back(ActiveScope{
      .node_index = node_index,
      .start = Clock::now(),
      .child_elapsed_ms = 0.0,
  });
}

void ProfilerSession::endScope() {
  if (!m_enabled) {
    return;
  }
  if (m_scope_stack.empty()) {
    throw std::logic_error("ProfilerSession::endScope called with empty scope stack");
  }

  ActiveScope scope = m_scope_stack.back();
  m_scope_stack.pop_back();

  const auto end = Clock::now();
  const double elapsed_ms =
      std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(end - scope.start).count();

  ProfileNode& node = m_nodes[scope.node_index];
  node.inclusive_ms += elapsed_ms;
  const double exclusive_ms = std::max(0.0, elapsed_ms - scope.child_elapsed_ms);
  node.exclusive_ms += exclusive_ms;

  if (!m_scope_stack.empty()) {
    m_scope_stack.back().child_elapsed_ms += elapsed_ms;
  }
}

void ProfilerSession::addBytesMoved(std::uint64_t bytes) {
  if (!m_enabled) {
    return;
  }
  const std::size_t node_index = m_scope_stack.empty() ? rootNodeIndex() : m_scope_stack.back().node_index;
  m_nodes[node_index].bytes_moved += bytes;
}

CounterRegistry& ProfilerSession::counters() noexcept { return m_counters; }

const CounterRegistry& ProfilerSession::counters() const noexcept { return m_counters; }

AllocatorStats& ProfilerSession::allocatorStats() noexcept { return m_allocator_stats; }

const AllocatorStats& ProfilerSession::allocatorStats() const noexcept { return m_allocator_stats; }

const std::vector<ProfileNode>& ProfilerSession::nodes() const noexcept { return m_nodes; }

std::size_t ProfilerSession::rootNodeIndex() const noexcept { return 0; }

void ProfilerSession::reset() {
  m_nodes.clear();
  m_nodes.push_back(ProfileNode{.name = "root", .call_count = 0});
  m_scope_stack.clear();
  m_counters.reset();
  m_allocator_stats.reset();
}

std::size_t ProfilerSession::findOrCreateChild(std::size_t parent_index, std::string_view phase_name) {
  const std::string phase_key(phase_name);
  const auto it = m_nodes[parent_index].child_lookup.find(phase_key);
  if (it != m_nodes[parent_index].child_lookup.end()) {
    return it->second;
  }

  const std::size_t child_index = m_nodes.size();
  m_nodes.push_back(ProfileNode{.name = phase_key, .call_count = 0});
  m_nodes[parent_index].children.push_back(child_index);
  m_nodes[parent_index].child_lookup.emplace(phase_key, child_index);
  return child_index;
}

ScopedProfile::ScopedProfile(ProfilerSession* session, std::string_view phase_name) : m_session(session) {
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

  out << "  \"counters\": {";
  bool first = true;
  for (const auto& [key, value] : session.counters().entries()) {
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

}  // namespace cosmosim::core
