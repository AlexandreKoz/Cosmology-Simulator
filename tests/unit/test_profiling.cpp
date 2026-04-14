#include <cassert>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>

#include "cosmosim/core/profiling.hpp"

namespace {

void testNestedTimerAndBytes() {
  cosmosim::core::ProfilerSession session(true);

  {
    cosmosim::core::ScopedProfile step_scope(&session, "step");
    session.addBytesMoved(128);
    {
      cosmosim::core::ScopedProfile gravity_scope(&session, "gravity");
      session.addBytesMoved(64);
    }
  }

  const auto& nodes = session.nodes();
  assert(nodes.size() >= 3);

  const cosmosim::core::ProfileNode& step = nodes.at(1);
  const cosmosim::core::ProfileNode& gravity = nodes.at(2);

  assert(step.name == "step");
  assert(step.call_count == 1);
  assert(step.inclusive_ms >= step.exclusive_ms);
  assert(step.bytes_moved == 128);

  assert(gravity.name == "gravity");
  assert(gravity.call_count == 1);
  assert(gravity.inclusive_ms >= gravity.exclusive_ms);
  assert(gravity.bytes_moved == 64);
}

void testCounterAndAllocatorAggregation() {
  cosmosim::core::ProfilerSession session(true);

  session.counters().addCount("tree.node_visits", 5);
  session.counters().addCount("tree.node_visits", 7);
  session.counters().setCount("hydro.face_fluxes", 11);

  session.allocatorStats().recordAllocate(256);
  session.allocatorStats().recordAllocate(128);
  session.allocatorStats().recordFree(64);

  assert(session.counters().count("tree.node_visits") == 12);
  assert(session.counters().count("hydro.face_fluxes") == 11);

  const auto allocator = session.allocatorStats().snapshot();
  assert(allocator.alloc_calls == 2);
  assert(allocator.free_calls == 1);
  assert(allocator.bytes_allocated == 384);
  assert(allocator.bytes_freed == 64);
  assert(allocator.peak_live_bytes == 384);
  assert(allocator.live_bytes == 320);
}

void testJsonAndCsvReportWriters() {
  cosmosim::core::ProfilerSession session(true);
  {
    cosmosim::core::ScopedProfile root_scope(&session, "source_terms");
    session.counters().addCount("refinement.operations", 9);
    session.allocatorStats().recordAllocate(512);
    session.addBytesMoved(1024);
  }

  const auto temp_dir = std::filesystem::temp_directory_path();
  const auto json_path = temp_dir / "cosmosim_profile_unit.json";
  const auto csv_path = temp_dir / "cosmosim_profile_unit.csv";

  cosmosim::core::writeProfilerReportJson(session, json_path);
  cosmosim::core::writeProfilerReportCsv(session, csv_path);

  assert(std::filesystem::exists(json_path));
  assert(std::filesystem::exists(csv_path));

  const std::uintmax_t json_size = std::filesystem::file_size(json_path);
  const std::uintmax_t csv_size = std::filesystem::file_size(csv_path);
  assert(json_size > 0);
  assert(csv_size > 0);

  std::filesystem::remove(json_path);
  std::filesystem::remove(csv_path);
}

void testOperationalEventReportWriter() {
  cosmosim::core::ProfilerSession session(true);
  session.recordEvent(cosmosim::core::RuntimeEvent{
      .event_kind = "config.freeze",
      .severity = cosmosim::core::RuntimeEventSeverity::kInfo,
      .subsystem = "core.config",
      .step_index = 4,
      .simulation_time_code = 0.125,
      .scale_factor = 0.5,
      .message = "configuration frozen",
      .payload = {{"schema_version", "1"}},
  });
  session.recordEvent(cosmosim::core::RuntimeEvent{
      .event_kind = "restart.write.failure",
      .severity = cosmosim::core::RuntimeEventSeverity::kError,
      .subsystem = "io.restart",
      .step_index = 4,
      .simulation_time_code = 0.125,
      .scale_factor = 0.5,
      .message = "restart write failed",
      .payload = {{"error", "disk full"}},
  });

  const auto path = std::filesystem::temp_directory_path() / "cosmosim_operational_events_unit.json";
  cosmosim::core::writeOperationalReportJson(session, path, "unit_profiling", "cafef00d");
  assert(std::filesystem::exists(path));

  std::ifstream in(path);
  const std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
  assert(text.find("\"run_label\": \"unit_profiling\"") != std::string::npos);
  assert(text.find("\"provenance_config_hash_hex\": \"cafef00d\"") != std::string::npos);
  assert(text.find("\"event_kind\": \"restart.write.failure\"") != std::string::npos);
  assert(text.find("\"status\": \"error\"") != std::string::npos);

  std::filesystem::remove(path);
}

}  // namespace

int main() {
  testNestedTimerAndBytes();
  testCounterAndAllocatorAggregation();
  testJsonAndCsvReportWriters();
  testOperationalEventReportWriter();
  return 0;
}
