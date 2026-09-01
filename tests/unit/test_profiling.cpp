#include <cassert>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <limits>
#include <string>
#include <system_error>

#include "cosmosim/core/profiling.hpp"
#include "cosmosim/core/memory_governor.hpp"
#include "cosmosim/core/build_config.hpp"
#include "cosmosim/core/openmp_runtime.hpp"
#include "../support/test_temp_workspace.hpp"

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

  auto scratch = cosmosim::test_support::TestTempWorkspace::createUniqueDirectory("profiling_unit_reports");
  const auto json_path = scratch.path("profile.json");
  const auto csv_path = scratch.path("profile.csv");

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

  const auto path = cosmosim::test_support::TestTempWorkspace::uniqueProcessLocalPath("cosmosim_operational_events_unit.json");
  cosmosim::core::writeOperationalReportJson(session, path, "unit_profiling", "cafef00d");
  assert(std::filesystem::exists(path));

  std::ifstream in(path);
  const std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
  in.close();
  assert(text.find("\"run_label\": \"unit_profiling\"") != std::string::npos);
  assert(text.find("\"provenance_config_hash_hex\": \"cafef00d\"") != std::string::npos);
  assert(text.find("\"event_kind\": \"restart.write.failure\"") != std::string::npos);
  assert(text.find("\"status\": \"error\"") != std::string::npos);

  std::error_code cleanup_error;
  std::filesystem::remove(path, cleanup_error);
}


void testOperationalReportIncludesMemoryAccounting() {
  cosmosim::core::ProfilerSession session(true);
  cosmosim::core::MemoryReport report;
  report.totals.persistent_total_bytes = 1024;
  report.totals.transient_total_bytes = 256;
  report.notes.push_back("unknown external allocations are not fully tracked");
  cosmosim::core::MemoryGovernor governor(cosmosim::core::MemoryGovernorPolicy{
      .hard_limit_bytes = 4096U,
      .external_runtime_reserve_bytes = 128U,
  });
  governor.setBaselineOwnedBytes(1024U);
  auto reservation = governor.reserve(
      cosmosim::core::MemoryClass::kScratchArena, 256U, "unit.profiler");
  reservation.commit();
  cosmosim::core::attachMemoryGovernorSnapshot(report, governor);
  cosmosim::core::attachProcessMemoryObservation(
      report,
      cosmosim::core::ProcessMemoryObservation{
          .current_rss_bytes = 2048U,
          .peak_rss_bytes = 3072U,
          .pss_bytes = 1536U,
      });
  report.distributed_process_memory.rank_count = 2;
  report.distributed_process_memory.current_rss =
      cosmosim::core::DistributedMemoryMetricSummary{
          .valid = true,
          .local_bytes = 2048U,
          .global_sum_bytes = 4608U,
          .rank_max_bytes = 2560U,
          .rank_mean_bytes = 2304.0,
          .max_to_mean_imbalance_ratio = 2560.0 / 2304.0,
      };
  session.setMemoryReport(std::move(report));

  const auto path = cosmosim::test_support::TestTempWorkspace::uniqueProcessLocalPath("cosmosim_operational_events_memory_unit.json");
  cosmosim::core::writeOperationalReportJson(session, path, "unit_memory", "beadfeed");

  std::ifstream in(path);
  const std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
  in.close();
  assert(text.find("\"memory_report\"") != std::string::npos);
  assert(text.find("\"persistent_total_bytes\": 1024") != std::string::npos);
  assert(text.find("\"transient_total_bytes\": 256") != std::string::npos);
  assert(text.find("\"governor\"") != std::string::npos);
  assert(text.find("\"hard_limit_bytes\": 4096") != std::string::npos);
  assert(text.find("\"committed_bytes\": 256") != std::string::npos);
  assert(text.find("\"pressure\": \"green\"") != std::string::npos);
  assert(text.find("\"reservation_rejection_count\": 0") != std::string::npos);
  assert(text.find("\"process_memory\"") != std::string::npos);
  assert(text.find("\"known_accounted_bytes\": 1408") != std::string::npos);
  assert(text.find("\"observed_rss_bytes\": 2048") != std::string::npos);
  assert(text.find("\"observed_peak_rss_bytes\": 3072") != std::string::npos);
  assert(text.find("\"observed_pss_bytes\": 1536") != std::string::npos);
  assert(text.find("\"unexplained_resident_bytes\": 640") != std::string::npos);
  assert(text.find("\"distributed_process_memory\"") != std::string::npos);
  assert(text.find("\"rank_max_bytes\": 2560") != std::string::npos);
  std::error_code cleanup_error;
  std::filesystem::remove(path, cleanup_error);
}


void testOpenMpThreadLocalAggregation() {
#if COSMOSIM_HAVE_OPENMP
  cosmosim::core::configureOpenMpThreads(2);
  const int workers = cosmosim::core::configuredOpenMpThreadCount();
  assert(workers >= 1);
  cosmosim::core::ProfilerSession session(true);
#pragma omp parallel
  {
    cosmosim::core::ScopedProfile scope(&session, "omp_worker");
    for (int i = 0; i < 100; ++i) {
      session.counters().addCount("omp.iterations", 1);
    }
    session.addBytesMoved(64);
  }
  assert(session.counters().count("omp.iterations") == static_cast<std::uint64_t>(workers) * 100U);
  const auto& nodes = session.nodes();
  bool found = false;
  for (const auto& node : nodes) {
    if (node.name == "omp_worker") {
      assert(node.call_count == static_cast<std::uint64_t>(workers));
      assert(node.bytes_moved == static_cast<std::uint64_t>(workers) * 64U);
      found = true;
    }
  }
  assert(found);
  cosmosim::core::configureOpenMpThreads(1);
#endif
}

void testScopeMacroAndEscapedOperationalFields() {
  cosmosim::core::ProfilerSession session(true);
  {
    COSMOSIM_PROFILE_SCOPE(&session, "macro_outer");
    COSMOSIM_PROFILE_SCOPE(&session, "macro_inner");
  }
  session.recordEvent(cosmosim::core::RuntimeEvent{
      .event_kind = "escaped\tkind",
      .severity = cosmosim::core::RuntimeEventSeverity::kWarning,
      .subsystem = "core.profiler",
      .step_index = 7,
      .simulation_time_code = std::numeric_limits<double>::quiet_NaN(),
      .scale_factor = std::numeric_limits<double>::infinity(),
      .message = "line1\nline2\r\tcontrol\x01",
      .payload = {{"zeta", "last"}, {"alpha", "first"}},
  });
  const auto path = cosmosim::test_support::TestTempWorkspace::uniqueProcessLocalPath("cosmosim_profiler_escape_unit.json");
  cosmosim::core::writeOperationalReportJson(session, path, "escape_test", "hash");
  std::ifstream in(path);
  const std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
  assert(text.find("\\t") != std::string::npos);
  assert(text.find("\\u0001") != std::string::npos);
  assert(text.find("\"simulation_time_code\": null") != std::string::npos);
  assert(text.find("\"scale_factor\": null") != std::string::npos);
  assert(text.find("\"alpha\": \"first\"") < text.find("\"zeta\": \"last\""));
  std::error_code cleanup_error;
  std::filesystem::remove(path, cleanup_error);
}

}  // namespace

int main() {
  testNestedTimerAndBytes();
  testCounterAndAllocatorAggregation();
  testJsonAndCsvReportWriters();
  testOperationalEventReportWriter();
  testOperationalReportIncludesMemoryAccounting();
  testOpenMpThreadLocalAggregation();
  testScopeMacroAndEscapedOperationalFields();
  return 0;
}
