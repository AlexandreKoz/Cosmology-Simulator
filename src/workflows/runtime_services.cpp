#include "cosmosim/workflows/runtime_services.hpp"

#include <algorithm>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <string>

#include "cosmosim/core/memory_accounting.hpp"
#include "cosmosim/parallel/distributed_memory.hpp"

namespace cosmosim::workflows {
namespace {

[[nodiscard]] core::DistributedMemoryMetricSummary reduceMemoryMetric(
    const parallel::MpiContext& mpi_context,
    const std::optional<std::uint64_t>& local_value) {
  const std::uint64_t available_rank_count = mpi_context.allreduceSumUint64(
      local_value.has_value() ? 1U : 0U);
  const int rank_count = std::max(mpi_context.worldSize(), 1);
  core::DistributedMemoryMetricSummary summary;
  if (available_rank_count != static_cast<std::uint64_t>(rank_count)) {
    return summary;
  }
  const std::uint64_t local = *local_value;
  const std::uint64_t global_sum = mpi_context.allreduceSumUint64(local);
  const std::uint64_t rank_max = mpi_context.allreduceMaxUint64(local);
  const double rank_mean = static_cast<double>(global_sum) /
      static_cast<double>(rank_count);
  summary.valid = true;
  summary.local_bytes = local;
  summary.global_sum_bytes = global_sum;
  summary.rank_max_bytes = rank_max;
  summary.rank_mean_bytes = rank_mean;
  summary.max_to_mean_imbalance_ratio = rank_mean > 0.0
      ? static_cast<double>(rank_max) / rank_mean
      : 0.0;
  return summary;
}

}  // namespace

void attachDistributedMemoryTelemetry(
    core::MemoryReport& report,
    const RuntimeServices& services) {
  const int rank_count = std::max(services.mpi_context.worldSize(), 1);
  const std::uint64_t local_owned = core::checkedMemoryBytesAdd(
      report.totals.persistent_total_bytes,
      report.totals.transient_total_bytes,
      "distributed memory owned capacity total");
  const std::uint64_t global_owned =
      services.mpi_context.allreduceSumUint64(local_owned);
  const std::uint64_t rank_max_owned =
      services.mpi_context.allreduceMaxUint64(local_owned);
  const double mean_owned = static_cast<double>(global_owned) /
      static_cast<double>(rank_count);
  report.distributed = core::DistributedMemorySummary{
      .valid = true,
      .rank_count = rank_count,
      .local_owned_bytes = local_owned,
      .global_sum_owned_bytes = global_owned,
      .rank_max_owned_bytes = rank_max_owned,
      .max_to_mean_imbalance_ratio = mean_owned > 0.0
          ? static_cast<double>(rank_max_owned) / mean_owned
          : 0.0,
  };

  report.distributed_process_memory.rank_count = rank_count;
  std::optional<std::uint64_t> known_accounted;
  std::optional<std::uint64_t> current_rss;
  std::optional<std::uint64_t> peak_rss;
  if (report.process_memory_reconciliation.has_value()) {
    known_accounted = report.process_memory_reconciliation->known_accounted_bytes;
    current_rss = report.process_memory_reconciliation->observed_rss_bytes;
    peak_rss = report.process_memory_reconciliation->observed_peak_rss_bytes;
  }
  report.distributed_process_memory.known_accounted =
      reduceMemoryMetric(services.mpi_context, known_accounted);
  report.distributed_process_memory.current_rss =
      reduceMemoryMetric(services.mpi_context, current_rss);
  report.distributed_process_memory.peak_rss =
      reduceMemoryMetric(services.mpi_context, peak_rss);
  report.distributed_process_memory.communication_high_water =
      reduceMemoryMetric(
          services.mpi_context,
          std::optional<std::uint64_t>(
              core::memoryReportCommunicationHighWaterBytes(report)));
}

FailureCoordinator::FailureCoordinator(const RuntimeServices& services) noexcept
    : m_services(services) {}

std::uint64_t FailureCoordinator::failedRankCount(bool local_failed) const {
  return m_services.mpi_context.allreduceSumUint64(local_failed ? 1ULL : 0ULL);
}

void FailureCoordinator::rethrowCollectiveFailure(
    const std::exception_ptr& local_failure,
    std::string_view phase_name) const {
  const std::uint64_t failed_rank_count =
      failedRankCount(local_failure != nullptr);
  if (failed_rank_count == 0U) {
    return;
  }
  if (local_failure != nullptr) {
    std::rethrow_exception(local_failure);
  }
  throw std::runtime_error(
      "collective phase '" + std::string(phase_name) + "' aborted because " +
      std::to_string(failed_rank_count) + " peer rank(s) reported failure");
}

}  // namespace cosmosim::workflows
