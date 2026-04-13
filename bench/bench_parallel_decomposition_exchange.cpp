#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <vector>

#include "bench/reporting/bench_report.hpp"
#include "cosmosim/parallel/distributed_memory.hpp"

int main() {
  constexpr std::size_t k_entity_count = 200000;
  constexpr int k_world_size = 4;

  const cosmosim::bench::BenchmarkExecutionConfig execution = cosmosim::bench::defaultExecutionConfig(3, 8);

  std::vector<cosmosim::parallel::DecompositionItem> items(k_entity_count);
  for (std::size_t i = 0; i < k_entity_count; ++i) {
    items[i].entity_id = static_cast<std::uint64_t>(i);
    items[i].kind = (i % 8 == 0) ? cosmosim::parallel::DecompositionEntityKind::kAmrPatch
                                 : cosmosim::parallel::DecompositionEntityKind::kParticle;
    items[i].x_comov = static_cast<double>(i % 512) / 512.0;
    items[i].y_comov = static_cast<double>((i / 64) % 512) / 512.0;
    items[i].z_comov = static_cast<double>((i / 4096) % 128) / 128.0;
    items[i].work_units = 1.0 + static_cast<double>(i % 11);
    items[i].memory_bytes = 64U + static_cast<std::uint64_t>((i % 13) * 16U);
  }

  cosmosim::parallel::DecompositionConfig config;
  config.world_size = k_world_size;
  config.work_weight = 1.0;
  config.memory_weight = 1.0 / 1024.0;

  double warmup_checksum = 0.0;
  for (std::size_t i = 0; i < execution.warmup_iterations; ++i) {
    const auto plan = cosmosim::parallel::buildMortonSfcDecomposition(items, config);
    warmup_checksum += plan.metrics.max_weighted_load;
  }

  const auto begin = std::chrono::steady_clock::now();
  std::uint64_t total_send_bytes = 0;
  std::uint64_t total_recv_bytes = 0;
  double checksum = warmup_checksum;

  for (std::size_t iter = 0; iter < execution.measurement_iterations; ++iter) {
    const auto plan = cosmosim::parallel::buildMortonSfcDecomposition(items, config);
    checksum += plan.metrics.weighted_imbalance_ratio;

    std::vector<int> ghost_owner_rank(items.size());
    for (std::size_t i = 0; i < items.size(); ++i) {
      ghost_owner_rank[i] = (plan.owning_rank_by_item[i] + 1) % k_world_size;
    }

    const auto ghost_plan = cosmosim::parallel::buildGhostExchangePlan(
        0,
        ghost_owner_rank,
        sizeof(std::uint64_t) + 3U * sizeof(double));
    total_send_bytes += ghost_plan.send_bytes;
    total_recv_bytes += ghost_plan.recv_bytes;
  }

  const auto end = std::chrono::steady_clock::now();
  const double elapsed_ms = std::chrono::duration<double, std::milli>(end - begin).count();

  cosmosim::bench::BenchmarkReporter reporter("bench_parallel_decomposition_exchange");
  cosmosim::bench::addExecutionFields(reporter, execution);
  reporter.addField("entity_count", k_entity_count);
  reporter.addField("world_size", k_world_size);
  reporter.addField("elapsed_ms", elapsed_ms);
  reporter.addField(
      "decomposition_per_iter_ms",
      elapsed_ms / static_cast<double>(execution.measurement_iterations));
  reporter.addField("ghost_send_bytes_total", total_send_bytes);
  reporter.addField("ghost_recv_bytes_total", total_recv_bytes);
  cosmosim::bench::addBandwidthFields(
      reporter,
      total_send_bytes + total_recv_bytes,
      elapsed_ms,
      "comm_bytes_total",
      "comm_effective_bandwidth_gb_s");
  reporter.addField("checksum", checksum);
  reporter.write(std::cout);

  return 0;
}
