#include <algorithm>
#include <chrono>
#include <cmath>
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
    const bool in_cluster = (i % 7U) < 5U;
    const double cluster_center = in_cluster ? 0.15 : 0.82;
    items[i].x_comov = std::fmod(cluster_center + 0.03 * std::sin(0.011 * static_cast<double>(i + 1U)), 1.0);
    if (items[i].x_comov < 0.0) {
      items[i].x_comov += 1.0;
    }
    items[i].y_comov = std::fmod(cluster_center + 0.04 * std::cos(0.017 * static_cast<double>(i + 3U)), 1.0);
    if (items[i].y_comov < 0.0) {
      items[i].y_comov += 1.0;
    }
    items[i].z_comov = std::fmod(cluster_center + 0.02 * std::sin(0.013 * static_cast<double>(i + 7U)), 1.0);
    if (items[i].z_comov < 0.0) {
      items[i].z_comov += 1.0;
    }
    items[i].active_target_count_recent = in_cluster ? 35U : 2U;
    items[i].remote_tree_interactions_recent = in_cluster ? 51U : 1U;
    items[i].work_units = 1.0 + static_cast<double>(i % 11);
    items[i].memory_bytes = 64U + static_cast<std::uint64_t>((i % 13) * 16U);
  }

  cosmosim::parallel::DecompositionConfig config;
  config.world_size = k_world_size;
  config.work_weight = 1.0;
  config.memory_weight = 1.0 / 1024.0;
  config.active_target_weight = 2.5;
  config.remote_tree_interaction_weight = 1.75;

  double warmup_checksum = 0.0;
  for (std::size_t i = 0; i < execution.warmup_iterations; ++i) {
    const auto plan = cosmosim::parallel::buildMortonSfcDecomposition(items, config);
    warmup_checksum += plan.metrics.max_weighted_load;
  }

  const auto begin = std::chrono::steady_clock::now();
  std::uint64_t total_send_bytes = 0;
  std::uint64_t total_recv_bytes = 0;
  double checksum = warmup_checksum;
  double weighted_imbalance_accum = 0.0;
  std::uint64_t remote_interactions_accum = 0;

  for (std::size_t iter = 0; iter < execution.measurement_iterations; ++iter) {
    const auto plan = cosmosim::parallel::buildMortonSfcDecomposition(items, config);
    checksum += plan.metrics.weighted_imbalance_ratio;
    weighted_imbalance_accum += plan.metrics.weighted_imbalance_ratio;
    for (const std::uint64_t per_rank : plan.metrics.remote_tree_interactions_by_rank) {
      remote_interactions_accum += per_rank;
    }

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
  reporter.addField(
      "mean_weighted_imbalance_ratio",
      weighted_imbalance_accum / static_cast<double>(execution.measurement_iterations));
  reporter.addField(
      "mean_remote_interactions_per_iter",
      static_cast<double>(remote_interactions_accum) / static_cast<double>(execution.measurement_iterations));
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
