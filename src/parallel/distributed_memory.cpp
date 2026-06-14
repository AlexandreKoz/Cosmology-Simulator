#include "cosmosim/parallel/distributed_memory.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <limits>
#include <numeric>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>

#include "cosmosim/core/build_config.hpp"

#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
#include <mpi.h>
#endif

namespace cosmosim::parallel {
namespace {

[[nodiscard]] double clampUnit(double value) {
  if (value <= 0.0) {
    return 0.0;
  }
  if (value >= 1.0) {
    return std::nextafter(1.0, 0.0);
  }
  return value;
}

[[nodiscard]] std::uint32_t quantize10bit(double coordinate, double min_coord, double max_coord) {
  const double extent = max_coord - min_coord;
  if (!(extent > 0.0)) {
    throw std::invalid_argument("decomposition domain extents must be positive");
  }
  const double normalized = clampUnit((coordinate - min_coord) / extent);
  constexpr double k_scale = 1024.0;
  const double scaled = std::floor(normalized * k_scale);
  const auto q = static_cast<std::uint32_t>(scaled);
  return std::min<std::uint32_t>(q, 1023U);
}

[[nodiscard]] std::uint64_t expandBits3d(std::uint32_t x) {
  std::uint64_t value = x & 0x3ffU;
  value = (value | (value << 16)) & 0x30000ffU;
  value = (value | (value << 8)) & 0x300f00fU;
  value = (value | (value << 4)) & 0x30c30c3U;
  value = (value | (value << 2)) & 0x9249249U;
  return value;
}

[[nodiscard]] std::uint64_t mortonKey3d(std::uint32_t x, std::uint32_t y, std::uint32_t z) {
  return (expandBits3d(z) << 2U) | (expandBits3d(y) << 1U) | expandBits3d(x);
}

[[nodiscard]] bool hasNonZeroComponentWeight(const DecompositionWeightCoefficients& weights) {
  return weights.particle_count != 0.0 || weights.gas_cell != 0.0 || weights.tree_interaction != 0.0 ||
      weights.pm_mesh != 0.0 || weights.amr_patch != 0.0 || weights.active_fraction != 0.0 ||
      weights.memory_pressure != 0.0 || weights.gpu_occupancy != 0.0 || weights.generic_work != 0.0;
}

[[nodiscard]] DecompositionWorkComponents effectiveWorkComponents(const DecompositionItem& item) {
  if (item.work_components.has_explicit_components) {
    DecompositionWorkComponents components = item.work_components;
    components.particle_count_cost = std::max(0.0, components.particle_count_cost);
    components.gas_cell_cost = std::max(0.0, components.gas_cell_cost);
    components.tree_interaction_cost = std::max(0.0, components.tree_interaction_cost);
    components.pm_mesh_cost = std::max(0.0, components.pm_mesh_cost);
    components.amr_patch_cost = std::max(0.0, components.amr_patch_cost);
    components.active_fraction_cost = std::max(0.0, components.active_fraction_cost);
    components.memory_pressure_cost = std::max(0.0, components.memory_pressure_cost);
    components.gpu_occupancy_cost = std::max(0.0, components.gpu_occupancy_cost);
    components.generic_work_cost = std::max(0.0, components.generic_work_cost);
    return components;
  }

  DecompositionWorkComponents components;
  components.particle_count_cost = (item.kind == DecompositionEntityKind::kParticle) ? 1.0 : 0.0;
  components.gas_cell_cost = (item.kind == DecompositionEntityKind::kHydroCell) ? 1.0 : 0.0;
  components.amr_patch_cost = (item.kind == DecompositionEntityKind::kAmrPatch) ? 1.0 : 0.0;
  components.pm_mesh_cost = (item.kind == DecompositionEntityKind::kPmMeshCell) ? 1.0 : 0.0;
  components.tree_interaction_cost = static_cast<double>(item.remote_tree_interactions_recent);
  components.active_fraction_cost = static_cast<double>(item.active_target_count_recent);
  components.memory_pressure_cost = static_cast<double>(item.memory_bytes);
  components.generic_work_cost = std::max(0.0, item.work_units);
  return components;
}

[[nodiscard]] double componentWeightedLoad(
    const DecompositionWorkComponents& components,
    const DecompositionWeightCoefficients& weights) {
  return weights.particle_count * components.particle_count_cost +
      weights.gas_cell * components.gas_cell_cost +
      weights.tree_interaction * components.tree_interaction_cost +
      weights.pm_mesh * components.pm_mesh_cost +
      weights.amr_patch * components.amr_patch_cost +
      weights.active_fraction * components.active_fraction_cost +
      weights.memory_pressure * components.memory_pressure_cost +
      weights.gpu_occupancy * components.gpu_occupancy_cost +
      weights.generic_work * components.generic_work_cost;
}

[[nodiscard]] double legacyWeightedLoad(const DecompositionItem& item, const DecompositionConfig& config) {
  const double owned_particle_term =
      (item.kind == DecompositionEntityKind::kParticle) ? 1.0 : 0.0;
  const double active_target_term = static_cast<double>(item.active_target_count_recent);
  const double remote_tree_term = static_cast<double>(item.remote_tree_interactions_recent);
  const double work_term = std::max(0.0, item.work_units);
  const double memory_term = static_cast<double>(item.memory_bytes);
  return config.owned_particle_weight * owned_particle_term +
      config.active_target_weight * active_target_term +
      config.remote_tree_interaction_weight * remote_tree_term +
      config.work_weight * work_term + config.memory_weight * memory_term;
}

[[nodiscard]] double weightedLoad(const DecompositionItem& item, const DecompositionConfig& config) {
  const DecompositionWorkComponents components = effectiveWorkComponents(item);
  if (config.prefer_component_work_model &&
      item.work_components.has_explicit_components && hasNonZeroComponentWeight(config.component_weights)) {
    return std::max(0.0, componentWeightedLoad(components, config.component_weights));
  }
  const double legacy = legacyWeightedLoad(item, config);
  if (legacy > 0.0) {
    return legacy;
  }
  const double component_fallback = componentWeightedLoad(components, config.component_weights);
  if (component_fallback > 0.0) {
    return component_fallback;
  }
  return (item.kind == DecompositionEntityKind::kParticle) ? 1.0 : std::max(1.0, components.rawTotal());
}

void validateComponentWeights(const DecompositionWeightCoefficients& weights) {
  if (weights.particle_count < 0.0 || weights.gas_cell < 0.0 || weights.tree_interaction < 0.0 ||
      weights.pm_mesh < 0.0 || weights.amr_patch < 0.0 || weights.active_fraction < 0.0 ||
      weights.memory_pressure < 0.0 || weights.gpu_occupancy < 0.0 || weights.generic_work < 0.0) {
    throw std::invalid_argument("decomposition component weights must be non-negative");
  }
}

void addWorkComponentsToMetrics(
    LoadBalanceMetrics& metrics,
    std::size_t rank,
    const DecompositionWorkComponents& components,
    double sign) {
  metrics.particle_count_cost_by_rank[rank] += sign * components.particle_count_cost;
  metrics.gas_cell_cost_by_rank[rank] += sign * components.gas_cell_cost;
  metrics.tree_interaction_cost_by_rank[rank] += sign * components.tree_interaction_cost;
  metrics.pm_mesh_cost_by_rank[rank] += sign * components.pm_mesh_cost;
  metrics.amr_patch_cost_by_rank[rank] += sign * components.amr_patch_cost;
  metrics.active_fraction_cost_by_rank[rank] += sign * components.active_fraction_cost;
  metrics.memory_pressure_cost_by_rank[rank] += sign * components.memory_pressure_cost;
  metrics.gpu_occupancy_cost_by_rank[rank] += sign * components.gpu_occupancy_cost;
  metrics.generic_work_cost_by_rank[rank] += sign * components.generic_work_cost;
}

template <typename T>
void appendPod(std::vector<std::uint8_t>& bytes, const T& value) {
  const std::size_t old_size = bytes.size();
  bytes.resize(old_size + sizeof(T));
  std::memcpy(bytes.data() + old_size, &value, sizeof(T));
}

template <typename T>
[[nodiscard]] T readPod(const std::vector<std::uint8_t>& bytes, std::size_t* offset) {
  if (*offset + sizeof(T) > bytes.size()) {
    throw std::runtime_error("ghost buffer decode overflow");
  }
  T value{};
  std::memcpy(&value, bytes.data() + *offset, sizeof(T));
  *offset += sizeof(T);
  return value;
}

[[nodiscard]] std::vector<std::string> splitLines(const std::string& text) {
  std::vector<std::string> lines;
  std::stringstream stream(text);
  std::string line;
  while (std::getline(stream, line)) {
    lines.push_back(line);
  }
  return lines;
}

[[nodiscard]] double absoluteValue(double value) {
  return (value < 0.0) ? -value : value;
}

[[nodiscard]] double stableRelativeError(double measured, double reference, double absolute_error) {
  const double denom = std::max(absoluteValue(reference), std::numeric_limits<double>::min());
  return absolute_error / denom;
}

[[nodiscard]] std::string rankConfigValueString(std::uint64_t value) {
  return std::to_string(value);
}

[[nodiscard]] std::string rankConfigValueString(int value) {
  return std::to_string(value);
}

[[nodiscard]] std::string rankConfigValueString(bool value) {
  return value ? "true" : "false";
}

void appendRankConfigMismatch(
    RankConfigConsensus* consensus,
    RankConfigMismatchProperty property,
    int baseline_rank,
    int rank,
    std::string baseline_value,
    std::string rank_value) {
  consensus->mismatches.push_back(RankConfigMismatch{
      .property = property,
      .baseline_rank = baseline_rank,
      .rank = rank,
      .baseline_value = std::move(baseline_value),
      .rank_value = std::move(rank_value),
  });
}

[[nodiscard]] constexpr std::size_t ghostExchangeRecordBytes() {
  // entity_id + position[3] + mass + density + velocity[3] + pressure + internal_energy.
  return sizeof(std::uint64_t) + sizeof(double) * 10U;
}

[[nodiscard]] bool laneIsPresentOrEmpty(std::size_t size, std::size_t expected) noexcept {
  return size == 0 || size == expected;
}

[[nodiscard]] double optionalLaneValue(const std::vector<double>& lane, std::size_t index) {
  return lane.empty() ? 0.0 : lane[index];
}

void resizeOptionalLaneForCommit(std::vector<double>* lane, std::size_t size) {
  if (lane->empty()) {
    lane->assign(size, 0.0);
  } else if (lane->size() != size) {
    throw std::invalid_argument("ghost optional lane size does not match storage size");
  }
}

void validateTransferDescriptor(
    const GhostTransferDescriptor& descriptor,
    GhostTransferRole expected_role,
    int expected_peer_rank,
    std::size_t expected_neighbor_slot,
    std::span<const std::uint32_t> expected_indices) {
  if (descriptor.role != expected_role) {
    throw std::invalid_argument("ghost transfer descriptor role does not match container");
  }
  if (descriptor.peer_rank != expected_peer_rank) {
    throw std::invalid_argument("ghost transfer descriptor peer_rank does not match neighbor slot");
  }
  if (descriptor.neighbor_slot != expected_neighbor_slot) {
    throw std::invalid_argument("ghost transfer descriptor neighbor_slot mismatch");
  }
  if (descriptor.local_indices.size() != expected_indices.size() ||
      !std::equal(descriptor.local_indices.begin(), descriptor.local_indices.end(), expected_indices.begin())) {
    throw std::invalid_argument("ghost transfer descriptor indices drift from canonical plan indices");
  }
  if (expected_role == GhostTransferRole::kOutboundSend &&
      descriptor.intent != GhostTransferIntent::kGhostRefreshRequest &&
      descriptor.intent != GhostTransferIntent::kOwnershipMigrationSend) {
    throw std::invalid_argument("outbound transfer intent must be ghost refresh request or migration send");
  }
  if (expected_role == GhostTransferRole::kInboundReceive &&
      descriptor.intent != GhostTransferIntent::kGhostRefreshReceiveStaging &&
      descriptor.intent != GhostTransferIntent::kOwnershipMigrationReceiveStaging) {
    throw std::invalid_argument("inbound transfer intent must be receive-staging intent");
  }
  if (descriptor.local_indices.empty() &&
      !(expected_role == GhostTransferRole::kOutboundSend &&
        descriptor.intent == GhostTransferIntent::kGhostRefreshRequest)) {
    throw std::invalid_argument("ghost transfer descriptor local_indices must be non-empty");
  }
  if (descriptor.intent == GhostTransferIntent::kGhostRefreshRequest ||
      descriptor.intent == GhostTransferIntent::kGhostRefreshReceiveStaging) {
    if (descriptor.expected_post_transfer_residency != LocalIndexResidency::kGhost) {
      throw std::invalid_argument("ghost refresh transfers must keep ghost post-transfer residency");
    }
  }
}

}  // namespace

bool GhostLayerEpoch::matches(const GhostLayerEpoch& expected) const noexcept {
  return decomposition_epoch == expected.decomposition_epoch && ghost_sync_epoch == expected.ghost_sync_epoch &&
      particle_index_generation == expected.particle_index_generation;
}

double DecompositionWorkComponents::rawTotal() const noexcept {
  return particle_count_cost + gas_cell_cost + tree_interaction_cost + pm_mesh_cost + amr_patch_cost +
      active_fraction_cost + memory_pressure_cost + gpu_occupancy_cost + generic_work_cost;
}

void validateOwnershipDescriptor(const OwnershipDescriptor& descriptor) {
  if (descriptor.owner_rank < 0 || descriptor.local_rank < 0) {
    throw std::invalid_argument("ownership descriptor ranks must be non-negative");
  }
  switch (descriptor.kind) {
    case ExchangeObjectKind::kLocalParticle:
      if (!descriptor.is_authoritative || !descriptor.is_mutable || descriptor.owner_rank != descriptor.local_rank) {
        throw std::invalid_argument("local particle descriptor must be authoritative and mutable only on owner rank");
      }
      break;
    case ExchangeObjectKind::kImportedGhostParticle:
      if (descriptor.is_authoritative || descriptor.is_mutable || descriptor.owner_rank == descriptor.local_rank) {
        throw std::invalid_argument("imported ghost particle descriptor must be non-authoritative read-only remote state");
      }
      break;
    case ExchangeObjectKind::kTreePseudoParticle:
      if (descriptor.is_authoritative || descriptor.is_mutable) {
        throw std::invalid_argument("tree pseudo-particle descriptor must be derived read-only exchange state");
      }
      break;
    case ExchangeObjectKind::kPmMeshCell:
      if (!descriptor.is_authoritative || descriptor.owner_rank != descriptor.local_rank) {
        throw std::invalid_argument("PM mesh cell descriptor must be authoritative on its owning mesh rank");
      }
      break;
    case ExchangeObjectKind::kHydroGhostCell:
      if (descriptor.is_authoritative || descriptor.is_mutable || descriptor.owner_rank == descriptor.local_rank) {
        throw std::invalid_argument("hydro ghost cell descriptor must be read-only boundary state on consumer rank");
      }
      break;
    case ExchangeObjectKind::kAmrPatchMetadata:
      if (descriptor.is_mutable && descriptor.owner_rank != descriptor.local_rank) {
        throw std::invalid_argument("remote AMR patch metadata cannot be mutable on non-owner rank");
      }
      break;
  }
}

DecompositionPlan buildMortonSfcDecomposition(std::span<const DecompositionItem> items, const DecompositionConfig& config) {
  if (config.world_size <= 0) {
    throw std::invalid_argument("world_size must be positive");
  }
  if (config.owned_particle_weight < 0.0 || config.active_target_weight < 0.0 ||
      config.remote_tree_interaction_weight < 0.0 || config.work_weight < 0.0 || config.memory_weight < 0.0) {
    throw std::invalid_argument("decomposition weights must be non-negative");
  }
  validateComponentWeights(config.component_weights);

  struct KeyedItem {
    std::size_t index = 0;
    std::uint64_t morton_key = 0;
    double weighted_load = 0.0;
    DecompositionWorkComponents components{};
  };

  std::vector<KeyedItem> keyed(items.size());
  for (std::size_t i = 0; i < items.size(); ++i) {
    const std::uint32_t qx = quantize10bit(items[i].x_comov, config.domain_x_min_comov, config.domain_x_max_comov);
    const std::uint32_t qy = quantize10bit(items[i].y_comov, config.domain_y_min_comov, config.domain_y_max_comov);
    const std::uint32_t qz = quantize10bit(items[i].z_comov, config.domain_z_min_comov, config.domain_z_max_comov);
    keyed[i] = KeyedItem{
        .index = i,
        .morton_key = mortonKey3d(qx, qy, qz),
        .weighted_load = weightedLoad(items[i], config),
        .components = effectiveWorkComponents(items[i]),
    };
  }

  std::stable_sort(keyed.begin(), keyed.end(), [&](const KeyedItem& a, const KeyedItem& b) {
    if (a.morton_key != b.morton_key) {
      return a.morton_key < b.morton_key;
    }
    return items[a.index].entity_id < items[b.index].entity_id;
  });

  DecompositionPlan plan;
  plan.owning_rank_by_item.assign(items.size(), 0);
  plan.sorted_indices.resize(items.size());
  plan.ranges_by_rank.assign(static_cast<std::size_t>(config.world_size), RankRange{});
  plan.metrics.weighted_load_by_rank.assign(static_cast<std::size_t>(config.world_size), 0.0);
  plan.metrics.memory_bytes_by_rank.assign(static_cast<std::size_t>(config.world_size), 0ULL);
  plan.metrics.owned_particles_by_rank.assign(static_cast<std::size_t>(config.world_size), 0ULL);
  plan.metrics.active_targets_by_rank.assign(static_cast<std::size_t>(config.world_size), 0ULL);
  plan.metrics.remote_tree_interactions_by_rank.assign(static_cast<std::size_t>(config.world_size), 0ULL);
  plan.metrics.particle_count_cost_by_rank.assign(static_cast<std::size_t>(config.world_size), 0.0);
  plan.metrics.gas_cell_cost_by_rank.assign(static_cast<std::size_t>(config.world_size), 0.0);
  plan.metrics.tree_interaction_cost_by_rank.assign(static_cast<std::size_t>(config.world_size), 0.0);
  plan.metrics.pm_mesh_cost_by_rank.assign(static_cast<std::size_t>(config.world_size), 0.0);
  plan.metrics.amr_patch_cost_by_rank.assign(static_cast<std::size_t>(config.world_size), 0.0);
  plan.metrics.active_fraction_cost_by_rank.assign(static_cast<std::size_t>(config.world_size), 0.0);
  plan.metrics.memory_pressure_cost_by_rank.assign(static_cast<std::size_t>(config.world_size), 0.0);
  plan.metrics.gpu_occupancy_cost_by_rank.assign(static_cast<std::size_t>(config.world_size), 0.0);
  plan.metrics.generic_work_cost_by_rank.assign(static_cast<std::size_t>(config.world_size), 0.0);

  for (std::size_t sorted_pos = 0; sorted_pos < keyed.size(); ++sorted_pos) {
    plan.sorted_indices[sorted_pos] = keyed[sorted_pos].index;
  }

  const double total_load = std::accumulate(
      keyed.begin(), keyed.end(), 0.0, [](double acc, const KeyedItem& entry) { return acc + entry.weighted_load; });

  if (!keyed.empty()) {
    const std::size_t active_rank_count =
        std::min<std::size_t>(static_cast<std::size_t>(config.world_size), keyed.size());
    const double target_per_active_rank =
        (active_rank_count > 0) ? total_load / static_cast<double>(active_rank_count) : 0.0;

    std::size_t current_rank = 0;
    std::size_t rank_begin = 0;
    double cumulative_load = 0.0;

    for (std::size_t sorted_pos = 0; sorted_pos < keyed.size(); ++sorted_pos) {
      cumulative_load += keyed[sorted_pos].weighted_load;
      if (current_rank + 1 >= active_rank_count) {
        continue;
      }

      const std::size_t items_remaining = keyed.size() - (sorted_pos + 1U);
      const std::size_t ranks_remaining = active_rank_count - (current_rank + 1U);
      const bool can_cut_after_current = (sorted_pos + 1U > rank_begin) && (items_remaining >= ranks_remaining);
      if (!can_cut_after_current) {
        continue;
      }

      const bool must_cut_to_keep_one_item_per_remaining_rank = items_remaining == ranks_remaining;
      const double next_target_prefix = target_per_active_rank * static_cast<double>(current_rank + 1U);
      const bool crossed_target = cumulative_load >= next_target_prefix;
      if (!must_cut_to_keep_one_item_per_remaining_rank && !crossed_target) {
        continue;
      }

      plan.ranges_by_rank[current_rank] = RankRange{
          .begin_sorted = rank_begin,
          .end_sorted = sorted_pos + 1U};
      ++current_rank;
      rank_begin = sorted_pos + 1U;
    }

    plan.ranges_by_rank[current_rank] = RankRange{.begin_sorted = rank_begin, .end_sorted = keyed.size()};
    for (std::size_t rank = current_rank + 1U; rank < static_cast<std::size_t>(config.world_size); ++rank) {
      plan.ranges_by_rank[rank] = RankRange{.begin_sorted = keyed.size(), .end_sorted = keyed.size()};
    }

    for (std::size_t rank = 0; rank < plan.ranges_by_rank.size(); ++rank) {
      const RankRange range = plan.ranges_by_rank[rank];
      for (std::size_t sorted_pos = range.begin_sorted; sorted_pos < range.end_sorted; ++sorted_pos) {
        const std::size_t original_index = keyed[sorted_pos].index;
        plan.owning_rank_by_item[original_index] = static_cast<int>(rank);
        plan.metrics.weighted_load_by_rank[rank] += keyed[sorted_pos].weighted_load;
        plan.metrics.memory_bytes_by_rank[rank] += items[original_index].memory_bytes;
        if (items[original_index].kind == DecompositionEntityKind::kParticle) {
          ++plan.metrics.owned_particles_by_rank[rank];
        }
        plan.metrics.active_targets_by_rank[rank] += items[original_index].active_target_count_recent;
        plan.metrics.remote_tree_interactions_by_rank[rank] += items[original_index].remote_tree_interactions_recent;
        addWorkComponentsToMetrics(plan.metrics, rank, keyed[sorted_pos].components, 1.0);
      }
    }
  }

  const auto max_load_it = std::max_element(plan.metrics.weighted_load_by_rank.begin(), plan.metrics.weighted_load_by_rank.end());
  plan.metrics.max_weighted_load =
      (max_load_it == plan.metrics.weighted_load_by_rank.end()) ? 0.0 : *max_load_it;
  plan.metrics.mean_weighted_load =
      plan.metrics.weighted_load_by_rank.empty()
          ? 0.0
          : (std::accumulate(plan.metrics.weighted_load_by_rank.begin(), plan.metrics.weighted_load_by_rank.end(), 0.0) /
             static_cast<double>(plan.metrics.weighted_load_by_rank.size()));
  plan.metrics.weighted_imbalance_ratio =
      (plan.metrics.mean_weighted_load > 0.0) ? (plan.metrics.max_weighted_load / plan.metrics.mean_weighted_load) : 0.0;

  plan.metrics.total_memory_bytes = std::accumulate(
      plan.metrics.memory_bytes_by_rank.begin(), plan.metrics.memory_bytes_by_rank.end(), 0ULL);
  const auto max_mem_it = std::max_element(plan.metrics.memory_bytes_by_rank.begin(), plan.metrics.memory_bytes_by_rank.end());
  plan.metrics.max_memory_bytes = (max_mem_it == plan.metrics.memory_bytes_by_rank.end()) ? 0ULL : *max_mem_it;
  const double mean_memory = plan.metrics.memory_bytes_by_rank.empty()
                                 ? 0.0
                                 : (static_cast<double>(plan.metrics.total_memory_bytes) /
                                    static_cast<double>(plan.metrics.memory_bytes_by_rank.size()));
  plan.metrics.memory_imbalance_ratio = (mean_memory > 0.0) ? (static_cast<double>(plan.metrics.max_memory_bytes) / mean_memory) : 0.0;

  return plan;
}

void applyRuntimeDecompositionFeedback(
    std::span<DecompositionItem> items,
    const DecompositionRuntimeMeasurements& measurements,
    const DecompositionFeedbackCoefficients& coefficients) {
  if (!measurements.has_measurements || items.empty()) {
    return;
  }
  if (coefficients.measured_tree_pair < 0.0 || coefficients.measured_pm_cell < 0.0 ||
      coefficients.measured_amr_cell < 0.0 || coefficients.measured_hydro_face < 0.0 ||
      coefficients.measured_wall_ms < 0.0) {
    throw std::invalid_argument("runtime decomposition feedback coefficients must be non-negative");
  }

  auto proxy_sum = [](std::span<DecompositionItem> entries, auto getter) {
    double sum = 0.0;
    for (const DecompositionItem& item : entries) {
      sum += std::max(0.0, getter(effectiveWorkComponents(item)));
    }
    return sum;
  };
  const double tree_proxy_sum = proxy_sum(items, [](const DecompositionWorkComponents& c) { return c.tree_interaction_cost; });
  const double pm_proxy_sum = proxy_sum(items, [](const DecompositionWorkComponents& c) { return c.pm_mesh_cost; });
  const double amr_proxy_sum = proxy_sum(items, [](const DecompositionWorkComponents& c) { return c.amr_patch_cost; });
  const double gas_proxy_sum = proxy_sum(items, [](const DecompositionWorkComponents& c) { return c.gas_cell_cost; });
  const double memory_proxy_sum = proxy_sum(items, [](const DecompositionWorkComponents& c) { return c.memory_pressure_cost; });
  const double generic_proxy_sum = proxy_sum(items, [](const DecompositionWorkComponents& c) { return std::max(1.0, c.generic_work_cost); });

  const double tree_total = coefficients.measured_tree_pair *
      static_cast<double>(measurements.tree_pair_evaluations_recent) +
      static_cast<double>(measurements.tree_remote_request_bytes_recent) / 1024.0;
  const double pm_total = coefficients.measured_pm_cell *
      static_cast<double>(measurements.pm_mesh_cells_touched_recent) +
      static_cast<double>(measurements.pm_fft_transpose_bytes_recent) / 1024.0;
  const double amr_total = coefficients.measured_amr_cell *
      static_cast<double>(measurements.amr_patch_cells_updated_recent);
  const double hydro_total = coefficients.measured_hydro_face *
      static_cast<double>(measurements.hydro_face_fluxes_recent);
  const double memory_total = static_cast<double>(measurements.ghost_exchange_bytes_recent);
  const double wall_total = coefficients.measured_wall_ms *
      (measurements.tree_wall_ms_recent + measurements.pm_wall_ms_recent + measurements.amr_wall_ms_recent +
       measurements.hydro_wall_ms_recent);

  auto distribute = [](double total, double proxy, double sum, std::size_t count) {
    if (!(total > 0.0)) {
      return 0.0;
    }
    if (sum > 0.0) {
      return total * std::max(0.0, proxy) / sum;
    }
    return total / static_cast<double>(std::max<std::size_t>(count, 1U));
  };

  for (DecompositionItem& item : items) {
    DecompositionWorkComponents components = effectiveWorkComponents(item);
    components.tree_interaction_cost += distribute(tree_total, components.tree_interaction_cost, tree_proxy_sum, items.size());
    components.pm_mesh_cost += distribute(pm_total, components.pm_mesh_cost, pm_proxy_sum, items.size());
    components.amr_patch_cost += distribute(amr_total, components.amr_patch_cost, amr_proxy_sum, items.size());
    components.gas_cell_cost += distribute(hydro_total, components.gas_cell_cost, gas_proxy_sum, items.size());
    components.memory_pressure_cost += distribute(memory_total, components.memory_pressure_cost, memory_proxy_sum, items.size());
    components.generic_work_cost += distribute(wall_total, std::max(1.0, components.generic_work_cost), generic_proxy_sum, items.size());
    components.gpu_occupancy_cost += std::max(0.0, measurements.gpu_kernel_ms_recent) *
        std::max(0.0, measurements.accelerator_occupancy_fraction_recent);
    components.has_explicit_components = true;
    item.work_components = components;
  }
}


LoadBalanceMetrics computeCurrentOwnershipLoadBalanceMetrics(
    std::span<const DecompositionItem> items,
    const DecompositionConfig& config) {
  if (config.world_size <= 0) {
    throw std::invalid_argument("current ownership metrics require positive world_size");
  }
  validateComponentWeights(config.component_weights);

  LoadBalanceMetrics metrics;
  const std::size_t rank_count = static_cast<std::size_t>(config.world_size);
  metrics.weighted_load_by_rank.assign(rank_count, 0.0);
  metrics.memory_bytes_by_rank.assign(rank_count, 0ULL);
  metrics.owned_particles_by_rank.assign(rank_count, 0ULL);
  metrics.active_targets_by_rank.assign(rank_count, 0ULL);
  metrics.remote_tree_interactions_by_rank.assign(rank_count, 0ULL);
  metrics.particle_count_cost_by_rank.assign(rank_count, 0.0);
  metrics.gas_cell_cost_by_rank.assign(rank_count, 0.0);
  metrics.tree_interaction_cost_by_rank.assign(rank_count, 0.0);
  metrics.pm_mesh_cost_by_rank.assign(rank_count, 0.0);
  metrics.amr_patch_cost_by_rank.assign(rank_count, 0.0);
  metrics.active_fraction_cost_by_rank.assign(rank_count, 0.0);
  metrics.memory_pressure_cost_by_rank.assign(rank_count, 0.0);
  metrics.gpu_occupancy_cost_by_rank.assign(rank_count, 0.0);
  metrics.generic_work_cost_by_rank.assign(rank_count, 0.0);

  for (const DecompositionItem& item : items) {
    if (item.current_owner_rank < 0 || item.current_owner_rank >= config.world_size) {
      throw std::invalid_argument("decomposition item current_owner_rank is outside runtime world size");
    }
    const std::size_t rank = static_cast<std::size_t>(item.current_owner_rank);
    metrics.weighted_load_by_rank[rank] += weightedLoad(item, config);
    metrics.memory_bytes_by_rank[rank] += item.memory_bytes;
    if (item.kind == DecompositionEntityKind::kParticle) {
      ++metrics.owned_particles_by_rank[rank];
    }
    metrics.active_targets_by_rank[rank] += item.active_target_count_recent;
    metrics.remote_tree_interactions_by_rank[rank] += item.remote_tree_interactions_recent;
    addWorkComponentsToMetrics(metrics, rank, effectiveWorkComponents(item), 1.0);
  }

  const auto max_load_it = std::max_element(metrics.weighted_load_by_rank.begin(), metrics.weighted_load_by_rank.end());
  metrics.max_weighted_load = (max_load_it == metrics.weighted_load_by_rank.end()) ? 0.0 : *max_load_it;
  metrics.mean_weighted_load = metrics.weighted_load_by_rank.empty()
      ? 0.0
      : (std::accumulate(metrics.weighted_load_by_rank.begin(), metrics.weighted_load_by_rank.end(), 0.0) /
         static_cast<double>(metrics.weighted_load_by_rank.size()));
  metrics.weighted_imbalance_ratio =
      (metrics.mean_weighted_load > 0.0) ? (metrics.max_weighted_load / metrics.mean_weighted_load) : 0.0;

  metrics.total_memory_bytes = std::accumulate(metrics.memory_bytes_by_rank.begin(), metrics.memory_bytes_by_rank.end(), 0ULL);
  const auto max_mem_it = std::max_element(metrics.memory_bytes_by_rank.begin(), metrics.memory_bytes_by_rank.end());
  metrics.max_memory_bytes = (max_mem_it == metrics.memory_bytes_by_rank.end()) ? 0ULL : *max_mem_it;
  const double mean_memory = metrics.memory_bytes_by_rank.empty()
      ? 0.0
      : (static_cast<double>(metrics.total_memory_bytes) / static_cast<double>(metrics.memory_bytes_by_rank.size()));
  metrics.memory_imbalance_ratio =
      (mean_memory > 0.0) ? (static_cast<double>(metrics.max_memory_bytes) / mean_memory) : 0.0;
  return metrics;
}

RuntimeRebalancePlan buildRuntimeRebalancePlan(
    std::span<const DecompositionItem> items,
    const DecompositionConfig& decomposition_config,
    const RuntimeRebalanceConfig& rebalance_config) {
  if (rebalance_config.world_size <= 0) {
    throw std::invalid_argument("runtime rebalance world_size must be positive");
  }
  if (rebalance_config.imbalance_trigger_ratio < 1.0 || rebalance_config.memory_trigger_ratio < 1.0 ||
      rebalance_config.max_migrated_load_fraction < 0.0 || rebalance_config.max_migrated_load_fraction > 1.0) {
    throw std::invalid_argument("runtime rebalance thresholds are invalid");
  }
  if (decomposition_config.world_size != rebalance_config.world_size) {
    throw std::invalid_argument("runtime rebalance config world_size must match decomposition world_size");
  }
  RuntimeRebalancePlan rebalance;
  rebalance.current_metrics = computeCurrentOwnershipLoadBalanceMetrics(items, decomposition_config);
  rebalance.target_decomposition = buildMortonSfcDecomposition(items, decomposition_config);
  if (items.empty()) {
    rebalance.reason = "empty_decomposition";
    return rebalance;
  }

  const bool load_imbalanced = rebalance.current_metrics.weighted_imbalance_ratio >=
      rebalance_config.imbalance_trigger_ratio;
  const bool memory_imbalanced = rebalance.current_metrics.memory_imbalance_ratio >=
      rebalance_config.memory_trigger_ratio;
  if (!load_imbalanced && !memory_imbalanced) {
    rebalance.reason = "below_rebalance_threshold";
    return rebalance;
  }

  const double total_load = std::accumulate(
      rebalance.current_metrics.weighted_load_by_rank.begin(),
      rebalance.current_metrics.weighted_load_by_rank.end(),
      0.0);
  const double max_migrated_load = rebalance_config.max_migrated_load_fraction * std::max(0.0, total_load);

  for (std::size_t item_index = 0; item_index < items.size(); ++item_index) {
    const int old_owner = items[item_index].current_owner_rank;
    const int new_owner = rebalance.target_decomposition.owning_rank_by_item[item_index];
    if (old_owner < 0 || old_owner == new_owner) {
      continue;
    }
    const double item_load = weightedLoad(items[item_index], decomposition_config);
    if (items[item_index].kind == DecompositionEntityKind::kParticle && rebalance_config.allow_particle_migration) {
      if (max_migrated_load > 0.0 && rebalance.migrated_load + item_load > max_migrated_load &&
          !rebalance.particle_migrations.empty()) {
        continue;
      }
      rebalance.particle_migrations.push_back(ParticleMigrationIntent{
          .particle_id = items[item_index].entity_id,
          .item_index = item_index,
          .old_owner_rank = old_owner,
          .new_owner_rank = new_owner,
          .work_units = item_load,
      });
      rebalance.migrated_load += item_load;
    } else if (items[item_index].kind == DecompositionEntityKind::kAmrPatch &&
               rebalance_config.allow_amr_patch_reassignment) {
      rebalance.amr_patch_ownership_updates.push_back(AmrPatchOwnershipUpdate{
          .patch_id = items[item_index].entity_id,
          .old_owner_rank = old_owner,
          .new_owner_rank = new_owner,
      });
      rebalance.migrated_load += item_load;
    }
  }

  rebalance.should_rebalance = !rebalance.particle_migrations.empty() || !rebalance.amr_patch_ownership_updates.empty();
  rebalance.migrated_load_fraction = (total_load > 0.0) ? (rebalance.migrated_load / total_load) : 0.0;
  rebalance.reason = load_imbalanced && memory_imbalanced ? "load_and_memory_imbalance" :
      (load_imbalanced ? "load_imbalance" : "memory_imbalance");
  return rebalance;
}

std::vector<DecompositionItem> gatherDecompositionItemsAcrossRanks(
    const MpiContext& mpi_context,
    std::span<const DecompositionItem> local_items) {
  static_assert(std::is_trivially_copyable_v<DecompositionItem>);
  if (mpi_context.worldSize() == 1) {
    return std::vector<DecompositionItem>(local_items.begin(), local_items.end());
  }
  if (!mpi_context.isEnabled()) {
    throw std::runtime_error("multi-rank decomposition item gather requires MPI to be enabled");
  }
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  const int local_bytes = static_cast<int>(local_items.size_bytes());
  std::vector<int> recv_counts(static_cast<std::size_t>(mpi_context.worldSize()), 0);
  MPI_Allgather(&local_bytes, 1, MPI_INT, recv_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);
  std::vector<int> displacements(recv_counts.size(), 0);
  int total_bytes = 0;
  for (std::size_t i = 0; i < recv_counts.size(); ++i) {
    if (recv_counts[i] % static_cast<int>(sizeof(DecompositionItem)) != 0) {
      throw std::runtime_error("decomposition item gather received a non-record-aligned byte count");
    }
    displacements[i] = total_bytes;
    total_bytes += recv_counts[i];
  }
  std::vector<std::uint8_t> recv_bytes(static_cast<std::size_t>(total_bytes));
  MPI_Allgatherv(
      local_items.data(),
      local_bytes,
      MPI_BYTE,
      recv_bytes.data(),
      recv_counts.data(),
      displacements.data(),
      MPI_BYTE,
      MPI_COMM_WORLD);
  std::vector<DecompositionItem> gathered(static_cast<std::size_t>(total_bytes) / sizeof(DecompositionItem));
  std::memcpy(gathered.data(), recv_bytes.data(), recv_bytes.size());
  return gathered;
#else
  throw std::runtime_error("multi-rank decomposition item gather requires an MPI build");
#endif
}

GhostExchangePlan buildGhostExchangePlan(
    int world_rank,
    std::span<const LocalGhostDescriptor> local_ghost_descriptors,
    std::size_t bytes_per_ghost) {
  if (bytes_per_ghost == 0) {
    throw std::invalid_argument("bytes_per_ghost must be positive");
  }
  GhostExchangePlan plan;
  std::vector<int> owners;
  owners.reserve(local_ghost_descriptors.size());

  bool has_epoch = false;
  GhostLayerEpoch common_epoch{};
  for (const LocalGhostDescriptor descriptor : local_ghost_descriptors) {
    if (!has_epoch) {
      common_epoch = descriptor.epoch;
      has_epoch = true;
    } else if (!descriptor.epoch.matches(common_epoch)) {
      throw std::invalid_argument("ghost descriptors in one exchange plan must share a common epoch");
    }
    if (descriptor.owning_rank < 0) {
      throw std::invalid_argument("ghost owner rank must be non-negative");
    }
    if (descriptor.residency == LocalIndexResidency::kOwned) {
      if (descriptor.owning_rank != world_rank) {
        throw std::invalid_argument("owned local index must have world_rank ownership");
      }
      continue;
    }
    if (descriptor.owning_rank == world_rank) {
      throw std::invalid_argument("ghost local index cannot be owned by world_rank");
    }
    owners.push_back(descriptor.owning_rank);
  }

  plan.epoch = common_epoch;
  plan.exchange_sequence = common_epoch.ghost_sync_epoch;
  std::sort(owners.begin(), owners.end());
  owners.erase(std::unique(owners.begin(), owners.end()), owners.end());

  plan.neighbor_ranks = owners;
  plan.send_local_indices_by_neighbor.assign(owners.size(), {});
  plan.recv_local_indices_by_neighbor.assign(owners.size(), {});
  plan.outbound_transfers.assign(owners.size(), {});
  plan.inbound_transfers.assign(owners.size(), {});

  for (std::uint32_t local_index = 0; local_index < local_ghost_descriptors.size(); ++local_index) {
    const LocalGhostDescriptor descriptor = local_ghost_descriptors[local_index];
    if (descriptor.residency == LocalIndexResidency::kOwned) {
      continue;
    }
    const auto it = std::lower_bound(owners.begin(), owners.end(), descriptor.owning_rank);
    if (it == owners.end() || *it != descriptor.owning_rank) {
      throw std::logic_error("owner rank map mismatch");
    }
    const std::size_t neighbor_slot = static_cast<std::size_t>(std::distance(owners.begin(), it));
    plan.recv_local_indices_by_neighbor[neighbor_slot].push_back(local_index);
  }

  for (std::size_t i = 0; i < owners.size(); ++i) {
    // Descriptor-only planning can identify local ghost import slots, but not the peer-owned
    // source rows that must be exported back. Keep outbound payload indices empty rather than
    // pretending the local ghost rows themselves are valid send sources.
    plan.send_local_indices_by_neighbor[i].clear();
    plan.outbound_transfers[i] = GhostTransferDescriptor{
        .role = GhostTransferRole::kOutboundSend,
        .intent = GhostTransferIntent::kGhostRefreshRequest,
        .peer_rank = owners[i],
        .neighbor_slot = i,
        .expected_post_transfer_residency = LocalIndexResidency::kGhost,
        .local_indices = plan.send_local_indices_by_neighbor[i],
    };
    plan.inbound_transfers[i] = GhostTransferDescriptor{
        .role = GhostTransferRole::kInboundReceive,
        .intent = GhostTransferIntent::kGhostRefreshReceiveStaging,
        .peer_rank = owners[i],
        .neighbor_slot = i,
        .expected_post_transfer_residency = LocalIndexResidency::kGhost,
        .local_indices = plan.recv_local_indices_by_neighbor[i],
    };
    plan.recv_bytes +=
        static_cast<std::uint64_t>(plan.recv_local_indices_by_neighbor[i].size()) * bytes_per_ghost;
    plan.send_bytes +=
        static_cast<std::uint64_t>(plan.send_local_indices_by_neighbor[i].size()) * bytes_per_ghost;
  }

  validateGhostExchangePlan(plan);
  return plan;
}

GhostExchangePlan buildGhostExchangePlan(
    int world_rank,
    std::span<const int> ghost_owner_rank_by_local_index,
    std::size_t bytes_per_ghost) {
  std::vector<LocalGhostDescriptor> descriptors;
  descriptors.reserve(ghost_owner_rank_by_local_index.size());
  for (const int owner_rank : ghost_owner_rank_by_local_index) {
    descriptors.push_back(LocalGhostDescriptor{
        .residency = (owner_rank == world_rank) ? LocalIndexResidency::kOwned : LocalIndexResidency::kGhost,
        .owning_rank = owner_rank,
    });
  }
  return buildGhostExchangePlan(world_rank, descriptors, bytes_per_ghost);
}

GhostExchangePlan buildExplicitGhostExchangePlan(
    int world_rank,
    std::span<const int> neighbor_ranks,
    std::span<const std::vector<std::uint32_t>> send_local_indices_by_neighbor,
    std::span<const std::vector<std::uint32_t>> recv_local_indices_by_neighbor,
    std::size_t bytes_per_ghost,
    const GhostLayerEpoch& epoch,
    bool enable_nonblocking_overlap) {
  if (world_rank < 0) {
    throw std::invalid_argument("world_rank must be non-negative");
  }
  if (bytes_per_ghost == 0) {
    throw std::invalid_argument("bytes_per_ghost must be positive");
  }
  if (neighbor_ranks.size() != send_local_indices_by_neighbor.size() ||
      neighbor_ranks.size() != recv_local_indices_by_neighbor.size()) {
    throw std::invalid_argument("explicit ghost exchange plan container sizes must match");
  }

  GhostExchangePlan plan;
  plan.neighbor_ranks.assign(neighbor_ranks.begin(), neighbor_ranks.end());
  plan.send_local_indices_by_neighbor.assign(
      send_local_indices_by_neighbor.begin(), send_local_indices_by_neighbor.end());
  plan.recv_local_indices_by_neighbor.assign(
      recv_local_indices_by_neighbor.begin(), recv_local_indices_by_neighbor.end());
  plan.outbound_transfers.resize(neighbor_ranks.size());
  plan.inbound_transfers.resize(neighbor_ranks.size());
  plan.epoch = epoch;
  plan.exchange_sequence = epoch.ghost_sync_epoch;
  plan.uses_blocking_exchange = true;
  plan.nonblocking_overlap_enabled = enable_nonblocking_overlap;

  for (std::size_t i = 0; i < neighbor_ranks.size(); ++i) {
    if (neighbor_ranks[i] < 0 || neighbor_ranks[i] == world_rank) {
      throw std::invalid_argument("explicit ghost exchange neighbor rank must be a remote non-negative rank");
    }
    plan.outbound_transfers[i] = GhostTransferDescriptor{
        .role = GhostTransferRole::kOutboundSend,
        .intent = GhostTransferIntent::kGhostRefreshRequest,
        .peer_rank = neighbor_ranks[i],
        .neighbor_slot = i,
        .expected_post_transfer_residency = LocalIndexResidency::kGhost,
        .local_indices = plan.send_local_indices_by_neighbor[i],
    };
    plan.inbound_transfers[i] = GhostTransferDescriptor{
        .role = GhostTransferRole::kInboundReceive,
        .intent = GhostTransferIntent::kGhostRefreshReceiveStaging,
        .peer_rank = neighbor_ranks[i],
        .neighbor_slot = i,
        .expected_post_transfer_residency = LocalIndexResidency::kGhost,
        .local_indices = plan.recv_local_indices_by_neighbor[i],
    };
    plan.send_bytes += static_cast<std::uint64_t>(plan.send_local_indices_by_neighbor[i].size()) * bytes_per_ghost;
    plan.recv_bytes += static_cast<std::uint64_t>(plan.recv_local_indices_by_neighbor[i].size()) * bytes_per_ghost;
  }

  validateGhostExchangePlan(plan);
  return plan;
}

void validateGhostExchangePlan(const GhostExchangePlan& plan) {
  if (!plan.uses_blocking_exchange && !plan.nonblocking_overlap_enabled) {
    throw std::invalid_argument("ghost exchange plan must expose either the default blocking path or an explicit overlap path");
  }
  if (plan.nonblocking_overlap_enabled && !plan.uses_blocking_exchange) {
    throw std::invalid_argument("nonblocking ghost exchange overlap must share the blocking ownership contract");
  }
  const std::size_t neighbor_count = plan.neighbor_ranks.size();
  if (plan.send_local_indices_by_neighbor.size() != neighbor_count ||
      plan.recv_local_indices_by_neighbor.size() != neighbor_count ||
      plan.outbound_transfers.size() != neighbor_count ||
      plan.inbound_transfers.size() != neighbor_count) {
    throw std::invalid_argument("ghost exchange plan containers must have matching neighbor counts");
  }
  for (std::size_t i = 0; i < neighbor_count; ++i) {
    if (i > 0 && plan.neighbor_ranks[i - 1] >= plan.neighbor_ranks[i]) {
      throw std::invalid_argument("ghost exchange plan neighbor_ranks must be strictly increasing");
    }
    validateTransferDescriptor(
        plan.outbound_transfers[i],
        GhostTransferRole::kOutboundSend,
        plan.neighbor_ranks[i],
        i,
        plan.send_local_indices_by_neighbor[i]);
    validateTransferDescriptor(
        plan.inbound_transfers[i],
        GhostTransferRole::kInboundReceive,
        plan.neighbor_ranks[i],
        i,
        plan.recv_local_indices_by_neighbor[i]);
  }
}

void validateGhostTransferAgainstResidency(
    const GhostTransferDescriptor& descriptor,
    std::span<const LocalGhostDescriptor> local_ghost_descriptors,
    int world_rank) {
  if (world_rank < 0) {
    throw std::invalid_argument("world_rank must be non-negative");
  }
  for (const std::uint32_t local_index : descriptor.local_indices) {
    if (local_index >= local_ghost_descriptors.size()) {
      throw std::out_of_range("ghost transfer descriptor local index out of residency table range");
    }
    const LocalGhostDescriptor local = local_ghost_descriptors[local_index];
    if (descriptor.role == GhostTransferRole::kOutboundSend) {
      if (local.residency != LocalIndexResidency::kOwned || local.owning_rank != world_rank) {
        throw std::invalid_argument("outbound ghost or migration payload must be packed from authoritative local state");
      }
    } else {
      if (descriptor.intent == GhostTransferIntent::kGhostRefreshReceiveStaging) {
        if (local.residency != LocalIndexResidency::kGhost || local.owning_rank == world_rank) {
          throw std::invalid_argument("ghost refresh receive staging must unpack into remote-owned ghost slots");
        }
      } else if (descriptor.intent == GhostTransferIntent::kOwnershipMigrationReceiveStaging) {
        if (descriptor.expected_post_transfer_residency != LocalIndexResidency::kOwned) {
          throw std::invalid_argument("ownership migration receive staging must produce owned local state");
        }
      }
    }
  }
}

void validateBlockingGhostExchangeContracts(
    const GhostExchangePlan& plan,
    std::span<const LocalGhostDescriptor> local_ghost_descriptors,
    int world_rank,
    const GhostLayerEpoch& expected_epoch) {
  validateGhostExchangePlan(plan);
  if (!plan.uses_blocking_exchange) {
    throw std::invalid_argument("default ghost exchange path must be blocking and correctness-first");
  }
  if (!plan.epoch.matches(expected_epoch)) {
    throw std::invalid_argument("ghost exchange plan epoch is stale for the current decomposition/sync generation");
  }
  for (const GhostTransferDescriptor& descriptor : plan.outbound_transfers) {
    validateGhostTransferAgainstResidency(descriptor, local_ghost_descriptors, world_rank);
  }
  for (const GhostTransferDescriptor& descriptor : plan.inbound_transfers) {
    validateGhostTransferAgainstResidency(descriptor, local_ghost_descriptors, world_rank);
  }
}


LocalOwnershipIdentitySummary summarizeLocalOwnedParticleIds(std::span<const std::uint64_t> local_particle_ids) {
  LocalOwnershipIdentitySummary summary;
  summary.local_owned_count = static_cast<std::uint64_t>(local_particle_ids.size());
  std::unordered_set<std::uint64_t> seen;
  seen.reserve(local_particle_ids.size());
  for (const std::uint64_t particle_id : local_particle_ids) {
    summary.local_particle_id_sum += particle_id;
    summary.local_particle_id_square_sum += particle_id * particle_id;
    summary.local_particle_id_xor ^= particle_id;
    if (!seen.insert(particle_id).second) {
      summary.local_particle_ids_unique = false;
    }
  }
  return summary;
}

ExactOwnershipPartitionReport validateExactGlobalOwnershipPartition(
    const MpiContext& mpi_context,
    std::span<const std::uint64_t> local_owned_particle_ids,
    std::span<const std::uint64_t> expected_global_particle_ids) {
  ExactOwnershipPartitionReport report;
  std::vector<std::uint64_t> local_sorted(local_owned_particle_ids.begin(), local_owned_particle_ids.end());
  std::sort(local_sorted.begin(), local_sorted.end());
  report.local_particle_ids_unique = std::adjacent_find(local_sorted.begin(), local_sorted.end()) == local_sorted.end();

  std::vector<std::uint64_t> global_ids;
  if (!mpi_context.isEnabled()) {
    global_ids = std::move(local_sorted);
  } else {
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
    const std::uint64_t local_count = static_cast<std::uint64_t>(local_owned_particle_ids.size());
    std::vector<std::uint64_t> counts64(static_cast<std::size_t>(mpi_context.worldSize()), 0U);
    MPI_Allgather(
        const_cast<std::uint64_t*>(&local_count),
        1,
        MPI_UINT64_T,
        counts64.data(),
        1,
        MPI_UINT64_T,
        MPI_COMM_WORLD);
    std::vector<int> recv_counts(static_cast<std::size_t>(mpi_context.worldSize()), 0);
    std::vector<int> recv_displs(static_cast<std::size_t>(mpi_context.worldSize()), 0);
    std::uint64_t total_count = 0;
    for (int rank = 0; rank < mpi_context.worldSize(); ++rank) {
      if (counts64[static_cast<std::size_t>(rank)] > static_cast<std::uint64_t>(std::numeric_limits<int>::max())) {
        throw std::overflow_error("exact ownership partition gather exceeds MPI int count limit");
      }
      recv_counts[static_cast<std::size_t>(rank)] = static_cast<int>(counts64[static_cast<std::size_t>(rank)]);
      if (rank > 0) {
        recv_displs[static_cast<std::size_t>(rank)] =
            recv_displs[static_cast<std::size_t>(rank - 1)] + recv_counts[static_cast<std::size_t>(rank - 1)];
      }
      total_count += counts64[static_cast<std::size_t>(rank)];
    }
    if (total_count > static_cast<std::uint64_t>(std::numeric_limits<int>::max())) {
      throw std::overflow_error("exact ownership partition global count exceeds MPI int count limit");
    }
    global_ids.resize(static_cast<std::size_t>(total_count));
    MPI_Allgatherv(
        const_cast<std::uint64_t*>(local_owned_particle_ids.data()),
        static_cast<int>(local_owned_particle_ids.size()),
        MPI_UINT64_T,
        global_ids.data(),
        recv_counts.data(),
        recv_displs.data(),
        MPI_UINT64_T,
        MPI_COMM_WORLD);
#else
    throw std::runtime_error("exact global ownership validation requires MPI support when MPI context is enabled");
#endif
  }

  std::sort(global_ids.begin(), global_ids.end());
  report.global_owned_count = static_cast<std::uint64_t>(global_ids.size());
  for (auto it = global_ids.begin(); it != global_ids.end();) {
    const auto range = std::equal_range(it, global_ids.end(), *it);
    if (std::distance(range.first, range.second) > 1) {
      report.duplicate_particle_ids.push_back(*it);
    }
    it = range.second;
  }
  report.globally_unique = report.duplicate_particle_ids.empty();

  std::vector<std::uint64_t> expected(expected_global_particle_ids.begin(), expected_global_particle_ids.end());
  std::sort(expected.begin(), expected.end());
  expected.erase(std::unique(expected.begin(), expected.end()), expected.end());
  std::vector<std::uint64_t> global_unique = global_ids;
  global_unique.erase(std::unique(global_unique.begin(), global_unique.end()), global_unique.end());
  std::set_difference(
      expected.begin(),
      expected.end(),
      global_unique.begin(),
      global_unique.end(),
      std::back_inserter(report.missing_expected_particle_ids));
  std::set_difference(
      global_unique.begin(),
      global_unique.end(),
      expected.begin(),
      expected.end(),
      std::back_inserter(report.extra_particle_ids));
  report.matches_expected_ids =
      report.missing_expected_particle_ids.empty() && report.extra_particle_ids.empty() &&
      global_unique.size() == expected.size();
  return report;
}

bool partitionIdentityMatchesGeneratedSet(
    const LocalOwnershipIdentitySummary& reduced_global_summary,
    std::uint64_t expected_global_count,
    std::uint64_t expected_particle_id_sum,
    std::uint64_t expected_particle_id_square_sum,
    std::uint64_t expected_particle_id_xor) {
  return reduced_global_summary.local_particle_ids_unique &&
      reduced_global_summary.local_owned_count == expected_global_count &&
      reduced_global_summary.local_particle_id_sum == expected_particle_id_sum &&
      reduced_global_summary.local_particle_id_square_sum == expected_particle_id_square_sum &&
      reduced_global_summary.local_particle_id_xor == expected_particle_id_xor;
}

bool partitionIdentityMatchesGeneratedSet(
    const LocalOwnershipIdentitySummary& reduced_global_summary,
    std::uint64_t expected_global_count,
    std::uint64_t expected_particle_id_sum,
    std::uint64_t expected_particle_id_xor) {
  return partitionIdentityMatchesGeneratedSet(
      reduced_global_summary,
      expected_global_count,
      expected_particle_id_sum,
      reduced_global_summary.local_particle_id_square_sum,
      expected_particle_id_xor);
}

double deterministicRankOrderedSum(std::span<const double> per_rank_values) {
  long double sum = 0.0;
  for (const double value : per_rank_values) {
    sum += static_cast<long double>(value);
  }
  return static_cast<double>(sum);
}

ReductionAgreement compareReductionAgreement(
    std::span<const double> per_rank_values,
    double measured_sum) {
  const double deterministic_baseline_sum = deterministicRankOrderedSum(per_rank_values);
  const double absolute_error = absoluteValue(measured_sum - deterministic_baseline_sum);
  return ReductionAgreement{
      .deterministic_baseline_sum = deterministic_baseline_sum,
      .measured_sum = measured_sum,
      .absolute_error = absolute_error,
      .relative_error = stableRelativeError(measured_sum, deterministic_baseline_sum, absolute_error),
  };
}

bool satisfiesReductionAgreement(
    const ReductionAgreement& agreement,
    const ReductionAgreementPolicy& policy) {
  if (policy.absolute_tolerance < 0.0 || policy.relative_tolerance < 0.0) {
    throw std::invalid_argument("reduction agreement tolerances must be non-negative");
  }
  const bool absolute_ok = agreement.absolute_error <= policy.absolute_tolerance;
  const bool relative_ok = agreement.relative_error <= policy.relative_tolerance;
  switch (policy.mode) {
    case ReductionAgreementMode::kAbsoluteOnly:
      return absolute_ok;
    case ReductionAgreementMode::kRelativeOnly:
      return relative_ok;
    case ReductionAgreementMode::kAbsoluteAndRelative:
      return absolute_ok && relative_ok;
    case ReductionAgreementMode::kAbsoluteOrRelative:
      return absolute_ok || relative_ok;
  }
  throw std::invalid_argument("unknown reduction agreement mode");
}

bool RankConfigConsensus::allConsistent() const noexcept {
  return normalized_config_hash_match && mpi_ranks_expected_match && deterministic_reduction_match;
}

RankConfigConsensus evaluateRankConfigConsensus(std::span<const RankConfigDigest> digests) {
  RankConfigConsensus consensus;
  if (digests.empty()) {
    return consensus;
  }

  const RankConfigDigest baseline = digests.front();
  for (const RankConfigDigest& digest : digests) {
    bool rank_matches = true;
    if (digest.normalized_config_hash != baseline.normalized_config_hash) {
      consensus.normalized_config_hash_match = false;
      rank_matches = false;
      appendRankConfigMismatch(
          &consensus,
          RankConfigMismatchProperty::kNormalizedConfigHash,
          baseline.world_rank,
          digest.world_rank,
          rankConfigValueString(baseline.normalized_config_hash),
          rankConfigValueString(digest.normalized_config_hash));
    }
    if (digest.mpi_ranks_expected != baseline.mpi_ranks_expected) {
      consensus.mpi_ranks_expected_match = false;
      rank_matches = false;
      appendRankConfigMismatch(
          &consensus,
          RankConfigMismatchProperty::kMpiRanksExpected,
          baseline.world_rank,
          digest.world_rank,
          rankConfigValueString(baseline.mpi_ranks_expected),
          rankConfigValueString(digest.mpi_ranks_expected));
    }
    if (digest.deterministic_reduction != baseline.deterministic_reduction) {
      consensus.deterministic_reduction_match = false;
      rank_matches = false;
      appendRankConfigMismatch(
          &consensus,
          RankConfigMismatchProperty::kDeterministicReduction,
          baseline.world_rank,
          digest.world_rank,
          rankConfigValueString(baseline.deterministic_reduction),
          rankConfigValueString(digest.deterministic_reduction));
    }
    if (!rank_matches) {
      consensus.mismatched_ranks.push_back(digest.world_rank);
    }
  }
  return consensus;
}

std::size_t ghostRefreshPayloadRecordBytes() noexcept {
  return ghostExchangeRecordBytes();
}

void validateGhostRefreshPayloadDescriptor(const GhostTransferDescriptor& descriptor) {
  if (descriptor.intent != GhostTransferIntent::kGhostRefreshRequest &&
      descriptor.intent != GhostTransferIntent::kGhostRefreshReceiveStaging) {
    throw std::invalid_argument(
        "GhostExchangeBuffer carries ghost-refresh payloads only; ownership migration must use ParticleMigrationRecord");
  }
  if (descriptor.expected_post_transfer_residency != LocalIndexResidency::kGhost) {
    throw std::invalid_argument("ghost-refresh payload descriptor must produce ghost residency");
  }
}

int ghostExchangePairStableTag(int tag_base, int local_rank, int peer_rank) {
  if (tag_base < 0 || local_rank < 0 || peer_rank < 0 || local_rank == peer_rank) {
    throw std::invalid_argument("ghostExchangePairStableTag: invalid rank pair or tag base");
  }
  // MPI receive matching constrains source and destination ranks. Keep this tag
  // independent of local neighbor-slot order; sequence separation is provided by
  // ghostExchangeSequencedTag for overlapping or repeated phases.
  return tag_base;
}

int ghostExchangeSequencedTag(
    int tag_base,
    int local_rank,
    int peer_rank,
    std::uint64_t exchange_sequence) {
  constexpr int k_sequence_stride = 16;
  constexpr int k_sequence_window = 64;
  if (tag_base < 0) {
    throw std::invalid_argument("ghostExchangeSequencedTag: tag_base must be non-negative");
  }
  const int phased_base = tag_base +
      static_cast<int>(exchange_sequence % static_cast<std::uint64_t>(k_sequence_window)) * k_sequence_stride;
  return ghostExchangePairStableTag(phased_base, local_rank, peer_rank);
}

bool GhostExchangeBufferSoA::isConsistent() const noexcept {
  const std::size_t n = entity_id.size();
  return laneIsPresentOrEmpty(position_x_comoving.size(), n) &&
      laneIsPresentOrEmpty(position_y_comoving.size(), n) &&
      laneIsPresentOrEmpty(position_z_comoving.size(), n) &&
      laneIsPresentOrEmpty(mass_code.size(), n) &&
      laneIsPresentOrEmpty(density_code.size(), n) &&
      laneIsPresentOrEmpty(velocity_x_code.size(), n) &&
      laneIsPresentOrEmpty(velocity_y_code.size(), n) &&
      laneIsPresentOrEmpty(velocity_z_code.size(), n) &&
      laneIsPresentOrEmpty(pressure_code.size(), n) &&
      laneIsPresentOrEmpty(internal_energy_code.size(), n);
}

std::size_t GhostExchangeBufferSoA::size() const noexcept { return entity_id.size(); }

bool GhostExchangeBufferSoA::hasGravityPayload() const noexcept {
  const std::size_t n = entity_id.size();
  return position_x_comoving.size() == n && position_y_comoving.size() == n &&
      position_z_comoving.size() == n && mass_code.size() == n;
}

bool GhostExchangeBufferSoA::hasHydroPayload() const noexcept {
  const std::size_t n = entity_id.size();
  return density_code.size() == n && velocity_x_code.size() == n && velocity_y_code.size() == n &&
      velocity_z_code.size() == n && pressure_code.size() == n && internal_energy_code.size() == n;
}

std::size_t ReadOnlyGhostExchangeView::size() const noexcept { return entity_id.size(); }

bool ReadOnlyGhostExchangeView::isConsistent() const noexcept {
  const std::size_t n = entity_id.size();
  return laneIsPresentOrEmpty(position_x_comoving.size(), n) &&
      laneIsPresentOrEmpty(position_y_comoving.size(), n) &&
      laneIsPresentOrEmpty(position_z_comoving.size(), n) &&
      laneIsPresentOrEmpty(mass_code.size(), n) &&
      laneIsPresentOrEmpty(density_code.size(), n) &&
      laneIsPresentOrEmpty(velocity_x_code.size(), n) &&
      laneIsPresentOrEmpty(velocity_y_code.size(), n) &&
      laneIsPresentOrEmpty(velocity_z_code.size(), n) &&
      laneIsPresentOrEmpty(pressure_code.size(), n) &&
      laneIsPresentOrEmpty(internal_energy_code.size(), n);
}

bool ReadOnlyGhostExchangeView::isFresh(const GhostLayerEpoch& expected_epoch) const noexcept {
  return epoch.matches(expected_epoch);
}

ReadOnlyGhostExchangeView makeReadOnlyGhostExchangeView(const GhostExchangeBufferSoA& storage) {
  if (!storage.isConsistent()) {
    throw std::invalid_argument("ghost storage must be component-consistent before building a read-only view");
  }
  return ReadOnlyGhostExchangeView{
      .epoch = storage.epoch,
      .entity_id = storage.entity_id,
      .position_x_comoving = storage.position_x_comoving,
      .position_y_comoving = storage.position_y_comoving,
      .position_z_comoving = storage.position_z_comoving,
      .mass_code = storage.mass_code,
      .density_code = storage.density_code,
      .velocity_x_code = storage.velocity_x_code,
      .velocity_y_code = storage.velocity_y_code,
      .velocity_z_code = storage.velocity_z_code,
      .pressure_code = storage.pressure_code,
      .internal_energy_code = storage.internal_energy_code,
  };
}

void requireFreshGhostExchangeView(
    const ReadOnlyGhostExchangeView& view,
    const GhostLayerEpoch& expected_epoch) {
  if (!view.isConsistent()) {
    throw std::invalid_argument("read-only ghost view component sizes are inconsistent");
  }
  if (!view.isFresh(expected_epoch)) {
    throw std::invalid_argument("read-only ghost view is stale for the current exchange epoch");
  }
}

void GhostExchangeBuffer::clear() { m_bytes.clear(); }

std::size_t GhostExchangeBuffer::byteSize() const noexcept { return m_bytes.size(); }

std::span<const std::uint8_t> GhostExchangeBuffer::encodedBytes() const noexcept { return m_bytes; }

void GhostExchangeBuffer::replaceEncodedBytes(std::vector<std::uint8_t> bytes) {
  m_bytes = std::move(bytes);
}

void GhostExchangeBuffer::packFrom(const GhostExchangeBufferSoA& source, std::span<const std::uint32_t> local_indices) {
  if (!source.isConsistent()) {
    throw std::invalid_argument("ghost source SoA fields must have matching sizes");
  }

  m_bytes.clear();
  appendPod<std::uint64_t>(m_bytes, static_cast<std::uint64_t>(local_indices.size()));

  for (const std::uint32_t index : local_indices) {
    if (index >= source.size()) {
      throw std::out_of_range("ghost pack local index out of range");
    }

    appendPod<std::uint64_t>(m_bytes, source.entity_id[index]);
    appendPod<double>(m_bytes, optionalLaneValue(source.position_x_comoving, index));
    appendPod<double>(m_bytes, optionalLaneValue(source.position_y_comoving, index));
    appendPod<double>(m_bytes, optionalLaneValue(source.position_z_comoving, index));
    appendPod<double>(m_bytes, optionalLaneValue(source.mass_code, index));
    appendPod<double>(m_bytes, optionalLaneValue(source.density_code, index));
    appendPod<double>(m_bytes, optionalLaneValue(source.velocity_x_code, index));
    appendPod<double>(m_bytes, optionalLaneValue(source.velocity_y_code, index));
    appendPod<double>(m_bytes, optionalLaneValue(source.velocity_z_code, index));
    appendPod<double>(m_bytes, optionalLaneValue(source.pressure_code, index));
    appendPod<double>(m_bytes, optionalLaneValue(source.internal_energy_code, index));
  }
}

void GhostExchangeBuffer::packFrom(
    const GhostTransferDescriptor& descriptor,
    const GhostExchangeBufferSoA& source,
    std::span<const std::uint32_t> local_indices) {
  validateGhostRefreshPayloadDescriptor(descriptor);
  if (descriptor.local_indices.size() != local_indices.size() ||
      !std::equal(descriptor.local_indices.begin(), descriptor.local_indices.end(), local_indices.begin())) {
    throw std::invalid_argument("ghost descriptor indices must match packed local indices");
  }
  packFrom(source, local_indices);
}

void GhostExchangeBuffer::unpackAppendTo(GhostExchangeBufferSoA& destination) const {
  if (!destination.isConsistent()) {
    throw std::invalid_argument("ghost destination SoA fields must have matching sizes");
  }
  if (m_bytes.size() < sizeof(std::uint64_t)) {
    throw std::runtime_error("ghost buffer is too small");
  }

  std::size_t offset = 0;
  const std::uint64_t count = readPod<std::uint64_t>(m_bytes, &offset);
  const std::uint64_t expected_payload_bytes =
      static_cast<std::uint64_t>(sizeof(std::uint64_t)) + count * static_cast<std::uint64_t>(ghostExchangeRecordBytes());
  if (expected_payload_bytes != static_cast<std::uint64_t>(m_bytes.size())) {
    throw std::runtime_error("ghost buffer payload shape does not match encoded count");
  }

  const std::size_t append_count = static_cast<std::size_t>(count);
  destination.entity_id.reserve(destination.entity_id.size() + append_count);
  destination.position_x_comoving.reserve(destination.position_x_comoving.size() + append_count);
  destination.position_y_comoving.reserve(destination.position_y_comoving.size() + append_count);
  destination.position_z_comoving.reserve(destination.position_z_comoving.size() + append_count);
  destination.mass_code.reserve(destination.mass_code.size() + append_count);
  destination.density_code.reserve(destination.density_code.size() + append_count);
  destination.velocity_x_code.reserve(destination.velocity_x_code.size() + append_count);
  destination.velocity_y_code.reserve(destination.velocity_y_code.size() + append_count);
  destination.velocity_z_code.reserve(destination.velocity_z_code.size() + append_count);
  destination.pressure_code.reserve(destination.pressure_code.size() + append_count);
  destination.internal_energy_code.reserve(destination.internal_energy_code.size() + append_count);

  for (std::uint64_t i = 0; i < count; ++i) {
    destination.entity_id.push_back(readPod<std::uint64_t>(m_bytes, &offset));
    destination.position_x_comoving.push_back(readPod<double>(m_bytes, &offset));
    destination.position_y_comoving.push_back(readPod<double>(m_bytes, &offset));
    destination.position_z_comoving.push_back(readPod<double>(m_bytes, &offset));
    destination.mass_code.push_back(readPod<double>(m_bytes, &offset));
    destination.density_code.push_back(readPod<double>(m_bytes, &offset));
    destination.velocity_x_code.push_back(readPod<double>(m_bytes, &offset));
    destination.velocity_y_code.push_back(readPod<double>(m_bytes, &offset));
    destination.velocity_z_code.push_back(readPod<double>(m_bytes, &offset));
    destination.pressure_code.push_back(readPod<double>(m_bytes, &offset));
    destination.internal_energy_code.push_back(readPod<double>(m_bytes, &offset));
  }

  if (offset != m_bytes.size()) {
    throw std::runtime_error("ghost buffer decode found trailing bytes");
  }
}

void GhostExchangeBuffer::unpackAppendTo(
    const GhostTransferDescriptor& descriptor,
    GhostExchangeBufferSoA& destination) const {
  validateGhostRefreshPayloadDescriptor(descriptor);
  if (m_bytes.size() < sizeof(std::uint64_t)) {
    throw std::runtime_error("ghost buffer is too small");
  }
  std::size_t offset = 0;
  const std::uint64_t encoded_count = readPod<std::uint64_t>(m_bytes, &offset);
  if (encoded_count != descriptor.local_indices.size()) {
    throw std::runtime_error("ghost buffer encoded count does not match receive descriptor slots");
  }
  unpackAppendTo(destination);
}

std::string DistributedRestartState::serialize() const {
  std::ostringstream stream;
  stream << std::setprecision(std::numeric_limits<double>::max_digits10);
  stream << "schema_version=" << schema_version << '\n';
  stream << "decomposition_epoch=" << decomposition_epoch << '\n';
  stream << "world_size=" << world_size << '\n';
  stream << "pm_grid_nx=" << pm_grid_nx << '\n';
  stream << "pm_grid_ny=" << pm_grid_ny << '\n';
  stream << "pm_grid_nz=" << pm_grid_nz << '\n';
  stream << "pm_decomposition_mode=" << pm_decomposition_mode << '\n';
  stream << "gravity_kick_opportunity=" << gravity_kick_opportunity << '\n';
  stream << "pm_update_cadence_steps=" << pm_update_cadence_steps << '\n';
  stream << "long_range_field_version=" << long_range_field_version << '\n';
  stream << "last_long_range_refresh_opportunity=" << last_long_range_refresh_opportunity << '\n';
  stream << "long_range_field_built_step_index=" << long_range_field_built_step_index << '\n';
  stream << "long_range_field_built_scale_factor=" << long_range_field_built_scale_factor << '\n';
  stream << "long_range_restart_policy=" << long_range_restart_policy << '\n';
  stream << "item_count=" << owning_rank_by_item.size() << '\n';
  for (std::size_t i = 0; i < owning_rank_by_item.size(); ++i) {
    stream << "rank[" << i << "]=" << owning_rank_by_item[i] << '\n';
  }
  stream << "pm_slab_rank_count=" << pm_slab_begin_x_by_rank.size() << '\n';
  for (std::size_t rank = 0; rank < pm_slab_begin_x_by_rank.size(); ++rank) {
    stream << "pm_slab_begin_x[" << rank << "]=" << pm_slab_begin_x_by_rank[rank] << '\n';
    stream << "pm_slab_end_x[" << rank << "]=" << pm_slab_end_x_by_rank[rank] << '\n';
  }
  return stream.str();
}

DistributedRestartState DistributedRestartState::deserialize(const std::string& encoded) {
  DistributedRestartState state;
  const std::vector<std::string> lines = splitLines(encoded);
  std::size_t expected_item_count = 0;
  std::vector<bool> seen_rank_entry;
  std::size_t expected_slab_rank_count = 0;
  std::vector<bool> seen_slab_begin;
  std::vector<bool> seen_slab_end;

  for (const std::string& line : lines) {
    if (line.empty()) {
      continue;
    }

    const std::size_t eq = line.find('=');
    if (eq == std::string::npos) {
      throw std::invalid_argument("invalid restart encoding line");
    }
    const std::string key = line.substr(0, eq);
    const std::string value = line.substr(eq + 1);

    if (key == "schema_version") {
      state.schema_version = static_cast<std::uint32_t>(std::stoul(value));
    } else if (key == "decomposition_epoch") {
      state.decomposition_epoch = std::stoull(value);
    } else if (key == "world_size") {
      state.world_size = std::stoi(value);
    } else if (key == "pm_grid_nx") {
      state.pm_grid_nx = static_cast<std::size_t>(std::stoull(value));
    } else if (key == "pm_grid_ny") {
      state.pm_grid_ny = static_cast<std::size_t>(std::stoull(value));
    } else if (key == "pm_grid_nz") {
      state.pm_grid_nz = static_cast<std::size_t>(std::stoull(value));
    } else if (key == "pm_decomposition_mode") {
      state.pm_decomposition_mode = value;
    } else if (key == "gravity_kick_opportunity") {
      state.gravity_kick_opportunity = std::stoull(value);
    } else if (key == "pm_update_cadence_steps") {
      state.pm_update_cadence_steps = std::stoull(value);
    } else if (key == "long_range_field_version") {
      state.long_range_field_version = std::stoull(value);
    } else if (key == "last_long_range_refresh_opportunity") {
      state.last_long_range_refresh_opportunity = std::stoull(value);
    } else if (key == "long_range_field_built_step_index") {
      state.long_range_field_built_step_index = std::stoull(value);
    } else if (key == "long_range_field_built_scale_factor") {
      state.long_range_field_built_scale_factor = std::stod(value);
    } else if (key == "long_range_restart_policy") {
      state.long_range_restart_policy = value;
    } else if (key == "item_count") {
      expected_item_count = static_cast<std::size_t>(std::stoull(value));
      state.owning_rank_by_item.assign(expected_item_count, 0);
      seen_rank_entry.assign(expected_item_count, false);
    } else if (key == "pm_slab_rank_count") {
      expected_slab_rank_count = static_cast<std::size_t>(std::stoull(value));
      state.pm_slab_begin_x_by_rank.assign(expected_slab_rank_count, 0);
      state.pm_slab_end_x_by_rank.assign(expected_slab_rank_count, 0);
      seen_slab_begin.assign(expected_slab_rank_count, false);
      seen_slab_end.assign(expected_slab_rank_count, false);
    } else if (key.rfind("rank[", 0) == 0) {
      const std::size_t open = key.find('[');
      const std::size_t close = key.find(']');
      if (open == std::string::npos || close == std::string::npos || close <= open + 1) {
        throw std::invalid_argument("invalid rank entry in restart encoding");
      }
      const std::size_t index = static_cast<std::size_t>(std::stoull(key.substr(open + 1, close - open - 1)));
      if (index >= state.owning_rank_by_item.size()) {
        throw std::out_of_range("restart rank index out of bounds");
      }
      if (seen_rank_entry[index]) {
        throw std::invalid_argument("duplicate restart rank entry");
      }
      state.owning_rank_by_item[index] = std::stoi(value);
      seen_rank_entry[index] = true;
    } else if (key.rfind("pm_slab_begin_x[", 0) == 0 || key.rfind("pm_slab_end_x[", 0) == 0) {
      const bool is_begin = key.rfind("pm_slab_begin_x[", 0) == 0;
      const std::size_t open = key.find('[');
      const std::size_t close = key.find(']');
      if (open == std::string::npos || close == std::string::npos || close <= open + 1) {
        throw std::invalid_argument("invalid PM slab entry in restart encoding");
      }
      const std::size_t rank_index = static_cast<std::size_t>(std::stoull(key.substr(open + 1, close - open - 1)));
      if (rank_index >= expected_slab_rank_count) {
        throw std::out_of_range("restart PM slab rank index out of bounds");
      }
      if (is_begin) {
        if (seen_slab_begin[rank_index]) {
          throw std::invalid_argument("duplicate PM slab begin entry");
        }
        state.pm_slab_begin_x_by_rank[rank_index] = static_cast<std::size_t>(std::stoull(value));
        seen_slab_begin[rank_index] = true;
      } else {
        if (seen_slab_end[rank_index]) {
          throw std::invalid_argument("duplicate PM slab end entry");
        }
        state.pm_slab_end_x_by_rank[rank_index] = static_cast<std::size_t>(std::stoull(value));
        seen_slab_end[rank_index] = true;
      }
    }
  }

  if (state.owning_rank_by_item.size() != expected_item_count) {
    throw std::runtime_error("restart decode item count mismatch");
  }
  if (state.world_size <= 0) {
    throw std::invalid_argument("restart world_size must be positive");
  }
  if (state.pm_update_cadence_steps == 0) {
    throw std::invalid_argument("restart PM cadence must be >= 1");
  }
  if (state.last_long_range_refresh_opportunity > state.gravity_kick_opportunity) {
    throw std::invalid_argument("restart cadence state is inconsistent: last refresh opportunity exceeds kick opportunity");
  }
  if (state.long_range_field_version == 0 && state.last_long_range_refresh_opportunity != 0) {
    throw std::invalid_argument("restart cadence state is inconsistent: non-zero refresh opportunity with zero field version");
  }
  if (state.long_range_restart_policy != "deterministic_rebuild") {
    throw std::invalid_argument("restart long-range policy is unsupported: " + state.long_range_restart_policy);
  }
  if (state.schema_version >= 2) {
    if (state.pm_grid_nx == 0 || state.pm_grid_ny == 0 || state.pm_grid_nz == 0) {
      throw std::invalid_argument("restart PM grid dimensions must be > 0 for schema_version >= 2");
    }
    if (state.pm_decomposition_mode.empty()) {
      throw std::invalid_argument("restart PM decomposition mode must be non-empty");
    }
    if (expected_slab_rank_count != static_cast<std::size_t>(state.world_size)) {
      throw std::invalid_argument("restart PM slab rank count must match world_size");
    }
    for (std::size_t rank = 0; rank < expected_slab_rank_count; ++rank) {
      if (!seen_slab_begin[rank] || !seen_slab_end[rank]) {
        throw std::runtime_error("restart decode missing PM slab ownership entry");
      }
      if (state.pm_slab_end_x_by_rank[rank] < state.pm_slab_begin_x_by_rank[rank]) {
        throw std::invalid_argument("restart PM slab end_x must be >= begin_x");
      }
    }
  }
  for (bool seen : seen_rank_entry) {
    if (!seen) {
      throw std::runtime_error("restart decode missing ownership entry");
    }
  }
  return state;
}

DistributedRestartCompatibilityReport evaluateDistributedRestartCompatibility(
    const DistributedRestartState& restart_state,
    const DistributedExecutionTopology& runtime_topology) {
  DistributedRestartCompatibilityReport report;
  if (restart_state.world_size != runtime_topology.world_size) {
    report.world_size_match = false;
    report.mismatch_messages.push_back(
        "world_size mismatch: restart=" + std::to_string(restart_state.world_size) +
        ", runtime=" + std::to_string(runtime_topology.world_size));
  }
  if (restart_state.pm_grid_nx != runtime_topology.pm_slab.global_nx ||
      restart_state.pm_grid_ny != runtime_topology.pm_slab.global_ny ||
      restart_state.pm_grid_nz != runtime_topology.pm_slab.global_nz) {
    report.pm_grid_shape_match = false;
    report.mismatch_messages.push_back(
        "PM grid mismatch: restart=(" + std::to_string(restart_state.pm_grid_nx) + "," +
        std::to_string(restart_state.pm_grid_ny) + "," + std::to_string(restart_state.pm_grid_nz) +
        "), runtime=(" + std::to_string(runtime_topology.pm_slab.global_nx) + "," +
        std::to_string(runtime_topology.pm_slab.global_ny) + "," +
        std::to_string(runtime_topology.pm_slab.global_nz) + ")");
  }
  if (restart_state.pm_decomposition_mode != runtime_topology.pm_decomposition_mode) {
    report.pm_decomposition_mode_match = false;
    report.mismatch_messages.push_back(
        "PM decomposition mode mismatch: restart=" + restart_state.pm_decomposition_mode +
        ", runtime=" + runtime_topology.pm_decomposition_mode);
  }
  if (restart_state.pm_update_cadence_steps == 0) {
    report.pm_cadence_steps_match = false;
    report.mismatch_messages.push_back(
        "PM cadence mismatch: restart cadence must be >= 1, got 0");
  }
  if (restart_state.last_long_range_refresh_opportunity > restart_state.gravity_kick_opportunity) {
    report.gravity_kick_state_match = false;
    report.mismatch_messages.push_back(
        "gravity kick mismatch: last refresh opportunity exceeds current kick opportunity");
  }
  if ((restart_state.long_range_field_version == 0) !=
      (restart_state.last_long_range_refresh_opportunity == 0)) {
    report.long_range_field_state_match = false;
    report.mismatch_messages.push_back(
        "long-range field mismatch: field version and refresh opportunity are inconsistent");
  }
  if (runtime_topology.world_rank < 0 ||
      runtime_topology.world_rank >= static_cast<int>(restart_state.pm_slab_begin_x_by_rank.size())) {
    report.pm_local_slab_match = false;
    report.mismatch_messages.push_back("runtime world_rank is outside restart PM slab ownership table");
    return report;
  }
  const std::size_t rank = static_cast<std::size_t>(runtime_topology.world_rank);
  const std::size_t begin_restart = restart_state.pm_slab_begin_x_by_rank[rank];
  const std::size_t end_restart = restart_state.pm_slab_end_x_by_rank[rank];
  if (begin_restart != runtime_topology.pm_slab.owned_x.begin_x ||
      end_restart != runtime_topology.pm_slab.owned_x.end_x) {
    report.pm_local_slab_match = false;
    report.mismatch_messages.push_back(
        "PM slab mismatch for rank " + std::to_string(runtime_topology.world_rank) +
        ": restart=[" + std::to_string(begin_restart) + "," + std::to_string(end_restart) +
        "), runtime=[" + std::to_string(runtime_topology.pm_slab.owned_x.begin_x) + "," +
        std::to_string(runtime_topology.pm_slab.owned_x.end_x) + ")");
  }
  return report;
}

PmMeshOwnershipDescriptor PmSlabLayout::ownershipDescriptor(
    std::uint64_t decomposition_epoch,
    std::string decomposition_mode) const {
  PmMeshOwnershipDescriptor descriptor{
      .decomposition_mode = std::move(decomposition_mode),
      .owner_rank = world_rank,
      .decomposition_epoch = decomposition_epoch,
      .global_nx = global_nx,
      .global_ny = global_ny,
      .global_nz = global_nz,
      .begin_x = owned_x.begin_x,
      .end_x = owned_x.end_x,
  };
  validatePmMeshOwnershipDescriptor(descriptor);
  return descriptor;
}

void validatePmMeshOwnershipDescriptor(const PmMeshOwnershipDescriptor& descriptor) {
  if (descriptor.decomposition_mode != "slab" && descriptor.decomposition_mode != "pencil") {
    throw std::invalid_argument("PM mesh ownership descriptor has unsupported decomposition mode");
  }
  if (descriptor.owner_rank < 0) {
    throw std::invalid_argument("PM mesh owner_rank must be non-negative");
  }
  if (descriptor.global_nx == 0 || descriptor.global_ny == 0 || descriptor.global_nz == 0) {
    throw std::invalid_argument("PM mesh descriptor global dimensions must be positive");
  }
  if (descriptor.begin_x > descriptor.end_x || descriptor.end_x > descriptor.global_nx) {
    throw std::invalid_argument("PM mesh x ownership range is invalid");
  }
}

void validateTreePseudoParticleDescriptor(const TreePseudoParticleDescriptor& descriptor) {
  if (descriptor.source_rank < 0) {
    throw std::invalid_argument("tree pseudo-particle source_rank must be non-negative");
  }
  if (!descriptor.derived_not_authoritative) {
    throw std::invalid_argument("tree pseudo-particles must be marked as derived non-authoritative state");
  }
}

void validateTreePseudoParticlePacket(const TreePseudoParticlePacket& packet) {
  validateTreePseudoParticleDescriptor(packet.descriptor);
  if (packet.source_count == 0 && packet.mass_code != 0.0) {
    throw std::invalid_argument("empty tree pseudo-particle packet cannot carry non-zero mass");
  }
  if (packet.mass_code < 0.0 || !std::isfinite(packet.mass_code)) {
    throw std::invalid_argument("tree pseudo-particle mass must be finite and non-negative");
  }
  const std::array values{packet.center_x_comoving, packet.center_y_comoving, packet.center_z_comoving,
                          packet.min_x_comoving, packet.max_x_comoving, packet.min_y_comoving,
                          packet.max_y_comoving, packet.min_z_comoving, packet.max_z_comoving};
  for (const double value : values) {
    if (!std::isfinite(value)) {
      throw std::invalid_argument("tree pseudo-particle packet contains non-finite geometry");
    }
  }
  if (packet.min_x_comoving > packet.max_x_comoving || packet.min_y_comoving > packet.max_y_comoving ||
      packet.min_z_comoving > packet.max_z_comoving) {
    throw std::invalid_argument("tree pseudo-particle packet bounds are invalid");
  }
  if (packet.child_count > 8U) {
    throw std::invalid_argument("tree pseudo-particle packet child count is invalid");
  }
  if (packet.is_leaf > 1U) {
    throw std::invalid_argument("tree pseudo-particle packet leaf flag is invalid");
  }
}

void validateHydroGhostCellDescriptor(const HydroGhostCellDescriptor& descriptor) {
  if (descriptor.owner_rank < 0 || descriptor.consumer_rank < 0) {
    throw std::invalid_argument("hydro ghost cell ranks must be non-negative");
  }
  if (descriptor.owner_rank == descriptor.consumer_rank) {
    throw std::invalid_argument("hydro ghost cell must be consumed on a non-owner rank");
  }
  if (!descriptor.boundary_state_only) {
    throw std::invalid_argument("hydro ghost cells are boundary exchange state, not authoritative conserved truth");
  }
}

void validateAmrPatchExchangeDescriptor(const AmrPatchExchangeDescriptor& descriptor) {
  if (descriptor.owner_rank < 0 || descriptor.peer_rank < 0) {
    throw std::invalid_argument("AMR patch exchange ranks must be non-negative");
  }
  if (!descriptor.metadata_only && descriptor.owner_rank != descriptor.peer_rank) {
    throw std::invalid_argument("remote AMR patch exchange cannot mutate authoritative patch metadata");
  }
}

void validateAmrPatchPayloadRecord(const AmrPatchPayloadRecord& record) {
  if (record.owner_rank < 0) {
    throw std::invalid_argument("AMR patch payload owner_rank must be non-negative");
  }
  if (record.patch_id == 0) {
    throw std::invalid_argument("AMR patch payload patch_id must be non-zero");
  }
  if (record.cell_count == 0) {
    throw std::invalid_argument("AMR patch payload cannot describe an empty patch");
  }
  if (!std::isfinite(record.cell_mass_sum_code) || !std::isfinite(record.gas_internal_energy_sum_code)) {
    throw std::invalid_argument("AMR patch payload contains non-finite cell sums");
  }
}

void validateAmrPatchCellPayloadRecord(const AmrPatchCellPayloadRecord& record) {
  if (record.owner_rank < 0) {
    throw std::invalid_argument("AMR patch cell payload owner_rank must be non-negative");
  }
  if (record.patch_id == 0) {
    throw std::invalid_argument("AMR patch cell payload patch_id must be non-zero");
  }
  if (record.gas_cell_id == 0) {
    throw std::invalid_argument("AMR patch cell payload must carry stable gas-cell identity");
  }
  if (!std::isfinite(record.center_x_comoving) || !std::isfinite(record.center_y_comoving) ||
      !std::isfinite(record.center_z_comoving) || !std::isfinite(record.mass_code) ||
      !std::isfinite(record.density_code) || !std::isfinite(record.pressure_code) ||
      !std::isfinite(record.internal_energy_code) || !std::isfinite(record.temperature_code) ||
      !std::isfinite(record.sound_speed_code)) {
    throw std::invalid_argument("AMR patch cell payload contains non-finite state");
  }
}

void validateHydroConservativeFluxCorrectionRecord(const HydroConservativeFluxCorrectionRecord& record) {
  if (record.gas_cell_id == 0) {
    throw std::invalid_argument("hydro conservative flux correction requires a non-zero gas_cell_id");
  }
  if (record.source_rank < 0 || record.owner_rank < 0) {
    throw std::invalid_argument("hydro conservative flux correction ranks must be non-negative");
  }
  if (!std::isfinite(record.delta_mass_density_comoving) ||
      !std::isfinite(record.delta_momentum_density_x_comoving) ||
      !std::isfinite(record.delta_momentum_density_y_comoving) ||
      !std::isfinite(record.delta_momentum_density_z_comoving) ||
      !std::isfinite(record.delta_total_energy_density_comoving)) {
    throw std::invalid_argument("hydro conservative flux correction contains non-finite state");
  }
}

void recordDistributedProfiling(
    core::ProfilerSession* profiler,
    const LoadBalanceMetrics& metrics,
    std::uint64_t ghost_exchange_send_bytes,
    std::uint64_t ghost_exchange_recv_bytes) {
  if (profiler == nullptr) {
    return;
  }

  profiler->counters().setCount("parallel.ghost_exchange_send_bytes", ghost_exchange_send_bytes);
  profiler->counters().setCount("parallel.ghost_exchange_recv_bytes", ghost_exchange_recv_bytes);
  profiler->counters().setCount("parallel.total_memory_bytes", metrics.total_memory_bytes);

  const std::uint64_t imbalance_ppm = (metrics.weighted_imbalance_ratio <= 0.0)
                                          ? 0ULL
                                          : static_cast<std::uint64_t>(std::llround(metrics.weighted_imbalance_ratio * 1.0e6));
  profiler->counters().setCount("parallel.weighted_imbalance_ratio_ppm", imbalance_ppm);
  auto record_component_total = [&](std::string_view name, const std::vector<double>& values) {
    const double total = std::accumulate(values.begin(), values.end(), 0.0);
    profiler->counters().setCount(std::string(name), static_cast<std::uint64_t>(std::llround(std::max(0.0, total))));
  };
  record_component_total("parallel.weight_component_particle_count", metrics.particle_count_cost_by_rank);
  record_component_total("parallel.weight_component_gas_cell", metrics.gas_cell_cost_by_rank);
  record_component_total("parallel.weight_component_tree_interaction", metrics.tree_interaction_cost_by_rank);
  record_component_total("parallel.weight_component_pm_mesh", metrics.pm_mesh_cost_by_rank);
  record_component_total("parallel.weight_component_amr_patch", metrics.amr_patch_cost_by_rank);
  record_component_total("parallel.weight_component_active_fraction", metrics.active_fraction_cost_by_rank);
  record_component_total("parallel.weight_component_memory_pressure", metrics.memory_pressure_cost_by_rank);
  record_component_total("parallel.weight_component_gpu_occupancy", metrics.gpu_occupancy_cost_by_rank);

  if (!metrics.weighted_load_by_rank.empty()) {
    const auto max_rank_it = std::max_element(metrics.weighted_load_by_rank.begin(), metrics.weighted_load_by_rank.end());
    const std::size_t max_rank = static_cast<std::size_t>(std::distance(metrics.weighted_load_by_rank.begin(), max_rank_it));
    profiler->counters().setCount("parallel.max_weighted_load_rank", static_cast<std::uint64_t>(max_rank));
    profiler->counters().setCount(
        "parallel.max_weighted_load",
        static_cast<std::uint64_t>(std::llround(std::max(0.0, *max_rank_it))));

    struct ComponentView {
      std::string_view name;
      const std::vector<double>* values;
    };
    const std::array<ComponentView, 9> components{{
        {"particle_count", &metrics.particle_count_cost_by_rank},
        {"gas_cell", &metrics.gas_cell_cost_by_rank},
        {"tree_interaction", &metrics.tree_interaction_cost_by_rank},
        {"pm_mesh", &metrics.pm_mesh_cost_by_rank},
        {"amr_patch", &metrics.amr_patch_cost_by_rank},
        {"active_fraction", &metrics.active_fraction_cost_by_rank},
        {"memory_pressure", &metrics.memory_pressure_cost_by_rank},
        {"gpu_occupancy", &metrics.gpu_occupancy_cost_by_rank},
        {"generic_work", &metrics.generic_work_cost_by_rank},
    }};
    std::string_view dominant_component = "none";
    double dominant_value = 0.0;
    for (const ComponentView& component : components) {
      if (component.values->size() <= max_rank) {
        continue;
      }
      const double value = (*component.values)[max_rank];
      if (value > dominant_value) {
        dominant_value = value;
        dominant_component = component.name;
      }
    }
    profiler->counters().setCount(
        "parallel.max_rank_dominant_component_cost",
        static_cast<std::uint64_t>(std::llround(std::max(0.0, dominant_value))));
    profiler->recordEvent(core::RuntimeEvent{
        .event_kind = "parallel.decomposition.hotspot",
        .severity = metrics.weighted_imbalance_ratio > 1.25 ? core::RuntimeEventSeverity::kWarning
                                                            : core::RuntimeEventSeverity::kInfo,
        .subsystem = "parallel.domain_decomposition",
        .message = "domain decomposition load hotspot attribution",
        .payload = {{"max_rank", std::to_string(max_rank)},
                    {"max_rank_load", std::to_string(*max_rank_it)},
                    {"mean_load", std::to_string(metrics.mean_weighted_load)},
                    {"weighted_imbalance_ratio", std::to_string(metrics.weighted_imbalance_ratio)},
                    {"memory_imbalance_ratio", std::to_string(metrics.memory_imbalance_ratio)},
                    {"dominant_component", std::string(dominant_component)},
                    {"dominant_component_cost", std::to_string(dominant_value)}}});
  }

  const std::uint64_t bytes_moved = ghost_exchange_send_bytes + ghost_exchange_recv_bytes;
  profiler->addBytesMoved(bytes_moved);
}

MpiContext::MpiContext() {
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  int initialized = 0;
  MPI_Initialized(&initialized);
  if (initialized != 0) {
    m_is_enabled = true;
    MPI_Comm_size(MPI_COMM_WORLD, &m_world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &m_world_rank);
  }
#endif
}

MpiContext::MpiContext(bool is_enabled, int world_size, int world_rank)
    : m_is_enabled(is_enabled), m_world_size(world_size), m_world_rank(world_rank) {
  if (world_size <= 0) {
    throw std::invalid_argument("MpiContext world_size must be positive");
  }
  if (world_rank < 0 || world_rank >= world_size) {
    throw std::invalid_argument("MpiContext world_rank must be within [0, world_size)");
  }
}

bool MpiContext::isEnabled() const noexcept { return m_is_enabled; }

bool MpiContext::isRoot() const noexcept { return m_world_rank == 0; }

int MpiContext::worldSize() const noexcept { return m_world_size; }

int MpiContext::worldRank() const noexcept { return m_world_rank; }

void MpiContext::validateExpectedWorldSizeOrThrow(int expected_world_size) const {
  if (expected_world_size <= 0) {
    throw std::invalid_argument("expected_world_size must be positive");
  }
  if (expected_world_size != m_world_size) {
    throw std::runtime_error(
        "parallel.mpi_ranks_expected does not match runtime world size: expected=" +
        std::to_string(expected_world_size) + ", runtime=" + std::to_string(m_world_size));
  }
}

double MpiContext::allreduceSumDouble(double local_value) const {
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  if (m_is_enabled) {
    double global = 0.0;
    MPI_Allreduce(&local_value, &global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return global;
  }
#endif
  return local_value;
}

std::uint64_t MpiContext::allreduceSumUint64(std::uint64_t local_value) const {
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  if (m_is_enabled) {
    std::uint64_t global = 0;
    MPI_Allreduce(&local_value, &global, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
    return global;
  }
#endif
  return local_value;
}

std::uint64_t MpiContext::allreduceXorUint64(std::uint64_t local_value) const {
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  if (m_is_enabled) {
    std::uint64_t global = 0;
    MPI_Allreduce(&local_value, &global, 1, MPI_UINT64_T, MPI_BXOR, MPI_COMM_WORLD);
    return global;
  }
#endif
  return local_value;
}

std::vector<TreePseudoParticlePacket> executeBlockingTreePseudoParticleExchange(
    const MpiContext& mpi_context,
    const TreePseudoParticlePacket& local_packet) {
  validateTreePseudoParticlePacket(local_packet);
  if (local_packet.descriptor.source_rank != mpi_context.worldRank()) {
    throw std::invalid_argument("tree pseudo-particle packet source rank does not match MPI context");
  }
  if (!mpi_context.isEnabled()) {
    return {local_packet};
  }
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  std::vector<TreePseudoParticlePacket> packets(static_cast<std::size_t>(mpi_context.worldSize()));
  MPI_Allgather(
      const_cast<TreePseudoParticlePacket*>(&local_packet),
      static_cast<int>(sizeof(TreePseudoParticlePacket)),
      MPI_BYTE,
      packets.data(),
      static_cast<int>(sizeof(TreePseudoParticlePacket)),
      MPI_BYTE,
      MPI_COMM_WORLD);
  for (int rank = 0; rank < mpi_context.worldSize(); ++rank) {
    const TreePseudoParticlePacket& packet = packets[static_cast<std::size_t>(rank)];
    validateTreePseudoParticlePacket(packet);
    if (packet.descriptor.source_rank != rank) {
      throw std::runtime_error("tree pseudo-particle exchange returned a packet with mismatched source rank");
    }
  }
  return packets;
#else
  throw std::runtime_error("tree pseudo-particle exchange requires MPI support when MPI context is enabled");
#endif
}

std::vector<TreePseudoParticlePacket> executeBlockingTreePseudoParticleHierarchyExchange(
    const MpiContext& mpi_context,
    std::span<const TreePseudoParticlePacket> local_packets,
    std::uint64_t exchange_sequence) {
  (void)exchange_sequence;
  for (const TreePseudoParticlePacket& packet : local_packets) {
    validateTreePseudoParticlePacket(packet);
    if (packet.descriptor.source_rank != mpi_context.worldRank()) {
      throw std::invalid_argument("tree pseudo hierarchy packet source rank does not match MPI context");
    }
  }
  if (!mpi_context.isEnabled()) {
    return std::vector<TreePseudoParticlePacket>(local_packets.begin(), local_packets.end());
  }
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  const int world_size = mpi_context.worldSize();
  const std::uint64_t local_count = static_cast<std::uint64_t>(local_packets.size());
  std::vector<std::uint64_t> counts64(static_cast<std::size_t>(world_size), 0U);
  MPI_Allgather(
      const_cast<std::uint64_t*>(&local_count),
      1,
      MPI_UINT64_T,
      counts64.data(),
      1,
      MPI_UINT64_T,
      MPI_COMM_WORLD);
  std::vector<int> recv_counts(static_cast<std::size_t>(world_size), 0);
  std::vector<int> recv_displs(static_cast<std::size_t>(world_size), 0);
  std::uint64_t total_packets = 0;
  for (int rank = 0; rank < world_size; ++rank) {
    const std::uint64_t bytes = counts64[static_cast<std::size_t>(rank)] * sizeof(TreePseudoParticlePacket);
    if (bytes > static_cast<std::uint64_t>(std::numeric_limits<int>::max())) {
      throw std::overflow_error("tree pseudo hierarchy exchange byte count exceeds MPI int limit");
    }
    recv_counts[static_cast<std::size_t>(rank)] = static_cast<int>(bytes);
    if (rank > 0) {
      recv_displs[static_cast<std::size_t>(rank)] =
          recv_displs[static_cast<std::size_t>(rank - 1)] + recv_counts[static_cast<std::size_t>(rank - 1)];
    }
    total_packets += counts64[static_cast<std::size_t>(rank)];
  }
  const std::uint64_t total_bytes = total_packets * sizeof(TreePseudoParticlePacket);
  if (total_bytes > static_cast<std::uint64_t>(std::numeric_limits<int>::max())) {
    throw std::overflow_error("tree pseudo hierarchy exchange total byte count exceeds MPI int limit");
  }
  std::vector<TreePseudoParticlePacket> result(static_cast<std::size_t>(total_packets));
  MPI_Allgatherv(
      const_cast<TreePseudoParticlePacket*>(local_packets.data()),
      static_cast<int>(local_packets.size() * sizeof(TreePseudoParticlePacket)),
      MPI_BYTE,
      result.data(),
      recv_counts.data(),
      recv_displs.data(),
      MPI_BYTE,
      MPI_COMM_WORLD);
  std::vector<std::uint64_t> per_rank_count(static_cast<std::size_t>(world_size), 0U);
  for (const TreePseudoParticlePacket& packet : result) {
    validateTreePseudoParticlePacket(packet);
    if (packet.descriptor.source_rank < 0 || packet.descriptor.source_rank >= world_size) {
      throw std::runtime_error("tree pseudo hierarchy exchange returned packet with invalid source rank");
    }
    ++per_rank_count[static_cast<std::size_t>(packet.descriptor.source_rank)];
  }
  for (int rank = 0; rank < world_size; ++rank) {
    if (per_rank_count[static_cast<std::size_t>(rank)] != counts64[static_cast<std::size_t>(rank)]) {
      throw std::runtime_error("tree pseudo hierarchy exchange source-rank coverage mismatch");
    }
  }
  return result;
#else
  throw std::runtime_error("tree pseudo hierarchy exchange requires MPI support when MPI context is enabled");
#endif
}

std::vector<AmrPatchPayloadRecord> executeBlockingAmrPatchPayloadExchange(
    const MpiContext& mpi_context,
    std::span<const AmrPatchPayloadRecord> local_records,
    std::uint64_t exchange_sequence) {
  (void)exchange_sequence;
  for (const AmrPatchPayloadRecord& record : local_records) {
    validateAmrPatchPayloadRecord(record);
    if (record.owner_rank != mpi_context.worldRank()) {
      throw std::invalid_argument("AMR patch payload exchange received a record not owned by this rank");
    }
  }
  if (!mpi_context.isEnabled()) {
    return std::vector<AmrPatchPayloadRecord>(local_records.begin(), local_records.end());
  }
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  const int world_size = mpi_context.worldSize();
  const std::uint64_t local_count = static_cast<std::uint64_t>(local_records.size());
  std::vector<std::uint64_t> counts64(static_cast<std::size_t>(world_size), 0U);
  MPI_Allgather(
      const_cast<std::uint64_t*>(&local_count),
      1,
      MPI_UINT64_T,
      counts64.data(),
      1,
      MPI_UINT64_T,
      MPI_COMM_WORLD);
  std::vector<int> recv_counts(static_cast<std::size_t>(world_size), 0);
  std::vector<int> recv_displs(static_cast<std::size_t>(world_size), 0);
  std::uint64_t total_records64 = 0;
  for (int rank = 0; rank < world_size; ++rank) {
    if (counts64[static_cast<std::size_t>(rank)] > static_cast<std::uint64_t>(std::numeric_limits<int>::max())) {
      throw std::overflow_error("AMR patch payload exchange count exceeds MPI int limit");
    }
    recv_counts[static_cast<std::size_t>(rank)] = static_cast<int>(counts64[static_cast<std::size_t>(rank)] * sizeof(AmrPatchPayloadRecord));
    if (rank > 0) {
      recv_displs[static_cast<std::size_t>(rank)] = recv_displs[static_cast<std::size_t>(rank - 1)] + recv_counts[static_cast<std::size_t>(rank - 1)];
    }
    total_records64 += counts64[static_cast<std::size_t>(rank)];
  }
  if (total_records64 > static_cast<std::uint64_t>(std::numeric_limits<int>::max())) {
    throw std::overflow_error("AMR patch payload exchange total count exceeds MPI int limit");
  }
  std::vector<AmrPatchPayloadRecord> result(static_cast<std::size_t>(total_records64));
  MPI_Allgatherv(
      const_cast<AmrPatchPayloadRecord*>(local_records.data()),
      static_cast<int>(local_records.size() * sizeof(AmrPatchPayloadRecord)),
      MPI_BYTE,
      result.data(),
      recv_counts.data(),
      recv_displs.data(),
      MPI_BYTE,
      MPI_COMM_WORLD);
  for (const AmrPatchPayloadRecord& record : result) {
    validateAmrPatchPayloadRecord(record);
  }
  return result;
#else
  throw std::runtime_error("AMR patch payload exchange requires MPI support when MPI context is enabled");
#endif
}


std::vector<AmrPatchCellPayloadRecord> executeBlockingAmrPatchCellPayloadExchange(
    const MpiContext& mpi_context,
    std::span<const AmrPatchCellPayloadRecord> local_records,
    std::uint64_t exchange_sequence) {
  (void)exchange_sequence;
  for (const AmrPatchCellPayloadRecord& record : local_records) {
    validateAmrPatchCellPayloadRecord(record);
    if (record.owner_rank != mpi_context.worldRank()) {
      throw std::invalid_argument("AMR patch cell payload exchange received a cell not owned by this rank");
    }
  }
  if (!mpi_context.isEnabled()) {
    return std::vector<AmrPatchCellPayloadRecord>(local_records.begin(), local_records.end());
  }
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  const int world_size = mpi_context.worldSize();
  const std::uint64_t local_count = static_cast<std::uint64_t>(local_records.size());
  std::vector<std::uint64_t> counts64(static_cast<std::size_t>(world_size), 0U);
  MPI_Allgather(
      const_cast<std::uint64_t*>(&local_count),
      1,
      MPI_UINT64_T,
      counts64.data(),
      1,
      MPI_UINT64_T,
      MPI_COMM_WORLD);
  std::vector<int> recv_counts(static_cast<std::size_t>(world_size), 0);
  std::vector<int> recv_displs(static_cast<std::size_t>(world_size), 0);
  std::uint64_t total_records64 = 0;
  for (int rank = 0; rank < world_size; ++rank) {
    const std::uint64_t bytes64 = counts64[static_cast<std::size_t>(rank)] * sizeof(AmrPatchCellPayloadRecord);
    if (bytes64 > static_cast<std::uint64_t>(std::numeric_limits<int>::max())) {
      throw std::overflow_error("AMR patch cell payload exchange byte count exceeds MPI int limit");
    }
    recv_counts[static_cast<std::size_t>(rank)] = static_cast<int>(bytes64);
    if (rank > 0) {
      recv_displs[static_cast<std::size_t>(rank)] =
          recv_displs[static_cast<std::size_t>(rank - 1)] + recv_counts[static_cast<std::size_t>(rank - 1)];
    }
    total_records64 += counts64[static_cast<std::size_t>(rank)];
  }
  const std::uint64_t total_bytes64 = total_records64 * sizeof(AmrPatchCellPayloadRecord);
  if (total_bytes64 > static_cast<std::uint64_t>(std::numeric_limits<int>::max())) {
    throw std::overflow_error("AMR patch cell payload exchange total byte count exceeds MPI int limit");
  }
  std::vector<AmrPatchCellPayloadRecord> result(static_cast<std::size_t>(total_records64));
  MPI_Allgatherv(
      const_cast<AmrPatchCellPayloadRecord*>(local_records.data()),
      static_cast<int>(local_records.size() * sizeof(AmrPatchCellPayloadRecord)),
      MPI_BYTE,
      result.data(),
      recv_counts.data(),
      recv_displs.data(),
      MPI_BYTE,
      MPI_COMM_WORLD);
  for (const AmrPatchCellPayloadRecord& record : result) {
    validateAmrPatchCellPayloadRecord(record);
  }
  return result;
#else
  throw std::runtime_error("AMR patch cell payload exchange requires MPI support when MPI context is enabled");
#endif
}

std::vector<HydroConservativeFluxCorrectionRecord> executeBlockingHydroConservativeFluxCorrectionExchange(
    const MpiContext& mpi_context,
    std::span<const HydroConservativeFluxCorrectionRecord> local_records,
    std::uint64_t exchange_sequence) {
  (void)exchange_sequence;
  for (const HydroConservativeFluxCorrectionRecord& record : local_records) {
    validateHydroConservativeFluxCorrectionRecord(record);
    if (record.source_rank != mpi_context.worldRank()) {
      throw std::invalid_argument("hydro conservative flux correction source rank does not match MPI context");
    }
  }
  if (!mpi_context.isEnabled()) {
    return std::vector<HydroConservativeFluxCorrectionRecord>(local_records.begin(), local_records.end());
  }
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  const int world_size = mpi_context.worldSize();
  const std::uint64_t local_count = static_cast<std::uint64_t>(local_records.size());
  std::vector<std::uint64_t> counts64(static_cast<std::size_t>(world_size), 0U);
  MPI_Allgather(
      const_cast<std::uint64_t*>(&local_count),
      1,
      MPI_UINT64_T,
      counts64.data(),
      1,
      MPI_UINT64_T,
      MPI_COMM_WORLD);
  std::vector<int> recv_counts(static_cast<std::size_t>(world_size), 0);
  std::vector<int> recv_displs(static_cast<std::size_t>(world_size), 0);
  std::uint64_t total_records64 = 0;
  for (int rank = 0; rank < world_size; ++rank) {
    const std::uint64_t bytes64 = counts64[static_cast<std::size_t>(rank)] * sizeof(HydroConservativeFluxCorrectionRecord);
    if (bytes64 > static_cast<std::uint64_t>(std::numeric_limits<int>::max())) {
      throw std::overflow_error("hydro conservative flux correction exchange byte count exceeds MPI int limit");
    }
    recv_counts[static_cast<std::size_t>(rank)] = static_cast<int>(bytes64);
    if (rank > 0) {
      recv_displs[static_cast<std::size_t>(rank)] =
          recv_displs[static_cast<std::size_t>(rank - 1)] + recv_counts[static_cast<std::size_t>(rank - 1)];
    }
    total_records64 += counts64[static_cast<std::size_t>(rank)];
  }
  const std::uint64_t total_bytes64 = total_records64 * sizeof(HydroConservativeFluxCorrectionRecord);
  if (total_bytes64 > static_cast<std::uint64_t>(std::numeric_limits<int>::max())) {
    throw std::overflow_error("hydro conservative flux correction exchange total byte count exceeds MPI int limit");
  }
  std::vector<HydroConservativeFluxCorrectionRecord> result(static_cast<std::size_t>(total_records64));
  MPI_Allgatherv(
      const_cast<HydroConservativeFluxCorrectionRecord*>(local_records.data()),
      static_cast<int>(local_records.size() * sizeof(HydroConservativeFluxCorrectionRecord)),
      MPI_BYTE,
      result.data(),
      recv_counts.data(),
      recv_displs.data(),
      MPI_BYTE,
      MPI_COMM_WORLD);
  for (const HydroConservativeFluxCorrectionRecord& record : result) {
    validateHydroConservativeFluxCorrectionRecord(record);
    if (record.source_rank < 0 || record.source_rank >= world_size ||
        record.owner_rank < 0 || record.owner_rank >= world_size) {
      throw std::runtime_error("hydro conservative flux correction exchange returned invalid rank metadata");
    }
  }
  return result;
#else
  throw std::runtime_error("hydro conservative flux correction exchange requires MPI support when MPI context is enabled");
#endif
}

PmSlabHaloExchangeResult executeBlockingPmSlabHaloExchange(
    const MpiContext& mpi_context,
    const PmSlabLayout& layout,
    std::span<const double> local_scalar_field,
    std::size_t halo_depth_x,
    bool periodic_x,
    std::uint64_t exchange_sequence) {
  if (!layout.isValid()) {
    throw std::invalid_argument("PM slab halo exchange requires a valid slab layout");
  }
  if (layout.world_size != mpi_context.worldSize() || layout.world_rank != mpi_context.worldRank()) {
    throw std::invalid_argument("PM slab halo exchange layout world metadata must match MPI context");
  }
  if (halo_depth_x == 0 || layout.world_size == 1) {
    return {};
  }
  const std::size_t plane_size = layout.global_ny * layout.global_nz;
  if (local_scalar_field.size() != layout.localCellCount()) {
    throw std::invalid_argument("PM slab halo exchange field size does not match local slab cell count");
  }
  const std::size_t depth = std::min(halo_depth_x, layout.local_nx());
  PmSlabHaloExchangeResult result;
  result.halo_depth_x = depth;
  const int left_peer = (layout.world_rank > 0) ? (layout.world_rank - 1) : (periodic_x ? layout.world_size - 1 : -1);
  const int right_peer = (layout.world_rank + 1 < layout.world_size) ? (layout.world_rank + 1) : (periodic_x ? 0 : -1);
  result.left_peer_rank = left_peer;
  result.right_peer_rank = right_peer;
  result.left_halo.assign(depth * plane_size, 0.0);
  result.right_halo.assign(depth * plane_size, 0.0);
  if (!mpi_context.isEnabled()) {
    throw std::runtime_error("PM slab halo exchange requires MPI for distributed layouts");
  }
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  std::vector<double> send_left(depth * plane_size, 0.0);
  std::vector<double> send_right(depth * plane_size, 0.0);
  std::copy_n(local_scalar_field.begin(), send_left.size(), send_left.begin());
  std::copy_n(local_scalar_field.end() - static_cast<std::ptrdiff_t>(send_right.size()), send_right.size(), send_right.begin());
  constexpr int k_pm_halo_tag_base = 8810;
  constexpr int k_send_left_side = 0;
  constexpr int k_send_right_side = 1;
  const auto edge_index = [&](int peer) {
    const int local = mpi_context.worldRank();
    return (std::abs(local - peer) == 1) ? std::min(local, peer) : (layout.world_size - 1);
  };
  const auto side_tag = [&](int peer, int side) {
    return k_pm_halo_tag_base + edge_index(peer) * 2 + side;
  };
  const auto sendrecv_planes = [&](int peer,
                                   const std::vector<double>& send,
                                   std::vector<double>& recv,
                                   int send_side,
                                   int recv_side) {
    if (peer < 0) {
      return;
    }
    MPI_Sendrecv(
        const_cast<double*>(send.data()),
        static_cast<int>(send.size()),
        MPI_DOUBLE,
        peer,
        ghostExchangeSequencedTag(side_tag(peer, send_side), mpi_context.worldRank(), peer, exchange_sequence),
        recv.data(),
        static_cast<int>(recv.size()),
        MPI_DOUBLE,
        peer,
        ghostExchangeSequencedTag(side_tag(peer, recv_side), mpi_context.worldRank(), peer, exchange_sequence),
        MPI_COMM_WORLD,
        MPI_STATUS_IGNORE);
    result.sent_bytes += static_cast<std::uint64_t>(send.size() * sizeof(double));
    result.received_bytes += static_cast<std::uint64_t>(recv.size() * sizeof(double));
  };
  // Send the left edge to the left peer and receive that peer's right edge as our left halo.
  sendrecv_planes(left_peer, send_left, result.left_halo, k_send_left_side, k_send_right_side);
  // Send the right edge to the right peer and receive that peer's left edge as our right halo.
  sendrecv_planes(right_peer, send_right, result.right_halo, k_send_right_side, k_send_left_side);
  return result;
#else
  throw std::runtime_error("PM slab halo exchange requires MPI support when MPI context is enabled");
#endif
}

GhostRefreshCommitReport commitBlockingGhostRefreshResult(
    GhostExchangeBufferSoA& ghost_storage,
    std::span<const LocalGhostDescriptor> local_ghost_descriptors,
    const GhostExchangePlan& plan,
    const BlockingGhostExchangeResult& result,
    const GhostLayerEpoch& expected_epoch) {
  validateGhostExchangePlan(plan);
  if (!plan.epoch.matches(expected_epoch)) {
    throw std::invalid_argument("commitBlockingGhostRefreshResult: ghost exchange plan epoch is stale");
  }
  if (!ghost_storage.isConsistent() || !result.received_ghosts.isConsistent()) {
    throw std::invalid_argument("commitBlockingGhostRefreshResult: ghost storage and result payloads must be component-consistent");
  }
  if (ghost_storage.size() < local_ghost_descriptors.size()) {
    throw std::invalid_argument("commitBlockingGhostRefreshResult: ghost storage must expose one slot per local descriptor");
  }
  if (!result.received_ghosts.epoch.matches(expected_epoch)) {
    throw std::invalid_argument("commitBlockingGhostRefreshResult: received ghost payload epoch is stale");
  }

  std::size_t expected_count = 0;
  for (const auto& indices : plan.recv_local_indices_by_neighbor) {
    expected_count += indices.size();
  }
  if (result.received_ghosts.size() != expected_count) {
    throw std::invalid_argument("commitBlockingGhostRefreshResult: received payload count does not match plan receive slots");
  }

  GhostRefreshCommitReport report;
  std::size_t result_row = 0;
  for (std::size_t slot = 0; slot < plan.recv_local_indices_by_neighbor.size(); ++slot) {
    for (const std::uint32_t local_index : plan.recv_local_indices_by_neighbor[slot]) {
      if (local_index >= local_ghost_descriptors.size() || local_index >= ghost_storage.size()) {
        throw std::out_of_range("commitBlockingGhostRefreshResult: receive slot index out of range");
      }
      const LocalGhostDescriptor descriptor = local_ghost_descriptors[local_index];
      if (descriptor.residency != LocalIndexResidency::kGhost || descriptor.owning_rank != plan.neighbor_ranks[slot]) {
        throw std::invalid_argument("commitBlockingGhostRefreshResult: receive slot is not a ghost owned by the exchange peer");
      }
      if (!descriptor.epoch.matches(expected_epoch)) {
        throw std::invalid_argument("commitBlockingGhostRefreshResult: local ghost descriptor is stale");
      }
      if (result.received_ghosts.entity_id[result_row] != descriptor.particle_id) {
        throw std::invalid_argument("commitBlockingGhostRefreshResult: received entity_id does not match ghost slot particle_id");
      }
      ghost_storage.entity_id[local_index] = result.received_ghosts.entity_id[result_row];
      const std::size_t storage_size = ghost_storage.size();
      resizeOptionalLaneForCommit(&ghost_storage.position_x_comoving, storage_size);
      resizeOptionalLaneForCommit(&ghost_storage.position_y_comoving, storage_size);
      resizeOptionalLaneForCommit(&ghost_storage.position_z_comoving, storage_size);
      resizeOptionalLaneForCommit(&ghost_storage.mass_code, storage_size);
      resizeOptionalLaneForCommit(&ghost_storage.density_code, storage_size);
      resizeOptionalLaneForCommit(&ghost_storage.velocity_x_code, storage_size);
      resizeOptionalLaneForCommit(&ghost_storage.velocity_y_code, storage_size);
      resizeOptionalLaneForCommit(&ghost_storage.velocity_z_code, storage_size);
      resizeOptionalLaneForCommit(&ghost_storage.pressure_code, storage_size);
      resizeOptionalLaneForCommit(&ghost_storage.internal_energy_code, storage_size);
      ghost_storage.position_x_comoving[local_index] = optionalLaneValue(result.received_ghosts.position_x_comoving, result_row);
      ghost_storage.position_y_comoving[local_index] = optionalLaneValue(result.received_ghosts.position_y_comoving, result_row);
      ghost_storage.position_z_comoving[local_index] = optionalLaneValue(result.received_ghosts.position_z_comoving, result_row);
      ghost_storage.mass_code[local_index] = optionalLaneValue(result.received_ghosts.mass_code, result_row);
      ghost_storage.density_code[local_index] = optionalLaneValue(result.received_ghosts.density_code, result_row);
      ghost_storage.velocity_x_code[local_index] = optionalLaneValue(result.received_ghosts.velocity_x_code, result_row);
      ghost_storage.velocity_y_code[local_index] = optionalLaneValue(result.received_ghosts.velocity_y_code, result_row);
      ghost_storage.velocity_z_code[local_index] = optionalLaneValue(result.received_ghosts.velocity_z_code, result_row);
      ghost_storage.pressure_code[local_index] = optionalLaneValue(result.received_ghosts.pressure_code, result_row);
      ghost_storage.internal_energy_code[local_index] = optionalLaneValue(result.received_ghosts.internal_energy_code, result_row);
      ++result_row;
      ++report.updated_ghost_slots;
    }
  }
  ghost_storage.epoch = expected_epoch;
  report.committed_payload_bytes = static_cast<std::uint64_t>(report.updated_ghost_slots) *
      static_cast<std::uint64_t>(ghostExchangeRecordBytes());
  return report;
}

void invalidateGhostCache(GhostCacheLifecycle& lifecycle, const GhostLayerEpoch& next_epoch) {
  lifecycle.epoch = next_epoch;
  lifecycle.valid = false;
  ++lifecycle.invalidation_count;
}

void markGhostCacheCommitted(GhostCacheLifecycle& lifecycle, const GhostLayerEpoch& committed_epoch) {
  lifecycle.epoch = committed_epoch;
  lifecycle.valid = true;
  ++lifecycle.refresh_count;
}

void requireValidGhostCache(
    const GhostCacheLifecycle& lifecycle,
    const GhostLayerEpoch& expected_epoch,
    std::string_view caller) {
  if (!lifecycle.valid || !lifecycle.epoch.matches(expected_epoch)) {
    throw std::runtime_error(std::string(caller) + ": stale or invalid ghost cache used by solver");
  }
}


BlockingGhostRefreshExchange executeBlockingGhostRefreshExchangeFromDescriptors(
    const MpiContext& mpi_context,
    std::span<const LocalGhostDescriptor> local_ghost_descriptors,
    const GhostExchangeBufferSoA& authoritative_local_state,
    const GhostLayerEpoch& expected_epoch) {
  if (!authoritative_local_state.isConsistent()) {
    throw std::invalid_argument(
        "executeBlockingGhostRefreshExchangeFromDescriptors: authoritative payload state is inconsistent");
  }
  if (!authoritative_local_state.epoch.matches(expected_epoch)) {
    throw std::invalid_argument(
        "executeBlockingGhostRefreshExchangeFromDescriptors: authoritative payload epoch is stale");
  }
  if (authoritative_local_state.size() < local_ghost_descriptors.size()) {
    throw std::invalid_argument(
        "executeBlockingGhostRefreshExchangeFromDescriptors: payload state must expose one row per local descriptor");
  }

  const int world_rank = mpi_context.worldRank();
  const int world_size = mpi_context.worldSize();
  if (world_rank < 0 || world_size <= 0 || world_rank >= world_size) {
    throw std::invalid_argument("executeBlockingGhostRefreshExchangeFromDescriptors: invalid MPI context");
  }

  std::unordered_map<std::uint64_t, std::uint32_t> owned_index_by_particle_id;
  owned_index_by_particle_id.reserve(local_ghost_descriptors.size());
  std::vector<std::vector<std::uint64_t>> requested_particle_ids_by_rank(static_cast<std::size_t>(world_size));
  std::vector<std::vector<std::uint32_t>> recv_indices_by_rank(static_cast<std::size_t>(world_size));

  for (std::uint32_t local_index = 0; local_index < local_ghost_descriptors.size(); ++local_index) {
    const LocalGhostDescriptor descriptor = local_ghost_descriptors[local_index];
    if (!descriptor.epoch.matches(expected_epoch)) {
      throw std::invalid_argument(
          "executeBlockingGhostRefreshExchangeFromDescriptors: stale local ghost descriptor epoch");
    }
    if (descriptor.owning_rank < 0 || descriptor.owning_rank >= world_size) {
      throw std::invalid_argument(
          "executeBlockingGhostRefreshExchangeFromDescriptors: descriptor owning_rank outside MPI world");
    }
    if (authoritative_local_state.entity_id[local_index] != descriptor.particle_id) {
      throw std::invalid_argument(
          "executeBlockingGhostRefreshExchangeFromDescriptors: payload entity_id does not match descriptor particle_id");
    }
    if (descriptor.residency == LocalIndexResidency::kOwned) {
      if (descriptor.owning_rank != world_rank) {
        throw std::invalid_argument(
            "executeBlockingGhostRefreshExchangeFromDescriptors: owned descriptor has nonlocal owner");
      }
      auto [_, inserted] = owned_index_by_particle_id.emplace(descriptor.particle_id, local_index);
      if (!inserted) {
        throw std::invalid_argument(
            "executeBlockingGhostRefreshExchangeFromDescriptors: duplicate owned particle_id in descriptor table");
      }
    } else {
      if (descriptor.owning_rank == world_rank) {
        throw std::invalid_argument(
            "executeBlockingGhostRefreshExchangeFromDescriptors: ghost descriptor is owned by local rank");
      }
      requested_particle_ids_by_rank[static_cast<std::size_t>(descriptor.owning_rank)].push_back(descriptor.particle_id);
      recv_indices_by_rank[static_cast<std::size_t>(descriptor.owning_rank)].push_back(local_index);
    }
  }

  const bool has_local_ghost_demands = std::any_of(
      recv_indices_by_rank.begin(), recv_indices_by_rank.end(), [](const auto& rows) { return !rows.empty(); });
  if (!mpi_context.isEnabled()) {
    if (has_local_ghost_demands) {
      throw std::runtime_error(
          "executeBlockingGhostRefreshExchangeFromDescriptors: non-empty ghost demand requires MPI");
    }
    return BlockingGhostRefreshExchange{
        .plan = buildExplicitGhostExchangePlan(
            world_rank,
            std::span<const int>{},
            std::span<const std::vector<std::uint32_t>>{},
            std::span<const std::vector<std::uint32_t>>{},
            ghostRefreshPayloadRecordBytes(),
            expected_epoch),
        .result = BlockingGhostExchangeResult{.received_ghosts = GhostExchangeBufferSoA{.epoch = expected_epoch}},
    };
  }

#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  constexpr int k_request_count_tag_base = 4810;
  constexpr int k_request_payload_tag_base = 5810;
  std::vector<std::vector<std::uint32_t>> send_indices_by_rank(static_cast<std::size_t>(world_size));
  for (int peer_rank = 0; peer_rank < world_size; ++peer_rank) {
    if (peer_rank == world_rank) {
      continue;
    }
    const auto& request_ids = requested_particle_ids_by_rank[static_cast<std::size_t>(peer_rank)];
    const std::uint64_t send_count = static_cast<std::uint64_t>(request_ids.size());
    std::uint64_t recv_count = 0;
    MPI_Sendrecv(
        &send_count,
        1,
        MPI_UINT64_T,
        peer_rank,
        ghostExchangeSequencedTag(k_request_count_tag_base, world_rank, peer_rank, expected_epoch.ghost_sync_epoch),
        &recv_count,
        1,
        MPI_UINT64_T,
        peer_rank,
        ghostExchangeSequencedTag(k_request_count_tag_base, world_rank, peer_rank, expected_epoch.ghost_sync_epoch),
        MPI_COMM_WORLD,
        MPI_STATUS_IGNORE);

    if (send_count > static_cast<std::uint64_t>(std::numeric_limits<int>::max()) ||
        recv_count > static_cast<std::uint64_t>(std::numeric_limits<int>::max())) {
      throw std::overflow_error(
          "executeBlockingGhostRefreshExchangeFromDescriptors: ghost request count exceeds MPI int count limit");
    }

    std::vector<std::uint64_t> received_request_ids(static_cast<std::size_t>(recv_count));
    MPI_Sendrecv(
        const_cast<std::uint64_t*>(request_ids.data()),
        static_cast<int>(request_ids.size()),
        MPI_UINT64_T,
        peer_rank,
        ghostExchangeSequencedTag(k_request_payload_tag_base, world_rank, peer_rank, expected_epoch.ghost_sync_epoch),
        received_request_ids.data(),
        static_cast<int>(received_request_ids.size()),
        MPI_UINT64_T,
        peer_rank,
        ghostExchangeSequencedTag(k_request_payload_tag_base, world_rank, peer_rank, expected_epoch.ghost_sync_epoch),
        MPI_COMM_WORLD,
        MPI_STATUS_IGNORE);

    auto& send_rows = send_indices_by_rank[static_cast<std::size_t>(peer_rank)];
    send_rows.reserve(received_request_ids.size());
    for (const std::uint64_t particle_id : received_request_ids) {
      const auto it = owned_index_by_particle_id.find(particle_id);
      if (it == owned_index_by_particle_id.end()) {
        throw std::runtime_error(
            "executeBlockingGhostRefreshExchangeFromDescriptors: peer requested a particle_id not owned by this rank");
      }
      send_rows.push_back(it->second);
    }
  }

  std::vector<int> neighbor_ranks;
  std::vector<std::vector<std::uint32_t>> send_indices_by_neighbor;
  std::vector<std::vector<std::uint32_t>> recv_indices_by_neighbor;
  for (int peer_rank = 0; peer_rank < world_size; ++peer_rank) {
    if (peer_rank == world_rank) {
      continue;
    }
    const auto& send_rows = send_indices_by_rank[static_cast<std::size_t>(peer_rank)];
    const auto& recv_rows = recv_indices_by_rank[static_cast<std::size_t>(peer_rank)];
    if (send_rows.empty() && recv_rows.empty()) {
      continue;
    }
    neighbor_ranks.push_back(peer_rank);
    send_indices_by_neighbor.push_back(send_rows);
    recv_indices_by_neighbor.push_back(recv_rows);
  }

  BlockingGhostRefreshExchange exchange;
  exchange.plan = buildExplicitGhostExchangePlan(
      world_rank,
      neighbor_ranks,
      send_indices_by_neighbor,
      recv_indices_by_neighbor,
      ghostRefreshPayloadRecordBytes(),
      expected_epoch);
  exchange.plan.exchange_sequence = expected_epoch.ghost_sync_epoch;
  exchange.result = executeBlockingGhostRefreshExchange(
      mpi_context,
      exchange.plan,
      local_ghost_descriptors,
      authoritative_local_state,
      expected_epoch);
  return exchange;
#else
  throw std::runtime_error(
      "executeBlockingGhostRefreshExchangeFromDescriptors: MPI support is not compiled in");
#endif
}

BlockingGhostExchangeResult executeBlockingGhostRefreshExchange(
    const MpiContext& mpi_context,
    const GhostExchangePlan& plan,
    std::span<const LocalGhostDescriptor> local_ghost_descriptors,
    const GhostExchangeBufferSoA& authoritative_local_state,
    const GhostLayerEpoch& expected_epoch) {
  validateBlockingGhostExchangeContracts(
      plan, local_ghost_descriptors, mpi_context.worldRank(), expected_epoch);
  if (!authoritative_local_state.isConsistent()) {
    throw std::invalid_argument("executeBlockingGhostRefreshExchange: authoritative local ghost payload state is inconsistent");
  }
  if (!authoritative_local_state.epoch.matches(expected_epoch)) {
    throw std::invalid_argument("executeBlockingGhostRefreshExchange: authoritative local payload epoch is stale");
  }

  BlockingGhostExchangeResult result;
  result.received_ghosts.epoch = expected_epoch;
  if (plan.neighbor_ranks.empty()) {
    return result;
  }
  if (!mpi_context.isEnabled()) {
    throw std::runtime_error(
        "executeBlockingGhostRefreshExchange: non-empty ghost exchange requires MPI; serial path must have no neighbors");
  }

#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  constexpr int k_size_tag_base = 6810;
  constexpr int k_payload_tag_base = 7810;
  for (std::size_t slot = 0; slot < plan.neighbor_ranks.size(); ++slot) {
    for (const std::uint32_t local_index : plan.send_local_indices_by_neighbor[slot]) {
      if (local_index >= authoritative_local_state.size() || local_index >= local_ghost_descriptors.size()) {
        throw std::out_of_range("executeBlockingGhostRefreshExchange: send descriptor index is outside local payload state");
      }
      if (authoritative_local_state.entity_id[local_index] != local_ghost_descriptors[local_index].particle_id) {
        throw std::invalid_argument("executeBlockingGhostRefreshExchange: send payload entity_id does not match owned descriptor particle_id");
      }
    }

    GhostExchangeBuffer send_buffer;
    send_buffer.packFrom(
        plan.outbound_transfers[slot],
        authoritative_local_state,
        plan.send_local_indices_by_neighbor[slot]);

    const std::uint64_t send_size = static_cast<std::uint64_t>(send_buffer.byteSize());
    std::uint64_t recv_size = 0;
    const int peer_rank = plan.neighbor_ranks[slot];
    MPI_Sendrecv(
        &send_size,
        1,
        MPI_UINT64_T,
        peer_rank,
        ghostExchangeSequencedTag(k_size_tag_base, mpi_context.worldRank(), peer_rank, plan.exchange_sequence),
        &recv_size,
        1,
        MPI_UINT64_T,
        peer_rank,
        ghostExchangeSequencedTag(k_size_tag_base, mpi_context.worldRank(), peer_rank, plan.exchange_sequence),
        MPI_COMM_WORLD,
        MPI_STATUS_IGNORE);

    if (send_size > static_cast<std::uint64_t>(std::numeric_limits<int>::max()) ||
        recv_size > static_cast<std::uint64_t>(std::numeric_limits<int>::max())) {
      throw std::overflow_error("executeBlockingGhostRefreshExchange: ghost payload exceeds MPI int count limit");
    }

    std::vector<std::uint8_t> recv_bytes(static_cast<std::size_t>(recv_size));
    const auto send_bytes = send_buffer.encodedBytes();
    MPI_Sendrecv(
        const_cast<std::uint8_t*>(send_bytes.data()),
        static_cast<int>(send_bytes.size()),
        MPI_BYTE,
        peer_rank,
        ghostExchangeSequencedTag(k_payload_tag_base, mpi_context.worldRank(), peer_rank, plan.exchange_sequence),
        recv_bytes.data(),
        static_cast<int>(recv_bytes.size()),
        MPI_BYTE,
        peer_rank,
        ghostExchangeSequencedTag(k_payload_tag_base, mpi_context.worldRank(), peer_rank, plan.exchange_sequence),
        MPI_COMM_WORLD,
        MPI_STATUS_IGNORE);

    GhostExchangeBuffer recv_buffer;
    recv_buffer.replaceEncodedBytes(std::move(recv_bytes));
    const std::size_t old_received_count = result.received_ghosts.size();
    recv_buffer.unpackAppendTo(plan.inbound_transfers[slot], result.received_ghosts);
    for (std::size_t i = 0; i < plan.recv_local_indices_by_neighbor[slot].size(); ++i) {
      const std::uint32_t local_index = plan.recv_local_indices_by_neighbor[slot][i];
      if (local_index >= local_ghost_descriptors.size()) {
        throw std::out_of_range("executeBlockingGhostRefreshExchange: receive descriptor index is outside residency table");
      }
      if (result.received_ghosts.entity_id[old_received_count + i] !=
          local_ghost_descriptors[local_index].particle_id) {
        throw std::invalid_argument("executeBlockingGhostRefreshExchange: received ghost entity_id does not match receive descriptor particle_id");
      }
    }
    result.sent_bytes += send_size;
    result.received_bytes += recv_size;
  }
  return result;
#else
  throw std::runtime_error("executeBlockingGhostRefreshExchange: MPI support is not compiled in");
#endif
}


bool RankDeviceAssignment::isValid() const noexcept {
  if (requested_device_count < 0 || visible_device_count < 0 || active_device_count < 0) {
    return false;
  }
  if (!uses_cuda) {
    return assigned_device_index == -1;
  }
  return active_device_count > 0 && assigned_device_index >= 0 && assigned_device_index < active_device_count &&
      visible_device_count >= active_device_count;
}

RankDeviceAssignment selectRankDeviceAssignment(
    int world_rank,
    int configured_gpu_devices,
    bool cuda_runtime_available,
    int visible_device_count) {
  if (world_rank < 0) {
    throw std::invalid_argument("world_rank must be non-negative");
  }
  if (configured_gpu_devices < 0) {
    throw std::invalid_argument("configured_gpu_devices must be >= 0");
  }
  if (visible_device_count < 0) {
    throw std::invalid_argument("visible_device_count must be >= 0");
  }

  RankDeviceAssignment assignment;
  assignment.requested_device_count = configured_gpu_devices;
  assignment.visible_device_count = visible_device_count;

  if (configured_gpu_devices == 0) {
    return assignment;
  }
  if (!cuda_runtime_available || visible_device_count == 0) {
    throw std::runtime_error(
        "parallel.gpu_devices requested CUDA PM execution, but no CUDA runtime devices are available");
  }
  if (configured_gpu_devices > visible_device_count) {
    throw std::runtime_error(
        "parallel.gpu_devices exceeds visible CUDA devices: requested=" + std::to_string(configured_gpu_devices) +
        ", visible=" + std::to_string(visible_device_count));
  }

  assignment.uses_cuda = true;
  assignment.active_device_count = configured_gpu_devices;
  assignment.assigned_device_index = world_rank % configured_gpu_devices;
  return assignment;
}

DistributedExecutionTopology buildDistributedExecutionTopology(
    std::size_t global_nx,
    std::size_t global_ny,
    std::size_t global_nz,
    const MpiContext& mpi_context,
    int mpi_ranks_expected,
    int configured_gpu_devices,
    bool cuda_runtime_available,
    int visible_device_count,
    std::string pm_decomposition_mode) {
  mpi_context.validateExpectedWorldSizeOrThrow(mpi_ranks_expected);

  DistributedExecutionTopology topology;
  topology.world_size = mpi_context.worldSize();
  topology.world_rank = mpi_context.worldRank();
  topology.mpi_enabled = mpi_context.isEnabled();
  topology.pm_decomposition_mode = std::move(pm_decomposition_mode);
  topology.pm_slab = makePmSlabLayout(global_nx, global_ny, global_nz, mpi_context.worldSize(), mpi_context.worldRank());
  topology.device_assignment =
      selectRankDeviceAssignment(mpi_context.worldRank(), configured_gpu_devices, cuda_runtime_available, visible_device_count);
  if (!topology.device_assignment.isValid()) {
    throw std::runtime_error("constructed distributed execution topology is invalid");
  }
  return topology;
}


}  // namespace cosmosim::parallel
