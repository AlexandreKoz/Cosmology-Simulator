#include "cosmosim/parallel/distributed_memory.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <iomanip>
#include <limits>
#include <numeric>
#include <optional>
#include <sstream>
#include <streambuf>
#include <unordered_map>
#include <unordered_set>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/core/memory_governor.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "parallel/internal/memory_constrained_sfc.hpp"
#include "parallel/internal/distributed_memory_support.hpp"

#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
#include <mpi.h>
#endif

namespace cosmosim::parallel {

// Core SFC decomposition and runtime-rebalance planning live here. MPI transport,
// TreePM, AMR, hydro, restart, PM, ghost-codec, diagnostics, and device-topology
// implementations are split into the sibling distributed_memory_*.cpp files.
namespace {

using internal::injectMpiTestFault;


[[nodiscard]] double clampUnit(double value) {
  if (value <= 0.0) {
    return 0.0;
  }
  if (value >= 1.0) {
    return std::nextafter(1.0, 0.0);
  }
  return value;
}

[[nodiscard]] std::uint64_t checkedUint64Add(
    std::uint64_t lhs,
    std::uint64_t rhs,
    std::string_view context) {
  if (rhs > std::numeric_limits<std::uint64_t>::max() - lhs) {
    throw std::overflow_error(std::string(context) + ": uint64 addition overflows");
  }
  return lhs + rhs;
}

[[nodiscard]] std::uint64_t checkedUint64Multiply(
    std::uint64_t lhs,
    std::uint64_t rhs,
    std::string_view context) {
  if (lhs != 0U && rhs > std::numeric_limits<std::uint64_t>::max() / lhs) {
    throw std::overflow_error(std::string(context) + ": uint64 multiplication overflows");
  }
  return lhs * rhs;
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

[[nodiscard]] DecompositionWorkComponents clampWorkComponents(
    DecompositionWorkComponents components) {
  components.particle_count_cost = std::max(0.0, components.particle_count_cost);
  components.gas_cell_cost = std::max(0.0, components.gas_cell_cost);
  components.tree_interaction_cost = std::max(0.0, components.tree_interaction_cost);
  components.pm_mesh_cost = std::max(0.0, components.pm_mesh_cost);
  components.amr_patch_cost = std::max(0.0, components.amr_patch_cost);
  components.active_fraction_cost = std::max(0.0, components.active_fraction_cost);
  components.memory_pressure_cost = std::max(0.0, components.memory_pressure_cost);
  components.transient_memory_cost = std::max(0.0, components.transient_memory_cost);
  components.source_event_cost = std::max(0.0, components.source_event_cost);
  components.communication_cost = std::max(0.0, components.communication_cost);
  components.gpu_occupancy_cost = std::max(0.0, components.gpu_occupancy_cost);
  components.generic_work_cost = std::max(0.0, components.generic_work_cost);
  return components;
}

[[nodiscard]] DecompositionWorkComponents effectiveWorkComponents(const DecompositionItem& item) {
  if (item.work_components.has_explicit_components) {
    return clampWorkComponents(item.work_components);
  }

  DecompositionWorkComponents components;
  components.particle_count_cost = (item.kind == DecompositionEntityKind::kParticle) ? 1.0 : 0.0;
  components.gas_cell_cost = (item.kind == DecompositionEntityKind::kHydroCell) ? 1.0 : 0.0;
  components.amr_patch_cost = (item.kind == DecompositionEntityKind::kAmrPatch) ? 1.0 : 0.0;
  components.pm_mesh_cost = (item.kind == DecompositionEntityKind::kPmMeshCell) ? 1.0 : 0.0;
  components.tree_interaction_cost = static_cast<double>(item.remote_tree_interactions_recent);
  components.active_fraction_cost = static_cast<double>(item.active_target_count_recent);
  components.memory_pressure_cost = static_cast<double>(item.memory_bytes);
  components.transient_memory_cost = 0.0;
  components.source_event_cost = 0.0;
  components.communication_cost = static_cast<double>(item.remote_tree_interactions_recent);
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
      weights.memory_pressure * (components.memory_pressure_cost + components.transient_memory_cost) +
      weights.gpu_occupancy * components.gpu_occupancy_cost +
      weights.generic_work * (components.generic_work_cost + components.source_event_cost + components.communication_cost);
}

// Legacy-field weighted load for one entity, then the component fallback
// chain shared by the item adapter and the explicit-components authority.
[[nodiscard]] double weightedLoadTail(
    const DecompositionWorkComponents& components,
    DecompositionEntityKind kind,
    double owned_particle_term,
    double active_target_term,
    double remote_tree_term,
    double work_term,
    double memory_term,
    const DecompositionConfig& config) {
  const double legacy = config.owned_particle_weight * owned_particle_term +
      config.active_target_weight * active_target_term +
      config.remote_tree_interaction_weight * remote_tree_term +
      config.work_weight * work_term + config.memory_weight * memory_term;
  if (legacy > 0.0) {
    return legacy;
  }
  const double component_fallback = componentWeightedLoad(components, config.component_weights);
  if (component_fallback > 0.0) {
    return component_fallback;
  }
  return (kind == DecompositionEntityKind::kParticle) ? 1.0 : std::max(1.0, components.rawTotal());
}

[[nodiscard]] double weightedLoad(const DecompositionItem& item, const DecompositionConfig& config) {
  const DecompositionWorkComponents components = effectiveWorkComponents(item);
  if (config.prefer_component_work_model &&
      item.work_components.has_explicit_components && hasNonZeroComponentWeight(config.component_weights)) {
    return std::max(0.0, componentWeightedLoad(components, config.component_weights));
  }
  return weightedLoadTail(
      components, item.kind,
      (item.kind == DecompositionEntityKind::kParticle) ? 1.0 : 0.0,
      static_cast<double>(item.active_target_count_recent),
      static_cast<double>(item.remote_tree_interactions_recent),
      std::max(0.0, item.work_units),
      static_cast<double>(item.memory_bytes),
      config);
}

[[nodiscard]] std::uint64_t sfcKeyForItem(const DecompositionItem& item, const DecompositionConfig& config) {
  return sfcKeyForPosition(item.x_comov, item.y_comov, item.z_comov, config);
}

using internal::SfcCutPoint;
using internal::lessSfcPoint;

[[nodiscard]] int ownerForSfcPoint(
    const SfcCutPoint point,
    std::span<const SfcCutPoint> cuts,
    int world_size) {
  const auto it = std::lower_bound(cuts.begin(), cuts.end(), point, lessSfcPoint);
  const auto rank = static_cast<int>(std::distance(cuts.begin(), it));
  return std::min(std::max(rank, 0), std::max(world_size - 1, 0));
}

void validateComponentWeights(const DecompositionWeightCoefficients& weights) {
  if (weights.particle_count < 0.0 || weights.gas_cell < 0.0 || weights.tree_interaction < 0.0 ||
      weights.pm_mesh < 0.0 || weights.amr_patch < 0.0 || weights.active_fraction < 0.0 ||
      weights.memory_pressure < 0.0 || weights.gpu_occupancy < 0.0 || weights.generic_work < 0.0) {
    throw std::invalid_argument("decomposition component weights must be non-negative");
  }
}

}  // namespace

std::uint64_t sfcKeyForPosition(
    double x_comov,
    double y_comov,
    double z_comov,
    const DecompositionConfig& config) {
  const std::uint32_t qx = quantize10bit(x_comov, config.domain_x_min_comov, config.domain_x_max_comov);
  const std::uint32_t qy = quantize10bit(y_comov, config.domain_y_min_comov, config.domain_y_max_comov);
  const std::uint32_t qz = quantize10bit(z_comov, config.domain_z_min_comov, config.domain_z_max_comov);
  return mortonKey3d(qx, qy, qz);
}

double weightedLoadFromExplicitComponents(
    const DecompositionWorkComponents& components,
    DecompositionEntityKind kind,
    std::uint64_t active_target_count_recent,
    std::uint64_t remote_tree_interactions_recent,
    double work_units,
    std::uint64_t memory_bytes,
    const DecompositionConfig& config) {
  const DecompositionWorkComponents clamped = clampWorkComponents(components);
  if (config.prefer_component_work_model && hasNonZeroComponentWeight(config.component_weights)) {
    return std::max(0.0, componentWeightedLoad(clamped, config.component_weights));
  }
  return weightedLoadTail(
      clamped, kind,
      (kind == DecompositionEntityKind::kParticle) ? 1.0 : 0.0,
      static_cast<double>(active_target_count_recent),
      static_cast<double>(remote_tree_interactions_recent),
      std::max(0.0, work_units),
      static_cast<double>(memory_bytes),
      config);
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
  metrics.transient_memory_cost_by_rank[rank] += sign * components.transient_memory_cost;
  metrics.source_event_cost_by_rank[rank] += sign * components.source_event_cost;
  metrics.communication_cost_by_rank[rank] += sign * components.communication_cost;
  metrics.gpu_occupancy_cost_by_rank[rank] += sign * components.gpu_occupancy_cost;
  metrics.generic_work_cost_by_rank[rank] += sign * components.generic_work_cost;
}

bool runtimeDecompositionHasGasIncidenceSource(const core::SimulationState& state) noexcept {
  return !state.gas_cell_identity.records().empty() && state.cells.size() != 0U;
}

RuntimeDecompositionSourceMemoryEstimate estimateRuntimeDecompositionSourceStorage(
    const core::SimulationState& state,
    std::span<const std::uint32_t> active_particle_indices) {
  RuntimeDecompositionSourceMemoryEstimate estimate;
  if (!active_particle_indices.empty()) {
    estimate.active_mask_bytes = core::checkedIntegralNarrow<std::uint64_t>(
        core::checkedSizeMultiply(
            state.particles.size(), sizeof(std::uint8_t),
            "runtime decomposition active mask estimate"),
        "runtime decomposition active mask estimate byte width");
  }
  estimate.compact_patch_index_bytes = core::checkedIntegralNarrow<std::uint64_t>(
      core::checkedSizeMultiply(
          state.patches.size(), sizeof(std::uint32_t),
          "runtime decomposition compact patch index estimate"),
      "runtime decomposition compact patch index estimate byte width");
  if (runtimeDecompositionHasGasIncidenceSource(state)) {
    estimate.gas_incidence_offset_bytes = core::checkedIntegralNarrow<std::uint64_t>(
        core::checkedSizeMultiply(
            core::checkedSizeAdd(
                state.particles.size(), std::size_t{1U},
                "runtime decomposition gas offset estimate"),
            sizeof(std::uint32_t),
            "runtime decomposition gas offset estimate"),
        "runtime decomposition gas offset estimate byte width");
    estimate.gas_incidence_index_bytes = core::checkedIntegralNarrow<std::uint64_t>(
        core::checkedSizeMultiply(
            state.cells.size(), sizeof(std::uint32_t),
            "runtime decomposition gas index estimate"),
        "runtime decomposition gas index estimate byte width");
    estimate.gas_incidence_construction_bytes = core::checkedIntegralNarrow<std::uint64_t>(
        core::checkedSizeMultiply(
            state.cells.size(), sizeof(std::uint64_t) * 2U,
            "runtime decomposition gas construction estimate"),
        "runtime decomposition gas construction estimate byte width");
  }
  estimate.other_source_owned_bytes = core::checkedIntegralNarrow<std::uint64_t>(
      core::checkedSizeMultiply(
          std::size_t{5U}, sizeof(std::uint64_t),
          "runtime decomposition species table estimate"),
      "runtime decomposition species table estimate byte width");
  std::uint64_t total_bytes = 0U;
  total_bytes = core::checkedMemoryBytesAdd(
      total_bytes, estimate.active_mask_bytes,
      "runtime decomposition source estimate total");
  total_bytes = core::checkedMemoryBytesAdd(
      total_bytes, estimate.compact_patch_index_bytes,
      "runtime decomposition source estimate total");
  total_bytes = core::checkedMemoryBytesAdd(
      total_bytes, estimate.gas_incidence_offset_bytes,
      "runtime decomposition source estimate total");
  total_bytes = core::checkedMemoryBytesAdd(
      total_bytes, estimate.gas_incidence_index_bytes,
      "runtime decomposition source estimate total");
  total_bytes = core::checkedMemoryBytesAdd(
      total_bytes, estimate.gas_incidence_construction_bytes,
      "runtime decomposition source estimate total");
  total_bytes = core::checkedMemoryBytesAdd(
      total_bytes, estimate.other_source_owned_bytes,
      "runtime decomposition source estimate total");
  estimate.total_bytes = total_bytes;
  return estimate;
}

CompactStartupPlannerMemoryEstimate estimateCompactStartupPlannerTransientBytes(
    std::size_t entity_upper_bound,
    std::size_t patch_count,
    std::size_t cell_count,
    std::size_t world_size,
    std::size_t pm_x_bins,
    std::size_t density_grid_cell_count) {
  CompactStartupPlannerMemoryEstimate estimate;
  // Compact records: startup records.reserve(particles + patches) capacity.
  estimate.compact_record_bytes = core::checkedIntegralNarrow<std::uint64_t>(
      core::checkedSizeMultiply(
          entity_upper_bound, sizeof(CompactRuntimeDecompositionRecord),
          "compact startup planner record bytes"),
      "compact startup planner record byte width");
  // MemoryGroup: the cut algorithm groups.reserve(entity_upper_bound) and the
  // worst case is one group per record (all distinct (key, id) points). The
  // mirrored field layout below is enforced against the actual local struct
  // by a static_assert inside buildMortonSfcDecompositionFromCompact.
  estimate.memory_group_bytes = core::checkedIntegralNarrow<std::uint64_t>(
      core::checkedSizeMultiply(
          entity_upper_bound,
          sizeof(std::size_t) * 2U + sizeof(std::uint64_t) + sizeof(double),
          "compact startup planner memory group bytes"),
      "compact startup planner memory group byte width");
  // Plan owner vector (sentinel-initialized, int per record) and sorted
  // index vector (size_t per record), sized to the record population.
  estimate.owning_rank_bytes = core::checkedIntegralNarrow<std::uint64_t>(
      core::checkedSizeMultiply(
          entity_upper_bound, sizeof(int),
          "compact startup planner owner vector bytes"),
      "compact startup planner owner vector byte width");
  estimate.sorted_index_bytes = core::checkedIntegralNarrow<std::uint64_t>(
      core::checkedSizeMultiply(
          entity_upper_bound, sizeof(std::size_t),
          "compact startup planner sorted index bytes"),
      "compact startup planner sorted index byte width");
  // Startup occupancy grids: density occupancy, active occupancy, and gas
  // occupancy (three uint32 lanes of density_grid_cell_count each) plus the
  // PM-x occupancy lane. Grid policy stays owned by the startup workflow;
  // the caller supplies its current dimensions.
  estimate.occupancy_bytes = core::checkedIntegralNarrow<std::uint64_t>(
      core::checkedSizeAdd(
          core::checkedSizeMultiply(
              core::checkedSizeAdd(
                  core::checkedSizeMultiply(
                      density_grid_cell_count, 3U,
                      "compact startup planner density grid lane count"),
                  pm_x_bins, "compact startup planner occupancy lane count"),
              sizeof(std::uint32_t),
              "compact startup planner occupancy bytes"),
          0U, "compact startup planner occupancy bytes"),
      "compact startup planner occupancy byte width");
  // Patch mapping scratch: patch_cell_count and the included-patch-row
  // write-back map, each reserved to patch_count (upper bound).
  estimate.patch_mapping_bytes = core::checkedIntegralNarrow<std::uint64_t>(
      core::checkedSizeMultiply(
          patch_count, sizeof(std::uint32_t) * 2U,
          "compact startup planner patch mapping bytes"),
      "compact startup planner patch mapping byte width");
  // Rank-scaled plan lanes: LoadBalanceMetrics holds one double weighted-load
  // lane, five uint64 lanes (memory, owned, active, remote, peak), and
  // thirteen double component lanes per rank, plus ranges_by_rank (two
  // size_t per rank), plus the mandatory_begin cut array ((world_size + 1)
  // size_t). Also the bounded per-entity temporaries that coexist with the
  // planner population: the worst-case rowsForParentParticleId lookup (all
  // cell rows) and the worst-case seen-patch list (all patch rows).
  const std::size_t rank_metric_lane = core::checkedSizeAdd(
      core::checkedSizeAdd(
          sizeof(double) * 14U,
          sizeof(std::uint64_t) * 5U,
          "compact startup planner rank metric lane"),
      sizeof(std::size_t) * 2U,
      "compact startup planner rank metric lane");
  std::size_t other_scratch = core::checkedSizeMultiply(
      world_size, rank_metric_lane,
      "compact startup planner rank-scaled scratch");
  other_scratch = core::checkedSizeAdd(
      other_scratch,
      core::checkedSizeMultiply(
          core::checkedSizeAdd(world_size, 1U,
              "compact startup planner mandatory lane count"),
          sizeof(std::size_t),
          "compact startup planner mandatory lane bytes"),
      "compact startup planner rank-scaled scratch");
  other_scratch = core::checkedSizeAdd(
      other_scratch,
      core::checkedSizeMultiply(
          cell_count, sizeof(std::uint32_t),
          "compact startup planner gas row lookup bytes"),
      "compact startup planner rank-scaled scratch");
  other_scratch = core::checkedSizeAdd(
      other_scratch,
      core::checkedSizeMultiply(
          patch_count, sizeof(std::uint32_t),
          "compact startup planner seen patch bytes"),
      "compact startup planner rank-scaled scratch");
  estimate.other_known_scratch_bytes =
      core::checkedIntegralNarrow<std::uint64_t>(
          other_scratch, "compact startup planner other scratch byte width");

  std::uint64_t total_bytes = 0U;
  total_bytes = core::checkedMemoryBytesAdd(
      total_bytes, estimate.compact_record_bytes,
      "compact startup planner transient total");
  total_bytes = core::checkedMemoryBytesAdd(
      total_bytes, estimate.memory_group_bytes,
      "compact startup planner transient total");
  total_bytes = core::checkedMemoryBytesAdd(
      total_bytes, estimate.owning_rank_bytes,
      "compact startup planner transient total");
  total_bytes = core::checkedMemoryBytesAdd(
      total_bytes, estimate.sorted_index_bytes,
      "compact startup planner transient total");
  total_bytes = core::checkedMemoryBytesAdd(
      total_bytes, estimate.occupancy_bytes,
      "compact startup planner transient total");
  total_bytes = core::checkedMemoryBytesAdd(
      total_bytes, estimate.patch_mapping_bytes,
      "compact startup planner transient total");
  total_bytes = core::checkedMemoryBytesAdd(
      total_bytes, estimate.other_known_scratch_bytes,
      "compact startup planner transient total");
  estimate.total_bytes = total_bytes;
  return estimate;
}

bool GhostLayerEpoch::matches(const GhostLayerEpoch& expected) const noexcept {
  return decomposition_epoch == expected.decomposition_epoch && ghost_sync_epoch == expected.ghost_sync_epoch &&
      particle_index_generation == expected.particle_index_generation;
}

double DecompositionWorkComponents::rawTotal() const noexcept {
  return particle_count_cost + gas_cell_cost + tree_interaction_cost + pm_mesh_cost + amr_patch_cost +
      active_fraction_cost + memory_pressure_cost + transient_memory_cost + source_event_cost +
      communication_cost + gpu_occupancy_cost + generic_work_cost;
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

std::uint64_t sfcKeyForCompactRuntimeRecord(
    const CompactRuntimeDecompositionRecord& record) noexcept {
  return record.sfc_key;
}

DecompositionPlan buildMortonSfcDecompositionFromCompact(
    std::span<CompactRuntimeDecompositionRecord> records,
    const DecompositionConfig& config,
    std::span<const DecompositionWorkComponents> components_by_local_index) {
  if (config.world_size <= 0) {
    throw std::invalid_argument("world_size must be positive");
  }
  if (config.owned_particle_weight < 0.0 || config.active_target_weight < 0.0 ||
      config.remote_tree_interaction_weight < 0.0 || config.work_weight < 0.0 || config.memory_weight < 0.0) {
    throw std::invalid_argument("decomposition weights must be non-negative");
  }
  validateComponentWeights(config.component_weights);
  if (!components_by_local_index.empty() && components_by_local_index.size() < records.size()) {
    throw std::invalid_argument("decomposition component span does not cover compact records");
  }
  for (const CompactRuntimeDecompositionRecord& record : records) {
    if (!std::isfinite(record.weighted_load)) {
      throw std::invalid_argument("decomposition item weighted load must be finite");
    }
    if (record.local_index >= records.size()) {
      throw std::invalid_argument(
          "compact decomposition record local_index must be dense within [0, record count)");
    }
  }

  std::sort(records.begin(), records.end(),
            [](const CompactRuntimeDecompositionRecord& lhs, const CompactRuntimeDecompositionRecord& rhs) {
              if (lhs.sfc_key != rhs.sfc_key) {
                return lhs.sfc_key < rhs.sfc_key;
              }
              if (lhs.entity_id != rhs.entity_id) {
                return lhs.entity_id < rhs.entity_id;
              }
              return lhs.local_index < rhs.local_index;
            });

  DecompositionPlan plan;
  // Sentinel owner: -1 is never a valid MPI rank, so a missing local_index
  // can never masquerade as rank 0 ownership. Duplicate writes are rejected
  // below and full coverage is verified after the fill pass.
  constexpr int k_unassigned_owner_rank = -1;
  plan.owning_rank_by_item.assign(records.size(), k_unassigned_owner_rank);
  plan.sorted_indices.resize(records.size());
  plan.ranges_by_rank.assign(static_cast<std::size_t>(config.world_size), RankRange{});
  plan.metrics.weighted_load_by_rank.assign(static_cast<std::size_t>(config.world_size), 0.0);
  plan.metrics.memory_bytes_by_rank.assign(static_cast<std::size_t>(config.world_size), 0ULL);
  plan.metrics.peak_memory_bytes_by_rank.assign(static_cast<std::size_t>(config.world_size), 0ULL);
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
  plan.metrics.transient_memory_cost_by_rank.assign(static_cast<std::size_t>(config.world_size), 0.0);
  plan.metrics.source_event_cost_by_rank.assign(static_cast<std::size_t>(config.world_size), 0.0);
  plan.metrics.communication_cost_by_rank.assign(static_cast<std::size_t>(config.world_size), 0.0);
  plan.metrics.gpu_occupancy_cost_by_rank.assign(static_cast<std::size_t>(config.world_size), 0.0);
  plan.metrics.generic_work_cost_by_rank.assign(static_cast<std::size_t>(config.world_size), 0.0);

  for (std::size_t sorted_pos = 0; sorted_pos < records.size(); ++sorted_pos) {
    plan.sorted_indices[sorted_pos] = records[sorted_pos].local_index;
  }

  const double total_load = std::accumulate(
      records.begin(), records.end(), 0.0,
      [](double acc, const CompactRuntimeDecompositionRecord& entry) { return acc + entry.weighted_load; });

  if (!records.empty()) {
    const bool enforce_rank_memory_limit = config.max_rank_memory_bytes != 0U;
    if (enforce_rank_memory_limit &&
        config.rank_transient_reserve_bytes >= config.max_rank_memory_bytes) {
      throw std::invalid_argument(
          "decomposition transient reserve must be smaller than the hard rank memory ceiling");
    }
    const std::uint64_t persistent_rank_limit = enforce_rank_memory_limit
        ? config.max_rank_memory_bytes - config.rank_transient_reserve_bytes
        : std::numeric_limits<std::uint64_t>::max();

    // An indivisible SFC point may represent several entities. Keep equal
    // (key, ID) points together, matching the distributed ownership contract.
    // Work is a secondary objective: a right-to-left greedy packing computes
    // the earliest prefix that may be left to each rank without making the
    // suffix infeasible. Greedy packing is exact for nonnegative integer bytes
    // on an ordered sequence; no work-optimal cut may violate this boundary.
    struct MemoryGroup {
      std::size_t begin = 0U;
      std::size_t end = 0U;
      std::uint64_t memory_bytes = 0U;
      double weighted_load = 0.0;
    };
    // estimateCompactStartupPlannerTransientBytes mirrors this layout; keep
    // the mirror honest so planner admission cannot drift from reality.
    static_assert(sizeof(MemoryGroup) ==
            sizeof(std::size_t) * 2U + sizeof(std::uint64_t) + sizeof(double),
        "MemoryGroup layout must match estimateCompactStartupPlannerTransientBytes");
    std::vector<MemoryGroup> groups;
    groups.reserve(records.size());
    for (std::size_t pos = 0U; pos < records.size(); ++pos) {
      const auto& record = records[pos];
      const std::uint64_t bytes = record.memory_bytes;
      if (groups.empty() ||
          records[groups.back().begin].sfc_key != record.sfc_key ||
          records[groups.back().begin].entity_id != record.entity_id) {
        groups.push_back(MemoryGroup{.begin = pos, .end = pos});
      }
      auto& group = groups.back();
      group.end = pos + 1U;
      group.memory_bytes = checkedUint64Add(group.memory_bytes, bytes, "SFC grouped memory");
      group.weighted_load += record.weighted_load;
      if (enforce_rank_memory_limit && group.memory_bytes > persistent_rank_limit) {
        throw std::runtime_error(
            "SFC decomposition cannot satisfy hard rank memory ceiling: one indivisible SFC group exceeds the persistent allowance");
      }
    }
    const std::size_t active_rank_count =
        std::min<std::size_t>(static_cast<std::size_t>(config.world_size), groups.size());
    const double target_per_active_rank = active_rank_count > 0U
        ? total_load / static_cast<double>(active_rank_count) : 0.0;
    std::vector<std::size_t> mandatory_begin(active_rank_count + 1U, 0U);
    mandatory_begin[0] = groups.size();
    if (enforce_rank_memory_limit) {
      std::size_t suffix_begin = groups.size();
      std::uint64_t suffix_bytes = 0U;
      for (std::size_t ranks = 1U; ranks < active_rank_count; ++ranks) {
        while (suffix_begin > 0U &&
               groups[suffix_begin - 1U].memory_bytes <= persistent_rank_limit - suffix_bytes) {
          --suffix_begin;
          suffix_bytes += groups[suffix_begin].memory_bytes;
        }
        mandatory_begin[ranks] = suffix_begin;
        suffix_bytes = 0U;
      }
    }
    std::size_t current_rank = 0U;
    std::size_t rank_begin = 0U;
    double cumulative_load = 0.0;
    std::uint64_t current_rank_memory_bytes = 0U;
    for (std::size_t group_pos = 0U; group_pos < groups.size(); ++group_pos) {
      const auto& group = groups[group_pos];
      const bool rank_has_items = group_pos > rank_begin;
      const bool would_exceed_memory = enforce_rank_memory_limit && rank_has_items &&
          group.memory_bytes > persistent_rank_limit - current_rank_memory_bytes;
      if (would_exceed_memory) {
        if (current_rank + 1U >= active_rank_count ||
            group_pos < mandatory_begin[active_rank_count - current_rank - 1U]) {
          throw std::runtime_error(
              "SFC decomposition cannot satisfy hard rank memory ceiling with the available ranks");
        }
        plan.ranges_by_rank[current_rank] = RankRange{
            .begin_sorted = groups[rank_begin].begin, .end_sorted = group.begin};
        ++current_rank;
        rank_begin = group_pos;
        current_rank_memory_bytes = 0U;
      }
      current_rank_memory_bytes = checkedUint64Add(
          current_rank_memory_bytes, group.memory_bytes, "SFC rank persistent memory accumulation");
      cumulative_load += group.weighted_load;
      if (current_rank + 1U >= active_rank_count) {
        continue;
      }
      const std::size_t groups_remaining = groups.size() - group_pos - 1U;
      const std::size_t ranks_remaining = active_rank_count - current_rank - 1U;
      if (groups_remaining < ranks_remaining) {
        continue;
      }
      const bool must_cut = groups_remaining == ranks_remaining;
      const bool crossed_target = cumulative_load >=
          target_per_active_rank * static_cast<double>(current_rank + 1U);
      const bool suffix_feasible = !enforce_rank_memory_limit ||
          group_pos + 1U >= mandatory_begin[ranks_remaining];
      if (!suffix_feasible || (!must_cut && !crossed_target)) {
        continue;
      }
      plan.ranges_by_rank[current_rank] = RankRange{
          .begin_sorted = groups[rank_begin].begin, .end_sorted = group.end};
      ++current_rank;
      rank_begin = group_pos + 1U;
      current_rank_memory_bytes = 0U;
    }
    if (!groups.empty()) {
      plan.ranges_by_rank[current_rank] = RankRange{
          .begin_sorted = groups[rank_begin].begin, .end_sorted = records.size()};
    }
    for (std::size_t rank = current_rank + 1U;
         rank < static_cast<std::size_t>(config.world_size); ++rank) {
      plan.ranges_by_rank[rank] = RankRange{.begin_sorted = records.size(), .end_sorted = records.size()};
    }

    for (std::size_t rank = 0; rank < plan.ranges_by_rank.size(); ++rank) {
      const RankRange range = plan.ranges_by_rank[rank];
      for (std::size_t sorted_pos = range.begin_sorted; sorted_pos < range.end_sorted; ++sorted_pos) {
        const CompactRuntimeDecompositionRecord& record = records[sorted_pos];
        int& owner_slot = plan.owning_rank_by_item[record.local_index];
        if (owner_slot != k_unassigned_owner_rank) {
          throw std::invalid_argument(
              "compact decomposition local_index must be unique within the record population");
        }
        owner_slot = static_cast<int>(rank);
        plan.metrics.weighted_load_by_rank[rank] += record.weighted_load;
        plan.metrics.memory_bytes_by_rank[rank] += record.memory_bytes;
        if (record.kind == DecompositionEntityKind::kParticle) {
          ++plan.metrics.owned_particles_by_rank[rank];
        }
        plan.metrics.active_targets_by_rank[rank] += record.active_target_count_recent;
        plan.metrics.remote_tree_interactions_by_rank[rank] += record.remote_tree_interactions_recent;
        if (!components_by_local_index.empty()) {
          addWorkComponentsToMetrics(
              plan.metrics, rank, components_by_local_index[record.local_index], 1.0);
        }
      }
    }
  }

  if (std::any_of(
          plan.owning_rank_by_item.begin(), plan.owning_rank_by_item.end(),
          [](int owner) { return owner == k_unassigned_owner_rank; })) {
    throw std::invalid_argument(
        "compact decomposition local_index must form a permutation of [0, record count)");
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
  std::uint64_t total_peak_memory = 0U;
  for (std::size_t rank = 0; rank < plan.metrics.memory_bytes_by_rank.size(); ++rank) {
    const bool rank_has_work = plan.ranges_by_rank[rank].begin_sorted != plan.ranges_by_rank[rank].end_sorted;
    plan.metrics.peak_memory_bytes_by_rank[rank] = rank_has_work
        ? checkedUint64Add(plan.metrics.memory_bytes_by_rank[rank], config.rank_transient_reserve_bytes,
                           "SFC rank peak memory metric")
        : 0U;
    total_peak_memory = checkedUint64Add(total_peak_memory, plan.metrics.peak_memory_bytes_by_rank[rank],
                                         "SFC total peak memory metric");
  }
  const auto max_peak_it = std::max_element(
      plan.metrics.peak_memory_bytes_by_rank.begin(), plan.metrics.peak_memory_bytes_by_rank.end());
  plan.metrics.max_peak_memory_bytes =
      (max_peak_it == plan.metrics.peak_memory_bytes_by_rank.end()) ? 0U : *max_peak_it;
  const double mean_peak_memory = plan.metrics.peak_memory_bytes_by_rank.empty()
      ? 0.0
      : static_cast<double>(total_peak_memory) / static_cast<double>(plan.metrics.peak_memory_bytes_by_rank.size());
  plan.metrics.peak_memory_imbalance_ratio = mean_peak_memory > 0.0
      ? static_cast<double>(plan.metrics.max_peak_memory_bytes) / mean_peak_memory
      : 0.0;
  if (config.max_rank_memory_bytes != 0U && plan.metrics.max_peak_memory_bytes > config.max_rank_memory_bytes) {
    throw std::logic_error("SFC decomposition produced a rank above its hard memory ceiling");
  }

  return plan;
}

DecompositionPlan buildMortonSfcDecomposition(
    std::span<const DecompositionItem> items, const DecompositionConfig& config) {
  if (config.world_size <= 0) {
    throw std::invalid_argument("world_size must be positive");
  }
  if (config.owned_particle_weight < 0.0 || config.active_target_weight < 0.0 ||
      config.remote_tree_interaction_weight < 0.0 || config.work_weight < 0.0 || config.memory_weight < 0.0) {
    throw std::invalid_argument("decomposition weights must be non-negative");
  }
  validateComponentWeights(config.component_weights);

  std::vector<CompactRuntimeDecompositionRecord> records;
  std::vector<DecompositionWorkComponents> components;
  records.reserve(items.size());
  components.reserve(items.size());
  for (std::size_t i = 0; i < items.size(); ++i) {
    const DecompositionItem& item = items[i];
    const double load = weightedLoad(item, config);
    if (!std::isfinite(load)) {
      throw std::invalid_argument("decomposition item weighted load must be finite");
    }
    records.push_back(CompactRuntimeDecompositionRecord{
        .entity_id = item.entity_id,
        .sfc_key = sfcKeyForItem(item, config),
        .memory_bytes = item.memory_bytes,
        .weighted_load = load,
        .local_index = i,
        .active_target_count_recent = item.active_target_count_recent,
        .remote_tree_interactions_recent = item.remote_tree_interactions_recent,
        .current_owner_rank = item.current_owner_rank,
        .kind = item.kind,
    });
    components.push_back(effectiveWorkComponents(item));
  }
  return buildMortonSfcDecompositionFromCompact(records, config, components);
}

std::vector<TopDomainLeaf> buildAuthoritativeTopDomainLeaves(
    std::span<const DecompositionItem> local_items,
    const DecompositionConfig& config,
    int owner_rank,
    std::uint64_t decomposition_epoch,
    std::size_t max_leaves_per_rank) {
  if (config.world_size <= 0 || owner_rank < 0 || owner_rank >= config.world_size) {
    throw std::invalid_argument("top-domain leaf builder received invalid rank metadata");
  }
  if (max_leaves_per_rank == 0U) {
    throw std::invalid_argument("top-domain leaf builder requires at least one leaf slot per rank");
  }

  struct KeyedLocalItem {
    const DecompositionItem* item = nullptr;
    std::uint64_t key = 0U;
  };
  std::vector<KeyedLocalItem> keyed;
  keyed.reserve(local_items.size());
  for (const DecompositionItem& item : local_items) {
    if (item.current_owner_rank != owner_rank) {
      continue;
    }
    const std::array values{
        item.x_comov, item.y_comov, item.z_comov,
        item.has_spatial_bounds ? item.min_x_comov : item.x_comov,
        item.has_spatial_bounds ? item.max_x_comov : item.x_comov,
        item.has_spatial_bounds ? item.min_y_comov : item.y_comov,
        item.has_spatial_bounds ? item.max_y_comov : item.y_comov,
        item.has_spatial_bounds ? item.min_z_comov : item.z_comov,
        item.has_spatial_bounds ? item.max_z_comov : item.z_comov};
    if (std::any_of(values.begin(), values.end(), [](double value) { return !std::isfinite(value); })) {
      throw std::invalid_argument("top-domain leaf builder found non-finite decomposition geometry");
    }
    if (item.has_spatial_bounds &&
        (item.min_x_comov > item.max_x_comov || item.min_y_comov > item.max_y_comov ||
         item.min_z_comov > item.max_z_comov)) {
      throw std::invalid_argument("top-domain leaf builder found inverted decomposition bounds");
    }
    keyed.push_back(KeyedLocalItem{.item = &item, .key = sfcKeyForItem(item, config)});
  }

  std::stable_sort(keyed.begin(), keyed.end(), [](const KeyedLocalItem& lhs, const KeyedLocalItem& rhs) {
    if (lhs.key != rhs.key) {
      return lhs.key < rhs.key;
    }
    return lhs.item->entity_id < rhs.item->entity_id;
  });
  if (keyed.empty()) {
    return {};
  }

  const std::size_t leaf_count = std::min(max_leaves_per_rank, keyed.size());
  std::vector<TopDomainLeaf> leaves;
  leaves.reserve(leaf_count);
  for (std::size_t leaf_ordinal = 0U; leaf_ordinal < leaf_count; ++leaf_ordinal) {
    const std::size_t begin = (leaf_ordinal * keyed.size()) / leaf_count;
    const std::size_t end = ((leaf_ordinal + 1U) * keyed.size()) / leaf_count;
    if (begin == end) {
      continue;
    }
    const DecompositionItem& first = *keyed[begin].item;
    const auto item_bounds = [](const DecompositionItem& item) {
      return std::array<double, 6>{
          item.has_spatial_bounds ? item.min_x_comov : item.x_comov,
          item.has_spatial_bounds ? item.max_x_comov : item.x_comov,
          item.has_spatial_bounds ? item.min_y_comov : item.y_comov,
          item.has_spatial_bounds ? item.max_y_comov : item.y_comov,
          item.has_spatial_bounds ? item.min_z_comov : item.z_comov,
          item.has_spatial_bounds ? item.max_z_comov : item.z_comov};
    };
    const auto first_bounds = item_bounds(first);
    TopDomainLeaf leaf{
        .owner_rank = owner_rank,
        .decomposition_epoch = decomposition_epoch,
        .sfc_key_begin = keyed[begin].key,
        .sfc_key_end = keyed[end - 1U].key,
        .min_x_comov = first_bounds[0],
        .max_x_comov = first_bounds[1],
        .min_y_comov = first_bounds[2],
        .max_y_comov = first_bounds[3],
        .min_z_comov = first_bounds[4],
        .max_z_comov = first_bounds[5],
        .periodic_geometry = true,
    };
    std::uint64_t id_hash = 1469598103934665603ULL;
    const auto mix_id = [&id_hash](std::uint64_t value) {
      id_hash ^= value;
      id_hash *= 1099511628211ULL;
    };
    mix_id(static_cast<std::uint64_t>(static_cast<std::uint32_t>(owner_rank)));
    // Multiple leaf groups may legitimately share identical begin/end Morton
    // keys when tied keys straddle a deterministic group boundary. Include the
    // stable within-rank group ordinal so those leaves remain distinct pseudo
    // nodes while preserving identity across geometry-only refits.
    mix_id(static_cast<std::uint64_t>(leaf_ordinal));
    mix_id(leaf.sfc_key_begin);
    mix_id(leaf.sfc_key_end);
    for (std::size_t slot = begin; slot < end; ++slot) {
      const DecompositionItem& item = *keyed[slot].item;
      const auto bounds = item_bounds(item);
      leaf.min_x_comov = std::min(leaf.min_x_comov, bounds[0]);
      leaf.max_x_comov = std::max(leaf.max_x_comov, bounds[1]);
      leaf.min_y_comov = std::min(leaf.min_y_comov, bounds[2]);
      leaf.max_y_comov = std::max(leaf.max_y_comov, bounds[3]);
      leaf.min_z_comov = std::min(leaf.min_z_comov, bounds[4]);
      leaf.max_z_comov = std::max(leaf.max_z_comov, bounds[5]);
      leaf.work_weight += weightedLoad(item, config);
      ++leaf.entity_count;
    }
    leaf.domain_leaf_id = id_hash;
    leaves.push_back(leaf);
  }
  return leaves;
}

std::vector<TopDomainLeaf> buildAuthoritativeTopDomainLeavesFromCompact(
    std::span<const DecompositionItem> local_items,
    const DecompositionConfig& config,
    int owner_rank,
    std::uint64_t decomposition_epoch,
    std::size_t max_leaves_per_rank) {
  if (config.world_size <= 0 || owner_rank < 0 || owner_rank >= config.world_size) {
    throw std::invalid_argument("compact top-domain leaf builder received invalid rank metadata");
  }
  if (max_leaves_per_rank == 0U) {
    throw std::invalid_argument("compact top-domain leaf builder requires at least one leaf slot per rank");
  }

  std::vector<CompactRuntimeDecompositionRecord> records;
  records.reserve(local_items.size());
  for (std::size_t i = 0; i < local_items.size(); ++i) {
    const DecompositionItem& item = local_items[i];
    if (item.current_owner_rank == owner_rank) {
      const std::array values{
          item.x_comov, item.y_comov, item.z_comov,
          item.has_spatial_bounds ? item.min_x_comov : item.x_comov,
          item.has_spatial_bounds ? item.max_x_comov : item.x_comov,
          item.has_spatial_bounds ? item.min_y_comov : item.y_comov,
          item.has_spatial_bounds ? item.max_y_comov : item.y_comov,
          item.has_spatial_bounds ? item.min_z_comov : item.z_comov,
          item.has_spatial_bounds ? item.max_z_comov : item.z_comov};
      if (std::any_of(values.begin(), values.end(), [](double value) { return !std::isfinite(value); })) {
        throw std::invalid_argument("compact top-domain leaf builder found non-finite decomposition geometry");
      }
      if (item.has_spatial_bounds &&
          (item.min_x_comov > item.max_x_comov || item.min_y_comov > item.max_y_comov ||
           item.min_z_comov > item.max_z_comov)) {
        throw std::invalid_argument("compact top-domain leaf builder found inverted decomposition bounds");
      }
      records.push_back(CompactRuntimeDecompositionRecord{
          .entity_id = item.entity_id,
          .sfc_key = sfcKeyForItem(item, config),
          .memory_bytes = item.memory_bytes,
          .weighted_load = weightedLoad(item, config),
          .local_index = i,
          .active_target_count_recent = item.active_target_count_recent,
          .remote_tree_interactions_recent = item.remote_tree_interactions_recent,
          .current_owner_rank = item.current_owner_rank,
          .kind = item.kind,
      });
    }
  }

  std::sort(records.begin(), records.end(),
            [](const CompactRuntimeDecompositionRecord& lhs,
               const CompactRuntimeDecompositionRecord& rhs) {
              if (lhs.sfc_key != rhs.sfc_key) {
                return lhs.sfc_key < rhs.sfc_key;
              }
              if (lhs.entity_id != rhs.entity_id) {
                return lhs.entity_id < rhs.entity_id;
              }
              return lhs.local_index < rhs.local_index;
            });
  if (records.empty()) {
    return {};
  }

  const std::size_t leaf_count = std::min(max_leaves_per_rank, records.size());
  std::vector<TopDomainLeaf> leaves;
  leaves.reserve(leaf_count);
  const auto item_bounds = [](const DecompositionItem& item) {
    return std::array<double, 6>{
        item.has_spatial_bounds ? item.min_x_comov : item.x_comov,
        item.has_spatial_bounds ? item.max_x_comov : item.x_comov,
        item.has_spatial_bounds ? item.min_y_comov : item.y_comov,
        item.has_spatial_bounds ? item.max_y_comov : item.y_comov,
        item.has_spatial_bounds ? item.min_z_comov : item.z_comov,
        item.has_spatial_bounds ? item.max_z_comov : item.z_comov};
  };
  for (std::size_t leaf_ordinal = 0U; leaf_ordinal < leaf_count; ++leaf_ordinal) {
    const std::size_t begin = (leaf_ordinal * records.size()) / leaf_count;
    const std::size_t end = ((leaf_ordinal + 1U) * records.size()) / leaf_count;
    if (begin == end) {
      continue;
    }
    const DecompositionItem& first = local_items[records[begin].local_index];
    const auto first_bounds = item_bounds(first);
    TopDomainLeaf leaf{
        .owner_rank = owner_rank,
        .decomposition_epoch = decomposition_epoch,
        .sfc_key_begin = records[begin].sfc_key,
        .sfc_key_end = records[end - 1U].sfc_key,
        .min_x_comov = first_bounds[0],
        .max_x_comov = first_bounds[1],
        .min_y_comov = first_bounds[2],
        .max_y_comov = first_bounds[3],
        .min_z_comov = first_bounds[4],
        .max_z_comov = first_bounds[5],
        .periodic_geometry = true,
    };
    std::uint64_t id_hash = 1469598103934665603ULL;
    const auto mix_id = [&id_hash](std::uint64_t value) {
      id_hash ^= value;
      id_hash *= 1099511628211ULL;
    };
    mix_id(static_cast<std::uint64_t>(static_cast<std::uint32_t>(owner_rank)));
    mix_id(static_cast<std::uint64_t>(leaf_ordinal));
    mix_id(leaf.sfc_key_begin);
    mix_id(leaf.sfc_key_end);
    for (std::size_t slot = begin; slot < end; ++slot) {
      const DecompositionItem& item = local_items[records[slot].local_index];
      const auto bounds = item_bounds(item);
      leaf.min_x_comov = std::min(leaf.min_x_comov, bounds[0]);
      leaf.max_x_comov = std::max(leaf.max_x_comov, bounds[1]);
      leaf.min_y_comov = std::min(leaf.min_y_comov, bounds[2]);
      leaf.max_y_comov = std::max(leaf.max_y_comov, bounds[3]);
      leaf.min_z_comov = std::min(leaf.min_z_comov, bounds[4]);
      leaf.max_z_comov = std::max(leaf.max_z_comov, bounds[5]);
      leaf.work_weight += records[slot].weighted_load;
      ++leaf.entity_count;
    }
    leaf.domain_leaf_id = id_hash;
    leaves.push_back(leaf);
  }
  return leaves;
}

struct RuntimeFeedbackNormalization {
  double tree_proxy_sum = 0.0;
  double pm_proxy_sum = 0.0;
  double amr_proxy_sum = 0.0;
  double gas_proxy_sum = 0.0;
  double memory_proxy_sum = 0.0;
  double generic_proxy_sum = 0.0;
  std::size_t entity_count = 0U;
};

struct RuntimeSourceComponentScratch {};

void validateRuntimeDecompositionSource(const RuntimeDecompositionSourceView& source);
[[nodiscard]] std::uint64_t sourceSfcKey(
    const RuntimeDecompositionSourceView& source,
    std::size_t local_index,
    DecompositionEntityKind kind,
    const DecompositionConfig& config);
[[nodiscard]] CompactRuntimeDecompositionRecord makeSourceRecord(
    const RuntimeDecompositionSourceView& source,
    std::size_t local_index,
    DecompositionEntityKind kind,
    const DecompositionConfig& config,
    const DecompositionRuntimeMeasurements& measurements,
    const DecompositionFeedbackCoefficients& coefficients,
    const RuntimeFeedbackNormalization& normalization,
    RuntimeSourceComponentScratch& scratch);

std::vector<TopDomainLeaf> buildAuthoritativeTopDomainLeavesFromSource(
    const RuntimeDecompositionSourceView& source,
    const DecompositionConfig& config,
    int owner_rank,
    std::uint64_t decomposition_epoch,
    std::size_t max_leaves_per_rank) {
  if (config.world_size <= 0 || owner_rank < 0 || owner_rank >= config.world_size ||
      source.world_rank < 0 || source.world_rank >= config.world_size) {
    throw std::invalid_argument("compact source top-domain leaf builder received invalid rank metadata");
  }
  if (max_leaves_per_rank == 0U) {
    throw std::invalid_argument("compact source top-domain leaf builder requires at least one leaf slot per rank");
  }
  validateRuntimeDecompositionSource(source);
  const DecompositionRuntimeMeasurements no_measurements{};
  const DecompositionFeedbackCoefficients default_coefficients{};
  const RuntimeFeedbackNormalization normalization{
      .entity_count = source.localEntityCount()};
  RuntimeSourceComponentScratch scratch;
  std::vector<CompactTopDomainSeedRecord> records;
  records.reserve(source.localEntityCount());
  const auto append_record = [&](std::size_t local_index, DecompositionEntityKind kind) {
    const int owner = kind == DecompositionEntityKind::kParticle
        ? (static_cast<int>(source.particleOwningRank(local_index)))
        : static_cast<int>(source.patch_owning_rank[local_index]);
    if (owner != owner_rank) {
      return;
    }
    const std::uint64_t entity_id = kind == DecompositionEntityKind::kParticle
        ? source.particle_ids[local_index]
        : source.patch_ids[local_index];
    records.push_back(CompactTopDomainSeedRecord{
        .sfc_key = sourceSfcKey(source, local_index, kind, config),
        .entity_id = entity_id,
        .local_index = local_index,
        .kind = kind,
    });
  };
  for (std::size_t particle = 0; particle < source.particle_count; ++particle) {
    append_record(particle, DecompositionEntityKind::kParticle);
  }
  for (const std::uint32_t patch : source.compact_patch_indices) {
    append_record(patch, DecompositionEntityKind::kAmrPatch);
  }
  std::sort(records.begin(), records.end(),
            [](const CompactTopDomainSeedRecord& lhs,
               const CompactTopDomainSeedRecord& rhs) {
              if (lhs.sfc_key != rhs.sfc_key) {
                return lhs.sfc_key < rhs.sfc_key;
              }
              if (lhs.entity_id != rhs.entity_id) {
                return lhs.entity_id < rhs.entity_id;
              }
              return lhs.local_index < rhs.local_index;
            });
  if (records.empty()) {
    return {};
  }
  const auto bounds_for = [&](const CompactTopDomainSeedRecord& record) {
    if (record.kind == DecompositionEntityKind::kParticle) {
      return std::array<double, 6>{
          source.particle_x_comoving[record.local_index],
          source.particle_x_comoving[record.local_index],
          source.particle_y_comoving[record.local_index],
          source.particle_y_comoving[record.local_index],
          source.particle_z_comoving[record.local_index],
          source.particle_z_comoving[record.local_index]};
    }
    const std::uint32_t first_cell = source.patch_first_cells[record.local_index];
    const std::uint32_t cell_count = source.patch_cell_counts[record.local_index];
    std::array<double, 6> bounds{
        source.cell_x_comoving[first_cell],
        source.cell_x_comoving[first_cell],
        source.cell_y_comoving[first_cell],
        source.cell_y_comoving[first_cell],
        source.cell_z_comoving[first_cell],
        source.cell_z_comoving[first_cell]};
    for (std::uint32_t offset = 0; offset < cell_count; ++offset) {
      const std::size_t cell = static_cast<std::size_t>(first_cell) + offset;
      bounds[0] = std::min(bounds[0], source.cell_x_comoving[cell]);
      bounds[1] = std::max(bounds[1], source.cell_x_comoving[cell]);
      bounds[2] = std::min(bounds[2], source.cell_y_comoving[cell]);
      bounds[3] = std::max(bounds[3], source.cell_y_comoving[cell]);
      bounds[4] = std::min(bounds[4], source.cell_z_comoving[cell]);
      bounds[5] = std::max(bounds[5], source.cell_z_comoving[cell]);
    }
    return bounds;
  };
  const auto work_for = [&](const CompactTopDomainSeedRecord& seed) {
    const CompactRuntimeDecompositionRecord record = makeSourceRecord(
        source, seed.local_index, seed.kind, config, no_measurements,
        default_coefficients, normalization, scratch);
    return record.weighted_load;
  };
  const std::size_t leaf_count = std::min(max_leaves_per_rank, records.size());
  std::vector<TopDomainLeaf> leaves;
  leaves.reserve(leaf_count);
  for (std::size_t leaf_ordinal = 0; leaf_ordinal < leaf_count; ++leaf_ordinal) {
    const std::size_t begin = (leaf_ordinal * records.size()) / leaf_count;
    const std::size_t end = ((leaf_ordinal + 1U) * records.size()) / leaf_count;
    if (begin == end) {
      continue;
    }
    const auto first_bounds = bounds_for(records[begin]);
    TopDomainLeaf leaf{
        .owner_rank = owner_rank,
        .decomposition_epoch = decomposition_epoch,
        .sfc_key_begin = records[begin].sfc_key,
        .sfc_key_end = records[end - 1U].sfc_key,
        .min_x_comov = first_bounds[0],
        .max_x_comov = first_bounds[1],
        .min_y_comov = first_bounds[2],
        .max_y_comov = first_bounds[3],
        .min_z_comov = first_bounds[4],
        .max_z_comov = first_bounds[5],
        .periodic_geometry = true,
    };
    std::uint64_t id_hash = 1469598103934665603ULL;
    const auto mix_id = [&id_hash](std::uint64_t value) {
      id_hash ^= value;
      id_hash *= 1099511628211ULL;
    };
    mix_id(static_cast<std::uint64_t>(static_cast<std::uint32_t>(owner_rank)));
    mix_id(static_cast<std::uint64_t>(leaf_ordinal));
    mix_id(leaf.sfc_key_begin);
    mix_id(leaf.sfc_key_end);
    for (std::size_t slot = begin; slot < end; ++slot) {
      const auto bounds = bounds_for(records[slot]);
      leaf.min_x_comov = std::min(leaf.min_x_comov, bounds[0]);
      leaf.max_x_comov = std::max(leaf.max_x_comov, bounds[1]);
      leaf.min_y_comov = std::min(leaf.min_y_comov, bounds[2]);
      leaf.max_y_comov = std::max(leaf.max_y_comov, bounds[3]);
      leaf.min_z_comov = std::min(leaf.min_z_comov, bounds[4]);
      leaf.max_z_comov = std::max(leaf.max_z_comov, bounds[5]);
      leaf.work_weight += work_for(records[slot]);
      ++leaf.entity_count;
    }
    leaf.domain_leaf_id = id_hash;
    leaves.push_back(leaf);
  }
  return leaves;
}

std::vector<TopDomainLeaf> refitAuthoritativeTopDomainLeaves(
    std::span<const TopDomainLeaf> seed_leaves,
    std::span<const double> pos_x_comoving,
    std::span<const double> pos_y_comoving,
    std::span<const double> pos_z_comoving,
    const DecompositionConfig& config,
    int owner_rank,
    std::uint64_t decomposition_epoch,
    TopDomainGeometryRefitDiagnostics* diagnostics) {
  if (config.world_size <= 0 || owner_rank < 0 || owner_rank >= config.world_size) {
    throw std::invalid_argument("top-domain leaf refit received invalid rank metadata");
  }
  if (pos_x_comoving.size() != pos_y_comoving.size() ||
      pos_x_comoving.size() != pos_z_comoving.size()) {
    throw std::invalid_argument("top-domain leaf refit source coordinate extents disagree");
  }
  if (diagnostics != nullptr) {
    *diagnostics = TopDomainGeometryRefitDiagnostics{};
    diagnostics->source_count = static_cast<std::uint64_t>(pos_x_comoving.size());
  }

  // Preserve seed ownership/epoch/SFC intervals; reset bounds for a fresh
  // O(N) AABB expansion. Seeds owned by another rank or stamped with a stale
  // epoch are rejected rather than silently repaired.
  std::vector<TopDomainLeaf> working;
  working.reserve(seed_leaves.size());
  for (const TopDomainLeaf& seed : seed_leaves) {
    if (seed.owner_rank != owner_rank) {
      throw std::invalid_argument("top-domain leaf refit seed is owned by another rank");
    }
    if (seed.decomposition_epoch != decomposition_epoch) {
      throw std::invalid_argument("top-domain leaf refit seed has a stale decomposition epoch");
    }
    TopDomainLeaf leaf = seed;
    leaf.min_x_comov = std::numeric_limits<double>::infinity();
    leaf.max_x_comov = -std::numeric_limits<double>::infinity();
    leaf.min_y_comov = std::numeric_limits<double>::infinity();
    leaf.max_y_comov = -std::numeric_limits<double>::infinity();
    leaf.min_z_comov = std::numeric_limits<double>::infinity();
    leaf.max_z_comov = -std::numeric_limits<double>::infinity();
    leaf.work_weight = 0.0;
    leaf.entity_count = 0U;

    working.push_back(leaf);
  }

  if (working.empty()) {
    if (diagnostics != nullptr && !pos_x_comoving.empty()) {
      // No seed geometry: leave diagnostics empty-leaf so the caller can
      // report no_geometry_installed rather than inventing intervals.
      diagnostics->empty_seed_leaf_count = 0U;
    }
    return {};
  }

  const auto sfc_key_for_source = [&](double x, double y, double z) {
    const std::uint32_t qx = quantize10bit(x, config.domain_x_min_comov, config.domain_x_max_comov);
    const std::uint32_t qy = quantize10bit(y, config.domain_y_min_comov, config.domain_y_max_comov);
    const std::uint32_t qz = quantize10bit(z, config.domain_z_min_comov, config.domain_z_max_comov);
    return mortonKey3d(qx, qy, qz);
  };
  const auto interval_distance =
      [](std::uint64_t key, std::uint64_t begin,
         std::uint64_t end) -> std::uint64_t {
    if (key < begin) {
      return begin - key;
    }
    if (key > end) {
      return key - end;
    }
    return std::uint64_t{0};
  };

  for (std::size_t source_index = 0; source_index < pos_x_comoving.size(); ++source_index) {
    const double x = pos_x_comoving[source_index];
    const double y = pos_y_comoving[source_index];
    const double z = pos_z_comoving[source_index];
    if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z)) {
      throw std::invalid_argument("top-domain leaf refit found non-finite source coordinates");
    }
    const std::uint64_t key = sfc_key_for_source(x, y, z);
    std::size_t selected = working.size();
    for (std::size_t i = 0; i < working.size(); ++i) {
      if (key >= working[i].sfc_key_begin && key <= working[i].sfc_key_end) {
        selected = i;
        break;
      }
    }
    if (selected == working.size()) {
      // Outside every seed interval: assign to the nearest interval so the
      // expanded AABB still covers the source without seam splitting.
      std::uint64_t best_distance = std::numeric_limits<std::uint64_t>::max();
      for (std::size_t i = 0; i < working.size(); ++i) {
        const std::uint64_t distance = interval_distance(
            key, working[i].sfc_key_begin, working[i].sfc_key_end);
        if (selected == working.size() || distance < best_distance ||
            (distance == best_distance &&
             working[i].domain_leaf_id < working[selected].domain_leaf_id)) {
          best_distance = distance;
          selected = i;
        }
      }
      if (selected == working.size()) {
        selected = 0U;
      }
      if (diagnostics != nullptr) {
        ++diagnostics->out_of_seed_range_source_count;
      }
    }
    TopDomainLeaf& leaf = working[selected];
    leaf.min_x_comov = std::min(leaf.min_x_comov, x);
    leaf.max_x_comov = std::max(leaf.max_x_comov, x);
    leaf.min_y_comov = std::min(leaf.min_y_comov, y);
    leaf.max_y_comov = std::max(leaf.max_y_comov, y);
    leaf.min_z_comov = std::min(leaf.min_z_comov, z);
    leaf.max_z_comov = std::max(leaf.max_z_comov, z);
    ++leaf.entity_count;
  }

  std::vector<TopDomainLeaf> refreshed;
  refreshed.reserve(working.size());
  for (TopDomainLeaf& leaf : working) {
    if (leaf.entity_count == 0U) {
      if (diagnostics != nullptr) {
        ++diagnostics->empty_seed_leaf_count;
      }
      // Omit empty leaves so published bounds are always finite and no
      // zero-entity routing packet is manufactured.
      continue;
    }
    refreshed.push_back(leaf);
    if (diagnostics != nullptr) {
      ++diagnostics->refreshed_leaf_count;
    }
  }
  return refreshed;
}

std::uint64_t topDomainGeometryFingerprint(
    std::span<const TopDomainLeaf> leaves) noexcept {
  std::uint64_t hash = 1469598103934665603ULL;
  const auto mix = [&hash](std::uint64_t value) {
    hash ^= value;
    hash *= 1099511628211ULL;
  };
  for (const TopDomainLeaf& leaf : leaves) {
    mix(leaf.domain_leaf_id);
    mix(static_cast<std::uint64_t>(static_cast<std::uint32_t>(std::max(leaf.owner_rank, 0))));
    mix(leaf.decomposition_epoch);
    mix(leaf.sfc_key_begin);
    mix(leaf.sfc_key_end);
    mix(std::bit_cast<std::uint64_t>(leaf.min_x_comov));
    mix(std::bit_cast<std::uint64_t>(leaf.max_x_comov));
    mix(std::bit_cast<std::uint64_t>(leaf.min_y_comov));
    mix(std::bit_cast<std::uint64_t>(leaf.max_y_comov));
    mix(std::bit_cast<std::uint64_t>(leaf.min_z_comov));
    mix(std::bit_cast<std::uint64_t>(leaf.max_z_comov));
    mix(leaf.entity_count);
  }
  return hash;
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

  // Pass 1: one streaming scan for every proxy sum (coefficients unchanged).
  double tree_proxy_sum = 0.0;
  double pm_proxy_sum = 0.0;
  double amr_proxy_sum = 0.0;
  double gas_proxy_sum = 0.0;
  double memory_proxy_sum = 0.0;
  double generic_proxy_sum = 0.0;
  for (const DecompositionItem& item : items) {
    const DecompositionWorkComponents components = effectiveWorkComponents(item);
    tree_proxy_sum += std::max(0.0, components.tree_interaction_cost);
    pm_proxy_sum += std::max(0.0, components.pm_mesh_cost);
    amr_proxy_sum += std::max(0.0, components.amr_patch_cost);
    gas_proxy_sum += std::max(0.0, components.gas_cell_cost);
    memory_proxy_sum += std::max(0.0, components.memory_pressure_cost);
    generic_proxy_sum += std::max(1.0, components.generic_work_cost);
  }

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

  // Pass 2: one streaming distribute/write scan over the same items.
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

void validateRuntimeDecompositionSource(const RuntimeDecompositionSourceView& source) {
  if (source.world_rank < 0) {
    throw std::invalid_argument("runtime decomposition source world_rank must be non-negative");
  }
  if (source.particle_ids.size() != source.particle_count ||
      source.particle_x_comoving.size() != source.particle_count ||
      source.particle_y_comoving.size() != source.particle_count ||
      source.particle_z_comoving.size() != source.particle_count ||
      (!source.particle_metadata_uniform &&
       source.particle_species_tag.size() != source.particle_count)) {
    throw std::invalid_argument("runtime decomposition source particle spans are inconsistent");
  }
  if (!source.particle_metadata_uniform && !source.particle_owning_rank.empty() &&
      source.particle_owning_rank.size() != source.particle_count) {
    throw std::invalid_argument("runtime decomposition source particle ownership span is inconsistent");
  }
  if (!source.active_particle_mask.empty() &&
      source.active_particle_mask.size() != source.particle_count) {
    throw std::invalid_argument("runtime decomposition source active mask span is inconsistent");
  }
  if (source.patch_count != 0U &&
      (source.patch_ids.size() != source.patch_count ||
       source.patch_levels.size() != source.patch_count ||
       source.patch_owning_rank.size() != source.patch_count ||
       source.patch_first_cells.size() != source.patch_count ||
       source.patch_cell_counts.size() != source.patch_count ||
       source.patch_cell_dim_x.size() != source.patch_count ||
       source.patch_cell_dim_y.size() != source.patch_count ||
       source.patch_cell_dim_z.size() != source.patch_count)) {
    throw std::invalid_argument("runtime decomposition source patch spans are inconsistent");
  }
  if (source.cell_x_comoving.size() != source.cell_y_comoving.size() ||
      source.cell_x_comoving.size() != source.cell_z_comoving.size() ||
      source.cell_x_comoving.size() != source.cell_patch_indices.size()) {
    throw std::invalid_argument("runtime decomposition source cell spans are inconsistent");
  }
  // Empty gas offsets + empty indices is the zero-allocation representation
  // for gas-empty/DMO states: every particle has zero associated gas patches.
  if (source.gas_patch_list_offsets.empty()) {
    if (!source.gas_patch_indices.empty()) {
      throw std::invalid_argument(
          "runtime decomposition source gas patch indices exist without offsets");
    }
  } else {
    if (source.gas_patch_list_offsets.size() != source.particle_count + 1U) {
      throw std::invalid_argument("runtime decomposition source gas patch offsets are inconsistent");
    }
    if (source.gas_patch_list_offsets.front() != 0U) {
      throw std::invalid_argument("runtime decomposition source gas patch offset origin is nonzero");
    }
    if (source.gas_patch_list_offsets.back() != source.gas_patch_indices.size()) {
      throw std::invalid_argument("runtime decomposition source gas patch offset endpoint is inconsistent");
    }
    for (std::size_t particle = 0; particle < source.particle_count; ++particle) {
      if (source.gas_patch_list_offsets[particle] >
          source.gas_patch_list_offsets[particle + 1U]) {
        throw std::invalid_argument("runtime decomposition source gas patch offsets are not monotonic");
      }
    }
  }
  for (std::size_t particle = 0; particle < source.particle_count; ++particle) {
    const std::uint32_t species_tag = source.particleSpeciesTag(particle);
    if (species_tag >= source.particle_memory_bytes_by_species.size()) {
      throw std::invalid_argument("runtime decomposition source species tag exceeds memory table");
    }
  }
  for (std::size_t patch_ordinal = 0;
       patch_ordinal < source.compact_patch_indices.size();
       ++patch_ordinal) {
    const std::uint32_t patch_index = source.compact_patch_indices[patch_ordinal];
    if (patch_index >= source.patch_count) {
      throw std::invalid_argument("runtime decomposition source compact patch index is out of range");
    }
    if (source.patch_cell_counts[patch_index] == 0U) {
      throw std::invalid_argument("runtime decomposition source compact patch has no cells");
    }
    const std::uint64_t patch_end =
        static_cast<std::uint64_t>(source.patch_first_cells[patch_index]) +
        source.patch_cell_counts[patch_index];
    if (patch_end > source.cell_x_comoving.size()) {
      throw std::invalid_argument("runtime decomposition source patch cell range is out of range");
    }
  }
  for (std::size_t cell = 0; cell < source.cell_patch_indices.size(); ++cell) {
    if (source.cell_patch_indices[cell] >= source.patch_count) {
      throw std::invalid_argument("runtime decomposition source cell patch index is out of range");
    }
  }
  for (const std::uint32_t patch_index : source.gas_patch_indices) {
    if (patch_index >= source.patch_count) {
      throw std::invalid_argument("runtime decomposition source gas patch incidence is out of range");
    }
  }
  for (const double coordinate : source.particle_x_comoving) {
    if (!std::isfinite(coordinate)) {
      throw std::invalid_argument("runtime decomposition source particle x coordinate is non-finite");
    }
  }
  for (const double coordinate : source.particle_y_comoving) {
    if (!std::isfinite(coordinate)) {
      throw std::invalid_argument("runtime decomposition source particle y coordinate is non-finite");
    }
  }
  for (const double coordinate : source.particle_z_comoving) {
    if (!std::isfinite(coordinate)) {
      throw std::invalid_argument("runtime decomposition source particle z coordinate is non-finite");
    }
  }
  for (const double coordinate : source.cell_x_comoving) {
    if (!std::isfinite(coordinate)) {
      throw std::invalid_argument("runtime decomposition source cell x coordinate is non-finite");
    }
  }
  for (const double coordinate : source.cell_y_comoving) {
    if (!std::isfinite(coordinate)) {
      throw std::invalid_argument("runtime decomposition source cell y coordinate is non-finite");
    }
  }
  for (const double coordinate : source.cell_z_comoving) {
    if (!std::isfinite(coordinate)) {
      throw std::invalid_argument("runtime decomposition source cell z coordinate is non-finite");
    }
  }
}

[[nodiscard]] double sourceAmrPatchCost(
    const RuntimeDecompositionSourceView& source,
    std::size_t particle_index,
    RuntimeSourceComponentScratch&) {
  const std::uint32_t species_tag = source.particleSpeciesTag(particle_index);
  if (species_tag != source.gas_species_tag) {
    return 0.0;
  }
  if (source.gas_patch_list_offsets.empty()) {
    // Empty incidence structure: zero associated gas patches (DMO/gas-empty).
    return 0.0;
  }
  const std::size_t begin = source.gas_patch_list_offsets[particle_index];
  const std::size_t end = source.gas_patch_list_offsets[particle_index + 1U];
  double cost = 0.0;
  for (std::size_t entry = begin; entry < end; ++entry) {
    const std::uint32_t patch_index = source.gas_patch_indices[entry];
    const std::int32_t level = std::max<std::int32_t>(source.patch_levels[patch_index], 0);
    cost += static_cast<double>(source.patch_cell_counts[patch_index]) *
        (1.0 + static_cast<double>(level));
  }
  return cost;
}

[[nodiscard]] DecompositionWorkComponents baseComponentsForSourceParticle(
    const RuntimeDecompositionSourceView& source,
    std::size_t particle_index,
    RuntimeSourceComponentScratch& scratch) {
  const std::uint32_t species_tag = source.particleSpeciesTag(particle_index);
  const std::uint64_t memory_bytes =
      source.particle_memory_bytes_by_species[species_tag];
  const bool active =
      !source.active_particle_mask.empty() &&
      source.active_particle_mask[particle_index] != 0U;
  DecompositionWorkComponents components;
  components.particle_count_cost = 1.0;
  components.gas_cell_cost = species_tag == source.gas_species_tag ? 1.0 : 0.0;
  components.tree_interaction_cost = 1.0;
  components.pm_mesh_cost = 1.0;
  components.amr_patch_cost =
      sourceAmrPatchCost(source, particle_index, scratch);
  components.active_fraction_cost = active ? 1.0 : 0.0;
  components.memory_pressure_cost = static_cast<double>(memory_bytes);
  components.transient_memory_cost =
      species_tag == source.gas_species_tag && active
          ? static_cast<double>(source.gas_transient_memory_bytes_per_cell)
          : 0.0;
  components.source_event_cost =
      (species_tag == source.star_species_tag ||
       species_tag == source.black_hole_species_tag)
          ? 1.0
          : 0.0;
  components.communication_cost = 2.0;
  components.gpu_occupancy_cost = 0.0;
  components.generic_work_cost = 1.0;
  components.has_explicit_components = true;
  return components;
}

[[nodiscard]] DecompositionWorkComponents baseComponentsForSourcePatch(
    const RuntimeDecompositionSourceView& source,
    std::size_t patch_index) {
  const std::uint32_t cell_count = source.patch_cell_counts[patch_index];
  const std::uint64_t memory_bytes = static_cast<std::uint64_t>(cell_count) *
      (sizeof(double) * 8U + sizeof(std::uint32_t) * 2U);
  const std::int32_t level = std::max<std::int32_t>(source.patch_levels[patch_index], 0);
  DecompositionWorkComponents components;
  components.amr_patch_cost = static_cast<double>(cell_count) *
      (1.0 + static_cast<double>(level));
  components.memory_pressure_cost = static_cast<double>(memory_bytes);
  components.transient_memory_cost = static_cast<double>(cell_count) *
      static_cast<double>(source.gas_transient_memory_bytes_per_cell);
  components.communication_cost = 2.0 * (
      static_cast<double>(source.patch_cell_dim_x[patch_index]) *
          source.patch_cell_dim_y[patch_index] +
      static_cast<double>(source.patch_cell_dim_x[patch_index]) *
          source.patch_cell_dim_z[patch_index] +
      static_cast<double>(source.patch_cell_dim_y[patch_index]) *
          source.patch_cell_dim_z[patch_index]);
  components.generic_work_cost = static_cast<double>(cell_count);
  components.has_explicit_components = true;
  return components;
}

[[nodiscard]] DecompositionWorkComponents applySourceMeasuredFeedback(
    const DecompositionWorkComponents& base,
    const DecompositionRuntimeMeasurements& measurements,
    const DecompositionFeedbackCoefficients& coefficients,
    const RuntimeFeedbackNormalization& normalization) {
  if (!measurements.has_measurements) {
    return base;
  }
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
      (measurements.tree_wall_ms_recent + measurements.pm_wall_ms_recent +
       measurements.amr_wall_ms_recent + measurements.hydro_wall_ms_recent);
  const auto entity_count =
      std::max<std::size_t>(normalization.entity_count, 1U);
  const auto distribute =
      [entity_count](double total, double proxy, double sum) {
    if (!(total > 0.0)) {
      return 0.0;
    }
    if (sum > 0.0) {
      return total * std::max(0.0, proxy) / sum;
    }
    return total / static_cast<double>(entity_count);
  };
  DecompositionWorkComponents components = base;
  components.tree_interaction_cost += distribute(
      tree_total, base.tree_interaction_cost, normalization.tree_proxy_sum);
  components.pm_mesh_cost += distribute(
      pm_total, base.pm_mesh_cost, normalization.pm_proxy_sum);
  components.amr_patch_cost += distribute(
      amr_total, base.amr_patch_cost, normalization.amr_proxy_sum);
  components.gas_cell_cost += distribute(
      hydro_total, base.gas_cell_cost, normalization.gas_proxy_sum);
  components.memory_pressure_cost += distribute(
      memory_total, base.memory_pressure_cost, normalization.memory_proxy_sum);
  components.generic_work_cost += distribute(
      wall_total, std::max(1.0, base.generic_work_cost), normalization.generic_proxy_sum);
  components.gpu_occupancy_cost += std::max(0.0, measurements.gpu_kernel_ms_recent) *
      std::max(0.0, measurements.accelerator_occupancy_fraction_recent);
  components.has_explicit_components = true;
  return components;
}

[[nodiscard]] RuntimeFeedbackNormalization computeSourceFeedbackNormalization(
    const RuntimeDecompositionSourceView& source,
    const DecompositionRuntimeMeasurements& measurements,
    const DecompositionFeedbackCoefficients& coefficients) {
  if (!measurements.has_measurements) {
    return RuntimeFeedbackNormalization{.entity_count = source.localEntityCount()};
  }
  if (coefficients.measured_tree_pair < 0.0 || coefficients.measured_pm_cell < 0.0 ||
      coefficients.measured_amr_cell < 0.0 || coefficients.measured_hydro_face < 0.0 ||
      coefficients.measured_wall_ms < 0.0) {
    throw std::invalid_argument("runtime decomposition feedback coefficients must be non-negative");
  }
  RuntimeFeedbackNormalization normalization{.entity_count = source.localEntityCount()};
  RuntimeSourceComponentScratch scratch;
  for (std::size_t particle = 0; particle < source.particle_count; ++particle) {
    const DecompositionWorkComponents components =
        baseComponentsForSourceParticle(source, particle, scratch);
    normalization.tree_proxy_sum += std::max(0.0, components.tree_interaction_cost);
    normalization.pm_proxy_sum += std::max(0.0, components.pm_mesh_cost);
    normalization.amr_proxy_sum += std::max(0.0, components.amr_patch_cost);
    normalization.gas_proxy_sum += std::max(0.0, components.gas_cell_cost);
    normalization.memory_proxy_sum += std::max(0.0, components.memory_pressure_cost);
    normalization.generic_proxy_sum += std::max(1.0, components.generic_work_cost);
  }
  for (const std::uint32_t patch_index : source.compact_patch_indices) {
    const DecompositionWorkComponents components =
        baseComponentsForSourcePatch(source, patch_index);
    normalization.tree_proxy_sum += std::max(0.0, components.tree_interaction_cost);
    normalization.pm_proxy_sum += std::max(0.0, components.pm_mesh_cost);
    normalization.amr_proxy_sum += std::max(0.0, components.amr_patch_cost);
    normalization.gas_proxy_sum += std::max(0.0, components.gas_cell_cost);
    normalization.memory_proxy_sum += std::max(0.0, components.memory_pressure_cost);
    normalization.generic_proxy_sum += std::max(1.0, components.generic_work_cost);
  }
  return normalization;
}

[[nodiscard]] std::uint64_t sourceSfcKey(
    const RuntimeDecompositionSourceView& source,
    std::size_t local_index,
    DecompositionEntityKind kind,
    const DecompositionConfig& config) {
  if (kind == DecompositionEntityKind::kParticle) {
    return sfcKeyForItem(
        DecompositionItem{
            .x_comov = source.particle_x_comoving[local_index],
            .y_comov = source.particle_y_comoving[local_index],
            .z_comov = source.particle_z_comoving[local_index],
        },
        config);
  }
  const std::uint32_t first_cell = source.patch_first_cells[local_index];
  return sfcKeyForItem(
      DecompositionItem{
          .x_comov = source.cell_x_comoving[first_cell],
          .y_comov = source.cell_y_comoving[first_cell],
          .z_comov = source.cell_z_comoving[first_cell],
      },
      config);
}

[[nodiscard]] double sourceWeightedLoad(
    const RuntimeDecompositionSourceView& source,
    const CompactRuntimeDecompositionRecord& record,
    const DecompositionWorkComponents& components,
    const DecompositionConfig& config) {
  if (config.prefer_component_work_model &&
      components.has_explicit_components &&
      hasNonZeroComponentWeight(config.component_weights)) {
    return std::max(0.0, componentWeightedLoad(components, config.component_weights));
  }
  const double owned_particle_term =
      record.kind == DecompositionEntityKind::kParticle ? 1.0 : 0.0;
  const double work_units = record.kind == DecompositionEntityKind::kParticle
      ? 1.0
      : static_cast<double>(source.patch_cell_counts[record.local_index]);
  const double legacy = config.owned_particle_weight * owned_particle_term +
      config.active_target_weight *
          static_cast<double>(record.active_target_count_recent) +
      config.remote_tree_interaction_weight *
          static_cast<double>(record.remote_tree_interactions_recent) +
      config.work_weight * work_units +
      config.memory_weight * static_cast<double>(record.memory_bytes);
  if (legacy > 0.0) {
    return legacy;
  }
  const double component_fallback = componentWeightedLoad(components, config.component_weights);
  if (component_fallback > 0.0) {
    return component_fallback;
  }
  return record.kind == DecompositionEntityKind::kParticle
      ? 1.0
      : std::max(1.0, components.rawTotal());
}

[[nodiscard]] CompactRuntimeDecompositionRecord makeSourceRecord(
    const RuntimeDecompositionSourceView& source,
    std::size_t local_index,
    DecompositionEntityKind kind,
    const DecompositionConfig& config,
    const DecompositionRuntimeMeasurements& measurements,
    const DecompositionFeedbackCoefficients& coefficients,
    const RuntimeFeedbackNormalization& normalization,
    RuntimeSourceComponentScratch& scratch) {
  const bool particle = kind == DecompositionEntityKind::kParticle;
  const std::uint32_t species_tag = particle
      ? source.particleSpeciesTag(local_index)
      : 0U;
  const DecompositionWorkComponents base = particle
      ? baseComponentsForSourceParticle(source, local_index, scratch)
      : baseComponentsForSourcePatch(source, local_index);
  const DecompositionWorkComponents components = applySourceMeasuredFeedback(
      base, measurements, coefficients, normalization);
  const int owner = particle
      ? (static_cast<int>(source.particleOwningRank(local_index)))
      : static_cast<int>(source.patch_owning_rank[local_index]);
  if (owner < 0) {
    throw std::invalid_argument("runtime decomposition source owner is negative");
  }
  const std::uint64_t active = particle &&
          !source.active_particle_mask.empty() &&
          source.active_particle_mask[local_index] != 0U
      ? 1U
      : 0U;
  CompactRuntimeDecompositionRecord record{
      .entity_id = particle ? source.particle_ids[local_index]
                            : source.patch_ids[local_index],
      .sfc_key = sourceSfcKey(source, local_index, kind, config),
      .memory_bytes = particle
          ? source.particle_memory_bytes_by_species[species_tag]
          : static_cast<std::uint64_t>(source.patch_cell_counts[local_index]) *
                (sizeof(double) * 8U + sizeof(std::uint32_t) * 2U),
      .weighted_load = 0.0,
      .local_index = local_index,
      .active_target_count_recent = active,
      .remote_tree_interactions_recent = particle ? 1U : 0U,
      .current_owner_rank = owner,
      .kind = kind,
  };
  record.weighted_load = sourceWeightedLoad(source, record, components, config);
  if (!std::isfinite(record.weighted_load)) {
    throw std::invalid_argument("runtime decomposition source weighted load must be finite");
  }
  return record;
}

[[nodiscard]] std::vector<CompactRuntimeDecompositionRecord> makeSourceRecords(
    const RuntimeDecompositionSourceView& source,
    const DecompositionConfig& config,
    const DecompositionRuntimeMeasurements& measurements,
    const DecompositionFeedbackCoefficients& coefficients) {
  validateRuntimeDecompositionSource(source);
  const RuntimeFeedbackNormalization normalization =
      computeSourceFeedbackNormalization(source, measurements, coefficients);
  RuntimeSourceComponentScratch scratch;
  std::vector<CompactRuntimeDecompositionRecord> records;
  records.reserve(source.localEntityCount());
  for (std::size_t particle = 0; particle < source.particle_count; ++particle) {
    records.push_back(makeSourceRecord(
        source, particle, DecompositionEntityKind::kParticle, config,
        measurements, coefficients, normalization, scratch));
  }
  for (const std::uint32_t patch : source.compact_patch_indices) {
    records.push_back(makeSourceRecord(
        source, patch, DecompositionEntityKind::kAmrPatch, config,
        measurements, coefficients, normalization, scratch));
  }
  return records;
}

[[nodiscard]] DecompositionWorkComponents sourceComponentsForRecord(
    const RuntimeDecompositionSourceView& source,
    const CompactRuntimeDecompositionRecord& record,
    const DecompositionRuntimeMeasurements& measurements,
    const DecompositionFeedbackCoefficients& coefficients,
    const RuntimeFeedbackNormalization& normalization,
    RuntimeSourceComponentScratch& scratch) {
  if (record.kind == DecompositionEntityKind::kParticle) {
    return applySourceMeasuredFeedback(
        baseComponentsForSourceParticle(source, record.local_index, scratch),
        measurements, coefficients, normalization);
  }
  if (record.kind == DecompositionEntityKind::kAmrPatch) {
    return applySourceMeasuredFeedback(
        baseComponentsForSourcePatch(source, record.local_index),
        measurements, coefficients, normalization);
  }
  throw std::invalid_argument("compact runtime source contains an unsupported entity kind");
}

std::vector<CompactRuntimeDecompositionRecord> makeCompactRuntimeDecompositionRecords(
    std::span<const DecompositionItem> items,
    const DecompositionConfig& config) {
  std::vector<CompactRuntimeDecompositionRecord> records;
  records.reserve(items.size());
  for (std::size_t i = 0; i < items.size(); ++i) {
    const DecompositionItem& item = items[i];
    if (item.current_owner_rank < 0 || item.current_owner_rank >= config.world_size) {
      throw std::invalid_argument("decomposition item current_owner_rank is outside runtime world size");
    }
    const double load = weightedLoad(item, config);
    if (!std::isfinite(load)) {
      throw std::invalid_argument("decomposition item weighted load must be finite");
    }
    records.push_back(CompactRuntimeDecompositionRecord{
        .entity_id = item.entity_id,
        .sfc_key = sfcKeyForItem(item, config),
        .memory_bytes = item.memory_bytes,
        .weighted_load = load,
        .local_index = i,
        .active_target_count_recent = item.active_target_count_recent,
        .remote_tree_interactions_recent = item.remote_tree_interactions_recent,
        .current_owner_rank = item.current_owner_rank,
        .kind = item.kind,
    });
  }
  return records;
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
  metrics.peak_memory_bytes_by_rank.assign(rank_count, 0ULL);
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
  metrics.transient_memory_cost_by_rank.assign(rank_count, 0.0);
  metrics.source_event_cost_by_rank.assign(rank_count, 0.0);
  metrics.communication_cost_by_rank.assign(rank_count, 0.0);
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
  std::uint64_t total_peak_memory = 0U;
  for (std::size_t rank = 0; rank < rank_count; ++rank) {
    const bool rank_has_work = metrics.weighted_load_by_rank[rank] > 0.0 || metrics.memory_bytes_by_rank[rank] != 0U;
    metrics.peak_memory_bytes_by_rank[rank] = rank_has_work
        ? checkedUint64Add(metrics.memory_bytes_by_rank[rank], config.rank_transient_reserve_bytes,
                           "current ownership rank peak memory")
        : 0U;
    total_peak_memory = checkedUint64Add(total_peak_memory, metrics.peak_memory_bytes_by_rank[rank],
                                         "current ownership total peak memory");
  }
  const auto max_peak_it = std::max_element(metrics.peak_memory_bytes_by_rank.begin(), metrics.peak_memory_bytes_by_rank.end());
  metrics.max_peak_memory_bytes =
      (max_peak_it == metrics.peak_memory_bytes_by_rank.end()) ? 0U : *max_peak_it;
  const double mean_peak_memory = rank_count == 0U
      ? 0.0
      : static_cast<double>(total_peak_memory) / static_cast<double>(rank_count);
  metrics.peak_memory_imbalance_ratio = mean_peak_memory > 0.0
      ? static_cast<double>(metrics.max_peak_memory_bytes) / mean_peak_memory
      : 0.0;
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
  const bool hard_memory_violated = decomposition_config.max_rank_memory_bytes != 0U &&
      rebalance.current_metrics.max_peak_memory_bytes > decomposition_config.max_rank_memory_bytes;
  if (!load_imbalanced && !memory_imbalanced && !hard_memory_violated) {
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
      if (!hard_memory_violated && max_migrated_load > 0.0 && rebalance.migrated_load + item_load > max_migrated_load &&
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
  rebalance.reason = hard_memory_violated ? "rank_memory_limit" :
      (load_imbalanced && memory_imbalanced ? "load_and_memory_imbalance" :
       (load_imbalanced ? "load_imbalance" : "memory_imbalance"));
  return rebalance;
}

RuntimeRebalancePlan buildDistributedRuntimeRebalancePlan(
    const MpiContext& mpi_context,
    std::span<const DecompositionItem> local_items,
    const DecompositionConfig& decomposition_config,
    const RuntimeRebalanceConfig& rebalance_config) {
  // All rank-local entry validation precedes the first collective. A malformed
  // configuration on one rank must not strand its peers in sample exchange.
  std::exception_ptr entry_failure;
  try {
    if (rebalance_config.world_size <= 0 ||
        decomposition_config.world_size != rebalance_config.world_size ||
        mpi_context.worldSize() != rebalance_config.world_size) {
      throw std::invalid_argument("distributed runtime rebalance world sizes must agree");
    }
    if (rebalance_config.imbalance_trigger_ratio < 1.0 ||
        rebalance_config.memory_trigger_ratio < 1.0 ||
        rebalance_config.max_migrated_load_fraction < 0.0 ||
        rebalance_config.max_migrated_load_fraction > 1.0) {
      throw std::invalid_argument("distributed runtime rebalance thresholds are invalid");
    }
    if (mpi_context.worldSize() > 1 && !mpi_context.isEnabled()) {
      throw std::runtime_error("distributed runtime rebalance requires MPI when world_size > 1");
    }
  } catch (...) {
    entry_failure = std::current_exception();
  }
  if (mpi_context.isEnabled() && mpi_context.worldSize() > 1) {
    mpi_context.rethrowCollectivePreparationFailure(entry_failure,
        "distributed rebalance entry validation");
  } else if (entry_failure != nullptr) {
    std::rethrow_exception(entry_failure);
  }
  if (mpi_context.worldSize() == 1) {
    RuntimeRebalancePlan serial = buildRuntimeRebalancePlan(local_items, decomposition_config, rebalance_config);
    serial.used_distributed_sfc_cuts = false;
    serial.local_entities_considered = static_cast<std::uint64_t>(local_items.size());
    serial.global_entities_considered = static_cast<std::uint64_t>(local_items.size());
    return serial;
  }

  struct LocalKeyedItem {
    std::size_t index = 0;
    SfcCutPoint point{};
    double weighted_load = 0.0;
    DecompositionWorkComponents components{};
  };
  struct CompactCutSample {
    std::uint64_t key = 0;
    std::uint64_t entity_id = 0;
    double represented_load = 0.0;
  };
  static_assert(std::is_trivially_copyable_v<CompactCutSample>);

  std::vector<LocalKeyedItem> keyed;
  std::vector<CompactCutSample> local_samples;
  std::exception_ptr local_preparation_failure;
  try {
    injectMpiTestFault(mpi_context, "sfc_local_preparation");
    keyed.resize(local_items.size());
  for (std::size_t i = 0; i < local_items.size(); ++i) {
    if (local_items[i].current_owner_rank < 0 || local_items[i].current_owner_rank >= decomposition_config.world_size) {
      throw std::invalid_argument("distributed decomposition item current_owner_rank is outside world size");
    }
    keyed[i] = LocalKeyedItem{
        .index = i,
        .point = SfcCutPoint{.key = sfcKeyForItem(local_items[i], decomposition_config),
                             .entity_id = local_items[i].entity_id},
        .weighted_load = weightedLoad(local_items[i], decomposition_config),
        .components = effectiveWorkComponents(local_items[i]),
    };
  }
  std::stable_sort(keyed.begin(), keyed.end(), [](const LocalKeyedItem& lhs, const LocalKeyedItem& rhs) {
    return lessSfcPoint(lhs.point, rhs.point);
  });

    constexpr std::size_t k_samples_per_rank = 256U;
  if (!keyed.empty()) {
    const std::size_t sample_count = std::min(k_samples_per_rank, keyed.size());
    local_samples.reserve(sample_count);
    for (std::size_t sample = 0; sample < sample_count; ++sample) {
      const std::size_t begin = sample * keyed.size() / sample_count;
      const std::size_t end = (sample + 1U) * keyed.size() / sample_count;
      double bucket_load = 0.0;
      for (std::size_t pos = begin; pos < end; ++pos) {
        bucket_load += keyed[pos].weighted_load;
      }
      const LocalKeyedItem& boundary = keyed[end - 1U];
      local_samples.push_back(CompactCutSample{
          .key = boundary.point.key,
          .entity_id = boundary.point.entity_id,
          .represented_load = bucket_load,
      });
    }
  }

  } catch (...) {
    local_preparation_failure = std::current_exception();
  }
  mpi_context.rethrowCollectivePreparationFailure(
      local_preparation_failure, "distributed rebalance local SFC preparation");

  std::vector<CompactCutSample> global_samples;
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  std::size_t local_sample_bytes = 0U;
  std::exception_ptr sample_preparation_failure;
  try {
    injectMpiTestFault(mpi_context, "sfc_cut_sample");
    local_sample_bytes = core::checkedSizeMultiply(
        local_samples.size(), sizeof(CompactCutSample),
        "distributed rebalance local cut sample byte count");
  } catch (...) {
    sample_preparation_failure = std::current_exception();
  }
  mpi_context.rethrowCollectivePreparationFailure(
      sample_preparation_failure,
      "distributed rebalance cut-sample local preparation");

  const auto local_sample_wire = std::span<const std::uint8_t>(
      reinterpret_cast<const std::uint8_t*>(local_samples.data()),
      local_sample_bytes);
  std::vector<std::uint8_t> recv_bytes =
      mpi_context.allgatherBytesBounded(local_sample_wire);

  std::exception_ptr sample_decode_failure;
  try {
    if (recv_bytes.size() % sizeof(CompactCutSample) != 0U) {
      throw std::runtime_error(
          "distributed rebalance cut sample exchange returned partial record bytes");
    }
    global_samples.resize(recv_bytes.size() / sizeof(CompactCutSample));
    if (!recv_bytes.empty()) {
      std::memcpy(global_samples.data(), recv_bytes.data(), recv_bytes.size());
    }
  } catch (...) {
    sample_decode_failure = std::current_exception();
  }
  mpi_context.rethrowCollectivePreparationFailure(
      sample_decode_failure,
      "distributed rebalance cut-sample reassembly");
#else
  throw std::runtime_error("distributed runtime rebalance requires an MPI-enabled build");
#endif

  std::sort(global_samples.begin(), global_samples.end(), [](const CompactCutSample& lhs, const CompactCutSample& rhs) {
    return lessSfcPoint(SfcCutPoint{.key = lhs.key, .entity_id = lhs.entity_id},
                        SfcCutPoint{.key = rhs.key, .entity_id = rhs.entity_id});
  });

  const double global_total_load = mpi_context.allreduceSumDouble(std::accumulate(
      keyed.begin(), keyed.end(), 0.0, [](double acc, const LocalKeyedItem& item) { return acc + item.weighted_load; }));
  const std::uint64_t global_entity_count =
      mpi_context.allreduceSumUint64(static_cast<std::uint64_t>(local_items.size()));
  std::vector<SfcCutPoint> cuts;
  if (!global_samples.empty() && global_total_load > 0.0) {
    const double target_per_rank = global_total_load / static_cast<double>(mpi_context.worldSize());
    double cumulative_sample_load = 0.0;
    std::size_t next_cut_rank = 1U;
    for (const CompactCutSample& sample : global_samples) {
      cumulative_sample_load += std::max(0.0, sample.represented_load);
      if (next_cut_rank < static_cast<std::size_t>(mpi_context.worldSize()) &&
          cumulative_sample_load >= target_per_rank * static_cast<double>(next_cut_rank)) {
        cuts.push_back(SfcCutPoint{.key = sample.key, .entity_id = sample.entity_id});
        ++next_cut_rank;
      }
    }
  }
  while (cuts.size() + 1U < static_cast<std::size_t>(mpi_context.worldSize())) {
    const SfcCutPoint final_point = global_samples.empty()
        ? SfcCutPoint{}
        : SfcCutPoint{.key = global_samples.back().key, .entity_id = global_samples.back().entity_id};
    cuts.push_back(final_point);
  }

  // The sampled work cuts are only proposals. Before materializing migration
  // intents, check their exact global persistent-memory totals. If a proposal
  // is unsafe, refine the cuts with distributed prefix queries rather than
  // gathering the population or incorrectly declaring a feasible case OOM.
  // The refinement uses O(N_local) existing sorted records, O(N_local/256)
  // prefix metadata, and O(P) control state. No per-particle global truth is
  // replicated. Equal (key, ID) records are an indivisible ownership group.
  std::vector<std::uint64_t> memory_block_prefix;
  if (decomposition_config.max_rank_memory_bytes != 0U) {
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
    if (decomposition_config.rank_transient_reserve_bytes >=
        decomposition_config.max_rank_memory_bytes) {
      throw std::invalid_argument(
          "distributed decomposition transient reserve must be smaller than the hard rank memory ceiling");
    }
    const std::uint64_t persistent_limit =
        decomposition_config.max_rank_memory_bytes - decomposition_config.rank_transient_reserve_bytes;
    const std::size_t rank_count = static_cast<std::size_t>(mpi_context.worldSize());
    std::vector<std::uint64_t> proposed_memory;
    std::vector<std::uint64_t> proposed_high;
    std::exception_ptr memory_preparation_failure;
    try {
      injectMpiTestFault(mpi_context, "sfc_memory_preflight");
      proposed_memory.assign(rank_count, 0U);
      proposed_high.assign(rank_count, 0U);
      for (const auto& entry : keyed) {
        const std::size_t owner = static_cast<std::size_t>(
            ownerForSfcPoint(entry.point, cuts, mpi_context.worldSize()));
        proposed_memory[owner] = checkedUint64Add(
            proposed_memory[owner], local_items[entry.index].memory_bytes,
            "distributed proposed rank persistent memory");
      }
    } catch (...) {
      memory_preparation_failure = std::current_exception();
    }
    mpi_context.rethrowCollectivePreparationFailure(
        memory_preparation_failure, "distributed memory cut preflight");
    // Detect overflow before any rank can mistake a wrapped sum for safety.
    for (std::size_t rank = 0U; rank < rank_count; ++rank) {
      proposed_high[rank] = proposed_memory[rank] >> 32U;
      proposed_memory[rank] &= 0xffffffffULL;
    }
    if (MPI_Allreduce(MPI_IN_PLACE, proposed_memory.data(), mpi_context.worldSize(),
                      MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD) != MPI_SUCCESS ||
        MPI_Allreduce(MPI_IN_PLACE, proposed_high.data(), mpi_context.worldSize(),
                      MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD) != MPI_SUCCESS) {
      throw std::runtime_error("distributed memory cut preflight Allreduce failed");
    }
    for (std::size_t rank = 0U; rank < rank_count; ++rank) {
      const std::uint64_t high = checkedUint64Add(
          proposed_high[rank], proposed_memory[rank] >> 32U, "distributed proposed memory carry");
      if (high > 0xffffffffULL) {
        throw std::overflow_error("distributed proposed rank memory exceeds uint64");
      }
      proposed_memory[rank] = (high << 32U) | (proposed_memory[rank] & 0xffffffffULL);
    }
    [[maybe_unused]] std::uint64_t checked_total_memory = 0U;
    for (const std::uint64_t bytes : proposed_memory) {
      checked_total_memory = checkedUint64Add(
          checked_total_memory, bytes, "distributed total persistent memory");
    }
    const bool needs_memory_repair = std::any_of(
        proposed_memory.begin(), proposed_memory.end(),
        [persistent_limit](std::uint64_t bytes) { return bytes > persistent_limit; });
    if (needs_memory_repair) {
      constexpr std::size_t k_memory_prefix_block_size = 256U;
      std::uint64_t local_total_memory = 0U;
      std::exception_ptr prefix_preparation_failure;
      try {
        injectMpiTestFault(mpi_context, "sfc_memory_prefix");
        const std::size_t block_count = keyed.size() / k_memory_prefix_block_size +
            (keyed.size() % k_memory_prefix_block_size != 0U ? 1U : 0U);
        memory_block_prefix.reserve(core::checkedSizeAdd(block_count, 1U, "memory prefix block count"));
        for (std::size_t pos = 0U; pos < keyed.size(); ++pos) {
          if (pos % k_memory_prefix_block_size == 0U) {
            memory_block_prefix.push_back(local_total_memory);
          }
          local_total_memory = checkedUint64Add(
              local_total_memory, local_items[keyed[pos].index].memory_bytes,
              "distributed local persistent memory total");
        }
        memory_block_prefix.push_back(local_total_memory);
      } catch (...) {
        prefix_preparation_failure = std::current_exception();
      }
      mpi_context.rethrowCollectivePreparationFailure(
          prefix_preparation_failure, "distributed memory prefix preparation");

      // Split each exact uint64 sum into 32-bit limbs. With MPI's signed-int
      // rank-count bound, each limb sum fits uint64; reconstruct with checked
      // carry propagation. A wrapped MPI_SUM must never certify a false fit.
      const auto exactGlobalSum = [&](std::uint64_t local_bytes) {
        const std::uint64_t low = mpi_context.allreduceSumUint64(local_bytes & 0xffffffffULL);
        const std::uint64_t high = mpi_context.allreduceSumUint64(local_bytes >> 32U);
        const std::uint64_t high_with_carry = checkedUint64Add(
            high, low >> 32U, "distributed memory sum carry");
        if (high_with_carry > 0xffffffffULL) {
          throw std::overflow_error("distributed total persistent memory exceeds uint64");
        }
        return (high_with_carry << 32U) | (low & 0xffffffffULL);
      };
      const std::uint64_t total_memory = exactGlobalSum(local_total_memory);
      const auto localPrefix = [&](SfcCutPoint point) {
        const auto it = std::upper_bound(
            keyed.begin(), keyed.end(), point,
            [](SfcCutPoint value, const LocalKeyedItem& entry) {
              return lessSfcPoint(value, entry.point);
            });
        const std::size_t end = static_cast<std::size_t>(std::distance(keyed.begin(), it));
        const std::size_t block = end / k_memory_prefix_block_size;
        std::uint64_t bytes = memory_block_prefix[block];
        for (std::size_t pos = block * k_memory_prefix_block_size; pos < end; ++pos) {
          bytes += local_items[keyed[pos].index].memory_bytes;
        }
        return bytes;
      };
      const auto globalPrefix = [&](SfcCutPoint point) {
        return exactGlobalSum(localPrefix(point));
      };
      // Select an actual SFC point whose inclusive prefix first reaches the
      // requested byte threshold. Morton keys occupy 30 bits; IDs occupy 64.
      // Both searches are monotone, fixed-iteration and population-independent.
      const auto firstPrefixAtLeast = [&](std::uint64_t target) {
        std::uint64_t key_low = 0U;
        std::uint64_t key_high = (1ULL << 30U) - 1U;
        while (key_low < key_high) {
          const std::uint64_t mid = key_low + (key_high - key_low) / 2U;
          if (globalPrefix(SfcCutPoint{mid, std::numeric_limits<std::uint64_t>::max()}) >= target) {
            key_high = mid;
          } else {
            key_low = mid + 1U;
          }
        }
        std::uint64_t id_low = 0U;
        std::uint64_t id_high = std::numeric_limits<std::uint64_t>::max();
        while (id_low < id_high) {
          const std::uint64_t mid = id_low + (id_high - id_low) / 2U;
          if (globalPrefix(SfcCutPoint{key_low, mid}) >= target) {
            id_high = mid;
          } else {
            id_low = mid + 1U;
          }
        }
        return SfcCutPoint{key_low, id_low};
      };
      // One bounded rank-level candidate per process. This also handles
      // empty ranks and points with zero memory without inventing a sentinel
      // that might coincide with a real particle ID.
      const auto globalNeighbor = [&](SfcCutPoint point, bool predecessor) {
        struct Candidate { std::uint64_t valid, key, entity_id; };
        Candidate local{};
        const auto it = std::lower_bound(
            keyed.begin(), keyed.end(), point,
            [](const LocalKeyedItem& entry, SfcCutPoint value) {
              return lessSfcPoint(entry.point, value);
            });
        if (predecessor ? it != keyed.begin() : it != keyed.end()) {
          const auto& entry = predecessor ? *std::prev(it) : *it;
          local = Candidate{1U, entry.point.key, entry.point.entity_id};
        }
        const auto wire = std::span<const std::uint8_t>(
            reinterpret_cast<const std::uint8_t*>(&local), sizeof(local));
        const std::vector<std::uint8_t> received = mpi_context.allgatherBytesBounded(wire);
        if (received.size() != core::checkedSizeMultiply(
                rank_count, sizeof(Candidate), "distributed SFC neighbor extent")) {
          throw std::runtime_error("distributed SFC neighbor exchange has an invalid size");
        }
        std::optional<SfcCutPoint> result;
        for (std::size_t rank = 0U; rank < rank_count; ++rank) {
          Candidate candidate{};
          std::memcpy(&candidate, received.data() + rank * sizeof(Candidate), sizeof(Candidate));
          if (candidate.valid == 0U) { continue; }
          const SfcCutPoint value{candidate.key, candidate.entity_id};
          if (!result.has_value() ||
              (predecessor ? lessSfcPoint(*result, value) : lessSfcPoint(value, *result))) {
            result = value;
          }
        }
        return result;
      };
      // Prepare all helper-owned O(P) metadata before the first refinement
      // collective. A rank-local allocation failure must be agreed upon before
      // any peer enters the fixed-sequence prefix reductions.
      std::vector<SfcCutPoint> mandatory_cuts;
      std::vector<SfcCutPoint> repaired_cuts;
      std::exception_ptr repair_preparation_failure;
      try {
        injectMpiTestFault(mpi_context, "sfc_memory_repair_metadata");
        mandatory_cuts.resize(rank_count - 1U);
        repaired_cuts.reserve(rank_count - 1U);
      } catch (...) {
        repair_preparation_failure = std::current_exception();
      }
      mpi_context.rethrowCollectivePreparationFailure(
          repair_preparation_failure, "distributed memory repair metadata preparation");
      const auto first_candidate = globalNeighbor(SfcCutPoint{}, false);
      if (!first_candidate.has_value()) {
        throw std::logic_error("distributed memory repair requires a nonempty population");
      }
      const SfcCutPoint first_point = *first_candidate;
      const SfcCutPoint last_point = global_samples.empty()
          ? first_point
          : SfcCutPoint{global_samples.back().key, global_samples.back().entity_id};
      cuts = internal::repairMemoryConstrainedSfcCuts(
          std::span<const SfcCutPoint>(cuts), persistent_limit, total_memory,
          first_point, last_point, globalPrefix, firstPrefixAtLeast, globalNeighbor,
          std::span<SfcCutPoint>(mandatory_cuts), repaired_cuts);
    }
#endif
  }

  RuntimeRebalancePlan rebalance;
  rebalance.used_distributed_sfc_cuts = true;
  rebalance.local_entities_considered = static_cast<std::uint64_t>(local_items.size());
  rebalance.global_entities_considered = global_entity_count;
  rebalance.target_decomposition.owning_rank_by_item.assign(local_items.size(), 0);
  rebalance.target_decomposition.sorted_indices.reserve(local_items.size());
  rebalance.target_decomposition.ranges_by_rank.assign(static_cast<std::size_t>(mpi_context.worldSize()), RankRange{});
  for (const SfcCutPoint cut : cuts) {
    rebalance.sfc_cut_keys.push_back(cut.key);
    rebalance.sfc_cut_entity_ids.push_back(cut.entity_id);
  }

  auto zero_metrics = [&]() {
    LoadBalanceMetrics metrics;
    const std::size_t rank_count = static_cast<std::size_t>(mpi_context.worldSize());
    metrics.weighted_load_by_rank.assign(rank_count, 0.0);
    metrics.memory_bytes_by_rank.assign(rank_count, 0ULL);
    metrics.peak_memory_bytes_by_rank.assign(rank_count, 0ULL);
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
    metrics.transient_memory_cost_by_rank.assign(rank_count, 0.0);
    metrics.source_event_cost_by_rank.assign(rank_count, 0.0);
    metrics.communication_cost_by_rank.assign(rank_count, 0.0);
    metrics.gpu_occupancy_cost_by_rank.assign(rank_count, 0.0);
    metrics.generic_work_cost_by_rank.assign(rank_count, 0.0);
    return metrics;
  };
  auto accumulate_item = [&](LoadBalanceMetrics& metrics, std::size_t rank, const DecompositionItem& item) {
    metrics.weighted_load_by_rank[rank] += weightedLoad(item, decomposition_config);
    metrics.memory_bytes_by_rank[rank] += item.memory_bytes;
    if (item.kind == DecompositionEntityKind::kParticle) {
      ++metrics.owned_particles_by_rank[rank];
    }
    metrics.active_targets_by_rank[rank] += item.active_target_count_recent;
    metrics.remote_tree_interactions_by_rank[rank] += item.remote_tree_interactions_recent;
    addWorkComponentsToMetrics(metrics, rank, effectiveWorkComponents(item), 1.0);
  };
  [[maybe_unused]] auto finalize_metrics = [&](LoadBalanceMetrics& metrics) {
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
        : static_cast<double>(metrics.total_memory_bytes) / static_cast<double>(metrics.memory_bytes_by_rank.size());
    metrics.memory_imbalance_ratio = (mean_memory > 0.0) ? static_cast<double>(metrics.max_memory_bytes) / mean_memory : 0.0;
    std::uint64_t total_peak_memory = 0U;
    for (std::size_t rank = 0; rank < metrics.memory_bytes_by_rank.size(); ++rank) {
      const bool rank_has_work = metrics.weighted_load_by_rank[rank] > 0.0 || metrics.memory_bytes_by_rank[rank] != 0U;
      metrics.peak_memory_bytes_by_rank[rank] = rank_has_work
          ? checkedUint64Add(metrics.memory_bytes_by_rank[rank], decomposition_config.rank_transient_reserve_bytes,
                             "distributed decomposition rank peak memory")
          : 0U;
      total_peak_memory = checkedUint64Add(total_peak_memory, metrics.peak_memory_bytes_by_rank[rank],
                                           "distributed decomposition total peak memory");
    }
    const auto max_peak_it = std::max_element(metrics.peak_memory_bytes_by_rank.begin(), metrics.peak_memory_bytes_by_rank.end());
    metrics.max_peak_memory_bytes =
        (max_peak_it == metrics.peak_memory_bytes_by_rank.end()) ? 0U : *max_peak_it;
    const double mean_peak_memory = metrics.peak_memory_bytes_by_rank.empty()
        ? 0.0
        : static_cast<double>(total_peak_memory) / static_cast<double>(metrics.peak_memory_bytes_by_rank.size());
    metrics.peak_memory_imbalance_ratio = mean_peak_memory > 0.0
        ? static_cast<double>(metrics.max_peak_memory_bytes) / mean_peak_memory
        : 0.0;
  };

  LoadBalanceMetrics local_current = zero_metrics();
  LoadBalanceMetrics local_target = zero_metrics();
  std::uint64_t local_moved_entities = 0;
  std::uint64_t local_moved_bytes = 0;
  double local_migrated_load = 0.0;
  for (std::size_t sorted_pos = 0; sorted_pos < keyed.size(); ++sorted_pos) {
    const LocalKeyedItem& entry = keyed[sorted_pos];
    const DecompositionItem& item = local_items[entry.index];
    rebalance.target_decomposition.sorted_indices.push_back(entry.index);
    const int new_owner = ownerForSfcPoint(entry.point, cuts, mpi_context.worldSize());
    rebalance.target_decomposition.owning_rank_by_item[entry.index] = new_owner;
    accumulate_item(local_current, static_cast<std::size_t>(item.current_owner_rank), item);
    accumulate_item(local_target, static_cast<std::size_t>(new_owner), item);
    if (item.current_owner_rank != new_owner) {
      ++local_moved_entities;
      local_moved_bytes += item.memory_bytes;
      local_migrated_load += entry.weighted_load;
      if (item.kind == DecompositionEntityKind::kParticle && rebalance_config.allow_particle_migration) {
        rebalance.particle_migrations.push_back(ParticleMigrationIntent{
            .particle_id = item.entity_id,
            .item_index = entry.index,
            .old_owner_rank = item.current_owner_rank,
            .new_owner_rank = new_owner,
            .work_units = entry.weighted_load,
        });
      } else if (item.kind == DecompositionEntityKind::kAmrPatch && rebalance_config.allow_amr_patch_reassignment) {
        rebalance.amr_patch_ownership_updates.push_back(AmrPatchOwnershipUpdate{
            .patch_id = item.entity_id,
            .old_owner_rank = item.current_owner_rank,
            .new_owner_rank = new_owner,
        });
      }
    }
  }

#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  const int metric_rank_count = mpi_context.worldSize();
  auto allreduce_double_vector = [&](std::vector<double>& values) {
    if (MPI_Allreduce(
            MPI_IN_PLACE, values.data(), metric_rank_count,
            MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD) != MPI_SUCCESS) {
      throw std::runtime_error(
          "distributed rebalance double metric Allreduce failed");
    }
  };
  auto allreduce_uint64_vector = [&](std::vector<std::uint64_t>& values) {
    if (MPI_Allreduce(
            MPI_IN_PLACE, values.data(), metric_rank_count,
            MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD) != MPI_SUCCESS) {
      throw std::runtime_error(
          "distributed rebalance uint64 metric Allreduce failed");
    }
  };
  auto reduce_metrics = [&](LoadBalanceMetrics& metrics) {
    allreduce_double_vector(metrics.weighted_load_by_rank);
    allreduce_uint64_vector(metrics.memory_bytes_by_rank);
    allreduce_uint64_vector(metrics.owned_particles_by_rank);
    allreduce_uint64_vector(metrics.active_targets_by_rank);
    allreduce_uint64_vector(metrics.remote_tree_interactions_by_rank);
    allreduce_double_vector(metrics.particle_count_cost_by_rank);
    allreduce_double_vector(metrics.gas_cell_cost_by_rank);
    allreduce_double_vector(metrics.tree_interaction_cost_by_rank);
    allreduce_double_vector(metrics.pm_mesh_cost_by_rank);
    allreduce_double_vector(metrics.amr_patch_cost_by_rank);
    allreduce_double_vector(metrics.active_fraction_cost_by_rank);
    allreduce_double_vector(metrics.memory_pressure_cost_by_rank);
    allreduce_double_vector(metrics.transient_memory_cost_by_rank);
    allreduce_double_vector(metrics.source_event_cost_by_rank);
    allreduce_double_vector(metrics.communication_cost_by_rank);
    allreduce_double_vector(metrics.gpu_occupancy_cost_by_rank);
    allreduce_double_vector(metrics.generic_work_cost_by_rank);
    finalize_metrics(metrics);
  };
  reduce_metrics(local_current);
  reduce_metrics(local_target);
#endif
  rebalance.current_metrics = std::move(local_current);
  rebalance.target_decomposition.metrics = std::move(local_target);

  rebalance.local_entities_moved = local_moved_entities;
  rebalance.global_entities_moved = mpi_context.allreduceSumUint64(local_moved_entities);
  rebalance.local_bytes_moved = local_moved_bytes;
  rebalance.global_bytes_moved = mpi_context.allreduceSumUint64(local_moved_bytes);
  rebalance.migrated_load = mpi_context.allreduceSumDouble(local_migrated_load);
  rebalance.migrated_load_fraction = (global_total_load > 0.0) ? rebalance.migrated_load / global_total_load : 0.0;
  rebalance.local_control_bytes = static_cast<std::uint64_t>(local_samples.size() * sizeof(CompactCutSample));
  rebalance.global_control_bytes = mpi_context.allreduceSumUint64(rebalance.local_control_bytes);
  rebalance.peak_temporary_bytes = static_cast<std::uint64_t>(
      keyed.capacity() * sizeof(LocalKeyedItem) + local_samples.capacity() * sizeof(CompactCutSample) +
      global_samples.capacity() * sizeof(CompactCutSample) +
      memory_block_prefix.capacity() * sizeof(std::uint64_t));
  rebalance.cut_displacement_fraction =
      (global_entity_count > 0U) ? static_cast<double>(rebalance.global_entities_moved) / static_cast<double>(global_entity_count) : 0.0;

  const bool load_imbalanced =
      rebalance.current_metrics.weighted_imbalance_ratio >= rebalance_config.imbalance_trigger_ratio;
  const bool memory_imbalanced =
      rebalance.current_metrics.memory_imbalance_ratio >= rebalance_config.memory_trigger_ratio;
  const bool hard_memory_violated = decomposition_config.max_rank_memory_bytes != 0U &&
      rebalance.current_metrics.max_peak_memory_bytes > decomposition_config.max_rank_memory_bytes;
  const bool target_memory_safe = decomposition_config.max_rank_memory_bytes == 0U ||
      rebalance.target_decomposition.metrics.max_peak_memory_bytes <= decomposition_config.max_rank_memory_bytes;
  if (!target_memory_safe) {
    if (hard_memory_violated) {
      throw std::runtime_error(
          "distributed rebalance cannot find a hard-memory-safe target for the current overloaded ranks");
    }
    rebalance.particle_migrations.clear();
    rebalance.amr_patch_ownership_updates.clear();
    rebalance.should_rebalance = false;
    rebalance.reason = "target_rank_memory_limit_exceeded";
    return rebalance;
  }
  const bool local_has_actionable_migration =
      !rebalance.particle_migrations.empty() ||
      !rebalance.amr_patch_ownership_updates.empty();
  // This vote is collective even on ranks that already have local work. Do
  // not place it behind a short-circuiting local predicate: ranks with no
  // migration would otherwise enter the Allreduce while actionable ranks
  // advance to the payload exchange.
  const std::uint64_t actionable_migration_rank_count =
      mpi_context.allreduceSumUint64(local_has_actionable_migration ? 1ULL : 0ULL);
  if (global_entity_count == 0U) {
    rebalance.reason = "empty_decomposition";
  } else if (!load_imbalanced && !memory_imbalanced && !hard_memory_violated) {
    rebalance.reason = "below_rebalance_threshold";
  } else if (!hard_memory_violated &&
             rebalance.migrated_load_fraction > rebalance_config.max_migrated_load_fraction &&
             rebalance_config.max_migrated_load_fraction < 1.0) {
    rebalance.reason = "migration_fraction_limited";
  } else {
    rebalance.should_rebalance =
        rebalance.global_entities_moved > 0U &&
        actionable_migration_rank_count > 0U;
    if (hard_memory_violated && !rebalance.should_rebalance) {
      throw std::runtime_error(
          "distributed rebalance cannot relieve the existing hard rank memory violation with the allowed migration operations");
    }
    rebalance.reason = hard_memory_violated ? "rank_memory_limit" :
           (load_imbalanced && memory_imbalanced ? "load_and_memory_imbalance" :
            (load_imbalanced ? "load_imbalance" : "memory_imbalance"));
  }
  return rebalance;
}

RuntimeRebalancePlan buildCompactDistributedRuntimeRebalancePlan(
    const MpiContext& mpi_context,
    std::span<const DecompositionItem> local_items,
    const DecompositionConfig& decomposition_config,
    const RuntimeRebalanceConfig& rebalance_config) {
  std::exception_ptr entry_failure;
  try {
    if (rebalance_config.world_size <= 0 ||
        decomposition_config.world_size != rebalance_config.world_size ||
        mpi_context.worldSize() != rebalance_config.world_size) {
      throw std::invalid_argument("distributed runtime rebalance world sizes must agree");
    }
    if (rebalance_config.imbalance_trigger_ratio < 1.0 ||
        rebalance_config.memory_trigger_ratio < 1.0 ||
        rebalance_config.max_migrated_load_fraction < 0.0 ||
        rebalance_config.max_migrated_load_fraction > 1.0) {
      throw std::invalid_argument("distributed runtime rebalance thresholds are invalid");
    }
    if (mpi_context.worldSize() > 1 && !mpi_context.isEnabled()) {
      throw std::runtime_error("distributed runtime rebalance requires MPI when world_size > 1");
    }
  } catch (...) {
    entry_failure = std::current_exception();
  }
  if (mpi_context.isEnabled() && mpi_context.worldSize() > 1) {
    mpi_context.rethrowCollectivePreparationFailure(entry_failure,
        "distributed rebalance entry validation");
  } else if (entry_failure != nullptr) {
    std::rethrow_exception(entry_failure);
  }

  struct CompactCutSample {
    std::uint64_t key = 0;
    std::uint64_t entity_id = 0;
    double represented_load = 0.0;
  };
  static_assert(std::is_trivially_copyable_v<CompactCutSample>);

  RuntimeRebalancePlan rebalance;
   rebalance.used_compact_planner = true;
   rebalance.local_entities_considered = static_cast<std::uint64_t>(local_items.size());
   rebalance.planner_local_entity_count = rebalance.local_entities_considered;


  std::vector<CompactRuntimeDecompositionRecord> records;
  std::vector<CompactCutSample> local_samples;
  std::exception_ptr local_preparation_failure;
  try {
    injectMpiTestFault(mpi_context, "sfc_local_preparation");
    records = makeCompactRuntimeDecompositionRecords(local_items, decomposition_config);
    std::sort(records.begin(), records.end(),
              [](const CompactRuntimeDecompositionRecord& lhs,
                 const CompactRuntimeDecompositionRecord& rhs) {
                if (lhs.sfc_key != rhs.sfc_key) {
                  return lhs.sfc_key < rhs.sfc_key;
                }
                if (lhs.entity_id != rhs.entity_id) {
                  return lhs.entity_id < rhs.entity_id;
                }
                return lhs.local_index < rhs.local_index;
              });
    constexpr std::size_t k_samples_per_rank = 256U;
    if (!records.empty()) {
      const std::size_t sample_count = std::min(k_samples_per_rank, records.size());
      local_samples.reserve(sample_count);
      for (std::size_t sample = 0; sample < sample_count; ++sample) {
        const std::size_t begin = sample * records.size() / sample_count;
        const std::size_t end = (sample + 1U) * records.size() / sample_count;
        double bucket_load = 0.0;
        for (std::size_t pos = begin; pos < end; ++pos) {
          bucket_load += records[pos].weighted_load;
        }
        const CompactRuntimeDecompositionRecord& boundary = records[end - 1U];
        local_samples.push_back(CompactCutSample{
            .key = boundary.sfc_key,
            .entity_id = boundary.entity_id,
            .represented_load = bucket_load,
        });
      }
    }
  } catch (...) {
    local_preparation_failure = std::current_exception();
  }
  if (mpi_context.isEnabled() && mpi_context.worldSize() > 1) {
    mpi_context.rethrowCollectivePreparationFailure(
        local_preparation_failure, "distributed rebalance local SFC preparation");
  } else if (local_preparation_failure != nullptr) {
    std::rethrow_exception(local_preparation_failure);
  }

  if (mpi_context.worldSize() == 1) {
    rebalance.target_decomposition = buildMortonSfcDecomposition(local_items, decomposition_config);
    rebalance.current_metrics = computeCurrentOwnershipLoadBalanceMetrics(local_items, decomposition_config);
    rebalance.global_entities_considered = rebalance.local_entities_considered;
    rebalance.used_distributed_sfc_cuts = false;
    rebalance.planner_record_bytes = static_cast<std::uint64_t>(
        records.capacity() * sizeof(CompactRuntimeDecompositionRecord));
    rebalance.planner_sample_bytes = static_cast<std::uint64_t>(
        local_samples.capacity() * sizeof(CompactCutSample));
    rebalance.planner_prefix_bytes = 0U;
    rebalance.planner_migration_intent_bytes = 0U;
     rebalance.planner_peak_temporary_bytes = rebalance.planner_record_bytes +
         rebalance.planner_sample_bytes;
     rebalance.planner_local_peak_temporary_bytes =
         rebalance.planner_peak_temporary_bytes;
     rebalance.planner_bytes_per_entity = rebalance.local_entities_considered == 0U
         ? 0.0
         : static_cast<double>(rebalance.planner_peak_temporary_bytes) /
               static_cast<double>(rebalance.local_entities_considered);

    rebalance.peak_temporary_bytes = rebalance.planner_peak_temporary_bytes;
    const double total_load = std::accumulate(
        rebalance.current_metrics.weighted_load_by_rank.begin(),
        rebalance.current_metrics.weighted_load_by_rank.end(), 0.0);
    const bool load_imbalanced = rebalance.current_metrics.weighted_imbalance_ratio >=
        rebalance_config.imbalance_trigger_ratio;
    const bool memory_imbalanced = rebalance.current_metrics.memory_imbalance_ratio >=
        rebalance_config.memory_trigger_ratio;
    const bool hard_memory_violated = decomposition_config.max_rank_memory_bytes != 0U &&
        rebalance.current_metrics.max_peak_memory_bytes > decomposition_config.max_rank_memory_bytes;
    if (local_items.empty()) {
      rebalance.reason = "empty_decomposition";
      return rebalance;
    }
    const double max_migrated_load = rebalance_config.max_migrated_load_fraction * std::max(0.0, total_load);
    for (std::size_t item_index = 0; item_index < local_items.size(); ++item_index) {
      const int old_owner = local_items[item_index].current_owner_rank;
      const int new_owner = rebalance.target_decomposition.owning_rank_by_item[item_index];
      if (old_owner < 0 || old_owner == new_owner) {
        continue;
      }
      const double item_load = weightedLoad(local_items[item_index], decomposition_config);
      if (local_items[item_index].kind == DecompositionEntityKind::kParticle &&
          rebalance_config.allow_particle_migration) {
        if (!hard_memory_violated && max_migrated_load > 0.0 &&
            rebalance.migrated_load + item_load > max_migrated_load &&
            !rebalance.particle_migrations.empty()) {
          continue;
        }
        rebalance.particle_migrations.push_back(ParticleMigrationIntent{
            .particle_id = local_items[item_index].entity_id,
            .item_index = item_index,
            .old_owner_rank = old_owner,
            .new_owner_rank = new_owner,
            .work_units = item_load,
        });
        rebalance.migrated_load += item_load;
      } else if (local_items[item_index].kind == DecompositionEntityKind::kAmrPatch &&
                 rebalance_config.allow_amr_patch_reassignment) {
        rebalance.amr_patch_ownership_updates.push_back(AmrPatchOwnershipUpdate{
            .patch_id = local_items[item_index].entity_id,
            .old_owner_rank = old_owner,
            .new_owner_rank = new_owner,
        });
        rebalance.migrated_load += item_load;
      }
    }
    rebalance.should_rebalance = !rebalance.particle_migrations.empty() ||
        !rebalance.amr_patch_ownership_updates.empty();
    rebalance.migrated_load_fraction = (total_load > 0.0) ? (rebalance.migrated_load / total_load) : 0.0;
    if (!load_imbalanced && !memory_imbalanced && !hard_memory_violated) {
      rebalance.reason = "below_rebalance_threshold";
      rebalance.should_rebalance = false;
      rebalance.particle_migrations.clear();
      rebalance.amr_patch_ownership_updates.clear();
      rebalance.migrated_load = 0.0;
      rebalance.migrated_load_fraction = 0.0;
    } else {
      rebalance.reason = hard_memory_violated ? "rank_memory_limit" :
          (load_imbalanced && memory_imbalanced ? "load_and_memory_imbalance" :
           (load_imbalanced ? "load_imbalance" : "memory_imbalance"));
    }
    rebalance.planner_migration_intent_bytes = static_cast<std::uint64_t>(
        rebalance.particle_migrations.capacity() * sizeof(ParticleMigrationIntent) +
        rebalance.amr_patch_ownership_updates.capacity() * sizeof(AmrPatchOwnershipUpdate));
     rebalance.planner_peak_temporary_bytes += rebalance.planner_migration_intent_bytes;
     rebalance.planner_local_peak_temporary_bytes =
         rebalance.planner_peak_temporary_bytes;
     rebalance.peak_temporary_bytes = rebalance.planner_peak_temporary_bytes;
     rebalance.planner_bytes_per_entity = rebalance.local_entities_considered == 0U
         ? 0.0
         : static_cast<double>(rebalance.planner_peak_temporary_bytes) /
               static_cast<double>(rebalance.local_entities_considered);

    return rebalance;
  }

  std::vector<CompactCutSample> global_samples;
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  std::size_t local_sample_bytes = 0U;
  std::exception_ptr sample_preparation_failure;
  try {
    injectMpiTestFault(mpi_context, "sfc_cut_sample");
    local_sample_bytes = core::checkedSizeMultiply(
        local_samples.size(), sizeof(CompactCutSample),
        "compact rebalance local cut sample byte count");
  } catch (...) {
    sample_preparation_failure = std::current_exception();
  }
  mpi_context.rethrowCollectivePreparationFailure(
      sample_preparation_failure,
      "distributed rebalance cut-sample local preparation");

  const auto local_sample_wire = std::span<const std::uint8_t>(
      reinterpret_cast<const std::uint8_t*>(local_samples.data()),
      local_sample_bytes);
  std::vector<std::uint8_t> recv_bytes =
      mpi_context.allgatherBytesBounded(local_sample_wire);

  std::exception_ptr sample_decode_failure;
  try {
    if (recv_bytes.size() % sizeof(CompactCutSample) != 0U) {
      throw std::runtime_error(
          "distributed rebalance cut sample exchange returned partial record bytes");
    }
    global_samples.resize(recv_bytes.size() / sizeof(CompactCutSample));
    if (!recv_bytes.empty()) {
      std::memcpy(global_samples.data(), recv_bytes.data(), recv_bytes.size());
    }
  } catch (...) {
    sample_decode_failure = std::current_exception();
  }
  mpi_context.rethrowCollectivePreparationFailure(
      sample_decode_failure,
      "distributed rebalance cut-sample reassembly");
#else
  throw std::runtime_error("distributed runtime rebalance requires an MPI-enabled build");
#endif

  std::sort(global_samples.begin(), global_samples.end(),
            [](const CompactCutSample& lhs, const CompactCutSample& rhs) {
              return lessSfcPoint(SfcCutPoint{.key = lhs.key, .entity_id = lhs.entity_id},
                                  SfcCutPoint{.key = rhs.key, .entity_id = rhs.entity_id});
            });

  const double global_total_load = mpi_context.allreduceSumDouble(std::accumulate(
      records.begin(), records.end(), 0.0,
      [](double acc, const CompactRuntimeDecompositionRecord& item) {
        return acc + item.weighted_load;
      }));
  const std::uint64_t global_entity_count =
      mpi_context.allreduceSumUint64(static_cast<std::uint64_t>(local_items.size()));
  rebalance.global_entities_considered = global_entity_count;

  std::vector<SfcCutPoint> cuts;
  if (!global_samples.empty() && global_total_load > 0.0) {
    const double target_per_rank = global_total_load / static_cast<double>(mpi_context.worldSize());
    double cumulative_sample_load = 0.0;
    std::size_t next_cut_rank = 1U;
    for (const CompactCutSample& sample : global_samples) {
      cumulative_sample_load += std::max(0.0, sample.represented_load);
      if (next_cut_rank < static_cast<std::size_t>(mpi_context.worldSize()) &&
          cumulative_sample_load >= target_per_rank * static_cast<double>(next_cut_rank)) {
        cuts.push_back(SfcCutPoint{.key = sample.key, .entity_id = sample.entity_id});
        ++next_cut_rank;
      }
    }
  }
  while (cuts.size() + 1U < static_cast<std::size_t>(mpi_context.worldSize())) {
    const SfcCutPoint final_point = global_samples.empty()
        ? SfcCutPoint{}
        : SfcCutPoint{.key = global_samples.back().key, .entity_id = global_samples.back().entity_id};
    cuts.push_back(final_point);
  }

  std::vector<std::uint64_t> memory_block_prefix;
  if (decomposition_config.max_rank_memory_bytes != 0U) {
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
    if (decomposition_config.rank_transient_reserve_bytes >=
        decomposition_config.max_rank_memory_bytes) {
      throw std::invalid_argument(
          "distributed decomposition transient reserve must be smaller than the hard rank memory ceiling");
    }
    const std::uint64_t persistent_limit =
        decomposition_config.max_rank_memory_bytes - decomposition_config.rank_transient_reserve_bytes;
    const std::size_t rank_count = static_cast<std::size_t>(mpi_context.worldSize());
    std::vector<std::uint64_t> proposed_memory;
    std::vector<std::uint64_t> proposed_high;
    std::exception_ptr memory_preparation_failure;
    try {
      injectMpiTestFault(mpi_context, "sfc_memory_preflight");
      proposed_memory.assign(rank_count, 0U);
      proposed_high.assign(rank_count, 0U);
      for (const auto& entry : records) {
        const SfcCutPoint point{.key = entry.sfc_key, .entity_id = entry.entity_id};
        const std::size_t owner = static_cast<std::size_t>(
            ownerForSfcPoint(point, cuts, mpi_context.worldSize()));
        proposed_memory[owner] = checkedUint64Add(
            proposed_memory[owner], entry.memory_bytes,
            "compact proposed rank persistent memory");
      }
    } catch (...) {
      memory_preparation_failure = std::current_exception();
    }
    mpi_context.rethrowCollectivePreparationFailure(
        memory_preparation_failure, "distributed memory cut preflight");
    for (std::size_t rank = 0U; rank < rank_count; ++rank) {
      proposed_high[rank] = proposed_memory[rank] >> 32U;
      proposed_memory[rank] &= 0xffffffffULL;
    }
    if (MPI_Allreduce(MPI_IN_PLACE, proposed_memory.data(), mpi_context.worldSize(),
                      MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD) != MPI_SUCCESS ||
        MPI_Allreduce(MPI_IN_PLACE, proposed_high.data(), mpi_context.worldSize(),
                      MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD) != MPI_SUCCESS) {
      throw std::runtime_error("distributed memory cut preflight Allreduce failed");
    }
    for (std::size_t rank = 0U; rank < rank_count; ++rank) {
      const std::uint64_t high = checkedUint64Add(
          proposed_high[rank], proposed_memory[rank] >> 32U, "compact proposed memory carry");
      if (high > 0xffffffffULL) {
        throw std::overflow_error("compact proposed rank memory exceeds uint64");
      }
      proposed_memory[rank] = (high << 32U) | (proposed_memory[rank] & 0xffffffffULL);
    }
    [[maybe_unused]] std::uint64_t checked_total_memory = 0U;
    for (const std::uint64_t bytes : proposed_memory) {
      checked_total_memory = checkedUint64Add(
          checked_total_memory, bytes, "compact total persistent memory");
    }
    const bool needs_memory_repair = std::any_of(
        proposed_memory.begin(), proposed_memory.end(),
        [persistent_limit](std::uint64_t bytes) { return bytes > persistent_limit; });
    if (needs_memory_repair) {
      constexpr std::size_t k_memory_prefix_block_size = 256U;
      std::uint64_t local_total_memory = 0U;
      std::exception_ptr prefix_preparation_failure;
      try {
        injectMpiTestFault(mpi_context, "sfc_memory_prefix");
        const std::size_t block_count = records.size() / k_memory_prefix_block_size +
            (records.size() % k_memory_prefix_block_size != 0U ? 1U : 0U);
        memory_block_prefix.reserve(core::checkedSizeAdd(block_count, 1U, "memory prefix block count"));
        for (std::size_t pos = 0U; pos < records.size(); ++pos) {
          if (pos % k_memory_prefix_block_size == 0U) {
            memory_block_prefix.push_back(local_total_memory);
          }
          local_total_memory = checkedUint64Add(
              local_total_memory, records[pos].memory_bytes,
              "compact local persistent memory total");
        }
        memory_block_prefix.push_back(local_total_memory);
      } catch (...) {
        prefix_preparation_failure = std::current_exception();
      }
      mpi_context.rethrowCollectivePreparationFailure(
          prefix_preparation_failure, "distributed memory prefix preparation");

      const auto exactGlobalSum = [&](std::uint64_t local_bytes) {
        const std::uint64_t low = mpi_context.allreduceSumUint64(local_bytes & 0xffffffffULL);
        const std::uint64_t high = mpi_context.allreduceSumUint64(local_bytes >> 32U);
        const std::uint64_t high_with_carry = checkedUint64Add(
            high, low >> 32U, "compact memory sum carry");
        if (high_with_carry > 0xffffffffULL) {
          throw std::overflow_error("compact total persistent memory exceeds uint64");
        }
        return (high_with_carry << 32U) | (low & 0xffffffffULL);
      };
      const std::uint64_t total_memory = exactGlobalSum(local_total_memory);
      const auto localPrefix = [&](SfcCutPoint point) {
        const auto it = std::upper_bound(
            records.begin(), records.end(), point,
            [](SfcCutPoint value, const CompactRuntimeDecompositionRecord& entry) {
              return lessSfcPoint(value, SfcCutPoint{entry.sfc_key, entry.entity_id});
            });
        const std::size_t end = static_cast<std::size_t>(std::distance(records.begin(), it));
        const std::size_t block = end / k_memory_prefix_block_size;
        std::uint64_t bytes = memory_block_prefix[block];
        for (std::size_t pos = block * k_memory_prefix_block_size; pos < end; ++pos) {
          bytes += records[pos].memory_bytes;
        }
        return bytes;
      };
      const auto globalPrefix = [&](SfcCutPoint point) {
        return exactGlobalSum(localPrefix(point));
      };
      const auto firstPrefixAtLeast = [&](std::uint64_t target) {
        std::uint64_t key_low = 0U;
        std::uint64_t key_high = (1ULL << 30U) - 1U;
        while (key_low < key_high) {
          const std::uint64_t mid = key_low + (key_high - key_low) / 2U;
          if (globalPrefix(SfcCutPoint{mid, std::numeric_limits<std::uint64_t>::max()}) >= target) {
            key_high = mid;
          } else {
            key_low = mid + 1U;
          }
        }
        std::uint64_t id_low = 0U;
        std::uint64_t id_high = std::numeric_limits<std::uint64_t>::max();
        while (id_low < id_high) {
          const std::uint64_t mid = id_low + (id_high - id_low) / 2U;
          if (globalPrefix(SfcCutPoint{key_low, mid}) >= target) {
            id_high = mid;
          } else {
            id_low = mid + 1U;
          }
        }
        return SfcCutPoint{key_low, id_low};
      };
      const auto globalNeighbor = [&](SfcCutPoint point, bool predecessor) {
        struct Candidate { std::uint64_t valid, key, entity_id; };
        Candidate local{};
        const auto it = std::lower_bound(
            records.begin(), records.end(), point,
            [](const CompactRuntimeDecompositionRecord& entry, SfcCutPoint value) {
              return lessSfcPoint(SfcCutPoint{entry.sfc_key, entry.entity_id}, value);
            });
        if (predecessor ? it != records.begin() : it != records.end()) {
          const auto& entry = predecessor ? *std::prev(it) : *it;
          local = Candidate{1U, entry.sfc_key, entry.entity_id};
        }
        const auto wire = std::span<const std::uint8_t>(
            reinterpret_cast<const std::uint8_t*>(&local), sizeof(local));
        const std::vector<std::uint8_t> received = mpi_context.allgatherBytesBounded(wire);
        if (received.size() != core::checkedSizeMultiply(
                rank_count, sizeof(Candidate), "compact SFC neighbor extent")) {
          throw std::runtime_error("distributed SFC neighbor exchange has an invalid size");
        }
        std::optional<SfcCutPoint> result;
        for (std::size_t rank = 0U; rank < rank_count; ++rank) {
          Candidate candidate{};
          std::memcpy(&candidate, received.data() + rank * sizeof(Candidate), sizeof(Candidate));
          if (candidate.valid == 0U) { continue; }
          const SfcCutPoint value{candidate.key, candidate.entity_id};
          if (!result.has_value() ||
              (predecessor ? lessSfcPoint(*result, value) : lessSfcPoint(value, *result))) {
            result = value;
          }
        }
        return result;
      };
      std::vector<SfcCutPoint> mandatory_cuts;
      std::vector<SfcCutPoint> repaired_cuts;
      std::exception_ptr repair_preparation_failure;
      try {
        injectMpiTestFault(mpi_context, "sfc_memory_repair_metadata");
        mandatory_cuts.resize(rank_count - 1U);
        repaired_cuts.reserve(rank_count - 1U);
      } catch (...) {
        repair_preparation_failure = std::current_exception();
      }
      mpi_context.rethrowCollectivePreparationFailure(
          repair_preparation_failure, "distributed memory repair metadata preparation");
      const auto first_candidate = globalNeighbor(SfcCutPoint{}, false);
      if (!first_candidate.has_value()) {
        throw std::logic_error("distributed memory repair requires a nonempty population");
      }
      const SfcCutPoint first_point = *first_candidate;
      const SfcCutPoint last_point = global_samples.empty()
          ? first_point
          : SfcCutPoint{global_samples.back().key, global_samples.back().entity_id};
      cuts = internal::repairMemoryConstrainedSfcCuts(
          std::span<const SfcCutPoint>(cuts), persistent_limit, total_memory,
          first_point, last_point, globalPrefix, firstPrefixAtLeast, globalNeighbor,
          std::span<SfcCutPoint>(mandatory_cuts), repaired_cuts);
    }
#endif
  }

  rebalance.used_distributed_sfc_cuts = true;
  for (const SfcCutPoint cut : cuts) {
    rebalance.sfc_cut_keys.push_back(cut.key);
    rebalance.sfc_cut_entity_ids.push_back(cut.entity_id);
  }

  auto zero_metrics = [&]() {
    LoadBalanceMetrics metrics;
    const std::size_t rank_count = static_cast<std::size_t>(mpi_context.worldSize());
    metrics.weighted_load_by_rank.assign(rank_count, 0.0);
    metrics.memory_bytes_by_rank.assign(rank_count, 0ULL);
    metrics.peak_memory_bytes_by_rank.assign(rank_count, 0ULL);
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
    metrics.transient_memory_cost_by_rank.assign(rank_count, 0.0);
    metrics.source_event_cost_by_rank.assign(rank_count, 0.0);
    metrics.communication_cost_by_rank.assign(rank_count, 0.0);
    metrics.gpu_occupancy_cost_by_rank.assign(rank_count, 0.0);
    metrics.generic_work_cost_by_rank.assign(rank_count, 0.0);
    return metrics;
  };
  auto accumulate_item = [&](LoadBalanceMetrics& metrics, std::size_t rank, const DecompositionItem& item) {
    metrics.weighted_load_by_rank[rank] += weightedLoad(item, decomposition_config);
    metrics.memory_bytes_by_rank[rank] += item.memory_bytes;
    if (item.kind == DecompositionEntityKind::kParticle) {
      ++metrics.owned_particles_by_rank[rank];
    }
    metrics.active_targets_by_rank[rank] += item.active_target_count_recent;
    metrics.remote_tree_interactions_by_rank[rank] += item.remote_tree_interactions_recent;
    addWorkComponentsToMetrics(metrics, rank, effectiveWorkComponents(item), 1.0);
  };
  [[maybe_unused]] auto finalize_metrics = [&](LoadBalanceMetrics& metrics) {
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
        : static_cast<double>(metrics.total_memory_bytes) / static_cast<double>(metrics.memory_bytes_by_rank.size());
    metrics.memory_imbalance_ratio = (mean_memory > 0.0) ? static_cast<double>(metrics.max_memory_bytes) / mean_memory : 0.0;
    std::uint64_t total_peak_memory = 0U;
    for (std::size_t rank = 0; rank < metrics.memory_bytes_by_rank.size(); ++rank) {
      const bool rank_has_work = metrics.weighted_load_by_rank[rank] > 0.0 || metrics.memory_bytes_by_rank[rank] != 0U;
      metrics.peak_memory_bytes_by_rank[rank] = rank_has_work
          ? checkedUint64Add(metrics.memory_bytes_by_rank[rank], decomposition_config.rank_transient_reserve_bytes,
                             "compact decomposition rank peak memory")
          : 0U;
      total_peak_memory = checkedUint64Add(total_peak_memory, metrics.peak_memory_bytes_by_rank[rank],
                                           "compact decomposition total peak memory");
    }
    const auto max_peak_it = std::max_element(metrics.peak_memory_bytes_by_rank.begin(), metrics.peak_memory_bytes_by_rank.end());
    metrics.max_peak_memory_bytes =
        (max_peak_it == metrics.peak_memory_bytes_by_rank.end()) ? 0U : *max_peak_it;
    const double mean_peak_memory = metrics.peak_memory_bytes_by_rank.empty()
        ? 0.0
        : static_cast<double>(total_peak_memory) / static_cast<double>(metrics.peak_memory_bytes_by_rank.size());
    metrics.peak_memory_imbalance_ratio = mean_peak_memory > 0.0
        ? static_cast<double>(metrics.max_peak_memory_bytes) / mean_peak_memory
        : 0.0;
  };

  LoadBalanceMetrics local_current = zero_metrics();
  LoadBalanceMetrics local_target = zero_metrics();
  std::uint64_t local_moved_entities = 0;
  std::uint64_t local_moved_bytes = 0;
  double local_migrated_load = 0.0;
  for (const CompactRuntimeDecompositionRecord& entry : records) {
    const DecompositionItem& item = local_items[entry.local_index];
    const SfcCutPoint point{.key = entry.sfc_key, .entity_id = entry.entity_id};
    const int new_owner = ownerForSfcPoint(point, cuts, mpi_context.worldSize());
    accumulate_item(local_current, static_cast<std::size_t>(item.current_owner_rank), item);
    accumulate_item(local_target, static_cast<std::size_t>(new_owner), item);
    if (item.current_owner_rank != new_owner) {
      ++local_moved_entities;
      local_moved_bytes += entry.memory_bytes;
      local_migrated_load += entry.weighted_load;
      if (item.kind == DecompositionEntityKind::kParticle && rebalance_config.allow_particle_migration) {
        rebalance.particle_migrations.push_back(ParticleMigrationIntent{
            .particle_id = item.entity_id,
            .item_index = entry.local_index,
            .old_owner_rank = item.current_owner_rank,
            .new_owner_rank = new_owner,
            .work_units = entry.weighted_load,
        });
      } else if (item.kind == DecompositionEntityKind::kAmrPatch && rebalance_config.allow_amr_patch_reassignment) {
        rebalance.amr_patch_ownership_updates.push_back(AmrPatchOwnershipUpdate{
            .patch_id = item.entity_id,
            .old_owner_rank = item.current_owner_rank,
            .new_owner_rank = new_owner,
        });
      }
    }
  }

#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  const int metric_rank_count = mpi_context.worldSize();
  auto allreduce_double_vector = [&](std::vector<double>& values) {
    if (MPI_Allreduce(
            MPI_IN_PLACE, values.data(), metric_rank_count,
            MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD) != MPI_SUCCESS) {
      throw std::runtime_error(
          "compact rebalance double metric Allreduce failed");
    }
  };
  auto allreduce_uint64_vector = [&](std::vector<std::uint64_t>& values) {
    if (MPI_Allreduce(
            MPI_IN_PLACE, values.data(), metric_rank_count,
            MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD) != MPI_SUCCESS) {
      throw std::runtime_error(
          "compact rebalance uint64 metric Allreduce failed");
    }
  };
  auto reduce_metrics = [&](LoadBalanceMetrics& metrics) {
    allreduce_double_vector(metrics.weighted_load_by_rank);
    allreduce_uint64_vector(metrics.memory_bytes_by_rank);
    allreduce_uint64_vector(metrics.owned_particles_by_rank);
    allreduce_uint64_vector(metrics.active_targets_by_rank);
    allreduce_uint64_vector(metrics.remote_tree_interactions_by_rank);
    allreduce_double_vector(metrics.particle_count_cost_by_rank);
    allreduce_double_vector(metrics.gas_cell_cost_by_rank);
    allreduce_double_vector(metrics.tree_interaction_cost_by_rank);
    allreduce_double_vector(metrics.pm_mesh_cost_by_rank);
    allreduce_double_vector(metrics.amr_patch_cost_by_rank);
    allreduce_double_vector(metrics.active_fraction_cost_by_rank);
    allreduce_double_vector(metrics.memory_pressure_cost_by_rank);
    allreduce_double_vector(metrics.transient_memory_cost_by_rank);
    allreduce_double_vector(metrics.source_event_cost_by_rank);
    allreduce_double_vector(metrics.communication_cost_by_rank);
    allreduce_double_vector(metrics.gpu_occupancy_cost_by_rank);
    allreduce_double_vector(metrics.generic_work_cost_by_rank);
    finalize_metrics(metrics);
  };
  reduce_metrics(local_current);
  reduce_metrics(local_target);
#endif
  rebalance.current_metrics = std::move(local_current);
  rebalance.target_decomposition.metrics = std::move(local_target);

  rebalance.local_entities_moved = local_moved_entities;
  rebalance.global_entities_moved = mpi_context.allreduceSumUint64(local_moved_entities);
  rebalance.local_bytes_moved = local_moved_bytes;
  rebalance.global_bytes_moved = mpi_context.allreduceSumUint64(local_moved_bytes);
  rebalance.migrated_load = mpi_context.allreduceSumDouble(local_migrated_load);
  rebalance.migrated_load_fraction = (global_total_load > 0.0) ? rebalance.migrated_load / global_total_load : 0.0;
  rebalance.local_control_bytes = static_cast<std::uint64_t>(local_samples.size() * sizeof(CompactCutSample));
  rebalance.global_control_bytes = mpi_context.allreduceSumUint64(rebalance.local_control_bytes);
  rebalance.planner_record_bytes = static_cast<std::uint64_t>(
      records.capacity() * sizeof(CompactRuntimeDecompositionRecord));
  rebalance.planner_sample_bytes = static_cast<std::uint64_t>(
      (local_samples.capacity() + global_samples.capacity()) * sizeof(CompactCutSample));
  rebalance.planner_prefix_bytes = static_cast<std::uint64_t>(
      memory_block_prefix.capacity() * sizeof(std::uint64_t));
  rebalance.planner_migration_intent_bytes = static_cast<std::uint64_t>(
      rebalance.particle_migrations.capacity() * sizeof(ParticleMigrationIntent) +
      rebalance.amr_patch_ownership_updates.capacity() * sizeof(AmrPatchOwnershipUpdate));
   rebalance.planner_peak_temporary_bytes = rebalance.planner_record_bytes +
       rebalance.planner_sample_bytes + rebalance.planner_prefix_bytes +
       rebalance.planner_migration_intent_bytes;
   rebalance.planner_local_peak_temporary_bytes =
       rebalance.planner_peak_temporary_bytes;
   rebalance.planner_bytes_per_entity = rebalance.local_entities_considered == 0U
       ? 0.0
       : static_cast<double>(rebalance.planner_peak_temporary_bytes) /
             static_cast<double>(rebalance.local_entities_considered);

  rebalance.peak_temporary_bytes = rebalance.planner_peak_temporary_bytes;
  rebalance.cut_displacement_fraction =
      (global_entity_count > 0U) ? static_cast<double>(rebalance.global_entities_moved) / static_cast<double>(global_entity_count) : 0.0;

  const bool load_imbalanced =
      rebalance.current_metrics.weighted_imbalance_ratio >= rebalance_config.imbalance_trigger_ratio;
  const bool memory_imbalanced =
      rebalance.current_metrics.memory_imbalance_ratio >= rebalance_config.memory_trigger_ratio;
  const bool hard_memory_violated = decomposition_config.max_rank_memory_bytes != 0U &&
      rebalance.current_metrics.max_peak_memory_bytes > decomposition_config.max_rank_memory_bytes;
  const bool target_memory_safe = decomposition_config.max_rank_memory_bytes == 0U ||
      rebalance.target_decomposition.metrics.max_peak_memory_bytes <= decomposition_config.max_rank_memory_bytes;
  if (!target_memory_safe) {
    if (hard_memory_violated) {
      throw std::runtime_error(
          "distributed rebalance cannot find a hard-memory-safe target for the current overloaded ranks");
    }
    rebalance.particle_migrations.clear();
    rebalance.amr_patch_ownership_updates.clear();
    rebalance.should_rebalance = false;
    rebalance.reason = "target_rank_memory_limit_exceeded";
    rebalance.planner_migration_intent_bytes = 0U;
     rebalance.planner_peak_temporary_bytes = rebalance.planner_record_bytes +
         rebalance.planner_sample_bytes + rebalance.planner_prefix_bytes +
         rebalance.planner_migration_intent_bytes;
     rebalance.planner_local_peak_temporary_bytes =
         rebalance.planner_peak_temporary_bytes;
     rebalance.planner_bytes_per_entity = rebalance.local_entities_considered == 0U
         ? 0.0
         : static_cast<double>(rebalance.planner_peak_temporary_bytes) /
               static_cast<double>(rebalance.local_entities_considered);

    rebalance.peak_temporary_bytes = rebalance.planner_peak_temporary_bytes;
    return rebalance;
  }
  const bool local_has_actionable_migration =
      !rebalance.particle_migrations.empty() ||
      !rebalance.amr_patch_ownership_updates.empty();
  const std::uint64_t actionable_migration_rank_count =
      mpi_context.allreduceSumUint64(local_has_actionable_migration ? 1ULL : 0ULL);
  if (global_entity_count == 0U) {
    rebalance.reason = "empty_decomposition";
  } else if (!load_imbalanced && !memory_imbalanced && !hard_memory_violated) {
    rebalance.reason = "below_rebalance_threshold";
  } else if (!hard_memory_violated &&
             rebalance.migrated_load_fraction > rebalance_config.max_migrated_load_fraction &&
             rebalance_config.max_migrated_load_fraction < 1.0) {
    rebalance.reason = "migration_fraction_limited";
  } else {
    rebalance.should_rebalance =
        rebalance.global_entities_moved > 0U &&
        actionable_migration_rank_count > 0U;
    if (hard_memory_violated && !rebalance.should_rebalance) {
      throw std::runtime_error(
          "distributed rebalance cannot relieve the existing hard rank memory violation with the allowed migration operations");
    }
    rebalance.reason = hard_memory_violated ? "rank_memory_limit" :
           (load_imbalanced && memory_imbalanced ? "load_and_memory_imbalance" :
            (load_imbalanced ? "load_imbalance" : "memory_imbalance"));
  }
  return rebalance;
}

void initializeCompactMetrics(
    LoadBalanceMetrics& metrics,
    std::size_t rank_count) {
  metrics.weighted_load_by_rank.assign(rank_count, 0.0);
  metrics.memory_bytes_by_rank.assign(rank_count, 0ULL);
  metrics.peak_memory_bytes_by_rank.assign(rank_count, 0ULL);
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
  metrics.transient_memory_cost_by_rank.assign(rank_count, 0.0);
  metrics.source_event_cost_by_rank.assign(rank_count, 0.0);
  metrics.communication_cost_by_rank.assign(rank_count, 0.0);
  metrics.gpu_occupancy_cost_by_rank.assign(rank_count, 0.0);
  metrics.generic_work_cost_by_rank.assign(rank_count, 0.0);
}

void finalizeCompactMetrics(
    LoadBalanceMetrics& metrics,
    const DecompositionConfig& config) {
  const auto max_load_it = std::max_element(
      metrics.weighted_load_by_rank.begin(), metrics.weighted_load_by_rank.end());
  metrics.max_weighted_load = max_load_it == metrics.weighted_load_by_rank.end()
      ? 0.0
      : *max_load_it;
  metrics.mean_weighted_load = metrics.weighted_load_by_rank.empty()
      ? 0.0
      : std::accumulate(
            metrics.weighted_load_by_rank.begin(),
            metrics.weighted_load_by_rank.end(),
            0.0) /
            static_cast<double>(metrics.weighted_load_by_rank.size());
  metrics.weighted_imbalance_ratio = metrics.mean_weighted_load > 0.0
      ? metrics.max_weighted_load / metrics.mean_weighted_load
      : 0.0;
  metrics.total_memory_bytes = std::accumulate(
      metrics.memory_bytes_by_rank.begin(), metrics.memory_bytes_by_rank.end(), 0ULL);
  const auto max_mem_it = std::max_element(
      metrics.memory_bytes_by_rank.begin(), metrics.memory_bytes_by_rank.end());
  metrics.max_memory_bytes = max_mem_it == metrics.memory_bytes_by_rank.end()
      ? 0ULL
      : *max_mem_it;
  const double mean_memory = metrics.memory_bytes_by_rank.empty()
      ? 0.0
      : static_cast<double>(metrics.total_memory_bytes) /
            static_cast<double>(metrics.memory_bytes_by_rank.size());
  metrics.memory_imbalance_ratio = mean_memory > 0.0
      ? static_cast<double>(metrics.max_memory_bytes) / mean_memory
      : 0.0;
  std::uint64_t total_peak_memory = 0U;
  for (std::size_t rank = 0; rank < metrics.memory_bytes_by_rank.size(); ++rank) {
    const bool rank_has_work =
        metrics.weighted_load_by_rank[rank] > 0.0 ||
        metrics.memory_bytes_by_rank[rank] != 0U;
    metrics.peak_memory_bytes_by_rank[rank] = rank_has_work
        ? checkedUint64Add(
              metrics.memory_bytes_by_rank[rank],
              config.rank_transient_reserve_bytes,
              "compact decomposition rank peak memory")
        : 0U;
    total_peak_memory = checkedUint64Add(
        total_peak_memory,
        metrics.peak_memory_bytes_by_rank[rank],
        "compact decomposition total peak memory");
  }
  const auto max_peak_it = std::max_element(
      metrics.peak_memory_bytes_by_rank.begin(),
      metrics.peak_memory_bytes_by_rank.end());
  metrics.max_peak_memory_bytes = max_peak_it == metrics.peak_memory_bytes_by_rank.end()
      ? 0ULL
      : *max_peak_it;
  const double mean_peak_memory = metrics.peak_memory_bytes_by_rank.empty()
      ? 0.0
      : static_cast<double>(total_peak_memory) /
            static_cast<double>(metrics.peak_memory_bytes_by_rank.size());
  metrics.peak_memory_imbalance_ratio = mean_peak_memory > 0.0
      ? static_cast<double>(metrics.max_peak_memory_bytes) / mean_peak_memory
      : 0.0;
}

void accumulateCompactSourceMetrics(
    LoadBalanceMetrics& metrics,
    std::size_t rank,
    const CompactRuntimeDecompositionRecord& record,
    const DecompositionWorkComponents& components) {
  metrics.weighted_load_by_rank[rank] += record.weighted_load;
  metrics.memory_bytes_by_rank[rank] += record.memory_bytes;
  if (record.kind == DecompositionEntityKind::kParticle) {
    ++metrics.owned_particles_by_rank[rank];
  }
  metrics.active_targets_by_rank[rank] += record.active_target_count_recent;
  metrics.remote_tree_interactions_by_rank[rank] +=
      record.remote_tree_interactions_recent;
  addWorkComponentsToMetrics(metrics, rank, components, 1.0);
}

RuntimeRebalancePlan buildCompactDistributedRuntimeRebalancePlan(
    const MpiContext& mpi_context,
    const RuntimeDecompositionSourceView& source,
    const DecompositionConfig& decomposition_config,
    const RuntimeRebalanceConfig& rebalance_config,
    const DecompositionRuntimeMeasurements& measurements,
    const DecompositionFeedbackCoefficients& feedback_coefficients) {
  std::exception_ptr entry_failure;
  try {
    if (rebalance_config.world_size <= 0 ||
        decomposition_config.world_size != rebalance_config.world_size ||
        mpi_context.worldSize() != rebalance_config.world_size) {
      throw std::invalid_argument("distributed runtime rebalance world sizes must agree");
    }
    if (rebalance_config.imbalance_trigger_ratio < 1.0 ||
        rebalance_config.memory_trigger_ratio < 1.0 ||
        rebalance_config.max_migrated_load_fraction < 0.0 ||
        rebalance_config.max_migrated_load_fraction > 1.0) {
      throw std::invalid_argument("runtime rebalance thresholds are invalid");
    }
    if (mpi_context.worldSize() > 1 && !mpi_context.isEnabled()) {
      throw std::runtime_error("distributed runtime rebalance requires MPI when world_size > 1");
    }
  } catch (...) {
    entry_failure = std::current_exception();
  }
  mpi_context.rethrowCollectivePreparationFailure(
      entry_failure, "compact source rebalance entry validation");

  struct CompactCutSample {
    std::uint64_t key = 0;
    std::uint64_t entity_id = 0;
    double represented_load = 0.0;
  };
  static_assert(std::is_trivially_copyable_v<CompactCutSample>);

  RuntimeRebalancePlan rebalance;
  rebalance.used_compact_planner = true;
  rebalance.local_entities_considered = static_cast<std::uint64_t>(source.localEntityCount());
  rebalance.planner_local_entity_count = rebalance.local_entities_considered;

  std::vector<CompactRuntimeDecompositionRecord> records;
  std::vector<CompactCutSample> local_samples;
  std::vector<CompactCutSample> global_samples;
  std::vector<std::uint8_t> recv_bytes;
  std::vector<std::uint64_t> memory_block_prefix;
  std::vector<SfcCutPoint> cuts;
  std::vector<std::uint64_t> proposed_memory;
  std::vector<std::uint64_t> proposed_high;
  std::vector<SfcCutPoint> mandatory_cuts;
  std::vector<SfcCutPoint> repaired_cuts;
  RuntimeSourceComponentScratch component_scratch;
  std::exception_ptr local_preparation_failure;
  try {
     validateRuntimeDecompositionSource(source);
     if (source.world_rank >= decomposition_config.world_size) {
       throw std::invalid_argument("runtime decomposition source world_rank is outside world size");
     }
     if (source.particle_metadata_uniform || !source.particle_owning_rank.empty()) {
      for (std::size_t particle = 0; particle < source.particle_count; ++particle) {
        const std::uint32_t owner = source.particleOwningRank(particle);
        if (static_cast<std::uint64_t>(owner) >=
            static_cast<std::uint64_t>(decomposition_config.world_size)) {
          throw std::invalid_argument("runtime decomposition source particle owner is outside world size");
        }
      }
    }
    for (const std::uint32_t owner : source.patch_owning_rank) {
      if (static_cast<std::uint64_t>(owner) >=
          static_cast<std::uint64_t>(decomposition_config.world_size)) {
        throw std::invalid_argument("runtime decomposition source patch owner is outside world size");
      }
    }
    records = makeSourceRecords(
        source, decomposition_config, measurements, feedback_coefficients);
    cuts.reserve(static_cast<std::size_t>(std::max(0, decomposition_config.world_size - 1)));
    rebalance.sfc_cut_keys.reserve(cuts.capacity());
    rebalance.sfc_cut_entity_ids.reserve(cuts.capacity());
    std::sort(records.begin(), records.end(),
              [](const CompactRuntimeDecompositionRecord& lhs,
                 const CompactRuntimeDecompositionRecord& rhs) {
                if (lhs.sfc_key != rhs.sfc_key) {
                  return lhs.sfc_key < rhs.sfc_key;
                }
                if (lhs.entity_id != rhs.entity_id) {
                  return lhs.entity_id < rhs.entity_id;
                }
                return lhs.local_index < rhs.local_index;
              });
    constexpr std::size_t k_samples_per_rank = 256U;
    if (!records.empty()) {
      const std::size_t sample_count = std::min(k_samples_per_rank, records.size());
      local_samples.reserve(sample_count);
      for (std::size_t sample = 0; sample < sample_count; ++sample) {
        const std::size_t begin = sample * records.size() / sample_count;
        const std::size_t end = (sample + 1U) * records.size() / sample_count;
        double bucket_load = 0.0;
        for (std::size_t pos = begin; pos < end; ++pos) {
          bucket_load += records[pos].weighted_load;
        }
        const CompactRuntimeDecompositionRecord& boundary = records[end - 1U];
        local_samples.push_back(CompactCutSample{
            .key = boundary.sfc_key,
            .entity_id = boundary.entity_id,
            .represented_load = bucket_load,
        });
      }
    }
  } catch (...) {
    local_preparation_failure = std::current_exception();
  }
  mpi_context.rethrowCollectivePreparationFailure(
      local_preparation_failure, "compact source rebalance local preparation");

  const auto set_planner_telemetry = [&] {
    const std::uint64_t record_bytes = static_cast<std::uint64_t>(
        records.capacity() * sizeof(CompactRuntimeDecompositionRecord));
    const std::uint64_t rank_count_u64 = checkedUint64Add(
        static_cast<std::uint64_t>(mpi_context.worldSize()), 0U,
        "compact source rank count");
    const std::size_t round_limit = std::max<std::size_t>(
        1U, mpiTransportRoundLimitBytes());
    const std::size_t sample_record_bytes = sizeof(CompactCutSample);
    const std::size_t sample_payload_bytes = core::checkedSizeMultiply(
        core::checkedSizeMultiply(
            static_cast<std::size_t>(rank_count_u64), 256U,
            "compact source bounded all-gather sample count"),
        sample_record_bytes,
        "compact source bounded all-gather sample bytes");
    const std::uint64_t round_payload_bytes = static_cast<std::uint64_t>(
        sample_payload_bytes);
    const std::size_t allgather_round_count = std::max<std::size_t>(
        1U,
        (sample_payload_bytes + round_limit - 1U) / round_limit);
    const std::uint64_t sample_bytes = checkedUint64Add(
        checkedUint64Add(
            checkedUint64Add(
                static_cast<std::uint64_t>(
                    local_samples.capacity() * sizeof(CompactCutSample)),
                static_cast<std::uint64_t>(
                    global_samples.capacity() * sizeof(CompactCutSample)),
                "compact source sample byte overflow"),
            static_cast<std::uint64_t>(recv_bytes.capacity()),
            "compact source receive sample byte overflow"),
        round_payload_bytes,
        "compact source bounded all-gather round sample byte overflow");
    const std::uint64_t gather_metadata_bytes = checkedUint64Add(
        checkedUint64Multiply(
            rank_count_u64,
            sizeof(std::uint64_t) * 4U +
                sizeof(int) * 2U + sizeof(std::size_t),
            "compact source bounded all-gather metadata bytes"),
        checkedUint64Multiply(
            static_cast<std::uint64_t>(allgather_round_count),
            checkedUint64Multiply(
                rank_count_u64,
                sizeof(int) * 2U + sizeof(std::size_t),
                "compact source bounded all-gather round metadata bytes"),
            "compact source bounded all-gather round count bytes"),
        "compact source bounded all-gather metadata byte overflow");
    const std::uint64_t candidate_bytes = checkedUint64Multiply(
        rank_count_u64, sizeof(std::uint64_t) * 3U,
        "compact source repair candidate bytes");
    const std::uint64_t prefix_bytes = static_cast<std::uint64_t>(
        memory_block_prefix.capacity() * sizeof(std::uint64_t));
    const std::uint64_t intent_bytes = static_cast<std::uint64_t>(
        rebalance.particle_migrations.capacity() * sizeof(ParticleMigrationIntent) +
        rebalance.amr_patch_ownership_updates.capacity() * sizeof(AmrPatchOwnershipUpdate));
    const std::uint64_t metric_lane_bytes = checkedUint64Add(
        checkedUint64Multiply(
            5U, sizeof(std::uint64_t),
            "compact source metric lane width"),
        checkedUint64Multiply(
            13U, sizeof(double),
            "compact source metric lane width"),
        "compact source metric lane width");
    const std::uint64_t metric_bytes = checkedUint64Multiply(
        checkedUint64Multiply(
            rank_count_u64, 2U, "compact source metric vector count"),
        metric_lane_bytes, "compact source metric bytes");
    const std::uint64_t other_bytes = checkedUint64Add(
        checkedUint64Add(
            source.source_scratch_bytes, metric_bytes,
            "compact source known scratch byte overflow"),
        checkedUint64Add(
            gather_metadata_bytes, candidate_bytes,
            "compact source known scratch byte overflow"),
        "compact source known scratch byte overflow");
    const std::uint64_t scratch_bytes = other_bytes;
    const std::uint64_t control_bytes = checkedUint64Add(
        checkedUint64Add(
            scratch_bytes,
            static_cast<std::uint64_t>(cuts.capacity() * sizeof(SfcCutPoint)),
            "compact source cut scratch byte overflow"),
        static_cast<std::uint64_t>(
            (rebalance.sfc_cut_keys.capacity() +
             rebalance.sfc_cut_entity_ids.capacity()) * sizeof(std::uint64_t)),
        "compact source cut output byte overflow");
    const std::uint64_t repair_bytes = checkedUint64Add(
        control_bytes,
        static_cast<std::uint64_t>(
            (proposed_memory.capacity() + proposed_high.capacity()) * sizeof(std::uint64_t) +
            (mandatory_cuts.capacity() + repaired_cuts.capacity()) * sizeof(SfcCutPoint)),
        "compact source repair scratch byte overflow");
    rebalance.planner_record_bytes = record_bytes;
    rebalance.planner_sample_bytes = sample_bytes;
    rebalance.planner_prefix_bytes = prefix_bytes;
    rebalance.planner_migration_intent_bytes = intent_bytes;
    rebalance.planner_other_known_scratch_bytes = repair_bytes;
    rebalance.planner_peak_temporary_bytes = record_bytes;
    rebalance.planner_peak_temporary_bytes = checkedUint64Add(
        rebalance.planner_peak_temporary_bytes, sample_bytes,
        "compact source planner peak byte overflow");
    rebalance.planner_peak_temporary_bytes = checkedUint64Add(
        rebalance.planner_peak_temporary_bytes, prefix_bytes,
        "compact source planner peak byte overflow");
    rebalance.planner_peak_temporary_bytes = checkedUint64Add(
        rebalance.planner_peak_temporary_bytes, intent_bytes,
        "compact source planner peak byte overflow");
    rebalance.planner_peak_temporary_bytes = checkedUint64Add(
        rebalance.planner_peak_temporary_bytes, repair_bytes,
        "compact source planner peak byte overflow");
    rebalance.planner_local_peak_temporary_bytes =
        rebalance.planner_peak_temporary_bytes;
    rebalance.peak_temporary_bytes = rebalance.planner_peak_temporary_bytes;
    rebalance.planner_bytes_per_entity = source.localEntityCount() == 0U
        ? 0.0
        : static_cast<double>(rebalance.planner_peak_temporary_bytes) /
              static_cast<double>(source.localEntityCount());
  };

  RuntimeFeedbackNormalization normalization;
  std::exception_ptr normalization_failure;
  try {
    normalization = computeSourceFeedbackNormalization(
        source, measurements, feedback_coefficients);
  } catch (...) {
    normalization_failure = std::current_exception();
  }
  mpi_context.rethrowCollectivePreparationFailure(
      normalization_failure, "compact source feedback normalization");
  const std::size_t rank_count = static_cast<std::size_t>(mpi_context.worldSize());
  const auto make_metrics = [&]() {
    LoadBalanceMetrics metrics;
    initializeCompactMetrics(metrics, rank_count);
    return metrics;
  };

  if (mpi_context.worldSize() == 1) {
    rebalance.global_entities_considered = rebalance.local_entities_considered;
    rebalance.used_distributed_sfc_cuts = false;
    rebalance.current_metrics = make_metrics();
    rebalance.target_decomposition.metrics = make_metrics();
    double serial_amr_migrated_load = 0.0;
    for (const CompactRuntimeDecompositionRecord& record : records) {
      const DecompositionWorkComponents components = sourceComponentsForRecord(
          source, record, measurements, feedback_coefficients, normalization,
          component_scratch);
      accumulateCompactSourceMetrics(
          rebalance.current_metrics,
          static_cast<std::size_t>(record.current_owner_rank),
          record,
          components);
      accumulateCompactSourceMetrics(
          rebalance.target_decomposition.metrics, 0U, record, components);
      if (record.current_owner_rank != 0 &&
          record.kind == DecompositionEntityKind::kParticle &&
          rebalance_config.allow_particle_migration) {
        rebalance.particle_migrations.push_back(ParticleMigrationIntent{
            .particle_id = record.entity_id,
            .item_index = record.local_index,
            .old_owner_rank = record.current_owner_rank,
            .new_owner_rank = 0,
            .work_units = record.weighted_load,
        });
      } else if (record.current_owner_rank != 0 &&
                 record.kind == DecompositionEntityKind::kAmrPatch &&
                 rebalance_config.allow_amr_patch_reassignment) {
        rebalance.amr_patch_ownership_updates.push_back(AmrPatchOwnershipUpdate{
            .patch_id = record.entity_id,
            .old_owner_rank = record.current_owner_rank,
            .new_owner_rank = 0,
        });
         serial_amr_migrated_load += record.weighted_load;

      }
    }
    finalizeCompactMetrics(rebalance.current_metrics, decomposition_config);
    finalizeCompactMetrics(rebalance.target_decomposition.metrics, decomposition_config);
    rebalance.local_entities_moved = static_cast<std::uint64_t>(
        std::count_if(
            records.begin(), records.end(),
            [](const CompactRuntimeDecompositionRecord& record) {
              return record.current_owner_rank != 0;
            }));
    rebalance.global_entities_moved = rebalance.local_entities_moved;
    rebalance.local_bytes_moved = std::accumulate(
        records.begin(), records.end(), std::uint64_t{0},
        [](std::uint64_t sum, const CompactRuntimeDecompositionRecord& record) {
          return checkedUint64Add(
              sum,
              record.current_owner_rank == 0 ? 0U : record.memory_bytes,
              "compact serial moved bytes");
        });
    rebalance.global_bytes_moved = rebalance.local_bytes_moved;
    rebalance.migrated_load = std::accumulate(
        rebalance.particle_migrations.begin(), rebalance.particle_migrations.end(), 0.0,
        [](double sum, const ParticleMigrationIntent& intent) {
          return sum + intent.work_units;
        }) + serial_amr_migrated_load;
    const double total_load = std::accumulate(
        records.begin(), records.end(), 0.0,
        [](double sum, const CompactRuntimeDecompositionRecord& record) {
          return sum + record.weighted_load;
        });
    rebalance.migrated_load_fraction = total_load > 0.0
        ? rebalance.migrated_load / total_load
        : 0.0;
    const bool load_imbalanced =
        rebalance.current_metrics.weighted_imbalance_ratio >=
        rebalance_config.imbalance_trigger_ratio;
    const bool memory_imbalanced =
        rebalance.current_metrics.memory_imbalance_ratio >=
        rebalance_config.memory_trigger_ratio;
    const bool hard_memory_violated =
        decomposition_config.max_rank_memory_bytes != 0U &&
        rebalance.current_metrics.max_peak_memory_bytes >
            decomposition_config.max_rank_memory_bytes;
    if (records.empty()) {
      rebalance.reason = "empty_decomposition";
    } else if (!load_imbalanced && !memory_imbalanced && !hard_memory_violated) {
      rebalance.reason = "below_rebalance_threshold";
      rebalance.particle_migrations.clear();
      rebalance.amr_patch_ownership_updates.clear();
      rebalance.migrated_load = 0.0;
      rebalance.migrated_load_fraction = 0.0;
    } else if (!hard_memory_violated &&
               rebalance.migrated_load_fraction >
                   rebalance_config.max_migrated_load_fraction &&
               rebalance_config.max_migrated_load_fraction < 1.0) {
      rebalance.reason = "migration_fraction_limited";
      rebalance.particle_migrations.clear();
      rebalance.amr_patch_ownership_updates.clear();
      rebalance.migrated_load = 0.0;
      rebalance.migrated_load_fraction = 0.0;
    } else {
      rebalance.should_rebalance = true;
      rebalance.reason = hard_memory_violated
          ? "rank_memory_limit"
          : (load_imbalanced && memory_imbalanced
                 ? "load_and_memory_imbalance"
                 : (load_imbalanced ? "load_imbalance" : "memory_imbalance"));
    }
    set_planner_telemetry();
    return rebalance;
  }

#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  {
    std::exception_ptr sample_failure;
    try {
      const std::size_t local_sample_bytes = core::checkedSizeMultiply(
          local_samples.size(), sizeof(CompactCutSample),
          "compact source local sample byte count");
      const auto local_sample_wire = std::span<const std::uint8_t>(
          reinterpret_cast<const std::uint8_t*>(local_samples.data()),
          local_sample_bytes);
      recv_bytes = mpi_context.allgatherBytesBounded(local_sample_wire);
      if (recv_bytes.size() % sizeof(CompactCutSample) != 0U) {
        throw std::runtime_error(
            "compact source cut sample exchange returned partial record bytes");
      }
      global_samples.resize(recv_bytes.size() / sizeof(CompactCutSample));
      if (!recv_bytes.empty()) {
        std::memcpy(global_samples.data(), recv_bytes.data(), recv_bytes.size());
      }
    } catch (...) {
      sample_failure = std::current_exception();
    }
    mpi_context.rethrowCollectivePreparationFailure(
        sample_failure, "compact source cut-sample exchange");
  }

  std::sort(global_samples.begin(), global_samples.end(),
            [](const CompactCutSample& lhs, const CompactCutSample& rhs) {
              return lessSfcPoint(
                  SfcCutPoint{.key = lhs.key, .entity_id = lhs.entity_id},
                  SfcCutPoint{.key = rhs.key, .entity_id = rhs.entity_id});
            });

  const double global_total_load = mpi_context.allreduceSumDouble(
      std::accumulate(
          records.begin(), records.end(), 0.0,
          [](double sum, const CompactRuntimeDecompositionRecord& record) {
            return sum + record.weighted_load;
          }));
  const std::uint64_t global_entity_count =
      mpi_context.allreduceSumUint64(rebalance.local_entities_considered);
  rebalance.global_entities_considered = global_entity_count;
  if (!global_samples.empty() && global_total_load > 0.0) {
    const double target_per_rank =
        global_total_load / static_cast<double>(mpi_context.worldSize());
    double cumulative_sample_load = 0.0;
    std::size_t next_cut_rank = 1U;
    for (const CompactCutSample& sample : global_samples) {
      cumulative_sample_load += std::max(0.0, sample.represented_load);
      if (next_cut_rank < static_cast<std::size_t>(mpi_context.worldSize()) &&
          cumulative_sample_load >= target_per_rank * static_cast<double>(next_cut_rank)) {
        cuts.push_back(SfcCutPoint{.key = sample.key, .entity_id = sample.entity_id});
        ++next_cut_rank;
      }
    }
  }
  while (cuts.size() + 1U < static_cast<std::size_t>(mpi_context.worldSize())) {
    cuts.push_back(global_samples.empty()
                       ? SfcCutPoint{}
                       : SfcCutPoint{
                             global_samples.back().key,
                             global_samples.back().entity_id});
  }

  if (decomposition_config.max_rank_memory_bytes != 0U) {
    // Declared outside the try block so the post-collective repair predicate
    // can capture it; it is only consumed when the collective preparation
    // phase above succeeded on every rank.
    std::uint64_t persistent_limit = 0U;
    std::exception_ptr memory_failure;
    try {
      if (decomposition_config.rank_transient_reserve_bytes >=
          decomposition_config.max_rank_memory_bytes) {
        throw std::invalid_argument(
            "distributed decomposition transient reserve must be smaller than the hard rank memory ceiling");
      }
      persistent_limit =
          decomposition_config.max_rank_memory_bytes -
          decomposition_config.rank_transient_reserve_bytes;
      proposed_memory.assign(rank_count, 0U);
      proposed_high.assign(rank_count, 0U);
      for (const CompactRuntimeDecompositionRecord& entry : records) {
        const SfcCutPoint point{.key = entry.sfc_key, .entity_id = entry.entity_id};
        const std::size_t owner = static_cast<std::size_t>(
            ownerForSfcPoint(point, cuts, mpi_context.worldSize()));
        proposed_memory[owner] = checkedUint64Add(
            proposed_memory[owner], entry.memory_bytes,
            "compact source proposed rank memory");
      }
      for (std::size_t rank = 0; rank < rank_count; ++rank) {
        proposed_high[rank] = proposed_memory[rank] >> 32U;
        proposed_memory[rank] &= 0xffffffffULL;
      }
    } catch (...) {
      memory_failure = std::current_exception();
    }
    mpi_context.rethrowCollectivePreparationFailure(
        memory_failure, "compact source memory preflight");
    if (MPI_Allreduce(
            MPI_IN_PLACE, proposed_memory.data(), mpi_context.worldSize(),
            MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD) != MPI_SUCCESS ||
        MPI_Allreduce(
            MPI_IN_PLACE, proposed_high.data(), mpi_context.worldSize(),
            MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD) != MPI_SUCCESS) {
      throw std::runtime_error("compact source memory preflight Allreduce failed");
    }
     for (std::size_t rank = 0; rank < rank_count; ++rank) {
       const std::uint64_t high = checkedUint64Add(
           proposed_high[rank], proposed_memory[rank] >> 32U,
           "compact source memory carry");
       if (high > 0xffffffffULL) {
         throw std::overflow_error("compact source proposed rank memory exceeds uint64");
       }
       proposed_memory[rank] = (high << 32U) |
           (proposed_memory[rank] & 0xffffffffULL);
     }

    const bool needs_memory_repair = std::any_of(
        proposed_memory.begin(), proposed_memory.end(),
        [persistent_limit](std::uint64_t bytes) {
          return bytes > persistent_limit;
        });
    if (needs_memory_repair) {
      constexpr std::size_t k_memory_prefix_block_size = 256U;
      std::exception_ptr prefix_failure;
      std::uint64_t local_total_memory = 0U;
      try {
        const std::size_t block_count = records.size() / k_memory_prefix_block_size +
            (records.size() % k_memory_prefix_block_size != 0U ? 1U : 0U);
        memory_block_prefix.reserve(core::checkedSizeAdd(
            block_count, 1U, "compact source memory prefix block count"));
        for (std::size_t pos = 0; pos < records.size(); ++pos) {
          if (pos % k_memory_prefix_block_size == 0U) {
            memory_block_prefix.push_back(local_total_memory);
          }
          local_total_memory = checkedUint64Add(
              local_total_memory, records[pos].memory_bytes,
              "compact source local persistent memory");
        }
        memory_block_prefix.push_back(local_total_memory);
      } catch (...) {
        prefix_failure = std::current_exception();
      }
      mpi_context.rethrowCollectivePreparationFailure(
          prefix_failure, "compact source memory prefix preparation");
      const auto exact_global_sum = [&](std::uint64_t local_bytes) {
        const std::uint64_t low =
            mpi_context.allreduceSumUint64(local_bytes & 0xffffffffULL);
        const std::uint64_t high =
            mpi_context.allreduceSumUint64(local_bytes >> 32U);
        const std::uint64_t carry = checkedUint64Add(
            high, low >> 32U, "compact source memory sum carry");
        if (carry > 0xffffffffULL) {
          throw std::overflow_error("compact source total memory exceeds uint64");
        }
        return (carry << 32U) | (low & 0xffffffffULL);
      };
      const auto local_prefix = [&](SfcCutPoint point) {
        const auto it = std::upper_bound(
            records.begin(), records.end(), point,
            [](SfcCutPoint value, const CompactRuntimeDecompositionRecord& entry) {
              return lessSfcPoint(
                  value, SfcCutPoint{entry.sfc_key, entry.entity_id});
            });
        const std::size_t end = static_cast<std::size_t>(
            std::distance(records.begin(), it));
        const std::size_t block = end / k_memory_prefix_block_size;
        std::uint64_t bytes = memory_block_prefix[block];
        for (std::size_t pos = block * k_memory_prefix_block_size; pos < end; ++pos) {
          bytes = checkedUint64Add(
              bytes, records[pos].memory_bytes,
              "compact source local prefix sum");
        }
        return bytes;
      };
      const auto global_prefix = [&](SfcCutPoint point) {
        return exact_global_sum(local_prefix(point));
      };
      const auto first_prefix_at_least = [&](std::uint64_t target) {
        std::uint64_t key_low = 0U;
        std::uint64_t key_high = (1ULL << 30U) - 1U;
        while (key_low < key_high) {
          const std::uint64_t mid = key_low + (key_high - key_low) / 2U;
          if (global_prefix(SfcCutPoint{mid, std::numeric_limits<std::uint64_t>::max()}) >= target) {
            key_high = mid;
          } else {
            key_low = mid + 1U;
          }
        }
        std::uint64_t id_low = 0U;
        std::uint64_t id_high = std::numeric_limits<std::uint64_t>::max();
        while (id_low < id_high) {
          const std::uint64_t mid = id_low + (id_high - id_low) / 2U;
          if (global_prefix(SfcCutPoint{key_low, mid}) >= target) {
            id_high = mid;
          } else {
            id_low = mid + 1U;
          }
        }
        return SfcCutPoint{key_low, id_low};
      };
      const auto global_neighbor = [&](SfcCutPoint point, bool predecessor) {
        struct Candidate {
          std::uint64_t valid;
          std::uint64_t key;
          std::uint64_t entity_id;
        };
        Candidate local{};
        const auto it = std::lower_bound(
            records.begin(), records.end(), point,
            [](const CompactRuntimeDecompositionRecord& entry, SfcCutPoint value) {
              return lessSfcPoint(
                  SfcCutPoint{entry.sfc_key, entry.entity_id}, value);
            });
        if (predecessor ? it != records.begin() : it != records.end()) {
          const auto& entry = predecessor ? *std::prev(it) : *it;
          local = Candidate{1U, entry.sfc_key, entry.entity_id};
        }
        const auto wire = std::span<const std::uint8_t>(
            reinterpret_cast<const std::uint8_t*>(&local), sizeof(local));
        const std::vector<std::uint8_t> received =
            mpi_context.allgatherBytesBounded(wire);
        if (received.size() != core::checkedSizeMultiply(
                rank_count, sizeof(Candidate),
                "compact source SFC neighbor extent")) {
          throw std::runtime_error("compact source SFC neighbor exchange has invalid size");
        }
        std::optional<SfcCutPoint> result;
        for (std::size_t rank = 0; rank < rank_count; ++rank) {
          Candidate candidate{};
          std::memcpy(
              &candidate, received.data() + rank * sizeof(Candidate),
              sizeof(Candidate));
          if (candidate.valid == 0U) {
            continue;
          }
          const SfcCutPoint value{candidate.key, candidate.entity_id};
          if (!result.has_value() ||
              (predecessor ? lessSfcPoint(*result, value) : lessSfcPoint(value, *result))) {
            result = value;
          }
        }
        return result;
      };
      std::exception_ptr repair_failure;
      try {
        mandatory_cuts.resize(rank_count - 1U);
        repaired_cuts.reserve(rank_count - 1U);
      } catch (...) {
        repair_failure = std::current_exception();
      }
      mpi_context.rethrowCollectivePreparationFailure(
          repair_failure, "compact source memory repair metadata preparation");
      const auto first_candidate = global_neighbor(SfcCutPoint{}, false);
      if (!first_candidate.has_value()) {
        throw std::logic_error("compact source memory repair requires a nonempty population");
      }
      const std::uint64_t total_memory = exact_global_sum(local_total_memory);
      const SfcCutPoint last_point = global_samples.empty()
          ? *first_candidate
          : SfcCutPoint{global_samples.back().key, global_samples.back().entity_id};
      cuts = internal::repairMemoryConstrainedSfcCuts(
          std::span<const SfcCutPoint>(cuts),
          decomposition_config.max_rank_memory_bytes -
              decomposition_config.rank_transient_reserve_bytes,
          total_memory,
          *first_candidate,
          last_point,
          global_prefix,
          first_prefix_at_least,
          global_neighbor,
          std::span<SfcCutPoint>(mandatory_cuts),
          repaired_cuts);
    }
  }

  LoadBalanceMetrics local_current = make_metrics();
  LoadBalanceMetrics local_target = make_metrics();
  std::uint64_t local_moved_entities = 0U;
  std::uint64_t local_moved_bytes = 0U;
  double local_migrated_load = 0.0;
  std::size_t particle_intent_count = 0U;
  std::size_t patch_intent_count = 0U;
  for (const CompactRuntimeDecompositionRecord& entry : records) {
    const int new_owner = ownerForSfcPoint(
        SfcCutPoint{.key = entry.sfc_key, .entity_id = entry.entity_id},
        cuts,
        mpi_context.worldSize());
    if (entry.current_owner_rank == new_owner) {
      continue;
    }
    if (entry.kind == DecompositionEntityKind::kParticle &&
        rebalance_config.allow_particle_migration) {
      ++particle_intent_count;
    } else if (entry.kind == DecompositionEntityKind::kAmrPatch &&
               rebalance_config.allow_amr_patch_reassignment) {
      ++patch_intent_count;
    }
  }
  std::exception_ptr intent_preparation_failure;
  try {
    rebalance.particle_migrations.reserve(particle_intent_count);
    rebalance.amr_patch_ownership_updates.reserve(patch_intent_count);
  } catch (...) {
    intent_preparation_failure = std::current_exception();
  }
  mpi_context.rethrowCollectivePreparationFailure(
      intent_preparation_failure, "compact source migration intent preparation");

  std::exception_ptr metric_preparation_failure;
  try {
    for (const CompactRuntimeDecompositionRecord& entry : records) {
      const DecompositionWorkComponents components = sourceComponentsForRecord(
          source, entry, measurements, feedback_coefficients, normalization,
          component_scratch);
      const int new_owner = ownerForSfcPoint(
          SfcCutPoint{.key = entry.sfc_key, .entity_id = entry.entity_id},
          cuts,
          mpi_context.worldSize());
      accumulateCompactSourceMetrics(
          local_current,
          static_cast<std::size_t>(entry.current_owner_rank),
          entry,
          components);
      accumulateCompactSourceMetrics(
          local_target, static_cast<std::size_t>(new_owner), entry, components);
      if (entry.current_owner_rank != new_owner) {
        ++local_moved_entities;
        local_moved_bytes = checkedUint64Add(
            local_moved_bytes, entry.memory_bytes,
            "compact source moved bytes");
        local_migrated_load += entry.weighted_load;
        if (entry.kind == DecompositionEntityKind::kParticle &&
            rebalance_config.allow_particle_migration) {
          rebalance.particle_migrations.push_back(ParticleMigrationIntent{
              .particle_id = entry.entity_id,
              .item_index = entry.local_index,
              .old_owner_rank = entry.current_owner_rank,
              .new_owner_rank = new_owner,
              .work_units = entry.weighted_load,
          });
        } else if (entry.kind == DecompositionEntityKind::kAmrPatch &&
                   rebalance_config.allow_amr_patch_reassignment) {
          rebalance.amr_patch_ownership_updates.push_back(AmrPatchOwnershipUpdate{
              .patch_id = entry.entity_id,
              .old_owner_rank = entry.current_owner_rank,
              .new_owner_rank = new_owner,
          });
        }
      }
    }
  } catch (...) {
    metric_preparation_failure = std::current_exception();
  }
  mpi_context.rethrowCollectivePreparationFailure(
      metric_preparation_failure, "compact source metric and intent reconstruction");

  const auto allreduce_double_vector = [&](std::vector<double>& values) {
    if (MPI_Allreduce(
            MPI_IN_PLACE, values.data(), mpi_context.worldSize(),
            MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD) != MPI_SUCCESS) {
      throw std::runtime_error("compact source double metric Allreduce failed");
    }
  };
  const auto allreduce_uint64_vector = [&](std::vector<std::uint64_t>& values) {
    if (MPI_Allreduce(
            MPI_IN_PLACE, values.data(), mpi_context.worldSize(),
            MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD) != MPI_SUCCESS) {
      throw std::runtime_error("compact source uint64 metric Allreduce failed");
    }
  };
  const auto reduce_metrics = [&](LoadBalanceMetrics& metrics) {
    allreduce_double_vector(metrics.weighted_load_by_rank);
    allreduce_uint64_vector(metrics.memory_bytes_by_rank);
    allreduce_uint64_vector(metrics.owned_particles_by_rank);
    allreduce_uint64_vector(metrics.active_targets_by_rank);
    allreduce_uint64_vector(metrics.remote_tree_interactions_by_rank);
    allreduce_double_vector(metrics.particle_count_cost_by_rank);
    allreduce_double_vector(metrics.gas_cell_cost_by_rank);
    allreduce_double_vector(metrics.tree_interaction_cost_by_rank);
    allreduce_double_vector(metrics.pm_mesh_cost_by_rank);
    allreduce_double_vector(metrics.amr_patch_cost_by_rank);
    allreduce_double_vector(metrics.active_fraction_cost_by_rank);
    allreduce_double_vector(metrics.memory_pressure_cost_by_rank);
    allreduce_double_vector(metrics.transient_memory_cost_by_rank);
    allreduce_double_vector(metrics.source_event_cost_by_rank);
    allreduce_double_vector(metrics.communication_cost_by_rank);
    allreduce_double_vector(metrics.gpu_occupancy_cost_by_rank);
    allreduce_double_vector(metrics.generic_work_cost_by_rank);
    finalizeCompactMetrics(metrics, decomposition_config);
  };
  reduce_metrics(local_current);
  reduce_metrics(local_target);
  rebalance.current_metrics = std::move(local_current);
  rebalance.target_decomposition.metrics = std::move(local_target);
  rebalance.local_entities_moved = local_moved_entities;
  rebalance.global_entities_moved = mpi_context.allreduceSumUint64(local_moved_entities);
  rebalance.local_bytes_moved = local_moved_bytes;
  rebalance.global_bytes_moved = mpi_context.allreduceSumUint64(local_moved_bytes);
  rebalance.migrated_load = mpi_context.allreduceSumDouble(local_migrated_load);
  rebalance.migrated_load_fraction = global_total_load > 0.0
      ? rebalance.migrated_load / global_total_load
      : 0.0;
  rebalance.local_control_bytes = static_cast<std::uint64_t>(
      local_samples.capacity() * sizeof(CompactCutSample));
  rebalance.global_control_bytes = mpi_context.allreduceSumUint64(
      rebalance.local_control_bytes);
  for (const SfcCutPoint& cut : cuts) {
    rebalance.sfc_cut_keys.push_back(cut.key);
    rebalance.sfc_cut_entity_ids.push_back(cut.entity_id);
  }
  std::exception_ptr telemetry_failure;
  try {
    set_planner_telemetry();
  } catch (...) {
    telemetry_failure = std::current_exception();
  }
  mpi_context.rethrowCollectivePreparationFailure(
      telemetry_failure, "compact source planner telemetry");
  rebalance.used_distributed_sfc_cuts = true;

  const bool load_imbalanced =
      rebalance.current_metrics.weighted_imbalance_ratio >=
      rebalance_config.imbalance_trigger_ratio;
  const bool memory_imbalanced =
      rebalance.current_metrics.memory_imbalance_ratio >=
      rebalance_config.memory_trigger_ratio;
  const bool hard_memory_violated =
      decomposition_config.max_rank_memory_bytes != 0U &&
      rebalance.current_metrics.max_peak_memory_bytes >
          decomposition_config.max_rank_memory_bytes;
  const bool target_memory_safe =
      decomposition_config.max_rank_memory_bytes == 0U ||
      rebalance.target_decomposition.metrics.max_peak_memory_bytes <=
          decomposition_config.max_rank_memory_bytes;
  if (!target_memory_safe) {
    if (hard_memory_violated) {
      throw std::runtime_error(
          "distributed source rebalance cannot find a hard-memory-safe target");
    }
    rebalance.particle_migrations.clear();
    rebalance.amr_patch_ownership_updates.clear();
    rebalance.should_rebalance = false;
    rebalance.reason = "target_rank_memory_limit_exceeded";
    set_planner_telemetry();
    return rebalance;
  }
  const bool local_has_actionable_migration =
      !rebalance.particle_migrations.empty() ||
      !rebalance.amr_patch_ownership_updates.empty();
  const std::uint64_t actionable_migration_rank_count =
      mpi_context.allreduceSumUint64(local_has_actionable_migration ? 1ULL : 0ULL);
  if (global_entity_count == 0U) {
    rebalance.reason = "empty_decomposition";
  } else if (!load_imbalanced && !memory_imbalanced && !hard_memory_violated) {
    rebalance.reason = "below_rebalance_threshold";
    rebalance.particle_migrations.clear();
    rebalance.amr_patch_ownership_updates.clear();
    rebalance.migrated_load = 0.0;
    rebalance.migrated_load_fraction = 0.0;
  } else if (!hard_memory_violated &&
             rebalance.migrated_load_fraction >
                 rebalance_config.max_migrated_load_fraction &&
             rebalance_config.max_migrated_load_fraction < 1.0) {
    rebalance.reason = "migration_fraction_limited";
  } else {
    rebalance.should_rebalance =
        rebalance.global_entities_moved > 0U &&
        actionable_migration_rank_count > 0U;
    if (hard_memory_violated && !rebalance.should_rebalance) {
      throw std::runtime_error(
          "distributed source rebalance cannot relieve the hard rank memory violation");
    }
    rebalance.reason = hard_memory_violated
        ? "rank_memory_limit"
        : (load_imbalanced && memory_imbalanced
               ? "load_and_memory_imbalance"
               : (load_imbalanced ? "load_imbalance" : "memory_imbalance"));
  }
  set_planner_telemetry();
  return rebalance;
#else
  throw std::runtime_error("distributed compact source rebalance requires an MPI-enabled build");
#endif
}


}  // namespace cosmosim::parallel
