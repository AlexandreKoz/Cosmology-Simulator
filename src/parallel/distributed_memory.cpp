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

#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
[[nodiscard]] bool queryActiveMpiWorld(int& world_size, int& world_rank) noexcept {
  world_size = 1;
  world_rank = 0;
  int initialized = 0;
  MPI_Initialized(&initialized);
  if (initialized == 0) {
    return false;
  }
  int finalized = 0;
  MPI_Finalized(&finalized);
  if (finalized != 0) {
    return false;
  }
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  return true;
}

[[nodiscard]] int queryNodeLocalRank(int fallback_world_rank) {
  MPI_Comm local_comm = MPI_COMM_NULL;
  const int split_result = MPI_Comm_split_type(
      MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &local_comm);
  if (split_result != MPI_SUCCESS || local_comm == MPI_COMM_NULL) {
    throw std::runtime_error("MPI_Comm_split_type(MPI_COMM_TYPE_SHARED) failed while determining node-local rank");
  }
  int local_rank = fallback_world_rank;
  const int rank_result = MPI_Comm_rank(local_comm, &local_rank);
  const int free_result = MPI_Comm_free(&local_comm);
  if (rank_result != MPI_SUCCESS || free_result != MPI_SUCCESS) {
    throw std::runtime_error("MPI node-local communicator query failed");
  }
  return local_rank;
}
#endif

void injectMpiTestFault(const MpiContext& mpi_context, std::string_view phase) {
#if COSMOSIM_ENABLE_TESTS
  const char* raw = std::getenv("COSMOSIM_MPI_TEST_FAULT");
  if (raw == nullptr || *raw == '\0') {
    return;
  }
  const std::string specification(raw);
  const std::size_t separator = specification.rfind(':');
  if (separator == std::string::npos) {
    return;
  }
  int configured_rank = -1;
  try {
    configured_rank = std::stoi(specification.substr(separator + 1U));
  } catch (...) {
    return;
  }
  if (configured_rank == mpi_context.worldRank() &&
      specification.substr(0U, separator) == phase) {
    throw std::runtime_error(
        "test-only injected MPI preparation failure at phase " + std::string(phase));
  }
#else
  static_cast<void>(mpi_context);
  static_cast<void>(phase);
#endif
}

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

[[nodiscard]] std::uint64_t sfcKeyForItem(const DecompositionItem& item, const DecompositionConfig& config) {
  const std::uint32_t qx = quantize10bit(item.x_comov, config.domain_x_min_comov, config.domain_x_max_comov);
  const std::uint32_t qy = quantize10bit(item.y_comov, config.domain_y_min_comov, config.domain_y_max_comov);
  const std::uint32_t qz = quantize10bit(item.z_comov, config.domain_z_min_comov, config.domain_z_max_comov);
  return mortonKey3d(qx, qy, qz);
}

struct SfcCutPoint {
  std::uint64_t key = 0;
  std::uint64_t entity_id = 0;
};

[[nodiscard]] bool lessSfcPoint(const SfcCutPoint& lhs, const SfcCutPoint& rhs) noexcept {
  if (lhs.key != rhs.key) {
    return lhs.key < rhs.key;
  }
  return lhs.entity_id < rhs.entity_id;
}

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
  (void)measured;
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

constexpr std::uint16_t kGhostLanePositionX = 1U << 0U;
constexpr std::uint16_t kGhostLanePositionY = 1U << 1U;
constexpr std::uint16_t kGhostLanePositionZ = 1U << 2U;
constexpr std::uint16_t kGhostLaneMass = 1U << 3U;
constexpr std::uint16_t kGhostLaneDensity = 1U << 4U;
constexpr std::uint16_t kGhostLaneVelocityX = 1U << 5U;
constexpr std::uint16_t kGhostLaneVelocityY = 1U << 6U;
constexpr std::uint16_t kGhostLaneVelocityZ = 1U << 7U;
constexpr std::uint16_t kGhostLanePressure = 1U << 8U;
constexpr std::uint16_t kGhostLaneInternalEnergy = 1U << 9U;
constexpr std::uint16_t kGhostKnownLaneMask = (1U << 10U) - 1U;

[[nodiscard]] std::uint16_t ghostOptionalLaneMask(const GhostExchangeBufferSoA& source) noexcept {
  std::uint16_t mask = 0U;
  const auto add_if_present = [&mask](const std::vector<double>& lane, std::uint16_t bit) {
    if (!lane.empty()) {
      mask = static_cast<std::uint16_t>(mask | bit);
    }
  };
  add_if_present(source.position_x_comoving, kGhostLanePositionX);
  add_if_present(source.position_y_comoving, kGhostLanePositionY);
  add_if_present(source.position_z_comoving, kGhostLanePositionZ);
  add_if_present(source.mass_code, kGhostLaneMass);
  add_if_present(source.density_code, kGhostLaneDensity);
  add_if_present(source.velocity_x_code, kGhostLaneVelocityX);
  add_if_present(source.velocity_y_code, kGhostLaneVelocityY);
  add_if_present(source.velocity_z_code, kGhostLaneVelocityZ);
  add_if_present(source.pressure_code, kGhostLanePressure);
  add_if_present(source.internal_energy_code, kGhostLaneInternalEnergy);
  return mask;
}

[[nodiscard]] bool laneIsPresentOrEmpty(std::size_t size, std::size_t expected) noexcept {
  return size == 0 || size == expected;
}

[[nodiscard]] double optionalLaneValue(const std::vector<double>& lane, std::size_t index) {
  return lane.empty() ? 0.0 : lane[index];
}

void commitOptionalGhostLane(
    std::vector<double>* destination,
    const std::vector<double>& source,
    std::size_t destination_size,
    std::size_t destination_index,
    std::size_t source_index) {
  if (source.empty()) {
    return;
  }
  if (source_index >= source.size()) {
    throw std::invalid_argument("ghost source optional lane does not cover committed row");
  }
  if (destination->empty()) {
    destination->assign(destination_size, 0.0);
  } else if (destination->size() != destination_size) {
    throw std::invalid_argument("ghost optional lane size does not match storage size");
  }
  (*destination)[destination_index] = source[source_index];
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
    keyed[i] = KeyedItem{
        .index = i,
        .morton_key = sfcKeyForItem(items[i], config),
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

RuntimeRebalancePlan buildDistributedRuntimeRebalancePlan(
    const MpiContext& mpi_context,
    std::span<const DecompositionItem> local_items,
    const DecompositionConfig& decomposition_config,
    const RuntimeRebalanceConfig& rebalance_config) {
  if (rebalance_config.world_size <= 0 || decomposition_config.world_size != rebalance_config.world_size ||
      mpi_context.worldSize() != rebalance_config.world_size) {
    throw std::invalid_argument("distributed runtime rebalance world sizes must agree");
  }
  if (rebalance_config.imbalance_trigger_ratio < 1.0 || rebalance_config.memory_trigger_ratio < 1.0 ||
      rebalance_config.max_migrated_load_fraction < 0.0 || rebalance_config.max_migrated_load_fraction > 1.0) {
    throw std::invalid_argument("distributed runtime rebalance thresholds are invalid");
  }
  if (mpi_context.worldSize() == 1) {
    RuntimeRebalancePlan serial = buildRuntimeRebalancePlan(local_items, decomposition_config, rebalance_config);
    serial.used_distributed_sfc_cuts = false;
    serial.local_entities_considered = static_cast<std::uint64_t>(local_items.size());
    serial.global_entities_considered = static_cast<std::uint64_t>(local_items.size());
    return serial;
  }
  if (!mpi_context.isEnabled()) {
    throw std::runtime_error("distributed runtime rebalance requires MPI when world_size > 1");
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

  std::vector<LocalKeyedItem> keyed(local_items.size());
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
  std::vector<CompactCutSample> local_samples;
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

  std::vector<CompactCutSample> global_samples;
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  std::size_t local_sample_bytes = 0U;
  std::exception_ptr sample_preparation_failure;
  try {
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
  [[maybe_unused]] auto finalize_metrics = [](LoadBalanceMetrics& metrics) {
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
      global_samples.capacity() * sizeof(CompactCutSample));
  rebalance.cut_displacement_fraction =
      (global_entity_count > 0U) ? static_cast<double>(rebalance.global_entities_moved) / static_cast<double>(global_entity_count) : 0.0;

  const bool load_imbalanced =
      rebalance.current_metrics.weighted_imbalance_ratio >= rebalance_config.imbalance_trigger_ratio;
  const bool memory_imbalanced =
      rebalance.current_metrics.memory_imbalance_ratio >= rebalance_config.memory_trigger_ratio;
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
  } else if (!load_imbalanced && !memory_imbalanced) {
    rebalance.reason = "below_rebalance_threshold";
  } else if (rebalance.migrated_load_fraction > rebalance_config.max_migrated_load_fraction &&
             rebalance_config.max_migrated_load_fraction < 1.0) {
    rebalance.reason = "migration_fraction_limited";
  } else {
    rebalance.reason = load_imbalanced && memory_imbalanced ? "load_and_memory_imbalance" :
        (load_imbalanced ? "load_imbalance" : "memory_imbalance");
    rebalance.should_rebalance =
        rebalance.global_entities_moved > 0U &&
        actionable_migration_rank_count > 0U;
  }
  return rebalance;
}

BoundedMpiTransferPlan planBoundedMpiTransferRounds(
    std::span<const std::size_t> logical_counts,
    std::size_t mpi_count_limit,
    std::size_t round_count_limit) {
  if (mpi_count_limit == 0U ||
      mpi_count_limit > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    throw std::invalid_argument(
        "bounded MPI planner count limit must be within classic MPI int range");
  }
  if (round_count_limit == 0U || round_count_limit > mpi_count_limit) {
    throw std::invalid_argument(
        "bounded MPI planner round limit must be positive and no larger than the MPI count limit");
  }

  BoundedMpiTransferPlan plan;
  plan.logical_counts.assign(logical_counts.begin(), logical_counts.end());
  plan.logical_displacements.resize(logical_counts.size(), 0U);
  if (logical_counts.empty()) {
    return plan;
  }
  if (logical_counts.size() > round_count_limit) {
    throw std::invalid_argument(
        "bounded MPI planner round limit is too small to address every participant safely");
  }

  std::size_t logical_total = 0U;
  for (std::size_t peer = 0; peer < logical_counts.size(); ++peer) {
    plan.logical_displacements[peer] = logical_total;
    logical_total = core::checkedSizeAdd(
        logical_total, logical_counts[peer], "bounded MPI logical prefix");
  }
  plan.logical_total_count = logical_total;

  const std::size_t fair_round_share = round_count_limit / logical_counts.size();
  plan.per_peer_count_limit = std::min(mpi_count_limit, fair_round_share);
  if (plan.per_peer_count_limit == 0U) {
    throw std::invalid_argument("bounded MPI planner computed a zero per-peer round limit");
  }

  std::size_t round_count = 0U;
  for (const std::size_t count : logical_counts) {
    const std::size_t peer_rounds = count / plan.per_peer_count_limit +
        (count % plan.per_peer_count_limit == 0U ? 0U : 1U);
    round_count = std::max(round_count, peer_rounds);
  }
  plan.rounds.reserve(round_count);
  std::vector<std::size_t> consumed(logical_counts.size(), 0U);
  for (std::size_t round_index = 0; round_index < round_count; ++round_index) {
    BoundedMpiRoundLayout round;
    round.counts.resize(logical_counts.size(), 0);
    round.displacements.resize(logical_counts.size(), 0);
    round.logical_offsets = consumed;

    std::size_t round_total = 0U;
    for (std::size_t peer = 0; peer < logical_counts.size(); ++peer) {
      const std::size_t remaining = logical_counts[peer] - consumed[peer];
      const std::size_t count = std::min(plan.per_peer_count_limit, remaining);
      round.displacements[peer] = core::checkedIntegralNarrow<int>(
          round_total, "bounded MPI round displacement");
      round.counts[peer] = core::checkedIntegralNarrow<int>(
          count, "bounded MPI round count");
      round_total = core::checkedSizeAdd(
          round_total, count, "bounded MPI round aggregate");
      if (round_total > round_count_limit || round_total > mpi_count_limit) {
        throw std::overflow_error("bounded MPI round aggregate exceeds representability limit");
      }
      consumed[peer] = core::checkedSizeAdd(
          consumed[peer], count, "bounded MPI logical coverage");
    }
    round.round_count = round_total;
    plan.rounds.push_back(std::move(round));
  }
  if (consumed != plan.logical_counts) {
    throw std::logic_error("bounded MPI planner did not cover the complete logical payload");
  }
  return plan;
}

std::size_t mpiTransportRoundLimitBytes() {
#if COSMOSIM_ENABLE_TESTS
  const char* raw = std::getenv("COSMOSIM_MPI_TEST_TRANSPORT_LIMIT_BYTES");
  if (raw != nullptr && *raw != '\0') {
    try {
      const unsigned long long parsed = std::stoull(raw);
      if (parsed > 0ULL && parsed <= static_cast<unsigned long long>(std::numeric_limits<int>::max())) {
        return static_cast<std::size_t>(parsed);
      }
    } catch (...) {
      // Invalid test-only overrides fall back to the production-safe default.
    }
  }
#endif
  return k_default_mpi_transport_round_bytes;
}

DirectedAmrPatchCellTransferPlan planDirectedAmrPatchCellTransfer(
    std::uint64_t logical_record_count,
    std::size_t transport_round_limit_bytes) {
  const std::size_t requested_limit = transport_round_limit_bytes == 0U
      ? mpiTransportRoundLimitBytes()
      : transport_round_limit_bytes;
  const std::size_t bounded_limit = std::min(
      requested_limit,
      static_cast<std::size_t>(std::numeric_limits<int>::max()));
  const std::size_t records_per_round = bounded_limit / sizeof(AmrPatchCellPayloadRecord);
  if (records_per_round == 0U) {
    throw std::invalid_argument(
        "directed AMR patch-cell transport round cannot hold one complete record");
  }
  const std::uint64_t round_count_u64 = logical_record_count / records_per_round +
      (logical_record_count % records_per_round == 0U ? 0U : 1U);
  return DirectedAmrPatchCellTransferPlan{
      .logical_record_count = logical_record_count,
      .records_per_round = records_per_round,
      .round_count = round_count_u64,
  };
}

std::vector<std::size_t> planDirectedAmrPatchCellTransferRounds(
    std::uint64_t logical_record_count,
    std::size_t transport_round_limit_bytes) {
  constexpr std::uint64_t k_max_materialized_round_count = 65'536U;
  const DirectedAmrPatchCellTransferPlan plan = planDirectedAmrPatchCellTransfer(
      logical_record_count, transport_round_limit_bytes);
  if (plan.round_count > k_max_materialized_round_count) {
    throw std::length_error(
        "directed AMR diagnostic round materialization exceeds bounded metadata cap; "
        "use the constant-space transfer plan");
  }
  const std::size_t round_count = core::checkedIntegralNarrow<std::size_t>(
      plan.round_count, "directed AMR patch-cell materialized round count");
  std::vector<std::size_t> rounds;
  rounds.reserve(round_count);
  std::uint64_t remaining = plan.logical_record_count;
  while (remaining != 0U) {
    const std::size_t records = static_cast<std::size_t>(std::min<std::uint64_t>(
        remaining, static_cast<std::uint64_t>(plan.records_per_round)));
    rounds.push_back(records);
    remaining -= records;
  }
  return rounds;
}

std::vector<std::vector<std::uint8_t>> exchangeBoundedAlltoallBytes(
    const MpiContext& mpi_context,
    const std::vector<std::vector<std::uint8_t>>& send_payloads) {
  const int world_size = mpi_context.worldSize();
  if (world_size <= 0) {
    throw std::invalid_argument("bounded byte all-to-all requires positive world size");
  }
  if (!mpi_context.isEnabled()) {
    if (world_size != 1 || send_payloads.size() != 1U) {
      throw std::runtime_error(
          "bounded byte all-to-all requires MPI when world size exceeds one");
    }
    return {send_payloads.front()};
  }
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  std::vector<std::uint64_t> send_counts64;
  std::vector<std::uint64_t> recv_counts64;
  std::exception_ptr local_failure;
  try {
    if (send_payloads.size() != static_cast<std::size_t>(world_size)) {
      throw std::invalid_argument(
          "bounded byte all-to-all payload rank extent must match MPI world size");
    }
    send_counts64.resize(static_cast<std::size_t>(world_size), 0U);
    recv_counts64.resize(static_cast<std::size_t>(world_size), 0U);
    for (std::size_t rank = 0; rank < send_payloads.size(); ++rank) {
      send_counts64[rank] = core::checkedIntegralNarrow<std::uint64_t>(
          send_payloads[rank].size(), "bounded byte all-to-all send count");
    }
  } catch (...) {
    local_failure = std::current_exception();
  }
  mpi_context.rethrowCollectivePreparationFailure(
      local_failure, "bounded byte all-to-all pre-count preparation");

  if (MPI_Alltoall(
          send_counts64.data(), 1, MPI_UINT64_T,
          recv_counts64.data(), 1, MPI_UINT64_T,
          MPI_COMM_WORLD) != MPI_SUCCESS) {
    throw std::runtime_error("bounded byte all-to-all count exchange failed");
  }

  BoundedMpiTransferPlan send_plan;
  BoundedMpiTransferPlan recv_plan;
  std::vector<std::vector<std::uint8_t>> recv_payloads;
  std::vector<std::uint8_t> send_round_buffer;
  std::vector<std::uint8_t> recv_round_buffer;
  BoundedMpiRoundLayout zero_round;
  local_failure = nullptr;
  try {
    std::vector<std::size_t> send_counts(send_counts64.size(), 0U);
    std::vector<std::size_t> recv_counts(recv_counts64.size(), 0U);
    for (std::size_t rank = 0; rank < send_counts.size(); ++rank) {
      send_counts[rank] = core::checkedIntegralNarrow<std::size_t>(
          send_counts64[rank], "bounded byte all-to-all logical send count");
      recv_counts[rank] = core::checkedIntegralNarrow<std::size_t>(
          recv_counts64[rank], "bounded byte all-to-all logical receive count");
    }
    const std::size_t round_limit = mpiTransportRoundLimitBytes();
    send_plan = planBoundedMpiTransferRounds(
        send_counts,
        static_cast<std::size_t>(std::numeric_limits<int>::max()),
        round_limit);
    recv_plan = planBoundedMpiTransferRounds(
        recv_counts,
        static_cast<std::size_t>(std::numeric_limits<int>::max()),
        round_limit);

    recv_payloads.resize(static_cast<std::size_t>(world_size));
    for (std::size_t rank = 0; rank < recv_payloads.size(); ++rank) {
      recv_payloads[rank].resize(recv_counts[rank]);
    }
    std::size_t maximum_send_round = 0U;
    for (const auto& round : send_plan.rounds) {
      maximum_send_round = std::max(maximum_send_round, round.round_count);
    }
    std::size_t maximum_recv_round = 0U;
    for (const auto& round : recv_plan.rounds) {
      maximum_recv_round = std::max(maximum_recv_round, round.round_count);
    }
    send_round_buffer.resize(maximum_send_round);
    recv_round_buffer.resize(maximum_recv_round);
    zero_round.counts.resize(static_cast<std::size_t>(world_size), 0);
    zero_round.displacements.resize(static_cast<std::size_t>(world_size), 0);
    zero_round.logical_offsets.resize(static_cast<std::size_t>(world_size), 0U);
    injectMpiTestFault(mpi_context, "alltoall_post_count");
  } catch (...) {
    local_failure = std::current_exception();
  }
  mpi_context.rethrowCollectivePreparationFailure(
      local_failure, "bounded byte all-to-all payload preparation");

  const std::uint64_t local_round_count = static_cast<std::uint64_t>(
      std::max(send_plan.rounds.size(), recv_plan.rounds.size()));
  const std::uint64_t global_round_count =
      mpi_context.allreduceMaxUint64(local_round_count);
  for (std::uint64_t round_index = 0U; round_index < global_round_count;
       ++round_index) {
    const BoundedMpiRoundLayout& send_round =
        round_index < send_plan.rounds.size()
            ? send_plan.rounds[static_cast<std::size_t>(round_index)]
            : zero_round;
    const BoundedMpiRoundLayout& recv_round =
        round_index < recv_plan.rounds.size()
            ? recv_plan.rounds[static_cast<std::size_t>(round_index)]
            : zero_round;

    for (std::size_t peer = 0; peer < send_round.counts.size(); ++peer) {
      const std::size_t count = static_cast<std::size_t>(send_round.counts[peer]);
      if (count == 0U) {
        continue;
      }
      std::memcpy(
          send_round_buffer.data() +
              static_cast<std::size_t>(send_round.displacements[peer]),
          send_payloads[peer].data() + send_round.logical_offsets[peer],
          count);
    }
    if (MPI_Alltoallv(
            send_round_buffer.empty() ? nullptr : send_round_buffer.data(),
            send_round.counts.data(), send_round.displacements.data(), MPI_BYTE,
            recv_round_buffer.empty() ? nullptr : recv_round_buffer.data(),
            recv_round.counts.data(), recv_round.displacements.data(), MPI_BYTE,
            MPI_COMM_WORLD) != MPI_SUCCESS) {
      throw std::runtime_error("bounded byte all-to-all payload exchange failed");
    }
    for (std::size_t peer = 0; peer < recv_round.counts.size(); ++peer) {
      const std::size_t count = static_cast<std::size_t>(recv_round.counts[peer]);
      if (count == 0U) {
        continue;
      }
      std::memcpy(
          recv_payloads[peer].data() + recv_round.logical_offsets[peer],
          recv_round_buffer.data() +
              static_cast<std::size_t>(recv_round.displacements[peer]),
          count);
    }
  }
  return recv_payloads;
#else
  throw std::runtime_error(
      "bounded byte all-to-all requires an MPI-enabled build");
#endif
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
  std::size_t local_byte_count = 0U;
  std::exception_ptr local_preparation_failure;
  try {
    local_byte_count = core::checkedSizeMultiply(
        local_items.size(), sizeof(DecompositionItem),
        "decomposition item gather local byte count");
  } catch (...) {
    local_preparation_failure = std::current_exception();
  }
  mpi_context.rethrowCollectivePreparationFailure(
      local_preparation_failure,
      "decomposition item gather local byte framing");

  const auto local_bytes = std::span<const std::uint8_t>(
      reinterpret_cast<const std::uint8_t*>(local_items.data()),
      local_byte_count);
  const std::vector<std::uint8_t> recv_bytes =
      mpi_context.allgatherBytesBounded(local_bytes);

  std::vector<DecompositionItem> gathered;
  std::exception_ptr local_reassembly_failure;
  try {
    if (recv_bytes.size() % sizeof(DecompositionItem) != 0U) {
      throw std::runtime_error(
          "decomposition item gather returned a non-record-aligned byte count");
    }
    gathered.resize(recv_bytes.size() / sizeof(DecompositionItem));
    if (!recv_bytes.empty()) {
      std::memcpy(gathered.data(), recv_bytes.data(), recv_bytes.size());
    }
  } catch (...) {
    local_reassembly_failure = std::current_exception();
  }
  mpi_context.rethrowCollectivePreparationFailure(
      local_reassembly_failure,
      "decomposition item gather receive reassembly");
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
    std::span<const std::uint64_t> expected_local_reference_particle_ids) {
  ExactOwnershipPartitionReport report;
  report.global_owned_count = mpi_context.allreduceSumUint64(
      static_cast<std::uint64_t>(local_owned_particle_ids.size()));

  const auto sorted_unique = [](std::span<const std::uint64_t> values) {
    std::vector<std::uint64_t> sorted(values.begin(), values.end());
    std::sort(sorted.begin(), sorted.end());
    const bool unique =
        std::adjacent_find(sorted.begin(), sorted.end()) == sorted.end();
    return std::pair{std::move(sorted), unique};
  };

  if (!mpi_context.isEnabled()) {
    auto [current, current_unique] = sorted_unique(local_owned_particle_ids);
    auto [expected, expected_unique] =
        sorted_unique(expected_local_reference_particle_ids);
    report.local_particle_ids_unique = current_unique;
    report.globally_unique = current_unique;
    for (auto it = current.begin(); it != current.end();) {
      const auto range = std::equal_range(it, current.end(), *it);
      if (std::distance(range.first, range.second) > 1) {
        report.duplicate_particle_ids.push_back(*it);
      }
      it = range.second;
    }
    current.erase(std::unique(current.begin(), current.end()), current.end());
    expected.erase(std::unique(expected.begin(), expected.end()), expected.end());
    std::set_difference(
        expected.begin(), expected.end(), current.begin(), current.end(),
        std::back_inserter(report.missing_expected_particle_ids));
    std::set_difference(
        current.begin(), current.end(), expected.begin(), expected.end(),
        std::back_inserter(report.extra_particle_ids));
    report.matches_expected_ids = expected_unique &&
        report.missing_expected_particle_ids.empty() &&
        report.extra_particle_ids.empty();
    return report;
  }

#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  const int world_size = mpi_context.worldSize();
  const auto exchange_by_hash = [&](std::span<const std::uint64_t> ids) {
    std::vector<std::vector<std::uint8_t>> send_payloads;
    std::exception_ptr local_preparation_failure;
    try {
      send_payloads.resize(static_cast<std::size_t>(world_size));
      for (const std::uint64_t id : ids) {
        const std::uint64_t mixed = id ^ (id >> 33U) ^ (id << 11U);
        const std::size_t owner = static_cast<std::size_t>(
            mixed % static_cast<std::uint64_t>(world_size));
        auto& payload = send_payloads[owner];
        const std::size_t old_size = payload.size();
        const std::size_t new_size = core::checkedSizeAdd(
            old_size, sizeof(id), "exact ownership hash bucket byte growth");
        payload.resize(new_size);
        std::memcpy(payload.data() + old_size, &id, sizeof(id));
      }
    } catch (...) {
      local_preparation_failure = std::current_exception();
    }
    mpi_context.rethrowCollectivePreparationFailure(
        local_preparation_failure,
        "exact ownership hash local payload preparation");

    std::vector<std::vector<std::uint8_t>> recv_payloads =
        exchangeBoundedAlltoallBytes(mpi_context, send_payloads);

    std::vector<std::uint64_t> recv;
    std::exception_ptr local_decode_failure;
    try {
      std::size_t total_ids = 0U;
      for (const auto& payload : recv_payloads) {
        if (payload.size() % sizeof(std::uint64_t) != 0U) {
          throw std::runtime_error(
              "exact ownership hash exchange returned partial uint64 record bytes");
        }
        total_ids = core::checkedSizeAdd(
            total_ids, payload.size() / sizeof(std::uint64_t),
            "exact ownership hash receive record total");
      }
      recv.resize(total_ids);
      std::size_t destination = 0U;
      for (const auto& payload : recv_payloads) {
        const std::size_t count = payload.size() / sizeof(std::uint64_t);
        if (!payload.empty()) {
          std::memcpy(
              recv.data() + destination, payload.data(), payload.size());
        }
        destination += count;
      }
    } catch (...) {
      local_decode_failure = std::current_exception();
    }
    mpi_context.rethrowCollectivePreparationFailure(
        local_decode_failure,
        "exact ownership hash receive reassembly");
    return recv;
  };

  std::vector<std::uint64_t> current_partition =
      exchange_by_hash(local_owned_particle_ids);
  std::vector<std::uint64_t> expected_partition =
      exchange_by_hash(expected_local_reference_particle_ids);
  std::sort(current_partition.begin(), current_partition.end());
  std::sort(expected_partition.begin(), expected_partition.end());

  for (auto it = current_partition.begin(); it != current_partition.end();) {
    const auto range = std::equal_range(it, current_partition.end(), *it);
    if (std::distance(range.first, range.second) > 1) {
      report.duplicate_particle_ids.push_back(*it);
    }
    it = range.second;
  }
  const bool local_unique = report.duplicate_particle_ids.empty();
  current_partition.erase(
      std::unique(current_partition.begin(), current_partition.end()),
      current_partition.end());
  expected_partition.erase(
      std::unique(expected_partition.begin(), expected_partition.end()),
      expected_partition.end());
  std::set_difference(
      expected_partition.begin(), expected_partition.end(),
      current_partition.begin(), current_partition.end(),
      std::back_inserter(report.missing_expected_particle_ids));
  std::set_difference(
      current_partition.begin(), current_partition.end(),
      expected_partition.begin(), expected_partition.end(),
      std::back_inserter(report.extra_particle_ids));

  const int local_duplicate = local_unique ? 0 : 1;
  const int local_mismatch =
      (report.missing_expected_particle_ids.empty() &&
       report.extra_particle_ids.empty()) ? 0 : 1;
  int any_duplicate = 0;
  int any_mismatch = 0;
  MPI_Allreduce(
      &local_duplicate, &any_duplicate, 1, MPI_INT, MPI_MAX,
      MPI_COMM_WORLD);
  MPI_Allreduce(
      &local_mismatch, &any_mismatch, 1, MPI_INT, MPI_MAX,
      MPI_COMM_WORLD);
  report.local_particle_ids_unique = any_duplicate == 0;
  report.globally_unique = any_duplicate == 0;
  report.matches_expected_ids = any_mismatch == 0;
  return report;
#else
  throw std::runtime_error(
      "exact global ownership validation requires MPI support when MPI context is enabled");
#endif
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
  return n != 0U && density_code.size() == n && velocity_x_code.size() == n && velocity_y_code.size() == n &&
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
  appendPod<std::uint16_t>(m_bytes, ghostOptionalLaneMask(source));

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
  constexpr std::size_t header_bytes = sizeof(std::uint64_t) + sizeof(std::uint16_t);
  if (m_bytes.size() < header_bytes) {
    throw std::runtime_error("ghost buffer is too small");
  }

  std::size_t offset = 0;
  const std::uint64_t count = readPod<std::uint64_t>(m_bytes, &offset);
  const std::uint16_t lane_mask = readPod<std::uint16_t>(m_bytes, &offset);
  if ((lane_mask & static_cast<std::uint16_t>(~kGhostKnownLaneMask)) != 0U) {
    throw std::runtime_error("ghost buffer contains an unknown optional-lane presence bit");
  }
  const std::uint64_t expected_payload_bytes =
      static_cast<std::uint64_t>(header_bytes) + count * static_cast<std::uint64_t>(ghostExchangeRecordBytes());
  if (expected_payload_bytes != static_cast<std::uint64_t>(m_bytes.size())) {
    throw std::runtime_error("ghost buffer payload shape does not match encoded count");
  }

  const std::size_t append_count = static_cast<std::size_t>(count);
  const std::size_t existing_count = destination.entity_id.size();
  const auto require_lane_schema = [existing_count, lane_mask](
                                       const std::vector<double>& lane,
                                       std::uint16_t bit) {
    if (existing_count == 0U) {
      return;
    }
    const bool existing_present = lane.size() == existing_count;
    const bool incoming_present = (lane_mask & bit) != 0U;
    if (existing_present != incoming_present) {
      throw std::runtime_error("ghost peer payloads disagree on optional-lane presence");
    }
  };
  require_lane_schema(destination.position_x_comoving, kGhostLanePositionX);
  require_lane_schema(destination.position_y_comoving, kGhostLanePositionY);
  require_lane_schema(destination.position_z_comoving, kGhostLanePositionZ);
  require_lane_schema(destination.mass_code, kGhostLaneMass);
  require_lane_schema(destination.density_code, kGhostLaneDensity);
  require_lane_schema(destination.velocity_x_code, kGhostLaneVelocityX);
  require_lane_schema(destination.velocity_y_code, kGhostLaneVelocityY);
  require_lane_schema(destination.velocity_z_code, kGhostLaneVelocityZ);
  require_lane_schema(destination.pressure_code, kGhostLanePressure);
  require_lane_schema(destination.internal_energy_code, kGhostLaneInternalEnergy);

  destination.entity_id.reserve(destination.entity_id.size() + append_count);
  const auto reserve_if_present = [append_count, lane_mask](
                                      std::vector<double>* lane,
                                      std::uint16_t bit) {
    if ((lane_mask & bit) != 0U) {
      lane->reserve(lane->size() + append_count);
    }
  };
  reserve_if_present(&destination.position_x_comoving, kGhostLanePositionX);
  reserve_if_present(&destination.position_y_comoving, kGhostLanePositionY);
  reserve_if_present(&destination.position_z_comoving, kGhostLanePositionZ);
  reserve_if_present(&destination.mass_code, kGhostLaneMass);
  reserve_if_present(&destination.density_code, kGhostLaneDensity);
  reserve_if_present(&destination.velocity_x_code, kGhostLaneVelocityX);
  reserve_if_present(&destination.velocity_y_code, kGhostLaneVelocityY);
  reserve_if_present(&destination.velocity_z_code, kGhostLaneVelocityZ);
  reserve_if_present(&destination.pressure_code, kGhostLanePressure);
  reserve_if_present(&destination.internal_energy_code, kGhostLaneInternalEnergy);

  for (std::uint64_t i = 0; i < count; ++i) {
    destination.entity_id.push_back(readPod<std::uint64_t>(m_bytes, &offset));
    const auto append_lane = [&](std::vector<double>* lane, std::uint16_t bit) {
      const double value = readPod<double>(m_bytes, &offset);
      if ((lane_mask & bit) != 0U) {
        lane->push_back(value);
      }
    };
    append_lane(&destination.position_x_comoving, kGhostLanePositionX);
    append_lane(&destination.position_y_comoving, kGhostLanePositionY);
    append_lane(&destination.position_z_comoving, kGhostLanePositionZ);
    append_lane(&destination.mass_code, kGhostLaneMass);
    append_lane(&destination.density_code, kGhostLaneDensity);
    append_lane(&destination.velocity_x_code, kGhostLaneVelocityX);
    append_lane(&destination.velocity_y_code, kGhostLaneVelocityY);
    append_lane(&destination.velocity_z_code, kGhostLaneVelocityZ);
    append_lane(&destination.pressure_code, kGhostLanePressure);
    append_lane(&destination.internal_energy_code, kGhostLaneInternalEnergy);
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
  for (std::size_t index = 0; index < state.owning_rank_by_item.size(); ++index) {
    const int rank = state.owning_rank_by_item[index];
    if (rank < 0 || rank >= state.world_size) {
      throw std::invalid_argument(
          "restart ownership entry rank is outside world_size at item " +
          std::to_string(index) + ": rank=" + std::to_string(rank) +
          ", world_size=" + std::to_string(state.world_size));
    }
  }
  return state;
}

DistributedRestartCompatibilityReport evaluateDistributedRestartCompatibility(
    const DistributedRestartState& restart_state,
    const DistributedExecutionTopology& runtime_topology) {
  DistributedRestartCompatibilityReport report;
  if (restart_state.schema_version != 2) {
    report.supported_schema_match = false;
    report.mismatch_messages.push_back(
        "distributed restart schema mismatch: expected=2, observed=" +
        std::to_string(restart_state.schema_version));
  }
  if (restart_state.world_size != runtime_topology.world_size) {
    report.world_size_match = false;
    report.mismatch_messages.push_back(
        "world_size mismatch: restart=" + std::to_string(restart_state.world_size) +
        ", runtime=" + std::to_string(runtime_topology.world_size));
  }
  if (restart_state.pm_slab_begin_x_by_rank.size() != restart_state.pm_slab_end_x_by_rank.size() ||
      restart_state.pm_slab_begin_x_by_rank.size() != static_cast<std::size_t>(std::max(restart_state.world_size, 0))) {
    report.pm_slab_table_shape_match = false;
    report.mismatch_messages.push_back(
        "PM slab table mismatch: begin_count=" +
        std::to_string(restart_state.pm_slab_begin_x_by_rank.size()) +
        ", end_count=" + std::to_string(restart_state.pm_slab_end_x_by_rank.size()) +
        ", restart_world_size=" + std::to_string(restart_state.world_size));
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
    report.mismatch_messages.push_back(
        "runtime world_rank is outside restart PM slab ownership table: rank=" +
        std::to_string(runtime_topology.world_rank) +
        ", table_size=" + std::to_string(restart_state.pm_slab_begin_x_by_rank.size()));
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
  if (descriptor.wire_version != 1U) {
    throw std::invalid_argument("tree pseudo-particle wire version is unsupported");
  }
  if (descriptor.source_rank < 0) {
    throw std::invalid_argument("tree pseudo-particle source_rank must be non-negative");
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
  if (packet.geometry_frame > 1U) {
    throw std::invalid_argument("tree pseudo-particle packet geometry frame is invalid");
  }
}

void validateHydroGhostCellDescriptor(const HydroGhostCellDescriptor& descriptor) {
  if (descriptor.gas_cell_id == 0) {
    throw std::invalid_argument("hydro ghost cell descriptor requires a non-zero gas_cell_id");
  }
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

void validateHydroGhostCellRequest(const HydroGhostCellRequest& request) {
  validateHydroGhostCellDescriptor(request.descriptor);
  if (request.face_key == 0) {
    throw std::invalid_argument("hydro ghost cell request requires a non-zero face key");
  }
  if (request.axis > 2U || request.side > 1U) {
    throw std::invalid_argument("hydro ghost cell request carries invalid face orientation metadata");
  }
}

void validateHydroGhostCellPayloadRecord(const HydroGhostCellPayloadRecord& record) {
  validateHydroGhostCellDescriptor(record.descriptor);
  if (record.face_key == 0) {
    throw std::invalid_argument("hydro ghost cell payload requires a non-zero face key");
  }
  if (!std::isfinite(record.mass_density_comoving) ||
      !std::isfinite(record.momentum_density_x_comoving) ||
      !std::isfinite(record.momentum_density_y_comoving) ||
      !std::isfinite(record.momentum_density_z_comoving) ||
      !std::isfinite(record.total_energy_density_comoving) ||
      !std::isfinite(record.metal_mass_density_comoving)) {
    throw std::invalid_argument("hydro ghost cell payload contains non-finite conserved state");
  }
  if (record.mass_density_comoving <= 0.0) {
    throw std::invalid_argument("hydro ghost cell payload requires positive mass density");
  }
  if (record.metal_mass_density_comoving < 0.0 ||
      record.metal_mass_density_comoving > record.mass_density_comoving * (1.0 + 1.0e-12)) {
    throw std::invalid_argument("hydro ghost cell payload contains invalid metal mass density");
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
  if (record.extent_x_comoving <= 0.0 || record.extent_y_comoving <= 0.0 || record.extent_z_comoving <= 0.0 ||
      record.cell_dim_x == 0U || record.cell_dim_y == 0U || record.cell_dim_z == 0U) {
    throw std::invalid_argument("AMR patch payload requires explicit positive patch geometry");
  }
  const std::uint64_t geometry_cells =
      static_cast<std::uint64_t>(record.cell_dim_x) *
      static_cast<std::uint64_t>(record.cell_dim_y) *
      static_cast<std::uint64_t>(record.cell_dim_z);
  if (geometry_cells != record.cell_count) {
    throw std::invalid_argument("AMR patch payload cell_count does not match explicit patch geometry");
  }
  if (!std::isfinite(record.origin_x_comoving) || !std::isfinite(record.origin_y_comoving) ||
      !std::isfinite(record.origin_z_comoving) || !std::isfinite(record.extent_x_comoving) ||
      !std::isfinite(record.extent_y_comoving) || !std::isfinite(record.extent_z_comoving)) {
    throw std::invalid_argument("AMR patch payload contains non-finite patch geometry");
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
      !std::isfinite(record.velocity_x_peculiar) || !std::isfinite(record.velocity_y_peculiar) ||
      !std::isfinite(record.velocity_z_peculiar) || !std::isfinite(record.density_code) ||
      !std::isfinite(record.pressure_code) ||
      !std::isfinite(record.internal_energy_code) || !std::isfinite(record.temperature_code) ||
      !std::isfinite(record.sound_speed_code) || !std::isfinite(record.metal_mass_code)) {
    throw std::invalid_argument("AMR patch cell payload contains non-finite state");
  }
  if (record.density_code <= 0.0 || record.pressure_code <= 0.0) {
    throw std::invalid_argument("AMR patch cell payload requires positive thermodynamic state");
  }
  if (record.metal_mass_code < 0.0 ||
      record.metal_mass_code > record.mass_code + 1.0e-12 * std::max(1.0, std::abs(record.mass_code))) {
    throw std::invalid_argument("AMR patch cell payload metal mass must lie within the gas-cell mass");
  }
}

void validateAmrFluxRegisterPayloadRecord(const AmrFluxRegisterPayloadRecord& record) {
  if (record.register_key == 0U || record.coarse_patch_id == 0U || record.coarse_gas_cell_id == 0U) {
    throw std::invalid_argument("AMR flux-register payload requires non-zero stable identity fields");
  }
  if (record.source_rank < 0 || record.owner_rank < 0) {
    throw std::invalid_argument("AMR flux-register payload ranks must be non-negative");
  }
  if (record.axis > 2U || record.orientation > 1U) {
    throw std::invalid_argument("AMR flux-register payload carries invalid face identity");
  }
  if (record.face_area_comov <= 0.0 || record.coarse_area_comov <= 0.0 ||
      record.fine_area_comov <= 0.0 || record.dt_code <= 0.0) {
    throw std::invalid_argument("AMR flux-register payload requires positive area and timestep metadata");
  }
  if (record.coarse_face_count == 0U && record.fine_face_count == 0U) {
    throw std::invalid_argument("AMR flux-register payload must carry at least one face contribution");
  }
  const std::array values{
      record.coarse_mass_flux_code,
      record.coarse_momentum_x_flux_code,
      record.coarse_momentum_y_flux_code,
      record.coarse_momentum_z_flux_code,
      record.coarse_total_energy_flux_code,
      record.coarse_metal_mass_flux_code,
      record.fine_mass_flux_code,
      record.fine_momentum_x_flux_code,
      record.fine_momentum_y_flux_code,
      record.fine_momentum_z_flux_code,
      record.fine_total_energy_flux_code,
      record.fine_metal_mass_flux_code,
      record.face_area_comov,
      record.coarse_area_comov,
      record.fine_area_comov,
      record.dt_code};
  for (const double value : values) {
    if (!std::isfinite(value)) {
      throw std::invalid_argument("AMR flux-register payload contains non-finite values");
    }
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
      !std::isfinite(record.delta_total_energy_density_comoving) ||
      !std::isfinite(record.delta_metal_mass_density_comoving)) {
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
        .step_index = std::nullopt,
        .simulation_time_code = std::nullopt,
        .scale_factor = std::nullopt,
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
  m_is_enabled = queryActiveMpiWorld(m_world_size, m_world_rank);
  m_local_rank = m_is_enabled ? queryNodeLocalRank(m_world_rank) : 0;
#endif
}

MpiContext::MpiContext(bool is_enabled, int world_size, int world_rank)
    : m_is_enabled(is_enabled), m_world_size(world_size), m_world_rank(world_rank), m_local_rank(world_rank) {
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

int MpiContext::localRank() const noexcept { return m_local_rank; }

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

std::uint64_t MpiContext::allreduceMaxUint64(std::uint64_t local_value) const {
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  if (m_is_enabled) {
    std::uint64_t global = 0;
    MPI_Allreduce(&local_value, &global, 1, MPI_UINT64_T, MPI_MAX, MPI_COMM_WORLD);
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


void MpiContext::rethrowCollectivePreparationFailure(
    const std::exception_ptr& local_failure,
    std::string_view phase_name) const {
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  if (m_is_enabled) {
    const int local_failed = local_failure != nullptr ? 1 : 0;
    int failed_count = 0;
    if (MPI_Allreduce(
            &local_failed, &failed_count, 1, MPI_INT, MPI_SUM,
            MPI_COMM_WORLD) != MPI_SUCCESS) {
      throw std::runtime_error(
          "MPI readiness Allreduce failed during collective preparation phase");
    }
    if (failed_count == 0) {
      return;
    }

    const int local_candidate = local_failure != nullptr ? m_world_rank : m_world_size;
    int failure_rank = m_world_size;
    if (MPI_Allreduce(
            &local_candidate, &failure_rank, 1, MPI_INT, MPI_MIN,
            MPI_COMM_WORLD) != MPI_SUCCESS) {
      throw std::runtime_error(
          "MPI failure-rank Allreduce failed during collective preparation phase");
    }

    constexpr std::size_t k_maximum_message_bytes = 2047U;
    std::array<char, k_maximum_message_bytes + 1U> message_buffer{};
    std::uint32_t message_length = 0U;
    if (m_world_rank == failure_rank) {
      const char* message = "unknown non-standard exception";
      try {
        std::rethrow_exception(local_failure);
      } catch (const std::exception& error) {
        message = error.what();
      } catch (...) {
      }
      const std::size_t raw_length = std::char_traits<char>::length(message);
      message_length = static_cast<std::uint32_t>(
          std::min(raw_length, k_maximum_message_bytes));
      std::copy_n(message, message_length, message_buffer.data());
    }
    if (MPI_Bcast(
            &message_length, 1, MPI_UINT32_T, failure_rank,
            MPI_COMM_WORLD) != MPI_SUCCESS ||
        MPI_Bcast(
            message_buffer.data(), static_cast<int>(message_buffer.size()), MPI_CHAR,
            failure_rank, MPI_COMM_WORLD) != MPI_SUCCESS) {
      throw std::runtime_error(
          "MPI diagnostic broadcast failed during collective preparation phase");
    }
    throw std::runtime_error(
        "collective preparation phase '" + std::string(phase_name) +
        "' failed on rank " + std::to_string(failure_rank) + ": " +
        std::string(message_buffer.data(),
                    message_buffer.data() + message_length));
  }
#endif
  if (local_failure != nullptr) {
    std::rethrow_exception(local_failure);
  }
}

std::vector<std::uint8_t> MpiContext::gatherBytesToRoot(
    std::span<const std::uint8_t> local_bytes,
    int root_rank) const {
  if (root_rank < 0 || root_rank >= m_world_size) {
    throw std::invalid_argument("MpiContext::gatherBytesToRoot invalid root rank");
  }
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  if (m_is_enabled) {
    std::uint64_t local_count64 = 0U;
    std::vector<std::uint64_t> counts64;
    std::exception_ptr local_failure;
    try {
      local_count64 = core::checkedIntegralNarrow<std::uint64_t>(
          local_bytes.size(), "gatherBytesToRoot local byte count");
      counts64.resize(static_cast<std::size_t>(m_world_size), 0U);
    } catch (...) {
      local_failure = std::current_exception();
    }
    rethrowCollectivePreparationFailure(
        local_failure, "gatherBytesToRoot count-buffer preparation");

    if (MPI_Allgather(
            &local_count64, 1, MPI_UINT64_T,
            counts64.data(), 1, MPI_UINT64_T,
            MPI_COMM_WORLD) != MPI_SUCCESS) {
      throw std::runtime_error("gatherBytesToRoot count Allgather failed");
    }

    BoundedMpiTransferPlan plan;
    std::vector<std::uint8_t> gathered;
    std::vector<std::uint8_t> round_receive_buffer;
    local_failure = nullptr;
    try {
      std::vector<std::size_t> logical_counts(counts64.size(), 0U);
      for (std::size_t rank = 0; rank < counts64.size(); ++rank) {
        logical_counts[rank] = core::checkedIntegralNarrow<std::size_t>(
            counts64[rank], "gatherBytesToRoot received byte count");
      }
      plan = planBoundedMpiTransferRounds(
          logical_counts,
          static_cast<std::size_t>(std::numeric_limits<int>::max()),
          mpiTransportRoundLimitBytes());
      if (m_world_rank == root_rank) {
        gathered.resize(plan.logical_total_count);
        std::size_t maximum_round_count = 0U;
        for (const auto& round : plan.rounds) {
          maximum_round_count = std::max(maximum_round_count, round.round_count);
        }
        round_receive_buffer.resize(maximum_round_count);
      }
      injectMpiTestFault(*this, "gather_post_count");
    } catch (...) {
      local_failure = std::current_exception();
    }
    rethrowCollectivePreparationFailure(
        local_failure, "gatherBytesToRoot payload preparation");

    for (const BoundedMpiRoundLayout& round : plan.rounds) {
      const std::size_t local_rank = static_cast<std::size_t>(m_world_rank);
      const int send_count = round.counts[local_rank];
      const std::size_t local_offset = round.logical_offsets[local_rank];
      const std::uint8_t* send_pointer = send_count == 0
          ? nullptr
          : local_bytes.data() + local_offset;
      const int status = MPI_Gatherv(
          const_cast<std::uint8_t*>(send_pointer), send_count, MPI_BYTE,
          m_world_rank == root_rank && !round_receive_buffer.empty()
              ? round_receive_buffer.data()
              : nullptr,
          m_world_rank == root_rank ? round.counts.data() : nullptr,
          m_world_rank == root_rank ? round.displacements.data() : nullptr,
          MPI_BYTE, root_rank, MPI_COMM_WORLD);
      if (status != MPI_SUCCESS) {
        throw std::runtime_error("gatherBytesToRoot bounded payload Gatherv failed");
      }
      if (m_world_rank == root_rank) {
        for (std::size_t rank = 0; rank < round.counts.size(); ++rank) {
          const std::size_t count = static_cast<std::size_t>(round.counts[rank]);
          if (count == 0U) {
            continue;
          }
          // The planner checked the logical prefix and every consumed peer offset
          // before any payload round, so this addition cannot overflow here.
          const std::size_t destination_offset =
              plan.logical_displacements[rank] + round.logical_offsets[rank];
          std::memcpy(
              gathered.data() + destination_offset,
              round_receive_buffer.data() +
                  static_cast<std::size_t>(round.displacements[rank]),
              count);
        }
      }
    }
    return gathered;
  }
#endif
  if (root_rank != 0) {
    throw std::invalid_argument("serial MpiContext only supports root rank 0");
  }
  return {local_bytes.begin(), local_bytes.end()};
}

std::vector<std::uint8_t> MpiContext::broadcastBytesFromRoot(
    std::span<const std::uint8_t> root_bytes,
    int root_rank) const {
  if (root_rank < 0 || root_rank >= m_world_size) {
    throw std::invalid_argument("MpiContext::broadcastBytesFromRoot invalid root rank");
  }
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  if (m_is_enabled) {
    std::uint64_t byte_count64 = 0U;
    std::exception_ptr local_failure;
    try {
      if (m_world_rank == root_rank) {
        byte_count64 = core::checkedIntegralNarrow<std::uint64_t>(
            root_bytes.size(), "broadcastBytesFromRoot root byte count");
      }
    } catch (...) {
      local_failure = std::current_exception();
    }
    rethrowCollectivePreparationFailure(
        local_failure, "broadcastBytesFromRoot size preparation");
    if (MPI_Bcast(
            &byte_count64, 1, MPI_UINT64_T, root_rank,
            MPI_COMM_WORLD) != MPI_SUCCESS) {
      throw std::runtime_error("broadcastBytesFromRoot size Bcast failed");
    }

    std::vector<std::uint8_t> bytes;
    local_failure = nullptr;
    try {
      const std::size_t byte_count = core::checkedIntegralNarrow<std::size_t>(
          byte_count64, "broadcastBytesFromRoot receive byte count");
      bytes.resize(byte_count);
      if (m_world_rank == root_rank && !root_bytes.empty()) {
        std::copy(root_bytes.begin(), root_bytes.end(), bytes.begin());
      }
      injectMpiTestFault(*this, "broadcast_post_count");
    } catch (...) {
      local_failure = std::current_exception();
    }
    rethrowCollectivePreparationFailure(
        local_failure, "broadcastBytesFromRoot payload preparation");

    const std::size_t round_limit = std::min(
        mpiTransportRoundLimitBytes(),
        static_cast<std::size_t>(std::numeric_limits<int>::max()));
    for (std::size_t offset = 0U; offset < bytes.size();) {
      const std::size_t chunk_size = std::min(round_limit, bytes.size() - offset);
      const int chunk_count = core::checkedIntegralNarrow<int>(
          chunk_size, "broadcastBytesFromRoot bounded chunk count");
      if (MPI_Bcast(
              bytes.data() + offset, chunk_count, MPI_BYTE,
              root_rank, MPI_COMM_WORLD) != MPI_SUCCESS) {
        throw std::runtime_error("broadcastBytesFromRoot bounded payload Bcast failed");
      }
      // chunk_size is bounded by bytes.size() - offset.
      offset += chunk_size;
    }
    return bytes;
  }
#endif
  if (root_rank != 0) {
    throw std::invalid_argument("serial MpiContext only supports root rank 0");
  }
  return {root_bytes.begin(), root_bytes.end()};
}

std::vector<std::uint8_t> MpiContext::allgatherBytesBounded(
    std::span<const std::uint8_t> local_bytes) const {
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  if (m_is_enabled) {
    std::uint64_t local_count64 = 0U;
    std::vector<std::uint64_t> counts64;
    std::exception_ptr local_failure;
    try {
      local_count64 = core::checkedIntegralNarrow<std::uint64_t>(
          local_bytes.size(), "allgatherBytesBounded local byte count");
      counts64.resize(static_cast<std::size_t>(m_world_size), 0U);
    } catch (...) {
      local_failure = std::current_exception();
    }
    rethrowCollectivePreparationFailure(
        local_failure, "allgatherBytesBounded count-buffer preparation");
    if (MPI_Allgather(
            &local_count64, 1, MPI_UINT64_T,
            counts64.data(), 1, MPI_UINT64_T,
            MPI_COMM_WORLD) != MPI_SUCCESS) {
      throw std::runtime_error("allgatherBytesBounded count Allgather failed");
    }

    BoundedMpiTransferPlan plan;
    std::vector<std::uint8_t> gathered;
    std::vector<std::uint8_t> round_receive_buffer;
    local_failure = nullptr;
    try {
      std::vector<std::size_t> logical_counts(counts64.size(), 0U);
      for (std::size_t rank = 0; rank < counts64.size(); ++rank) {
        logical_counts[rank] = core::checkedIntegralNarrow<std::size_t>(
            counts64[rank], "allgatherBytesBounded received byte count");
      }
      plan = planBoundedMpiTransferRounds(
          logical_counts,
          static_cast<std::size_t>(std::numeric_limits<int>::max()),
          mpiTransportRoundLimitBytes());
      gathered.resize(plan.logical_total_count);
      std::size_t maximum_round_count = 0U;
      for (const auto& round : plan.rounds) {
        maximum_round_count = std::max(maximum_round_count, round.round_count);
      }
      round_receive_buffer.resize(maximum_round_count);
      injectMpiTestFault(*this, "allgather_post_count");
    } catch (...) {
      local_failure = std::current_exception();
    }
    rethrowCollectivePreparationFailure(
        local_failure, "allgatherBytesBounded payload preparation");

    for (const BoundedMpiRoundLayout& round : plan.rounds) {
      const std::size_t local_rank = static_cast<std::size_t>(m_world_rank);
      const int send_count = round.counts[local_rank];
      const std::uint8_t* send_pointer = send_count == 0
          ? nullptr
          : local_bytes.data() + round.logical_offsets[local_rank];
      if (MPI_Allgatherv(
              const_cast<std::uint8_t*>(send_pointer), send_count, MPI_BYTE,
              round_receive_buffer.empty() ? nullptr : round_receive_buffer.data(),
              round.counts.data(), round.displacements.data(), MPI_BYTE,
              MPI_COMM_WORLD) != MPI_SUCCESS) {
        throw std::runtime_error("allgatherBytesBounded bounded payload Allgatherv failed");
      }
      for (std::size_t rank = 0; rank < round.counts.size(); ++rank) {
        const std::size_t count = static_cast<std::size_t>(round.counts[rank]);
        if (count == 0U) {
          continue;
        }
        // The planner checked the logical prefix and every consumed peer offset
        // before any payload round, so this addition cannot overflow here.
        const std::size_t destination_offset =
            plan.logical_displacements[rank] + round.logical_offsets[rank];
        std::memcpy(
            gathered.data() + destination_offset,
            round_receive_buffer.data() +
                static_cast<std::size_t>(round.displacements[rank]),
            count);
      }
    }
    return gathered;
  }
#endif
  return {local_bytes.begin(), local_bytes.end()};
}

namespace {

constexpr std::uint32_t k_tree_pseudo_wire_version = 1U;
constexpr std::size_t k_tree_pseudo_wire_bytes = 152U;

void appendTreeWireU32(std::vector<std::uint8_t>& bytes, std::uint32_t value) {
  for (unsigned shift = 0U; shift < 32U; shift += 8U) {
    bytes.push_back(static_cast<std::uint8_t>((value >> shift) & 0xffU));
  }
}

void appendTreeWireU64(std::vector<std::uint8_t>& bytes, std::uint64_t value) {
  for (unsigned shift = 0U; shift < 64U; shift += 8U) {
    bytes.push_back(static_cast<std::uint8_t>((value >> shift) & 0xffU));
  }
}

void appendTreeWireDouble(std::vector<std::uint8_t>& bytes, double value) {
  appendTreeWireU64(bytes, std::bit_cast<std::uint64_t>(value));
}

[[nodiscard]] std::uint32_t readTreeWireU32(std::span<const std::uint8_t> bytes, std::size_t& offset) {
  if (offset > bytes.size() || bytes.size() - offset < sizeof(std::uint32_t)) {
    throw std::runtime_error("tree pseudo-particle wire record is truncated");
  }
  std::uint32_t value = 0U;
  for (unsigned shift = 0U; shift < 32U; shift += 8U) {
    value |= static_cast<std::uint32_t>(bytes[offset++]) << shift;
  }
  return value;
}

[[nodiscard]] std::uint64_t readTreeWireU64(std::span<const std::uint8_t> bytes, std::size_t& offset) {
  if (offset > bytes.size() || bytes.size() - offset < sizeof(std::uint64_t)) {
    throw std::runtime_error("tree pseudo-particle wire record is truncated");
  }
  std::uint64_t value = 0U;
  for (unsigned shift = 0U; shift < 64U; shift += 8U) {
    value |= static_cast<std::uint64_t>(bytes[offset++]) << shift;
  }
  return value;
}

[[nodiscard]] double readTreeWireDouble(std::span<const std::uint8_t> bytes, std::size_t& offset) {
  return std::bit_cast<double>(readTreeWireU64(bytes, offset));
}

[[nodiscard, maybe_unused]] std::vector<std::uint8_t> encodeTreePseudoPackets(
    std::span<const TreePseudoParticlePacket> packets) {
  if (packets.size() > std::numeric_limits<std::size_t>::max() / k_tree_pseudo_wire_bytes) {
    throw std::overflow_error("tree pseudo-particle wire payload size overflows size_t");
  }
  std::vector<std::uint8_t> bytes;
  bytes.reserve(packets.size() * k_tree_pseudo_wire_bytes);
  for (const TreePseudoParticlePacket& packet : packets) {
    appendTreeWireU32(bytes, packet.descriptor.wire_version);
    appendTreeWireU64(bytes, packet.descriptor.pseudo_particle_id);
    appendTreeWireU32(bytes, static_cast<std::uint32_t>(packet.descriptor.source_rank));
    appendTreeWireU64(bytes, packet.descriptor.decomposition_epoch);
    appendTreeWireU32(bytes, packet.descriptor.derived_not_authoritative ? 1U : 0U);
    appendTreeWireU64(bytes, packet.descriptor.force_epoch);
    appendTreeWireU64(bytes, packet.descriptor.exchange_sequence);
    appendTreeWireU32(bytes, packet.geometry_frame);
    appendTreeWireDouble(bytes, packet.mass_code);
    appendTreeWireDouble(bytes, packet.center_x_comoving);
    appendTreeWireDouble(bytes, packet.center_y_comoving);
    appendTreeWireDouble(bytes, packet.center_z_comoving);
    appendTreeWireDouble(bytes, packet.min_x_comoving);
    appendTreeWireDouble(bytes, packet.max_x_comoving);
    appendTreeWireDouble(bytes, packet.min_y_comoving);
    appendTreeWireDouble(bytes, packet.max_y_comoving);
    appendTreeWireDouble(bytes, packet.min_z_comoving);
    appendTreeWireDouble(bytes, packet.max_z_comoving);
    appendTreeWireU64(bytes, packet.source_count);
    appendTreeWireU32(bytes, packet.hierarchy_level);
    appendTreeWireU32(bytes, packet.local_node_index);
    appendTreeWireU32(bytes, packet.child_count);
    appendTreeWireU32(bytes, packet.is_leaf);
  }
  if (bytes.size() != packets.size() * k_tree_pseudo_wire_bytes) {
    throw std::logic_error("tree pseudo-particle wire encoder size contract failed");
  }
  return bytes;
}

[[nodiscard, maybe_unused]] std::vector<TreePseudoParticlePacket> decodeTreePseudoPackets(
    std::span<const std::uint8_t> bytes) {
  if (bytes.size() % k_tree_pseudo_wire_bytes != 0U) {
    throw std::runtime_error("tree pseudo-particle wire payload is misaligned");
  }
  std::vector<TreePseudoParticlePacket> packets;
  packets.reserve(bytes.size() / k_tree_pseudo_wire_bytes);
  std::size_t offset = 0U;
  while (offset < bytes.size()) {
    TreePseudoParticlePacket packet;
    packet.descriptor.wire_version = readTreeWireU32(bytes, offset);
    packet.descriptor.pseudo_particle_id = readTreeWireU64(bytes, offset);
    packet.descriptor.source_rank = static_cast<int>(readTreeWireU32(bytes, offset));
    packet.descriptor.decomposition_epoch = readTreeWireU64(bytes, offset);
    packet.descriptor.derived_not_authoritative = readTreeWireU32(bytes, offset) != 0U;
    packet.descriptor.force_epoch = readTreeWireU64(bytes, offset);
    packet.descriptor.exchange_sequence = readTreeWireU64(bytes, offset);
    packet.geometry_frame = static_cast<std::uint8_t>(readTreeWireU32(bytes, offset));
    packet.mass_code = readTreeWireDouble(bytes, offset);
    packet.center_x_comoving = readTreeWireDouble(bytes, offset);
    packet.center_y_comoving = readTreeWireDouble(bytes, offset);
    packet.center_z_comoving = readTreeWireDouble(bytes, offset);
    packet.min_x_comoving = readTreeWireDouble(bytes, offset);
    packet.max_x_comoving = readTreeWireDouble(bytes, offset);
    packet.min_y_comoving = readTreeWireDouble(bytes, offset);
    packet.max_y_comoving = readTreeWireDouble(bytes, offset);
    packet.min_z_comoving = readTreeWireDouble(bytes, offset);
    packet.max_z_comoving = readTreeWireDouble(bytes, offset);
    packet.source_count = readTreeWireU64(bytes, offset);
    packet.hierarchy_level = readTreeWireU32(bytes, offset);
    packet.local_node_index = readTreeWireU32(bytes, offset);
    packet.child_count = static_cast<std::uint8_t>(readTreeWireU32(bytes, offset));
    packet.is_leaf = static_cast<std::uint8_t>(readTreeWireU32(bytes, offset));
    packets.push_back(packet);
  }
  return packets;
}

}  // namespace

std::vector<TreePseudoParticlePacket> executeBlockingTreePseudoParticleExchange(
    const MpiContext& mpi_context,
    const TreePseudoParticlePacket& local_packet) {
  std::exception_ptr local_validation_failure;
  try {
    validateTreePseudoParticlePacket(local_packet);
    if (local_packet.descriptor.source_rank != mpi_context.worldRank()) {
      throw std::invalid_argument("tree pseudo-particle packet source rank does not match MPI context");
    }
  } catch (...) {
    local_validation_failure = std::current_exception();
  }
  if (!mpi_context.isEnabled()) {
    if (local_validation_failure != nullptr) {
      std::rethrow_exception(local_validation_failure);
    }
    return {local_packet};
  }
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  std::uint64_t local_failure_vote = local_validation_failure != nullptr ? 1U : 0U;
  std::uint64_t global_failure_votes = 0U;
  MPI_Allreduce(
      &local_failure_vote, &global_failure_votes, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
  if (global_failure_votes != 0U) {
    throw std::runtime_error(
        "tree pseudo-particle exchange rejected invalid local input on one or more ranks");
  }
  std::vector<std::uint8_t> local_wire;
  std::exception_ptr local_encode_failure;
  try {
    local_wire = encodeTreePseudoPackets(
        std::span<const TreePseudoParticlePacket>(&local_packet, 1U));
    if (local_wire.size() != k_tree_pseudo_wire_bytes) {
      throw std::logic_error("tree pseudo-particle single-record wire size mismatch");
    }
  } catch (...) {
    local_encode_failure = std::current_exception();
  }
  local_failure_vote = local_encode_failure != nullptr ? 1U : 0U;
  global_failure_votes = 0U;
  MPI_Allreduce(
      &local_failure_vote, &global_failure_votes, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
  if (global_failure_votes != 0U) {
    throw std::runtime_error(
        "tree pseudo-particle exchange failed to encode local input on one or more ranks");
  }
  std::vector<std::uint8_t> gathered_wire(
      static_cast<std::size_t>(mpi_context.worldSize()) * k_tree_pseudo_wire_bytes,
      0U);
  MPI_Allgather(
      const_cast<std::uint8_t*>(local_wire.data()),
      static_cast<int>(k_tree_pseudo_wire_bytes),
      MPI_BYTE,
      gathered_wire.data(),
      static_cast<int>(k_tree_pseudo_wire_bytes),
      MPI_BYTE,
      MPI_COMM_WORLD);
  std::vector<TreePseudoParticlePacket> packets = decodeTreePseudoPackets(gathered_wire);
  for (int rank = 0; rank < mpi_context.worldSize(); ++rank) {
    const TreePseudoParticlePacket& packet = packets[static_cast<std::size_t>(rank)];
    validateTreePseudoParticlePacket(packet);
    if (packet.descriptor.source_rank != rank) {
      throw std::runtime_error("tree pseudo-particle exchange returned a packet with mismatched source rank");
    }
    if (packet.descriptor.exchange_sequence != local_packet.descriptor.exchange_sequence ||
        packet.descriptor.decomposition_epoch != local_packet.descriptor.decomposition_epoch ||
        packet.descriptor.force_epoch != local_packet.descriptor.force_epoch) {
      throw std::runtime_error("tree pseudo-particle exchange returned mixed protocol epochs");
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
  std::exception_ptr local_validation_failure;
  try {
    for (const TreePseudoParticlePacket& packet : local_packets) {
      validateTreePseudoParticlePacket(packet);
      if (packet.descriptor.source_rank != mpi_context.worldRank()) {
        throw std::invalid_argument("tree pseudo hierarchy packet source rank does not match MPI context");
      }
      if (packet.descriptor.exchange_sequence != exchange_sequence) {
        throw std::invalid_argument("tree pseudo hierarchy packet exchange sequence does not match exchange call");
      }
    }
  } catch (...) {
    local_validation_failure = std::current_exception();
  }
  if (!mpi_context.isEnabled()) {
    if (local_validation_failure != nullptr) {
      std::rethrow_exception(local_validation_failure);
    }
    return std::vector<TreePseudoParticlePacket>(local_packets.begin(), local_packets.end());
  }
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  int communicator_world_size = 1;
  int communicator_world_rank = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &communicator_world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &communicator_world_rank);
  const auto coordinate_failure = [](std::exception_ptr local_failure, std::string_view phase) {
    const std::uint64_t local_failure_vote = local_failure != nullptr ? 1U : 0U;
    std::uint64_t global_failure_votes = 0U;
    MPI_Allreduce(
        &local_failure_vote, &global_failure_votes, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
    if (global_failure_votes == 0U) {
      return;
    }
    if (local_failure != nullptr) {
      std::rethrow_exception(local_failure);
    }
    throw std::runtime_error(
        "tree pseudo hierarchy exchange failed during " + std::string(phase) +
        " on one or more peer ranks");
  };

  std::exception_ptr communicator_validation_failure;
  try {
    if (mpi_context.worldSize() != communicator_world_size ||
        mpi_context.worldRank() != communicator_world_rank) {
      throw std::invalid_argument(
          "tree pseudo hierarchy MPI context does not match MPI_COMM_WORLD");
    }
  } catch (...) {
    communicator_validation_failure = std::current_exception();
  }
  coordinate_failure(communicator_validation_failure, "communicator validation");
  coordinate_failure(local_validation_failure, "local input validation");

  const int world_size = communicator_world_size;
  std::uint64_t min_exchange_sequence = 0U;
  std::uint64_t max_exchange_sequence = 0U;
  MPI_Allreduce(
      &exchange_sequence, &min_exchange_sequence, 1, MPI_UINT64_T, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(
      &exchange_sequence, &max_exchange_sequence, 1, MPI_UINT64_T, MPI_MAX, MPI_COMM_WORLD);
  if (min_exchange_sequence != max_exchange_sequence) {
    throw std::runtime_error("tree pseudo hierarchy ranks disagree on exchange sequence");
  }
  std::uint64_t local_count = 0U;
  std::vector<std::uint64_t> counts64;
  std::exception_ptr count_buffer_failure;
  try {
    if (local_packets.size() > static_cast<std::size_t>(std::numeric_limits<std::uint64_t>::max())) {
      throw std::overflow_error("tree pseudo hierarchy local packet count exceeds uint64_t");
    }
    local_count = static_cast<std::uint64_t>(local_packets.size());
    counts64.assign(static_cast<std::size_t>(world_size), 0U);
  } catch (...) {
    count_buffer_failure = std::current_exception();
  }
  coordinate_failure(count_buffer_failure, "count-buffer preparation");
  MPI_Allgather(
      &local_count,
      1,
      MPI_UINT64_T,
      counts64.data(),
      1,
      MPI_UINT64_T,
      MPI_COMM_WORLD);
  std::vector<std::uint8_t> local_wire;
  std::exception_ptr local_encode_failure;
  try {
    local_wire = encodeTreePseudoPackets(local_packets);
    const std::size_t expected_local_bytes = core::checkedSizeMultiply(
        core::checkedIntegralNarrow<std::size_t>(
            counts64[static_cast<std::size_t>(communicator_world_rank)],
            "tree pseudo hierarchy local packet count"),
        k_tree_pseudo_wire_bytes,
        "tree pseudo hierarchy encoded byte count");
    if (local_wire.size() != expected_local_bytes) {
      throw std::logic_error(
          "tree pseudo hierarchy encoded byte count disagrees with gathered packet count");
    }
  } catch (...) {
    local_encode_failure = std::current_exception();
  }
  coordinate_failure(local_encode_failure, "local wire encoding");

  std::vector<std::uint8_t> result_wire =
      mpi_context.allgatherBytesBounded(local_wire);
  std::vector<TreePseudoParticlePacket> result;
  std::exception_ptr received_wire_validation_failure;
  try {
    result = decodeTreePseudoPackets(result_wire);
    std::vector<std::uint64_t> per_rank_count(static_cast<std::size_t>(world_size), 0U);
    std::vector<std::uint32_t> per_rank_root_count(static_cast<std::size_t>(world_size), 0U);
    std::vector<std::uint8_t> per_rank_authoritative(static_cast<std::size_t>(world_size), 0U);
    std::vector<std::uint8_t> per_rank_derived(static_cast<std::size_t>(world_size), 0U);
    std::vector<std::unordered_set<std::uint64_t>> per_rank_ids(static_cast<std::size_t>(world_size));
    bool have_epoch_contract = false;
    std::uint64_t expected_decomposition_epoch = 0U;
    std::uint64_t expected_force_epoch = 0U;
    std::uint8_t expected_geometry_frame = 0U;
    for (const TreePseudoParticlePacket& packet : result) {
      validateTreePseudoParticlePacket(packet);
      if (packet.descriptor.source_rank < 0 || packet.descriptor.source_rank >= world_size) {
        throw std::runtime_error("tree pseudo hierarchy exchange returned packet with invalid source rank");
      }
      if (packet.descriptor.exchange_sequence != exchange_sequence) {
        throw std::runtime_error("tree pseudo hierarchy exchange returned a stale exchange sequence");
      }
      if (!have_epoch_contract) {
        expected_decomposition_epoch = packet.descriptor.decomposition_epoch;
        expected_force_epoch = packet.descriptor.force_epoch;
        expected_geometry_frame = packet.geometry_frame;
        have_epoch_contract = true;
      } else if (packet.descriptor.decomposition_epoch != expected_decomposition_epoch ||
                 packet.descriptor.force_epoch != expected_force_epoch ||
                 packet.geometry_frame != expected_geometry_frame) {
        throw std::runtime_error("tree pseudo hierarchy exchange returned mixed epochs or geometry frames");
      }
      const std::size_t source_rank = static_cast<std::size_t>(packet.descriptor.source_rank);
      if (packet.descriptor.derived_not_authoritative) {
        per_rank_derived[source_rank] = 1U;
      } else {
        per_rank_authoritative[source_rank] = 1U;
      }
      if (!per_rank_ids[source_rank].insert(packet.descriptor.pseudo_particle_id).second) {
        throw std::runtime_error("tree pseudo hierarchy exchange returned duplicate node identity");
      }
      ++per_rank_count[source_rank];
      if (packet.hierarchy_level == 0U) {
        ++per_rank_root_count[source_rank];
      }
    }
    for (int rank = 0; rank < world_size; ++rank) {
      if (per_rank_count[static_cast<std::size_t>(rank)] != counts64[static_cast<std::size_t>(rank)]) {
        throw std::runtime_error("tree pseudo hierarchy exchange source-rank coverage mismatch");
      }
      if (per_rank_authoritative[static_cast<std::size_t>(rank)] != 0U &&
          per_rank_derived[static_cast<std::size_t>(rank)] != 0U) {
        throw std::runtime_error(
            "tree top-domain exchange cannot mix authoritative and derived geometry for one rank");
      }
      if (per_rank_derived[static_cast<std::size_t>(rank)] != 0U &&
          per_rank_root_count[static_cast<std::size_t>(rank)] != 1U) {
        throw std::runtime_error(
            "tree pseudo hierarchy exchange requires exactly one root descriptor per participating rank");
      }
    }
  } catch (...) {
    received_wire_validation_failure = std::current_exception();
  }
  coordinate_failure(received_wire_validation_failure, "received wire decoding and validation");
  return result;
#else
  throw std::runtime_error("tree pseudo hierarchy exchange requires MPI support when MPI context is enabled");
#endif
}


namespace {

struct AmrRankEnvelopeRecord {
  int rank = -1;
  std::uint8_t has_patches = 0;
  std::uint64_t decomposition_epoch = 0;
  double min_x_comoving = 0.0;
  double max_x_comoving = 0.0;
  double min_y_comoving = 0.0;
  double max_y_comoving = 0.0;
  double min_z_comoving = 0.0;
  double max_z_comoving = 0.0;
  double max_cell_width_comoving = 0.0;
};

[[nodiscard]] double amrPatchMaxX(const AmrPatchPayloadRecord& record) noexcept {
  return record.origin_x_comoving + record.extent_x_comoving;
}

[[nodiscard]] double amrPatchMaxY(const AmrPatchPayloadRecord& record) noexcept {
  return record.origin_y_comoving + record.extent_y_comoving;
}

[[nodiscard]] double amrPatchMaxZ(const AmrPatchPayloadRecord& record) noexcept {
  return record.origin_z_comoving + record.extent_z_comoving;
}

[[nodiscard]] double amrPatchMaxCellWidth(const AmrPatchPayloadRecord& record) noexcept {
  const double dx = record.extent_x_comoving / static_cast<double>(std::max<std::uint16_t>(record.cell_dim_x, 1U));
  const double dy = record.extent_y_comoving / static_cast<double>(std::max<std::uint16_t>(record.cell_dim_y, 1U));
  const double dz = record.extent_z_comoving / static_cast<double>(std::max<std::uint16_t>(record.cell_dim_z, 1U));
  return std::max({dx, dy, dz});
}

[[nodiscard, maybe_unused]] AmrRankEnvelopeRecord buildLocalAmrEnvelope(
    std::span<const AmrPatchPayloadRecord> records,
    int world_rank) {
  AmrRankEnvelopeRecord envelope;
  envelope.rank = world_rank;
  if (records.empty()) {
    return envelope;
  }
  envelope.has_patches = 1U;
  envelope.decomposition_epoch = records.front().decomposition_epoch;
  envelope.min_x_comoving = std::numeric_limits<double>::infinity();
  envelope.min_y_comoving = std::numeric_limits<double>::infinity();
  envelope.min_z_comoving = std::numeric_limits<double>::infinity();
  envelope.max_x_comoving = -std::numeric_limits<double>::infinity();
  envelope.max_y_comoving = -std::numeric_limits<double>::infinity();
  envelope.max_z_comoving = -std::numeric_limits<double>::infinity();
  for (const AmrPatchPayloadRecord& record : records) {
    validateAmrPatchPayloadRecord(record);
    if (record.owner_rank != world_rank) {
      throw std::invalid_argument("directed AMR envelope can only be built from local authoritative patches");
    }
    if (record.decomposition_epoch != envelope.decomposition_epoch) {
      throw std::invalid_argument("directed AMR envelope found mixed decomposition epochs in one exchange");
    }
    envelope.min_x_comoving = std::min(envelope.min_x_comoving, record.origin_x_comoving);
    envelope.min_y_comoving = std::min(envelope.min_y_comoving, record.origin_y_comoving);
    envelope.min_z_comoving = std::min(envelope.min_z_comoving, record.origin_z_comoving);
    envelope.max_x_comoving = std::max(envelope.max_x_comoving, amrPatchMaxX(record));
    envelope.max_y_comoving = std::max(envelope.max_y_comoving, amrPatchMaxY(record));
    envelope.max_z_comoving = std::max(envelope.max_z_comoving, amrPatchMaxZ(record));
    envelope.max_cell_width_comoving = std::max(envelope.max_cell_width_comoving, amrPatchMaxCellWidth(record));
  }
  return envelope;
}

[[nodiscard]] bool amrIntervalsOverlap(double amin, double amax, double bmin, double bmax, double reach) noexcept {
  return (amin - reach) <= bmax && (bmin - reach) <= amax;
}

[[nodiscard, maybe_unused]] bool amrEnvelopeMayNeedPeer(
    const AmrRankEnvelopeRecord& local,
    const AmrRankEnvelopeRecord& peer) noexcept {
  if (local.has_patches == 0U || peer.has_patches == 0U || local.rank == peer.rank) {
    return false;
  }
  const double reach = std::max(local.max_cell_width_comoving, peer.max_cell_width_comoving);
  return amrIntervalsOverlap(local.min_x_comoving, local.max_x_comoving, peer.min_x_comoving, peer.max_x_comoving, reach) &&
      amrIntervalsOverlap(local.min_y_comoving, local.max_y_comoving, peer.min_y_comoving, peer.max_y_comoving, reach) &&
      amrIntervalsOverlap(local.min_z_comoving, local.max_z_comoving, peer.min_z_comoving, peer.max_z_comoving, reach);
}

[[nodiscard]] bool amrPatchesMayShareInterface(
    const AmrPatchPayloadRecord& lhs,
    const AmrPatchPayloadRecord& rhs) noexcept {
  const double reach = std::max(amrPatchMaxCellWidth(lhs), amrPatchMaxCellWidth(rhs));
  return amrIntervalsOverlap(lhs.origin_x_comoving, amrPatchMaxX(lhs), rhs.origin_x_comoving, amrPatchMaxX(rhs), reach) &&
      amrIntervalsOverlap(lhs.origin_y_comoving, amrPatchMaxY(lhs), rhs.origin_y_comoving, amrPatchMaxY(rhs), reach) &&
      amrIntervalsOverlap(lhs.origin_z_comoving, amrPatchMaxZ(lhs), rhs.origin_z_comoving, amrPatchMaxZ(rhs), reach);
}

[[nodiscard, maybe_unused]] std::vector<int> discoverCandidateAmrPeers(
    const MpiContext& mpi_context,
    std::span<const AmrPatchPayloadRecord> local_patch_records,
    DirectedAmrExchangeDiagnostics* diagnostics) {
  const int world_size = mpi_context.worldSize();
#if !defined(COSMOSIM_ENABLE_MPI) || !COSMOSIM_ENABLE_MPI
  (void)local_patch_records;
#endif
  if (!mpi_context.isEnabled() || world_size <= 1) {
    if (diagnostics != nullptr) {
      diagnostics->control_plane_bytes += sizeof(AmrRankEnvelopeRecord);
    }
    return {};
  }
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  const int world_rank = mpi_context.worldRank();
  const AmrRankEnvelopeRecord local_envelope = buildLocalAmrEnvelope(local_patch_records, world_rank);
  static_assert(std::is_trivially_copyable_v<AmrRankEnvelopeRecord>);
  std::vector<AmrRankEnvelopeRecord> envelopes(static_cast<std::size_t>(world_size));
  MPI_Allgather(
      const_cast<AmrRankEnvelopeRecord*>(&local_envelope),
      static_cast<int>(sizeof(AmrRankEnvelopeRecord)),
      MPI_BYTE,
      envelopes.data(),
      static_cast<int>(sizeof(AmrRankEnvelopeRecord)),
      MPI_BYTE,
      MPI_COMM_WORLD);
  if (diagnostics != nullptr) {
    diagnostics->control_plane_bytes += static_cast<std::uint64_t>(world_size) * sizeof(AmrRankEnvelopeRecord);
  }
  std::vector<int> peers;
  for (const AmrRankEnvelopeRecord& envelope : envelopes) {
    if (envelope.rank < 0 || envelope.rank >= world_size) {
      throw std::runtime_error("directed AMR control-plane envelope returned invalid rank metadata");
    }
    if (amrEnvelopeMayNeedPeer(local_envelope, envelope)) {
      peers.push_back(envelope.rank);
    }
  }
  std::sort(peers.begin(), peers.end());
  peers.erase(std::unique(peers.begin(), peers.end()), peers.end());
  if (diagnostics != nullptr) {
    diagnostics->candidate_peer_count = static_cast<std::uint64_t>(peers.size());
  }
  return peers;
#else
  throw std::runtime_error("directed AMR peer discovery requires MPI support when MPI context is enabled");
#endif
}

template <typename T>
[[nodiscard]] std::vector<T> exchangePodRecordsWithPeer(
    const MpiContext& mpi_context,
    int peer_rank,
    std::span<const T> local_records,
    int count_tag_base,
    int payload_tag_base,
    std::uint64_t exchange_sequence,
    const char* caller,
    std::uint64_t* transport_round_count = nullptr) {
  static_assert(std::is_trivially_copyable_v<T>);
  if (peer_rank < 0 || peer_rank >= mpi_context.worldSize() || peer_rank == mpi_context.worldRank()) {
    throw std::invalid_argument(std::string(caller) + ": invalid peer rank");
  }
  if (!mpi_context.isEnabled()) {
    if (!local_records.empty()) {
      throw std::runtime_error(std::string(caller) + ": non-empty peer exchange requires MPI");
    }
    return {};
  }
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  const std::uint64_t send_count = static_cast<std::uint64_t>(local_records.size());
  std::uint64_t recv_count = 0;
  if (MPI_Sendrecv(
          const_cast<std::uint64_t*>(&send_count),
          1,
          MPI_UINT64_T,
          peer_rank,
          ghostExchangeSequencedTag(count_tag_base, mpi_context.worldRank(), peer_rank, exchange_sequence),
          &recv_count,
          1,
          MPI_UINT64_T,
          peer_rank,
          ghostExchangeSequencedTag(count_tag_base, mpi_context.worldRank(), peer_rank, exchange_sequence),
          MPI_COMM_WORLD,
          MPI_STATUS_IGNORE) != MPI_SUCCESS) {
    throw std::runtime_error(std::string(caller) + ": directed AMR record-count Sendrecv failed");
  }

  std::vector<T> received;
  std::size_t receive_count = 0U;
  std::size_t records_per_round = 0U;
  std::exception_ptr local_preparation_failure;
  try {
    receive_count = core::checkedIntegralNarrow<std::size_t>(
        recv_count, std::string(caller) + " receive record count");
    received.resize(receive_count);
    const std::size_t round_limit_bytes = std::min(
        mpiTransportRoundLimitBytes(),
        static_cast<std::size_t>(std::numeric_limits<int>::max()));
    records_per_round = round_limit_bytes / sizeof(T);
    if (records_per_round == 0U) {
      throw std::overflow_error(
          std::string(caller) + ": one record exceeds the bounded MPI transport round");
    }
  } catch (...) {
    local_preparation_failure = std::current_exception();
  }

  int local_ready = local_preparation_failure == nullptr ? 1 : 0;
  int peer_ready = 0;
  if (MPI_Sendrecv(
          &local_ready,
          1,
          MPI_INT,
          peer_rank,
          ghostExchangeSequencedTag(
              count_tag_base + 100, mpi_context.worldRank(), peer_rank, exchange_sequence),
          &peer_ready,
          1,
          MPI_INT,
          peer_rank,
          ghostExchangeSequencedTag(
              count_tag_base + 100, mpi_context.worldRank(), peer_rank, exchange_sequence),
          MPI_COMM_WORLD,
          MPI_STATUS_IGNORE) != MPI_SUCCESS) {
    throw std::runtime_error(std::string(caller) + ": directed AMR preparation-readiness Sendrecv failed");
  }
  if (local_ready == 0) {
    std::rethrow_exception(local_preparation_failure);
  }
  if (peer_ready == 0) {
    throw std::runtime_error(std::string(caller) + ": peer rejected directed AMR receive preparation");
  }

  std::size_t send_offset = 0U;
  std::size_t recv_offset = 0U;
  std::uint64_t rounds = 0U;
  while (send_offset < local_records.size() || recv_offset < received.size()) {
    const std::size_t send_records = std::min(records_per_round, local_records.size() - send_offset);
    const std::size_t recv_records = std::min(records_per_round, received.size() - recv_offset);
    const std::size_t send_bytes = core::checkedSizeMultiply(
        send_records, sizeof(T), std::string(caller) + " bounded send bytes");
    const std::size_t recv_bytes = core::checkedSizeMultiply(
        recv_records, sizeof(T), std::string(caller) + " bounded receive bytes");
    const int send_count_bytes = core::checkedIntegralNarrow<int>(
        send_bytes, std::string(caller) + " bounded send count");
    const int recv_count_bytes = core::checkedIntegralNarrow<int>(
        recv_bytes, std::string(caller) + " bounded receive count");
    T* recv_pointer = recv_records == 0U ? nullptr : received.data() + recv_offset;
    const T* send_pointer = send_records == 0U ? nullptr : local_records.data() + send_offset;
    if (MPI_Sendrecv(
            const_cast<T*>(send_pointer),
            send_count_bytes,
            MPI_BYTE,
            peer_rank,
            ghostExchangeSequencedTag(payload_tag_base, mpi_context.worldRank(), peer_rank, exchange_sequence),
            recv_pointer,
            recv_count_bytes,
            MPI_BYTE,
            peer_rank,
            ghostExchangeSequencedTag(payload_tag_base, mpi_context.worldRank(), peer_rank, exchange_sequence),
            MPI_COMM_WORLD,
            MPI_STATUS_IGNORE) != MPI_SUCCESS) {
      throw std::runtime_error(std::string(caller) + ": bounded directed AMR payload Sendrecv failed");
    }
    send_offset += send_records;
    recv_offset += recv_records;
    ++rounds;
  }
  if (transport_round_count != nullptr) {
    *transport_round_count += rounds;
  }
  return received;
#else
  throw std::runtime_error(std::string(caller) + ": MPI support is not compiled in");
#endif
}

void exchangeAmrPatchCellRecordStreamWithPeer(
    const MpiContext& mpi_context,
    int peer_rank,
    std::uint64_t logical_send_count,
    std::span<const AmrPatchBoundaryCellRequest> local_requests,
    std::span<const AmrPatchPayloadRecord> remote_patches,
    std::span<const AmrPatchBoundaryCellRequest> remote_requests,
    const DirectedAmrPatchCellPayloadProvider& provider,
    const DirectedAmrPatchCellPayloadConsumer& consumer,
    const DirectedAmrPatchCellAdmission& admission,
    int count_tag_base,
    int payload_tag_base,
    std::uint64_t exchange_sequence,
    std::size_t local_transport_round_limit_bytes,
    DirectedAmrExchangeDiagnostics& diagnostics) {
  if (!provider || !consumer) {
    throw std::invalid_argument(
        "directed AMR patch-cell streaming requires both producer and consumer callbacks");
  }
  if (!mpi_context.isEnabled()) {
    if (logical_send_count != 0U) {
      throw std::runtime_error(
          "directed AMR patch-cell streaming requires MPI for non-empty payloads");
    }
    return;
  }
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  std::uint64_t logical_receive_count = 0U;
  if (MPI_Sendrecv(
          &logical_send_count,
          1,
          MPI_UINT64_T,
          peer_rank,
          ghostExchangeSequencedTag(count_tag_base, mpi_context.worldRank(), peer_rank, exchange_sequence),
          &logical_receive_count,
          1,
          MPI_UINT64_T,
          peer_rank,
          ghostExchangeSequencedTag(count_tag_base, mpi_context.worldRank(), peer_rank, exchange_sequence),
          MPI_COMM_WORLD,
          MPI_STATUS_IGNORE) != MPI_SUCCESS) {
    throw std::runtime_error("directed AMR patch-cell record-count Sendrecv failed");
  }

  const std::size_t local_limit = std::min(
      local_transport_round_limit_bytes == 0U
          ? mpiTransportRoundLimitBytes()
          : local_transport_round_limit_bytes,
      static_cast<std::size_t>(std::numeric_limits<int>::max()));
  std::uint64_t local_limit_u64 = static_cast<std::uint64_t>(local_limit);
  std::uint64_t peer_limit_u64 = 0U;
  if (MPI_Sendrecv(
          &local_limit_u64,
          1,
          MPI_UINT64_T,
          peer_rank,
          ghostExchangeSequencedTag(count_tag_base + 1, mpi_context.worldRank(), peer_rank, exchange_sequence),
          &peer_limit_u64,
          1,
          MPI_UINT64_T,
          peer_rank,
          ghostExchangeSequencedTag(count_tag_base + 1, mpi_context.worldRank(), peer_rank, exchange_sequence),
          MPI_COMM_WORLD,
          MPI_STATUS_IGNORE) != MPI_SUCCESS) {
    throw std::runtime_error("directed AMR patch-cell transport-limit Sendrecv failed");
  }
  const std::size_t agreed_limit = std::min(
      local_limit,
      core::checkedIntegralNarrow<std::size_t>(
          peer_limit_u64, "directed AMR peer transport limit"));

  DirectedAmrPatchCellTransferPlan send_plan;
  DirectedAmrPatchCellTransferPlan receive_plan;
  std::exception_ptr local_admission_failure;
  try {
    send_plan = planDirectedAmrPatchCellTransfer(
        logical_send_count, agreed_limit);
    receive_plan = planDirectedAmrPatchCellTransfer(
        logical_receive_count, agreed_limit);
    if (admission) {
      admission(peer_rank, remote_patches, remote_requests, logical_receive_count);
    }
  } catch (...) {
    local_admission_failure = std::current_exception();
  }
  int local_admitted = local_admission_failure == nullptr ? 1 : 0;
  int peer_admitted = 0;
  if (MPI_Sendrecv(
          &local_admitted,
          1,
          MPI_INT,
          peer_rank,
          ghostExchangeSequencedTag(count_tag_base + 2, mpi_context.worldRank(), peer_rank, exchange_sequence),
          &peer_admitted,
          1,
          MPI_INT,
          peer_rank,
          ghostExchangeSequencedTag(count_tag_base + 2, mpi_context.worldRank(), peer_rank, exchange_sequence),
          MPI_COMM_WORLD,
          MPI_STATUS_IGNORE) != MPI_SUCCESS) {
    throw std::runtime_error("directed AMR patch-cell admission-readiness Sendrecv failed");
  }
  if (local_admitted == 0) {
    std::rethrow_exception(local_admission_failure);
  }
  if (peer_admitted == 0) {
    throw std::runtime_error("directed AMR peer rejected patch-cell memory admission");
  }

  const std::uint64_t round_count = std::max(send_plan.round_count, receive_plan.round_count);
  std::vector<AmrPatchCellPayloadRecord> send_chunk;
  std::vector<AmrPatchCellPayloadRecord> receive_chunk;
  const std::size_t max_records = send_plan.records_per_round;
  send_chunk.reserve(max_records);
  receive_chunk.reserve(max_records);
  diagnostics.patch_cell_send_capacity_high_water_bytes = std::max(
      diagnostics.patch_cell_send_capacity_high_water_bytes,
      static_cast<std::uint64_t>(core::checkedSizeMultiply(
          send_chunk.capacity(), sizeof(AmrPatchCellPayloadRecord),
          "directed AMR streamed send capacity")));
  diagnostics.patch_cell_receive_capacity_high_water_bytes = std::max(
      diagnostics.patch_cell_receive_capacity_high_water_bytes,
      static_cast<std::uint64_t>(core::checkedSizeMultiply(
          receive_chunk.capacity(), sizeof(AmrPatchCellPayloadRecord),
          "directed AMR streamed receive capacity")));
  diagnostics.communication_workspace_high_water_bytes = std::max(
      diagnostics.communication_workspace_high_water_bytes,
      core::checkedMemoryBytesAdd(
          diagnostics.patch_cell_send_capacity_high_water_bytes,
          diagnostics.patch_cell_receive_capacity_high_water_bytes,
          "directed AMR streamed simultaneous transport capacity"));

  std::uint64_t send_offset = 0U;
  std::uint64_t receive_offset = 0U;
  for (std::uint64_t round = 0U; round < round_count; ++round) {
    const std::size_t send_records = round < send_plan.round_count
        ? static_cast<std::size_t>(std::min<std::uint64_t>(
              send_plan.logical_record_count - send_offset,
              static_cast<std::uint64_t>(send_plan.records_per_round)))
        : 0U;
    const std::size_t receive_records = round < receive_plan.round_count
        ? static_cast<std::size_t>(std::min<std::uint64_t>(
              receive_plan.logical_record_count - receive_offset,
              static_cast<std::uint64_t>(receive_plan.records_per_round)))
        : 0U;

    std::exception_ptr local_producer_failure;
    try {
      send_chunk.clear();
      if (send_records != 0U) {
        provider(local_requests, send_offset, send_records, send_chunk);
      }
      if (send_chunk.size() != send_records) {
        throw std::runtime_error(
            "directed AMR patch-cell producer did not provide the requested streamed record count");
      }
    } catch (...) {
      local_producer_failure = std::current_exception();
    }
    int local_ready = local_producer_failure == nullptr ? 1 : 0;
    int peer_ready = 0;
    if (MPI_Sendrecv(
            &local_ready,
            1,
            MPI_INT,
            peer_rank,
            ghostExchangeSequencedTag(count_tag_base + 3, mpi_context.worldRank(), peer_rank, exchange_sequence),
            &peer_ready,
            1,
            MPI_INT,
            peer_rank,
            ghostExchangeSequencedTag(count_tag_base + 3, mpi_context.worldRank(), peer_rank, exchange_sequence),
            MPI_COMM_WORLD,
            MPI_STATUS_IGNORE) != MPI_SUCCESS) {
      throw std::runtime_error("directed AMR patch-cell producer-readiness Sendrecv failed");
    }
    if (local_ready == 0) {
      std::rethrow_exception(local_producer_failure);
    }
    if (peer_ready == 0) {
      throw std::runtime_error("directed AMR peer rejected patch-cell producer preparation");
    }

    receive_chunk.resize(receive_records);
    const std::size_t send_bytes = core::checkedSizeMultiply(
        send_records, sizeof(AmrPatchCellPayloadRecord),
        "directed AMR streamed send bytes");
    const std::size_t receive_bytes = core::checkedSizeMultiply(
        receive_records, sizeof(AmrPatchCellPayloadRecord),
        "directed AMR streamed receive bytes");
    const int send_count_bytes = core::checkedIntegralNarrow<int>(
        send_bytes, "directed AMR streamed send count");
    const int receive_count_bytes = core::checkedIntegralNarrow<int>(
        receive_bytes, "directed AMR streamed receive count");
    if (MPI_Sendrecv(
            send_records == 0U ? nullptr : send_chunk.data(),
            send_count_bytes,
            MPI_BYTE,
            peer_rank,
            ghostExchangeSequencedTag(payload_tag_base, mpi_context.worldRank(), peer_rank, exchange_sequence),
            receive_records == 0U ? nullptr : receive_chunk.data(),
            receive_count_bytes,
            MPI_BYTE,
            peer_rank,
            ghostExchangeSequencedTag(payload_tag_base, mpi_context.worldRank(), peer_rank, exchange_sequence),
            MPI_COMM_WORLD,
            MPI_STATUS_IGNORE) != MPI_SUCCESS) {
      throw std::runtime_error("directed AMR bounded streamed patch-cell Sendrecv failed");
    }

    std::exception_ptr local_consumer_failure;
    try {
      if (!receive_chunk.empty()) {
        consumer(peer_rank, receive_chunk);
      }
    } catch (...) {
      local_consumer_failure = std::current_exception();
    }
    int local_consumed = local_consumer_failure == nullptr ? 1 : 0;
    int peer_consumed = 0;
    if (MPI_Sendrecv(
            &local_consumed,
            1,
            MPI_INT,
            peer_rank,
            ghostExchangeSequencedTag(count_tag_base + 4, mpi_context.worldRank(), peer_rank, exchange_sequence),
            &peer_consumed,
            1,
            MPI_INT,
            peer_rank,
            ghostExchangeSequencedTag(count_tag_base + 4, mpi_context.worldRank(), peer_rank, exchange_sequence),
            MPI_COMM_WORLD,
            MPI_STATUS_IGNORE) != MPI_SUCCESS) {
      throw std::runtime_error("directed AMR patch-cell consumer-readiness Sendrecv failed");
    }
    if (local_consumed == 0) {
      std::rethrow_exception(local_consumer_failure);
    }
    if (peer_consumed == 0) {
      throw std::runtime_error("directed AMR peer rejected streamed patch-cell consumption");
    }

    send_offset += static_cast<std::uint64_t>(send_records);
    receive_offset += static_cast<std::uint64_t>(receive_records);
    ++diagnostics.patch_cell_transport_round_count;
  }
  if (send_offset != logical_send_count) {
    throw std::logic_error("directed AMR streamed patch-cell sender did not cover logical payload");
  }
  if (receive_offset != logical_receive_count) {
    throw std::logic_error("directed AMR streamed patch-cell receiver did not cover logical payload");
  }
#else
  (void)peer_rank;
  (void)logical_send_count;
  (void)local_requests;
  (void)remote_patches;
  (void)remote_requests;
  (void)provider;
  (void)consumer;
  (void)admission;
  (void)count_tag_base;
  (void)payload_tag_base;
  (void)exchange_sequence;
  (void)local_transport_round_limit_bytes;
  (void)diagnostics;
  throw std::runtime_error("directed AMR patch-cell streaming requires MPI support");
#endif
}

[[nodiscard]] std::uint8_t amrBoundaryFaceBit(AmrPatchBoundaryFace face) noexcept {
  return static_cast<std::uint8_t>(face);
}

[[nodiscard]] double amrInterfaceTolerance(
    const AmrPatchPayloadRecord& lhs,
    const AmrPatchPayloadRecord& rhs) noexcept {
  const double scale = std::max({
      1.0,
      std::abs(lhs.origin_x_comoving), std::abs(amrPatchMaxX(lhs)),
      std::abs(lhs.origin_y_comoving), std::abs(amrPatchMaxY(lhs)),
      std::abs(lhs.origin_z_comoving), std::abs(amrPatchMaxZ(lhs)),
      std::abs(rhs.origin_x_comoving), std::abs(amrPatchMaxX(rhs)),
      std::abs(rhs.origin_y_comoving), std::abs(amrPatchMaxY(rhs)),
      std::abs(rhs.origin_z_comoving), std::abs(amrPatchMaxZ(rhs))});
  return 1.0e-10 * scale;
}

[[nodiscard]] bool amrIntervalsShareArea(
    double lhs_min,
    double lhs_max,
    double rhs_min,
    double rhs_max,
    double tolerance) noexcept {
  return std::min(lhs_max, rhs_max) - std::max(lhs_min, rhs_min) > tolerance;
}

[[nodiscard]] std::uint8_t amrPatchBoundaryMaskForPeer(
    const AmrPatchPayloadRecord& local,
    const AmrPatchPayloadRecord& remote) noexcept {
  if (local.owner_rank == remote.owner_rank) {
    return 0U;
  }
  const std::array<double, 3> local_min{
      local.origin_x_comoving, local.origin_y_comoving, local.origin_z_comoving};
  const std::array<double, 3> local_max{
      amrPatchMaxX(local), amrPatchMaxY(local), amrPatchMaxZ(local)};
  const std::array<double, 3> remote_min{
      remote.origin_x_comoving, remote.origin_y_comoving, remote.origin_z_comoving};
  const std::array<double, 3> remote_max{
      amrPatchMaxX(remote), amrPatchMaxY(remote), amrPatchMaxZ(remote)};
  const double tolerance = amrInterfaceTolerance(local, remote);
  std::uint8_t mask = 0U;
  for (std::size_t axis = 0; axis < 3U; ++axis) {
    const std::size_t transverse_a = (axis + 1U) % 3U;
    const std::size_t transverse_b = (axis + 2U) % 3U;
    if (!amrIntervalsShareArea(
            local_min[transverse_a], local_max[transverse_a],
            remote_min[transverse_a], remote_max[transverse_a], tolerance) ||
        !amrIntervalsShareArea(
            local_min[transverse_b], local_max[transverse_b],
            remote_min[transverse_b], remote_max[transverse_b], tolerance)) {
      continue;
    }
    if (std::abs(local_min[axis] - remote_max[axis]) <= tolerance) {
      mask |= amrBoundaryFaceBit(axis == 0U ? AmrPatchBoundaryFace::kXLower :
          (axis == 1U ? AmrPatchBoundaryFace::kYLower : AmrPatchBoundaryFace::kZLower));
    }
    if (std::abs(local_max[axis] - remote_min[axis]) <= tolerance) {
      mask |= amrBoundaryFaceBit(axis == 0U ? AmrPatchBoundaryFace::kXUpper :
          (axis == 1U ? AmrPatchBoundaryFace::kYUpper : AmrPatchBoundaryFace::kZUpper));
    }
  }
  return mask;
}

[[nodiscard]] std::size_t amrBoundaryFaceIndex(AmrPatchBoundaryFace face) noexcept {
  switch (face) {
    case AmrPatchBoundaryFace::kXLower: return 0U;
    case AmrPatchBoundaryFace::kXUpper: return 1U;
    case AmrPatchBoundaryFace::kYLower: return 2U;
    case AmrPatchBoundaryFace::kYUpper: return 3U;
    case AmrPatchBoundaryFace::kZLower: return 4U;
    case AmrPatchBoundaryFace::kZUpper: return 5U;
  }
  return 0U;
}

[[nodiscard]] std::size_t amrBoundaryFaceAxis(AmrPatchBoundaryFace face) noexcept {
  switch (face) {
    case AmrPatchBoundaryFace::kXLower:
    case AmrPatchBoundaryFace::kXUpper:
      return 0U;
    case AmrPatchBoundaryFace::kYLower:
    case AmrPatchBoundaryFace::kYUpper:
      return 1U;
    case AmrPatchBoundaryFace::kZLower:
    case AmrPatchBoundaryFace::kZUpper:
      return 2U;
  }
  return 0U;
}

[[nodiscard]] std::uint16_t amrPatchDimension(
    const AmrPatchPayloadRecord& patch,
    std::size_t axis) {
  switch (axis) {
    case 0U: return patch.cell_dim_x;
    case 1U: return patch.cell_dim_y;
    case 2U: return patch.cell_dim_z;
    default: throw std::out_of_range("directed AMR boundary axis out of range");
  }
}

[[nodiscard]] double amrPatchExtent(
    const AmrPatchPayloadRecord& patch,
    std::size_t axis) {
  switch (axis) {
    case 0U: return patch.extent_x_comoving;
    case 1U: return patch.extent_y_comoving;
    case 2U: return patch.extent_z_comoving;
    default: throw std::out_of_range("directed AMR boundary axis out of range");
  }
}

[[nodiscard]] std::uint16_t requiredSourceBoundaryDepth(
    const AmrPatchPayloadRecord& source,
    const AmrPatchPayloadRecord& target,
    std::size_t axis) {
  const std::uint16_t source_dim = amrPatchDimension(source, axis);
  const std::uint16_t target_dim = amrPatchDimension(target, axis);
  if (source_dim == 0U || target_dim == 0U) {
    throw std::invalid_argument("directed AMR boundary depth requires positive patch dimensions");
  }
  const double source_width = amrPatchExtent(source, axis) / static_cast<double>(source_dim);
  const double target_width = amrPatchExtent(target, axis) / static_cast<double>(target_dim);
  if (!(source_width > 0.0) || !(target_width > 0.0) ||
      !std::isfinite(source_width) || !std::isfinite(target_width)) {
    throw std::invalid_argument("directed AMR boundary depth requires finite positive cell widths");
  }

  const double ratio = target_width / source_width;
  const double ratio_tolerance = 1.0e-10 * std::max(1.0, std::abs(ratio));
  if (ratio <= 1.0 + ratio_tolerance) {
    return 1U;
  }
  const double rounded = std::round(ratio);
  if (std::abs(ratio - rounded) > ratio_tolerance || rounded < 1.0 ||
      rounded > static_cast<double>(std::numeric_limits<std::uint16_t>::max())) {
    throw std::runtime_error(
        "directed AMR fine-to-coarse boundary depth is not an aligned integer cell-width ratio");
  }
  const auto depth = static_cast<std::uint16_t>(rounded);
  if (depth > source_dim) {
    throw std::runtime_error(
        "directed AMR fine-to-coarse boundary depth exceeds source patch dimension");
  }
  return depth;
}

[[nodiscard]] AmrPatchBoundaryCellRequest amrPatchBoundaryRequestForPeer(
    const AmrPatchPayloadRecord& source,
    const AmrPatchPayloadRecord& target) {
  AmrPatchBoundaryCellRequest request;
  request.patch_id = source.patch_id;
  request.boundary_face_mask = amrPatchBoundaryMaskForPeer(source, target);
  constexpr std::array<AmrPatchBoundaryFace, 6> k_faces{
      AmrPatchBoundaryFace::kXLower,
      AmrPatchBoundaryFace::kXUpper,
      AmrPatchBoundaryFace::kYLower,
      AmrPatchBoundaryFace::kYUpper,
      AmrPatchBoundaryFace::kZLower,
      AmrPatchBoundaryFace::kZUpper};
  for (const AmrPatchBoundaryFace face : k_faces) {
    if ((request.boundary_face_mask & amrBoundaryFaceBit(face)) == 0U) {
      continue;
    }
    request.boundary_face_depths[amrBoundaryFaceIndex(face)] =
        requiredSourceBoundaryDepth(source, target, amrBoundaryFaceAxis(face));
  }
  return request;
}

void mergeBoundaryRequest(
    AmrPatchBoundaryCellRequest& destination,
    const AmrPatchBoundaryCellRequest& source) {
  if (destination.patch_id == 0U) {
    destination.patch_id = source.patch_id;
  }
  if (destination.patch_id != source.patch_id) {
    throw std::invalid_argument("cannot merge directed AMR boundary requests for different patches");
  }
  destination.boundary_face_mask |= source.boundary_face_mask;
  for (std::size_t face = 0; face < destination.boundary_face_depths.size(); ++face) {
    destination.boundary_face_depths[face] = std::max(
        destination.boundary_face_depths[face], source.boundary_face_depths[face]);
  }
}

[[nodiscard]] std::uint16_t requestFaceDepth(
    const AmrPatchBoundaryCellRequest& request,
    AmrPatchBoundaryFace face) {
  const bool selected = (request.boundary_face_mask & amrBoundaryFaceBit(face)) != 0U;
  const std::uint16_t depth = request.boundary_face_depths[amrBoundaryFaceIndex(face)];
  if (selected && depth == 0U) {
    throw std::invalid_argument("directed AMR selected boundary face has zero source depth");
  }
  if (!selected && depth != 0U) {
    throw std::invalid_argument("directed AMR unselected boundary face has non-zero source depth");
  }
  return selected ? depth : 0U;
}

[[nodiscard]] std::size_t requestedBoundaryCellCount(
    const AmrPatchPayloadRecord& patch,
    const AmrPatchBoundaryCellRequest& request) {
  if (request.patch_id != patch.patch_id || request.boundary_face_mask == 0U) {
    throw std::invalid_argument("directed AMR boundary count request does not match patch metadata");
  }
  const auto selected_positions = [&request](
      std::uint16_t dim,
      AmrPatchBoundaryFace lower,
      AmrPatchBoundaryFace upper) -> std::size_t {
    const std::size_t lower_depth = requestFaceDepth(request, lower);
    const std::size_t upper_depth = requestFaceDepth(request, upper);
    if (lower_depth > dim || upper_depth > dim) {
      throw std::invalid_argument("directed AMR boundary depth exceeds patch dimension");
    }
    return std::min<std::size_t>(
        dim, core::checkedSizeAdd(lower_depth, upper_depth,
                                  "directed AMR selected boundary depth"));
  };
  const std::size_t nx = patch.cell_dim_x;
  const std::size_t ny = patch.cell_dim_y;
  const std::size_t nz = patch.cell_dim_z;
  const std::size_t total = core::checkedSizeProduct3(
      nx, ny, nz, "directed AMR boundary cell count");
  const std::size_t interior_x = nx - selected_positions(
      patch.cell_dim_x, AmrPatchBoundaryFace::kXLower, AmrPatchBoundaryFace::kXUpper);
  const std::size_t interior_y = ny - selected_positions(
      patch.cell_dim_y, AmrPatchBoundaryFace::kYLower, AmrPatchBoundaryFace::kYUpper);
  const std::size_t interior_z = nz - selected_positions(
      patch.cell_dim_z, AmrPatchBoundaryFace::kZLower, AmrPatchBoundaryFace::kZUpper);
  const std::size_t unselected = core::checkedSizeProduct3(
      interior_x, interior_y, interior_z, "directed AMR unselected interior cell count");
  return total - unselected;
}

[[nodiscard]] AmrPatchBoundaryCellRequest oneLayerBoundaryRequest(
    const AmrPatchPayloadRecord& patch,
    std::uint8_t mask) {
  AmrPatchBoundaryCellRequest request;
  request.patch_id = patch.patch_id;
  request.boundary_face_mask = mask;
  constexpr std::array<AmrPatchBoundaryFace, 6> k_faces{
      AmrPatchBoundaryFace::kXLower,
      AmrPatchBoundaryFace::kXUpper,
      AmrPatchBoundaryFace::kYLower,
      AmrPatchBoundaryFace::kYUpper,
      AmrPatchBoundaryFace::kZLower,
      AmrPatchBoundaryFace::kZUpper};
  for (const AmrPatchBoundaryFace face : k_faces) {
    if ((mask & amrBoundaryFaceBit(face)) != 0U) {
      request.boundary_face_depths[amrBoundaryFaceIndex(face)] = 1U;
    }
  }
  return request;
}

[[nodiscard]] bool patchCellOffsetMatchesBoundaryRequest(
    const AmrPatchPayloadRecord& patch,
    std::uint32_t offset,
    const AmrPatchBoundaryCellRequest& request) {
  if (request.patch_id != patch.patch_id || patch.cell_dim_x == 0U ||
      patch.cell_dim_y == 0U || patch.cell_dim_z == 0U || offset >= patch.cell_count) {
    return false;
  }
  const std::size_t nx = patch.cell_dim_x;
  const std::size_t ny = patch.cell_dim_y;
  const std::size_t nz = patch.cell_dim_z;
  const std::size_t plane = core::checkedSizeMultiply(nx, ny, "directed AMR patch plane");
  const std::size_t i = offset % nx;
  const std::size_t j = (offset / nx) % ny;
  const std::size_t k = offset / plane;
  const std::size_t x_lower = requestFaceDepth(request, AmrPatchBoundaryFace::kXLower);
  const std::size_t x_upper = requestFaceDepth(request, AmrPatchBoundaryFace::kXUpper);
  const std::size_t y_lower = requestFaceDepth(request, AmrPatchBoundaryFace::kYLower);
  const std::size_t y_upper = requestFaceDepth(request, AmrPatchBoundaryFace::kYUpper);
  const std::size_t z_lower = requestFaceDepth(request, AmrPatchBoundaryFace::kZLower);
  const std::size_t z_upper = requestFaceDepth(request, AmrPatchBoundaryFace::kZUpper);
  if (x_lower > nx || x_upper > nx || y_lower > ny || y_upper > ny ||
      z_lower > nz || z_upper > nz) {
    throw std::invalid_argument("directed AMR boundary depth exceeds patch dimension");
  }
  return (x_lower != 0U && i < x_lower) ||
      (x_upper != 0U && i >= nx - x_upper) ||
      (y_lower != 0U && j < y_lower) ||
      (y_upper != 0U && j >= ny - y_upper) ||
      (z_lower != 0U && k < z_lower) ||
      (z_upper != 0U && k >= nz - z_upper);
}

[[nodiscard]] std::uint64_t checkedU64AddLocal(
    std::uint64_t lhs,
    std::uint64_t rhs,
    std::string_view context) {
  if (rhs > std::numeric_limits<std::uint64_t>::max() - lhs) {
    throw std::overflow_error(std::string(context) + ": uint64 addition overflow");
  }
  return lhs + rhs;
}

[[nodiscard]] std::uint64_t checkedRecordTrafficBytes(
    std::size_t sent_records,
    std::size_t received_records,
    std::size_t record_bytes,
    std::string_view context) {
  const std::size_t record_count = core::checkedSizeAdd(
      sent_records, received_records, std::string(context) + " record count");
  const std::size_t bytes = core::checkedSizeMultiply(
      record_count, record_bytes, std::string(context) + " byte count");
  return core::checkedIntegralNarrow<std::uint64_t>(
      bytes, std::string(context) + " uint64 byte count");
}

[[nodiscard]] std::unordered_map<std::uint64_t, AmrPatchBoundaryCellRequest>
requestByPatchId(std::span<const AmrPatchBoundaryCellRequest> requests) {
  std::unordered_map<std::uint64_t, AmrPatchBoundaryCellRequest> result;
  result.reserve(requests.size());
  for (const AmrPatchBoundaryCellRequest& request : requests) {
    if (request.patch_id == 0U || request.boundary_face_mask == 0U) {
      throw std::invalid_argument("directed AMR boundary request is empty or malformed");
    }
    auto [it, inserted] = result.emplace(request.patch_id, request);
    if (!inserted) {
      mergeBoundaryRequest(it->second, request);
    }
  }
  return result;
}

}  // namespace

std::vector<AmrPatchBoundaryCellRequest> planDirectedAmrPatchBoundaryCellRequests(
    std::span<const AmrPatchPayloadRecord> local_patch_records,
    std::span<const AmrPatchPayloadRecord> remote_patch_records) {
  std::vector<AmrPatchBoundaryCellRequest> requests;
  requests.reserve(local_patch_records.size());
  for (const AmrPatchPayloadRecord& local_record : local_patch_records) {
    validateAmrPatchPayloadRecord(local_record);
    AmrPatchBoundaryCellRequest combined;
    combined.patch_id = local_record.patch_id;
    for (const AmrPatchPayloadRecord& remote_record : remote_patch_records) {
      validateAmrPatchPayloadRecord(remote_record);
      const AmrPatchBoundaryCellRequest one =
          amrPatchBoundaryRequestForPeer(local_record, remote_record);
      if (one.boundary_face_mask != 0U) {
        mergeBoundaryRequest(combined, one);
      }
    }
    if (combined.boundary_face_mask != 0U) {
      requests.push_back(combined);
    }
  }
  return requests;
}

std::size_t directedAmrPatchBoundaryCellCount(
    const AmrPatchPayloadRecord& patch,
    const AmrPatchBoundaryCellRequest& request) {
  validateAmrPatchPayloadRecord(patch);
  if (request.boundary_face_mask == 0U) {
    return 0U;
  }
  return requestedBoundaryCellCount(patch, request);
}

std::size_t directedAmrPatchBoundaryCellCount(
    const AmrPatchPayloadRecord& patch,
    std::uint8_t boundary_face_mask) {
  validateAmrPatchPayloadRecord(patch);
  if (boundary_face_mask == 0U) {
    return 0U;
  }
  return requestedBoundaryCellCount(
      patch, oneLayerBoundaryRequest(patch, boundary_face_mask));
}

DirectedAmrPatchPayloadExchange executeBlockingDirectedAmrPatchPayloadExchange(
    const MpiContext& mpi_context,
    std::span<const AmrPatchPayloadRecord> local_patch_records,
    const DirectedAmrPatchCellPayloadProvider& cell_payload_provider,
    const DirectedAmrPatchCellPayloadConsumer& cell_payload_consumer,
    const DirectedAmrPatchCellAdmission& cell_payload_admission,
    std::size_t transport_round_limit_bytes,
    std::uint64_t exchange_sequence) {
#if !defined(COSMOSIM_ENABLE_MPI) || !COSMOSIM_ENABLE_MPI
  (void)exchange_sequence;
#endif
  const int world_rank = mpi_context.worldRank();
  std::unordered_map<std::uint64_t, AmrPatchPayloadRecord> local_patch_by_id;
  local_patch_by_id.reserve(local_patch_records.size());
  for (const AmrPatchPayloadRecord& record : local_patch_records) {
    validateAmrPatchPayloadRecord(record);
    if (record.owner_rank != world_rank) {
      throw std::invalid_argument("directed AMR patch exchange received non-local authoritative patch metadata");
    }
    const auto [it, inserted] = local_patch_by_id.emplace(record.patch_id, record);
    if (!inserted) {
      throw std::invalid_argument("directed AMR patch exchange found duplicate local patch metadata");
    }
  }
  if ((!cell_payload_provider || !cell_payload_consumer) && !local_patch_records.empty()) {
    throw std::invalid_argument(
        "directed AMR patch exchange requires boundary-cell producer and consumer callbacks");
  }

  DirectedAmrPatchPayloadExchange result;
  if (!mpi_context.isEnabled() || mpi_context.worldSize() <= 1) {
    return result;
  }

#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  constexpr int k_patch_count_tag_base = 8910;
  constexpr int k_patch_payload_tag_base = 9910;
  constexpr int k_cell_count_tag_base = 10910;
  constexpr int k_cell_payload_tag_base = 11910;
  const std::vector<int> candidate_peers = discoverCandidateAmrPeers(
      mpi_context, local_patch_records, &result.diagnostics);
  for (const int peer_rank : candidate_peers) {
    std::vector<AmrPatchPayloadRecord> remote_patches = exchangePodRecordsWithPeer(
        mpi_context,
        peer_rank,
        local_patch_records,
        k_patch_count_tag_base,
        k_patch_payload_tag_base,
        exchange_sequence,
        "executeBlockingDirectedAmrPatchPayloadExchange/patch");
    for (const AmrPatchPayloadRecord& record : remote_patches) {
      validateAmrPatchPayloadRecord(record);
      if (record.owner_rank != peer_rank) {
        throw std::runtime_error("directed AMR patch exchange received metadata from a rank that is not the owner");
      }
    }
    result.diagnostics.directed_patch_descriptor_records_sent +=
        static_cast<std::uint64_t>(local_patch_records.size());
    result.diagnostics.directed_patch_descriptor_records_received +=
        static_cast<std::uint64_t>(remote_patches.size());
    result.diagnostics.patch_descriptor_bytes = checkedU64AddLocal(
        result.diagnostics.patch_descriptor_bytes,
        checkedRecordTrafficBytes(
            local_patch_records.size(), remote_patches.size(),
            sizeof(AmrPatchPayloadRecord), "directed AMR patch descriptor traffic"),
        "directed AMR accumulated patch descriptor traffic");

    const std::vector<AmrPatchBoundaryCellRequest> local_requests =
        planDirectedAmrPatchBoundaryCellRequests(local_patch_records, remote_patches);
    const std::vector<AmrPatchBoundaryCellRequest> remote_requests =
        planDirectedAmrPatchBoundaryCellRequests(remote_patches, local_patch_records);
    if (local_requests.empty() && remote_requests.empty()) {
      continue;
    }
    ++result.diagnostics.neighbor_peer_count;

    const auto local_request_by_id = requestByPatchId(local_requests);
    const auto remote_request_by_id = requestByPatchId(remote_requests);
    std::uint64_t peer_interface_count = 0U;
    for (const AmrPatchPayloadRecord& local_record : local_patch_records) {
      for (const AmrPatchPayloadRecord& remote_record : remote_patches) {
        if (amrPatchBoundaryMaskForPeer(local_record, remote_record) != 0U) {
          peer_interface_count = checkedU64AddLocal(
              peer_interface_count, 1U, "directed AMR peer interface count");
        }
      }
    }
    result.diagnostics.remote_interface_count = checkedU64AddLocal(
        result.diagnostics.remote_interface_count,
        peer_interface_count,
        "directed AMR accumulated remote interface count");

    for (const AmrPatchPayloadRecord& remote_patch : remote_patches) {
      if (remote_request_by_id.contains(remote_patch.patch_id)) {
        result.patch_payloads_received.push_back(remote_patch);
      }
    }

    std::uint64_t logical_send_count = 0U;
    for (const AmrPatchBoundaryCellRequest& request : local_requests) {
      const AmrPatchPayloadRecord& patch = local_patch_by_id.at(request.patch_id);
      logical_send_count = checkedU64AddLocal(
          logical_send_count,
          static_cast<std::uint64_t>(requestedBoundaryCellCount(
              patch, request)),
          "directed AMR logical send boundary-cell count");
    }

    std::unordered_map<std::uint64_t, AmrPatchPayloadRecord> remote_patch_by_id;
    remote_patch_by_id.reserve(remote_patches.size());
    for (const AmrPatchPayloadRecord& patch : remote_patches) {
      remote_patch_by_id.emplace(patch.patch_id, patch);
    }
    std::unordered_map<std::uint64_t, std::size_t> observed_local_count;
    observed_local_count.reserve(local_requests.size());
    std::unordered_map<std::uint64_t, std::size_t> observed_remote_count;
    observed_remote_count.reserve(remote_requests.size());

    const DirectedAmrPatchCellPayloadProvider validated_provider =
        [&](std::span<const AmrPatchBoundaryCellRequest> requests,
            std::uint64_t first_record,
            std::size_t max_records,
            std::vector<AmrPatchCellPayloadRecord>& output) {
          cell_payload_provider(requests, first_record, max_records, output);
          for (const AmrPatchCellPayloadRecord& record : output) {
            validateAmrPatchCellPayloadRecord(record);
            if (record.owner_rank != world_rank) {
              throw std::invalid_argument(
                  "directed AMR boundary payload producer returned a non-local cell");
            }
            const auto patch_it = local_patch_by_id.find(record.patch_id);
            const auto request_it = local_request_by_id.find(record.patch_id);
            if (patch_it == local_patch_by_id.end() || request_it == local_request_by_id.end() ||
                !patchCellOffsetMatchesBoundaryRequest(
                    patch_it->second, record.local_cell_offset, request_it->second)) {
              throw std::invalid_argument(
                  "directed AMR boundary payload producer returned a non-interface cell");
            }
            ++observed_local_count[record.patch_id];
          }
        };
    const DirectedAmrPatchCellPayloadConsumer validated_consumer =
        [&](int source_rank, std::span<const AmrPatchCellPayloadRecord> records) {
          if (source_rank != peer_rank) {
            throw std::runtime_error("directed AMR patch-cell consumer received wrong peer rank");
          }
          for (const AmrPatchCellPayloadRecord& record : records) {
            validateAmrPatchCellPayloadRecord(record);
            if (record.owner_rank != peer_rank) {
              throw std::runtime_error(
                  "directed AMR patch-cell stream received payload from a rank that is not the owner");
            }
            const auto patch_it = remote_patch_by_id.find(record.patch_id);
            const auto request_it = remote_request_by_id.find(record.patch_id);
            if (patch_it == remote_patch_by_id.end() || request_it == remote_request_by_id.end() ||
                !patchCellOffsetMatchesBoundaryRequest(
                    patch_it->second, record.local_cell_offset, request_it->second)) {
              throw std::runtime_error(
                  "directed AMR patch-cell stream received a non-interface cell payload");
            }
            ++observed_remote_count[record.patch_id];
          }
          cell_payload_consumer(source_rank, records);
        };

    exchangeAmrPatchCellRecordStreamWithPeer(
        mpi_context,
        peer_rank,
        logical_send_count,
        local_requests,
        remote_patches,
        remote_requests,
        validated_provider,
        validated_consumer,
        cell_payload_admission,
        k_cell_count_tag_base,
        k_cell_payload_tag_base,
        exchange_sequence,
        transport_round_limit_bytes,
        result.diagnostics);

    for (const AmrPatchBoundaryCellRequest& request : local_requests) {
      const AmrPatchPayloadRecord& patch = local_patch_by_id.at(request.patch_id);
      if (observed_local_count[request.patch_id] !=
          requestedBoundaryCellCount(patch, request)) {
        throw std::runtime_error(
            "directed AMR streamed producer did not cover requested local boundary cells");
      }
    }
    std::uint64_t logical_receive_count = 0U;
    for (const AmrPatchBoundaryCellRequest& request : remote_requests) {
      const AmrPatchPayloadRecord& patch = remote_patch_by_id.at(request.patch_id);
      const std::size_t expected = requestedBoundaryCellCount(
          patch, request);
      if (observed_remote_count[request.patch_id] != expected) {
        throw std::runtime_error(
            "directed AMR streamed remote boundary coverage mismatch");
      }
      logical_receive_count = checkedU64AddLocal(
          logical_receive_count,
          static_cast<std::uint64_t>(expected),
          "directed AMR logical receive boundary-cell count");
    }
    result.diagnostics.directed_patch_cell_records_sent = checkedU64AddLocal(
        result.diagnostics.directed_patch_cell_records_sent,
        logical_send_count,
        "directed AMR accumulated sent boundary-cell records");
    result.diagnostics.directed_patch_cell_records_received = checkedU64AddLocal(
        result.diagnostics.directed_patch_cell_records_received,
        logical_receive_count,
        "directed AMR accumulated received boundary-cell records");
    const std::uint64_t traffic_records = checkedU64AddLocal(
        logical_send_count,
        logical_receive_count,
        "directed AMR streamed traffic record count");
    const std::uint64_t traffic_bytes = core::checkedMemoryBytesAdd(
        0U,
        core::checkedIntegralNarrow<std::uint64_t>(
            core::checkedSizeMultiply(
                core::checkedIntegralNarrow<std::size_t>(
                    traffic_records,
                    "directed AMR streamed traffic records size_t"),
                sizeof(AmrPatchCellPayloadRecord),
                "directed AMR streamed traffic bytes"),
            "directed AMR streamed traffic bytes uint64"),
        "directed AMR streamed traffic byte accounting");
    result.diagnostics.patch_cell_payload_bytes = checkedU64AddLocal(
        result.diagnostics.patch_cell_payload_bytes,
        traffic_bytes,
        "directed AMR accumulated streamed patch-cell traffic");
  }
  result.diagnostics.remote_patch_ghost_count =
      static_cast<std::uint64_t>(result.patch_payloads_received.size());
  return result;
#else
  throw std::runtime_error("directed AMR patch payload exchange requires MPI support when MPI context is enabled");
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
      throw std::invalid_argument("AMR patch payload compatibility exchange received non-local authoritative patch metadata");
    }
  }
  return std::vector<AmrPatchPayloadRecord>(local_records.begin(), local_records.end());
}

std::vector<AmrPatchCellPayloadRecord> executeBlockingAmrPatchCellPayloadExchange(
    const MpiContext& mpi_context,
    std::span<const AmrPatchCellPayloadRecord> local_records,
    std::uint64_t exchange_sequence) {
  for (const AmrPatchCellPayloadRecord& record : local_records) {
    validateAmrPatchCellPayloadRecord(record);
    if (record.owner_rank != mpi_context.worldRank()) {
      throw std::invalid_argument("directed AMR patch-cell compatibility exchange received non-local authoritative cell payload");
    }
  }
  (void)exchange_sequence;
  return std::vector<AmrPatchCellPayloadRecord>(local_records.begin(), local_records.end());
}

std::vector<AmrFluxRegisterPayloadRecord> executeBlockingAmrFluxRegisterPayloadExchange(
    const MpiContext& mpi_context,
    std::span<const AmrFluxRegisterPayloadRecord> local_records,
    std::uint64_t exchange_sequence) {
#if !defined(COSMOSIM_ENABLE_MPI) || !COSMOSIM_ENABLE_MPI
  (void)exchange_sequence;
#endif
  const int world_rank = mpi_context.worldRank();
  const int world_size = mpi_context.worldSize();
  for (const AmrFluxRegisterPayloadRecord& record : local_records) {
    validateAmrFluxRegisterPayloadRecord(record);
    if (record.source_rank != world_rank) {
      throw std::invalid_argument("AMR flux-register payload source rank does not match MPI context");
    }
    if (record.owner_rank < 0 || record.owner_rank >= world_size) {
      throw std::invalid_argument("AMR flux-register payload owner rank is outside MPI world");
    }
  }
  if (!mpi_context.isEnabled() || world_size <= 1) {
    return std::vector<AmrFluxRegisterPayloadRecord>(local_records.begin(), local_records.end());
  }
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  constexpr int k_flux_count_tag_base = 12910;
  constexpr int k_flux_payload_tag_base = 13910;
  std::vector<AmrFluxRegisterPayloadRecord> inbound_records;
  for (int peer_rank = 0; peer_rank < world_size; ++peer_rank) {
    if (peer_rank == world_rank) {
      continue;
    }
    std::vector<AmrFluxRegisterPayloadRecord> outbound_to_peer;
    outbound_to_peer.reserve(local_records.size());
    for (const AmrFluxRegisterPayloadRecord& record : local_records) {
      if (record.owner_rank == peer_rank) {
        outbound_to_peer.push_back(record);
      }
    }
    std::vector<AmrFluxRegisterPayloadRecord> received_from_peer = exchangePodRecordsWithPeer(
        mpi_context,
        peer_rank,
        std::span<const AmrFluxRegisterPayloadRecord>(outbound_to_peer),
        k_flux_count_tag_base,
        k_flux_payload_tag_base,
        exchange_sequence,
        "executeBlockingAmrFluxRegisterPayloadExchange");
    for (const AmrFluxRegisterPayloadRecord& record : received_from_peer) {
      validateAmrFluxRegisterPayloadRecord(record);
      if (record.owner_rank != world_rank || record.source_rank != peer_rank) {
        throw std::runtime_error("AMR flux-register directed exchange returned stale source/owner rank metadata");
      }
    }
    inbound_records.insert(inbound_records.end(), received_from_peer.begin(), received_from_peer.end());
  }
  for (const AmrFluxRegisterPayloadRecord& record : local_records) {
    if (record.owner_rank == world_rank) {
      inbound_records.push_back(record);
    }
  }
  return inbound_records;
#else
  throw std::runtime_error("AMR flux-register payload exchange requires MPI support when MPI context is enabled");
#endif
}

namespace {

template <typename Record>
[[nodiscard]] std::vector<Record> allgatherTrivialRecordsBounded(
    const MpiContext& mpi_context,
    std::span<const Record> local_records,
    std::string_view phase) {
  static_assert(std::is_trivially_copyable_v<Record>);

  std::size_t local_byte_count = 0U;
  std::exception_ptr local_preparation_failure;
  try {
    local_byte_count = core::checkedSizeMultiply(
        local_records.size(), sizeof(Record),
        std::string(phase) + " local byte count");
  } catch (...) {
    local_preparation_failure = std::current_exception();
  }
  mpi_context.rethrowCollectivePreparationFailure(
      local_preparation_failure,
      std::string(phase) + " local wire preparation");

  const auto local_wire = std::span<const std::uint8_t>(
      reinterpret_cast<const std::uint8_t*>(local_records.data()),
      local_byte_count);
  std::vector<std::uint8_t> gathered_wire =
      mpi_context.allgatherBytesBounded(local_wire);

  std::vector<Record> result;
  std::exception_ptr local_reassembly_failure;
  try {
    if (gathered_wire.size() % sizeof(Record) != 0U) {
      throw std::runtime_error(
          std::string(phase) + " returned partial record bytes");
    }
    result.resize(gathered_wire.size() / sizeof(Record));
    if (!gathered_wire.empty()) {
      std::memcpy(result.data(), gathered_wire.data(), gathered_wire.size());
    }
  } catch (...) {
    local_reassembly_failure = std::current_exception();
  }
  mpi_context.rethrowCollectivePreparationFailure(
      local_reassembly_failure,
      std::string(phase) + " receive reassembly");
  return result;
}

}  // namespace

std::vector<HydroConservativeFluxCorrectionRecord> executeBlockingHydroConservativeFluxCorrectionExchange(
    const MpiContext& mpi_context,
    std::span<const HydroConservativeFluxCorrectionRecord> local_records,
    std::uint64_t exchange_sequence) {
  (void)exchange_sequence;
  std::exception_ptr local_validation_failure;
  try {
    for (const HydroConservativeFluxCorrectionRecord& record : local_records) {
      validateHydroConservativeFluxCorrectionRecord(record);
      if (record.source_rank != mpi_context.worldRank()) {
        throw std::invalid_argument(
            "hydro conservative flux correction source rank does not match MPI context");
      }
    }
  } catch (...) {
    local_validation_failure = std::current_exception();
  }
  if (mpi_context.isEnabled()) {
    mpi_context.rethrowCollectivePreparationFailure(
        local_validation_failure,
        "hydro conservative flux correction local validation");
  } else if (local_validation_failure != nullptr) {
    std::rethrow_exception(local_validation_failure);
  }
  if (!mpi_context.isEnabled()) {
    return std::vector<HydroConservativeFluxCorrectionRecord>(
        local_records.begin(), local_records.end());
  }
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  const int world_size = mpi_context.worldSize();
  std::vector<HydroConservativeFluxCorrectionRecord> result =
      allgatherTrivialRecordsBounded(
          mpi_context, local_records,
          "hydro conservative flux correction exchange");

  std::exception_ptr local_result_failure;
  try {
    for (const HydroConservativeFluxCorrectionRecord& record : result) {
      validateHydroConservativeFluxCorrectionRecord(record);
      if (record.source_rank < 0 || record.source_rank >= world_size ||
          record.owner_rank < 0 || record.owner_rank >= world_size) {
        throw std::runtime_error(
            "hydro conservative flux correction exchange returned invalid rank metadata");
      }
    }
  } catch (...) {
    local_result_failure = std::current_exception();
  }
  mpi_context.rethrowCollectivePreparationFailure(
      local_result_failure,
      "hydro conservative flux correction receive validation");
  return result;
#else
  throw std::runtime_error(
      "hydro conservative flux correction exchange requires MPI support when MPI context is enabled");
#endif
}

std::vector<HydroGhostCellRequest> executeBlockingHydroGhostCellRequestExchange(
    const MpiContext& mpi_context,
    std::span<const HydroGhostCellRequest> local_requests,
    std::uint64_t exchange_sequence) {
  (void)exchange_sequence;
  std::exception_ptr local_validation_failure;
  try {
    for (const HydroGhostCellRequest& request : local_requests) {
      validateHydroGhostCellRequest(request);
      if (request.descriptor.consumer_rank != mpi_context.worldRank()) {
        throw std::invalid_argument(
            "hydro ghost cell request consumer rank does not match MPI context");
      }
    }
  } catch (...) {
    local_validation_failure = std::current_exception();
  }
  if (mpi_context.isEnabled()) {
    mpi_context.rethrowCollectivePreparationFailure(
        local_validation_failure,
        "hydro ghost cell request local validation");
  } else if (local_validation_failure != nullptr) {
    std::rethrow_exception(local_validation_failure);
  }
  if (!mpi_context.isEnabled()) {
    return std::vector<HydroGhostCellRequest>(
        local_requests.begin(), local_requests.end());
  }
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  const int world_size = mpi_context.worldSize();
  std::vector<HydroGhostCellRequest> result =
      allgatherTrivialRecordsBounded(
          mpi_context, local_requests,
          "hydro ghost cell request exchange");

  std::exception_ptr local_result_failure;
  try {
    for (const HydroGhostCellRequest& request : result) {
      validateHydroGhostCellRequest(request);
      if (request.descriptor.owner_rank < 0 ||
          request.descriptor.owner_rank >= world_size ||
          request.descriptor.consumer_rank < 0 ||
          request.descriptor.consumer_rank >= world_size) {
        throw std::runtime_error(
            "hydro ghost cell request exchange returned invalid rank metadata");
      }
    }
  } catch (...) {
    local_result_failure = std::current_exception();
  }
  mpi_context.rethrowCollectivePreparationFailure(
      local_result_failure,
      "hydro ghost cell request receive validation");
  return result;
#else
  throw std::runtime_error(
      "hydro ghost cell request exchange requires MPI support when MPI context is enabled");
#endif
}

std::vector<HydroGhostCellPayloadRecord> executeBlockingHydroGhostCellPayloadExchange(
    const MpiContext& mpi_context,
    std::span<const HydroGhostCellPayloadRecord> local_records,
    std::uint64_t exchange_sequence) {
  (void)exchange_sequence;
  std::exception_ptr local_validation_failure;
  try {
    for (const HydroGhostCellPayloadRecord& record : local_records) {
      validateHydroGhostCellPayloadRecord(record);
      if (record.descriptor.owner_rank != mpi_context.worldRank()) {
        throw std::invalid_argument(
            "hydro ghost cell payload owner rank does not match MPI context");
      }
    }
  } catch (...) {
    local_validation_failure = std::current_exception();
  }
  if (mpi_context.isEnabled()) {
    mpi_context.rethrowCollectivePreparationFailure(
        local_validation_failure,
        "hydro ghost cell payload local validation");
  } else if (local_validation_failure != nullptr) {
    std::rethrow_exception(local_validation_failure);
  }
  if (!mpi_context.isEnabled()) {
    return std::vector<HydroGhostCellPayloadRecord>(
        local_records.begin(), local_records.end());
  }
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  const int world_size = mpi_context.worldSize();
  std::vector<HydroGhostCellPayloadRecord> result =
      allgatherTrivialRecordsBounded(
          mpi_context, local_records,
          "hydro ghost cell payload exchange");

  std::exception_ptr local_result_failure;
  try {
    for (const HydroGhostCellPayloadRecord& record : result) {
      validateHydroGhostCellPayloadRecord(record);
      if (record.descriptor.owner_rank < 0 ||
          record.descriptor.owner_rank >= world_size ||
          record.descriptor.consumer_rank < 0 ||
          record.descriptor.consumer_rank >= world_size) {
        throw std::runtime_error(
            "hydro ghost cell payload exchange returned invalid rank metadata");
      }
    }
  } catch (...) {
    local_result_failure = std::current_exception();
  }
  mpi_context.rethrowCollectivePreparationFailure(
      local_result_failure,
      "hydro ghost cell payload receive validation");
  return result;
#else
  throw std::runtime_error(
      "hydro ghost cell payload exchange requires MPI support when MPI context is enabled");
#endif
}

PmSlabHaloExchangeResult executeBlockingPmSlabHaloExchange(
    const MpiContext& mpi_context,
    const PmSlabLayout& layout,
    std::span<const double> local_scalar_field,
    std::size_t halo_depth_x,
    bool periodic_x,
    std::uint64_t exchange_sequence) {
#if !defined(COSMOSIM_ENABLE_MPI) || !COSMOSIM_ENABLE_MPI
  (void)exchange_sequence;
#endif
  PmSlabHaloExchangeResult result;
  std::vector<double> send_left;
  std::vector<double> send_right;
  std::size_t halo_value_count = 0U;
  std::uint64_t payload_bytes = 0U;
  int left_peer = -1;
  int right_peer = -1;
  bool no_exchange = false;
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  int communicator_world_size = 1;
  int communicator_world_rank = 0;
  const bool communicator_mpi_active =
      queryActiveMpiWorld(communicator_world_size, communicator_world_rank);
#endif

  std::exception_ptr local_preparation_failure;
  try {
    if (!layout.isValid()) {
      throw std::invalid_argument("PM slab halo exchange requires a valid slab layout");
    }
    if (layout.world_size != mpi_context.worldSize() ||
        layout.world_rank != mpi_context.worldRank()) {
      throw std::invalid_argument(
          "PM slab halo exchange layout world metadata must match MPI context");
    }
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
    if (!communicator_mpi_active &&
        (layout.world_size > 1 || mpi_context.isEnabled())) {
      throw std::invalid_argument(
          "PM slab halo exchange requires an active MPI_COMM_WORLD for an enabled or distributed context");
    }
    if (communicator_mpi_active &&
        (layout.world_size != communicator_world_size ||
         layout.world_rank != communicator_world_rank)) {
      throw std::invalid_argument(
          "PM slab halo exchange layout world metadata must match MPI_COMM_WORLD");
    }
#endif
    if (layout.global_ny >
        std::numeric_limits<std::size_t>::max() / layout.global_nz) {
      throw std::overflow_error("PM slab halo exchange plane size overflows size_t");
    }
    const std::size_t plane_size = layout.global_ny * layout.global_nz;
    if (layout.local_nx() >
        std::numeric_limits<std::size_t>::max() / plane_size) {
      throw std::overflow_error("PM slab halo exchange local field size overflows size_t");
    }
    const std::size_t expected_local_values = layout.local_nx() * plane_size;
    if (local_scalar_field.size() != expected_local_values) {
      throw std::invalid_argument(
          "PM slab halo exchange field size does not match local slab cell count");
    }
    if (halo_depth_x == 0 || layout.world_size == 1 || layout.local_nx() == 0) {
      no_exchange = true;
    } else {
      if (!mpi_context.isEnabled()) {
        throw std::runtime_error(
            "PM slab halo exchange requires MPI for distributed layouts");
      }
      // A single-neighbor payload cannot span a second slab. Use the smallest
      // non-empty slab extent so every communicating rank posts matching counts.
      std::size_t minimum_nonempty_slab_nx =
          std::numeric_limits<std::size_t>::max();
      for (int rank = 0; rank < layout.world_size; ++rank) {
        const PmSlabRange owned =
            pmOwnedXRangeForRank(layout.global_nx, layout.world_size, rank);
        if (owned.extentX() > 0U) {
          minimum_nonempty_slab_nx =
              std::min(minimum_nonempty_slab_nx, owned.extentX());
        }
      }
      if (minimum_nonempty_slab_nx ==
          std::numeric_limits<std::size_t>::max()) {
        throw std::logic_error("PM slab halo exchange layout has no non-empty owner");
      }
      const std::size_t depth =
          std::min(halo_depth_x, minimum_nonempty_slab_nx);
      if (depth > std::numeric_limits<std::size_t>::max() / plane_size) {
        throw std::overflow_error("PM slab halo exchange payload size overflows size_t");
      }
      halo_value_count = depth * plane_size;
      if (halo_value_count >
          static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw std::overflow_error(
            "PM slab halo exchange payload count exceeds MPI int limit");
      }
      if (halo_value_count >
          std::numeric_limits<std::uint64_t>::max() / sizeof(double)) {
        throw std::overflow_error(
            "PM slab halo exchange byte diagnostics overflow uint64_t");
      }
      payload_bytes =
          static_cast<std::uint64_t>(halo_value_count) * sizeof(double);
      result.halo_depth_x = depth;
      if (layout.owned_x.begin_x > 0) {
        left_peer = pmOwnerRankForGlobalX(
            layout.global_nx,
            layout.world_size,
            layout.owned_x.begin_x - 1U);
      } else if (periodic_x) {
        left_peer = pmOwnerRankForGlobalX(
            layout.global_nx,
            layout.world_size,
            layout.global_nx - 1U);
      }
      if (layout.owned_x.end_x < layout.global_nx) {
        right_peer = pmOwnerRankForGlobalX(
            layout.global_nx,
            layout.world_size,
            layout.owned_x.end_x);
      } else if (periodic_x) {
        right_peer = pmOwnerRankForGlobalX(
            layout.global_nx, layout.world_size, 0U);
      }
      result.left_peer_rank = left_peer;
      result.right_peer_rank = right_peer;
      result.left_halo.assign(halo_value_count, 0.0);
      result.right_halo.assign(halo_value_count, 0.0);
      send_left.assign(halo_value_count, 0.0);
      send_right.assign(halo_value_count, 0.0);
      const std::span<const double> left_source =
          local_scalar_field.first(halo_value_count);
      const std::span<const double> right_source =
          local_scalar_field.last(halo_value_count);
      std::copy(left_source.begin(), left_source.end(), send_left.begin());
      std::copy(right_source.begin(), right_source.end(), send_right.begin());

      const int local_rank = mpi_context.worldRank();
      const auto is_remote_peer = [&](int peer) {
        return peer >= 0 && peer != local_rank;
      };
      const std::uint64_t remote_side_count =
          static_cast<std::uint64_t>(is_remote_peer(left_peer)) +
          static_cast<std::uint64_t>(is_remote_peer(right_peer));
      if (remote_side_count > 0 &&
          payload_bytes >
              std::numeric_limits<std::uint64_t>::max() / remote_side_count) {
        throw std::overflow_error(
            "PM slab halo exchange aggregate byte diagnostics overflow uint64_t");
      }
      result.sent_bytes = payload_bytes * remote_side_count;
      result.received_bytes = payload_bytes * remote_side_count;
    }
  } catch (...) {
    local_preparation_failure = std::current_exception();
  }

#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  if (communicator_mpi_active && communicator_world_size > 1) {
    const std::uint64_t local_failure_vote =
        local_preparation_failure ? 1U : 0U;
    std::uint64_t failure_count = 0U;
    MPI_Allreduce(
        &local_failure_vote,
        &failure_count,
        1,
        MPI_UINT64_T,
        MPI_SUM,
        MPI_COMM_WORLD);
    if (failure_count != 0U) {
      if (local_preparation_failure) {
        std::rethrow_exception(local_preparation_failure);
      }
      throw std::runtime_error(
          "PM slab halo exchange peer rejected protocol preparation");
    }
    const std::array<std::uint64_t, 7> local_protocol_identity{
        static_cast<std::uint64_t>(layout.global_nx),
        static_cast<std::uint64_t>(layout.global_ny),
        static_cast<std::uint64_t>(layout.global_nz),
        static_cast<std::uint64_t>(halo_depth_x),
        periodic_x ? 1U : 0U,
        exchange_sequence,
        static_cast<std::uint64_t>(communicator_world_size),
    };
    std::array<std::uint64_t, 7> minimum_protocol_identity{};
    std::array<std::uint64_t, 7> maximum_protocol_identity{};
    MPI_Allreduce(
        local_protocol_identity.data(),
        minimum_protocol_identity.data(),
        static_cast<int>(local_protocol_identity.size()),
        MPI_UINT64_T,
        MPI_MIN,
        MPI_COMM_WORLD);
    MPI_Allreduce(
        local_protocol_identity.data(),
        maximum_protocol_identity.data(),
        static_cast<int>(local_protocol_identity.size()),
        MPI_UINT64_T,
        MPI_MAX,
        MPI_COMM_WORLD);
    if (minimum_protocol_identity != maximum_protocol_identity) {
      throw std::runtime_error(
          "PM slab halo exchange ranks disagree on global shape, halo depth, "
          "boundary mode, or exchange sequence");
    }
  }
#endif
  if (local_preparation_failure) {
    std::rethrow_exception(local_preparation_failure);
  }
  if (no_exchange) {
    return result;
  }

#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
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

  const int local_rank = mpi_context.worldRank();
  const auto is_remote_peer = [&](int peer) {
    return peer >= 0 && peer != local_rank;
  };
  if (left_peer == local_rank) {
    std::copy(send_right.begin(), send_right.end(), result.left_halo.begin());
  }
  if (right_peer == local_rank) {
    std::copy(send_left.begin(), send_left.end(), result.right_halo.begin());
  }

  std::array<MPI_Request, 4> requests{};
  int request_count = 0;
  const int mpi_value_count = static_cast<int>(halo_value_count);
  const auto post_receive = [&](int peer, std::vector<double>& receive, int sender_side) {
    if (!is_remote_peer(peer)) {
      return;
    }
    MPI_Irecv(
        receive.data(),
        mpi_value_count,
        MPI_DOUBLE,
        peer,
        ghostExchangeSequencedTag(side_tag(peer, sender_side), local_rank, peer, exchange_sequence),
        MPI_COMM_WORLD,
        &requests[static_cast<std::size_t>(request_count++)]);
  };
  const auto post_send = [&](int peer, const std::vector<double>& send, int sender_side) {
    if (!is_remote_peer(peer)) {
      return;
    }
    MPI_Isend(
        const_cast<double*>(send.data()),
        mpi_value_count,
        MPI_DOUBLE,
        peer,
        ghostExchangeSequencedTag(side_tag(peer, sender_side), local_rank, peer, exchange_sequence),
        MPI_COMM_WORLD,
        &requests[static_cast<std::size_t>(request_count++)]);
  };

  // Post both receives before either send. Distinct side tags keep two-sided
  // traffic unambiguous when the periodic left and right owner are the same rank.
  post_receive(left_peer, result.left_halo, k_send_right_side);
  post_receive(right_peer, result.right_halo, k_send_left_side);
  post_send(left_peer, send_left, k_send_left_side);
  post_send(right_peer, send_right, k_send_right_side);
  if (request_count > 0) {
    MPI_Waitall(request_count, requests.data(), MPI_STATUSES_IGNORE);
  }
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
      commitOptionalGhostLane(&ghost_storage.position_x_comoving, result.received_ghosts.position_x_comoving, storage_size, local_index, result_row);
      commitOptionalGhostLane(&ghost_storage.position_y_comoving, result.received_ghosts.position_y_comoving, storage_size, local_index, result_row);
      commitOptionalGhostLane(&ghost_storage.position_z_comoving, result.received_ghosts.position_z_comoving, storage_size, local_index, result_row);
      commitOptionalGhostLane(&ghost_storage.mass_code, result.received_ghosts.mass_code, storage_size, local_index, result_row);
      commitOptionalGhostLane(&ghost_storage.density_code, result.received_ghosts.density_code, storage_size, local_index, result_row);
      commitOptionalGhostLane(&ghost_storage.velocity_x_code, result.received_ghosts.velocity_x_code, storage_size, local_index, result_row);
      commitOptionalGhostLane(&ghost_storage.velocity_y_code, result.received_ghosts.velocity_y_code, storage_size, local_index, result_row);
      commitOptionalGhostLane(&ghost_storage.velocity_z_code, result.received_ghosts.velocity_z_code, storage_size, local_index, result_row);
      commitOptionalGhostLane(&ghost_storage.pressure_code, result.received_ghosts.pressure_code, storage_size, local_index, result_row);
      commitOptionalGhostLane(&ghost_storage.internal_energy_code, result.received_ghosts.internal_energy_code, storage_size, local_index, result_row);
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
  const int world_rank = mpi_context.worldRank();
  const int world_size = mpi_context.worldSize();
  std::unordered_map<std::uint64_t, std::uint32_t> owned_index_by_particle_id;
  std::vector<std::vector<std::uint64_t>> requested_particle_ids_by_rank;
  std::vector<std::vector<std::uint32_t>> recv_indices_by_rank;
  std::exception_ptr local_preparation_failure;
  try {
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
    if (world_rank < 0 || world_size <= 0 || world_rank >= world_size) {
      throw std::invalid_argument(
          "executeBlockingGhostRefreshExchangeFromDescriptors: invalid MPI context");
    }

    owned_index_by_particle_id.reserve(local_ghost_descriptors.size());
    requested_particle_ids_by_rank.resize(static_cast<std::size_t>(world_size));
    recv_indices_by_rank.resize(static_cast<std::size_t>(world_size));
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
        requested_particle_ids_by_rank[static_cast<std::size_t>(descriptor.owning_rank)].push_back(
            descriptor.particle_id);
        recv_indices_by_rank[static_cast<std::size_t>(descriptor.owning_rank)].push_back(local_index);
      }
    }
  } catch (...) {
    local_preparation_failure = std::current_exception();
  }
  mpi_context.rethrowCollectivePreparationFailure(
      local_preparation_failure,
      "ghost refresh descriptor/request preparation");
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
        .result = BlockingGhostExchangeResult{.received_ghosts = GhostExchangeBufferSoA{
            .epoch = expected_epoch,
            .entity_id = {},
            .position_x_comoving = {},
            .position_y_comoving = {},
            .position_z_comoving = {},
            .mass_code = {},
            .density_code = {},
            .velocity_x_code = {},
            .velocity_y_code = {},
            .velocity_z_code = {},
            .pressure_code = {},
            .internal_energy_code = {},
        }},
    };
  }

#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  constexpr int k_request_count_tag_base = 4810;
  constexpr int k_request_ready_tag_base = 5310;
  constexpr int k_request_payload_tag_base = 5810;
  constexpr int k_request_mapping_ready_tag_base = 6310;
  std::vector<std::vector<std::uint32_t>> send_indices_by_rank;
  local_preparation_failure = nullptr;
  try {
    send_indices_by_rank.resize(static_cast<std::size_t>(world_size));
  } catch (...) {
    local_preparation_failure = std::current_exception();
  }
  mpi_context.rethrowCollectivePreparationFailure(
      local_preparation_failure,
      "ghost refresh request exchange storage preparation");

  std::exception_ptr first_protocol_failure;
  const auto record_peer_failure = [&](int peer_rank, std::string_view phase) {
    if (first_protocol_failure != nullptr) {
      return;
    }
    try {
      throw std::runtime_error(
          "ghost refresh peer " + std::to_string(peer_rank) +
          " rejected " + std::string(phase));
    } catch (...) {
      first_protocol_failure = std::current_exception();
    }
  };

  for (int peer_rank = 0; peer_rank < world_size; ++peer_rank) {
    if (peer_rank == world_rank) {
      continue;
    }
    const auto& request_ids =
        requested_particle_ids_by_rank[static_cast<std::size_t>(peer_rank)];
    const std::uint64_t send_count = core::checkedIntegralNarrow<std::uint64_t>(
        request_ids.size(), "ghost refresh request send count");
    std::uint64_t recv_count = 0U;
    if (MPI_Sendrecv(
            const_cast<std::uint64_t*>(&send_count), 1, MPI_UINT64_T,
            peer_rank,
            ghostExchangeSequencedTag(
                k_request_count_tag_base, world_rank, peer_rank,
                expected_epoch.ghost_sync_epoch),
            &recv_count, 1, MPI_UINT64_T,
            peer_rank,
            ghostExchangeSequencedTag(
                k_request_count_tag_base, world_rank, peer_rank,
                expected_epoch.ghost_sync_epoch),
            MPI_COMM_WORLD, MPI_STATUS_IGNORE) != MPI_SUCCESS) {
      throw std::runtime_error(
          "executeBlockingGhostRefreshExchangeFromDescriptors: request-count Sendrecv failed");
    }

    std::vector<std::uint64_t> received_request_ids;
    int request_send_count = 0;
    int request_recv_count = 0;
    std::exception_ptr peer_local_failure = first_protocol_failure;
    if (peer_local_failure == nullptr) {
      try {
        request_send_count = core::checkedIntegralNarrow<int>(
            send_count, "ghost refresh request send count");
        const std::size_t checked_recv_count = core::checkedIntegralNarrow<std::size_t>(
            recv_count, "ghost refresh request receive count");
        request_recv_count = core::checkedIntegralNarrow<int>(
            checked_recv_count, "ghost refresh request receive MPI count");
        received_request_ids.resize(checked_recv_count);
        injectMpiTestFault(mpi_context, "ghost_request_post_metadata");
      } catch (...) {
        peer_local_failure = std::current_exception();
        if (first_protocol_failure == nullptr) {
          first_protocol_failure = peer_local_failure;
        }
      }
    }

    int local_ready = peer_local_failure == nullptr ? 1 : 0;
    int peer_ready = 0;
    if (MPI_Sendrecv(
            &local_ready, 1, MPI_INT, peer_rank,
            ghostExchangeSequencedTag(
                k_request_ready_tag_base, world_rank, peer_rank,
                expected_epoch.ghost_sync_epoch),
            &peer_ready, 1, MPI_INT, peer_rank,
            ghostExchangeSequencedTag(
                k_request_ready_tag_base, world_rank, peer_rank,
                expected_epoch.ghost_sync_epoch),
            MPI_COMM_WORLD, MPI_STATUS_IGNORE) != MPI_SUCCESS) {
      throw std::runtime_error(
          "executeBlockingGhostRefreshExchangeFromDescriptors: request-readiness Sendrecv failed");
    }
    if (local_ready == 0 || peer_ready == 0) {
      if (peer_ready == 0) {
        record_peer_failure(peer_rank, "request payload preparation");
      }
      continue;
    }

    if (MPI_Sendrecv(
            request_ids.empty() ? nullptr : const_cast<std::uint64_t*>(request_ids.data()),
            request_send_count,
            MPI_UINT64_T, peer_rank,
            ghostExchangeSequencedTag(
                k_request_payload_tag_base, world_rank, peer_rank,
                expected_epoch.ghost_sync_epoch),
            received_request_ids.empty() ? nullptr : received_request_ids.data(),
            request_recv_count,
            MPI_UINT64_T, peer_rank,
            ghostExchangeSequencedTag(
                k_request_payload_tag_base, world_rank, peer_rank,
                expected_epoch.ghost_sync_epoch),
            MPI_COMM_WORLD, MPI_STATUS_IGNORE) != MPI_SUCCESS) {
      throw std::runtime_error(
          "executeBlockingGhostRefreshExchangeFromDescriptors: request-payload Sendrecv failed");
    }

    peer_local_failure = nullptr;
    try {
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
      injectMpiTestFault(mpi_context, "ghost_request_post_payload");
    } catch (...) {
      peer_local_failure = std::current_exception();
      if (first_protocol_failure == nullptr) {
        first_protocol_failure = peer_local_failure;
      }
    }

    local_ready = peer_local_failure == nullptr ? 1 : 0;
    peer_ready = 0;
    if (MPI_Sendrecv(
            &local_ready, 1, MPI_INT, peer_rank,
            ghostExchangeSequencedTag(
                k_request_mapping_ready_tag_base, world_rank, peer_rank,
                expected_epoch.ghost_sync_epoch),
            &peer_ready, 1, MPI_INT, peer_rank,
            ghostExchangeSequencedTag(
                k_request_mapping_ready_tag_base, world_rank, peer_rank,
                expected_epoch.ghost_sync_epoch),
            MPI_COMM_WORLD, MPI_STATUS_IGNORE) != MPI_SUCCESS) {
      throw std::runtime_error(
          "executeBlockingGhostRefreshExchangeFromDescriptors: request-mapping readiness Sendrecv failed");
    }
    if (peer_ready == 0) {
      record_peer_failure(peer_rank, "request mapping validation");
    }
  }

  mpi_context.rethrowCollectivePreparationFailure(
      first_protocol_failure,
      "ghost refresh request discovery");
  std::vector<int> neighbor_ranks;
  std::vector<std::vector<std::uint32_t>> send_indices_by_neighbor;
  std::vector<std::vector<std::uint32_t>> recv_indices_by_neighbor;
  BlockingGhostRefreshExchange exchange;
  local_preparation_failure = nullptr;
  try {
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
    exchange.plan = buildExplicitGhostExchangePlan(
        world_rank,
        neighbor_ranks,
        send_indices_by_neighbor,
        recv_indices_by_neighbor,
        ghostRefreshPayloadRecordBytes(),
        expected_epoch);
    exchange.plan.exchange_sequence = expected_epoch.ghost_sync_epoch;
  } catch (...) {
    local_preparation_failure = std::current_exception();
  }
  mpi_context.rethrowCollectivePreparationFailure(
      local_preparation_failure,
      "ghost refresh payload-plan preparation");
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
  BlockingGhostExchangeResult result;
  result.received_ghosts.epoch = expected_epoch;

  std::vector<GhostExchangeBuffer> send_buffers;
  std::exception_ptr local_preparation_failure;
  try {
    validateBlockingGhostExchangeContracts(
        plan, local_ghost_descriptors, mpi_context.worldRank(), expected_epoch);
    if (!authoritative_local_state.isConsistent()) {
      throw std::invalid_argument(
          "executeBlockingGhostRefreshExchange: authoritative local ghost payload state is inconsistent");
    }
    if (!authoritative_local_state.epoch.matches(expected_epoch)) {
      throw std::invalid_argument(
          "executeBlockingGhostRefreshExchange: authoritative local payload epoch is stale");
    }
    if (!plan.neighbor_ranks.empty() && !mpi_context.isEnabled()) {
      throw std::runtime_error(
          "executeBlockingGhostRefreshExchange: non-empty ghost exchange requires MPI; serial path must have no neighbors");
    }

    send_buffers.resize(plan.neighbor_ranks.size());
    std::size_t total_receive_records = 0U;
    for (std::size_t slot = 0; slot < plan.neighbor_ranks.size(); ++slot) {
      for (const std::uint32_t local_index : plan.send_local_indices_by_neighbor[slot]) {
        if (local_index >= authoritative_local_state.size() ||
            local_index >= local_ghost_descriptors.size()) {
          throw std::out_of_range(
              "executeBlockingGhostRefreshExchange: send descriptor index is outside local payload state");
        }
        if (authoritative_local_state.entity_id[local_index] !=
            local_ghost_descriptors[local_index].particle_id) {
          throw std::invalid_argument(
              "executeBlockingGhostRefreshExchange: send payload entity_id does not match owned descriptor particle_id");
        }
      }
      send_buffers[slot].packFrom(
          plan.outbound_transfers[slot],
          authoritative_local_state,
          plan.send_local_indices_by_neighbor[slot]);
      total_receive_records = core::checkedSizeAdd(
          total_receive_records,
          plan.recv_local_indices_by_neighbor[slot].size(),
          "ghost refresh total receive record count");
    }

    // Reserve only lanes that are semantically present in the authoritative
    // payload. Optional-lane absence is part of the wire contract: in
    // particular, generic DMO particle ghosts must not allocate or advertise
    // gas/hydro lanes merely because an exchange occurs.
    result.received_ghosts.entity_id.reserve(total_receive_records);
    const auto reserve_if_source_present = [total_receive_records](
                                                const std::vector<double>& source,
                                                std::vector<double>* destination) {
      if (!source.empty()) {
        destination->reserve(total_receive_records);
      }
    };
    reserve_if_source_present(authoritative_local_state.position_x_comoving, &result.received_ghosts.position_x_comoving);
    reserve_if_source_present(authoritative_local_state.position_y_comoving, &result.received_ghosts.position_y_comoving);
    reserve_if_source_present(authoritative_local_state.position_z_comoving, &result.received_ghosts.position_z_comoving);
    reserve_if_source_present(authoritative_local_state.mass_code, &result.received_ghosts.mass_code);
    reserve_if_source_present(authoritative_local_state.density_code, &result.received_ghosts.density_code);
    reserve_if_source_present(authoritative_local_state.velocity_x_code, &result.received_ghosts.velocity_x_code);
    reserve_if_source_present(authoritative_local_state.velocity_y_code, &result.received_ghosts.velocity_y_code);
    reserve_if_source_present(authoritative_local_state.velocity_z_code, &result.received_ghosts.velocity_z_code);
    reserve_if_source_present(authoritative_local_state.pressure_code, &result.received_ghosts.pressure_code);
    reserve_if_source_present(authoritative_local_state.internal_energy_code, &result.received_ghosts.internal_energy_code);
  } catch (...) {
    local_preparation_failure = std::current_exception();
  }

  if (mpi_context.isEnabled()) {
    mpi_context.rethrowCollectivePreparationFailure(
        local_preparation_failure,
        "ghost refresh payload pre-communication preparation");
  } else if (local_preparation_failure != nullptr) {
    std::rethrow_exception(local_preparation_failure);
  }
  // A rank with no local peer payload still remains a participant in the
  // MPI-world protocol: MPI-enabled ranks must reach the final distributed
  // failure agreement below even when this peer loop is empty. The serial
  // no-neighbor path can return locally because it has no world collective.
  if (!mpi_context.isEnabled() && plan.neighbor_ranks.empty()) {
    return result;
  }

#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  constexpr int k_size_tag_base = 6810;
  constexpr int k_ready_tag_base = 7310;
  constexpr int k_payload_tag_base = 7810;
  constexpr int k_decode_ready_tag_base = 8310;
  std::exception_ptr first_protocol_failure;
  const auto record_peer_failure = [&](int peer_rank, std::string_view phase) {
    if (first_protocol_failure != nullptr) {
      return;
    }
    try {
      throw std::runtime_error(
          "ghost payload peer " + std::to_string(peer_rank) +
          " rejected " + std::string(phase));
    } catch (...) {
      first_protocol_failure = std::current_exception();
    }
  };

  for (std::size_t slot = 0; slot < plan.neighbor_ranks.size(); ++slot) {
    const auto send_bytes = send_buffers[slot].encodedBytes();
    const std::uint64_t send_size = core::checkedIntegralNarrow<std::uint64_t>(
        send_bytes.size(), "ghost refresh payload send size");
    std::uint64_t recv_size = 0U;
    const int peer_rank = plan.neighbor_ranks[slot];
    if (MPI_Sendrecv(
            const_cast<std::uint64_t*>(&send_size), 1, MPI_UINT64_T,
            peer_rank,
            ghostExchangeSequencedTag(
                k_size_tag_base, mpi_context.worldRank(), peer_rank,
                plan.exchange_sequence),
            &recv_size, 1, MPI_UINT64_T,
            peer_rank,
            ghostExchangeSequencedTag(
                k_size_tag_base, mpi_context.worldRank(), peer_rank,
                plan.exchange_sequence),
            MPI_COMM_WORLD, MPI_STATUS_IGNORE) != MPI_SUCCESS) {
      throw std::runtime_error(
          "executeBlockingGhostRefreshExchange: payload-size Sendrecv failed");
    }

    std::vector<std::uint8_t> recv_bytes;
    int payload_send_count = 0;
    int payload_recv_count = 0;
    std::exception_ptr peer_local_failure = first_protocol_failure;
    if (peer_local_failure == nullptr) {
      try {
        payload_send_count = core::checkedIntegralNarrow<int>(
            send_size, "ghost refresh payload send MPI count");
        const std::size_t checked_recv_size = core::checkedIntegralNarrow<std::size_t>(
            recv_size, "ghost refresh payload receive size");
        payload_recv_count = core::checkedIntegralNarrow<int>(
            checked_recv_size, "ghost refresh payload receive MPI count");
        recv_bytes.resize(checked_recv_size);
        injectMpiTestFault(mpi_context, "ghost_payload_post_metadata");
      } catch (...) {
        peer_local_failure = std::current_exception();
        if (first_protocol_failure == nullptr) {
          first_protocol_failure = peer_local_failure;
        }
      }
    }

    int local_ready = peer_local_failure == nullptr ? 1 : 0;
    int peer_ready = 0;
    if (MPI_Sendrecv(
            &local_ready, 1, MPI_INT, peer_rank,
            ghostExchangeSequencedTag(
                k_ready_tag_base, mpi_context.worldRank(), peer_rank,
                plan.exchange_sequence),
            &peer_ready, 1, MPI_INT, peer_rank,
            ghostExchangeSequencedTag(
                k_ready_tag_base, mpi_context.worldRank(), peer_rank,
                plan.exchange_sequence),
            MPI_COMM_WORLD, MPI_STATUS_IGNORE) != MPI_SUCCESS) {
      throw std::runtime_error(
          "executeBlockingGhostRefreshExchange: payload-readiness Sendrecv failed");
    }
    if (local_ready == 0 || peer_ready == 0) {
      if (peer_ready == 0) {
        record_peer_failure(peer_rank, "payload preparation");
      }
      continue;
    }

    if (MPI_Sendrecv(
            send_bytes.empty() ? nullptr : const_cast<std::uint8_t*>(send_bytes.data()),
            payload_send_count,
            MPI_BYTE, peer_rank,
            ghostExchangeSequencedTag(
                k_payload_tag_base, mpi_context.worldRank(), peer_rank,
                plan.exchange_sequence),
            recv_bytes.empty() ? nullptr : recv_bytes.data(),
            payload_recv_count,
            MPI_BYTE, peer_rank,
            ghostExchangeSequencedTag(
                k_payload_tag_base, mpi_context.worldRank(), peer_rank,
                plan.exchange_sequence),
            MPI_COMM_WORLD, MPI_STATUS_IGNORE) != MPI_SUCCESS) {
      throw std::runtime_error(
          "executeBlockingGhostRefreshExchange: ghost payload Sendrecv failed");
    }

    peer_local_failure = nullptr;
    try {
      GhostExchangeBuffer recv_buffer;
      recv_buffer.replaceEncodedBytes(std::move(recv_bytes));
      const std::size_t old_received_count = result.received_ghosts.size();
      recv_buffer.unpackAppendTo(
          plan.inbound_transfers[slot], result.received_ghosts);
      for (std::size_t i = 0;
           i < plan.recv_local_indices_by_neighbor[slot].size(); ++i) {
        const std::uint32_t local_index =
            plan.recv_local_indices_by_neighbor[slot][i];
        if (local_index >= local_ghost_descriptors.size()) {
          throw std::out_of_range(
              "executeBlockingGhostRefreshExchange: receive descriptor index is outside residency table");
        }
        if (result.received_ghosts.entity_id[old_received_count + i] !=
            local_ghost_descriptors[local_index].particle_id) {
          throw std::invalid_argument(
              "executeBlockingGhostRefreshExchange: received ghost entity_id does not match receive descriptor particle_id");
        }
      }
      injectMpiTestFault(mpi_context, "ghost_payload_post_payload");
      result.sent_bytes = core::checkedSizeAdd(
          core::checkedIntegralNarrow<std::size_t>(
              result.sent_bytes, "ghost refresh sent byte accumulator"),
          core::checkedIntegralNarrow<std::size_t>(
              send_size, "ghost refresh sent byte count"),
          "ghost refresh sent byte aggregate");
      result.received_bytes = core::checkedSizeAdd(
          core::checkedIntegralNarrow<std::size_t>(
              result.received_bytes, "ghost refresh received byte accumulator"),
          core::checkedIntegralNarrow<std::size_t>(
              recv_size, "ghost refresh received byte count"),
          "ghost refresh received byte aggregate");
    } catch (...) {
      peer_local_failure = std::current_exception();
      if (first_protocol_failure == nullptr) {
        first_protocol_failure = peer_local_failure;
      }
    }

    local_ready = peer_local_failure == nullptr ? 1 : 0;
    peer_ready = 0;
    if (MPI_Sendrecv(
            &local_ready, 1, MPI_INT, peer_rank,
            ghostExchangeSequencedTag(
                k_decode_ready_tag_base, mpi_context.worldRank(), peer_rank,
                plan.exchange_sequence),
            &peer_ready, 1, MPI_INT, peer_rank,
            ghostExchangeSequencedTag(
                k_decode_ready_tag_base, mpi_context.worldRank(), peer_rank,
                plan.exchange_sequence),
            MPI_COMM_WORLD, MPI_STATUS_IGNORE) != MPI_SUCCESS) {
      throw std::runtime_error(
          "executeBlockingGhostRefreshExchange: payload-decode readiness Sendrecv failed");
    }
    if (peer_ready == 0) {
      record_peer_failure(peer_rank, "payload decode/validation");
    }
  }

  mpi_context.rethrowCollectivePreparationFailure(
      first_protocol_failure,
      "ghost refresh sequential peer payload exchange");
  return result;
#else
  throw std::runtime_error(
      "executeBlockingGhostRefreshExchange: MPI support is not compiled in");
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
    int local_rank,
    int configured_gpu_devices,
    bool cuda_runtime_available,
    int visible_device_count) {
  if (local_rank < 0) {
    throw std::invalid_argument("local_rank must be non-negative");
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
  assignment.assigned_device_index = local_rank % configured_gpu_devices;
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
  topology.local_rank = mpi_context.localRank();
  topology.mpi_enabled = mpi_context.isEnabled();
  topology.pm_decomposition_mode = std::move(pm_decomposition_mode);
  topology.pm_slab = makePmSlabLayout(global_nx, global_ny, global_nz, mpi_context.worldSize(), mpi_context.worldRank());
  topology.device_assignment =
      selectRankDeviceAssignment(mpi_context.localRank(), configured_gpu_devices, cuda_runtime_available, visible_device_count);
  if (!topology.device_assignment.isValid()) {
    throw std::runtime_error("constructed distributed execution topology is invalid");
  }
  return topology;
}


}  // namespace cosmosim::parallel
