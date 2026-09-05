#include "cosmosim/core/simulation_state.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <new>
#include <numeric>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>
#include <unordered_map>
#include <unordered_set>

namespace cosmosim::core {


void GasCellIdentityMap::assign(std::vector<GasCellIdentityRecord> records) {
  GasCellIdentityMap candidate;
  candidate.m_records = std::move(records);
  if (!candidate.rebuildLookupTables()) {
    throw std::invalid_argument(
        "GasCellIdentityMap.assign: gas_cell_id must be nonzero, present parent_particle_id must be nonzero, and gas_cell_id/local_cell_row must be unique");
  }

  m_records = std::move(candidate.m_records);
  m_index_by_gas_cell_id = std::move(candidate.m_index_by_gas_cell_id);
  m_index_by_local_row = std::move(candidate.m_index_by_local_row);
  ++m_generation;
}

void GasCellIdentityMap::assignWithGeneration(
    std::vector<GasCellIdentityRecord> records,
    std::uint64_t generation) {
  GasCellIdentityMap candidate;
  candidate.m_records = std::move(records);
  if (!candidate.rebuildLookupTables()) {
    throw std::invalid_argument(
        "GasCellIdentityMap.assignWithGeneration: gas_cell_id must be nonzero, present parent_particle_id must be nonzero, and gas_cell_id/local_cell_row must be unique");
  }

  m_records = std::move(candidate.m_records);
  m_index_by_gas_cell_id = std::move(candidate.m_index_by_gas_cell_id);
  m_index_by_local_row = std::move(candidate.m_index_by_local_row);
  m_generation = generation;
}

void GasCellIdentityMap::clear() noexcept {
  m_records.clear();
  m_index_by_gas_cell_id.clear();
  m_index_by_local_row.clear();
  ++m_generation;
}

std::size_t GasCellIdentityMap::size() const noexcept { return m_records.size(); }

bool GasCellIdentityMap::empty() const noexcept { return m_records.empty(); }

std::uint64_t GasCellIdentityMap::generation() const noexcept { return m_generation; }

std::uint64_t GasCellIdentityMap::recordsCurrentSizeBytes() const {
  return checkedIntegralNarrow<std::uint64_t>(
      checkedSizeMultiply(
          m_records.size(), sizeof(GasCellIdentityRecord),
          "GasCellIdentityMap record current bytes"),
      "GasCellIdentityMap record current byte width");
}

std::uint64_t GasCellIdentityMap::recordsOwnedCapacityBytes() const {
  return checkedIntegralNarrow<std::uint64_t>(
      checkedSizeMultiply(
          m_records.capacity(), sizeof(GasCellIdentityRecord),
          "GasCellIdentityMap record capacity bytes"),
      "GasCellIdentityMap record capacity byte width");
}

std::uint64_t GasCellIdentityMap::lookupOwnedBytesConservative() const {
  const auto map_bytes = [](const auto& map, std::string_view context) {
    using Map = std::remove_cvref_t<decltype(map)>;
    // libstdc++/libc++ do not expose node allocation size. Account the bucket
    // array exactly by bucket_count and conservatively model every live node
    // as value_type plus two links. This is deliberately an upper-oriented
    // baseline model; allocator bookkeeping remains inside the process safety
    // margin and RSS/PSS reconciliation.
    const std::size_t bucket_bytes = checkedSizeMultiply(
        map.bucket_count(), sizeof(void*), context);
    const std::size_t node_bytes_each = checkedSizeAdd(
        sizeof(typename Map::value_type), 2U * sizeof(void*), context);
    const std::size_t node_bytes = checkedSizeMultiply(map.size(), node_bytes_each, context);
    return checkedSizeAdd(bucket_bytes, node_bytes, context);
  };
  const std::size_t gas_lookup = map_bytes(
      m_index_by_gas_cell_id, "GasCellIdentityMap gas-id lookup bytes");
  const std::size_t row_lookup = map_bytes(
      m_index_by_local_row, "GasCellIdentityMap row lookup bytes");
  return checkedIntegralNarrow<std::uint64_t>(
      checkedSizeAdd(gas_lookup, row_lookup, "GasCellIdentityMap lookup total bytes"),
      "GasCellIdentityMap lookup total byte width");
}

bool GasCellIdentityMap::rebuildLookupTables() {
  std::unordered_map<std::uint64_t, std::size_t> index_by_gas_cell_id;
  std::unordered_map<std::uint32_t, std::size_t> index_by_local_row;
  index_by_gas_cell_id.reserve(m_records.size());
  index_by_local_row.reserve(m_records.size());

  for (std::size_t index = 0; index < m_records.size(); ++index) {
    const auto& record = m_records[index];
    if (record.gas_cell_id == 0U ||
        (record.parent_particle_id.has_value() && *record.parent_particle_id == 0U)) {
      return false;
    }
    if (!index_by_gas_cell_id.emplace(record.gas_cell_id, index).second) {
      return false;
    }
    if (!index_by_local_row.emplace(record.local_cell_row, index).second) {
      return false;
    }
  }

  m_index_by_gas_cell_id.swap(index_by_gas_cell_id);
  m_index_by_local_row.swap(index_by_local_row);
  return true;
}

bool GasCellIdentityMap::isConsistent() const noexcept {
  if (m_index_by_gas_cell_id.size() != m_records.size() || m_index_by_local_row.size() != m_records.size()) {
    return false;
  }
  for (std::size_t index = 0; index < m_records.size(); ++index) {
    const auto& record = m_records[index];
    if (record.gas_cell_id == 0U ||
        (record.parent_particle_id.has_value() && *record.parent_particle_id == 0U)) {
      return false;
    }
    const auto gas_it = m_index_by_gas_cell_id.find(record.gas_cell_id);
    if (gas_it == m_index_by_gas_cell_id.end() || gas_it->second != index) {
      return false;
    }
    const auto row_it = m_index_by_local_row.find(record.local_cell_row);
    if (row_it == m_index_by_local_row.end() || row_it->second != index) {
      return false;
    }
  }
  return true;
}

bool GasCellIdentityMap::coversDenseLocalRows(std::size_t cell_count) const noexcept {
  if (m_records.size() != cell_count) {
    return false;
  }
  for (std::size_t row = 0; row < cell_count; ++row) {
    if (m_index_by_local_row.find(checkedLocalCellRow(row, "GasCellIdentityMap::coversDenseLocalRows")) == m_index_by_local_row.end()) {
      return false;
    }
  }
  return true;
}

void GasCellIdentityMap::requireCoversDenseLocalRows(std::size_t cell_count, std::string_view caller) const {
  if (!isConsistent()) {
    throw std::runtime_error(std::string(caller) + ": GasCellIdentityMap is internally inconsistent");
  }
  if (!coversDenseLocalRows(cell_count)) {
    throw std::runtime_error(
        std::string(caller) + ": GasCellIdentityMap does not cover exactly the dense local cell rows [0, cell_count)");
  }
}

std::span<const GasCellIdentityRecord> GasCellIdentityMap::records() const noexcept { return m_records; }

const GasCellIdentityRecord* GasCellIdentityMap::findByGasCellId(std::uint64_t gas_cell_id) const noexcept {
  const auto it = m_index_by_gas_cell_id.find(gas_cell_id);
  if (it == m_index_by_gas_cell_id.end()) {
    return nullptr;
  }
  return &m_records[it->second];
}

const GasCellIdentityRecord* GasCellIdentityMap::findByLocalRow(std::uint32_t local_cell_row) const noexcept {
  const auto it = m_index_by_local_row.find(local_cell_row);
  if (it == m_index_by_local_row.end()) {
    return nullptr;
  }
  return &m_records[it->second];
}

std::optional<std::uint32_t> GasCellIdentityMap::rowForGasCellId(std::uint64_t gas_cell_id) const noexcept {
  const auto* record = findByGasCellId(gas_cell_id);
  if (record == nullptr) {
    return std::nullopt;
  }
  return record->local_cell_row;
}

std::optional<std::uint64_t> GasCellIdentityMap::gasCellIdForLocalRow(std::uint32_t local_cell_row) const noexcept {
  const auto* record = findByLocalRow(local_cell_row);
  if (record == nullptr) {
    return std::nullopt;
  }
  return record->gas_cell_id;
}

std::optional<std::uint64_t> GasCellIdentityMap::parentParticleIdForGasCellId(
    std::uint64_t gas_cell_id) const noexcept {
  const auto* record = findByGasCellId(gas_cell_id);
  if (record == nullptr) {
    return std::nullopt;
  }
  return record->parent_particle_id;
}

std::optional<std::uint64_t> GasCellIdentityMap::owningPatchIdForGasCellId(
    std::uint64_t gas_cell_id) const noexcept {
  const auto* record = findByGasCellId(gas_cell_id);
  if (record == nullptr) {
    return std::nullopt;
  }
  return record->owning_patch_id;
}

std::vector<std::uint32_t> GasCellIdentityMap::rowsForParentParticleId(std::uint64_t parent_particle_id) const {
  std::vector<std::uint32_t> rows;
  for (const auto& record : m_records) {
    if (record.parent_particle_id.has_value() && *record.parent_particle_id == parent_particle_id) {
      rows.push_back(record.local_cell_row);
    }
  }
  std::sort(rows.begin(), rows.end());
  return rows;
}

std::vector<std::uint32_t> GasCellIdentityMap::rowsForPatch(std::uint64_t owning_patch_id) const {
  std::vector<std::uint32_t> rows;
  for (const auto& record : m_records) {
    if (record.owning_patch_id == owning_patch_id) {
      rows.push_back(record.local_cell_row);
    }
  }
  std::sort(rows.begin(), rows.end());
  return rows;
}

std::vector<std::uint32_t> buildGasCellNewToOldRowMap(
    const GasCellIdentityMap& old_map,
    const GasCellIdentityMap& new_map) {
  if (!old_map.isConsistent() || !new_map.isConsistent()) {
    throw std::runtime_error("buildGasCellNewToOldRowMap: inconsistent gas-cell identity map");
  }

  std::vector<std::uint32_t> new_to_old(new_map.size(), kInvalidGasCellRow);
  for (const auto& new_record : new_map.records()) {
    if (new_record.local_cell_row >= new_to_old.size()) {
      throw std::runtime_error("buildGasCellNewToOldRowMap: new map has non-dense local rows");
    }
    if (const auto old_row = old_map.rowForGasCellId(new_record.gas_cell_id); old_row.has_value()) {
      new_to_old[new_record.local_cell_row] = *old_row;
    }
  }
  return new_to_old;
}

bool ParticleReorderMap::isConsistent(std::size_t particle_count) const {
  if (old_to_new_index.size() != particle_count || new_to_old_index.size() != particle_count) {
    return false;
  }

  std::vector<std::uint8_t> visited(particle_count, 0U);
  for (std::size_t i = 0; i < particle_count; ++i) {
    const auto mapped = old_to_new_index[i];
    if (mapped >= particle_count) {
      return false;
    }
    if (visited[mapped] != 0U) {
      return false;
    }
    visited[mapped] = 1U;
    if (new_to_old_index[mapped] != i) {
      return false;
    }
  }
  return true;
}

std::uint64_t PendingFluxRegisterStore::currentSizeBytes() const {
  return checkedIntegralNarrow<std::uint64_t>(
      checkedSizeMultiply(
          m_records.size(), sizeof(PendingFluxRegisterRecord),
          "PendingFluxRegisterStore current bytes"),
      "PendingFluxRegisterStore current byte width");
}

std::uint64_t PendingFluxRegisterStore::ownedCapacityBytes() const {
  return checkedIntegralNarrow<std::uint64_t>(
      checkedSizeMultiply(
          m_records.capacity(), sizeof(PendingFluxRegisterRecord),
          "PendingFluxRegisterStore capacity bytes"),
      "PendingFluxRegisterStore capacity byte width");
}

std::uint64_t AmrTemporalBoundaryHistoryStore::currentSizeBytes() const {
  std::size_t bytes = checkedSizeMultiply(
      m_records.size(), sizeof(AmrTemporalBoundaryHistoryRecord),
      "AmrTemporalBoundaryHistoryStore outer current bytes");
  for (const AmrTemporalBoundaryHistoryRecord& record : m_records) {
    bytes = checkedSizeAdd(
        bytes,
        checkedSizeMultiply(
            record.cells.size(), sizeof(AmrTemporalBoundaryHistoryCellRecord),
            "AmrTemporalBoundaryHistoryStore nested current bytes"),
        "AmrTemporalBoundaryHistoryStore current total bytes");
  }
  return checkedIntegralNarrow<std::uint64_t>(
      bytes, "AmrTemporalBoundaryHistoryStore current byte width");
}

std::uint64_t AmrTemporalBoundaryHistoryStore::ownedCapacityBytes() const {
  std::size_t bytes = checkedSizeMultiply(
      m_records.capacity(), sizeof(AmrTemporalBoundaryHistoryRecord),
      "AmrTemporalBoundaryHistoryStore outer capacity bytes");
  for (const AmrTemporalBoundaryHistoryRecord& record : m_records) {
    bytes = checkedSizeAdd(
        bytes,
        checkedSizeMultiply(
            record.cells.capacity(), sizeof(AmrTemporalBoundaryHistoryCellRecord),
            "AmrTemporalBoundaryHistoryStore nested capacity bytes"),
        "AmrTemporalBoundaryHistoryStore capacity total bytes");
  }
  return checkedIntegralNarrow<std::uint64_t>(
      bytes, "AmrTemporalBoundaryHistoryStore capacity byte width");
}

MonotonicScratchAllocator::MonotonicScratchAllocator(std::size_t initial_capacity_bytes)
    : MonotonicScratchAllocator(nullptr, initial_capacity_bytes) {}

MonotonicScratchAllocator::MonotonicScratchAllocator(
    MemoryGovernor* memory_governor,
    std::size_t initial_capacity_bytes)
    : m_memory_governor(memory_governor) {
  if (initial_capacity_bytes > 0U) {
    m_next_block_capacity_bytes = std::max<std::size_t>(1024U, initial_capacity_bytes);
    (void)appendBlock(initial_capacity_bytes);
  }
}

MonotonicScratchAllocator::Block& MonotonicScratchAllocator::appendBlock(
    std::size_t minimum_capacity_bytes) {
  const std::size_t capacity = std::max(minimum_capacity_bytes, m_next_block_capacity_bytes);
  const std::size_t new_total_capacity = checkedSizeAdd(
      m_total_capacity_bytes, capacity, "MonotonicScratchAllocator.appendBlock total capacity");

  MemoryReservation reservation;
  if (m_memory_governor != nullptr) {
    reservation = m_memory_governor->reserve(
        MemoryClass::kScratchArena,
        static_cast<std::uint64_t>(capacity),
        "core.transient_step_scratch");
  }

  Block block;
  // Reserve-before-allocate: if the backing allocation throws, the pending
  // reservation remains local and its destructor returns the reservation.
  block.storage = std::make_unique<std::byte[]>(capacity);
  block.capacity_bytes = capacity;
  block.offset_bytes = 0U;
  if (reservation.valid()) {
    reservation.commit();
    block.reservation = std::move(reservation);
  }
  // If vector growth throws, destruction of the local block releases both the
  // backing storage and its committed governor reservation.
  m_blocks.push_back(std::move(block));
  m_total_capacity_bytes = new_total_capacity;
  m_capacity_high_water_bytes = std::max(
      m_capacity_high_water_bytes, m_total_capacity_bytes);

  if (capacity <= std::numeric_limits<std::size_t>::max() / 2U) {
    m_next_block_capacity_bytes = capacity * 2U;
  } else {
    m_next_block_capacity_bytes = std::numeric_limits<std::size_t>::max();
  }
  return m_blocks.back();
}

std::byte* MonotonicScratchAllocator::allocateBytes(std::size_t bytes, std::size_t alignment) {
  if (alignment == 0U || (alignment & (alignment - 1U)) != 0U) {
    throw std::invalid_argument(
        "MonotonicScratchAllocator.allocateBytes: alignment must be a nonzero power of two");
  }
  if (bytes == 0U) {
    return nullptr;
  }

  const std::size_t worst_case_capacity = checkedSizeAdd(
      bytes, alignment - 1U, "MonotonicScratchAllocator.allocateBytes block capacity");

  const auto try_allocate_from_block = [&](Block& block) -> std::byte* {
    if (block.offset_bytes > block.capacity_bytes) {
      throw std::logic_error("MonotonicScratchAllocator block offset exceeds capacity");
    }
    const std::size_t previous_offset = block.offset_bytes;
    void* candidate = block.storage.get() + block.offset_bytes;
    std::size_t available = block.capacity_bytes - block.offset_bytes;
    void* aligned = std::align(alignment, bytes, candidate, available);
    if (aligned == nullptr) {
      return nullptr;
    }
    auto* aligned_bytes = static_cast<std::byte*>(aligned);
    const auto aligned_offset = static_cast<std::size_t>(aligned_bytes - block.storage.get());
    block.offset_bytes = checkedSizeAdd(
        aligned_offset, bytes, "MonotonicScratchAllocator.allocateBytes");
    const std::size_t consumed_bytes = block.offset_bytes - previous_offset;
    m_current_used_bytes = checkedSizeAdd(
        m_current_used_bytes, consumed_bytes,
        "MonotonicScratchAllocator.allocateBytes current used");
    m_logical_high_water_bytes = std::max(
        m_logical_high_water_bytes, m_current_used_bytes);
    return aligned_bytes;
  };

  for (std::size_t block_index = m_current_block; block_index < m_blocks.size(); ++block_index) {
    Block& block = m_blocks[block_index];
    if (std::byte* allocation = try_allocate_from_block(block); allocation != nullptr) {
      m_current_block = block_index;
      return allocation;
    }
  }

  Block& block = appendBlock(worst_case_capacity);
  m_current_block = m_blocks.size() - 1U;
  if (std::byte* allocation = try_allocate_from_block(block); allocation != nullptr) {
    return allocation;
  }
  throw std::logic_error("MonotonicScratchAllocator internal block sizing error");
}

void MonotonicScratchAllocator::reset() {
  for (Block& block : m_blocks) {
    block.offset_bytes = 0U;
  }
  m_current_block = 0U;
  m_current_used_bytes = 0U;
}

std::size_t MonotonicScratchAllocator::usedBytes() const noexcept {
  return m_current_used_bytes;
}

std::size_t MonotonicScratchAllocator::capacityBytes() const noexcept {
  return m_total_capacity_bytes;
}

std::size_t MonotonicScratchAllocator::logicalHighWaterBytes() const noexcept {
  return m_logical_high_water_bytes;
}

std::size_t MonotonicScratchAllocator::capacityHighWaterBytes() const noexcept {
  return m_capacity_high_water_bytes;
}

bool MonotonicScratchAllocator::governed() const noexcept {
  return m_memory_governor != nullptr;
}

ParticleReorderMap buildParticleReorderMap(const SimulationState& state, ParticleReorderMode mode) {
  ParticleReorderMap reorder_map;
  reorder_map.new_to_old_index.resize(state.particles.size());
  std::iota(reorder_map.new_to_old_index.begin(), reorder_map.new_to_old_index.end(), 0U);

  const auto key_comp = [&](std::uint32_t lhs, std::uint32_t rhs) {
    if (mode == ParticleReorderMode::kByTimeBin) {
      throw std::invalid_argument(
          "buildParticleReorderMap(kByTimeBin): time-bin ordering requires scheduler authority; "
          "use buildParticleReorderMapByScheduler(state, scheduler) instead of derived time_bin mirrors");
    }
    if (mode == ParticleReorderMode::kBySfcKey) {
      const auto lhs_key = state.particle_sidecar.sfc_key[lhs];
      const auto rhs_key = state.particle_sidecar.sfc_key[rhs];
      return std::tuple{lhs_key, lhs} < std::tuple{rhs_key, rhs};
    }
    const auto lhs_key = state.particle_sidecar.species_tag[lhs];
    const auto rhs_key = state.particle_sidecar.species_tag[rhs];
    return std::tuple{lhs_key, lhs} < std::tuple{rhs_key, rhs};
  };

  std::stable_sort(reorder_map.new_to_old_index.begin(), reorder_map.new_to_old_index.end(), key_comp);

  reorder_map.old_to_new_index.resize(state.particles.size());
  for (std::size_t new_index = 0; new_index < reorder_map.new_to_old_index.size(); ++new_index) {
    const auto old_index = reorder_map.new_to_old_index[new_index];
    reorder_map.old_to_new_index[old_index] = checkedLocalParticleRow(new_index, "buildParticleReorderMap new particle row");
  }

  return reorder_map;
}

template <typename T>
void reorderAlignedVector(AlignedVector<T>& values, std::span<const std::uint32_t> new_to_old_index) {
  if (values.size() != new_to_old_index.size()) {
    throw std::invalid_argument("reorderAlignedVector: permutation length must match lane length");
  }
  AlignedVector<T> reordered(values.size());
  for (std::size_t i = 0; i < new_to_old_index.size(); ++i) {
    reordered[i] = values[new_to_old_index[i]];
  }
  values.swap(reordered);
}

template <typename T>
void reorderOptionalParentLane(AlignedVector<T>& values, std::span<const std::uint32_t> new_to_old_index) {
  if (!values.empty()) {
    reorderAlignedVector(values, new_to_old_index);
  }
}

std::vector<std::uint32_t> buildSidecarRowOrderByParent(
    std::span<const std::uint32_t> particle_index,
    std::span<const std::uint32_t> old_to_new_index) {
  for (const auto parent_index : particle_index) {
    if (parent_index >= old_to_new_index.size()) {
      throw std::out_of_range("reorderParticles: sidecar particle index out of range");
    }
  }
  std::vector<std::uint32_t> row_order(particle_index.size());
  std::iota(row_order.begin(), row_order.end(), 0U);
  std::stable_sort(
      row_order.begin(),
      row_order.end(),
      [&](std::uint32_t lhs, std::uint32_t rhs) {
        const auto lhs_particle = particle_index[lhs];
        const auto rhs_particle = particle_index[rhs];
        return std::tuple{old_to_new_index[lhs_particle], lhs} <
               std::tuple{old_to_new_index[rhs_particle], rhs};
      });
  return row_order;
}

template <typename Visitor>
void visitStarSidecarLanes(StarParticleSidecar& sidecar, Visitor&& visitor) {
  visitor(sidecar.particle_index);
  visitor(sidecar.formation_scale_factor);
  visitor(sidecar.birth_mass_code);
  visitor(sidecar.metallicity_mass_fraction);
  visitor(sidecar.stellar_age_years_last);
  visitor(sidecar.stellar_returned_mass_cumulative_code);
  visitor(sidecar.stellar_returned_metals_cumulative_code);
  visitor(sidecar.stellar_newly_synthesized_metals_cumulative_code);
  visitor(sidecar.stellar_feedback_energy_cumulative_erg);
  visitor(sidecar.enrichment_carry_mass_code);
  visitor(sidecar.enrichment_carry_metals_code);
  visitor(sidecar.enrichment_carry_feedback_energy_erg);
  visitor(sidecar.enrichment_carry_momentum_code);
  visitor(sidecar.stellar_deposited_mass_cumulative_code);
  visitor(sidecar.stellar_deposited_metals_cumulative_code);
  visitor(sidecar.stellar_deposited_feedback_energy_cumulative_erg);
  for (std::size_t channel = 0; channel < sidecar.stellar_returned_mass_channel_cumulative_code.size(); ++channel) {
    visitor(sidecar.stellar_returned_mass_channel_cumulative_code[channel]);
    visitor(sidecar.stellar_returned_metals_channel_cumulative_code[channel]);
    visitor(sidecar.stellar_feedback_energy_channel_cumulative_erg[channel]);
  }
}

template <typename Visitor>
void visitBlackHoleSidecarLanes(BlackHoleParticleSidecar& sidecar, Visitor&& visitor) {
  visitor(sidecar.particle_index);
  visitor(sidecar.host_cell_index);
  visitor(sidecar.subgrid_mass_code);
  visitor(sidecar.accretion_rate_code);
  visitor(sidecar.feedback_energy_code);
  visitor(sidecar.eddington_ratio);
  visitor(sidecar.cumulative_accreted_mass_code);
  visitor(sidecar.cumulative_feedback_energy_code);
  visitor(sidecar.duty_cycle_active_time_code);
  visitor(sidecar.duty_cycle_total_time_code);
}

template <typename Visitor>
void visitTracerSidecarLanes(TracerParticleSidecar& sidecar, Visitor&& visitor) {
  visitor(sidecar.particle_index);
  visitor(sidecar.parent_particle_id);
  visitor(sidecar.injection_step);
  visitor(sidecar.host_cell_index);
  visitor(sidecar.mass_fraction_of_host);
  visitor(sidecar.last_host_mass_code);
  visitor(sidecar.cumulative_exchanged_mass_code);
}

template <typename Sidecar, typename LaneVisitor>
void moveSidecarRowsWithParent(
    Sidecar& sidecar,
    std::span<const std::uint32_t> old_to_new_index,
    LaneVisitor&& visit_lanes) {
  const auto row_order = buildSidecarRowOrderByParent(sidecar.particle_index, old_to_new_index);
  visit_lanes(sidecar, [&](auto& lane) { reorderAlignedVector(lane, row_order); });
}

void remapSidecarParticleIndex(
    AlignedVector<std::uint32_t>& particle_index,
    std::span<const std::uint32_t> old_to_new_index) {
  for (auto& index : particle_index) {
    if (index >= old_to_new_index.size()) {
      throw std::out_of_range("reorderParticles: sidecar particle index out of range");
    }
    index = old_to_new_index[index];
  }
}

void syncStarSidecarRows(
    StarParticleSidecar& sidecar,
    std::span<const std::uint32_t> old_to_new_index,
    SidecarSyncMode mode) {
  if (mode == SidecarSyncMode::kMoveWithParent) {
    moveSidecarRowsWithParent(sidecar, old_to_new_index, [](StarParticleSidecar& lanes, auto&& visitor) {
      visitStarSidecarLanes(lanes, visitor);
    });
  }
  remapSidecarParticleIndex(sidecar.particle_index, old_to_new_index);
}

void syncBlackHoleSidecarRows(
    BlackHoleParticleSidecar& sidecar,
    std::span<const std::uint32_t> old_to_new_index,
    SidecarSyncMode mode) {
  if (mode == SidecarSyncMode::kMoveWithParent) {
    moveSidecarRowsWithParent(sidecar, old_to_new_index, [](BlackHoleParticleSidecar& lanes, auto&& visitor) {
      visitBlackHoleSidecarLanes(lanes, visitor);
    });
  }
  remapSidecarParticleIndex(sidecar.particle_index, old_to_new_index);
}

void syncTracerSidecarRows(
    TracerParticleSidecar& sidecar,
    std::span<const std::uint32_t> old_to_new_index,
    SidecarSyncMode mode) {
  if (mode == SidecarSyncMode::kMoveWithParent) {
    moveSidecarRowsWithParent(sidecar, old_to_new_index, [](TracerParticleSidecar& lanes, auto&& visitor) {
      visitTracerSidecarLanes(lanes, visitor);
    });
  }
  remapSidecarParticleIndex(sidecar.particle_index, old_to_new_index);
}

void reorderParticles(
    SimulationState& state,
    const ParticleReorderMap& reorder_map,
    const SidecarSyncPolicy& sync_policy) {
  if (!reorder_map.isConsistent(state.particles.size())) {
    throw std::invalid_argument("reorderParticles: inconsistent reorder map");
  }
  if (state.cells.size() > 0) {
    requireGasCellIdentityMapCoversDenseRows(state, "reorderParticles");
    if (!gasCellIdentityMapMatchesSidecarLanes(state)) {
      throw std::invalid_argument("reorderParticles: GasCellIdentityMap must match gas-cell sidecar mirror lanes");
    }
    if (gasCellIdentityMapMatchesParticleBoundState(state)) {
      const auto gas_globals = state.particle_species_index.globalIndices(ParticleSpecies::kGas);
      std::vector<std::uint64_t> gas_particle_ids_old;
      gas_particle_ids_old.reserve(gas_globals.size());
      for (const auto gas_global : gas_globals) {
        gas_particle_ids_old.push_back(state.particle_sidecar.particle_id[gas_global]);
      }
      std::vector<std::uint64_t> gas_particle_ids_new;
      gas_particle_ids_new.reserve(gas_globals.size());
      for (const auto old_index : reorder_map.new_to_old_index) {
        if (state.particle_sidecar.species_tag[old_index] == static_cast<std::uint32_t>(ParticleSpecies::kGas)) {
          gas_particle_ids_new.push_back(state.particle_sidecar.particle_id[old_index]);
        }
      }
      if (gas_particle_ids_new != gas_particle_ids_old) {
        throw std::invalid_argument(
            "reorderParticles: legacy particle-bound gas-cell layout forbids gas-particle relative reorder without "
            "gas-cell rebuild");
      }
    }
  }

  const std::span<const std::uint32_t> new_to_old_index = reorder_map.new_to_old_index;

  reorderAlignedVector(state.particles.position_x_comoving, new_to_old_index);
  reorderAlignedVector(state.particles.position_y_comoving, new_to_old_index);
  reorderAlignedVector(state.particles.position_z_comoving, new_to_old_index);
  reorderAlignedVector(state.particles.velocity_x_peculiar, new_to_old_index);
  reorderAlignedVector(state.particles.velocity_y_peculiar, new_to_old_index);
  reorderAlignedVector(state.particles.velocity_z_peculiar, new_to_old_index);
  reorderAlignedVector(state.particles.mass_code, new_to_old_index);
  reorderAlignedVector(state.particles.time_bin, new_to_old_index);

  reorderAlignedVector(state.particle_sidecar.particle_id, new_to_old_index);
  reorderAlignedVector(state.particle_sidecar.sfc_key, new_to_old_index);
  reorderAlignedVector(state.particle_sidecar.species_tag, new_to_old_index);
  reorderAlignedVector(state.particle_sidecar.particle_flags, new_to_old_index);
  reorderAlignedVector(state.particle_sidecar.owning_rank, new_to_old_index);
  reorderAlignedVector(state.particle_sidecar.last_drift_time_code, new_to_old_index);
  reorderAlignedVector(state.particle_sidecar.last_drift_scale_factor, new_to_old_index);
  reorderOptionalParentLane(state.particle_sidecar.gravity_softening_comoving, new_to_old_index);
  reorderOptionalParentLane(state.particle_sidecar.has_gravity_softening_override, new_to_old_index);

  syncStarSidecarRows(state.star_particles, reorder_map.old_to_new_index, sync_policy.star_particles);
  syncBlackHoleSidecarRows(state.black_holes, reorder_map.old_to_new_index, sync_policy.black_holes);
  syncTracerSidecarRows(state.tracers, reorder_map.old_to_new_index, sync_policy.tracers);

  state.rebuildSpeciesIndex();
  debugAssertSpeciesSidecarOwnershipInvariants(state);
  state.bumpParticleIndexGeneration();
}

void debugAssertNoStaleParticleIndices(const SimulationState& state) {
  auto check_indices = [&](std::span<const std::uint32_t> indices, const char* name) {
    for (const auto index : indices) {
      if (index >= state.particles.size()) {
        throw std::runtime_error(std::string("debugAssertNoStaleParticleIndices: stale index in ") + name);
      }
    }
  };

  check_indices(state.star_particles.particle_index, "star_particles");
  check_indices(state.black_holes.particle_index, "black_holes");
  check_indices(state.tracers.particle_index, "tracers");
}

void debugAssertSpeciesSidecarOwnershipInvariants(const SimulationState& state) {
  if (!state.validateOwnershipInvariants()) {
    throw std::runtime_error(
        "debugAssertSpeciesSidecarOwnershipInvariants: species sidecars must be one-to-one with eligible parents");
  }
}

void refreshGasCellIdentityFromParticleOrder(SimulationState& state) {
  const auto gas_globals = state.particle_species_index.globalIndices(ParticleSpecies::kGas);
  if (state.cells.size() != gas_globals.size() || state.gas_cells.size() != gas_globals.size()) {
    throw std::runtime_error(
        "refreshGasCellIdentityFromParticleOrder: gas-cell count must match canonical gas-particle count");
  }
  state.gas_cells.gas_cell_id.resize(gas_globals.size());
  state.gas_cells.parent_particle_id.resize(gas_globals.size());
  for (std::size_t cell_index = 0; cell_index < gas_globals.size(); ++cell_index) {
    const std::uint32_t particle_index = gas_globals[cell_index];
    const std::uint64_t particle_id = state.particle_sidecar.particle_id[particle_index];
    state.gas_cells.gas_cell_id[cell_index] = particle_id;
    state.gas_cells.parent_particle_id[cell_index] = particle_id;
    state.gas_cells.velocity_x_peculiar[cell_index] = state.particles.velocity_x_peculiar[particle_index];
    state.gas_cells.velocity_y_peculiar[cell_index] = state.particles.velocity_y_peculiar[particle_index];
    state.gas_cells.velocity_z_peculiar[cell_index] = state.particles.velocity_z_peculiar[particle_index];
  }
  refreshGasCellIdentityMapFromParticleBoundState(state);
  state.bumpCellIndexGeneration();
}

bool gasCellIdentityMatchesParticleOrder(const SimulationState& state) {
  const auto gas_globals = state.particle_species_index.globalIndices(ParticleSpecies::kGas);
  if (state.cells.size() != gas_globals.size() || state.gas_cells.size() != gas_globals.size()) {
    return false;
  }
  if (state.gas_cells.gas_cell_id.size() != gas_globals.size() ||
      state.gas_cells.parent_particle_id.size() != gas_globals.size()) {
    return false;
  }
  std::unordered_set<std::uint64_t> seen_ids;
  seen_ids.reserve(gas_globals.size());
  for (std::size_t cell_index = 0; cell_index < gas_globals.size(); ++cell_index) {
    const std::uint32_t particle_index = gas_globals[cell_index];
    if (particle_index >= state.particles.size()) {
      return false;
    }
    const std::uint64_t expected_particle_id = state.particle_sidecar.particle_id[particle_index];
    if (state.gas_cells.parent_particle_id[cell_index] != expected_particle_id ||
        state.gas_cells.gas_cell_id[cell_index] != expected_particle_id) {
      return false;
    }
    if (!seen_ids.insert(state.gas_cells.gas_cell_id[cell_index]).second) {
      return false;
    }
  }
  return true;
}

void refreshGasCellIdentityMapFromParticleBoundState(SimulationState& state) {
  if (!gasCellIdentityMatchesParticleOrder(state)) {
    throw std::runtime_error(
        "refreshGasCellIdentityMapFromParticleBoundState: legacy gas-cell lanes must match particle-bound order");
  }

  std::vector<GasCellIdentityRecord> records;
  records.reserve(state.cells.size());
  for (std::size_t cell_index = 0; cell_index < state.cells.size(); ++cell_index) {
    const std::uint64_t gas_cell_id = state.gas_cells.gas_cell_id[cell_index];
    const std::uint64_t parent_particle_id = state.gas_cells.parent_particle_id[cell_index];
    if (gas_cell_id == 0U || parent_particle_id == 0U) {
      throw std::runtime_error(
          "refreshGasCellIdentityMapFromParticleBoundState: legacy particle-bound gas-cell IDs must be nonzero");
    }
    std::uint64_t owning_patch_id = 0;
    if (!state.patches.patch_id.empty()) {
      const std::uint32_t patch_index = state.cells.patch_index[cell_index];
      if (patch_index >= state.patches.size()) {
        throw std::runtime_error(
            "refreshGasCellIdentityMapFromParticleBoundState: cell patch_index is outside PatchSoa");
      }
      owning_patch_id = state.patches.patch_id[patch_index];
    }
    records.push_back(GasCellIdentityRecord{
        .gas_cell_id = gas_cell_id,
        .parent_particle_id = parent_particle_id,
        .owning_patch_id = owning_patch_id,
        .local_cell_row = checkedLocalCellRow(cell_index, "SimulationState gas-cell row"),
    });
  }
  state.gas_cell_identity.assign(std::move(records));
}

void refreshGasCellIdentityMapFromSidecarLanes(SimulationState& state) {
  if (state.cells.size() != state.gas_cells.size() ||
      state.gas_cells.gas_cell_id.size() != state.cells.size() ||
      state.gas_cells.parent_particle_id.size() != state.cells.size()) {
    throw std::runtime_error(
        "refreshGasCellIdentityMapFromSidecarLanes: gas-cell sidecar lanes must cover dense CellSoa rows");
  }
  if (!state.patches.patch_id.empty() && state.cells.patch_index.size() != state.cells.size()) {
    throw std::runtime_error(
        "refreshGasCellIdentityMapFromSidecarLanes: cells.patch_index must cover dense CellSoa rows when patches exist");
  }

  std::vector<GasCellIdentityRecord> records;
  records.reserve(state.cells.size());
  for (std::size_t cell_index = 0; cell_index < state.cells.size(); ++cell_index) {
    const std::uint64_t gas_cell_id = state.gas_cells.gas_cell_id[cell_index];
    if (gas_cell_id == 0U) {
      throw std::runtime_error(
          "refreshGasCellIdentityMapFromSidecarLanes: gas_cell_id lanes must be nonzero");
    }

    std::uint64_t owning_patch_id = 0U;
    if (!state.patches.patch_id.empty()) {
      const std::uint32_t patch_index = state.cells.patch_index[cell_index];
      if (patch_index >= state.patches.size()) {
        throw std::runtime_error(
            "refreshGasCellIdentityMapFromSidecarLanes: cell patch_index is outside PatchSoa");
      }
      owning_patch_id = state.patches.patch_id[patch_index];
    }

    const std::uint64_t mirrored_parent = state.gas_cells.parent_particle_id[cell_index];
    records.push_back(GasCellIdentityRecord{
        .gas_cell_id = gas_cell_id,
        .parent_particle_id = mirrored_parent == 0U
            ? std::optional<std::uint64_t>{}
            : std::optional<std::uint64_t>{mirrored_parent},
        .owning_patch_id = owning_patch_id,
        .local_cell_row = checkedLocalCellRow(cell_index, "SimulationState gas-cell row"),
    });
  }
  state.gas_cell_identity.assign(std::move(records));
}

bool gasCellIdentityMapMatchesParticleBoundState(const SimulationState& state) {
  if (!gasCellIdentityMatchesParticleOrder(state)) {
    return false;
  }
  if (!state.gas_cell_identity.isConsistent() ||
      !state.gas_cell_identity.coversDenseLocalRows(state.cells.size())) {
    return false;
  }
  for (std::size_t cell_index = 0; cell_index < state.cells.size(); ++cell_index) {
    const auto* record = state.gas_cell_identity.findByLocalRow(checkedLocalCellRow(cell_index, "SimulationState gas-cell row"));
    if (record == nullptr) {
      return false;
    }
    if (record->gas_cell_id != state.gas_cells.gas_cell_id[cell_index]) {
      return false;
    }
    if (!record->parent_particle_id.has_value() ||
        *record->parent_particle_id != state.gas_cells.parent_particle_id[cell_index]) {
      return false;
    }
    if (record->parent_particle_id.value() == 0U || record->gas_cell_id == 0U) {
      return false;
    }
    if (!state.patches.patch_id.empty()) {
      const std::uint32_t patch_index = state.cells.patch_index[cell_index];
      if (patch_index >= state.patches.size() || record->owning_patch_id != state.patches.patch_id[patch_index]) {
        return false;
      }
    }
  }
  return true;
}

bool gasCellIdentityMapMatchesSidecarLanes(const SimulationState& state) {
  if (!state.gas_cell_identity.isConsistent() ||
      !state.gas_cell_identity.coversDenseLocalRows(state.cells.size())) {
    return false;
  }
  if (state.gas_cells.gas_cell_id.size() != state.cells.size() ||
      state.gas_cells.parent_particle_id.size() != state.cells.size()) {
    return false;
  }
  for (std::size_t cell_index = 0; cell_index < state.cells.size(); ++cell_index) {
    const auto* record = state.gas_cell_identity.findByLocalRow(checkedLocalCellRow(cell_index, "SimulationState gas-cell row"));
    if (record == nullptr || record->gas_cell_id == 0U) {
      return false;
    }
    if (state.gas_cells.gas_cell_id[cell_index] != record->gas_cell_id) {
      return false;
    }
    const std::uint64_t mirrored_parent =
        record->parent_particle_id.has_value() ? *record->parent_particle_id : 0U;
    if (state.gas_cells.parent_particle_id[cell_index] != mirrored_parent) {
      return false;
    }
    if (!state.patches.patch_id.empty()) {
      const std::uint32_t patch_index = state.cells.patch_index[cell_index];
      if (patch_index >= state.patches.size() || record->owning_patch_id != state.patches.patch_id[patch_index]) {
        return false;
      }
    }
  }
  return true;
}

void SimulationState::refreshGasCellIdentityFromParticleOrder() {
  cosmosim::core::refreshGasCellIdentityFromParticleOrder(*this);
}

bool SimulationState::gasCellIdentityMatchesParticleOrder() const {
  return cosmosim::core::gasCellIdentityMatchesParticleOrder(*this);
}

void SimulationState::refreshGasCellIdentityMapFromParticleBoundState() {
  cosmosim::core::refreshGasCellIdentityMapFromParticleBoundState(*this);
}

void SimulationState::refreshGasCellIdentityMapFromSidecarLanes() {
  cosmosim::core::refreshGasCellIdentityMapFromSidecarLanes(*this);
}

void SimulationState::synchronizeGasCellIdentityCompatibilityMirrors() {
  gas_cell_identity.requireCoversDenseLocalRows(cells.size(),
      "synchronizeGasCellIdentityCompatibilityMirrors");
  if (gas_cells.gas_cell_id.size() != cells.size() ||
      gas_cells.parent_particle_id.size() != cells.size()) {
    throw std::runtime_error(
        "synchronizeGasCellIdentityCompatibilityMirrors: gas-cell sidecar lanes must cover dense CellSoa rows");
  }

  std::unordered_map<std::uint64_t, std::uint32_t> patch_row_by_id;
  if (!patches.patch_id.empty()) {
    if (!patches.isConsistent() || cells.patch_index.size() != cells.size()) {
      throw std::runtime_error(
          "synchronizeGasCellIdentityCompatibilityMirrors: PatchSoa and CellSoa patch-index lane must be consistent");
    }
    patch_row_by_id.reserve(patches.size());
    for (std::uint32_t patch_row = 0; patch_row < patches.size(); ++patch_row) {
      if (!patch_row_by_id.emplace(patches.patch_id[patch_row], patch_row).second) {
        throw std::runtime_error(
            "synchronizeGasCellIdentityCompatibilityMirrors: duplicate PatchSoa patch_id");
      }
    }
  }

  for (const GasCellIdentityRecord& record : gas_cell_identity.records()) {
    if (record.local_cell_row >= cells.size()) {
      throw std::runtime_error(
          "synchronizeGasCellIdentityCompatibilityMirrors: identity local row is outside CellSoa");
    }
    const std::uint32_t row = record.local_cell_row;
    gas_cells.gas_cell_id[row] = record.gas_cell_id;
    gas_cells.parent_particle_id[row] = record.parent_particle_id.value_or(0U);
    if (!patches.patch_id.empty()) {
      const auto patch_it = patch_row_by_id.find(record.owning_patch_id);
      if (patch_it == patch_row_by_id.end()) {
        throw std::runtime_error(
            "synchronizeGasCellIdentityCompatibilityMirrors: identity owning_patch_id is not present in PatchSoa");
      }
      cells.patch_index[row] = patch_it->second;
    }
  }

  if (!gasCellIdentityMapMatchesSidecarLanes()) {
    throw std::runtime_error(
        "synchronizeGasCellIdentityCompatibilityMirrors: map-to-mirror synchronization failed validation");
  }
}

void SimulationState::replaceGasCellIdentityRecords(std::vector<GasCellIdentityRecord> records) {
  if (records.size() != cells.size()) {
    throw std::invalid_argument(
        "replaceGasCellIdentityRecords: record count must match dense CellSoa rows");
  }
  gas_cell_identity.assign(std::move(records));
  synchronizeGasCellIdentityCompatibilityMirrors();
  bumpCellIndexGeneration();
}

void SimulationState::restoreGasCellIdentityRecords(
    std::vector<GasCellIdentityRecord> records,
    std::uint64_t generation) {
  if (records.size() != cells.size()) {
    throw std::invalid_argument(
        "restoreGasCellIdentityRecords: record count must match dense CellSoa rows");
  }
  gas_cell_identity.assignWithGeneration(std::move(records), generation);
  synchronizeGasCellIdentityCompatibilityMirrors();
  bumpCellIndexGeneration();
}

bool SimulationState::gasCellIdentityMapMatchesParticleBoundState() const {
  return cosmosim::core::gasCellIdentityMapMatchesParticleBoundState(*this);
}

bool SimulationState::gasCellIdentityMapMatchesSidecarLanes() const {
  return cosmosim::core::gasCellIdentityMapMatchesSidecarLanes(*this);
}

std::uint64_t SimulationState::gasCellIdentityGeneration() const noexcept {
  return gas_cell_identity.generation();
}

void requireGasCellIdentityMapCoversDenseRows(const SimulationState& state, std::string_view caller) {
  state.gas_cell_identity.requireCoversDenseLocalRows(state.cells.size(), caller);
  if (!gasCellIdentityMapMatchesSidecarLanes(state)) {
    throw std::runtime_error(std::string(caller) +
        ": GasCellIdentityMap and compatibility gas-cell mirror lanes diverged; "
        "use SimulationState::replaceGasCellIdentityRecords() or "
        "synchronizeGasCellIdentityCompatibilityMirrors()");
  }
}

void SimulationState::requireGasCellIdentityMapCoversDenseRows(std::string_view caller) const {
  cosmosim::core::requireGasCellIdentityMapCoversDenseRows(*this, caller);
}

void requireGasCellIdentityMapFresh(
    const SimulationState& state,
    std::uint64_t expected_generation,
    std::string_view caller) {
  if (state.gas_cell_identity.generation() != expected_generation) {
    throw std::runtime_error(std::string(caller) + ": stale GasCellIdentityMap generation");
  }
  requireGasCellIdentityMapCoversDenseRows(state, caller);
}

void SimulationState::requireGasCellIdentityMapFresh(
    std::uint64_t expected_generation,
    std::string_view caller) const {
  cosmosim::core::requireGasCellIdentityMapFresh(*this, expected_generation, caller);
}

std::optional<std::uint32_t> SimulationState::rowForGasCellId(std::uint64_t gas_cell_id) const noexcept {
  return gas_cell_identity.rowForGasCellId(gas_cell_id);
}

std::optional<std::uint64_t> SimulationState::gasCellIdForLocalRow(std::uint32_t local_cell_row) const noexcept {
  return gas_cell_identity.gasCellIdForLocalRow(local_cell_row);
}

std::optional<std::uint64_t> SimulationState::parentParticleIdForGasCellId(
    std::uint64_t gas_cell_id) const noexcept {
  return gas_cell_identity.parentParticleIdForGasCellId(gas_cell_id);
}

std::optional<std::uint64_t> SimulationState::owningPatchIdForGasCellId(std::uint64_t gas_cell_id) const noexcept {
  return gas_cell_identity.owningPatchIdForGasCellId(gas_cell_id);
}

void legacyRequireParticleBoundGasCellContract(const SimulationState& state, std::string_view caller) {
  const std::size_t gas_particle_count = state.particle_species_index.count(ParticleSpecies::kGas);
  const std::size_t cell_count = state.cells.size();
  const std::size_t gas_cell_count = state.gas_cells.size();
  if (cell_count == 0 && gas_cell_count == 0) {
    return;
  }
  if (cell_count != gas_cell_count) {
    throw std::runtime_error(std::string(caller) +
        ": legacy/import particle-bound gas-cell compatibility contract violated: CellSoa row count does not match "
        "GasCellSidecar row count");
  }
  if (gas_particle_count != cell_count) {
    throw std::runtime_error(std::string(caller) +
        ": legacy/import particle-bound gas-cell compatibility contract violated: local gas particle count must equal "
        "local gas cell row count; production hydro must use GasCellIdentityMap instead");
  }
  if (!gasCellIdentityMatchesParticleOrder(state)) {
    throw std::runtime_error(std::string(caller) +
        ": legacy/import particle-bound gas-cell compatibility contract violated: gas_cell_id and parent_particle_id "
        "lanes must equal the canonical gas particle IDs in local gas-cell row order");
  }
  if (!gasCellIdentityMapMatchesParticleBoundState(state)) {
    throw std::runtime_error(std::string(caller) +
        ": legacy/import particle-bound gas-cell compatibility contract violated: SimulationState::gas_cell_identity "
        "must match the compatibility gas_cell_id/parent_particle_id mirror lanes and cover dense local rows");
  }
}

void requireParticleBoundGasCellContract(const SimulationState& state, std::string_view caller) {
  legacyRequireParticleBoundGasCellContract(state, caller);
}

std::optional<std::uint64_t> parentParticleIdForGasCellRow(const SimulationState& state, std::uint32_t cell_index) {
  if (cell_index >= state.cells.size()) {
    throw std::out_of_range("parentParticleIdForGasCellRow: cell index out of range");
  }
  state.requireGasCellIdentityMapCoversDenseRows("parentParticleIdForGasCellRow");
  const auto* record = state.gas_cell_identity.findByLocalRow(cell_index);
  if (record == nullptr) {
    throw std::runtime_error("parentParticleIdForGasCellRow: GasCellIdentityMap is missing local row");
  }
  return record->parent_particle_id;
}

std::uint32_t gasCellRowForParticleId(const SimulationState& state, std::uint64_t particle_id) {
  legacyRequireParticleBoundGasCellContract(state, "gasCellRowForParticleId");
  for (std::uint32_t cell_index = 0; cell_index < state.gas_cells.parent_particle_id.size(); ++cell_index) {
    if (state.gas_cells.parent_particle_id[cell_index] == particle_id) {
      return cell_index;
    }
  }
  throw std::out_of_range("gasCellRowForParticleId: gas particle ID is not attached to a local gas-cell row");
}

std::uint32_t gasParticleIndexForCellRow(const SimulationState& state, std::uint32_t cell_index) {
  legacyRequireParticleBoundGasCellContract(state, "gasParticleIndexForCellRow");
  if (cell_index >= state.cells.size()) {
    throw std::out_of_range("gasParticleIndexForCellRow: cell index out of range");
  }
  const auto gas_globals = state.particle_species_index.globalIndices(ParticleSpecies::kGas);
  return gas_globals[cell_index];
}

void SimulationState::requireParticleBoundGasCellContract(std::string_view caller) const {
  cosmosim::core::requireParticleBoundGasCellContract(*this, caller);
}

void SimulationState::legacyRequireParticleBoundGasCellContract(std::string_view caller) const {
  cosmosim::core::legacyRequireParticleBoundGasCellContract(*this, caller);
}

std::optional<std::uint64_t> SimulationState::parentParticleIdForGasCellRow(std::uint32_t cell_index) const {
  return cosmosim::core::parentParticleIdForGasCellRow(*this, cell_index);
}

std::uint32_t SimulationState::gasCellRowForParticleId(std::uint64_t particle_id) const {
  return cosmosim::core::gasCellRowForParticleId(*this, particle_id);
}

std::uint32_t SimulationState::gasParticleIndexForCellRow(std::uint32_t cell_index) const {
  return cosmosim::core::gasParticleIndexForCellRow(*this, cell_index);
}

void debugAssertGasCellIdentityContract(const SimulationState& state) {
  legacyRequireParticleBoundGasCellContract(state, "debugAssertGasCellIdentityContract");
}

}  // namespace cosmosim::core
