#include "cosmosim/core/simulation_state.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <new>
#include <numeric>
#include <stdexcept>
#include <tuple>

namespace cosmosim::core {

bool ParticleReorderMap::isConsistent(std::size_t particle_count) const noexcept {
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

MonotonicScratchAllocator::MonotonicScratchAllocator(std::size_t initial_capacity_bytes)
    : m_storage(initial_capacity_bytes), m_offset_bytes(0) {}

std::byte* MonotonicScratchAllocator::allocateBytes(std::size_t bytes, std::size_t alignment) {
  if (alignment == 0 || (alignment & (alignment - 1U)) != 0) {
    throw std::invalid_argument("MonotonicScratchAllocator.allocateBytes: alignment must be power-of-two");
  }

  if (bytes == 0) {
    return m_storage.data() + m_offset_bytes;
  }

  const std::size_t aligned_offset = (m_offset_bytes + alignment - 1U) & ~(alignment - 1U);
  const std::size_t required_size = aligned_offset + bytes;

  if (required_size > m_storage.size()) {
    const std::size_t grow_size = std::max(required_size, std::max<std::size_t>(1024, m_storage.size() * 2));
    m_storage.resize(grow_size);
  }

  auto* ptr = m_storage.data() + aligned_offset;
  m_offset_bytes = required_size;
  return ptr;
}

void MonotonicScratchAllocator::reset() { m_offset_bytes = 0; }

std::size_t MonotonicScratchAllocator::capacityBytes() const noexcept { return m_storage.size(); }

ParticleReorderMap buildParticleReorderMap(const SimulationState& state, ParticleReorderMode mode) {
  ParticleReorderMap reorder_map;
  reorder_map.new_to_old_index.resize(state.particles.size());
  std::iota(reorder_map.new_to_old_index.begin(), reorder_map.new_to_old_index.end(), 0U);

  const auto key_comp = [&](std::uint32_t lhs, std::uint32_t rhs) {
    if (mode == ParticleReorderMode::kByTimeBin) {
      const auto lhs_key = state.particles.time_bin[lhs];
      const auto rhs_key = state.particles.time_bin[rhs];
      return std::tuple{lhs_key, lhs} < std::tuple{rhs_key, rhs};
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
    reorder_map.old_to_new_index[old_index] = static_cast<std::uint32_t>(new_index);
  }

  return reorder_map;
}

template <typename T>
void reorderAlignedVector(AlignedVector<T>& values, std::span<const std::uint32_t> new_to_old_index) {
  AlignedVector<T> reordered(values.size());
  for (std::size_t i = 0; i < new_to_old_index.size(); ++i) {
    reordered[i] = values[new_to_old_index[i]];
  }
  values.swap(reordered);
}

void reorderParticles(
    SimulationState& state,
    const ParticleReorderMap& reorder_map,
    const SidecarSyncPolicy& sync_policy) {
  if (!reorder_map.isConsistent(state.particles.size())) {
    throw std::invalid_argument("reorderParticles: inconsistent reorder map");
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

  auto remap_sidecar_index = [&](AlignedVector<std::uint32_t>& particle_index, SidecarSyncMode mode) {
    if (mode == SidecarSyncMode::kMoveWithParent) {
      reorderAlignedVector(particle_index, new_to_old_index);
      return;
    }
    for (auto& index : particle_index) {
      if (index >= reorder_map.old_to_new_index.size()) {
        throw std::out_of_range("reorderParticles: sidecar particle index out of range");
      }
      index = reorder_map.old_to_new_index[index];
    }
  };

  remap_sidecar_index(state.star_particles.particle_index, sync_policy.star_particles);
  remap_sidecar_index(state.black_holes.particle_index, sync_policy.black_holes);
  remap_sidecar_index(state.tracers.particle_index, sync_policy.tracers);

  state.rebuildSpeciesIndex();
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

}  // namespace cosmosim::core
