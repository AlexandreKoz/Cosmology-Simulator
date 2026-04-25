#include "cosmosim/core/simulation_state.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <new>
#include <numeric>
#include <stdexcept>
#include <tuple>
#include <vector>

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
  if (!state.particle_sidecar.gravity_softening_comoving.empty()) {
    reorderAlignedVector(state.particle_sidecar.gravity_softening_comoving, new_to_old_index);
  }

  auto remap_sidecar_index = [&](AlignedVector<std::uint32_t>& particle_index) {
    for (auto& index : particle_index) {
      if (index >= reorder_map.old_to_new_index.size()) {
        throw std::out_of_range("reorderParticles: sidecar particle index out of range");
      }
      index = reorder_map.old_to_new_index[index];
    }
  };

  if (sync_policy.star_particles == SidecarSyncMode::kMoveWithParent) {
    const auto row_order =
        buildSidecarRowOrderByParent(state.star_particles.particle_index, reorder_map.old_to_new_index);
    reorderAlignedVector(state.star_particles.particle_index, row_order);
    reorderAlignedVector(state.star_particles.formation_scale_factor, row_order);
    reorderAlignedVector(state.star_particles.birth_mass_code, row_order);
    reorderAlignedVector(state.star_particles.metallicity_mass_fraction, row_order);
    reorderAlignedVector(state.star_particles.stellar_age_years_last, row_order);
    reorderAlignedVector(state.star_particles.stellar_returned_mass_cumulative_code, row_order);
    reorderAlignedVector(state.star_particles.stellar_returned_metals_cumulative_code, row_order);
    reorderAlignedVector(state.star_particles.stellar_feedback_energy_cumulative_erg, row_order);
    for (std::size_t channel = 0; channel < state.star_particles.stellar_returned_mass_channel_cumulative_code.size();
         ++channel) {
      reorderAlignedVector(state.star_particles.stellar_returned_mass_channel_cumulative_code[channel], row_order);
      reorderAlignedVector(state.star_particles.stellar_returned_metals_channel_cumulative_code[channel], row_order);
      reorderAlignedVector(state.star_particles.stellar_feedback_energy_channel_cumulative_erg[channel], row_order);
    }
    remap_sidecar_index(state.star_particles.particle_index);
  } else {
    remap_sidecar_index(state.star_particles.particle_index);
  }
  if (sync_policy.black_holes == SidecarSyncMode::kMoveWithParent) {
    const auto row_order =
        buildSidecarRowOrderByParent(state.black_holes.particle_index, reorder_map.old_to_new_index);
    reorderAlignedVector(state.black_holes.particle_index, row_order);
    reorderAlignedVector(state.black_holes.host_cell_index, row_order);
    reorderAlignedVector(state.black_holes.subgrid_mass_code, row_order);
    reorderAlignedVector(state.black_holes.accretion_rate_code, row_order);
    reorderAlignedVector(state.black_holes.feedback_energy_code, row_order);
    reorderAlignedVector(state.black_holes.eddington_ratio, row_order);
    reorderAlignedVector(state.black_holes.cumulative_accreted_mass_code, row_order);
    reorderAlignedVector(state.black_holes.cumulative_feedback_energy_code, row_order);
    reorderAlignedVector(state.black_holes.duty_cycle_active_time_code, row_order);
    reorderAlignedVector(state.black_holes.duty_cycle_total_time_code, row_order);
    remap_sidecar_index(state.black_holes.particle_index);
  } else {
    remap_sidecar_index(state.black_holes.particle_index);
  }
  if (sync_policy.tracers == SidecarSyncMode::kMoveWithParent) {
    const auto row_order =
        buildSidecarRowOrderByParent(state.tracers.particle_index, reorder_map.old_to_new_index);
    reorderAlignedVector(state.tracers.particle_index, row_order);
    reorderAlignedVector(state.tracers.parent_particle_id, row_order);
    reorderAlignedVector(state.tracers.injection_step, row_order);
    reorderAlignedVector(state.tracers.host_cell_index, row_order);
    reorderAlignedVector(state.tracers.mass_fraction_of_host, row_order);
    reorderAlignedVector(state.tracers.last_host_mass_code, row_order);
    reorderAlignedVector(state.tracers.cumulative_exchanged_mass_code, row_order);
    remap_sidecar_index(state.tracers.particle_index);
  } else {
    remap_sidecar_index(state.tracers.particle_index);
  }

  state.rebuildSpeciesIndex();
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

}  // namespace cosmosim::core
