#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <vector>

#include "cosmosim/core/checked_arithmetic.hpp"
#include "cosmosim/core/memory_governor.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/hydro/hydro_core_solver.hpp"
#include "cosmosim/parallel/distributed_memory.hpp"

namespace cosmosim::workflows::internal {

// Conservative per-particle persistent footprint used by runtime decomposition
// memory accounting (species lanes plus required module sidecar rows). This is
// the runtime module-sidecar-aware estimator; the startup initial-placement
// weighting keeps its own historical formula so startup weights do not shift.
[[nodiscard]] inline std::uint64_t estimateParticleMemoryBytesForDecomposition(
    const core::SimulationState& state,
    std::uint32_t species_tag) {
  std::uint64_t bytes = sizeof(double) * 7U + sizeof(std::uint64_t) * 2U + sizeof(std::uint32_t) * 3U;
  if (!state.particle_sidecar.gravity_softening_comoving.empty()) {
    bytes += sizeof(double);
  }
  if (!state.particle_sidecar.has_gravity_softening_override.empty()) {
    bytes += sizeof(std::uint8_t);
  }
  if (species_tag == static_cast<std::uint32_t>(core::ParticleSpecies::kGas)) {
    bytes += sizeof(double) * 8U + sizeof(std::uint64_t) * 2U + sizeof(std::uint32_t) * 2U;
  } else if (species_tag == static_cast<std::uint32_t>(core::ParticleSpecies::kStar)) {
    bytes += sizeof(std::uint32_t) + sizeof(double) * 13U;
  } else if (species_tag == static_cast<std::uint32_t>(core::ParticleSpecies::kBlackHole)) {
    bytes += sizeof(std::uint32_t) * 2U + sizeof(double) * 8U;
  } else if (species_tag == static_cast<std::uint32_t>(core::ParticleSpecies::kTracer)) {
    bytes += sizeof(std::uint64_t) * 2U + sizeof(std::uint32_t) * 2U + sizeof(double) * 3U;
  }
  for (const core::ModuleSidecarBlock* block_ptr : state.sidecars.blocksSortedByName()) {
    const core::ModuleSidecarBlock& block = *block_ptr;
    if (!block.particle_indexed || block.row_stride_bytes == 0U) {
      continue;
    }
    const bool species_mask_requires_row = (block.required_species_mask & (1U << species_tag)) != 0U ||
        (block.requirement.kind == core::ModuleSidecarRequirementKind::kSpeciesMask &&
         (block.requirement.species_mask & (1U << species_tag)) != 0U);
    const bool predicate_may_require_row =
        (block.requirement.kind == core::ModuleSidecarRequirementKind::kGasDensityAtLeast &&
         species_tag == static_cast<std::uint32_t>(core::ParticleSpecies::kGas)) ||
        (block.requirement.kind == core::ModuleSidecarRequirementKind::kBlackHoleAccretionAtLeast &&
         species_tag == static_cast<std::uint32_t>(core::ParticleSpecies::kBlackHole)) ||
        block.requirement.kind == core::ModuleSidecarRequirementKind::kParticleFlagMask;
    if (species_mask_requires_row || predicate_may_require_row) {
      bytes += block.row_stride_bytes;
    }
  }
  return bytes;
}

// Concrete storage behind parallel::RuntimeDecompositionSourceView. Builds the
// derived active mask, gas particle-to-patch incidence, and compact nonempty
// patch index from the canonical SimulationState; all other view spans alias
// authoritative state lanes. Gas-empty/DMO states allocate empty offset/index
// spans through parallel::runtimeDecompositionHasGasIncidenceSource — the same
// condition parallel::estimateRuntimeDecompositionSourceStorage charges for
// MemoryGovernor admission.
class RuntimeDecompositionSourceStorage {
 public:
  RuntimeDecompositionSourceStorage(
      const core::SimulationState& state,
      int world_rank,
      std::span<const std::uint32_t> active_particle_indices)
      : m_state(state), m_world_rank(world_rank) {
    if (!state.patches.isConsistent()) {
      throw std::runtime_error(
          "runtime decomposition source requires consistent AMR patch metadata");
    }
    if (!active_particle_indices.empty()) {
      m_active_mask.assign(state.particles.size(), 0U);
      for (const std::uint32_t particle_index : active_particle_indices) {
        if (particle_index >= m_active_mask.size()) {
          throw std::out_of_range(
              "runtime decomposition active particle index is outside SimulationState");
        }
        m_active_mask[particle_index] = 1U;
      }
    }
    struct GasPatchEntry {
      std::uint64_t particle_id = 0;
      std::uint32_t patch_index = 0;
    };
    std::vector<GasPatchEntry> gas_patch_entries;
    const bool has_gas_incidence_source =
        parallel::runtimeDecompositionHasGasIncidenceSource(state);
    if (has_gas_incidence_source) {
      state.gas_cell_identity.requireCoversDenseLocalRows(
          state.cells.size(), "runtime decomposition source gas identity");
      gas_patch_entries.reserve(state.cells.size());
      for (std::size_t cell = 0; cell < state.cells.size(); ++cell) {
        const core::GasCellIdentityRecord* record =
            state.gas_cell_identity.findByLocalRow(
                core::checkedLocalCellRow(
                    cell, "runtime decomposition source gas identity row"));
        if (record == nullptr || !record->parent_particle_id.has_value() ||
            record->parent_particle_id.value() == 0U) {
          continue;
        }
        gas_patch_entries.push_back(GasPatchEntry{
            .particle_id = record->parent_particle_id.value(),
            .patch_index = state.cells.patch_index[cell],
        });
      }
    }
    std::sort(gas_patch_entries.begin(), gas_patch_entries.end(),
              [](const GasPatchEntry& lhs, const GasPatchEntry& rhs) {
                if (lhs.particle_id != rhs.particle_id) {
                  return lhs.particle_id < rhs.particle_id;
                }
                return lhs.patch_index < rhs.patch_index;
              });
    gas_patch_entries.erase(
        std::unique(
            gas_patch_entries.begin(), gas_patch_entries.end(),
            [](const GasPatchEntry& lhs, const GasPatchEntry& rhs) {
              return lhs.particle_id == rhs.particle_id &&
                  lhs.patch_index == rhs.patch_index;
            }),
        gas_patch_entries.end());
    // Empty gas incidence (pure DMO / gas-empty state) is represented by
    // empty offset/index spans; no O(N_particles) all-zero offset array.
    if (!gas_patch_entries.empty()) {
      m_gas_patch_list_offsets.assign(state.particles.size() + 1U, 0U);
      m_gas_patch_indices.reserve(gas_patch_entries.size());
      for (std::size_t particle = 0; particle < state.particles.size(); ++particle) {
        const std::uint64_t particle_id = state.particle_sidecar.particle_id[particle];
        const auto begin = std::lower_bound(
            gas_patch_entries.begin(), gas_patch_entries.end(), particle_id,
            [](const GasPatchEntry& entry, std::uint64_t value) {
              return entry.particle_id < value;
            });
        for (auto it = begin;
             it != gas_patch_entries.end() && it->particle_id == particle_id;
             ++it) {
          m_gas_patch_indices.push_back(it->patch_index);
        }
        m_gas_patch_list_offsets[particle + 1U] =
            static_cast<std::uint32_t>(m_gas_patch_indices.size());
      }
    }
    m_compact_patch_indices.reserve(state.patches.size());
    for (std::size_t patch = 0; patch < state.patches.size(); ++patch) {
      if (state.patches.cell_count[patch] == 0U) {
        continue;
      }
      if (state.patches.first_cell[patch] >= state.cells.size()) {
        throw std::out_of_range(
            "runtime decomposition nonempty AMR patch has no authoritative cell range");
      }
      m_compact_patch_indices.push_back(static_cast<std::uint32_t>(patch));
    }
    for (std::size_t species = 0; species < core::k_particle_species_count; ++species) {
      m_particle_memory_bytes_by_species[species] =
          estimateParticleMemoryBytesForDecomposition(
              state, static_cast<std::uint32_t>(species));
    }
    m_source_scratch_bytes = 0U;
    m_source_scratch_bytes = core::checkedMemoryBytesAdd(
        m_source_scratch_bytes,
        core::checkedSizeMultiply(
            m_active_mask.capacity(), sizeof(std::uint8_t),
            "runtime decomposition active mask capacity"),
        "runtime decomposition source scratch byte overflow");
    m_source_scratch_bytes = core::checkedMemoryBytesAdd(
        m_source_scratch_bytes,
        core::checkedSizeMultiply(
            m_gas_patch_list_offsets.capacity(), sizeof(std::uint32_t),
            "runtime decomposition gas patch offset capacity"),
        "runtime decomposition source scratch byte overflow");
    m_source_scratch_bytes = core::checkedMemoryBytesAdd(
        m_source_scratch_bytes,
        core::checkedSizeMultiply(
            m_gas_patch_indices.capacity(), sizeof(std::uint32_t),
            "runtime decomposition gas patch incidence capacity"),
        "runtime decomposition source scratch byte overflow");
     m_source_scratch_bytes = core::checkedMemoryBytesAdd(
         m_source_scratch_bytes,
         core::checkedSizeMultiply(
             m_compact_patch_indices.capacity(), sizeof(std::uint32_t),
             "runtime decomposition compact patch capacity"),
         "runtime decomposition source scratch byte overflow");
     if (has_gas_incidence_source) {
       m_source_scratch_bytes = core::checkedMemoryBytesAdd(
           m_source_scratch_bytes,
           core::checkedSizeMultiply(
               m_state.cells.size(), sizeof(std::uint64_t) * 2U,
               "runtime decomposition gas incidence construction capacity"),
           "runtime decomposition source scratch byte overflow");
     }

  }

  [[nodiscard]] parallel::RuntimeDecompositionSourceView view() const noexcept {
    parallel::RuntimeDecompositionSourceView source;
    source.world_rank = m_world_rank;
    source.particle_count = m_state.particles.size();
    source.patch_count = m_state.patches.size();
    source.particle_ids = m_state.particle_sidecar.particle_id;
    source.particle_x_comoving = m_state.particles.position_x_comoving;
    source.particle_y_comoving = m_state.particles.position_y_comoving;
    source.particle_z_comoving = m_state.particles.position_z_comoving;
    source.particle_species_tag = m_state.particle_sidecar.species_tag;
    source.particle_owning_rank = m_state.particle_sidecar.owning_rank;
    source.active_particle_mask = m_active_mask;
    source.patch_ids = m_state.patches.patch_id;
    source.patch_levels = m_state.patches.level;
    source.patch_owning_rank = m_state.patches.owning_rank;
    source.patch_first_cells = m_state.patches.first_cell;
    source.patch_cell_counts = m_state.patches.cell_count;
    source.compact_patch_indices = m_compact_patch_indices;
    source.patch_cell_dim_x = m_state.patches.cell_dim_x;
    source.patch_cell_dim_y = m_state.patches.cell_dim_y;
    source.patch_cell_dim_z = m_state.patches.cell_dim_z;
    source.cell_x_comoving = m_state.cells.center_x_comoving;
    source.cell_y_comoving = m_state.cells.center_y_comoving;
    source.cell_z_comoving = m_state.cells.center_z_comoving;
    source.cell_patch_indices = m_state.cells.patch_index;
    source.gas_patch_list_offsets = m_gas_patch_list_offsets;
    source.gas_patch_indices = m_gas_patch_indices;
    source.particle_memory_bytes_by_species = m_particle_memory_bytes_by_species;
    source.gas_species_tag = static_cast<std::uint32_t>(core::ParticleSpecies::kGas);
    source.star_species_tag = static_cast<std::uint32_t>(core::ParticleSpecies::kStar);
    source.black_hole_species_tag =
        static_cast<std::uint32_t>(core::ParticleSpecies::kBlackHole);
    source.gas_transient_memory_bytes_per_cell =
        hydro::k_hydro_runtime_batch_scratch_budget_bytes_per_cell;
    source.source_scratch_bytes = m_source_scratch_bytes;
    return source;
  }

 private:
  const core::SimulationState& m_state;
  int m_world_rank = 0;
   std::vector<std::uint8_t> m_active_mask;
   std::vector<std::uint32_t> m_gas_patch_list_offsets;
   std::vector<std::uint32_t> m_gas_patch_indices;

  std::vector<std::uint32_t> m_compact_patch_indices;
   std::array<std::uint64_t, 5U> m_particle_memory_bytes_by_species{};

  std::uint64_t m_source_scratch_bytes = 0U;
};

}  // namespace cosmosim::workflows::internal
