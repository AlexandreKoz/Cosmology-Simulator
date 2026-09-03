#include "cosmosim/core/simulation_state.hpp"

#include <limits>
#include <stdexcept>

namespace cosmosim::core {

void SimulationState::resizeParticles(std::size_t count) {
  (void)checkedLocalCount(count, kMaxLocalParticleCount, "particle", "SimulationState::resizeParticles");
  particles.resize(count);
  particle_sidecar.resize(count);
  bumpParticleIndexGeneration();
}

void SimulationState::resizeCells(std::size_t count) {
  (void)checkedLocalCount(count, kMaxLocalCellCount, "gas-cell", "SimulationState::resizeCells");
  cells.resize(count);
  gas_cells.resize(count);
  gas_cell_identity.clear();
  bumpCellIndexGeneration();
}

void SimulationState::resizePatches(std::size_t count) {
  (void)checkedLocalCount(count, kMaxLocalPatchCount, "AMR patch", "SimulationState::resizePatches");
  patches.resize(count);
}

bool SimulationState::validateOwnershipInvariants() const {
  if (!particles.isConsistent() || !particle_sidecar.isConsistent() || !cells.isConsistent() ||
      !gas_cells.isConsistent() || !patches.isConsistent() || !star_particles.isConsistent() ||
      !black_holes.isConsistent() || !tracers.isConsistent()) {
    return false;
  }

  if (particles.size() != particle_sidecar.size()) {
    return false;
  }

  if (!validateUniqueParticleIds()) {
    return false;
  }

  if (cells.size() != gas_cells.size()) {
    return false;
  }

  if (!gas_cell_identity.isConsistent() || !gas_cell_identity.coversDenseLocalRows(cells.size()) ||
      !gasCellIdentityMapMatchesSidecarLanes()) {
    return false;
  }

  if (!species.isConsistentWith(particle_sidecar)) {
    return false;
  }

  std::vector<std::uint8_t> star_rows_by_particle(particles.size(), 0);
  std::vector<std::uint8_t> bh_rows_by_particle(particles.size(), 0);
  std::vector<std::uint8_t> tracer_rows_by_particle(particles.size(), 0);

  if (patches.size() != 0U) {
    std::vector<std::uint32_t> cell_owner(cells.size(), std::numeric_limits<std::uint32_t>::max());
    for (std::size_t patch = 0; patch < patches.size(); ++patch) {
      const std::uint64_t begin = patches.first_cell[patch];
      const std::uint64_t count = patches.cell_count[patch];
      if (begin > cells.size() || count > static_cast<std::uint64_t>(cells.size()) - begin) {
        return false;
      }
      const std::uint64_t end = begin + count;
      for (std::uint64_t cell = begin; cell < end; ++cell) {
        const auto cell_index = static_cast<std::size_t>(cell);
        if (cell_owner[cell_index] != std::numeric_limits<std::uint32_t>::max()) {
          return false;
        }
        cell_owner[cell_index] = checkedLocalPatchRow(patch, "SimulationState::validateOwnershipInvariants patch row");
      }
    }
    for (std::size_t cell = 0; cell < cells.size(); ++cell) {
      if (cell_owner[cell] == std::numeric_limits<std::uint32_t>::max() ||
          cells.patch_index[cell] != cell_owner[cell]) {
        return false;
      }
    }
  } else if (cells.size() != 0U) {
    // Non-AMR hydro states are allowed to carry cells without PatchSoa ownership.
    // In that mode patch_index is not authoritative and must remain the neutral zero value.
    for (const std::uint32_t patch_index : cells.patch_index) {
      if (patch_index != 0U) {
        return false;
      }
    }
  }

  for (std::size_t i = 0; i < star_particles.size(); ++i) {
    const auto index = star_particles.particle_index[i];
    if (index >= particles.size()) {
      return false;
    }
    if (particle_sidecar.species_tag[index] != static_cast<std::uint32_t>(ParticleSpecies::kStar)) {
      return false;
    }
    if (++star_rows_by_particle[index] != 1) {
      return false;
    }
  }

  for (std::size_t i = 0; i < black_holes.size(); ++i) {
    const auto index = black_holes.particle_index[i];
    if (index >= particles.size()) {
      return false;
    }
    if (black_holes.host_cell_index[i] != kInvalidGasCellRow &&
        black_holes.host_cell_index[i] >= cells.size()) {
      return false;
    }
    if (particle_sidecar.species_tag[index] != static_cast<std::uint32_t>(ParticleSpecies::kBlackHole)) {
      return false;
    }
    if (++bh_rows_by_particle[index] != 1) {
      return false;
    }
  }

  for (std::size_t i = 0; i < tracers.size(); ++i) {
    const auto index = tracers.particle_index[i];
    if (index >= particles.size()) {
      return false;
    }
    if (tracers.host_cell_index[i] != kInvalidGasCellRow &&
        tracers.host_cell_index[i] >= cells.size()) {
      return false;
    }
    if (tracers.mass_fraction_of_host[i] < 0.0 || tracers.last_host_mass_code[i] < 0.0) {
      return false;
    }
    if (particle_sidecar.species_tag[index] != static_cast<std::uint32_t>(ParticleSpecies::kTracer)) {
      return false;
    }
    if (++tracer_rows_by_particle[index] != 1) {
      return false;
    }
  }

  for (std::size_t particle_index = 0; particle_index < particles.size(); ++particle_index) {
    const auto species_tag = particle_sidecar.species_tag[particle_index];
    const bool has_star_row = star_rows_by_particle[particle_index] == 1;
    const bool has_bh_row = bh_rows_by_particle[particle_index] == 1;
    const bool has_tracer_row = tracer_rows_by_particle[particle_index] == 1;

    if (species_tag == static_cast<std::uint32_t>(ParticleSpecies::kStar)) {
      if (!has_star_row || has_bh_row || has_tracer_row) {
        return false;
      }
    } else if (species_tag == static_cast<std::uint32_t>(ParticleSpecies::kBlackHole)) {
      if (has_star_row || !has_bh_row || has_tracer_row) {
        return false;
      }
    } else if (species_tag == static_cast<std::uint32_t>(ParticleSpecies::kTracer)) {
      if (has_star_row || has_bh_row || !has_tracer_row) {
        return false;
      }
    } else if (has_star_row || has_bh_row || has_tracer_row) {
      return false;
    }
  }

  return true;
}

std::uint64_t SimulationState::particleIndexGeneration() const noexcept {
  return m_particle_index_generation;
}

std::uint64_t SimulationState::cellIndexGeneration() const noexcept {
  return m_cell_index_generation;
}

std::uint64_t SimulationState::gravitySourceGeneration() const noexcept {
  return m_gravity_source_generation;
}

void SimulationState::bumpGravitySourceGeneration() noexcept {
  if (m_gravity_source_generation != std::numeric_limits<std::uint64_t>::max()) {
    ++m_gravity_source_generation;
  }
}

void SimulationState::bumpParticleIndexGeneration() noexcept {
  ++m_particle_index_generation;
}

void SimulationState::bumpCellIndexGeneration() noexcept {
  ++m_cell_index_generation;
}

}  // namespace cosmosim::core
