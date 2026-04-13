#include "cosmosim/core/simulation_state.hpp"

namespace cosmosim::core {

void SimulationState::resizeParticles(std::size_t count) {
  particles.resize(count);
  particle_sidecar.resize(count);
}

void SimulationState::resizeCells(std::size_t count) {
  cells.resize(count);
  gas_cells.resize(count);
}

void SimulationState::resizePatches(std::size_t count) { patches.resize(count); }

bool SimulationState::validateOwnershipInvariants() const {
  if (!particles.isConsistent() || !particle_sidecar.isConsistent() || !cells.isConsistent() ||
      !gas_cells.isConsistent() || !patches.isConsistent() || !star_particles.isConsistent() ||
      !black_holes.isConsistent() || !tracers.isConsistent()) {
    return false;
  }

  if (particles.size() != particle_sidecar.size()) {
    return false;
  }

  if (cells.size() != gas_cells.size()) {
    return false;
  }

  if (!species.isConsistentWith(particle_sidecar)) {
    return false;
  }

  for (std::size_t i = 0; i < cells.patch_index.size(); ++i) {
    if (cells.patch_index[i] >= patches.size() && patches.size() > 0) {
      return false;
    }
  }

  for (std::size_t patch = 0; patch < patches.size(); ++patch) {
    const std::uint64_t begin = patches.first_cell[patch];
    const std::uint64_t count = patches.cell_count[patch];
    if (begin + count > cells.size()) {
      return false;
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
  }

  for (std::size_t i = 0; i < black_holes.size(); ++i) {
    const auto index = black_holes.particle_index[i];
    if (index >= particles.size()) {
      return false;
    }
    if (particle_sidecar.species_tag[index] != static_cast<std::uint32_t>(ParticleSpecies::kBlackHole)) {
      return false;
    }
  }

  for (std::size_t i = 0; i < tracers.size(); ++i) {
    const auto index = tracers.particle_index[i];
    if (index >= particles.size()) {
      return false;
    }
    if (tracers.host_cell_index[i] >= cells.size() && cells.size() > 0) {
      return false;
    }
    if (tracers.mass_fraction_of_host[i] < 0.0 || tracers.last_host_mass_code[i] < 0.0) {
      return false;
    }
    if (particle_sidecar.species_tag[index] != static_cast<std::uint32_t>(ParticleSpecies::kTracer)) {
      return false;
    }
  }

  return true;
}

}  // namespace cosmosim::core
