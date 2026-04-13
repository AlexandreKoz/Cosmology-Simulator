#include "cosmosim/core/simulation_state.hpp"

#include <array>
#include <stdexcept>
#include <unordered_set>

namespace cosmosim::core {
namespace {
constexpr std::size_t k_species_count = 5;

[[nodiscard]] bool isValidSpeciesTag(std::uint32_t value) {
  return value <= static_cast<std::uint32_t>(ParticleSpecies::kTracer);
}
}  // namespace

std::uint64_t SpeciesContainer::totalCount() const noexcept {
  std::uint64_t total = 0;
  for (const auto count : count_by_species) {
    total += count;
  }
  return total;
}

bool SpeciesContainer::isConsistentWith(const ParticleSidecar& sidecar) const noexcept {
  std::array<std::uint64_t, k_species_count> measured{};
  for (const auto tag : sidecar.species_tag) {
    if (!isValidSpeciesTag(tag)) {
      return false;
    }
    ++measured.at(tag);
  }
  return measured == count_by_species;
}

void ParticleSpeciesIndex::rebuild(const ParticleSidecar& sidecar) {
  for (auto& indices : global_index_by_species) {
    indices.clear();
  }

  local_index_by_global.resize(sidecar.size());
  for (std::size_t global_index = 0; global_index < sidecar.size(); ++global_index) {
    const auto tag = sidecar.species_tag[global_index];
    if (!isValidSpeciesTag(tag)) {
      throw std::invalid_argument("ParticleSpeciesIndex.rebuild: invalid species tag");
    }
    auto& species_indices = global_index_by_species[tag];
    local_index_by_global[global_index] = static_cast<std::uint32_t>(species_indices.size());
    species_indices.push_back(static_cast<std::uint32_t>(global_index));
  }
}

std::size_t ParticleSpeciesIndex::count(ParticleSpecies species) const noexcept {
  return global_index_by_species[static_cast<std::uint32_t>(species)].size();
}

std::span<const std::uint32_t> ParticleSpeciesIndex::globalIndices(ParticleSpecies species) const noexcept {
  return global_index_by_species[static_cast<std::uint32_t>(species)];
}

std::uint32_t ParticleSpeciesIndex::localIndex(std::uint32_t global_index) const {
  if (global_index >= local_index_by_global.size()) {
    throw std::out_of_range("ParticleSpeciesIndex.localIndex: global index out of range");
  }
  return local_index_by_global[global_index];
}

std::uint32_t ParticleSpeciesIndex::globalIndex(ParticleSpecies species, std::uint32_t local_index) const {
  const auto& species_indices = global_index_by_species[static_cast<std::uint32_t>(species)];
  if (local_index >= species_indices.size()) {
    throw std::out_of_range("ParticleSpeciesIndex.globalIndex: local index out of range");
  }
  return species_indices[local_index];
}

void SimulationState::rebuildSpeciesIndex() { particle_species_index.rebuild(particle_sidecar); }

bool SimulationState::validateUniqueParticleIds() const {
  std::unordered_set<std::uint64_t> ids;
  ids.reserve(particle_sidecar.particle_id.size());
  for (const auto id : particle_sidecar.particle_id) {
    if (!ids.insert(id).second) {
      return false;
    }
  }
  return true;
}

ParticleTransferPacket SimulationState::packSpeciesTransferPacket(ParticleSpecies species_tag) const {
  const auto indices = particle_species_index.globalIndices(species_tag);

  ParticleTransferPacket packet;
  packet.species = species_tag;
  packet.particle_id.resize(indices.size());
  packet.position_x_comoving.resize(indices.size());
  packet.position_y_comoving.resize(indices.size());
  packet.position_z_comoving.resize(indices.size());
  packet.velocity_x_peculiar.resize(indices.size());
  packet.velocity_y_peculiar.resize(indices.size());
  packet.velocity_z_peculiar.resize(indices.size());
  packet.mass_code.resize(indices.size());
  packet.time_bin.resize(indices.size());
  packet.owning_rank.resize(indices.size());

  for (std::size_t i = 0; i < indices.size(); ++i) {
    const auto source = indices[i];
    packet.particle_id[i] = particle_sidecar.particle_id[source];
    packet.position_x_comoving[i] = particles.position_x_comoving[source];
    packet.position_y_comoving[i] = particles.position_y_comoving[source];
    packet.position_z_comoving[i] = particles.position_z_comoving[source];
    packet.velocity_x_peculiar[i] = particles.velocity_x_peculiar[source];
    packet.velocity_y_peculiar[i] = particles.velocity_y_peculiar[source];
    packet.velocity_z_peculiar[i] = particles.velocity_z_peculiar[source];
    packet.mass_code[i] = particles.mass_code[source];
    packet.time_bin[i] = particles.time_bin[source];
    packet.owning_rank[i] = particle_sidecar.owning_rank[source];
  }

  return packet;
}

}  // namespace cosmosim::core
