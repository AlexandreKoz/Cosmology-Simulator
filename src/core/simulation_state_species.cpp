#include "cosmosim/core/simulation_state.hpp"

#include <array>
#include <algorithm>
#include <limits>
#include <stdexcept>
#include <unordered_set>

namespace cosmosim::core {
namespace {
constexpr std::size_t k_species_count = 5;

[[nodiscard]] bool isValidSpeciesTag(std::uint32_t value) {
  return value <= static_cast<std::uint32_t>(ParticleSpecies::kTracer);
}

[[nodiscard]] std::size_t findRequiredSidecarRow(
    std::span<const std::uint32_t> particle_indices,
    std::uint32_t particle_index,
    const char* sidecar_name) {
  for (std::size_t row = 0; row < particle_indices.size(); ++row) {
    if (particle_indices[row] == particle_index) {
      return row;
    }
  }
  throw std::runtime_error(std::string("missing required sidecar row for particle migration: ") + sidecar_name);
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

std::vector<ParticleMigrationRecord> SimulationState::packParticleMigrationRecords(
    std::span<const std::uint32_t> local_indices) const {
  std::vector<ParticleMigrationRecord> records;
  records.reserve(local_indices.size());
  for (const auto index : local_indices) {
    if (index >= particles.size()) {
      throw std::out_of_range("packParticleMigrationRecords: local index out of range");
    }
  }

  for (const auto index : local_indices) {
    ParticleMigrationRecord record;
    record.particle_id = particle_sidecar.particle_id[index];
    record.sfc_key = particle_sidecar.sfc_key[index];
    record.species_tag = particle_sidecar.species_tag[index];
    record.particle_flags = particle_sidecar.particle_flags[index];
    record.owning_rank = particle_sidecar.owning_rank[index];
    record.position_x_comoving = particles.position_x_comoving[index];
    record.position_y_comoving = particles.position_y_comoving[index];
    record.position_z_comoving = particles.position_z_comoving[index];
    record.velocity_x_peculiar = particles.velocity_x_peculiar[index];
    record.velocity_y_peculiar = particles.velocity_y_peculiar[index];
    record.velocity_z_peculiar = particles.velocity_z_peculiar[index];
    record.mass_code = particles.mass_code[index];
    record.time_bin = particles.time_bin[index];

    if (record.species_tag == static_cast<std::uint32_t>(ParticleSpecies::kStar)) {
      const std::size_t row =
          findRequiredSidecarRow(star_particles.particle_index, index, "star_particles");
      record.has_star_fields = true;
      record.star_fields.formation_scale_factor = star_particles.formation_scale_factor[row];
      record.star_fields.birth_mass_code = star_particles.birth_mass_code[row];
      record.star_fields.metallicity_mass_fraction = star_particles.metallicity_mass_fraction[row];
      record.star_fields.stellar_age_years_last = star_particles.stellar_age_years_last[row];
      record.star_fields.stellar_returned_mass_cumulative_code =
          star_particles.stellar_returned_mass_cumulative_code[row];
      record.star_fields.stellar_returned_metals_cumulative_code =
          star_particles.stellar_returned_metals_cumulative_code[row];
      record.star_fields.stellar_feedback_energy_cumulative_erg =
          star_particles.stellar_feedback_energy_cumulative_erg[row];
      for (std::size_t channel = 0; channel < record.star_fields.stellar_returned_mass_channel_cumulative_code.size();
           ++channel) {
        record.star_fields.stellar_returned_mass_channel_cumulative_code[channel] =
            star_particles.stellar_returned_mass_channel_cumulative_code[channel][row];
        record.star_fields.stellar_returned_metals_channel_cumulative_code[channel] =
            star_particles.stellar_returned_metals_channel_cumulative_code[channel][row];
        record.star_fields.stellar_feedback_energy_channel_cumulative_erg[channel] =
            star_particles.stellar_feedback_energy_channel_cumulative_erg[channel][row];
      }
    } else if (record.species_tag == static_cast<std::uint32_t>(ParticleSpecies::kBlackHole)) {
      const std::size_t row =
          findRequiredSidecarRow(black_holes.particle_index, index, "black_holes");
      record.has_black_hole_fields = true;
      record.black_hole_fields.host_cell_index = black_holes.host_cell_index[row];
      record.black_hole_fields.subgrid_mass_code = black_holes.subgrid_mass_code[row];
      record.black_hole_fields.accretion_rate_code = black_holes.accretion_rate_code[row];
      record.black_hole_fields.feedback_energy_code = black_holes.feedback_energy_code[row];
      record.black_hole_fields.eddington_ratio = black_holes.eddington_ratio[row];
      record.black_hole_fields.cumulative_accreted_mass_code = black_holes.cumulative_accreted_mass_code[row];
      record.black_hole_fields.cumulative_feedback_energy_code = black_holes.cumulative_feedback_energy_code[row];
      record.black_hole_fields.duty_cycle_active_time_code = black_holes.duty_cycle_active_time_code[row];
      record.black_hole_fields.duty_cycle_total_time_code = black_holes.duty_cycle_total_time_code[row];
    } else if (record.species_tag == static_cast<std::uint32_t>(ParticleSpecies::kTracer)) {
      const std::size_t row = findRequiredSidecarRow(tracers.particle_index, index, "tracers");
      record.has_tracer_fields = true;
      record.tracer_fields.parent_particle_id = tracers.parent_particle_id[row];
      record.tracer_fields.injection_step = tracers.injection_step[row];
      record.tracer_fields.host_cell_index = tracers.host_cell_index[row];
      record.tracer_fields.mass_fraction_of_host = tracers.mass_fraction_of_host[row];
      record.tracer_fields.last_host_mass_code = tracers.last_host_mass_code[row];
      record.tracer_fields.cumulative_exchanged_mass_code = tracers.cumulative_exchanged_mass_code[row];
    }

    records.push_back(record);
  }
  return records;
}

void SimulationState::commitParticleMigration(const ParticleMigrationCommit& commit) {
  if (commit.world_rank < 0) {
    throw std::invalid_argument("commitParticleMigration: world_rank must be non-negative");
  }
  const std::size_t particle_count = particles.size();
  std::vector<std::uint8_t> outbound_mask(particle_count, 0U);
  for (const auto index : commit.outbound_local_indices) {
    if (index >= particle_count) {
      throw std::out_of_range("commitParticleMigration: outbound index out of range");
    }
    if (outbound_mask[index] != 0U) {
      throw std::invalid_argument("commitParticleMigration: duplicate outbound index");
    }
    if (particle_sidecar.owning_rank[index] != static_cast<std::uint32_t>(commit.world_rank)) {
      throw std::invalid_argument("commitParticleMigration: outbound index is not owned by committing rank");
    }
    outbound_mask[index] = 1U;
  }

  std::vector<std::uint8_t> stale_ghost_mask(particle_count, 0U);
  for (const auto index : commit.stale_local_ghost_indices) {
    if (index >= particle_count) {
      throw std::out_of_range("commitParticleMigration: stale ghost index out of range");
    }
    stale_ghost_mask[index] = 1U;
  }

  std::vector<std::uint8_t> remove_mask = outbound_mask;
  for (std::size_t i = 0; i < remove_mask.size(); ++i) {
    if (stale_ghost_mask[i] != 0U) {
      remove_mask[i] = 1U;
    }
  }

  std::vector<std::uint32_t> old_to_new(particle_count, std::numeric_limits<std::uint32_t>::max());
  std::size_t kept_count = 0;
  for (std::size_t i = 0; i < particle_count; ++i) {
    if (remove_mask[i] == 0U) {
      old_to_new[i] = static_cast<std::uint32_t>(kept_count++);
    }
  }

  const std::size_t final_count = kept_count + commit.inbound_records.size();
  ParticleSoa new_particles;
  new_particles.resize(final_count);
  ParticleSidecar new_sidecar;
  new_sidecar.resize(final_count);

  std::size_t write_index = 0;
  for (std::size_t i = 0; i < particle_count; ++i) {
    if (remove_mask[i] != 0U) {
      continue;
    }
    new_particles.position_x_comoving[write_index] = particles.position_x_comoving[i];
    new_particles.position_y_comoving[write_index] = particles.position_y_comoving[i];
    new_particles.position_z_comoving[write_index] = particles.position_z_comoving[i];
    new_particles.velocity_x_peculiar[write_index] = particles.velocity_x_peculiar[i];
    new_particles.velocity_y_peculiar[write_index] = particles.velocity_y_peculiar[i];
    new_particles.velocity_z_peculiar[write_index] = particles.velocity_z_peculiar[i];
    new_particles.mass_code[write_index] = particles.mass_code[i];
    new_particles.time_bin[write_index] = particles.time_bin[i];
    new_sidecar.particle_id[write_index] = particle_sidecar.particle_id[i];
    new_sidecar.sfc_key[write_index] = particle_sidecar.sfc_key[i];
    new_sidecar.species_tag[write_index] = particle_sidecar.species_tag[i];
    new_sidecar.particle_flags[write_index] = particle_sidecar.particle_flags[i];
    new_sidecar.owning_rank[write_index] = particle_sidecar.owning_rank[i];
    ++write_index;
  }

  for (const auto& inbound : commit.inbound_records) {
    if (!isValidSpeciesTag(inbound.species_tag)) {
      throw std::invalid_argument("commitParticleMigration: inbound record has invalid species tag");
    }
    if (inbound.owning_rank != static_cast<std::uint32_t>(commit.world_rank)) {
      throw std::invalid_argument("commitParticleMigration: inbound record ownership must equal commit world rank");
    }
    new_particles.position_x_comoving[write_index] = inbound.position_x_comoving;
    new_particles.position_y_comoving[write_index] = inbound.position_y_comoving;
    new_particles.position_z_comoving[write_index] = inbound.position_z_comoving;
    new_particles.velocity_x_peculiar[write_index] = inbound.velocity_x_peculiar;
    new_particles.velocity_y_peculiar[write_index] = inbound.velocity_y_peculiar;
    new_particles.velocity_z_peculiar[write_index] = inbound.velocity_z_peculiar;
    new_particles.mass_code[write_index] = inbound.mass_code;
    new_particles.time_bin[write_index] = inbound.time_bin;
    new_sidecar.particle_id[write_index] = inbound.particle_id;
    new_sidecar.sfc_key[write_index] = inbound.sfc_key;
    new_sidecar.species_tag[write_index] = inbound.species_tag;
    new_sidecar.particle_flags[write_index] = inbound.particle_flags;
    new_sidecar.owning_rank[write_index] = inbound.owning_rank;
    ++write_index;
  }

  auto rebuildStarSidecar = [&](StarParticleSidecar* destination) {
    std::vector<std::uint32_t> kept_rows;
    for (std::size_t row = 0; row < star_particles.size(); ++row) {
      const auto old_particle = star_particles.particle_index[row];
      if (old_particle >= remove_mask.size()) {
        throw std::runtime_error("commitParticleMigration: stale star sidecar index");
      }
      if (remove_mask[old_particle] == 0U) {
        kept_rows.push_back(static_cast<std::uint32_t>(row));
      }
    }
    destination->resize(kept_rows.size());
    for (std::size_t row = 0; row < kept_rows.size(); ++row) {
      const std::size_t source = kept_rows[row];
      const auto old_particle = star_particles.particle_index[source];
      destination->particle_index[row] = old_to_new[old_particle];
      destination->formation_scale_factor[row] = star_particles.formation_scale_factor[source];
      destination->birth_mass_code[row] = star_particles.birth_mass_code[source];
      destination->metallicity_mass_fraction[row] = star_particles.metallicity_mass_fraction[source];
      destination->stellar_age_years_last[row] = star_particles.stellar_age_years_last[source];
      destination->stellar_returned_mass_cumulative_code[row] =
          star_particles.stellar_returned_mass_cumulative_code[source];
      destination->stellar_returned_metals_cumulative_code[row] =
          star_particles.stellar_returned_metals_cumulative_code[source];
      destination->stellar_feedback_energy_cumulative_erg[row] =
          star_particles.stellar_feedback_energy_cumulative_erg[source];
      for (std::size_t channel = 0; channel < destination->stellar_returned_mass_channel_cumulative_code.size(); ++channel) {
        destination->stellar_returned_mass_channel_cumulative_code[channel][row] =
            star_particles.stellar_returned_mass_channel_cumulative_code[channel][source];
        destination->stellar_returned_metals_channel_cumulative_code[channel][row] =
            star_particles.stellar_returned_metals_channel_cumulative_code[channel][source];
        destination->stellar_feedback_energy_channel_cumulative_erg[channel][row] =
            star_particles.stellar_feedback_energy_channel_cumulative_erg[channel][source];
      }
    }
    for (const auto& inbound : commit.inbound_records) {
      if (!inbound.has_star_fields) {
        continue;
      }
      const std::size_t row = destination->size();
      destination->resize(row + 1);
      destination->particle_index[row] = static_cast<std::uint32_t>(kept_count + (&inbound - commit.inbound_records.data()));
      destination->formation_scale_factor[row] = inbound.star_fields.formation_scale_factor;
      destination->birth_mass_code[row] = inbound.star_fields.birth_mass_code;
      destination->metallicity_mass_fraction[row] = inbound.star_fields.metallicity_mass_fraction;
      destination->stellar_age_years_last[row] = inbound.star_fields.stellar_age_years_last;
      destination->stellar_returned_mass_cumulative_code[row] =
          inbound.star_fields.stellar_returned_mass_cumulative_code;
      destination->stellar_returned_metals_cumulative_code[row] =
          inbound.star_fields.stellar_returned_metals_cumulative_code;
      destination->stellar_feedback_energy_cumulative_erg[row] =
          inbound.star_fields.stellar_feedback_energy_cumulative_erg;
      for (std::size_t channel = 0; channel < inbound.star_fields.stellar_returned_mass_channel_cumulative_code.size();
           ++channel) {
        destination->stellar_returned_mass_channel_cumulative_code[channel][row] =
            inbound.star_fields.stellar_returned_mass_channel_cumulative_code[channel];
        destination->stellar_returned_metals_channel_cumulative_code[channel][row] =
            inbound.star_fields.stellar_returned_metals_channel_cumulative_code[channel];
        destination->stellar_feedback_energy_channel_cumulative_erg[channel][row] =
            inbound.star_fields.stellar_feedback_energy_channel_cumulative_erg[channel];
      }
    }
  };

  auto rebuildBlackHoleSidecar = [&](BlackHoleParticleSidecar* destination) {
    std::vector<std::uint32_t> kept_rows;
    for (std::size_t row = 0; row < black_holes.size(); ++row) {
      const auto old_particle = black_holes.particle_index[row];
      if (old_particle >= remove_mask.size()) {
        throw std::runtime_error("commitParticleMigration: stale black-hole sidecar index");
      }
      if (remove_mask[old_particle] == 0U) {
        kept_rows.push_back(static_cast<std::uint32_t>(row));
      }
    }
    destination->resize(kept_rows.size());
    for (std::size_t row = 0; row < kept_rows.size(); ++row) {
      const std::size_t source = kept_rows[row];
      const auto old_particle = black_holes.particle_index[source];
      destination->particle_index[row] = old_to_new[old_particle];
      destination->host_cell_index[row] = black_holes.host_cell_index[source];
      destination->subgrid_mass_code[row] = black_holes.subgrid_mass_code[source];
      destination->accretion_rate_code[row] = black_holes.accretion_rate_code[source];
      destination->feedback_energy_code[row] = black_holes.feedback_energy_code[source];
      destination->eddington_ratio[row] = black_holes.eddington_ratio[source];
      destination->cumulative_accreted_mass_code[row] = black_holes.cumulative_accreted_mass_code[source];
      destination->cumulative_feedback_energy_code[row] = black_holes.cumulative_feedback_energy_code[source];
      destination->duty_cycle_active_time_code[row] = black_holes.duty_cycle_active_time_code[source];
      destination->duty_cycle_total_time_code[row] = black_holes.duty_cycle_total_time_code[source];
    }
    for (const auto& inbound : commit.inbound_records) {
      if (!inbound.has_black_hole_fields) {
        continue;
      }
      const std::size_t row = destination->size();
      destination->resize(row + 1);
      destination->particle_index[row] = static_cast<std::uint32_t>(kept_count + (&inbound - commit.inbound_records.data()));
      destination->host_cell_index[row] = inbound.black_hole_fields.host_cell_index;
      destination->subgrid_mass_code[row] = inbound.black_hole_fields.subgrid_mass_code;
      destination->accretion_rate_code[row] = inbound.black_hole_fields.accretion_rate_code;
      destination->feedback_energy_code[row] = inbound.black_hole_fields.feedback_energy_code;
      destination->eddington_ratio[row] = inbound.black_hole_fields.eddington_ratio;
      destination->cumulative_accreted_mass_code[row] = inbound.black_hole_fields.cumulative_accreted_mass_code;
      destination->cumulative_feedback_energy_code[row] =
          inbound.black_hole_fields.cumulative_feedback_energy_code;
      destination->duty_cycle_active_time_code[row] = inbound.black_hole_fields.duty_cycle_active_time_code;
      destination->duty_cycle_total_time_code[row] = inbound.black_hole_fields.duty_cycle_total_time_code;
    }
  };

  auto rebuildTracerSidecar = [&](TracerParticleSidecar* destination) {
    std::vector<std::uint32_t> kept_rows;
    for (std::size_t row = 0; row < tracers.size(); ++row) {
      const auto old_particle = tracers.particle_index[row];
      if (old_particle >= remove_mask.size()) {
        throw std::runtime_error("commitParticleMigration: stale tracer sidecar index");
      }
      if (remove_mask[old_particle] == 0U) {
        kept_rows.push_back(static_cast<std::uint32_t>(row));
      }
    }
    destination->resize(kept_rows.size());
    for (std::size_t row = 0; row < kept_rows.size(); ++row) {
      const std::size_t source = kept_rows[row];
      const auto old_particle = tracers.particle_index[source];
      destination->particle_index[row] = old_to_new[old_particle];
      destination->parent_particle_id[row] = tracers.parent_particle_id[source];
      destination->injection_step[row] = tracers.injection_step[source];
      destination->host_cell_index[row] = tracers.host_cell_index[source];
      destination->mass_fraction_of_host[row] = tracers.mass_fraction_of_host[source];
      destination->last_host_mass_code[row] = tracers.last_host_mass_code[source];
      destination->cumulative_exchanged_mass_code[row] = tracers.cumulative_exchanged_mass_code[source];
    }
    for (const auto& inbound : commit.inbound_records) {
      if (!inbound.has_tracer_fields) {
        continue;
      }
      const std::size_t row = destination->size();
      destination->resize(row + 1);
      destination->particle_index[row] = static_cast<std::uint32_t>(kept_count + (&inbound - commit.inbound_records.data()));
      destination->parent_particle_id[row] = inbound.tracer_fields.parent_particle_id;
      destination->injection_step[row] = inbound.tracer_fields.injection_step;
      destination->host_cell_index[row] = inbound.tracer_fields.host_cell_index;
      destination->mass_fraction_of_host[row] = inbound.tracer_fields.mass_fraction_of_host;
      destination->last_host_mass_code[row] = inbound.tracer_fields.last_host_mass_code;
      destination->cumulative_exchanged_mass_code[row] = inbound.tracer_fields.cumulative_exchanged_mass_code;
    }
  };

  StarParticleSidecar rebuilt_star;
  BlackHoleParticleSidecar rebuilt_black_holes;
  TracerParticleSidecar rebuilt_tracers;
  rebuildStarSidecar(&rebuilt_star);
  rebuildBlackHoleSidecar(&rebuilt_black_holes);
  rebuildTracerSidecar(&rebuilt_tracers);

  particles = std::move(new_particles);
  particle_sidecar = std::move(new_sidecar);
  star_particles = std::move(rebuilt_star);
  black_holes = std::move(rebuilt_black_holes);
  tracers = std::move(rebuilt_tracers);

  species.count_by_species.fill(0);
  for (const auto tag : particle_sidecar.species_tag) {
    if (!isValidSpeciesTag(tag)) {
      throw std::invalid_argument("commitParticleMigration: invalid species tag in rebuilt sidecar");
    }
    ++species.count_by_species[tag];
  }
  rebuildSpeciesIndex();
  bumpParticleIndexGeneration();
}

}  // namespace cosmosim::core
