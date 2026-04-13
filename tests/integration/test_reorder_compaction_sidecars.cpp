#include <array>
#include <cassert>
#include <cstdint>

#include "cosmosim/core/simulation_state.hpp"

namespace {

void compact_front(cosmosim::core::SimulationState& state, std::size_t keep_count) {
  state.particles.position_x_comoving.resize(keep_count);
  state.particles.position_y_comoving.resize(keep_count);
  state.particles.position_z_comoving.resize(keep_count);
  state.particles.velocity_x_peculiar.resize(keep_count);
  state.particles.velocity_y_peculiar.resize(keep_count);
  state.particles.velocity_z_peculiar.resize(keep_count);
  state.particles.mass_code.resize(keep_count);
  state.particles.time_bin.resize(keep_count);

  state.particle_sidecar.particle_id.resize(keep_count);
  state.particle_sidecar.sfc_key.resize(keep_count);
  state.particle_sidecar.species_tag.resize(keep_count);
  state.particle_sidecar.particle_flags.resize(keep_count);
  state.particle_sidecar.owning_rank.resize(keep_count);

  auto compact_index_sidecar = [keep_count](cosmosim::core::AlignedVector<std::uint32_t>& indices) {
    std::size_t write_pos = 0;
    for (const auto index : indices) {
      if (index < keep_count) {
        indices[write_pos] = index;
        ++write_pos;
      }
    }
    indices.resize(write_pos);
  };

  compact_index_sidecar(state.star_particles.particle_index);
  state.star_particles.formation_scale_factor.resize(state.star_particles.particle_index.size());
  state.star_particles.birth_mass_code.resize(state.star_particles.particle_index.size());
  state.star_particles.metallicity_mass_fraction.resize(state.star_particles.particle_index.size());

  compact_index_sidecar(state.black_holes.particle_index);
  state.black_holes.host_cell_index.resize(state.black_holes.particle_index.size());
  state.black_holes.subgrid_mass_code.resize(state.black_holes.particle_index.size());
  state.black_holes.accretion_rate_code.resize(state.black_holes.particle_index.size());
  state.black_holes.feedback_energy_code.resize(state.black_holes.particle_index.size());
  state.black_holes.eddington_ratio.resize(state.black_holes.particle_index.size());
  state.black_holes.cumulative_accreted_mass_code.resize(state.black_holes.particle_index.size());
  state.black_holes.cumulative_feedback_energy_code.resize(state.black_holes.particle_index.size());
  state.black_holes.duty_cycle_active_time_code.resize(state.black_holes.particle_index.size());
  state.black_holes.duty_cycle_total_time_code.resize(state.black_holes.particle_index.size());

  compact_index_sidecar(state.tracers.particle_index);
  state.tracers.parent_particle_id.resize(state.tracers.particle_index.size());
  state.tracers.injection_step.resize(state.tracers.particle_index.size());

  std::array<std::uint64_t, 5> counts{};
  for (const auto tag : state.particle_sidecar.species_tag) {
    ++counts[tag];
  }
  state.species.count_by_species = counts;
  state.rebuildSpeciesIndex();
}

}  // namespace

int main() {
  cosmosim::core::SimulationState state;
  state.resizeParticles(8);
  const std::array<std::uint32_t, 8> species{2, 0, 4, 2, 1, 3, 0, 2};
  const std::array<std::uint8_t, 8> time_bins{1, 0, 2, 1, 0, 2, 1, 0};
  const std::array<std::uint64_t, 8> sfc{55, 11, 77, 33, 22, 88, 44, 66};

  for (std::size_t i = 0; i < 8; ++i) {
    state.particle_sidecar.particle_id[i] = 500 + i;
    state.particle_sidecar.species_tag[i] = species[i];
    state.particle_sidecar.sfc_key[i] = sfc[i];
    state.particles.time_bin[i] = time_bins[i];
    state.particles.position_x_comoving[i] = static_cast<double>(i);
    state.particles.position_y_comoving[i] = static_cast<double>(i);
    state.particles.position_z_comoving[i] = static_cast<double>(i);
    state.particles.velocity_x_peculiar[i] = 0.0;
    state.particles.velocity_y_peculiar[i] = 0.0;
    state.particles.velocity_z_peculiar[i] = 0.0;
    state.particles.mass_code[i] = 1.0;
  }

  state.star_particles.resize(3);
  state.star_particles.particle_index = {0, 3, 7};
  state.black_holes.resize(1);
  state.black_holes.particle_index[0] = 5;
  state.tracers.resize(1);
  state.tracers.particle_index[0] = 2;

  std::array<std::uint64_t, 5> counts{};
  for (const auto tag : species) {
    ++counts[tag];
  }
  state.species.count_by_species = counts;
  state.rebuildSpeciesIndex();

  const auto reorder_map = cosmosim::core::buildParticleReorderMap(
      state,
      cosmosim::core::ParticleReorderMode::kBySpecies);
  cosmosim::core::reorderParticles(state, reorder_map);
  cosmosim::core::debugAssertNoStaleParticleIndices(state);
  assert(state.validateOwnershipInvariants());

  compact_front(state, 6);
  cosmosim::core::debugAssertNoStaleParticleIndices(state);
  assert(state.validateOwnershipInvariants());

  for (const auto index : state.star_particles.particle_index) {
    assert(state.particle_sidecar.species_tag[index] == static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kStar));
  }

  return 0;
}
