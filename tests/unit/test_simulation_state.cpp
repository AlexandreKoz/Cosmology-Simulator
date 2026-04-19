#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include "cosmosim/core/simulation_state.hpp"

int main() {
  static_assert(
      sizeof(cosmosim::core::GravityParticleKernelView) ==
          sizeof(std::span<std::uint32_t>) + (7 * sizeof(std::span<double>)),
      "GravityParticleKernelView hot contract changed: unexpected extra field(s)");
  static_assert(
      sizeof(cosmosim::core::HydroCellKernelView) ==
          sizeof(std::span<std::uint32_t>) + (6 * sizeof(std::span<double>)),
      "HydroCellKernelView hot contract changed: unexpected extra field(s)");

  cosmosim::core::SimulationState state;
  state.resizeParticles(5);
  state.resizeCells(3);
  state.resizePatches(1);

  for (std::size_t i = 0; i < 5; ++i) {
    state.particle_sidecar.particle_id[i] = 1000 + i;
    state.particle_sidecar.species_tag[i] =
        static_cast<std::uint32_t>((i < 2) ? cosmosim::core::ParticleSpecies::kDarkMatter
                                           : cosmosim::core::ParticleSpecies::kStar);
    state.particles.position_x_comoving[i] = static_cast<double>(i);
    state.particles.position_y_comoving[i] = static_cast<double>(i + 1);
    state.particles.position_z_comoving[i] = static_cast<double>(i + 2);
    state.particles.velocity_x_peculiar[i] = 0.1;
    state.particles.velocity_y_peculiar[i] = 0.2;
    state.particles.velocity_z_peculiar[i] = 0.3;
    state.particles.mass_code[i] = 1.0;
    state.particles.time_bin[i] = static_cast<std::uint8_t>(i % 3);
  }
  state.species.count_by_species = {2, 0, 3, 0, 0};

  state.star_particles.resize(3);
  for (std::size_t i = 0; i < 3; ++i) {
    state.star_particles.particle_index[i] = static_cast<std::uint32_t>(i + 2);
    state.star_particles.formation_scale_factor[i] = 0.5;
    state.star_particles.birth_mass_code[i] = 1.1;
    state.star_particles.metallicity_mass_fraction[i] = 0.02;
  }

  state.patches.patch_id[0] = 42;
  state.patches.level[0] = 0;
  state.patches.first_cell[0] = 0;
  state.patches.cell_count[0] = 3;

  for (std::size_t i = 0; i < 3; ++i) {
    state.cells.center_x_comoving[i] = static_cast<double>(i);
    state.cells.center_y_comoving[i] = static_cast<double>(i) + 0.5;
    state.cells.center_z_comoving[i] = static_cast<double>(i) + 1.0;
    state.cells.mass_code[i] = 3.0;
    state.cells.time_bin[i] = 1;
    state.cells.patch_index[i] = 0;

    state.gas_cells.density_code[i] = 10.0 + static_cast<double>(i);
    state.gas_cells.pressure_code[i] = 1.0 + static_cast<double>(i);
    state.gas_cells.internal_energy_code[i] = 0.1 + static_cast<double>(i);
    state.gas_cells.temperature_code[i] = 100.0 + static_cast<double>(i);
    state.gas_cells.sound_speed_code[i] = 2.0;
    state.gas_cells.recon_gradient_x[i] = 0.01;
    state.gas_cells.recon_gradient_y[i] = 0.02;
    state.gas_cells.recon_gradient_z[i] = 0.03;
  }

  state.rebuildSpeciesIndex();
  assert(state.particle_species_index.count(cosmosim::core::ParticleSpecies::kDarkMatter) == 2);
  assert(state.particle_species_index.count(cosmosim::core::ParticleSpecies::kStar) == 3);
  assert(state.particle_species_index.localIndex(4) == 2);
  assert(state.particle_species_index.globalIndex(cosmosim::core::ParticleSpecies::kStar, 1) == 3);

  assert(state.validateOwnershipInvariants());
  assert(state.validateUniqueParticleIds());

  state.metadata.run_name = "unit_state";
  state.metadata.normalized_config_hash = 1234;
  state.metadata.normalized_config_hash_hex = "0x4d2";
  state.metadata.step_index = 17;
  state.metadata.scale_factor = 0.5;

  const std::string serialized = state.metadata.serialize();
  const auto parsed = cosmosim::core::StateMetadata::deserialize(serialized);
  assert(parsed.run_name == state.metadata.run_name);
  assert(parsed.normalized_config_hash == state.metadata.normalized_config_hash);
  assert(parsed.step_index == state.metadata.step_index);
  assert(parsed.scale_factor == state.metadata.scale_factor);

  const auto packet = state.packSpeciesTransferPacket(cosmosim::core::ParticleSpecies::kStar);
  assert(packet.particle_id.size() == 3);
  assert(packet.particle_id[0] == 1002);
  for (std::size_t i = 0; i < packet.particle_id.size(); ++i) {
    const auto global = state.particle_species_index.globalIndex(
        cosmosim::core::ParticleSpecies::kStar,
        static_cast<std::uint32_t>(i));
    assert(packet.particle_id[i] == state.particle_sidecar.particle_id[global]);
    assert(packet.owning_rank[i] == state.particle_sidecar.owning_rank[global]);
    assert(packet.mass_code[i] == state.particles.mass_code[global]);
  }

  cosmosim::core::TransientStepWorkspace workspace;
  const std::array<std::uint32_t, 2> particle_indices{1, 3};
  auto particle_view = cosmosim::core::buildParticleActiveView(state, particle_indices, workspace);
  assert(particle_view.size() == 2);
  assert(particle_view.particle_id[0] == 1001);
  assert(particle_view.position_x_comoving[1] == 3.0);

  const std::array<std::uint32_t, 1> cell_indices{2};
  auto cell_view = cosmosim::core::buildCellActiveView(state, cell_indices, workspace);
  assert(cell_view.size() == 1);
  assert(cell_view.density_code[0] == 12.0);

  const std::array<std::uint32_t, 2> kernel_particle_indices{0, 2};
  const auto particle_id_before = state.particle_sidecar.particle_id;
  auto gravity_view =
      cosmosim::core::buildGravityParticleKernelView(state, kernel_particle_indices, workspace);
  gravity_view.position_x_comoving[0] = 9.0;
  gravity_view.mass_code[1] = 8.0;
  cosmosim::core::scatterGravityParticleKernelView(gravity_view, state);
  assert(state.particles.position_x_comoving[0] == 9.0);
  assert(state.particles.mass_code[2] == 8.0);
  assert(state.particle_sidecar.particle_id == particle_id_before);

  const std::array<std::uint32_t, 2> kernel_cell_indices{0, 1};
  auto hydro_view = cosmosim::core::buildHydroCellKernelView(state, kernel_cell_indices, workspace);
  hydro_view.center_x_comoving[1] = 44.0;
  hydro_view.density_code[0] = 55.0;
  hydro_view.pressure_code[1] = 66.0;
  cosmosim::core::scatterHydroCellKernelView(hydro_view, state);
  assert(state.cells.center_x_comoving[1] == 44.0);
  assert(state.gas_cells.density_code[0] == 55.0);
  assert(state.gas_cells.pressure_code[1] == 66.0);

  bool threw = false;
  try {
    const std::array<std::uint32_t, 1> bad_indices{8};
    (void)cosmosim::core::buildParticleActiveView(state, bad_indices, workspace);
  } catch (const std::out_of_range&) {
    threw = true;
  }
  assert(threw);

  state.particle_sidecar.particle_id[1] = state.particle_sidecar.particle_id[0];
  assert(!state.validateUniqueParticleIds());
  state.particle_sidecar.particle_id[1] = 1001;
  assert(state.validateUniqueParticleIds());

  state.species.count_by_species = {3, 0, 2, 0, 0};
  assert(!state.validateOwnershipInvariants());
  state.species.count_by_species = {2, 0, 3, 0, 0};
  assert(state.validateOwnershipInvariants());

  state.star_particles.resize(0);
  assert(!state.validateOwnershipInvariants());
  state.star_particles.resize(3);
  for (std::size_t i = 0; i < 3; ++i) {
    state.star_particles.particle_index[i] = static_cast<std::uint32_t>(i + 2);
    state.star_particles.formation_scale_factor[i] = 0.5;
    state.star_particles.birth_mass_code[i] = 1.1;
    state.star_particles.metallicity_mass_fraction[i] = 0.02;
  }
  assert(state.validateOwnershipInvariants());

  state.star_particles.particle_index[1] = 2;
  assert(!state.validateOwnershipInvariants());

  state.star_particles.particle_index[1] = 3;
  assert(state.validateOwnershipInvariants());

  const std::vector<std::uint32_t> outbound_indices = {3};
  auto outbound_records = state.packParticleMigrationRecords(outbound_indices);
  assert(outbound_records.size() == 1);
  assert(outbound_records[0].particle_id == 1003);
  outbound_records[0].owning_rank = 0;

  cosmosim::core::ParticleMigrationRecord inbound_record;
  inbound_record.particle_id = 9090;
  inbound_record.sfc_key = 777;
  inbound_record.species_tag = static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kStar);
  inbound_record.particle_flags = 42;
  inbound_record.owning_rank = 0;
  inbound_record.position_x_comoving = 11.0;
  inbound_record.position_y_comoving = 12.0;
  inbound_record.position_z_comoving = 13.0;
  inbound_record.velocity_x_peculiar = 1.0;
  inbound_record.velocity_y_peculiar = 2.0;
  inbound_record.velocity_z_peculiar = 3.0;
  inbound_record.mass_code = 9.5;
  inbound_record.time_bin = 2;
  inbound_record.has_star_fields = true;
  inbound_record.star_fields.formation_scale_factor = 0.75;
  inbound_record.star_fields.birth_mass_code = 5.0;
  inbound_record.star_fields.metallicity_mass_fraction = 0.04;
  inbound_record.star_fields.stellar_age_years_last = 1.0e7;

  cosmosim::core::ParticleMigrationCommit commit;
  commit.world_rank = 0;
  commit.outbound_local_indices = outbound_indices;
  commit.inbound_records = {inbound_record};
  commit.stale_local_ghost_indices = {1};
  state.commitParticleMigration(commit);

  assert(state.validateOwnershipInvariants());
  assert(state.particles.size() == 4);
  assert(state.particle_sidecar.owning_rank[3] == 0);
  assert(state.particle_sidecar.particle_id[3] == 9090);
  assert(state.star_particles.size() == 3);
  assert(state.species.count_by_species[static_cast<std::size_t>(cosmosim::core::ParticleSpecies::kStar)] == 3);

  bool found_migrated_star = false;
  for (std::size_t row = 0; row < state.star_particles.size(); ++row) {
    const auto particle_index = state.star_particles.particle_index[row];
    if (state.particle_sidecar.particle_id[particle_index] == 9090) {
      found_migrated_star = true;
      assert(state.star_particles.formation_scale_factor[row] == 0.75);
      assert(state.star_particles.birth_mass_code[row] == 5.0);
    }
  }
  assert(found_migrated_star);

  return 0;
}
