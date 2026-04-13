#include <cassert>
#include <cmath>
#include <filesystem>
#include <stdexcept>
#include <vector>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/core/provenance.hpp"
#include "cosmosim/core/time_integration.hpp"
#include "cosmosim/io/restart_checkpoint.hpp"

#if COSMOSIM_ENABLE_HDF5
#include <hdf5.h>
#endif

namespace {

void populateState(cosmosim::core::SimulationState& state) {
  state.resizeParticles(5);
  state.resizeCells(3);
  state.resizePatches(1);

  state.species.count_by_species = {2, 1, 1, 1, 0};
  for (std::size_t i = 0; i < state.particles.size(); ++i) {
    state.particle_sidecar.particle_id[i] = 100 + i;
    state.particle_sidecar.sfc_key[i] = 200 + i;
    state.particle_sidecar.owning_rank[i] = 0;
    state.particle_sidecar.particle_flags[i] = static_cast<std::uint32_t>(i);
    state.particles.position_x_comoving[i] = 0.1 * static_cast<double>(i + 1);
    state.particles.position_y_comoving[i] = 0.2 * static_cast<double>(i + 1);
    state.particles.position_z_comoving[i] = 0.3 * static_cast<double>(i + 1);
    state.particles.velocity_x_peculiar[i] = 1.0 + static_cast<double>(i);
    state.particles.velocity_y_peculiar[i] = 2.0 + static_cast<double>(i);
    state.particles.velocity_z_peculiar[i] = 3.0 + static_cast<double>(i);
    state.particles.mass_code[i] = 10.0 + static_cast<double>(i);
    state.particles.time_bin[i] = static_cast<std::uint8_t>(i % 2);
  }
  state.particle_sidecar.species_tag = {
      static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kDarkMatter),
      static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kDarkMatter),
      static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kGas),
      static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kStar),
      static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kBlackHole)};

  state.star_particles.resize(1);
  state.star_particles.particle_index[0] = 3;
  state.star_particles.formation_scale_factor[0] = 0.4;
  state.star_particles.birth_mass_code[0] = 9.0;
  state.star_particles.metallicity_mass_fraction[0] = 0.02;

  state.black_holes.resize(1);
  state.black_holes.particle_index[0] = 4;
  state.black_holes.host_cell_index[0] = 3;
  state.black_holes.subgrid_mass_code[0] = 12.0;
  state.black_holes.accretion_rate_code[0] = 0.1;
  state.black_holes.feedback_energy_code[0] = 0.8;
  state.black_holes.eddington_ratio[0] = 0.2;
  state.black_holes.cumulative_accreted_mass_code[0] = 0.3;
  state.black_holes.cumulative_feedback_energy_code[0] = 1.2;
  state.black_holes.duty_cycle_active_time_code[0] = 4.0;
  state.black_holes.duty_cycle_total_time_code[0] = 5.0;

  state.tracers.resize(0);

  for (std::size_t i = 0; i < state.cells.size(); ++i) {
    state.cells.center_x_comoving[i] = 0.5 * static_cast<double>(i + 1);
    state.cells.center_y_comoving[i] = 0.6 * static_cast<double>(i + 1);
    state.cells.center_z_comoving[i] = 0.7 * static_cast<double>(i + 1);
    state.cells.mass_code[i] = 20.0 + static_cast<double>(i);
    state.cells.time_bin[i] = 0;
    state.cells.patch_index[i] = 0;

    state.gas_cells.density_code[i] = 30.0 + static_cast<double>(i);
    state.gas_cells.pressure_code[i] = 40.0 + static_cast<double>(i);
    state.gas_cells.internal_energy_code[i] = 50.0 + static_cast<double>(i);
    state.gas_cells.temperature_code[i] = 60.0 + static_cast<double>(i);
    state.gas_cells.sound_speed_code[i] = 70.0 + static_cast<double>(i);
    state.gas_cells.recon_gradient_x[i] = 0.0;
    state.gas_cells.recon_gradient_y[i] = 0.0;
    state.gas_cells.recon_gradient_z[i] = 0.0;
  }

  state.patches.patch_id[0] = 11;
  state.patches.level[0] = 0;
  state.patches.first_cell[0] = 0;
  state.patches.cell_count[0] = 3;

  state.metadata.run_name = "restart_integration";
  state.metadata.step_index = 77;
  state.metadata.scale_factor = 0.25;
  state.metadata.normalized_config_hash_hex = "beadfeedbeadfeed";

  cosmosim::core::ModuleSidecarBlock module_sidecar;
  module_sidecar.module_name = "hydro";
  module_sidecar.schema_version = 3;
  module_sidecar.payload = {std::byte{0x01}, std::byte{0x02}, std::byte{0x03}};
  state.sidecars.upsert(std::move(module_sidecar));

  state.rebuildSpeciesIndex();
  assert(state.validateOwnershipInvariants());
}

void testRestartRoundtrip() {
#if COSMOSIM_ENABLE_HDF5
  cosmosim::core::SimulationState state;
  populateState(state);

  cosmosim::core::IntegratorState integrator_state;
  integrator_state.current_time_code = 0.123;
  integrator_state.current_scale_factor = 0.25;
  integrator_state.dt_time_code = 0.001;
  integrator_state.step_index = 77;
  integrator_state.time_bins.hierarchical_enabled = true;
  integrator_state.time_bins.active_bin = 1;
  integrator_state.time_bins.max_bin = 3;

  cosmosim::core::HierarchicalTimeBinScheduler scheduler(3);
  scheduler.reset(5, 2, 8);
  scheduler.setElementBin(1, 1, scheduler.currentTick());
  scheduler.requestBinTransition(4, 3);
  scheduler.beginSubstep();
  scheduler.endSubstep();

  cosmosim::io::RestartWritePayload payload;
  payload.state = &state;
  payload.integrator_state = &integrator_state;
  payload.scheduler = &scheduler;
  payload.normalized_config_text = "schema_version = 1\nmode = zoom_in\n";
  payload.normalized_config_hash_hex = cosmosim::core::stableConfigHashHex(payload.normalized_config_text);
  payload.provenance = cosmosim::core::makeProvenanceRecord(payload.normalized_config_hash_hex, "deadbeef");

  const std::filesystem::path checkpoint_path =
      std::filesystem::temp_directory_path() / "cosmosim_restart_roundtrip.hdf5";

  cosmosim::io::writeRestartCheckpointHdf5(checkpoint_path, payload);
  const cosmosim::io::RestartReadResult restored = cosmosim::io::readRestartCheckpointHdf5(checkpoint_path);

  assert(restored.state.validateOwnershipInvariants());
  assert(restored.state.metadata.run_name == state.metadata.run_name);
  assert(restored.state.metadata.schema_version == state.metadata.schema_version);
  assert(restored.state.metadata.step_index == state.metadata.step_index);
  assert(std::abs(restored.state.metadata.scale_factor - state.metadata.scale_factor) < 1.0e-15);
  assert(restored.state.metadata.normalized_config_hash_hex == state.metadata.normalized_config_hash_hex);
  assert(restored.state.particle_sidecar.particle_id == state.particle_sidecar.particle_id);
  assert(restored.state.particle_sidecar.sfc_key == state.particle_sidecar.sfc_key);
  assert(restored.state.particle_sidecar.species_tag == state.particle_sidecar.species_tag);
  assert(restored.state.particle_sidecar.particle_flags == state.particle_sidecar.particle_flags);
  assert(restored.state.particle_sidecar.owning_rank == state.particle_sidecar.owning_rank);
  assert(restored.state.star_particles.particle_index == state.star_particles.particle_index);
  assert(restored.state.star_particles.formation_scale_factor == state.star_particles.formation_scale_factor);
  assert(restored.state.black_holes.particle_index == state.black_holes.particle_index);
  assert(restored.state.black_holes.host_cell_index == state.black_holes.host_cell_index);
  assert(restored.state.black_holes.subgrid_mass_code == state.black_holes.subgrid_mass_code);
  assert(restored.state.black_holes.accretion_rate_code == state.black_holes.accretion_rate_code);
  assert(restored.state.black_holes.feedback_energy_code == state.black_holes.feedback_energy_code);
  assert(restored.state.black_holes.eddington_ratio == state.black_holes.eddington_ratio);
  assert(
      restored.state.black_holes.cumulative_accreted_mass_code ==
      state.black_holes.cumulative_accreted_mass_code);
  assert(
      restored.state.black_holes.cumulative_feedback_energy_code ==
      state.black_holes.cumulative_feedback_energy_code);
  assert(
      restored.state.black_holes.duty_cycle_active_time_code ==
      state.black_holes.duty_cycle_active_time_code);
  assert(
      restored.state.black_holes.duty_cycle_total_time_code ==
      state.black_holes.duty_cycle_total_time_code);
  assert(restored.state.species.count_by_species == state.species.count_by_species);
  assert(restored.integrator_state.step_index == integrator_state.step_index);
  assert(std::abs(restored.integrator_state.current_time_code - integrator_state.current_time_code) < 1.0e-15);
  assert(restored.integrator_state.scheme == integrator_state.scheme);
  assert(restored.integrator_state.time_bins.hierarchical_enabled == integrator_state.time_bins.hierarchical_enabled);
  assert(restored.integrator_state.time_bins.active_bin == integrator_state.time_bins.active_bin);
  assert(restored.integrator_state.time_bins.max_bin == integrator_state.time_bins.max_bin);
  assert(restored.scheduler_state.max_bin == scheduler.maxBin());
  assert(restored.scheduler_state.current_tick == scheduler.currentTick());
  assert(restored.scheduler_state.bin_index.size() == scheduler.elementCount());
  const auto original_scheduler_state = scheduler.exportPersistentState();
  assert(restored.scheduler_state.bin_index == original_scheduler_state.bin_index);
  assert(restored.scheduler_state.next_activation_tick == original_scheduler_state.next_activation_tick);
  assert(restored.scheduler_state.active_flag == original_scheduler_state.active_flag);
  assert(restored.scheduler_state.pending_bin_index == original_scheduler_state.pending_bin_index);
  assert(restored.normalized_config_hash_hex == payload.normalized_config_hash_hex);
  assert(restored.normalized_config_text == payload.normalized_config_text);
  assert(restored.provenance.config_hash_hex == payload.provenance.config_hash_hex);
  assert(restored.state.sidecars.find("hydro") != nullptr);
  const cosmosim::core::ModuleSidecarBlock* hydro_sidecar = restored.state.sidecars.find("hydro");
  assert(hydro_sidecar->schema_version == 3);
  assert(hydro_sidecar->payload.size() == 3);
  assert(hydro_sidecar->payload[0] == std::byte{0x01});
  assert(hydro_sidecar->payload[1] == std::byte{0x02});
  assert(hydro_sidecar->payload[2] == std::byte{0x03});
  assert(restored.payload_hash == cosmosim::io::restartPayloadIntegrityHash(payload));
  assert(restored.payload_hash_hex == cosmosim::io::restartPayloadIntegrityHashHex(payload));

  cosmosim::core::HierarchicalTimeBinScheduler resumed_scheduler(restored.scheduler_state.max_bin);
  resumed_scheduler.importPersistentState(restored.scheduler_state);
  assert(resumed_scheduler.currentTick() == scheduler.currentTick());

  hid_t tamper_file = H5Fopen(checkpoint_path.string().c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  assert(tamper_file >= 0);
  hid_t tamper_attr = H5Aopen(tamper_file, "payload_integrity_hash", H5P_DEFAULT);
  assert(tamper_attr >= 0);
  std::uint64_t bad_hash = 0;
  assert(H5Awrite(tamper_attr, H5T_NATIVE_UINT64, &bad_hash) >= 0);
  H5Aclose(tamper_attr);
  H5Fclose(tamper_file);

  bool integrity_threw = false;
  try {
    (void)cosmosim::io::readRestartCheckpointHdf5(checkpoint_path);
  } catch (const std::runtime_error&) {
    integrity_threw = true;
  }
  assert(integrity_threw);

  std::filesystem::remove(checkpoint_path);
#else
  bool threw = false;
  try {
    cosmosim::io::RestartWritePayload payload;
    cosmosim::io::writeRestartCheckpointHdf5("unused.hdf5", payload);
  } catch (const std::runtime_error&) {
    threw = true;
  }
  assert(threw);
#endif
}

}  // namespace

int main() {
  testRestartRoundtrip();
  return 0;
}
