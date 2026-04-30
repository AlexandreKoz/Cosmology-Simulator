#include <cassert>
#include <cmath>
#include <filesystem>
#include <span>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/core/provenance.hpp"
#include "cosmosim/core/time_integration.hpp"
#include "cosmosim/gravity/tree_softening.hpp"
#include "cosmosim/io/restart_checkpoint.hpp"

#if COSMOSIM_ENABLE_HDF5
#include <hdf5.h>
#endif

namespace {

void populateState(cosmosim::core::SimulationState& state) {
  state.resizeParticles(6);
  state.resizeCells(3);
  state.resizePatches(1);

  state.species.count_by_species = {1, 3, 1, 1, 0};
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
  state.particle_sidecar.gravity_softening_comoving = {0.001, 0.002, 0.003, 0.004, 0.005, 0.006};
  state.particle_sidecar.has_gravity_softening_override = {1U, 1U, 1U, 1U, 1U, 1U};
  state.particle_sidecar.species_tag = {
      static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kDarkMatter),
      static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kGas),
      static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kGas),
      static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kGas),
      static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kStar),
      static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kBlackHole)};

  state.star_particles.resize(1);
  state.star_particles.particle_index[0] = 4;
  state.star_particles.formation_scale_factor[0] = 0.4;
  state.star_particles.birth_mass_code[0] = 9.0;
  state.star_particles.metallicity_mass_fraction[0] = 0.02;
  state.star_particles.stellar_age_years_last[0] = 1.25e7;
  state.star_particles.stellar_returned_mass_cumulative_code[0] = 0.4;
  state.star_particles.stellar_returned_metals_cumulative_code[0] = 0.03;
  state.star_particles.stellar_feedback_energy_cumulative_erg[0] = 5.0e49;
  for (std::size_t channel = 0; channel < state.star_particles.stellar_returned_mass_channel_cumulative_code.size(); ++channel) {
    state.star_particles.stellar_returned_mass_channel_cumulative_code[channel][0] =
        0.01 * static_cast<double>(channel + 1);
    state.star_particles.stellar_returned_metals_channel_cumulative_code[channel][0] =
        0.001 * static_cast<double>(channel + 1);
    state.star_particles.stellar_feedback_energy_channel_cumulative_erg[channel][0] =
        1.0e48 * static_cast<double>(channel + 1);
  }

  state.black_holes.resize(1);
  state.black_holes.particle_index[0] = 5;
  state.black_holes.host_cell_index[0] = 2;
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

std::unordered_map<std::uint64_t, double> gasDensityByParticleId(const cosmosim::core::SimulationState& state) {
  std::unordered_map<std::uint64_t, double> by_id;
  const auto gas_globals = state.particle_species_index.globalIndices(cosmosim::core::ParticleSpecies::kGas);
  assert(gas_globals.size() == state.cells.size());
  for (std::size_t cell = 0; cell < state.cells.size(); ++cell) {
    by_id.emplace(state.particle_sidecar.particle_id[gas_globals[cell]], state.gas_cells.density_code[cell]);
  }
  return by_id;
}

std::vector<std::uint64_t> particleIdsForIndices(
    const cosmosim::core::SimulationState& state,
    std::span<const std::uint32_t> indices) {
  std::vector<std::uint64_t> ids;
  ids.reserve(indices.size());
  for (const std::uint32_t index : indices) {
    ids.push_back(state.particle_sidecar.particle_id[index]);
  }
  return ids;
}

std::vector<std::uint64_t> schedulerActiveIdsFromPersistentState(
    const cosmosim::core::SimulationState& state,
    const cosmosim::core::TimeBinPersistentState& persistent_state) {
  cosmosim::core::HierarchicalTimeBinScheduler scheduler_probe(persistent_state.max_bin);
  scheduler_probe.importPersistentState(persistent_state);
  const auto active = scheduler_probe.beginSubstep();
  return particleIdsForIndices(state, active);
}

void assertSofteningPriorityResolution(const cosmosim::core::SimulationState& state) {
  cosmosim::gravity::TreeSofteningPolicy global_softening;
  global_softening.epsilon_comoving = 0.123;

  cosmosim::gravity::TreeSofteningSpeciesPolicy species_softening;
  species_softening.enabled = true;
  species_softening.epsilon_comoving_by_species = {0.011, 0.022, 0.033, 0.044, 0.055};

  const cosmosim::gravity::TreeSofteningView softening_view{
      .source_species_tag = std::span<const std::uint32_t>(
          state.particle_sidecar.species_tag.data(),
          state.particle_sidecar.species_tag.size()),
      .source_particle_epsilon_comoving = std::span<const double>(
          state.particle_sidecar.gravity_softening_comoving.data(),
          state.particle_sidecar.gravity_softening_comoving.size()),
      .source_particle_epsilon_override_mask = std::span<const std::uint8_t>(
          state.particle_sidecar.has_gravity_softening_override.data(),
          state.particle_sidecar.has_gravity_softening_override.size()),
      .target_particle_epsilon_comoving = {},
      .species_policy = species_softening,
  };

  assert(std::abs(cosmosim::gravity::resolveSourceSofteningEpsilon(0, global_softening, softening_view) - 0.001) < 1.0e-15);
  assert(std::abs(cosmosim::gravity::resolveSourceSofteningEpsilon(2, global_softening, softening_view) - 0.003) < 1.0e-15);
}

void testRestartRoundtrip() {
#if COSMOSIM_ENABLE_HDF5
  cosmosim::core::SimulationState state;
  populateState(state);
  const auto gas_density_before = gasDensityByParticleId(state);

  cosmosim::core::IntegratorState integrator_state;
  integrator_state.current_time_code = 0.123;
  integrator_state.current_scale_factor = 0.25;
  integrator_state.dt_time_code = 0.001;
  integrator_state.step_index = 77;
  integrator_state.time_bins.hierarchical_enabled = true;
  integrator_state.time_bins.active_bin = 1;
  integrator_state.time_bins.max_bin = 3;

  cosmosim::core::HierarchicalTimeBinScheduler scheduler(3);
  scheduler.reset(6, 2, 8);
  scheduler.setElementBin(1, 1, scheduler.currentTick());
  scheduler.requestBinTransition(4, 3);
  scheduler.beginSubstep();
  scheduler.endSubstep();
  const std::vector<std::uint64_t> expected_active_particle_ids =
      schedulerActiveIdsFromPersistentState(state, scheduler.exportPersistentState());

  cosmosim::io::RestartWritePayload payload;
  payload.state = &state;
  payload.integrator_state = &integrator_state;
  payload.scheduler = &scheduler;
  payload.normalized_config_text = "schema_version = 1\nmode = zoom_in\n";
  payload.normalized_config_hash_hex = cosmosim::core::stableConfigHashHex(payload.normalized_config_text);
  payload.provenance = cosmosim::core::makeProvenanceRecord(payload.normalized_config_hash_hex, "deadbeef");
  payload.provenance.gravity_treepm_pm_grid = 32;
  payload.provenance.gravity_treepm_pm_grid_nx = 32;
  payload.provenance.gravity_treepm_pm_grid_ny = 24;
  payload.provenance.gravity_treepm_pm_grid_nz = 16;
  payload.provenance.gravity_treepm_assignment_scheme = "tsc";
  payload.provenance.gravity_treepm_window_deconvolution = true;
  payload.provenance.gravity_treepm_asmth_cells = 1.25;
  payload.provenance.gravity_treepm_rcut_cells = 5.5;
  payload.provenance.gravity_treepm_mesh_spacing_mpc_comoving = std::cbrt(0.25 * 0.5 * 0.75);
  payload.provenance.gravity_treepm_mesh_spacing_x_mpc_comoving = 0.25;
  payload.provenance.gravity_treepm_mesh_spacing_y_mpc_comoving = 0.5;
  payload.provenance.gravity_treepm_mesh_spacing_z_mpc_comoving = 0.75;
  payload.provenance.gravity_treepm_split_scale_mpc_comoving = 0.3125;
  payload.provenance.gravity_treepm_cutoff_radius_mpc_comoving = 1.375;
  payload.provenance.gravity_treepm_update_cadence_steps = 3;
  payload.provenance.gravity_softening_policy = "comoving_fixed";
  payload.provenance.gravity_softening_kernel = "plummer";
  payload.provenance.gravity_softening_epsilon_kpc_comoving = 2.0;
  payload.provenance.gravity_pm_fft_backend = "fftw3";
  payload.distributed_gravity_state.schema_version = 2;
  payload.distributed_gravity_state.decomposition_epoch = integrator_state.step_index;
  payload.distributed_gravity_state.world_size = 1;
  payload.distributed_gravity_state.pm_grid_nx = 32;
  payload.distributed_gravity_state.pm_grid_ny = 32;
  payload.distributed_gravity_state.pm_grid_nz = 32;
  payload.distributed_gravity_state.pm_decomposition_mode = "slab";
  payload.distributed_gravity_state.gravity_kick_opportunity = 9;
  payload.distributed_gravity_state.pm_update_cadence_steps = 3;
  payload.distributed_gravity_state.long_range_field_version = 4;
  payload.distributed_gravity_state.last_long_range_refresh_opportunity = 9;
  payload.distributed_gravity_state.long_range_field_built_step_index = integrator_state.step_index;
  payload.distributed_gravity_state.long_range_field_built_scale_factor = integrator_state.current_scale_factor;
  payload.distributed_gravity_state.long_range_restart_policy = "deterministic_rebuild";
  payload.distributed_gravity_state.owning_rank_by_item.assign(state.particles.size(), 0);
  payload.distributed_gravity_state.pm_slab_begin_x_by_rank = {0};
  payload.distributed_gravity_state.pm_slab_end_x_by_rank = {32};

  const std::filesystem::path checkpoint_path =
      std::filesystem::temp_directory_path() / "cosmosim_restart_roundtrip.hdf5";

  cosmosim::io::writeRestartCheckpointHdf5(checkpoint_path, payload);
  const cosmosim::io::RestartReadResult restored = cosmosim::io::readRestartCheckpointHdf5(checkpoint_path);

  assert(restored.state.validateOwnershipInvariants());
  cosmosim::core::debugAssertGasCellIdentityContract(restored.state);
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
  assert(restored.state.particle_sidecar.gravity_softening_comoving == state.particle_sidecar.gravity_softening_comoving);
  assert(restored.state.star_particles.particle_index == state.star_particles.particle_index);
  assert(restored.state.star_particles.formation_scale_factor == state.star_particles.formation_scale_factor);
  assert(restored.state.star_particles.birth_mass_code == state.star_particles.birth_mass_code);
  assert(restored.state.star_particles.metallicity_mass_fraction == state.star_particles.metallicity_mass_fraction);
  assert(restored.state.star_particles.stellar_age_years_last == state.star_particles.stellar_age_years_last);
  assert(
      restored.state.star_particles.stellar_returned_mass_cumulative_code ==
      state.star_particles.stellar_returned_mass_cumulative_code);
  assert(
      restored.state.star_particles.stellar_returned_metals_cumulative_code ==
      state.star_particles.stellar_returned_metals_cumulative_code);
  assert(
      restored.state.star_particles.stellar_feedback_energy_cumulative_erg ==
      state.star_particles.stellar_feedback_energy_cumulative_erg);
  for (std::size_t channel = 0; channel < state.star_particles.stellar_returned_mass_channel_cumulative_code.size(); ++channel) {
    assert(
        restored.state.star_particles.stellar_returned_mass_channel_cumulative_code[channel] ==
        state.star_particles.stellar_returned_mass_channel_cumulative_code[channel]);
    assert(
        restored.state.star_particles.stellar_returned_metals_channel_cumulative_code[channel] ==
        state.star_particles.stellar_returned_metals_channel_cumulative_code[channel]);
    assert(
        restored.state.star_particles.stellar_feedback_energy_channel_cumulative_erg[channel] ==
        state.star_particles.stellar_feedback_energy_channel_cumulative_erg[channel]);
  }
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
  assert(restored.state.particle_species_index.globalIndices(cosmosim::core::ParticleSpecies::kDarkMatter).size() == 1);
  assert(restored.state.particle_species_index.globalIndices(cosmosim::core::ParticleSpecies::kGas).size() == 3);
  assert(restored.state.particle_species_index.globalIndices(cosmosim::core::ParticleSpecies::kStar).size() == 1);
  assert(restored.state.particle_species_index.globalIndices(cosmosim::core::ParticleSpecies::kBlackHole).size() == 1);
  assert(gasDensityByParticleId(restored.state) == gas_density_before);
  assert(restored.state.gasCellIdentityMatchesParticleOrder());
  assert(restored.state.gas_cells.gas_cell_id == state.gas_cells.gas_cell_id);
  assert(restored.state.gas_cells.parent_particle_id == state.gas_cells.parent_particle_id);
  assertSofteningPriorityResolution(restored.state);
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
  assert(restored.provenance.gravity_treepm_pm_grid == payload.provenance.gravity_treepm_pm_grid);
  assert(
      restored.provenance.gravity_treepm_assignment_scheme ==
      payload.provenance.gravity_treepm_assignment_scheme);
  assert(
      restored.provenance.gravity_treepm_window_deconvolution ==
      payload.provenance.gravity_treepm_window_deconvolution);
  assert(restored.provenance.gravity_treepm_asmth_cells == payload.provenance.gravity_treepm_asmth_cells);
  assert(restored.provenance.gravity_treepm_rcut_cells == payload.provenance.gravity_treepm_rcut_cells);
  assert(
      restored.provenance.gravity_treepm_mesh_spacing_mpc_comoving ==
      payload.provenance.gravity_treepm_mesh_spacing_mpc_comoving);
  assert(
      restored.provenance.gravity_treepm_mesh_spacing_x_mpc_comoving ==
      payload.provenance.gravity_treepm_mesh_spacing_x_mpc_comoving);
  assert(
      restored.provenance.gravity_treepm_mesh_spacing_y_mpc_comoving ==
      payload.provenance.gravity_treepm_mesh_spacing_y_mpc_comoving);
  assert(
      restored.provenance.gravity_treepm_mesh_spacing_z_mpc_comoving ==
      payload.provenance.gravity_treepm_mesh_spacing_z_mpc_comoving);
  assert(
      restored.provenance.gravity_treepm_split_scale_mpc_comoving ==
      payload.provenance.gravity_treepm_split_scale_mpc_comoving);
  assert(
      restored.provenance.gravity_treepm_cutoff_radius_mpc_comoving ==
      payload.provenance.gravity_treepm_cutoff_radius_mpc_comoving);
  assert(
      restored.provenance.gravity_treepm_update_cadence_steps ==
      payload.provenance.gravity_treepm_update_cadence_steps);
  assert(restored.provenance.gravity_softening_policy == payload.provenance.gravity_softening_policy);
  assert(restored.provenance.gravity_softening_kernel == payload.provenance.gravity_softening_kernel);
  assert(
      restored.provenance.gravity_softening_epsilon_kpc_comoving ==
      payload.provenance.gravity_softening_epsilon_kpc_comoving);
  assert(restored.provenance.gravity_pm_fft_backend == payload.provenance.gravity_pm_fft_backend);
  assert(restored.distributed_gravity_state.schema_version == payload.distributed_gravity_state.schema_version);
  assert(restored.distributed_gravity_state.decomposition_epoch == payload.distributed_gravity_state.decomposition_epoch);
  assert(restored.distributed_gravity_state.world_size == payload.distributed_gravity_state.world_size);
  assert(restored.distributed_gravity_state.pm_grid_nx == payload.distributed_gravity_state.pm_grid_nx);
  assert(restored.distributed_gravity_state.pm_grid_ny == payload.distributed_gravity_state.pm_grid_ny);
  assert(restored.distributed_gravity_state.pm_grid_nz == payload.distributed_gravity_state.pm_grid_nz);
  assert(restored.distributed_gravity_state.owning_rank_by_item == payload.distributed_gravity_state.owning_rank_by_item);
  assert(restored.distributed_gravity_state.pm_slab_begin_x_by_rank == payload.distributed_gravity_state.pm_slab_begin_x_by_rank);
  assert(restored.distributed_gravity_state.pm_slab_end_x_by_rank == payload.distributed_gravity_state.pm_slab_end_x_by_rank);
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
  const auto resumed_active = resumed_scheduler.beginSubstep();
  assert(particleIdsForIndices(restored.state, resumed_active) == expected_active_particle_ids);

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

  const std::filesystem::path schema_mismatch_path =
      std::filesystem::temp_directory_path() / "cosmosim_restart_schema_mismatch.hdf5";
  cosmosim::io::writeRestartCheckpointHdf5(schema_mismatch_path, payload);
  hid_t schema_file = H5Fopen(schema_mismatch_path.string().c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  assert(schema_file >= 0);
  hid_t schema_attr = H5Aopen(schema_file, "restart_schema_version", H5P_DEFAULT);
  assert(schema_attr >= 0);
  std::uint32_t bad_schema = 2;
  assert(H5Awrite(schema_attr, H5T_NATIVE_UINT32, &bad_schema) >= 0);
  H5Aclose(schema_attr);
  H5Fclose(schema_file);

  bool schema_threw = false;
  try {
    (void)cosmosim::io::readRestartCheckpointHdf5(schema_mismatch_path);
  } catch (const std::runtime_error&) {
    schema_threw = true;
  }
  assert(schema_threw);
  std::filesystem::remove(schema_mismatch_path);

  const std::filesystem::path missing_required_path =
      std::filesystem::temp_directory_path() / "cosmosim_restart_missing_required.hdf5";
  cosmosim::io::writeRestartCheckpointHdf5(missing_required_path, payload);
  hid_t missing_file = H5Fopen(missing_required_path.string().c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  assert(missing_file >= 0);
  assert(H5Ldelete(missing_file, "/scheduler/pending_bin_index", H5P_DEFAULT) >= 0);
  H5Fclose(missing_file);
  bool missing_required_threw = false;
  std::string missing_required_error;
  try {
    (void)cosmosim::io::readRestartCheckpointHdf5(missing_required_path);
  } catch (const std::runtime_error& error) {
    missing_required_threw = true;
    missing_required_error = error.what();
  }
  assert(missing_required_threw);
  const bool missing_field_named = missing_required_error.find("pending_bin_index") != std::string::npos;
  const bool integrity_guard_rejected =
      missing_required_error.find("integrity hash mismatch") != std::string::npos;
  assert(missing_field_named || integrity_guard_rejected);
  std::filesystem::remove(missing_required_path);

  const std::filesystem::path legacy_softening_path =
      std::filesystem::temp_directory_path() / "cosmosim_restart_legacy_softening_compat.hdf5";
  cosmosim::core::SimulationState legacy_softening_state = state;
  legacy_softening_state.particle_sidecar.gravity_softening_comoving.clear();
  legacy_softening_state.particle_sidecar.has_gravity_softening_override.clear();
  cosmosim::io::RestartWritePayload legacy_softening_payload = payload;
  legacy_softening_payload.state = &legacy_softening_state;
  legacy_softening_payload.provenance = cosmosim::core::makeProvenanceRecord(payload.normalized_config_hash_hex, "deadbeef");
  legacy_softening_payload.provenance.gravity_softening_policy = payload.provenance.gravity_softening_policy;
  legacy_softening_payload.provenance.gravity_softening_kernel = payload.provenance.gravity_softening_kernel;
  legacy_softening_payload.provenance.gravity_softening_epsilon_kpc_comoving =
      payload.provenance.gravity_softening_epsilon_kpc_comoving;
  cosmosim::io::writeRestartCheckpointHdf5(legacy_softening_path, legacy_softening_payload);
  const cosmosim::io::RestartReadResult legacy_softening_result =
      cosmosim::io::readRestartCheckpointHdf5(legacy_softening_path);
  assert(legacy_softening_result.state.particle_sidecar.gravity_softening_comoving.empty());
  std::filesystem::remove(legacy_softening_path);

  const std::filesystem::path finalize_failure_dir =
      std::filesystem::temp_directory_path() / "cosmosim_restart_finalize_failure_target";
  std::filesystem::remove_all(finalize_failure_dir);
  std::filesystem::create_directories(finalize_failure_dir);

  cosmosim::io::RestartWritePolicy failure_policy;
  failure_policy.temporary_suffix = ".partial";

  bool finalize_threw = false;
  try {
    cosmosim::io::writeRestartCheckpointHdf5(finalize_failure_dir, payload, failure_policy);
  } catch (const std::runtime_error&) {
    finalize_threw = true;
  }
  assert(finalize_threw);
  assert(std::filesystem::is_directory(finalize_failure_dir));
  const std::filesystem::path partial_path = finalize_failure_dir.string() + failure_policy.temporary_suffix;
  assert(std::filesystem::exists(partial_path));
  std::filesystem::remove(partial_path);
  std::filesystem::remove_all(finalize_failure_dir);
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
