#include <array>
#include <cassert>
#include <cmath>
#include <filesystem>
#include <optional>
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

cosmosim::core::HierarchicalTimeBinScheduler makeMigrationScheduler(
    const cosmosim::core::SimulationState& state) {
  std::uint8_t max_bin = 0;
  for (const std::uint8_t bin : state.particles.time_bin) {
    if (bin > max_bin) {
      max_bin = bin;
    }
  }
  cosmosim::core::HierarchicalTimeBinScheduler scheduler(max_bin);
  scheduler.reset(static_cast<std::uint32_t>(state.particles.size()), 0U, 0U);
  for (std::uint32_t i = 0; i < state.particles.size(); ++i) {
    scheduler.setElementBin(i, state.particles.time_bin[i], scheduler.currentTick());
  }
  return scheduler;
}


void populateState(cosmosim::core::SimulationState& state) {
  state.resizeParticles(7);
  state.resizeCells(3);
  state.resizePatches(1);

  state.species.count_by_species = {1, 3, 1, 1, 1};
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
  state.particle_sidecar.gravity_softening_comoving = {0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007};
  state.particle_sidecar.has_gravity_softening_override = {1U, 1U, 1U, 1U, 1U, 1U, 1U};
  state.particle_sidecar.species_tag = {
      static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kDarkMatter),
      static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kGas),
      static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kGas),
      static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kGas),
      static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kStar),
      static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kBlackHole),
      static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kTracer)};

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

  state.tracers.resize(1);
  state.tracers.particle_index[0] = 6;
  state.tracers.parent_particle_id[0] = state.particle_sidecar.particle_id[1];
  state.tracers.injection_step[0] = 1234;
  state.tracers.host_cell_index[0] = 1;
  state.tracers.mass_fraction_of_host[0] = 0.125;
  state.tracers.last_host_mass_code[0] = 22.5;
  state.tracers.cumulative_exchanged_mass_code[0] = 0.75;

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
    state.gas_cells.velocity_x_peculiar[i] = 80.0 + static_cast<double>(i);
    state.gas_cells.velocity_y_peculiar[i] = 90.0 + static_cast<double>(i);
    state.gas_cells.velocity_z_peculiar[i] = 100.0 + static_cast<double>(i);
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

  cosmosim::core::ModuleSidecarBlock indexed_module_sidecar;
  indexed_module_sidecar.module_name = "chemistry";
  indexed_module_sidecar.schema_version = 7;
  indexed_module_sidecar.particle_indexed = true;
  indexed_module_sidecar.row_stride_bytes = 2;
  indexed_module_sidecar.required_species_mask =
      1U << static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kGas);
  indexed_module_sidecar.particle_id_by_row = {101, 102, 103};
  indexed_module_sidecar.payload = {
      std::byte{0x10},
      std::byte{0x11},
      std::byte{0x12},
      std::byte{0x13},
      std::byte{0x14},
      std::byte{0x15}};
  state.sidecars.upsert(std::move(indexed_module_sidecar));

  state.rebuildSpeciesIndex();
  state.refreshGasCellIdentityFromParticleOrder();
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


void assertParticleTimeBinsMatchScheduler(
    const cosmosim::core::SimulationState& state,
    const cosmosim::core::TimeBinPersistentState& scheduler_state) {
  assert(state.particles.time_bin.size() == scheduler_state.bin_index.size());
  for (std::size_t i = 0; i < state.particles.time_bin.size(); ++i) {
    assert(state.particles.time_bin[i] == scheduler_state.bin_index[i]);
  }
}

void assertGasCellTimeBinsMatchScheduler(
    const cosmosim::core::SimulationState& state,
    const cosmosim::core::TimeBinPersistentState& scheduler_state) {
  assert(state.cells.time_bin.size() == scheduler_state.bin_index.size());
  for (std::size_t row = 0; row < state.cells.time_bin.size(); ++row) {
    assert(state.cells.time_bin[row] == scheduler_state.bin_index[row]);
  }
}

std::vector<std::uint64_t> gasCellIdsByDenseRow(const cosmosim::core::SimulationState& state) {
  state.requireGasCellIdentityMapCoversDenseRows("gasCellIdsByDenseRow");
  std::vector<std::uint64_t> ids;
  ids.reserve(state.cells.size());
  for (std::uint32_t row = 0; row < state.cells.size(); ++row) {
    const auto gas_cell_id = state.gas_cell_identity.gasCellIdForLocalRow(row);
    assert(gas_cell_id.has_value());
    ids.push_back(*gas_cell_id);
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

void fillRestartPayload(
    cosmosim::io::RestartWritePayload& payload,
    cosmosim::core::SimulationState& state,
    cosmosim::core::IntegratorState& integrator_state,
    cosmosim::core::HierarchicalTimeBinScheduler& scheduler) {
  payload.persistent_state.simulation_state = &state;
  payload.integrator_state = &integrator_state;
  payload.scheduler = &scheduler;
  payload.normalized_config_text = "schema_version = 1\nmode = zoom_in\n";
  payload.normalized_config_hash_hex = cosmosim::core::stableConfigHashHex(payload.normalized_config_text);
  payload.provenance = cosmosim::core::makeProvenanceRecord(payload.normalized_config_hash_hex, "deadbeef");
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
  integrator_state.current_redshift = 3.0;
  integrator_state.current_hubble_rate_code = 0.071;
  integrator_state.time_si_per_code = 3.0856775814913673e16;
  integrator_state.dt_time_code = 0.001;
  integrator_state.last_drift_factor_code = 0.25;
  integrator_state.last_first_kick_factor_code = 0.125;
  integrator_state.last_second_kick_factor_code = 0.126;
  integrator_state.last_first_hubble_drag_factor = 0.997;
  integrator_state.last_second_hubble_drag_factor = 0.996;
  integrator_state.current_boundary_kind = cosmosim::core::StepBoundaryKind::kCheckpointPoint;
  integrator_state.last_completed_boundary_kind = cosmosim::core::StepBoundaryKind::kCheckpointPoint;
  integrator_state.last_completed_restart_safe = true;
  integrator_state.step_index = 77;
  integrator_state.time_bins.hierarchical_enabled = true;
  integrator_state.time_bins.active_bin = 1;
  integrator_state.time_bins.max_bin = 3;

  cosmosim::core::HierarchicalTimeBinScheduler scheduler(3);
  scheduler.reset(static_cast<std::uint32_t>(state.particles.size()), 2, 8);
  scheduler.setElementBin(1, 1, scheduler.currentTick());
  scheduler.submitCandidateBin(4, 3, cosmosim::core::TimeStepCandidateSource::kUserClamp);
  scheduler.beginSubstep();
  scheduler.endSubstep();
  cosmosim::core::syncTimeBinMirrorsFromScheduler(
      scheduler, state, cosmosim::core::TimeBinMirrorDomain::kParticles);
  const std::vector<std::uint64_t> expected_active_particle_ids =
      schedulerActiveIdsFromPersistentState(state, scheduler.exportPersistentState());

  // This scheduler deliberately differs from the particle scheduler.  It is
  // keyed by stable gas_cell_id and persists its own bins/activation state;
  // parent particle rows must not be used to reconstruct it.
  cosmosim::core::HierarchicalTimeBinScheduler gas_cell_scheduler(3);
  gas_cell_scheduler.reset(static_cast<std::uint32_t>(state.cells.size()), 1, 8);
  gas_cell_scheduler.setElementBin(0, 3, gas_cell_scheduler.currentTick());
  gas_cell_scheduler.setElementBin(1, 0, gas_cell_scheduler.currentTick());
  gas_cell_scheduler.submitCandidateBin(2, 2, cosmosim::core::TimeStepCandidateSource::kUserClamp);
  gas_cell_scheduler.beginSubstep();
  gas_cell_scheduler.endSubstep();
  assert(gas_cell_scheduler.currentTick() == scheduler.currentTick());
  cosmosim::core::syncGasCellTimeBinMirrorsFromGasCellScheduler(gas_cell_scheduler, state);
  const auto original_gas_cell_scheduler_state = gas_cell_scheduler.exportPersistentState();
  const auto original_gas_cell_scheduler_ids = gasCellIdsByDenseRow(state);

  cosmosim::io::RestartWritePayload payload;
  payload.persistent_state.simulation_state = &state;
  payload.integrator_state = &integrator_state;
  payload.scheduler = &scheduler;
  payload.gas_cell_scheduler = &gas_cell_scheduler;
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
  payload.provenance.gravity_treepm_tree_opening_criterion = "relative_force_error";
  payload.provenance.gravity_treepm_tree_opening_theta = 0.61;
  payload.provenance.gravity_treepm_tree_relative_force_tolerance = 0.0045;
  payload.provenance.gravity_treepm_tree_relative_force_acceleration_floor = 1.0e-23;
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
  payload.output_cadence_state.output_enabled = true;
  payload.output_cadence_state.write_restarts = true;
  payload.output_cadence_state.snapshot_due = false;
  payload.output_cadence_state.checkpoint_due = false;
  payload.output_cadence_state.last_completed_step_index = integrator_state.step_index;
  payload.output_cadence_state.snapshot_interval_steps = 5;
  payload.output_cadence_state.next_snapshot_step_index = 80;
  payload.output_cadence_state.snapshot_stem = "snapshot";
  payload.output_cadence_state.restart_stem = "restart";
  payload.stochastic_state.modules.push_back(cosmosim::io::StochasticModulePersistentState{
      .module_name = "star_formation",
      .schema_version = 1,
      .rng_policy = "stateless_splitmix64(seed,step_index,cell_index,rank_local_seed_offset)",
      .random_seed = 123456789ull,
      .rank_local_seed_offset = 0,
      .last_committed_step_index = integrator_state.step_index,
      .deterministic_from_serialized_inputs = true,
  });
  payload.stochastic_state.modules.push_back(cosmosim::io::StochasticModulePersistentState{
      .module_name = "stellar_feedback",
      .schema_version = 1,
      .rng_policy = "stateless_splitmix64(seed,step_index,star_index)",
      .random_seed = 42424242ull,
      .rank_local_seed_offset = 0,
      .last_committed_step_index = integrator_state.step_index,
      .deterministic_from_serialized_inputs = true,
  });

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
  assert(restored.state.particle_sidecar.has_gravity_softening_override == state.particle_sidecar.has_gravity_softening_override);
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
  assert(restored.state.tracers.particle_index == state.tracers.particle_index);
  assert(restored.state.tracers.parent_particle_id == state.tracers.parent_particle_id);
  assert(restored.state.tracers.injection_step == state.tracers.injection_step);
  assert(restored.state.tracers.host_cell_index == state.tracers.host_cell_index);
  assert(restored.state.tracers.mass_fraction_of_host == state.tracers.mass_fraction_of_host);
  assert(restored.state.tracers.last_host_mass_code == state.tracers.last_host_mass_code);
  assert(restored.state.tracers.cumulative_exchanged_mass_code == state.tracers.cumulative_exchanged_mass_code);
  assert(restored.state.species.count_by_species == state.species.count_by_species);
  assert(restored.state.particle_species_index.globalIndices(cosmosim::core::ParticleSpecies::kDarkMatter).size() == 1);
  assert(restored.state.particle_species_index.globalIndices(cosmosim::core::ParticleSpecies::kGas).size() == 3);
  assert(restored.state.particle_species_index.globalIndices(cosmosim::core::ParticleSpecies::kStar).size() == 1);
  assert(restored.state.particle_species_index.globalIndices(cosmosim::core::ParticleSpecies::kBlackHole).size() == 1);
  assert(restored.state.particle_species_index.globalIndices(cosmosim::core::ParticleSpecies::kTracer).size() == 1);
  assert(gasDensityByParticleId(restored.state) == gas_density_before);
  assert(restored.state.gasCellIdentityMatchesParticleOrder());
  assert(restored.state.gas_cells.gas_cell_id == state.gas_cells.gas_cell_id);
  assert(restored.state.gas_cells.parent_particle_id == state.gas_cells.parent_particle_id);
  assert(restored.state.gas_cell_identity.size() == state.gas_cell_identity.size());
  assert(restored.state.gasCellIdentityMapMatchesSidecarLanes());
  for (std::uint32_t row = 0; row < restored.state.cells.size(); ++row) {
    const auto* restored_record = restored.state.gas_cell_identity.findByLocalRow(row);
    const auto* original_record = state.gas_cell_identity.findByLocalRow(row);
    assert(restored_record != nullptr);
    assert(original_record != nullptr);
    assert(restored_record->gas_cell_id == original_record->gas_cell_id);
    assert(restored_record->parent_particle_id == original_record->parent_particle_id);
    assert(restored_record->owning_patch_id == original_record->owning_patch_id);
    assert(restored_record->local_cell_row == original_record->local_cell_row);
  }
  assert(restored.state.gas_cells.velocity_x_peculiar == state.gas_cells.velocity_x_peculiar);
  assert(restored.state.gas_cells.velocity_y_peculiar == state.gas_cells.velocity_y_peculiar);
  assert(restored.state.gas_cells.velocity_z_peculiar == state.gas_cells.velocity_z_peculiar);
  assertSofteningPriorityResolution(restored.state);
  assert(restored.integrator_state.step_index == integrator_state.step_index);
  assert(std::abs(restored.integrator_state.current_time_code - integrator_state.current_time_code) < 1.0e-15);
  assert(std::abs(restored.integrator_state.current_scale_factor - integrator_state.current_scale_factor) < 1.0e-15);
  assert(std::abs(restored.integrator_state.current_redshift - integrator_state.current_redshift) < 1.0e-15);
  assert(std::abs(restored.integrator_state.current_hubble_rate_code - integrator_state.current_hubble_rate_code) < 1.0e-15);
  assert(std::abs(restored.integrator_state.time_si_per_code - integrator_state.time_si_per_code) < 1.0e-15 * integrator_state.time_si_per_code);
  assert(std::abs(restored.integrator_state.last_drift_factor_code - integrator_state.last_drift_factor_code) < 1.0e-15);
  assert(std::abs(restored.integrator_state.last_first_kick_factor_code - integrator_state.last_first_kick_factor_code) < 1.0e-15);
  assert(std::abs(restored.integrator_state.last_second_kick_factor_code - integrator_state.last_second_kick_factor_code) < 1.0e-15);
  assert(std::abs(restored.integrator_state.last_first_hubble_drag_factor - integrator_state.last_first_hubble_drag_factor) < 1.0e-15);
  assert(std::abs(restored.integrator_state.last_second_hubble_drag_factor - integrator_state.last_second_hubble_drag_factor) < 1.0e-15);
  assert(restored.integrator_state.last_completed_boundary_kind == integrator_state.last_completed_boundary_kind);
  assert(restored.integrator_state.last_completed_restart_safe == integrator_state.last_completed_restart_safe);
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
  assertParticleTimeBinsMatchScheduler(restored.state, restored.scheduler_state);
  assert(restored.gas_cell_scheduler_state.max_bin == gas_cell_scheduler.maxBin());
  assert(restored.gas_cell_scheduler_state.current_tick == gas_cell_scheduler.currentTick());
  assert(restored.gas_cell_scheduler_ids == original_gas_cell_scheduler_ids);
  assert(restored.gas_cell_scheduler_state.bin_index == original_gas_cell_scheduler_state.bin_index);
  assert(
      restored.gas_cell_scheduler_state.next_activation_tick ==
      original_gas_cell_scheduler_state.next_activation_tick);
  assert(restored.gas_cell_scheduler_state.active_flag == original_gas_cell_scheduler_state.active_flag);
  assert(
      restored.gas_cell_scheduler_state.pending_bin_index ==
      original_gas_cell_scheduler_state.pending_bin_index);
  assertGasCellTimeBinsMatchScheduler(restored.state, restored.gas_cell_scheduler_state);
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
  assert(
      restored.provenance.gravity_treepm_tree_opening_criterion ==
      payload.provenance.gravity_treepm_tree_opening_criterion);
  assert(
      restored.provenance.gravity_treepm_tree_opening_theta ==
      payload.provenance.gravity_treepm_tree_opening_theta);
  assert(
      restored.provenance.gravity_treepm_tree_relative_force_tolerance ==
      payload.provenance.gravity_treepm_tree_relative_force_tolerance);
  assert(
      restored.provenance.gravity_treepm_tree_relative_force_acceleration_floor ==
      payload.provenance.gravity_treepm_tree_relative_force_acceleration_floor);
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
  assert(restored.output_cadence_state.output_enabled == payload.output_cadence_state.output_enabled);
  assert(restored.output_cadence_state.write_restarts == payload.output_cadence_state.write_restarts);
  assert(restored.output_cadence_state.snapshot_due == payload.output_cadence_state.snapshot_due);
  assert(restored.output_cadence_state.checkpoint_due == payload.output_cadence_state.checkpoint_due);
  assert(restored.output_cadence_state.last_completed_step_index == payload.output_cadence_state.last_completed_step_index);
  assert(restored.output_cadence_state.snapshot_interval_steps == payload.output_cadence_state.snapshot_interval_steps);
  assert(restored.output_cadence_state.next_snapshot_step_index == payload.output_cadence_state.next_snapshot_step_index);
  assert(restored.output_cadence_state.snapshot_stem == payload.output_cadence_state.snapshot_stem);
  assert(restored.output_cadence_state.restart_stem == payload.output_cadence_state.restart_stem);
  assert(restored.diagnostics.restart_schema_name == cosmosim::io::restartSchema().name);
  assert(restored.diagnostics.restart_schema_version == cosmosim::io::restartSchema().version);
  assert(restored.diagnostics.current_boundary_kind ==
         std::string(cosmosim::core::stepBoundaryKindName(integrator_state.current_boundary_kind)));
  assert(restored.diagnostics.last_completed_boundary_kind ==
         std::string(cosmosim::core::stepBoundaryKindName(integrator_state.last_completed_boundary_kind)));
  assert(restored.diagnostics.restart_safe);
  assert(restored.diagnostics.step_index == integrator_state.step_index);
  assert(restored.diagnostics.scheduler_current_tick == scheduler.currentTick());
  assert(restored.diagnostics.scheduler_max_bin == scheduler.maxBin());
  assert(restored.diagnostics.scheduler_element_count == state.particles.size());
  assert(restored.diagnostics.gas_cell_scheduler_current_tick == gas_cell_scheduler.currentTick());
  assert(restored.diagnostics.gas_cell_scheduler_max_bin == gas_cell_scheduler.maxBin());
  assert(restored.diagnostics.gas_cell_scheduler_element_count == state.cells.size());
  assert(restored.diagnostics.pm_cadence_steps == integrator_state.pm_sync_state.cadenceSteps());
  assert(restored.diagnostics.pm_gravity_kick_opportunity == integrator_state.pm_sync_state.gravityKickOpportunity());
  assert(restored.diagnostics.pm_field_version == integrator_state.pm_sync_state.fieldVersion());
  assert(restored.diagnostics.pm_long_range_field_valid == integrator_state.pm_long_range_field_valid);
  assert(restored.diagnostics.output_enabled == payload.output_cadence_state.output_enabled);
  assert(restored.diagnostics.output_next_snapshot_step_index == payload.output_cadence_state.next_snapshot_step_index);
  assert(restored.diagnostics.stochastic_module_count == payload.stochastic_state.modules.size());
  assert(restored.stochastic_state.modules.size() == payload.stochastic_state.modules.size());
  assert(restored.stochastic_state.modules[0].module_name == "star_formation");
  assert(restored.stochastic_state.modules[0].schema_version == 1);
  assert(restored.stochastic_state.modules[0].rng_policy == "stateless_splitmix64(seed,step_index,cell_index,rank_local_seed_offset)");
  assert(restored.stochastic_state.modules[0].random_seed == 123456789ull);
  assert(restored.stochastic_state.modules[0].rank_local_seed_offset == 0);
  assert(restored.stochastic_state.modules[0].last_committed_step_index == integrator_state.step_index);
  assert(restored.stochastic_state.modules[0].deterministic_from_serialized_inputs);
  assert(restored.stochastic_state.modules[1].module_name == "stellar_feedback");
  assert(restored.stochastic_state.modules[1].schema_version == 1);
  assert(restored.stochastic_state.modules[1].rng_policy == "stateless_splitmix64(seed,step_index,star_index)");
  assert(restored.stochastic_state.modules[1].random_seed == 42424242ull);
  assert(restored.stochastic_state.modules[1].rank_local_seed_offset == 0);
  assert(restored.stochastic_state.modules[1].last_committed_step_index == integrator_state.step_index);
  assert(restored.stochastic_state.modules[1].deterministic_from_serialized_inputs);
  assert(restored.state.sidecars.find("hydro") != nullptr);
  const cosmosim::core::ModuleSidecarBlock* hydro_sidecar = restored.state.sidecars.find("hydro");
  assert(hydro_sidecar->schema_version == 3);
  assert(hydro_sidecar->payload.size() == 3);
  assert(hydro_sidecar->payload[0] == std::byte{0x01});
  assert(hydro_sidecar->payload[1] == std::byte{0x02});
  assert(hydro_sidecar->payload[2] == std::byte{0x03});
  const cosmosim::core::ModuleSidecarBlock* chemistry_sidecar = restored.state.sidecars.find("chemistry");
  assert(chemistry_sidecar != nullptr);
  assert(chemistry_sidecar->schema_version == 7);
  assert(chemistry_sidecar->particle_indexed);
  assert(chemistry_sidecar->row_stride_bytes == 2);
  assert(chemistry_sidecar->required_species_mask ==
         (1U << static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kGas)));
  assert((chemistry_sidecar->particle_id_by_row == std::vector<std::uint64_t>{101, 102, 103}));
  assert(chemistry_sidecar->payload.size() == 6);
  assert(chemistry_sidecar->payload[0] == std::byte{0x10});
  assert(chemistry_sidecar->payload[5] == std::byte{0x15});
  assert(restored.payload_hash == cosmosim::io::restartPayloadIntegrityHash(payload));
  assert(restored.payload_hash_hex == cosmosim::io::restartPayloadIntegrityHashHex(payload));

  cosmosim::core::HierarchicalTimeBinScheduler resumed_scheduler(restored.scheduler_state.max_bin);
  resumed_scheduler.importPersistentState(restored.scheduler_state);
  assert(resumed_scheduler.currentTick() == scheduler.currentTick());
  const auto resumed_active = resumed_scheduler.beginSubstep();
  assert(particleIdsForIndices(restored.state, resumed_active) == expected_active_particle_ids);

  cosmosim::core::SimulationState stale_state = state;
  stale_state.particles.time_bin[0] = static_cast<std::uint8_t>(stale_state.particles.time_bin[0] ^ 1U);
  cosmosim::io::RestartWritePayload stale_payload = payload;
  stale_payload.persistent_state.simulation_state = &stale_state;
  bool stale_writer_threw = false;
  try {
    (void)cosmosim::io::restartPayloadIntegrityHash(stale_payload);
  } catch (const std::invalid_argument&) {
    stale_writer_threw = true;
  }
  assert(stale_writer_threw);

  const std::filesystem::path stale_mirror_path =
      std::filesystem::temp_directory_path() / "cosmosim_restart_stale_timebin_mirror.hdf5";
  cosmosim::io::writeRestartCheckpointHdf5(stale_mirror_path, payload);
  hid_t stale_file = H5Fopen(stale_mirror_path.string().c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  assert(stale_file >= 0);
  hid_t stale_dataset = H5Dopen2(stale_file, "/state/particles/time_bin", H5P_DEFAULT);
  assert(stale_dataset >= 0);
  std::vector<std::uint8_t> stale_bins(state.particles.time_bin.begin(), state.particles.time_bin.end());
  stale_bins[0] = static_cast<std::uint8_t>(stale_bins[0] ^ 1U);
  assert(H5Dwrite(stale_dataset, H5T_NATIVE_UINT8, H5S_ALL, H5S_ALL, H5P_DEFAULT, stale_bins.data()) >= 0);
  H5Dclose(stale_dataset);
  H5Fclose(stale_file);
  bool stale_reader_threw = false;
  try {
    (void)cosmosim::io::readRestartCheckpointHdf5(stale_mirror_path);
  } catch (const std::invalid_argument&) {
    stale_reader_threw = true;
  }
  assert(stale_reader_threw);
  std::filesystem::remove(stale_mirror_path);

  const std::filesystem::path invalid_scheduler_path =
      std::filesystem::temp_directory_path() / "cosmosim_restart_invalid_scheduler_active_flag.hdf5";
  cosmosim::io::writeRestartCheckpointHdf5(invalid_scheduler_path, payload);
  hid_t invalid_scheduler_file = H5Fopen(invalid_scheduler_path.string().c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  assert(invalid_scheduler_file >= 0);
  hid_t invalid_active_dataset = H5Dopen2(invalid_scheduler_file, "/scheduler/active_flag", H5P_DEFAULT);
  assert(invalid_active_dataset >= 0);
  const auto scheduler_state_for_invalid_file = scheduler.exportPersistentState();
  std::vector<std::uint8_t> invalid_active_flags(scheduler_state_for_invalid_file.active_flag.begin(),
                                                 scheduler_state_for_invalid_file.active_flag.end());
  invalid_active_flags[0] = 2;
  assert(H5Dwrite(
             invalid_active_dataset, H5T_NATIVE_UINT8, H5S_ALL, H5S_ALL, H5P_DEFAULT, invalid_active_flags.data()) >=
         0);
  H5Dclose(invalid_active_dataset);
  H5Fclose(invalid_scheduler_file);
  bool invalid_scheduler_reader_threw = false;
  try {
    (void)cosmosim::io::readRestartCheckpointHdf5(invalid_scheduler_path);
  } catch (const std::invalid_argument&) {
    invalid_scheduler_reader_threw = true;
  }
  assert(invalid_scheduler_reader_threw);
  std::filesystem::remove(invalid_scheduler_path);

  const std::filesystem::path invalid_patch_path =
      std::filesystem::temp_directory_path() / "cosmosim_restart_invalid_hydro_patch_geometry.hdf5";
  cosmosim::io::writeRestartCheckpointHdf5(invalid_patch_path, payload);
  hid_t invalid_patch_file = H5Fopen(invalid_patch_path.string().c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  assert(invalid_patch_file >= 0);
  hid_t invalid_patch_dataset = H5Dopen2(invalid_patch_file, "/state/cells/patch_index", H5P_DEFAULT);
  assert(invalid_patch_dataset >= 0);
  std::vector<std::uint32_t> invalid_patch_index(state.cells.patch_index.begin(), state.cells.patch_index.end());
  invalid_patch_index.back() = 1U;
  assert(H5Dwrite(
             invalid_patch_dataset, H5T_NATIVE_UINT32, H5S_ALL, H5S_ALL, H5P_DEFAULT,
             invalid_patch_index.data()) >= 0);
  H5Dclose(invalid_patch_dataset);
  H5Fclose(invalid_patch_file);
  bool invalid_patch_reader_threw = false;
  try {
    (void)cosmosim::io::readRestartCheckpointHdf5(invalid_patch_path);
  } catch (const std::exception& ex) {
    const std::string message = ex.what();
    invalid_patch_reader_threw =
        message.find("/state/cells/patch_index") != std::string::npos ||
        message.find("/state/patches") != std::string::npos ||
        message.find("payload integrity hash mismatch") != std::string::npos;
  }
  assert(invalid_patch_reader_threw);
  std::filesystem::remove(invalid_patch_path);

  const std::filesystem::path missing_parent_flag_path =
      std::filesystem::temp_directory_path() / "cosmosim_restart_missing_gas_identity_has_parent.hdf5";
  cosmosim::io::writeRestartCheckpointHdf5(missing_parent_flag_path, payload);
  hid_t missing_parent_flag_file = H5Fopen(missing_parent_flag_path.string().c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  assert(missing_parent_flag_file >= 0);
  assert(H5Ldelete(missing_parent_flag_file, "/state/gas_cell_identity/has_parent_particle", H5P_DEFAULT) >= 0);
  H5Fclose(missing_parent_flag_file);
  bool missing_parent_flag_threw = false;
  try {
    (void)cosmosim::io::readRestartCheckpointHdf5(missing_parent_flag_path);
  } catch (const std::runtime_error& error) {
    missing_parent_flag_threw =
        std::string(error.what()).find("/state/gas_cell_identity/has_parent_particle") != std::string::npos;
  }
  assert(missing_parent_flag_threw);
  std::filesystem::remove(missing_parent_flag_path);

  const std::filesystem::path duplicate_gas_id_path =
      std::filesystem::temp_directory_path() / "cosmosim_restart_duplicate_gas_cell_id.hdf5";
  cosmosim::io::writeRestartCheckpointHdf5(duplicate_gas_id_path, payload);
  hid_t duplicate_gas_id_file = H5Fopen(duplicate_gas_id_path.string().c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  assert(duplicate_gas_id_file >= 0);
  hid_t duplicate_gas_id_dataset =
      H5Dopen2(duplicate_gas_id_file, "/state/gas_cell_identity/gas_cell_id", H5P_DEFAULT);
  assert(duplicate_gas_id_dataset >= 0);
  std::vector<std::uint64_t> duplicate_gas_ids(state.gas_cells.gas_cell_id.begin(), state.gas_cells.gas_cell_id.end());
  duplicate_gas_ids[1] = duplicate_gas_ids[0];
  assert(H5Dwrite(
             duplicate_gas_id_dataset, H5T_NATIVE_UINT64, H5S_ALL, H5S_ALL, H5P_DEFAULT,
             duplicate_gas_ids.data()) >= 0);
  H5Dclose(duplicate_gas_id_dataset);
  H5Fclose(duplicate_gas_id_file);
  bool duplicate_gas_id_threw = false;
  try {
    (void)cosmosim::io::readRestartCheckpointHdf5(duplicate_gas_id_path);
  } catch (const std::exception& error) {
    duplicate_gas_id_threw =
        std::string(error.what()).find("gas_cell_id") != std::string::npos ||
        std::string(error.what()).find("payload integrity hash mismatch") != std::string::npos;
  }
  assert(duplicate_gas_id_threw);
  std::filesystem::remove(duplicate_gas_id_path);

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
  const bool missing_field_named =
      missing_required_error.find("/scheduler/pending_bin_index") != std::string::npos;
  assert(missing_field_named);
  std::filesystem::remove(missing_required_path);

  const std::filesystem::path missing_softening_mask_path =
      std::filesystem::temp_directory_path() / "cosmosim_restart_missing_softening_mask.hdf5";
  cosmosim::io::writeRestartCheckpointHdf5(missing_softening_mask_path, payload);
  hid_t missing_softening_file = H5Fopen(missing_softening_mask_path.string().c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  assert(missing_softening_file >= 0);
  assert(H5Ldelete(missing_softening_file, "/state/particle_sidecar/has_gravity_softening_override", H5P_DEFAULT) >= 0);
  H5Fclose(missing_softening_file);
  bool missing_softening_threw = false;
  try {
    (void)cosmosim::io::readRestartCheckpointHdf5(missing_softening_mask_path);
  } catch (const std::runtime_error& error) {
    missing_softening_threw = std::string(error.what()).find("has_gravity_softening_override") != std::string::npos;
  }
  assert(missing_softening_threw);
  std::filesystem::remove(missing_softening_mask_path);

  const std::filesystem::path legacy_softening_path =
      std::filesystem::temp_directory_path() / "cosmosim_restart_legacy_softening_compat.hdf5";
  cosmosim::core::SimulationState legacy_softening_state = state;
  legacy_softening_state.particle_sidecar.gravity_softening_comoving.clear();
  legacy_softening_state.particle_sidecar.has_gravity_softening_override.clear();
  cosmosim::io::RestartWritePayload legacy_softening_payload = payload;
  legacy_softening_payload.persistent_state.simulation_state = &legacy_softening_state;
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

void testRestartAfterReorderAndMigration() {
#if COSMOSIM_ENABLE_HDF5
  cosmosim::core::SimulationState state;
  populateState(state);

  state.particle_sidecar.sfc_key = {10, 20, 21, 22, 1, 2, 23};
  const auto reorder = cosmosim::core::buildParticleReorderMap(state, cosmosim::core::ParticleReorderMode::kBySfcKey);
  cosmosim::core::reorderParticles(state, reorder);
  assert(state.validateOwnershipInvariants());

  cosmosim::core::IntegratorState integrator_state;
  integrator_state.current_time_code = 0.25;
  integrator_state.current_scale_factor = 0.5;
  integrator_state.dt_time_code = 0.002;
  integrator_state.step_index = 88;
  cosmosim::core::HierarchicalTimeBinScheduler scheduler(3);
  scheduler.reset(static_cast<std::uint32_t>(state.particles.size()), 2, 10);
  cosmosim::core::syncTimeBinMirrorsFromScheduler(
      scheduler, state, cosmosim::core::TimeBinMirrorDomain::kParticlesAndCells);

  cosmosim::io::RestartWritePayload reorder_payload;
  fillRestartPayload(reorder_payload, state, integrator_state, scheduler);
  const auto expected_density_by_id = gasDensityByParticleId(state);
  const std::filesystem::path reorder_path =
      std::filesystem::temp_directory_path() / "cosmosim_restart_after_reorder.hdf5";
  cosmosim::io::writeRestartCheckpointHdf5(reorder_path, reorder_payload);
  const auto reordered_restore = cosmosim::io::readRestartCheckpointHdf5(reorder_path);
  assertParticleTimeBinsMatchScheduler(reordered_restore.state, reordered_restore.scheduler_state);
  assert(reordered_restore.state.particle_sidecar.particle_id == state.particle_sidecar.particle_id);
  assert(reordered_restore.state.star_particles.particle_index == state.star_particles.particle_index);
  assert(reordered_restore.state.black_holes.particle_index == state.black_holes.particle_index);
  assert(gasDensityByParticleId(reordered_restore.state) == expected_density_by_id);
  assert(reordered_restore.state.particle_sidecar.has_gravity_softening_override == state.particle_sidecar.has_gravity_softening_override);
  std::filesystem::remove(reorder_path);

  const std::uint32_t black_hole_index = state.black_holes.particle_index[0];
  auto migrated_records = state.packParticleMigrationRecords(std::array<std::uint32_t, 1>{black_hole_index}, makeMigrationScheduler(state));
  assert(migrated_records.size() == 1);
  migrated_records[0].owning_rank = 0;
  cosmosim::core::ParticleMigrationCommit commit;
  commit.world_rank = 0;
  commit.outbound_local_indices = {black_hole_index};
  commit.inbound_records = migrated_records;
  state.commitParticleMigration(commit);
  assert(state.validateOwnershipInvariants());
  scheduler.reset(static_cast<std::uint32_t>(state.particles.size()), 2, 12);
  cosmosim::core::syncTimeBinMirrorsFromScheduler(
      scheduler, state, cosmosim::core::TimeBinMirrorDomain::kParticlesAndCells);

  cosmosim::io::RestartWritePayload migration_payload;
  fillRestartPayload(migration_payload, state, integrator_state, scheduler);
  const std::filesystem::path migration_path =
      std::filesystem::temp_directory_path() / "cosmosim_restart_after_migration.hdf5";
  cosmosim::io::writeRestartCheckpointHdf5(migration_path, migration_payload);
  const auto migration_restore = cosmosim::io::readRestartCheckpointHdf5(migration_path);
  assert(migration_restore.state.validateOwnershipInvariants());
  assertParticleTimeBinsMatchScheduler(migration_restore.state, migration_restore.scheduler_state);
  assert(migration_restore.state.particle_sidecar.particle_id == state.particle_sidecar.particle_id);
  assert(migration_restore.state.black_holes.particle_index == state.black_holes.particle_index);
  assert(migration_restore.state.particle_sidecar.gravity_softening_comoving == state.particle_sidecar.gravity_softening_comoving);
  assert(migration_restore.state.particle_sidecar.has_gravity_softening_override == state.particle_sidecar.has_gravity_softening_override);
  assert(migration_restore.state.tracers.particle_index == state.tracers.particle_index);
  assert(migration_restore.state.tracers.parent_particle_id == state.tracers.parent_particle_id);
  assert(migration_restore.state.tracers.host_cell_index == state.tracers.host_cell_index);
  assert(migration_restore.distributed_gravity_state.owning_rank_by_item == migration_payload.distributed_gravity_state.owning_rank_by_item);
  std::filesystem::remove(migration_path);
#endif
}

void testParentlessGasCellRestartRoundtrip() {
#if COSMOSIM_ENABLE_HDF5
  cosmosim::core::SimulationState state;
  state.resizeParticles(0);
  state.resizeCells(2);
  state.resizePatches(1);
  state.patches.patch_id[0] = 9100;
  state.patches.level[0] = 2;
  state.patches.first_cell[0] = 0;
  state.patches.cell_count[0] = 2;
  state.patches.owning_rank[0] = 0;
  state.cells.patch_index[0] = 0;
  state.cells.patch_index[1] = 0;
  state.cells.center_x_comoving = {0.25, 0.75};
  state.cells.center_y_comoving = {0.5, 0.5};
  state.cells.center_z_comoving = {0.125, 0.875};
  state.cells.mass_code = {2.0, 3.0};
  state.gas_cell_identity.assign({
      {.gas_cell_id = 8101, .parent_particle_id = std::nullopt, .owning_patch_id = 9100, .local_cell_row = 0},
      {.gas_cell_id = 8102, .parent_particle_id = std::nullopt, .owning_patch_id = 9100, .local_cell_row = 1},
  });
  state.gas_cells.gas_cell_id = {8101, 8102};
  state.gas_cells.parent_particle_id = {0, 0};
  state.gas_cells.velocity_x_peculiar = {1.25, -0.5};
  state.gas_cells.velocity_y_peculiar = {0.75, 2.5};
  state.gas_cells.velocity_z_peculiar = {-1.0, 0.125};
  state.gas_cells.density_code = {4.0, 5.0};
  state.gas_cells.pressure_code = {6.0, 7.0};
  state.gas_cells.internal_energy_code = {8.0, 9.0};
  state.gas_cells.temperature_code = {10.0, 11.0};
  state.gas_cells.sound_speed_code = {12.0, 13.0};
  state.rebuildSpeciesIndex();
  assert(state.validateOwnershipInvariants());

  cosmosim::core::IntegratorState integrator_state;
  integrator_state.current_boundary_kind = cosmosim::core::StepBoundaryKind::kCheckpointPoint;
  integrator_state.last_completed_boundary_kind = cosmosim::core::StepBoundaryKind::kCheckpointPoint;
  integrator_state.last_completed_restart_safe = true;
  integrator_state.step_index = 9;
  cosmosim::core::HierarchicalTimeBinScheduler scheduler(0);
  scheduler.reset(0, 0, 0);
  cosmosim::io::RestartWritePayload payload;
  fillRestartPayload(payload, state, integrator_state, scheduler);

  const std::filesystem::path checkpoint_path =
      std::filesystem::temp_directory_path() / "cosmosim_restart_parentless_gas_cells.hdf5";
  cosmosim::io::writeRestartCheckpointHdf5(checkpoint_path, payload);
  const auto restored = cosmosim::io::readRestartCheckpointHdf5(checkpoint_path);
  assert(restored.state.validateOwnershipInvariants());
  assert(restored.state.gas_cells.gas_cell_id == state.gas_cells.gas_cell_id);
  assert(restored.state.gas_cells.parent_particle_id == state.gas_cells.parent_particle_id);
  assert(restored.state.gas_cells.velocity_x_peculiar == state.gas_cells.velocity_x_peculiar);
  assert(restored.state.gas_cells.velocity_y_peculiar == state.gas_cells.velocity_y_peculiar);
  assert(restored.state.gas_cells.velocity_z_peculiar == state.gas_cells.velocity_z_peculiar);
  for (std::uint32_t row = 0; row < restored.state.cells.size(); ++row) {
    const auto* record = restored.state.gas_cell_identity.findByLocalRow(row);
    assert(record != nullptr);
    assert(!record->parent_particle_id.has_value());
    assert(record->owning_patch_id == 9100);
  }
  std::filesystem::remove(checkpoint_path);
#endif
}

}  // namespace

int main() {
  testRestartRoundtrip();
  testRestartAfterReorderAndMigration();
  testParentlessGasCellRestartRoundtrip();
  return 0;
}
