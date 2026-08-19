#include <cassert>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/core/provenance.hpp"
#include "cosmosim/io/snapshot_hdf5.hpp"
#include "cosmosim/io/io_contract.hpp"
#include "cosmosim/io/restart_checkpoint.hpp"
#include "../support/test_temp_workspace.hpp"

#if COSMOSIM_ENABLE_HDF5
#include <hdf5.h>
#endif

namespace {

[[nodiscard]] bool containsString(
    const std::vector<std::string>& values,
    std::string_view expected) {
  for (const std::string& value : values) {
    if (value == expected) {
      return true;
    }
  }
  return false;
}

void fillMixedSpeciesState(cosmosim::core::SimulationState& state) {
  state.resizeParticles(8);
  state.resizeCells(2);
  state.cells.mass_code[0] = 10.0;
  state.cells.mass_code[1] = 15.0;
  for (std::size_t i = 0; i < state.particles.size(); ++i) {
    state.particles.position_x_comoving[i] = static_cast<double>(i) * 0.1;
    state.particles.position_y_comoving[i] = static_cast<double>(i) * 0.2;
    state.particles.position_z_comoving[i] = static_cast<double>(i) * 0.3;
    state.particles.velocity_x_peculiar[i] = static_cast<double>(i) * 1.0;
    state.particles.velocity_y_peculiar[i] = static_cast<double>(i) * 2.0;
    state.particles.velocity_z_peculiar[i] = static_cast<double>(i) * 3.0;
    state.particles.mass_code[i] = 100.0 + static_cast<double>(i);
    state.particle_sidecar.particle_id[i] = 1000 + i;
    state.particle_sidecar.owning_rank[i] = 0;
  }
  state.particle_sidecar.gravity_softening_comoving = {0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06};

  state.particle_sidecar.species_tag[0] = static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kDarkMatter);
  state.particle_sidecar.species_tag[1] = static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kDarkMatter);
  state.particle_sidecar.species_tag[2] = static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kGas);
  state.particle_sidecar.species_tag[3] = static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kGas);
  state.particle_sidecar.species_tag[4] = static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kStar);
  state.particle_sidecar.species_tag[5] = static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kStar);
  state.particle_sidecar.species_tag[6] = static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kTracer);
  state.particle_sidecar.species_tag[7] = static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kBlackHole);

  state.particles.mass_code[2] = state.cells.mass_code[0];
  state.particles.mass_code[3] = state.cells.mass_code[1];
  for (std::size_t gas_row = 0; gas_row < 2; ++gas_row) {
    const std::size_t particle_index = gas_row + 2;
    state.cells.center_x_comoving[gas_row] = state.particles.position_x_comoving[particle_index];
    state.cells.center_y_comoving[gas_row] = state.particles.position_y_comoving[particle_index];
    state.cells.center_z_comoving[gas_row] = state.particles.position_z_comoving[particle_index];
    state.gas_cells.gas_cell_id[gas_row] = 9000 + gas_row;
    state.gas_cells.parent_particle_id[gas_row] = state.particle_sidecar.particle_id[particle_index];
    state.gas_cells.velocity_x_peculiar[gas_row] = state.particles.velocity_x_peculiar[particle_index];
    state.gas_cells.velocity_y_peculiar[gas_row] = state.particles.velocity_y_peculiar[particle_index];
    state.gas_cells.velocity_z_peculiar[gas_row] = state.particles.velocity_z_peculiar[particle_index];
    state.gas_cells.internal_energy_code[gas_row] = 2.5 + static_cast<double>(gas_row);
    state.gas_cells.density_code[gas_row] = 4.0 + static_cast<double>(gas_row);
    state.gas_cells.pressure_code[gas_row] = 1.5 + 0.25 * static_cast<double>(gas_row);
    state.gas_cells.temperature_code[gas_row] = 1000.0 + 100.0 * static_cast<double>(gas_row);
    state.gas_cells.sound_speed_code[gas_row] = 0.5 + 0.1 * static_cast<double>(gas_row);
    state.gas_cells.metal_mass_code[gas_row] = state.cells.mass_code[gas_row] * (0.01 + 0.01 * gas_row);
  }
  state.restoreGasCellIdentityRecords(
      {
          {.gas_cell_id = 9000, .parent_particle_id = 1002, .owning_patch_id = 0, .local_cell_row = 0},
          {.gas_cell_id = 9001, .parent_particle_id = 1003, .owning_patch_id = 0, .local_cell_row = 1},
      },
      1U);

  state.star_particles.resize(2);
  for (std::size_t star_row = 0; star_row < 2; ++star_row) {
    const std::size_t particle_index = star_row + 4;
    state.star_particles.particle_index[star_row] = static_cast<std::uint32_t>(particle_index);
    state.star_particles.formation_scale_factor[star_row] = 0.25 + 0.05 * star_row;
    state.star_particles.birth_mass_code[star_row] = state.particles.mass_code[particle_index];
    state.star_particles.metallicity_mass_fraction[star_row] = 0.015 + 0.005 * star_row;
    state.star_particles.birth_key[star_row] = 0xabc000ULL + star_row;
    state.star_particles.parent_gas_cell_id[star_row] = 9000 + star_row;
    state.star_particles.birth_tick[star_row] = 77 + star_row;
    state.star_particles.birth_ordinal[star_row] = static_cast<std::uint32_t>(star_row);
  }

  state.tracers.resize(1);
  state.tracers.particle_index[0] = 6;
  state.tracers.parent_particle_id[0] = 1005;
  state.tracers.injection_step[0] = 11;
  state.tracers.host_cell_index[0] = 1;
  state.tracers.mass_fraction_of_host[0] = 0.25;
  state.tracers.last_host_mass_code[0] = state.cells.mass_code[1];
  state.tracers.cumulative_exchanged_mass_code[0] = 0.1;

  state.black_holes.resize(1);
  state.black_holes.particle_index[0] = 7;
  state.black_holes.host_cell_index[0] = 1;
  state.black_holes.subgrid_mass_code[0] = 250.0;
  state.black_holes.accretion_rate_code[0] = 0.75;
  state.black_holes.feedback_energy_code[0] = 12.5;
  state.black_holes.eddington_ratio[0] = 0.125;
  state.black_holes.cumulative_accreted_mass_code[0] = 20.0;
  state.black_holes.cumulative_feedback_energy_code[0] = 88.0;
  state.black_holes.duty_cycle_active_time_code[0] = 3.0;
  state.black_holes.duty_cycle_total_time_code[0] = 9.0;

  state.metadata.scale_factor = 0.5;
  state.metadata.run_name = "snapshot_roundtrip";
  state.species.count_by_species[static_cast<std::size_t>(cosmosim::core::ParticleSpecies::kDarkMatter)] = 2;
  state.species.count_by_species[static_cast<std::size_t>(cosmosim::core::ParticleSpecies::kGas)] = 2;
  state.species.count_by_species[static_cast<std::size_t>(cosmosim::core::ParticleSpecies::kStar)] = 2;
  state.species.count_by_species[static_cast<std::size_t>(cosmosim::core::ParticleSpecies::kTracer)] = 1;
  state.species.count_by_species[static_cast<std::size_t>(cosmosim::core::ParticleSpecies::kBlackHole)] = 1;
  state.rebuildSpeciesIndex();
}

void testDecoupledGasIdentityRoundtrip() {
#if COSMOSIM_ENABLE_HDF5
  auto config = cosmosim::core::makeUnvalidatedSimulationConfigForTests();
  config.output.run_name = "snapshot_decoupled_gas_identity";

  cosmosim::core::SimulationState state;
  state.resizeParticles(2U);
  state.resizeCells(3U);
  for (std::size_t particle_index = 0; particle_index < 2U; ++particle_index) {
    state.particle_sidecar.particle_id[particle_index] = 5001U + particle_index;
    state.particle_sidecar.species_tag[particle_index] =
        static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kGas);
    state.particle_sidecar.owning_rank[particle_index] = 0U;
    state.particles.mass_code[particle_index] = 1.0;
  }
  for (std::size_t row = 0; row < 3U; ++row) {
    state.cells.center_x_comoving[row] = 0.1 + 0.2 * static_cast<double>(row);
    state.cells.center_y_comoving[row] = 0.2 + 0.1 * static_cast<double>(row);
    state.cells.center_z_comoving[row] = 0.3 + 0.05 * static_cast<double>(row);
    state.cells.mass_code[row] = 2.0 + static_cast<double>(row);
    state.cells.patch_index[row] = 0U;
    state.gas_cells.velocity_x_peculiar[row] = 0.01 * static_cast<double>(row + 1U);
    state.gas_cells.velocity_y_peculiar[row] = -0.02 * static_cast<double>(row + 1U);
    state.gas_cells.velocity_z_peculiar[row] = 0.03 * static_cast<double>(row + 1U);
    state.gas_cells.density_code[row] = 4.0 + static_cast<double>(row);
    state.gas_cells.internal_energy_code[row] = 5.0 + static_cast<double>(row);
    state.gas_cells.metal_mass_code[row] = 0.01 * state.cells.mass_code[row];
  }
  state.replaceGasCellIdentityRecords({
      {.gas_cell_id = 7001U,
       .parent_particle_id = std::nullopt,
       .owning_patch_id = 42U,
       .local_cell_row = 0U},
      {.gas_cell_id = 7002U,
       .parent_particle_id = 5001U,
       .owning_patch_id = 42U,
       .local_cell_row = 1U},
      {.gas_cell_id = 7003U,
       .parent_particle_id = 5001U,
       .owning_patch_id = 43U,
       .local_cell_row = 2U},
  });
  state.rebuildSpeciesIndex();

  const std::filesystem::path path =
      cosmosim::test_support::TestTempWorkspace::uniqueProcessLocalPath("cosmosim_snapshot_decoupled_gas_identity.hdf5");
  cosmosim::io::SnapshotWritePayload payload;
  payload.state = &state;
  payload.config = &config;
  payload.normalized_config_text = "schema_version = 1\n";
  payload.provenance = cosmosim::core::makeProvenanceRecord(
      "gas_identity", "gas_identity", 0, payload.normalized_config_text);
  cosmosim::io::writeGadgetArepoSnapshotHdf5(path, payload);

  const auto roundtrip =
      cosmosim::io::readGadgetArepoSnapshotHdf5(path, config);
  assert(roundtrip.report.schema_version == 6U);
  assert(roundtrip.state.cells.size() == 3U);
  assert(roundtrip.state.gas_cell_identity.size() == 3U);
  assert(!roundtrip.state.parentParticleIdForGasCellId(7001U).has_value());
  assert(roundtrip.state.parentParticleIdForGasCellId(7002U).value() == 5001U);
  assert(roundtrip.state.parentParticleIdForGasCellId(7003U).value() == 5001U);
  assert(roundtrip.state.gas_cell_identity.rowsForParentParticleId(5001U).size() == 2U);
  assert(roundtrip.state.owningPatchIdForGasCellId(7001U).value() == 42U);
  assert(roundtrip.state.owningPatchIdForGasCellId(7003U).value() == 43U);
  assert(std::abs(roundtrip.state.cells.center_x_comoving[0] - state.cells.center_x_comoving[0]) < 1.0e-14);
  assert(std::abs(roundtrip.state.cells.center_x_comoving[2] - state.cells.center_x_comoving[2]) < 1.0e-14);
  assert(std::abs(roundtrip.state.cells.mass_code[1] - state.cells.mass_code[1]) < 1.0e-14);
  assert(roundtrip.report.analysis_ready);
  assert(!roundtrip.report.evolution_ready);
  assert(!roundtrip.report.evolution_readiness_reasons.empty());
  std::filesystem::remove(path);
#endif
}

void testRoundtripMixedSpeciesSnapshot() {
  auto config = cosmosim::core::makeUnvalidatedSimulationConfigForTests();
  config.output.run_name = "snapshot_roundtrip";
  config.cosmology.box_size_x_mpc_comoving = 10.0;
  config.cosmology.box_size_y_mpc_comoving = 8.0;
  config.cosmology.box_size_z_mpc_comoving = 6.0;
  config.cosmology.box_size_mpc_comoving = config.cosmology.box_size_x_mpc_comoving;

  cosmosim::core::SimulationState state;
  fillMixedSpeciesState(state);

  const std::filesystem::path snapshot_directory =
      cosmosim::test_support::TestTempWorkspace::uniqueProcessLocalPath("cosmosim_snapshot_roundtrip_set");
  std::filesystem::remove_all(snapshot_directory);
  std::filesystem::create_directories(snapshot_directory);
  const std::filesystem::path snapshot_path = snapshot_directory / "snap_000.hdf5";

#if COSMOSIM_ENABLE_HDF5
  cosmosim::io::SnapshotWritePayload payload;
  payload.state = &state;
  payload.config = &config;
  payload.normalized_config_text = "schema_version=1\nmode=zoom_in\n";
  payload.provenance = cosmosim::core::makeProvenanceRecord("abc123", "deadbeef", 0, payload.normalized_config_text);
  payload.provenance.gravity_treepm_pm_grid = 9;
  payload.provenance.gravity_treepm_pm_grid_nx = 9;
  payload.provenance.gravity_treepm_pm_grid_ny = 7;
  payload.provenance.gravity_treepm_pm_grid_nz = 5;
  payload.provenance.gravity_treepm_assignment_scheme = "cic";
  payload.provenance.gravity_treepm_window_deconvolution = true;
  payload.provenance.gravity_treepm_asmth_cells = 1.75;
  payload.provenance.gravity_treepm_rcut_cells = 6.0;
  payload.provenance.gravity_treepm_mesh_spacing_mpc_comoving = std::cbrt((10.0 / 9.0) * (8.0 / 7.0) * (6.0 / 5.0));
  payload.provenance.gravity_treepm_mesh_spacing_x_mpc_comoving = 10.0 / 9.0;
  payload.provenance.gravity_treepm_mesh_spacing_y_mpc_comoving = 8.0 / 7.0;
  payload.provenance.gravity_treepm_mesh_spacing_z_mpc_comoving = 6.0 / 5.0;
  payload.provenance.gravity_treepm_split_scale_mpc_comoving = 1.75 * payload.provenance.gravity_treepm_mesh_spacing_mpc_comoving;
  payload.provenance.gravity_treepm_cutoff_radius_mpc_comoving = 6.0 / 9.0;
  payload.provenance.gravity_treepm_update_cadence_steps = 1;
  payload.provenance.gravity_treepm_tree_opening_criterion = "relative_force_error";
  payload.provenance.gravity_treepm_tree_opening_theta = 0.58;
  payload.provenance.gravity_treepm_tree_relative_force_tolerance = 0.003;
  payload.provenance.gravity_treepm_tree_relative_force_acceleration_floor = 1.0e-22;
  payload.provenance.gravity_softening_policy = "comoving_fixed";
  payload.provenance.gravity_softening_kernel = "plummer";
  payload.provenance.gravity_softening_epsilon_kpc_comoving = 1.5;
  payload.provenance.gravity_pm_fft_backend = "fftw3";

  cosmosim::io::SnapshotIoPolicy policy;
  policy.enable_compression = true;
  policy.compression_level = 1;
  policy.chunk_particle_count = 2;
  cosmosim::io::writeGadgetArepoSnapshotHdf5(snapshot_path, payload, policy);
  const auto pending_single = cosmosim::io::inspectSnapshotSet(snapshot_path);
  assert(!pending_single.complete);
  assert(pending_single.num_files_per_snapshot == 1U);
  assert(!pending_single.generation_id.empty());
  const std::filesystem::path single_integrity_path =
      std::filesystem::path(snapshot_path.string() + ".member.integrity");
  assert(std::filesystem::is_regular_file(single_integrity_path));
  std::filesystem::remove(single_integrity_path);
  bool missing_integrity_rejected = false;
  try {
    cosmosim::io::writeSnapshotSetCompletionMarker(
        snapshot_directory, pending_single.generation_id, 1U,
        pending_single.global_part_count, false);
  } catch (const std::runtime_error&) {
    missing_integrity_rejected = true;
  }
  assert(missing_integrity_rejected);
  cosmosim::io::writeGadgetArepoSnapshotHdf5(snapshot_path, payload, policy);
  assert(std::filesystem::is_regular_file(single_integrity_path));
  cosmosim::io::writeSnapshotSetCompletionMarker(
      snapshot_directory, pending_single.generation_id, 1U,
      pending_single.global_part_count, false);
  const auto complete_single = cosmosim::io::inspectSnapshotSet(snapshot_directory);
  assert(complete_single.complete);
  assert(std::filesystem::is_regular_file(
      snapshot_directory / (pending_single.generation_id + ".complete")));
  cosmosim::io::validateSnapshotSetHdf5(snapshot_directory).requireValid();

  hid_t inspect_file = H5Fopen(snapshot_path.string().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  assert(inspect_file >= 0);
  assert(H5Aexists(inspect_file, "cosmosim_file_kind") > 0);
  assert(H5Lexists(inspect_file, "/PartType0/StarFormationRate", H5P_DEFAULT) > 0);
  assert(H5Lexists(inspect_file, "/PartType0/ParentParticleIDs", H5P_DEFAULT) > 0);
  assert(H5Lexists(inspect_file, "/PartType0/HasParentParticle", H5P_DEFAULT) > 0);
  assert(H5Lexists(inspect_file, "/PartType0/OwningPatchIDs", H5P_DEFAULT) > 0);
  assert(H5Lexists(inspect_file, "/PartType0/ColdCloudMassFraction", H5P_DEFAULT) > 0);
  assert(H5Lexists(inspect_file, "/PartType0/EffectivePressure", H5P_DEFAULT) > 0);
  assert(H5Lexists(inspect_file, "/PartType0/EffectiveInternalEnergy", H5P_DEFAULT) > 0);
  assert(H5Lexists(inspect_file, "/PartType0/IsOnEffectiveEos", H5P_DEFAULT) > 0);
  assert(H5Lexists(inspect_file, "/PartType0/Pressure", H5P_DEFAULT) > 0);
  assert(H5Lexists(inspect_file, "/PartType0/CHUI_TemperatureCode", H5P_DEFAULT) > 0);
  assert(H5Lexists(inspect_file, "/PartType0/CHUI_SoundSpeedCode", H5P_DEFAULT) > 0);
  assert(H5Lexists(inspect_file, "/PartType5/CHUI_BHSubgridMass", H5P_DEFAULT) > 0);
  assert(H5Lexists(inspect_file, "/PartType5/CHUI_BHAccretionRateMsunPerYr", H5P_DEFAULT) > 0);
  assert(H5Lexists(inspect_file, "/PartType5/CHUI_BHEddingtonRatio", H5P_DEFAULT) > 0);
  hid_t inspect_header = H5Gopen2(inspect_file, "/Header", H5P_DEFAULT);
  assert(inspect_header >= 0);
  double box_x = 0.0;
  double box_y = 0.0;
  double box_z = 0.0;
  hid_t attr = H5Aopen(inspect_header, "CHUIBoxSizeX_MpcComoving", H5P_DEFAULT);
  assert(attr >= 0);
  assert(H5Aread(attr, H5T_NATIVE_DOUBLE, &box_x) >= 0);
  H5Aclose(attr);
  attr = H5Aopen(inspect_header, "CHUIBoxSizeY_MpcComoving", H5P_DEFAULT);
  assert(attr >= 0);
  assert(H5Aread(attr, H5T_NATIVE_DOUBLE, &box_y) >= 0);
  H5Aclose(attr);
  attr = H5Aopen(inspect_header, "CHUIBoxSizeZ_MpcComoving", H5P_DEFAULT);
  assert(attr >= 0);
  assert(H5Aread(attr, H5T_NATIVE_DOUBLE, &box_z) >= 0);
  H5Aclose(attr);
  assert(std::abs(box_x - config.cosmology.box_size_x_mpc_comoving) < 1.0e-12);
  assert(std::abs(box_y - config.cosmology.box_size_y_mpc_comoving) < 1.0e-12);
  assert(std::abs(box_z - config.cosmology.box_size_z_mpc_comoving) < 1.0e-12);
  H5Gclose(inspect_header);
  H5Fclose(inspect_file);

  const cosmosim::io::SnapshotReadResult roundtrip =
      cosmosim::io::readGadgetArepoSnapshotHdf5(snapshot_path, config);
  const auto& schema = cosmosim::io::gadgetArepoSchemaMap();

  assert(roundtrip.state.particles.size() == state.particles.size());
  assert(roundtrip.state.validateUniqueParticleIds());
  assert(roundtrip.state.metadata.scale_factor == state.metadata.scale_factor);
  assert(roundtrip.normalized_config_text == payload.normalized_config_text);
  const auto& shared_contract = cosmosim::io::sharedIoContractNames();
  assert(roundtrip.report.schema_name == schema.schema_name);
  assert(roundtrip.report.schema_version == schema.schema_version);
  assert(roundtrip.report.file_kind == shared_contract.science_snapshot_file_kind);
  assert(!roundtrip.report.restart_compatible);
  assert(roundtrip.report.header_time == state.metadata.scale_factor);
  assert(roundtrip.report.header_redshift == 1.0 / state.metadata.scale_factor - 1.0);
  assert(roundtrip.report.header_box_size_x == config.cosmology.box_size_x_mpc_comoving);
  assert(roundtrip.report.header_box_size_y == config.cosmology.box_size_y_mpc_comoving);
  assert(roundtrip.report.header_box_size_z == config.cosmology.box_size_z_mpc_comoving);
  assert(roundtrip.report.header_omega_matter == config.cosmology.omega_matter);
  assert(roundtrip.report.header_omega_lambda == config.cosmology.omega_lambda);
  assert(roundtrip.report.header_omega_baryon == config.cosmology.omega_baryon);
  assert(roundtrip.report.header_hubble_param == config.cosmology.hubble_param);
  assert(roundtrip.provenance.schema_version == payload.provenance.schema_version);
  assert(roundtrip.provenance.git_sha == payload.provenance.git_sha);
  assert(roundtrip.provenance.config_hash_hex == payload.provenance.config_hash_hex);
  assert(roundtrip.provenance.integrity_digest_algorithm == "sha256");
  assert(
      roundtrip.provenance.normalized_config_sha256_hex ==
      cosmosim::core::strongConfigHashSha256Hex(payload.normalized_config_text));
  assert(roundtrip.provenance.enabled_features == payload.provenance.enabled_features);
  assert(roundtrip.provenance.gravity_treepm_pm_grid == payload.provenance.gravity_treepm_pm_grid);
  assert(roundtrip.provenance.gravity_treepm_pm_grid_nx == payload.provenance.gravity_treepm_pm_grid_nx);
  assert(roundtrip.provenance.gravity_treepm_pm_grid_ny == payload.provenance.gravity_treepm_pm_grid_ny);
  assert(roundtrip.provenance.gravity_treepm_pm_grid_nz == payload.provenance.gravity_treepm_pm_grid_nz);
  assert(roundtrip.provenance.gravity_treepm_assignment_scheme == payload.provenance.gravity_treepm_assignment_scheme);
  assert(
      roundtrip.provenance.gravity_treepm_window_deconvolution ==
      payload.provenance.gravity_treepm_window_deconvolution);
  assert(roundtrip.provenance.gravity_treepm_asmth_cells == payload.provenance.gravity_treepm_asmth_cells);
  assert(roundtrip.provenance.gravity_treepm_rcut_cells == payload.provenance.gravity_treepm_rcut_cells);
  assert(
      roundtrip.provenance.gravity_treepm_mesh_spacing_mpc_comoving ==
      payload.provenance.gravity_treepm_mesh_spacing_mpc_comoving);
  assert(
      roundtrip.provenance.gravity_treepm_mesh_spacing_x_mpc_comoving ==
      payload.provenance.gravity_treepm_mesh_spacing_x_mpc_comoving);
  assert(
      roundtrip.provenance.gravity_treepm_mesh_spacing_y_mpc_comoving ==
      payload.provenance.gravity_treepm_mesh_spacing_y_mpc_comoving);
  assert(
      roundtrip.provenance.gravity_treepm_mesh_spacing_z_mpc_comoving ==
      payload.provenance.gravity_treepm_mesh_spacing_z_mpc_comoving);
  assert(
      roundtrip.provenance.gravity_treepm_split_scale_mpc_comoving ==
      payload.provenance.gravity_treepm_split_scale_mpc_comoving);
  assert(
      roundtrip.provenance.gravity_treepm_cutoff_radius_mpc_comoving ==
      payload.provenance.gravity_treepm_cutoff_radius_mpc_comoving);
  assert(
      roundtrip.provenance.gravity_treepm_update_cadence_steps ==
      payload.provenance.gravity_treepm_update_cadence_steps);
  assert(
      roundtrip.provenance.gravity_treepm_tree_opening_criterion ==
      payload.provenance.gravity_treepm_tree_opening_criterion);
  assert(
      roundtrip.provenance.gravity_treepm_tree_opening_theta ==
      payload.provenance.gravity_treepm_tree_opening_theta);
  assert(
      roundtrip.provenance.gravity_treepm_tree_relative_force_tolerance ==
      payload.provenance.gravity_treepm_tree_relative_force_tolerance);
  assert(
      roundtrip.provenance.gravity_treepm_tree_relative_force_acceleration_floor ==
      payload.provenance.gravity_treepm_tree_relative_force_acceleration_floor);
  assert(roundtrip.provenance.gravity_softening_policy == payload.provenance.gravity_softening_policy);
  assert(roundtrip.provenance.gravity_softening_kernel == payload.provenance.gravity_softening_kernel);
  assert(
      roundtrip.provenance.gravity_softening_epsilon_kpc_comoving ==
      payload.provenance.gravity_softening_epsilon_kpc_comoving);
  assert(roundtrip.provenance.gravity_pm_fft_backend == payload.provenance.gravity_pm_fft_backend);
  assert(containsString(roundtrip.report.present_aliases, "/PartType0/Coordinates=Coordinates"));
  assert(containsString(roundtrip.report.present_aliases, "/PartType1/Coordinates=Coordinates"));
  assert(containsString(roundtrip.report.present_aliases, "/PartType3/Coordinates=Coordinates"));
  assert(containsString(roundtrip.report.present_aliases, "/PartType4/Coordinates=Coordinates"));
  assert(containsString(roundtrip.report.present_aliases, "/PartType5/Coordinates=Coordinates"));
  assert(roundtrip.state.tracers.size() == 1);
  assert(roundtrip.state.tracers.parent_particle_id[0] == 1005);
  assert(roundtrip.state.tracers.injection_step[0] == 11);
  assert(roundtrip.state.tracers.host_cell_index[0] == 1);
  assert(std::abs(roundtrip.state.tracers.mass_fraction_of_host[0] - 0.25) < 1.0e-12);
  assert(std::abs(roundtrip.state.tracers.cumulative_exchanged_mass_code[0] - 0.1) < 1.0e-12);
  assert(roundtrip.state.cells.size() == 2);
  assert(roundtrip.state.gas_cells.gas_cell_id[0] == 9000);
  assert(roundtrip.state.gas_cells.gas_cell_id[1] == 9001);
  assert(roundtrip.state.parentParticleIdForGasCellId(9000).value() == 1002);
  assert(roundtrip.state.parentParticleIdForGasCellId(9001).value() == 1003);
  assert(std::abs(roundtrip.state.gas_cells.metal_mass_code[0] - 0.1) < 1.0e-12);
  assert(std::abs(roundtrip.state.gas_cells.metal_mass_code[1] - 0.3) < 1.0e-12);
  for (std::size_t gas_row = 0; gas_row < 2; ++gas_row) {
    assert(std::abs(roundtrip.state.gas_cells.pressure_code[gas_row] - state.gas_cells.pressure_code[gas_row]) < 1.0e-12);
    assert(std::abs(roundtrip.state.gas_cells.temperature_code[gas_row] - state.gas_cells.temperature_code[gas_row]) < 1.0e-12);
    assert(std::abs(roundtrip.state.gas_cells.sound_speed_code[gas_row] - state.gas_cells.sound_speed_code[gas_row]) < 1.0e-12);
  }
  assert(roundtrip.state.black_holes.size() == 1);
  assert(roundtrip.state.black_holes.particle_index[0] < roundtrip.state.particles.size());
  assert(roundtrip.state.black_holes.host_cell_index[0] == 1);
  assert(std::abs(roundtrip.state.black_holes.subgrid_mass_code[0] - state.black_holes.subgrid_mass_code[0]) < 1.0e-12);
  assert(std::abs(roundtrip.state.black_holes.accretion_rate_code[0] - state.black_holes.accretion_rate_code[0]) < 1.0e-12);
  assert(std::abs(roundtrip.state.black_holes.feedback_energy_code[0] - state.black_holes.feedback_energy_code[0]) < 1.0e-12);
  assert(std::abs(roundtrip.state.black_holes.eddington_ratio[0] - state.black_holes.eddington_ratio[0]) < 1.0e-12);
  assert(std::abs(roundtrip.state.black_holes.cumulative_accreted_mass_code[0] - state.black_holes.cumulative_accreted_mass_code[0]) < 1.0e-12);
  assert(std::abs(roundtrip.state.black_holes.cumulative_feedback_energy_code[0] - state.black_holes.cumulative_feedback_energy_code[0]) < 1.0e-12);
  assert(roundtrip.state.star_particles.size() == 2);
  for (std::size_t star_row = 0; star_row < 2; ++star_row) {
    assert(roundtrip.state.star_particles.birth_key[star_row] == state.star_particles.birth_key[star_row]);
    assert(roundtrip.state.star_particles.parent_gas_cell_id[star_row] ==
           state.star_particles.parent_gas_cell_id[star_row]);
    assert(roundtrip.state.star_particles.birth_tick[star_row] == state.star_particles.birth_tick[star_row]);
    assert(roundtrip.state.star_particles.birth_ordinal[star_row] == state.star_particles.birth_ordinal[star_row]);
    assert(std::abs(roundtrip.state.star_particles.birth_mass_code[star_row] -
                    state.star_particles.birth_mass_code[star_row]) < 1.0e-12);
    assert(std::abs(roundtrip.state.star_particles.metallicity_mass_fraction[star_row] -
                    state.star_particles.metallicity_mass_fraction[star_row]) < 1.0e-12);
  }
  assert(roundtrip.state.particle_sidecar.gravity_softening_comoving.size() ==
         roundtrip.state.particles.size());
  const auto assert_softening_for_id = [&](
      const std::uint64_t particle_id, const double expected_softening) {
    std::size_t roundtrip_index = roundtrip.state.particles.size();
    for (std::size_t j = 0; j < roundtrip.state.particles.size(); ++j) {
      if (roundtrip.state.particle_sidecar.particle_id[j] == particle_id) {
        roundtrip_index = j;
        break;
      }
    }
    assert(roundtrip_index < roundtrip.state.particles.size());
    assert(std::abs(
               roundtrip.state.particle_sidecar.gravity_softening_comoving[roundtrip_index] -
               expected_softening) < 1.0e-12);
  };
  for (std::size_t i = 0; i < state.particles.size(); ++i) {
    const auto species = static_cast<cosmosim::core::ParticleSpecies>(
        state.particle_sidecar.species_tag[i]);
    if (species != cosmosim::core::ParticleSpecies::kGas) {
      assert_softening_for_id(
          state.particle_sidecar.particle_id[i],
          state.particle_sidecar.gravity_softening_comoving[i]);
    }
  }
  // Gas cells are authoritative PartType0 rows and are no longer required to
  // inherit a particle-side softening value through optional parent lineage.
  // The collisionless softening diagnostic remains round-tripped exactly.

  double checksum_in = 0.0;
  double checksum_out = 0.0;
  for (std::size_t i = 0; i < state.particles.size(); ++i) {
    checksum_in += state.particles.position_x_comoving[i] + state.particles.mass_code[i] * 0.001;
    checksum_out += roundtrip.state.particles.position_x_comoving[i] +
                    roundtrip.state.particles.mass_code[i] * 0.001;
  }
  assert(checksum_in == checksum_out);

  bool restart_reader_rejected_snapshot = false;
  std::string restart_reader_error;
  try {
    (void)cosmosim::io::readRestartCheckpointHdf5(snapshot_path);
  } catch (const std::runtime_error& ex) {
    restart_reader_rejected_snapshot = true;
    restart_reader_error = ex.what();
  }
  assert(restart_reader_rejected_snapshot);
  assert(restart_reader_error.find("science_snapshot") != std::string::npos);

  hid_t legacy_file = H5Fopen(snapshot_path.string().c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  assert(legacy_file >= 0);
  hid_t legacy_provenance = H5Gopen2(legacy_file, "/Provenance", H5P_DEFAULT);
  assert(legacy_provenance >= 0);
  for (const char* attribute_name : {
           "gravity_treepm_tree_opening_criterion",
           "gravity_treepm_tree_opening_theta",
           "gravity_treepm_tree_relative_force_tolerance",
           "gravity_treepm_tree_relative_force_acceleration_floor"}) {
    assert(H5Adelete(legacy_provenance, attribute_name) >= 0);
  }
  H5Gclose(legacy_provenance);
  H5Fclose(legacy_file);

  const cosmosim::io::SnapshotReadResult legacy_roundtrip =
      cosmosim::io::readGadgetArepoSnapshotHdf5(snapshot_path, config);
  assert(legacy_roundtrip.provenance.gravity_treepm_tree_opening_criterion == "com_distance");
  assert(legacy_roundtrip.provenance.gravity_treepm_tree_opening_theta == 0.7);
  assert(legacy_roundtrip.provenance.gravity_treepm_tree_relative_force_tolerance == 0.005);
  assert(
      legacy_roundtrip.provenance.gravity_treepm_tree_relative_force_acceleration_floor ==
      1.0e-30);

  std::filesystem::remove_all(snapshot_directory);
#else
  bool threw = false;
  std::string error_message;
  try {
    cosmosim::io::SnapshotWritePayload payload;
    payload.state = &state;
    payload.config = &config;
    cosmosim::io::writeGadgetArepoSnapshotHdf5(snapshot_path, payload);
  } catch (const std::runtime_error& ex) {
    threw = true;
    error_message = ex.what();
  }
  assert(threw);
  assert(error_message.find("COSMOSIM_ENABLE_HDF5=OFF") != std::string::npos);

  threw = false;
  error_message.clear();
  try {
    static_cast<void>(cosmosim::io::readGadgetArepoSnapshotHdf5(snapshot_path, config));
  } catch (const std::runtime_error& ex) {
    threw = true;
    error_message = ex.what();
  }
  assert(threw);
  assert(error_message.find("COSMOSIM_ENABLE_HDF5=OFF") != std::string::npos);
#endif
}

void testMassTableFallbackSnapshotImport() {
  auto config = cosmosim::core::makeUnvalidatedSimulationConfigForTests();
#if COSMOSIM_ENABLE_HDF5
  const std::filesystem::path snapshot_path =
      cosmosim::test_support::TestTempWorkspace::uniqueProcessLocalPath("cosmosim_snapshot_mass_table_only.hdf5");

  hid_t file = H5Fcreate(snapshot_path.string().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  assert(file >= 0);
  hid_t header = H5Gcreate2(file, "/Header", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  assert(header >= 0);
  hid_t part1 = H5Gcreate2(file, "/PartType1", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  assert(part1 >= 0);

  const std::array<std::uint32_t, 6> counts = {0, 2, 0, 0, 0, 0};
  const std::array<std::uint32_t, 6> zeros = {0, 0, 0, 0, 0, 0};
  const std::array<double, 6> mass_table = {0.0, 5.0, 0.0, 0.0, 0.0, 0.0};
  hsize_t vec_dims[1] = {6};
  hid_t vec_space = H5Screate_simple(1, vec_dims, nullptr);
  assert(vec_space >= 0);
  hid_t attr = H5Acreate2(header, "NumPart_ThisFile", H5T_STD_U32LE, vec_space, H5P_DEFAULT, H5P_DEFAULT);
  assert(attr >= 0);
  assert(H5Awrite(attr, H5T_NATIVE_UINT32, counts.data()) >= 0);
  H5Aclose(attr);
  attr = H5Acreate2(header, "NumPart_Total", H5T_STD_U32LE, vec_space, H5P_DEFAULT, H5P_DEFAULT);
  assert(attr >= 0);
  assert(H5Awrite(attr, H5T_NATIVE_UINT32, counts.data()) >= 0);
  H5Aclose(attr);
  attr = H5Acreate2(header, "NumPart_Total_HighWord", H5T_STD_U32LE, vec_space, H5P_DEFAULT, H5P_DEFAULT);
  assert(attr >= 0);
  assert(H5Awrite(attr, H5T_NATIVE_UINT32, zeros.data()) >= 0);
  H5Aclose(attr);
  attr = H5Acreate2(header, "MassTable", H5T_IEEE_F64LE, vec_space, H5P_DEFAULT, H5P_DEFAULT);
  assert(attr >= 0);
  assert(H5Awrite(attr, H5T_NATIVE_DOUBLE, mass_table.data()) >= 0);
  H5Aclose(attr);
  H5Sclose(vec_space);

  const std::array<double, 6> coords = {0.0, 0.1, 0.2, 1.0, 1.1, 1.2};
  const std::array<double, 6> vels = {10.0, 11.0, 12.0, 20.0, 21.0, 22.0};
  const std::array<std::uint64_t, 2> ids = {101, 102};
  hsize_t coords_dims[2] = {2, 3};
  hid_t coords_space = H5Screate_simple(2, coords_dims, nullptr);
  assert(coords_space >= 0);
  hid_t dataset = H5Dcreate2(part1, "Coordinates", H5T_IEEE_F64LE, coords_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  assert(dataset >= 0);
  assert(H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, coords.data()) >= 0);
  H5Dclose(dataset);
  dataset = H5Dcreate2(part1, "Velocities", H5T_IEEE_F64LE, coords_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  assert(dataset >= 0);
  assert(H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, vels.data()) >= 0);
  H5Dclose(dataset);
  H5Sclose(coords_space);

  hsize_t id_dims[1] = {2};
  hid_t id_space = H5Screate_simple(1, id_dims, nullptr);
  assert(id_space >= 0);
  dataset = H5Dcreate2(part1, "ParticleIDs", H5T_STD_U64LE, id_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  assert(dataset >= 0);
  assert(H5Dwrite(dataset, H5T_NATIVE_UINT64, H5S_ALL, H5S_ALL, H5P_DEFAULT, ids.data()) >= 0);
  H5Dclose(dataset);
  H5Sclose(id_space);

  H5Gclose(part1);
  H5Gclose(header);
  H5Fclose(file);

  cosmosim::io::SnapshotReadOptions import_options;
  // This synthetic fixture deliberately stores internal/native code values;
  // external imports must now state that semantic contract explicitly.
  import_options.dialect = cosmosim::io::SnapshotDialect::kChuiNative;
  const auto imported = cosmosim::io::importExternalSnapshotHdf5(
      snapshot_path, config, import_options);
  assert(imported.state.particles.size() == 2);
  assert(imported.state.particles.mass_code[0] == 5.0);
  assert(imported.state.particles.mass_code[1] == 5.0);
  assert(containsString(imported.report.defaulted_fields, "/PartType1/Masses=MassTable"));
  std::filesystem::remove(snapshot_path);
#else
  (void)config;
#endif
}


void testPersistentIdAndMissingFieldContracts() {
#if COSMOSIM_ENABLE_HDF5
  auto config = cosmosim::core::makeUnvalidatedSimulationConfigForTests();
  config.output.run_name = "snapshot_missing_contracts";
  config.cosmology.box_size_x_mpc_comoving = 10.0;
  config.cosmology.box_size_y_mpc_comoving = 10.0;
  config.cosmology.box_size_z_mpc_comoving = 10.0;
  config.cosmology.box_size_mpc_comoving = 10.0;

  cosmosim::core::SimulationState invalid_ids;
  invalid_ids.resizeParticles(1U);
  invalid_ids.metadata.scale_factor = 0.5;
  invalid_ids.particles.mass_code[0] = 1.0;
  invalid_ids.particle_sidecar.particle_id[0] = 0U;
  invalid_ids.particle_sidecar.species_tag[0] =
      static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kDarkMatter);
  invalid_ids.species.count_by_species[
      static_cast<std::size_t>(cosmosim::core::ParticleSpecies::kDarkMatter)] = 1U;
  invalid_ids.rebuildSpeciesIndex();
  cosmosim::io::SnapshotWritePayload invalid_payload;
  invalid_payload.state = &invalid_ids;
  invalid_payload.config = &config;
  invalid_payload.normalized_config_text = "schema_version=1\n";
  invalid_payload.provenance = cosmosim::core::makeProvenanceRecord(
      "zero_id", "zero_id", 0, invalid_payload.normalized_config_text);
  const auto zero_path = cosmosim::test_support::TestTempWorkspace::uniqueProcessLocalPath("cosmosim_zero_id_snapshot.hdf5");
  bool zero_rejected = false;
  try {
    cosmosim::io::writeScienceSnapshotHdf5(zero_path, invalid_payload);
  } catch (const std::exception&) {
    zero_rejected = true;
  }
  assert(zero_rejected);
  std::filesystem::remove(zero_path);

  cosmosim::core::SimulationState state;
  fillMixedSpeciesState(state);
  const auto path = cosmosim::test_support::TestTempWorkspace::uniqueProcessLocalPath("cosmosim_missing_star_tracer_fields.hdf5");
  cosmosim::io::SnapshotWritePayload payload;
  payload.state = &state;
  payload.config = &config;
  payload.normalized_config_text = "schema_version=1\n";
  payload.provenance = cosmosim::core::makeProvenanceRecord(
      "missing_fields", "missing_fields", 0, payload.normalized_config_text);
  cosmosim::io::writeScienceSnapshotHdf5(path, payload);
  hid_t file = H5Fopen(path.string().c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  assert(file >= 0);
  assert(H5Ldelete(file, "/PartType4/Metallicity", H5P_DEFAULT) >= 0);
  assert(H5Ldelete(file, "/PartType3/TracerParentParticleID", H5P_DEFAULT) >= 0);
  assert(H5Fclose(file) >= 0);

  bool reject_missing = false;
  try {
    static_cast<void>(cosmosim::io::readGadgetArepoSnapshotHdf5(path, config));
  } catch (const std::runtime_error&) {
    reject_missing = true;
  }
  assert(reject_missing);
  assert(!cosmosim::io::validateSnapshotSetHdf5(path).valid);

  cosmosim::io::SnapshotReadOptions unavailable_options;
  unavailable_options.missing_field_policy =
      cosmosim::io::SnapshotMissingFieldPolicy::kMarkUnavailable;
  const auto unavailable =
      cosmosim::io::readGadgetArepoSnapshotHdf5(path, config, unavailable_options);
  assert(containsString(unavailable.report.unavailable_fields, "/PartType4/Metallicity"));
  assert(containsString(unavailable.report.unavailable_fields, "/PartType3/TracerParentParticleID"));
  assert(unavailable.report.analysis_ready);
  assert(!unavailable.report.evolution_ready);
  std::filesystem::remove(path);
#endif
}

void testSnapshotSetCompletionContract() {
#if COSMOSIM_ENABLE_HDF5
  auto config = cosmosim::core::makeUnvalidatedSimulationConfigForTests();
  config.output.run_name = "snapshot_set_contract";
  config.cosmology.box_size_x_mpc_comoving = 10.0;
  config.cosmology.box_size_y_mpc_comoving = 8.0;
  config.cosmology.box_size_z_mpc_comoving = 6.0;
  config.cosmology.box_size_mpc_comoving = 10.0;

  const std::array<std::uint64_t, 6> global_counts = {0U, 4U, 0U, 0U, 0U, 0U};
  const std::filesystem::path directory =
      cosmosim::test_support::TestTempWorkspace::uniqueProcessLocalPath("cosmosim_snapshot_set_contract");
  std::filesystem::remove_all(directory);
  std::filesystem::create_directories(directory);

  for (std::uint32_t member_index = 0; member_index < 2U; ++member_index) {
    cosmosim::core::SimulationState state;
    state.resizeParticles(2U);
    state.metadata.scale_factor = 0.5;
    for (std::size_t row = 0; row < 2U; ++row) {
      state.particles.position_x_comoving[row] = 0.5 + member_index + 0.1 * row;
      state.particles.position_y_comoving[row] = 0.5 + 0.1 * row;
      state.particles.position_z_comoving[row] = 0.5 + 0.1 * row;
      state.particles.mass_code[row] = 1.0;
      state.particle_sidecar.particle_id[row] =
          1U + static_cast<std::uint64_t>(member_index) * 100U + row;
      state.particle_sidecar.species_tag[row] =
          static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kDarkMatter);
      state.particle_sidecar.owning_rank[row] = member_index;
    }
    state.species.count_by_species[
        static_cast<std::size_t>(cosmosim::core::ParticleSpecies::kDarkMatter)] = 2U;
    state.rebuildSpeciesIndex();

    cosmosim::io::SnapshotWritePayload payload;
    payload.state = &state;
    payload.config = &config;
    payload.normalized_config_text = "schema_version=1\nmode=cosmo_cube\n";
    payload.provenance = cosmosim::core::makeProvenanceRecord(
        "snapshot_set_hash", "snapshot_set_sha", static_cast<int>(member_index),
        payload.normalized_config_text);
    payload.set_member.member_index = member_index;
    payload.set_member.num_files_per_snapshot = 2U;
    payload.set_member.global_part_count = global_counts;
    payload.set_member.has_global_part_count = true;
    payload.set_member.generation_id = "snapshot_set_generation";
    cosmosim::io::writeScienceSnapshotHdf5(
        directory / ("snap_042." + std::to_string(member_index) + ".hdf5"),
        payload);
  }

  auto before_marker = cosmosim::io::inspectSnapshotSet(directory);
  assert(before_marker.num_files_per_snapshot == 2U);
  assert(before_marker.global_part_count == global_counts);
  assert(!before_marker.complete);

  cosmosim::io::writeSnapshotSetCompletionMarker(
      directory, "snapshot_set_generation", 2U, global_counts, false);
  const auto complete = cosmosim::io::inspectSnapshotSet(directory);
  assert(complete.complete);
  assert(complete.member_paths.size() == 2U);
  cosmosim::io::validateSnapshotSetHdf5(directory).requireValid();
  const auto merged = cosmosim::io::readCosmoSimScienceSnapshotHdf5(directory, config);
  assert(merged.state.particles.size() == 4U);
  assert(merged.report.analysis_ready);
  assert(merged.report.evolution_ready);
  assert(merged.state.validatePersistentParticleIds());

  cosmosim::io::SnapshotReadOptions constrained;
  constrained.budget.max_particles = 3U;
  bool budget_rejected = false;
  try {
    static_cast<void>(cosmosim::io::readCosmoSimScienceSnapshotHdf5(
        directory, config, constrained));
  } catch (const std::length_error&) {
    budget_rejected = true;
  }
  assert(budget_rejected);

  cosmosim::io::SnapshotReadOptions byte_constrained;
  byte_constrained.budget.max_materialized_bytes = 1U;
  bool byte_budget_rejected = false;
  try {
    static_cast<void>(cosmosim::io::readCosmoSimScienceSnapshotHdf5(
        directory, config, byte_constrained));
  } catch (const std::length_error&) {
    byte_budget_rejected = true;
  }
  assert(byte_budget_rejected);

  // Marker presence alone is insufficient: it must bind the set metadata.
  const std::filesystem::path marker = directory / "snapshot_set_generation.complete";
  std::ifstream original_marker_stream(marker);
  const std::string original_marker(
      (std::istreambuf_iterator<char>(original_marker_stream)),
      std::istreambuf_iterator<char>());
  assert(!original_marker.empty());
  {
    std::ofstream stream(marker, std::ios::trunc);
    stream << "schema=chui_snapshot_set_v1\n"
           << "generation_id=snapshot_set_generation\n"
           << "num_files_per_snapshot=2\n"
           << "global_part_count=1,1,0,0,0,0\n";
  }
  assert(!cosmosim::io::inspectSnapshotSet(directory).complete);

  // Restore the valid strong manifest, then prove member-level scientific
  // disagreement is rejected independently of manifest presence.
  {
    std::ofstream stream(marker, std::ios::trunc);
    stream << original_marker;
  }
  const auto member_one = directory / "snap_042.1.hdf5";
  hid_t mismatch_file = H5Fopen(member_one.string().c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  assert(mismatch_file >= 0);
  hid_t mismatch_header = H5Gopen2(mismatch_file, "/Header", H5P_DEFAULT);
  assert(mismatch_header >= 0);
  hid_t omega_attr = H5Aopen(mismatch_header, "OmegaBaryon", H5P_DEFAULT);
  assert(omega_attr >= 0);
  double omega_baryon = 0.123456789;
  assert(H5Awrite(omega_attr, H5T_NATIVE_DOUBLE, &omega_baryon) >= 0);
  H5Aclose(omega_attr);
  H5Gclose(mismatch_header);
  H5Fclose(mismatch_file);
  bool member_mismatch_rejected = false;
  try {
    static_cast<void>(cosmosim::io::inspectSnapshotSet(directory));
  } catch (const std::runtime_error&) {
    member_mismatch_rejected = true;
  }
  assert(member_mismatch_rejected);
  std::filesystem::remove_all(directory);
#endif
}

}  // namespace

int main() {
  testDecoupledGasIdentityRoundtrip();
  testRoundtripMixedSpeciesSnapshot();
  testMassTableFallbackSnapshotImport();
  testPersistentIdAndMissingFieldContracts();
  testSnapshotSetCompletionContract();
  return 0;
}
