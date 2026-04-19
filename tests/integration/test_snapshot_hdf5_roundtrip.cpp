#include <cassert>
#include <array>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/core/provenance.hpp"
#include "cosmosim/io/snapshot_hdf5.hpp"

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
  state.resizeParticles(7);
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

  state.particle_sidecar.species_tag[0] = static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kDarkMatter);
  state.particle_sidecar.species_tag[1] = static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kDarkMatter);
  state.particle_sidecar.species_tag[2] = static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kGas);
  state.particle_sidecar.species_tag[3] = static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kGas);
  state.particle_sidecar.species_tag[4] = static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kStar);
  state.particle_sidecar.species_tag[5] = static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kStar);
  state.particle_sidecar.species_tag[6] = static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kTracer);
  state.tracers.resize(1);
  state.tracers.particle_index[0] = 6;
  state.tracers.parent_particle_id[0] = 1005;
  state.tracers.injection_step[0] = 11;
  state.tracers.host_cell_index[0] = 1;
  state.tracers.mass_fraction_of_host[0] = 0.25;
  state.tracers.last_host_mass_code[0] = state.cells.mass_code[1];
  state.tracers.cumulative_exchanged_mass_code[0] = 0.1;

  state.metadata.scale_factor = 0.5;
  state.metadata.run_name = "snapshot_roundtrip";
  state.species.count_by_species[static_cast<std::size_t>(cosmosim::core::ParticleSpecies::kDarkMatter)] = 2;
  state.species.count_by_species[static_cast<std::size_t>(cosmosim::core::ParticleSpecies::kGas)] = 2;
  state.species.count_by_species[static_cast<std::size_t>(cosmosim::core::ParticleSpecies::kStar)] = 2;
  state.species.count_by_species[static_cast<std::size_t>(cosmosim::core::ParticleSpecies::kTracer)] = 1;
  state.rebuildSpeciesIndex();
}

void testRoundtripMixedSpeciesSnapshot() {
  cosmosim::core::SimulationConfig config;
  config.output.run_name = "snapshot_roundtrip";

  cosmosim::core::SimulationState state;
  fillMixedSpeciesState(state);

  const std::filesystem::path snapshot_path =
      std::filesystem::temp_directory_path() / "cosmosim_snapshot_roundtrip.hdf5";

#if COSMOSIM_ENABLE_HDF5
  cosmosim::io::SnapshotWritePayload payload;
  payload.state = &state;
  payload.config = &config;
  payload.normalized_config_text = "schema_version=1\nmode=zoom_in\n";
  payload.provenance = cosmosim::core::makeProvenanceRecord("abc123", "deadbeef", 0);
  payload.provenance.gravity_treepm_pm_grid = 9;
  payload.provenance.gravity_treepm_assignment_scheme = "cic";
  payload.provenance.gravity_treepm_window_deconvolution = true;
  payload.provenance.gravity_treepm_asmth_cells = 1.75;
  payload.provenance.gravity_treepm_rcut_cells = 6.0;
  payload.provenance.gravity_treepm_mesh_spacing_mpc_comoving = 1.0 / 9.0;
  payload.provenance.gravity_treepm_split_scale_mpc_comoving = 1.75 / 9.0;
  payload.provenance.gravity_treepm_cutoff_radius_mpc_comoving = 6.0 / 9.0;
  payload.provenance.gravity_treepm_update_cadence_steps = 1;
  payload.provenance.gravity_softening_policy = "comoving_fixed";
  payload.provenance.gravity_softening_kernel = "plummer";
  payload.provenance.gravity_softening_epsilon_kpc_comoving = 1.5;
  payload.provenance.gravity_pm_fft_backend = "fftw3";

  cosmosim::io::SnapshotIoPolicy policy;
  policy.enable_compression = true;
  policy.compression_level = 1;
  policy.chunk_particle_count = 2;
  cosmosim::io::writeGadgetArepoSnapshotHdf5(snapshot_path, payload, policy);

  const cosmosim::io::SnapshotReadResult roundtrip =
      cosmosim::io::readGadgetArepoSnapshotHdf5(snapshot_path, config);
  const auto& schema = cosmosim::io::gadgetArepoSchemaMap();

  assert(roundtrip.state.particles.size() == state.particles.size());
  assert(roundtrip.state.validateUniqueParticleIds());
  assert(roundtrip.state.metadata.scale_factor == state.metadata.scale_factor);
  assert(roundtrip.normalized_config_text == payload.normalized_config_text);
  assert(roundtrip.report.schema_name == schema.schema_name);
  assert(roundtrip.report.schema_version == schema.schema_version);
  assert(roundtrip.provenance.schema_version == payload.provenance.schema_version);
  assert(roundtrip.provenance.git_sha == payload.provenance.git_sha);
  assert(roundtrip.provenance.config_hash_hex == payload.provenance.config_hash_hex);
  assert(roundtrip.provenance.enabled_features == payload.provenance.enabled_features);
  assert(roundtrip.provenance.gravity_treepm_pm_grid == payload.provenance.gravity_treepm_pm_grid);
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
      roundtrip.provenance.gravity_treepm_split_scale_mpc_comoving ==
      payload.provenance.gravity_treepm_split_scale_mpc_comoving);
  assert(
      roundtrip.provenance.gravity_treepm_cutoff_radius_mpc_comoving ==
      payload.provenance.gravity_treepm_cutoff_radius_mpc_comoving);
  assert(
      roundtrip.provenance.gravity_treepm_update_cadence_steps ==
      payload.provenance.gravity_treepm_update_cadence_steps);
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
  assert(roundtrip.state.tracers.size() == 1);
  assert(roundtrip.state.tracers.parent_particle_id[0] == 1005);
  assert(roundtrip.state.tracers.injection_step[0] == 11);
  assert(roundtrip.state.tracers.host_cell_index[0] == 1);
  assert(std::abs(roundtrip.state.tracers.mass_fraction_of_host[0] - 0.25) < 1.0e-12);
  assert(std::abs(roundtrip.state.tracers.cumulative_exchanged_mass_code[0] - 0.1) < 1.0e-12);

  double checksum_in = 0.0;
  double checksum_out = 0.0;
  for (std::size_t i = 0; i < state.particles.size(); ++i) {
    checksum_in += state.particles.position_x_comoving[i] + state.particles.mass_code[i] * 0.001;
    checksum_out += roundtrip.state.particles.position_x_comoving[i] +
                    roundtrip.state.particles.mass_code[i] * 0.001;
  }
  assert(checksum_in == checksum_out);

  std::filesystem::remove(snapshot_path);
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
  cosmosim::core::SimulationConfig config;
#if COSMOSIM_ENABLE_HDF5
  const std::filesystem::path snapshot_path =
      std::filesystem::temp_directory_path() / "cosmosim_snapshot_mass_table_only.hdf5";

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

  const auto imported = cosmosim::io::readGadgetArepoSnapshotHdf5(snapshot_path, config);
  assert(imported.state.particles.size() == 2);
  assert(imported.state.particles.mass_code[0] == 5.0);
  assert(imported.state.particles.mass_code[1] == 5.0);
  assert(containsString(imported.report.defaulted_fields, "/PartType1/Masses=MassTable"));
  std::filesystem::remove(snapshot_path);
#else
  (void)config;
#endif
}

}  // namespace

int main() {
  testRoundtripMixedSpeciesSnapshot();
  testMassTableFallbackSnapshotImport();
  return 0;
}
