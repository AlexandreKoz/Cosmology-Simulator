#include <cassert>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/core/provenance.hpp"
#include "cosmosim/io/snapshot_hdf5.hpp"

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

}  // namespace

int main() {
  testRoundtripMixedSpeciesSnapshot();
  return 0;
}
