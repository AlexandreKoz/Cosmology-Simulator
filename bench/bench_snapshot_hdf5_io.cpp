#include <chrono>
#include <cstddef>
#include <filesystem>
#include <iostream>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/core/provenance.hpp"
#include "cosmosim/io/snapshot_hdf5.hpp"

namespace {

void buildSyntheticState(cosmosim::core::SimulationState& state, std::size_t particle_count) {
  state.resizeParticles(particle_count);
  for (std::size_t i = 0; i < particle_count; ++i) {
    state.particles.position_x_comoving[i] = static_cast<double>(i) * 0.001;
    state.particles.position_y_comoving[i] = static_cast<double>(i) * 0.002;
    state.particles.position_z_comoving[i] = static_cast<double>(i) * 0.003;
    state.particles.velocity_x_peculiar[i] = static_cast<double>(i) * 0.01;
    state.particles.velocity_y_peculiar[i] = static_cast<double>(i) * 0.02;
    state.particles.velocity_z_peculiar[i] = static_cast<double>(i) * 0.03;
    state.particles.mass_code[i] = 1.0 + static_cast<double>(i % 64);
    state.particle_sidecar.particle_id[i] = 100000 + i;
    state.particle_sidecar.owning_rank[i] = 0;
    state.particle_sidecar.species_tag[i] =
        (i % 4 == 0) ? static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kGas)
                     : static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kDarkMatter);
  }
  state.metadata.scale_factor = 1.0;
  state.rebuildSpeciesIndex();
}

void runCase(bool use_compression) {
  cosmosim::core::SimulationConfig config;
  cosmosim::core::SimulationState state;
  constexpr std::size_t particle_count = 1u << 18;
  buildSyntheticState(state, particle_count);

  cosmosim::io::SnapshotWritePayload payload;
  payload.state = &state;
  payload.config = &config;
  payload.normalized_config_text = "schema_version=1\nmode=zoom_in\n";
  payload.provenance = cosmosim::core::makeProvenanceRecord("bench_hash", "bench_sha", 0);

  cosmosim::io::SnapshotIoPolicy policy;
  policy.enable_compression = use_compression;
  policy.compression_level = use_compression ? 2 : 0;
  policy.chunk_particle_count = 1u << 14;

  const std::filesystem::path out_path = std::filesystem::temp_directory_path() /
                                         (use_compression ? "cosmosim_bench_snapshot_deflate.hdf5"
                                                          : "cosmosim_bench_snapshot_plain.hdf5");

  const auto write_start = std::chrono::steady_clock::now();
  cosmosim::io::writeGadgetArepoSnapshotHdf5(out_path, payload, policy);
  const auto write_end = std::chrono::steady_clock::now();

  const auto read_start = std::chrono::steady_clock::now();
  const auto read_result = cosmosim::io::readGadgetArepoSnapshotHdf5(out_path, config);
  const auto read_end = std::chrono::steady_clock::now();

  const auto write_us =
      std::chrono::duration_cast<std::chrono::microseconds>(write_end - write_start).count();
  const auto read_us =
      std::chrono::duration_cast<std::chrono::microseconds>(read_end - read_start).count();

  const double approx_bytes = static_cast<double>(particle_count) * (3 + 3 + 1) * sizeof(double);
  const double write_gib_s = (approx_bytes / (1024.0 * 1024.0 * 1024.0)) / (static_cast<double>(write_us) * 1.0e-6);
  const double read_gib_s = (approx_bytes / (1024.0 * 1024.0 * 1024.0)) / (static_cast<double>(read_us) * 1.0e-6);

  std::cout << "bench_snapshot_hdf5_io"
            << " build_type=" << COSMOSIM_BUILD_TYPE
            << " threads=1"
            << " compression=" << (use_compression ? "deflate" : "none")
            << " particle_count=" << particle_count
            << " write_us=" << write_us
            << " read_us=" << read_us
            << " write_gib_s=" << write_gib_s
            << " read_gib_s=" << read_gib_s
            << " checksum=" << read_result.state.particles.size() << '\n';

  std::filesystem::remove(out_path);
}

}  // namespace

int main() {
#if COSMOSIM_ENABLE_HDF5
  runCase(false);
  runCase(true);
#else
  std::cout << "bench_snapshot_hdf5_io skipped: COSMOSIM_ENABLE_HDF5=OFF\n";
#endif
  return 0;
}
