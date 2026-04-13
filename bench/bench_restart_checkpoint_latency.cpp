#include <chrono>
#include <cstddef>
#include <filesystem>
#include <iostream>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/core/provenance.hpp"
#include "cosmosim/io/restart_checkpoint.hpp"

namespace {

cosmosim::core::SimulationState makeBenchmarkState(std::size_t particle_count) {
  cosmosim::core::SimulationState state;
  state.resizeParticles(particle_count);
  state.resizeCells(particle_count / 4 + 1);
  state.resizePatches(1);

  state.species.count_by_species = {particle_count, 0, 0, 0, 0};
  for (std::size_t i = 0; i < particle_count; ++i) {
    state.particle_sidecar.particle_id[i] = i + 1;
    state.particle_sidecar.species_tag[i] = static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kDarkMatter);
    state.particle_sidecar.owning_rank[i] = 0;
    state.particle_sidecar.sfc_key[i] = i;
    state.particle_sidecar.particle_flags[i] = 0;
    state.particles.position_x_comoving[i] = static_cast<double>(i) * 1.0e-3;
    state.particles.position_y_comoving[i] = static_cast<double>(i) * 2.0e-3;
    state.particles.position_z_comoving[i] = static_cast<double>(i) * 3.0e-3;
    state.particles.velocity_x_peculiar[i] = 0.0;
    state.particles.velocity_y_peculiar[i] = 0.0;
    state.particles.velocity_z_peculiar[i] = 0.0;
    state.particles.mass_code[i] = 1.0;
    state.particles.time_bin[i] = 0;
  }

  state.patches.patch_id[0] = 1;
  state.patches.level[0] = 0;
  state.patches.first_cell[0] = 0;
  state.patches.cell_count[0] = static_cast<std::uint32_t>(state.cells.size());
  for (std::size_t i = 0; i < state.cells.size(); ++i) {
    state.cells.patch_index[i] = 0;
  }

  state.rebuildSpeciesIndex();
  return state;
}

}  // namespace

int main() {
#if COSMOSIM_ENABLE_HDF5
  const std::size_t particle_count = 200000;
  cosmosim::core::SimulationState state = makeBenchmarkState(particle_count);

  cosmosim::core::IntegratorState integrator_state;
  integrator_state.dt_time_code = 0.001;

  cosmosim::core::HierarchicalTimeBinScheduler scheduler(0);
  scheduler.reset(static_cast<std::uint32_t>(particle_count), 0);

  cosmosim::io::RestartWritePayload payload;
  payload.state = &state;
  payload.integrator_state = &integrator_state;
  payload.scheduler = &scheduler;
  payload.normalized_config_text = "schema_version = 1\nmode = cosmo_cube\n";
  payload.normalized_config_hash_hex = cosmosim::core::stableConfigHashHex(payload.normalized_config_text);
  payload.provenance = cosmosim::core::makeProvenanceRecord(payload.normalized_config_hash_hex, "benchmark");

  const std::filesystem::path checkpoint_path =
      std::filesystem::temp_directory_path() / "cosmosim_restart_bench.hdf5";

  const auto setup_begin = std::chrono::steady_clock::now();
  const std::uint64_t payload_hash = cosmosim::io::restartPayloadIntegrityHash(payload);
  const auto setup_end = std::chrono::steady_clock::now();

  const auto write_begin = std::chrono::steady_clock::now();
  cosmosim::io::writeRestartCheckpointHdf5(checkpoint_path, payload);
  const auto write_end = std::chrono::steady_clock::now();

  const std::uintmax_t bytes = std::filesystem::file_size(checkpoint_path);
  const std::chrono::duration<double> setup_seconds = setup_end - setup_begin;
  const std::chrono::duration<double> write_seconds = write_end - write_begin;
  const double bandwidth_mib_s = (static_cast<double>(bytes) / (1024.0 * 1024.0)) / write_seconds.count();

  std::cout << "bench=restart_checkpoint_latency\n";
  std::cout << "build_type=" << COSMOSIM_BUILD_TYPE << "\n";
  std::cout << "threads_hint=1\n";
  std::cout << "features=hdf5=" << COSMOSIM_ENABLE_HDF5 << "\n";
  std::cout << "particle_count=" << particle_count << "\n";
  std::cout << "payload_hash=" << payload_hash << "\n";
  std::cout << "setup_seconds=" << setup_seconds.count() << "\n";
  std::cout << "steady_state_write_seconds=" << write_seconds.count() << "\n";
  std::cout << "checkpoint_size_bytes=" << bytes << "\n";
  std::cout << "effective_write_bandwidth_mib_s=" << bandwidth_mib_s << "\n";

  std::filesystem::remove(checkpoint_path);
#else
  std::cout << "bench=restart_checkpoint_latency skipped (COSMOSIM_ENABLE_HDF5=OFF)\n";
#endif

  return 0;
}
