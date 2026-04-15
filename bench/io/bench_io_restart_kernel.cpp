#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <vector>

#include "bench/reporting/bench_report.hpp"
#include "cosmosim/core/provenance.hpp"
#include "cosmosim/core/time_integration.hpp"
#include "cosmosim/io/restart_checkpoint.hpp"

namespace {

void seedState(cosmosim::core::SimulationState& state, std::size_t particle_count) {
  state.resizeParticles(particle_count);
  state.resizeCells(particle_count / 4 + 1);
  state.resizePatches(1);

  state.species.count_by_species = {particle_count, 0, 0, 0, 0};
  for (std::size_t i = 0; i < particle_count; ++i) {
    state.particles.position_x_comoving[i] = static_cast<double>(i) * 1.0e-4;
    state.particles.position_y_comoving[i] = static_cast<double>(i) * 2.0e-4;
    state.particles.position_z_comoving[i] = static_cast<double>(i) * 3.0e-4;
    state.particles.velocity_x_peculiar[i] = 0.01 * static_cast<double>(i % 7U);
    state.particles.velocity_y_peculiar[i] = 0.02 * static_cast<double>(i % 5U);
    state.particles.velocity_z_peculiar[i] = 0.03 * static_cast<double>(i % 3U);
    state.particles.mass_code[i] = 1.0;
    state.particles.time_bin[i] = 0;
    state.particle_sidecar.particle_id[i] = 1000000U + i;
    state.particle_sidecar.sfc_key[i] = i;
    state.particle_sidecar.species_tag[i] = static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kDarkMatter);
    state.particle_sidecar.particle_flags[i] = 0;
    state.particle_sidecar.owning_rank[i] = 0;
  }

  state.patches.patch_id[0] = 1;
  state.patches.level[0] = 0;
  state.patches.first_cell[0] = 0;
  state.patches.cell_count[0] = static_cast<std::uint32_t>(state.cells.size());
  for (std::size_t i = 0; i < state.cells.size(); ++i) {
    state.cells.patch_index[i] = 0;
  }

  state.rebuildSpeciesIndex();
}

}  // namespace

int main() {
#if COSMOSIM_ENABLE_HDF5
  constexpr std::size_t k_particle_count = 1U << 16;

  const auto execution = cosmosim::bench::defaultExecutionConfig(1, 6);

  cosmosim::core::SimulationState state;
  seedState(state, k_particle_count);

  cosmosim::core::IntegratorState integrator{};
  integrator.current_time_code = 0.4;
  integrator.step_index = 7;

  cosmosim::core::HierarchicalTimeBinScheduler scheduler(2);
  scheduler.reset(static_cast<std::uint32_t>(k_particle_count), 0, 0);

  cosmosim::io::RestartWritePayload payload;
  payload.state = &state;
  payload.integrator_state = &integrator;
  payload.scheduler = &scheduler;
  payload.normalized_config_text = "mode=zoom_in\n";
  payload.normalized_config_hash_hex = cosmosim::core::stableConfigHashHex(payload.normalized_config_text);
  state.metadata.normalized_config_hash = cosmosim::core::stableConfigHash(payload.normalized_config_text);
  state.metadata.normalized_config_hash_hex = payload.normalized_config_hash_hex;
  state.metadata.step_index = integrator.step_index;
  state.metadata.scale_factor = 1.0;
  payload.provenance = cosmosim::core::makeProvenanceRecord(payload.normalized_config_hash_hex, "bench_restart", 0);

  const std::filesystem::path checkpoint_path =
      std::filesystem::temp_directory_path() / "cosmosim_bench_restart_kernel.h5";

  for (std::size_t iter = 0; iter < execution.warmup_iterations; ++iter) {
    cosmosim::io::writeRestartCheckpointHdf5(checkpoint_path, payload);
    const auto read_result = cosmosim::io::readRestartCheckpointHdf5(checkpoint_path);
    (void)read_result;
  }

  double write_ms_accum = 0.0;
  double read_ms_accum = 0.0;
  for (std::size_t iter = 0; iter < execution.measurement_iterations; ++iter) {
    const auto write_begin = cosmosim::bench::BenchmarkClock::now();
    cosmosim::io::writeRestartCheckpointHdf5(checkpoint_path, payload);
    const auto write_end = cosmosim::bench::BenchmarkClock::now();
    write_ms_accum += cosmosim::bench::BenchmarkClock::millisecondsBetween(write_begin, write_end);

    const auto read_begin = cosmosim::bench::BenchmarkClock::now();
    const auto read_result = cosmosim::io::readRestartCheckpointHdf5(checkpoint_path);
    const auto read_end = cosmosim::bench::BenchmarkClock::now();
    (void)read_result;
    read_ms_accum += cosmosim::bench::BenchmarkClock::millisecondsBetween(read_begin, read_end);
  }

  const std::uint64_t bytes_per_iteration =
      static_cast<std::uint64_t>(k_particle_count) * (7U * sizeof(double) + sizeof(std::uint64_t) + sizeof(std::uint32_t));

  cosmosim::bench::BenchmarkReporter reporter("bench_io_restart_kernel");
  cosmosim::bench::addExecutionFields(reporter, execution);
  reporter.addField("particle_count", k_particle_count);
  reporter.addField("write_ms", write_ms_accum);
  reporter.addField("read_ms", read_ms_accum);
  reporter.addField("iterations", execution.measurement_iterations);
  cosmosim::bench::addBandwidthFields(
      reporter,
      bytes_per_iteration * execution.measurement_iterations,
      write_ms_accum,
      "write_bytes_proxy",
      "write_bandwidth_proxy_gb_s");
  cosmosim::bench::addBandwidthFields(
      reporter,
      bytes_per_iteration * execution.measurement_iterations,
      read_ms_accum,
      "read_bytes_proxy",
      "read_bandwidth_proxy_gb_s");
  reporter.write();

  std::filesystem::remove(checkpoint_path);
#else
  cosmosim::bench::BenchmarkReporter reporter("bench_io_restart_kernel");
  reporter.addField("skipped", "COSMOSIM_ENABLE_HDF5=OFF");
  reporter.write();
#endif
  return 0;
}
