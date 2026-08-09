#include <array>
#include <cassert>
#include <cstdint>
#include <filesystem>
#include <stdexcept>
#include <string>

#include <mpi.h>

#include "cosmosim/core/config.hpp"
#include "cosmosim/core/provenance.hpp"
#include "cosmosim/io/snapshot_hdf5.hpp"

namespace {

cosmosim::core::SimulationState makeLocalState(int rank) {
  cosmosim::core::SimulationState state;
  state.resizeParticles(2);
  state.metadata.scale_factor = 0.5;
  state.metadata.run_name = "snapshot_set_mpi";
  for (std::size_t i = 0; i < state.particles.size(); ++i) {
    const std::uint64_t id = 1U + static_cast<std::uint64_t>(rank) * 100U + i;
    state.particles.position_x_comoving[i] = 0.25 + static_cast<double>(rank) + 0.1 * static_cast<double>(i);
    state.particles.position_y_comoving[i] = 0.5 + 0.05 * static_cast<double>(i);
    state.particles.position_z_comoving[i] = 0.75 + 0.05 * static_cast<double>(i);
    state.particles.velocity_x_peculiar[i] = 10.0 + static_cast<double>(rank);
    state.particles.velocity_y_peculiar[i] = 20.0 + static_cast<double>(i);
    state.particles.velocity_z_peculiar[i] = -5.0;
    state.particles.mass_code[i] = 2.0;
    state.particle_sidecar.particle_id[i] = id;
    state.particle_sidecar.species_tag[i] =
        static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kDarkMatter);
    state.particle_sidecar.owning_rank[i] = rank;
  }
  state.species.count_by_species[
      static_cast<std::size_t>(cosmosim::core::ParticleSpecies::kDarkMatter)] = 2U;
  state.rebuildSpeciesIndex();
  return state;
}

}  // namespace

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  int rank = 0;
  int size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 2 && size != 4) {
    MPI_Finalize();
    return 0;
  }

  const std::filesystem::path directory =
      std::filesystem::temp_directory_path() /
      ("cosmosim_snapshot_set_mpi_" + std::to_string(size));
  if (rank == 0) {
    std::filesystem::remove_all(directory);
    std::filesystem::create_directories(directory);
  }
  MPI_Barrier(MPI_COMM_WORLD);

  auto config = cosmosim::core::makeUnvalidatedSimulationConfigForTests();
  config.output.run_name = "snapshot_set_mpi";
  config.cosmology.box_size_mpc_comoving = 32.0;
  config.cosmology.box_size_x_mpc_comoving = 32.0;
  config.cosmology.box_size_y_mpc_comoving = 32.0;
  config.cosmology.box_size_z_mpc_comoving = 32.0;
  cosmosim::core::SimulationState state = makeLocalState(rank);

  const std::array<std::uint64_t, 6> global_counts{
      0U, static_cast<std::uint64_t>(2 * size), 0U, 0U, 0U, 0U};
  const std::string normalized = "schema_version=1\nmode=cosmo_cube\n";
  const std::string generation = "snapshot_set_mpi_generation_" + std::to_string(size);

  cosmosim::io::SnapshotWritePayload payload;
  payload.state = &state;
  payload.config = &config;
  payload.normalized_config_text = normalized;
  payload.provenance = cosmosim::core::makeProvenanceRecord(
      "snapshot_set_mpi_config_hash", "snapshot_set_mpi_sha", rank, normalized);
  payload.set_member.member_index = static_cast<std::uint32_t>(rank);
  payload.set_member.num_files_per_snapshot = static_cast<std::uint32_t>(size);
  payload.set_member.global_part_count = global_counts;
  payload.set_member.has_global_part_count = true;
  payload.set_member.generation_id = generation;

  const std::filesystem::path member =
      directory / ("snap_007." + std::to_string(rank) + ".hdf5");
  cosmosim::io::writeScienceSnapshotHdf5(member, payload);
  MPI_Barrier(MPI_COMM_WORLD);

  int root_ok = 1;
  if (rank == 0) {
    try {
      cosmosim::io::writeSnapshotSetCompletionMarker(
          directory, generation, static_cast<std::uint32_t>(size), global_counts, false);
      const auto inspection = cosmosim::io::inspectSnapshotSet(directory);
      assert(inspection.complete);
      assert(inspection.num_files_per_snapshot == static_cast<std::uint32_t>(size));
      assert(inspection.global_part_count == global_counts);
      assert(inspection.member_paths.size() == static_cast<std::size_t>(size));
      cosmosim::io::validateSnapshotSetHdf5(directory).requireValid();

      const auto merged = cosmosim::io::readCosmoSimScienceSnapshotHdf5(directory, config);
      assert(merged.state.particles.size() == static_cast<std::size_t>(2 * size));
      assert(merged.state.validatePersistentParticleIds());

      // Incomplete sets must fail closed even if a completion manifest remains.
      std::filesystem::remove(directory / "snap_007.1.hdf5");
      bool rejected = false;
      try {
        static_cast<void>(cosmosim::io::inspectSnapshotSet(directory));
      } catch (const std::exception&) {
        rejected = true;
      }
      assert(rejected);
    } catch (...) {
      root_ok = 0;
    }
  }
  MPI_Bcast(&root_ok, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  if (rank == 0) std::filesystem::remove_all(directory);
  MPI_Finalize();
  return root_ok == 1 ? 0 : 1;
}
