#include <array>
#include <cassert>
#include <cstdint>
#include <filesystem>
#include <stdexcept>
#include <string>

#include <mpi.h>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/core/config.hpp"
#include "cosmosim/core/provenance.hpp"
#include "cosmosim/io/snapshot_hdf5.hpp"

#include "../support/mpi_test_workspace.hpp"

namespace {

cosmosim::core::SimulationState makeLocalState(int rank, int size) {
  cosmosim::core::SimulationState state;
  // Exercise uneven ownership and an empty participant in the collective writer.
  const std::size_t local_count = rank == size - 1
      ? 0U
      : static_cast<std::size_t>(rank % 3 + 1);
  state.resizeParticles(local_count);
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
      static_cast<std::size_t>(cosmosim::core::ParticleSpecies::kDarkMatter)] =
      static_cast<std::uint64_t>(local_count);
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
  if (size != 2 && size != 3 && size != 4 && size != 8) {
    MPI_Finalize();
    return 1;
  }

  auto workspace = cosmosim::test_support::createMpiSharedWorkspace(
      "cosmosim_snapshot_set_mpi_" + std::to_string(size));
  const std::filesystem::path& directory = workspace.root();
  MPI_Barrier(MPI_COMM_WORLD);

  auto config = cosmosim::core::makeUnvalidatedSimulationConfigForTests();
  config.output.run_name = "snapshot_set_mpi";
  config.cosmology.box_size_mpc_comoving = 32.0;
  config.cosmology.box_size_x_mpc_comoving = 32.0;
  config.cosmology.box_size_y_mpc_comoving = 32.0;
  config.cosmology.box_size_z_mpc_comoving = 32.0;
  cosmosim::core::SimulationState state = makeLocalState(rank, size);

  const std::uint64_t local_particle_count =
      static_cast<std::uint64_t>(state.particles.size());
  std::uint64_t global_particle_count = 0U;
  MPI_Allreduce(
      &local_particle_count, &global_particle_count, 1, MPI_UINT64_T, MPI_SUM,
      MPI_COMM_WORLD);
  const std::array<std::uint64_t, 6> global_counts{
      0U, global_particle_count, 0U, 0U, 0U, 0U};
  const std::string normalized = "schema_version=1\nmode=cosmo_cube\n";
  const std::string generation = "snapshot_set_mpi_generation_" + std::to_string(size);

  cosmosim::io::SnapshotWritePayload payload;
  payload.state = &state;
  payload.config = &config;
  payload.normalized_config_text = normalized;
  // Collective metadata must be identical on every rank. Production obtains
  // this by broadcasting rank-zero provenance before HDF5 creation.
  payload.provenance = cosmosim::core::makeProvenanceRecord(
      "snapshot_set_mpi_config_hash", "snapshot_set_mpi_sha", 0, normalized);
  payload.set_member.global_part_count = global_counts;
  payload.set_member.has_global_part_count = true;
  payload.set_member.generation_id = generation;

  int root_ok = 1;
#if COSMOSIM_HDF5_PARALLEL
  payload.set_member.member_index = 0U;
  payload.set_member.num_files_per_snapshot = 1U;
  payload.collective_single_file = true;
  std::uint64_t row_offset = 0U;
  MPI_Exscan(
      &local_particle_count, &row_offset, 1, MPI_UINT64_T, MPI_SUM,
      MPI_COMM_WORLD);
  if (rank == 0) row_offset = 0U;
  payload.collective_file_row_offset[1] = row_offset;

  cosmosim::io::SnapshotIoPolicy policy;
  policy.caller_managed_publication = true;
  const std::filesystem::path final_path = directory / "snap_007.hdf5";
  const std::filesystem::path partial_path = directory / "snap_007.hdf5.partial";
  cosmosim::io::writeScienceSnapshotHdf5(partial_path, payload, policy);
  cosmosim::io::verifySingleFileScienceSnapshotPartitionHdf5(
      partial_path, payload, policy);
  MPI_Barrier(MPI_COMM_WORLD);

  if (rank == 0) {
    try {
      cosmosim::io::publishSingleFileScienceSnapshot(
          partial_path, final_path, payload.set_member, false);
      cosmosim::io::writeSnapshotSetCompletionMarker(
          final_path, generation, 1U, global_counts, false,
          cosmosim::io::SnapshotCompletionIntegrityMode::kDistributedScienceReadback);
      const std::filesystem::path completion = directory / "snap_007.complete";
      const auto inspection = cosmosim::io::inspectSnapshotSet(completion);
      assert(inspection.complete);
      assert(inspection.num_files_per_snapshot == 1U);
      assert(inspection.global_part_count == global_counts);
      assert(inspection.member_paths.size() == 1U);
      cosmosim::io::validateSnapshotSetHdf5(completion).requireValid();

      const auto direct = cosmosim::io::readGadgetArepoSnapshotHdf5(
          final_path, config, cosmosim::io::SnapshotReadOptions{
              .require_complete_chui_set = false});
      assert(direct.state.particles.size() ==
             static_cast<std::size_t>(global_particle_count));
      assert(direct.state.validatePersistentParticleIds());
    } catch (...) {
      root_ok = 0;
    }
  }
#else
  // Serial-HDF5 MPI builds retain an explicit compatibility test for the
  // legacy sharded topology. The production runtime's auto/single policy
  // fails closed instead of silently selecting this branch.
  payload.set_member.member_index = static_cast<std::uint32_t>(rank);
  payload.set_member.num_files_per_snapshot = static_cast<std::uint32_t>(size);
  const std::filesystem::path member =
      directory / ("snap_007." + std::to_string(rank) + ".hdf5");
  cosmosim::io::writeScienceSnapshotHdf5(member, payload);
  MPI_Barrier(MPI_COMM_WORLD);

  if (rank == 0) {
    try {
      cosmosim::io::writeSnapshotSetCompletionMarker(
          member, generation, static_cast<std::uint32_t>(size), global_counts, false);
      const std::filesystem::path completion = directory / "snap_007.complete";
      const auto inspection = cosmosim::io::inspectSnapshotSet(completion);
      assert(inspection.complete);
      assert(inspection.num_files_per_snapshot == static_cast<std::uint32_t>(size));
      assert(inspection.global_part_count == global_counts);
      assert(inspection.member_paths.size() == static_cast<std::size_t>(size));
      cosmosim::io::validateSnapshotSetHdf5(completion).requireValid();

      const auto merged = cosmosim::io::readCosmoSimScienceSnapshotHdf5(completion, config);
      assert(merged.state.particles.size() ==
             static_cast<std::size_t>(global_particle_count));
      assert(merged.state.validatePersistentParticleIds());

      std::filesystem::remove(directory / "snap_007.1.hdf5");
      bool rejected = false;
      try {
        static_cast<void>(cosmosim::io::inspectSnapshotSet(completion));
      } catch (const std::exception&) {
        rejected = true;
      }
      assert(rejected);
    } catch (...) {
      root_ok = 0;
    }
  }
#endif
  MPI_Bcast(&root_ok, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  return root_ok == 1 ? 0 : 1;
}
