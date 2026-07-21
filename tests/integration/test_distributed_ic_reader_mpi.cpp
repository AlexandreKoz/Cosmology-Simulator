#include <array>
#include <cassert>
#include <cstdint>
#include <filesystem>
#include <numeric>
#include <string>
#include <vector>

#include <hdf5.h>
#include <mpi.h>

#include "cosmosim/core/config.hpp"
#include "cosmosim/io/ic_reader.hpp"
#include "cosmosim/parallel/distributed_memory.hpp"

namespace {

class Hdf5Handle {
 public:
  explicit Hdf5Handle(hid_t value = -1) : m_value(value) {}
  Hdf5Handle(const Hdf5Handle&) = delete;
  Hdf5Handle& operator=(const Hdf5Handle&) = delete;
  ~Hdf5Handle() {
    if (m_value < 0) return;
    switch (H5Iget_type(m_value)) {
      case H5I_FILE: H5Fclose(m_value); break;
      case H5I_GROUP: H5Gclose(m_value); break;
      case H5I_DATASET: H5Dclose(m_value); break;
      case H5I_DATASPACE: H5Sclose(m_value); break;
      case H5I_ATTR: H5Aclose(m_value); break;
      default: break;
    }
  }
  [[nodiscard]] hid_t get() const noexcept { return m_value; }
 private:
  hid_t m_value = -1;
};

template <typename T>
void writeAttribute(hid_t group, const char* name, hid_t type, const T& value) {
  Hdf5Handle space(H5Screate(H5S_SCALAR));
  Hdf5Handle attr(H5Acreate2(group, name, type, space.get(), H5P_DEFAULT, H5P_DEFAULT));
  assert(attr.get() >= 0);
  assert(H5Awrite(attr.get(), type, &value) >= 0);
}

template <typename T, std::size_t N>
void writeArrayAttribute(hid_t group, const char* name, hid_t type, const std::array<T, N>& values) {
  hsize_t dims[1] = {N};
  Hdf5Handle space(H5Screate_simple(1, dims, nullptr));
  Hdf5Handle attr(H5Acreate2(group, name, type, space.get(), H5P_DEFAULT, H5P_DEFAULT));
  assert(attr.get() >= 0);
  assert(H5Awrite(attr.get(), type, values.data()) >= 0);
}

void writeDataset1d(hid_t group, const char* name, hid_t file_type, hid_t memory_type,
                    const void* values, std::size_t count) {
  hsize_t dims[1] = {count};
  Hdf5Handle space(H5Screate_simple(1, dims, nullptr));
  Hdf5Handle dataset(H5Dcreate2(group, name, file_type, space.get(), H5P_DEFAULT,
                                H5P_DEFAULT, H5P_DEFAULT));
  assert(dataset.get() >= 0);
  assert(H5Dwrite(dataset.get(), memory_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, values) >= 0);
}

void writeDatasetVec3(hid_t group, const char* name, const std::vector<float>& values) {
  hsize_t dims[2] = {values.size() / 3U, 3U};
  Hdf5Handle space(H5Screate_simple(2, dims, nullptr));
  Hdf5Handle dataset(H5Dcreate2(group, name, H5T_IEEE_F32LE, space.get(), H5P_DEFAULT,
                                H5P_DEFAULT, H5P_DEFAULT));
  assert(dataset.get() >= 0);
  assert(H5Dwrite(dataset.get(), H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                  values.data()) >= 0);
}

void writeMember(const std::filesystem::path& path, std::uint32_t member,
                 bool duplicate_ids) {
  Hdf5Handle file(H5Fcreate(path.string().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT));
  Hdf5Handle header(H5Gcreate2(file.get(), "/Header", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
  const std::array<std::uint32_t, 6> local{0U, 4U, 0U, 0U, 0U, 0U};
  const std::array<std::uint32_t, 6> total{0U, 8U, 0U, 0U, 0U, 0U};
  const std::array<std::uint32_t, 6> high{};
  const std::array<double, 6> mass_table{};
  writeArrayAttribute(header.get(), "NumPart_ThisFile", H5T_NATIVE_UINT32, local);
  writeArrayAttribute(header.get(), "NumPart_Total", H5T_NATIVE_UINT32, total);
  writeArrayAttribute(header.get(), "NumPart_Total_HighWord", H5T_NATIVE_UINT32, high);
  writeArrayAttribute(header.get(), "MassTable", H5T_NATIVE_DOUBLE, mass_table);
  writeAttribute(header.get(), "Time", H5T_NATIVE_DOUBLE, 1.0);
  writeAttribute(header.get(), "Redshift", H5T_NATIVE_DOUBLE, 0.0);
  writeAttribute(header.get(), "BoxSize", H5T_NATIVE_DOUBLE, 50000.0);
  writeAttribute(header.get(), "Omega0", H5T_NATIVE_DOUBLE, 0.315);
  writeAttribute(header.get(), "OmegaLambda", H5T_NATIVE_DOUBLE, 0.685);
  writeAttribute(header.get(), "HubbleParam", H5T_NATIVE_DOUBLE, 0.674);
  writeAttribute(header.get(), "NumFilesPerSnapshot", H5T_NATIVE_UINT32, 2U);

  Hdf5Handle dm(H5Gcreate2(file.get(), "/PartType1", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
  std::vector<float> coordinates;
  std::vector<float> velocities(12U, 0.0F);
  std::vector<float> masses(4U, 2.0F);
  std::vector<std::uint64_t> ids(4U);
  for (std::uint32_t local_index = 0U; local_index < 4U; ++local_index) {
    const std::uint32_t global_index = member * 4U + local_index;
    coordinates.insert(coordinates.end(), {
        3125.0F + static_cast<float>(global_index) * 6250.0F, 100.0F, 100.0F});
    ids[local_index] = duplicate_ids && global_index == 7U ? 1U : global_index + 1U;
  }
  writeDatasetVec3(dm.get(), "Coordinates", coordinates);
  writeDatasetVec3(dm.get(), "Velocities", velocities);
  writeDataset1d(dm.get(), "Masses", H5T_IEEE_F32LE, H5T_NATIVE_FLOAT,
                 masses.data(), masses.size());
  writeDataset1d(dm.get(), "ParticleIDs", H5T_STD_U64LE, H5T_NATIVE_UINT64,
                 ids.data(), ids.size());
}

}  // namespace

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  int world_rank = 0;
  int world_size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  const bool duplicate_ids = argc > 1 && std::string(argv[1]) == "duplicate";
  const auto base = std::filesystem::temp_directory_path() /
      (duplicate_ids ? "cosmosim_distributed_ic_duplicate" : "cosmosim_distributed_ic");
  const auto first = std::filesystem::path(base.string() + ".0.hdf5");
  const auto second = std::filesystem::path(base.string() + ".1.hdf5");
  if (world_rank == 0) {
    writeMember(first, 0U, duplicate_ids);
    writeMember(second, 1U, duplicate_ids);
  }
  MPI_Barrier(MPI_COMM_WORLD);

  auto config = cosmosim::core::makeUnvalidatedSimulationConfigForTests();
  config.mode.ic_convention =
      cosmosim::core::InitialConditionConvention::kGadgetArepoBridgeV1;
  config.mode.ic_staging_particle_count = 2U;
  const cosmosim::parallel::MpiContext mpi_context(true, world_size, world_rank);

  bool rejected_duplicate = false;
  try {
    const auto result = cosmosim::io::readDistributedGadgetArepoHdf5Ic(
        first, config, mpi_context,
        cosmosim::io::IcImportOptions{.chunk_particle_count = 4U});
    assert(!duplicate_ids);
    assert(result.report.already_partitioned);
    assert(result.report.counters.peak_staging_bytes < 4096U);
    assert(mpi_context.allreduceSumUint64(result.report.counters.chunks_assigned) ==
           (world_size == 1 ? 2U : 4U));
    assert(mpi_context.allreduceSumUint64(result.report.counters.records_read) == 8U);
    assert(mpi_context.allreduceSumUint64(result.report.counters.records_converted) == 8U);
    assert(mpi_context.allreduceSumUint64(result.report.counters.records_routed) == 8U);
    assert(result.state.validateOwnershipInvariants());
    for (std::uint32_t owner : result.state.particle_sidecar.owning_rank) {
      assert(owner == static_cast<std::uint32_t>(world_rank));
    }
    const std::uint64_t local_count = result.state.particles.size();
    const std::uint64_t global_count = mpi_context.allreduceSumUint64(local_count);
    assert(global_count == 8U);
    const std::uint64_t local_id_sum = std::accumulate(
        result.state.particle_sidecar.particle_id.begin(),
        result.state.particle_sidecar.particle_id.end(), std::uint64_t{0});
    assert(mpi_context.allreduceSumUint64(local_id_sum) == 36U);
    if (world_size > 1) {
      assert(local_count < global_count);
    }
  } catch (const std::runtime_error&) {
    rejected_duplicate = duplicate_ids;
  }
  if (duplicate_ids) assert(rejected_duplicate);

  MPI_Barrier(MPI_COMM_WORLD);
  if (world_rank == 0) {
    std::filesystem::remove(first);
    std::filesystem::remove(second);
  }
  MPI_Finalize();
  return 0;
}
