#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>

#include <hdf5.h>

#include "cosmosim/core/config.hpp"
#include "cosmosim/core/profiling.hpp"
#include "cosmosim/io/ic_reader.hpp"
#include "cosmosim/parallel/distributed_memory.hpp"
#include "cosmosim/workflows/runtime_services.hpp"
#include "workflows/internal/initial_condition_runtime.hpp"

namespace {

class Hdf5Handle {
 public:
  explicit Hdf5Handle(hid_t value = -1) : m_value(value) {}
  Hdf5Handle(const Hdf5Handle&) = delete;
  Hdf5Handle& operator=(const Hdf5Handle&) = delete;
  Hdf5Handle(Hdf5Handle&& other) noexcept : m_value(other.m_value) {
    other.m_value = -1;
  }
  ~Hdf5Handle() {
    if (m_value < 0) return;
    const H5I_type_t type = H5Iget_type(m_value);
    if (type == H5I_FILE) H5Fclose(m_value);
    else if (type == H5I_GROUP) H5Gclose(m_value);
    else if (type == H5I_DATASET) H5Dclose(m_value);
    else if (type == H5I_DATASPACE) H5Sclose(m_value);
    else if (type == H5I_ATTR) H5Aclose(m_value);
  }
  [[nodiscard]] hid_t get() const noexcept { return m_value; }

 private:
  hid_t m_value = -1;
};

template <typename T>
void writeAttributeScalar(
    hid_t group, const char* name, hid_t file_type, hid_t memory_type,
    const T& value) {
  Hdf5Handle space(H5Screate(H5S_SCALAR));
  Hdf5Handle attribute(H5Acreate2(
      group, name, file_type, space.get(), H5P_DEFAULT, H5P_DEFAULT));
  assert(attribute.get() >= 0);
  assert(H5Awrite(attribute.get(), memory_type, &value) >= 0);
}

template <typename T>
void writeAttributeArray6(
    hid_t group, const char* name, hid_t file_type, hid_t memory_type,
    const std::array<T, 6>& values) {
  hsize_t dimensions[1]{6U};
  Hdf5Handle space(H5Screate_simple(1, dimensions, nullptr));
  Hdf5Handle attribute(H5Acreate2(
      group, name, file_type, space.get(), H5P_DEFAULT, H5P_DEFAULT));
  assert(attribute.get() >= 0);
  assert(H5Awrite(attribute.get(), memory_type, values.data()) >= 0);
}

template <typename T>
void writeDataset1d(
    hid_t group, const char* name, hid_t file_type, hid_t memory_type,
    const T* values, std::size_t count) {
  hsize_t dimensions[1]{static_cast<hsize_t>(count)};
  Hdf5Handle space(H5Screate_simple(1, dimensions, nullptr));
  Hdf5Handle dataset(H5Dcreate2(
      group, name, file_type, space.get(), H5P_DEFAULT, H5P_DEFAULT,
      H5P_DEFAULT));
  assert(dataset.get() >= 0);
  assert(H5Dwrite(
      dataset.get(), memory_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, values) >=
      0);
}

void writeGasSource(const std::filesystem::path& path) {
  Hdf5Handle file(H5Fcreate(
      path.string().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT));
  Hdf5Handle header(H5Gcreate2(
      file.get(), "/Header", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
  std::array<std::int32_t, 6> local{1, 0, 0, 0, 0, 0};
  std::array<std::uint32_t, 6> total{1, 0, 0, 0, 0, 0};
  std::array<std::uint32_t, 6> high{};
  std::array<double, 6> mass_table{};
  writeAttributeArray6(
      header.get(), "NumPart_ThisFile", H5T_STD_I32LE, H5T_NATIVE_INT32,
      local);
  writeAttributeArray6(
      header.get(), "NumPart_Total", H5T_STD_U32LE, H5T_NATIVE_UINT32,
      total);
  writeAttributeArray6(
      header.get(), "NumPart_Total_HighWord", H5T_STD_U32LE,
      H5T_NATIVE_UINT32, high);
  writeAttributeArray6(
      header.get(), "MassTable", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE,
      mass_table);
  const double time = 1.0;
  const double redshift = 0.0;
  const double box_size = 50.0;
  const double omega_matter = 0.315;
  const double omega_lambda = 0.685;
  const double hubble_param = 0.674;
  const std::int32_t file_count = 1;
  writeAttributeScalar(
      header.get(), "Time", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, time);
  writeAttributeScalar(
      header.get(), "Redshift", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE,
      redshift);
  writeAttributeScalar(
      header.get(), "BoxSize", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE,
      box_size);
  writeAttributeScalar(
      header.get(), "Omega0", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE,
      omega_matter);
  writeAttributeScalar(
      header.get(), "OmegaLambda", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE,
      omega_lambda);
  writeAttributeScalar(
      header.get(), "HubbleParam", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE,
      hubble_param);
  writeAttributeScalar(
      header.get(), "NumFilesPerSnapshot", H5T_STD_I32LE,
      H5T_NATIVE_INT32, file_count);

  Hdf5Handle gas(H5Gcreate2(
      file.get(), "/PartType0", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
  const std::array<double, 3> coordinates{1.0, 2.0, 3.0};
  const std::array<double, 3> velocities{4.0, 0.0, 0.0};
  hsize_t vector_dimensions[2]{1U, 3U};
  Hdf5Handle vector_space(H5Screate_simple(2, vector_dimensions, nullptr));
  Hdf5Handle position(H5Dcreate2(
      gas.get(), "Coordinates", H5T_IEEE_F64LE, vector_space.get(),
      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
  Hdf5Handle velocity(H5Dcreate2(
      gas.get(), "Velocities", H5T_IEEE_F64LE, vector_space.get(),
      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
  assert(H5Dwrite(
      position.get(), H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
      coordinates.data()) >= 0);
  assert(H5Dwrite(
      velocity.get(), H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
      velocities.data()) >= 0);
  const double mass = 2.0;
  const double internal_energy = 5.0;
  const std::uint64_t id = 101U;
  writeDataset1d(
      gas.get(), "Masses", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, &mass, 1U);
  writeDataset1d(
      gas.get(), "ParticleIDs", H5T_STD_U64LE, H5T_NATIVE_UINT64, &id, 1U);
  writeDataset1d(
      gas.get(), "InternalEnergy", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE,
      &internal_energy, 1U);
}

std::string directConfigText(const std::filesystem::path& source) {
  return
      "[mode]\nmode = zoom_in\nic_file = " + source.string() +
      "\nic_convention = gadget_arepo_bridge_v1\n"
      "ic_bridge_source_length_unit_to_si = 3.0856775814913673e22\n"
      "ic_bridge_source_mass_unit_to_si = 1.98847e30\n"
      "ic_bridge_source_velocity_unit_to_si = 1000\n"
      "ic_bridge_coordinate_frame = comoving\n"
      "ic_bridge_velocity_convention = physical_peculiar\n"
      "ic_bridge_length_hubble_exponent = 0\n"
      "ic_bridge_length_scale_factor_exponent = 0\n"
      "ic_bridge_mass_hubble_exponent = 0\n"
      "ic_bridge_mass_scale_factor_exponent = 0\n"
      "ic_bridge_velocity_hubble_exponent = 0\n"
      "ic_bridge_velocity_scale_factor_exponent = 0\n"
      "ic_gas_density_policy = use_config_value\n"
      "ic_gas_density_value_code = 7.5\n";
}

void writeManifestRuntimeConfig(
    const std::filesystem::path& path,
    std::string_view optional_ic_file = {}) {
  std::ofstream output(path);
  output << "[mode]\nmode = zoom_in\n";
  if (!optional_ic_file.empty()) {
    output << "ic_file = " << optional_ic_file << "\n";
  }
  output << "ic_convention = manifest_v1\n"
         << "ic_manifest_file = source.audit.json\n"
         << "\n[output]\nrun_name = manifest_runtime_test\n"
         << "output_directory = outputs\n";
}

}  // namespace

int main() {
  const std::filesystem::path root =
      std::filesystem::current_path() / "manifest_runtime_test_tmp";
  std::filesystem::remove_all(root);
  std::filesystem::create_directories(root);
  const std::filesystem::path source = root / "source.hdf5";
  const std::filesystem::path manifest_path = root / "source.audit.json";
  writeGasSource(source);

  const auto direct_frozen = cosmosim::core::loadFrozenConfigFromString(
      directConfigText(source), "direct_manifest_seed");
  const auto direct = cosmosim::io::readGadgetArepoHdf5Ic(
      source, direct_frozen.config);
  assert(direct.report.manifest.has_value());
  auto manifest = *direct.report.manifest;
  manifest.source_files = {source.filename()};
  cosmosim::io::writeIcManifestJson(manifest, manifest_path);

  const std::filesystem::path runtime_config = root / "runtime.param.txt";
  writeManifestRuntimeConfig(runtime_config);
  const auto frozen = cosmosim::core::loadFrozenConfigFromFile(runtime_config);
  assert(
      frozen.config.mode.ic_gas_density_policy ==
      cosmosim::core::InitialConditionMissingFieldPolicy::kReject);
  assert(frozen.config.mode.ic_file == "generated");

  cosmosim::parallel::MpiContext mpi_context(false, 1, 0);
  cosmosim::core::ProfilerSession profiler(false);
  const cosmosim::workflows::RuntimeServices services{
      .mpi_context = mpi_context,
      .profiler = profiler,
      .deterministic_execution = true};
  const cosmosim::workflows::internal::InitialConditionRuntime runtime(
      frozen, services);
  const std::filesystem::path runtime_run_directory = root / "run";
  std::filesystem::create_directories(runtime_run_directory);
  const auto started = runtime.materialize(
      cosmosim::workflows::ReferenceWorkflowOptions{}, runtime_run_directory);
  assert(started.state.particles.size() == 1U);
  assert(started.state.gas_cells.size() == 1U);
  assert(
      started.import_report.provenance_authority ==
      "supplied_manifest_v1");
  assert(started.state.gas_cells.density_code[0] == 7.5);
  assert(
      started.state.gas_cells.density_code ==
      direct.state.gas_cells.density_code);
  assert(started.import_report.manifest.has_value());
  const auto contract = std::find_if(
      started.import_report.manifest->missing_field_contracts.begin(),
      started.import_report.manifest->missing_field_contracts.end(),
      [](const cosmosim::io::IcMissingFieldContract& value) {
        return value.field_path == "/PartType0/Density";
      });
  assert(contract !=
      started.import_report.manifest->missing_field_contracts.end());
  assert(
      contract->policy ==
      cosmosim::io::IcMissingFieldPolicy::kUseConfigValue);
  assert(contract->configured_value_code == 7.5);

  const std::filesystem::path conflicting_config =
      root / "runtime_conflicting.param.txt";
  writeManifestRuntimeConfig(conflicting_config, "different.hdf5");
  const auto conflicting =
      cosmosim::core::loadFrozenConfigFromFile(conflicting_config);
  const cosmosim::workflows::internal::InitialConditionRuntime
      conflicting_runtime(conflicting, services);
  bool conflict_rejected = false;
  try {
    const std::filesystem::path conflicting_run_directory =
        root / "conflicting_run";
    std::filesystem::create_directories(conflicting_run_directory);
    (void)conflicting_runtime.materialize(
        cosmosim::workflows::ReferenceWorkflowOptions{},
        conflicting_run_directory);
  } catch (const std::invalid_argument& error) {
    conflict_rejected =
        std::string(error.what()).find("conflicting mode.ic_file") !=
        std::string::npos;
  }
  assert(conflict_rejected);

  std::filesystem::remove_all(root);
  return 0;
}
