#include <cassert>
#include <array>
#include <cstdint>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <vector>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/io/ic_reader.hpp"

#if COSMOSIM_ENABLE_HDF5
#include <hdf5.h>
#endif

namespace {

void testGeneratedIsolatedIcSpeciesAndOwnership() {
  cosmosim::core::SimulationConfig config;
  config.output.run_name = "unit_ic_reader";

  const cosmosim::io::IcReadResult result =
      cosmosim::io::buildGeneratedIsolatedIc(config, 5, 3, 1000);

  assert(result.state.particles.size() == 8);
  assert(result.state.species.count_by_species[static_cast<std::size_t>(cosmosim::core::ParticleSpecies::kDarkMatter)] ==
         5);
  assert(result.state.species.count_by_species[static_cast<std::size_t>(cosmosim::core::ParticleSpecies::kGas)] == 3);
  assert(result.state.validateOwnershipInvariants());
}

void testGeneratedConverterDefaultAudit() {
  cosmosim::core::SimulationConfig config;
  const cosmosim::io::IcReadResult result = cosmosim::io::convertGeneratedIsolatedIcToState(config, 4);

  assert(result.state.particles.size() == 20);
  assert(!result.report.defaulted_fields.empty());
}

void testHdf5GateBehavior() {
#if COSMOSIM_ENABLE_HDF5
  cosmosim::core::SimulationConfig config;
  bool threw = false;
  try {
    const auto result = cosmosim::io::readGadgetArepoHdf5Ic("/definitely/missing/file.hdf5", config);
    (void)result;
  } catch (const std::runtime_error&) {
    threw = true;
  }
  assert(threw);
#else
  cosmosim::core::SimulationConfig config;
  bool threw = false;
  try {
    const auto result = cosmosim::io::readGadgetArepoHdf5Ic("ics.hdf5", config);
    (void)result;
  } catch (const std::runtime_error&) {
    threw = true;
  }
  assert(threw);
#endif
}

#if COSMOSIM_ENABLE_HDF5
class Hdf5Handle {
 public:
  explicit Hdf5Handle(hid_t handle = -1) : handle_(handle) {}
  Hdf5Handle(const Hdf5Handle&) = delete;
  Hdf5Handle& operator=(const Hdf5Handle&) = delete;
  Hdf5Handle(Hdf5Handle&& other) noexcept : handle_(other.handle_) { other.handle_ = -1; }
  ~Hdf5Handle() {
    if (handle_ >= 0) {
      const H5I_type_t type = H5Iget_type(handle_);
      if (type == H5I_FILE) {
        H5Fclose(handle_);
      } else if (type == H5I_GROUP) {
        H5Gclose(handle_);
      } else if (type == H5I_DATASPACE) {
        H5Sclose(handle_);
      } else if (type == H5I_DATASET) {
        H5Dclose(handle_);
      } else if (type == H5I_ATTR) {
        H5Aclose(handle_);
      }
    }
  }
  [[nodiscard]] hid_t get() const { return handle_; }

 private:
  hid_t handle_ = -1;
};

void writeHeaderAttributeU32x6(hid_t header_group, const char* name, const std::array<std::uint32_t, 6>& values) {
  hsize_t dims[1] = {6};
  Hdf5Handle space(H5Screate_simple(1, dims, nullptr));
  Hdf5Handle attr(H5Acreate2(header_group, name, H5T_NATIVE_UINT32, space.get(), H5P_DEFAULT, H5P_DEFAULT));
  assert(attr.get() >= 0);
  assert(H5Awrite(attr.get(), H5T_NATIVE_UINT32, values.data()) >= 0);
}

void writeHeaderAttributeF64x6(hid_t header_group, const char* name, const std::array<double, 6>& values) {
  hsize_t dims[1] = {6};
  Hdf5Handle space(H5Screate_simple(1, dims, nullptr));
  Hdf5Handle attr(H5Acreate2(header_group, name, H5T_NATIVE_DOUBLE, space.get(), H5P_DEFAULT, H5P_DEFAULT));
  assert(attr.get() >= 0);
  assert(H5Awrite(attr.get(), H5T_NATIVE_DOUBLE, values.data()) >= 0);
}

void writeHeaderAttributeF64(hid_t header_group, const char* name, double value) {
  hsize_t dims[1] = {1};
  Hdf5Handle space(H5Screate_simple(1, dims, nullptr));
  Hdf5Handle attr(H5Acreate2(header_group, name, H5T_NATIVE_DOUBLE, space.get(), H5P_DEFAULT, H5P_DEFAULT));
  assert(attr.get() >= 0);
  assert(H5Awrite(attr.get(), H5T_NATIVE_DOUBLE, &value) >= 0);
}

void writeDataset1d(hid_t group, const char* name, const std::vector<double>& values) {
  hsize_t dims[1] = {values.size()};
  Hdf5Handle space(H5Screate_simple(1, dims, nullptr));
  Hdf5Handle dataset(H5Dcreate2(group, name, H5T_IEEE_F64LE, space.get(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
  assert(dataset.get() >= 0);
  assert(H5Dwrite(dataset.get(), H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, values.data()) >= 0);
}

void writeDataset1dIds(hid_t group, const char* name, const std::vector<std::uint64_t>& values) {
  hsize_t dims[1] = {values.size()};
  Hdf5Handle space(H5Screate_simple(1, dims, nullptr));
  Hdf5Handle dataset(H5Dcreate2(group, name, H5T_STD_U64LE, space.get(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
  assert(dataset.get() >= 0);
  assert(H5Dwrite(dataset.get(), H5T_NATIVE_UINT64, H5S_ALL, H5S_ALL, H5P_DEFAULT, values.data()) >= 0);
}

void writeDataset2dVec3(hid_t group, const char* name, const std::vector<double>& values) {
  assert(values.size() % 3 == 0);
  hsize_t dims[2] = {values.size() / 3, 3};
  Hdf5Handle space(H5Screate_simple(2, dims, nullptr));
  Hdf5Handle dataset(H5Dcreate2(group, name, H5T_IEEE_F64LE, space.get(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
  assert(dataset.get() >= 0);
  assert(H5Dwrite(dataset.get(), H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, values.data()) >= 0);
}

std::filesystem::path writeMinimalIcFile(bool include_density) {
  const std::filesystem::path path =
      std::filesystem::temp_directory_path() /
      (include_density ? "cosmosim_ic_reader_gas_present.hdf5" : "cosmosim_ic_reader_gas_missing_density.hdf5");
  Hdf5Handle file(H5Fcreate(path.string().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT));
  Hdf5Handle header(H5Gcreate2(file.get(), "/Header", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
  writeHeaderAttributeU32x6(header.get(), "NumPart_ThisFile", {2, 0, 0, 0, 0, 0});
  writeHeaderAttributeF64x6(header.get(), "MassTable", {0.0, 0.0, 0.0, 0.0, 0.0, 0.0});
  writeHeaderAttributeF64(header.get(), "Time", 1.0);

  Hdf5Handle gas(H5Gcreate2(file.get(), "/PartType0", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
  writeDataset2dVec3(gas.get(), "Coordinates", {1.0, 2.0, 3.0, 4.0, 5.0, 6.0});
  writeDataset2dVec3(gas.get(), "Velocities", {10.0, 0.0, 0.0, 20.0, 0.0, 0.0});
  writeDataset1d(gas.get(), "Masses", {5.0, 6.0});
  writeDataset1dIds(gas.get(), "ParticleIDs", {101, 102});
  writeDataset1d(gas.get(), "InternalEnergy", {100.0, 200.0});
  if (include_density) {
    writeDataset1d(gas.get(), "Density", {1.5, 2.5});
  }
  writeDataset1d(gas.get(), "Metallicity", {0.02, 0.03});
  return path;
}

void testHdf5GasThermoMapping() {
  cosmosim::core::SimulationConfig config;
  config.output.run_name = "ic_reader_hdf5_gas_mapping";
  const std::filesystem::path path = writeMinimalIcFile(true);

  const cosmosim::io::IcReadResult result = cosmosim::io::readGadgetArepoHdf5Ic(path, config);
  assert(result.state.cells.size() == 2);
  assert(result.state.gas_cells.internal_energy_code[0] == 100.0);
  assert(result.state.gas_cells.internal_energy_code[1] == 200.0);
  assert(result.state.gas_cells.density_code[0] == 1.5);
  assert(result.state.gas_cells.density_code[1] == 2.5);
  assert(result.state.cells.center_x_comoving[0] == result.state.particles.position_x_comoving[0]);

  bool found_thermo_bypass = false;
  bool found_metallicity_unsupported = false;
  for (const std::string& value : result.report.unsupported_fields) {
    if (value.find("thermodynamic fields currently bypassed") != std::string::npos) {
      found_thermo_bypass = true;
    }
    if (value.find("PartType0/Metallicity") != std::string::npos) {
      found_metallicity_unsupported = true;
    }
  }
  assert(!found_thermo_bypass);
  assert(found_metallicity_unsupported);
  std::filesystem::remove(path);
}

void testHdf5GasOptionalDensityMissingBehavior() {
  cosmosim::core::SimulationConfig config;
  config.output.run_name = "ic_reader_hdf5_optional";
  const std::filesystem::path path = writeMinimalIcFile(false);

  const cosmosim::io::IcReadResult result = cosmosim::io::readGadgetArepoHdf5Ic(path, config);
  assert(result.state.gas_cells.density_code[0] == 0.0);
  assert(result.state.gas_cells.density_code[1] == 0.0);

  bool recorded_missing_density = false;
  bool recorded_defaulted_density = false;
  for (const std::string& value : result.report.missing_optional_fields) {
    if (value == "/PartType0/Density") {
      recorded_missing_density = true;
      break;
    }
  }
  for (const std::string& value : result.report.defaulted_fields) {
    if (value == "/PartType0/Density=zero") {
      recorded_defaulted_density = true;
      break;
    }
  }
  assert(recorded_missing_density);
  assert(recorded_defaulted_density);
  std::filesystem::remove(path);
}
#endif

}  // namespace

int main() {
  testGeneratedIsolatedIcSpeciesAndOwnership();
  testGeneratedConverterDefaultAudit();
  testHdf5GateBehavior();
#if COSMOSIM_ENABLE_HDF5
  testHdf5GasThermoMapping();
  testHdf5GasOptionalDensityMissingBehavior();
#endif
  return 0;
}
