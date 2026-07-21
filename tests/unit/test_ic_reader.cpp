#include <cassert>
#include <algorithm>
#include <array>
#include <cmath>
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

cosmosim::io::IcManifest makeValidManifest() {
  cosmosim::io::IcManifest manifest;
  manifest.source_files = {"ics.0.hdf5", "ics.1.hdf5"};
  manifest.source_sha256 = {std::string(64, 'a'), std::string(64, 'b')};
  manifest.source_provenance_ids = {
      "sha256:" + manifest.source_sha256[0],
      "sha256:" + manifest.source_sha256[1]};
  manifest.source_file_sizes_bytes = {123U, 456U};
  manifest.original_header_attributes = {"Header=file0", "Header=file1"};
  manifest.num_files_per_snapshot = 2U;
  manifest.num_part_this_file = {
      std::array<std::uint64_t, 6>{2U, 3U, 0U, 0U, 0U, 0U},
      std::array<std::uint64_t, 6>{1U, 4294967300ULL, 0U, 0U, 0U, 0U}};
  manifest.num_part_total = {3U, 4294967303ULL, 0U, 0U, 0U, 0U};
  manifest.num_part_total_high_word = {0U, 1U, 0U, 0U, 0U, 0U};
  manifest.mass_table = {0.0, 2.0, 0.0, 0.0, 0.0, 0.0};
  manifest.box_size = 100.0;
  manifest.scale_factor = 0.5;
  manifest.redshift = 1.0;
  manifest.omega_matter = 0.3;
  manifest.omega_lambda = 0.7;
  manifest.hubble_param = 0.5;
  manifest.fields.push_back(cosmosim::io::IcFieldManifest{
      .source_file_index = 0U,
      .dataset_path = "/PartType1/Coordinates",
      .selected_alias = "Coordinates",
      .scalar_type = "float64",
      .scalar_class = cosmosim::io::IcScalarClass::kFloatingPoint,
      .byte_width = 8U,
      .is_signed = true,
      .byte_order = cosmosim::io::IcByteOrder::kLittleEndian,
      .rank = 2U,
      .dimensions = {7U, 3U},
      .record_count = 7U,
      .base_unit_to_si = 10.0,
      .hubble_exponent = -1.0,
      .scale_factor_exponent = 0.0,
      .coordinate_frame = cosmosim::io::IcCoordinateFrame::kComoving,
      .velocity_convention = cosmosim::io::IcVelocityConvention::kNotVelocity,
      .semantics = cosmosim::io::IcFieldSemantics::kCoordinate,
      .disposition = cosmosim::io::IcFieldDisposition::kConverted,
      .source_unit = "source_length",
      .target_unit = "runtime_length",
      .conversion_equation =
          "target = stored * base_unit_to_si * h^hubble_exponent * "
          "a^scale_factor_exponent / target_si_per_code"});
  return manifest;
}

void testManifestValidationAndConversions() {
  cosmosim::io::IcManifest manifest = makeValidManifest();
  cosmosim::io::validateIcManifest(manifest);
  assert(cosmosim::io::icStoredToSiMultiplier(
             manifest.fields.front(), manifest.hubble_param,
             manifest.scale_factor) == 20.0);
  assert(cosmosim::io::icVelocityConventionMultiplier(
             cosmosim::io::IcVelocityConvention::kPhysicalPeculiar, 0.25) ==
         1.0);
  assert(cosmosim::io::icVelocityConventionMultiplier(
             cosmosim::io::IcVelocityConvention::kSqrtAScaledPeculiar,
             0.25) == 2.0);
  assert(cosmosim::io::icVelocityConventionMultiplier(
             cosmosim::io::IcVelocityConvention::kComovingCoordinateRate,
             0.25) == 0.25);
  const std::string json = cosmosim::io::serializeIcManifestJson(manifest);
  assert(json.find("\"dialect\": \"gadget_arepo_bridge_v1\"") !=
         std::string::npos);
  assert(json.find("\"num_files_per_snapshot\": 2") !=
         std::string::npos);
  const cosmosim::io::IcManifest decoded =
      cosmosim::io::deserializeIcManifestJson(json);
  assert(decoded.source_sha256 == manifest.source_sha256);
  assert(decoded.fields.front().scalar_type == "float64");
  assert(decoded.fields.front().dimensions ==
         manifest.fields.front().dimensions);
  const std::filesystem::path manifest_path =
      std::filesystem::current_path() / "chui_ic_manifest_test.json";
  cosmosim::io::writeIcManifestJson(manifest, manifest_path);
  assert(std::filesystem::exists(manifest_path));
  const cosmosim::io::IcManifest loaded =
      cosmosim::io::readIcManifestJson(manifest_path);
  assert(loaded.source_provenance_ids == manifest.source_provenance_ids);
  std::filesystem::remove(manifest_path);

  manifest.num_part_this_file[0][2] = 1U;
  manifest.num_part_total[2] = 1U;
  bool rejected_implicit_family_mapping = false;
  try {
    cosmosim::io::validateIcManifest(manifest);
  } catch (const std::invalid_argument&) {
    rejected_implicit_family_mapping = true;
  }
  assert(rejected_implicit_family_mapping);
  manifest.species_policy[2] =
      cosmosim::io::IcSpeciesPolicy::kCollisionlessFamily2AsDarkMatter;
  cosmosim::io::validateIcManifest(manifest);
}

void testGeneratedIsolatedIcSpeciesAndOwnership() {
  auto config = cosmosim::core::makeUnvalidatedSimulationConfigForTests();
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
  auto config = cosmosim::core::makeUnvalidatedSimulationConfigForTests();
  const cosmosim::io::IcReadResult result = cosmosim::io::convertGeneratedIsolatedIcToState(config, 4);

  assert(result.state.particles.size() == 20);
  assert(!result.report.defaulted_fields.empty());
}

void testHdf5GateBehavior() {
#if COSMOSIM_ENABLE_HDF5
  auto config = cosmosim::core::makeUnvalidatedSimulationConfigForTests();
  bool threw = false;
  try {
    const auto result = cosmosim::io::readGadgetArepoHdf5Ic("/definitely/missing/file.hdf5", config);
    (void)result;
  } catch (const std::runtime_error&) {
    threw = true;
  }
  assert(threw);
#else
  auto config = cosmosim::core::makeUnvalidatedSimulationConfigForTests();
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
      } else if (type == H5I_DATATYPE) {
        H5Tclose(handle_);
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
  Hdf5Handle space(H5Screate(H5S_SCALAR));
  Hdf5Handle attr(H5Acreate2(header_group, name, H5T_NATIVE_DOUBLE, space.get(), H5P_DEFAULT, H5P_DEFAULT));
  assert(attr.get() >= 0);
  assert(H5Awrite(attr.get(), H5T_NATIVE_DOUBLE, &value) >= 0);
}

void writeHeaderAttributeU32(hid_t header_group, const char* name, std::uint32_t value) {
  Hdf5Handle space(H5Screate(H5S_SCALAR));
  Hdf5Handle attr(H5Acreate2(header_group, name, H5T_NATIVE_UINT32, space.get(), H5P_DEFAULT, H5P_DEFAULT));
  assert(attr.get() >= 0);
  assert(H5Awrite(attr.get(), H5T_NATIVE_UINT32, &value) >= 0);
}

void writeHeaderAttributeString(
    hid_t header_group, const char* name, const std::string& value) {
  Hdf5Handle space(H5Screate(H5S_SCALAR));
  Hdf5Handle type(H5Tcopy(H5T_C_S1));
  assert(type.get() >= 0);
  assert(H5Tset_size(type.get(), value.size() + 1U) >= 0);
  assert(H5Tset_strpad(type.get(), H5T_STR_NULLTERM) >= 0);
  Hdf5Handle attr(H5Acreate2(
      header_group, name, type.get(), space.get(), H5P_DEFAULT, H5P_DEFAULT));
  assert(attr.get() >= 0);
  assert(H5Awrite(attr.get(), type.get(), value.c_str()) >= 0);
}

void writeCanonicalHeaderAttributes(hid_t header_group, bool valid_schema) {
  writeHeaderAttributeString(
      header_group, "ChuiIcSchemaName",
      valid_schema ? "chui_canonical_v1" : "unknown_canonical_schema");
  writeHeaderAttributeU32(header_group, "ChuiIcSchemaVersion", 1U);
  writeHeaderAttributeF64(
      header_group, "ChuiLengthUnitToSI", 3.0856775814913673e19);
  writeHeaderAttributeF64(
      header_group, "ChuiMassUnitToSI", 1.98847e30);
  writeHeaderAttributeF64(
      header_group, "ChuiVelocityUnitToSI", 1.0e3);
  writeHeaderAttributeString(header_group, "ChuiCoordinateFrame", "comoving");
  writeHeaderAttributeString(
      header_group, "ChuiVelocityConvention", "physical_peculiar");
  writeHeaderAttributeString(
      header_group, "ConversionManifestSha256", std::string(64U, 'a'));
}

void writeRequiredHeaderWithTotals(
    hid_t header_group,
    const std::array<std::uint32_t, 6>& counts,
    const std::array<std::uint32_t, 6>& totals,
    std::uint32_t num_files,
    double scale_factor = 1.0,
    double box_size = 50000.0) {
  writeHeaderAttributeU32x6(header_group, "NumPart_ThisFile", counts);
  writeHeaderAttributeU32x6(header_group, "NumPart_Total", totals);
  writeHeaderAttributeU32x6(
      header_group, "NumPart_Total_HighWord", {0, 0, 0, 0, 0, 0});
  writeHeaderAttributeF64x6(
      header_group, "MassTable", {0.0, 0.0, 0.0, 0.0, 0.0, 0.0});
  writeHeaderAttributeF64(header_group, "Time", scale_factor);
  writeHeaderAttributeF64(
      header_group, "Redshift", 1.0 / scale_factor - 1.0);
  writeHeaderAttributeF64(header_group, "BoxSize", box_size);
  writeHeaderAttributeF64(header_group, "Omega0", 0.315);
  writeHeaderAttributeF64(header_group, "OmegaLambda", 0.685);
  writeHeaderAttributeF64(header_group, "HubbleParam", 0.674);
  writeHeaderAttributeU32(header_group, "NumFilesPerSnapshot", num_files);
}

void writeRequiredHeader(
    hid_t header_group,
    const std::array<std::uint32_t, 6>& counts,
    double scale_factor = 1.0) {
  writeRequiredHeaderWithTotals(
      header_group, counts, counts, 1U, scale_factor);
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

void writeDataset1dF32(
    hid_t group, const char* name, const std::vector<float>& values) {
  hsize_t dims[1] = {values.size()};
  Hdf5Handle space(H5Screate_simple(1, dims, nullptr));
  Hdf5Handle dataset(H5Dcreate2(
      group, name, H5T_IEEE_F32LE, space.get(), H5P_DEFAULT, H5P_DEFAULT,
      H5P_DEFAULT));
  assert(dataset.get() >= 0);
  assert(H5Dwrite(
             dataset.get(), H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,
             H5P_DEFAULT, values.data()) >= 0);
}

void writeDataset2dVec3F32(
    hid_t group, const char* name, const std::vector<float>& values) {
  assert(values.size() % 3 == 0);
  hsize_t dims[2] = {values.size() / 3, 3};
  Hdf5Handle space(H5Screate_simple(2, dims, nullptr));
  Hdf5Handle dataset(H5Dcreate2(
      group, name, H5T_IEEE_F32LE, space.get(), H5P_DEFAULT, H5P_DEFAULT,
      H5P_DEFAULT));
  assert(dataset.get() >= 0);
  assert(H5Dwrite(
             dataset.get(), H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,
             H5P_DEFAULT, values.data()) >= 0);
}

std::filesystem::path writeMinimalIcFile(
    bool include_density,
    bool duplicate_ids = false) {
  const std::filesystem::path path =
      std::filesystem::temp_directory_path() /
      (duplicate_ids
           ? "cosmosim_ic_reader_duplicate_ids.hdf5"
           : (include_density ? "cosmosim_ic_reader_gas_present.hdf5" : "cosmosim_ic_reader_gas_missing_density.hdf5"));
  Hdf5Handle file(H5Fcreate(path.string().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT));
  Hdf5Handle header(H5Gcreate2(file.get(), "/Header", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
  writeRequiredHeader(header.get(), {2, 0, 0, 0, 0, 0});

  Hdf5Handle gas(H5Gcreate2(file.get(), "/PartType0", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
  writeDataset2dVec3(gas.get(), "Coordinates", {1.0, 2.0, 3.0, 4.0, 5.0, 6.0});
  writeDataset2dVec3(gas.get(), "Velocities", {10.0, 0.0, 0.0, 20.0, 0.0, 0.0});
  writeDataset1d(gas.get(), "Masses", {5.0, 6.0});
  writeDataset1dIds(
      gas.get(), "ParticleIDs", duplicate_ids
          ? std::vector<std::uint64_t>{101, 101}
          : std::vector<std::uint64_t>{101, 102});
  writeDataset1d(gas.get(), "InternalEnergy", {100.0, 200.0});
  if (include_density) {
    writeDataset1d(gas.get(), "Density", {1.5, 2.5});
  }
  writeDataset1d(gas.get(), "Metallicity", {0.02, 0.03});
  return path;
}

std::filesystem::path writeCanonicalDmIcFile(bool valid_schema) {
  const std::filesystem::path path =
      std::filesystem::temp_directory_path() /
      (valid_schema ? "cosmosim_ic_reader_canonical.hdf5"
                    : "cosmosim_ic_reader_bad_canonical.hdf5");
  Hdf5Handle file(H5Fcreate(
      path.string().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT));
  Hdf5Handle header(H5Gcreate2(
      file.get(), "/Header", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
  writeRequiredHeader(header.get(), {0, 1, 0, 0, 0, 0});
  writeCanonicalHeaderAttributes(header.get(), valid_schema);
  Hdf5Handle dm(H5Gcreate2(
      file.get(), "/PartType1", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
  writeDataset2dVec3(dm.get(), "Coordinates", {10.0, 20.0, 30.0});
  writeDataset2dVec3(dm.get(), "Velocities", {2.0, 3.0, 4.0});
  writeDataset1d(dm.get(), "Masses", {5.0});
  writeDataset1dIds(dm.get(), "ParticleIDs", {301});
  return path;
}

std::filesystem::path writeMinimalBlackHoleIcFile() {
  const std::filesystem::path path =
      std::filesystem::temp_directory_path() /
      "cosmosim_ic_reader_black_hole.hdf5";
  Hdf5Handle file(H5Fcreate(
      path.string().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT));
  Hdf5Handle header(H5Gcreate2(
      file.get(), "/Header", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
  writeRequiredHeader(header.get(), {0, 0, 0, 0, 0, 1});
  Hdf5Handle black_hole(H5Gcreate2(
      file.get(), "/PartType5", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
  writeDataset2dVec3(black_hole.get(), "Coordinates", {10.0, 20.0, 30.0});
  writeDataset2dVec3(black_hole.get(), "Velocities", {2.0, 3.0, 4.0});
  writeDataset1d(black_hole.get(), "Masses", {5.0});
  writeDataset1dIds(black_hole.get(), "ParticleIDs", {501});
  writeDataset1d(black_hole.get(), "BH_Mass", {7.0});
  writeDataset1d(black_hole.get(), "BH_Mdot", {0.25});
  return path;
}

std::filesystem::path writeMinimalFamily2IcFile() {
  const std::filesystem::path path =
      std::filesystem::temp_directory_path() /
      "cosmosim_ic_reader_family2.hdf5";
  Hdf5Handle file(H5Fcreate(
      path.string().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT));
  Hdf5Handle header(H5Gcreate2(
      file.get(), "/Header", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
  writeRequiredHeader(header.get(), {0, 0, 1, 0, 0, 0});
  Hdf5Handle family(H5Gcreate2(
      file.get(), "/PartType2", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
  writeDataset2dVec3(family.get(), "Coordinates", {10.0, 20.0, 30.0});
  writeDataset2dVec3(family.get(), "Velocities", {2.0, 3.0, 4.0});
  writeDataset1d(family.get(), "Masses", {5.0});
  writeDataset1dIds(family.get(), "ParticleIDs", {201});
  return path;
}


std::filesystem::path writeMinimalStarIcFile() {
  const std::filesystem::path path =
      std::filesystem::temp_directory_path() /
      "cosmosim_ic_reader_star.hdf5";
  Hdf5Handle file(H5Fcreate(
      path.string().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT));
  Hdf5Handle header(H5Gcreate2(
      file.get(), "/Header", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
  writeRequiredHeader(header.get(), {0, 0, 0, 0, 1, 0});
  Hdf5Handle star(H5Gcreate2(
      file.get(), "/PartType4", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
  writeDataset2dVec3(star.get(), "Coordinates", {2.0, 3.0, 4.0});
  writeDataset2dVec3(star.get(), "Velocities", {5.0, 6.0, 7.0});
  writeDataset1d(star.get(), "Masses", {8.0});
  writeDataset1dIds(star.get(), "ParticleIDs", {401});
  writeDataset1d(star.get(), "GFM_StellarFormationTime", {0.25});
  writeDataset1d(star.get(), "GFM_InitialMass", {9.0});
  writeDataset1d(star.get(), "GFM_Metallicity", {0.02});
  return path;
}

std::vector<std::filesystem::path> writeMultifileDmSet(
    const std::string& stem,
    bool duplicate_ids = false,
    bool inconsistent_box = false,
    bool inconsistent_schema = false,
    bool omit_second_file = false) {
  std::vector<std::filesystem::path> paths;
  for (std::uint32_t file_index = 0U; file_index < 2U; ++file_index) {
    const auto path = std::filesystem::temp_directory_path() /
        (stem + "." + std::to_string(file_index) + ".hdf5");
    paths.push_back(path);
    if (omit_second_file && file_index == 1U) {
      std::filesystem::remove(path);
      continue;
    }
    Hdf5Handle file(H5Fcreate(
        path.string().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT));
    Hdf5Handle header(H5Gcreate2(
        file.get(), "/Header", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
    writeRequiredHeaderWithTotals(
        header.get(), {0, 1, 0, 0, 0, 0}, {0, 2, 0, 0, 0, 0}, 2U,
        1.0, inconsistent_box && file_index == 1U ? 49000.0 : 50000.0);
    Hdf5Handle dm(H5Gcreate2(
        file.get(), "/PartType1", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
    if (inconsistent_schema && file_index == 1U) {
      writeDataset2dVec3(dm.get(), "Coordinates", {30.0, 1.0, 1.0});
      writeDataset2dVec3(dm.get(), "Velocities", {0.0, 0.0, 0.0});
      writeDataset1d(dm.get(), "Masses", {2.0});
    } else {
      writeDataset2dVec3F32(
          dm.get(), "Coordinates",
          {file_index == 0U ? 10.0F : 30.0F, 1.0F, 1.0F});
      writeDataset2dVec3F32(
          dm.get(), "Velocities", {0.0F, 0.0F, 0.0F});
      writeDataset1dF32(dm.get(), "Masses", {2.0F});
    }
    writeDataset1dIds(
        dm.get(), "ParticleIDs",
        {duplicate_ids ? 601U : 601U + file_index});
  }
  return paths;
}

void removePaths(const std::vector<std::filesystem::path>& paths) {
  for (const auto& path : paths) {
    std::filesystem::remove(path);
  }
}

void testCanonicalHeaderContract() {
  auto config = cosmosim::core::makeUnvalidatedSimulationConfigForTests();
  config.mode.ic_convention =
      cosmosim::core::InitialConditionConvention::kChuiCanonicalV1;
  config.units.length_unit = "kpc";
  const auto canonical_path = writeCanonicalDmIcFile(true);
  const auto result =
      cosmosim::io::readGadgetArepoHdf5Ic(canonical_path, config);
  assert(result.report.manifest.has_value());
  assert(result.report.manifest->dialect ==
         cosmosim::io::IcDialect::kChuiCanonicalV1);
  assert(result.state.particles.size() == 1U);
  assert(result.state.particles.position_x_comoving[0] == 10.0);
  assert(result.state.particles.velocity_x_peculiar[0] == 2.0);
  std::filesystem::remove(canonical_path);

  const auto invalid_path = writeCanonicalDmIcFile(false);
  bool invalid_rejected = false;
  try {
    (void)cosmosim::io::readGadgetArepoHdf5Ic(invalid_path, config);
  } catch (const std::runtime_error& error) {
    invalid_rejected =
        std::string(error.what()).find("unsupported schema") !=
        std::string::npos;
  }
  assert(invalid_rejected);
  std::filesystem::remove(invalid_path);
}

void testHdf5StarSidecarAndMultifileSchema() {
  auto config = cosmosim::core::makeUnvalidatedSimulationConfigForTests();
  const auto star_path = writeMinimalStarIcFile();
  const auto star_result =
      cosmosim::io::readGadgetArepoHdf5Ic(star_path, config);
  assert(star_result.state.star_particles.size() == 1U);
  assert(star_result.state.star_particles.particle_index[0] == 0U);
  assert(star_result.state.star_particles.formation_scale_factor[0] == 0.25);
  assert(star_result.state.star_particles.birth_mass_code[0] == 9.0);
  assert(star_result.state.star_particles.metallicity_mass_fraction[0] == 0.02);
  std::filesystem::remove(star_path);

  const auto paths = writeMultifileDmSet("cosmosim_ic_reader_multifile");
  const auto result = cosmosim::io::readGadgetArepoHdf5Ic(
      paths.front(), config,
      cosmosim::io::IcImportOptions{.chunk_particle_count = 1U});
  assert(result.state.particles.size() == 2U);
  assert(result.report.manifest->num_files_per_snapshot == 2U);
  assert(result.report.counters.files_assigned == 2U);
  assert(result.report.counters.chunks_assigned == 2U);
  bool observed_float32 = false;
  bool observed_scalar_attribute = false;
  for (const auto& field : result.report.manifest->fields) {
    if (field.source_file_index == 0U &&
        field.dataset_path == "/PartType1/Coordinates") {
      observed_float32 = field.scalar_type == "float32" &&
                         field.byte_width == 4U && field.rank == 2U;
    }
    if (field.source_file_index == 0U &&
        field.dataset_path == "/Header/BoxSize") {
      observed_scalar_attribute = field.rank == 0U &&
                                  field.dimensions.empty() &&
                                  field.record_count == 1U;
    }
  }
  assert(observed_float32);
  assert(observed_scalar_attribute);
  removePaths(paths);

  const auto duplicate_paths = writeMultifileDmSet(
      "cosmosim_ic_reader_multifile_duplicate", true);
  bool duplicate_rejected = false;
  try {
    (void)cosmosim::io::readGadgetArepoHdf5Ic(
        duplicate_paths.front(), config,
        cosmosim::io::IcImportOptions{.chunk_particle_count = 1U});
  } catch (const std::runtime_error& error) {
    duplicate_rejected =
        std::string(error.what()).find("duplicate particle IDs") !=
        std::string::npos;
  }
  assert(duplicate_rejected);
  removePaths(duplicate_paths);

  const auto inconsistent_paths = writeMultifileDmSet(
      "cosmosim_ic_reader_multifile_inconsistent", false, true);
  bool header_rejected = false;
  try {
    (void)cosmosim::io::readGadgetArepoHdf5Ic(
        inconsistent_paths.front(), config);
  } catch (const std::runtime_error& error) {
    header_rejected =
        std::string(error.what()).find("inconsistent cosmology") !=
        std::string::npos;
  }
  assert(header_rejected);
  removePaths(inconsistent_paths);

  const auto schema_paths = writeMultifileDmSet(
      "cosmosim_ic_reader_multifile_schema", false, false, true);
  bool schema_rejected = false;
  try {
    (void)cosmosim::io::readGadgetArepoHdf5Ic(
        schema_paths.front(), config);
  } catch (const std::runtime_error& error) {
    schema_rejected =
        std::string(error.what()).find("inconsistent source schema") !=
        std::string::npos;
  }
  assert(schema_rejected);
  removePaths(schema_paths);

  const auto missing_paths = writeMultifileDmSet(
      "cosmosim_ic_reader_multifile_missing", false, false, false, true);
  bool missing_rejected = false;
  try {
    (void)cosmosim::io::readGadgetArepoHdf5Ic(
        missing_paths.front(), config);
  } catch (const std::runtime_error& error) {
    missing_rejected =
        std::string(error.what()).find("missing multifile IC member") !=
        std::string::npos;
  }
  assert(missing_rejected);
  removePaths(missing_paths);
}

void testHdf5GasThermoMapping() {
  auto config = cosmosim::core::makeUnvalidatedSimulationConfigForTests();
  config.output.run_name = "ic_reader_hdf5_gas_mapping";
  const std::filesystem::path path = writeMinimalIcFile(true);

  const cosmosim::io::IcReadResult result = cosmosim::io::readGadgetArepoHdf5Ic(path, config);
  assert(result.state.cells.size() == 2);
  assert(result.state.gas_cells.internal_energy_code[0] == 100.0);
  assert(result.state.gas_cells.internal_energy_code[1] == 200.0);
  assert(std::abs(result.state.gas_cells.density_code[0] - 1.5e9) < 1.0e-6);
  assert(std::abs(result.state.gas_cells.density_code[1] - 2.5e9) < 1.0e-6);
  assert(result.report.manifest.has_value());
  assert(result.report.manifest->dialect ==
         cosmosim::io::IcDialect::kGadgetArepoBridgeV1);
  assert(result.state.cells.center_x_comoving[0] == result.state.particles.position_x_comoving[0]);

  auto mismatched_manifest = *result.report.manifest;
  const auto coordinate_field = std::find_if(
      mismatched_manifest.fields.begin(), mismatched_manifest.fields.end(),
      [](const cosmosim::io::IcFieldManifest& field) {
        return field.dataset_path == "/PartType0/Coordinates";
      });
  assert(coordinate_field != mismatched_manifest.fields.end());
  coordinate_field->byte_width = 4U;
  bool supplied_schema_rejected = false;
  try {
    (void)cosmosim::io::readGadgetArepoHdf5Ic(
        path, config,
        cosmosim::io::IcImportOptions{.manifest = &mismatched_manifest});
  } catch (const std::runtime_error& error) {
    supplied_schema_rejected =
        std::string(error.what()).find("schema does not match") !=
        std::string::npos;
  }
  assert(supplied_schema_rejected);

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
  auto config = cosmosim::core::makeUnvalidatedSimulationConfigForTests();
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

void testHdf5BlackHoleAndValidationFailures() {
  auto config = cosmosim::core::makeUnvalidatedSimulationConfigForTests();
  const std::filesystem::path black_hole_path = writeMinimalBlackHoleIcFile();
  const cosmosim::io::IcReadResult black_hole_result =
      cosmosim::io::readGadgetArepoHdf5Ic(black_hole_path, config);
  assert(black_hole_result.state.particles.size() == 1U);
  assert(black_hole_result.state.particle_sidecar.species_tag[0] ==
         static_cast<std::uint32_t>(
             cosmosim::core::ParticleSpecies::kBlackHole));
  assert(black_hole_result.state.black_holes.size() == 1U);
  assert(black_hole_result.state.black_holes.particle_index[0] == 0U);
  assert(black_hole_result.state.black_holes.subgrid_mass_code[0] == 7.0);
  assert(black_hole_result.state.black_holes.accretion_rate_code[0] > 0.0);
  std::filesystem::remove(black_hole_path);

  const std::filesystem::path family2_path = writeMinimalFamily2IcFile();
  bool implicit_family2_rejected = false;
  try {
    (void)cosmosim::io::readGadgetArepoHdf5Ic(family2_path, config);
  } catch (const std::invalid_argument& error) {
    implicit_family2_rejected =
        std::string(error.what()).find("explicit reject policy") !=
        std::string::npos;
  }
  assert(implicit_family2_rejected);
  cosmosim::io::IcSchemaSummary family2_schema;
  family2_schema.count_by_type = {0, 0, 1, 0, 0, 0};
  family2_schema.total_count_by_type = family2_schema.count_by_type;
  family2_schema.mass_table = {0, 0, 0, 0, 0, 0};
  family2_schema.num_files_per_snapshot = 1U;
  family2_schema.box_size = 50000.0;
  family2_schema.scale_factor = 1.0;
  family2_schema.redshift = 0.0;
  family2_schema.omega_matter = 0.315;
  family2_schema.omega_lambda = 0.685;
  family2_schema.hubble_param = 0.674;
  cosmosim::io::IcManifest family2_manifest =
      cosmosim::io::makeGadgetArepoBridgeV1Manifest(
          family2_path, family2_schema);
  family2_manifest.species_policy[2] =
      cosmosim::io::IcSpeciesPolicy::kCollisionlessFamily2AsDarkMatter;
  const cosmosim::io::IcReadResult family2_result =
      cosmosim::io::readGadgetArepoHdf5Ic(
          family2_path,
          config,
          cosmosim::io::IcImportOptions{.manifest = &family2_manifest});
  assert(family2_result.state.particle_sidecar.species_tag[0] ==
         static_cast<std::uint32_t>(
             cosmosim::core::ParticleSpecies::kDarkMatter));
  std::filesystem::remove(family2_path);

  const std::filesystem::path duplicate_path = writeMinimalIcFile(true, true);
  bool duplicate_rejected = false;
  try {
    (void)cosmosim::io::readGadgetArepoHdf5Ic(duplicate_path, config);
  } catch (const std::runtime_error& error) {
    duplicate_rejected =
        std::string(error.what()).find("duplicate particle IDs") !=
        std::string::npos;
  }
  assert(duplicate_rejected);
  std::filesystem::remove(duplicate_path);

  const std::filesystem::path malformed_path = writeMinimalIcFile(true);
  {
    Hdf5Handle file(H5Fopen(
        malformed_path.string().c_str(), H5F_ACC_RDWR, H5P_DEFAULT));
    assert(H5Ldelete(file.get(), "/PartType0/Coordinates", H5P_DEFAULT) >= 0);
    Hdf5Handle gas(H5Gopen2(file.get(), "/PartType0", H5P_DEFAULT));
    writeDataset1d(gas.get(), "Coordinates", {1, 2, 3, 4, 5, 6});
  }
  bool malformed_rejected = false;
  try {
    (void)cosmosim::io::readGadgetArepoHdf5Ic(malformed_path, config);
  } catch (const std::exception& error) {
    const std::string message = error.what();
    malformed_rejected =
        message.find("dimensions [N,3]") != std::string::npos ||
        message.find("record count disagrees") != std::string::npos;
  }
  assert(malformed_rejected);
  std::filesystem::remove(malformed_path);

  const std::filesystem::path mismatch_path = writeMinimalIcFile(true);
  config.cosmology.omega_matter = 0.4;
  bool cosmology_rejected = false;
  try {
    (void)cosmosim::io::readGadgetArepoHdf5Ic(mismatch_path, config);
  } catch (const std::runtime_error& error) {
    cosmology_rejected =
        std::string(error.what()).find("cosmology/BoxSize/start epoch") !=
        std::string::npos;
  }
  assert(cosmology_rejected);
  std::filesystem::remove(mismatch_path);
}
#endif

}  // namespace

int main() {
  testManifestValidationAndConversions();
  testGeneratedIsolatedIcSpeciesAndOwnership();
  testGeneratedConverterDefaultAudit();
  testHdf5GateBehavior();
#if COSMOSIM_ENABLE_HDF5
  testCanonicalHeaderContract();
  testHdf5StarSidecarAndMultifileSchema();
  testHdf5GasThermoMapping();
  testHdf5GasOptionalDensityMissingBehavior();
  testHdf5BlackHoleAndValidationFailures();
#endif
  return 0;
}
