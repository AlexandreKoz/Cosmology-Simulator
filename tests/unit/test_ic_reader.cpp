#include <cassert>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/io/ic_reader.hpp"
#include "io/internal/ic_canonical_limits.hpp"
#include "io/internal/ic_reader_session.hpp"

#if COSMOSIM_ENABLE_HDF5
#include <hdf5.h>
#endif

namespace {


void testCanonicalSingleFileCountLimit() {
  std::array<std::uint64_t, 6> counts{};
  counts[0] = static_cast<std::uint64_t>(UINT32_MAX);
  cosmosim::io::internal::validateCanonicalSingleFileCounts(counts);

  counts[0] = static_cast<std::uint64_t>(UINT32_MAX) + 1ULL;
  bool rejected = false;
  try {
    cosmosim::io::internal::validateCanonicalSingleFileCounts(counts);
  } catch (const std::length_error& error) {
    rejected = std::string(error.what()).find("PartType0") !=
        std::string::npos;
  }
  assert(rejected);

  counts = {};
  counts[1] = static_cast<std::uint64_t>(UINT32_MAX);
  counts[4] = static_cast<std::uint64_t>(UINT32_MAX);
  cosmosim::io::internal::validateCanonicalSingleFileCounts(counts);
  counts[4] += 1ULL;
  rejected = false;
  try {
    cosmosim::io::internal::validateCanonicalSingleFileCounts(counts);
  } catch (const std::length_error& error) {
    rejected = std::string(error.what()).find("PartType4") !=
        std::string::npos;
  }
  assert(rejected);
}

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
      .length_power = 1,
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

  manifest.fields.front().length_power = 0;
  bool rejected_dimension_drift = false;
  try {
    cosmosim::io::validateIcManifest(manifest);
  } catch (const std::invalid_argument&) {
    rejected_dimension_drift = true;
  }
  assert(rejected_dimension_drift);
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
  Hdf5Handle& operator=(Hdf5Handle&& other) noexcept {
    if (this != &other) {
      this->~Hdf5Handle();
      handle_ = other.handle_;
      other.handle_ = -1;
    }
    return *this;
  }
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

void writeHeaderAttributeI32x6(
    hid_t header_group, const char* name,
    const std::array<std::int32_t, 6>& values) {
  hsize_t dims[1] = {6};
  Hdf5Handle space(H5Screate_simple(1, dims, nullptr));
  Hdf5Handle attr(H5Acreate2(
      header_group, name, H5T_STD_I32LE, space.get(), H5P_DEFAULT,
      H5P_DEFAULT));
  assert(attr.get() >= 0);
  assert(H5Awrite(attr.get(), H5T_NATIVE_INT32, values.data()) >= 0);
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

void writeHeaderAttributeI32(
    hid_t header_group, const char* name, std::int32_t value) {
  Hdf5Handle space(H5Screate(H5S_SCALAR));
  Hdf5Handle attr(H5Acreate2(
      header_group, name, H5T_STD_I32LE, space.get(), H5P_DEFAULT,
      H5P_DEFAULT));
  assert(attr.get() >= 0);
  assert(H5Awrite(attr.get(), H5T_NATIVE_INT32, &value) >= 0);
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

void writeCanonicalHeaderAttributes(
    hid_t header_group, bool valid_schema, std::string_view digest,
    std::string_view sidecar_name, std::string_view marker_name) {
  writeHeaderAttributeString(
      header_group, "ChuiIcSchemaName",
      valid_schema ? "chui_canonical_v1" : "unknown_canonical_schema");
  writeHeaderAttributeU32(header_group, "ChuiIcSchemaVersion", 2U);
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
      header_group, "ConversionManifestSha256", std::string(digest));
  writeHeaderAttributeString(
      header_group, "ConversionManifestSidecar", std::string(sidecar_name));
  writeHeaderAttributeString(
      header_group, "ConversionBundleMarker", std::string(marker_name));
}

void writeRequiredHeaderWithTotals(
    hid_t header_group,
    const std::array<std::uint32_t, 6>& counts,
    const std::array<std::uint32_t, 6>& totals,
    std::uint32_t num_files,
    double scale_factor = 1.0,
    double box_size = 50000.0,
    bool signed_arepo_counts = true) {
  if (signed_arepo_counts) {
    std::array<std::int32_t, 6> signed_counts{};
    for (std::size_t i = 0; i < counts.size(); ++i) {
      assert(counts[i] <= static_cast<std::uint32_t>(INT32_MAX));
      signed_counts[i] = static_cast<std::int32_t>(counts[i]);
    }
    writeHeaderAttributeI32x6(
        header_group, "NumPart_ThisFile", signed_counts);
  } else {
    writeHeaderAttributeU32x6(header_group, "NumPart_ThisFile", counts);
  }
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
  if (signed_arepo_counts) {
    assert(num_files <= static_cast<std::uint32_t>(INT32_MAX));
    writeHeaderAttributeI32(
        header_group, "NumFilesPerSnapshot",
        static_cast<std::int32_t>(num_files));
  } else {
    writeHeaderAttributeU32(header_group, "NumFilesPerSnapshot", num_files);
  }
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

void writeDatasetBytes(
    hid_t group, const char* name, const std::string& values) {
  hsize_t dims[1] = {values.size()};
  Hdf5Handle space(H5Screate_simple(1, dims, nullptr));
  Hdf5Handle dataset(H5Dcreate2(
      group, name, H5T_STD_U8LE, space.get(), H5P_DEFAULT, H5P_DEFAULT,
      H5P_DEFAULT));
  assert(dataset.get() >= 0);
  assert(H5Dwrite(
             dataset.get(), H5T_NATIVE_UINT8, H5S_ALL, H5S_ALL,
             H5P_DEFAULT, values.data()) >= 0);
}

void writeTextFile(
    const std::filesystem::path& path, const std::string& value) {
  std::ofstream output(path, std::ios::binary | std::ios::trunc);
  assert(output);
  output.write(value.data(), static_cast<std::streamsize>(value.size()));
  output.flush();
  assert(output);
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
    bool duplicate_ids = false,
    bool signed_arepo_counts = true) {
  const std::filesystem::path path =
      std::filesystem::temp_directory_path() /
      (duplicate_ids
           ? "cosmosim_ic_reader_duplicate_ids.hdf5"
           : (include_density
                  ? (signed_arepo_counts
                         ? "cosmosim_ic_reader_gas_present_signed.hdf5"
                         : "cosmosim_ic_reader_gas_present_unsigned.hdf5")
                  : "cosmosim_ic_reader_gas_missing_density.hdf5"));
  Hdf5Handle file(H5Fcreate(path.string().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT));
  Hdf5Handle header(H5Gcreate2(file.get(), "/Header", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
  writeRequiredHeaderWithTotals(
      header.get(), {2, 0, 0, 0, 0, 0}, {2, 0, 0, 0, 0, 0}, 1U,
      1.0, 50000.0, signed_arepo_counts);

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


std::filesystem::path writeDimensionalGasBlackHoleIcFile() {
  const auto path = std::filesystem::temp_directory_path() /
      "cosmosim_ic_dimensional_contract.hdf5";
  Hdf5Handle file(H5Fcreate(
      path.string().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT));
  Hdf5Handle header(H5Gcreate2(
      file.get(), "/Header", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
  writeRequiredHeaderWithTotals(
      header.get(), {1, 0, 0, 0, 0, 1}, {1, 0, 0, 0, 0, 1}, 1U,
      0.5, 12.5);
  {
    Hdf5Handle hubble(H5Aopen(header.get(), "HubbleParam", H5P_DEFAULT));
    const double value = 0.5;
    assert(hubble.get() >= 0);
    assert(H5Awrite(hubble.get(), H5T_NATIVE_DOUBLE, &value) >= 0);
  }

  Hdf5Handle gas(H5Gcreate2(
      file.get(), "/PartType0", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
  writeDataset2dVec3(gas.get(), "Coordinates", {1.0, 2.0, 3.0});
  writeDataset2dVec3(gas.get(), "Velocities", {1.0, 0.0, 0.0});
  writeDataset1d(gas.get(), "Masses", {1.0});
  writeDataset1dIds(gas.get(), "ParticleIDs", {1001U});
  writeDataset1d(gas.get(), "InternalEnergy", {4.0});
  writeDataset1d(gas.get(), "Density", {2.0});

  Hdf5Handle black_hole(H5Gcreate2(
      file.get(), "/PartType5", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
  writeDataset2dVec3(black_hole.get(), "Coordinates", {2.0, 2.0, 2.0});
  writeDataset2dVec3(black_hole.get(), "Velocities", {1.0, 0.0, 0.0});
  writeDataset1d(black_hole.get(), "Masses", {1.0});
  writeDataset1dIds(black_hole.get(), "ParticleIDs", {5001U});
  writeDataset1d(black_hole.get(), "BH_Mass", {3.0});
  writeDataset1d(black_hole.get(), "BH_Mdot", {0.25});
  return path;
}

std::filesystem::path canonicalSidecarPath(
    const std::filesystem::path& canonical_path) {
  return canonical_path.string() + ".manifest.json";
}

std::filesystem::path canonicalMarkerPath(
    const std::filesystem::path& canonical_path) {
  return canonical_path.string() + ".complete";
}

void removeCanonicalBundle(const std::filesystem::path& canonical_path) {
  std::filesystem::remove(canonical_path);
  std::filesystem::remove(canonicalSidecarPath(canonical_path));
  std::filesystem::remove(canonicalMarkerPath(canonical_path));
}

std::filesystem::path writeCanonicalDmIcFile(bool valid_schema) {
  const std::filesystem::path path =
      std::filesystem::temp_directory_path() /
      (valid_schema ? "cosmosim_ic_reader_canonical.hdf5"
                    : "cosmosim_ic_reader_bad_canonical.hdf5");
  removeCanonicalBundle(path);
  const std::string manifest_json =
      cosmosim::io::serializeIcManifestJson(makeValidManifest());
  const std::string digest = cosmosim::io::icSha256Hex(manifest_json);
  const auto sidecar_path = canonicalSidecarPath(path);
  const auto marker_path = canonicalMarkerPath(path);

  Hdf5Handle file(H5Fcreate(
      path.string().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT));
  Hdf5Handle header(H5Gcreate2(
      file.get(), "/Header", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
  writeRequiredHeader(header.get(), {0, 1, 0, 0, 0, 0});
  writeCanonicalHeaderAttributes(
      header.get(), valid_schema, digest, sidecar_path.filename().string(),
      marker_path.filename().string());
  Hdf5Handle provenance(H5Gcreate2(
      file.get(), "/Provenance", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
  writeDatasetBytes(provenance.get(), "ConversionManifestJson", manifest_json);
  Hdf5Handle dm(H5Gcreate2(
      file.get(), "/PartType1", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
  writeDataset2dVec3(dm.get(), "Coordinates", {10.0, 20.0, 30.0});
  writeDataset2dVec3(dm.get(), "Velocities", {2.0, 3.0, 4.0});
  writeDataset1d(dm.get(), "Masses", {5.0});
  writeDataset1dIds(dm.get(), "ParticleIDs", {301});
  file = Hdf5Handle{};
  writeTextFile(sidecar_path, manifest_json);
  writeTextFile(
      marker_path,
      "chui_ic_bundle_v1\nsha256=" + digest +
          "\ncanonical=" + path.filename().string() +
          "\nmanifest=" + sidecar_path.filename().string() + "\n");
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
        1.0, inconsistent_box && file_index == 1U ? 49000.0 : 50000.0,
        file_index == 0U);
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


cosmosim::core::SimulationConfig makeExplicitBridgeConfig() {
  auto config = cosmosim::core::makeUnvalidatedSimulationConfigForTests();
  config.mode.ic_convention =
      cosmosim::core::InitialConditionConvention::kGadgetArepoBridgeV1;
  config.mode.ic_bridge_source_length_unit_to_si = 3.0856775814913673e19;
  config.mode.ic_bridge_source_mass_unit_to_si = 1.98847e30;
  config.mode.ic_bridge_source_velocity_unit_to_si = 1.0e3;
  config.mode.ic_bridge_coordinate_frame =
      cosmosim::core::InitialConditionCoordinateFrame::kComoving;
  config.mode.ic_bridge_velocity_convention =
      cosmosim::core::InitialConditionVelocityConvention::kPhysicalPeculiar;
  config.mode.ic_bridge_length_hubble_exponent = 0.0;
  config.mode.ic_bridge_length_scale_factor_exponent = 0.0;
  config.mode.ic_bridge_mass_hubble_exponent = 0.0;
  config.mode.ic_bridge_mass_scale_factor_exponent = 0.0;
  config.mode.ic_bridge_velocity_hubble_exponent = 0.0;
  config.mode.ic_bridge_velocity_scale_factor_exponent = 0.0;
  return config;
}

void removeDataset(
    const std::filesystem::path& path,
    const char* group_name,
    const char* dataset_name) {
  Hdf5Handle file(H5Fopen(
      path.string().c_str(), H5F_ACC_RDWR, H5P_DEFAULT));
  Hdf5Handle group(H5Gopen2(file.get(), group_name, H5P_DEFAULT));
  assert(file.get() >= 0 && group.get() >= 0);
  assert(H5Ldelete(group.get(), dataset_name, H5P_DEFAULT) >= 0);
}

void replaceDataset1d(
    const std::filesystem::path& path,
    const char* group_name,
    const char* dataset_name,
    const std::vector<double>& values) {
  Hdf5Handle file(H5Fopen(
      path.string().c_str(), H5F_ACC_RDWR, H5P_DEFAULT));
  Hdf5Handle group(H5Gopen2(file.get(), group_name, H5P_DEFAULT));
  assert(file.get() >= 0 && group.get() >= 0);
  assert(H5Ldelete(group.get(), dataset_name, H5P_DEFAULT) >= 0);
  writeDataset1d(group.get(), dataset_name, values);
}

void expectIcReadFailure(
    const std::filesystem::path& path,
    const cosmosim::core::SimulationConfig& config,
    std::string_view expected_text);

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
  assert(std::abs(result.state.particles.velocity_x_peculiar[0] - 2.0) < 1.0e-12);
  assert(result.report.manifest_verified);
  assert(result.report.verified_manifest_sha256.size() == 64U);
  removeCanonicalBundle(canonical_path);

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
  removeCanonicalBundle(invalid_path);

  const auto tampered_path = writeCanonicalDmIcFile(true);
  writeTextFile(canonicalSidecarPath(tampered_path), "{}\n");
  bool tamper_rejected = false;
  try {
    (void)cosmosim::io::readGadgetArepoHdf5Ic(tampered_path, config);
  } catch (const std::runtime_error& error) {
    tamper_rejected =
        std::string(error.what()).find("sidecar manifest SHA-256") !=
        std::string::npos;
  }
  assert(tamper_rejected);
  removeCanonicalBundle(tampered_path);
}

void testHdf5StarSidecarAndMultifileSchema() {
  auto config = makeExplicitBridgeConfig();
  const auto star_path = writeMinimalStarIcFile();
  const auto star_result =
      cosmosim::io::readGadgetArepoHdf5Ic(star_path, config);
  assert(star_result.state.star_particles.size() == 1U);
  assert(star_result.state.star_particles.particle_index[0] == 0U);
  assert(star_result.state.star_particles.formation_scale_factor[0] == 0.25);
  assert(star_result.state.star_particles.birth_mass_code[0] == 9.0);
  assert(star_result.state.star_particles.metallicity_mass_fraction[0] == 0.02);
  std::filesystem::remove(star_path);

  auto missing_star_path = writeMinimalStarIcFile();
  removeDataset(
      missing_star_path, "/PartType4", "GFM_StellarFormationTime");
  expectIcReadFailure(
      missing_star_path, config,
      "/PartType4/StellarFormationTime is missing and its normalized missing-field policy is reject");
  missing_star_path = writeMinimalStarIcFile();
  removeDataset(
      missing_star_path, "/PartType4", "GFM_StellarFormationTime");
  config.mode.ic_star_formation_time_policy =
      cosmosim::core::InitialConditionMissingFieldPolicy::kDialectDefinedDefault;
  const auto defaulted_star =
      cosmosim::io::readGadgetArepoHdf5Ic(missing_star_path, config);
  assert(defaulted_star.state.star_particles.formation_scale_factor[0] == 1.0);
  assert(!defaulted_star.report.defaulted_fields.empty());
  std::filesystem::remove(missing_star_path);
  config.mode.ic_star_formation_time_policy =
      cosmosim::core::InitialConditionMissingFieldPolicy::kReject;

  const auto unsigned_path = writeMinimalIcFile(true, false, false);
  const auto unsigned_result = cosmosim::io::readGadgetArepoHdf5Ic(
      unsigned_path, config,
      cosmosim::io::IcImportOptions{.chunk_particle_count = 1U});
  assert(unsigned_result.state.particles.size() == 2U);
  assert(unsigned_result.report.counters.chunks_assigned == 2U);
  assert(unsigned_result.report.counters.source_file_open_count == 1U);
  assert(unsigned_result.report.counters.source_dataset_open_count > 0U);
  assert(unsigned_result.report.counters.source_dataset_open_count < 14U);
  const auto unsigned_header = std::find_if(
      unsigned_result.report.manifest->fields.begin(),
      unsigned_result.report.manifest->fields.end(),
      [](const cosmosim::io::IcFieldManifest& field) {
        return field.dataset_path == "/Header/NumPart_ThisFile";
      });
  assert(unsigned_header != unsigned_result.report.manifest->fields.end());
  assert(!unsigned_header->is_signed && unsigned_header->byte_width == 4U);
  std::filesystem::remove(unsigned_path);

  const auto paths = writeMultifileDmSet("cosmosim_ic_reader_multifile");
  const auto result = cosmosim::io::readGadgetArepoHdf5Ic(
      paths.front(), config,
      cosmosim::io::IcImportOptions{.chunk_particle_count = 1U});
  assert(result.state.particles.size() == 2U);
  assert(result.report.manifest->num_files_per_snapshot == 2U);
  assert(result.report.counters.files_assigned == 2U);
  assert(result.report.counters.chunks_assigned == 2U);
  bool observed_signed_count = false;
  bool observed_unsigned_count = false;
  for (const auto& field : result.report.manifest->fields) {
    if (field.dataset_path == "/Header/NumPart_ThisFile") {
      observed_signed_count = observed_signed_count || field.is_signed;
      observed_unsigned_count = observed_unsigned_count || !field.is_signed;
    }
  }
  assert(observed_signed_count && observed_unsigned_count);
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
  auto config = makeExplicitBridgeConfig();
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
  assert(result.report.counters.source_file_open_count == 1U);
  assert(result.report.counters.full_file_hash_pass_count == 3U);
  assert(result.report.counters.source_identity_validation_count == 2U);
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
  } catch (const std::exception& error) {
    const std::string message = error.what();
    supplied_schema_rejected =
        message.find("schema does not match") != std::string::npos ||
        message.find("datatype/rank/dimensions/count are inconsistent") !=
            std::string::npos;
  }
  assert(supplied_schema_rejected);

  auto hash_mismatched_manifest = *result.report.manifest;
  hash_mismatched_manifest.source_sha256.front() = std::string(64U, '0');
  hash_mismatched_manifest.source_provenance_ids.front() =
      "sha256:" + hash_mismatched_manifest.source_sha256.front();
  bool supplied_hash_rejected = false;
  try {
    (void)cosmosim::io::readGadgetArepoHdf5Ic(
        path, config,
        cosmosim::io::IcImportOptions{.manifest = &hash_mismatched_manifest});
  } catch (const std::runtime_error& error) {
    supplied_hash_rejected =
        std::string(error.what()).find("SHA-256") != std::string::npos ||
        std::string(error.what()).find("provenance") != std::string::npos;
  }
  assert(supplied_hash_rejected);

  const auto expect_invalid_supplied_manifest = [&](
      const cosmosim::io::IcManifest& manifest,
      std::string_view expected_text) {
    bool rejected = false;
    try {
      (void)cosmosim::io::readGadgetArepoHdf5Ic(
          path, config,
          cosmosim::io::IcImportOptions{.manifest = &manifest});
    } catch (const std::exception& error) {
      rejected = std::string(error.what()).find(expected_text) !=
          std::string::npos;
    }
    assert(rejected);
  };

  auto empty_hash_manifest = *result.report.manifest;
  empty_hash_manifest.source_sha256.clear();
  expect_invalid_supplied_manifest(
      empty_hash_manifest, "one count/hash/size/header record");

  auto incomplete_source_manifest = *result.report.manifest;
  incomplete_source_manifest.original_header_attributes.clear();
  expect_invalid_supplied_manifest(
      incomplete_source_manifest, "one count/hash/size/header record");

  auto invalid_dialect_manifest = *result.report.manifest;
  invalid_dialect_manifest.dialect_version = "99";
  expect_invalid_supplied_manifest(
      invalid_dialect_manifest, "unsupported IC dialect version");

  auto mismatched_vector_manifest = *result.report.manifest;
  mismatched_vector_manifest.source_file_sizes_bytes.push_back(1U);
  expect_invalid_supplied_manifest(
      mismatched_vector_manifest, "one count/hash/size/header record");

  assert(std::abs(result.state.gas_cells.metal_mass_code[0] - 0.10) <
         1.0e-12);
  assert(std::abs(result.state.gas_cells.metal_mass_code[1] - 0.18) <
         1.0e-12);
  const auto metallicity_field = std::find_if(
      result.report.manifest->fields.begin(),
      result.report.manifest->fields.end(),
      [](const cosmosim::io::IcFieldManifest& field) {
        return field.dataset_path == "/PartType0/Metallicity";
      });
  assert(metallicity_field != result.report.manifest->fields.end());
  assert(metallicity_field->disposition ==
         cosmosim::io::IcFieldDisposition::kConverted);
  assert(metallicity_field->semantics ==
         cosmosim::io::IcFieldSemantics::kIntensive);
  assert(metallicity_field->conversion_equation.find("target = stored") !=
         std::string::npos);
  assert(std::find(
             result.report.manifest->converted_fields.begin(),
             result.report.manifest->converted_fields.end(),
             "/PartType0/Metallicity") !=
         result.report.manifest->converted_fields.end());
  assert(std::find(
             result.report.manifest->conversion_equations.begin(),
             result.report.manifest->conversion_equations.end(),
             metallicity_field->conversion_equation) !=
         result.report.manifest->conversion_equations.end());
  for (const std::string& value : result.report.unsupported_fields) {
    assert(value.find("thermodynamic fields currently bypassed") ==
           std::string::npos);
    assert(value.find("PartType0/Metallicity") == std::string::npos);
  }
  std::filesystem::remove(path);

  for (const std::vector<double>& invalid_values : {
           std::vector<double>{-0.01, 0.03},
           std::vector<double>{1.01, 0.03},
           std::vector<double>{
               std::numeric_limits<double>::quiet_NaN(), 0.03}}) {
    const auto invalid_path = writeMinimalIcFile(true);
    replaceDataset1d(
        invalid_path, "/PartType0", "Metallicity", invalid_values);
    expectIcReadFailure(
        invalid_path, config,
        "density, internal energy, and metallicity must be finite and physical");
  }
}

void expectIcReadFailure(
    const std::filesystem::path& path,
    const cosmosim::core::SimulationConfig& config,
    std::string_view expected_text);

void testHdf5GasOptionalDensityMissingBehavior() {
  auto config = makeExplicitBridgeConfig();
  config.output.run_name = "ic_reader_hdf5_optional";

  auto path = writeMinimalIcFile(true);
  removeDataset(path, "/PartType0", "Metallicity");
  const auto zero_metallicity =
      cosmosim::io::readGadgetArepoHdf5Ic(path, config);
  assert(zero_metallicity.state.gas_cells.metal_mass_code[0] == 0.0);
  assert(zero_metallicity.state.gas_cells.metal_mass_code[1] == 0.0);
  assert(zero_metallicity.report.manifest.has_value());
  const auto metallicity_contract = std::find_if(
      zero_metallicity.report.manifest->missing_field_contracts.begin(),
      zero_metallicity.report.manifest->missing_field_contracts.end(),
      [](const cosmosim::io::IcMissingFieldContract& value) {
        return value.field_path == "/PartType0/Metallicity";
      });
  assert(metallicity_contract !=
         zero_metallicity.report.manifest->missing_field_contracts.end());
  assert(metallicity_contract->policy ==
         cosmosim::io::IcMissingFieldPolicy::kDialectDefinedDefault);
  assert(metallicity_contract->resolution.find("zero gas metallicity") !=
         std::string::npos);
  std::filesystem::remove(path);

  path = writeMinimalIcFile(true);
  removeDataset(path, "/PartType0", "InternalEnergy");
  expectIcReadFailure(
      path, config,
      "/PartType0/InternalEnergy is missing and its normalized missing-field policy is reject");

  path = writeMinimalIcFile(false);
  expectIcReadFailure(
      path, config,
      "/PartType0/Density is missing and its normalized missing-field policy is reject");

  path = writeMinimalIcFile(false);
  config.mode.ic_gas_density_policy =
      cosmosim::core::InitialConditionMissingFieldPolicy::kUseConfigValue;
  config.mode.ic_gas_density_value_code = 3.25;
  const auto result = cosmosim::io::readGadgetArepoHdf5Ic(path, config);
  assert(result.state.gas_cells.density_code[0] == 3.25);
  assert(result.state.gas_cells.density_code[1] == 3.25);
  assert(result.report.manifest.has_value());
  const auto contract = std::find_if(
      result.report.manifest->missing_field_contracts.begin(),
      result.report.manifest->missing_field_contracts.end(),
      [](const cosmosim::io::IcMissingFieldContract& value) {
        return value.field_path == "/PartType0/Density";
      });
  assert(contract != result.report.manifest->missing_field_contracts.end());
  assert(contract->policy == cosmosim::io::IcMissingFieldPolicy::kUseConfigValue);
  assert(contract->configured_value_code == 3.25);

  auto manifest_authority_config = config;
  manifest_authority_config.mode.ic_gas_density_policy =
      cosmosim::core::InitialConditionMissingFieldPolicy::kReject;
  manifest_authority_config.mode.ic_gas_density_value_code = 0.0;
  const std::filesystem::path ignored_requested_path =
      path.parent_path() / "manifest_authority_ignored_source.hdf5";
  const auto manifest_authoritative =
      cosmosim::io::readGadgetArepoHdf5Ic(
          ignored_requested_path, manifest_authority_config,
          cosmosim::io::IcImportOptions{
              .manifest = &*result.report.manifest});
  assert(manifest_authoritative.state.gas_cells.density_code[0] == 3.25);
  assert(manifest_authoritative.state.gas_cells.density_code[1] == 3.25);
  assert(
      manifest_authoritative.report.provenance_authority ==
      "supplied_manifest_v1");
  std::filesystem::remove(path);
}


void testSharedDimensionalConversionContract() {
  auto config = makeExplicitBridgeConfig();
  config.mode.ic_bridge_coordinate_frame =
      cosmosim::core::InitialConditionCoordinateFrame::kPhysical;
  config.mode.ic_bridge_velocity_convention =
      cosmosim::core::InitialConditionVelocityConvention::kSqrtAScaledPeculiar;
  config.mode.ic_bridge_length_hubble_exponent = -1.0;
  config.mode.ic_bridge_mass_hubble_exponent = -1.0;
  config.units.length_unit = "kpc";
  config.cosmology.box_size_x_mpc_comoving = 0.05;
  config.cosmology.box_size_y_mpc_comoving = 0.05;
  config.cosmology.box_size_z_mpc_comoving = 0.05;
  config.cosmology.box_size_mpc_comoving = 0.05;
  config.cosmology.hubble_param = 0.5;
  config.numerics.a_begin = 0.5;
  const auto path = writeDimensionalGasBlackHoleIcFile();

  const auto direct = cosmosim::io::readGadgetArepoHdf5Ic(
      path, config, cosmosim::io::IcImportOptions{
          .validate_runtime_cosmology = false});
  assert(direct.report.manifest.has_value());
  assert(direct.state.particles.size() == 2U);
  assert(std::abs(direct.state.particles.position_x_comoving[0] - 4.0) < 1.0e-12);
  assert(std::abs(
      direct.state.particles.velocity_x_peculiar[0] -
      1.0 / std::sqrt(0.5)) < 1.0e-12);
  assert(std::abs(direct.state.particles.mass_code[0] - 2.0) < 1.0e-12);
  assert(std::abs(direct.state.gas_cells.internal_energy_code[0] - 4.0) < 1.0e-12);
  assert(std::abs(direct.state.gas_cells.density_code[0] - 0.0625) < 1.0e-12);
  assert(std::abs(direct.state.black_holes.subgrid_mass_code[0] - 6.0) < 1.0e-12);
  assert(std::abs(
      direct.state.black_holes.accretion_rate_code[0] - 0.25) < 1.0e-12);

  const cosmosim::io::IcManifest manifest = *direct.report.manifest;
  const auto manifest_driven = cosmosim::io::readGadgetArepoHdf5Ic(
      path, config, cosmosim::io::IcImportOptions{
          .validate_runtime_cosmology = false, .manifest = &manifest});
  assert(manifest_driven.state.particles.position_x_comoving ==
         direct.state.particles.position_x_comoving);
  assert(manifest_driven.state.particles.velocity_x_peculiar ==
         direct.state.particles.velocity_x_peculiar);
  assert(manifest_driven.state.particles.mass_code ==
         direct.state.particles.mass_code);
  assert(manifest_driven.state.gas_cells.internal_energy_code ==
         direct.state.gas_cells.internal_energy_code);
  assert(manifest_driven.state.gas_cells.density_code ==
         direct.state.gas_cells.density_code);
  assert(manifest_driven.state.black_holes.accretion_rate_code ==
         direct.state.black_holes.accretion_rate_code);

  auto physical_velocity_config = config;
  physical_velocity_config.mode.ic_bridge_velocity_convention =
      cosmosim::core::InitialConditionVelocityConvention::kPhysicalPeculiar;
  const auto physical_velocity = cosmosim::io::readGadgetArepoHdf5Ic(
      path, physical_velocity_config, cosmosim::io::IcImportOptions{
          .validate_runtime_cosmology = false});
  assert(std::abs(
      physical_velocity.state.particles.velocity_x_peculiar[0] - 1.0) <
      1.0e-12);
  assert(physical_velocity.state.gas_cells.internal_energy_code ==
         direct.state.gas_cells.internal_energy_code);
  assert(physical_velocity.state.black_holes.accretion_rate_code ==
         direct.state.black_holes.accretion_rate_code);
  std::filesystem::remove(path);
}

void testHdf5BlackHoleAndValidationFailures() {
  auto config = makeExplicitBridgeConfig();
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

  auto missing_mdot_path = writeMinimalBlackHoleIcFile();
  removeDataset(missing_mdot_path, "/PartType5", "BH_Mdot");
  expectIcReadFailure(
      missing_mdot_path, config,
      "/PartType5/BH_Mdot is missing and its normalized missing-field policy is reject");
  missing_mdot_path = writeMinimalBlackHoleIcFile();
  removeDataset(missing_mdot_path, "/PartType5", "BH_Mdot");
  config.mode.ic_bh_mdot_policy =
      cosmosim::core::InitialConditionMissingFieldPolicy::kDialectDefinedDefault;
  const auto defaulted_bh =
      cosmosim::io::readGadgetArepoHdf5Ic(missing_mdot_path, config);
  assert(defaulted_bh.state.black_holes.accretion_rate_code[0] == 0.0);
  assert(!defaulted_bh.report.defaulted_fields.empty());
  assert(defaulted_bh.report.manifest.has_value());
  config.mode.ic_bh_mdot_policy =
      cosmosim::core::InitialConditionMissingFieldPolicy::kReject;
  const auto manifest_defaulted_bh =
      cosmosim::io::readGadgetArepoHdf5Ic(
          missing_mdot_path, config,
          cosmosim::io::IcImportOptions{
              .manifest = &*defaulted_bh.report.manifest});
  assert(
      manifest_defaulted_bh.state.black_holes.accretion_rate_code[0] == 0.0);
  std::filesystem::remove(missing_mdot_path);

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
  auto family2_config = config;
  family2_config.mode.ic_part_type2_policy =
      cosmosim::core::InitialConditionSpeciesPolicy::kDarkMatter;
  const auto family2_inspection = cosmosim::io::readGadgetArepoHdf5Ic(
      family2_path, family2_config);
  assert(family2_inspection.report.manifest.has_value());
  cosmosim::io::IcManifest family2_manifest =
      *family2_inspection.report.manifest;
  assert(
      family2_manifest.species_policy[2] ==
      cosmosim::io::IcSpeciesPolicy::kCollisionlessFamily2AsDarkMatter);
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

void expectIcReadFailure(
    const std::filesystem::path& path,
    const cosmosim::core::SimulationConfig& config,
    std::string_view expected_text) {
  bool rejected = false;
  try {
    (void)cosmosim::io::readGadgetArepoHdf5Ic(path, config);
  } catch (const std::exception& error) {
    rejected = expected_text.empty() ||
        std::string(error.what()).find(expected_text) != std::string::npos;
  }
  assert(rejected);
  std::filesystem::remove(path);
}

void replaceNumFilesWithSignedValue(
    const std::filesystem::path& path, std::int64_t value) {
  Hdf5Handle file(H5Fopen(
      path.string().c_str(), H5F_ACC_RDWR, H5P_DEFAULT));
  Hdf5Handle header(H5Gopen2(file.get(), "/Header", H5P_DEFAULT));
  assert(H5Adelete(header.get(), "NumFilesPerSnapshot") >= 0);
  Hdf5Handle space(H5Screate(H5S_SCALAR));
  Hdf5Handle attribute(H5Acreate2(
      header.get(), "NumFilesPerSnapshot", H5T_STD_I64LE, space.get(),
      H5P_DEFAULT, H5P_DEFAULT));
  assert(H5Awrite(attribute.get(), H5T_NATIVE_INT64, &value) >= 0);
}

void replaceNumFilesWithUnsignedValue(
    const std::filesystem::path& path, std::uint64_t value) {
  Hdf5Handle file(H5Fopen(
      path.string().c_str(), H5F_ACC_RDWR, H5P_DEFAULT));
  Hdf5Handle header(H5Gopen2(file.get(), "/Header", H5P_DEFAULT));
  assert(H5Adelete(header.get(), "NumFilesPerSnapshot") >= 0);
  Hdf5Handle space(H5Screate(H5S_SCALAR));
  Hdf5Handle attribute(H5Acreate2(
      header.get(), "NumFilesPerSnapshot", H5T_STD_U64LE, space.get(),
      H5P_DEFAULT, H5P_DEFAULT));
  assert(H5Awrite(attribute.get(), H5T_NATIVE_UINT64, &value) >= 0);
}

void replaceNumFilesWithFloatingValue(
    const std::filesystem::path& path, double value) {
  Hdf5Handle file(H5Fopen(
      path.string().c_str(), H5F_ACC_RDWR, H5P_DEFAULT));
  Hdf5Handle header(H5Gopen2(file.get(), "/Header", H5P_DEFAULT));
  assert(H5Adelete(header.get(), "NumFilesPerSnapshot") >= 0);
  Hdf5Handle space(H5Screate(H5S_SCALAR));
  Hdf5Handle attribute(H5Acreate2(
      header.get(), "NumFilesPerSnapshot", H5T_IEEE_F64LE, space.get(),
      H5P_DEFAULT, H5P_DEFAULT));
  assert(H5Awrite(attribute.get(), H5T_NATIVE_DOUBLE, &value) >= 0);
}

void replaceThisFileWithSignedCounts(
    const std::filesystem::path& path,
    const std::array<std::int32_t, 6>& values) {
  Hdf5Handle file(H5Fopen(
      path.string().c_str(), H5F_ACC_RDWR, H5P_DEFAULT));
  Hdf5Handle header(H5Gopen2(file.get(), "/Header", H5P_DEFAULT));
  assert(H5Adelete(header.get(), "NumPart_ThisFile") >= 0);
  writeHeaderAttributeI32x6(header.get(), "NumPart_ThisFile", values);
}

void replaceHeaderCountAttributeWithShape(
    const std::filesystem::path& path,
    const char* name,
    std::size_t count) {
  Hdf5Handle file(
      H5Fopen(path.string().c_str(), H5F_ACC_RDWR, H5P_DEFAULT));
  Hdf5Handle header(H5Gopen2(file.get(), "/Header", H5P_DEFAULT));
  assert(H5Adelete(header.get(), name) >= 0);
  const hsize_t dims[1]{static_cast<hsize_t>(count)};
  Hdf5Handle space(H5Screate_simple(1, dims, nullptr));
  Hdf5Handle attribute(H5Acreate2(
      header.get(), name, H5T_STD_U32LE, space.get(), H5P_DEFAULT,
      H5P_DEFAULT));
  std::vector<std::uint32_t> values(count, 0U);
  assert(H5Awrite(attribute.get(), H5T_NATIVE_UINT32, values.data()) >= 0);
}

void replaceDatasetWithFloatingIds(
    const std::filesystem::path& path,
    bool signed_negative) {
  Hdf5Handle file(
      H5Fopen(path.string().c_str(), H5F_ACC_RDWR, H5P_DEFAULT));
  Hdf5Handle gas(H5Gopen2(file.get(), "/PartType0", H5P_DEFAULT));
  assert(H5Ldelete(gas.get(), "ParticleIDs", H5P_DEFAULT) >= 0);
  hsize_t dims[1]{2U};
  Hdf5Handle space(H5Screate_simple(1, dims, nullptr));
  if (signed_negative) {
    Hdf5Handle dataset(H5Dcreate2(
        gas.get(), "ParticleIDs", H5T_STD_I64LE, space.get(), H5P_DEFAULT,
        H5P_DEFAULT, H5P_DEFAULT));
    const std::array<std::int64_t, 2> values{-1, 2};
    assert(H5Dwrite(
               dataset.get(), H5T_NATIVE_INT64, H5S_ALL, H5S_ALL,
               H5P_DEFAULT, values.data()) >= 0);
  } else {
    Hdf5Handle dataset(H5Dcreate2(
        gas.get(), "ParticleIDs", H5T_IEEE_F64LE, space.get(), H5P_DEFAULT,
        H5P_DEFAULT, H5P_DEFAULT));
    const std::array<double, 2> values{1.0, 2.0};
    assert(H5Dwrite(
               dataset.get(), H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
               H5P_DEFAULT, values.data()) >= 0);
  }
}


std::filesystem::path writeTracerPolicyIcFile(bool valid_parent) {
  const auto path = std::filesystem::temp_directory_path() /
      (valid_parent ? "cosmosim_ic_tracer_parent_valid.hdf5"
                    : "cosmosim_ic_tracer_parent_invalid.hdf5");
  Hdf5Handle file(H5Fcreate(
      path.string().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT));
  Hdf5Handle header(H5Gcreate2(
      file.get(), "/Header", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
  writeRequiredHeader(header.get(), {1, 0, 1, 0, 0, 0});

  Hdf5Handle gas(H5Gcreate2(
      file.get(), "/PartType0", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
  writeDataset2dVec3(gas.get(), "Coordinates", {1.0, 2.0, 3.0});
  writeDataset2dVec3(gas.get(), "Velocities", {0.0, 0.0, 0.0});
  writeDataset1d(gas.get(), "Masses", {2.0});
  writeDataset1dIds(gas.get(), "ParticleIDs", {8001U});
  writeDataset1d(gas.get(), "InternalEnergy", {1.0});
  writeDataset1d(gas.get(), "Density", {1.0});

  Hdf5Handle tracer(H5Gcreate2(
      file.get(), "/PartType2", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
  writeDataset2dVec3(tracer.get(), "Coordinates", {1.0, 2.0, 3.0});
  writeDataset2dVec3(tracer.get(), "Velocities", {0.0, 0.0, 0.0});
  writeDataset1d(tracer.get(), "Masses", {0.1});
  writeDataset1dIds(tracer.get(), "ParticleIDs", {9001U});
  writeDataset1dIds(
      tracer.get(), "ParentParticleIDs",
      {valid_parent ? 8001U : 9999U});
  writeDataset1dIds(tracer.get(), "InjectionStep", {7U});
  writeDataset1dIds(tracer.get(), "HostCellIndex", {42U});
  writeDataset1d(tracer.get(), "MassFractionOfHost", {0.05});
  writeDataset1d(tracer.get(), "LastHostMass", {2.0});
  writeDataset1d(tracer.get(), "CumulativeExchangedMass", {0.0});
  return path;
}

void testTracerParentAndHostRemap() {
  auto config = makeExplicitBridgeConfig();
  config.mode.ic_part_type2_policy =
      cosmosim::core::InitialConditionSpeciesPolicy::kTracer;

  const auto valid_path = writeTracerPolicyIcFile(true);
  const auto result = cosmosim::io::readGadgetArepoHdf5Ic(valid_path, config);
  assert(result.state.tracers.size() == 1U);
  assert(result.state.cells.size() == 1U);
  assert(result.state.tracers.parent_particle_id[0] == 8001U);
  assert(result.state.tracers.host_cell_index[0] == 0U);
  assert(result.report.manifest.has_value());
  assert(std::find_if(
             result.report.manifest->dropped_fields.begin(),
             result.report.manifest->dropped_fields.end(),
             [](const std::string& value) {
               return value.find("HostCellIndex") != std::string::npos;
             }) != result.report.manifest->dropped_fields.end());
  std::filesystem::remove(valid_path);

  const auto invalid_path = writeTracerPolicyIcFile(false);
  expectIcReadFailure(
      invalid_path, config,
      "tracer parent must resolve to a gas cell on the same final owner rank");
}


void testReaderSessionSourceIdentity() {
  using cosmosim::io::IcImportCounters;
  using cosmosim::io::internal::IcReaderSession;

  auto path = writeMinimalIcFile(true);
  const std::uint64_t size = std::filesystem::file_size(path);
  const std::string sha = cosmosim::io::icSha256FileHex(path);
  {
    IcImportCounters counters;
    IcReaderSession session(path, size, sha, counters);
    session.revalidateSourceIdentity(counters);
    assert(counters.source_file_open_count == 1U);
    assert(counters.full_file_hash_pass_count == 2U);
    assert(counters.source_identity_validation_count == 2U);
    assert(counters.hash_bytes_read == 2U * size);
  }

  bool hash_mismatch_rejected = false;
  try {
    IcImportCounters counters;
    IcReaderSession session(path, size, std::string(64U, '0'), counters);
    static_cast<void>(session);
  } catch (const std::runtime_error& error) {
    hash_mismatch_rejected =
        std::string(error.what()).find("SHA-256 mismatch") !=
        std::string::npos;
  }
  assert(hash_mismatch_rejected);

  {
    IcImportCounters counters;
    IcReaderSession session(path, size, sha, counters);
    std::ofstream append(path, std::ios::binary | std::ios::app);
    append.put('\0');
    append.close();
    bool size_change_rejected = false;
    try {
      session.revalidateSourceIdentity(counters);
    } catch (const std::runtime_error& error) {
      size_change_rejected =
          std::string(error.what()).find("identity changed") !=
              std::string::npos ||
          std::string(error.what()).find("size changed") !=
              std::string::npos;
    }
    assert(size_change_rejected);
  }
  std::filesystem::remove(path);

  path = writeMinimalIcFile(true);
  const auto replacement = std::filesystem::path(path.string() + ".replacement");
  const auto displaced = std::filesystem::path(path.string() + ".displaced");
  std::filesystem::copy_file(path, replacement);
  const std::uint64_t replacement_size = std::filesystem::file_size(path);
  const std::string replacement_sha = cosmosim::io::icSha256FileHex(path);
  {
    IcImportCounters counters;
    IcReaderSession session(
        path, replacement_size, replacement_sha, counters);
    std::filesystem::rename(path, displaced);
    std::filesystem::rename(replacement, path);
    bool replacement_rejected = false;
    try {
      session.revalidateSourceIdentity(counters);
    } catch (const std::runtime_error& error) {
      replacement_rejected =
          std::string(error.what()).find("identity changed") !=
          std::string::npos;
    }
    assert(replacement_rejected);
  }
  std::filesystem::remove(path);
  std::filesystem::remove(displaced);

  path = writeMinimalIcFile(true);
  const std::uint64_t mutation_size = std::filesystem::file_size(path);
  const std::string mutation_sha = cosmosim::io::icSha256FileHex(path);
  {
    IcImportCounters counters;
    IcReaderSession session(path, mutation_size, mutation_sha, counters);
    std::fstream mutation(
        path, std::ios::binary | std::ios::in | std::ios::out);
    mutation.seekg(-1, std::ios::end);
    char value = 0;
    mutation.read(&value, 1);
    value = static_cast<char>(value ^ 0x01);
    mutation.seekp(-1, std::ios::end);
    mutation.write(&value, 1);
    mutation.close();
    bool mutation_rejected = false;
    try {
      session.revalidateSourceIdentity(counters);
    } catch (const std::runtime_error& error) {
      mutation_rejected =
          std::string(error.what()).find("identity changed") !=
              std::string::npos ||
          std::string(error.what()).find("content changed") !=
              std::string::npos;
    }
    assert(mutation_rejected);
  }
  std::filesystem::remove(path);
}

void testHdf5MalformedSchemaSafety() {
  const auto config = makeExplicitBridgeConfig();

  auto path = writeMinimalIcFile(true);
  replaceHeaderCountAttributeWithShape(path, "NumPart_ThisFile", 7U);
  expectIcReadFailure(path, config, "expected [6]");

  path = writeMinimalIcFile(true);
  replaceHeaderCountAttributeWithShape(path, "NumPart_ThisFile", 5U);
  expectIcReadFailure(path, config, "expected [6]");

  path = writeMinimalIcFile(true);
  replaceThisFileWithSignedCounts(path, {-1, 0, 0, 0, 0, 0});
  expectIcReadFailure(path, config, "negative particle count");

  path = writeMinimalIcFile(true);
  replaceNumFilesWithSignedValue(path, -1);
  expectIcReadFailure(path, config, "must be non-negative");

  path = writeMinimalIcFile(true);
  replaceNumFilesWithUnsignedValue(
      path, static_cast<std::uint64_t>(UINT32_MAX) + 1ULL);
  expectIcReadFailure(path, config, "exceeds uint32 range");

  path = writeMinimalIcFile(true);
  replaceNumFilesWithFloatingValue(path, 1.0);
  expectIcReadFailure(path, config, "invalid datatype class");

  path = writeMinimalIcFile(true);
  {
    Hdf5Handle file(
        H5Fopen(path.string().c_str(), H5F_ACC_RDWR, H5P_DEFAULT));
    Hdf5Handle header(H5Gopen2(file.get(), "/Header", H5P_DEFAULT));
    assert(H5Adelete(header.get(), "Time") >= 0);
    hsize_t dims[1]{1U};
    Hdf5Handle space(H5Screate_simple(1, dims, nullptr));
    Hdf5Handle attribute(H5Acreate2(
        header.get(), "Time", H5T_IEEE_F64LE, space.get(), H5P_DEFAULT,
        H5P_DEFAULT));
    const double value = 1.0;
    assert(H5Awrite(attribute.get(), H5T_NATIVE_DOUBLE, &value) >= 0);
  }
  expectIcReadFailure(path, config, "expected []");

  path = writeMinimalIcFile(true);
  replaceDatasetWithFloatingIds(path, false);
  expectIcReadFailure(path, config, "unsigned integer datatype");

  path = writeMinimalIcFile(true);
  replaceDatasetWithFloatingIds(path, true);
  expectIcReadFailure(path, config, "unsigned integer datatype");

  path = writeMinimalIcFile(true);
  {
    Hdf5Handle file(
        H5Fopen(path.string().c_str(), H5F_ACC_RDWR, H5P_DEFAULT));
    Hdf5Handle group(H5Gcreate2(
        file.get(), "/PartType6", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
    writeDataset1d(group.get(), "Masses", {1.0});
  }
  expectIcReadFailure(path, config, "unsupported populated particle-family");

  path = writeMinimalIcFile(true);
  {
    Hdf5Handle file(
        H5Fopen(path.string().c_str(), H5F_ACC_RDWR, H5P_DEFAULT));
    Hdf5Handle gas(H5Gopen2(file.get(), "/PartType0", H5P_DEFAULT));
    writeDataset2dVec3(
        gas.get(), "Position", {1.0, 2.0, 3.0, 4.0, 5.0, 6.0});
  }
  expectIcReadFailure(path, config, "ambiguous aliases");

  path = writeMinimalIcFile(true);
  {
    Hdf5Handle file(
        H5Fopen(path.string().c_str(), H5F_ACC_RDWR, H5P_DEFAULT));
    Hdf5Handle gas(H5Gopen2(file.get(), "/PartType0", H5P_DEFAULT));
    writeDataset1d(gas.get(), "UnrecognizedScalar", {3.0, 4.0});
  }
  const auto result = cosmosim::io::readGadgetArepoHdf5Ic(path, config);
  assert(result.report.manifest.has_value());
  const auto& fields = result.report.manifest->fields;
  const auto dropped = std::find_if(
      fields.begin(), fields.end(), [](const cosmosim::io::IcFieldManifest& field) {
        return field.dataset_path == "/PartType0/UnrecognizedScalar";
      });
  assert(dropped != fields.end());
  assert(dropped->disposition == cosmosim::io::IcFieldDisposition::kDropped);
  assert(std::find_if(
             result.report.manifest->dropped_fields.begin(),
             result.report.manifest->dropped_fields.end(),
             [](const std::string& value) {
               return value.find("UnrecognizedScalar") != std::string::npos;
             }) != result.report.manifest->dropped_fields.end());
  std::filesystem::remove(path);
}

#endif

}  // namespace

int main() {
  testCanonicalSingleFileCountLimit();
  testManifestValidationAndConversions();
  testGeneratedIsolatedIcSpeciesAndOwnership();
  testGeneratedConverterDefaultAudit();
  testHdf5GateBehavior();
#if COSMOSIM_ENABLE_HDF5
  testCanonicalHeaderContract();
  testHdf5StarSidecarAndMultifileSchema();
  testHdf5GasThermoMapping();
  testHdf5GasOptionalDensityMissingBehavior();
  testSharedDimensionalConversionContract();
  testHdf5BlackHoleAndValidationFailures();
  testTracerParentAndHostRemap();
  testReaderSessionSourceIdentity();
  testHdf5MalformedSchemaSafety();
#endif
  return 0;
}
