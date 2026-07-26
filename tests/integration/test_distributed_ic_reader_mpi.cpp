#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <limits>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include <hdf5.h>
#include <mpi.h>

#include "cosmosim/core/config.hpp"
#include "cosmosim/io/ic_reader.hpp"
#include "cosmosim/parallel/distributed_memory.hpp"

namespace {

constexpr double K_MPC_TO_SI = 3.0856775814913673e22;
constexpr double K_SOLAR_MASS_TO_SI = 1.98847e30;
constexpr double K_KM_S_TO_SI = 1000.0;
constexpr std::uint32_t K_MEMBER_COUNT = 2U;
constexpr std::array<std::uint32_t, 6> K_LOCAL_COUNTS{
    8U, 16U, 0U, 0U, 4U, 2U};
constexpr std::array<std::uint32_t, 6> K_TOTAL_COUNTS{
    16U, 32U, 0U, 0U, 8U, 4U};
constexpr std::uint64_t K_GLOBAL_PARTICLE_COUNT = 60U;
constexpr double K_BOX_SIZE = 64.0;

class Hdf5Handle {
 public:
  explicit Hdf5Handle(hid_t value = -1) : m_value(value) {}
  Hdf5Handle(const Hdf5Handle&) = delete;
  Hdf5Handle& operator=(const Hdf5Handle&) = delete;
  ~Hdf5Handle() {
    if (m_value < 0) {
      return;
    }
    switch (H5Iget_type(m_value)) {
      case H5I_FILE:
        H5Fclose(m_value);
        break;
      case H5I_GROUP:
        H5Gclose(m_value);
        break;
      case H5I_DATASET:
        H5Dclose(m_value);
        break;
      case H5I_DATASPACE:
        H5Sclose(m_value);
        break;
      case H5I_ATTR:
        H5Aclose(m_value);
        break;
      default:
        break;
    }
  }
  [[nodiscard]] hid_t get() const noexcept { return m_value; }

 private:
  hid_t m_value = -1;
};

template <typename T>
void writeAttribute(hid_t group, const char* name, hid_t type, const T& value) {
  Hdf5Handle space(H5Screate(H5S_SCALAR));
  Hdf5Handle attr(
      H5Acreate2(group, name, type, space.get(), H5P_DEFAULT, H5P_DEFAULT));
  assert(attr.get() >= 0);
  assert(H5Awrite(attr.get(), type, &value) >= 0);
}

template <typename T, std::size_t N>
void writeArrayAttribute(
    hid_t group,
    const char* name,
    hid_t type,
    const std::array<T, N>& values) {
  const hsize_t dims[1] = {N};
  Hdf5Handle space(H5Screate_simple(1, dims, nullptr));
  Hdf5Handle attr(
      H5Acreate2(group, name, type, space.get(), H5P_DEFAULT, H5P_DEFAULT));
  assert(attr.get() >= 0);
  assert(H5Awrite(attr.get(), type, values.data()) >= 0);
}

void writeStringAttribute(
    hid_t group,
    const char* name,
    std::string_view value) {
  Hdf5Handle type(H5Tcopy(H5T_C_S1));
  assert(type.get() >= 0);
  assert(H5Tset_size(type.get(), value.size() + 1U) >= 0);
  Hdf5Handle space(H5Screate(H5S_SCALAR));
  Hdf5Handle attr(H5Acreate2(
      group, name, type.get(), space.get(), H5P_DEFAULT, H5P_DEFAULT));
  assert(attr.get() >= 0);
  std::string storage(value);
  assert(H5Awrite(attr.get(), type.get(), storage.c_str()) >= 0);
}

void writeDataset1d(
    hid_t group,
    const char* name,
    hid_t file_type,
    hid_t memory_type,
    const void* values,
    std::size_t count) {
  const hsize_t dims[1] = {count};
  Hdf5Handle space(H5Screate_simple(1, dims, nullptr));
  Hdf5Handle dataset(H5Dcreate2(
      group, name, file_type, space.get(), H5P_DEFAULT, H5P_DEFAULT,
      H5P_DEFAULT));
  assert(dataset.get() >= 0);
  assert(H5Dwrite(
             dataset.get(), memory_type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
             values) >= 0);
}

void writeByteDataset(
    hid_t group, const char* name, const std::string& values) {
  const hsize_t dims[1] = {values.size()};
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

void overwriteStringAttribute(
    const std::filesystem::path& path,
    const char* attribute_name,
    std::string_view value) {
  Hdf5Handle file(
      H5Fopen(path.string().c_str(), H5F_ACC_RDWR, H5P_DEFAULT));
  assert(file.get() >= 0);
  Hdf5Handle header(H5Gopen2(file.get(), "/Header", H5P_DEFAULT));
  assert(header.get() >= 0);
  Hdf5Handle attribute(H5Aopen(header.get(), attribute_name, H5P_DEFAULT));
  assert(attribute.get() >= 0);
  Hdf5Handle type(H5Aget_type(attribute.get()));
  assert(type.get() >= 0);
  std::string storage(value);
  assert(H5Awrite(attribute.get(), type.get(), storage.c_str()) >= 0);
}

[[nodiscard]] std::string canonicalAuditJson() {
  cosmosim::io::IcManifest manifest;
  manifest.source_files = {"clean_room_arepo_source.hdf5"};
  manifest.source_sha256 = {std::string(64U, 'a')};
  manifest.source_provenance_ids = {
      "sha256:" + manifest.source_sha256.front()};
  manifest.source_file_sizes_bytes = {1U};
  manifest.original_header_attributes = {"clean-room canonical fixture"};
  manifest.num_files_per_snapshot = 1U;
  manifest.num_part_this_file = {std::array<std::uint64_t, 6>{}};
  manifest.num_part_total = {};
  manifest.num_part_total_high_word = {};
  manifest.mass_table = {};
  manifest.box_size = K_BOX_SIZE;
  manifest.scale_factor = 1.0;
  manifest.redshift = 0.0;
  manifest.omega_matter = 0.315;
  manifest.omega_lambda = 0.685;
  manifest.hubble_param = 0.674;
  return cosmosim::io::serializeIcManifestJson(manifest);
}

void writeDatasetVec3(
    hid_t group,
    const char* name,
    const std::vector<float>& values) {
  const hsize_t dims[2] = {values.size() / 3U, 3U};
  Hdf5Handle space(H5Screate_simple(2, dims, nullptr));
  Hdf5Handle dataset(H5Dcreate2(
      group, name, H5T_IEEE_F32LE, space.get(), H5P_DEFAULT, H5P_DEFAULT,
      H5P_DEFAULT));
  assert(dataset.get() >= 0);
  assert(H5Dwrite(
             dataset.get(), H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,
             H5P_DEFAULT, values.data()) >= 0);
}

[[nodiscard]] float positionX(std::uint32_t global_slot) {
  return static_cast<float>(
      (static_cast<double>(global_slot) + 0.5) *
      (K_BOX_SIZE / static_cast<double>(K_GLOBAL_PARTICLE_COUNT)));
}

void writeCommonParticleFields(
    hid_t group,
    std::uint32_t member,
    std::uint32_t local_count,
    std::uint32_t global_slot_begin,
    std::uint64_t id_begin,
    float mass,
    bool duplicate_last_id = false) {
  std::vector<float> coordinates;
  std::vector<float> velocities;
  std::vector<float> masses(local_count, mass);
  std::vector<std::uint64_t> ids(local_count);
  coordinates.reserve(static_cast<std::size_t>(local_count) * 3U);
  velocities.reserve(static_cast<std::size_t>(local_count) * 3U);
  for (std::uint32_t local_index = 0U; local_index < local_count;
       ++local_index) {
    const std::uint32_t member_offset = member * local_count;
    const std::uint32_t type_index = member_offset + local_index;
    coordinates.insert(
        coordinates.end(),
        {positionX(global_slot_begin + type_index),
         1.0F + static_cast<float>(member),
         2.0F + static_cast<float>(local_index % 7U)});
    velocities.insert(
        velocities.end(),
        {0.1F * static_cast<float>(local_index + 1U), 0.0F, 0.0F});
    ids[local_index] = id_begin + type_index;
  }
  if (duplicate_last_id && !ids.empty()) {
    ids.back() = 1000U;
  }
  writeDatasetVec3(group, "Coordinates", coordinates);
  writeDatasetVec3(group, "Velocities", velocities);
  writeDataset1d(
      group, "Masses", H5T_IEEE_F32LE, H5T_NATIVE_FLOAT, masses.data(),
      masses.size());
  writeDataset1d(
      group, "ParticleIDs", H5T_STD_U64LE, H5T_NATIVE_UINT64, ids.data(),
      ids.size());
}

void writeMember(
    const std::filesystem::path& path,
    std::uint32_t member,
    bool duplicate_ids,
    bool canonical_header = false,
    bool alternate_convention = false) {
  Hdf5Handle file(
      H5Fcreate(path.string().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT));
  Hdf5Handle header(
      H5Gcreate2(file.get(), "/Header", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
  const std::array<std::uint32_t, 6> high{};
  const std::array<double, 6> mass_table{};
  std::array<std::int32_t, 6> signed_local_counts{};
  std::transform(
      K_LOCAL_COUNTS.begin(), K_LOCAL_COUNTS.end(),
      signed_local_counts.begin(),
      [](std::uint32_t value) { return static_cast<std::int32_t>(value); });
  writeArrayAttribute(
      header.get(), "NumPart_ThisFile", H5T_NATIVE_INT32,
      signed_local_counts);
  writeArrayAttribute(
      header.get(), "NumPart_Total", H5T_NATIVE_UINT32, K_TOTAL_COUNTS);
  writeArrayAttribute(
      header.get(), "NumPart_Total_HighWord", H5T_NATIVE_UINT32, high);
  writeArrayAttribute(
      header.get(), "MassTable", H5T_NATIVE_DOUBLE, mass_table);
  const double scale_factor = alternate_convention ? 0.5 : 1.0;
  const double redshift = alternate_convention ? 1.0 : 0.0;
  const double hubble_param = alternate_convention ? 0.5 : 0.674;
  writeAttribute(header.get(), "Time", H5T_NATIVE_DOUBLE, scale_factor);
  writeAttribute(header.get(), "Redshift", H5T_NATIVE_DOUBLE, redshift);
  writeAttribute(header.get(), "BoxSize", H5T_NATIVE_DOUBLE, K_BOX_SIZE);
  writeAttribute(header.get(), "Omega0", H5T_NATIVE_DOUBLE, 0.315);
  writeAttribute(header.get(), "OmegaLambda", H5T_NATIVE_DOUBLE, 0.685);
  writeAttribute(header.get(), "HubbleParam", H5T_NATIVE_DOUBLE, hubble_param);
  writeAttribute(
      header.get(), "NumFilesPerSnapshot", H5T_NATIVE_INT32,
      static_cast<std::int32_t>(K_MEMBER_COUNT));
  if (canonical_header) {
    writeStringAttribute(header.get(), "ChuiIcSchemaName", "chui_canonical_v1");
    const std::string audit_json = canonicalAuditJson();
    const std::string audit_digest = cosmosim::io::icSha256Hex(audit_json);
    const std::filesystem::path sidecar_path =
        path.string() + ".manifest.json";
    const std::filesystem::path marker_path = path.string() + ".complete";
    writeAttribute(
        header.get(), "ChuiIcSchemaVersion", H5T_NATIVE_UINT32,
        std::uint32_t{2U});
    writeStringAttribute(header.get(), "ChuiCoordinateFrame", "comoving");
    writeStringAttribute(
        header.get(), "ChuiVelocityConvention", "physical_peculiar");
    writeStringAttribute(
        header.get(), "ConversionManifestSha256", audit_digest);
    writeStringAttribute(
        header.get(), "ConversionManifestSidecar",
        sidecar_path.filename().string());
    writeStringAttribute(
        header.get(), "ConversionBundleMarker",
        marker_path.filename().string());
    Hdf5Handle provenance(H5Gcreate2(
        file.get(), "/Provenance", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
    writeByteDataset(provenance.get(), "ConversionManifestJson", audit_json);
    writeTextFile(sidecar_path, audit_json);
    writeTextFile(
        marker_path,
        "chui_ic_bundle_v1\nsha256=" + audit_digest +
            "\ncanonical=" + path.filename().string() +
            "\nmanifest=" + sidecar_path.filename().string() + "\n");
    writeAttribute(
        header.get(), "ChuiLengthUnitToSI", H5T_NATIVE_DOUBLE, K_MPC_TO_SI);
    writeAttribute(
        header.get(), "ChuiMassUnitToSI", H5T_NATIVE_DOUBLE,
        K_SOLAR_MASS_TO_SI);
    writeAttribute(
        header.get(), "ChuiVelocityUnitToSI", H5T_NATIVE_DOUBLE,
        K_KM_S_TO_SI);
  }

  Hdf5Handle gas(H5Gcreate2(
      file.get(), "/PartType0", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
  writeCommonParticleFields(
      gas.get(), member, K_LOCAL_COUNTS[0], 0U, 1000U, 1.0F,
      duplicate_ids && member == 1U);
  std::vector<float> internal_energy(K_LOCAL_COUNTS[0]);
  std::vector<float> density(K_LOCAL_COUNTS[0]);
  std::vector<float> metallicity(K_LOCAL_COUNTS[0], 0.02F);
  for (std::size_t i = 0U; i < internal_energy.size(); ++i) {
    internal_energy[i] = 4.0F + static_cast<float>(i);
    density[i] = 2.0F + 0.1F * static_cast<float>(i);
  }
  writeDataset1d(
      gas.get(), "InternalEnergy", H5T_IEEE_F32LE, H5T_NATIVE_FLOAT,
      internal_energy.data(), internal_energy.size());
  writeDataset1d(
      gas.get(), "Density", H5T_IEEE_F32LE, H5T_NATIVE_FLOAT, density.data(),
      density.size());
  writeDataset1d(
      gas.get(), "Metallicity", H5T_IEEE_F32LE, H5T_NATIVE_FLOAT,
      metallicity.data(), metallicity.size());

  Hdf5Handle dm(H5Gcreate2(
      file.get(), "/PartType1", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
  writeCommonParticleFields(
      dm.get(), member, K_LOCAL_COUNTS[1], 16U, 2000U, 2.0F);

  Hdf5Handle star(H5Gcreate2(
      file.get(), "/PartType4", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
  writeCommonParticleFields(
      star.get(), member, K_LOCAL_COUNTS[4], 48U, 4000U, 3.0F);
  std::vector<float> formation(K_LOCAL_COUNTS[4], 0.5F);
  std::vector<float> initial_mass(K_LOCAL_COUNTS[4], 3.5F);
  std::vector<float> star_metallicity(K_LOCAL_COUNTS[4], 0.015F);
  writeDataset1d(
      star.get(), "GFM_StellarFormationTime", H5T_IEEE_F32LE,
      H5T_NATIVE_FLOAT, formation.data(), formation.size());
  writeDataset1d(
      star.get(), "GFM_InitialMass", H5T_IEEE_F32LE, H5T_NATIVE_FLOAT,
      initial_mass.data(), initial_mass.size());
  writeDataset1d(
      star.get(), "GFM_Metallicity", H5T_IEEE_F32LE, H5T_NATIVE_FLOAT,
      star_metallicity.data(), star_metallicity.size());

  Hdf5Handle black_hole(H5Gcreate2(
      file.get(), "/PartType5", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
  writeCommonParticleFields(
      black_hole.get(), member, K_LOCAL_COUNTS[5], 56U, 5000U, 4.0F);
  std::vector<float> black_hole_mass(K_LOCAL_COUNTS[5], 5.0F);
  std::vector<float> black_hole_mdot(K_LOCAL_COUNTS[5], 0.25F);
  writeDataset1d(
      black_hole.get(), "BH_Mass", H5T_IEEE_F32LE, H5T_NATIVE_FLOAT,
      black_hole_mass.data(), black_hole_mass.size());
  writeDataset1d(
      black_hole.get(), "BH_Mdot", H5T_IEEE_F32LE, H5T_NATIVE_FLOAT,
      black_hole_mdot.data(), black_hole_mdot.size());
}

void writeScalingMember(
    const std::filesystem::path& path,
    std::uint32_t member,
    std::uint32_t local_dm_count) {
  Hdf5Handle file(
      H5Fcreate(path.string().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT));
  Hdf5Handle header(
      H5Gcreate2(file.get(), "/Header", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
  std::array<std::uint32_t, 6> local{};
  std::array<std::uint32_t, 6> total{};
  local[1] = local_dm_count;
  total[1] = local_dm_count * K_MEMBER_COUNT;
  const std::array<std::uint32_t, 6> high{};
  const std::array<double, 6> mass_table{};
  std::array<std::int32_t, 6> signed_local{};
  std::transform(
      local.begin(), local.end(), signed_local.begin(),
      [](std::uint32_t value) { return static_cast<std::int32_t>(value); });
  writeArrayAttribute(
      header.get(), "NumPart_ThisFile", H5T_NATIVE_INT32, signed_local);
  writeArrayAttribute(
      header.get(), "NumPart_Total", H5T_NATIVE_UINT32, total);
  writeArrayAttribute(
      header.get(), "NumPart_Total_HighWord", H5T_NATIVE_UINT32, high);
  writeArrayAttribute(
      header.get(), "MassTable", H5T_NATIVE_DOUBLE, mass_table);
  writeAttribute(header.get(), "Time", H5T_NATIVE_DOUBLE, 1.0);
  writeAttribute(header.get(), "Redshift", H5T_NATIVE_DOUBLE, 0.0);
  writeAttribute(header.get(), "BoxSize", H5T_NATIVE_DOUBLE, K_BOX_SIZE);
  writeAttribute(header.get(), "Omega0", H5T_NATIVE_DOUBLE, 0.315);
  writeAttribute(header.get(), "OmegaLambda", H5T_NATIVE_DOUBLE, 0.685);
  writeAttribute(header.get(), "HubbleParam", H5T_NATIVE_DOUBLE, 0.674);
  writeAttribute(
      header.get(), "NumFilesPerSnapshot", H5T_NATIVE_INT32,
      static_cast<std::int32_t>(K_MEMBER_COUNT));

  Hdf5Handle dm(H5Gcreate2(
      file.get(), "/PartType1", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
  const std::uint64_t global_count =
      static_cast<std::uint64_t>(local_dm_count) * K_MEMBER_COUNT;
  std::vector<float> coordinates;
  std::vector<float> velocities(
      static_cast<std::size_t>(local_dm_count) * 3U, 0.0F);
  std::vector<float> masses(local_dm_count, 2.0F);
  std::vector<std::uint64_t> ids(local_dm_count);
  coordinates.reserve(static_cast<std::size_t>(local_dm_count) * 3U);
  for (std::uint32_t local_index = 0U; local_index < local_dm_count;
       ++local_index) {
    const std::uint64_t global_index =
        static_cast<std::uint64_t>(member) * local_dm_count + local_index;
    const float x = static_cast<float>(
        (static_cast<double>(global_index) + 0.5) *
        (K_BOX_SIZE / static_cast<double>(global_count)));
    coordinates.insert(coordinates.end(), {x, 1.0F, 2.0F});
    ids[local_index] = 100000U + global_index;
  }
  writeDatasetVec3(dm.get(), "Coordinates", coordinates);
  writeDatasetVec3(dm.get(), "Velocities", velocities);
  writeDataset1d(
      dm.get(), "Masses", H5T_IEEE_F32LE, H5T_NATIVE_FLOAT, masses.data(),
      masses.size());
  writeDataset1d(
      dm.get(), "ParticleIDs", H5T_STD_U64LE, H5T_NATIVE_UINT64, ids.data(),
      ids.size());
}

[[nodiscard]] cosmosim::io::IcReadResult runScalingImport(
    const std::filesystem::path& base,
    std::uint32_t local_dm_count,
    const cosmosim::core::SimulationConfig& config,
    const cosmosim::parallel::MpiContext& mpi_context) {
  const auto first = std::filesystem::path(base.string() + ".0.hdf5");
  const auto second = std::filesystem::path(base.string() + ".1.hdf5");
  if (mpi_context.isRoot()) {
    writeScalingMember(first, 0U, local_dm_count);
    writeScalingMember(second, 1U, local_dm_count);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  auto result = cosmosim::io::readDistributedGadgetArepoHdf5Ic(
      first, config, mpi_context,
      cosmosim::io::IcImportOptions{.chunk_particle_count = 7U});
  MPI_Barrier(MPI_COMM_WORLD);
  if (mpi_context.isRoot()) {
    std::filesystem::remove(first);
    std::filesystem::remove(second);
    std::filesystem::remove(first.string() + ".manifest.json");
    std::filesystem::remove(second.string() + ".manifest.json");
    std::filesystem::remove(first.string() + ".complete");
    std::filesystem::remove(second.string() + ".complete");
  }
  return result;
}


void writePolicyMember(
    const std::filesystem::path& path,
    std::uint32_t member) {
  Hdf5Handle file(
      H5Fcreate(path.string().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT));
  Hdf5Handle header(
      H5Gcreate2(file.get(), "/Header", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
  std::array<std::uint32_t, 6> local{};
  std::array<std::uint32_t, 6> total{};
  local[0] = 3U;
  local[2] = 3U;
  local[3] = 2U;
  total[0] = 6U;
  total[2] = 6U;
  total[3] = 4U;
  const std::array<std::uint32_t, 6> high{};
  const std::array<double, 6> mass_table{};
  std::array<std::int32_t, 6> signed_local{};
  std::transform(
      local.begin(), local.end(), signed_local.begin(),
      [](std::uint32_t value) { return static_cast<std::int32_t>(value); });
  writeArrayAttribute(
      header.get(), "NumPart_ThisFile", H5T_NATIVE_INT32, signed_local);
  writeArrayAttribute(
      header.get(), "NumPart_Total", H5T_NATIVE_UINT32, total);
  writeArrayAttribute(
      header.get(), "NumPart_Total_HighWord", H5T_NATIVE_UINT32, high);
  writeArrayAttribute(
      header.get(), "MassTable", H5T_NATIVE_DOUBLE, mass_table);
  writeAttribute(header.get(), "Time", H5T_NATIVE_DOUBLE, 1.0);
  writeAttribute(header.get(), "Redshift", H5T_NATIVE_DOUBLE, 0.0);
  writeAttribute(header.get(), "BoxSize", H5T_NATIVE_DOUBLE, K_BOX_SIZE);
  writeAttribute(header.get(), "Omega0", H5T_NATIVE_DOUBLE, 0.315);
  writeAttribute(header.get(), "OmegaLambda", H5T_NATIVE_DOUBLE, 0.685);
  writeAttribute(header.get(), "HubbleParam", H5T_NATIVE_DOUBLE, 0.674);
  writeAttribute(
      header.get(), "NumFilesPerSnapshot", H5T_NATIVE_INT32,
      static_cast<std::int32_t>(K_MEMBER_COUNT));

  Hdf5Handle gas(H5Gcreate2(
      file.get(), "/PartType0", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
  writeCommonParticleFields(gas.get(), member, 3U, 0U, 8000U, 3.0F);
  std::vector<float> internal_energy(3U, 2.0F);
  std::vector<float> density(3U, 1.0F);
  writeDataset1d(
      gas.get(), "InternalEnergy", H5T_IEEE_F32LE, H5T_NATIVE_FLOAT,
      internal_energy.data(), internal_energy.size());
  writeDataset1d(
      gas.get(), "Density", H5T_IEEE_F32LE, H5T_NATIVE_FLOAT,
      density.data(), density.size());

  Hdf5Handle tracer(H5Gcreate2(
      file.get(), "/PartType2", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
  writeCommonParticleFields(tracer.get(), member, 3U, 0U, 6000U, 1.0F);
  std::vector<std::uint64_t> parent(3U);
  std::vector<std::uint64_t> injection(3U);
  std::vector<std::uint64_t> host(3U);
  std::vector<float> fraction(3U, 0.25F);
  std::vector<float> last_mass(3U, 4.0F);
  std::vector<float> exchanged(3U, 0.0F);
  for (std::size_t index = 0U; index < parent.size(); ++index) {
    parent[index] = 8000U + member * 3U + index;
    injection[index] = 10U + index;
    host[index] = index;
  }
  writeDataset1d(
      tracer.get(), "ParentParticleIDs", H5T_STD_U64LE, H5T_NATIVE_UINT64,
      parent.data(), parent.size());
  writeDataset1d(
      tracer.get(), "InjectionStep", H5T_STD_U64LE, H5T_NATIVE_UINT64,
      injection.data(), injection.size());
  writeDataset1d(
      tracer.get(), "HostCellIndex", H5T_STD_U64LE, H5T_NATIVE_UINT64,
      host.data(), host.size());
  writeDataset1d(
      tracer.get(), "MassFractionOfHost", H5T_IEEE_F32LE, H5T_NATIVE_FLOAT,
      fraction.data(), fraction.size());
  writeDataset1d(
      tracer.get(), "LastHostMass", H5T_IEEE_F32LE, H5T_NATIVE_FLOAT,
      last_mass.data(), last_mass.size());
  writeDataset1d(
      tracer.get(), "CumulativeExchangedMass", H5T_IEEE_F32LE,
      H5T_NATIVE_FLOAT, exchanged.data(), exchanged.size());

  Hdf5Handle star(H5Gcreate2(
      file.get(), "/PartType3", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
  writeCommonParticleFields(star.get(), member, 2U, 6U, 7000U, 2.0F);
}

[[nodiscard]] std::string broadcastTestString(
    std::string value,
    int world_rank) {
  std::uint64_t size = world_rank == 0 ? value.size() : 0U;
  MPI_Bcast(&size, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
  value.resize(static_cast<std::size_t>(size));
  if (size > 0U) {
    assert(size <= static_cast<std::uint64_t>(std::numeric_limits<int>::max()));
    MPI_Bcast(
        value.data(), static_cast<int>(size), MPI_CHAR, 0, MPI_COMM_WORLD);
  }
  return value;
}

[[nodiscard]] cosmosim::core::SimulationConfig makeBridgeConfig() {
  auto config = cosmosim::core::makeUnvalidatedSimulationConfigForTests();
  config.mode.ic_convention =
      cosmosim::core::InitialConditionConvention::kGadgetArepoBridgeV1;
  config.mode.ic_bridge_source_length_unit_to_si = K_MPC_TO_SI;
  config.mode.ic_bridge_source_mass_unit_to_si = K_SOLAR_MASS_TO_SI;
  config.mode.ic_bridge_source_velocity_unit_to_si = K_KM_S_TO_SI;
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
  config.mode.ic_staging_particle_count = 21U;
  return config;
}


[[nodiscard]] cosmosim::core::SimulationConfig makeCanonicalConfig() {
  auto config = makeBridgeConfig();
  config.mode.ic_convention =
      cosmosim::core::InitialConditionConvention::kChuiCanonicalV1;
  return config;
}

[[nodiscard]] cosmosim::core::SimulationConfig makeAlternateConfig() {
  auto config = makeBridgeConfig();
  config.mode.ic_bridge_coordinate_frame =
      cosmosim::core::InitialConditionCoordinateFrame::kPhysical;
  config.mode.ic_bridge_velocity_convention =
      cosmosim::core::InitialConditionVelocityConvention::kSqrtAScaledPeculiar;
  config.mode.ic_bridge_length_hubble_exponent = -1.0;
  config.mode.ic_bridge_mass_hubble_exponent = -1.0;
  config.numerics.a_begin = 0.5;
  config.numerics.z_begin = 1.0;
  config.cosmology.hubble_param = 0.5;
  config.cosmology.box_size_x_mpc_comoving = 256.0;
  config.cosmology.box_size_y_mpc_comoving = 256.0;
  config.cosmology.box_size_z_mpc_comoving = 256.0;
  config.cosmology.box_size_mpc_comoving = 256.0;
  return config;
}

[[nodiscard]] cosmosim::core::SimulationConfig makePolicyConfig() {
  auto config = makeBridgeConfig();
  config.mode.ic_part_type2_policy =
      cosmosim::core::InitialConditionSpeciesPolicy::kTracer;
  config.mode.ic_part_type3_policy =
      cosmosim::core::InitialConditionSpeciesPolicy::kStar;
  return config;
}

[[nodiscard]] bool startsWith(std::string_view value, std::string_view prefix) {
  return value.size() >= prefix.size() &&
      value.substr(0U, prefix.size()) == prefix;
}

void assertCounterAgreesAcrossRanks(std::uint64_t local_value) {
  std::uint64_t minimum = 0U;
  std::uint64_t maximum = 0U;
  MPI_Allreduce(
      &local_value, &minimum, 1, MPI_UINT64_T, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(
      &local_value, &maximum, 1, MPI_UINT64_T, MPI_MAX, MPI_COMM_WORLD);
  assert(minimum == maximum);
}

void assertEquivalentImportedState(
    const cosmosim::core::SimulationState& direct,
    const cosmosim::core::SimulationState& supplied) {
  assert(direct.particles.position_x_comoving ==
         supplied.particles.position_x_comoving);
  assert(direct.particles.position_y_comoving ==
         supplied.particles.position_y_comoving);
  assert(direct.particles.position_z_comoving ==
         supplied.particles.position_z_comoving);
  assert(direct.particles.velocity_x_peculiar ==
         supplied.particles.velocity_x_peculiar);
  assert(direct.particles.velocity_y_peculiar ==
         supplied.particles.velocity_y_peculiar);
  assert(direct.particles.velocity_z_peculiar ==
         supplied.particles.velocity_z_peculiar);
  assert(direct.particles.mass_code == supplied.particles.mass_code);

  assert(direct.particle_sidecar.particle_id ==
         supplied.particle_sidecar.particle_id);
  assert(direct.particle_sidecar.species_tag ==
         supplied.particle_sidecar.species_tag);
  assert(direct.particle_sidecar.owning_rank ==
         supplied.particle_sidecar.owning_rank);

  assert(direct.cells.center_x_comoving == supplied.cells.center_x_comoving);
  assert(direct.cells.center_y_comoving == supplied.cells.center_y_comoving);
  assert(direct.cells.center_z_comoving == supplied.cells.center_z_comoving);
  assert(direct.cells.mass_code == supplied.cells.mass_code);

  assert(direct.gas_cells.gas_cell_id == supplied.gas_cells.gas_cell_id);
  assert(direct.gas_cells.parent_particle_id ==
         supplied.gas_cells.parent_particle_id);
  assert(direct.gas_cells.density_code == supplied.gas_cells.density_code);
  assert(direct.gas_cells.internal_energy_code ==
         supplied.gas_cells.internal_energy_code);
  assert(direct.gas_cells.temperature_code ==
         supplied.gas_cells.temperature_code);

  assert(direct.star_particles.particle_index ==
         supplied.star_particles.particle_index);
  assert(direct.star_particles.formation_scale_factor ==
         supplied.star_particles.formation_scale_factor);
  assert(direct.star_particles.birth_mass_code ==
         supplied.star_particles.birth_mass_code);
  assert(direct.star_particles.metallicity_mass_fraction ==
         supplied.star_particles.metallicity_mass_fraction);

  assert(direct.black_holes.particle_index ==
         supplied.black_holes.particle_index);
  assert(direct.black_holes.subgrid_mass_code ==
         supplied.black_holes.subgrid_mass_code);
  assert(direct.black_holes.accretion_rate_code ==
         supplied.black_holes.accretion_rate_code);
  assert(direct.black_holes.feedback_energy_code ==
         supplied.black_holes.feedback_energy_code);

  assert(direct.tracers.particle_index == supplied.tracers.particle_index);
  assert(direct.metadata.scale_factor == supplied.metadata.scale_factor);
}

void setFaultInjection(std::string_view mode) {
  if (!startsWith(mode, "fault_")) {
    return;
  }
  const std::string phase(mode.substr(6U));
  const std::string value = phase + ":1";
  assert(setenv("COSMOSIM_IC_TEST_FAULT", value.c_str(), 1) == 0);
}

}  // namespace

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  int world_rank = 0;
  int world_size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  const std::string mode = argc > 1 ? argv[1] : "normal";
  const bool duplicate_ids = mode == "duplicate";
  const bool fault_mode = startsWith(mode, "fault_");
  const bool route_mutation =
      mode == "route_loss" || mode == "route_duplicate";
  const bool canonical_manifest_mode =
      startsWith(mode, "canonical_manifest");
  const bool canonical_mode =
      mode == "canonical" || canonical_manifest_mode;
  const bool alternate_mode = mode == "alternate";
  const bool bridge_manifest_canonical_dialect =
      mode == "bridge_manifest_canonical_dialect";
  const bool invalid_manifest_mode = mode == "invalid_manifest";
  const bool manifest_mode = mode == "manifest" || canonical_manifest_mode ||
      invalid_manifest_mode || bridge_manifest_canonical_dialect;
  const bool policy_mode = mode == "type_policy";
  if ((fault_mode || route_mutation) && world_size < 2) {
    MPI_Finalize();
    return 2;
  }
  setFaultInjection(mode);
  if (mode == "route_loss") {
    assert(setenv("COSMOSIM_IC_TEST_ROUTE_MUTATION", "drop:1", 1) == 0);
  } else if (mode == "route_duplicate") {
    assert(setenv("COSMOSIM_IC_TEST_ROUTE_MUTATION", "duplicate:1", 1) == 0);
  }

  auto config = policy_mode ? makePolicyConfig()
      : alternate_mode ? makeAlternateConfig()
      : canonical_mode ? makeCanonicalConfig()
      : makeBridgeConfig();
  const cosmosim::parallel::MpiContext mpi_context(
      true, world_size, world_rank);

  if (mode == "scaling") {
    const auto scaling_root = std::filesystem::temp_directory_path();
    const auto small = runScalingImport(
        scaling_root / "cosmosim_distributed_ic_scaling_small", 32U, config,
        mpi_context);
    const auto large = runScalingImport(
        scaling_root / "cosmosim_distributed_ic_scaling_large", 128U, config,
        mpi_context);
    const std::uint64_t small_global = mpi_context.allreduceSumUint64(
        small.state.particles.size());
    const std::uint64_t large_global = mpi_context.allreduceSumUint64(
        large.state.particles.size());
    assert(small_global == 64U);
    assert(large_global == 256U);
    assert(large_global == 4U * small_global);
    const std::uint64_t small_batches = mpi_context.allreduceSumUint64(
        small.report.counters.routing_batch_count);
    const std::uint64_t large_batches = mpi_context.allreduceSumUint64(
        large.report.counters.routing_batch_count);
    const std::uint64_t small_chunks = mpi_context.allreduceSumUint64(
        small.report.counters.chunks_assigned);
    const std::uint64_t large_chunks = mpi_context.allreduceSumUint64(
        large.report.counters.chunks_assigned);
    const std::uint64_t small_main_exchanges = mpi_context.allreduceSumUint64(
        small.report.counters.main_exchange_count);
    const std::uint64_t large_main_exchanges = mpi_context.allreduceSumUint64(
        large.report.counters.main_exchange_count);
    assert(small_batches < small_chunks);
    assert(large_batches < large_chunks);
    assert(small_main_exchanges == small_batches);
    assert(large_main_exchanges == large_batches);
    assert(large_batches > small_batches);
    assert(
        large.report.counters.peak_staging_bytes <=
        small.report.counters.peak_staging_bytes + 64U * 1024U);
    if (world_size > 1) {
      assert(
          small.state.particles.position_x_comoving.capacity() <
          small_global);
      assert(
          large.state.particles.position_x_comoving.capacity() <
          large_global);
    }
    MPI_Finalize();
    return 0;
  }

  const auto base = std::filesystem::temp_directory_path() /
      ("cosmosim_distributed_ic_acceptance_" + mode);
  const auto first = std::filesystem::path(base.string() + ".0.hdf5");
  const auto second = std::filesystem::path(base.string() + ".1.hdf5");
  if (world_rank == 0) {
    if (policy_mode) {
      writePolicyMember(first, 0U);
      writePolicyMember(second, 1U);
    } else {
      writeMember(
          first, 0U, duplicate_ids, canonical_mode, alternate_mode);
      writeMember(
          second, 1U, duplicate_ids, canonical_mode, alternate_mode);
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);

  std::optional<cosmosim::io::IcManifest> supplied_manifest;
  if (manifest_mode) {
    std::string manifest_json;
    if (world_rank == 0) {
      const auto serial = cosmosim::io::readGadgetArepoHdf5Ic(
          first, config,
          cosmosim::io::IcImportOptions{.chunk_particle_count = 7U});
      assert(serial.report.manifest.has_value());
      manifest_json =
          cosmosim::io::serializeIcManifestJson(*serial.report.manifest);
    }
    manifest_json = broadcastTestString(std::move(manifest_json), world_rank);
    supplied_manifest =
        cosmosim::io::deserializeIcManifestJson(manifest_json);
    config.mode.ic_convention =
        cosmosim::core::InitialConditionConvention::kManifestV1;
    config.mode.ic_manifest_file = "in_memory_mpi_acceptance_manifest.json";

    if (invalid_manifest_mode) {
      supplied_manifest->source_sha256.clear();
    } else if (mode == "canonical_manifest_source_hash_mismatch") {
      supplied_manifest->source_sha256.front() = std::string(64U, '0');
      supplied_manifest->source_provenance_ids.front() =
          "sha256:" + supplied_manifest->source_sha256.front();
    } else if (mode == "canonical_manifest_source_path_mismatch") {
      supplied_manifest->source_files.front() =
          std::filesystem::path(base.string() + ".missing.hdf5");
    } else if (mode == "canonical_manifest_bridge_dialect") {
      supplied_manifest->dialect =
          cosmosim::io::IcDialect::kGadgetArepoBridgeV1;
      supplied_manifest->canonical_source_manifest_verified = false;
      supplied_manifest->canonical_source_manifest_sha256.clear();
    } else if (bridge_manifest_canonical_dialect) {
      supplied_manifest->dialect = cosmosim::io::IcDialect::kChuiCanonicalV1;
      supplied_manifest->canonical_source_manifest_verified = false;
      supplied_manifest->canonical_source_manifest_sha256.clear();
    }
  }

  if (world_rank == 0 && mode == "canonical_manifest_tampered_sidecar") {
    writeTextFile(first.string() + ".manifest.json", "tampered manifest");
  } else if (world_rank == 0 &&
             mode == "canonical_manifest_missing_marker") {
    std::filesystem::remove(first.string() + ".complete");
  } else if (world_rank == 0 &&
             mode == "canonical_manifest_bad_digest") {
    overwriteStringAttribute(
        first, "ConversionManifestSha256", std::string(64U, '0'));
  }
  MPI_Barrier(MPI_COMM_WORLD);

  const bool manifest_rejection_expected =
      invalid_manifest_mode ||
      mode == "canonical_manifest_tampered_sidecar" ||
      mode == "canonical_manifest_missing_marker" ||
      mode == "canonical_manifest_bad_digest" ||
      mode == "canonical_manifest_source_hash_mismatch" ||
      mode == "canonical_manifest_source_path_mismatch" ||
      mode == "canonical_manifest_bridge_dialect" ||
      bridge_manifest_canonical_dialect;
  const bool rejection_expected =
      duplicate_ids || fault_mode || route_mutation ||
      manifest_rejection_expected;
  bool rejected_as_expected = false;
  try {
    cosmosim::io::IcImportOptions import_options;
    import_options.chunk_particle_count = 7U;
    if (supplied_manifest.has_value()) {
      import_options.manifest = &*supplied_manifest;
    }
    std::optional<cosmosim::io::IcReadResult> direct_canonical_result;
    if (mode == "canonical_manifest") {
      auto direct_config = makeCanonicalConfig();
      direct_canonical_result =
          cosmosim::io::readDistributedGadgetArepoHdf5Ic(
              first, direct_config, mpi_context,
              cosmosim::io::IcImportOptions{.chunk_particle_count = 7U});
    }
    const auto result = cosmosim::io::readDistributedGadgetArepoHdf5Ic(
        first, config, mpi_context, import_options);
    assert(!rejection_expected);
    if (direct_canonical_result.has_value()) {
      assertEquivalentImportedState(
          direct_canonical_result->state, result.state);
    }
    assert(result.report.already_partitioned);
    assert(result.state.validateOwnershipInvariants());
    if (canonical_mode) {
      assert(result.report.manifest_verified);
      assert(!result.report.verified_manifest_sha256.empty());
    }
    if (canonical_manifest_mode) {
      assert(result.report.provenance_authority == "supplied_manifest_v1");
      assert(result.report.manifest.has_value());
      assert(
          result.report.manifest->dialect ==
          cosmosim::io::IcDialect::kChuiCanonicalV1);
    }

    const std::uint64_t expected_global = policy_mode
        ? 16U
        : K_GLOBAL_PARTICLE_COUNT;
    const std::uint64_t expected_chunks = policy_mode ? 6U : 14U;
    assert(
        mpi_context.allreduceSumUint64(
            result.report.counters.files_assigned) == K_MEMBER_COUNT);
    assert(
        mpi_context.allreduceSumUint64(
            result.report.counters.chunks_assigned) == expected_chunks);
    const std::uint64_t global_batches = mpi_context.allreduceSumUint64(
        result.report.counters.routing_batch_count);
    const std::uint64_t global_file_opens = mpi_context.allreduceSumUint64(
        result.report.counters.source_file_open_count);
    const std::uint64_t global_dataset_opens = mpi_context.allreduceSumUint64(
        result.report.counters.source_dataset_open_count);
    const std::uint64_t global_hash_passes = mpi_context.allreduceSumUint64(
        result.report.counters.full_file_hash_pass_count);
    const std::uint64_t global_reader_batches = mpi_context.allreduceSumUint64(
        result.report.counters.reader_batches_assigned);
    const std::uint64_t global_reader_records = mpi_context.allreduceSumUint64(
        result.report.counters.reader_records_assigned);
    const std::uint64_t global_reader_imbalance = mpi_context.allreduceSumUint64(
        result.report.counters.reader_record_imbalance);
    const std::uint64_t global_main_exchanges = mpi_context.allreduceSumUint64(
        result.report.counters.main_exchange_count);
    const std::uint64_t global_exact_audit_exchanges =
        mpi_context.allreduceSumUint64(
            result.report.counters.exact_audit_exchange_count);
    const std::uint64_t global_routing_collective_phases =
        mpi_context.allreduceSumUint64(
            result.report.counters.routing_collective_phase_count);
    assert(global_batches > 0U && global_batches <= expected_chunks);
    if (!policy_mode) assert(global_batches < expected_chunks);
    assert(global_file_opens == K_MEMBER_COUNT);
    assert(global_hash_passes <= 3U * K_MEMBER_COUNT);
    assert(global_dataset_opens > 0U);
    assert(global_reader_batches == global_batches);
    assert(global_reader_records == expected_global);
    assert(global_reader_imbalance <= expected_global);
    assert(global_main_exchanges == global_batches);
    assert(global_exact_audit_exchanges >= global_batches);
    assert(global_routing_collective_phases > global_batches);
    assert(global_routing_collective_phases <= 32U * global_batches);
    const std::uint64_t global_collective_phases =
        mpi_context.allreduceSumUint64(
            result.report.counters.collective_phase_count);
    assert(global_collective_phases > global_batches);

    if (world_size > 1) {
      const std::uint64_t local_actual_collectives =
          result.report.counters.mpi_collective_call_count;
      assertCounterAgreesAcrossRanks(local_actual_collectives);
      assertCounterAgreesAcrossRanks(
          result.report.counters.routing_mpi_collective_call_count);
      assertCounterAgreesAcrossRanks(
          result.report.counters.nonrouting_mpi_collective_call_count);
      assertCounterAgreesAcrossRanks(
          result.report.counters.mpi_allreduce_call_count);
      assertCounterAgreesAcrossRanks(
          result.report.counters.mpi_bcast_call_count);
      assertCounterAgreesAcrossRanks(
          result.report.counters.mpi_gather_call_count);
      assertCounterAgreesAcrossRanks(
          result.report.counters.mpi_gatherv_call_count);
      assertCounterAgreesAcrossRanks(
          result.report.counters.mpi_alltoall_call_count);
      assertCounterAgreesAcrossRanks(
          result.report.counters.mpi_alltoallv_call_count);
      assertCounterAgreesAcrossRanks(
          result.report.counters.logical_consensus_phase_count);
      assertCounterAgreesAcrossRanks(
          result.report.counters.routing_logical_consensus_phase_count);
      assertCounterAgreesAcrossRanks(
          result.report.counters.distributed_id_audit_round_count);
      assert(local_actual_collectives > 0U);
      assert(
          result.report.counters.routing_mpi_collective_call_count ==
          cosmosim::io::kIcRoutingMpiCollectiveCallsPerBatchV1 *
              global_batches);
      assert(
          result.report.counters.nonrouting_mpi_collective_call_count +
              result.report.counters.routing_mpi_collective_call_count ==
          local_actual_collectives);
      const std::uint64_t expected_nonrouting_collectives =
          cosmosim::io::kIcNonroutingMpiCollectiveFixedCallsV1 +
          (import_options.validate_runtime_cosmology ? 1U : 0U) +
          cosmosim::io::kIcNonroutingMpiCollectiveCallsPerSourceFileV1 *
              K_MEMBER_COUNT +
          cosmosim::io::kIcNonroutingMpiCollectiveCallsPerIdAuditRoundV1 *
              result.report.counters.distributed_id_audit_round_count +
          result.report.counters.mpi_bcast_call_count;
      assert(
          result.report.counters.nonrouting_mpi_collective_call_count ==
          expected_nonrouting_collectives);
      assert(
          local_actual_collectives ==
          cosmosim::io::kIcRoutingMpiCollectiveCallsPerBatchV1 *
                  global_batches +
              expected_nonrouting_collectives);
      assert(
          result.report.counters.mpi_allreduce_call_count +
              result.report.counters.mpi_bcast_call_count +
              result.report.counters.mpi_gather_call_count +
              result.report.counters.mpi_gatherv_call_count +
              result.report.counters.mpi_alltoall_call_count +
              result.report.counters.mpi_alltoallv_call_count ==
          local_actual_collectives);
      assert(
          result.report.counters.logical_consensus_phase_count ==
          result.report.counters.collective_phase_count);
      assert(
          result.report.counters.routing_logical_consensus_phase_count ==
          result.report.counters.routing_collective_phase_count);
      const double expected_collectives_per_million =
          static_cast<double>(local_actual_collectives) * 1.0e6 /
          static_cast<double>(expected_global);
      assert(std::abs(
          result.report.counters.collectives_per_million_records -
          expected_collectives_per_million) < 1.0e-9);
    } else {
      assert(result.report.counters.mpi_collective_call_count == 0U);
      assert(result.report.counters.collectives_per_million_records == 0.0);
    }
    assert(
        mpi_context.allreduceSumUint64(
            result.report.counters.records_read) == expected_global);
    assert(
        mpi_context.allreduceSumUint64(
            result.report.counters.records_converted) == expected_global);
    assert(
        mpi_context.allreduceSumUint64(
            result.report.counters.records_routed) == expected_global);

    const std::uint64_t local_count = result.state.particles.size();
    const std::uint64_t global_count =
        mpi_context.allreduceSumUint64(local_count);
    assert(global_count == expected_global);
    if (policy_mode) {
      assert(
          mpi_context.allreduceSumUint64(result.state.tracers.size()) == 6U);
      assert(
          mpi_context.allreduceSumUint64(
              result.state.star_particles.size()) == 4U);
      assert(
          mpi_context.allreduceSumUint64(result.state.cells.size()) == 6U);
      assert(
          mpi_context.allreduceSumUint64(result.state.black_holes.size()) ==
          0U);
    } else {
      assert(
          mpi_context.allreduceSumUint64(result.state.cells.size()) ==
          K_TOTAL_COUNTS[0]);
      assert(
          mpi_context.allreduceSumUint64(
              result.state.star_particles.size()) == K_TOTAL_COUNTS[4]);
      assert(
          mpi_context.allreduceSumUint64(result.state.black_holes.size()) ==
          K_TOTAL_COUNTS[5]);
    }

    for (std::uint32_t owner : result.state.particle_sidecar.owning_rank) {
      assert(owner == static_cast<std::uint32_t>(world_rank));
    }
    for (double position_y : result.state.particles.position_y_comoving) {
      assert(std::isfinite(position_y) && position_y >= 0.0);
    }
    for (double position_z : result.state.particles.position_z_comoving) {
      assert(std::isfinite(position_z) && position_z >= 0.0);
    }

    if (alternate_mode) {
      int local_gas_match = 0;
      int local_bh_match = 0;
      double local_error = 0.0;
      for (std::size_t row = 0U; row < result.state.gas_cells.size(); ++row) {
        if (result.state.gas_cells.parent_particle_id[row] == 1000U) {
          ++local_gas_match;
          local_error = std::max(
              local_error,
              std::abs(result.state.gas_cells.internal_energy_code[row] -
                       4.0));
          local_error = std::max(
              local_error,
              std::abs(result.state.gas_cells.density_code[row] - 0.0625));
        }
      }
      for (std::size_t row = 0U; row < result.state.black_holes.size(); ++row) {
        const std::uint32_t particle_index =
            result.state.black_holes.particle_index[row];
        if (result.state.particle_sidecar.particle_id[particle_index] ==
            5000U) {
          ++local_bh_match;
          local_error = std::max(
              local_error,
              std::abs(result.state.black_holes.subgrid_mass_code[row] -
                       10.0));
          local_error = std::max(
              local_error,
              std::abs(result.state.black_holes.accretion_rate_code[row] -
                       0.25));
        }
      }
      int global_gas_match = 0;
      int global_bh_match = 0;
      double global_error = 0.0;
      MPI_Allreduce(
          &local_gas_match, &global_gas_match, 1, MPI_INT, MPI_SUM,
          MPI_COMM_WORLD);
      MPI_Allreduce(
          &local_bh_match, &global_bh_match, 1, MPI_INT, MPI_SUM,
          MPI_COMM_WORLD);
      MPI_Allreduce(
          &local_error, &global_error, 1, MPI_DOUBLE, MPI_MAX,
          MPI_COMM_WORLD);
      assert(global_gas_match == 1);
      assert(global_bh_match == 1);
      assert(global_error < 1.0e-12);
    }

    if (world_size > 1) {
      assert(local_count < global_count);
      assert(
          result.state.particles.position_x_comoving.capacity() <
          global_count);
    }
    assert(result.report.counters.peak_staging_bytes < 512U * 1024U);
    assert(result.report.counters.hash_bytes_read > 0U || world_rank != 0);
    assert(result.report.counters.payload_bytes_read > 0U ||
           result.report.counters.chunks_assigned == 0U);
  } catch (const std::runtime_error&) {
    rejected_as_expected = rejection_expected;
  }

  if (rejection_expected) {
    assert(rejected_as_expected);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  if (world_rank == 0) {
    std::filesystem::remove(first);
    std::filesystem::remove(second);
    std::filesystem::remove(first.string() + ".manifest.json");
    std::filesystem::remove(second.string() + ".manifest.json");
    std::filesystem::remove(first.string() + ".complete");
    std::filesystem::remove(second.string() + ".complete");
  }
  unsetenv("COSMOSIM_IC_TEST_FAULT");
  unsetenv("COSMOSIM_IC_TEST_ROUTE_MUTATION");
  MPI_Finalize();
  return 0;
}

