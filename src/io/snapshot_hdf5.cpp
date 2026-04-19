#include "cosmosim/io/snapshot_hdf5.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/core/version.hpp"

#if COSMOSIM_ENABLE_HDF5
#include <hdf5.h>
#endif

namespace cosmosim::io {
namespace {

constexpr std::uint32_t k_species_dark_matter =
    static_cast<std::uint32_t>(core::ParticleSpecies::kDarkMatter);
constexpr std::uint32_t k_species_gas = static_cast<std::uint32_t>(core::ParticleSpecies::kGas);
constexpr std::uint32_t k_species_star = static_cast<std::uint32_t>(core::ParticleSpecies::kStar);
constexpr std::uint32_t k_species_black_hole =
    static_cast<std::uint32_t>(core::ParticleSpecies::kBlackHole);
constexpr std::uint32_t k_species_tracer = static_cast<std::uint32_t>(core::ParticleSpecies::kTracer);

[[nodiscard]] std::size_t mapSpeciesTagToPartType(std::uint32_t species_tag) {
  if (species_tag == k_species_gas) {
    return 0;
  }
  if (species_tag == k_species_star) {
    return 4;
  }
  if (species_tag == k_species_tracer) {
    return 3;
  }
  if (species_tag == k_species_black_hole) {
    return 5;
  }
  return 1;
}

[[nodiscard]] std::uint32_t mapPartTypeToSpeciesTag(std::size_t part_type) {
  if (part_type == 0) {
    return k_species_gas;
  }
  if (part_type == 4) {
    return k_species_star;
  }
  if (part_type == 3) {
    return k_species_tracer;
  }
  if (part_type == 5) {
    return k_species_black_hole;
  }
  return k_species_dark_matter;
}

[[nodiscard]] std::string toTypeAliasPath(std::size_t type_index) {
  return "/ParticleType" + std::to_string(type_index);
}

#if COSMOSIM_ENABLE_HDF5

class Hdf5Handle {
 public:
  explicit Hdf5Handle(hid_t handle = -1) : m_handle(handle) {}
  Hdf5Handle(const Hdf5Handle&) = delete;
  Hdf5Handle& operator=(const Hdf5Handle&) = delete;
  Hdf5Handle(Hdf5Handle&& other) noexcept : m_handle(other.m_handle) { other.m_handle = -1; }
  Hdf5Handle& operator=(Hdf5Handle&& other) noexcept {
    if (this != &other) {
      close();
      m_handle = other.m_handle;
      other.m_handle = -1;
    }
    return *this;
  }
  ~Hdf5Handle() { close(); }

  [[nodiscard]] hid_t get() const { return m_handle; }
  [[nodiscard]] bool valid() const { return m_handle >= 0; }

 private:
  void close() {
    if (m_handle >= 0) {
      const H5I_type_t type = H5Iget_type(m_handle);
      if (type == H5I_FILE) {
        H5Fclose(m_handle);
      } else if (type == H5I_GROUP) {
        H5Gclose(m_handle);
      } else if (type == H5I_DATASET) {
        H5Dclose(m_handle);
      } else if (type == H5I_DATASPACE) {
        H5Sclose(m_handle);
      } else if (type == H5I_ATTR) {
        H5Aclose(m_handle);
      } else if (type == H5I_DATATYPE) {
        H5Tclose(m_handle);
      } else if (type == H5I_GENPROP_LST) {
        H5Pclose(m_handle);
      }
      m_handle = -1;
    }
  }

  hid_t m_handle = -1;
};

[[nodiscard]] bool hdf5PathExists(hid_t parent, const std::string& path) {
  return H5Lexists(parent, path.c_str(), H5P_DEFAULT) > 0;
}

void writeScalarStringAttribute(hid_t location, const std::string& key, std::string_view value) {
  Hdf5Handle scalar_space(H5Screate(H5S_SCALAR));
  if (!scalar_space.valid()) {
    throw std::runtime_error("failed to allocate scalar dataspace for attribute: " + key);
  }

  Hdf5Handle string_type(H5Tcopy(H5T_C_S1));
  if (!string_type.valid()) {
    throw std::runtime_error("failed to allocate HDF5 string datatype");
  }
  if (H5Tset_size(string_type.get(), value.size()) < 0 || H5Tset_strpad(string_type.get(), H5T_STR_NULLTERM) < 0) {
    throw std::runtime_error("failed to configure HDF5 string datatype");
  }

  Hdf5Handle attr(H5Acreate2(location, key.c_str(), string_type.get(), scalar_space.get(), H5P_DEFAULT, H5P_DEFAULT));
  if (!attr.valid()) {
    throw std::runtime_error("failed to create attribute: " + key);
  }

  std::string local(value);
  if (H5Awrite(attr.get(), string_type.get(), local.c_str()) < 0) {
    throw std::runtime_error("failed to write attribute: " + key);
  }
}

void writeScalarDoubleAttribute(hid_t location, const std::string& key, double value) {
  Hdf5Handle scalar_space(H5Screate(H5S_SCALAR));
  Hdf5Handle attr(H5Acreate2(location, key.c_str(), H5T_IEEE_F64LE, scalar_space.get(), H5P_DEFAULT, H5P_DEFAULT));
  if (!attr.valid() || H5Awrite(attr.get(), H5T_NATIVE_DOUBLE, &value) < 0) {
    throw std::runtime_error("failed to write double attribute: " + key);
  }
}

void writeScalarUint32Attribute(hid_t location, const std::string& key, std::uint32_t value) {
  Hdf5Handle scalar_space(H5Screate(H5S_SCALAR));
  Hdf5Handle attr(H5Acreate2(location, key.c_str(), H5T_STD_U32LE, scalar_space.get(), H5P_DEFAULT, H5P_DEFAULT));
  if (!attr.valid() || H5Awrite(attr.get(), H5T_NATIVE_UINT32, &value) < 0) {
    throw std::runtime_error("failed to write uint32 attribute: " + key);
  }
}

void writeHeaderArrays(
    hid_t header_group,
    const std::array<std::uint64_t, 6>& part_count,
    const std::array<double, 6>& mass_table,
    const core::SimulationConfig& config,
    const core::SimulationState& state) {
  std::array<std::uint32_t, 6> part_count_u32{};
  std::array<std::uint32_t, 6> total_part_count_u32{};
  for (std::size_t i = 0; i < part_count.size(); ++i) {
    if (part_count[i] > static_cast<std::uint64_t>(std::numeric_limits<std::uint32_t>::max())) {
      throw std::runtime_error("snapshot writer currently supports <= uint32 particle counts per type");
    }
    part_count_u32[i] = static_cast<std::uint32_t>(part_count[i]);
    total_part_count_u32[i] = part_count_u32[i];
  }

  hsize_t vector_dims[1] = {6};
  Hdf5Handle vector_space(H5Screate_simple(1, vector_dims, nullptr));
  if (!vector_space.valid()) {
    throw std::runtime_error("failed to create vector dataspace for Header attributes");
  }

  Hdf5Handle attr_num_this(
      H5Acreate2(header_group, "NumPart_ThisFile", H5T_STD_U32LE, vector_space.get(), H5P_DEFAULT, H5P_DEFAULT));
  Hdf5Handle attr_num_total(
      H5Acreate2(header_group, "NumPart_Total", H5T_STD_U32LE, vector_space.get(), H5P_DEFAULT, H5P_DEFAULT));
  Hdf5Handle attr_num_total_high(H5Acreate2(
      header_group,
      "NumPart_Total_HighWord",
      H5T_STD_U32LE,
      vector_space.get(),
      H5P_DEFAULT,
      H5P_DEFAULT));
  Hdf5Handle attr_mass(
      H5Acreate2(header_group, "MassTable", H5T_IEEE_F64LE, vector_space.get(), H5P_DEFAULT, H5P_DEFAULT));

  const std::array<std::uint32_t, 6> zeros{};
  if (!attr_num_this.valid() || H5Awrite(attr_num_this.get(), H5T_NATIVE_UINT32, part_count_u32.data()) < 0 ||
      !attr_num_total.valid() || H5Awrite(attr_num_total.get(), H5T_NATIVE_UINT32, total_part_count_u32.data()) < 0 ||
      !attr_num_total_high.valid() || H5Awrite(attr_num_total_high.get(), H5T_NATIVE_UINT32, zeros.data()) < 0 ||
      !attr_mass.valid() || H5Awrite(attr_mass.get(), H5T_NATIVE_DOUBLE, mass_table.data()) < 0) {
    throw std::runtime_error("failed to write required Header vector attributes");
  }

  writeScalarDoubleAttribute(header_group, "Time", state.metadata.scale_factor);
  writeScalarDoubleAttribute(header_group, "Redshift", (1.0 / state.metadata.scale_factor) - 1.0);
  writeScalarDoubleAttribute(header_group, "BoxSize", config.cosmology.box_size_mpc_comoving);
  writeScalarDoubleAttribute(header_group, "Omega0", config.cosmology.omega_matter);
  writeScalarDoubleAttribute(header_group, "OmegaLambda", config.cosmology.omega_lambda);
  writeScalarDoubleAttribute(header_group, "OmegaBaryon", config.cosmology.omega_baryon);
  writeScalarDoubleAttribute(header_group, "HubbleParam", config.cosmology.hubble_param);

  writeScalarUint32Attribute(header_group, "NumFilesPerSnapshot", 1);
  writeScalarUint32Attribute(header_group, "Flag_Sfr", config.physics.enable_star_formation ? 1u : 0u);
  writeScalarUint32Attribute(header_group, "Flag_Cooling", config.physics.enable_cooling ? 1u : 0u);
  writeScalarUint32Attribute(header_group, "Flag_StellarAge", 1u);
  writeScalarUint32Attribute(header_group, "Flag_Metals", 1u);
  writeScalarUint32Attribute(header_group, "Flag_Feedback", config.physics.enable_feedback ? 1u : 0u);
}

Hdf5Handle createDatasetProperties(std::size_t element_count, std::size_t element_width, const SnapshotIoPolicy& policy) {
  Hdf5Handle properties(H5Pcreate(H5P_DATASET_CREATE));
  if (!properties.valid()) {
    throw std::runtime_error("failed to create HDF5 dataset creation properties");
  }

  const hsize_t chunk_rows = std::max<hsize_t>(1, std::min<hsize_t>(static_cast<hsize_t>(policy.chunk_particle_count), static_cast<hsize_t>(element_count)));
  if (element_width == 1) {
    hsize_t chunk_dims[1] = {chunk_rows};
    if (H5Pset_chunk(properties.get(), 1, chunk_dims) < 0) {
      throw std::runtime_error("failed to set chunking for 1D dataset");
    }
  } else {
    hsize_t chunk_dims[2] = {chunk_rows, static_cast<hsize_t>(element_width)};
    if (H5Pset_chunk(properties.get(), 2, chunk_dims) < 0) {
      throw std::runtime_error("failed to set chunking for 2D dataset");
    }
  }

  if (policy.enable_compression) {
    if (H5Pset_deflate(properties.get(), policy.compression_level) < 0) {
      throw std::runtime_error("failed to set deflate compression level");
    }
  }

  return properties;
}

void writeDataset1d(
    hid_t group,
    std::string_view name,
    const double* values,
    std::size_t count,
    const SnapshotIoPolicy& policy) {
  hsize_t dims[1] = {static_cast<hsize_t>(count)};
  Hdf5Handle dataspace(H5Screate_simple(1, dims, nullptr));
  Hdf5Handle properties = createDatasetProperties(count, 1, policy);
  Hdf5Handle dataset(
      H5Dcreate2(group, std::string(name).c_str(), H5T_IEEE_F64LE, dataspace.get(), H5P_DEFAULT, properties.get(), H5P_DEFAULT));
  if (!dataset.valid() || H5Dwrite(dataset.get(), H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, values) < 0) {
    throw std::runtime_error("failed to write dataset: " + std::string(name));
  }
}

void writeDataset1dU64(
    hid_t group,
    std::string_view name,
    const std::uint64_t* values,
    std::size_t count,
    const SnapshotIoPolicy& policy) {
  hsize_t dims[1] = {static_cast<hsize_t>(count)};
  Hdf5Handle dataspace(H5Screate_simple(1, dims, nullptr));
  Hdf5Handle properties = createDatasetProperties(count, 1, policy);
  Hdf5Handle dataset(
      H5Dcreate2(group, std::string(name).c_str(), H5T_STD_U64LE, dataspace.get(), H5P_DEFAULT, properties.get(), H5P_DEFAULT));
  if (!dataset.valid() || H5Dwrite(dataset.get(), H5T_NATIVE_UINT64, H5S_ALL, H5S_ALL, H5P_DEFAULT, values) < 0) {
    throw std::runtime_error("failed to write dataset: " + std::string(name));
  }
}

void writeDataset2d3(
    hid_t group,
    std::string_view name,
    const std::vector<double>& packed_values,
    std::size_t count,
    const SnapshotIoPolicy& policy) {
  hsize_t dims[2] = {static_cast<hsize_t>(count), 3};
  Hdf5Handle dataspace(H5Screate_simple(2, dims, nullptr));
  Hdf5Handle properties = createDatasetProperties(count, 3, policy);
  Hdf5Handle dataset(
      H5Dcreate2(group, std::string(name).c_str(), H5T_IEEE_F64LE, dataspace.get(), H5P_DEFAULT, properties.get(), H5P_DEFAULT));
  if (!dataset.valid() ||
      H5Dwrite(dataset.get(), H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, packed_values.data()) < 0) {
    throw std::runtime_error("failed to write dataset: " + std::string(name));
  }
}

[[nodiscard]] std::vector<std::uint32_t> collectGlobalIndicesForPartType(const core::SimulationState& state, std::size_t part_type) {
  std::vector<std::uint32_t> indices;
  indices.reserve(state.particles.size());
  for (std::size_t i = 0; i < state.particles.size(); ++i) {
    if (mapSpeciesTagToPartType(state.particle_sidecar.species_tag[i]) == part_type) {
      indices.push_back(static_cast<std::uint32_t>(i));
    }
  }
  return indices;
}

void packCoords(
    const core::SimulationState& state,
    const std::vector<std::uint32_t>& indices,
    std::vector<double>& out_coords,
    std::vector<double>& out_velocities,
    std::vector<double>& out_masses,
    std::vector<std::uint64_t>& out_ids) {
  out_coords.resize(indices.size() * 3);
  out_velocities.resize(indices.size() * 3);
  out_masses.resize(indices.size());
  out_ids.resize(indices.size());
  for (std::size_t i = 0; i < indices.size(); ++i) {
    const std::size_t global_i = indices[i];
    out_coords[i * 3 + 0] = state.particles.position_x_comoving[global_i];
    out_coords[i * 3 + 1] = state.particles.position_y_comoving[global_i];
    out_coords[i * 3 + 2] = state.particles.position_z_comoving[global_i];
    out_velocities[i * 3 + 0] = state.particles.velocity_x_peculiar[global_i];
    out_velocities[i * 3 + 1] = state.particles.velocity_y_peculiar[global_i];
    out_velocities[i * 3 + 2] = state.particles.velocity_z_peculiar[global_i];
    out_masses[i] = state.particles.mass_code[global_i];
    out_ids[i] = state.particle_sidecar.particle_id[global_i];
  }
}

void readDatasetChunk1d(
    hid_t group,
    const std::string& dataset_name,
    std::size_t start,
    std::size_t count,
    std::vector<double>& out) {
  Hdf5Handle dataset(H5Dopen2(group, dataset_name.c_str(), H5P_DEFAULT));
  Hdf5Handle file_space(H5Dget_space(dataset.get()));
  hsize_t file_offset[1] = {static_cast<hsize_t>(start)};
  hsize_t file_count[1] = {static_cast<hsize_t>(count)};
  H5Sselect_hyperslab(file_space.get(), H5S_SELECT_SET, file_offset, nullptr, file_count, nullptr);
  hsize_t mem_dims[1] = {static_cast<hsize_t>(count)};
  Hdf5Handle mem_space(H5Screate_simple(1, mem_dims, nullptr));
  out.resize(count);
  if (!dataset.valid() || H5Dread(dataset.get(), H5T_NATIVE_DOUBLE, mem_space.get(), file_space.get(), H5P_DEFAULT, out.data()) < 0) {
    throw std::runtime_error("failed to read 1D dataset: " + dataset_name);
  }
}

void readDatasetChunk2d(
    hid_t group,
    const std::string& dataset_name,
    std::size_t start,
    std::size_t count,
    std::vector<double>& out) {
  Hdf5Handle dataset(H5Dopen2(group, dataset_name.c_str(), H5P_DEFAULT));
  Hdf5Handle file_space(H5Dget_space(dataset.get()));
  hsize_t file_offset[2] = {static_cast<hsize_t>(start), 0};
  hsize_t file_count[2] = {static_cast<hsize_t>(count), 3};
  H5Sselect_hyperslab(file_space.get(), H5S_SELECT_SET, file_offset, nullptr, file_count, nullptr);
  hsize_t mem_dims[2] = {static_cast<hsize_t>(count), 3};
  Hdf5Handle mem_space(H5Screate_simple(2, mem_dims, nullptr));
  out.resize(count * 3);
  if (!dataset.valid() || H5Dread(dataset.get(), H5T_NATIVE_DOUBLE, mem_space.get(), file_space.get(), H5P_DEFAULT, out.data()) < 0) {
    throw std::runtime_error("failed to read 2D dataset: " + dataset_name);
  }
}

void readDatasetChunkIds(
    hid_t group,
    const std::string& dataset_name,
    std::size_t start,
    std::size_t count,
    std::vector<std::uint64_t>& out) {
  Hdf5Handle dataset(H5Dopen2(group, dataset_name.c_str(), H5P_DEFAULT));
  Hdf5Handle file_space(H5Dget_space(dataset.get()));
  hsize_t file_offset[1] = {static_cast<hsize_t>(start)};
  hsize_t file_count[1] = {static_cast<hsize_t>(count)};
  H5Sselect_hyperslab(file_space.get(), H5S_SELECT_SET, file_offset, nullptr, file_count, nullptr);
  hsize_t mem_dims[1] = {static_cast<hsize_t>(count)};
  Hdf5Handle mem_space(H5Screate_simple(1, mem_dims, nullptr));
  out.resize(count);
  if (!dataset.valid() || H5Dread(dataset.get(), H5T_NATIVE_UINT64, mem_space.get(), file_space.get(), H5P_DEFAULT, out.data()) < 0) {
    throw std::runtime_error("failed to read ids dataset: " + dataset_name);
  }
}

void readDatasetChunkU32(
    hid_t group,
    const std::string& dataset_name,
    std::size_t start,
    std::size_t count,
    std::vector<std::uint32_t>& out) {
  Hdf5Handle dataset(H5Dopen2(group, dataset_name.c_str(), H5P_DEFAULT));
  Hdf5Handle file_space(H5Dget_space(dataset.get()));
  hsize_t file_offset[1] = {static_cast<hsize_t>(start)};
  hsize_t file_count[1] = {static_cast<hsize_t>(count)};
  H5Sselect_hyperslab(file_space.get(), H5S_SELECT_SET, file_offset, nullptr, file_count, nullptr);
  hsize_t mem_dims[1] = {static_cast<hsize_t>(count)};
  Hdf5Handle mem_space(H5Screate_simple(1, mem_dims, nullptr));
  out.resize(count);
  if (!dataset.valid() || H5Dread(dataset.get(), H5T_NATIVE_UINT32, mem_space.get(), file_space.get(), H5P_DEFAULT, out.data()) < 0) {
    throw std::runtime_error("failed to read uint32 dataset: " + dataset_name);
  }
}

[[nodiscard]] std::string pickAlias(
    hid_t group,
    const GadgetArepoFieldAliases& aliases,
    SnapshotIoReport& report,
    const std::string& context_path,
    bool required) {
  for (std::string_view alias : aliases.read_aliases) {
    if (hdf5PathExists(group, std::string(alias))) {
      report.present_aliases.push_back(context_path + "=" + std::string(alias));
      return std::string(alias);
    }
  }
  if (required) {
    throw std::runtime_error("required dataset missing under: " + context_path);
  }
  return {};
}

void readScalarStringAttribute(hid_t location, const std::string& key, std::string& out_value) {
  Hdf5Handle attr(H5Aopen(location, key.c_str(), H5P_DEFAULT));
  if (!attr.valid()) {
    return;
  }
  Hdf5Handle type(H5Aget_type(attr.get()));
  const std::size_t length = H5Tget_size(type.get());
  std::string buffer(length, '\0');
  if (H5Aread(attr.get(), type.get(), buffer.data()) < 0) {
    throw std::runtime_error("failed reading attribute: " + key);
  }
  if (!buffer.empty() && buffer.back() == '\0') {
    buffer.pop_back();
  }
  out_value = buffer;
}

[[nodiscard]] bool readScalarUint32Attribute(hid_t location, const std::string& key, std::uint32_t& out_value) {
  Hdf5Handle attr(H5Aopen(location, key.c_str(), H5P_DEFAULT));
  if (!attr.valid()) {
    return false;
  }
  if (H5Aread(attr.get(), H5T_NATIVE_UINT32, &out_value) < 0) {
    throw std::runtime_error("failed reading attribute: " + key);
  }
  return true;
}

[[nodiscard]] bool readScalarDoubleAttribute(hid_t location, const std::string& key, double& out_value) {
  Hdf5Handle attr(H5Aopen(location, key.c_str(), H5P_DEFAULT));
  if (!attr.valid()) {
    return false;
  }
  if (H5Aread(attr.get(), H5T_NATIVE_DOUBLE, &out_value) < 0) {
    throw std::runtime_error("failed reading attribute: " + key);
  }
  return true;
}

void readHeaderCounts(hid_t header_group, std::array<std::uint64_t, 6>& out_counts) {
  Hdf5Handle attr_count(H5Aopen(header_group, "NumPart_ThisFile", H5P_DEFAULT));
  std::array<std::uint32_t, 6> raw_count{};
  if (!attr_count.valid() || H5Aread(attr_count.get(), H5T_NATIVE_UINT32, raw_count.data()) < 0) {
    throw std::runtime_error("Header/NumPart_ThisFile missing or unreadable");
  }
  for (std::size_t i = 0; i < out_counts.size(); ++i) {
    out_counts[i] = raw_count[i];
  }
}

void readOptionalHeaderMassTable(
    hid_t header_group,
    std::array<double, 6>& out_mass_table) {
  out_mass_table.fill(0.0);
  Hdf5Handle attr_mass(H5Aopen(header_group, "MassTable", H5P_DEFAULT));
  if (attr_mass.valid()) {
    if (H5Aread(attr_mass.get(), H5T_NATIVE_DOUBLE, out_mass_table.data()) < 0) {
      throw std::runtime_error("Header/MassTable unreadable");
    }
  }
}

void readOptionalHeaderDouble(hid_t header_group, const std::string& key, double default_value, double& out_value) {
  out_value = default_value;
  Hdf5Handle attr(H5Aopen(header_group, key.c_str(), H5P_DEFAULT));
  if (attr.valid()) {
    H5Aread(attr.get(), H5T_NATIVE_DOUBLE, &out_value);
  }
}

#endif

}  // namespace

const GadgetArepoSchemaMap& gadgetArepoSchemaMap() {
  static const GadgetArepoSchemaMap k_schema{};
  return k_schema;
}

void writeGadgetArepoSnapshotHdf5(
    const std::filesystem::path& output_path,
    const SnapshotWritePayload& payload,
    const SnapshotIoPolicy& policy) {
#if !COSMOSIM_ENABLE_HDF5
  static_cast<void>(output_path);
  static_cast<void>(payload);
  static_cast<void>(policy);
  throw std::runtime_error("COSMOSIM_ENABLE_HDF5=OFF: snapshot writer unavailable");
#else
  if (payload.state == nullptr || payload.config == nullptr) {
    throw std::runtime_error("snapshot writer requires non-null state and config");
  }
  if (policy.chunk_particle_count == 0) {
    throw std::runtime_error("chunk_particle_count must be > 0");
  }
  if (policy.enable_compression && (policy.compression_level < 0 || policy.compression_level > 9)) {
    throw std::runtime_error("compression_level must be in [0, 9]");
  }

  const core::SimulationState& state = *payload.state;
  const core::SimulationConfig& config = *payload.config;
  std::vector<std::int64_t> tracer_row_by_particle(state.particles.size(), -1);
  for (std::size_t tracer_row = 0; tracer_row < state.tracers.size(); ++tracer_row) {
    const std::uint32_t particle_index = state.tracers.particle_index[tracer_row];
    if (particle_index < tracer_row_by_particle.size()) {
      tracer_row_by_particle[particle_index] = static_cast<std::int64_t>(tracer_row);
    }
  }

  std::array<std::uint64_t, 6> count_by_type{};
  std::array<double, 6> mass_table{};

  for (std::size_t i = 0; i < state.particles.size(); ++i) {
    const std::size_t type_index = mapSpeciesTagToPartType(state.particle_sidecar.species_tag[i]);
    ++count_by_type[type_index];
  }

  const std::filesystem::path parent_dir = output_path.parent_path();
  if (!parent_dir.empty()) {
    std::filesystem::create_directories(parent_dir);
  }
  const std::filesystem::path temp_path = output_path.string() + ".tmp";

  Hdf5Handle file(H5Fcreate(temp_path.string().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT));
  if (!file.valid()) {
    throw std::runtime_error("failed creating snapshot file: " + temp_path.string());
  }

  const auto& schema = gadgetArepoSchemaMap();
  Hdf5Handle header_group(H5Gcreate2(file.get(), std::string(schema.header_group).c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
  if (!header_group.valid()) {
    throw std::runtime_error("failed creating /Header group");
  }

  writeHeaderArrays(header_group.get(), count_by_type, mass_table, config, state);
  writeScalarStringAttribute(header_group.get(), "CosmoSimSchemaName", schema.schema_name);
  writeScalarUint32Attribute(header_group.get(), "CosmoSimSchemaVersion", schema.schema_version);
  writeScalarStringAttribute(header_group.get(), "CosmoSimBuild", core::buildProvenance());

  Hdf5Handle config_group(
      H5Gcreate2(file.get(), std::string(schema.config_group).c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
  writeScalarStringAttribute(config_group.get(), std::string(schema.config_normalized_attribute), payload.normalized_config_text);

  Hdf5Handle provenance_group(
      H5Gcreate2(file.get(), std::string(schema.provenance_group).c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
  writeScalarStringAttribute(provenance_group.get(), "schema_version", payload.provenance.schema_version);
  writeScalarStringAttribute(provenance_group.get(), "git_sha", payload.provenance.git_sha.empty() ? payload.git_sha : payload.provenance.git_sha);
  writeScalarStringAttribute(provenance_group.get(), "compiler_id", payload.provenance.compiler_id);
  writeScalarStringAttribute(provenance_group.get(), "compiler_version", payload.provenance.compiler_version);
  writeScalarStringAttribute(provenance_group.get(), "build_preset", payload.provenance.build_preset);
  writeScalarStringAttribute(provenance_group.get(), "enabled_features", payload.provenance.enabled_features);
  writeScalarStringAttribute(provenance_group.get(), "config_hash_hex", payload.provenance.config_hash_hex);
  writeScalarStringAttribute(provenance_group.get(), "timestamp_utc", payload.provenance.timestamp_utc);
  writeScalarStringAttribute(provenance_group.get(), "hardware_summary", payload.provenance.hardware_summary);
  writeScalarUint32Attribute(
      provenance_group.get(),
      "gravity_treepm_pm_grid",
      static_cast<std::uint32_t>(std::max(payload.provenance.gravity_treepm_pm_grid, 0)));
  writeScalarStringAttribute(
      provenance_group.get(),
      "gravity_treepm_assignment_scheme",
      payload.provenance.gravity_treepm_assignment_scheme);
  writeScalarStringAttribute(
      provenance_group.get(),
      "gravity_treepm_window_deconvolution",
      payload.provenance.gravity_treepm_window_deconvolution ? "true" : "false");
  writeScalarDoubleAttribute(
      provenance_group.get(),
      "gravity_treepm_asmth_cells",
      payload.provenance.gravity_treepm_asmth_cells);
  writeScalarDoubleAttribute(
      provenance_group.get(),
      "gravity_treepm_rcut_cells",
      payload.provenance.gravity_treepm_rcut_cells);
  writeScalarDoubleAttribute(
      provenance_group.get(),
      "gravity_treepm_mesh_spacing_mpc_comoving",
      payload.provenance.gravity_treepm_mesh_spacing_mpc_comoving);
  writeScalarDoubleAttribute(
      provenance_group.get(),
      "gravity_treepm_split_scale_mpc_comoving",
      payload.provenance.gravity_treepm_split_scale_mpc_comoving);
  writeScalarDoubleAttribute(
      provenance_group.get(),
      "gravity_treepm_cutoff_radius_mpc_comoving",
      payload.provenance.gravity_treepm_cutoff_radius_mpc_comoving);
  writeScalarUint32Attribute(
      provenance_group.get(),
      "gravity_treepm_update_cadence_steps",
      static_cast<std::uint32_t>(std::max(payload.provenance.gravity_treepm_update_cadence_steps, 0)));
  writeScalarStringAttribute(
      provenance_group.get(),
      "gravity_softening_policy",
      payload.provenance.gravity_softening_policy);
  writeScalarStringAttribute(
      provenance_group.get(),
      "gravity_softening_kernel",
      payload.provenance.gravity_softening_kernel);
  writeScalarDoubleAttribute(
      provenance_group.get(),
      "gravity_softening_epsilon_kpc_comoving",
      payload.provenance.gravity_softening_epsilon_kpc_comoving);
  writeScalarStringAttribute(
      provenance_group.get(),
      "gravity_pm_fft_backend",
      payload.provenance.gravity_pm_fft_backend);

  for (std::size_t type_index = 0; type_index < schema.part_type_group.size(); ++type_index) {
    const std::vector<std::uint32_t> indices = collectGlobalIndicesForPartType(state, type_index);
    if (indices.empty()) {
      continue;
    }

    std::vector<double> coords;
    std::vector<double> velocities;
    std::vector<double> masses;
    std::vector<std::uint64_t> ids;
    packCoords(state, indices, coords, velocities, masses, ids);

    Hdf5Handle type_group(H5Gcreate2(
        file.get(),
        std::string(schema.part_type_group[type_index]).c_str(),
        H5P_DEFAULT,
        H5P_DEFAULT,
        H5P_DEFAULT));
    if (!type_group.valid()) {
      throw std::runtime_error("failed creating part type group");
    }

    writeDataset2d3(type_group.get(), schema.coordinates.canonical_name, coords, indices.size(), policy);
    writeDataset2d3(type_group.get(), schema.velocities.canonical_name, velocities, indices.size(), policy);
    writeDataset1d(type_group.get(), schema.masses.canonical_name, masses.data(), indices.size(), policy);
    writeDataset1dU64(type_group.get(), schema.particle_ids.canonical_name, ids.data(), indices.size(), policy);
    if (type_index == 3) {
      std::vector<std::uint64_t> parent_particle_id(indices.size(), 0);
      std::vector<std::uint64_t> injection_step(indices.size(), 0);
      std::vector<std::uint32_t> host_cell_index(indices.size(), 0);
      std::vector<double> mass_fraction_of_host(indices.size(), 0.0);
      std::vector<double> cumulative_exchanged_mass_code(indices.size(), 0.0);
      for (std::size_t i = 0; i < indices.size(); ++i) {
        const std::int64_t tracer_row = tracer_row_by_particle[indices[i]];
        if (tracer_row < 0) {
          continue;
        }
        parent_particle_id[i] = state.tracers.parent_particle_id[tracer_row];
        injection_step[i] = state.tracers.injection_step[tracer_row];
        host_cell_index[i] = state.tracers.host_cell_index[tracer_row];
        mass_fraction_of_host[i] = state.tracers.mass_fraction_of_host[tracer_row];
        cumulative_exchanged_mass_code[i] =
            state.tracers.cumulative_exchanged_mass_code[tracer_row];
      }
      writeDataset1dU64(type_group.get(), "TracerParentParticleID", parent_particle_id.data(), indices.size(), policy);
      writeDataset1dU64(type_group.get(), "TracerInjectionStep", injection_step.data(), indices.size(), policy);
      hsize_t dims[1] = {static_cast<hsize_t>(indices.size())};
      Hdf5Handle dataspace(H5Screate_simple(1, dims, nullptr));
      Hdf5Handle properties = createDatasetProperties(indices.size(), 1, policy);
      Hdf5Handle host_ds(H5Dcreate2(
          type_group.get(), "TracerHostCellIndex", H5T_STD_U32LE, dataspace.get(), H5P_DEFAULT, properties.get(), H5P_DEFAULT));
      Hdf5Handle frac_ds(H5Dcreate2(
          type_group.get(), "TracerMassFractionOfHost", H5T_IEEE_F64LE, dataspace.get(), H5P_DEFAULT, properties.get(), H5P_DEFAULT));
      Hdf5Handle exchange_ds(H5Dcreate2(
          type_group.get(), "TracerCumulativeExchangedMassCode", H5T_IEEE_F64LE, dataspace.get(), H5P_DEFAULT, properties.get(), H5P_DEFAULT));
      if (!host_ds.valid() ||
          H5Dwrite(host_ds.get(), H5T_NATIVE_UINT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, host_cell_index.data()) < 0 ||
          !frac_ds.valid() ||
          H5Dwrite(frac_ds.get(), H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, mass_fraction_of_host.data()) < 0 ||
          !exchange_ds.valid() ||
          H5Dwrite(exchange_ds.get(), H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, cumulative_exchanged_mass_code.data()) < 0) {
        throw std::runtime_error("failed to write tracer sidecar datasets");
      }
    }

    if (policy.write_particle_type_alias_groups) {
      const std::string alias_group_path = toTypeAliasPath(type_index);
      H5Lcreate_hard(file.get(), std::string(schema.part_type_group[type_index]).c_str(), file.get(), alias_group_path.c_str(), H5P_DEFAULT, H5P_DEFAULT);
    }
  }

  if (H5Fflush(file.get(), H5F_SCOPE_GLOBAL) < 0) {
    throw std::runtime_error("failed to flush snapshot file");
  }

  if (std::filesystem::exists(output_path)) {
    std::filesystem::remove(output_path);
  }
  std::filesystem::rename(temp_path, output_path);
#endif
}

SnapshotReadResult readGadgetArepoSnapshotHdf5(
    const std::filesystem::path& input_path,
    const core::SimulationConfig& config,
    const SnapshotReadOptions& options) {
#if !COSMOSIM_ENABLE_HDF5
  static_cast<void>(input_path);
  static_cast<void>(config);
  static_cast<void>(options);
  throw std::runtime_error("COSMOSIM_ENABLE_HDF5=OFF: snapshot reader unavailable");
#else
  SnapshotReadResult result;
  Hdf5Handle file(H5Fopen(input_path.string().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
  if (!file.valid()) {
    throw std::runtime_error("failed opening snapshot file: " + input_path.string());
  }

  const auto& schema = gadgetArepoSchemaMap();
  Hdf5Handle header_group(H5Gopen2(file.get(), std::string(schema.header_group).c_str(), H5P_DEFAULT));
  if (!header_group.valid()) {
    throw std::runtime_error("snapshot missing " + std::string(schema.header_group));
  }

  readOptionalHeaderDouble(header_group.get(), "Time", 1.0, result.state.metadata.scale_factor);
  result.report.schema_name = std::string(schema.schema_name);
  readScalarStringAttribute(header_group.get(), "CosmoSimSchemaName", result.report.schema_name);

  std::uint32_t schema_version = schema.schema_version;
  {
    Hdf5Handle attr(H5Aopen(header_group.get(), "CosmoSimSchemaVersion", H5P_DEFAULT));
    if (attr.valid()) {
      H5Aread(attr.get(), H5T_NATIVE_UINT32, &schema_version);
    }
  }
  result.report.schema_version = schema_version;

  std::array<std::uint64_t, 6> header_counts{};
  std::array<double, 6> header_mass_table{};
  readHeaderCounts(header_group.get(), header_counts);
  readOptionalHeaderMassTable(header_group.get(), header_mass_table);
  std::size_t total_count = 0;
  for (std::uint64_t count : header_counts) {
    total_count += static_cast<std::size_t>(count);
  }
  result.state.resizeParticles(total_count);
  std::vector<std::uint32_t> tracer_particle_index;
  std::vector<std::uint64_t> tracer_parent_particle_id;
  std::vector<std::uint64_t> tracer_injection_step;
  std::vector<std::uint32_t> tracer_host_cell_index;
  std::vector<double> tracer_mass_fraction_of_host;
  std::vector<double> tracer_cumulative_exchanged_mass_code;

  std::size_t global_offset = 0;
  for (std::size_t type_index = 0; type_index < header_counts.size(); ++type_index) {
    const std::size_t local_count = static_cast<std::size_t>(header_counts[type_index]);
    if (local_count == 0) {
      continue;
    }

    Hdf5Handle group(H5Gopen2(file.get(), std::string(schema.part_type_group[type_index]).c_str(), H5P_DEFAULT));
    if (!group.valid()) {
      const std::string alias_group = toTypeAliasPath(type_index);
      group = Hdf5Handle(H5Gopen2(file.get(), alias_group.c_str(), H5P_DEFAULT));
      if (!group.valid()) {
        throw std::runtime_error("missing particle group for type " + std::to_string(type_index));
      }
      result.report.present_aliases.push_back(std::string(schema.part_type_group[type_index]) + "=" + alias_group);
    }

    const std::string coordinates_name = pickAlias(
        group.get(),
        schema.coordinates,
        result.report,
        std::string(schema.part_type_group[type_index]) + "/Coordinates",
        true);
    const std::string velocities_name = pickAlias(
        group.get(),
        schema.velocities,
        result.report,
        std::string(schema.part_type_group[type_index]) + "/Velocities",
        options.require_velocities);
    const std::string ids_name = pickAlias(
        group.get(),
        schema.particle_ids,
        result.report,
        std::string(schema.part_type_group[type_index]) + "/ParticleIDs",
        options.require_ids);
    const std::string masses_name = pickAlias(
        group.get(),
        schema.masses,
        result.report,
        std::string(schema.part_type_group[type_index]) + "/Masses",
        false);

    std::vector<double> coords_chunk;
    std::vector<double> vel_chunk;
    std::vector<double> mass_chunk;
    std::vector<std::uint64_t> ids_chunk;
    std::vector<std::uint64_t> tracer_parent_chunk;
    std::vector<std::uint64_t> tracer_step_chunk;
    std::vector<std::uint32_t> tracer_host_chunk;
    std::vector<double> tracer_fraction_chunk;
    std::vector<double> tracer_exchange_chunk;

    readDatasetChunk2d(group.get(), coordinates_name, 0, local_count, coords_chunk);
    if (!velocities_name.empty()) {
      readDatasetChunk2d(group.get(), velocities_name, 0, local_count, vel_chunk);
    } else {
      vel_chunk.assign(local_count * 3, 0.0);
      result.report.defaulted_fields.push_back(std::string(schema.part_type_group[type_index]) + "/Velocities=zero");
    }

    if (!ids_name.empty()) {
      readDatasetChunkIds(group.get(), ids_name, 0, local_count, ids_chunk);
    } else {
      ids_chunk.resize(local_count);
      for (std::size_t i = 0; i < local_count; ++i) {
        ids_chunk[i] = static_cast<std::uint64_t>(global_offset + i + 1);
      }
      result.report.defaulted_fields.push_back(std::string(schema.part_type_group[type_index]) + "/ParticleIDs=generated");
    }

    if (!masses_name.empty()) {
      readDatasetChunk1d(group.get(), masses_name, 0, local_count, mass_chunk);
    } else if (options.allow_mass_table_fallback) {
      const double constant_mass = header_mass_table[type_index];
      if (!(constant_mass > 0.0)) {
        throw std::runtime_error(
            "Masses missing and Header/MassTable fallback is unavailable for type " +
            std::to_string(type_index));
      }
      mass_chunk.assign(local_count, constant_mass);
      result.report.defaulted_fields.push_back(
          std::string(schema.part_type_group[type_index]) + "/Masses=MassTable");
    } else {
      throw std::runtime_error("Masses missing and fallback disabled");
    }
    if (type_index == 3) {
      if (hdf5PathExists(group.get(), "TracerParentParticleID")) {
        readDatasetChunkIds(group.get(), "TracerParentParticleID", 0, local_count, tracer_parent_chunk);
      } else {
        tracer_parent_chunk.assign(local_count, 0);
      }
      if (hdf5PathExists(group.get(), "TracerInjectionStep")) {
        readDatasetChunkIds(group.get(), "TracerInjectionStep", 0, local_count, tracer_step_chunk);
      } else {
        tracer_step_chunk.assign(local_count, 0);
      }
      if (hdf5PathExists(group.get(), "TracerHostCellIndex")) {
        readDatasetChunkU32(group.get(), "TracerHostCellIndex", 0, local_count, tracer_host_chunk);
      } else {
        tracer_host_chunk.assign(local_count, 0);
      }
      if (hdf5PathExists(group.get(), "TracerMassFractionOfHost")) {
        readDatasetChunk1d(group.get(), "TracerMassFractionOfHost", 0, local_count, tracer_fraction_chunk);
      } else {
        tracer_fraction_chunk.assign(local_count, 0.0);
      }
      if (hdf5PathExists(group.get(), "TracerCumulativeExchangedMassCode")) {
        readDatasetChunk1d(group.get(), "TracerCumulativeExchangedMassCode", 0, local_count, tracer_exchange_chunk);
      } else {
        tracer_exchange_chunk.assign(local_count, 0.0);
      }
    }

    for (std::size_t i = 0; i < local_count; ++i) {
      const std::size_t global_i = global_offset + i;
      result.state.particles.position_x_comoving[global_i] = coords_chunk[i * 3 + 0];
      result.state.particles.position_y_comoving[global_i] = coords_chunk[i * 3 + 1];
      result.state.particles.position_z_comoving[global_i] = coords_chunk[i * 3 + 2];
      result.state.particles.velocity_x_peculiar[global_i] = vel_chunk[i * 3 + 0];
      result.state.particles.velocity_y_peculiar[global_i] = vel_chunk[i * 3 + 1];
      result.state.particles.velocity_z_peculiar[global_i] = vel_chunk[i * 3 + 2];
      result.state.particles.mass_code[global_i] = mass_chunk[i];
      result.state.particles.time_bin[global_i] = 0;
      result.state.particle_sidecar.particle_id[global_i] = ids_chunk[i];
      result.state.particle_sidecar.species_tag[global_i] = mapPartTypeToSpeciesTag(type_index);
      result.state.particle_sidecar.owning_rank[global_i] = 0;
      result.state.species.count_by_species[result.state.particle_sidecar.species_tag[global_i]] += 1;
      if (type_index == 3) {
        tracer_particle_index.push_back(static_cast<std::uint32_t>(global_i));
        tracer_parent_particle_id.push_back(tracer_parent_chunk[i]);
        tracer_injection_step.push_back(tracer_step_chunk[i]);
        tracer_host_cell_index.push_back(tracer_host_chunk[i]);
        tracer_mass_fraction_of_host.push_back(tracer_fraction_chunk[i]);
        tracer_cumulative_exchanged_mass_code.push_back(tracer_exchange_chunk[i]);
      }
    }

    global_offset += local_count;
  }

  result.state.metadata.run_name = config.output.run_name;
  result.state.tracers.resize(tracer_particle_index.size());
  for (std::size_t i = 0; i < tracer_particle_index.size(); ++i) {
    result.state.tracers.particle_index[i] = tracer_particle_index[i];
    result.state.tracers.parent_particle_id[i] = tracer_parent_particle_id[i];
    result.state.tracers.injection_step[i] = tracer_injection_step[i];
    result.state.tracers.host_cell_index[i] = tracer_host_cell_index[i];
    result.state.tracers.mass_fraction_of_host[i] = tracer_mass_fraction_of_host[i];
    result.state.tracers.last_host_mass_code[i] = 0.0;
    result.state.tracers.cumulative_exchanged_mass_code[i] = tracer_cumulative_exchanged_mass_code[i];
  }
  result.state.rebuildSpeciesIndex();

  Hdf5Handle config_group(H5Gopen2(file.get(), std::string(schema.config_group).c_str(), H5P_DEFAULT));
  if (config_group.valid()) {
    readScalarStringAttribute(
        config_group.get(),
        std::string(schema.config_normalized_attribute),
        result.normalized_config_text);
  }

  Hdf5Handle provenance_group(H5Gopen2(file.get(), std::string(schema.provenance_group).c_str(), H5P_DEFAULT));
  if (provenance_group.valid()) {
    readScalarStringAttribute(provenance_group.get(), "schema_version", result.provenance.schema_version);
    readScalarStringAttribute(provenance_group.get(), "git_sha", result.provenance.git_sha);
    readScalarStringAttribute(provenance_group.get(), "compiler_id", result.provenance.compiler_id);
    readScalarStringAttribute(provenance_group.get(), "compiler_version", result.provenance.compiler_version);
    readScalarStringAttribute(provenance_group.get(), "build_preset", result.provenance.build_preset);
    readScalarStringAttribute(provenance_group.get(), "enabled_features", result.provenance.enabled_features);
    readScalarStringAttribute(provenance_group.get(), "config_hash_hex", result.provenance.config_hash_hex);
    readScalarStringAttribute(provenance_group.get(), "timestamp_utc", result.provenance.timestamp_utc);
    readScalarStringAttribute(provenance_group.get(), "hardware_summary", result.provenance.hardware_summary);
    std::uint32_t gravity_pm_grid = 0;
    if (readScalarUint32Attribute(provenance_group.get(), "gravity_treepm_pm_grid", gravity_pm_grid)) {
      result.provenance.gravity_treepm_pm_grid = static_cast<int>(gravity_pm_grid);
    }
    readScalarStringAttribute(
        provenance_group.get(),
        "gravity_treepm_assignment_scheme",
        result.provenance.gravity_treepm_assignment_scheme);
    std::string gravity_deconvolution;
    readScalarStringAttribute(
        provenance_group.get(), "gravity_treepm_window_deconvolution", gravity_deconvolution);
    if (!gravity_deconvolution.empty()) {
      result.provenance.gravity_treepm_window_deconvolution = gravity_deconvolution == "true";
    }
    static_cast<void>(readScalarDoubleAttribute(
        provenance_group.get(),
        "gravity_treepm_asmth_cells",
        result.provenance.gravity_treepm_asmth_cells));
    static_cast<void>(readScalarDoubleAttribute(
        provenance_group.get(),
        "gravity_treepm_rcut_cells",
        result.provenance.gravity_treepm_rcut_cells));
    static_cast<void>(readScalarDoubleAttribute(
        provenance_group.get(),
        "gravity_treepm_mesh_spacing_mpc_comoving",
        result.provenance.gravity_treepm_mesh_spacing_mpc_comoving));
    static_cast<void>(readScalarDoubleAttribute(
        provenance_group.get(),
        "gravity_treepm_split_scale_mpc_comoving",
        result.provenance.gravity_treepm_split_scale_mpc_comoving));
    static_cast<void>(readScalarDoubleAttribute(
        provenance_group.get(),
        "gravity_treepm_cutoff_radius_mpc_comoving",
        result.provenance.gravity_treepm_cutoff_radius_mpc_comoving));
    std::uint32_t gravity_pm_cadence = 0;
    if (readScalarUint32Attribute(
            provenance_group.get(), "gravity_treepm_update_cadence_steps", gravity_pm_cadence)) {
      result.provenance.gravity_treepm_update_cadence_steps = static_cast<int>(gravity_pm_cadence);
    }
    readScalarStringAttribute(
        provenance_group.get(),
        "gravity_softening_policy",
        result.provenance.gravity_softening_policy);
    readScalarStringAttribute(
        provenance_group.get(),
        "gravity_softening_kernel",
        result.provenance.gravity_softening_kernel);
    static_cast<void>(readScalarDoubleAttribute(
        provenance_group.get(),
        "gravity_softening_epsilon_kpc_comoving",
        result.provenance.gravity_softening_epsilon_kpc_comoving));
    readScalarStringAttribute(
        provenance_group.get(),
        "gravity_pm_fft_backend",
        result.provenance.gravity_pm_fft_backend);
  }

  result.report.chunk_particle_count = 0;
  result.report.compression_enabled = false;
  result.report.compression_level = 0;
  return result;
#endif
}

}  // namespace cosmosim::io
