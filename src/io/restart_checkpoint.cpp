#include "cosmosim/io/restart_checkpoint.hpp"

#include <algorithm>
#include <bit>
#include <array>
#include <cerrno>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>
#include <span>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/io/io_contract.hpp"

#if COSMOSIM_ENABLE_HDF5
#include <fcntl.h>
#include <hdf5.h>
#include <unistd.h>
#endif

namespace cosmosim::io {
namespace {

constexpr std::uint64_t k_offset_basis = 14695981039346656037ull;
constexpr std::uint64_t k_prime = 1099511628211ull;

[[nodiscard]] std::string hexU64(std::uint64_t value) {
  constexpr char k_hex[] = "0123456789abcdef";
  std::string text(16, '0');
  for (int i = 15; i >= 0; --i) {
    text[static_cast<std::size_t>(i)] = k_hex[value & 0xF];
    value >>= 4;
  }
  return text;
}

[[nodiscard]] std::uint64_t fnv1aAppend(std::uint64_t seed, std::span<const std::byte> bytes) {
  std::uint64_t hash = seed;
  for (const std::byte value : bytes) {
    hash ^= static_cast<std::uint64_t>(std::to_integer<unsigned char>(value));
    hash *= k_prime;
  }
  return hash;
}

template <typename VectorLike>
[[nodiscard]] std::span<const std::byte> asBytesSpan(const VectorLike& values) {
  using ValueType = typename VectorLike::value_type;
  return {reinterpret_cast<const std::byte*>(values.data()), values.size() * sizeof(ValueType)};
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
    if (m_handle < 0) {
      return;
    }
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
    }
    m_handle = -1;
  }

  hid_t m_handle = -1;
};

void writeScalarStringAttribute(hid_t location, std::string_view key, const std::string& value) {
  Hdf5Handle scalar_space(H5Screate(H5S_SCALAR));
  Hdf5Handle string_type(H5Tcopy(H5T_C_S1));
  if (!scalar_space.valid() || !string_type.valid()) {
    throw std::runtime_error("failed to create string attribute type");
  }
  if (H5Tset_size(string_type.get(), value.size()) < 0 || H5Tset_strpad(string_type.get(), H5T_STR_NULLTERM) < 0) {
    throw std::runtime_error("failed to configure string attribute type");
  }
  Hdf5Handle attr(H5Acreate2(location, std::string(key).c_str(), string_type.get(), scalar_space.get(), H5P_DEFAULT, H5P_DEFAULT));
  if (!attr.valid() || H5Awrite(attr.get(), string_type.get(), value.c_str()) < 0) {
    throw std::runtime_error("failed to write string attribute: " + std::string(key));
  }
}

[[nodiscard]] std::string readScalarStringAttribute(hid_t location, std::string_view key) {
  Hdf5Handle attr(H5Aopen(location, std::string(key).c_str(), H5P_DEFAULT));
  if (!attr.valid()) {
    throw std::runtime_error("missing string attribute: " + std::string(key));
  }
  Hdf5Handle attr_type(H5Aget_type(attr.get()));
  const std::size_t length = H5Tget_size(attr_type.get());
  std::string value(length, '\0');
  if (H5Aread(attr.get(), attr_type.get(), value.data()) < 0) {
    throw std::runtime_error("failed reading string attribute: " + std::string(key));
  }
  while (!value.empty() && value.back() == '\0') {
    value.pop_back();
  }
  return value;
}

void writeScalarU32Attribute(hid_t location, std::string_view key, std::uint32_t value) {
  Hdf5Handle scalar_space(H5Screate(H5S_SCALAR));
  Hdf5Handle attr(H5Acreate2(location, std::string(key).c_str(), H5T_STD_U32LE, scalar_space.get(), H5P_DEFAULT, H5P_DEFAULT));
  if (!attr.valid() || H5Awrite(attr.get(), H5T_NATIVE_UINT32, &value) < 0) {
    throw std::runtime_error("failed writing u32 attribute: " + std::string(key));
  }
}

[[nodiscard]] std::uint32_t readScalarU32Attribute(hid_t location, std::string_view key) {
  Hdf5Handle attr(H5Aopen(location, std::string(key).c_str(), H5P_DEFAULT));
  std::uint32_t value = 0;
  if (!attr.valid() || H5Aread(attr.get(), H5T_NATIVE_UINT32, &value) < 0) {
    throw std::runtime_error("failed reading u32 attribute: " + std::string(key));
  }
  return value;
}

void writeScalarU64Attribute(hid_t location, std::string_view key, std::uint64_t value) {
  Hdf5Handle scalar_space(H5Screate(H5S_SCALAR));
  Hdf5Handle attr(H5Acreate2(location, std::string(key).c_str(), H5T_STD_U64LE, scalar_space.get(), H5P_DEFAULT, H5P_DEFAULT));
  if (!attr.valid() || H5Awrite(attr.get(), H5T_NATIVE_UINT64, &value) < 0) {
    throw std::runtime_error("failed writing u64 attribute: " + std::string(key));
  }
}

[[nodiscard]] std::uint64_t readScalarU64Attribute(hid_t location, std::string_view key) {
  Hdf5Handle attr(H5Aopen(location, std::string(key).c_str(), H5P_DEFAULT));
  std::uint64_t value = 0;
  if (!attr.valid() || H5Aread(attr.get(), H5T_NATIVE_UINT64, &value) < 0) {
    throw std::runtime_error("failed reading u64 attribute: " + std::string(key));
  }
  return value;
}

void writeScalarF64Attribute(hid_t location, std::string_view key, double value) {
  Hdf5Handle scalar_space(H5Screate(H5S_SCALAR));
  Hdf5Handle attr(H5Acreate2(location, std::string(key).c_str(), H5T_IEEE_F64LE, scalar_space.get(), H5P_DEFAULT, H5P_DEFAULT));
  if (!attr.valid() || H5Awrite(attr.get(), H5T_NATIVE_DOUBLE, &value) < 0) {
    throw std::runtime_error("failed writing f64 attribute: " + std::string(key));
  }
}

[[nodiscard]] double readScalarF64Attribute(hid_t location, std::string_view key) {
  Hdf5Handle attr(H5Aopen(location, std::string(key).c_str(), H5P_DEFAULT));
  double value = 0.0;
  if (!attr.valid() || H5Aread(attr.get(), H5T_NATIVE_DOUBLE, &value) < 0) {
    throw std::runtime_error("failed reading f64 attribute: " + std::string(key));
  }
  return value;
}

template <typename T>
void writeDataset1d(
    hid_t group,
    std::string_view name,
    hid_t file_type,
    hid_t memory_type,
    std::span<const T> values) {
  hsize_t dims[1] = {static_cast<hsize_t>(values.size())};
  Hdf5Handle space(H5Screate_simple(1, dims, nullptr));
  Hdf5Handle dataset(H5Dcreate2(group, std::string(name).c_str(), file_type, space.get(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
  if (!dataset.valid()) {
    throw std::runtime_error("failed creating dataset: " + std::string(name));
  }
  if (!values.empty() && H5Dwrite(dataset.get(), memory_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, values.data()) < 0) {
    throw std::runtime_error("failed writing dataset: " + std::string(name));
  }
}

template <typename T, typename Allocator>
void writeDataset1d(
    hid_t group,
    std::string_view name,
    hid_t file_type,
    hid_t memory_type,
    const std::vector<T, Allocator>& values) {
  writeDataset1d<T>(group, name, file_type, memory_type, std::span<const T>(values.data(), values.size()));
}

template <typename T>
[[nodiscard]] std::vector<T> readDataset1d(hid_t group, std::string_view name, hid_t memory_type) {
  Hdf5Handle dataset(H5Dopen2(group, std::string(name).c_str(), H5P_DEFAULT));
  if (!dataset.valid()) {
    throw std::runtime_error("missing dataset: " + std::string(name));
  }
  Hdf5Handle space(H5Dget_space(dataset.get()));
  hsize_t dims[1] = {0};
  if (H5Sget_simple_extent_ndims(space.get()) != 1 || H5Sget_simple_extent_dims(space.get(), dims, nullptr) != 1) {
    throw std::runtime_error("unexpected rank for dataset: " + std::string(name));
  }
  std::vector<T> values(static_cast<std::size_t>(dims[0]));
  if (!values.empty() && H5Dread(dataset.get(), memory_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, values.data()) < 0) {
    throw std::runtime_error("failed reading dataset: " + std::string(name));
  }
  return values;
}

template <typename T>
[[nodiscard]] core::AlignedVector<T> toAlignedVector(const std::vector<T>& values) {
  core::AlignedVector<T> aligned(values.size());
  std::copy(values.begin(), values.end(), aligned.begin());
  return aligned;
}

template <typename T>
[[nodiscard]] core::AlignedVector<T> readDataset1dAligned(hid_t group, std::string_view name, hid_t memory_type) {
  return toAlignedVector(readDataset1d<T>(group, name, memory_type));
}

void writeStringDataset(hid_t group, std::string_view name, const std::string& value) {
  const std::vector<std::uint8_t> bytes(value.begin(), value.end());
  writeDataset1d(group, name, H5T_STD_U8LE, H5T_NATIVE_UINT8, bytes);
}

[[nodiscard]] std::string readStringDataset(hid_t group, std::string_view name) {
  const auto bytes = readDataset1d<std::uint8_t>(group, name, H5T_NATIVE_UINT8);
  return std::string(bytes.begin(), bytes.end());
}

[[nodiscard]] hid_t openOrCreateGroup(hid_t parent, const std::string& name) {
  if (H5Lexists(parent, name.c_str(), H5P_DEFAULT) > 0) {
    return H5Gopen2(parent, name.c_str(), H5P_DEFAULT);
  }
  return H5Gcreate2(parent, name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
}

void maybeFsync(const std::filesystem::path& file_path, bool enabled) {
  if (!enabled) {
    return;
  }
#if defined(_WIN32)
  (void)file_path;
  throw std::runtime_error("restart fsync finalize is not implemented on this platform");
#else
  const int fd = ::open(file_path.c_str(), O_RDONLY);
  if (fd < 0) {
    throw std::runtime_error(
        "failed to open restart temporary file for fsync: " + file_path.string() + ": " + std::strerror(errno));
  }
  if (::fsync(fd) != 0) {
    const std::string message = std::strerror(errno);
    ::close(fd);
    throw std::runtime_error("failed to fsync restart temporary file: " + file_path.string() + ": " + message);
  }
  if (::close(fd) != 0) {
    throw std::runtime_error(
        "failed to close restart temporary file after fsync: " + file_path.string() + ": " + std::strerror(errno));
  }
#endif
}

void writeStateGroup(hid_t root, const core::SimulationState& state) {
  Hdf5Handle state_group(openOrCreateGroup(root, "/state"));
  Hdf5Handle particles_group(openOrCreateGroup(state_group.get(), "particles"));
  writeDataset1d(particles_group.get(), "position_x_comoving", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.particles.position_x_comoving);
  writeDataset1d(particles_group.get(), "position_y_comoving", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.particles.position_y_comoving);
  writeDataset1d(particles_group.get(), "position_z_comoving", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.particles.position_z_comoving);
  writeDataset1d(particles_group.get(), "velocity_x_peculiar", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.particles.velocity_x_peculiar);
  writeDataset1d(particles_group.get(), "velocity_y_peculiar", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.particles.velocity_y_peculiar);
  writeDataset1d(particles_group.get(), "velocity_z_peculiar", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.particles.velocity_z_peculiar);
  writeDataset1d(particles_group.get(), "mass_code", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.particles.mass_code);
  writeDataset1d(particles_group.get(), "time_bin", H5T_STD_U8LE, H5T_NATIVE_UINT8, state.particles.time_bin);

  Hdf5Handle particle_sidecar_group(openOrCreateGroup(state_group.get(), "particle_sidecar"));
  writeDataset1d(particle_sidecar_group.get(), "particle_id", H5T_STD_U64LE, H5T_NATIVE_UINT64, state.particle_sidecar.particle_id);
  writeDataset1d(particle_sidecar_group.get(), "sfc_key", H5T_STD_U64LE, H5T_NATIVE_UINT64, state.particle_sidecar.sfc_key);
  writeDataset1d(particle_sidecar_group.get(), "species_tag", H5T_STD_U32LE, H5T_NATIVE_UINT32, state.particle_sidecar.species_tag);
  writeDataset1d(particle_sidecar_group.get(), "particle_flags", H5T_STD_U32LE, H5T_NATIVE_UINT32, state.particle_sidecar.particle_flags);
  writeDataset1d(particle_sidecar_group.get(), "owning_rank", H5T_STD_U32LE, H5T_NATIVE_UINT32, state.particle_sidecar.owning_rank);

  Hdf5Handle cells_group(openOrCreateGroup(state_group.get(), "cells"));
  writeDataset1d(cells_group.get(), "center_x_comoving", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.cells.center_x_comoving);
  writeDataset1d(cells_group.get(), "center_y_comoving", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.cells.center_y_comoving);
  writeDataset1d(cells_group.get(), "center_z_comoving", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.cells.center_z_comoving);
  writeDataset1d(cells_group.get(), "mass_code", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.cells.mass_code);
  writeDataset1d(cells_group.get(), "time_bin", H5T_STD_U8LE, H5T_NATIVE_UINT8, state.cells.time_bin);
  writeDataset1d(cells_group.get(), "patch_index", H5T_STD_U32LE, H5T_NATIVE_UINT32, state.cells.patch_index);

  Hdf5Handle gas_group(openOrCreateGroup(state_group.get(), "gas_cells"));
  writeDataset1d(gas_group.get(), "density_code", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.gas_cells.density_code);
  writeDataset1d(gas_group.get(), "pressure_code", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.gas_cells.pressure_code);
  writeDataset1d(gas_group.get(), "internal_energy_code", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.gas_cells.internal_energy_code);
  writeDataset1d(gas_group.get(), "temperature_code", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.gas_cells.temperature_code);
  writeDataset1d(gas_group.get(), "sound_speed_code", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.gas_cells.sound_speed_code);
  writeDataset1d(gas_group.get(), "recon_gradient_x", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.gas_cells.recon_gradient_x);
  writeDataset1d(gas_group.get(), "recon_gradient_y", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.gas_cells.recon_gradient_y);
  writeDataset1d(gas_group.get(), "recon_gradient_z", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.gas_cells.recon_gradient_z);

  Hdf5Handle patches_group(openOrCreateGroup(state_group.get(), "patches"));
  writeDataset1d(patches_group.get(), "patch_id", H5T_STD_U64LE, H5T_NATIVE_UINT64, state.patches.patch_id);
  writeDataset1d(patches_group.get(), "level", H5T_STD_I32LE, H5T_NATIVE_INT32, state.patches.level);
  writeDataset1d(patches_group.get(), "first_cell", H5T_STD_U32LE, H5T_NATIVE_UINT32, state.patches.first_cell);
  writeDataset1d(patches_group.get(), "cell_count", H5T_STD_U32LE, H5T_NATIVE_UINT32, state.patches.cell_count);

  writeDataset1d(state_group.get(), "species_count_by_species", H5T_STD_U64LE, H5T_NATIVE_UINT64, std::vector<std::uint64_t>(state.species.count_by_species.begin(), state.species.count_by_species.end()));

  Hdf5Handle star_group(openOrCreateGroup(state_group.get(), "star_particles"));
  writeDataset1d(star_group.get(), "particle_index", H5T_STD_U32LE, H5T_NATIVE_UINT32, state.star_particles.particle_index);
  writeDataset1d(star_group.get(), "formation_scale_factor", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.star_particles.formation_scale_factor);
  writeDataset1d(star_group.get(), "birth_mass_code", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.star_particles.birth_mass_code);
  writeDataset1d(star_group.get(), "metallicity_mass_fraction", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.star_particles.metallicity_mass_fraction);

  Hdf5Handle bh_group(openOrCreateGroup(state_group.get(), "black_holes"));
  writeDataset1d(bh_group.get(), "particle_index", H5T_STD_U32LE, H5T_NATIVE_UINT32, state.black_holes.particle_index);
  writeDataset1d(bh_group.get(), "host_cell_index", H5T_STD_U32LE, H5T_NATIVE_UINT32, state.black_holes.host_cell_index);
  writeDataset1d(bh_group.get(), "subgrid_mass_code", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.black_holes.subgrid_mass_code);
  writeDataset1d(bh_group.get(), "accretion_rate_code", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.black_holes.accretion_rate_code);
  writeDataset1d(bh_group.get(), "feedback_energy_code", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.black_holes.feedback_energy_code);
  writeDataset1d(bh_group.get(), "eddington_ratio", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.black_holes.eddington_ratio);
  writeDataset1d(
      bh_group.get(),
      "cumulative_accreted_mass_code",
      H5T_IEEE_F64LE,
      H5T_NATIVE_DOUBLE,
      state.black_holes.cumulative_accreted_mass_code);
  writeDataset1d(
      bh_group.get(),
      "cumulative_feedback_energy_code",
      H5T_IEEE_F64LE,
      H5T_NATIVE_DOUBLE,
      state.black_holes.cumulative_feedback_energy_code);
  writeDataset1d(
      bh_group.get(),
      "duty_cycle_active_time_code",
      H5T_IEEE_F64LE,
      H5T_NATIVE_DOUBLE,
      state.black_holes.duty_cycle_active_time_code);
  writeDataset1d(
      bh_group.get(),
      "duty_cycle_total_time_code",
      H5T_IEEE_F64LE,
      H5T_NATIVE_DOUBLE,
      state.black_holes.duty_cycle_total_time_code);

  Hdf5Handle tracer_group(openOrCreateGroup(state_group.get(), "tracers"));
  writeDataset1d(tracer_group.get(), "particle_index", H5T_STD_U32LE, H5T_NATIVE_UINT32, state.tracers.particle_index);
  writeDataset1d(tracer_group.get(), "parent_particle_id", H5T_STD_U64LE, H5T_NATIVE_UINT64, state.tracers.parent_particle_id);
  writeDataset1d(tracer_group.get(), "injection_step", H5T_STD_U64LE, H5T_NATIVE_UINT64, state.tracers.injection_step);
  writeDataset1d(tracer_group.get(), "host_cell_index", H5T_STD_U32LE, H5T_NATIVE_UINT32, state.tracers.host_cell_index);
  writeDataset1d(tracer_group.get(), "mass_fraction_of_host", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.tracers.mass_fraction_of_host);
  writeDataset1d(tracer_group.get(), "last_host_mass_code", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.tracers.last_host_mass_code);
  writeDataset1d(
      tracer_group.get(),
      "cumulative_exchanged_mass_code",
      H5T_IEEE_F64LE,
      H5T_NATIVE_DOUBLE,
      state.tracers.cumulative_exchanged_mass_code);

  writeStringDataset(state_group.get(), "state_metadata", state.metadata.serialize());

  Hdf5Handle sidecars_group(openOrCreateGroup(state_group.get(), "module_sidecars"));
  const auto sidecar_blocks = state.sidecars.blocksSortedByName();
  std::string module_names_text;
  for (const core::ModuleSidecarBlock* block : sidecar_blocks) {
    Hdf5Handle module_group(openOrCreateGroup(sidecars_group.get(), block->module_name));
    writeScalarU32Attribute(module_group.get(), "schema_version", block->schema_version);
    std::vector<std::uint8_t> payload(block->payload.size());
    for (std::size_t i = 0; i < block->payload.size(); ++i) {
      payload[i] = std::to_integer<std::uint8_t>(block->payload[i]);
    }
    writeDataset1d(module_group.get(), "payload", H5T_STD_U8LE, H5T_NATIVE_UINT8, payload);
    module_names_text += block->module_name;
    module_names_text.push_back('\n');
  }
  writeStringDataset(state_group.get(), "module_sidecar_names", module_names_text);
}

void readStateGroup(hid_t root, core::SimulationState& state) {
  Hdf5Handle state_group(H5Gopen2(root, "/state", H5P_DEFAULT));
  Hdf5Handle particles_group(H5Gopen2(state_group.get(), "particles", H5P_DEFAULT));
  state.particles.position_x_comoving = readDataset1dAligned<double>(particles_group.get(), "position_x_comoving", H5T_NATIVE_DOUBLE);
  state.particles.position_y_comoving = readDataset1dAligned<double>(particles_group.get(), "position_y_comoving", H5T_NATIVE_DOUBLE);
  state.particles.position_z_comoving = readDataset1dAligned<double>(particles_group.get(), "position_z_comoving", H5T_NATIVE_DOUBLE);
  state.particles.velocity_x_peculiar = readDataset1dAligned<double>(particles_group.get(), "velocity_x_peculiar", H5T_NATIVE_DOUBLE);
  state.particles.velocity_y_peculiar = readDataset1dAligned<double>(particles_group.get(), "velocity_y_peculiar", H5T_NATIVE_DOUBLE);
  state.particles.velocity_z_peculiar = readDataset1dAligned<double>(particles_group.get(), "velocity_z_peculiar", H5T_NATIVE_DOUBLE);
  state.particles.mass_code = readDataset1dAligned<double>(particles_group.get(), "mass_code", H5T_NATIVE_DOUBLE);
  state.particles.time_bin = readDataset1dAligned<std::uint8_t>(particles_group.get(), "time_bin", H5T_NATIVE_UINT8);

  Hdf5Handle particle_sidecar_group(H5Gopen2(state_group.get(), "particle_sidecar", H5P_DEFAULT));
  state.particle_sidecar.particle_id =
      readDataset1dAligned<std::uint64_t>(particle_sidecar_group.get(), "particle_id", H5T_NATIVE_UINT64);
  state.particle_sidecar.sfc_key =
      readDataset1dAligned<std::uint64_t>(particle_sidecar_group.get(), "sfc_key", H5T_NATIVE_UINT64);
  state.particle_sidecar.species_tag =
      readDataset1dAligned<std::uint32_t>(particle_sidecar_group.get(), "species_tag", H5T_NATIVE_UINT32);
  state.particle_sidecar.particle_flags =
      readDataset1dAligned<std::uint32_t>(particle_sidecar_group.get(), "particle_flags", H5T_NATIVE_UINT32);
  state.particle_sidecar.owning_rank =
      readDataset1dAligned<std::uint32_t>(particle_sidecar_group.get(), "owning_rank", H5T_NATIVE_UINT32);

  Hdf5Handle cells_group(H5Gopen2(state_group.get(), "cells", H5P_DEFAULT));
  state.cells.center_x_comoving = readDataset1dAligned<double>(cells_group.get(), "center_x_comoving", H5T_NATIVE_DOUBLE);
  state.cells.center_y_comoving = readDataset1dAligned<double>(cells_group.get(), "center_y_comoving", H5T_NATIVE_DOUBLE);
  state.cells.center_z_comoving = readDataset1dAligned<double>(cells_group.get(), "center_z_comoving", H5T_NATIVE_DOUBLE);
  state.cells.mass_code = readDataset1dAligned<double>(cells_group.get(), "mass_code", H5T_NATIVE_DOUBLE);
  state.cells.time_bin = readDataset1dAligned<std::uint8_t>(cells_group.get(), "time_bin", H5T_NATIVE_UINT8);
  state.cells.patch_index = readDataset1dAligned<std::uint32_t>(cells_group.get(), "patch_index", H5T_NATIVE_UINT32);

  Hdf5Handle gas_group(H5Gopen2(state_group.get(), "gas_cells", H5P_DEFAULT));
  state.gas_cells.density_code = readDataset1dAligned<double>(gas_group.get(), "density_code", H5T_NATIVE_DOUBLE);
  state.gas_cells.pressure_code = readDataset1dAligned<double>(gas_group.get(), "pressure_code", H5T_NATIVE_DOUBLE);
  state.gas_cells.internal_energy_code = readDataset1dAligned<double>(gas_group.get(), "internal_energy_code", H5T_NATIVE_DOUBLE);
  state.gas_cells.temperature_code = readDataset1dAligned<double>(gas_group.get(), "temperature_code", H5T_NATIVE_DOUBLE);
  state.gas_cells.sound_speed_code = readDataset1dAligned<double>(gas_group.get(), "sound_speed_code", H5T_NATIVE_DOUBLE);
  state.gas_cells.recon_gradient_x = readDataset1dAligned<double>(gas_group.get(), "recon_gradient_x", H5T_NATIVE_DOUBLE);
  state.gas_cells.recon_gradient_y = readDataset1dAligned<double>(gas_group.get(), "recon_gradient_y", H5T_NATIVE_DOUBLE);
  state.gas_cells.recon_gradient_z = readDataset1dAligned<double>(gas_group.get(), "recon_gradient_z", H5T_NATIVE_DOUBLE);

  Hdf5Handle patches_group(H5Gopen2(state_group.get(), "patches", H5P_DEFAULT));
  state.patches.patch_id = readDataset1dAligned<std::uint64_t>(patches_group.get(), "patch_id", H5T_NATIVE_UINT64);
  state.patches.level = readDataset1dAligned<std::int32_t>(patches_group.get(), "level", H5T_NATIVE_INT32);
  state.patches.first_cell = readDataset1dAligned<std::uint32_t>(patches_group.get(), "first_cell", H5T_NATIVE_UINT32);
  state.patches.cell_count = readDataset1dAligned<std::uint32_t>(patches_group.get(), "cell_count", H5T_NATIVE_UINT32);

  const auto species_count = readDataset1d<std::uint64_t>(state_group.get(), "species_count_by_species", H5T_NATIVE_UINT64);
  if (species_count.size() != state.species.count_by_species.size()) {
    throw std::runtime_error("invalid species_count_by_species size in restart");
  }
  std::copy(species_count.begin(), species_count.end(), state.species.count_by_species.begin());

  Hdf5Handle star_group(H5Gopen2(state_group.get(), "star_particles", H5P_DEFAULT));
  state.star_particles.particle_index =
      readDataset1dAligned<std::uint32_t>(star_group.get(), "particle_index", H5T_NATIVE_UINT32);
  state.star_particles.formation_scale_factor =
      readDataset1dAligned<double>(star_group.get(), "formation_scale_factor", H5T_NATIVE_DOUBLE);
  state.star_particles.birth_mass_code = readDataset1dAligned<double>(star_group.get(), "birth_mass_code", H5T_NATIVE_DOUBLE);
  state.star_particles.metallicity_mass_fraction =
      readDataset1dAligned<double>(star_group.get(), "metallicity_mass_fraction", H5T_NATIVE_DOUBLE);

  Hdf5Handle bh_group(H5Gopen2(state_group.get(), "black_holes", H5P_DEFAULT));
  state.black_holes.particle_index = readDataset1dAligned<std::uint32_t>(bh_group.get(), "particle_index", H5T_NATIVE_UINT32);
  state.black_holes.host_cell_index = readDataset1dAligned<std::uint32_t>(bh_group.get(), "host_cell_index", H5T_NATIVE_UINT32);
  state.black_holes.subgrid_mass_code = readDataset1dAligned<double>(bh_group.get(), "subgrid_mass_code", H5T_NATIVE_DOUBLE);
  state.black_holes.accretion_rate_code =
      readDataset1dAligned<double>(bh_group.get(), "accretion_rate_code", H5T_NATIVE_DOUBLE);
  state.black_holes.feedback_energy_code =
      readDataset1dAligned<double>(bh_group.get(), "feedback_energy_code", H5T_NATIVE_DOUBLE);
  state.black_holes.eddington_ratio = readDataset1dAligned<double>(bh_group.get(), "eddington_ratio", H5T_NATIVE_DOUBLE);
  state.black_holes.cumulative_accreted_mass_code =
      readDataset1dAligned<double>(bh_group.get(), "cumulative_accreted_mass_code", H5T_NATIVE_DOUBLE);
  state.black_holes.cumulative_feedback_energy_code =
      readDataset1dAligned<double>(bh_group.get(), "cumulative_feedback_energy_code", H5T_NATIVE_DOUBLE);
  state.black_holes.duty_cycle_active_time_code =
      readDataset1dAligned<double>(bh_group.get(), "duty_cycle_active_time_code", H5T_NATIVE_DOUBLE);
  state.black_holes.duty_cycle_total_time_code =
      readDataset1dAligned<double>(bh_group.get(), "duty_cycle_total_time_code", H5T_NATIVE_DOUBLE);

  Hdf5Handle tracer_group(H5Gopen2(state_group.get(), "tracers", H5P_DEFAULT));
  state.tracers.particle_index = readDataset1dAligned<std::uint32_t>(tracer_group.get(), "particle_index", H5T_NATIVE_UINT32);
  state.tracers.parent_particle_id =
      readDataset1dAligned<std::uint64_t>(tracer_group.get(), "parent_particle_id", H5T_NATIVE_UINT64);
  state.tracers.injection_step =
      readDataset1dAligned<std::uint64_t>(tracer_group.get(), "injection_step", H5T_NATIVE_UINT64);
  state.tracers.host_cell_index =
      readDataset1dAligned<std::uint32_t>(tracer_group.get(), "host_cell_index", H5T_NATIVE_UINT32);
  state.tracers.mass_fraction_of_host =
      readDataset1dAligned<double>(tracer_group.get(), "mass_fraction_of_host", H5T_NATIVE_DOUBLE);
  state.tracers.last_host_mass_code =
      readDataset1dAligned<double>(tracer_group.get(), "last_host_mass_code", H5T_NATIVE_DOUBLE);
  state.tracers.cumulative_exchanged_mass_code =
      readDataset1dAligned<double>(tracer_group.get(), "cumulative_exchanged_mass_code", H5T_NATIVE_DOUBLE);

  state.metadata = core::StateMetadata::deserialize(readStringDataset(state_group.get(), "state_metadata"));
  state.rebuildSpeciesIndex();

  Hdf5Handle sidecars_group(H5Gopen2(state_group.get(), "module_sidecars", H5P_DEFAULT));
  const std::string module_names_text = readStringDataset(state_group.get(), "module_sidecar_names");
  std::istringstream module_names_stream(module_names_text);
  std::string module_name;
  while (std::getline(module_names_stream, module_name)) {
    if (module_name.empty()) {
      continue;
    }
    Hdf5Handle module_group(H5Gopen2(sidecars_group.get(), module_name.c_str(), H5P_DEFAULT));
    core::ModuleSidecarBlock block;
    block.module_name = module_name;
    block.schema_version = readScalarU32Attribute(module_group.get(), "schema_version");
    const auto payload_u8 = readDataset1d<std::uint8_t>(module_group.get(), "payload", H5T_NATIVE_UINT8);
    block.payload.resize(payload_u8.size());
    for (std::size_t i = 0; i < payload_u8.size(); ++i) {
      block.payload[i] = static_cast<std::byte>(payload_u8[i]);
    }
    state.sidecars.upsert(std::move(block));
  }

  if (!state.validateOwnershipInvariants()) {
    throw std::runtime_error("restart state failed ownership invariant validation");
  }
}
#endif

}  // namespace

const RestartSchema& restartSchema() {
  static const RestartSchema schema{};
  return schema;
}

bool isRestartSchemaCompatible(std::uint32_t file_schema_version) {
  return file_schema_version == restartSchema().version;
}

const std::vector<std::string_view>& exactRestartCompletenessChecklist() {
  static const std::vector<std::string_view> checklist = {
      "simulation_state_lanes_and_metadata",
      "module_sidecars_with_schema_versions",
      "integrator_state",
      "scheduler_persistent_state",
      "normalized_config_text_and_hash",
      "provenance_record",
      "payload_integrity_hash_and_hex"};
  return checklist;
}

std::uint64_t restartPayloadIntegrityHash(const RestartWritePayload& payload) {
  if (payload.state == nullptr || payload.integrator_state == nullptr || payload.scheduler == nullptr) {
    throw std::invalid_argument("restart payload must provide state, integrator_state, and scheduler");
  }
  validateContinuationMetadata(
      payload.normalized_config_text,
      payload.normalized_config_hash_hex,
      payload.provenance,
      "restart payload");

  std::uint64_t hash = k_offset_basis;

  const auto append_string = [&hash](const std::string& value) {
    hash = fnv1aAppend(hash, {reinterpret_cast<const std::byte*>(value.data()), value.size()});
  };

  const auto append_u64 = [&hash](std::uint64_t value) {
    const std::array<std::byte, sizeof(std::uint64_t)> bytes =
        std::bit_cast<std::array<std::byte, sizeof(std::uint64_t)>>(value);
    hash = fnv1aAppend(hash, bytes);
  };

  append_string(payload.normalized_config_text);
  append_string(payload.normalized_config_hash_hex);
  append_string(core::serializeProvenanceRecord(payload.provenance));
  append_string(payload.state->metadata.serialize());

  const auto append_any_vec = [&hash](const auto& values) { hash = fnv1aAppend(hash, asBytesSpan(values)); };

  const core::SimulationState& state = *payload.state;
  append_any_vec(state.particles.position_x_comoving);
  append_any_vec(state.particles.position_y_comoving);
  append_any_vec(state.particles.position_z_comoving);
  append_any_vec(state.particles.velocity_x_peculiar);
  append_any_vec(state.particles.velocity_y_peculiar);
  append_any_vec(state.particles.velocity_z_peculiar);
  append_any_vec(state.particles.mass_code);
  append_any_vec(state.particles.time_bin);
  append_any_vec(state.particle_sidecar.particle_id);
  append_any_vec(state.particle_sidecar.sfc_key);
  append_any_vec(state.particle_sidecar.species_tag);
  append_any_vec(state.particle_sidecar.particle_flags);
  append_any_vec(state.particle_sidecar.owning_rank);

  append_any_vec(state.cells.center_x_comoving);
  append_any_vec(state.cells.center_y_comoving);
  append_any_vec(state.cells.center_z_comoving);
  append_any_vec(state.cells.mass_code);
  append_any_vec(state.cells.time_bin);
  append_any_vec(state.cells.patch_index);

  append_any_vec(state.gas_cells.density_code);
  append_any_vec(state.gas_cells.pressure_code);
  append_any_vec(state.gas_cells.internal_energy_code);
  append_any_vec(state.gas_cells.temperature_code);
  append_any_vec(state.gas_cells.sound_speed_code);
  append_any_vec(state.gas_cells.recon_gradient_x);
  append_any_vec(state.gas_cells.recon_gradient_y);
  append_any_vec(state.gas_cells.recon_gradient_z);

  append_any_vec(state.patches.patch_id);
  append_any_vec(state.patches.level);
  append_any_vec(state.patches.first_cell);
  append_any_vec(state.patches.cell_count);
  append_any_vec(state.star_particles.particle_index);
  append_any_vec(state.star_particles.formation_scale_factor);
  append_any_vec(state.star_particles.birth_mass_code);
  append_any_vec(state.star_particles.metallicity_mass_fraction);
  append_any_vec(state.black_holes.particle_index);
  append_any_vec(state.black_holes.host_cell_index);
  append_any_vec(state.black_holes.subgrid_mass_code);
  append_any_vec(state.black_holes.accretion_rate_code);
  append_any_vec(state.black_holes.feedback_energy_code);
  append_any_vec(state.black_holes.eddington_ratio);
  append_any_vec(state.black_holes.cumulative_accreted_mass_code);
  append_any_vec(state.black_holes.cumulative_feedback_energy_code);
  append_any_vec(state.black_holes.duty_cycle_active_time_code);
  append_any_vec(state.black_holes.duty_cycle_total_time_code);
  append_any_vec(state.tracers.particle_index);
  append_any_vec(state.tracers.parent_particle_id);
  append_any_vec(state.tracers.injection_step);
  append_any_vec(state.tracers.host_cell_index);
  append_any_vec(state.tracers.mass_fraction_of_host);
  append_any_vec(state.tracers.last_host_mass_code);
  append_any_vec(state.tracers.cumulative_exchanged_mass_code);
  append_any_vec(state.species.count_by_species);

  const auto ordered_sidecars = state.sidecars.blocksSortedByName();
  for (const core::ModuleSidecarBlock* block : ordered_sidecars) {
    append_string(block->module_name);
    append_u64(static_cast<std::uint64_t>(block->schema_version));
    hash = fnv1aAppend(hash, std::span<const std::byte>(block->payload.data(), block->payload.size()));
  }

  append_u64(payload.integrator_state->step_index);
  append_u64(std::bit_cast<std::uint64_t>(payload.integrator_state->current_time_code));
  append_u64(std::bit_cast<std::uint64_t>(payload.integrator_state->current_scale_factor));
  append_u64(std::bit_cast<std::uint64_t>(payload.integrator_state->dt_time_code));
  append_u64(static_cast<std::uint64_t>(payload.integrator_state->scheme));
  append_u64(payload.integrator_state->time_bins.hierarchical_enabled ? 1ull : 0ull);
  append_u64(static_cast<std::uint64_t>(payload.integrator_state->time_bins.active_bin));
  append_u64(static_cast<std::uint64_t>(payload.integrator_state->time_bins.max_bin));

  const core::TimeBinPersistentState scheduler_state = payload.scheduler->exportPersistentState();
  append_u64(scheduler_state.current_tick);
  append_u64(static_cast<std::uint64_t>(scheduler_state.max_bin));
  append_any_vec(scheduler_state.bin_index);
  append_any_vec(scheduler_state.next_activation_tick);
  append_any_vec(scheduler_state.active_flag);
  append_any_vec(scheduler_state.pending_bin_index);

  return hash;
}

std::string restartPayloadIntegrityHashHex(const RestartWritePayload& payload) {
  return hexU64(restartPayloadIntegrityHash(payload));
}

void writeRestartCheckpointHdf5(
    const std::filesystem::path& output_path,
    const RestartWritePayload& payload,
    const RestartWritePolicy& policy) {
#if !COSMOSIM_ENABLE_HDF5
  (void)output_path;
  (void)payload;
  (void)policy;
  throw std::runtime_error("restart checkpoint requires COSMOSIM_ENABLE_HDF5=ON");
#else
  if (payload.state == nullptr || payload.integrator_state == nullptr || payload.scheduler == nullptr) {
    throw std::invalid_argument("restart write payload must include state, integrator_state, and scheduler");
  }
  validateContinuationMetadata(
      payload.normalized_config_text,
      payload.normalized_config_hash_hex,
      payload.provenance,
      "restart writer");
  if (!payload.state->validateOwnershipInvariants()) {
    throw std::invalid_argument("cannot checkpoint invalid simulation state");
  }

  std::filesystem::create_directories(output_path.parent_path());
  const std::filesystem::path temporary_path = output_path.string() + policy.temporary_suffix;

  Hdf5Handle file(H5Fcreate(temporary_path.string().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT));
  if (!file.valid()) {
    throw std::runtime_error("failed to create temporary restart file: " + temporary_path.string());
  }

  writeScalarStringAttribute(file.get(), "restart_schema_name", restartSchema().name);
  writeScalarU32Attribute(file.get(), "restart_schema_version", restartSchema().version);
  writeScalarStringAttribute(file.get(), "normalized_config_hash_hex", payload.normalized_config_hash_hex);
  writeScalarStringAttribute(file.get(), "payload_integrity_hash_hex", restartPayloadIntegrityHashHex(payload));
  writeScalarU64Attribute(file.get(), "payload_integrity_hash", restartPayloadIntegrityHash(payload));

  const auto& shared_names = sharedIoContractNames();
  writeStringDataset(
      file.get(),
      std::string(shared_names.normalized_config_text_dataset),
      payload.normalized_config_text);
  writeStringDataset(
      file.get(),
      std::string(shared_names.provenance_record_dataset),
      core::serializeProvenanceRecord(payload.provenance));

  writeStateGroup(file.get(), *payload.state);

  Hdf5Handle integrator_group(openOrCreateGroup(file.get(), "/integrator"));
  writeScalarF64Attribute(integrator_group.get(), "current_time_code", payload.integrator_state->current_time_code);
  writeScalarF64Attribute(integrator_group.get(), "current_scale_factor", payload.integrator_state->current_scale_factor);
  writeScalarF64Attribute(integrator_group.get(), "dt_time_code", payload.integrator_state->dt_time_code);
  writeScalarU64Attribute(integrator_group.get(), "step_index", payload.integrator_state->step_index);
  writeScalarU32Attribute(integrator_group.get(), "scheme", static_cast<std::uint32_t>(payload.integrator_state->scheme));
  writeScalarU32Attribute(integrator_group.get(), "time_bins_hierarchical", payload.integrator_state->time_bins.hierarchical_enabled ? 1U : 0U);
  writeScalarU32Attribute(integrator_group.get(), "time_bins_active_bin", payload.integrator_state->time_bins.active_bin);
  writeScalarU32Attribute(integrator_group.get(), "time_bins_max_bin", payload.integrator_state->time_bins.max_bin);

  const core::TimeBinPersistentState scheduler_state = payload.scheduler->exportPersistentState();
  Hdf5Handle scheduler_group(openOrCreateGroup(file.get(), "/scheduler"));
  writeScalarU64Attribute(scheduler_group.get(), "current_tick", scheduler_state.current_tick);
  writeScalarU32Attribute(scheduler_group.get(), "max_bin", scheduler_state.max_bin);
  writeDataset1d(scheduler_group.get(), "bin_index", H5T_STD_U8LE, H5T_NATIVE_UINT8, scheduler_state.bin_index);
  writeDataset1d(
      scheduler_group.get(),
      "next_activation_tick",
      H5T_STD_U64LE,
      H5T_NATIVE_UINT64,
      scheduler_state.next_activation_tick);
  writeDataset1d(scheduler_group.get(), "active_flag", H5T_STD_U8LE, H5T_NATIVE_UINT8, scheduler_state.active_flag);
  writeDataset1d(
      scheduler_group.get(),
      "pending_bin_index",
      H5T_STD_U8LE,
      H5T_NATIVE_UINT8,
      scheduler_state.pending_bin_index);

  if (H5Fflush(file.get(), H5F_SCOPE_GLOBAL) < 0) {
    throw std::runtime_error("failed to flush temporary restart file");
  }
  file = Hdf5Handle();

  maybeFsync(temporary_path, policy.enable_fsync_finalize);

  std::error_code rename_error;
  std::filesystem::rename(temporary_path, output_path, rename_error);
  if (rename_error) {
    throw std::runtime_error(
        "failed to atomically finalize restart checkpoint from '" + temporary_path.string() + "' to '" +
        output_path.string() + "': " + rename_error.message());
  }
#endif
}

RestartReadResult readRestartCheckpointHdf5(const std::filesystem::path& input_path) {
#if !COSMOSIM_ENABLE_HDF5
  (void)input_path;
  throw std::runtime_error("restart checkpoint requires COSMOSIM_ENABLE_HDF5=ON");
#else
  Hdf5Handle file(H5Fopen(input_path.string().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
  if (!file.valid()) {
    throw std::runtime_error("failed to open restart checkpoint: " + input_path.string());
  }

  RestartReadResult result;
  const std::string schema_name = readScalarStringAttribute(file.get(), "restart_schema_name");
  const std::uint32_t schema_version = readScalarU32Attribute(file.get(), "restart_schema_version");
  if (schema_name != restartSchema().name || !isRestartSchemaCompatible(schema_version)) {
    throw std::runtime_error(
        "restart schema is not compatible: file='" + schema_name + "' v" + std::to_string(schema_version) +
        ", expected='" + restartSchema().name + "' v" + std::to_string(restartSchema().version));
  }

  result.normalized_config_hash_hex = readScalarStringAttribute(file.get(), "normalized_config_hash_hex");
  result.payload_hash_hex = readScalarStringAttribute(file.get(), "payload_integrity_hash_hex");
  result.payload_hash = readScalarU64Attribute(file.get(), "payload_integrity_hash");

  const auto& shared_names = sharedIoContractNames();
  result.normalized_config_text =
      readStringDataset(file.get(), std::string(shared_names.normalized_config_text_dataset));
  result.provenance = core::deserializeProvenanceRecord(
      readStringDataset(file.get(), std::string(shared_names.provenance_record_dataset)));

  readStateGroup(file.get(), result.state);

  Hdf5Handle integrator_group(H5Gopen2(file.get(), "/integrator", H5P_DEFAULT));
  result.integrator_state.current_time_code = readScalarF64Attribute(integrator_group.get(), "current_time_code");
  result.integrator_state.current_scale_factor = readScalarF64Attribute(integrator_group.get(), "current_scale_factor");
  result.integrator_state.dt_time_code = readScalarF64Attribute(integrator_group.get(), "dt_time_code");
  result.integrator_state.step_index = readScalarU64Attribute(integrator_group.get(), "step_index");
  result.integrator_state.scheme =
      static_cast<core::TimeStepScheme>(readScalarU32Attribute(integrator_group.get(), "scheme"));
  result.integrator_state.time_bins.hierarchical_enabled =
      readScalarU32Attribute(integrator_group.get(), "time_bins_hierarchical") != 0;
  result.integrator_state.time_bins.active_bin =
      static_cast<std::uint8_t>(readScalarU32Attribute(integrator_group.get(), "time_bins_active_bin"));
  result.integrator_state.time_bins.max_bin =
      static_cast<std::uint8_t>(readScalarU32Attribute(integrator_group.get(), "time_bins_max_bin"));

  Hdf5Handle scheduler_group(H5Gopen2(file.get(), "/scheduler", H5P_DEFAULT));
  result.scheduler_state.current_tick = readScalarU64Attribute(scheduler_group.get(), "current_tick");
  result.scheduler_state.max_bin = static_cast<std::uint8_t>(readScalarU32Attribute(scheduler_group.get(), "max_bin"));
  result.scheduler_state.bin_index = readDataset1d<std::uint8_t>(scheduler_group.get(), "bin_index", H5T_NATIVE_UINT8);
  result.scheduler_state.next_activation_tick =
      readDataset1d<std::uint64_t>(scheduler_group.get(), "next_activation_tick", H5T_NATIVE_UINT64);
  result.scheduler_state.active_flag = readDataset1d<std::uint8_t>(scheduler_group.get(), "active_flag", H5T_NATIVE_UINT8);
  result.scheduler_state.pending_bin_index =
      readDataset1d<std::uint8_t>(scheduler_group.get(), "pending_bin_index", H5T_NATIVE_UINT8);

  RestartWritePayload verify_payload;
  verify_payload.state = &result.state;
  verify_payload.integrator_state = &result.integrator_state;
  core::HierarchicalTimeBinScheduler verify_scheduler(result.scheduler_state.max_bin);
  verify_scheduler.importPersistentState(result.scheduler_state);
  verify_payload.scheduler = &verify_scheduler;
  verify_payload.normalized_config_hash_hex = result.normalized_config_hash_hex;
  verify_payload.normalized_config_text = result.normalized_config_text;
  verify_payload.provenance = result.provenance;

  const std::uint64_t computed_hash = restartPayloadIntegrityHash(verify_payload);
  if (computed_hash != result.payload_hash || hexU64(computed_hash) != result.payload_hash_hex) {
    throw std::runtime_error("restart payload integrity hash mismatch");
  }

  return result;
#endif
}

}  // namespace cosmosim::io
