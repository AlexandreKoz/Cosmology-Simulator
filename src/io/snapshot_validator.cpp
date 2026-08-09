#include "cosmosim/io/snapshot_hdf5.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <limits>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/core/checked_arithmetic.hpp"
#include "internal/snapshot_field_contract.hpp"

#if COSMOSIM_ENABLE_HDF5
#include <hdf5.h>
#endif

namespace cosmosim::io {
namespace {

#if COSMOSIM_ENABLE_HDF5
class H5Handle {
 public:
  explicit H5Handle(hid_t value = -1) : m_value(value) {}
  H5Handle(const H5Handle&) = delete;
  H5Handle& operator=(const H5Handle&) = delete;
  H5Handle(H5Handle&& other) noexcept : m_value(other.m_value) { other.m_value = -1; }
  H5Handle& operator=(H5Handle&& other) noexcept {
    if (this != &other) {
      close();
      m_value = other.m_value;
      other.m_value = -1;
    }
    return *this;
  }
  ~H5Handle() { close(); }
  [[nodiscard]] hid_t get() const noexcept { return m_value; }
  [[nodiscard]] bool valid() const noexcept { return m_value >= 0; }
 private:
  void close() noexcept {
    if (m_value < 0) return;
    switch (H5Iget_type(m_value)) {
      case H5I_FILE: H5Fclose(m_value); break;
      case H5I_GROUP: H5Gclose(m_value); break;
      case H5I_DATASET: H5Dclose(m_value); break;
      case H5I_DATASPACE: H5Sclose(m_value); break;
      case H5I_DATATYPE: H5Tclose(m_value); break;
      case H5I_ATTR: H5Aclose(m_value); break;
      default: break;
    }
    m_value = -1;
  }
  hid_t m_value = -1;
};

[[nodiscard]] std::uint64_t checkedAdd(std::uint64_t lhs, std::uint64_t rhs, std::string_view label) {
  if (rhs > std::numeric_limits<std::uint64_t>::max() - lhs) {
    throw std::overflow_error("snapshot validator: overflow while accounting " + std::string(label));
  }
  return lhs + rhs;
}

[[nodiscard]] std::uint64_t checkedMul(std::uint64_t lhs, std::uint64_t rhs, std::string_view label) {
  if (lhs != 0U && rhs > std::numeric_limits<std::uint64_t>::max() / lhs) {
    throw std::overflow_error("snapshot validator: overflow while accounting " + std::string(label));
  }
  return lhs * rhs;
}

[[nodiscard]] std::array<std::uint64_t, 6> readCountArray(hid_t header, const char* name) {
  H5Handle attr(H5Aopen(header, name, H5P_DEFAULT));
  if (!attr.valid()) throw std::runtime_error(std::string("snapshot validator: missing /Header/") + name);
  H5Handle space(H5Aget_space(attr.get()));
  H5Handle type(H5Aget_type(attr.get()));
  if (!space.valid() || !type.valid() || H5Sget_simple_extent_ndims(space.get()) != 1 ||
      H5Tget_class(type.get()) != H5T_INTEGER) {
    throw std::runtime_error(std::string("snapshot validator: invalid count attribute ") + name);
  }
  hsize_t dims[1] = {0U};
  if (H5Sget_simple_extent_dims(space.get(), dims, nullptr) < 0 || dims[0] != 6U) {
    throw std::runtime_error(std::string("snapshot validator: count attribute must have six entries: ") + name);
  }
  std::array<std::uint64_t, 6> result{};
  if (H5Aread(attr.get(), H5T_NATIVE_UINT64, result.data()) < 0) {
    throw std::runtime_error(std::string("snapshot validator: failed reading count attribute ") + name);
  }
  return result;
}

[[nodiscard]] std::array<double, 6> readMassTable(hid_t header) {
  std::array<double, 6> result{};
  if (H5Aexists(header, "MassTable") <= 0) return result;
  H5Handle attr(H5Aopen(header, "MassTable", H5P_DEFAULT));
  H5Handle space(attr.valid() ? H5Aget_space(attr.get()) : -1);
  if (!attr.valid() || !space.valid() || H5Sget_simple_extent_ndims(space.get()) != 1) {
    throw std::runtime_error("snapshot validator: invalid MassTable attribute");
  }
  hsize_t dims[1] = {0U};
  if (H5Sget_simple_extent_dims(space.get(), dims, nullptr) < 0 || dims[0] != 6U ||
      H5Aread(attr.get(), H5T_NATIVE_DOUBLE, result.data()) < 0) {
    throw std::runtime_error("snapshot validator: invalid MassTable attribute contents");
  }
  return result;
}

[[nodiscard]] H5Handle openPartType(hid_t file, std::size_t type_index) {
  const std::string canonical = "/PartType" + std::to_string(type_index);
  if (H5Lexists(file, canonical.c_str(), H5P_DEFAULT) > 0) {
    return H5Handle(H5Gopen2(file, canonical.c_str(), H5P_DEFAULT));
  }
  const std::string alias = "/ParticleType" + std::to_string(type_index);
  if (H5Lexists(file, alias.c_str(), H5P_DEFAULT) > 0) {
    return H5Handle(H5Gopen2(file, alias.c_str(), H5P_DEFAULT));
  }
  return H5Handle{};
}

[[nodiscard]] std::size_t requireDataset(
    hid_t group,
    const char* name,
    std::uint64_t expected_rows,
    H5T_class_t expected_class,
    std::size_t expected_width_bytes,
    std::uint32_t vector_width,
    const SnapshotReadBudget& budget,
    SnapshotValidationReport& report) {
  H5Handle dataset(H5Dopen2(group, name, H5P_DEFAULT));
  if (!dataset.valid()) throw std::runtime_error(std::string("snapshot validator: missing dataset ") + name);
  H5Handle space(H5Dget_space(dataset.get()));
  H5Handle type(H5Dget_type(dataset.get()));
  if (!space.valid() || !type.valid() || H5Tget_class(type.get()) != expected_class) {
    throw std::runtime_error(std::string("snapshot validator: dataset type mismatch: ") + name);
  }
  const std::size_t actual_width_bytes = H5Tget_size(type.get());
  const bool width_matches = expected_width_bytes != 0U
      ? actual_width_bytes == expected_width_bytes
      : (actual_width_bytes == 4U || actual_width_bytes == 8U);
  if (!width_matches) {
    throw std::runtime_error(std::string("snapshot validator: dataset scalar width mismatch: ") + name);
  }
  const int rank = H5Sget_simple_extent_ndims(space.get());
  hsize_t dims[2] = {0U, 0U};
  if (rank < 1 || rank > 2 || H5Sget_simple_extent_dims(space.get(), dims, nullptr) < 0) {
    throw std::runtime_error(std::string("snapshot validator: invalid dataset rank: ") + name);
  }
  if (dims[0] != expected_rows ||
      (vector_width == 1U && rank != 1) ||
      (vector_width != 1U && (rank != 2 || dims[1] != vector_width))) {
    throw std::runtime_error(std::string("snapshot validator: dataset shape mismatch: ") + name);
  }
  const std::uint64_t scalars = checkedMul(expected_rows, vector_width, name);
  const std::uint64_t bytes = checkedMul(scalars, actual_width_bytes, name);
  if (bytes > budget.max_dataset_bytes) {
    throw std::length_error(std::string("snapshot validator: dataset exceeds max_dataset_bytes: ") + name);
  }
  report.datasets_checked = checkedAdd(report.datasets_checked, 1U, "dataset count");
  report.scalar_values_checked = checkedAdd(report.scalar_values_checked, scalars, "scalar values");
  return actual_width_bytes;
}

[[nodiscard]] bool datasetExists(hid_t group, const char* name) {
  return H5Lexists(group, name, H5P_DEFAULT) > 0;
}

void streamFiniteDoubleDataset(
    hid_t group,
    const char* name,
    std::uint64_t rows,
    std::uint32_t width,
    bool require_positive,
    bool check_metallicity,
    std::optional<std::array<double, 3>> coordinate_upper_bounds) {
  H5Handle dataset(H5Dopen2(group, name, H5P_DEFAULT));
  if (!dataset.valid()) throw std::runtime_error(std::string("snapshot validator: cannot open dataset ") + name);
  constexpr std::uint64_t kChunkRows = 65536U;
  const std::uint64_t chunk_rows = std::min(rows, kChunkRows);
  std::vector<double> buffer(static_cast<std::size_t>(checkedMul(chunk_rows, width, name)));
  for (std::uint64_t offset = 0U; offset < rows; offset += kChunkRows) {
    const std::uint64_t count = std::min(kChunkRows, rows - offset);
    H5Handle file_space(H5Dget_space(dataset.get()));
    if (!file_space.valid()) throw std::runtime_error("snapshot validator: failed obtaining file dataspace");
    hsize_t start[2] = {offset, 0U};
    hsize_t dims[2] = {count, width};
    const int rank = width == 1U ? 1 : 2;
    if (H5Sselect_hyperslab(file_space.get(), H5S_SELECT_SET, start, nullptr, dims, nullptr) < 0) {
      throw std::runtime_error("snapshot validator: failed selecting dataset hyperslab");
    }
    H5Handle memory_space(H5Screate_simple(rank, dims, nullptr));
    if (!memory_space.valid() || H5Dread(dataset.get(), H5T_NATIVE_DOUBLE, memory_space.get(),
                                         file_space.get(), H5P_DEFAULT, buffer.data()) < 0) {
      throw std::runtime_error(std::string("snapshot validator: failed streaming dataset ") + name);
    }
    const std::size_t scalar_count = static_cast<std::size_t>(count * width);
    for (std::size_t i = 0; i < scalar_count; ++i) {
      const double value = buffer[i];
      if (!std::isfinite(value)) {
        throw std::runtime_error(std::string("snapshot validator: non-finite value in ") + name);
      }
      if (require_positive && !(value > 0.0)) {
        throw std::runtime_error(std::string("snapshot validator: non-positive value in ") + name);
      }
      if (check_metallicity && (value < 0.0 || value > 1.0)) {
        throw std::runtime_error(std::string("snapshot validator: metallicity outside [0,1] in ") + name);
      }
      if (coordinate_upper_bounds.has_value()) {
        const std::size_t component = width == 1U ? 0U : (i % width);
        const double upper_bound = coordinate_upper_bounds->at(component);
        if (!(upper_bound > 0.0) || !std::isfinite(upper_bound)) {
          throw std::runtime_error("snapshot validator: invalid coordinate box bound");
        }
        const double tolerance = 128.0 * std::numeric_limits<double>::epsilon() *
                                 std::max(1.0, upper_bound);
        if (value < -tolerance || value > upper_bound + tolerance) {
          throw std::runtime_error(std::string("snapshot validator: coordinate outside box in ") + name);
        }
      }
    }
  }
}

void streamParticleIds(
    hid_t group,
    std::uint64_t rows,
    std::vector<std::uint64_t>& ids,
    const SnapshotReadBudget& budget,
    SnapshotValidationReport& report) {
  H5Handle dataset(H5Dopen2(group, "ParticleIDs", H5P_DEFAULT));
  if (!dataset.valid()) throw std::runtime_error("snapshot validator: cannot open ParticleIDs");
  constexpr std::uint64_t kChunkRows = 65536U;
  std::vector<std::uint64_t> buffer(static_cast<std::size_t>(std::min(rows, kChunkRows)));
  for (std::uint64_t offset = 0U; offset < rows; offset += kChunkRows) {
    const std::uint64_t count = std::min(kChunkRows, rows - offset);
    H5Handle file_space(H5Dget_space(dataset.get()));
    hsize_t start[1] = {offset};
    hsize_t dims[1] = {count};
    if (!file_space.valid() ||
        H5Sselect_hyperslab(file_space.get(), H5S_SELECT_SET, start, nullptr, dims, nullptr) < 0) {
      throw std::runtime_error("snapshot validator: failed selecting ParticleIDs hyperslab");
    }
    H5Handle memory_space(H5Screate_simple(1, dims, nullptr));
    if (!memory_space.valid() || H5Dread(dataset.get(), H5T_NATIVE_UINT64, memory_space.get(),
                                         file_space.get(), H5P_DEFAULT, buffer.data()) < 0) {
      throw std::runtime_error("snapshot validator: failed streaming ParticleIDs");
    }
    for (std::size_t i = 0; i < static_cast<std::size_t>(count); ++i) {
      if (buffer[i] == 0U) throw std::runtime_error("snapshot validator: ParticleIDs contains zero");
      ids.push_back(buffer[i]);
    }
    report.particle_ids_checked = checkedAdd(report.particle_ids_checked, count, "particle IDs");
    const std::uint64_t id_bytes = checkedMul(ids.size(), sizeof(std::uint64_t), "validation IDs");
    report.validation_peak_id_bytes = std::max(report.validation_peak_id_bytes, id_bytes);
    if (id_bytes > budget.max_validation_id_bytes) {
      throw std::length_error("snapshot validator: global ID uniqueness proof exceeds max_validation_id_bytes");
    }
  }
}

void validateChuiExtensions(
    hid_t group,
    std::size_t type_index,
    std::uint64_t rows,
    const SnapshotReadBudget& budget,
    SnapshotValidationReport& report) {
  internal::forEachRequiredChuiSnapshotField(type_index, [&](const auto& field) {
    H5T_class_t type_class = H5T_NO_CLASS;
    std::size_t width_bytes = 0U;
    switch (field.storage) {
      case internal::SnapshotFieldStorage::kFloat64:
        type_class = H5T_FLOAT; width_bytes = 8U; break;
      case internal::SnapshotFieldStorage::kUint64:
        type_class = H5T_INTEGER; width_bytes = 8U; break;
      case internal::SnapshotFieldStorage::kUint32:
        type_class = H5T_INTEGER; width_bytes = 4U; break;
      case internal::SnapshotFieldStorage::kUint8:
        type_class = H5T_INTEGER; width_bytes = 1U; break;
    }
    static_cast<void>(requireDataset(
        group, field.name, rows, type_class, width_bytes, 1U, budget, report));
    if (field.require_finite) {
      streamFiniteDoubleDataset(
          group, field.name, rows, 1U, false, field.metallicity_fraction, std::nullopt);
    }
  });
}

void validateMember(
    const std::filesystem::path& path,
    const SnapshotSetInspection& inspection,
    const SnapshotReadOptions& options,
    std::vector<std::uint64_t>& ids,
    SnapshotValidationReport& report) {
  H5Handle file(H5Fopen(path.string().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
  if (!file.valid()) throw std::runtime_error("snapshot validator: cannot open member " + path.string());
  H5Handle header(H5Gopen2(file.get(), "/Header", H5P_DEFAULT));
  if (!header.valid()) throw std::runtime_error("snapshot validator: member missing /Header");
  const auto local_counts = readCountArray(header.get(), "NumPart_ThisFile");
  const auto mass_table = readMassTable(header.get());
  const bool chui_native = inspection.dialect == SnapshotDialect::kChuiNative ||
                           inspection.schema_name.rfind("chui_science_snapshot_", 0U) == 0U;
  for (std::size_t type_index = 0; type_index < local_counts.size(); ++type_index) {
    const std::uint64_t rows = local_counts[type_index];
    if (rows == 0U) continue;
    H5Handle group = openPartType(file.get(), type_index);
    if (!group.valid()) {
      throw std::runtime_error("snapshot validator: missing particle group for populated type " +
                               std::to_string(type_index));
    }
    static_cast<void>(requireDataset(
        group.get(), "Coordinates", rows, H5T_FLOAT, chui_native ? 8U : 0U, 3U, options.budget, report));
    std::optional<std::array<double, 3>> coordinate_bounds;
    if (chui_native) {
      double storage_scale = 1.0;
      if (inspection.dialect == SnapshotDialect::kArepoFormat3 ||
          inspection.dialect == SnapshotDialect::kGadget4Hdf5) {
        storage_scale = inspection.hubble_param;
      }
      coordinate_bounds = std::array<double, 3>{
          inspection.box_size_x * storage_scale,
          inspection.box_size_y * storage_scale,
          inspection.box_size_z * storage_scale};
    }
    streamFiniteDoubleDataset(
        group.get(), "Coordinates", rows, 3U, false, false, coordinate_bounds);
    if (options.require_velocities || chui_native) {
      static_cast<void>(requireDataset(
          group.get(), "Velocities", rows, H5T_FLOAT, chui_native ? 8U : 0U, 3U, options.budget, report));
      streamFiniteDoubleDataset(group.get(), "Velocities", rows, 3U, false, false, std::nullopt);
    }
    if (options.require_ids || chui_native) {
      static_cast<void>(requireDataset(
          group.get(), "ParticleIDs", rows, H5T_INTEGER, chui_native ? 8U : 0U, 1U, options.budget, report));
      streamParticleIds(group.get(), rows, ids, options.budget, report);
    }
    if (datasetExists(group.get(), "Masses")) {
      static_cast<void>(requireDataset(
          group.get(), "Masses", rows, H5T_FLOAT, chui_native ? 8U : 0U, 1U, options.budget, report));
      streamFiniteDoubleDataset(group.get(), "Masses", rows, 1U, true, false, std::nullopt);
    } else if (!(mass_table[type_index] > 0.0 && std::isfinite(mass_table[type_index]))) {
      throw std::runtime_error("snapshot validator: populated type has neither Masses nor positive MassTable");
    }
    if (chui_native) validateChuiExtensions(group.get(), type_index, rows, options.budget, report);
  }
}
#endif

}  // namespace

void SnapshotValidationReport::requireValid() const {
  if (valid) return;
  std::ostringstream message;
  message << "snapshot validation failed";
  for (const auto& diagnostic : diagnostics) message << "; " << diagnostic;
  throw std::runtime_error(message.str());
}

SnapshotValidationReport validateSnapshotSetHdf5(
    const std::filesystem::path& input_path,
    const SnapshotReadOptions& options) {
  SnapshotValidationReport report;
#if !COSMOSIM_ENABLE_HDF5
  static_cast<void>(input_path);
  static_cast<void>(options);
  report.diagnostics.push_back("COSMOSIM_ENABLE_HDF5=OFF: snapshot validation unavailable");
  return report;
#else
  try {
    report.inspection = inspectSnapshotSet(input_path, options);
    if (!report.inspection.complete) {
      throw std::runtime_error("logical snapshot set is incomplete");
    }
    std::uint64_t expected_particles = 0U;
    for (const auto count : report.inspection.global_part_count) {
      expected_particles = checkedAdd(expected_particles, count, "global particle count");
    }
    if (expected_particles > options.budget.max_particles) {
      throw std::length_error("snapshot validator: global particle count exceeds max_particles");
    }
    const bool validate_ids = options.require_ids ||
        report.inspection.dialect == SnapshotDialect::kChuiNative ||
        report.inspection.schema_name.rfind("chui_science_snapshot_", 0U) == 0U;
    std::vector<std::uint64_t> ids;
    if (validate_ids) {
      const std::uint64_t id_bytes = checkedMul(expected_particles, sizeof(std::uint64_t), "global ID buffer");
      if (id_bytes > options.budget.max_validation_id_bytes) {
        throw std::length_error("snapshot validator: uniqueness proof exceeds max_validation_id_bytes");
      }
      ids.reserve(static_cast<std::size_t>(expected_particles));
    }
    for (const auto& member : report.inspection.member_paths) {
      validateMember(member, report.inspection, options, ids, report);
    }
    if (validate_ids) {
      if (ids.size() != expected_particles) {
        throw std::runtime_error(
            "snapshot validator: validated ParticleIDs count disagrees with global header count");
      }
      std::sort(ids.begin(), ids.end());
      if (std::adjacent_find(ids.begin(), ids.end()) != ids.end()) {
        throw std::runtime_error("snapshot validator: duplicate ParticleIDs across logical snapshot set");
      }
    }
    report.valid = true;
    report.diagnostics.push_back("direct HDF5 structural/scientific validation passed");
  } catch (const std::exception& error) {
    report.valid = false;
    report.diagnostics.push_back(error.what());
  }
  return report;
#endif
}

}  // namespace cosmosim::io
