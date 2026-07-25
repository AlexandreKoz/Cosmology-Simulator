#include "io/internal/ic_reader_session.hpp"

#if COSMOSIM_ENABLE_HDF5
#include <array>
#include <limits>
#include <stdexcept>
#include <string>

namespace cosmosim::io::internal {
namespace {

void checkedCounterAdd(
    std::uint64_t& destination,
    std::uint64_t value,
    std::string_view counter_name) {
  if (destination > std::numeric_limits<std::uint64_t>::max() - value) {
    throw std::overflow_error(
        "IC " + std::string(counter_name) + " counter overflow");
  }
  destination += value;
}

}  // namespace

IcReaderSession::IcReaderSession(
    const std::filesystem::path& path,
    std::uint64_t expected_size_bytes,
    std::string_view expected_sha256,
    IcImportCounters& counters)
    : m_path(path),
      m_file(H5Fopen(path.string().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT)) {
  if (!m_file.valid()) {
    throw std::runtime_error(
        "failed to open persistent IC reader session for " + path.string());
  }
  checkedCounterAdd(
      counters.source_file_open_count, 1U, "source_file_open_count");
  const std::uint64_t observed_size = std::filesystem::file_size(path);
  if (observed_size != expected_size_bytes) {
    throw std::runtime_error(
        "IC source identity changed before payload read: file size mismatch for " +
        path.string());
  }
  const std::string observed_sha256 = icSha256FileHex(path);
  checkedCounterAdd(counters.hash_bytes_read, observed_size, "hash_bytes_read");
  if (observed_sha256 != expected_sha256) {
    throw std::runtime_error(
        "IC source identity changed before payload read: SHA-256 mismatch for " +
        path.string());
  }
}

hid_t IcReaderSession::dataset(
    const IcFieldManifest& field,
    IcImportCounters& counters) {
  auto existing = m_datasets.find(field.dataset_path);
  if (existing != m_datasets.end()) {
    return existing->second.get();
  }
  const std::size_t separator = field.dataset_path.rfind('/');
  if (separator == std::string::npos) {
    throw std::logic_error(
        "IC dataset manifest path is not absolute: " + field.dataset_path);
  }
  const std::string source_path =
      field.dataset_path.substr(0U, separator + 1U) + field.selected_alias;
  Hdf5Handle handle(H5Dopen2(m_file.get(), source_path.c_str(), H5P_DEFAULT));
  if (!handle.valid()) {
    throw std::runtime_error(
        "failed to open persistent IC dataset " + source_path + " in " +
        m_path.string());
  }
  checkedCounterAdd(
      counters.source_dataset_open_count, 1U, "source_dataset_open_count");
  auto [inserted, ok] =
      m_datasets.emplace(field.dataset_path, std::move(handle));
  if (!ok) {
    throw std::logic_error("failed to cache persistent IC dataset handle");
  }
  return inserted->second.get();
}

void readChunkDouble(
    IcReaderSession& session,
    const IcFieldManifest& field,
    std::size_t start,
    std::size_t count,
    std::size_t components,
    std::vector<double>& out,
    IcImportCounters& counters) {
  const hid_t dataset = session.dataset(field, counters);
  Hdf5Handle file_space(H5Dget_space(dataset));
  if (!file_space.valid()) {
    throw std::runtime_error("failed to inspect " + field.dataset_path);
  }
  const int rank = components == 1U ? 1 : 2;
  std::array<hsize_t, 2> offset{static_cast<hsize_t>(start), 0U};
  std::array<hsize_t, 2> extent{
      static_cast<hsize_t>(count), static_cast<hsize_t>(components)};
  if (H5Sselect_hyperslab(
          file_space.get(), H5S_SELECT_SET, offset.data(), nullptr,
          extent.data(), nullptr) < 0) {
    throw std::runtime_error("failed hyperslab for " + field.dataset_path);
  }
  Hdf5Handle memory_space(H5Screate_simple(rank, extent.data(), nullptr));
  if (!memory_space.valid()) {
    throw std::runtime_error(
        "failed to create memory dataspace for " + field.dataset_path);
  }
  out.resize(count * components);
  if (H5Dread(
          dataset, H5T_NATIVE_DOUBLE, memory_space.get(), file_space.get(),
          H5P_DEFAULT, out.data()) < 0) {
    throw std::runtime_error("failed chunk read for " + field.dataset_path);
  }
  const std::uint64_t value_count =
      static_cast<std::uint64_t>(count) * components;
  checkedCounterAdd(
      counters.payload_bytes_read, value_count * field.byte_width,
      "payload_bytes_read");
  checkedCounterAdd(
      counters.converted_payload_bytes, value_count * sizeof(double),
      "converted_payload_bytes");
}

void readChunkU64(
    IcReaderSession& session,
    const IcFieldManifest& field,
    std::size_t start,
    std::size_t count,
    std::vector<std::uint64_t>& out,
    IcImportCounters& counters) {
  const hid_t dataset = session.dataset(field, counters);
  Hdf5Handle file_space(H5Dget_space(dataset));
  if (!file_space.valid()) {
    throw std::runtime_error("failed to inspect " + field.dataset_path);
  }
  hsize_t offset[1]{static_cast<hsize_t>(start)};
  hsize_t extent[1]{static_cast<hsize_t>(count)};
  if (H5Sselect_hyperslab(
          file_space.get(), H5S_SELECT_SET, offset, nullptr, extent,
          nullptr) < 0) {
    throw std::runtime_error("failed hyperslab for " + field.dataset_path);
  }
  Hdf5Handle memory_space(H5Screate_simple(1, extent, nullptr));
  if (!memory_space.valid()) {
    throw std::runtime_error(
        "failed to create memory dataspace for " + field.dataset_path);
  }
  out.resize(count);
  if (H5Dread(
          dataset, H5T_NATIVE_UINT64, memory_space.get(), file_space.get(),
          H5P_DEFAULT, out.data()) < 0) {
    throw std::runtime_error("failed ID chunk read for " + field.dataset_path);
  }
  checkedCounterAdd(
      counters.payload_bytes_read,
      static_cast<std::uint64_t>(count) * field.byte_width,
      "payload_bytes_read");
  checkedCounterAdd(
      counters.converted_payload_bytes,
      static_cast<std::uint64_t>(count) * sizeof(std::uint64_t),
      "converted_payload_bytes");
}

}  // namespace cosmosim::io::internal
#endif
