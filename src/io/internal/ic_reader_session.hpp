#pragma once

#include "cosmosim/core/build_config.hpp"

#if COSMOSIM_ENABLE_HDF5
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

#include <hdf5.h>

#include "cosmosim/io/ic_reader.hpp"
#include "io/internal/ic_hdf5_handle.hpp"

namespace cosmosim::io::internal {

class IcReaderSession {
 public:
  IcReaderSession(
      const std::filesystem::path& path,
      std::uint64_t expected_size_bytes,
      std::string_view expected_sha256,
      IcImportCounters& counters);

  IcReaderSession(const IcReaderSession&) = delete;
  IcReaderSession& operator=(const IcReaderSession&) = delete;
  IcReaderSession(IcReaderSession&&) noexcept = default;
  IcReaderSession& operator=(IcReaderSession&&) noexcept = default;

  [[nodiscard]] hid_t dataset(
      const IcFieldManifest& field,
      IcImportCounters& counters);

 private:
  std::filesystem::path m_path;
  Hdf5Handle m_file;
  std::unordered_map<std::string, Hdf5Handle> m_datasets;
};

void readChunkDouble(
    IcReaderSession& session,
    const IcFieldManifest& field,
    std::size_t start,
    std::size_t count,
    std::size_t components,
    std::vector<double>& out,
    IcImportCounters& counters);

void readChunkU64(
    IcReaderSession& session,
    const IcFieldManifest& field,
    std::size_t start,
    std::size_t count,
    std::vector<std::uint64_t>& out,
    IcImportCounters& counters);

}  // namespace cosmosim::io::internal
#endif
