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
#include <chrono>

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
      IcIntegrityMode integrity_mode,
      IcImportCounters& counters);
  // Compatibility/internal-test overload preserves the historical strict
  // completion rehash contract. Production call sites pass an explicit mode.
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

  // Revalidate source identity after the final payload read. Strict mode also
  // recomputes the complete-file digest. This is explicit because destructors
  // must not throw.
  void revalidateSourceIdentity(IcImportCounters& counters) const;

 private:
  struct SourceIdentity {
    std::uint64_t size_bytes = 0U;
    std::filesystem::file_time_type last_write_time{};
    std::uint64_t device = 0U;
    std::uint64_t inode = 0U;
    bool has_native_identity = false;
  };

  [[nodiscard]] static SourceIdentity captureSourceIdentity(
      const std::filesystem::path& path);
  static void requireSameIdentity(
      const SourceIdentity& expected,
      const SourceIdentity& observed,
      const std::filesystem::path& path,
      std::string_view phase);

  std::filesystem::path m_path;
  std::uint64_t m_expected_size_bytes = 0U;
  std::string m_expected_sha256;
  IcIntegrityMode m_integrity_mode = IcIntegrityMode::kVerifiedIdentity;
  SourceIdentity m_open_identity;
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
