#include "io/internal/transactional_file.hpp"

#include <atomic>
#include <chrono>
#include <cerrno>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <random>
#include <stdexcept>
#include <string>
#include <system_error>

#if defined(_WIN32)
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#else
#include <fcntl.h>
#include <unistd.h>
#endif

namespace cosmosim::io::internal {
namespace {

std::atomic<std::uint64_t> g_temp_counter{0};

#if defined(_WIN32)
[[nodiscard]] std::string platformErrorMessage(const char* prefix) {
  return std::string(prefix) + " (Win32 error " +
      std::to_string(static_cast<unsigned long>(GetLastError())) + ")";
}
#endif

[[nodiscard]] std::filesystem::path effectiveParent(
    const std::filesystem::path& path) {
  const auto parent = path.parent_path();
  return parent.empty() ? std::filesystem::path(".") : parent;
}

}  // namespace

void ensureParentDirectory(const std::filesystem::path& output_path) {
  const std::filesystem::path parent = output_path.parent_path();
  if (parent.empty()) {
    return;
  }
  std::error_code error;
  std::filesystem::create_directories(parent, error);
  if (error) {
    throw std::runtime_error(
        "failed to create output parent directory " + parent.string() +
        ": " + error.message());
  }
}

void validateTemporarySuffix(std::string_view suffix) {
  if (suffix.empty()) {
    throw std::invalid_argument("transactional temporary suffix must not be empty");
  }
  if (suffix.find('/') != std::string_view::npos ||
      suffix.find('\\') != std::string_view::npos) {
    throw std::invalid_argument(
        "transactional temporary suffix must not contain path separators");
  }
  if (suffix == "." || suffix == "..") {
    throw std::invalid_argument("transactional temporary suffix is invalid");
  }
}

std::filesystem::path makeUniqueTemporarySibling(
    const std::filesystem::path& output_path,
    std::string_view suffix) {
  validateTemporarySuffix(suffix);
  ensureParentDirectory(output_path);
  const auto parent = effectiveParent(output_path);
  const std::string base = output_path.filename().string();
  if (base.empty()) {
    throw std::invalid_argument("transactional output path must name a file");
  }

  std::random_device random_device;
  for (int attempt = 0; attempt < 64; ++attempt) {
    const auto now = static_cast<std::uint64_t>(
        std::chrono::steady_clock::now().time_since_epoch().count());
    const std::uint64_t counter = g_temp_counter.fetch_add(1, std::memory_order_relaxed);
    const std::uint64_t random_bits =
        (static_cast<std::uint64_t>(random_device()) << 32U) ^
        static_cast<std::uint64_t>(random_device());
    const std::string candidate_name =
        base + std::string(suffix) + "." + std::to_string(now) + "." +
        std::to_string(counter) + "." + std::to_string(random_bits);
    const auto candidate = parent / candidate_name;
    if (candidate == output_path) {
      continue;
    }
    std::error_code error;
    const bool exists = std::filesystem::exists(candidate, error);
    if (!error && !exists) {
      return candidate;
    }
  }
  throw std::runtime_error(
      "failed to allocate a collision-free temporary sibling for " +
      output_path.string());
}

void syncFileForDurability(const std::filesystem::path& path) {
#if defined(_WIN32)
  HANDLE handle = CreateFileW(
      path.c_str(), GENERIC_READ, FILE_SHARE_READ | FILE_SHARE_WRITE | FILE_SHARE_DELETE,
      nullptr, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, nullptr);
  if (handle == INVALID_HANDLE_VALUE) {
    throw std::runtime_error(platformErrorMessage("failed to open file for durable flush"));
  }
  if (!FlushFileBuffers(handle)) {
    const std::string message = platformErrorMessage("failed to FlushFileBuffers");
    CloseHandle(handle);
    throw std::runtime_error(message);
  }
  if (!CloseHandle(handle)) {
    throw std::runtime_error(platformErrorMessage("failed to close durable-flush handle"));
  }
#else
  const int fd = ::open(path.c_str(), O_RDONLY);
  if (fd < 0) {
    throw std::runtime_error(
        "failed to open file for fsync " + path.string() + ": " +
        std::strerror(errno));
  }
  if (::fsync(fd) != 0) {
    const std::string message = std::strerror(errno);
    ::close(fd);
    throw std::runtime_error(
        "failed to fsync file " + path.string() + ": " + message);
  }
  if (::close(fd) != 0) {
    throw std::runtime_error(
        "failed to close file after fsync " + path.string() + ": " +
        std::strerror(errno));
  }
#endif
}

void syncDirectoryForDurability(const std::filesystem::path& directory) {
#if defined(_WIN32)
  // ReplaceFileW / MoveFileExW with MOVEFILE_WRITE_THROUGH provide the Windows
  // publication durability primitive. Directory handles are not flushed here.
  static_cast<void>(directory);
#else
  const std::filesystem::path target = directory.empty() ? std::filesystem::path(".") : directory;
  const int fd = ::open(target.c_str(), O_RDONLY | O_DIRECTORY);
  if (fd < 0) {
    throw std::runtime_error(
        "failed to open directory for fsync " + target.string() + ": " +
        std::strerror(errno));
  }
  if (::fsync(fd) != 0) {
    const std::string message = std::strerror(errno);
    ::close(fd);
    throw std::runtime_error(
        "failed to fsync directory " + target.string() + ": " + message);
  }
  if (::close(fd) != 0) {
    throw std::runtime_error(
        "failed to close directory after fsync " + target.string() + ": " +
        std::strerror(errno));
  }
#endif
}

void atomicReplaceFile(
    const std::filesystem::path& temporary_path,
    const std::filesystem::path& output_path,
    FileDurability durability) {
  if (temporary_path == output_path) {
    throw std::invalid_argument("transactional temporary path equals final path");
  }
  ensureParentDirectory(output_path);
  if (durability == FileDurability::kDurablePublication) {
    syncFileForDurability(temporary_path);
  }

#if defined(_WIN32)
  const bool destination_exists = std::filesystem::exists(output_path);
  BOOL ok = FALSE;
  if (destination_exists) {
    ok = ReplaceFileW(
        output_path.c_str(), temporary_path.c_str(), nullptr,
        REPLACEFILE_WRITE_THROUGH, nullptr, nullptr);
    if (!ok && GetLastError() == ERROR_FILE_NOT_FOUND) {
      ok = MoveFileExW(
          temporary_path.c_str(), output_path.c_str(),
          MOVEFILE_REPLACE_EXISTING |
              (durability == FileDurability::kDurablePublication
                   ? MOVEFILE_WRITE_THROUGH
                   : 0));
    }
  } else {
    ok = MoveFileExW(
        temporary_path.c_str(), output_path.c_str(),
        MOVEFILE_REPLACE_EXISTING |
            (durability == FileDurability::kDurablePublication
                 ? MOVEFILE_WRITE_THROUGH
                 : 0));
  }
  if (!ok) {
    throw std::runtime_error(platformErrorMessage("failed to atomically publish file"));
  }
#else
  if (::rename(temporary_path.c_str(), output_path.c_str()) != 0) {
    throw std::runtime_error(
        "failed to atomically publish " + temporary_path.string() + " as " +
        output_path.string() + ": " + std::strerror(errno));
  }
  if (durability == FileDurability::kDurablePublication) {
    syncDirectoryForDurability(effectiveParent(output_path));
  }
#endif
}

TransactionalFileTarget::TransactionalFileTarget(
    std::filesystem::path output_path,
    std::string_view temporary_suffix,
    FileDurability durability)
    : m_output_path(std::move(output_path)),
      m_temporary_path(makeUniqueTemporarySibling(m_output_path, temporary_suffix)),
      m_durability(durability) {}

TransactionalFileTarget::~TransactionalFileTarget() {
  if (m_published || m_temporary_path.empty()) {
    return;
  }
  std::error_code ignored;
  std::filesystem::remove(m_temporary_path, ignored);
}

const std::filesystem::path& TransactionalFileTarget::outputPath() const noexcept {
  return m_output_path;
}

const std::filesystem::path& TransactionalFileTarget::temporaryPath() const noexcept {
  return m_temporary_path;
}

FileDurability TransactionalFileTarget::durability() const noexcept {
  return m_durability;
}

void TransactionalFileTarget::publish() {
  if (m_published) {
    throw std::logic_error("transactional file target was published twice");
  }
  atomicReplaceFile(m_temporary_path, m_output_path, m_durability);
  m_published = true;
}

void writeTextFileTransactionally(
    const std::filesystem::path& output_path,
    std::string_view contents,
    FileDurability durability) {
  TransactionalFileTarget target(output_path, ".part", durability);
  {
    std::ofstream output(target.temporaryPath(), std::ios::binary | std::ios::trunc);
    if (!output) {
      throw std::runtime_error(
          "failed to create transactional text file " +
          target.temporaryPath().string());
    }
    output.write(contents.data(), static_cast<std::streamsize>(contents.size()));
    output.flush();
    if (!output) {
      throw std::runtime_error(
          "failed to write transactional text file " +
          target.temporaryPath().string());
    }
  }
  target.publish();
}

}  // namespace cosmosim::io::internal
