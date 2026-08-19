#pragma once

#include <atomic>
#include <chrono>
#include <cctype>
#include <cstdint>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>

#if defined(_WIN32)
#include <process.h>
#else
#include <unistd.h>
#endif

namespace cosmosim::test_support {

// Collision-resistant test scratch namespace.  Every TestTempWorkspace instance
// owns exactly one directory and removes only that directory at destruction.
// Tests that need one path shared by MPI ranks should broadcast a nonce/token
// from rank zero and call createMpiShared() with the same token on every rank.
// Exactly one rank (normally rank zero) owns cleanup; all other ranks are
// non-owning users of the same directory.
class TestTempWorkspace {
 public:
  TestTempWorkspace() = default;
  TestTempWorkspace(const TestTempWorkspace&) = delete;
  TestTempWorkspace& operator=(const TestTempWorkspace&) = delete;

  TestTempWorkspace(TestTempWorkspace&& other) noexcept
      : m_root(std::move(other.m_root)), m_owns(other.m_owns) {
    other.m_owns = false;
  }

  TestTempWorkspace& operator=(TestTempWorkspace&& other) noexcept {
    if (this != &other) {
      cleanup();
      m_root = std::move(other.m_root);
      m_owns = other.m_owns;
      other.m_owns = false;
    }
    return *this;
  }

  ~TestTempWorkspace() { cleanup(); }

  [[nodiscard]] static TestTempWorkspace createProcessLocal(std::string_view test_stem) {
    return createOwned(uniqueProcessLocalPath(test_stem));
  }

  [[nodiscard]] static TestTempWorkspace createUniqueDirectory(std::string_view test_stem) {
    return createProcessLocal(test_stem);
  }

  [[nodiscard]] static TestTempWorkspace createMpiShared(
      std::string_view test_stem,
      std::uint64_t shared_run_token,
      bool owns_cleanup) {
    const auto root = scratchRoot() /
        (sanitize(test_stem) + "_mpi_" + std::to_string(shared_run_token));
    std::error_code error;
    const bool created = std::filesystem::create_directories(root, error);
    if (error || (!created && !std::filesystem::is_directory(root))) {
      throw std::runtime_error(
          "failed to create MPI-shared test temporary workspace '" + root.string() +
          "': " + error.message());
    }
    return TestTempWorkspace(root, owns_cleanup);
  }

  [[nodiscard]] static std::uint64_t uniqueRunToken() noexcept {
    const std::uint64_t sequence = s_sequence.fetch_add(1U, std::memory_order_relaxed);
    const auto ticks = static_cast<std::uint64_t>(
        std::chrono::steady_clock::now().time_since_epoch().count());
    // Mix process identity, monotonic time, and a process-local sequence.
    return ticks ^ (processId() * 0x9E3779B97F4A7C15ULL) ^
        (sequence * 0xBF58476D1CE4E5B9ULL);
  }

  [[nodiscard]] static std::filesystem::path uniqueProcessLocalFile(
      std::string_view stem,
      std::string_view extension) {
    auto path = uniqueProcessLocalPath(stem);
    if (!extension.empty()) {
      const std::string normalized_extension = extension.front() == '.'
          ? std::string(extension)
          : "." + std::string(extension);
      path += normalized_extension;
    }
    return path;
  }

  [[nodiscard]] static std::filesystem::path uniqueProcessLocalPath(std::string_view stem) {
    std::error_code error;
    std::filesystem::create_directories(scratchRoot(), error);
    if (error) {
      throw std::runtime_error(
          "failed to create CosmoSim test scratch root: " + error.message());
    }
    const std::uint64_t sequence = s_sequence.fetch_add(1U, std::memory_order_relaxed);
    const auto ticks = static_cast<std::uint64_t>(
        std::chrono::steady_clock::now().time_since_epoch().count());
    return scratchRoot() /
        (sanitize(stem) + "_p" + std::to_string(processId()) + "_" +
         std::to_string(ticks) + "_" + std::to_string(sequence));
  }

  [[nodiscard]] const std::filesystem::path& root() const noexcept { return m_root; }

  [[nodiscard]] std::filesystem::path path(std::string_view relative_name) const {
    if (m_root.empty()) {
      throw std::logic_error("TestTempWorkspace has no owned root");
    }
    return m_root / std::filesystem::path(relative_name);
  }

 private:
  explicit TestTempWorkspace(std::filesystem::path root, bool owns)
      : m_root(std::move(root)), m_owns(owns) {}

  [[nodiscard]] static TestTempWorkspace createOwned(const std::filesystem::path& root) {
    std::error_code error;
    const bool created = std::filesystem::create_directories(root, error);
    if (error || (!created && !std::filesystem::is_directory(root))) {
      throw std::runtime_error(
          "failed to create test temporary workspace '" + root.string() + "': " +
          error.message());
    }
    return TestTempWorkspace(root, true);
  }

  [[nodiscard]] static std::filesystem::path scratchRoot() {
    return std::filesystem::temp_directory_path() / "cosmosim_tests";
  }

  [[nodiscard]] static std::string sanitize(std::string_view value) {
    std::string result;
    result.reserve(value.size());
    for (const char ch : value) {
      const unsigned char byte = static_cast<unsigned char>(ch);
      if (std::isalnum(byte) != 0 || ch == '_' || ch == '-') {
        result.push_back(ch);
      } else {
        result.push_back('_');
      }
    }
    return result.empty() ? std::string("test") : result;
  }

  [[nodiscard]] static std::uint64_t processId() noexcept {
#if defined(_WIN32)
    return static_cast<std::uint64_t>(_getpid());
#else
    return static_cast<std::uint64_t>(getpid());
#endif
  }

  void cleanup() noexcept {
    if (!m_owns || m_root.empty()) {
      return;
    }
    std::error_code ignored;
    std::filesystem::remove_all(m_root, ignored);
    m_owns = false;
  }

  inline static std::atomic<std::uint64_t> s_sequence{0U};
  std::filesystem::path m_root;
  bool m_owns = false;
};

}  // namespace cosmosim::test_support
