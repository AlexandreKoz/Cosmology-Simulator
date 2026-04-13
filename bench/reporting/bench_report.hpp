#pragma once

#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "cosmosim/core/build_config.hpp"

namespace cosmosim::bench {

struct BenchmarkExecutionConfig {
  std::size_t warmup_iterations = 2;
  std::size_t measurement_iterations = 10;
  int threads = 1;
  int mpi_ranks = 1;
  std::string device = "cpu";
};

struct BenchmarkClock {
  using clock_type = std::chrono::steady_clock;
  using time_point = clock_type::time_point;

  [[nodiscard]] static time_point now() { return clock_type::now(); }

  [[nodiscard]] static double millisecondsBetween(time_point begin, time_point end) {
    return std::chrono::duration<double, std::milli>(end - begin).count();
  }
};

class BenchmarkReporter {
 public:
  explicit BenchmarkReporter(std::string benchmark_name) : m_benchmark_name(std::move(benchmark_name)) {}

  void addField(std::string key, std::string value) {
    m_fields.emplace_back(std::move(key), std::move(value));
  }

  void addField(std::string key, const char* value) { m_fields.emplace_back(std::move(key), std::string(value)); }

  template <typename Numeric>
  void addField(std::string key, Numeric value) {
    m_fields.emplace_back(std::move(key), std::to_string(value));
  }

  void write(std::ostream& out = std::cout) const {
    out << m_benchmark_name;
    for (const auto& [key, value] : m_fields) {
      out << ' ' << key << '=' << value;
    }
    out << '\n';
  }

 private:
  std::string m_benchmark_name;
  std::vector<std::pair<std::string, std::string>> m_fields;
};

[[nodiscard]] inline int parsePositiveEnv(std::string_view name, int fallback) {
  const char* raw_value = std::getenv(std::string(name).c_str());
  if (raw_value == nullptr) {
    return fallback;
  }
  char* end_ptr = nullptr;
  const long parsed = std::strtol(raw_value, &end_ptr, 10);
  if (end_ptr == raw_value || parsed <= 0) {
    return fallback;
  }
  return static_cast<int>(parsed);
}

[[nodiscard]] inline BenchmarkExecutionConfig defaultExecutionConfig(
    std::size_t warmup_iterations,
    std::size_t measurement_iterations) {
  BenchmarkExecutionConfig config;
  config.warmup_iterations = warmup_iterations;
  config.measurement_iterations = measurement_iterations;
  config.threads = parsePositiveEnv("COSMOSIM_BENCH_THREADS", 1);
  config.mpi_ranks = parsePositiveEnv("COSMOSIM_BENCH_MPI_RANKS", 1);
  if (const char* device = std::getenv("COSMOSIM_BENCH_DEVICE"); device != nullptr && *device != '\0') {
    config.device = device;
  }
  return config;
}

inline void addExecutionFields(BenchmarkReporter& reporter, const BenchmarkExecutionConfig& execution) {
  reporter.addField("build_type", COSMOSIM_BUILD_TYPE);
  reporter.addField("threads", execution.threads);
  reporter.addField("mpi_ranks", execution.mpi_ranks);
  reporter.addField("device", execution.device);
  reporter.addField("warmup_iterations", execution.warmup_iterations);
  reporter.addField("measurement_iterations", execution.measurement_iterations);
}

inline void addBandwidthFields(
    BenchmarkReporter& reporter,
    std::uint64_t bytes_moved,
    double elapsed_ms,
    std::string_view bytes_key = "bytes_moved",
    std::string_view bandwidth_key = "effective_bandwidth_gb_s") {
  reporter.addField(std::string(bytes_key), bytes_moved);
  const double elapsed_s = elapsed_ms * 1.0e-3;
  const double bandwidth_gb_s = elapsed_s > 0.0 ? static_cast<double>(bytes_moved) / elapsed_s * 1.0e-9 : 0.0;
  reporter.addField(std::string(bandwidth_key), bandwidth_gb_s);
}

}  // namespace cosmosim::bench
