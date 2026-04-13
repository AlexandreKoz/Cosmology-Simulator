#include <chrono>
#include <cstdint>
#include <iostream>

#include "cosmosim/core/config.hpp"
#include "cosmosim/core/simulation_mode.hpp"

namespace {

void runDispatchLoop(const cosmosim::core::ModePolicy& policy, std::uint64_t iterations) {
  volatile std::uint64_t sink = 0;
  for (std::uint64_t i = 0; i < iterations; ++i) {
    switch (policy.gravity_boundary) {
      case cosmosim::core::GravityBoundaryModel::kPeriodicPoisson:
        sink += 1;
        break;
      case cosmosim::core::GravityBoundaryModel::kIsolatedMonopoleDirichlet:
        sink += 2;
        break;
    }

    switch (policy.hydro_boundary) {
      case cosmosim::core::BoundaryCondition::kPeriodic:
        sink += 4;
        break;
      case cosmosim::core::BoundaryCondition::kOpen:
        sink += 8;
        break;
      case cosmosim::core::BoundaryCondition::kReflective:
        sink += 16;
        break;
    }
  }
  if (sink == 0) {
    std::cout << "sink=" << sink << '\n';
  }
}

}  // namespace

int main() {
  constexpr std::uint64_t iterations = 100000000;

  cosmosim::core::ModeConfig mode_config;
  mode_config.mode = cosmosim::core::SimulationMode::kZoomIn;
  const cosmosim::core::ModePolicy periodic_policy = cosmosim::core::buildModePolicy(mode_config);

  mode_config.mode = cosmosim::core::SimulationMode::kIsolatedCluster;
  const cosmosim::core::ModePolicy isolated_policy = cosmosim::core::buildModePolicy(mode_config);

  const auto setup_begin = std::chrono::steady_clock::now();
  const auto setup_end = std::chrono::steady_clock::now();

  const auto periodic_begin = std::chrono::steady_clock::now();
  runDispatchLoop(periodic_policy, iterations);
  const auto periodic_end = std::chrono::steady_clock::now();

  const auto isolated_begin = std::chrono::steady_clock::now();
  runDispatchLoop(isolated_policy, iterations);
  const auto isolated_end = std::chrono::steady_clock::now();

  const double setup_ms = std::chrono::duration<double, std::milli>(setup_end - setup_begin).count();
  const double periodic_ms =
      std::chrono::duration<double, std::milli>(periodic_end - periodic_begin).count();
  const double isolated_ms =
      std::chrono::duration<double, std::milli>(isolated_end - isolated_begin).count();

  const double periodic_ns_per_iter = periodic_ms * 1.0e6 / static_cast<double>(iterations);
  const double isolated_ns_per_iter = isolated_ms * 1.0e6 / static_cast<double>(iterations);

  std::cout << "bench_mode_dispatch\n";
  std::cout << "build_type="
#ifdef NDEBUG
            << "Release";
#else
            << "Debug";
#endif
  std::cout << " features=core-only\n";
  std::cout << "threads=1 (single-thread loop)\n";
  std::cout << "setup_ms=" << setup_ms << '\n';
  std::cout << "periodic_stage_dispatch_ms=" << periodic_ms << " ns_per_iter=" << periodic_ns_per_iter
            << '\n';
  std::cout << "isolated_stage_dispatch_ms=" << isolated_ms << " ns_per_iter=" << isolated_ns_per_iter
            << '\n';
  return 0;
}
