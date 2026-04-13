#include <chrono>
#include <cstdint>
#include <iostream>
#include <string>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/core/profiling.hpp"

namespace {

void runScopes(cosmosim::core::ProfilerSession& session, int iterations) {
  for (int i = 0; i < iterations; ++i) {
    cosmosim::core::ScopedProfile step_scope(&session, "step");
    session.addBytesMoved(256);
    {
      cosmosim::core::ScopedProfile gravity_scope(&session, "gravity.tree");
      session.counters().addCount("tree.node_visits", 32);
      session.addBytesMoved(512);
    }
    {
      cosmosim::core::ScopedProfile hydro_scope(&session, "hydro.face_flux");
      session.counters().addCount("hydro.face_fluxes", 64);
      session.addBytesMoved(768);
    }
  }
}

}  // namespace

int main() {
  constexpr int k_iterations = 150000;
  constexpr int k_warmup = 1000;

  cosmosim::core::ProfilerSession disabled_session(false);
  cosmosim::core::ProfilerSession enabled_session(true);

  runScopes(disabled_session, k_warmup);
  runScopes(enabled_session, k_warmup);

  const auto setup_start = std::chrono::steady_clock::now();
  disabled_session.reset();
  enabled_session.reset();
  const auto setup_end = std::chrono::steady_clock::now();

  const auto disabled_start = std::chrono::steady_clock::now();
  runScopes(disabled_session, k_iterations);
  const auto disabled_end = std::chrono::steady_clock::now();

  const auto enabled_start = std::chrono::steady_clock::now();
  runScopes(enabled_session, k_iterations);
  const auto enabled_end = std::chrono::steady_clock::now();

  const double setup_ms =
      std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(setup_end - setup_start).count();
  const double disabled_ms =
      std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(disabled_end - disabled_start).count();
  const double enabled_ms =
      std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(enabled_end - enabled_start).count();

  const double overhead_ms = enabled_ms - disabled_ms;
  const double overhead_pct = (disabled_ms > 0.0) ? (100.0 * overhead_ms / disabled_ms) : 0.0;
  const double scopes = static_cast<double>(k_iterations) * 3.0;
  const double enabled_scope_rate_mscope_s = scopes / (enabled_ms * 1.0e-3) * 1.0e-6;

  std::cout << "benchmark=profiling_overhead"
            << " build_type=" << COSMOSIM_BUILD_TYPE
            << " profiling_compile_time_enabled=" << COSMOSIM_ENABLE_PROFILING
            << " threads=1"
            << " gpu_events=disabled"
            << " iterations=" << k_iterations
            << " setup_ms=" << setup_ms
            << " disabled_ms=" << disabled_ms
            << " enabled_ms=" << enabled_ms
            << " overhead_ms=" << overhead_ms
            << " overhead_percent=" << overhead_pct
            << " enabled_scope_rate_mscope_s=" << enabled_scope_rate_mscope_s
            << " tree_node_visits=" << enabled_session.counters().count("tree.node_visits")
            << " hydro_face_fluxes=" << enabled_session.counters().count("hydro.face_fluxes")
            << '\n';

  return 0;
}
