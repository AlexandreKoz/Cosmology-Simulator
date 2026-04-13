#include <chrono>
#include <cstdint>
#include <iostream>

#include "cosmosim/core/time_integration.hpp"

int main() {
  constexpr std::uint32_t k_element_count = 500000;
  constexpr std::uint8_t k_max_bin = 6;
  constexpr std::uint32_t k_substeps = 4096;

  cosmosim::core::HierarchicalTimeBinScheduler scheduler(k_max_bin);
  scheduler.reset(k_element_count, k_max_bin, 0);

  for (std::uint32_t i = 0; i < k_element_count; ++i) {
    scheduler.setElementBin(i, static_cast<std::uint8_t>(i % (k_max_bin + 1)), 0);
  }

  std::uint64_t active_checksum = 0;
  std::uint64_t touch_checksum = 0;

  const auto begin = std::chrono::steady_clock::now();
  for (std::uint32_t step = 0; step < k_substeps; ++step) {
    const auto active = scheduler.beginSubstep();
    active_checksum += active.size();

    for (const std::uint32_t element : active) {
      touch_checksum += element;
      if ((element & 15U) == 0U) {
        const auto current_bin = scheduler.hotMetadata().bin_index[element];
        const auto target = static_cast<std::uint8_t>((current_bin == 0U) ? 1U : current_bin - 1U);
        scheduler.requestBinTransition(element, target);
      }
    }

    scheduler.endSubstep();
  }
  const auto end = std::chrono::steady_clock::now();

  const auto total_us = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
  const double active_per_us = static_cast<double>(active_checksum) / static_cast<double>(total_us);

  const auto& diagnostics = scheduler.diagnostics();

  std::cout << "build_type="
#ifdef NDEBUG
            << "release";
#else
            << "debug";
#endif
  std::cout << '\n';
  std::cout << "elements=" << k_element_count << '\n';
  std::cout << "max_bin=" << static_cast<int>(k_max_bin) << '\n';
  std::cout << "substeps=" << k_substeps << '\n';
  std::cout << "elapsed_us=" << total_us << '\n';
  std::cout << "active_checksum=" << active_checksum << '\n';
  std::cout << "touch_checksum=" << touch_checksum << '\n';
  std::cout << "active_per_us=" << active_per_us << '\n';
  std::cout << "promoted=" << diagnostics.promoted_elements << '\n';
  std::cout << "demoted=" << diagnostics.demoted_elements << '\n';
  std::cout << "illegal_transitions=" << diagnostics.illegal_transition_attempts << '\n';

  return 0;
}
