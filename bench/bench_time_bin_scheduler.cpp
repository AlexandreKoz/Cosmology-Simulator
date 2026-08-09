#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <string_view>
#include <vector>

#include "cosmosim/core/time_integration.hpp"

namespace {

using Clock = std::chrono::steady_clock;
using Scheduler = cosmosim::core::HierarchicalTimeBinScheduler;

struct Scenario {
  std::string_view name;
  std::uint32_t elements;
  std::uint8_t max_bin;
  std::uint32_t substeps;
  bool dynamic_transitions;
};

std::vector<std::vector<std::uint32_t>> buildMirrorBins(
    std::uint32_t element_count,
    std::uint8_t max_bin) {
  std::vector<std::vector<std::uint32_t>> bins(static_cast<std::size_t>(max_bin) + 1U);
  bins[max_bin].reserve(element_count);
  for (std::uint32_t element = 0; element < element_count; ++element) {
    bins[max_bin].push_back(element);
  }
  std::vector<std::size_t> position(element_count);
  std::iota(position.begin(), position.end(), 0U);

  // Reproduce reset(max_bin) followed by setElementBin(i, i % bins), including
  // the scheduler's swap-with-last erase ordering. This gives the exact pre-sort
  // member order for the static benchmark scenario.
  for (std::uint32_t element = 0; element < element_count; ++element) {
    const auto target = static_cast<std::uint8_t>(element % (static_cast<std::uint32_t>(max_bin) + 1U));
    if (target == max_bin) {
      continue;
    }
    auto& source = bins[max_bin];
    const std::size_t remove_pos = position[element];
    const std::uint32_t last = source.back();
    source[remove_pos] = last;
    position[last] = remove_pos;
    source.pop_back();
    position[element] = bins[target].size();
    bins[target].push_back(element);
  }
  return bins;
}

std::vector<std::uint32_t> preSortActiveForTick(
    const std::vector<std::vector<std::uint32_t>>& bins,
    std::uint64_t tick) {
  std::vector<std::uint32_t> active;
  std::size_t expected = 0;
  for (std::size_t bin = 0; bin < bins.size(); ++bin) {
    const std::uint64_t period = std::uint64_t{1} << bin;
    if (tick % period == 0U) {
      expected += bins[bin].size();
    }
  }
  active.reserve(expected);
  for (std::size_t bin = 0; bin < bins.size(); ++bin) {
    const std::uint64_t period = std::uint64_t{1} << bin;
    if (tick % period != 0U) {
      continue;
    }
    active.insert(active.end(), bins[bin].begin(), bins[bin].end());
  }
  return active;
}

void runScenario(const Scenario& scenario) {
  Scheduler scheduler(scenario.max_bin);
  scheduler.reset(scenario.elements, scenario.max_bin, 0U);
  for (std::uint32_t i = 0; i < scenario.elements; ++i) {
    scheduler.setElementBin(
        i,
        static_cast<std::uint8_t>(i % (static_cast<std::uint32_t>(scenario.max_bin) + 1U)),
        0U);
  }

  const auto mirror_bins = scenario.dynamic_transitions
      ? std::vector<std::vector<std::uint32_t>>{}
      : buildMirrorBins(scenario.elements, scenario.max_bin);

  std::uint64_t active_checksum = 0;
  std::uint64_t touch_checksum = 0;
  std::uint64_t begin_ns = 0;
  std::uint64_t transition_ns = 0;
  std::uint64_t isolated_sort_ns = 0;
  std::uint64_t isolated_sort_elements = 0;

  for (std::uint32_t step = 0; step < scenario.substeps; ++step) {
    const std::uint64_t tick = scheduler.currentTick();
    if (!scenario.dynamic_transitions) {
      auto unsorted = preSortActiveForTick(mirror_bins, tick);
      isolated_sort_elements += unsorted.size();
      const auto sort_begin = Clock::now();
      std::sort(unsorted.begin(), unsorted.end());
      const auto sort_end = Clock::now();
      isolated_sort_ns += static_cast<std::uint64_t>(
          std::chrono::duration_cast<std::chrono::nanoseconds>(sort_end - sort_begin).count());
    }

    const auto begin_start = Clock::now();
    const auto active = scheduler.beginSubstep();
    const auto begin_end = Clock::now();
    begin_ns += static_cast<std::uint64_t>(
        std::chrono::duration_cast<std::chrono::nanoseconds>(begin_end - begin_start).count());
    active_checksum += active.size();

    const auto transition_begin = Clock::now();
    for (const std::uint32_t element : active) {
      touch_checksum += element;
      if (scenario.dynamic_transitions && (element & 63U) == 0U) {
        const auto current_bin = scheduler.hotMetadata().bin_index[element];
        const auto target = static_cast<std::uint8_t>(
            (current_bin == 0U) ? std::min<std::uint8_t>(1U, scenario.max_bin) : current_bin - 1U);
        scheduler.submitCandidateBin(
            element,
            target,
            cosmosim::core::TimeStepCandidateSource::kUserClamp);
      }
    }
    scheduler.endSubstep();
    const auto transition_end = Clock::now();
    transition_ns += static_cast<std::uint64_t>(
        std::chrono::duration_cast<std::chrono::nanoseconds>(transition_end - transition_begin).count());
  }

  const double begin_ms = static_cast<double>(begin_ns) / 1.0e6;
  const double transition_ms = static_cast<double>(transition_ns) / 1.0e6;
  const double sort_ms = static_cast<double>(isolated_sort_ns) / 1.0e6;
  std::cout << "benchmark=time_bin_scheduler"
            << " scenario=" << scenario.name
            << " build_type="
#ifdef NDEBUG
            << "release"
#else
            << "debug"
#endif
            << " elements=" << scenario.elements
            << " max_bin=" << static_cast<int>(scenario.max_bin)
            << " substeps=" << scenario.substeps
            << " active_checksum=" << active_checksum
            << " begin_substep_ms=" << begin_ms
            << " transition_end_ms=" << transition_ms
            << " comparison_sort_reference_ms=" << sort_ms
            << " comparison_sort_reference_elements=" << isolated_sort_elements
            << " comparison_sort_to_radix_begin_ratio=" << (begin_ns > 0U ? static_cast<double>(isolated_sort_ns) / static_cast<double>(begin_ns) : 0.0)
            << " touch_checksum=" << touch_checksum
            << '\n';
}

}  // namespace

int main() {
  const Scenario scenarios[] = {
      {"workstation_static_100k", 100000U, 8U, 256U, false},
      {"workstation_static_500k", 500000U, 8U, 128U, false},
      {"late_time_fine_bins_500k", 500000U, 12U, 128U, false},
      {"dynamic_transitions_250k", 250000U, 8U, 256U, true},
  };
  for (const Scenario& scenario : scenarios) {
    runScenario(scenario);
  }
  return 0;
}
