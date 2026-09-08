#pragma once

#include <cstdint>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

#include "cosmosim/core/memory_accounting.hpp"
#include "cosmosim/core/memory_governor.hpp"

namespace cosmosim::core {

// Plans only physical replacement allocations, not logical population growth.
// The caller supplies the existing baseline owner and must exclude governed
// commitments. This is a phase-boundary utility, never a hot-loop allocator.
class RetainedCapacityTransaction final {
 public:
  using BaselineMeasure = std::function<std::uint64_t()>;

  explicit RetainedCapacityTransaction(BaselineMeasure measure)
      : m_measure(std::move(measure)) {}

  template <class TContainer>
  void add(TContainer& values, std::size_t target_capacity) {
    static_assert(std::is_nothrow_swappable_v<TContainer>,
                  "retained replacement requires nonthrowing container swap");
    if (m_consumed) throw std::logic_error("retained capacity plan is already consumed");
    if (target_capacity <= values.capacity()) return;
    if (target_capacity > values.max_size()) {
      throw std::length_error("retained capacity target exceeds container max_size");
    }
    const std::uint64_t bytes = checkedBytes(target_capacity,
        sizeof(typename TContainer::value_type));
    m_replacement_bytes = checkedMemoryBytesAdd(
        m_replacement_bytes, bytes, "retained capacity replacement total");
    m_actions.emplace_back([&values, target_capacity]() {
      // Construct a replacement before touching the original. Numeric SoA
      // lanes have a nonthrowing swap. Reject an allocator that grants more
      // than the declared capacity before that replacement becomes live.
      TContainer replacement;
      replacement.reserve(target_capacity);
      if (replacement.capacity() > target_capacity) {
        throw std::length_error("retained capacity allocator exceeded declared replacement bound");
      }
      replacement.insert(replacement.end(), values.begin(), values.end());
      values.swap(replacement);
    });
  }

  [[nodiscard]] std::uint64_t replacementBytes() const noexcept {
    return m_replacement_bytes;
  }

  void execute(MemoryGovernor* governor, MemoryClass memory_class,
               std::string_view owner) {
    if (m_consumed) throw std::logic_error("retained capacity plan is already consumed");
    m_consumed = true;
    if (m_actions.empty()) return;
    if (governor == nullptr) {
      for (auto& action : m_actions) action();
      return;
    }
    const std::uint64_t owned_before = m_measure();
    const auto snapshot = governor->snapshot();
    if (snapshot.baseline_owned_bytes < owned_before) {
      throw std::logic_error("retained capacity transaction has stale governor baseline");
    }
    auto reservation = governor->reserve(memory_class, m_replacement_bytes, owner);
    reservation.commit();
    const auto reconcile = [&]() {
      const std::uint64_t owned_after = m_measure();
      std::uint64_t baseline_after = snapshot.baseline_owned_bytes;
      if (owned_after >= owned_before) {
        baseline_after = checkedMemoryBytesAdd(baseline_after,
            owned_after - owned_before, "retained capacity baseline growth");
      } else {
        const std::uint64_t released = owned_before - owned_after;
        if (released > baseline_after) {
          throw std::logic_error("retained capacity baseline release exceeds process baseline");
        }
        baseline_after -= released;
      }
      reservation.reconcileBaselineOwnedAndRelease(baseline_after);
    };
    try {
      for (auto& action : m_actions) action();
    } catch (...) {
      // Earlier lanes may already have grown. Retain their physical accounting
      // even when a later allocation fails; never release their bytes as scratch.
      reconcile();
      throw;
    }
    reconcile();
  }

 private:
  [[nodiscard]] static std::uint64_t checkedBytes(std::size_t count,
                                                   std::size_t width) {
    if (width != 0U && count > std::numeric_limits<std::uint64_t>::max() / width) {
      throw std::overflow_error("retained capacity byte multiplication overflow");
    }
    return static_cast<std::uint64_t>(count) * static_cast<std::uint64_t>(width);
  }
  BaselineMeasure m_measure;
  std::uint64_t m_replacement_bytes = 0U;
  bool m_consumed = false;
  std::vector<std::function<void()>> m_actions;
};

}  // namespace cosmosim::core
