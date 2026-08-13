#pragma once

#include <cstdint>
#include <limits>
#include <stdexcept>

namespace cosmosim::gravity {

// Strong gravity-state identities deliberately do not implicitly convert to or
// from uint64_t.  Source mutation, tree rebuilds, PM field refreshes, force
// evaluation, and MPI ownership changes are independent state transitions even
// when the rung-zero production scheduler currently advances some of them at
// the same boundary.
template <typename TTag>
struct GravityIdentity {
  std::uint64_t value = 0U;

  [[nodiscard]] constexpr bool valid() const noexcept { return value != 0U; }
  [[nodiscard]] constexpr bool operator==(const GravityIdentity&) const noexcept = default;
};

struct GravitySourceGenerationTag;
struct TreeBuildGenerationTag;
struct PmFieldVersionTag;
struct DecompositionEpochTag;

using GravitySourceGeneration = GravityIdentity<GravitySourceGenerationTag>;
using TreeBuildGeneration = GravityIdentity<TreeBuildGenerationTag>;
using PmFieldVersion = GravityIdentity<PmFieldVersionTag>;
using DecompositionEpoch = GravityIdentity<DecompositionEpochTag>;

struct ForceEvaluationEpoch {
  std::uint64_t sequence = 0U;
  double scale_factor = 1.0;

  [[nodiscard]] constexpr bool valid() const noexcept {
    return sequence != 0U && scale_factor > 0.0;
  }
  [[nodiscard]] constexpr bool operator==(const ForceEvaluationEpoch&) const noexcept = default;
};

// Advance helpers centralize overflow behavior for identities owned by a local
// component.  Scheduler-owned identities are normally constructed from their
// authoritative scheduler counters instead.
template <typename TIdentity>
[[nodiscard]] TIdentity nextGravityIdentity(TIdentity current, const char* context) {
  if (current.value == std::numeric_limits<std::uint64_t>::max()) {
    throw std::overflow_error(context);
  }
  return TIdentity{current.value + 1U};
}

}  // namespace cosmosim::gravity
