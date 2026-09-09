#include <array>
#include <cassert>
#include <cstdint>
#include <limits>
#include <stdexcept>

#include "cosmosim/core/governed_scratch_arena.hpp"
#include "../../src/io/internal/sidecar_row_lookup.hpp"

int main() {
  using namespace cosmosim;
  const auto bytes = io::internal::sidecarLookupWorkspaceBytes(3U, 0U, 0U);
  assert(bytes == 3U * sizeof(std::uint64_t) + 256U);
  core::MemoryGovernor governor(core::MemoryGovernorPolicy{.hard_limit_bytes = bytes});
  {
    core::GovernedScratchArena arena(&governor, core::MemoryClass::kDiagnostic,
                                     bytes, "unit.sidecar_lookup");
    const std::array<std::uint32_t,3> indices{9U, 1U, 5U};
    io::internal::SidecarRowLookup lookup(indices, arena.resource(), "stellar");
    assert(lookup.rowFor(9U) == 0U);
    assert(lookup.rowFor(1U) == 1U);
    assert(lookup.rowFor(5U) == 2U);
    bool missing = false;
    try { (void)lookup.rowFor(7U); }
    catch (const std::runtime_error&) { missing = true; }
    assert(missing);
  }
  assert(governor.snapshot().committed_bytes == 0U);
  bool rejected = false;
  try {
    core::GovernedScratchArena arena(&governor, core::MemoryClass::kDiagnostic,
                                     bytes + 1U, "unit.over_limit");
  } catch (const core::MemoryAdmissionError&) { rejected = true; }
  assert(rejected);
  {
    core::GovernedScratchArena arena(nullptr, core::MemoryClass::kDiagnostic,
                                     bytes, "unit.duplicate");
    bool duplicate = false;
    try {
      io::internal::SidecarRowLookup lookup(
          std::array<std::uint32_t,3>{9U,1U,9U}, arena.resource(), "stellar");
    } catch (const std::runtime_error&) { duplicate = true; }
    assert(duplicate);
  }
  bool overflow = false;
  try {
    (void)io::internal::sidecarLookupWorkspaceBytes(
        std::numeric_limits<std::uint64_t>::max(), 1U, 0U);
  } catch (const std::overflow_error&) { overflow = true; }
  assert(overflow);
  assert(io::internal::sidecarLookupWorkspaceBytes(0U,0U,0U) == 0U);
}
