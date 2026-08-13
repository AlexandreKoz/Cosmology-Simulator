#pragma once

#include <cstdint>
#include <string_view>
#include <vector>

#include "cosmosim/core/simulation_state.hpp"

namespace cosmosim::workflows::internal {

struct AuthoritativeGravitySourceRows {
  std::vector<std::uint32_t> particle_rows;
  std::vector<std::uint32_t> gas_cell_rows;
};

// Selects only authoritative local gravity source rows. Collisionless and
// compact-object particles remain particle-owned. Generic gas-tagged particles
// are compatibility/lineage mirrors and are excluded; owned leaf gas cells are
// the gas source of truth. Covered coarse cells are omitted.
[[nodiscard]] AuthoritativeGravitySourceRows selectAuthoritativeGravitySourceRows(
    const core::SimulationState& state,
    std::uint32_t world_rank,
    std::uint32_t world_size,
    std::string_view caller);

}  // namespace cosmosim::workflows::internal
