#pragma once

#include <cstdint>
#include <optional>
#include <string_view>
#include <unordered_map>

#include "cosmosim/core/simulation_state.hpp"

namespace cosmosim::workflows::internal {

using ParticleRowById = std::unordered_map<std::uint64_t, std::uint32_t>;

[[nodiscard]] ParticleRowById buildParticleRowById(
    const core::SimulationState& state);

[[nodiscard]] std::optional<std::uint32_t> parentParticleRowForGasCellRow(
    const core::SimulationState& state,
    std::uint32_t cell_row,
    const ParticleRowById& particle_row_by_id,
    std::string_view caller);

[[nodiscard]] const core::GasCellIdentityRecord&
gasCellIdentityRecordForLocalRow(
    const core::SimulationState& state,
    std::uint32_t cell_row,
    std::string_view caller);

[[nodiscard]] std::uint32_t gasCellOwnerRankForLocalRow(
    const core::SimulationState& state,
    std::uint32_t cell_row,
    const ParticleRowById& particle_row_by_id,
    std::string_view caller);

// Gas-cell state is authoritative. This derives deterministic parent-particle
// compatibility mirrors by stable gas_cell_id after hydro/AMR mutations.
void synchronizeParentParticleCompatibilityMirrors(
    core::SimulationState& state,
    std::uint32_t world_rank,
    std::string_view caller);

}  // namespace cosmosim::workflows::internal
