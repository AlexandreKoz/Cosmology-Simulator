#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "cosmosim/core/config.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/hydro/hydro_cartesian_patch.hpp"

namespace cosmosim::workflows::internal {

// Fixed-patch compatibility layout. Geometry order is physical Cartesian order;
// dense CellSoa rows are transient storage locations only.
struct CartesianGasCellRowLayout {
  hydro::HydroCartesianPatchSpec spec;
  std::vector<std::uint32_t> dense_row_by_geometry_row;
  std::vector<std::uint32_t> geometry_row_by_dense_row;
  std::uint64_t mapping_signature = 0;
};

struct CartesianGasCellLayoutBuildResult {
  CartesianGasCellRowLayout layout;
  std::string diagnostic;

  [[nodiscard]] bool ok() const noexcept { return diagnostic.empty(); }
};

// Construct a fixed-patch layout only when physical Cartesian topology can be
// proven from complete explicit patch metadata or the full gas-cell center set.
// This routine deliberately never factors cell_count into an invented mesh.
[[nodiscard]] CartesianGasCellLayoutBuildResult buildCartesianGasCellRowLayout(
    const core::SimulationState& state,
    const core::SimulationConfig& config);

}  // namespace cosmosim::workflows::internal
