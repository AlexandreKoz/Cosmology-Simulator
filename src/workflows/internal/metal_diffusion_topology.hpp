#pragma once

#include <cstdint>
#include <span>
#include <vector>

#include "cosmosim/core/simulation_mode.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/physics/metal_diffusion.hpp"

namespace cosmosim::workflows::internal {

struct MetalDiffusionTopologyResult {
  std::vector<physics::MetalDiffusionFace> faces;
  std::uint64_t interior_patch_face_count = 0U;
  std::uint64_t cross_patch_face_count = 0U;
  std::uint64_t periodic_wrap_face_count = 0U;
};

// Construct the authoritative local diffusion interface graph from AMR patch
// geometry. Faces are oriented from the negative-axis cell to the positive-axis
// cell, including periodic wraps. The same graph is used to reconstruct the
// velocity gradients that set the Smagorinsky diffusivity, so patch tiling does
// not silently change the local closure.
[[nodiscard]] MetalDiffusionTopologyResult buildMetalDiffusionTopology(
    const core::SimulationState& state,
    std::span<const std::uint8_t> owned_leaf_mask,
    std::uint32_t world_rank,
    double scale_factor,
    core::BoundaryCondition hydro_boundary,
    std::span<physics::MetalDiffusionCell> diffusion_cells);

}  // namespace cosmosim::workflows::internal
