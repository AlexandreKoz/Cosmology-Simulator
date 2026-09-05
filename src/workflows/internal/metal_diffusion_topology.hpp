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
  std::vector<double> strain_magnitude_code;
  std::uint64_t construction_scratch_high_water_bytes = 0U;
  std::uint64_t interior_patch_face_count = 0U;
  std::uint64_t cross_patch_face_count = 0U;
  std::uint64_t periodic_wrap_face_count = 0U;
};

// Construct the authoritative local diffusion interface graph from AMR patch
// geometry. Faces are oriented from the negative-axis cell to the positive-axis
// cell, including periodic wraps. The same graph is used to reconstruct the
// trace-free strain magnitude that sets the Smagorinsky diffusivity. Gradient
// accumulation is cell-local scratch; no 3x3 gradient is retained per gas cell.
[[nodiscard]] MetalDiffusionTopologyResult buildMetalDiffusionTopology(
    const core::SimulationState& state,
    std::span<const std::uint8_t> owned_leaf_mask,
    std::uint32_t world_rank,
    double scale_factor,
    core::BoundaryCondition hydro_boundary);

}  // namespace cosmosim::workflows::internal
