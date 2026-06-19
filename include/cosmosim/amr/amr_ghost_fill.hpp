#pragma once

#include <cstddef>
#include <cstdint>
#include <span>

#include "cosmosim/amr/amr_hydro_geometry.hpp"
#include "cosmosim/hydro/hydro_core_solver.hpp"

namespace cosmosim::amr {

enum class AmrGhostSourceResidency : std::uint8_t {
  kLocal,
  kRemoteReadOnly,
};

struct AmrHydroGhostFillPatch {
  AmrHydroPatchGeometry* geometry = nullptr;
  hydro::HydroConservedStateSoa* conserved = nullptr;
  AmrGhostSourceResidency residency = AmrGhostSourceResidency::kLocal;
  std::uint64_t ghost_hydro_epoch = 0;
  std::uint64_t expected_ghost_hydro_epoch = 0;
  // The finite-volume solver consumes ghost states before MUSCL-Hancock
  // reconstruction.  ghost_fill_time_code is therefore the step-start time
  // for the patch update, not an implicit predictor half-time.
  double target_state_time_code = 0.0;
  double ghost_fill_time_code = 0.0;
  double source_current_state_time_code = 0.0;
  const core::AmrTemporalBoundaryHistoryStore* temporal_boundary_history = nullptr;
  bool enable_temporal_coarse_to_fine = false;
  bool requires_ghost_fill = true;
};

struct AmrHydroGhostFillDiagnostics {
  std::size_t physical_ghosts_filled = 0;
  std::size_t same_level_ghosts_filled = 0;
  std::size_t coarse_to_fine_ghosts_filled = 0;
  std::size_t fine_to_coarse_ghosts_filled = 0;
  std::size_t temporal_coarse_to_fine_ghosts_filled = 0;
  std::size_t temporal_endpoint_ghosts_filled = 0;
  std::size_t skipped_remote_ghosts = 0;
  std::size_t stale_epoch_rejections = 0;
  std::size_t temporal_same_level_mismatch_rejections = 0;
  std::size_t temporal_history_missing_rejections = 0;
  std::size_t temporal_history_invalid_rejections = 0;
  std::size_t temporal_time_out_of_range_rejections = 0;
  std::size_t temporal_fine_to_coarse_misalignment_rejections = 0;
  std::size_t temporal_geometry_mismatch_rejections = 0;
  std::size_t temporal_identity_mismatch_rejections = 0;
  std::size_t missing_source_records = 0;
  std::size_t unresolved_ghosts = 0;
};

[[nodiscard]] std::uint64_t amrPatchGeometryFingerprint(const PatchDescriptor& patch);

// Capture the start/end conserved states for coarse patches participating in a
// local AMR subcycle.  The history is persistent state: callers may checkpoint
// after the coarse advance and resume fine substeps using the same records.
void captureAmrTemporalBoundaryHistoryStart(
    core::SimulationState& state,
    std::span<const PatchDescriptor> patches,
    double interval_start_code,
    double adiabatic_index);
void captureAmrTemporalBoundaryHistoryEnd(
    core::SimulationState& state,
    std::span<const PatchDescriptor> patches,
    double interval_end_code,
    double adiabatic_index);
void retireAmrTemporalBoundaryHistory(core::SimulationState& state);


[[nodiscard]] AmrHydroGhostFillDiagnostics fillAmrHydroGhostCells(
    std::span<AmrHydroGhostFillPatch> patches,
    double adiabatic_index);

}  // namespace cosmosim::amr
