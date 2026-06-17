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
};

struct AmrHydroGhostFillDiagnostics {
  std::size_t physical_ghosts_filled = 0;
  std::size_t same_level_ghosts_filled = 0;
  std::size_t coarse_to_fine_ghosts_filled = 0;
  std::size_t fine_to_coarse_ghosts_filled = 0;
  std::size_t skipped_remote_ghosts = 0;
  std::size_t stale_epoch_rejections = 0;
  std::size_t missing_source_records = 0;
  std::size_t unresolved_ghosts = 0;
};

[[nodiscard]] AmrHydroGhostFillDiagnostics fillAmrHydroGhostCells(
    std::span<AmrHydroGhostFillPatch> patches,
    double adiabatic_index);

}  // namespace cosmosim::amr
