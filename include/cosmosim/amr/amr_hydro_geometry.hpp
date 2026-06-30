#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include "cosmosim/amr/amr_framework.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/hydro/hydro_core_solver.hpp"

namespace cosmosim::amr {

enum class AmrHydroBoundaryClass : std::uint8_t {
  kPhysical,
  kSameLevel,
  kCoarseFine,
};

enum class AmrHydroGhostFillStatus : std::uint8_t {
  kUnfilledPhysicalBoundary,
  kUnfilledSameLevel,
  kUnfilledCoarseFine,
  kFilledPhysicalBoundary,
  kFilledSameLevel,
  kFilledCoarseToFine,
  kFilledCoarseToFineTemporal,
  kFilledFineToCoarse,
  kSkippedRemote,
  kRejectedStaleRemote,
  kRejectedTemporalSameLevelMismatch,
  kRejectedTemporalHistoryMissing,
  kRejectedTemporalHistoryInvalid,
  kRejectedTemporalTimeOutOfRange,
  kRejectedTemporalFineToCoarseMismatch,
  kMissingSource,
};

struct AmrHydroGeometryOptions {
  hydro::HydroBoundaryKind physical_boundary_kind = hydro::HydroBoundaryKind::kOpen;
  std::array<AmrHydroBoundaryClass, 6> boundary_classes = {
      AmrHydroBoundaryClass::kPhysical,
      AmrHydroBoundaryClass::kPhysical,
      AmrHydroBoundaryClass::kPhysical,
      AmrHydroBoundaryClass::kPhysical,
      AmrHydroBoundaryClass::kPhysical,
      AmrHydroBoundaryClass::kPhysical,
  };
};

struct AmrHydroCellDescriptor {
  std::uint64_t patch_id = 0;
  std::uint64_t gas_cell_id = 0;
  std::uint32_t local_cell_row = core::kInvalidGasCellRow;
  std::size_t patch_local_cell = hydro::k_invalid_cell_index;
};

struct AmrHydroGhostDescriptor {
  std::uint64_t patch_id = 0;
  std::uint64_t source_patch_id = 0;
  std::uint64_t source_gas_cell_id = 0;
  std::uint64_t owner_gas_cell_id = 0;
  std::uint32_t owner_local_cell_row = core::kInvalidGasCellRow;
  std::size_t ghost_slot = hydro::k_invalid_ghost_cell_slot;
  std::size_t ghost_cell = hydro::k_invalid_cell_index;
  AmrHydroBoundaryClass boundary_class = AmrHydroBoundaryClass::kPhysical;
  AmrHydroGhostFillStatus fill_status = AmrHydroGhostFillStatus::kUnfilledPhysicalBoundary;
};

struct AmrHydroFaceDescriptor {
  std::uint64_t patch_id = 0;
  std::uint64_t face_id = 0;
  std::uint64_t owner_gas_cell_id = 0;
  std::uint64_t neighbor_gas_cell_id = 0;
  std::uint64_t ghost_source_gas_cell_id = 0;
  std::size_t owner_patch_cell = hydro::k_invalid_cell_index;
  std::size_t neighbor_patch_cell = hydro::k_invalid_cell_index;
  std::size_t ghost_slot = hydro::k_invalid_ghost_cell_slot;
  double area_comoving = 0.0;
  double normal_x = 0.0;
  double normal_y = 0.0;
  double normal_z = 0.0;
  hydro::HydroFaceAxis axis = hydro::HydroFaceAxis::kX;
  AmrHydroBoundaryClass boundary_class = AmrHydroBoundaryClass::kPhysical;
};

struct AmrHydroPatchGeometry {
  PatchDescriptor patch;
  hydro::HydroPatchGeometry geometry;
  std::vector<AmrHydroCellDescriptor> real_cells;
  std::vector<std::uint64_t> gas_cell_ids;
  std::vector<std::uint32_t> local_cell_rows;
  std::vector<AmrHydroGhostDescriptor> ghosts;
  std::vector<AmrHydroFaceDescriptor> faces;
  std::uint64_t source_gas_cell_identity_generation = 0;

  [[nodiscard]] std::span<const std::uint64_t> gasCellIds() const noexcept;
  [[nodiscard]] std::vector<std::size_t> internalFaceIndices() const;
};

[[nodiscard]] AmrHydroPatchGeometry buildAmrHydroPatchGeometry(
    const core::SimulationState& state,
    const PatchDescriptor& patch,
    const AmrHydroGeometryOptions& options = {});

[[nodiscard]] AmrHydroPatchGeometry buildRemoteAmrHydroPatchGeometry(
    const PatchDescriptor& patch,
    std::span<const std::uint64_t> gas_cell_ids,
    std::uint64_t source_gas_cell_identity_generation,
    const AmrHydroGeometryOptions& options = {});

[[nodiscard]] hydro::HydroConservedStateSoa loadAmrHydroConservedState(
    const core::SimulationState& state,
    const AmrHydroPatchGeometry& patch_geometry,
    double adiabatic_index);

}  // namespace cosmosim::amr
