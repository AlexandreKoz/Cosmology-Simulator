#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include "cosmosim/hydro/hydro_reconstruction.hpp"
#include "cosmosim/hydro/hydro_riemann.hpp"

namespace cosmosim::hydro {

constexpr std::size_t k_invalid_cell_index = static_cast<std::size_t>(-1);
constexpr std::size_t k_invalid_ghost_cell_slot = static_cast<std::size_t>(-1);

// Cell-centered primitive variables in comoving coordinates.
// rho_comoving: comoving mass density.
// vel_*_peculiar: peculiar velocity components.
// pressure_comoving: comoving-frame thermal pressure.
struct HydroPrimitiveState {
  double rho_comoving = 0.0;
  double vel_x_peculiar = 0.0;
  double vel_y_peculiar = 0.0;
  double vel_z_peculiar = 0.0;
  double pressure_comoving = 0.0;
};

// Cell-centered conserved variables in comoving coordinates.
// mass_density_comoving: rho.
// momentum_density_*_comoving: rho * u.
// total_energy_density_comoving: rho * (e_internal + 0.5 |u|^2).
struct HydroConservedState {
  double mass_density_comoving = 0.0;
  double momentum_density_x_comoving = 0.0;
  double momentum_density_y_comoving = 0.0;
  double momentum_density_z_comoving = 0.0;
  double total_energy_density_comoving = 0.0;

  HydroConservedState& operator+=(const HydroConservedState& rhs);
  HydroConservedState& operator-=(const HydroConservedState& rhs);
};

[[nodiscard]] HydroConservedState operator+(HydroConservedState lhs, const HydroConservedState& rhs);
[[nodiscard]] HydroConservedState operator-(HydroConservedState lhs, const HydroConservedState& rhs);
[[nodiscard]] HydroConservedState operator*(double scalar, HydroConservedState state);

class HydroConservedStateSoa {
 public:
  explicit HydroConservedStateSoa(std::size_t cell_count = 0);

  void resize(std::size_t cell_count);
  [[nodiscard]] std::size_t size() const;

  [[nodiscard]] HydroConservedState loadCell(std::size_t cell_index) const;
  void storeCell(std::size_t cell_index, const HydroConservedState& cell_state);

  [[nodiscard]] std::span<double> massDensityComoving();
  [[nodiscard]] std::span<const double> massDensityComoving() const;

  [[nodiscard]] std::span<double> momentumDensityXComoving();
  [[nodiscard]] std::span<const double> momentumDensityXComoving() const;

  [[nodiscard]] std::span<double> momentumDensityYComoving();
  [[nodiscard]] std::span<const double> momentumDensityYComoving() const;

  [[nodiscard]] std::span<double> momentumDensityZComoving();
  [[nodiscard]] std::span<const double> momentumDensityZComoving() const;

  [[nodiscard]] std::span<double> totalEnergyDensityComoving();
  [[nodiscard]] std::span<const double> totalEnergyDensityComoving() const;

 private:
  std::vector<double> m_mass_density_comoving;
  std::vector<double> m_momentum_density_x_comoving;
  std::vector<double> m_momentum_density_y_comoving;
  std::vector<double> m_momentum_density_z_comoving;
  std::vector<double> m_total_energy_density_comoving;
};

class HydroPrimitiveCacheSoa {
 public:
  explicit HydroPrimitiveCacheSoa(std::size_t cell_count = 0);

  void resize(std::size_t cell_count);
  [[nodiscard]] std::size_t size() const;

  [[nodiscard]] HydroPrimitiveState loadCell(std::size_t cell_index) const;
  void storeCell(std::size_t cell_index, const HydroPrimitiveState& primitive_state);

 private:
  std::vector<double> m_rho_comoving;
  std::vector<double> m_vel_x_peculiar;
  std::vector<double> m_vel_y_peculiar;
  std::vector<double> m_vel_z_peculiar;
  std::vector<double> m_pressure_comoving;
};

struct HydroPatchColdData {
  std::uint64_t patch_id = 0;
  int refinement_level = 0;
};

enum class HydroFaceAxis {
  kX,
  kY,
  kZ,
};

enum class HydroFaceSide {
  kLower,
  kUpper,
};

enum class HydroFluxRegisterFaceRole {
  kNone,
  kCoarse,
  kFine,
};

struct HydroFluxRegisterFace {
  HydroFluxRegisterFaceRole role = HydroFluxRegisterFaceRole::kNone;
  std::uint64_t register_key = 0;
  std::uint64_t coarse_patch_id = 0;
  std::uint64_t coarse_gas_cell_id = 0;
  std::size_t coarse_cell_index = k_invalid_cell_index;
  int level = 0;
  HydroFaceAxis axis = HydroFaceAxis::kX;
  HydroFaceSide orientation = HydroFaceSide::kLower;
  double coarse_orientation_sign = 1.0;
};

enum class HydroBoundaryKind {
  kInterior,
  kPeriodic,
  kOpen,
  kReflective,
  kImportedMpi,
};

enum class HydroGhostMutationRights {
  kReadOnlyImported,
  kWritablePhysicalBoundaryScratch,
};

struct HydroGhostCell {
  std::size_t owner_real_cell = k_invalid_cell_index;
  std::size_t source_real_cell = k_invalid_cell_index;
  std::size_t ghost_cell = k_invalid_cell_index;
  std::size_t ghost_slot = k_invalid_ghost_cell_slot;
  HydroBoundaryKind boundary_kind = HydroBoundaryKind::kInterior;
  HydroFaceAxis axis = HydroFaceAxis::kX;
  HydroFaceSide side = HydroFaceSide::kLower;
  HydroGhostMutationRights mutation_rights = HydroGhostMutationRights::kWritablePhysicalBoundaryScratch;
};

struct HydroFace {
  std::size_t owner_cell = k_invalid_cell_index;
  std::size_t neighbor_cell = k_invalid_cell_index;
  std::size_t owner_minus_cell = k_invalid_cell_index;
  std::size_t neighbor_plus_cell = k_invalid_cell_index;
  std::size_t ghost_cell_slot = k_invalid_ghost_cell_slot;
  double area_comoving = 0.0;
  double normal_x = 0.0;
  double normal_y = 0.0;
  double normal_z = 0.0;
  HydroFaceAxis axis = HydroFaceAxis::kX;
};

struct HydroPatchGeometry {
  std::size_t nx = 0;
  std::size_t ny = 0;
  std::size_t nz = 0;
  double origin_x_comoving = 0.0;
  double origin_y_comoving = 0.0;
  double origin_z_comoving = 0.0;
  double cell_width_x_comoving = 0.0;
  double cell_width_y_comoving = 0.0;
  double cell_width_z_comoving = 0.0;
  double cell_volume_comoving = 0.0;
  std::vector<HydroFace> faces;
  std::vector<HydroGhostCell> ghost_cells;
  std::vector<HydroFluxRegisterFace> flux_register_faces;

  [[nodiscard]] std::size_t cellCount() const noexcept;
  [[nodiscard]] std::size_t totalCellStorageCount() const noexcept;
  [[nodiscard]] std::size_t linearCellIndex(std::size_t i, std::size_t j, std::size_t k) const;
  [[nodiscard]] std::array<std::size_t, 3> cellIjk(std::size_t row) const;
  [[nodiscard]] std::size_t neighborCell(
      std::size_t row,
      int di,
      int dj,
      int dk) const noexcept;
};

struct HydroFluxRegisterRecord {
  HydroFluxRegisterFaceRole role = HydroFluxRegisterFaceRole::kNone;
  std::uint64_t register_key = 0;
  std::uint64_t coarse_patch_id = 0;
  std::uint64_t coarse_gas_cell_id = 0;
  std::size_t coarse_cell_index = k_invalid_cell_index;
  int level = 0;
  HydroFaceAxis axis = HydroFaceAxis::kX;
  HydroFaceSide orientation = HydroFaceSide::kLower;
  double face_area_comoving = 0.0;
  double dt_code = 0.0;
  HydroConservedState flux_code;
};

class HydroFluxRegisterSink {
 public:
  virtual ~HydroFluxRegisterSink() = default;

  virtual void recordFaceFlux(const HydroFluxRegisterRecord& record) = 0;
};

struct HydroActiveSetView {
  std::span<const std::size_t> active_cells;
  std::span<const std::size_t> active_faces;
};

struct HydroUpdateContext {
  double dt_code = 0.0;
  double scale_factor = 1.0;
  double hubble_rate_code = 0.0;
};

struct HydroSourceContext {
  HydroUpdateContext update;
  std::span<const double> gravity_accel_x_peculiar;
  std::span<const double> gravity_accel_y_peculiar;
  std::span<const double> gravity_accel_z_peculiar;
  std::span<const double> hydrogen_number_density_cgs;
  std::span<const double> metallicity_mass_fraction;
  std::span<const double> temperature_k;
  double redshift = 0.0;
};

// Volume-integrated conserved diagnostics over a selected set of real cells.
// Density-like conserved lanes are multiplied by HydroPatchGeometry::cell_volume_comoving.
// internal_energy is derived as total energy minus kinetic energy, not stored as an
// independent conserved authority.
struct HydroConservationTotals {
  double mass = 0.0;
  double momentum_x = 0.0;
  double momentum_y = 0.0;
  double momentum_z = 0.0;
  double total_energy = 0.0;
  double internal_energy = 0.0;

  HydroConservationTotals& operator+=(const HydroConservationTotals& rhs);
  HydroConservationTotals& operator-=(const HydroConservationTotals& rhs);
};

[[nodiscard]] HydroConservationTotals operator+(
    HydroConservationTotals lhs,
    const HydroConservationTotals& rhs);
[[nodiscard]] HydroConservationTotals operator-(
    HydroConservationTotals lhs,
    const HydroConservationTotals& rhs);

struct HydroConservationReport {
  HydroConservationTotals before;
  HydroConservationTotals after;
  HydroConservationTotals flux_delta;
  HydroConservationTotals source_delta;
  HydroConservationTotals floor_delta;
  HydroConservationTotals residual;
  std::uint64_t cell_count = 0;
  std::uint64_t internal_energy_floor_count = 0;
};

class HydroSourceTerm {
 public:
  virtual ~HydroSourceTerm() = default;

  [[nodiscard]] virtual HydroConservedState sourceForCell(
      std::size_t cell_index,
      const HydroConservedState& conserved,
      const HydroPrimitiveState& primitive,
      const HydroSourceContext& context) const = 0;
};

class ComovingGravityExpansionSource final : public HydroSourceTerm {
 public:
  [[nodiscard]] HydroConservedState sourceForCell(
      std::size_t cell_index,
      const HydroConservedState& conserved,
      const HydroPrimitiveState& primitive,
      const HydroSourceContext& context) const override;
};

struct HydroProfileEvent {
  std::uint64_t bytes_moved = 0;
  std::uint64_t face_count = 0;
  std::uint64_t limiter_clip_count = 0;
  std::uint64_t positivity_fallback_count = 0;
  std::uint64_t riemann_fallback_count = 0;
  double reconstruct_ms = 0.0;
  double riemann_ms = 0.0;
  double accumulate_ms = 0.0;
  double source_ms = 0.0;
  double total_ms = 0.0;
  HydroConservationReport conservation;
};

struct HydroScratchBuffers {
  // Transient per-step workspaces. None of these arrays are restart-authoritative.
  std::vector<HydroConservedState> cell_delta;
  std::vector<HydroPrimitiveState> left_states;
  std::vector<HydroPrimitiveState> right_states;
  std::vector<HydroConservedState> fluxes;
  std::vector<std::size_t> touched_cells;
  std::vector<std::size_t> full_active_cells;
  std::vector<std::size_t> full_active_faces;

  void resize(std::size_t cell_count, std::size_t active_face_count);
};

class HydroCoreSolver {
 public:
  explicit HydroCoreSolver(double adiabatic_index);

  [[nodiscard]] double adiabaticIndex() const;

  [[nodiscard]] static HydroConservedState conservedFromPrimitive(
      const HydroPrimitiveState& primitive,
      double adiabatic_index);

  [[nodiscard]] static HydroPrimitiveState primitiveFromConserved(
      const HydroConservedState& conserved,
      double adiabatic_index);

  [[nodiscard]] static HydroConservationTotals conservationTotals(
      const HydroConservedStateSoa& conserved,
      const HydroPatchGeometry& geometry,
      std::span<const std::size_t> cells);

  void advancePatch(
      HydroConservedStateSoa& conserved,
      const HydroPatchGeometry& geometry,
      const HydroUpdateContext& update,
      const HydroReconstruction& reconstruction,
      const HydroRiemannSolver& riemann_solver,
      std::span<const HydroSourceTerm* const> source_terms,
      const HydroSourceContext& source_context,
      HydroProfileEvent* profile = nullptr,
      HydroFluxRegisterSink* flux_register_sink = nullptr) const;

  void advancePatchWithScratch(
      HydroConservedStateSoa& conserved,
      const HydroPatchGeometry& geometry,
      const HydroUpdateContext& update,
      const HydroReconstruction& reconstruction,
      const HydroRiemannSolver& riemann_solver,
      std::span<const HydroSourceTerm* const> source_terms,
      const HydroSourceContext& source_context,
      HydroScratchBuffers& scratch,
      HydroPrimitiveCacheSoa* primitive_cache,
      HydroProfileEvent* profile = nullptr,
      HydroFluxRegisterSink* flux_register_sink = nullptr) const;

  void advancePatchActiveSet(
      HydroConservedStateSoa& conserved,
      const HydroPatchGeometry& geometry,
      const HydroActiveSetView& active_set,
      const HydroUpdateContext& update,
      const HydroReconstruction& reconstruction,
      const HydroRiemannSolver& riemann_solver,
      std::span<const HydroSourceTerm* const> source_terms,
      const HydroSourceContext& source_context,
      HydroProfileEvent* profile = nullptr,
      HydroFluxRegisterSink* flux_register_sink = nullptr) const;

  void advancePatchActiveSetWithScratch(
      HydroConservedStateSoa& conserved,
      const HydroPatchGeometry& geometry,
      const HydroActiveSetView& active_set,
      const HydroUpdateContext& update,
      const HydroReconstruction& reconstruction,
      const HydroRiemannSolver& riemann_solver,
      std::span<const HydroSourceTerm* const> source_terms,
      const HydroSourceContext& source_context,
      HydroScratchBuffers& scratch,
      HydroPrimitiveCacheSoa* primitive_cache,
      HydroProfileEvent* profile = nullptr,
      HydroFluxRegisterSink* flux_register_sink = nullptr) const;

 private:
  double m_adiabatic_index = 5.0 / 3.0;
};

}  // namespace cosmosim::hydro
