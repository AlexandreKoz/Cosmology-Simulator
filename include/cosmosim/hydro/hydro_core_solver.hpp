#pragma once

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include "cosmosim/hydro/hydro_reconstruction.hpp"
#include "cosmosim/hydro/hydro_riemann.hpp"

namespace cosmosim::hydro {

constexpr std::size_t k_invalid_cell_index = static_cast<std::size_t>(-1);

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

struct HydroFace {
  std::size_t owner_cell = k_invalid_cell_index;
  std::size_t neighbor_cell = k_invalid_cell_index;
  double area_comoving = 0.0;
  double normal_x = 0.0;
  double normal_y = 0.0;
  double normal_z = 0.0;
};

struct HydroPatchGeometry {
  double cell_volume_comoving = 0.0;
  std::vector<HydroFace> faces;
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
};

struct HydroScratchBuffers {
  std::vector<HydroConservedState> cell_delta;
  std::vector<HydroPrimitiveState> left_states;
  std::vector<HydroPrimitiveState> right_states;
  std::vector<HydroConservedState> fluxes;

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

  void advancePatch(
      HydroConservedStateSoa& conserved,
      const HydroPatchGeometry& geometry,
      const HydroUpdateContext& update,
      const HydroReconstruction& reconstruction,
      const HydroRiemannSolver& riemann_solver,
      std::span<const HydroSourceTerm* const> source_terms,
      const HydroSourceContext& source_context,
      HydroProfileEvent* profile = nullptr) const;

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
      HydroProfileEvent* profile = nullptr) const;

  void advancePatchActiveSet(
      HydroConservedStateSoa& conserved,
      const HydroPatchGeometry& geometry,
      const HydroActiveSetView& active_set,
      const HydroUpdateContext& update,
      const HydroReconstruction& reconstruction,
      const HydroRiemannSolver& riemann_solver,
      std::span<const HydroSourceTerm* const> source_terms,
      const HydroSourceContext& source_context,
      HydroProfileEvent* profile = nullptr) const;

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
      HydroProfileEvent* profile = nullptr) const;

 private:
  double m_adiabatic_index = 5.0 / 3.0;
};

}  // namespace cosmosim::hydro
