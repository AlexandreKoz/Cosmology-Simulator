#include "cosmosim/hydro/hydro_core_solver.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace cosmosim::hydro {
namespace {

constexpr double k_small = 1.0e-14;

[[nodiscard]] double dot3(double ax, double ay, double az, double bx, double by, double bz) {
  return ax * bx + ay * by + az * bz;
}

[[nodiscard]] double norm3(double x, double y, double z) {
  return std::sqrt(dot3(x, y, z, x, y, z));
}

[[nodiscard]] std::array<double, 3> gravityAtCell(
    std::size_t cell_index,
    std::span<const double> ax,
    std::span<const double> ay,
    std::span<const double> az) {
  const double gx = (cell_index < ax.size()) ? ax[cell_index] : 0.0;
  const double gy = (cell_index < ay.size()) ? ay[cell_index] : 0.0;
  const double gz = (cell_index < az.size()) ? az[cell_index] : 0.0;
  return {gx, gy, gz};
}

void validateAdvanceInputs(
    const HydroConservedStateSoa& conserved,
    const HydroPatchGeometry& geometry,
    const HydroUpdateContext& update,
    const HydroSourceContext& source_context,
    double adiabatic_index) {
  if (conserved.size() == 0) {
    throw std::invalid_argument("Hydro advance requires at least one cell");
  }
  if (geometry.cell_volume_comoving <= 0.0) {
    throw std::invalid_argument("Hydro advance requires positive cell_volume_comoving");
  }
  if (update.dt_code <= 0.0) {
    throw std::invalid_argument("Hydro advance requires dt_code > 0");
  }
  if (update.scale_factor <= 0.0) {
    throw std::invalid_argument("Hydro advance requires scale_factor > 0");
  }
  if (adiabatic_index <= 1.0) {
    throw std::invalid_argument("Hydro solver requires adiabatic_index > 1");
  }
  if (std::abs(source_context.update.dt_code - update.dt_code) > 1.0e-16 ||
      std::abs(source_context.update.scale_factor - update.scale_factor) > 1.0e-16 ||
      std::abs(source_context.update.hubble_rate_code - update.hubble_rate_code) > 1.0e-16) {
    throw std::invalid_argument("Hydro source_context.update must match update passed to advancePatch");
  }

  const std::size_t real_cell_count = geometry.cellCount() == 0 ? conserved.size() : geometry.cellCount();
  for (const HydroFace& face : geometry.faces) {
    if (face.owner_cell >= conserved.size()) {
      throw std::invalid_argument("Hydro face owner index is out of range");
    }
    if (face.owner_cell >= real_cell_count) {
      throw std::invalid_argument("Hydro face owner must be an authoritative real cell");
    }
    if (face.neighbor_cell != k_invalid_cell_index && face.neighbor_cell >= conserved.size()) {
      throw std::invalid_argument("Hydro face neighbor index is out of range");
    }
    if (face.owner_minus_cell != k_invalid_cell_index && face.owner_minus_cell >= conserved.size()) {
      throw std::invalid_argument("Hydro face owner-minus stencil index is out of range");
    }
    if (face.neighbor_plus_cell != k_invalid_cell_index && face.neighbor_plus_cell >= conserved.size()) {
      throw std::invalid_argument("Hydro face neighbor-plus stencil index is out of range");
    }
    if (face.ghost_cell_slot != k_invalid_ghost_cell_slot) {
      if (face.ghost_cell_slot >= geometry.ghost_cells.size()) {
        throw std::invalid_argument("Hydro face references an invalid ghost cell slot");
      }
      const HydroGhostCell& ghost = geometry.ghost_cells[face.ghost_cell_slot];
      if (face.neighbor_cell != ghost.ghost_cell || face.owner_cell != ghost.owner_real_cell) {
        throw std::invalid_argument("Hydro face ghost metadata does not match face endpoints");
      }
    }
    if (face.area_comoving <= 0.0) {
      throw std::invalid_argument("Hydro face area_comoving must be positive");
    }
    const double normal_norm = norm3(face.normal_x, face.normal_y, face.normal_z);
    if (std::abs(normal_norm - 1.0) > 1.0e-10) {
      throw std::invalid_argument("Hydro face normal must be unit length");
    }
  }
  if (!geometry.flux_register_faces.empty() &&
      geometry.flux_register_faces.size() != geometry.faces.size()) {
    throw std::invalid_argument("Hydro flux-register metadata must be empty or match face count");
  }
  for (const HydroFluxRegisterFace& register_face : geometry.flux_register_faces) {
    if (register_face.role == HydroFluxRegisterFaceRole::kNone) {
      continue;
    }
    if (register_face.register_key == 0U) {
      throw std::invalid_argument("Hydro flux-register face requires nonzero register_key");
    }
    if (register_face.coarse_patch_id == 0U) {
      throw std::invalid_argument("Hydro flux-register face requires nonzero coarse_patch_id");
    }
    if (register_face.coarse_cell_index == k_invalid_cell_index) {
      throw std::invalid_argument("Hydro flux-register face requires valid coarse_cell_index");
    }
    if (register_face.coarse_orientation_sign != 1.0 &&
        register_face.coarse_orientation_sign != -1.0) {
      throw std::invalid_argument("Hydro flux-register coarse_orientation_sign must be +1 or -1");
    }
  }
}

void validateActiveSet(const HydroActiveSetView& active_set, const HydroConservedStateSoa& conserved, const HydroPatchGeometry& geometry) {
  if (active_set.active_cells.empty()) {
    throw std::invalid_argument("Hydro active set requires at least one active cell");
  }

  const std::size_t real_cell_count = geometry.cellCount() == 0 ? conserved.size() : geometry.cellCount();
  for (std::size_t cell_index : active_set.active_cells) {
    if (cell_index >= conserved.size()) {
      throw std::invalid_argument("Hydro active cell index is out of range");
    }
    if (cell_index >= real_cell_count) {
      throw std::invalid_argument("Hydro active cells must be authoritative real cells");
    }
  }
  for (std::size_t face_index : active_set.active_faces) {
    if (face_index >= geometry.faces.size()) {
      throw std::invalid_argument("Hydro active face index is out of range");
    }
  }
}

void fillPrimitiveCache(
    const HydroConservedStateSoa& conserved,
    const HydroActiveSetView& active_set,
    const HydroPatchGeometry& geometry,
    double adiabatic_index,
    HydroPrimitiveCacheSoa& primitive_cache) {
  if (primitive_cache.size() != conserved.size()) {
    primitive_cache.resize(conserved.size());
  }

  for (std::size_t cell_index : active_set.active_cells) {
    primitive_cache.storeCell(
        cell_index,
        HydroCoreSolver::primitiveFromConserved(conserved.loadCell(cell_index), adiabatic_index));
  }
  for (std::size_t face_index : active_set.active_faces) {
    const HydroFace& face = geometry.faces[face_index];
    primitive_cache.storeCell(
        face.owner_cell,
        HydroCoreSolver::primitiveFromConserved(conserved.loadCell(face.owner_cell), adiabatic_index));
    if (face.neighbor_cell != k_invalid_cell_index) {
      primitive_cache.storeCell(
          face.neighbor_cell,
          HydroCoreSolver::primitiveFromConserved(conserved.loadCell(face.neighbor_cell), adiabatic_index));
    }
    if (face.owner_minus_cell != k_invalid_cell_index) {
      primitive_cache.storeCell(
          face.owner_minus_cell,
          HydroCoreSolver::primitiveFromConserved(conserved.loadCell(face.owner_minus_cell), adiabatic_index));
    }
    if (face.neighbor_plus_cell != k_invalid_cell_index) {
      primitive_cache.storeCell(
          face.neighbor_plus_cell,
          HydroCoreSolver::primitiveFromConserved(conserved.loadCell(face.neighbor_plus_cell), adiabatic_index));
    }
  }
}

[[nodiscard]] HydroConservationTotals totalsFromDelta(
    const HydroConservedState& delta,
    const HydroConservedState& reference_state,
    double cell_volume_comoving) {
  const double rho = std::max(reference_state.mass_density_comoving, k_small);
  const double old_kinetic_density = 0.5 *
      (reference_state.momentum_density_x_comoving * reference_state.momentum_density_x_comoving +
       reference_state.momentum_density_y_comoving * reference_state.momentum_density_y_comoving +
       reference_state.momentum_density_z_comoving * reference_state.momentum_density_z_comoving) /
      rho;
  const HydroConservedState updated = reference_state + delta;
  const double updated_rho = std::max(updated.mass_density_comoving, k_small);
  const double new_kinetic_density = 0.5 *
      (updated.momentum_density_x_comoving * updated.momentum_density_x_comoving +
       updated.momentum_density_y_comoving * updated.momentum_density_y_comoving +
       updated.momentum_density_z_comoving * updated.momentum_density_z_comoving) /
      updated_rho;

  HydroConservationTotals totals;
  totals.mass = delta.mass_density_comoving * cell_volume_comoving;
  totals.momentum_x = delta.momentum_density_x_comoving * cell_volume_comoving;
  totals.momentum_y = delta.momentum_density_y_comoving * cell_volume_comoving;
  totals.momentum_z = delta.momentum_density_z_comoving * cell_volume_comoving;
  totals.total_energy = delta.total_energy_density_comoving * cell_volume_comoving;
  totals.internal_energy =
      ((updated.total_energy_density_comoving - new_kinetic_density) -
       (reference_state.total_energy_density_comoving - old_kinetic_density)) *
      cell_volume_comoving;
  return totals;
}

[[nodiscard]] std::size_t fluxNeighborTargetCell(
    const HydroPatchGeometry& geometry,
    const HydroFace& face) {
  if (face.ghost_cell_slot == k_invalid_ghost_cell_slot) {
    return face.neighbor_cell;
  }
  if (face.ghost_cell_slot >= geometry.ghost_cells.size()) {
    throw std::invalid_argument("Hydro face references an invalid ghost cell slot");
  }

  const HydroGhostCell& ghost = geometry.ghost_cells[face.ghost_cell_slot];
  if (ghost.boundary_kind == HydroBoundaryKind::kPeriodic) {
    // Periodic ghosts are only reconstruction scratch.  The conservative
    // update must land on the wrapped authoritative real cell, otherwise a
    // production periodic patch leaks flux into non-authoritative ghost storage.
    return ghost.source_real_cell;
  }

  // Open and reflective physical-boundary ghosts are not authoritative cells.
  // Imported MPI ghosts are scratch rows that intentionally accumulate the
  // conservative opposite-side face delta; the workflow restores the row and
  // routes that delta back to the authoritative owner by stable gas_cell_id.
  if (ghost.boundary_kind == HydroBoundaryKind::kImportedMpi) {
    return ghost.ghost_cell;
  }
  return k_invalid_cell_index;
}

void maybeRecordFluxRegister(
    const HydroPatchGeometry& geometry,
    std::size_t face_index,
    const HydroUpdateContext& update,
    const HydroConservedState& flux_code,
    HydroFluxRegisterSink* flux_register_sink) {
  if (flux_register_sink == nullptr || geometry.flux_register_faces.empty()) {
    return;
  }
  const HydroFluxRegisterFace& register_face = geometry.flux_register_faces[face_index];
  if (register_face.role == HydroFluxRegisterFaceRole::kNone) {
    return;
  }
  const HydroConservedState oriented_flux = register_face.coarse_orientation_sign * flux_code;
  flux_register_sink->recordFaceFlux(HydroFluxRegisterRecord{
      .role = register_face.role,
      .register_key = register_face.register_key,
      .coarse_patch_id = register_face.coarse_patch_id,
      .coarse_gas_cell_id = register_face.coarse_gas_cell_id,
      .coarse_cell_index = register_face.coarse_cell_index,
      .level = register_face.level,
      .axis = register_face.axis,
      .orientation = register_face.orientation,
      .face_area_comoving = geometry.faces[face_index].area_comoving,
      .dt_code = update.dt_code,
      .flux_code = oriented_flux});
}

}  // namespace

HydroConservedState& HydroConservedState::operator+=(const HydroConservedState& rhs) {
  mass_density_comoving += rhs.mass_density_comoving;
  momentum_density_x_comoving += rhs.momentum_density_x_comoving;
  momentum_density_y_comoving += rhs.momentum_density_y_comoving;
  momentum_density_z_comoving += rhs.momentum_density_z_comoving;
  total_energy_density_comoving += rhs.total_energy_density_comoving;
  return *this;
}

HydroConservedState& HydroConservedState::operator-=(const HydroConservedState& rhs) {
  mass_density_comoving -= rhs.mass_density_comoving;
  momentum_density_x_comoving -= rhs.momentum_density_x_comoving;
  momentum_density_y_comoving -= rhs.momentum_density_y_comoving;
  momentum_density_z_comoving -= rhs.momentum_density_z_comoving;
  total_energy_density_comoving -= rhs.total_energy_density_comoving;
  return *this;
}

HydroConservedState operator+(HydroConservedState lhs, const HydroConservedState& rhs) {
  lhs += rhs;
  return lhs;
}

HydroConservedState operator-(HydroConservedState lhs, const HydroConservedState& rhs) {
  lhs -= rhs;
  return lhs;
}

HydroConservedState operator*(double scalar, HydroConservedState state) {
  state.mass_density_comoving *= scalar;
  state.momentum_density_x_comoving *= scalar;
  state.momentum_density_y_comoving *= scalar;
  state.momentum_density_z_comoving *= scalar;
  state.total_energy_density_comoving *= scalar;
  return state;
}

HydroConservationTotals& HydroConservationTotals::operator+=(const HydroConservationTotals& rhs) {
  mass += rhs.mass;
  momentum_x += rhs.momentum_x;
  momentum_y += rhs.momentum_y;
  momentum_z += rhs.momentum_z;
  total_energy += rhs.total_energy;
  internal_energy += rhs.internal_energy;
  return *this;
}

HydroConservationTotals& HydroConservationTotals::operator-=(const HydroConservationTotals& rhs) {
  mass -= rhs.mass;
  momentum_x -= rhs.momentum_x;
  momentum_y -= rhs.momentum_y;
  momentum_z -= rhs.momentum_z;
  total_energy -= rhs.total_energy;
  internal_energy -= rhs.internal_energy;
  return *this;
}

HydroConservationTotals operator+(HydroConservationTotals lhs, const HydroConservationTotals& rhs) {
  lhs += rhs;
  return lhs;
}

HydroConservationTotals operator-(HydroConservationTotals lhs, const HydroConservationTotals& rhs) {
  lhs -= rhs;
  return lhs;
}

HydroConservedStateSoa::HydroConservedStateSoa(std::size_t cell_count)
    : m_mass_density_comoving(cell_count, 0.0),
      m_momentum_density_x_comoving(cell_count, 0.0),
      m_momentum_density_y_comoving(cell_count, 0.0),
      m_momentum_density_z_comoving(cell_count, 0.0),
      m_total_energy_density_comoving(cell_count, 0.0) {}

void HydroConservedStateSoa::resize(std::size_t cell_count) {
  m_mass_density_comoving.resize(cell_count, 0.0);
  m_momentum_density_x_comoving.resize(cell_count, 0.0);
  m_momentum_density_y_comoving.resize(cell_count, 0.0);
  m_momentum_density_z_comoving.resize(cell_count, 0.0);
  m_total_energy_density_comoving.resize(cell_count, 0.0);
}

std::size_t HydroConservedStateSoa::size() const { return m_mass_density_comoving.size(); }

HydroConservedState HydroConservedStateSoa::loadCell(std::size_t cell_index) const {
  return HydroConservedState{
      .mass_density_comoving = m_mass_density_comoving.at(cell_index),
      .momentum_density_x_comoving = m_momentum_density_x_comoving.at(cell_index),
      .momentum_density_y_comoving = m_momentum_density_y_comoving.at(cell_index),
      .momentum_density_z_comoving = m_momentum_density_z_comoving.at(cell_index),
      .total_energy_density_comoving = m_total_energy_density_comoving.at(cell_index)};
}

void HydroConservedStateSoa::storeCell(std::size_t cell_index, const HydroConservedState& cell_state) {
  m_mass_density_comoving.at(cell_index) = cell_state.mass_density_comoving;
  m_momentum_density_x_comoving.at(cell_index) = cell_state.momentum_density_x_comoving;
  m_momentum_density_y_comoving.at(cell_index) = cell_state.momentum_density_y_comoving;
  m_momentum_density_z_comoving.at(cell_index) = cell_state.momentum_density_z_comoving;
  m_total_energy_density_comoving.at(cell_index) = cell_state.total_energy_density_comoving;
}

std::span<double> HydroConservedStateSoa::massDensityComoving() { return m_mass_density_comoving; }
std::span<const double> HydroConservedStateSoa::massDensityComoving() const { return m_mass_density_comoving; }
std::span<double> HydroConservedStateSoa::momentumDensityXComoving() { return m_momentum_density_x_comoving; }
std::span<const double> HydroConservedStateSoa::momentumDensityXComoving() const { return m_momentum_density_x_comoving; }
std::span<double> HydroConservedStateSoa::momentumDensityYComoving() { return m_momentum_density_y_comoving; }
std::span<const double> HydroConservedStateSoa::momentumDensityYComoving() const { return m_momentum_density_y_comoving; }
std::span<double> HydroConservedStateSoa::momentumDensityZComoving() { return m_momentum_density_z_comoving; }
std::span<const double> HydroConservedStateSoa::momentumDensityZComoving() const { return m_momentum_density_z_comoving; }
std::span<double> HydroConservedStateSoa::totalEnergyDensityComoving() { return m_total_energy_density_comoving; }
std::span<const double> HydroConservedStateSoa::totalEnergyDensityComoving() const { return m_total_energy_density_comoving; }

HydroPrimitiveCacheSoa::HydroPrimitiveCacheSoa(std::size_t cell_count)
    : m_rho_comoving(cell_count, 0.0),
      m_vel_x_peculiar(cell_count, 0.0),
      m_vel_y_peculiar(cell_count, 0.0),
      m_vel_z_peculiar(cell_count, 0.0),
      m_pressure_comoving(cell_count, 0.0) {}

void HydroPrimitiveCacheSoa::resize(std::size_t cell_count) {
  m_rho_comoving.resize(cell_count, 0.0);
  m_vel_x_peculiar.resize(cell_count, 0.0);
  m_vel_y_peculiar.resize(cell_count, 0.0);
  m_vel_z_peculiar.resize(cell_count, 0.0);
  m_pressure_comoving.resize(cell_count, 0.0);
}

std::size_t HydroPrimitiveCacheSoa::size() const { return m_rho_comoving.size(); }

HydroPrimitiveState HydroPrimitiveCacheSoa::loadCell(std::size_t cell_index) const {
  return HydroPrimitiveState{
      .rho_comoving = m_rho_comoving.at(cell_index),
      .vel_x_peculiar = m_vel_x_peculiar.at(cell_index),
      .vel_y_peculiar = m_vel_y_peculiar.at(cell_index),
      .vel_z_peculiar = m_vel_z_peculiar.at(cell_index),
      .pressure_comoving = m_pressure_comoving.at(cell_index)};
}

void HydroPrimitiveCacheSoa::storeCell(std::size_t cell_index, const HydroPrimitiveState& primitive_state) {
  m_rho_comoving.at(cell_index) = primitive_state.rho_comoving;
  m_vel_x_peculiar.at(cell_index) = primitive_state.vel_x_peculiar;
  m_vel_y_peculiar.at(cell_index) = primitive_state.vel_y_peculiar;
  m_vel_z_peculiar.at(cell_index) = primitive_state.vel_z_peculiar;
  m_pressure_comoving.at(cell_index) = primitive_state.pressure_comoving;
}

HydroConservedState ComovingGravityExpansionSource::sourceForCell(
    std::size_t cell_index,
    const HydroConservedState& conserved,
    const HydroPrimitiveState& primitive,
    const HydroSourceContext& context) const {
  const std::array<double, 3> g = gravityAtCell(
      cell_index,
      context.gravity_accel_x_peculiar,
      context.gravity_accel_y_peculiar,
      context.gravity_accel_z_peculiar);

  if (!std::isfinite(context.update.scale_factor) || context.update.scale_factor <= 0.0) {
    throw std::invalid_argument("comoving gravity source requires finite scale_factor > 0");
  }
  if (!std::isfinite(context.update.hubble_rate_code) || context.update.hubble_rate_code < 0.0) {
    throw std::invalid_argument("comoving gravity source requires finite hubble_rate_code >= 0");
  }
  // TreePM returns the scale-free comoving particle kernel A. The conserved
  // gas momentum stores rho_c u with physical peculiar velocity u, whose
  // equation is du/dt + H u = A/a^2. Keep the scale-factor response here so
  // particles and gas share one TreePM normalization without double counting.
  const double inverse_scale_factor_squared =
      1.0 / (context.update.scale_factor * context.update.scale_factor);
  const std::array<double, 3> peculiar_acceleration{
      g[0] * inverse_scale_factor_squared,
      g[1] * inverse_scale_factor_squared,
      g[2] * inverse_scale_factor_squared,
  };
  const double hubble_rate = context.update.hubble_rate_code;
  const double internal_energy_density = std::max(
      conserved.total_energy_density_comoving -
          0.5 * (conserved.momentum_density_x_comoving * conserved.momentum_density_x_comoving +
                 conserved.momentum_density_y_comoving * conserved.momentum_density_y_comoving +
                 conserved.momentum_density_z_comoving * conserved.momentum_density_z_comoving) /
              std::max(conserved.mass_density_comoving, k_small),
      0.0);
  const double kinetic_energy_density = std::max(
      conserved.total_energy_density_comoving - internal_energy_density,
      0.0);

  HydroConservedState source;
  source.mass_density_comoving = 0.0;
  source.momentum_density_x_comoving =
      primitive.rho_comoving * peculiar_acceleration[0] -
      hubble_rate * conserved.momentum_density_x_comoving;
  source.momentum_density_y_comoving =
      primitive.rho_comoving * peculiar_acceleration[1] -
      hubble_rate * conserved.momentum_density_y_comoving;
  source.momentum_density_z_comoving =
      primitive.rho_comoving * peculiar_acceleration[2] -
      hubble_rate * conserved.momentum_density_z_comoving;

  const double work_gravity = primitive.rho_comoving *
      (primitive.vel_x_peculiar * peculiar_acceleration[0] +
       primitive.vel_y_peculiar * peculiar_acceleration[1] +
       primitive.vel_z_peculiar * peculiar_acceleration[2]);
  // For rho_c, rho_c u, and rho_c(e + u^2/2), homogeneous expansion gives
  // -H rho_c u^2 - 3 H P. This is -2 H E for gamma=5/3, while the pressure
  // form remains correct for every gamma supported by HydroCoreSolver.
  source.total_energy_density_comoving =
      work_gravity - hubble_rate *
          (2.0 * kinetic_energy_density + 3.0 * primitive.pressure_comoving);
  return source;
}

HydroCoreSolver::HydroCoreSolver(double adiabatic_index) : m_adiabatic_index(adiabatic_index) {
  if (m_adiabatic_index <= 1.0) {
    throw std::invalid_argument("HydroCoreSolver requires adiabatic_index > 1");
  }
}

double HydroCoreSolver::adiabaticIndex() const { return m_adiabatic_index; }

void HydroScratchBuffers::resize(std::size_t cell_count, std::size_t active_face_count) {
  cell_delta.resize(cell_count);
  left_states.resize(active_face_count);
  right_states.resize(active_face_count);
  fluxes.resize(active_face_count);
}

HydroConservedState HydroCoreSolver::conservedFromPrimitive(
    const HydroPrimitiveState& primitive,
    double adiabatic_index) {
  if (primitive.rho_comoving <= 0.0) {
    throw std::invalid_argument("Primitive state requires rho_comoving > 0");
  }
  if (primitive.pressure_comoving <= 0.0) {
    throw std::invalid_argument("Primitive state requires pressure_comoving > 0");
  }

  const double velocity_squared =
      primitive.vel_x_peculiar * primitive.vel_x_peculiar +
      primitive.vel_y_peculiar * primitive.vel_y_peculiar +
      primitive.vel_z_peculiar * primitive.vel_z_peculiar;

  HydroConservedState conserved;
  conserved.mass_density_comoving = primitive.rho_comoving;
  conserved.momentum_density_x_comoving = primitive.rho_comoving * primitive.vel_x_peculiar;
  conserved.momentum_density_y_comoving = primitive.rho_comoving * primitive.vel_y_peculiar;
  conserved.momentum_density_z_comoving = primitive.rho_comoving * primitive.vel_z_peculiar;
  conserved.total_energy_density_comoving =
      primitive.pressure_comoving / (adiabatic_index - 1.0) + 0.5 * primitive.rho_comoving * velocity_squared;
  return conserved;
}

HydroPrimitiveState HydroCoreSolver::primitiveFromConserved(
    const HydroConservedState& conserved,
    double adiabatic_index) {
  if (conserved.mass_density_comoving <= 0.0) {
    throw std::invalid_argument("Conserved state requires mass_density_comoving > 0");
  }

  const double inv_rho = 1.0 / conserved.mass_density_comoving;
  HydroPrimitiveState primitive;
  primitive.rho_comoving = conserved.mass_density_comoving;
  primitive.vel_x_peculiar = conserved.momentum_density_x_comoving * inv_rho;
  primitive.vel_y_peculiar = conserved.momentum_density_y_comoving * inv_rho;
  primitive.vel_z_peculiar = conserved.momentum_density_z_comoving * inv_rho;

  const double kinetic_density = 0.5 *
      (conserved.momentum_density_x_comoving * conserved.momentum_density_x_comoving +
       conserved.momentum_density_y_comoving * conserved.momentum_density_y_comoving +
       conserved.momentum_density_z_comoving * conserved.momentum_density_z_comoving) *
      inv_rho;
  const double internal_density = conserved.total_energy_density_comoving - kinetic_density;
  primitive.pressure_comoving = std::max((adiabatic_index - 1.0) * internal_density, k_small);
  return primitive;
}

HydroConservationTotals HydroCoreSolver::conservationTotals(
    const HydroConservedStateSoa& conserved,
    const HydroPatchGeometry& geometry,
    std::span<const std::size_t> cells) {
  if (geometry.cell_volume_comoving <= 0.0) {
    throw std::invalid_argument("Hydro conservation totals require positive cell_volume_comoving");
  }

  HydroConservationTotals totals;
  const std::size_t real_cell_count = geometry.cellCount() == 0 ? conserved.size() : geometry.cellCount();
  for (std::size_t cell_index : cells) {
    if (cell_index >= conserved.size()) {
      throw std::out_of_range("Hydro conservation totals cell index is out of range");
    }
    if (cell_index >= real_cell_count) {
      throw std::invalid_argument("Hydro conservation totals require authoritative real cells");
    }
    const HydroConservedState state = conserved.loadCell(cell_index);
    const double rho = std::max(state.mass_density_comoving, k_small);
    const double kinetic_density = 0.5 *
        (state.momentum_density_x_comoving * state.momentum_density_x_comoving +
         state.momentum_density_y_comoving * state.momentum_density_y_comoving +
         state.momentum_density_z_comoving * state.momentum_density_z_comoving) /
        rho;
    totals.mass += state.mass_density_comoving * geometry.cell_volume_comoving;
    totals.momentum_x += state.momentum_density_x_comoving * geometry.cell_volume_comoving;
    totals.momentum_y += state.momentum_density_y_comoving * geometry.cell_volume_comoving;
    totals.momentum_z += state.momentum_density_z_comoving * geometry.cell_volume_comoving;
    totals.total_energy += state.total_energy_density_comoving * geometry.cell_volume_comoving;
    totals.internal_energy +=
        (state.total_energy_density_comoving - kinetic_density) * geometry.cell_volume_comoving;
  }
  return totals;
}

void HydroCoreSolver::advancePatch(
    HydroConservedStateSoa& conserved,
    const HydroPatchGeometry& geometry,
    const HydroUpdateContext& update,
    const HydroReconstruction& reconstruction,
    const HydroRiemannSolver& riemann_solver,
    std::span<const HydroSourceTerm* const> source_terms,
    const HydroSourceContext& source_context,
    HydroProfileEvent* profile,
    HydroFluxRegisterSink* flux_register_sink) const {
  HydroScratchBuffers scratch;
  HydroPrimitiveCacheSoa primitive_cache(conserved.size());
  advancePatchWithScratch(
      conserved,
      geometry,
      update,
      reconstruction,
      riemann_solver,
      source_terms,
      source_context,
      scratch,
      &primitive_cache,
      profile,
      flux_register_sink);
}

void HydroCoreSolver::advancePatchWithScratch(
    HydroConservedStateSoa& conserved,
    const HydroPatchGeometry& geometry,
    const HydroUpdateContext& update,
    const HydroReconstruction& reconstruction,
    const HydroRiemannSolver& riemann_solver,
    std::span<const HydroSourceTerm* const> source_terms,
    const HydroSourceContext& source_context,
    HydroScratchBuffers& scratch,
    HydroPrimitiveCacheSoa* primitive_cache,
    HydroProfileEvent* profile,
    HydroFluxRegisterSink* flux_register_sink) const {
  scratch.full_active_cells.resize(conserved.size());
  for (std::size_t i = 0; i < scratch.full_active_cells.size(); ++i) {
    scratch.full_active_cells[i] = i;
  }
  scratch.full_active_faces.resize(geometry.faces.size());
  for (std::size_t i = 0; i < scratch.full_active_faces.size(); ++i) {
    scratch.full_active_faces[i] = i;
  }
  advancePatchActiveSetWithScratch(
      conserved,
      geometry,
      HydroActiveSetView{.active_cells = scratch.full_active_cells, .active_faces = scratch.full_active_faces},
      update,
      reconstruction,
      riemann_solver,
      source_terms,
      source_context,
      scratch,
      primitive_cache,
      profile,
      flux_register_sink);
}

void HydroCoreSolver::advancePatchActiveSet(
    HydroConservedStateSoa& conserved,
    const HydroPatchGeometry& geometry,
    const HydroActiveSetView& active_set,
    const HydroUpdateContext& update,
    const HydroReconstruction& reconstruction,
    const HydroRiemannSolver& riemann_solver,
    std::span<const HydroSourceTerm* const> source_terms,
    const HydroSourceContext& source_context,
    HydroProfileEvent* profile,
    HydroFluxRegisterSink* flux_register_sink) const {
  HydroScratchBuffers scratch;
  HydroPrimitiveCacheSoa primitive_cache(conserved.size());
  advancePatchActiveSetWithScratch(
      conserved,
      geometry,
      active_set,
      update,
      reconstruction,
      riemann_solver,
      source_terms,
      source_context,
      scratch,
      &primitive_cache,
      profile,
      flux_register_sink);
}

void HydroCoreSolver::advancePatchActiveSetWithScratch(
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
    HydroProfileEvent* profile,
    HydroFluxRegisterSink* flux_register_sink) const {
  const std::uint64_t limiter_before = reconstruction.limiterClipCount();
  const std::uint64_t positivity_before = reconstruction.positivityFallbackCount();
  const std::uint64_t riemann_fallback_before = riemann_solver.fallbackCount();

  validateAdvanceInputs(conserved, geometry, update, source_context, m_adiabatic_index);
  validateActiveSet(active_set, conserved, geometry);
  scratch.resize(conserved.size(), active_set.active_faces.size());
  HydroConservationReport conservation_report;

  std::vector<unsigned char> active_cell_mask(conserved.size(), 0U);
  scratch.touched_cells.clear();
  scratch.touched_cells.reserve(active_set.active_faces.size() * 2U + active_set.active_cells.size());
  for (std::size_t cell_index : active_set.active_cells) {
    active_cell_mask[cell_index] = 1U;
    scratch.cell_delta[cell_index] = HydroConservedState{};
    scratch.touched_cells.push_back(cell_index);
  }
  for (std::size_t face_index : active_set.active_faces) {
    const HydroFace& face = geometry.faces[face_index];
    scratch.cell_delta[face.owner_cell] = HydroConservedState{};
    scratch.touched_cells.push_back(face.owner_cell);
    const std::size_t neighbor_target = fluxNeighborTargetCell(geometry, face);
    if (neighbor_target != k_invalid_cell_index) {
      scratch.cell_delta[neighbor_target] = HydroConservedState{};
      scratch.touched_cells.push_back(neighbor_target);
    }
  }
  std::sort(scratch.touched_cells.begin(), scratch.touched_cells.end());
  scratch.touched_cells.erase(
      std::unique(scratch.touched_cells.begin(), scratch.touched_cells.end()),
      scratch.touched_cells.end());
  const std::size_t real_cell_count = geometry.cellCount() == 0 ? conserved.size() : geometry.cellCount();
  scratch.full_active_cells.clear();
  scratch.full_active_cells.reserve(scratch.touched_cells.size());
  for (const std::size_t cell_index : scratch.touched_cells) {
    if (cell_index < real_cell_count) {
      scratch.full_active_cells.push_back(cell_index);
    }
  }
  conservation_report.before = conservationTotals(conserved, geometry, scratch.full_active_cells);
  conservation_report.cell_count = static_cast<std::uint64_t>(scratch.full_active_cells.size());
  if (primitive_cache != nullptr) {
    fillPrimitiveCache(conserved, active_set, geometry, m_adiabatic_index, *primitive_cache);
  }

  const auto total_start = std::chrono::steady_clock::now();
  const auto reconstruct_start = std::chrono::steady_clock::now();
  for (std::size_t active_face_slot = 0; active_face_slot < active_set.active_faces.size(); ++active_face_slot) {
    const std::size_t face_index = active_set.active_faces[active_face_slot];
    const bool consumed_cache =
        primitive_cache != nullptr &&
        reconstruction.reconstructFaceFromCache(
            *primitive_cache,
            geometry.faces[face_index],
            scratch.left_states[active_face_slot],
            scratch.right_states[active_face_slot]);
    if (!consumed_cache) {
      reconstruction.reconstructFace(
          conserved,
          geometry.faces[face_index],
          scratch.left_states[active_face_slot],
          scratch.right_states[active_face_slot],
          m_adiabatic_index);
    }
  }
  const auto reconstruct_stop = std::chrono::steady_clock::now();

  const auto riemann_start = std::chrono::steady_clock::now();
  for (std::size_t active_face_slot = 0; active_face_slot < active_set.active_faces.size(); ++active_face_slot) {
    const std::size_t face_index = active_set.active_faces[active_face_slot];
    scratch.fluxes[active_face_slot] = riemann_solver.computeFlux(
        scratch.left_states[active_face_slot],
        scratch.right_states[active_face_slot],
        geometry.faces[face_index],
        m_adiabatic_index);
    maybeRecordFluxRegister(
        geometry,
        face_index,
        update,
        scratch.fluxes[active_face_slot],
        flux_register_sink);
  }
  const auto riemann_stop = std::chrono::steady_clock::now();

  const auto accumulate_start = std::chrono::steady_clock::now();
  const double flux_scale = update.dt_code / (update.scale_factor * geometry.cell_volume_comoving);
  for (std::size_t active_face_slot = 0; active_face_slot < active_set.active_faces.size(); ++active_face_slot) {
    const std::size_t face_index = active_set.active_faces[active_face_slot];
    const HydroFace& face = geometry.faces[face_index];
    const HydroConservedState face_delta = (flux_scale * face.area_comoving) * scratch.fluxes[active_face_slot];
    scratch.cell_delta[face.owner_cell] -= face_delta;
    const std::size_t neighbor_target = fluxNeighborTargetCell(geometry, face);
    if (neighbor_target != k_invalid_cell_index) {
      scratch.cell_delta[neighbor_target] += face_delta;
    }
  }
  const auto accumulate_stop = std::chrono::steady_clock::now();

  const auto source_start = std::chrono::steady_clock::now();
  for (std::size_t cell_index : scratch.touched_cells) {
    const HydroConservedState old_cell = conserved.loadCell(cell_index);
    HydroPrimitiveState primitive = primitive_cache != nullptr
        ? primitive_cache->loadCell(cell_index)
        : primitiveFromConserved(old_cell, m_adiabatic_index);

    HydroConservedState source_total;
    if (active_cell_mask[cell_index] != 0U) {
      for (const HydroSourceTerm* source : source_terms) {
        if (source == nullptr) {
          continue;
        }
        source_total += source->sourceForCell(cell_index, old_cell, primitive, source_context);
      }
    }

    const HydroConservedState flux_delta = scratch.cell_delta[cell_index];
    const HydroConservedState source_delta = update.dt_code * source_total;
    conservation_report.flux_delta += totalsFromDelta(
        flux_delta,
        old_cell,
        geometry.cell_volume_comoving);
    conservation_report.source_delta += totalsFromDelta(
        source_delta,
        old_cell + flux_delta,
        geometry.cell_volume_comoving);

    HydroConservedState updated = old_cell + flux_delta + source_delta;
    const HydroConservedState before_floor = updated;
    updated.mass_density_comoving = std::max(updated.mass_density_comoving, k_small);
    const double kinetic_density = 0.5 *
        (updated.momentum_density_x_comoving * updated.momentum_density_x_comoving +
         updated.momentum_density_y_comoving * updated.momentum_density_y_comoving +
         updated.momentum_density_z_comoving * updated.momentum_density_z_comoving) /
        std::max(updated.mass_density_comoving, k_small);
    const double minimum_total_energy_density = kinetic_density + k_small;
    const bool applied_internal_energy_floor = updated.total_energy_density_comoving < minimum_total_energy_density;
    if (applied_internal_energy_floor) {
      updated.total_energy_density_comoving = minimum_total_energy_density;
      ++conservation_report.internal_energy_floor_count;
    }
    if (updated.mass_density_comoving != before_floor.mass_density_comoving ||
        updated.total_energy_density_comoving != before_floor.total_energy_density_comoving) {
      conservation_report.floor_delta += totalsFromDelta(
          updated - before_floor,
          before_floor,
          geometry.cell_volume_comoving);
    }
    conserved.storeCell(cell_index, updated);
    if (primitive_cache != nullptr) {
      primitive_cache->storeCell(cell_index, primitiveFromConserved(updated, m_adiabatic_index));
    }
  }
  const auto source_stop = std::chrono::steady_clock::now();
  const auto total_stop = std::chrono::steady_clock::now();
  conservation_report.after = conservationTotals(conserved, geometry, scratch.full_active_cells);
  conservation_report.residual = conservation_report.after - conservation_report.before -
      conservation_report.flux_delta - conservation_report.source_delta - conservation_report.floor_delta;

  if (profile != nullptr) {
    profile->face_count += static_cast<std::uint64_t>(active_set.active_faces.size());
    const std::uint64_t limiter_after = reconstruction.limiterClipCount();
    const std::uint64_t positivity_after = reconstruction.positivityFallbackCount();
    const std::uint64_t riemann_fallback_after = riemann_solver.fallbackCount();
    profile->limiter_clip_count += (limiter_after >= limiter_before) ? (limiter_after - limiter_before) : 0U;
    profile->positivity_fallback_count +=
        (positivity_after >= positivity_before) ? (positivity_after - positivity_before) : 0U;
    profile->riemann_fallback_count +=
        (riemann_fallback_after >= riemann_fallback_before) ? (riemann_fallback_after - riemann_fallback_before) : 0U;
    profile->reconstruct_ms += std::chrono::duration<double, std::milli>(reconstruct_stop - reconstruct_start).count();
    profile->riemann_ms += std::chrono::duration<double, std::milli>(riemann_stop - riemann_start).count();
    profile->accumulate_ms += std::chrono::duration<double, std::milli>(accumulate_stop - accumulate_start).count();
    profile->source_ms += std::chrono::duration<double, std::milli>(source_stop - source_start).count();
    profile->total_ms += std::chrono::duration<double, std::milli>(total_stop - total_start).count();
    profile->bytes_moved += static_cast<std::uint64_t>(
        (6U * scratch.touched_cells.size() + 10U * active_set.active_faces.size()) * sizeof(double));
    profile->conservation = conservation_report;
  }
}

}  // namespace cosmosim::hydro
