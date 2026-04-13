#include "cosmosim/hydro/hydro_riemann.hpp"

#include <algorithm>
#include <cmath>

#include "cosmosim/hydro/hydro_core_solver.hpp"

namespace cosmosim::hydro {
namespace {

constexpr double k_small = 1.0e-14;

[[nodiscard]] double dot3(double ax, double ay, double az, double bx, double by, double bz) {
  return ax * bx + ay * by + az * bz;
}

[[nodiscard]] double soundSpeed(double rho_comoving, double pressure_comoving, double adiabatic_index) {
  return std::sqrt(adiabatic_index * std::max(pressure_comoving, k_small) / std::max(rho_comoving, k_small));
}

[[nodiscard]] HydroConservedState eulerPhysicalFlux(
    const HydroPrimitiveState& primitive,
    const HydroFace& face,
    double adiabatic_index) {
  const double mass_density = primitive.rho_comoving;
  const double vel_n = dot3(
      primitive.vel_x_peculiar,
      primitive.vel_y_peculiar,
      primitive.vel_z_peculiar,
      face.normal_x,
      face.normal_y,
      face.normal_z);

  HydroConservedState flux;
  flux.mass_density_comoving = mass_density * vel_n;
  flux.momentum_density_x_comoving = mass_density * primitive.vel_x_peculiar * vel_n + primitive.pressure_comoving * face.normal_x;
  flux.momentum_density_y_comoving = mass_density * primitive.vel_y_peculiar * vel_n + primitive.pressure_comoving * face.normal_y;
  flux.momentum_density_z_comoving = mass_density * primitive.vel_z_peculiar * vel_n + primitive.pressure_comoving * face.normal_z;

  const HydroConservedState conserved = HydroCoreSolver::conservedFromPrimitive(primitive, adiabatic_index);
  flux.total_energy_density_comoving = (conserved.total_energy_density_comoving + primitive.pressure_comoving) * vel_n;
  return flux;
}

[[nodiscard]] HydroConservedState hllFlux(
    const HydroPrimitiveState& left_state,
    const HydroPrimitiveState& right_state,
    const HydroFace& face,
    double adiabatic_index) {
  const HydroConservedState u_left = HydroCoreSolver::conservedFromPrimitive(left_state, adiabatic_index);
  const HydroConservedState u_right = HydroCoreSolver::conservedFromPrimitive(right_state, adiabatic_index);

  const HydroConservedState f_left = eulerPhysicalFlux(left_state, face, adiabatic_index);
  const HydroConservedState f_right = eulerPhysicalFlux(right_state, face, adiabatic_index);

  const double vel_n_left = dot3(left_state.vel_x_peculiar, left_state.vel_y_peculiar, left_state.vel_z_peculiar, face.normal_x, face.normal_y, face.normal_z);
  const double vel_n_right = dot3(right_state.vel_x_peculiar, right_state.vel_y_peculiar, right_state.vel_z_peculiar, face.normal_x, face.normal_y, face.normal_z);

  const double c_left = soundSpeed(left_state.rho_comoving, left_state.pressure_comoving, adiabatic_index);
  const double c_right = soundSpeed(right_state.rho_comoving, right_state.pressure_comoving, adiabatic_index);

  const double s_left = std::min(vel_n_left - c_left, vel_n_right - c_right);
  const double s_right = std::max(vel_n_left + c_left, vel_n_right + c_right);

  if (s_left >= 0.0) {
    return f_left;
  }
  if (s_right <= 0.0) {
    return f_right;
  }

  const double inv_speed_delta = 1.0 / std::max(s_right - s_left, k_small);
  HydroConservedState flux;
  flux.mass_density_comoving =
      (s_right * f_left.mass_density_comoving - s_left * f_right.mass_density_comoving +
       s_left * s_right * (u_right.mass_density_comoving - u_left.mass_density_comoving)) *
      inv_speed_delta;
  flux.momentum_density_x_comoving =
      (s_right * f_left.momentum_density_x_comoving - s_left * f_right.momentum_density_x_comoving +
       s_left * s_right * (u_right.momentum_density_x_comoving - u_left.momentum_density_x_comoving)) *
      inv_speed_delta;
  flux.momentum_density_y_comoving =
      (s_right * f_left.momentum_density_y_comoving - s_left * f_right.momentum_density_y_comoving +
       s_left * s_right * (u_right.momentum_density_y_comoving - u_left.momentum_density_y_comoving)) *
      inv_speed_delta;
  flux.momentum_density_z_comoving =
      (s_right * f_left.momentum_density_z_comoving - s_left * f_right.momentum_density_z_comoving +
       s_left * s_right * (u_right.momentum_density_z_comoving - u_left.momentum_density_z_comoving)) *
      inv_speed_delta;
  flux.total_energy_density_comoving =
      (s_right * f_left.total_energy_density_comoving - s_left * f_right.total_energy_density_comoving +
       s_left * s_right * (u_right.total_energy_density_comoving - u_left.total_energy_density_comoving)) *
      inv_speed_delta;
  return flux;
}

}  // namespace

HydroConservedState HlleRiemannSolver::computeFlux(
    const HydroPrimitiveState& left_state,
    const HydroPrimitiveState& right_state,
    const HydroFace& face,
    double adiabatic_index) const {
  return hllFlux(left_state, right_state, face, adiabatic_index);
}

HydroConservedState HllcRiemannSolver::computeFlux(
    const HydroPrimitiveState& left_state,
    const HydroPrimitiveState& right_state,
    const HydroFace& face,
    double adiabatic_index) const {
  const double un_left = dot3(left_state.vel_x_peculiar, left_state.vel_y_peculiar, left_state.vel_z_peculiar,
      face.normal_x, face.normal_y, face.normal_z);
  const double un_right = dot3(right_state.vel_x_peculiar, right_state.vel_y_peculiar, right_state.vel_z_peculiar,
      face.normal_x, face.normal_y, face.normal_z);

  const double c_left = soundSpeed(left_state.rho_comoving, left_state.pressure_comoving, adiabatic_index);
  const double c_right = soundSpeed(right_state.rho_comoving, right_state.pressure_comoving, adiabatic_index);

  const double s_left = std::min(un_left - c_left, un_right - c_right);
  const double s_right = std::max(un_left + c_left, un_right + c_right);

  if (s_left >= 0.0) {
    return eulerPhysicalFlux(left_state, face, adiabatic_index);
  }
  if (s_right <= 0.0) {
    return eulerPhysicalFlux(right_state, face, adiabatic_index);
  }

  const HydroConservedState u_left = HydroCoreSolver::conservedFromPrimitive(left_state, adiabatic_index);
  const HydroConservedState u_right = HydroCoreSolver::conservedFromPrimitive(right_state, adiabatic_index);
  const HydroConservedState f_left = eulerPhysicalFlux(left_state, face, adiabatic_index);
  const HydroConservedState f_right = eulerPhysicalFlux(right_state, face, adiabatic_index);

  const double numerator = right_state.pressure_comoving - left_state.pressure_comoving +
      left_state.rho_comoving * un_left * (s_left - un_left) -
      right_state.rho_comoving * un_right * (s_right - un_right);
  const double denominator = left_state.rho_comoving * (s_left - un_left) -
      right_state.rho_comoving * (s_right - un_right);
  if (std::abs(denominator) < k_small) {
    ++m_fallback_count;
    return hllFlux(left_state, right_state, face, adiabatic_index);
  }
  const double s_star = numerator / denominator;

  const double p_star_left = left_state.pressure_comoving + left_state.rho_comoving * (s_left - un_left) * (s_star - un_left);
  const double p_star_right = right_state.pressure_comoving + right_state.rho_comoving * (s_right - un_right) * (s_star - un_right);
  const double p_star = 0.5 * (p_star_left + p_star_right);
  if (p_star <= 0.0) {
    ++m_fallback_count;
    return hllFlux(left_state, right_state, face, adiabatic_index);
  }

  auto starState = [&](const HydroPrimitiveState& prim, const HydroConservedState& cons, double s_k, double u_k_n) {
    const double denom = s_k - s_star;
    if (std::abs(denom) < k_small) {
      return HydroConservedState{};
    }
    const double rho_star = prim.rho_comoving * (s_k - u_k_n) / denom;
    if (!std::isfinite(rho_star) || rho_star <= 0.0) {
      return HydroConservedState{};
    }
    HydroConservedState star;
    star.mass_density_comoving = rho_star;

    const double velocity_projection = dot3(prim.vel_x_peculiar, prim.vel_y_peculiar, prim.vel_z_peculiar,
        face.normal_x, face.normal_y, face.normal_z);
    const double delta_normal = s_star - velocity_projection;
    star.momentum_density_x_comoving = rho_star * (prim.vel_x_peculiar + delta_normal * face.normal_x);
    star.momentum_density_y_comoving = rho_star * (prim.vel_y_peculiar + delta_normal * face.normal_y);
    star.momentum_density_z_comoving = rho_star * (prim.vel_z_peculiar + delta_normal * face.normal_z);

    const double denom_energy = prim.rho_comoving * (s_k - u_k_n);
    if (std::abs(denom_energy) < k_small) {
      return HydroConservedState{};
    }
    const double energy_num =
        (s_k - u_k_n) * cons.total_energy_density_comoving - prim.pressure_comoving * u_k_n + p_star * s_star;
    star.total_energy_density_comoving = rho_star * energy_num / denom_energy;
    if (!std::isfinite(star.total_energy_density_comoving) || star.total_energy_density_comoving <= 0.0) {
      return HydroConservedState{};
    }
    return star;
  };

  const HydroConservedState u_star_left = starState(left_state, u_left, s_left, un_left);
  const HydroConservedState u_star_right = starState(right_state, u_right, s_right, un_right);

  if (u_star_left.mass_density_comoving <= 0.0 || u_star_right.mass_density_comoving <= 0.0 ||
      !std::isfinite(u_star_left.total_energy_density_comoving) || !std::isfinite(u_star_right.total_energy_density_comoving)) {
    ++m_fallback_count;
    return hllFlux(left_state, right_state, face, adiabatic_index);
  }

  if (s_star >= 0.0) {
    return f_left + s_left * (u_star_left - u_left);
  }
  return f_right + s_right * (u_star_right - u_right);
}

std::uint64_t HllcRiemannSolver::fallbackCount() const {
  return m_fallback_count;
}

}  // namespace cosmosim::hydro
