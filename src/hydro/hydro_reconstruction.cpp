#include "cosmosim/hydro/hydro_reconstruction.hpp"

#include <cmath>
#include <stdexcept>
#include <string>

#include "cosmosim/hydro/hydro_core_solver.hpp"

namespace cosmosim::hydro {
namespace {

constexpr double k_small = 1.0e-14;

[[nodiscard]] bool enforcePrimitiveFloors(HydroPrimitiveState& state, double rho_floor, double pressure_floor) {
  bool changed = false;
  if (state.rho_comoving < rho_floor) {
    state.rho_comoving = rho_floor;
    changed = true;
  }
  if (state.pressure_comoving < pressure_floor) {
    state.pressure_comoving = pressure_floor;
    changed = true;
  }
  return changed;
}

[[nodiscard]] double minmodLimiter(double a, double b) {
  if (a * b <= 0.0) {
    return 0.0;
  }
  return (std::abs(a) < std::abs(b)) ? a : b;
}

[[nodiscard]] std::size_t axisIndex(HydroFaceAxis axis) {
  switch (axis) {
    case HydroFaceAxis::kX:
      return 0U;
    case HydroFaceAxis::kY:
      return 1U;
    case HydroFaceAxis::kZ:
      return 2U;
  }
  return 0U;
}

[[nodiscard]] double normalComponent(const HydroPrimitiveState& state, HydroFaceAxis axis) {
  switch (axis) {
    case HydroFaceAxis::kX:
      return state.vel_x_peculiar;
    case HydroFaceAxis::kY:
      return state.vel_y_peculiar;
    case HydroFaceAxis::kZ:
      return state.vel_z_peculiar;
  }
  return state.vel_x_peculiar;
}

void addScaledSlope(HydroPrimitiveState& state, const HydroPrimitiveState& slope, double scale) {
  state.rho_comoving += scale * slope.rho_comoving;
  state.vel_x_peculiar += scale * slope.vel_x_peculiar;
  state.vel_y_peculiar += scale * slope.vel_y_peculiar;
  state.vel_z_peculiar += scale * slope.vel_z_peculiar;
  state.pressure_comoving += scale * slope.pressure_comoving;
}

[[nodiscard]] double dtOverCellWidth(const HydroReconstructionPolicy& policy, HydroFaceAxis axis) {
  const double axis_value = policy.dt_over_cell_width_code[axisIndex(axis)];
  return axis_value > 0.0 ? axis_value : policy.dt_over_dx_code;
}

[[nodiscard]] double faceNormalSign(const HydroFace& face) {
  switch (face.axis) {
    case HydroFaceAxis::kX:
      return face.normal_x >= 0.0 ? 1.0 : -1.0;
    case HydroFaceAxis::kY:
      return face.normal_y >= 0.0 ? 1.0 : -1.0;
    case HydroFaceAxis::kZ:
      return face.normal_z >= 0.0 ? 1.0 : -1.0;
  }
  return 1.0;
}

}  // namespace

std::string_view hydroSlopeLimiterToString(HydroSlopeLimiter limiter) {
  switch (limiter) {
    case HydroSlopeLimiter::kMinmod:
      return "minmod";
    case HydroSlopeLimiter::kMonotonizedCentral:
      return "mc";
    case HydroSlopeLimiter::kVanLeer:
      return "van_leer";
  }
  return "unknown";
}

HydroSlopeLimiter hydroSlopeLimiterFromString(std::string_view name) {
  if (name == "minmod") {
    return HydroSlopeLimiter::kMinmod;
  }
  if (name == "mc" || name == "monotonized_central") {
    return HydroSlopeLimiter::kMonotonizedCentral;
  }
  if (name == "van_leer") {
    return HydroSlopeLimiter::kVanLeer;
  }
  throw std::invalid_argument("Unsupported hydro slope limiter: " + std::string(name));
}

double applyHydroSlopeLimiter(HydroSlopeLimiter limiter, double delta_minus, double delta_plus) {
  switch (limiter) {
    case HydroSlopeLimiter::kMinmod:
      return minmodLimiter(delta_minus, delta_plus);
    case HydroSlopeLimiter::kMonotonizedCentral: {
      const double mm = minmodLimiter(delta_minus, delta_plus);
      return minmodLimiter(0.5 * (delta_minus + delta_plus), 2.0 * mm);
    }
    case HydroSlopeLimiter::kVanLeer: {
      const double product = delta_minus * delta_plus;
      if (product <= 0.0) {
        return 0.0;
      }
      return 2.0 * product / (delta_minus + delta_plus);
    }
  }
  return 0.0;
}

bool HydroReconstruction::reconstructFaceFromCache(
    const HydroPrimitiveCacheSoa& primitive_cache,
    const HydroFace& face,
    HydroPrimitiveState& left_state,
    HydroPrimitiveState& right_state) const {
  (void)primitive_cache;
  (void)face;
  (void)left_state;
  (void)right_state;
  return false;
}

void PiecewiseConstantReconstruction::reconstructFace(
    const HydroConservedStateSoa& conserved,
    const HydroFace& face,
    HydroPrimitiveState& left_state,
    HydroPrimitiveState& right_state,
    double adiabatic_index) const {
  left_state = HydroCoreSolver::primitiveFromConserved(conserved.loadCell(face.owner_cell), adiabatic_index);
  if (face.neighbor_cell == k_invalid_cell_index) {
    right_state = left_state;
    return;
  }
  right_state = HydroCoreSolver::primitiveFromConserved(conserved.loadCell(face.neighbor_cell), adiabatic_index);
}

bool PiecewiseConstantReconstruction::reconstructFaceFromCache(
    const HydroPrimitiveCacheSoa& primitive_cache,
    const HydroFace& face,
    HydroPrimitiveState& left_state,
    HydroPrimitiveState& right_state) const {
  left_state = primitive_cache.loadCell(face.owner_cell);
  if (face.neighbor_cell == k_invalid_cell_index) {
    right_state = left_state;
    return true;
  }
  right_state = primitive_cache.loadCell(face.neighbor_cell);
  return true;
}

MusclHancockReconstruction::MusclHancockReconstruction(HydroReconstructionPolicy policy) : m_policy(policy) {
  if (m_policy.rho_floor <= 0.0 || m_policy.pressure_floor <= 0.0) {
    throw std::invalid_argument("MUSCL-Hancock reconstruction requires positive floors");
  }
  if (m_policy.dt_over_dx_code < 0.0 ||
      m_policy.dt_over_cell_width_code[0] < 0.0 ||
      m_policy.dt_over_cell_width_code[1] < 0.0 ||
      m_policy.dt_over_cell_width_code[2] < 0.0) {
    throw std::invalid_argument("MUSCL-Hancock reconstruction requires non-negative CFL policy values");
  }
}

bool MusclHancockReconstruction::reconstructFaceFromCache(
    const HydroPrimitiveCacheSoa& primitive_cache,
    const HydroFace& face,
    HydroPrimitiveState& left_state,
    HydroPrimitiveState& right_state) const {
  left_state = primitive_cache.loadCell(face.owner_cell);
  if (face.neighbor_cell == k_invalid_cell_index) {
    right_state = left_state;
    return true;
  }
  right_state = primitive_cache.loadCell(face.neighbor_cell);

  const HydroPrimitiveState left_minus = face.owner_minus_cell == k_invalid_cell_index
      ? left_state
      : primitive_cache.loadCell(face.owner_minus_cell);
  const HydroPrimitiveState right_plus = face.neighbor_plus_cell == k_invalid_cell_index
      ? right_state
      : primitive_cache.loadCell(face.neighbor_plus_cell);

  auto limited = [&](double qm, double q, double qp) {
    const double dm = q - qm;
    const double dp = qp - q;
    const double slope = applyHydroSlopeLimiter(m_policy.limiter, dm, dp);
    if (std::abs(slope) < std::abs(0.5 * (dm + dp))) {
      ++m_limiter_clip_count;
    }
    return slope;
  };

  const HydroPrimitiveState slope_left{
      .rho_comoving = limited(left_minus.rho_comoving, left_state.rho_comoving, right_state.rho_comoving),
      .vel_x_peculiar = limited(left_minus.vel_x_peculiar, left_state.vel_x_peculiar, right_state.vel_x_peculiar),
      .vel_y_peculiar = limited(left_minus.vel_y_peculiar, left_state.vel_y_peculiar, right_state.vel_y_peculiar),
      .vel_z_peculiar = limited(left_minus.vel_z_peculiar, left_state.vel_z_peculiar, right_state.vel_z_peculiar),
      .pressure_comoving = limited(left_minus.pressure_comoving, left_state.pressure_comoving, right_state.pressure_comoving)};
  const HydroPrimitiveState slope_right{
      .rho_comoving = limited(left_state.rho_comoving, right_state.rho_comoving, right_plus.rho_comoving),
      .vel_x_peculiar = limited(left_state.vel_x_peculiar, right_state.vel_x_peculiar, right_plus.vel_x_peculiar),
      .vel_y_peculiar = limited(left_state.vel_y_peculiar, right_state.vel_y_peculiar, right_plus.vel_y_peculiar),
      .vel_z_peculiar = limited(left_state.vel_z_peculiar, right_state.vel_z_peculiar, right_plus.vel_z_peculiar),
      .pressure_comoving = limited(left_state.pressure_comoving, right_state.pressure_comoving, right_plus.pressure_comoving)};

  const double orientation = faceNormalSign(face);
  addScaledSlope(left_state, slope_left, 0.5 * orientation);
  addScaledSlope(right_state, slope_right, -0.5 * orientation);

  const double cfl = dtOverCellWidth(m_policy, face.axis);
  if (m_policy.enable_muscl_hancock_predictor && cfl > 0.0) {
    auto predict = [&](HydroPrimitiveState& state, const HydroPrimitiveState& slope) {
      const double u_normal = normalComponent(state, face.axis);
      const double slope_u_normal = normalComponent(slope, face.axis);
      const double inv_rho = 1.0 / std::max(state.rho_comoving, k_small);

      const double drho_dt = -u_normal * slope.rho_comoving - state.rho_comoving * slope_u_normal;
      const double dux_dt = -u_normal * slope.vel_x_peculiar;
      const double duy_dt = -u_normal * slope.vel_y_peculiar;
      const double duz_dt = -u_normal * slope.vel_z_peculiar;
      const double dp_dt = -u_normal * slope.pressure_comoving -
          m_policy.adiabatic_index * state.pressure_comoving * slope_u_normal;

      state.rho_comoving += -0.5 * cfl * drho_dt;
      state.vel_x_peculiar += -0.5 * cfl * (face.axis == HydroFaceAxis::kX
          ? dux_dt - slope.pressure_comoving * inv_rho
          : dux_dt);
      state.vel_y_peculiar += -0.5 * cfl * (face.axis == HydroFaceAxis::kY
          ? duy_dt - slope.pressure_comoving * inv_rho
          : duy_dt);
      state.vel_z_peculiar += -0.5 * cfl * (face.axis == HydroFaceAxis::kZ
          ? duz_dt - slope.pressure_comoving * inv_rho
          : duz_dt);
      state.pressure_comoving += -0.5 * cfl * dp_dt;
    };

    predict(left_state, slope_left);
    predict(right_state, slope_right);
  }

  if (enforcePrimitiveFloors(left_state, m_policy.rho_floor, m_policy.pressure_floor)) {
    ++m_positivity_fallback_count;
  }
  if (enforcePrimitiveFloors(right_state, m_policy.rho_floor, m_policy.pressure_floor)) {
    ++m_positivity_fallback_count;
  }
  return true;
}

void MusclHancockReconstruction::reconstructFace(
    const HydroConservedStateSoa& conserved,
    const HydroFace& face,
    HydroPrimitiveState& left_state,
    HydroPrimitiveState& right_state,
    double adiabatic_index) const {
  left_state = HydroCoreSolver::primitiveFromConserved(conserved.loadCell(face.owner_cell), adiabatic_index);
  if (face.neighbor_cell == k_invalid_cell_index) {
    right_state = left_state;
    return;
  }
  right_state = HydroCoreSolver::primitiveFromConserved(conserved.loadCell(face.neighbor_cell), adiabatic_index);
}

HydroReconstructionPolicy MusclHancockReconstruction::policy() const { return m_policy; }
std::uint64_t MusclHancockReconstruction::limiterClipCount() const { return m_limiter_clip_count; }
std::uint64_t MusclHancockReconstruction::positivityFallbackCount() const { return m_positivity_fallback_count; }

}  // namespace cosmosim::hydro
