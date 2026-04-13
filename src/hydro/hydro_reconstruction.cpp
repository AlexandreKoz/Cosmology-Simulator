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

[[nodiscard]] bool isContiguousPair(std::size_t left_index, std::size_t right_index, std::size_t cell_count) {
  return right_index == left_index + 1U || (left_index + 1U == cell_count && right_index == 0U);
}

[[nodiscard]] std::size_t wrappedPlusOne(std::size_t cell_index, std::size_t cell_count) {
  return (cell_index + 1U) % cell_count;
}

[[nodiscard]] std::size_t wrappedMinusOne(std::size_t cell_index, std::size_t cell_count) {
  return (cell_index + cell_count - 1U) % cell_count;
}

[[nodiscard]] double minmodLimiter(double a, double b) {
  if (a * b <= 0.0) {
    return 0.0;
  }
  return (std::abs(a) < std::abs(b)) ? a : b;
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

  const std::size_t cell_count = primitive_cache.size();
  if (!isContiguousPair(face.owner_cell, face.neighbor_cell, cell_count)) {
    return true;
  }

  const std::size_t left_minus_index = wrappedMinusOne(face.owner_cell, cell_count);
  const std::size_t right_plus_index = wrappedPlusOne(face.neighbor_cell, cell_count);

  const HydroPrimitiveState left_minus = primitive_cache.loadCell(left_minus_index);
  const HydroPrimitiveState right_plus = primitive_cache.loadCell(right_plus_index);

  auto limited = [&](double qm, double q, double qp) {
    const double dm = q - qm;
    const double dp = qp - q;
    const double slope = applyHydroSlopeLimiter(m_policy.limiter, dm, dp);
    if (std::abs(slope) < std::abs(0.5 * (dm + dp))) {
      ++m_limiter_clip_count;
    }
    return slope;
  };

  const double slope_rho_left = limited(left_minus.rho_comoving, left_state.rho_comoving, right_state.rho_comoving);
  const double slope_u_left = limited(left_minus.vel_x_peculiar, left_state.vel_x_peculiar, right_state.vel_x_peculiar);
  const double slope_p_left = limited(left_minus.pressure_comoving, left_state.pressure_comoving, right_state.pressure_comoving);

  const double slope_rho_right = limited(left_state.rho_comoving, right_state.rho_comoving, right_plus.rho_comoving);
  const double slope_u_right = limited(left_state.vel_x_peculiar, right_state.vel_x_peculiar, right_plus.vel_x_peculiar);
  const double slope_p_right = limited(left_state.pressure_comoving, right_state.pressure_comoving, right_plus.pressure_comoving);

  left_state.rho_comoving += 0.5 * slope_rho_left;
  left_state.vel_x_peculiar += 0.5 * slope_u_left;
  left_state.pressure_comoving += 0.5 * slope_p_left;

  right_state.rho_comoving -= 0.5 * slope_rho_right;
  right_state.vel_x_peculiar -= 0.5 * slope_u_right;
  right_state.pressure_comoving -= 0.5 * slope_p_right;

  if (m_policy.enable_muscl_hancock_predictor && m_policy.dt_over_dx_code > 0.0) {
    const double cfl = m_policy.dt_over_dx_code;

    const double drho_dt_left = -left_state.vel_x_peculiar * slope_rho_left - left_state.rho_comoving * slope_u_left;
    const double du_dt_left = -left_state.vel_x_peculiar * slope_u_left - slope_p_left / std::max(left_state.rho_comoving, k_small);
    const double dp_dt_left = -left_state.vel_x_peculiar * slope_p_left - m_policy.adiabatic_index * left_state.pressure_comoving * slope_u_left;

    const double drho_dt_right = -right_state.vel_x_peculiar * slope_rho_right - right_state.rho_comoving * slope_u_right;
    const double du_dt_right = -right_state.vel_x_peculiar * slope_u_right - slope_p_right / std::max(right_state.rho_comoving, k_small);
    const double dp_dt_right = -right_state.vel_x_peculiar * slope_p_right - m_policy.adiabatic_index * right_state.pressure_comoving * slope_u_right;

    left_state.rho_comoving += -0.5 * cfl * drho_dt_left;
    left_state.vel_x_peculiar += -0.5 * cfl * du_dt_left;
    left_state.pressure_comoving += -0.5 * cfl * dp_dt_left;

    right_state.rho_comoving += -0.5 * cfl * drho_dt_right;
    right_state.vel_x_peculiar += -0.5 * cfl * du_dt_right;
    right_state.pressure_comoving += -0.5 * cfl * dp_dt_right;
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
