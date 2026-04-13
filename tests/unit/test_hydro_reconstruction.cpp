#include <cassert>
#include <cmath>

#include "cosmosim/hydro/hydro_core_solver.hpp"
#include "cosmosim/hydro/hydro_reconstruction.hpp"

namespace {

constexpr double k_tol = 1.0e-12;

void testLimiterStringRoundTrip() {
  using cosmosim::hydro::HydroSlopeLimiter;
  using cosmosim::hydro::hydroSlopeLimiterFromString;
  using cosmosim::hydro::hydroSlopeLimiterToString;

  assert(hydroSlopeLimiterFromString(hydroSlopeLimiterToString(HydroSlopeLimiter::kMinmod)) == HydroSlopeLimiter::kMinmod);
  assert(hydroSlopeLimiterFromString(hydroSlopeLimiterToString(HydroSlopeLimiter::kMonotonizedCentral)) ==
         HydroSlopeLimiter::kMonotonizedCentral);
  assert(hydroSlopeLimiterFromString(hydroSlopeLimiterToString(HydroSlopeLimiter::kVanLeer)) == HydroSlopeLimiter::kVanLeer);
}

void testMusclPositivityFallbackCount() {
  constexpr double gamma = 1.4;
  cosmosim::hydro::HydroPrimitiveCacheSoa cache(3);
  cache.storeCell(0, {.rho_comoving = 1.0, .vel_x_peculiar = 0.0, .pressure_comoving = 1.0});
  cache.storeCell(1, {.rho_comoving = 1.0e-14, .vel_x_peculiar = 0.0, .pressure_comoving = 1.0e-14});
  cache.storeCell(2, {.rho_comoving = 1.0, .vel_x_peculiar = 0.0, .pressure_comoving = 1.0});

  cosmosim::hydro::MusclHancockReconstruction reconstruction(cosmosim::hydro::HydroReconstructionPolicy{
      .limiter = cosmosim::hydro::HydroSlopeLimiter::kMonotonizedCentral,
      .dt_over_dx_code = 0.0,
      .rho_floor = 1.0e-8,
      .pressure_floor = 1.0e-8,
      .enable_muscl_hancock_predictor = false,
      .adiabatic_index = gamma});

  const cosmosim::hydro::HydroFace face{.owner_cell = 0, .neighbor_cell = 1, .area_comoving = 1.0, .normal_x = 1.0};
  cosmosim::hydro::HydroPrimitiveState left;
  cosmosim::hydro::HydroPrimitiveState right;
  const bool consumed = reconstruction.reconstructFaceFromCache(cache, face, left, right);
  assert(consumed);
  assert(right.rho_comoving >= 1.0e-8 - k_tol);
  assert(right.pressure_comoving >= 1.0e-8 - k_tol);
  assert(reconstruction.positivityFallbackCount() >= 1U);
}

}  // namespace

int main() {
  testLimiterStringRoundTrip();
  testMusclPositivityFallbackCount();
  return 0;
}
