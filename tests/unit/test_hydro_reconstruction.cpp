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

void testMusclUsesExplicitYAxisStencilForAllPrimitives() {
  cosmosim::hydro::HydroPrimitiveCacheSoa cache(4);
  for (std::size_t i = 0; i < 4; ++i) {
    const double row = static_cast<double>(i);
    cache.storeCell(i, {
        .rho_comoving = 1.0 + row,
        .vel_x_peculiar = 10.0 + row,
        .vel_y_peculiar = 20.0 + 2.0 * row,
        .vel_z_peculiar = -5.0 + 0.5 * row,
        .pressure_comoving = 2.0 + row});
  }

  cosmosim::hydro::MusclHancockReconstruction reconstruction(cosmosim::hydro::HydroReconstructionPolicy{
      .limiter = cosmosim::hydro::HydroSlopeLimiter::kMonotonizedCentral,
      .dt_over_dx_code = 0.0,
      .rho_floor = 1.0e-10,
      .pressure_floor = 1.0e-10,
      .enable_muscl_hancock_predictor = false});

  const cosmosim::hydro::HydroFace face{
      .owner_cell = 1,
      .neighbor_cell = 2,
      .owner_minus_cell = 0,
      .neighbor_plus_cell = 3,
      .area_comoving = 1.0,
      .normal_y = 1.0,
      .axis = cosmosim::hydro::HydroFaceAxis::kY};
  cosmosim::hydro::HydroPrimitiveState left;
  cosmosim::hydro::HydroPrimitiveState right;
  assert(reconstruction.reconstructFaceFromCache(cache, face, left, right));

  assert(std::abs(left.rho_comoving - 2.5) < k_tol);
  assert(std::abs(right.rho_comoving - 2.5) < k_tol);
  assert(std::abs(left.vel_x_peculiar - 11.5) < k_tol);
  assert(std::abs(right.vel_x_peculiar - 11.5) < k_tol);
  assert(std::abs(left.vel_y_peculiar - 23.0) < k_tol);
  assert(std::abs(right.vel_y_peculiar - 23.0) < k_tol);
  assert(std::abs(left.vel_z_peculiar + 4.25) < k_tol);
  assert(std::abs(right.vel_z_peculiar + 4.25) < k_tol);
  assert(std::abs(left.pressure_comoving - 3.5) < k_tol);
  assert(std::abs(right.pressure_comoving - 3.5) < k_tol);
}

void testMusclZAxisPredictorUsesZCflAndNormalVelocity() {
  cosmosim::hydro::HydroPrimitiveCacheSoa cache(4);
  cache.storeCell(0, {.rho_comoving = 1.0, .vel_x_peculiar = 0.1, .vel_y_peculiar = 0.2, .vel_z_peculiar = 0.3, .pressure_comoving = 1.0});
  cache.storeCell(1, {.rho_comoving = 1.2, .vel_x_peculiar = 0.2, .vel_y_peculiar = 0.4, .vel_z_peculiar = 0.6, .pressure_comoving = 1.1});
  cache.storeCell(2, {.rho_comoving = 1.4, .vel_x_peculiar = 0.3, .vel_y_peculiar = 0.6, .vel_z_peculiar = 0.9, .pressure_comoving = 1.2});
  cache.storeCell(3, {.rho_comoving = 1.6, .vel_x_peculiar = 0.4, .vel_y_peculiar = 0.8, .vel_z_peculiar = 1.2, .pressure_comoving = 1.3});

  cosmosim::hydro::MusclHancockReconstruction no_predictor(cosmosim::hydro::HydroReconstructionPolicy{
      .limiter = cosmosim::hydro::HydroSlopeLimiter::kMonotonizedCentral,
      .dt_over_dx_code = 0.0,
      .rho_floor = 1.0e-10,
      .pressure_floor = 1.0e-10,
      .enable_muscl_hancock_predictor = false});
  cosmosim::hydro::MusclHancockReconstruction z_predictor(cosmosim::hydro::HydroReconstructionPolicy{
      .limiter = cosmosim::hydro::HydroSlopeLimiter::kMonotonizedCentral,
      .dt_over_dx_code = 0.0,
      .dt_over_cell_width_code = {0.0, 0.0, 0.2},
      .rho_floor = 1.0e-10,
      .pressure_floor = 1.0e-10,
      .enable_muscl_hancock_predictor = true});

  const cosmosim::hydro::HydroFace face{
      .owner_cell = 1,
      .neighbor_cell = 2,
      .owner_minus_cell = 0,
      .neighbor_plus_cell = 3,
      .area_comoving = 1.0,
      .normal_z = 1.0,
      .axis = cosmosim::hydro::HydroFaceAxis::kZ};
  cosmosim::hydro::HydroPrimitiveState left_no;
  cosmosim::hydro::HydroPrimitiveState right_no;
  cosmosim::hydro::HydroPrimitiveState left_pred;
  cosmosim::hydro::HydroPrimitiveState right_pred;
  assert(no_predictor.reconstructFaceFromCache(cache, face, left_no, right_no));
  assert(z_predictor.reconstructFaceFromCache(cache, face, left_pred, right_pred));

  assert(std::abs(left_pred.rho_comoving - left_no.rho_comoving) > 1.0e-4);
  assert(std::abs(left_pred.vel_z_peculiar - left_no.vel_z_peculiar) > 1.0e-4);
  assert(std::abs(left_pred.vel_x_peculiar - left_no.vel_x_peculiar) > 1.0e-4);
  assert(std::abs(right_pred.pressure_comoving - right_no.pressure_comoving) > 1.0e-4);
}

}  // namespace

int main() {
  testLimiterStringRoundTrip();
  testMusclPositivityFallbackCount();
  testMusclUsesExplicitYAxisStencilForAllPrimitives();
  testMusclZAxisPredictorUsesZCflAndNormalVelocity();
  return 0;
}
