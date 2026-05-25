#include <cassert>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "cosmosim/hydro/hydro_core_solver.hpp"

namespace {

constexpr double k_tol = 1.0e-10;

void testLimiterLibrary() {
  using cosmosim::hydro::HydroSlopeLimiter;
  using cosmosim::hydro::applyHydroSlopeLimiter;

  const double mm = applyHydroSlopeLimiter(HydroSlopeLimiter::kMinmod, 2.0, 1.0);
  const double mc = applyHydroSlopeLimiter(HydroSlopeLimiter::kMonotonizedCentral, 2.0, 1.0);
  const double vl = applyHydroSlopeLimiter(HydroSlopeLimiter::kVanLeer, 2.0, 1.0);
  assert(std::abs(mm - 1.0) < k_tol);
  assert(mc >= mm && mc <= 1.5);
  assert(vl > mm && vl < 2.0);
  assert(std::abs(applyHydroSlopeLimiter(HydroSlopeLimiter::kVanLeer, -1.0, 2.0)) < k_tol);
}

void testPrimitiveConservedRoundTrip() {
  cosmosim::hydro::HydroPrimitiveState primitive;
  primitive.rho_comoving = 2.5;
  primitive.vel_x_peculiar = 1.2;
  primitive.vel_y_peculiar = -0.5;
  primitive.vel_z_peculiar = 0.75;
  primitive.pressure_comoving = 3.0;

  const double gamma = 5.0 / 3.0;
  const cosmosim::hydro::HydroConservedState conserved =
      cosmosim::hydro::HydroCoreSolver::conservedFromPrimitive(primitive, gamma);
  const cosmosim::hydro::HydroPrimitiveState round_trip =
      cosmosim::hydro::HydroCoreSolver::primitiveFromConserved(conserved, gamma);

  assert(std::abs(round_trip.rho_comoving - primitive.rho_comoving) < k_tol);
  assert(std::abs(round_trip.vel_x_peculiar - primitive.vel_x_peculiar) < k_tol);
  assert(std::abs(round_trip.vel_y_peculiar - primitive.vel_y_peculiar) < k_tol);
  assert(std::abs(round_trip.vel_z_peculiar - primitive.vel_z_peculiar) < k_tol);
  assert(std::abs(round_trip.pressure_comoving - primitive.pressure_comoving) < k_tol);
}

void testMusclReconstructionProducesFiniteStates() {
  constexpr double gamma = 1.4;
  cosmosim::hydro::HydroConservedStateSoa conserved(5);
  for (std::size_t i = 0; i < 5; ++i) {
    cosmosim::hydro::HydroPrimitiveState primitive;
    primitive.rho_comoving = 1.0 + 0.2 * static_cast<double>(i);
    primitive.vel_x_peculiar = 0.05 * static_cast<double>(i);
    primitive.pressure_comoving = 1.0 + 0.1 * static_cast<double>(i);
    conserved.storeCell(i, cosmosim::hydro::HydroCoreSolver::conservedFromPrimitive(primitive, gamma));
  }

  cosmosim::hydro::HydroPrimitiveCacheSoa cache(5);
  for (std::size_t i = 0; i < 5; ++i) {
    cache.storeCell(i, cosmosim::hydro::HydroCoreSolver::primitiveFromConserved(conserved.loadCell(i), gamma));
  }

  cosmosim::hydro::MusclHancockReconstruction reconstruction(cosmosim::hydro::HydroReconstructionPolicy{
      .limiter = cosmosim::hydro::HydroSlopeLimiter::kMonotonizedCentral,
      .dt_over_dx_code = 0.1,
      .rho_floor = 1.0e-8,
      .pressure_floor = 1.0e-8,
      .enable_muscl_hancock_predictor = true});

  cosmosim::hydro::HydroFace face{.owner_cell = 2, .neighbor_cell = 3, .area_comoving = 1.0, .normal_x = 1.0};
  cosmosim::hydro::HydroPrimitiveState left;
  cosmosim::hydro::HydroPrimitiveState right;
  const bool consumed = reconstruction.reconstructFaceFromCache(cache, face, left, right);
  assert(consumed);
  assert(left.rho_comoving > 0.0 && right.rho_comoving > 0.0);
  assert(left.pressure_comoving > 0.0 && right.pressure_comoving > 0.0);
}

void testComovingSourceTermSanity() {
  cosmosim::hydro::HydroPrimitiveState primitive;
  primitive.rho_comoving = 1.5;
  primitive.vel_x_peculiar = 2.0;
  primitive.vel_y_peculiar = -1.0;
  primitive.vel_z_peculiar = 0.5;
  primitive.pressure_comoving = 1.2;

  const double gamma = 5.0 / 3.0;
  const cosmosim::hydro::HydroConservedState conserved =
      cosmosim::hydro::HydroCoreSolver::conservedFromPrimitive(primitive, gamma);

  const std::vector<double> gravity_x{0.2};
  const std::vector<double> gravity_y{-0.1};
  const std::vector<double> gravity_z{0.05};

  cosmosim::hydro::HydroSourceContext context;
  context.update.dt_code = 0.01;
  context.update.scale_factor = 0.8;
  context.update.hubble_rate_code = 0.4;
  context.gravity_accel_x_peculiar = gravity_x;
  context.gravity_accel_y_peculiar = gravity_y;
  context.gravity_accel_z_peculiar = gravity_z;

  cosmosim::hydro::ComovingGravityExpansionSource source;
  const cosmosim::hydro::HydroConservedState source_state =
      source.sourceForCell(0, conserved, primitive, context);

  assert(std::abs(source_state.mass_density_comoving) < k_tol);
  assert(source_state.momentum_density_x_comoving < primitive.rho_comoving * gravity_x[0]);
}

void testRiemannSymmetryRegression() {
  const cosmosim::hydro::HydroPrimitiveState left{.rho_comoving = 1.0, .pressure_comoving = 1.0};
  const cosmosim::hydro::HydroPrimitiveState right = left;
  const cosmosim::hydro::HydroFace face{.owner_cell = 0, .neighbor_cell = 1, .area_comoving = 1.0, .normal_x = 1.0};
  cosmosim::hydro::HllcRiemannSolver solver;
  const auto flux = solver.computeFlux(left, right, face, 1.4);
  assert(std::abs(flux.mass_density_comoving) < k_tol);
  assert(std::abs(flux.momentum_density_x_comoving - 1.0) < 1.0e-8);
}

void testProfileFallbackCountersAreStepLocal() {
  constexpr double gamma = 1.4;
  cosmosim::hydro::HydroConservedStateSoa conserved(16);
  for (std::size_t i = 0; i < conserved.size(); ++i) {
    cosmosim::hydro::HydroPrimitiveState primitive;
    primitive.rho_comoving = (i < 8) ? 1.0 : 0.125;
    primitive.pressure_comoving = (i < 8) ? 1.0 : 0.1;
    conserved.storeCell(i, cosmosim::hydro::HydroCoreSolver::conservedFromPrimitive(primitive, gamma));
  }

  cosmosim::hydro::HydroPatchGeometry geometry;
  geometry.cell_volume_comoving = 1.0;
  for (std::size_t i = 0; i < conserved.size(); ++i) {
    geometry.faces.push_back(cosmosim::hydro::HydroFace{
        .owner_cell = i,
        .neighbor_cell = (i + 1U) % conserved.size(),
        .area_comoving = 1.0,
        .normal_x = 1.0});
  }

  cosmosim::hydro::HydroUpdateContext update{.dt_code = 1.0e-3, .scale_factor = 1.0, .hubble_rate_code = 0.0};
  cosmosim::hydro::HydroSourceContext source_context{.update = update};
  cosmosim::hydro::HydroCoreSolver solver(gamma);
  cosmosim::hydro::MusclHancockReconstruction reconstruction;
  cosmosim::hydro::HllcRiemannSolver riemann;
  cosmosim::hydro::HydroProfileEvent profile;

  solver.advancePatch(conserved, geometry, update, reconstruction, riemann, {}, source_context, &profile);
  const std::uint64_t first_limiter = profile.limiter_clip_count;

  solver.advancePatch(conserved, geometry, update, reconstruction, riemann, {}, source_context, &profile);
  const std::uint64_t second_increment = profile.limiter_clip_count - first_limiter;
  assert(second_increment < 10U * static_cast<std::uint64_t>(geometry.faces.size()));
}

void testActiveSetFullCoverageMatchesFullAdvance() {
  constexpr double gamma = 1.4;
  cosmosim::hydro::HydroConservedStateSoa full_state(6);
  cosmosim::hydro::HydroConservedStateSoa active_state(6);
  for (std::size_t i = 0; i < 6; ++i) {
    cosmosim::hydro::HydroPrimitiveState primitive;
    primitive.rho_comoving = 1.0 + 0.1 * static_cast<double>(i);
    primitive.vel_x_peculiar = 0.01 * static_cast<double>(i);
    primitive.pressure_comoving = 1.0 + 0.05 * static_cast<double>(i);
    const auto conserved = cosmosim::hydro::HydroCoreSolver::conservedFromPrimitive(primitive, gamma);
    full_state.storeCell(i, conserved);
    active_state.storeCell(i, conserved);
  }

  cosmosim::hydro::HydroPatchGeometry geometry;
  geometry.cell_volume_comoving = 1.0;
  for (std::size_t i = 0; i < full_state.size(); ++i) {
    geometry.faces.push_back(cosmosim::hydro::HydroFace{
        .owner_cell = i,
        .neighbor_cell = (i + 1U) % full_state.size(),
        .area_comoving = 1.0,
        .normal_x = 1.0});
  }

  cosmosim::hydro::HydroUpdateContext update{.dt_code = 5.0e-4, .scale_factor = 1.0, .hubble_rate_code = 0.0};
  cosmosim::hydro::HydroSourceContext source_context{.update = update};
  cosmosim::hydro::HydroCoreSolver solver(gamma);
  cosmosim::hydro::MusclHancockReconstruction reconstruction;
  cosmosim::hydro::HllcRiemannSolver riemann;

  solver.advancePatch(full_state, geometry, update, reconstruction, riemann, {}, source_context, nullptr);

  std::vector<std::size_t> active_cells(full_state.size());
  std::vector<std::size_t> active_faces(geometry.faces.size());
  for (std::size_t i = 0; i < active_cells.size(); ++i) active_cells[i] = i;
  for (std::size_t i = 0; i < active_faces.size(); ++i) active_faces[i] = i;

  cosmosim::hydro::HydroScratchBuffers scratch;
  cosmosim::hydro::HydroPrimitiveCacheSoa cache(active_state.size());
  solver.advancePatchActiveSetWithScratch(
      active_state,
      geometry,
      cosmosim::hydro::HydroActiveSetView{.active_cells = active_cells, .active_faces = active_faces},
      update,
      reconstruction,
      riemann,
      {},
      source_context,
      scratch,
      &cache,
      nullptr);

  for (std::size_t i = 0; i < full_state.size(); ++i) {
    const auto a = full_state.loadCell(i);
    const auto b = active_state.loadCell(i);
    assert(std::abs(a.mass_density_comoving - b.mass_density_comoving) < 1.0e-12);
    assert(std::abs(a.momentum_density_x_comoving - b.momentum_density_x_comoving) < 1.0e-12);
    assert(std::abs(a.total_energy_density_comoving - b.total_energy_density_comoving) < 1.0e-12);
  }
}

}  // namespace

int main() {
  testLimiterLibrary();
  testPrimitiveConservedRoundTrip();
  testMusclReconstructionProducesFiniteStates();
  testComovingSourceTermSanity();
  testRiemannSymmetryRegression();
  testProfileFallbackCountersAreStepLocal();
  testActiveSetFullCoverageMatchesFullAdvance();
  return 0;
}
