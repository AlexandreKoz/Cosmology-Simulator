#include <array>
#include <cassert>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "cosmosim/hydro/hydro_boundary_conditions.hpp"
#include "cosmosim/hydro/hydro_cartesian_patch.hpp"
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
      .dt_over_cell_width_code = {0.1, 0.1, 0.1},
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

  const double inverse_a2 = 1.0 / (context.update.scale_factor * context.update.scale_factor);
  const double kinetic_energy_density = 0.5 * primitive.rho_comoving *
      (primitive.vel_x_peculiar * primitive.vel_x_peculiar +
       primitive.vel_y_peculiar * primitive.vel_y_peculiar +
       primitive.vel_z_peculiar * primitive.vel_z_peculiar);
  const double expected_momentum_x =
      primitive.rho_comoving * gravity_x[0] * inverse_a2 -
      context.update.hubble_rate_code * conserved.momentum_density_x_comoving;
  const double expected_momentum_y =
      primitive.rho_comoving * gravity_y[0] * inverse_a2 -
      context.update.hubble_rate_code * conserved.momentum_density_y_comoving;
  const double expected_momentum_z =
      primitive.rho_comoving * gravity_z[0] * inverse_a2 -
      context.update.hubble_rate_code * conserved.momentum_density_z_comoving;
  const double expected_gravity_work = primitive.rho_comoving * inverse_a2 *
      (primitive.vel_x_peculiar * gravity_x[0] +
       primitive.vel_y_peculiar * gravity_y[0] +
       primitive.vel_z_peculiar * gravity_z[0]);
  const double expected_energy = expected_gravity_work - context.update.hubble_rate_code *
      (2.0 * kinetic_energy_density + 3.0 * primitive.pressure_comoving);
  assert(std::abs(source_state.mass_density_comoving) < k_tol);
  assert(std::abs(source_state.momentum_density_x_comoving - expected_momentum_x) < k_tol);
  assert(std::abs(source_state.momentum_density_y_comoving - expected_momentum_y) < k_tol);
  assert(std::abs(source_state.momentum_density_z_comoving - expected_momentum_z) < k_tol);
  assert(std::abs(source_state.total_energy_density_comoving - expected_energy) < k_tol);

  const std::vector<double> zero_gravity{0.0};
  context.gravity_accel_x_peculiar = zero_gravity;
  context.gravity_accel_y_peculiar = zero_gravity;
  context.gravity_accel_z_peculiar = zero_gravity;
  const cosmosim::hydro::HydroConservedState homogeneous_source =
      source.sourceForCell(0, conserved, primitive, context);
  assert(std::abs(
      homogeneous_source.total_energy_density_comoving +
      context.update.hubble_rate_code *
          (2.0 * kinetic_energy_density + 3.0 * primitive.pressure_comoving)) < k_tol);
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
        .owner_minus_cell = (i + conserved.size() - 1U) % conserved.size(),
        .neighbor_plus_cell = (i + 2U) % conserved.size(),
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
        .owner_minus_cell = (i + full_state.size() - 1U) % full_state.size(),
        .neighbor_plus_cell = (i + 2U) % full_state.size(),
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


void testActivePeriodicGhostFluxUpdatesWrappedRealCell() {
  constexpr double gamma = 1.4;
  cosmosim::hydro::HydroPatchGeometry geometry = cosmosim::hydro::makeCartesianPatchGeometry(
      cosmosim::hydro::HydroCartesianPatchSpec{
          .nx = 4,
          .ny = 1,
          .nz = 1,
          .origin_x_comoving = 0.0,
          .origin_y_comoving = 0.0,
          .origin_z_comoving = 0.0,
          .cell_width_x_comoving = 0.25,
          .cell_width_y_comoving = 1.0,
          .cell_width_z_comoving = 1.0});
  cosmosim::hydro::appendCartesianBoundaryGhostFaces(geometry, cosmosim::hydro::HydroBoundaryKind::kPeriodic);

  cosmosim::hydro::HydroConservedStateSoa conserved(geometry.totalCellStorageCount());
  for (std::size_t cell = 0; cell < geometry.cellCount(); ++cell) {
    cosmosim::hydro::HydroPrimitiveState primitive;
    primitive.rho_comoving = 1.0 + 0.2 * static_cast<double>(cell);
    primitive.vel_x_peculiar = 0.15 - 0.04 * static_cast<double>(cell);
    primitive.vel_y_peculiar = 0.01 * static_cast<double>(cell);
    primitive.vel_z_peculiar = -0.005 * static_cast<double>(cell);
    primitive.pressure_comoving = 1.0 + 0.1 * static_cast<double>(cell);
    conserved.storeCell(cell, cosmosim::hydro::HydroCoreSolver::conservedFromPrimitive(primitive, gamma));
  }
  cosmosim::hydro::fillHydroBoundaryGhostCells(conserved, geometry, gamma);

  const std::size_t lower_x_face = [&geometry]() {
    for (std::size_t face_index = 0; face_index < geometry.faces.size(); ++face_index) {
      const cosmosim::hydro::HydroFace& face = geometry.faces[face_index];
      if (face.ghost_cell_slot == cosmosim::hydro::k_invalid_ghost_cell_slot) {
        continue;
      }
      const cosmosim::hydro::HydroGhostCell& ghost = geometry.ghost_cells[face.ghost_cell_slot];
      if (ghost.boundary_kind == cosmosim::hydro::HydroBoundaryKind::kPeriodic &&
          ghost.axis == cosmosim::hydro::HydroFaceAxis::kX &&
          ghost.side == cosmosim::hydro::HydroFaceSide::kLower &&
          ghost.owner_real_cell == 0U &&
          ghost.source_real_cell == 3U) {
        return face_index;
      }
    }
    throw std::runtime_error("missing lower-x periodic ghost face");
  }();

  const std::vector<std::size_t> measured_cells{0U, 3U};
  const cosmosim::hydro::HydroConservationTotals before =
      cosmosim::hydro::HydroCoreSolver::conservationTotals(conserved, geometry, measured_cells);
  const cosmosim::hydro::HydroConservedState wrapped_before = conserved.loadCell(3U);

  cosmosim::hydro::HydroUpdateContext update{.dt_code = 1.0e-4, .scale_factor = 1.0, .hubble_rate_code = 0.0};
  cosmosim::hydro::HydroSourceContext source_context{.update = update};
  cosmosim::hydro::HydroCoreSolver solver(gamma);
  cosmosim::hydro::MusclHancockReconstruction reconstruction(cosmosim::hydro::HydroReconstructionPolicy{
      .limiter = cosmosim::hydro::HydroSlopeLimiter::kMonotonizedCentral,
      .dt_over_dx_code = update.dt_code / geometry.cell_width_x_comoving,
      .dt_over_cell_width_code = {
          update.dt_code / geometry.cell_width_x_comoving,
          update.dt_code / geometry.cell_width_y_comoving,
          update.dt_code / geometry.cell_width_z_comoving},
      .rho_floor = 1.0e-12,
      .pressure_floor = 1.0e-12,
      .enable_muscl_hancock_predictor = true,
      .adiabatic_index = gamma});
  cosmosim::hydro::HllcRiemannSolver riemann;
  cosmosim::hydro::HydroScratchBuffers scratch;
  cosmosim::hydro::HydroPrimitiveCacheSoa cache(conserved.size());
  cosmosim::hydro::HydroProfileEvent profile;
  const std::vector<std::size_t> active_cells{0U};
  const std::vector<std::size_t> active_faces{lower_x_face};

  solver.advancePatchActiveSetWithScratch(
      conserved,
      geometry,
      cosmosim::hydro::HydroActiveSetView{.active_cells = active_cells, .active_faces = active_faces},
      update,
      reconstruction,
      riemann,
      {},
      source_context,
      scratch,
      &cache,
      &profile);

  const cosmosim::hydro::HydroConservedState wrapped_after = conserved.loadCell(3U);
  assert(std::abs(wrapped_after.mass_density_comoving - wrapped_before.mass_density_comoving) > 1.0e-12);
  const cosmosim::hydro::HydroConservationTotals after =
      cosmosim::hydro::HydroCoreSolver::conservationTotals(conserved, geometry, measured_cells);
  assert(std::abs(after.mass - before.mass) < 1.0e-12);
  assert(std::abs(after.momentum_x - before.momentum_x) < 1.0e-12);
  assert(std::abs(after.momentum_y - before.momentum_y) < 1.0e-12);
  assert(std::abs(after.momentum_z - before.momentum_z) < 1.0e-12);
  assert(std::abs(after.total_energy - before.total_energy) < 1.0e-12);
  assert(profile.conservation.cell_count == 2U);
  assert(std::abs(profile.conservation.residual.mass) < 1.0e-12);
}

void testConservationReportSeparatesSourceTerms() {
  constexpr double gamma = 1.4;
  cosmosim::hydro::HydroConservedStateSoa conserved(2);
  for (std::size_t i = 0; i < conserved.size(); ++i) {
    cosmosim::hydro::HydroPrimitiveState primitive;
    primitive.rho_comoving = 1.0;
    primitive.vel_x_peculiar = 0.1;
    primitive.pressure_comoving = 1.0;
    conserved.storeCell(i, cosmosim::hydro::HydroCoreSolver::conservedFromPrimitive(primitive, gamma));
  }

  cosmosim::hydro::HydroPatchGeometry geometry;
  geometry.cell_volume_comoving = 2.0;
  geometry.faces.push_back(cosmosim::hydro::HydroFace{
      .owner_cell = 0,
      .neighbor_cell = 1,
      .area_comoving = 1.0,
      .normal_x = 1.0});

  cosmosim::hydro::HydroUpdateContext update{.dt_code = 1.0e-3, .scale_factor = 1.0, .hubble_rate_code = 0.0};
  std::vector<double> gravity_x{0.2, 0.2};
  std::vector<double> gravity_y{0.0, 0.0};
  std::vector<double> gravity_z{0.0, 0.0};
  cosmosim::hydro::HydroSourceContext source_context{
      .update = update,
      .gravity_accel_x_peculiar = gravity_x,
      .gravity_accel_y_peculiar = gravity_y,
      .gravity_accel_z_peculiar = gravity_z};

  cosmosim::hydro::HydroCoreSolver solver(gamma);
  cosmosim::hydro::MusclHancockReconstruction reconstruction;
  cosmosim::hydro::HllcRiemannSolver riemann;
  cosmosim::hydro::ComovingGravityExpansionSource source;
  std::array<const cosmosim::hydro::HydroSourceTerm*, 1> sources{&source};
  cosmosim::hydro::HydroProfileEvent profile;

  solver.advancePatch(conserved, geometry, update, reconstruction, riemann, sources, source_context, &profile);

  assert(profile.conservation.cell_count == 2U);
  assert(std::abs(profile.conservation.source_delta.momentum_x) > 0.0);
  assert(std::abs(profile.conservation.residual.mass) < 1.0e-12);
  assert(std::abs(profile.conservation.residual.momentum_x) < 1.0e-12);
  assert(std::abs(profile.conservation.residual.total_energy) < 1.0e-12);
}

void testInternalEnergyFloorIsReported() {
  constexpr double gamma = 1.4;
  cosmosim::hydro::HydroConservedStateSoa conserved(1);
  conserved.storeCell(0, cosmosim::hydro::HydroConservedState{
      .mass_density_comoving = 1.0,
      .momentum_density_x_comoving = 10.0,
      .momentum_density_y_comoving = 0.0,
      .momentum_density_z_comoving = 0.0,
      .total_energy_density_comoving = 1.0});

  cosmosim::hydro::HydroPatchGeometry geometry;
  geometry.cell_volume_comoving = 1.0;
  cosmosim::hydro::HydroUpdateContext update{.dt_code = 1.0e-3, .scale_factor = 1.0, .hubble_rate_code = 0.0};
  cosmosim::hydro::HydroSourceContext source_context{.update = update};
  cosmosim::hydro::HydroCoreSolver solver(gamma);
  cosmosim::hydro::MusclHancockReconstruction reconstruction;
  cosmosim::hydro::HllcRiemannSolver riemann;
  cosmosim::hydro::HydroProfileEvent profile;

  solver.advancePatch(conserved, geometry, update, reconstruction, riemann, {}, source_context, &profile);

  assert(profile.conservation.internal_energy_floor_count == 1U);
  assert(profile.conservation.floor_delta.total_energy > 0.0);
  const auto primitive = cosmosim::hydro::HydroCoreSolver::primitiveFromConserved(conserved.loadCell(0), gamma);
  assert(primitive.pressure_comoving > 0.0);
}

void testCartesianPatchGeometry2x2x2() {
  const cosmosim::hydro::HydroPatchGeometry geometry = cosmosim::hydro::makeCartesianPatchGeometry(
      cosmosim::hydro::HydroCartesianPatchSpec{
          .nx = 2,
          .ny = 2,
          .nz = 2,
          .origin_x_comoving = 1.0,
          .origin_y_comoving = 2.0,
          .origin_z_comoving = 3.0,
          .cell_width_x_comoving = 0.5,
          .cell_width_y_comoving = 0.25,
          .cell_width_z_comoving = 0.125});

  assert(geometry.cellCount() == 8U);
  assert(geometry.faces.size() == 12U);
  assert(std::abs(geometry.cell_volume_comoving - 0.015625) < k_tol);
  assert(geometry.linearCellIndex(1, 0, 0) == 1U);
  assert(geometry.linearCellIndex(0, 1, 0) == 2U);
  assert(geometry.linearCellIndex(0, 0, 1) == 4U);
  const auto ijk = geometry.cellIjk(6U);
  assert(ijk[0] == 0U && ijk[1] == 1U && ijk[2] == 1U);
  assert(geometry.neighborCell(0U, 1, 0, 0) == 1U);
  assert(geometry.neighborCell(0U, 0, 1, 0) == 2U);
  assert(geometry.neighborCell(0U, 0, 0, 1) == 4U);
  assert(geometry.neighborCell(0U, -1, 0, 0) == cosmosim::hydro::k_invalid_cell_index);

  std::size_t x_faces = 0;
  std::size_t y_faces = 0;
  std::size_t z_faces = 0;
  for (const cosmosim::hydro::HydroFace& face : geometry.faces) {
    if (face.axis == cosmosim::hydro::HydroFaceAxis::kX) {
      ++x_faces;
      assert(std::abs(face.area_comoving - 0.03125) < k_tol);
      assert(face.normal_x == 1.0 && face.normal_y == 0.0 && face.normal_z == 0.0);
    } else if (face.axis == cosmosim::hydro::HydroFaceAxis::kY) {
      ++y_faces;
      assert(std::abs(face.area_comoving - 0.0625) < k_tol);
      assert(face.normal_x == 0.0 && face.normal_y == 1.0 && face.normal_z == 0.0);
    } else {
      ++z_faces;
      assert(std::abs(face.area_comoving - 0.125) < k_tol);
      assert(face.normal_x == 0.0 && face.normal_y == 0.0 && face.normal_z == 1.0);
    }
    assert(face.neighbor_cell != cosmosim::hydro::k_invalid_cell_index);
  }
  assert(x_faces == 4U);
  assert(y_faces == 4U);
  assert(z_faces == 4U);

  const std::array<std::size_t, 3> factors = cosmosim::hydro::chooseNearCubicCartesianFactors(24U);
  assert(factors[0] == 4U && factors[1] == 3U && factors[2] == 2U);
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
  testActivePeriodicGhostFluxUpdatesWrappedRealCell();
  testConservationReportSeparatesSourceTerms();
  testInternalEnergyFloorIsReported();
  testCartesianPatchGeometry2x2x2();
  return 0;
}
