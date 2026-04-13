#include <chrono>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <vector>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/hydro/hydro_core_solver.hpp"

namespace {

cosmosim::hydro::HydroPatchGeometry makePeriodic1dGeometry(std::size_t cell_count) {
  cosmosim::hydro::HydroPatchGeometry geometry;
  geometry.cell_volume_comoving = 1.0;
  geometry.faces.reserve(cell_count);
  for (std::size_t i = 0; i < cell_count; ++i) {
    geometry.faces.push_back(cosmosim::hydro::HydroFace{
        .owner_cell = i,
        .neighbor_cell = (i + 1U) % cell_count,
        .area_comoving = 1.0,
        .normal_x = 1.0,
        .normal_y = 0.0,
        .normal_z = 0.0});
  }
  return geometry;
}

void fillInitialState(cosmosim::hydro::HydroConservedStateSoa& conserved, double gamma) {
  for (std::size_t i = 0; i < conserved.size(); ++i) {
    cosmosim::hydro::HydroPrimitiveState primitive;
    primitive.rho_comoving = (i % 32U < 16U) ? 1.0 : 0.5;
    primitive.vel_x_peculiar = (i % 8U) * 0.01;
    primitive.vel_y_peculiar = 0.0;
    primitive.vel_z_peculiar = 0.0;
    primitive.pressure_comoving = 1.0;
    conserved.storeCell(i, cosmosim::hydro::HydroCoreSolver::conservedFromPrimitive(primitive, gamma));
  }
}

}  // namespace

int main() {
  constexpr double k_gamma = 1.4;
  constexpr std::size_t k_cell_count = 131072;
  constexpr std::size_t k_iterations = 20;

  cosmosim::hydro::HydroConservedStateSoa conserved(k_cell_count);
  fillInitialState(conserved, k_gamma);

  const cosmosim::hydro::HydroPatchGeometry geometry = makePeriodic1dGeometry(k_cell_count);

  cosmosim::hydro::HydroUpdateContext update;
  update.dt_code = 2.5e-4;
  update.scale_factor = 1.0;
  update.hubble_rate_code = 0.0;

  cosmosim::hydro::HydroSourceContext source_context;
  source_context.update = update;

  cosmosim::hydro::HydroCoreSolver solver(k_gamma);
  cosmosim::hydro::MusclHancockReconstruction reconstruction(cosmosim::hydro::HydroReconstructionPolicy{
      .limiter = cosmosim::hydro::HydroSlopeLimiter::kMonotonizedCentral,
      .dt_over_dx_code = update.dt_code,
      .rho_floor = 1.0e-10,
      .pressure_floor = 1.0e-10,
      .enable_muscl_hancock_predictor = true});
  cosmosim::hydro::HllcRiemannSolver riemann_solver;
  cosmosim::hydro::HydroScratchBuffers scratch;
  cosmosim::hydro::HydroPrimitiveCacheSoa primitive_cache(k_cell_count);

  cosmosim::hydro::HydroProfileEvent profile;

  const auto setup_start = std::chrono::steady_clock::now();
  solver.advancePatchWithScratch(
      conserved,
      geometry,
      update,
      reconstruction,
      riemann_solver,
      {},
      source_context,
      scratch,
      &primitive_cache,
      &profile);
  const auto setup_stop = std::chrono::steady_clock::now();

  profile = {};
  const auto steady_start = std::chrono::steady_clock::now();
  for (std::size_t iter = 0; iter < k_iterations; ++iter) {
    solver.advancePatchWithScratch(
        conserved,
        geometry,
        update,
        reconstruction,
        riemann_solver,
        {},
        source_context,
        scratch,
        &primitive_cache,
        &profile);
  }
  const auto steady_stop = std::chrono::steady_clock::now();

  const double setup_ms = std::chrono::duration<double, std::milli>(setup_stop - setup_start).count();
  const double steady_ms = std::chrono::duration<double, std::milli>(steady_stop - steady_start).count();
  const double steady_seconds = steady_ms * 1.0e-3;

  const std::uint64_t total_faces = profile.face_count;
  const double face_throughput_mface_s = static_cast<double>(total_faces) / steady_seconds * 1.0e-6;
  const double effective_bandwidth_gb_s = static_cast<double>(profile.bytes_moved) / steady_seconds * 1.0e-9;

  std::cout << "bench_hydro_face_sweep"
            << " build_type=" << COSMOSIM_BUILD_TYPE
            << " threads=1"
            << " features=hydro_core_solver+muscl_hancock+hllc"
            << " cache=primitive+scratch_reuse"
            << " setup_ms=" << setup_ms
            << " steady_ms=" << steady_ms
            << " iterations=" << k_iterations
            << " cells=" << k_cell_count
            << " faces=" << total_faces
            << " reconstruct_ms=" << profile.reconstruct_ms
            << " riemann_ms=" << profile.riemann_ms
            << " accumulate_ms=" << profile.accumulate_ms
            << " source_ms=" << profile.source_ms
            << " total_ms=" << profile.total_ms
            << " face_throughput_mface_s=" << face_throughput_mface_s
            << " bytes_moved=" << profile.bytes_moved
            << " effective_bandwidth_gb_s=" << effective_bandwidth_gb_s
            << " limiter_clips=" << profile.limiter_clip_count
            << " reconstruction_positivity_fallbacks=" << profile.positivity_fallback_count
            << " riemann_hlle_fallbacks=" << profile.riemann_fallback_count
            << '\n';

  return 0;
}
