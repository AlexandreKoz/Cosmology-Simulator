#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <vector>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/gravity/tree_gravity.hpp"
#include "cosmosim/hydro/hydro_core_solver.hpp"

namespace {

cosmosim::hydro::HydroPatchGeometry makePeriodic1dGeometry(std::size_t cell_count) {
  cosmosim::hydro::HydroPatchGeometry geometry;
  geometry.cell_volume_comoving = 1.0 / static_cast<double>(cell_count);
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

}  // namespace

int main() {
  constexpr std::size_t k_particle_count = 8192;
  constexpr std::size_t k_active_count = 2048;

  std::vector<double> pos_x(k_particle_count, 0.0);
  std::vector<double> pos_y(k_particle_count, 0.0);
  std::vector<double> pos_z(k_particle_count, 0.0);
  std::vector<double> mass(k_particle_count, 1.0);
  for (std::size_t i = 0; i < k_particle_count; ++i) {
    pos_x[i] = std::fmod(static_cast<double>((37U * i + 13U) % 104729U) * 7.13e-5, 1.0);
    pos_y[i] = std::fmod(static_cast<double>((73U * i + 17U) % 130363U) * 5.11e-5, 1.0);
    pos_z[i] = std::fmod(static_cast<double>((97U * i + 19U) % 156007U) * 4.73e-5, 1.0);
  }
  std::vector<std::uint32_t> active(k_active_count, 0);
  for (std::size_t i = 0; i < k_active_count; ++i) {
    active[i] = static_cast<std::uint32_t>((17U * i) % k_particle_count);
  }

  std::vector<double> accel_x(k_active_count, 0.0);
  std::vector<double> accel_y(k_active_count, 0.0);
  std::vector<double> accel_z(k_active_count, 0.0);

  cosmosim::gravity::TreeGravitySolver tree_solver;
  cosmosim::gravity::TreeGravityOptions tree_options;
  tree_options.opening_theta = 0.7;
  tree_options.max_leaf_size = 16;
  tree_options.softening.epsilon_comoving = 3.0e-4;

  cosmosim::gravity::TreeGravityProfile tree_profile;
  const auto tree_start = std::chrono::steady_clock::now();
  tree_solver.build(pos_x, pos_y, pos_z, mass, tree_options, &tree_profile);
  tree_solver.evaluateActiveSet(pos_x, pos_y, pos_z, mass, active, accel_x, accel_y, accel_z, tree_options, &tree_profile);
  const auto tree_stop = std::chrono::steady_clock::now();

  constexpr double k_gamma = 1.4;
  constexpr std::size_t k_cells = 4096;
  cosmosim::hydro::HydroConservedStateSoa conserved(k_cells);
  for (std::size_t i = 0; i < k_cells; ++i) {
    cosmosim::hydro::HydroPrimitiveState primitive;
    primitive.rho_comoving = (i % 64U < 32U) ? 1.0 : 0.5;
    primitive.vel_x_peculiar = 0.1;
    primitive.pressure_comoving = 1.0;
    conserved.storeCell(i, cosmosim::hydro::HydroCoreSolver::conservedFromPrimitive(primitive, k_gamma));
  }

  cosmosim::hydro::HydroCoreSolver hydro_solver(k_gamma);
  const auto geometry = makePeriodic1dGeometry(k_cells);
  cosmosim::hydro::HydroUpdateContext update{.dt_code = 8.0e-4, .scale_factor = 1.0, .hubble_rate_code = 0.0};
  cosmosim::hydro::HydroSourceContext source_context{.update = update};
  cosmosim::hydro::MusclHancockReconstruction reconstruction(cosmosim::hydro::HydroReconstructionPolicy{
      .limiter = cosmosim::hydro::HydroSlopeLimiter::kMonotonizedCentral,
      .dt_over_dx_code = update.dt_code * static_cast<double>(k_cells),
      .rho_floor = 1.0e-10,
      .pressure_floor = 1.0e-10,
      .enable_muscl_hancock_predictor = true});
  cosmosim::hydro::HllcRiemannSolver riemann;
  cosmosim::hydro::HydroProfileEvent hydro_profile;

  const auto hydro_setup_start = std::chrono::steady_clock::now();
  hydro_solver.advancePatch(conserved, geometry, update, reconstruction, riemann, {}, source_context, &hydro_profile);
  const auto hydro_setup_stop = std::chrono::steady_clock::now();

  hydro_profile = {};
  const auto hydro_steady_start = std::chrono::steady_clock::now();
  for (int iter = 0; iter < 8; ++iter) {
    hydro_solver.advancePatch(conserved, geometry, update, reconstruction, riemann, {}, source_context, &hydro_profile);
  }
  const auto hydro_steady_stop = std::chrono::steady_clock::now();

  const double tree_ms = std::chrono::duration<double, std::milli>(tree_stop - tree_start).count();
  const double hydro_setup_ms = std::chrono::duration<double, std::milli>(hydro_setup_stop - hydro_setup_start).count();
  const double hydro_steady_ms = std::chrono::duration<double, std::milli>(hydro_steady_stop - hydro_steady_start).count();
  const double hydro_seconds = hydro_steady_ms * 1.0e-3;
  const double hydro_face_rate_mface_s = static_cast<double>(hydro_profile.face_count) / hydro_seconds * 1.0e-6;

  std::cout << "bench_validation_ladder"
            << " build_type=" << COSMOSIM_BUILD_TYPE
            << " threads=1"
            << " features=tree_gravity+hydro_core"
            << " hardware=unspecified"
            << " tree_particles=" << k_particle_count
            << " tree_active=" << k_active_count
            << " tree_total_ms=" << tree_ms
            << " tree_visited_nodes=" << tree_profile.visited_nodes
            << " tree_accepted_nodes=" << tree_profile.accepted_nodes
            << " hydro_cells=" << k_cells
            << " hydro_setup_ms=" << hydro_setup_ms
            << " hydro_steady_ms=" << hydro_steady_ms
            << " hydro_faces=" << hydro_profile.face_count
            << " hydro_face_rate_mface_s=" << hydro_face_rate_mface_s
            << " hydro_bytes_moved=" << hydro_profile.bytes_moved
            << '\n';

  return 0;
}
