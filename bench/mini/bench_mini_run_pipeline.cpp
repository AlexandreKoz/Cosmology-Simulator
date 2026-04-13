#include <cstddef>
#include <cstdint>
#include <sstream>
#include <string>
#include <vector>

#include "bench/reporting/bench_report.hpp"
#include "cosmosim/amr/amr_framework.hpp"
#include "cosmosim/core/config.hpp"
#include "cosmosim/core/simulation_mode.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/gravity/tree_gravity.hpp"
#include "cosmosim/hydro/hydro_core_solver.hpp"

namespace {

[[nodiscard]] cosmosim::core::FrozenConfig makeMiniConfig() {
  std::stringstream stream;
  stream << "schema_version = 1\n\n";
  stream << "[mode]\n";
  stream << "mode = zoom_in\n";
  stream << "ic_file = ics.hdf5\n";
  stream << "zoom_high_res_region = false\n\n";
  stream << "[numerics]\n";
  stream << "max_global_steps = 4\n";
  return cosmosim::core::loadFrozenConfigFromString(stream.str(), "bench_mini_run_pipeline");
}

void seedState(cosmosim::core::SimulationState& state, std::size_t particle_count) {
  state.resizeParticles(particle_count);
  for (std::size_t i = 0; i < particle_count; ++i) {
    state.particles.position_x_comoving[i] = static_cast<double>((i * 13U) % 1024U) / 1024.0;
    state.particles.position_y_comoving[i] = static_cast<double>((i * 17U) % 1024U) / 1024.0;
    state.particles.position_z_comoving[i] = static_cast<double>((i * 19U) % 1024U) / 1024.0;
    state.particles.mass_code[i] = 1.0;
  }
}

}  // namespace

int main() {
  constexpr std::size_t k_particle_count = 32768;
  constexpr std::size_t k_hydro_cells = 4096;

  const auto execution = cosmosim::bench::defaultExecutionConfig(1, 4);

  const cosmosim::core::FrozenConfig frozen_config = makeMiniConfig();
  const cosmosim::core::SimulationConfig& config = frozen_config.config;
  const cosmosim::core::ModePolicy mode_policy = cosmosim::core::buildModePolicy(config.mode);

  cosmosim::core::SimulationState state;
  seedState(state, k_particle_count);

  std::vector<std::uint32_t> active(k_particle_count);
  for (std::size_t i = 0; i < k_particle_count; ++i) {
    active[i] = static_cast<std::uint32_t>(i);
  }

  std::vector<double> accel_x(k_particle_count, 0.0);
  std::vector<double> accel_y(k_particle_count, 0.0);
  std::vector<double> accel_z(k_particle_count, 0.0);

  cosmosim::gravity::TreeGravitySolver gravity_solver;
  cosmosim::gravity::TreeGravityOptions gravity_options;
  gravity_options.opening_theta = 0.8;

  cosmosim::hydro::HydroConservedStateSoa hydro_state(k_hydro_cells);
  cosmosim::hydro::HydroPatchGeometry geometry;
  geometry.cell_volume_comoving = 1.0;
  geometry.faces.reserve(k_hydro_cells);
  for (std::size_t i = 0; i < k_hydro_cells; ++i) {
    geometry.faces.push_back(cosmosim::hydro::HydroFace{.owner_cell = i,
                                                        .neighbor_cell = (i + 1U) % k_hydro_cells,
                                                        .area_comoving = 1.0,
                                                        .normal_x = 1.0,
                                                        .normal_y = 0.0,
                                                        .normal_z = 0.0});
    cosmosim::hydro::HydroPrimitiveState primitive;
    primitive.rho_comoving = 1.0;
    primitive.pressure_comoving = 1.0;
    hydro_state.storeCell(i, cosmosim::hydro::HydroCoreSolver::conservedFromPrimitive(primitive, 1.4));
  }

  cosmosim::hydro::HydroCoreSolver hydro_solver(1.4);
  cosmosim::hydro::MusclHancockReconstruction reconstruction(cosmosim::hydro::HydroReconstructionPolicy{});
  cosmosim::hydro::HllcRiemannSolver riemann_solver;
  cosmosim::hydro::HydroUpdateContext update{.dt_code = 1.0e-3, .scale_factor = 1.0, .hubble_rate_code = 0.0};
  std::vector<double> zero(k_hydro_cells, 0.0);
  cosmosim::hydro::HydroSourceContext source_context{.update = update,
                                                     .gravity_accel_x_peculiar = zero,
                                                     .gravity_accel_y_peculiar = zero,
                                                     .gravity_accel_z_peculiar = zero,
                                                     .hydrogen_number_density_cgs = zero,
                                                     .metallicity_mass_fraction = zero,
                                                     .temperature_k = zero,
                                                     .redshift = 0.0};

  cosmosim::amr::PatchHierarchy hierarchy;
  cosmosim::amr::PatchDescriptor root;
  root.cell_dims = {8, 8, 8};
  const std::uint64_t root_id = hierarchy.createRootPatch(root);
  (void)root_id;

  for (std::size_t iter = 0; iter < execution.warmup_iterations; ++iter) {
    gravity_solver.build(
        state.particles.position_x_comoving,
        state.particles.position_y_comoving,
        state.particles.position_z_comoving,
        state.particles.mass_code,
        gravity_options,
        nullptr);
    gravity_solver.evaluateActiveSet(
        state.particles.position_x_comoving,
        state.particles.position_y_comoving,
        state.particles.position_z_comoving,
        state.particles.mass_code,
        active,
        accel_x,
        accel_y,
        accel_z,
        gravity_options,
        nullptr);

    hydro_solver.advancePatch(hydro_state, geometry, update, reconstruction, riemann_solver, {}, source_context, nullptr);
  }

  cosmosim::gravity::TreeGravityProfile gravity_profile{};
  cosmosim::hydro::HydroProfileEvent hydro_profile{};

  const auto begin = cosmosim::bench::BenchmarkClock::now();
  for (std::size_t iter = 0; iter < execution.measurement_iterations; ++iter) {
    gravity_solver.build(
        state.particles.position_x_comoving,
        state.particles.position_y_comoving,
        state.particles.position_z_comoving,
        state.particles.mass_code,
        gravity_options,
        &gravity_profile);
    gravity_solver.evaluateActiveSet(
        state.particles.position_x_comoving,
        state.particles.position_y_comoving,
        state.particles.position_z_comoving,
        state.particles.mass_code,
        active,
        accel_x,
        accel_y,
        accel_z,
        gravity_options,
        &gravity_profile);

    hydro_solver.advancePatch(hydro_state, geometry, update, reconstruction, riemann_solver, {}, source_context, &hydro_profile);
  }
  const auto end = cosmosim::bench::BenchmarkClock::now();

  const double measurement_ms = cosmosim::bench::BenchmarkClock::millisecondsBetween(begin, end);

  cosmosim::bench::BenchmarkReporter reporter("bench_mini_run_pipeline");
  cosmosim::bench::addExecutionFields(reporter, execution);
  reporter.addField("mode", cosmosim::core::modeToString(config.mode.mode));
  reporter.addField("hydro_boundary", cosmosim::core::boundaryConditionToString(mode_policy.hydro_boundary));
  reporter.addField("gravity_boundary", cosmosim::core::gravityBoundaryModelToString(mode_policy.gravity_boundary));
  reporter.addField("particle_count", k_particle_count);
  reporter.addField("hydro_cells", k_hydro_cells);
  reporter.addField("measurement_ms", measurement_ms);
  reporter.addField("gravity_build_ms", gravity_profile.build_ms);
  reporter.addField("gravity_traversal_ms", gravity_profile.traversal_ms);
  reporter.addField("hydro_total_ms", hydro_profile.total_ms);
  reporter.addField("hydro_face_count", hydro_profile.face_count);
  reporter.write();

  return 0;
}
