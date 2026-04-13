#include <cstddef>
#include <cstdint>
#include <vector>

#include "bench/reporting/bench_report.hpp"
#include "cosmosim/hydro/hydro_core_solver.hpp"

namespace {

[[nodiscard]] cosmosim::hydro::HydroPatchGeometry makePeriodic1dGeometry(std::size_t cell_count) {
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

void fillConserved(cosmosim::hydro::HydroConservedStateSoa& conserved, double adiabatic_index) {
  for (std::size_t i = 0; i < conserved.size(); ++i) {
    cosmosim::hydro::HydroPrimitiveState primitive;
    primitive.rho_comoving = (i % 32U < 16U) ? 1.0 : 0.5;
    primitive.vel_x_peculiar = (i % 16U) * 0.02;
    primitive.vel_y_peculiar = 0.01;
    primitive.vel_z_peculiar = 0.0;
    primitive.pressure_comoving = 1.0;
    conserved.storeCell(i, cosmosim::hydro::HydroCoreSolver::conservedFromPrimitive(primitive, adiabatic_index));
  }
}

}  // namespace

int main() {
  constexpr double k_adiabatic_index = 1.4;
  constexpr std::size_t k_cell_count = 65536;

  const auto execution = cosmosim::bench::defaultExecutionConfig(2, 16);

  cosmosim::hydro::HydroConservedStateSoa conserved(k_cell_count);
  fillConserved(conserved, k_adiabatic_index);
  const auto geometry = makePeriodic1dGeometry(k_cell_count);

  std::vector<double> gravity_x(k_cell_count, 0.0);
  std::vector<double> gravity_y(k_cell_count, 0.0);
  std::vector<double> gravity_z(k_cell_count, 0.0);
  std::vector<double> hydrogen_nh_cgs(k_cell_count, 0.1);
  std::vector<double> metallicity(k_cell_count, 0.02);
  std::vector<double> temperature_k(k_cell_count, 1.0e4);

  cosmosim::hydro::HydroUpdateContext update;
  update.dt_code = 2.5e-4;
  update.scale_factor = 1.0;

  cosmosim::hydro::HydroSourceContext source_context;
  source_context.update = update;
  source_context.gravity_accel_x_peculiar = gravity_x;
  source_context.gravity_accel_y_peculiar = gravity_y;
  source_context.gravity_accel_z_peculiar = gravity_z;
  source_context.hydrogen_number_density_cgs = hydrogen_nh_cgs;
  source_context.metallicity_mass_fraction = metallicity;
  source_context.temperature_k = temperature_k;

  cosmosim::hydro::HydroCoreSolver solver(k_adiabatic_index);
  cosmosim::hydro::MusclHancockReconstruction reconstruction(cosmosim::hydro::HydroReconstructionPolicy{
      .limiter = cosmosim::hydro::HydroSlopeLimiter::kMonotonizedCentral,
      .dt_over_dx_code = update.dt_code,
      .rho_floor = 1.0e-10,
      .pressure_floor = 1.0e-10,
      .enable_muscl_hancock_predictor = true});
  cosmosim::hydro::HllcRiemannSolver riemann_solver;

  cosmosim::hydro::HydroScratchBuffers scratch;
  cosmosim::hydro::HydroPrimitiveCacheSoa primitive_cache(k_cell_count);
  cosmosim::hydro::HydroProfileEvent profile{};

  for (std::size_t iter = 0; iter < execution.warmup_iterations; ++iter) {
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
        nullptr);
  }

  const auto begin = cosmosim::bench::BenchmarkClock::now();
  for (std::size_t iter = 0; iter < execution.measurement_iterations; ++iter) {
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
  const auto end = cosmosim::bench::BenchmarkClock::now();

  const double measurement_ms = cosmosim::bench::BenchmarkClock::millisecondsBetween(begin, end);
  const double face_updates_per_second = measurement_ms > 0.0
                                             ? static_cast<double>(profile.face_count) / (measurement_ms * 1.0e-3)
                                             : 0.0;

  cosmosim::bench::BenchmarkReporter reporter("bench_hydro_face_kernel");
  cosmosim::bench::addExecutionFields(reporter, execution);
  reporter.addField("cell_count", k_cell_count);
  reporter.addField("measurement_ms", measurement_ms);
  reporter.addField("face_count", profile.face_count);
  reporter.addField("reconstruct_ms", profile.reconstruct_ms);
  reporter.addField("riemann_ms", profile.riemann_ms);
  reporter.addField("accumulate_ms", profile.accumulate_ms);
  reporter.addField("source_ms", profile.source_ms);
  reporter.addField("total_ms", profile.total_ms);
  reporter.addField("face_updates_per_second", face_updates_per_second);
  reporter.addField("limiter_clip_count", profile.limiter_clip_count);
  reporter.addField("positivity_fallback_count", profile.positivity_fallback_count);
  reporter.addField("riemann_fallback_count", profile.riemann_fallback_count);
  cosmosim::bench::addBandwidthFields(reporter, profile.bytes_moved, measurement_ms);
  reporter.write();

  return 0;
}
