#include <algorithm>
#include <cassert>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "cosmosim/gravity/pm_solver.hpp"
#include "cosmosim/hydro/hydro_core_solver.hpp"
#include "cosmosim/core/config.hpp"
#include "cosmosim/physics/cooling_heating.hpp"
#include "validation_tolerance.hpp"

namespace {

constexpr double k_pi = 3.141592653589793238462643383279502884;

void requireOrThrow(bool condition, const std::string& message) {
  if (!condition) {
    throw std::runtime_error(message);
  }
}

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

void fillSodLikeInitialState(cosmosim::hydro::HydroConservedStateSoa& conserved, double gamma) {
  for (std::size_t i = 0; i < conserved.size(); ++i) {
    cosmosim::hydro::HydroPrimitiveState primitive;
    if (i < conserved.size() / 2U) {
      primitive.rho_comoving = 1.0;
      primitive.pressure_comoving = 1.0;
    } else {
      primitive.rho_comoving = 0.125;
      primitive.pressure_comoving = 0.1;
    }
    conserved.storeCell(i, cosmosim::hydro::HydroCoreSolver::conservedFromPrimitive(primitive, gamma));
  }
}

void testPmSingleMode(const cosmosim::validation::ValidationToleranceTable& tolerances) {
  const cosmosim::gravity::PmGridShape shape{32, 8, 8};
  cosmosim::gravity::PmGridStorage grid(shape);
  cosmosim::gravity::PmSolver solver(shape);

  cosmosim::gravity::PmSolveOptions options;
  options.box_size_mpc_comoving = 1.0;
  options.scale_factor = 1.0;
  options.gravitational_constant_code = 1.0;

  const std::size_t particle_count = 32;
  std::vector<double> pos_x(particle_count);
  std::vector<double> pos_y(particle_count, 0.25);
  std::vector<double> pos_z(particle_count, 0.75);
  std::vector<double> mass(particle_count, 1.0);
  for (std::size_t i = 0; i < particle_count; ++i) {
    pos_x[i] = (static_cast<double>(i) + 0.5) / static_cast<double>(particle_count) * options.box_size_mpc_comoving;
    mass[i] = 1.0 + 0.1 * std::sin(2.0 * k_pi * pos_x[i] / options.box_size_mpc_comoving);
  }

  std::vector<double> accel_x(particle_count, 0.0);
  std::vector<double> accel_y(particle_count, 0.0);
  std::vector<double> accel_z(particle_count, 0.0);
  std::vector<double> phi(particle_count, 0.0);

  solver.assignDensity(grid, pos_x, pos_y, pos_z, mass, options, nullptr);
  solver.solvePoissonPeriodic(grid, options, nullptr);
  solver.interpolateForces(grid, pos_x, pos_y, pos_z, accel_x, accel_y, accel_z, options, nullptr);
  solver.interpolatePotential(grid, pos_x, pos_y, pos_z, phi, options, nullptr);

  const double kx = 2.0 * k_pi / options.box_size_mpc_comoving;
  const double expected_amp = 4.0 * k_pi * options.gravitational_constant_code *
      options.scale_factor * options.scale_factor * 0.1 / kx;
  const double expected_phi_amp = -4.0 * k_pi * options.gravitational_constant_code *
      options.scale_factor * options.scale_factor * 0.1 / (kx * kx);

  double corr = 0.0;
  double norm_expected = 0.0;
  double norm_got = 0.0;
  double corr_phi = 0.0;
  double norm_phi_expected = 0.0;
  double norm_phi_got = 0.0;
  for (std::size_t i = 0; i < particle_count; ++i) {
    const double expected = expected_amp * std::cos(kx * pos_x[i]);
    const double expected_phi = expected_phi_amp * std::sin(kx * pos_x[i]);
    corr += expected * accel_x[i];
    norm_expected += expected * expected;
    norm_got += accel_x[i] * accel_x[i];
    corr_phi += expected_phi * phi[i];
    norm_phi_expected += expected_phi * expected_phi;
    norm_phi_got += phi[i] * phi[i];
  }

  const double cosine_similarity = corr / std::sqrt(std::max(norm_expected * norm_got, 1.0e-20));
  const double cosine_similarity_phi = corr_phi / std::sqrt(std::max(norm_phi_expected * norm_phi_got, 1.0e-20));
  requireOrThrow(
      std::abs(cosine_similarity) >= tolerances.require("gravity_pm_single_mode.min_cosine_similarity"),
      "gravity_pm_single_mode failed: cosine similarity below tolerance");
  requireOrThrow(
      std::abs(cosine_similarity_phi) >= tolerances.require("gravity_pm_single_mode.min_cosine_similarity"),
      "gravity_pm_single_mode failed: potential cosine similarity below tolerance");
}

void testHydroSodMassConservation(const cosmosim::validation::ValidationToleranceTable& tolerances) {
  constexpr std::size_t k_cell_count = 64;
  constexpr std::size_t k_step_count = 24;
  constexpr double k_gamma = 1.4;

  cosmosim::hydro::HydroConservedStateSoa conserved(k_cell_count);
  fillSodLikeInitialState(conserved, k_gamma);

  const auto geometry = makePeriodic1dGeometry(k_cell_count);

  cosmosim::hydro::HydroUpdateContext update{};
  update.dt_code = 6.0e-4;
  update.scale_factor = 1.0;

  cosmosim::hydro::HydroSourceContext source_context{};
  source_context.update = update;

  cosmosim::hydro::HydroCoreSolver solver(k_gamma);
  cosmosim::hydro::MusclHancockReconstruction reconstruction(cosmosim::hydro::HydroReconstructionPolicy{
      .limiter = cosmosim::hydro::HydroSlopeLimiter::kMonotonizedCentral,
      .dt_over_dx_code = update.dt_code * static_cast<double>(k_cell_count),
      .rho_floor = 1.0e-10,
      .pressure_floor = 1.0e-10,
      .enable_muscl_hancock_predictor = true});
  cosmosim::hydro::HllcRiemannSolver riemann;

  double initial_mass = 0.0;
  for (double rho : conserved.massDensityComoving()) {
    initial_mass += rho * geometry.cell_volume_comoving;
  }

  for (std::size_t step = 0; step < k_step_count; ++step) {
    solver.advancePatch(conserved, geometry, update, reconstruction, riemann, {}, source_context, nullptr);
  }

  double final_mass = 0.0;
  for (double rho : conserved.massDensityComoving()) {
    final_mass += rho * geometry.cell_volume_comoving;
  }

  const double mass_error = std::abs(final_mass - initial_mass);
  requireOrThrow(
      mass_error <= tolerances.require("hydro_sod_like.mass_abs_error"),
      "hydro_sod_like failed: mass conservation error above tolerance");
}

void testCoolingEnergyMonotonicity(const cosmosim::validation::ValidationToleranceTable& tolerances) {
  constexpr double k_gamma = 5.0 / 3.0;

  cosmosim::hydro::HydroPrimitiveState primitive{};
  primitive.rho_comoving = 1.0;
  primitive.pressure_comoving = 1.0;
  auto conserved = cosmosim::hydro::HydroCoreSolver::conservedFromPrimitive(primitive, k_gamma);

  cosmosim::physics::CoolingModelConfig model{};
  model.temperature_floor_k = 100.0;
  model.max_subcycles = 32;
  model.uv_background_model = cosmosim::core::UvBackgroundModel::kNone;

  cosmosim::physics::CoolingRateProvider provider(model);
  cosmosim::physics::CoolingSourceIntegrator integrator(1.0e-8);
  cosmosim::physics::CoolingHeatingSource cooling(provider, integrator);

  std::vector<double> nh(1, 5.0);
  std::vector<double> metal(1, 0.02);
  std::vector<double> temp(1, 2.0e6);
  cosmosim::hydro::HydroSourceContext source_context{};
  source_context.update.dt_code = 0.005;
  source_context.update.scale_factor = 1.0;
  source_context.hydrogen_number_density_cgs = nh;
  source_context.metallicity_mass_fraction = metal;
  source_context.temperature_k = temp;

  const double initial_total_energy = conserved.total_energy_density_comoving;
  for (int step = 0; step < 16; ++step) {
    const auto source = cooling.sourceForCell(0, conserved, primitive, source_context);
    conserved.total_energy_density_comoving += source.total_energy_density_comoving * source_context.update.dt_code;
  }

  const double allowed_increase = tolerances.require("cooling_box.total_energy_nonincrease_eps");
  requireOrThrow(
      conserved.total_energy_density_comoving <= initial_total_energy + allowed_increase,
      "cooling_box failed: thermal energy increased beyond allowed epsilon");
}

}  // namespace

int main() {
  const auto tolerances = cosmosim::validation::ValidationToleranceTable::loadFromFile(
      std::string(COSMOSIM_SOURCE_DIR) + "/validation/reference/validation_tolerances_v1.txt");

  testPmSingleMode(tolerances);
  testHydroSodMassConservation(tolerances);
  testCoolingEnergyMonotonicity(tolerances);
  return 0;
}
