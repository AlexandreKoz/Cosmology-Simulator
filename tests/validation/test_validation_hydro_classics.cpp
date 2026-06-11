#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <span>
#include <string>
#include <vector>

#include "cosmosim/hydro/hydro_core_solver.hpp"
#include "hydro_validation_helpers.hpp"

namespace {

namespace hv = cosmosim::tests::hydro_validation;

struct FixedGravitySource final : public cosmosim::hydro::HydroSourceTerm {
  std::span<const double> accel_x_code;
  std::span<const double> accel_y_code;
  std::span<const double> accel_z_code;

  [[nodiscard]] cosmosim::hydro::HydroConservedState sourceForCell(
      std::size_t cell_index,
      const cosmosim::hydro::HydroConservedState& conserved,
      const cosmosim::hydro::HydroPrimitiveState& primitive,
      const cosmosim::hydro::HydroSourceContext& context) const override {
    const double gx = cell_index < accel_x_code.size() ? accel_x_code[cell_index] : 0.0;
    const double gy = cell_index < accel_y_code.size() ? accel_y_code[cell_index] : 0.0;
    const double gz = cell_index < accel_z_code.size() ? accel_z_code[cell_index] : 0.0;

    cosmosim::hydro::HydroConservedState source{};
    source.momentum_density_x_comoving = primitive.rho_comoving * gx;
    source.momentum_density_y_comoving = primitive.rho_comoving * gy;
    source.momentum_density_z_comoving = primitive.rho_comoving * gz;
    source.total_energy_density_comoving = primitive.rho_comoving *
        (primitive.vel_x_peculiar * gx + primitive.vel_y_peculiar * gy + primitive.vel_z_peculiar * gz);
    (void)conserved;
    (void)context;
    return source;
  }
};

[[nodiscard]] double relativeDifference(double lhs, double rhs) {
  return std::abs(lhs - rhs) / std::max({std::abs(lhs), std::abs(rhs), 1.0e-14});
}

[[nodiscard]] double radiusFromCenter(const hv::CellCenter& center) {
  const double dx = center.x_comoving - 0.5;
  const double dy = center.y_comoving - 0.5;
  const double dz = center.z_comoving - 0.5;
  return std::sqrt(dx * dx + dy * dy + dz * dz);
}

void testSedovBlastGuard() {
  constexpr double gamma = 5.0 / 3.0;
  cosmosim::hydro::HydroPatchGeometry geometry = hv::makePeriodicCartesianPatch(8U, 8U, 8U);
  cosmosim::hydro::HydroConservedStateSoa conserved(geometry.cellCount());

  for (std::size_t cell = 0; cell < geometry.cellCount(); ++cell) {
    const hv::CellCenter center = hv::cellCenter(geometry, cell);
    const double radius = radiusFromCenter(center);
    cosmosim::hydro::HydroPrimitiveState primitive{};
    primitive.rho_comoving = 1.0;
    primitive.pressure_comoving = radius < 0.16 ? 8.0 : 1.0e-3;
    hv::storePrimitive(conserved, cell, primitive, gamma);
  }

  const auto initial = hv::conservedTotals(conserved, geometry);
  hv::advanceHydroSteps(conserved, geometry, hv::HydroStepConfig{
      .gamma = gamma,
      .dt_code = 6.0e-5,
      .step_count = 24,
      .rho_floor = 1.0e-9,
      .pressure_floor = 1.0e-9});

  hv::requireFinitePositiveState(conserved, gamma, "Sedov");
  const auto final = hv::conservedTotals(conserved, geometry);
  hv::requireOrThrow(
      relativeDifference(final.total_energy, initial.total_energy) < 2.0e-10,
      "Sedov: total energy drift in closed periodic patch");

  const std::array<double, 5> edges{0.0, 0.14, 0.24, 0.36, 0.70};
  const std::vector<hv::RadialBin> bins = hv::radialBins(conserved, geometry, gamma, edges);
  hv::requireOrThrow(bins[0].count > 0U && bins[1].count > 0U, "Sedov: missing radial diagnostics");
  hv::requireOrThrow(
      bins[1].mean_radial_velocity > 1.0e-3,
      "Sedov: blast did not produce outward radial motion in the shocked shell");
  hv::requireOrThrow(
      bins[0].mean_pressure > bins[3].mean_pressure,
      "Sedov: central pressure is not above far-field pressure");
  hv::requireOrThrow(
      bins[1].density_variance / (bins[1].mean_density * bins[1].mean_density) < 0.5,
      "Sedov: radial shell density variation is too large for CI symmetry guard");
}

void testNohConvergingInflowGuard() {
  constexpr double gamma = 5.0 / 3.0;
  cosmosim::hydro::HydroPatchGeometry geometry = hv::makePeriodicCartesianPatch(8U, 8U, 8U);
  cosmosim::hydro::HydroConservedStateSoa conserved(geometry.cellCount());

  for (std::size_t cell = 0; cell < geometry.cellCount(); ++cell) {
    const hv::CellCenter center = hv::cellCenter(geometry, cell);
    const double dx = center.x_comoving - 0.5;
    const double dy = center.y_comoving - 0.5;
    const double dz = center.z_comoving - 0.5;
    const double radius = std::sqrt(dx * dx + dy * dy + dz * dz);
    const double inv_radius = radius > 1.0e-14 ? 1.0 / radius : 0.0;
    cosmosim::hydro::HydroPrimitiveState primitive{};
    primitive.rho_comoving = 1.0;
    primitive.vel_x_peculiar = -0.35 * dx * inv_radius;
    primitive.vel_y_peculiar = -0.35 * dy * inv_radius;
    primitive.vel_z_peculiar = -0.35 * dz * inv_radius;
    primitive.pressure_comoving = 1.0e-3;
    hv::storePrimitive(conserved, cell, primitive, gamma);
  }

  hv::advanceHydroSteps(conserved, geometry, hv::HydroStepConfig{
      .gamma = gamma,
      .dt_code = 4.0e-5,
      .step_count = 30,
      .rho_floor = 1.0e-9,
      .pressure_floor = 1.0e-9,
      .limiter = cosmosim::hydro::HydroSlopeLimiter::kMinmod});

  hv::requireFinitePositiveState(conserved, gamma, "Noh");
  const std::array<double, 5> edges{0.0, 0.16, 0.28, 0.42, 0.70};
  const std::vector<hv::RadialBin> bins = hv::radialBins(conserved, geometry, gamma, edges);
  hv::requireOrThrow(
      bins[0].mean_density > 1.01,
      "Noh: converging inflow did not compress the central region");
  hv::requireOrThrow(
      bins[0].mean_density > bins[3].mean_density,
      "Noh: central density is not above the outer inflow density");
  hv::requireOrThrow(
      bins[1].density_variance / (bins[1].mean_density * bins[1].mean_density) < 0.4,
      "Noh: radial compression shell is too asymmetric");
}

void testGreshoVortexGuard() {
  constexpr double gamma = 1.4;
  cosmosim::hydro::HydroPatchGeometry geometry = hv::makePeriodicCartesianPatch(16U, 16U, 1U);
  cosmosim::hydro::HydroConservedStateSoa conserved(geometry.cellCount());

  double initial_lz = 0.0;
  double initial_vtheta_inner = 0.0;
  std::size_t inner_count = 0;
  for (std::size_t cell = 0; cell < geometry.cellCount(); ++cell) {
    const hv::CellCenter center = hv::cellCenter(geometry, cell);
    const double dx = center.x_comoving - 0.5;
    const double dy = center.y_comoving - 0.5;
    const double radius = std::sqrt(dx * dx + dy * dy);
    double vtheta = 0.0;
    double pressure = 3.0 + 4.0 * std::log(2.0);
    if (radius < 0.2) {
      vtheta = 5.0 * radius;
      pressure = 5.0 + 12.5 * radius * radius;
    } else if (radius < 0.4) {
      vtheta = 2.0 - 5.0 * radius;
      pressure = 9.0 + 12.5 * radius * radius - 20.0 * radius + 4.0 * std::log(5.0 * radius);
    }
    cosmosim::hydro::HydroPrimitiveState primitive{};
    primitive.rho_comoving = 1.0;
    if (radius > 1.0e-14) {
      primitive.vel_x_peculiar = -vtheta * dy / radius;
      primitive.vel_y_peculiar = vtheta * dx / radius;
    }
    primitive.pressure_comoving = pressure;
    hv::storePrimitive(conserved, cell, primitive, gamma);

    initial_lz += (dx * primitive.vel_y_peculiar - dy * primitive.vel_x_peculiar) *
        primitive.rho_comoving * geometry.cell_volume_comoving;
    if (radius > 0.08 && radius < 0.18) {
      initial_vtheta_inner += vtheta;
      ++inner_count;
    }
  }

  hv::advanceHydroSteps(conserved, geometry, hv::HydroStepConfig{
      .gamma = gamma,
      .dt_code = 5.0e-5,
      .step_count = 24,
      .rho_floor = 1.0e-9,
      .pressure_floor = 1.0e-9,
      .limiter = cosmosim::hydro::HydroSlopeLimiter::kVanLeer});

  hv::requireFinitePositiveState(conserved, gamma, "Gresho");
  double final_lz = 0.0;
  double final_vtheta_inner = 0.0;
  for (std::size_t cell = 0; cell < geometry.cellCount(); ++cell) {
    const hv::CellCenter center = hv::cellCenter(geometry, cell);
    const double dx = center.x_comoving - 0.5;
    const double dy = center.y_comoving - 0.5;
    const double radius = std::sqrt(dx * dx + dy * dy);
    const cosmosim::hydro::HydroPrimitiveState primitive =
        cosmosim::hydro::HydroCoreSolver::primitiveFromConserved(conserved.loadCell(cell), gamma);
    final_lz += (dx * primitive.vel_y_peculiar - dy * primitive.vel_x_peculiar) *
        primitive.rho_comoving * geometry.cell_volume_comoving;
    if (radius > 0.08 && radius < 0.18) {
      final_vtheta_inner +=
          (-primitive.vel_x_peculiar * dy + primitive.vel_y_peculiar * dx) / std::max(radius, 1.0e-14);
    }
  }
  initial_vtheta_inner /= static_cast<double>(inner_count);
  final_vtheta_inner /= static_cast<double>(inner_count);

  hv::requireOrThrow(
      relativeDifference(final_lz, initial_lz) < 0.15,
      "Gresho: angular momentum drift exceeds CI guard");
  hv::requireOrThrow(
      final_vtheta_inner > 0.75 * initial_vtheta_inner,
      "Gresho: inner vortex profile diffused catastrophically");
}

void testKelvinHelmholtzGuard() {
  constexpr double gamma = 1.4;
  cosmosim::hydro::HydroPatchGeometry geometry = hv::makePeriodicCartesianPatch(24U, 12U, 1U);
  cosmosim::hydro::HydroConservedStateSoa conserved(geometry.cellCount());

  double initial_transverse_ke = 0.0;
  for (std::size_t cell = 0; cell < geometry.cellCount(); ++cell) {
    const hv::CellCenter center = hv::cellCenter(geometry, cell);
    const double y_centered = center.y_comoving - 0.5;
    const bool dense_layer = std::abs(y_centered) < 0.25;
    const double perturb = 0.025 * std::sin(4.0 * hv::k_pi * center.x_comoving) *
        (std::exp(-std::pow((y_centered - 0.25) / 0.07, 2.0)) +
         std::exp(-std::pow((y_centered + 0.25) / 0.07, 2.0)));
    cosmosim::hydro::HydroPrimitiveState primitive{};
    primitive.rho_comoving = dense_layer ? 2.0 : 1.0;
    primitive.vel_x_peculiar = dense_layer ? 0.25 : -0.25;
    primitive.vel_y_peculiar = perturb;
    primitive.pressure_comoving = 2.5;
    hv::storePrimitive(conserved, cell, primitive, gamma);
    initial_transverse_ke +=
        0.5 * primitive.rho_comoving * primitive.vel_y_peculiar * primitive.vel_y_peculiar *
        geometry.cell_volume_comoving;
  }

  hv::advanceHydroSteps(conserved, geometry, hv::HydroStepConfig{
      .gamma = gamma,
      .dt_code = 7.5e-5,
      .step_count = 50,
      .rho_floor = 1.0e-9,
      .pressure_floor = 1.0e-9,
      .limiter = cosmosim::hydro::HydroSlopeLimiter::kVanLeer});

  hv::requireFinitePositiveState(conserved, gamma, "Kelvin-Helmholtz");
  double final_transverse_ke = 0.0;
  double interface_density_min = 10.0;
  double interface_density_max = 0.0;
  for (std::size_t cell = 0; cell < geometry.cellCount(); ++cell) {
    const hv::CellCenter center = hv::cellCenter(geometry, cell);
    const double y_abs = std::abs(center.y_comoving - 0.5);
    const cosmosim::hydro::HydroPrimitiveState primitive =
        cosmosim::hydro::HydroCoreSolver::primitiveFromConserved(conserved.loadCell(cell), gamma);
    final_transverse_ke +=
        0.5 * primitive.rho_comoving * primitive.vel_y_peculiar * primitive.vel_y_peculiar *
        geometry.cell_volume_comoving;
    if (std::abs(y_abs - 0.25) < 0.13) {
      interface_density_min = std::min(interface_density_min, primitive.rho_comoving);
      interface_density_max = std::max(interface_density_max, primitive.rho_comoving);
    }
  }

  hv::requireOrThrow(
      final_transverse_ke > 0.55 * initial_transverse_ke,
      "Kelvin-Helmholtz: transverse perturbation was erased catastrophically");
  hv::requireOrThrow(
      interface_density_min < 1.8 && interface_density_max > 1.2,
      "Kelvin-Helmholtz: shear-layer density contrast vanished or inverted");
}

void testEvrardCollapseToyGuard() {
  constexpr double gamma = 5.0 / 3.0;
  cosmosim::hydro::HydroPatchGeometry geometry = hv::makePeriodicCartesianPatch(8U, 8U, 8U);
  cosmosim::hydro::HydroConservedStateSoa conserved(geometry.cellCount());
  std::vector<double> accel_x(geometry.cellCount(), 0.0);
  std::vector<double> accel_y(geometry.cellCount(), 0.0);
  std::vector<double> accel_z(geometry.cellCount(), 0.0);

  for (std::size_t cell = 0; cell < geometry.cellCount(); ++cell) {
    const hv::CellCenter center = hv::cellCenter(geometry, cell);
    const double dx = center.x_comoving - 0.5;
    const double dy = center.y_comoving - 0.5;
    const double dz = center.z_comoving - 0.5;
    const double radius = std::sqrt(dx * dx + dy * dy + dz * dz);
    const double softened_r3 = std::pow(radius * radius + 0.04 * 0.04, 1.5);
    cosmosim::hydro::HydroPrimitiveState primitive{};
    primitive.rho_comoving = radius < 0.48 ? 1.0 / std::max(radius, 0.08) : 0.05;
    primitive.pressure_comoving = 0.015;
    hv::storePrimitive(conserved, cell, primitive, gamma);
    accel_x[cell] = -0.08 * dx / softened_r3;
    accel_y[cell] = -0.08 * dy / softened_r3;
    accel_z[cell] = -0.08 * dz / softened_r3;
  }

  const auto initial = hv::conservedTotals(conserved, geometry);
  FixedGravitySource source;
  source.accel_x_code = accel_x;
  source.accel_y_code = accel_y;
  source.accel_z_code = accel_z;
  const std::array<const cosmosim::hydro::HydroSourceTerm*, 1> sources{&source};
  hv::advanceHydroSteps(
      conserved,
      geometry,
      hv::HydroStepConfig{
          .gamma = gamma,
          .dt_code = 4.0e-5,
          .step_count = 36,
          .rho_floor = 1.0e-9,
          .pressure_floor = 1.0e-9,
          .limiter = cosmosim::hydro::HydroSlopeLimiter::kMinmod},
      sources);

  hv::requireFinitePositiveState(conserved, gamma, "Evrard toy");
  const auto final = hv::conservedTotals(conserved, geometry);
  hv::requireOrThrow(
      final.total_energy > 0.5 * initial.total_energy && final.total_energy < 2.5 * initial.total_energy,
      "Evrard toy: source-coupled energy left the bounded CI envelope");

  const std::array<double, 4> edges{0.0, 0.18, 0.34, 0.70};
  const std::vector<hv::RadialBin> bins = hv::radialBins(conserved, geometry, gamma, edges);
  hv::requireOrThrow(
      bins[1].mean_radial_velocity < -1.0e-4,
      "Evrard toy: analytic gravity source did not produce inward collapse trend");
  hv::requireOrThrow(
      bins[0].mean_density > bins[2].mean_density,
      "Evrard toy: central concentration was not preserved under collapse guard");
}

}  // namespace

int main() {
  testSedovBlastGuard();
  testNohConvergingInflowGuard();
  testGreshoVortexGuard();
  testKelvinHelmholtzGuard();
  testEvrardCollapseToyGuard();
  return 0;
}
