#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <vector>

#include "cosmosim/physics/metal_diffusion.hpp"

namespace {

using cosmosim::physics::MetalDiffusionBoundaryKind;
using cosmosim::physics::MetalDiffusionCell;
using cosmosim::physics::MetalDiffusionConfig;
using cosmosim::physics::MetalDiffusionFace;
using cosmosim::physics::MetalDiffusionModel;

constexpr double kPi = 3.141592653589793238462643383279502884;

MetalDiffusionConfig makeConfig(cosmosim::core::MetalDiffusionTimeIntegrator integrator) {
  MetalDiffusionConfig config;
  config.enabled = true;
  config.model = cosmosim::core::MetalDiffusionModel::kSmagorinsky;
  config.time_integrator = integrator;
  config.smagorinsky_coefficient = 0.02;
  config.parabolic_cfl = 0.4;
  config.max_subcycles = 1024;
  config.max_rkl_stages = 128;
  return config;
}

struct PeriodicLine {
  std::vector<MetalDiffusionCell> cells;
  std::vector<MetalDiffusionFace> faces;
};

PeriodicLine makePeriodicLine(std::size_t count, double amplitude) {
  PeriodicLine line;
  line.cells.resize(count);
  line.faces.resize(count);
  const double dx = 1.0 / static_cast<double>(count);
  for (std::size_t i = 0; i < count; ++i) {
    const double x = (static_cast<double>(i) + 0.5) * dx;
    auto& cell = line.cells[i];
    cell.gas_cell_id = static_cast<std::uint64_t>(i + 1U);
    cell.gas_mass_code = dx;
    cell.density_code = 1.0;
    cell.volume_code = dx;
    cell.filter_length_code = 1.0;
    cell.metal_mass_code = dx * (0.1 + amplitude * std::sin(2.0 * kPi * x));
    // diag(+1/2,-1/2,0) has |S*|=1, therefore kappa=C_mix.
    cell.velocity_gradient.grad[0][0] = 0.5;
    cell.velocity_gradient.grad[1][1] = -0.5;
    line.faces[i] = MetalDiffusionFace{
        .left_cell = static_cast<std::uint32_t>(i),
        .right_cell = static_cast<std::uint32_t>((i + 1U) % count),
        .area_code = 1.0,
        .center_distance_code = dx,
        .boundary_kind = i + 1U == count
            ? MetalDiffusionBoundaryKind::kPeriodic
            : MetalDiffusionBoundaryKind::kInternal,
    };
  }
  return line;
}

double totalMetals(const PeriodicLine& line) {
  return std::accumulate(
      line.cells.begin(), line.cells.end(), 0.0,
      [](double sum, const MetalDiffusionCell& cell) {
        return sum + cell.metal_mass_code;
      });
}

double sineAmplitude(const PeriodicLine& line) {
  const double dx = 1.0 / static_cast<double>(line.cells.size());
  double projection = 0.0;
  for (std::size_t i = 0; i < line.cells.size(); ++i) {
    const double x = (static_cast<double>(i) + 0.5) * dx;
    const double z = line.cells[i].metal_mass_code / line.cells[i].gas_mass_code;
    projection += (z - 0.1) * std::sin(2.0 * kPi * x);
  }
  return 2.0 * projection / static_cast<double>(line.cells.size());
}

void testConstantAndZeroCoefficientIdentity() {
  auto line = makePeriodicLine(16, 0.0);
  auto config = makeConfig(cosmosim::core::MetalDiffusionTimeIntegrator::kExplicitSubcycling);
  MetalDiffusionModel model(config);
  const auto before = totalMetals(line);
  const auto report = model.advance(line.cells, line.faces, 0.1);
  assert(report.limited_faces == 0U);
  assert(std::abs(totalMetals(line) - before) < 1.0e-14);
  for (const auto& cell : line.cells) {
    assert(std::abs(cell.metal_mass_code / cell.gas_mass_code - 0.1) < 1.0e-14);
  }

  auto zero_line = makePeriodicLine(16, 0.03);
  config.smagorinsky_coefficient = 0.0;
  MetalDiffusionModel zero_model(config);
  std::vector<double> initial;
  for (const auto& cell : zero_line.cells) initial.push_back(cell.metal_mass_code);
  const auto zero_report = zero_model.advance(zero_line.cells, zero_line.faces, 1.0);
  assert(zero_report.faces_evaluated == 0U);
  for (std::size_t i = 0; i < initial.size(); ++i) {
    assert(zero_line.cells[i].metal_mass_code == initial[i]);
  }
}

void testSineDecayConservationAndPositivity() {
  auto line = makePeriodicLine(64, 0.02);
  const MetalDiffusionModel model(
      makeConfig(cosmosim::core::MetalDiffusionTimeIntegrator::kExplicitSubcycling));
  const double before = totalMetals(line);
  const double amplitude_before = sineAmplitude(line);
  constexpr double dt = 5.0e-3;
  const auto report = model.advance(line.cells, line.faces, dt);
  const double amplitude_after = sineAmplitude(line);
  const double expected = amplitude_before *
      std::exp(-0.02 * std::pow(2.0 * kPi, 2.0) * dt);
  assert(std::abs(amplitude_after - expected) < 2.0e-5);
  assert(std::abs(totalMetals(line) - before) < 2.0e-14);
  assert(std::abs(report.conservation_residual_code) < 2.0e-14);
  for (const auto& cell : line.cells) {
    assert(cell.metal_mass_code >= 0.0);
    assert(cell.metal_mass_code <= cell.gas_mass_code);
  }
}

void testExplicitAndRkl2Agreement() {
  auto explicit_line = makePeriodicLine(48, 0.025);
  auto rkl_line = explicit_line;
  const MetalDiffusionModel explicit_model(
      makeConfig(cosmosim::core::MetalDiffusionTimeIntegrator::kExplicitSubcycling));
  const MetalDiffusionModel rkl_model(
      makeConfig(cosmosim::core::MetalDiffusionTimeIntegrator::kRkl2));
  const double stable = explicit_model.stableExplicitTimestepCode(
      explicit_line.cells, explicit_line.faces);
  const double dt = 6.0 * stable;
  const auto explicit_report = explicit_model.advance(explicit_line.cells, explicit_line.faces, dt);
  assert(explicit_report.subcycles >= 2U);
  const auto report = rkl_model.advance(rkl_line.cells, rkl_line.faces, dt);
  assert(report.rkl_stages >= 2U);
  for (std::size_t i = 0; i < explicit_line.cells.size(); ++i) {
    assert(std::abs(
        explicit_line.cells[i].metal_mass_code - rkl_line.cells[i].metal_mass_code) < 2.0e-7);
  }
}

void testSolidBodyRotationHasNoSmagorinskyDiffusion() {
  cosmosim::physics::MetalDiffusionVelocityGradient gradient;
  gradient.grad[0][1] = -3.0;
  gradient.grad[1][0] = 3.0;
  assert(cosmosim::physics::traceFreeStrainMagnitude(gradient) < 1.0e-15);

  MetalDiffusionCell cell;
  cell.filter_length_code = 2.0;
  cell.velocity_gradient = gradient;
  const auto config = makeConfig(
      cosmosim::core::MetalDiffusionTimeIntegrator::kExplicitSubcycling);
  assert(cosmosim::physics::smagorinskyMetalDiffusivityCode(config, cell) == 0.0);
}

void testStepProfileRemainsBounded() {
  auto line = makePeriodicLine(32, 0.0);
  for (std::size_t i = 0; i < line.cells.size(); ++i) {
    const double z = i < line.cells.size() / 2U ? 0.0 : 0.9;
    line.cells[i].metal_mass_code = line.cells[i].gas_mass_code * z;
  }
  const double before = totalMetals(line);
  const MetalDiffusionModel model(
      makeConfig(cosmosim::core::MetalDiffusionTimeIntegrator::kExplicitSubcycling));
  const double dt = 10.0 * model.stableExplicitTimestepCode(line.cells, line.faces);
  const auto report = model.advance(line.cells, line.faces, dt);
  assert(report.subcycles >= 10U);
  assert(std::abs(totalMetals(line) - before) < 2.0e-14);
  for (const auto& cell : line.cells) {
    const double z = cell.metal_mass_code / cell.gas_mass_code;
    assert(z >= 0.0 && z <= 1.0);
  }
}

void testCompactViewMatchesLegacyAndWorkspacePlateaus() {
  auto legacy_line = makePeriodicLine(64, 0.02);
  auto compact_line = legacy_line;
  const auto config = makeConfig(
      cosmosim::core::MetalDiffusionTimeIntegrator::kExplicitSubcycling);
  const MetalDiffusionModel model(config);
  const double stable = model.stableExplicitTimestepCode(
      legacy_line.cells, legacy_line.faces);
  const double dt = 3.0 * stable;
  const auto legacy_report = model.advance(
      legacy_line.cells, legacy_line.faces, dt);

  std::vector<double> gas_mass(compact_line.cells.size());
  std::vector<double> metal_mass(compact_line.cells.size());
  std::vector<double> rho_kappa(compact_line.cells.size());
  for (std::size_t i = 0U; i < compact_line.cells.size(); ++i) {
    gas_mass[i] = compact_line.cells[i].gas_mass_code;
    metal_mass[i] = compact_line.cells[i].metal_mass_code;
    rho_kappa[i] = compact_line.cells[i].density_code *
        cosmosim::physics::smagorinskyMetalDiffusivityCode(
            config, compact_line.cells[i]);
  }
  cosmosim::physics::MetalDiffusionWorkspace workspace;
  const auto compact_report = model.advanceFromView(
      cosmosim::physics::MetalDiffusionFieldView{
          .gas_mass_code = gas_mass,
          .metal_mass_code = metal_mass,
          .rho_kappa_code = rho_kappa,
      },
      compact_line.faces,
      dt,
      workspace);
  assert(std::abs(
      compact_report.conservation_residual_code -
      legacy_report.conservation_residual_code) < 1.0e-14);
  for (std::size_t i = 0U; i < metal_mass.size(); ++i) {
    assert(std::abs(
        metal_mass[i] - legacy_line.cells[i].metal_mass_code) < 1.0e-14);
  }
  const std::uint64_t required = model.requiredWorkspaceBytes(metal_mass.size());
  assert(required == 5U * metal_mass.size() * sizeof(double));
  assert(workspace.ownedCapacityBytes() >= required);
  const std::uint64_t retained_before = workspace.ownedCapacityBytes();
  const std::uint64_t high_water_before = workspace.highWaterBytes();
  (void)model.advanceFromView(
      cosmosim::physics::MetalDiffusionFieldView{
          .gas_mass_code = gas_mass,
          .metal_mass_code = metal_mass,
          .rho_kappa_code = rho_kappa,
      },
      compact_line.faces,
      dt,
      workspace);
  assert(workspace.ownedCapacityBytes() == retained_before);
  assert(workspace.highWaterBytes() == high_water_before);
}

}  // namespace

int main() {
  testConstantAndZeroCoefficientIdentity();
  testSineDecayConservationAndPositivity();
  testExplicitAndRkl2Agreement();
  testSolidBodyRotationHasNoSmagorinskyDiffusion();
  testStepProfileRemainsBounded();
  testCompactViewMatchesLegacyAndWorkspacePlateaus();
  return 0;
}
