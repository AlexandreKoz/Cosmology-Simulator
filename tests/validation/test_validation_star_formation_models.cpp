#include <cassert>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <memory>
#include <string>

#include "cosmosim/core/config.hpp"
#include "cosmosim/core/units.hpp"
#include "cosmosim/physics/effective_multiphase_ism.hpp"
#include "cosmosim/physics/star_formation.hpp"

namespace {

using cosmosim::physics::EffectiveMultiphaseEosTable;

std::shared_ptr<const EffectiveMultiphaseEosTable> makeEffectiveTable() {
  cosmosim::core::PhysicsConfig physics;
  physics.star_formation_model =
      cosmosim::core::StarFormationModelKind::kEffectiveMultiphaseTngLike;
  physics.sf_effective_eos_table_bins = 256U;
  physics.sf_effective_eos_max_density_ratio = 1.0e6;
  physics.sf_effective_q_eos = 0.3;
  const auto units = cosmosim::core::makeUnitSystem("kpc", "msun", "km_s");
  return std::make_shared<const EffectiveMultiphaseEosTable>(
      cosmosim::physics::makeEffectiveMultiphaseEosConfig(physics),
      units,
      cosmosim::physics::makeEffectiveIsmReferenceCoolingProvider(physics));
}

cosmosim::physics::StarFormationCellInput makeAdaptiveCell(
    std::uint64_t id,
    double density,
    double scale) {
  cosmosim::physics::StarFormationCellInput cell;
  cell.gas_cell_id = id;
  cell.gas_mass_code = density * scale * scale * scale;
  cell.gas_density_code = density;
  cell.cell_volume_code = scale * scale * scale;
  cell.gas_temperature_k = 100.0;
  cell.gas_sound_speed_code = 0.02;
  cell.velocity_divergence_code = -0.5;
  cell.velocity_gradient_frobenius_sq_code = 0.01;
  cell.gas_specific_internal_energy_code = 1.0e-3;
  cell.is_active = true;
  cell.is_owned = true;
  cell.is_leaf = true;
  return cell;
}

void writeEffectiveCurve(
    const EffectiveMultiphaseEosTable& table,
    const std::filesystem::path& output) {
  std::ofstream stream(output);
  stream << std::setprecision(17);
  stream << "density_ratio,density_phys_code,cold_mass_fraction,specific_internal_energy_eff_code,"
            "pressure_phys_code,signal_speed_code,star_formation_timescale_code,specific_sfr_code\n";
  for (const auto& entry : table.entries()) {
    const double specific_sfr = entry.cold_mass_fraction / entry.star_formation_timescale_code;
    stream << entry.density_ratio << ',' << entry.density_phys_code << ','
           << entry.cold_mass_fraction << ',' << entry.specific_internal_energy_eff_code << ','
           << entry.pressure_phys_code << ',' << std::sqrt(entry.signal_speed_squared_code) << ','
           << entry.star_formation_timescale_code << ',' << specific_sfr << '\n';
  }
}

void writeAdaptiveResolution(
    const std::filesystem::path& output) {
  cosmosim::physics::StarFormationConfig config;
  config.model = cosmosim::core::StarFormationModelKind::kAdaptiveBoundJeans;
  config.enabled = true;
  config.epsilon_ff = 1.0;
  config.bound_alpha_vir_max = 1.0;
  config.require_converging_flow = true;
  config.jeans_mass_floor_code = 0.0;
  config.newton_g_code = 1.0;
  config.density_is_comoving = false;
  config.geometry_is_comoving = false;
  config.stochastic_spawning = false;
  config.target_star_particle_mass_code = 1.0e-6;
  config.min_star_particle_mass_code = 1.0e-12;
  config.max_star_particle_mass_code = 1.0e30;
  config.max_fractional_mass_conversion = 0.25;
  cosmosim::physics::StarFormationModel model(config);

  std::ofstream stream(output);
  stream << std::setprecision(17);
  stream << "density_code,cell_scale_code,eligible,virial_parameter,jeans_mass_code,"
            "free_fall_time_code,expected_spawn_mass_code\n";
  std::uint64_t id = 10000U;
  for (double density : {0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0}) {
    for (double scale : {0.25, 0.5, 1.0, 2.0}) {
      const auto outcome = model.sampleCellOutcome(
          makeAdaptiveCell(id++, density, scale), 1.0e-3, 7U, 0U, 1.0);
      stream << density << ',' << scale << ',' << (outcome.eligible ? 1 : 0) << ','
             << outcome.virial_parameter << ',' << outcome.jeans_mass_code << ','
             << outcome.free_fall_time_code << ',' << outcome.expected_spawn_mass_code << '\n';
    }
  }
}

void writeStochasticEnsemble(
    const EffectiveMultiphaseEosTable& table,
    const std::filesystem::path& output) {
  cosmosim::physics::StarFormationConfig config;
  config.model = cosmosim::core::StarFormationModelKind::kEffectiveMultiphaseTngLike;
  config.enabled = true;
  config.stochastic_spawning = true;
  config.random_seed = 919191U;
  config.target_star_particle_mass_code = 1.0;
  config.min_star_particle_mass_code = 1.0;
  config.max_star_particle_mass_code = 1.0;
  config.max_spawn_particles_per_cell_step = 1U;
  config.max_fractional_mass_conversion = 0.25;
  config.min_remaining_gas_fraction = 0.0;
  config.min_remaining_gas_mass_code = 0.0;
  auto shared = std::make_shared<const EffectiveMultiphaseEosTable>(table);
  cosmosim::physics::StarFormationModel model(config, shared);

  constexpr std::uint64_t sample_count = 20000U;
  double expected_sum = 0.0;
  double realized_sum = 0.0;
  const double rho = 10.0 * table.thresholdDensityPhysCode();
  const auto eq = table.lookup(rho);
  assert(eq.valid);
  for (std::uint64_t i = 0; i < sample_count; ++i) {
    cosmosim::physics::StarFormationCellInput cell;
    cell.gas_cell_id = 500000U + i;
    cell.gas_mass_code = 100.0;
    cell.gas_density_code = rho;
    cell.cell_volume_code = cell.gas_mass_code / rho;
    cell.gas_temperature_k = 1.0e4;
    cell.gas_sound_speed_code = std::sqrt(eq.entry.signal_speed_squared_code);
    cell.gas_specific_internal_energy_code = eq.entry.specific_internal_energy_eff_code;
    cell.is_active = true;
    cell.is_owned = true;
    cell.is_leaf = true;
    const auto outcome = model.sampleCellOutcome(cell, 1.0e-4, 11U, 0U, 1.0);
    expected_sum += outcome.expected_spawn_mass_code;
    realized_sum += outcome.spawned_mass_code;
  }
  const double expected_mean = expected_sum / static_cast<double>(sample_count);
  const double realized_mean = realized_sum / static_cast<double>(sample_count);
  const double standard_error = std::sqrt(
      std::max(expected_mean * (1.0 - std::min(expected_mean, 1.0)), 1.0e-12) /
      static_cast<double>(sample_count));
  assert(std::abs(realized_mean - expected_mean) < 6.0 * standard_error + 5.0e-3);
  std::ofstream stream(output);
  stream << std::setprecision(17);
  stream << "sample_count,expected_mean_mass_code,realized_mean_mass_code,standard_error_code\n";
  stream << sample_count << ',' << expected_mean << ',' << realized_mean << ',' << standard_error << '\n';
}

void writeTimestepSubdivision(const std::filesystem::path& output) {
  std::ofstream stream(output);
  stream << std::setprecision(17);
  stream << "substeps,analytic_fraction,composed_fraction,absolute_error\n";
  const double rate = 0.37;
  const double total_dt = 0.8;
  const double analytic = -std::expm1(-rate * total_dt);
  for (int n : {1, 2, 4, 8, 16, 32, 64}) {
    const double per_step = -std::expm1(-rate * total_dt / static_cast<double>(n));
    const double composed = 1.0 - std::pow(1.0 - per_step, n);
    const double error = std::abs(composed - analytic);
    assert(error < 1.0e-14);
    stream << n << ',' << analytic << ',' << composed << ',' << error << '\n';
  }
}

void validateInterpolation(const EffectiveMultiphaseEosTable& table) {
  for (double ratio : {1.01, 1.5, 3.0, 10.0, 100.0, 1.0e4}) {
    const double rho = ratio * table.thresholdDensityPhysCode();
    const auto direct = table.evaluateDirect(rho);
    const auto interpolated = table.lookup(rho);
    assert(direct.valid && interpolated.valid);
    const auto raw_threshold = table.evaluateDirect(table.thresholdDensityPhysCode());
    const double continuity_offset = table.entries().front().specific_internal_energy_eff_code -
        raw_threshold.specific_internal_energy_eff_code;
    const double normalized_direct_pressure = (table.config().adiabatic_index - 1.0) * rho *
        (direct.specific_internal_energy_eff_code + continuity_offset);
    const double pressure_rel = std::abs(interpolated.entry.pressure_phys_code - normalized_direct_pressure) /
        std::max(std::abs(normalized_direct_pressure), 1.0e-30);
    assert(pressure_rel < 0.03);
  }
}

}  // namespace

int main(int argc, char** argv) {
  const auto table = makeEffectiveTable();
  validateInterpolation(*table);
  if (argc > 1) {
    const std::filesystem::path output_dir = argv[1];
    std::filesystem::create_directories(output_dir);
    writeEffectiveCurve(*table, output_dir / "effective_multiphase_eos_curve.csv");
    writeAdaptiveResolution(output_dir / "adaptive_resolution_behavior.csv");
    writeStochasticEnsemble(*table, output_dir / "stochastic_ensemble.csv");
    writeTimestepSubdivision(output_dir / "timestep_subdivision.csv");
    std::ofstream metadata(output_dir / "validation_metadata.txt");
    metadata << "validation_scope=analytic_and_controlled_model_validation\n";
    metadata << "observational_calibration=false\n";
    metadata << "isolated_disk_long_run=false\n";
    metadata << "effective_eos_table_hash=" << table->tableHashHex() << "\n";
    metadata << "parameter_set=" << table->config().parameter_set << "\n";
  }
  return 0;
}
