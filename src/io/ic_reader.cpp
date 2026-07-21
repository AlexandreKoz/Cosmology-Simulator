#include "cosmosim/io/ic_reader.hpp"

#include <cstddef>
#include <cstdint>

#include "cosmosim/core/units.hpp"

namespace cosmosim::io {
namespace {

constexpr std::uint32_t k_species_dark_matter =
    static_cast<std::uint32_t>(core::ParticleSpecies::kDarkMatter);
constexpr std::uint32_t k_species_gas =
    static_cast<std::uint32_t>(core::ParticleSpecies::kGas);

void fillSpeciesLedger(core::SimulationState& state) {
  state.species.count_by_species = {};
  for (std::uint32_t species_tag : state.particle_sidecar.species_tag) {
    if (species_tag < state.species.count_by_species.size()) {
      state.species.count_by_species[species_tag] += 1;
    }
  }
}

}  // namespace

IcManifest makeGadgetArepoBridgeV1Manifest(
    const std::filesystem::path& ic_path,
    const IcSchemaSummary& schema) {
  IcManifest manifest;
  manifest.dialect = IcDialect::kGadgetArepoBridgeV1;
  manifest.dialect_version = "1";
  manifest.converter_version = "runtime_inspection_required";
  manifest.source_files = {ic_path.lexically_normal()};
  manifest.num_part_this_file = {schema.count_by_type};
  manifest.num_part_total = schema.total_count_by_type;
  manifest.num_part_total_high_word = schema.total_count_high_word;
  manifest.num_files_per_snapshot = schema.num_files_per_snapshot;
  manifest.mass_table = schema.mass_table;
  manifest.box_size = schema.box_size;
  manifest.scale_factor = schema.scale_factor;
  manifest.redshift = schema.redshift;
  manifest.omega_matter = schema.omega_matter;
  manifest.omega_lambda = schema.omega_lambda;
  manifest.hubble_param = schema.hubble_param;
  // This is deliberately a policy template, not an observed audit manifest.
  // readGadgetArepoHdf5Ic() inspects every real dataset and fills datatype,
  // dimensions, hashes, aliases, and conversion equations before validation.
  return manifest;
}

IcReadResult buildGeneratedIsolatedIc(
    const core::SimulationConfig& config,
    std::size_t dark_matter_particle_count,
    std::size_t gas_particle_count,
    std::uint64_t particle_id_seed) {
  IcReadResult result;
  const std::size_t total_particle_count = dark_matter_particle_count + gas_particle_count;
  result.state.resizeParticles(total_particle_count);
  result.state.resizeCells(gas_particle_count);
  result.state.metadata.scale_factor = 1.0;
  result.state.metadata.run_name = config.output.run_name;
  result.state.metadata.snapshot_stem = config.output.output_stem;
  result.state.metadata.restart_stem = config.output.restart_stem;

  const core::UnitSystem target_units =
      core::makeUnitSystem(config.units.length_unit, config.units.mass_unit, config.units.velocity_unit);

  std::size_t global_index = 0;
  for (std::size_t i = 0; i < gas_particle_count; ++i, ++global_index) {
    const double x = static_cast<double>(i) * 0.01;
    result.state.particle_sidecar.particle_id[global_index] = particle_id_seed + global_index;
    result.state.particle_sidecar.species_tag[global_index] = k_species_gas;
    result.state.particle_sidecar.owning_rank[global_index] = 0;
    result.state.particles.position_x_comoving[global_index] = x;
    result.state.particles.position_y_comoving[global_index] = 0.5;
    result.state.particles.position_z_comoving[global_index] = 0.5;
    result.state.particles.velocity_x_peculiar[global_index] = 10.0;
    result.state.particles.velocity_y_peculiar[global_index] = 0.0;
    result.state.particles.velocity_z_peculiar[global_index] = 0.0;
    result.state.particles.mass_code[global_index] = target_units.massSiToCode(1.0e34);
    result.state.particles.time_bin[global_index] = 0;

    result.state.cells.center_x_comoving[i] = x;
    result.state.cells.center_y_comoving[i] = 0.5;
    result.state.cells.center_z_comoving[i] = 0.5;
    result.state.cells.mass_code[i] = result.state.particles.mass_code[global_index];
    result.state.cells.time_bin[i] = 0;
    result.state.cells.patch_index[i] = 0;
    result.state.gas_cells.density_code[i] = 1.0;
    result.state.gas_cells.pressure_code[i] = 1.0;
    result.state.gas_cells.internal_energy_code[i] = 1.5;
    result.state.gas_cells.temperature_code[i] = 1.0e4;
    result.state.gas_cells.sound_speed_code[i] = 1.0;
  }

  for (std::size_t i = 0; i < dark_matter_particle_count; ++i, ++global_index) {
    result.state.particle_sidecar.particle_id[global_index] = particle_id_seed + global_index;
    result.state.particle_sidecar.species_tag[global_index] = k_species_dark_matter;
    result.state.particle_sidecar.owning_rank[global_index] = 0;
    result.state.particles.position_x_comoving[global_index] = static_cast<double>(i) * 0.01;
    result.state.particles.position_y_comoving[global_index] = 0.0;
    result.state.particles.position_z_comoving[global_index] = 0.0;
    result.state.particles.velocity_x_peculiar[global_index] = 0.0;
    result.state.particles.velocity_y_peculiar[global_index] = 0.0;
    result.state.particles.velocity_z_peculiar[global_index] = 0.0;
    result.state.particles.mass_code[global_index] = target_units.massSiToCode(1.0e36);
    result.state.particles.time_bin[global_index] = 0;
  }

  fillSpeciesLedger(result.state);
  result.state.rebuildSpeciesIndex();
  if (gas_particle_count > 0) {
    result.state.refreshGasCellIdentityFromParticleOrder();
  }
  return result;
}

IcReadResult convertGeneratedIsolatedIcToState(
    const core::SimulationConfig& config,
    std::size_t particle_count_per_axis) {
  const std::size_t dm_count = particle_count_per_axis * particle_count_per_axis;
  const std::size_t gas_count = particle_count_per_axis;
  IcReadResult result = buildGeneratedIsolatedIc(config, dm_count, gas_count);
  result.report.defaulted_fields.push_back("gas/internal_energy_code=uniform_default");
  return result;
}


}  // namespace cosmosim::io
