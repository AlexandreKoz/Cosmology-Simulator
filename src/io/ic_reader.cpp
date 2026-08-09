#include "cosmosim/io/ic_reader.hpp"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>

#include "cosmosim/core/checked_arithmetic.hpp"
#include "cosmosim/core/constants.hpp"
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
  result.report.provenance_authority = "generated_runtime_config";
  const std::size_t total_particle_count = core::checkedSizeAdd(
      dark_matter_particle_count, gas_particle_count,
      "generated isolated IC particle count");
  if (total_particle_count != 0U && particle_id_seed == 0U) {
    throw std::invalid_argument(
        "generated isolated IC particle_id_seed must be nonzero");
  }
  if (total_particle_count > 0U &&
      particle_id_seed > std::numeric_limits<std::uint64_t>::max() -
          static_cast<std::uint64_t>(total_particle_count - 1U)) {
    throw std::overflow_error(
        "generated isolated IC particle IDs overflow uint64");
  }
  result.state.resizeParticles(total_particle_count);
  result.state.resizeCells(gas_particle_count);
  result.state.metadata.scale_factor = 1.0;
  result.state.metadata.run_name = config.output.run_name;
  result.state.metadata.snapshot_stem = config.output.output_stem;
  result.state.metadata.restart_stem = config.output.restart_stem;

  const core::UnitSystem target_units =
      core::makeUnitSystem(config.units.length_unit, config.units.mass_unit, config.units.velocity_unit);
  const double box_x_code = target_units.lengthSiToCode(
      config.cosmology.box_size_x_mpc_comoving * core::constants::k_megaparsec_si);
  const double box_y_code = target_units.lengthSiToCode(
      config.cosmology.box_size_y_mpc_comoving * core::constants::k_megaparsec_si);
  const double box_z_code = target_units.lengthSiToCode(
      config.cosmology.box_size_z_mpc_comoving * core::constants::k_megaparsec_si);
  if (!(box_x_code > 0.0) || !(box_y_code > 0.0) || !(box_z_code > 0.0)) {
    throw std::invalid_argument(
        "generated isolated IC requires positive configured box extents");
  }
  const auto centeredCoordinate = [](std::size_t index, std::size_t count,
                                     double extent) {
    if (count == 0U) return 0.5 * extent;
    return extent * (static_cast<double>(index) + 0.5) /
        static_cast<double>(count);
  };

  std::size_t global_index = 0;
  for (std::size_t i = 0; i < gas_particle_count; ++i, ++global_index) {
    const double x = centeredCoordinate(i, gas_particle_count, box_x_code);
    result.state.particle_sidecar.particle_id[global_index] = particle_id_seed + global_index;
    result.state.particle_sidecar.species_tag[global_index] = k_species_gas;
    result.state.particle_sidecar.owning_rank[global_index] = 0;
    result.state.particles.position_x_comoving[global_index] = x;
    result.state.particles.position_y_comoving[global_index] = 0.5 * box_y_code;
    result.state.particles.position_z_comoving[global_index] = 0.5 * box_z_code;
    result.state.particles.velocity_x_peculiar[global_index] = 10.0;
    result.state.particles.velocity_y_peculiar[global_index] = 0.0;
    result.state.particles.velocity_z_peculiar[global_index] = 0.0;
    result.state.particles.mass_code[global_index] = target_units.massSiToCode(1.0e34);
    result.state.particles.time_bin[global_index] = 0;

    result.state.cells.center_x_comoving[i] = x;
    result.state.cells.center_y_comoving[i] = 0.5 * box_y_code;
    result.state.cells.center_z_comoving[i] = 0.5 * box_z_code;
    result.state.cells.mass_code[i] = result.state.particles.mass_code[global_index];
    result.state.cells.time_bin[i] = 0;
    result.state.cells.patch_index[i] = 0;
    result.state.gas_cells.density_code[i] = 1.0;
    result.state.gas_cells.pressure_code[i] = 1.0;
    result.state.gas_cells.internal_energy_code[i] = 1.5;
    result.state.gas_cells.metal_mass_code[i] = 0.0;
    result.state.gas_cells.temperature_code[i] = 1.0e4;
    result.state.gas_cells.sound_speed_code[i] = 1.0;
  }

  for (std::size_t i = 0; i < dark_matter_particle_count; ++i, ++global_index) {
    result.state.particle_sidecar.particle_id[global_index] = particle_id_seed + global_index;
    result.state.particle_sidecar.species_tag[global_index] = k_species_dark_matter;
    result.state.particle_sidecar.owning_rank[global_index] = 0;
    result.state.particles.position_x_comoving[global_index] =
        centeredCoordinate(i, dark_matter_particle_count, box_x_code);
    result.state.particles.position_y_comoving[global_index] = 0.25 * box_y_code;
    result.state.particles.position_z_comoving[global_index] = 0.25 * box_z_code;
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
  const std::size_t dm_count = core::checkedSizeMultiply(
      particle_count_per_axis, particle_count_per_axis,
      "generated isolated IC particle_count_per_axis^2");
  const std::size_t gas_count = particle_count_per_axis;
  IcReadResult result = buildGeneratedIsolatedIc(config, dm_count, gas_count);
  result.report.defaulted_fields.push_back("gas/internal_energy_code=uniform_default");
  return result;
}


}  // namespace cosmosim::io
