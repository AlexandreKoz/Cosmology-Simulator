#include "cosmosim/io/ic_reader.hpp"
#include "io/internal/ic_file_set_common.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <cstdlib>
#include <filesystem>
#include <functional>
#include <fstream>
#include <iterator>
#include <limits>
#include <numeric>
#include <optional>
#include <set>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "cosmosim/core/build_config.hpp"
#include "io/internal/ic_byte_codec.hpp"
#include "io/internal/ic_canonical_bundle.hpp"
#include "io/internal/ic_conversion_catalog.hpp"
#include "io/internal/ic_hdf5_handle.hpp"
#include "io/internal/ic_reader_session.hpp"
#include "io/internal/ic_record_codec.hpp"
#include "io/internal/ic_stream_ingestion.hpp"
#include "cosmosim/core/units.hpp"
#include "cosmosim/parallel/distributed_memory.hpp"

#if COSMOSIM_ENABLE_HDF5
#include <hdf5.h>
#endif
#if COSMOSIM_ENABLE_MPI
#include <mpi.h>
#endif

namespace cosmosim::io::file_set_internal {

#if COSMOSIM_ENABLE_HDF5
using internal::IcReaderSession;

[[nodiscard]] const IcFieldManifest* findField(const IcManifest& manifest,std::size_t file_index,std::string_view path){const auto it=std::find_if(manifest.fields.begin(),manifest.fields.end(),[&](const IcFieldManifest& field){return field.source_file_index==file_index&&field.dataset_path==path;});return it==manifest.fields.end()?nullptr:&*it;}
[[nodiscard]] const IcFieldManifest& requireField(const IcManifest& manifest,std::size_t file_index,std::string_view path){const auto* field=findField(manifest,file_index,path);if(field==nullptr)throw std::runtime_error("manifest lacks inspected field "+std::string(path)+" for file "+std::to_string(file_index));return *field;}
[[nodiscard]] const IcMissingFieldContract& requireMissingFieldContract(
    const IcManifest& manifest,
    std::size_t file_index,
    std::string_view path) {
  const auto it = std::find_if(
      manifest.missing_field_contracts.begin(),
      manifest.missing_field_contracts.end(),
      [&](const IcMissingFieldContract& contract) {
        return contract.source_file_index == file_index &&
            contract.field_path == path;
      });
  if (it == manifest.missing_field_contracts.end()) {
    throw std::runtime_error(
        "manifest lacks explicit missing-field contract " +
        std::string(path) + " for file " + std::to_string(file_index));
  }
  return *it;
}
[[nodiscard]] double resolveMissingScalarValue(
    const IcMissingFieldContract& contract,
    double dialect_default) {
  if (contract.policy == IcMissingFieldPolicy::kUseConfigValue) {
    return contract.configured_value_code;
  }
  if (contract.policy == IcMissingFieldPolicy::kDialectDefinedDefault) {
    return dialect_default;
  }
  throw std::runtime_error(
      "unsupported resolved missing-field policy for " + contract.field_path);
}

void checkedCounterAdd(
    std::uint64_t& destination,
    std::uint64_t value,
    std::string_view counter_name) {
  if (destination > std::numeric_limits<std::uint64_t>::max() - value) {
    throw std::overflow_error(
        "IC " + std::string(counter_name) + " counter overflow");
  }
  destination += value;
}


void convertValues(
    std::vector<double>& values,
    const IcFieldManifest& field,
    const IcManifest& manifest,
    const core::UnitSystem& target_units) {
  const double factor =
      icFieldConversionMultiplier(field, manifest, target_units);
  for (double& value : values) {
    value *= factor;
    if (!std::isfinite(value)) {
      throw std::runtime_error(
          "non-finite converted value in " + field.dataset_path);
    }
  }
}

[[nodiscard]] std::uint64_t precedingRecordCount(
    const IcManifest& manifest,
    std::size_t file_index,
    std::size_t type_index) {
  std::uint64_t count = 0;
  for (std::size_t file = 0; file < manifest.num_part_this_file.size();
       ++file) {
    for (std::size_t type = 0; type < 6; ++type) {
      if (file == file_index && type == type_index) {
        return count;
      }
      if (count > std::numeric_limits<std::uint64_t>::max() -
                      manifest.num_part_this_file[file][type]) {
        throw std::overflow_error("generated IC ID offset overflow");
      }
      count += manifest.num_part_this_file[file][type];
    }
  }
  throw std::logic_error("invalid file/type offset");
}

[[nodiscard]] double convertedHeaderFieldCode(
    const IcManifest& manifest,
    const IcFieldManifest& field,
    double stored_value,
    const core::UnitSystem& target_units) {
  const double value = stored_value *
      icFieldConversionMultiplier(field, manifest, target_units);
  if (!std::isfinite(value)) {
    throw std::runtime_error(
        "non-finite converted Header value for " + field.dataset_path);
  }
  return value;
}

[[nodiscard]] double normalizePeriodicCoordinate(
    double value,
    double box_size,
    std::string_view axis_name) {
  if (!std::isfinite(value) || !(box_size > 0.0)) {
    throw std::runtime_error(
        "IC " + std::string(axis_name) + " coordinate is non-finite");
  }
  const double tolerance = 1.0e-12 * std::max(1.0, box_size);
  if (value < -tolerance || value > box_size + tolerance) {
    throw std::runtime_error(
        "IC " + std::string(axis_name) +
        " coordinate is outside the periodic box");
  }
  if (value < 0.0 || value >= box_size) {
    return 0.0;
  }
  return value;
}

void validateRecordScientificState(
    ParticleRecord& record,
    IcSpeciesPolicy policy,
    double box_size) {
  record.x = normalizePeriodicCoordinate(record.x, box_size, "x");
  record.y = normalizePeriodicCoordinate(record.y, box_size, "y");
  record.z = normalizePeriodicCoordinate(record.z, box_size, "z");
  if (!std::isfinite(record.vx) || !std::isfinite(record.vy) ||
      !std::isfinite(record.vz)) {
    throw std::runtime_error("IC velocity components must be finite");
  }
  if (record.id == 0U || !(record.mass > 0.0) ||
      !std::isfinite(record.mass)) {
    throw std::runtime_error(
        "IC particle IDs must be nonzero and masses finite/positive");
  }
  if (policy == IcSpeciesPolicy::kGas) {
    if (!std::isfinite(record.gas_density) || record.gas_density < 0.0 ||
        !std::isfinite(record.gas_internal_energy) ||
        record.gas_internal_energy < 0.0) {
      throw std::runtime_error(
          "gas density and internal energy must be finite and non-negative");
    }
  } else if (policy == IcSpeciesPolicy::kStar) {
    if (!std::isfinite(record.star_formation)) {
      throw std::runtime_error("stellar formation time must be finite");
    }
    if (record.star_formation < 0.0) {
      throw std::runtime_error(
          "negative stellar formation time identifies an AREPO wind particle; "
          "gadget_arepo_bridge_v1 does not silently reinterpret wind particles as ordinary stars");
    }
    if (
        !std::isfinite(record.star_birth_mass) ||
        !(record.star_birth_mass > 0.0) ||
        !std::isfinite(record.star_metallicity) ||
        record.star_metallicity < 0.0) {
      throw std::runtime_error("stellar IC sidecar values are invalid");
    }
  } else if (policy == IcSpeciesPolicy::kBlackHole) {
    if (!std::isfinite(record.bh_mass) || !(record.bh_mass > 0.0) ||
        !std::isfinite(record.bh_mdot) || record.bh_mdot < 0.0) {
      throw std::runtime_error("black-hole sidecar values must be physical");
    }
  } else if (policy == IcSpeciesPolicy::kTracer) {
    if (record.tracer_parent == 0U ||
        !std::isfinite(record.tracer_fraction) ||
        record.tracer_fraction < 0.0 || record.tracer_fraction > 1.0 ||
        !std::isfinite(record.tracer_last_host_mass) ||
        record.tracer_last_host_mass < 0.0 ||
        !std::isfinite(record.tracer_exchanged_mass)) {
      throw std::runtime_error("tracer IC sidecar values are invalid");
    }
  }
}

[[nodiscard]] std::vector<ParticleRecord> readRecordChunk(
    IcReaderSession& session,
    const Inspection& inspection,
    std::size_t file_index,
    std::size_t type_index,
    std::size_t start,
    std::size_t count,
    const core::SimulationConfig& config,
    const IcImportOptions& options,
    IcImportCounters& counters) {
  static_cast<void>(options);
  const IcManifest& manifest = inspection.manifest;
  const std::string prefix =
      "/PartType" + std::to_string(type_index) + "/";
  const core::UnitSystem target = core::makeUnitSystem(
      config.units.length_unit, config.units.mass_unit,
      config.units.velocity_unit);
  const double box_size = convertedHeaderFieldCode(
      manifest, requireField(manifest, file_index, "/Header/BoxSize"),
      manifest.box_size, target);

  std::vector<double> pos;
  std::vector<double> vel;
  std::vector<double> mass;
  std::vector<double> gas_u;
  std::vector<double> gas_rho;
  std::vector<double> star_time;
  std::vector<double> star_birth;
  std::vector<double> star_metal;
  std::vector<double> bh_mass;
  std::vector<double> bh_mdot;
  std::vector<double> tracer_fraction;
  std::vector<double> tracer_last;
  std::vector<double> tracer_exchange;
  std::vector<std::uint64_t> ids;
  std::vector<std::uint64_t> tracer_parent;
  std::vector<std::uint64_t> tracer_injection;
  std::vector<std::uint64_t> tracer_host64;

  const IcFieldManifest& position_field =
      requireField(manifest, file_index, prefix + "Coordinates");
  readChunkDouble(
      session, position_field, start, count, 3U, pos, counters);
  convertValues(pos, position_field, manifest, target);

  if (const auto* field =
          findField(manifest, file_index, prefix + "Velocities")) {
    readChunkDouble(session, *field, start, count, 3U, vel, counters);
    convertValues(vel, *field, manifest, target);
  } else {
    vel.assign(count * 3U, 0.0);
  }

  if (const auto* field =
          findField(manifest, file_index, prefix + "ParticleIDs")) {
    readChunkU64(session, *field, start, count, ids, counters);
  } else {
    ids.resize(count);
    const std::uint64_t base =
        precedingRecordCount(manifest, file_index, type_index) + start + 1U;
    for (std::size_t i = 0; i < count; ++i) {
      ids[i] = base + i;
    }
  }

  if (const auto* field =
          findField(manifest, file_index, prefix + "Masses")) {
    readChunkDouble(session, *field, start, count, 1U, mass, counters);
    convertValues(mass, *field, manifest, target);
  } else {
    mass.assign(count, manifest.mass_table[type_index]);
    convertValues(
        mass, requireField(manifest, file_index, "/Header/MassTable"),
        manifest, target);
  }

  const IcSpeciesPolicy policy = manifest.species_policy[type_index];
  if (policy == IcSpeciesPolicy::kGas) {
    if (const auto* field =
            findField(manifest, file_index, prefix + "InternalEnergy")) {
      readChunkDouble(
          session, *field, start, count, 1U, gas_u, counters);
      convertValues(gas_u, *field, manifest, target);
    } else {
      const auto& contract = requireMissingFieldContract(
          manifest, file_index, prefix + "InternalEnergy");
      gas_u.assign(count, resolveMissingScalarValue(contract, 0.0));
    }
    if (const auto* field =
            findField(manifest, file_index, prefix + "Density")) {
      readChunkDouble(
          session, *field, start, count, 1U, gas_rho, counters);
      convertValues(gas_rho, *field, manifest, target);
    } else {
      const auto& contract = requireMissingFieldContract(
          manifest, file_index, prefix + "Density");
      gas_rho.assign(count, resolveMissingScalarValue(contract, 0.0));
    }
  }

  if (policy == IcSpeciesPolicy::kStar) {
    if (const auto* field = findField(
            manifest, file_index, prefix + "StellarFormationTime")) {
      readChunkDouble(
          session, *field, start, count, 1U, star_time, counters);
    } else {
      const auto& contract = requireMissingFieldContract(
          manifest, file_index, prefix + "StellarFormationTime");
      star_time.assign(
          count, resolveMissingScalarValue(contract, manifest.scale_factor));
    }
    if (const auto* field =
            findField(manifest, file_index, prefix + "InitialMass")) {
      readChunkDouble(
          session, *field, start, count, 1U, star_birth, counters);
      convertValues(star_birth, *field, manifest, target);
    } else {
      const auto& contract = requireMissingFieldContract(
          manifest, file_index, prefix + "InitialMass");
      if (contract.policy == IcMissingFieldPolicy::kDialectDefinedDefault) {
        star_birth = mass;
      } else {
        star_birth.assign(
            count, resolveMissingScalarValue(contract, 0.0));
      }
    }
    if (const auto* field =
            findField(manifest, file_index, prefix + "Metallicity")) {
      readChunkDouble(
          session, *field, start, count, 1U, star_metal, counters);
    } else {
      const auto& contract = requireMissingFieldContract(
          manifest, file_index, prefix + "Metallicity");
      star_metal.assign(count, resolveMissingScalarValue(contract, 0.0));
    }
  }

  if (policy == IcSpeciesPolicy::kBlackHole) {
    const IcFieldManifest& field =
        requireField(manifest, file_index, prefix + "BH_Mass");
    readChunkDouble(
        session, field, start, count, 1U, bh_mass, counters);
    convertValues(bh_mass, field, manifest, target);
    if (const auto* mdot =
            findField(manifest, file_index, prefix + "BH_Mdot")) {
      readChunkDouble(
          session, *mdot, start, count, 1U, bh_mdot, counters);
      convertValues(bh_mdot, *mdot, manifest, target);
    } else {
      const auto& contract = requireMissingFieldContract(
          manifest, file_index, prefix + "BH_Mdot");
      bh_mdot.assign(count, resolveMissingScalarValue(contract, 0.0));
    }
  }

  if (policy == IcSpeciesPolicy::kTracer) {
    readChunkU64(
        session,
        requireField(manifest, file_index, prefix + "ParentParticleIDs"),
        start, count, tracer_parent, counters);
    readChunkU64(
        session,
        requireField(manifest, file_index, prefix + "InjectionStep"),
        start, count, tracer_injection, counters);
    tracer_host64.assign(count, kInvalidIndex);
    const auto& fraction_field =
        requireField(manifest, file_index, prefix + "MassFractionOfHost");
    readChunkDouble(
        session, fraction_field, start, count, 1U, tracer_fraction,
        counters);
    const auto& last_field =
        requireField(manifest, file_index, prefix + "LastHostMass");
    readChunkDouble(
        session, last_field, start, count, 1U, tracer_last, counters);
    convertValues(tracer_last, last_field, manifest, target);
    const auto& exchange_field = requireField(
        manifest, file_index, prefix + "CumulativeExchangedMass");
    readChunkDouble(
        session, exchange_field, start, count, 1U, tracer_exchange,
        counters);
    convertValues(tracer_exchange, exchange_field, manifest, target);
  }

  std::vector<ParticleRecord> records(count);
  for (std::size_t i = 0; i < count; ++i) {
    ParticleRecord& record = records[i];
    record.id = ids[i];
    record.species = speciesTag(policy);
    record.x = pos[i * 3U];
    record.y = pos[i * 3U + 1U];
    record.z = pos[i * 3U + 2U];
    record.vx = vel[i * 3U];
    record.vy = vel[i * 3U + 1U];
    record.vz = vel[i * 3U + 2U];
    record.mass = mass[i];
    if (policy == IcSpeciesPolicy::kGas) {
      record.gas_density = gas_rho[i];
      record.gas_internal_energy = gas_u[i];
    } else if (policy == IcSpeciesPolicy::kStar) {
      record.star_formation = star_time[i];
      record.star_birth_mass = star_birth[i];
      record.star_metallicity = star_metal[i];
    } else if (policy == IcSpeciesPolicy::kBlackHole) {
      record.bh_mass = bh_mass[i];
      record.bh_mdot = bh_mdot[i];
    } else if (policy == IcSpeciesPolicy::kTracer) {
      record.tracer_parent = tracer_parent[i];
      record.tracer_injection = tracer_injection[i];
      if (tracer_host64[i] > std::numeric_limits<std::uint32_t>::max()) {
        throw std::runtime_error(
            "tracer host cell index exceeds uint32 range");
      }
      record.tracer_host = static_cast<std::uint32_t>(tracer_host64[i]);
      record.tracer_fraction = tracer_fraction[i];
      record.tracer_last_host_mass = tracer_last[i];
      record.tracer_exchanged_mass = tracer_exchange[i];
    }
    validateRecordScientificState(record, policy, box_size);
  }

  std::uint64_t staging_bytes = vectorCapacityBytes(records);
  const auto account_vector = [&](const auto& values) {
    checkedCounterAdd(
        staging_bytes, vectorCapacityBytes(values), "peak_staging_bytes");
  };
  account_vector(pos);
  account_vector(vel);
  account_vector(mass);
  account_vector(gas_u);
  account_vector(gas_rho);
  account_vector(star_time);
  account_vector(star_birth);
  account_vector(star_metal);
  account_vector(bh_mass);
  account_vector(bh_mdot);
  account_vector(tracer_fraction);
  account_vector(tracer_last);
  account_vector(tracer_exchange);
  account_vector(ids);
  account_vector(tracer_parent);
  account_vector(tracer_injection);
  account_vector(tracer_host64);
  counters.peak_staging_bytes =
      std::max(counters.peak_staging_bytes, staging_bytes);
  checkedCounterAdd(counters.records_read, count, "records_read");
  checkedCounterAdd(counters.records_converted, count, "records_converted");
  counters.bytes_read = counters.metadata_bytes_read +
      counters.hash_bytes_read + counters.payload_bytes_read;
  return records;
}

void appendRecords(core::SimulationState& state,const std::vector<ParticleRecord>& records,std::uint32_t owner_rank){const std::size_t old_particles=state.particles.size();std::size_t gas_count=0,star_count=0,bh_count=0,tracer_count=0;for(const auto& r:records){if(r.species==static_cast<std::uint32_t>(core::ParticleSpecies::kGas))++gas_count;else if(r.species==static_cast<std::uint32_t>(core::ParticleSpecies::kStar))++star_count;else if(r.species==static_cast<std::uint32_t>(core::ParticleSpecies::kBlackHole))++bh_count;else if(r.species==static_cast<std::uint32_t>(core::ParticleSpecies::kTracer))++tracer_count;}
  const std::size_t old_gas=state.cells.size(),old_stars=state.star_particles.size(),old_bh=state.black_holes.size(),old_tracers=state.tracers.size();state.resizeParticles(old_particles+records.size());state.resizeCells(old_gas+gas_count);state.star_particles.resize(old_stars+star_count);state.black_holes.resize(old_bh+bh_count);state.tracers.resize(old_tracers+tracer_count);std::size_t gas_row=old_gas,star_row=old_stars,bh_row=old_bh,tracer_row=old_tracers;
  for(std::size_t i=0;i<records.size();++i){const auto& r=records[i];const std::size_t p=old_particles+i;state.particles.position_x_comoving[p]=r.x;state.particles.position_y_comoving[p]=r.y;state.particles.position_z_comoving[p]=r.z;state.particles.velocity_x_peculiar[p]=r.vx;state.particles.velocity_y_peculiar[p]=r.vy;state.particles.velocity_z_peculiar[p]=r.vz;state.particles.mass_code[p]=r.mass;state.particles.time_bin[p]=0U;state.particle_sidecar.particle_id[p]=r.id;state.particle_sidecar.species_tag[p]=r.species;state.particle_sidecar.owning_rank[p]=owner_rank;state.particle_sidecar.sfc_key[p]=0U;state.particle_sidecar.particle_flags[p]=0U;state.particle_sidecar.last_drift_time_code[p]=0.0;state.particle_sidecar.last_drift_scale_factor[p]=1.0;
    if(r.species==static_cast<std::uint32_t>(core::ParticleSpecies::kGas)){state.cells.center_x_comoving[gas_row]=r.x;state.cells.center_y_comoving[gas_row]=r.y;state.cells.center_z_comoving[gas_row]=r.z;state.cells.mass_code[gas_row]=r.mass;state.cells.time_bin[gas_row]=0U;state.cells.patch_index[gas_row]=0U;state.gas_cells.gas_cell_id[gas_row]=r.id;state.gas_cells.parent_particle_id[gas_row]=r.id;state.gas_cells.velocity_x_peculiar[gas_row]=r.vx;state.gas_cells.velocity_y_peculiar[gas_row]=r.vy;state.gas_cells.velocity_z_peculiar[gas_row]=r.vz;state.gas_cells.density_code[gas_row]=r.gas_density;state.gas_cells.internal_energy_code[gas_row]=r.gas_internal_energy;state.gas_cells.pressure_code[gas_row]=0.0;state.gas_cells.temperature_code[gas_row]=0.0;state.gas_cells.sound_speed_code[gas_row]=0.0;++gas_row;}
    else if(r.species==static_cast<std::uint32_t>(core::ParticleSpecies::kStar)){state.star_particles.particle_index[star_row]=static_cast<std::uint32_t>(p);state.star_particles.formation_scale_factor[star_row]=r.star_formation;state.star_particles.birth_mass_code[star_row]=r.star_birth_mass;state.star_particles.metallicity_mass_fraction[star_row]=r.star_metallicity;state.star_particles.stellar_age_years_last[star_row]=0.0;state.star_particles.stellar_returned_mass_cumulative_code[star_row]=0.0;state.star_particles.stellar_returned_metals_cumulative_code[star_row]=0.0;state.star_particles.stellar_feedback_energy_cumulative_erg[star_row]=0.0;for(std::size_t c=0;c<3U;++c){state.star_particles.stellar_returned_mass_channel_cumulative_code[c][star_row]=0.0;state.star_particles.stellar_returned_metals_channel_cumulative_code[c][star_row]=0.0;state.star_particles.stellar_feedback_energy_channel_cumulative_erg[c][star_row]=0.0;}++star_row;}
    else if(r.species==static_cast<std::uint32_t>(core::ParticleSpecies::kBlackHole)){state.black_holes.particle_index[bh_row]=static_cast<std::uint32_t>(p);state.black_holes.host_cell_index[bh_row]=kInvalidIndex;state.black_holes.subgrid_mass_code[bh_row]=r.bh_mass;state.black_holes.accretion_rate_code[bh_row]=r.bh_mdot;state.black_holes.feedback_energy_code[bh_row]=0.0;state.black_holes.eddington_ratio[bh_row]=0.0;state.black_holes.cumulative_accreted_mass_code[bh_row]=0.0;state.black_holes.cumulative_feedback_energy_code[bh_row]=0.0;state.black_holes.duty_cycle_active_time_code[bh_row]=0.0;state.black_holes.duty_cycle_total_time_code[bh_row]=0.0;++bh_row;}
    else if(r.species==static_cast<std::uint32_t>(core::ParticleSpecies::kTracer)){state.tracers.particle_index[tracer_row]=static_cast<std::uint32_t>(p);state.tracers.parent_particle_id[tracer_row]=r.tracer_parent;state.tracers.injection_step[tracer_row]=r.tracer_injection;state.tracers.host_cell_index[tracer_row]=r.tracer_host;state.tracers.mass_fraction_of_host[tracer_row]=r.tracer_fraction;state.tracers.last_host_mass_code[tracer_row]=r.tracer_last_host_mass;state.tracers.cumulative_exchanged_mass_code[tracer_row]=r.tracer_exchanged_mass;++tracer_row;}
  }
}

void finalizeImportedState(
    core::SimulationState& state,
    const IcManifest& manifest,
    const core::SimulationConfig& config) {
  state.metadata.run_name = config.output.run_name;
  state.metadata.scale_factor = manifest.scale_factor;
  state.species.count_by_species = {};
  for (const std::uint32_t species_tag :
       state.particle_sidecar.species_tag) {
    if (species_tag >= state.species.count_by_species.size()) {
      throw std::runtime_error(
          "IC import produced an out-of-range species tag");
    }
    ++state.species.count_by_species[species_tag];
  }
  state.rebuildSpeciesIndex();
  if (state.cells.size() > 0U) {
    state.refreshGasCellIdentityFromParticleOrder();
  }

  if (state.tracers.size() > 0U) {
    std::unordered_map<std::uint64_t, std::uint32_t> gas_row_by_parent;
    gas_row_by_parent.reserve(state.gas_cells.size());
    for (std::size_t row = 0U; row < state.gas_cells.size(); ++row) {
      const std::uint64_t parent_id =
          state.gas_cells.parent_particle_id[row];
      if (parent_id == 0U ||
          !gas_row_by_parent
               .emplace(parent_id, static_cast<std::uint32_t>(row))
               .second) {
        throw std::runtime_error(
            "IC gas parent identities must be nonzero and unique before "
            "tracer attachment");
      }
    }
    for (std::size_t row = 0U; row < state.tracers.size(); ++row) {
      const auto host = gas_row_by_parent.find(
          state.tracers.parent_particle_id[row]);
      if (host == gas_row_by_parent.end()) {
        throw std::runtime_error(
            "IC tracer parent must resolve to a gas cell on the same final "
            "owner rank");
      }
      state.tracers.host_cell_index[row] = host->second;
    }
  }

  if (!state.validateOwnershipInvariants()) {
    throw std::runtime_error(
        "IC import produced invalid species/sidecar/ownership invariants");
  }
}

[[nodiscard]] double convertedBoxSizeCode(
    const IcManifest& manifest,
    const core::SimulationConfig& config) {
  const core::UnitSystem target = core::makeUnitSystem(
      config.units.length_unit, config.units.mass_unit,
      config.units.velocity_unit);
  return convertedHeaderFieldCode(
      manifest, requireField(manifest, 0U, "/Header/BoxSize"),
      manifest.box_size, target);
}

void validateRuntimeCosmology(const IcManifest& manifest,const core::SimulationConfig& config){const double box_code=convertedBoxSizeCode(manifest,config);const core::UnitSystem target=core::makeUnitSystem(config.units.length_unit,config.units.mass_unit,config.units.velocity_unit);const core::UnitSystem mpc=core::makeUnitSystem("mpc","msun","km_s");const double box_mpc=box_code*target.length_si_per_code/mpc.length_si_per_code;if(!nearlyEqual(manifest.scale_factor,config.numerics.a_begin)||!nearlyEqual(manifest.omega_matter,config.cosmology.omega_matter)||!nearlyEqual(manifest.omega_lambda,config.cosmology.omega_lambda)||!nearlyEqual(manifest.hubble_param,config.cosmology.hubble_param)||!nearlyEqual(box_mpc,config.cosmology.box_size_x_mpc_comoving)||!nearlyEqual(box_mpc,config.cosmology.box_size_y_mpc_comoving)||!nearlyEqual(box_mpc,config.cosmology.box_size_z_mpc_comoving))throw std::runtime_error("IC manifest cosmology/BoxSize/start epoch does not match frozen runtime configuration");}

void validateSerialCountsAndIds(const core::SimulationState& state,const IcManifest& manifest){std::uint64_t expected=0;for(auto count:manifest.num_part_total){if(expected>std::numeric_limits<std::uint64_t>::max()-count)throw std::overflow_error("global particle count overflow");expected+=count;}if(state.particles.size()!=expected)throw std::runtime_error("IC import particle count mismatch");std::vector<std::uint64_t> ids(state.particle_sidecar.particle_id.begin(),state.particle_sidecar.particle_id.end());std::sort(ids.begin(),ids.end());if(std::adjacent_find(ids.begin(),ids.end())!=ids.end())throw std::runtime_error("IC import contains duplicate particle IDs");}

#endif  // COSMOSIM_ENABLE_HDF5

}  // namespace cosmosim::io::file_set_internal
