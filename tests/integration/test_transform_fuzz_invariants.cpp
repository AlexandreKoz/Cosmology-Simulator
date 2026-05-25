#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "cosmosim/core/provenance.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/core/time_integration.hpp"
#include "cosmosim/io/restart_checkpoint.hpp"

namespace {

using cosmosim::core::ParticleMigrationCommit;
using cosmosim::core::ParticleMigrationRecord;
using cosmosim::core::ParticleSpecies;
using cosmosim::core::SimulationState;

constexpr std::uint32_t speciesTag(const ParticleSpecies species) {
  return static_cast<std::uint32_t>(species);
}

[[noreturn]] void fail(const std::string& context, const std::string& message) {
  throw std::runtime_error(context + ": " + message);
}

void require(bool condition, const std::string& context, const std::string& message) {
  if (!condition) {
    fail(context, message);
  }
}

template <typename T>
void requireEqual(const T& actual, const T& expected, const std::string& context, const std::string& field) {
  if (!(actual == expected)) {
    std::ostringstream out;
    out << field << " mismatch";
    fail(context, out.str());
  }
}

struct ParticleOracleRecord {
  std::uint64_t sfc_key = 0;
  std::uint32_t species_tag = 0;
  std::uint32_t particle_flags = 0;
  std::uint32_t owning_rank = 0;
  std::uint8_t time_bin = 0;
  double position_x_comoving = 0.0;
  double position_y_comoving = 0.0;
  double position_z_comoving = 0.0;
  double velocity_x_peculiar = 0.0;
  double velocity_y_peculiar = 0.0;
  double velocity_z_peculiar = 0.0;
  double mass_code = 0.0;
  double gravity_softening_comoving = 0.0;
  bool has_gravity_softening_override = false;
};

struct StarOracleRecord {
  double formation_scale_factor = 0.0;
  double birth_mass_code = 0.0;
  double metallicity_mass_fraction = 0.0;
  double stellar_age_years_last = 0.0;
  double stellar_returned_mass_cumulative_code = 0.0;
  double stellar_returned_metals_cumulative_code = 0.0;
  double stellar_feedback_energy_cumulative_erg = 0.0;
  std::array<double, 3> stellar_returned_mass_channel_cumulative_code{};
  std::array<double, 3> stellar_returned_metals_channel_cumulative_code{};
  std::array<double, 3> stellar_feedback_energy_channel_cumulative_erg{};
};

struct BlackHoleOracleRecord {
  std::uint32_t host_cell_index = 0;
  double subgrid_mass_code = 0.0;
  double accretion_rate_code = 0.0;
  double feedback_energy_code = 0.0;
  double eddington_ratio = 0.0;
  double cumulative_accreted_mass_code = 0.0;
  double cumulative_feedback_energy_code = 0.0;
  double duty_cycle_active_time_code = 0.0;
  double duty_cycle_total_time_code = 0.0;
};

struct TracerOracleRecord {
  std::uint64_t parent_particle_id = 0;
  std::uint64_t injection_step = 0;
  std::uint32_t host_cell_index = 0;
  double mass_fraction_of_host = 0.0;
  double last_host_mass_code = 0.0;
  double cumulative_exchanged_mass_code = 0.0;
};

struct TransformOracle {
  std::unordered_map<std::uint64_t, ParticleOracleRecord> particles;
  std::unordered_map<std::uint64_t, StarOracleRecord> stars;
  std::unordered_map<std::uint64_t, BlackHoleOracleRecord> black_holes;
  std::unordered_map<std::uint64_t, TracerOracleRecord> tracers;
};

std::uint32_t findParticleIndexById(const SimulationState& state, const std::uint64_t particle_id) {
  for (std::uint32_t i = 0; i < state.particles.size(); ++i) {
    if (state.particle_sidecar.particle_id[i] == particle_id) {
      return i;
    }
  }
  throw std::runtime_error("particle ID not found in fuzz state");
}

void rebuildSpeciesLedger(SimulationState& state) {
  state.species.count_by_species.fill(0);
  for (const std::uint32_t tag : state.particle_sidecar.species_tag) {
    ++state.species.count_by_species[tag];
  }
  state.rebuildSpeciesIndex();
}

TransformOracle captureOracle(const SimulationState& state) {
  TransformOracle oracle;
  oracle.particles.reserve(state.particles.size());
  for (std::size_t i = 0; i < state.particles.size(); ++i) {
    oracle.particles.emplace(
        state.particle_sidecar.particle_id[i],
        ParticleOracleRecord{
            .sfc_key = state.particle_sidecar.sfc_key[i],
            .species_tag = state.particle_sidecar.species_tag[i],
            .particle_flags = state.particle_sidecar.particle_flags[i],
            .owning_rank = state.particle_sidecar.owning_rank[i],
            .time_bin = state.particles.time_bin[i],
            .position_x_comoving = state.particles.position_x_comoving[i],
            .position_y_comoving = state.particles.position_y_comoving[i],
            .position_z_comoving = state.particles.position_z_comoving[i],
            .velocity_x_peculiar = state.particles.velocity_x_peculiar[i],
            .velocity_y_peculiar = state.particles.velocity_y_peculiar[i],
            .velocity_z_peculiar = state.particles.velocity_z_peculiar[i],
            .mass_code = state.particles.mass_code[i],
            .gravity_softening_comoving = state.particle_sidecar.gravity_softening_comoving.empty()
                                               ? 0.0
                                               : state.particle_sidecar.gravity_softening_comoving[i],
            .has_gravity_softening_override = state.particle_sidecar.hasGravitySofteningOverride(i),
        });
  }

  for (std::size_t row = 0; row < state.star_particles.size(); ++row) {
    const std::uint64_t id = state.particle_sidecar.particle_id[state.star_particles.particle_index[row]];
    StarOracleRecord record{
        .formation_scale_factor = state.star_particles.formation_scale_factor[row],
        .birth_mass_code = state.star_particles.birth_mass_code[row],
        .metallicity_mass_fraction = state.star_particles.metallicity_mass_fraction[row],
        .stellar_age_years_last = state.star_particles.stellar_age_years_last[row],
        .stellar_returned_mass_cumulative_code = state.star_particles.stellar_returned_mass_cumulative_code[row],
        .stellar_returned_metals_cumulative_code = state.star_particles.stellar_returned_metals_cumulative_code[row],
        .stellar_feedback_energy_cumulative_erg = state.star_particles.stellar_feedback_energy_cumulative_erg[row],
    };
    for (std::size_t channel = 0; channel < 3; ++channel) {
      record.stellar_returned_mass_channel_cumulative_code[channel] =
          state.star_particles.stellar_returned_mass_channel_cumulative_code[channel][row];
      record.stellar_returned_metals_channel_cumulative_code[channel] =
          state.star_particles.stellar_returned_metals_channel_cumulative_code[channel][row];
      record.stellar_feedback_energy_channel_cumulative_erg[channel] =
          state.star_particles.stellar_feedback_energy_channel_cumulative_erg[channel][row];
    }
    oracle.stars.emplace(id, record);
  }

  for (std::size_t row = 0; row < state.black_holes.size(); ++row) {
    const std::uint64_t id = state.particle_sidecar.particle_id[state.black_holes.particle_index[row]];
    oracle.black_holes.emplace(
        id,
        BlackHoleOracleRecord{
            .host_cell_index = state.black_holes.host_cell_index[row],
            .subgrid_mass_code = state.black_holes.subgrid_mass_code[row],
            .accretion_rate_code = state.black_holes.accretion_rate_code[row],
            .feedback_energy_code = state.black_holes.feedback_energy_code[row],
            .eddington_ratio = state.black_holes.eddington_ratio[row],
            .cumulative_accreted_mass_code = state.black_holes.cumulative_accreted_mass_code[row],
            .cumulative_feedback_energy_code = state.black_holes.cumulative_feedback_energy_code[row],
            .duty_cycle_active_time_code = state.black_holes.duty_cycle_active_time_code[row],
            .duty_cycle_total_time_code = state.black_holes.duty_cycle_total_time_code[row],
        });
  }

  for (std::size_t row = 0; row < state.tracers.size(); ++row) {
    const std::uint64_t id = state.particle_sidecar.particle_id[state.tracers.particle_index[row]];
    oracle.tracers.emplace(
        id,
        TracerOracleRecord{
            .parent_particle_id = state.tracers.parent_particle_id[row],
            .injection_step = state.tracers.injection_step[row],
            .host_cell_index = state.tracers.host_cell_index[row],
            .mass_fraction_of_host = state.tracers.mass_fraction_of_host[row],
            .last_host_mass_code = state.tracers.last_host_mass_code[row],
            .cumulative_exchanged_mass_code = state.tracers.cumulative_exchanged_mass_code[row],
        });
  }
  return oracle;
}

void verifyOracle(const SimulationState& state, const TransformOracle& oracle, const std::string& context) {
  requireEqual(state.particles.size(), oracle.particles.size(), context, "particle count");
  require(state.validateOwnershipInvariants(), context, "ownership invariants failed");

  for (std::size_t i = 0; i < state.particles.size(); ++i) {
    const std::uint64_t id = state.particle_sidecar.particle_id[i];
    const auto it = oracle.particles.find(id);
    require(it != oracle.particles.end(), context, "unexpected particle ID " + std::to_string(id));
    const auto& record = it->second;
    requireEqual(state.particle_sidecar.sfc_key[i], record.sfc_key, context, "sfc_key for ID " + std::to_string(id));
    requireEqual(state.particle_sidecar.species_tag[i], record.species_tag, context, "species for ID " + std::to_string(id));
    requireEqual(state.particle_sidecar.particle_flags[i], record.particle_flags, context, "flags for ID " + std::to_string(id));
    requireEqual(state.particle_sidecar.owning_rank[i], record.owning_rank, context, "owner for ID " + std::to_string(id));
    requireEqual(state.particles.time_bin[i], record.time_bin, context, "time_bin for ID " + std::to_string(id));
    requireEqual(state.particles.position_x_comoving[i], record.position_x_comoving, context, "x for ID " + std::to_string(id));
    requireEqual(state.particles.position_y_comoving[i], record.position_y_comoving, context, "y for ID " + std::to_string(id));
    requireEqual(state.particles.position_z_comoving[i], record.position_z_comoving, context, "z for ID " + std::to_string(id));
    requireEqual(state.particles.velocity_x_peculiar[i], record.velocity_x_peculiar, context, "vx for ID " + std::to_string(id));
    requireEqual(state.particles.velocity_y_peculiar[i], record.velocity_y_peculiar, context, "vy for ID " + std::to_string(id));
    requireEqual(state.particles.velocity_z_peculiar[i], record.velocity_z_peculiar, context, "vz for ID " + std::to_string(id));
    requireEqual(state.particles.mass_code[i], record.mass_code, context, "mass for ID " + std::to_string(id));
    require(!state.particle_sidecar.gravity_softening_comoving.empty(), context, "softening lane absent");
    requireEqual(
        state.particle_sidecar.gravity_softening_comoving[i],
        record.gravity_softening_comoving,
        context,
        "softening for ID " + std::to_string(id));
    requireEqual(
        state.particle_sidecar.hasGravitySofteningOverride(i),
        record.has_gravity_softening_override,
        context,
        "softening override mask for ID " + std::to_string(id));
  }

  requireEqual(state.star_particles.size(), oracle.stars.size(), context, "star sidecar row count");
  for (std::size_t row = 0; row < state.star_particles.size(); ++row) {
    const std::uint64_t id = state.particle_sidecar.particle_id[state.star_particles.particle_index[row]];
    const auto it = oracle.stars.find(id);
    require(it != oracle.stars.end(), context, "unexpected star sidecar ID " + std::to_string(id));
    const auto& record = it->second;
    requireEqual(state.star_particles.formation_scale_factor[row], record.formation_scale_factor, context, "star formation scale");
    requireEqual(state.star_particles.birth_mass_code[row], record.birth_mass_code, context, "star birth mass");
    requireEqual(state.star_particles.metallicity_mass_fraction[row], record.metallicity_mass_fraction, context, "star metallicity");
    requireEqual(state.star_particles.stellar_age_years_last[row], record.stellar_age_years_last, context, "star age");
    requireEqual(
        state.star_particles.stellar_returned_mass_cumulative_code[row],
        record.stellar_returned_mass_cumulative_code,
        context,
        "star returned mass");
    requireEqual(
        state.star_particles.stellar_returned_metals_cumulative_code[row],
        record.stellar_returned_metals_cumulative_code,
        context,
        "star returned metals");
    requireEqual(
        state.star_particles.stellar_feedback_energy_cumulative_erg[row],
        record.stellar_feedback_energy_cumulative_erg,
        context,
        "star feedback energy");
    for (std::size_t channel = 0; channel < 3; ++channel) {
      requireEqual(
          state.star_particles.stellar_returned_mass_channel_cumulative_code[channel][row],
          record.stellar_returned_mass_channel_cumulative_code[channel],
          context,
          "star returned mass channel");
      requireEqual(
          state.star_particles.stellar_returned_metals_channel_cumulative_code[channel][row],
          record.stellar_returned_metals_channel_cumulative_code[channel],
          context,
          "star returned metals channel");
      requireEqual(
          state.star_particles.stellar_feedback_energy_channel_cumulative_erg[channel][row],
          record.stellar_feedback_energy_channel_cumulative_erg[channel],
          context,
          "star feedback channel");
    }
  }

  requireEqual(state.black_holes.size(), oracle.black_holes.size(), context, "black-hole sidecar row count");
  for (std::size_t row = 0; row < state.black_holes.size(); ++row) {
    const std::uint64_t id = state.particle_sidecar.particle_id[state.black_holes.particle_index[row]];
    const auto it = oracle.black_holes.find(id);
    require(it != oracle.black_holes.end(), context, "unexpected black-hole sidecar ID " + std::to_string(id));
    const auto& record = it->second;
    requireEqual(state.black_holes.host_cell_index[row], record.host_cell_index, context, "bh host cell");
    requireEqual(state.black_holes.subgrid_mass_code[row], record.subgrid_mass_code, context, "bh subgrid mass");
    requireEqual(state.black_holes.accretion_rate_code[row], record.accretion_rate_code, context, "bh accretion rate");
    requireEqual(state.black_holes.feedback_energy_code[row], record.feedback_energy_code, context, "bh feedback energy");
    requireEqual(state.black_holes.eddington_ratio[row], record.eddington_ratio, context, "bh eddington ratio");
    requireEqual(
        state.black_holes.cumulative_accreted_mass_code[row],
        record.cumulative_accreted_mass_code,
        context,
        "bh cumulative accreted mass");
    requireEqual(
        state.black_holes.cumulative_feedback_energy_code[row],
        record.cumulative_feedback_energy_code,
        context,
        "bh cumulative feedback");
    requireEqual(
        state.black_holes.duty_cycle_active_time_code[row],
        record.duty_cycle_active_time_code,
        context,
        "bh active time");
    requireEqual(
        state.black_holes.duty_cycle_total_time_code[row],
        record.duty_cycle_total_time_code,
        context,
        "bh total time");
  }

  requireEqual(state.tracers.size(), oracle.tracers.size(), context, "tracer sidecar row count");
  for (std::size_t row = 0; row < state.tracers.size(); ++row) {
    const std::uint64_t id = state.particle_sidecar.particle_id[state.tracers.particle_index[row]];
    const auto it = oracle.tracers.find(id);
    require(it != oracle.tracers.end(), context, "unexpected tracer sidecar ID " + std::to_string(id));
    const auto& record = it->second;
    requireEqual(state.tracers.parent_particle_id[row], record.parent_particle_id, context, "tracer parent");
    requireEqual(state.tracers.injection_step[row], record.injection_step, context, "tracer injection step");
    requireEqual(state.tracers.host_cell_index[row], record.host_cell_index, context, "tracer host cell");
    requireEqual(state.tracers.mass_fraction_of_host[row], record.mass_fraction_of_host, context, "tracer mass fraction");
    requireEqual(state.tracers.last_host_mass_code[row], record.last_host_mass_code, context, "tracer last host mass");
    requireEqual(
        state.tracers.cumulative_exchanged_mass_code[row],
        record.cumulative_exchanged_mass_code,
        context,
        "tracer cumulative exchanged mass");
  }
}

SimulationState generateFuzzState(const std::uint32_t seed, const std::size_t particle_count) {
  std::mt19937 rng(seed);
  SimulationState state;
  state.resizeParticles(particle_count);
  state.particle_sidecar.gravity_softening_comoving.resize(particle_count, 0.0);
  state.particle_sidecar.has_gravity_softening_override.resize(particle_count, 0U);

  const std::array<ParticleSpecies, 5> cycle{
      ParticleSpecies::kDarkMatter,
      ParticleSpecies::kGas,
      ParticleSpecies::kStar,
      ParticleSpecies::kBlackHole,
      ParticleSpecies::kTracer};
  std::vector<std::uint32_t> star_indices;
  std::vector<std::uint32_t> black_hole_indices;
  std::vector<std::uint32_t> tracer_indices;

  for (std::uint32_t i = 0; i < particle_count; ++i) {
    const auto species = cycle[(i * 3U + 1U) % cycle.size()];
    const std::uint64_t nonmonotonic_id = 900'000U + static_cast<std::uint64_t>((particle_count - i) * 37U) +
                                          static_cast<std::uint64_t>(rng() % 23U);
    state.particle_sidecar.particle_id[i] = nonmonotonic_id * 10U + static_cast<std::uint64_t>(i);
    state.particle_sidecar.sfc_key[i] = (static_cast<std::uint64_t>(rng()) << 16U) ^ static_cast<std::uint64_t>(i * 7919U);
    state.particle_sidecar.species_tag[i] = speciesTag(species);
    state.particle_sidecar.particle_flags[i] = 0xA500U + i * 17U;
    state.particle_sidecar.owning_rank[i] = 0;
    state.particles.time_bin[i] = static_cast<std::uint8_t>((rng() + i * 5U) % 8U);
    state.particles.position_x_comoving[i] = 0.01 * static_cast<double>(static_cast<int>(rng() % 200U) - 100);
    state.particles.position_y_comoving[i] = 0.02 * static_cast<double>(static_cast<int>(rng() % 200U) - 100);
    state.particles.position_z_comoving[i] = 0.03 * static_cast<double>(static_cast<int>(rng() % 200U) - 100);
    state.particles.velocity_x_peculiar[i] = -0.11 * static_cast<double>(i + 1U);
    state.particles.velocity_y_peculiar[i] = 0.13 * static_cast<double>((rng() % 31U) + i);
    state.particles.velocity_z_peculiar[i] = -0.17 * static_cast<double>((rng() % 29U) + 2U * i);
    state.particles.mass_code[i] = 1.0 + 0.25 * static_cast<double>(i) + 0.001 * static_cast<double>(rng() % 17U);
    state.particle_sidecar.gravity_softening_comoving[i] = 0.002 + 0.0001 * static_cast<double>((i * 19U) % 41U);
    if (((i + seed) % 5U) == 0U || ((i * 7U + seed) % 11U) == 3U) {
      state.particle_sidecar.setGravitySofteningOverride(i, 0.02 + 0.0003 * static_cast<double>(i));
    }

    if (species == ParticleSpecies::kStar) {
      star_indices.push_back(i);
    } else if (species == ParticleSpecies::kBlackHole) {
      black_hole_indices.push_back(i);
    } else if (species == ParticleSpecies::kTracer) {
      tracer_indices.push_back(i);
    }
  }

  state.star_particles.resize(star_indices.size());
  for (std::size_t row = 0; row < star_indices.size(); ++row) {
    state.star_particles.particle_index[row] = star_indices[row];
    state.star_particles.formation_scale_factor[row] = 0.35 + 0.01 * static_cast<double>(row);
    state.star_particles.birth_mass_code[row] = 2.0 + static_cast<double>(row);
    state.star_particles.metallicity_mass_fraction[row] = 0.001 * static_cast<double>(row + 1U);
    state.star_particles.stellar_age_years_last[row] = 1.0e6 * static_cast<double>(row + 1U);
    state.star_particles.stellar_returned_mass_cumulative_code[row] = 0.05 * static_cast<double>(row + 1U);
    state.star_particles.stellar_returned_metals_cumulative_code[row] = 0.005 * static_cast<double>(row + 1U);
    state.star_particles.stellar_feedback_energy_cumulative_erg[row] = 1.0e49 * static_cast<double>(row + 1U);
    for (std::size_t channel = 0; channel < 3; ++channel) {
      state.star_particles.stellar_returned_mass_channel_cumulative_code[channel][row] =
          0.1 * static_cast<double>((channel + 1U) * (row + 1U));
      state.star_particles.stellar_returned_metals_channel_cumulative_code[channel][row] =
          0.01 * static_cast<double>((channel + 1U) * (row + 1U));
      state.star_particles.stellar_feedback_energy_channel_cumulative_erg[channel][row] =
          1.0e48 * static_cast<double>((channel + 1U) * (row + 1U));
    }
  }

  state.black_holes.resize(black_hole_indices.size());
  for (std::size_t row = 0; row < black_hole_indices.size(); ++row) {
    state.black_holes.particle_index[row] = black_hole_indices[row];
    state.black_holes.host_cell_index[row] = 0;
    state.black_holes.subgrid_mass_code[row] = 10.0 + static_cast<double>(row);
    state.black_holes.accretion_rate_code[row] = 0.01 * static_cast<double>(row + 1U);
    state.black_holes.feedback_energy_code[row] = 0.02 * static_cast<double>(row + 1U);
    state.black_holes.eddington_ratio[row] = 0.1 * static_cast<double>(row + 1U);
    state.black_holes.cumulative_accreted_mass_code[row] = 0.3 * static_cast<double>(row + 1U);
    state.black_holes.cumulative_feedback_energy_code[row] = 0.4 * static_cast<double>(row + 1U);
    state.black_holes.duty_cycle_active_time_code[row] = 0.5 * static_cast<double>(row + 1U);
    state.black_holes.duty_cycle_total_time_code[row] = 0.6 * static_cast<double>(row + 1U);
  }

  state.tracers.resize(tracer_indices.size());
  for (std::size_t row = 0; row < tracer_indices.size(); ++row) {
    state.tracers.particle_index[row] = tracer_indices[row];
    state.tracers.parent_particle_id[row] = state.particle_sidecar.particle_id[(tracer_indices[row] + 1U) % particle_count];
    state.tracers.injection_step[row] = 1000U + row;
    state.tracers.host_cell_index[row] = 0;
    state.tracers.mass_fraction_of_host[row] = 0.01 * static_cast<double>(row + 1U);
    state.tracers.last_host_mass_code[row] = 1.5 + static_cast<double>(row);
    state.tracers.cumulative_exchanged_mass_code[row] = 0.2 * static_cast<double>(row + 1U);
  }

  rebuildSpeciesLedger(state);
  return state;
}

void compactSidecarsAfterFrontResize(SimulationState& state) {
  auto compact_star = [&]() {
    std::size_t write = 0;
    for (std::size_t row = 0; row < state.star_particles.size(); ++row) {
      if (state.star_particles.particle_index[row] >= state.particles.size()) {
        continue;
      }
      state.star_particles.particle_index[write] = state.star_particles.particle_index[row];
      state.star_particles.formation_scale_factor[write] = state.star_particles.formation_scale_factor[row];
      state.star_particles.birth_mass_code[write] = state.star_particles.birth_mass_code[row];
      state.star_particles.metallicity_mass_fraction[write] = state.star_particles.metallicity_mass_fraction[row];
      state.star_particles.stellar_age_years_last[write] = state.star_particles.stellar_age_years_last[row];
      state.star_particles.stellar_returned_mass_cumulative_code[write] =
          state.star_particles.stellar_returned_mass_cumulative_code[row];
      state.star_particles.stellar_returned_metals_cumulative_code[write] =
          state.star_particles.stellar_returned_metals_cumulative_code[row];
      state.star_particles.stellar_feedback_energy_cumulative_erg[write] =
          state.star_particles.stellar_feedback_energy_cumulative_erg[row];
      for (std::size_t channel = 0; channel < 3; ++channel) {
        state.star_particles.stellar_returned_mass_channel_cumulative_code[channel][write] =
            state.star_particles.stellar_returned_mass_channel_cumulative_code[channel][row];
        state.star_particles.stellar_returned_metals_channel_cumulative_code[channel][write] =
            state.star_particles.stellar_returned_metals_channel_cumulative_code[channel][row];
        state.star_particles.stellar_feedback_energy_channel_cumulative_erg[channel][write] =
            state.star_particles.stellar_feedback_energy_channel_cumulative_erg[channel][row];
      }
      ++write;
    }
    state.star_particles.resize(write);
  };

  auto compact_black_holes = [&]() {
    std::size_t write = 0;
    for (std::size_t row = 0; row < state.black_holes.size(); ++row) {
      if (state.black_holes.particle_index[row] >= state.particles.size()) {
        continue;
      }
      state.black_holes.particle_index[write] = state.black_holes.particle_index[row];
      state.black_holes.host_cell_index[write] = state.black_holes.host_cell_index[row];
      state.black_holes.subgrid_mass_code[write] = state.black_holes.subgrid_mass_code[row];
      state.black_holes.accretion_rate_code[write] = state.black_holes.accretion_rate_code[row];
      state.black_holes.feedback_energy_code[write] = state.black_holes.feedback_energy_code[row];
      state.black_holes.eddington_ratio[write] = state.black_holes.eddington_ratio[row];
      state.black_holes.cumulative_accreted_mass_code[write] = state.black_holes.cumulative_accreted_mass_code[row];
      state.black_holes.cumulative_feedback_energy_code[write] = state.black_holes.cumulative_feedback_energy_code[row];
      state.black_holes.duty_cycle_active_time_code[write] = state.black_holes.duty_cycle_active_time_code[row];
      state.black_holes.duty_cycle_total_time_code[write] = state.black_holes.duty_cycle_total_time_code[row];
      ++write;
    }
    state.black_holes.resize(write);
  };

  auto compact_tracers = [&]() {
    std::size_t write = 0;
    for (std::size_t row = 0; row < state.tracers.size(); ++row) {
      if (state.tracers.particle_index[row] >= state.particles.size()) {
        continue;
      }
      state.tracers.particle_index[write] = state.tracers.particle_index[row];
      state.tracers.parent_particle_id[write] = state.tracers.parent_particle_id[row];
      state.tracers.injection_step[write] = state.tracers.injection_step[row];
      state.tracers.host_cell_index[write] = state.tracers.host_cell_index[row];
      state.tracers.mass_fraction_of_host[write] = state.tracers.mass_fraction_of_host[row];
      state.tracers.last_host_mass_code[write] = state.tracers.last_host_mass_code[row];
      state.tracers.cumulative_exchanged_mass_code[write] = state.tracers.cumulative_exchanged_mass_code[row];
      ++write;
    }
    state.tracers.resize(write);
  };

  compact_star();
  compact_black_holes();
  compact_tracers();
  rebuildSpeciesLedger(state);
}

void reorderByPermutation(SimulationState& state, const std::vector<std::uint32_t>& new_to_old) {
  cosmosim::core::ParticleReorderMap map;
  map.new_to_old_index = new_to_old;
  map.old_to_new_index.resize(new_to_old.size());
  for (std::uint32_t new_index = 0; new_index < new_to_old.size(); ++new_index) {
    map.old_to_new_index[new_to_old[new_index]] = new_index;
  }
  cosmosim::core::reorderParticles(state, map);
}

void eraseMissingParticlesFromOracle(TransformOracle& oracle, const std::unordered_set<std::uint64_t>& kept_ids) {
  for (auto it = oracle.particles.begin(); it != oracle.particles.end();) {
    if (kept_ids.count(it->first) == 0U) {
      oracle.stars.erase(it->first);
      oracle.black_holes.erase(it->first);
      oracle.tracers.erase(it->first);
      it = oracle.particles.erase(it);
    } else {
      ++it;
    }
  }
}

void compactByDeterministicMask(SimulationState& state, TransformOracle& oracle, const std::uint32_t seed) {
  std::vector<std::uint32_t> kept;
  std::vector<std::uint32_t> dropped;
  kept.reserve(state.particles.size());
  dropped.reserve(state.particles.size());
  for (std::uint32_t i = 0; i < state.particles.size(); ++i) {
    const bool keep = ((state.particles.time_bin[i] + i + seed) % 4U) != 0U;
    (keep ? kept : dropped).push_back(i);
  }
  if (kept.size() < 8U) {
    kept.insert(kept.end(), dropped.begin(), dropped.begin() + static_cast<std::ptrdiff_t>(8U - kept.size()));
    dropped.erase(dropped.begin(), dropped.begin() + static_cast<std::ptrdiff_t>(8U - kept.size()));
  }

  std::vector<std::uint32_t> new_to_old = kept;
  new_to_old.insert(new_to_old.end(), dropped.begin(), dropped.end());
  reorderByPermutation(state, new_to_old);

  std::unordered_set<std::uint64_t> kept_ids;
  kept_ids.reserve(kept.size());
  for (std::size_t i = 0; i < kept.size(); ++i) {
    kept_ids.insert(state.particle_sidecar.particle_id[i]);
  }
  eraseMissingParticlesFromOracle(oracle, kept_ids);

  state.resizeParticles(kept.size());
  compactSidecarsAfterFrontResize(state);
}

std::vector<std::uint32_t> firstIndicesBySpecies(const SimulationState& state, const ParticleSpecies species, const std::size_t limit) {
  std::vector<std::uint32_t> indices;
  for (std::uint32_t i = 0; i < state.particles.size() && indices.size() < limit; ++i) {
    if (state.particle_sidecar.species_tag[i] == speciesTag(species)) {
      indices.push_back(i);
    }
  }
  return indices;
}

void migrateRoundTripSparseSidecars(SimulationState& state) {
  std::vector<std::uint32_t> outbound;
  for (const ParticleSpecies species : {ParticleSpecies::kStar, ParticleSpecies::kBlackHole, ParticleSpecies::kTracer}) {
    const auto found = firstIndicesBySpecies(state, species, 1U);
    outbound.insert(outbound.end(), found.begin(), found.end());
  }
  if (outbound.empty()) {
    outbound.push_back(0U);
  }
  std::sort(outbound.begin(), outbound.end());
  outbound.erase(std::unique(outbound.begin(), outbound.end()), outbound.end());

  std::vector<ParticleMigrationRecord> records = state.packParticleMigrationRecords(outbound);
  ParticleMigrationCommit commit;
  commit.world_rank = 0;
  commit.outbound_local_indices = outbound;
  commit.inbound_records = records;
  state.commitParticleMigration(commit);
}

void rebuildGasStateByParticleId(SimulationState& state, const std::string& context) {
  const auto gas_indices = state.particle_species_index.globalIndices(ParticleSpecies::kGas);
  state.resizeCells(gas_indices.size());
  state.refreshGasCellIdentityFromParticleOrder();
  for (std::size_t cell = 0; cell < state.cells.size(); ++cell) {
    const std::uint64_t parent_id = state.gas_cells.parent_particle_id[cell];
    state.cells.center_x_comoving[cell] = 0.001 * static_cast<double>(parent_id % 997U);
    state.cells.center_y_comoving[cell] = 0.002 * static_cast<double>(parent_id % 991U);
    state.cells.center_z_comoving[cell] = 0.003 * static_cast<double>(parent_id % 983U);
    state.cells.mass_code[cell] = 10.0 + static_cast<double>(parent_id % 17U);
    state.cells.time_bin[cell] = state.particles.time_bin[gas_indices[cell]];
    state.cells.patch_index[cell] = 0;
    state.gas_cells.density_code[cell] = 100.0 + static_cast<double>(parent_id % 31U);
    state.gas_cells.pressure_code[cell] = 200.0 + static_cast<double>(parent_id % 29U);
    state.gas_cells.internal_energy_code[cell] = 300.0 + static_cast<double>(parent_id % 23U);
    state.gas_cells.temperature_code[cell] = 400.0 + static_cast<double>(parent_id % 19U);
    state.gas_cells.sound_speed_code[cell] = 500.0 + static_cast<double>(parent_id % 13U);
  }
  require(state.gasCellIdentityMatchesParticleOrder(), context, "gas identity does not match particle order");
  require(state.validateOwnershipInvariants(), context, "ownership failed after gas rebuild");
}

void verifyRestartRoundTripIfEnabled(SimulationState state, const TransformOracle& oracle, const std::string& context) {
#if COSMOSIM_ENABLE_HDF5
  cosmosim::core::IntegratorState integrator_state;
  std::uint8_t scheduler_max_bin = 0;
  for (const auto bin : state.particles.time_bin) {
    scheduler_max_bin = std::max(scheduler_max_bin, bin);
  }
  cosmosim::core::HierarchicalTimeBinScheduler scheduler(scheduler_max_bin);
  scheduler.reset(static_cast<std::uint32_t>(state.particles.size()), 0, 0);
  for (std::uint32_t i = 0; i < static_cast<std::uint32_t>(state.particles.size()); ++i) {
    scheduler.setElementBin(i, state.particles.time_bin[i], scheduler.currentTick());
  }
  cosmosim::core::syncTimeBinMirrorsFromScheduler(scheduler, state);

  cosmosim::io::RestartWritePayload payload;
  payload.persistent_state.simulation_state = &state;
  payload.integrator_state = &integrator_state;
  payload.scheduler = &scheduler;
  payload.normalized_config_text = "schema_version = 1\n[mode]\nmode = zoom_in\n[fuzz]\nseed = fixed\n";
  payload.normalized_config_hash_hex = cosmosim::core::stableConfigHashHex(payload.normalized_config_text);
  state.metadata.run_name = "transform_fuzz_invariants";
  state.metadata.normalized_config_hash_hex = payload.normalized_config_hash_hex;
  payload.provenance = cosmosim::core::makeProvenanceRecord(payload.normalized_config_hash_hex, "transform-fuzz");
  payload.distributed_gravity_state.schema_version = 2;
  payload.distributed_gravity_state.world_size = 1;
  payload.distributed_gravity_state.pm_grid_nx = 4;
  payload.distributed_gravity_state.pm_grid_ny = 4;
  payload.distributed_gravity_state.pm_grid_nz = 4;
  payload.distributed_gravity_state.owning_rank_by_item.assign(state.particles.size(), 0U);
  payload.distributed_gravity_state.pm_slab_begin_x_by_rank = {0U};
  payload.distributed_gravity_state.pm_slab_end_x_by_rank = {4U};

  std::string path_token = context;
  for (char& ch : path_token) {
    if (ch == '/' || ch == '=' || ch == ':') {
      ch = '_';
    }
  }
  const auto path = std::filesystem::temp_directory_path() / ("cosmosim_transform_fuzz_" + path_token + ".hdf5");
  cosmosim::io::writeRestartCheckpointHdf5(path, payload);
  const auto restored = cosmosim::io::readRestartCheckpointHdf5(path);
  verifyOracle(restored.state, oracle, context + "/restart");
  require(restored.state.gasCellIdentityMatchesParticleOrder(), context, "restart gas identity mismatch");
  std::filesystem::remove(path);
#else
  (void)state;
  (void)oracle;
  (void)context;
#endif
}

void requireStaleScatterThrows(
    const cosmosim::core::GravityParticleKernelView& view,
    SimulationState& state,
    const std::string& context) {
  bool threw = false;
  try {
    cosmosim::core::scatterGravityParticleKernelView(view, state);
  } catch (const std::runtime_error&) {
    threw = true;
  }
  require(threw, context, "stale active particle scatter did not fail");
}

void runTransformFuzzHarness(const std::size_t trials, const std::uint32_t base_seed) {
  for (std::size_t trial = 0; trial < trials; ++trial) {
    const std::uint32_t seed = base_seed + static_cast<std::uint32_t>(trial * 97U);
    const std::string context = "seed=" + std::to_string(seed) + "/trial=" + std::to_string(trial);
    SimulationState state = generateFuzzState(seed, 18U + (trial % 5U));
    TransformOracle oracle = captureOracle(state);
    verifyOracle(state, oracle, context + "/initial");

    cosmosim::core::TransientStepWorkspace workspace;
    const std::vector<std::uint32_t> active_particles{0U, 2U, static_cast<std::uint32_t>(state.particles.size() - 1U)};
    auto stale_view = cosmosim::core::buildGravityParticleKernelView(state, active_particles, workspace);

    std::vector<std::uint32_t> permutation(state.particles.size());
    std::iota(permutation.begin(), permutation.end(), 0U);
    std::mt19937 rng(seed ^ 0xC05A2026U);
    std::shuffle(permutation.begin(), permutation.end(), rng);
    reorderByPermutation(state, permutation);
    verifyOracle(state, oracle, context + "/after-random-reorder");
    requireStaleScatterThrows(stale_view, state, context + "/stale-after-reorder");

    stale_view = cosmosim::core::buildGravityParticleKernelView(state, active_particles, workspace);
    compactByDeterministicMask(state, oracle, seed);
    verifyOracle(state, oracle, context + "/after-compact");
    requireStaleScatterThrows(stale_view, state, context + "/stale-after-compact");

    const auto sfc_reorder = cosmosim::core::buildParticleReorderMap(state, cosmosim::core::ParticleReorderMode::kBySfcKey);
    cosmosim::core::reorderParticles(state, sfc_reorder);
    verifyOracle(state, oracle, context + "/after-sfc-reorder");

    migrateRoundTripSparseSidecars(state);
    verifyOracle(state, oracle, context + "/after-migration-commit");

    rebuildGasStateByParticleId(state, context + "/gas-rebuild");
    verifyOracle(state, oracle, context + "/after-gas-rebuild");
    verifyRestartRoundTripIfEnabled(state, oracle, context);
  }
}

}  // namespace

int main(int argc, char** argv) {
  const bool long_validation = argc > 1 && std::string(argv[1]) == "--long";
  runTransformFuzzHarness(long_validation ? 32U : 6U, 0xC05A'2026U);
  return EXIT_SUCCESS;
}
