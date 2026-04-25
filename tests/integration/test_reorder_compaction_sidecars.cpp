#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#include "cosmosim/core/simulation_state.hpp"

namespace {

using cosmosim::core::ParticleSpecies;
using cosmosim::core::SimulationState;

struct ParticleIdentityRecord {
  std::uint32_t species_tag = 0;
  std::uint8_t time_bin = 0;
  double position_x = 0.0;
  double mass = 0.0;
  double softening = 0.0;
};

[[nodiscard]] std::unordered_map<std::uint64_t, ParticleIdentityRecord> captureParticleIdentity(const SimulationState& state) {
  std::unordered_map<std::uint64_t, ParticleIdentityRecord> map;
  for (std::size_t i = 0; i < state.particles.size(); ++i) {
    map[state.particle_sidecar.particle_id[i]] = ParticleIdentityRecord{
        .species_tag = state.particle_sidecar.species_tag[i],
        .time_bin = state.particles.time_bin[i],
        .position_x = state.particles.position_x_comoving[i],
        .mass = state.particles.mass_code[i],
        .softening = state.particle_sidecar.gravity_softening_comoving.empty()
                         ? 0.0
                         : state.particle_sidecar.gravity_softening_comoving[i],
    };
  }
  return map;
}

void rebuildSpeciesLedger(SimulationState& state) {
  state.species.count_by_species.fill(0);
  for (const auto tag : state.particle_sidecar.species_tag) {
    ++state.species.count_by_species[tag];
  }
  state.rebuildSpeciesIndex();
}

void seedSyntheticState(SimulationState& state) {
  state.resizeParticles(7);
  state.particle_sidecar.gravity_softening_comoving.resize(7);

  const std::array<std::uint64_t, 7> ids{701, 702, 703, 704, 705, 706, 707};
  const std::array<std::uint32_t, 7> species{
      static_cast<std::uint32_t>(ParticleSpecies::kDarkMatter),
      static_cast<std::uint32_t>(ParticleSpecies::kGas),
      static_cast<std::uint32_t>(ParticleSpecies::kStar),
      static_cast<std::uint32_t>(ParticleSpecies::kBlackHole),
      static_cast<std::uint32_t>(ParticleSpecies::kTracer),
      static_cast<std::uint32_t>(ParticleSpecies::kGas),
      static_cast<std::uint32_t>(ParticleSpecies::kStar)};

  for (std::size_t i = 0; i < ids.size(); ++i) {
    state.particle_sidecar.particle_id[i] = ids[i];
    state.particle_sidecar.species_tag[i] = species[i];
    state.particle_sidecar.particle_flags[i] = static_cast<std::uint32_t>(10 + i);
    state.particle_sidecar.owning_rank[i] = 0;
    state.particle_sidecar.sfc_key[i] = static_cast<std::uint64_t>(900 - 17 * i);

    state.particles.position_x_comoving[i] = static_cast<double>(i) + 0.25;
    state.particles.position_y_comoving[i] = static_cast<double>(i) + 0.5;
    state.particles.position_z_comoving[i] = static_cast<double>(i) + 0.75;
    state.particles.velocity_x_peculiar[i] = -static_cast<double>(i);
    state.particles.velocity_y_peculiar[i] = static_cast<double>(i) * 2.0;
    state.particles.velocity_z_peculiar[i] = static_cast<double>(i) * 3.0;
    state.particles.mass_code[i] = 100.0 + static_cast<double>(i);
    state.particles.time_bin[i] = static_cast<std::uint8_t>((i * 3U) % 4U);
    state.particle_sidecar.gravity_softening_comoving[i] = 0.005 * static_cast<double>(i + 1U);
  }

  state.star_particles.resize(2);
  state.star_particles.particle_index = {2U, 6U};
  for (std::size_t row = 0; row < state.star_particles.size(); ++row) {
    state.star_particles.formation_scale_factor[row] = 0.7 + 0.1 * static_cast<double>(row);
    state.star_particles.birth_mass_code[row] = 50.0 + static_cast<double>(row);
    state.star_particles.metallicity_mass_fraction[row] = 0.01 + 0.02 * static_cast<double>(row);
    state.star_particles.stellar_age_years_last[row] = 1.0e6 * (1.0 + static_cast<double>(row));
    state.star_particles.stellar_returned_mass_cumulative_code[row] = 0.2 + static_cast<double>(row);
    state.star_particles.stellar_returned_metals_cumulative_code[row] = 0.3 + static_cast<double>(row);
    state.star_particles.stellar_feedback_energy_cumulative_erg[row] = 1.0e50 * (2.0 + static_cast<double>(row));
    for (std::size_t channel = 0; channel < 3; ++channel) {
      state.star_particles.stellar_returned_mass_channel_cumulative_code[channel][row] =
          10.0 * static_cast<double>(channel + 1) + static_cast<double>(row);
      state.star_particles.stellar_returned_metals_channel_cumulative_code[channel][row] =
          20.0 * static_cast<double>(channel + 1) + static_cast<double>(row);
      state.star_particles.stellar_feedback_energy_channel_cumulative_erg[channel][row] =
          30.0 * static_cast<double>(channel + 1) + static_cast<double>(row);
    }
  }

  state.black_holes.resize(1);
  state.black_holes.particle_index[0] = 3;
  state.black_holes.host_cell_index[0] = 42;
  state.black_holes.subgrid_mass_code[0] = 9.0;
  state.black_holes.accretion_rate_code[0] = 8.0;
  state.black_holes.feedback_energy_code[0] = 7.0;
  state.black_holes.eddington_ratio[0] = 0.6;
  state.black_holes.cumulative_accreted_mass_code[0] = 5.0;
  state.black_holes.cumulative_feedback_energy_code[0] = 4.0;
  state.black_holes.duty_cycle_active_time_code[0] = 3.0;
  state.black_holes.duty_cycle_total_time_code[0] = 2.0;

  state.tracers.resize(1);
  state.tracers.particle_index[0] = 4;
  state.tracers.parent_particle_id[0] = 702;
  state.tracers.injection_step[0] = 99;
  state.tracers.host_cell_index[0] = 7;
  state.tracers.mass_fraction_of_host[0] = 0.123;
  state.tracers.last_host_mass_code[0] = 11.0;
  state.tracers.cumulative_exchanged_mass_code[0] = 1.25;

  rebuildSpeciesLedger(state);
  cosmosim::core::debugAssertNoStaleParticleIndices(state);
  assert(state.validateOwnershipInvariants());
}

void compactFront(SimulationState& state, std::size_t keep_count) {
  state.resizeParticles(keep_count);

  auto compact_star = [&]() {
    std::size_t write = 0;
    for (std::size_t row = 0; row < state.star_particles.size(); ++row) {
      if (state.star_particles.particle_index[row] >= keep_count) {
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
      if (state.black_holes.particle_index[row] >= keep_count) {
        continue;
      }
      state.black_holes.particle_index[write] = state.black_holes.particle_index[row];
      state.black_holes.host_cell_index[write] = state.black_holes.host_cell_index[row];
      state.black_holes.subgrid_mass_code[write] = state.black_holes.subgrid_mass_code[row];
      state.black_holes.accretion_rate_code[write] = state.black_holes.accretion_rate_code[row];
      state.black_holes.feedback_energy_code[write] = state.black_holes.feedback_energy_code[row];
      state.black_holes.eddington_ratio[write] = state.black_holes.eddington_ratio[row];
      state.black_holes.cumulative_accreted_mass_code[write] = state.black_holes.cumulative_accreted_mass_code[row];
      state.black_holes.cumulative_feedback_energy_code[write] =
          state.black_holes.cumulative_feedback_energy_code[row];
      state.black_holes.duty_cycle_active_time_code[write] = state.black_holes.duty_cycle_active_time_code[row];
      state.black_holes.duty_cycle_total_time_code[write] = state.black_holes.duty_cycle_total_time_code[row];
      ++write;
    }
    state.black_holes.resize(write);
  };

  auto compact_tracers = [&]() {
    std::size_t write = 0;
    for (std::size_t row = 0; row < state.tracers.size(); ++row) {
      if (state.tracers.particle_index[row] >= keep_count) {
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
  cosmosim::core::debugAssertNoStaleParticleIndices(state);
}

void assertParticleIdentityInvariant(
    const SimulationState& state,
    const std::unordered_map<std::uint64_t, ParticleIdentityRecord>& expected) {
  for (std::size_t i = 0; i < state.particles.size(); ++i) {
    const auto id = state.particle_sidecar.particle_id[i];
    const auto it = expected.find(id);
    if (it == expected.end()) {
      continue;
    }
    assert(state.particle_sidecar.species_tag[i] == it->second.species_tag);
    assert(state.particles.time_bin[i] == it->second.time_bin);
    assert(state.particles.position_x_comoving[i] == it->second.position_x);
    assert(state.particles.mass_code[i] == it->second.mass);
    assert(!state.particle_sidecar.gravity_softening_comoving.empty());
    assert(state.particle_sidecar.gravity_softening_comoving[i] == it->second.softening);
  }

  cosmosim::core::debugAssertNoStaleParticleIndices(state);
  assert(state.validateOwnershipInvariants());
}

void test_particle_resize_identity_invariants() {
  SimulationState state;
  seedSyntheticState(state);
  const auto baseline = captureParticleIdentity(state);

  state.resizeParticles(10);
  for (std::size_t i = 7; i < 10; ++i) {
    state.particle_sidecar.particle_id[i] = 900 + i;
    state.particle_sidecar.species_tag[i] = static_cast<std::uint32_t>(ParticleSpecies::kGas);
    state.particles.time_bin[i] = static_cast<std::uint8_t>(i % 4U);
    state.particles.position_x_comoving[i] = 50.0 + static_cast<double>(i);
    state.particles.mass_code[i] = 500.0 + static_cast<double>(i);
    state.particle_sidecar.gravity_softening_comoving[i] = 0.5 + static_cast<double>(i);
  }
  rebuildSpeciesLedger(state);

  assertParticleIdentityInvariant(state, baseline);

  compactFront(state, 5);
  assert(state.particles.size() == 5);
  assert(state.particle_sidecar.gravity_softening_comoving.size() == 5);
  assertParticleIdentityInvariant(state, baseline);
}

void test_particle_reorder_sidecar_invariants() {
  SimulationState state;
  seedSyntheticState(state);
  const auto baseline = captureParticleIdentity(state);

  const auto canonical = cosmosim::core::buildParticleReorderMap(state, cosmosim::core::ParticleReorderMode::kBySpecies);
  cosmosim::core::reorderParticles(state, canonical);
  assertParticleIdentityInvariant(state, baseline);

  const std::size_t n = state.particles.size();
  cosmosim::core::ParticleReorderMap reverse;
  reverse.new_to_old_index.resize(n);
  reverse.old_to_new_index.resize(n);
  for (std::size_t new_index = 0; new_index < n; ++new_index) {
    const auto old_index = static_cast<std::uint32_t>(n - 1U - new_index);
    reverse.new_to_old_index[new_index] = old_index;
    reverse.old_to_new_index[old_index] = static_cast<std::uint32_t>(new_index);
  }

  cosmosim::core::SidecarSyncPolicy move_sidecars;
  move_sidecars.star_particles = cosmosim::core::SidecarSyncMode::kMoveWithParent;
  move_sidecars.black_holes = cosmosim::core::SidecarSyncMode::kMoveWithParent;
  move_sidecars.tracers = cosmosim::core::SidecarSyncMode::kMoveWithParent;
  cosmosim::core::reorderParticles(state, reverse, move_sidecars);

  assertParticleIdentityInvariant(state, baseline);

  for (std::size_t row = 0; row < state.star_particles.size(); ++row) {
    const auto parent = state.star_particles.particle_index[row];
    const auto id = state.particle_sidecar.particle_id[parent];
    if (id == 703) {
      assert(state.star_particles.birth_mass_code[row] == 50.0);
      assert(state.star_particles.stellar_returned_metals_channel_cumulative_code[2][row] == 60.0);
    }
    if (id == 707) {
      assert(state.star_particles.birth_mass_code[row] == 51.0);
      assert(state.star_particles.stellar_returned_metals_channel_cumulative_code[2][row] == 61.0);
    }
  }

  compactFront(state, 4);
  assertParticleIdentityInvariant(state, baseline);
}

void test_particle_stale_view_invalidation() {
  SimulationState state;
  seedSyntheticState(state);

  cosmosim::core::TransientStepWorkspace workspace;
  std::vector<std::uint32_t> active{0U, 1U, 2U};
  auto gravity_view = cosmosim::core::buildGravityParticleKernelView(state, active, workspace);

  const auto reorder = cosmosim::core::buildParticleReorderMap(state, cosmosim::core::ParticleReorderMode::kBySfcKey);
  cosmosim::core::reorderParticles(state, reorder);

  bool stale_particle_view_threw = false;
  try {
    cosmosim::core::scatterGravityParticleKernelView(gravity_view, state);
  } catch (const std::runtime_error&) {
    stale_particle_view_threw = true;
  }
  assert(stale_particle_view_threw);

  state.resizeCells(3);
  for (std::size_t i = 0; i < state.cells.size(); ++i) {
    state.cells.center_x_comoving[i] = static_cast<double>(i);
    state.cells.center_y_comoving[i] = static_cast<double>(i);
    state.cells.center_z_comoving[i] = static_cast<double>(i);
    state.cells.mass_code[i] = 1.0;
    state.gas_cells.density_code[i] = 2.0;
    state.gas_cells.pressure_code[i] = 3.0;
  }

  std::vector<std::uint32_t> active_cells{0U, 2U};
  auto hydro_view = cosmosim::core::buildHydroCellKernelView(state, active_cells, workspace);
  state.resizeCells(2);

  bool stale_cell_view_threw = false;
  try {
    cosmosim::core::scatterHydroCellKernelView(hydro_view, state);
  } catch (const std::runtime_error&) {
    stale_cell_view_threw = true;
  }
  assert(stale_cell_view_threw);
}

}  // namespace

int main() {
  test_particle_resize_identity_invariants();
  test_particle_reorder_sidecar_invariants();
  test_particle_stale_view_invalidation();
  return 0;
}
