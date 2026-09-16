#include <cassert>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include "cosmosim/core/memory_governor.hpp"
#include "cosmosim/physics/stellar_evolution.hpp"

namespace {

[[nodiscard]] cosmosim::physics::StellarEvolutionTable testTable() {
  return cosmosim::physics::StellarEvolutionTable::loadFromTextFile(
      std::string(COSMOSIM_SOURCE_DIR) + "/resources/stellar_evolution/test_synthetic_v2.txt");
}

[[nodiscard]] cosmosim::core::SimulationState makeStarState(std::size_t star_count) {
  cosmosim::core::SimulationState state;
  state.resizeCells(1);
  state.cells.center_x_comoving[0] = 0.5;
  state.cells.center_y_comoving[0] = 0.5;
  state.cells.center_z_comoving[0] = 0.5;
  state.cells.mass_code[0] = 1.0;
  state.cells.patch_index[0] = 0;
  state.resizePatches(0);
  state.particles.mass_code.assign(star_count, 1.0);
  state.particles.position_x_comoving.assign(star_count, 0.5);
  state.particles.position_y_comoving.assign(star_count, 0.5);
  state.particles.position_z_comoving.assign(star_count, 0.5);
  state.star_particles.particle_index.resize(star_count);
  state.star_particles.birth_mass_code.assign(star_count, 1.0);
  state.star_particles.formation_scale_factor.assign(star_count, 0.5);
  state.star_particles.metallicity_mass_fraction.assign(star_count, 0.01);
  state.star_particles.stellar_age_years_last.assign(star_count, 1.0e7);
  state.star_particles.stellar_returned_mass_cumulative_code.assign(star_count, 0.0);
  state.star_particles.stellar_returned_metals_cumulative_code.assign(star_count, 0.0);
  state.star_particles.stellar_newly_synthesized_metals_cumulative_code.assign(star_count, 0.0);
  state.star_particles.stellar_feedback_energy_cumulative_erg.assign(star_count, 0.0);
  for (std::size_t c = 0; c < 3U; ++c) {
    state.star_particles.stellar_returned_mass_channel_cumulative_code[c].assign(star_count, 0.0);
    state.star_particles.stellar_returned_metals_channel_cumulative_code[c].assign(star_count, 0.0);
    state.star_particles.stellar_feedback_energy_channel_cumulative_erg[c].assign(star_count, 0.0);
  }
  for (std::size_t i = 0; i < star_count; ++i) {
    state.star_particles.particle_index[i] = static_cast<std::uint32_t>(i % 1U);
  }
  // Particle masses must cover referenced particle rows.
  if (!state.particles.mass_code.empty()) {
    state.particles.mass_code[0] = 10.0;
  }
  return state;
}

void testZeroBatch() {
  cosmosim::physics::StellarEvolutionBookkeeper keeper(
      cosmosim::physics::StellarEvolutionConfig{}, testTable());
  cosmosim::core::SimulationState state = makeStarState(0);
  cosmosim::physics::StellarEvolutionBatchWorkspace workspace;
  cosmosim::core::MemoryGovernor governor;
  const auto counters = keeper.evaluateElapsedYearsGoverned(
      state, std::span<const std::uint32_t>{}, 1.0e8, &governor, workspace);
  assert(counters.scanned_stars == 0U);
  assert(workspace.budgets.empty());
}

void testRetainedReuseGrowthOnlyAfterAdmission() {
  cosmosim::physics::StellarEvolutionBookkeeper keeper(
      cosmosim::physics::StellarEvolutionConfig{}, testTable());
  auto state = makeStarState(4);
  std::vector<std::uint32_t> active = {0U, 1U, 2U, 3U};
  cosmosim::physics::StellarEvolutionBatchWorkspace workspace;
  cosmosim::core::MemoryGovernor governor;
  const auto first = keeper.evaluateElapsedYearsGoverned(
      state, std::span<const std::uint32_t>(active.data(), active.size()), 1.0e8,
      &governor, workspace);
  const auto cap_after_first = workspace.budgets.capacity();
  assert(cap_after_first >= 4U);
  assert(workspace.reservation.committed());
  // Second batch of same size reuses retained capacity (no growth).
  const auto second = keeper.evaluateElapsedYearsGoverned(
      state, std::span<const std::uint32_t>(active.data(), active.size()), 1.0e8,
      &governor, workspace);
  assert(workspace.budgets.capacity() == cap_after_first);
  assert(first.evolved_stars == second.evolved_stars);
  assert(first.returned_mass_code == second.returned_mass_code);
}

void testTightHeadroomRejectionAndRetrySameEvents() {
  cosmosim::physics::StellarEvolutionBookkeeper keeper(
      cosmosim::physics::StellarEvolutionConfig{}, testTable());
  auto state = makeStarState(4);
  std::vector<std::uint32_t> active = {0U, 1U, 2U, 3U};
  const std::uint64_t need =
      cosmosim::physics::stellarEvolutionBatchStagingBytes(active.size());
  cosmosim::core::MemoryGovernor tight({.hard_limit_bytes = need - 1U});
  tight.setBaselineOwnedBytes(0U);
  cosmosim::physics::StellarEvolutionBatchWorkspace workspace;
  bool rejected = false;
  try {
    (void)keeper.evaluateElapsedYearsGoverned(
        state, std::span<const std::uint32_t>(active.data(), active.size()), 1.0e8,
        &tight, workspace);
  } catch (const cosmosim::core::MemoryAdmissionError&) {
    rejected = true;
  }
  assert(rejected);
  assert(workspace.budgets.capacity() == 0U);
  cosmosim::core::MemoryGovernor ample({.hard_limit_bytes = need + 1024U});
  ample.setBaselineOwnedBytes(0U);
  const auto counters = keeper.evaluateElapsedYearsGoverned(
      state, std::span<const std::uint32_t>(active.data(), active.size()), 1.0e8,
      &ample, workspace);
  // Diagnostic path produces identical counters/IDs.
  const auto report = keeper.evaluateElapsedYears(
      state, std::span<const std::uint32_t>(active.data(), active.size()), 1.0e8);
  assert(counters.evolved_stars == report.counters.evolved_stars);
  assert(std::abs(counters.returned_mass_code - report.counters.returned_mass_code) < 1.0e-12);
  assert(std::abs(counters.feedback_energy_erg - report.counters.feedback_energy_erg) < 1.0e-6);
  assert(workspace.budgets.size() == report.budgets.size());
  for (std::size_t i = 0; i < workspace.budgets.size(); ++i) {
    assert(workspace.budgets[i].star_index == report.budgets[i].star_index);
  }
}

}  // namespace

int main() {
  testZeroBatch();
  testRetainedReuseGrowthOnlyAfterAdmission();
  testTightHeadroomRejectionAndRetrySameEvents();
  return 0;
}
