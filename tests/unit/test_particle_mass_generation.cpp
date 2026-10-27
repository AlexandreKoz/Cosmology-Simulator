// SimulationState::particleMassGeneration contract (2026-10-27; AUTHORED,
// NOT RUN). Mass-derived caches (the hierarchical encounter criterion's local
// maximum mass) key on this generation, so every path that can change which
// mass a row holds must advance it, and paths that cannot must not be
// required to.
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "../support/test_require.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/core/simulation_state_runtime.hpp"

namespace {
using namespace cosmosim::core;

SimulationState makeState(std::size_t count) {
  SimulationState state;
  state.resizeParticles(count);
  for (std::size_t row = 0; row < count; ++row) {
    state.particle_sidecar.particle_id[row] = 10U + row;
    state.particle_sidecar.species_tag[row] = static_cast<std::uint32_t>(ParticleSpecies::kDarkMatter);
    state.particle_sidecar.owning_rank[row] = 0U;
    state.particles.mass_code[row] = 1.0 + static_cast<double>(row);  // mixed masses
  }
  state.rebuildSpeciesIndex();
  return state;
}

void testGenerationAdvances() {
  SimulationState empty;
  const auto initial = empty.particleMassGeneration();
  auto state = makeState(8U);
  // Population replacement (resize) advances the mass generation.
  CHUI_REQUIRE(state.particleMassGeneration() > initial);

  // An in-place mass write must be notified; it advances the mass generation
  // without touching the row-layout generation.
  const auto index_before = state.particleIndexGeneration();
  const auto mass_before = state.particleMassGeneration();
  state.particles.mass_code[3] = 42.0;
  state.notifyParticleMassMutation();
  CHUI_REQUIRE(state.particleMassGeneration() == mass_before + 1U);
  CHUI_REQUIRE(state.particleIndexGeneration() == index_before);

  // Any dense-row change (reorder, migration, compaction) also advances it.
  const auto mass_before_bump = state.particleMassGeneration();
  state.bumpParticleIndexGeneration();
  CHUI_REQUIRE(state.particleMassGeneration() > mass_before_bump);

  // A real reorder moves masses between rows and advances both generations.
  ParticleReorderMap map;
  for (std::uint32_t row = 0; row < 8U; ++row) {
    map.new_to_old_index.push_back(7U - row);
    map.old_to_new_index.push_back(7U - row);
  }
  const auto mass_before_reorder = state.particleMassGeneration();
  const auto index_before_reorder = state.particleIndexGeneration();
  reorderParticles(state, map);
  CHUI_REQUIRE(state.particleIndexGeneration() > index_before_reorder);
  CHUI_REQUIRE(state.particleMassGeneration() > mass_before_reorder);
  CHUI_REQUIRE(state.particles.mass_code[4] == 42.0);  // old row 3 -> new row 4

  // Positions and velocities are not mass state: drifting them must not
  // invalidate mass-derived caches (no hidden rescans on every block).
  const auto mass_before_drift = state.particleMassGeneration();
  for (auto& x : state.particles.position_x_comoving) x += 0.25;
  state.bumpGravitySourceGeneration();
  CHUI_REQUIRE(state.particleMassGeneration() == mass_before_drift);

  // Copies carry the generation (rollback snapshots never copy whole states,
  // but restart/test fixtures do).
  const SimulationState copy = state;
  CHUI_REQUIRE(copy.particleMassGeneration() == state.particleMassGeneration());
}

// Emulates the coordinator's cache key: (row count, mass generation, lane
// address). A key hit must imply the cached maximum equals a fresh scan.
void testCacheKeyEquivalence() {
  auto state = makeState(16U);
  struct Key {
    std::size_t rows;
    std::uint64_t generation;
    const double* data;
    bool operator==(const Key&) const = default;
  };
  const auto key = [&] {
    return Key{state.particles.mass_code.size(), state.particleMassGeneration(), state.particles.mass_code.data()};
  };
  const auto scan = [&] {
    double maximum = 0.0;
    for (const double mass : state.particles.mass_code) maximum = std::max(maximum, mass);
    return maximum;
  };
  Key cached_key = key();
  double cached = scan();
  // Unchanged population: legitimate hit, identical value.
  CHUI_REQUIRE(key() == cached_key && cached == scan());
  // In-place growth of one mass (e.g. black-hole accretion): key changes.
  state.particles.mass_code[5] = 1.0e3;
  state.notifyParticleMassMutation();
  CHUI_REQUIRE(!(key() == cached_key));
  cached_key = key();
  cached = scan();
  CHUI_REQUIRE(cached == 1.0e3);
  // Population change: key changes (count and generation).
  state.resizeParticles(20U);
  CHUI_REQUIRE(!(key() == cached_key));
}

}  // namespace

int main() {
  try {
    testGenerationAdvances();
    testCacheKeyEquivalence();
  } catch (const std::exception& error) {
    std::cerr << "particle mass generation: " << error.what() << '\n';
    return 1;
  }
  return 0;
}
