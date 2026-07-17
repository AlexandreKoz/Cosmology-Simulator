#include "cosmosim/workflows/runtime_resources.hpp"

#include <cassert>
#include <stdexcept>
#include <vector>

#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/core/time_integration.hpp"

namespace {

class MutableEpochSource final
    : public cosmosim::workflows::RuntimeResourceEpochSource {
 public:
  [[nodiscard]] cosmosim::workflows::RuntimeResourceEpoch
  currentRuntimeEpoch() const noexcept override {
    return epoch;
  }

  cosmosim::workflows::RuntimeResourceEpoch epoch{};
};

template <class View>
concept ExposesSimulationState = requires(View& view) {
  view.state;
};

template <class View>
concept ExposesEpochMutation = requires(View& view) {
  view.bumpParticleIndexGeneration();
};

template <class View>
concept ExposesOwnerContext = requires(View& view) {
  view.ownerContext();
};

template <class View>
concept ExposesParticleSidecar = requires(View& view) {
  view.particle_sidecar;
};

void expectStale(auto&& callback) {
  bool threw = false;
  try {
    callback();
  } catch (const std::runtime_error&) {
    threw = true;
  }
  assert(threw);
}

}  // namespace

int main() {
  using cosmosim::workflows::RuntimeEpochField;
  using cosmosim::workflows::RuntimeResourceLease;

  static_assert(
      !ExposesSimulationState<cosmosim::workflows::GravityStageView>);
  static_assert(
      !ExposesSimulationState<cosmosim::workflows::HydroAmrStageView>);
  static_assert(
      !ExposesSimulationState<cosmosim::workflows::SourceMutationStageView>);
  static_assert(
      !ExposesEpochMutation<cosmosim::workflows::MigrationOwnershipView>);
  static_assert(
      !ExposesOwnerContext<cosmosim::workflows::GravityStageView>);
  static_assert(
      !ExposesParticleSidecar<cosmosim::workflows::AnalysisStageView>);

  MutableEpochSource source;
  source.epoch = {
      .particle_index_generation = 4,
      .cell_index_generation = 7,
      .gas_identity_generation = 9,
      .scheduler_tick = 12,
      .step_index = 3,
  };

  cosmosim::workflows::GravityStageView gravity_view(RuntimeResourceLease(
      source,
      RuntimeEpochField::kParticleIndex |
          RuntimeEpochField::kSchedulerTick |
          RuntimeEpochField::kStepIndex));
  gravity_view.requireFresh();

  // A gravity view does not guard unrelated cell generations.
  ++source.epoch.cell_index_generation;
  gravity_view.requireFresh();

  ++source.epoch.scheduler_tick;
  expectStale([&]() { gravity_view.requireFresh(); });

  cosmosim::workflows::MigrationOwnershipView migration_view(
      RuntimeResourceLease(
          source,
          RuntimeEpochField::kParticleIndex |
              RuntimeEpochField::kCellIndex |
              RuntimeEpochField::kGasIdentity));
  migration_view.requireFresh();
  ++source.epoch.gas_identity_generation;
  expectStale([&]() { migration_view.requireFresh(); });

  cosmosim::workflows::AnalysisStageView analysis_view(RuntimeResourceLease(
      source,
      RuntimeEpochField::kStepIndex));
  analysis_view.requireFresh();
  ++source.epoch.step_index;
  expectStale([&]() { analysis_view.requireFresh(); });

  // The production epoch adapter observes the existing state/scheduler truth
  // rather than maintaining a second generation authority.
  cosmosim::core::SimulationState state;
  cosmosim::core::HierarchicalTimeBinScheduler scheduler(0);
  cosmosim::core::IntegratorState integrator_state;
  cosmosim::workflows::SimulationRuntimeEpochSource production_source(
      state, scheduler, integrator_state);
  cosmosim::workflows::MigrationOwnershipView production_view(
      RuntimeResourceLease(
          production_source,
          RuntimeEpochField::kParticleIndex |
              RuntimeEpochField::kCellIndex |
              RuntimeEpochField::kGasIdentity));
  production_view.requireFresh();
  state.bumpParticleIndexGeneration();
  expectStale([&]() { production_view.requireFresh(); });

  // A real dense-row reorder/compaction boundary bumps the same authoritative
  // generation observed by workflow leases; no manual lease invalidation is
  // required from the migration caller.
  state.resizeParticles(2U);
  state.particle_sidecar.particle_id = {2U, 1U};
  state.particle_sidecar.sfc_key = {20U, 10U};
  state.particle_sidecar.species_tag = {
      static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kDarkMatter),
      static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kDarkMatter)};
  state.species.count_by_species.fill(0U);
  state.species.count_by_species[static_cast<std::size_t>(
      cosmosim::core::ParticleSpecies::kDarkMatter)] = 2U;
  state.rebuildSpeciesIndex();
  cosmosim::workflows::MigrationOwnershipView pre_compaction_view(
      RuntimeResourceLease(
          production_source,
          RuntimeEpochField::kParticleIndex));
  pre_compaction_view.requireFresh();
  const cosmosim::core::ParticleReorderMap reorder =
      cosmosim::core::buildParticleReorderMap(
          state, cosmosim::core::ParticleReorderMode::kBySfcKey);
  cosmosim::core::reorderParticles(state, reorder);
  expectStale([&]() { pre_compaction_view.requireFresh(); });
  return 0;
}
