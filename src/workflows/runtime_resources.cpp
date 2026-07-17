#include "cosmosim/workflows/runtime_resources.hpp"

#include <stdexcept>
#include <string>
#include <utility>

#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/core/time_integration.hpp"

namespace cosmosim::workflows {
namespace {

[[nodiscard]] bool guards(
    RuntimeEpochField guarded_fields,
    RuntimeEpochField field) noexcept {
  return (guarded_fields & field) != RuntimeEpochField::kNone;
}

[[nodiscard]] std::string staleField(
    const RuntimeResourceEpoch& captured,
    const RuntimeResourceEpoch& current,
    RuntimeEpochField guarded_fields) {
  if (guards(guarded_fields, RuntimeEpochField::kParticleIndex) &&
      captured.particle_index_generation != current.particle_index_generation) {
    return "particle_index_generation";
  }
  if (guards(guarded_fields, RuntimeEpochField::kCellIndex) &&
      captured.cell_index_generation != current.cell_index_generation) {
    return "cell_index_generation";
  }
  if (guards(guarded_fields, RuntimeEpochField::kGasIdentity) &&
      captured.gas_identity_generation != current.gas_identity_generation) {
    return "gas_identity_generation";
  }
  if (guards(guarded_fields, RuntimeEpochField::kSchedulerTick) &&
      captured.scheduler_tick != current.scheduler_tick) {
    return "scheduler_tick";
  }
  if (guards(guarded_fields, RuntimeEpochField::kStepIndex) &&
      captured.step_index != current.step_index) {
    return "step_index";
  }
  return {};
}

}  // namespace

SimulationRuntimeEpochSource::SimulationRuntimeEpochSource(
    const core::SimulationState& state,
    const core::HierarchicalTimeBinScheduler& scheduler,
    const core::IntegratorState& integrator_state) noexcept
    : m_state(&state),
      m_scheduler(&scheduler),
      m_integrator_state(&integrator_state) {}

RuntimeResourceEpoch
SimulationRuntimeEpochSource::currentRuntimeEpoch() const noexcept {
  return RuntimeResourceEpoch{
      .particle_index_generation = m_state->particleIndexGeneration(),
      .cell_index_generation = m_state->cellIndexGeneration(),
      .gas_identity_generation = m_state->gasCellIdentityGeneration(),
      .scheduler_tick = m_scheduler->currentTick(),
      .step_index = m_integrator_state->step_index,
  };
}

RuntimeResourceLease::RuntimeResourceLease(
    const RuntimeResourceEpochSource& source,
    RuntimeEpochField guarded_fields) noexcept
    : m_source(&source),
      m_captured_epoch(source.currentRuntimeEpoch()),
      m_guarded_fields(guarded_fields) {}

const RuntimeResourceEpoch& RuntimeResourceLease::capturedEpoch() const noexcept {
  return m_captured_epoch;
}

RuntimeEpochField RuntimeResourceLease::guardedFields() const noexcept {
  return m_guarded_fields;
}

bool RuntimeResourceLease::isFresh() const noexcept {
  if (m_source == nullptr) {
    return false;
  }
  return staleField(
             m_captured_epoch,
             m_source->currentRuntimeEpoch(),
             m_guarded_fields)
      .empty();
}

void RuntimeResourceLease::requireFresh(std::string_view view_name) const {
  if (m_source == nullptr) {
    throw std::runtime_error(
        "runtime resource view '" + std::string(view_name) +
        "' has no epoch authority");
  }
  const std::string stale_field = staleField(
      m_captured_epoch,
      m_source->currentRuntimeEpoch(),
      m_guarded_fields);
  if (!stale_field.empty()) {
    throw std::runtime_error(
        "runtime resource view '" + std::string(view_name) +
        "' is stale after " + stale_field + " changed");
  }
}

#define COSMOSIM_DEFINE_RUNTIME_STAGE_VIEW(Type, Name)                    \
  Type::Type(RuntimeResourceLease lease) noexcept                        \
      : m_lease(std::move(lease)) {}                                     \
  Type::Type(                                                            \
      RuntimeResourceLease lease, core::StepContext& context) noexcept   \
      : m_lease(std::move(lease)), m_context(&context) {}                \
  void Type::requireFresh() const { m_lease.requireFresh(Name); }         \
  core::StepContext& Type::ownerContext() const {                         \
    if (m_context == nullptr) {                                           \
      throw std::logic_error(                                             \
          "runtime resource view '" Name "' has no stage context");      \
    }                                                                     \
    return *m_context;                                                     \
  }

COSMOSIM_DEFINE_RUNTIME_STAGE_VIEW(DriftParticleStageView, "drift_particles")
COSMOSIM_DEFINE_RUNTIME_STAGE_VIEW(GravityStageView, "gravity")
COSMOSIM_DEFINE_RUNTIME_STAGE_VIEW(HydroAmrStageView, "hydro_amr")
COSMOSIM_DEFINE_RUNTIME_STAGE_VIEW(SourceMutationStageView, "source_mutation")
COSMOSIM_DEFINE_RUNTIME_STAGE_VIEW(AnalysisStageView, "analysis")
COSMOSIM_DEFINE_RUNTIME_STAGE_VIEW(OutputRestartStageView, "output_restart")

#undef COSMOSIM_DEFINE_RUNTIME_STAGE_VIEW

MigrationOwnershipView::MigrationOwnershipView(
    RuntimeResourceLease lease) noexcept
    : m_lease(std::move(lease)) {}

void MigrationOwnershipView::requireFresh() const {
  m_lease.requireFresh("migration_ownership");
}

TimeCriteriaStageView::TimeCriteriaStageView(
    RuntimeResourceLease lease) noexcept
    : m_lease(std::move(lease)) {}

void TimeCriteriaStageView::requireFresh() const {
  m_lease.requireFresh("time_criteria");
}

}  // namespace cosmosim::workflows
