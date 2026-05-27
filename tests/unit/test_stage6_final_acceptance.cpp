#include <array>
#include <concepts>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include "cosmosim/analysis/diagnostics.hpp"
#include "cosmosim/analysis/halo_workflow.hpp"
#include "cosmosim/core/memory_accounting.hpp"
#include "cosmosim/core/provenance.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/core/time_integration.hpp"
#include "cosmosim/io/restart_checkpoint.hpp"
#include "cosmosim/gravity/tree_pm_coupling.hpp"
#include "cosmosim/physics/black_hole_agn.hpp"
#include "cosmosim/physics/star_formation.hpp"
#include "cosmosim/physics/stellar_evolution.hpp"
#include "cosmosim/physics/stellar_feedback.hpp"
#include "cosmosim/physics/tracer_support.hpp"

namespace {

template <typename T>
concept ExposesParticleIdentity = requires(T view) { view.particle_id; };

template <typename T>
concept ExposesParticleSpeciesTag = requires(T view) { view.species_tag; };

template <typename T>
concept ExposesOwningRank = requires(T view) { view.owning_rank; };

template <typename T>
concept ExposesGasPatchIndex = requires(T view) { view.patch_index; };

template <typename T>
concept ExposesGasTemperature = requires(T view) { view.temperature_code; };

template <typename T>
concept ExposesHydroReconstructionGradient = requires(T gas) {
  gas.recon_gradient_x;
  gas.recon_gradient_y;
  gas.recon_gradient_z;
};

template <typename T>
concept ExposesPersistentAccelerationLane = requires(T particles) {
  particles.acceleration_x_comoving;
  particles.acceleration_y_comoving;
  particles.acceleration_z_comoving;
};

template <typename T>
concept ExposesTransientWorkspaceInRestartPayload = requires(T payload) { payload.transient_workspace; };

template <typename T>
concept ExposesHydroScratchInRestartPayload = requires(T payload) { payload.hydro_scratch; };

template <typename T>
concept ExposesPmWorkspaceInRestartPayload = requires(T payload) { payload.pm_workspace; };

template <typename T>
concept ExposesMpiBuffersInRestartPayload = requires(T payload) { payload.mpi_buffers; };

template <typename T>
concept ExposesOutputBuffersInRestartPayload = requires(T payload) { payload.output_buffers; };

template <typename T>
concept StarFormationHasNarrowViewApi = requires(
    const T model,
    cosmosim::core::SimulationState& state,
    cosmosim::physics::StarFormationRuntimeView view) {
  { model.applyFromView(state, view, 0.001, 1.0, 7, 0) } -> std::same_as<cosmosim::physics::StarFormationStepReport>;
};

template <typename T>
concept StellarFeedbackHasNarrowViewApi = requires(
    const T model,
    cosmosim::core::SimulationState& state,
    cosmosim::physics::StellarFeedbackModuleState& module_state,
    cosmosim::physics::StellarFeedbackGeometryView geometry,
    cosmosim::physics::StellarFeedbackDepositionView deposition,
    std::span<const std::uint32_t> active_stars,
    std::span<const double> returned_mass,
    std::span<const double> returned_metals) {
  { model.applyWithViews(state, module_state, geometry, deposition, active_stars, returned_mass, returned_metals, 0.001) }
      -> std::same_as<cosmosim::physics::StellarFeedbackStepReport>;
};

template <typename T>
concept BlackHoleAgnHasNarrowViewApi = requires(
    const T model,
    cosmosim::physics::BlackHoleAgnAccretionView view) {
  { model.applyAccretionFromView(view, 0.001) } -> std::same_as<cosmosim::physics::BlackHoleAgnCounters>;
};

template <typename T>
concept TracerHasNarrowViewApi = requires(
    T model,
    cosmosim::physics::TracerHostMassView view) {
  { model.updateMassFromHostCellsView(view) } -> std::same_as<cosmosim::physics::TracerUpdateCounters>;
};

template <typename T>
concept StellarEvolutionHasNarrowViewApi = requires(
    T model,
    cosmosim::physics::StellarEvolutionRuntimeView view,
    double scale_factor,
    double dt_code) {
  { model.applyFromView(view, scale_factor, dt_code) } -> std::same_as<cosmosim::physics::StellarEvolutionStepReport>;
};

template <typename T>
concept ExposesAdaptiveTimeStepCriteriaView = requires(T view) {
  view.particles.velocity_x_peculiar;
  view.particles.accel_x_comoving;
  view.gas_cells.gas_particle_index_by_cell;
  view.gas_cells.sound_speed_code;
};

template <typename T>
concept FofHasNarrowHaloViewApi = requires(
    const T finder,
    cosmosim::analysis::HaloParticleView view,
    const cosmosim::core::SimulationConfig& config,
    cosmosim::analysis::FofProfilingCounters* counters) {
  { finder.buildCatalogFromView(view, config, 1, 1.0, counters) } -> std::same_as<cosmosim::analysis::HaloCatalog>;
};

void requireStage6CompileTimeContracts() {
  static_assert(cosmosim::core::k_is_canonical_particle_state_owner_v<cosmosim::core::SimulationState>);
  static_assert(!cosmosim::core::k_is_canonical_particle_state_owner_v<cosmosim::core::ParticleSoaStorage>);
  static_assert(!cosmosim::core::ParticleSoaStorage::k_owns_persistent_particle_truth);
  static_assert(!cosmosim::core::ParticleSoaStorage::k_is_restart_serializable);
  static_assert(!ExposesPersistentAccelerationLane<cosmosim::core::ParticleSoa>);
  static_assert(!ExposesHydroReconstructionGradient<cosmosim::core::GasCellSidecar>);
  static_assert(!ExposesParticleIdentity<cosmosim::core::GravityParticleKernelView>);
  static_assert(!ExposesParticleSpeciesTag<cosmosim::core::GravityParticleKernelView>);
  static_assert(!ExposesOwningRank<cosmosim::core::GravityParticleKernelView>);
  static_assert(!ExposesGasPatchIndex<cosmosim::core::HydroCellKernelView>);
  static_assert(!ExposesGasTemperature<cosmosim::core::HydroCellKernelView>);
  static_assert(!ExposesHydroReconstructionGradient<cosmosim::core::HydroCellKernelView>);
  static_assert(!ExposesTransientWorkspaceInRestartPayload<cosmosim::io::RestartWritePayload>);
  static_assert(!ExposesHydroScratchInRestartPayload<cosmosim::io::RestartWritePayload>);
  static_assert(!ExposesPmWorkspaceInRestartPayload<cosmosim::io::RestartWritePayload>);
  static_assert(!ExposesMpiBuffersInRestartPayload<cosmosim::io::RestartWritePayload>);
  static_assert(!ExposesOutputBuffersInRestartPayload<cosmosim::io::RestartWritePayload>);
  static_assert(StarFormationHasNarrowViewApi<cosmosim::physics::StarFormationModel>);
  static_assert(StellarFeedbackHasNarrowViewApi<cosmosim::physics::StellarFeedbackModel>);
  static_assert(BlackHoleAgnHasNarrowViewApi<cosmosim::physics::BlackHoleAgnModel>);
  static_assert(TracerHasNarrowViewApi<cosmosim::physics::TracerModel>);
  static_assert(StellarEvolutionHasNarrowViewApi<cosmosim::physics::StellarEvolutionBookkeeper>);
  static_assert(ExposesAdaptiveTimeStepCriteriaView<cosmosim::core::AdaptiveTimeStepCriteriaView>);
  static_assert(FofHasNarrowHaloViewApi<cosmosim::analysis::FofHaloFinder>);
}

cosmosim::core::SimulationState makeParticleState() {
  cosmosim::core::SimulationState state;
  state.resizeParticles(8);
  for (std::size_t i = 0; i < state.particles.size(); ++i) {
    state.particle_sidecar.particle_id[i] = 1000 + i;
    state.particle_sidecar.species_tag[i] = static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kDarkMatter);
    state.particle_sidecar.owning_rank[i] = 0;
    state.particles.position_x_comoving[i] = static_cast<double>(i);
    state.particles.position_y_comoving[i] = static_cast<double>(i) * 2.0;
    state.particles.position_z_comoving[i] = static_cast<double>(i) * 3.0;
    state.particles.velocity_x_peculiar[i] = 0.1;
    state.particles.velocity_y_peculiar[i] = 0.2;
    state.particles.velocity_z_peculiar[i] = 0.3;
    state.particles.mass_code[i] = 1.0 + static_cast<double>(i);
  }
  state.species.count_by_species[0] = state.particles.size();
  state.rebuildSpeciesIndex();
  return state;
}

cosmosim::core::SimulationState makeGasState() {
  cosmosim::core::SimulationState state;
  state.resizeParticles(5);
  state.resizeCells(3);
  state.species.count_by_species = {2, 3, 0, 0, 0};
  for (std::size_t i = 0; i < state.particles.size(); ++i) {
    state.particle_sidecar.particle_id[i] = 2000 + i;
    state.particle_sidecar.owning_rank[i] = 0;
    state.particle_sidecar.species_tag[i] = i < 2 ? static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kDarkMatter)
                                                  : static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kGas);
    state.particles.mass_code[i] = 1.0;
  }
  for (std::size_t i = 0; i < state.cells.size(); ++i) {
    state.cells.center_x_comoving[i] = static_cast<double>(i);
    state.cells.center_y_comoving[i] = 0.25 * static_cast<double>(i);
    state.cells.center_z_comoving[i] = 0.5 * static_cast<double>(i);
    state.cells.mass_code[i] = 10.0 + static_cast<double>(i);
    state.gas_cells.density_code[i] = 2.0 + static_cast<double>(i);
    state.gas_cells.pressure_code[i] = 1.0 + 0.5 * static_cast<double>(i);
    state.gas_cells.internal_energy_code[i] = 0.75 + 0.1 * static_cast<double>(i);
    state.gas_cells.temperature_code[i] = 100.0 + static_cast<double>(i);
  }
  state.rebuildSpeciesIndex();
  state.refreshGasCellIdentityFromParticleOrder();
  return state;
}

void testActiveKernelViewsRejectStaleGenerations() {
  cosmosim::core::SimulationState particle_state = makeParticleState();
  cosmosim::core::TransientStepWorkspace workspace;
  const std::array<std::uint32_t, 3> active_particles{1, 3, 6};
  auto particle_view = cosmosim::core::buildGravityParticleKernelView(particle_state, active_particles, workspace);
  particle_view.mass_code[1] = 123.0;
  particle_state.bumpParticleIndexGeneration();

  bool stale_particle_generation_threw = false;
  try {
    cosmosim::core::scatterGravityParticleKernelView(particle_view, particle_state);
  } catch (const std::runtime_error&) {
    stale_particle_generation_threw = true;
  }
  assert(stale_particle_generation_threw);

  cosmosim::core::SimulationState gas_state = makeGasState();
  cosmosim::core::TransientStepWorkspace hydro_workspace;
  const std::array<std::uint32_t, 2> active_cells{0, 2};
  auto cell_view = cosmosim::core::buildHydroCellKernelView(gas_state, active_cells, hydro_workspace);
  cell_view.mass_code[0] = 88.0;
  gas_state.bumpCellIndexGeneration();

  bool stale_cell_generation_threw = false;
  try {
    cosmosim::core::scatterHydroCellKernelView(cell_view, gas_state);
  } catch (const std::runtime_error&) {
    stale_cell_generation_threw = true;
  }
  assert(stale_cell_generation_threw);
}

void testMemoryReportsCoverStage6CategoriesWithoutCountingViewsAsOwners() {
  cosmosim::core::SimulationState state = makeGasState();
  state.star_particles.resize(1);
  state.black_holes.resize(1);
  state.tracers.resize(1);
  state.sidecars.upsert(cosmosim::core::ModuleSidecarBlock{
      .module_name = "stage6_probe",
      .schema_version = 1,
      .payload = {std::byte{0x01}, std::byte{0x02}, std::byte{0x03}}});

  cosmosim::core::TransientStepWorkspace workspace;
  workspace.particle_position_x_comoving.reserve(16);
  workspace.gravity_particle_index.reserve(16);
  workspace.hydro_cell_index.reserve(8);
  workspace.hydro_recon_gradient_x.reserve(8);
  static_cast<void>(workspace.scratch.allocateBytes(256, alignof(double)));

  const cosmosim::core::MemoryReport runtime_report = cosmosim::core::collectSimulationMemoryReport(state, &workspace);
  assert(runtime_report.totals.persistent_total_bytes > 0);
  assert(runtime_report.totals.transient_total_bytes > 0);
  assert(runtime_report.totals.persistent_by_subsystem[cosmosim::core::memorySubsystemIndex(cosmosim::core::MemorySubsystem::kParticles)] > 0);
  assert(runtime_report.totals.persistent_by_subsystem[cosmosim::core::memorySubsystemIndex(cosmosim::core::MemorySubsystem::kGasHydro)] > 0);
  assert(runtime_report.totals.persistent_by_subsystem[cosmosim::core::memorySubsystemIndex(cosmosim::core::MemorySubsystem::kSidecars)] > 0);
  assert(runtime_report.totals.transient_by_subsystem[cosmosim::core::memorySubsystemIndex(cosmosim::core::MemorySubsystem::kActiveSets)] > 0);
  assert(runtime_report.totals.transient_by_subsystem[cosmosim::core::memorySubsystemIndex(cosmosim::core::MemorySubsystem::kScratch)] >= 256);

  bool saw_nonowning_view_note = false;
  for (const std::string& note : runtime_report.notes) {
    saw_nonowning_view_note = saw_nonowning_view_note || note.find("spans/views") != std::string::npos;
  }
  assert(saw_nonowning_view_note);

  cosmosim::core::MemoryBudgetEstimateInput estimate_input;
  estimate_input.particle_capacity = 1024;
  estimate_input.gas_cell_capacity = 512;
  estimate_input.star_capacity = 128;
  estimate_input.black_hole_capacity = 8;
  estimate_input.tracer_capacity = 64;
  estimate_input.active_particle_capacity = 256;
  estimate_input.active_cell_capacity = 128;
  estimate_input.tree_node_capacity = 2048;
  estimate_input.pm_grid_cells = 4096;
  estimate_input.mpi_exchange_particle_capacity = 128;
  estimate_input.output_buffer_bytes = 8192;
  const cosmosim::core::MemoryReport estimate = cosmosim::core::estimatePreRunMemoryBudget(estimate_input);
  for (std::size_t i = 0; i < static_cast<std::size_t>(cosmosim::core::MemorySubsystem::kCount); ++i) {
    const auto subsystem = static_cast<cosmosim::core::MemorySubsystem>(i);
    const std::uint64_t total = estimate.totals.persistent_by_subsystem[i] + estimate.totals.transient_by_subsystem[i] +
                                estimate.totals.unknown_by_subsystem[i];
    assert(total > 0 || subsystem == cosmosim::core::MemorySubsystem::kScratch);
  }
}

void testDiagnosticsStateWrappersDelegateToNarrowViews() {
  cosmosim::core::SimulationState state = makeGasState();
  state.resizeParticles(6);
  state.species.count_by_species = {2, 3, 1, 0, 0};
  for (std::size_t i = 0; i < state.particles.size(); ++i) {
    state.particle_sidecar.particle_id[i] = 3000 + i;
    state.particle_sidecar.owning_rank[i] = 0;
    state.particle_sidecar.species_tag[i] = i < 2 ? static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kDarkMatter)
                                                  : (i < 5 ? static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kGas)
                                                           : static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kStar));
    state.particles.position_x_comoving[i] = static_cast<double>(i) * 0.1;
    state.particles.position_y_comoving[i] = static_cast<double>(i) * 0.2;
    state.particles.position_z_comoving[i] = static_cast<double>(i) * 0.3;
    state.particles.velocity_x_peculiar[i] = 0.1 * static_cast<double>(i);
    state.particles.velocity_y_peculiar[i] = 0.2 * static_cast<double>(i);
    state.particles.velocity_z_peculiar[i] = 0.3 * static_cast<double>(i);
    state.particles.mass_code[i] = 1.0 + static_cast<double>(i);
  }
  state.star_particles.resize(1);
  state.star_particles.particle_index[0] = 5;
  state.star_particles.formation_scale_factor[0] = 0.5;
  state.star_particles.birth_mass_code[0] = 2.0;
  state.rebuildSpeciesIndex();
  state.refreshGasCellIdentityFromParticleOrder();

  cosmosim::core::SimulationConfig config = cosmosim::core::makeUnvalidatedSimulationConfigForTests();
  cosmosim::analysis::DiagnosticsEngine engine(config);
  const cosmosim::analysis::DiagnosticsStateView view = cosmosim::analysis::buildDiagnosticsStateView(state);
  const auto health_from_state = engine.computeRunHealth(state);
  const auto health_from_view = engine.computeRunHealth(view);
  assert(health_from_state.particle_count == health_from_view.particle_count);
  assert(health_from_state.cell_count == health_from_view.cell_count);
  assert(health_from_state.non_finite_particles == health_from_view.non_finite_particles);

  const auto angular_from_state = engine.computeAngularMomentumBudget(state);
  const auto angular_from_view = engine.computeAngularMomentumBudget(view.particles);
  for (std::size_t component = 0; component < 3; ++component) {
    assert(std::abs(angular_from_state.total_l_code[component] - angular_from_view.total_l_code[component]) < 1.0e-14);
  }

  const auto sfh_from_state = engine.computeStarFormationHistory(state, 4);
  const auto sfh_from_view = engine.computeStarFormationHistory(view.stars, 4);
  assert(sfh_from_state.size() == sfh_from_view.size());
  for (std::size_t i = 0; i < sfh_from_state.size(); ++i) {
    assert(sfh_from_state[i].formed_mass_code == sfh_from_view[i].formed_mass_code);
  }
}


void testTreePmRuntimeMemoryReportAccountsLiveSolverBuffers() {
  cosmosim::gravity::TreePmCoordinator coordinator(cosmosim::gravity::PmGridShape{4, 4, 4});
  const cosmosim::core::MemoryReport report = coordinator.memoryReport();
  bool saw_pm_mesh = false;
  bool saw_tree_category = false;
  bool saw_mpi_buffers = false;
  for (const auto& entry : report.entries) {
    saw_pm_mesh = saw_pm_mesh || entry.label == "pm_mesh.density";
    saw_tree_category = saw_tree_category || entry.subsystem == cosmosim::core::MemorySubsystem::kTree;
    saw_mpi_buffers = saw_mpi_buffers || entry.label == "treepm.exchange.send_payload";
  }
  assert(saw_pm_mesh);
  assert(saw_tree_category);
  assert(saw_mpi_buffers);
  assert(report.totals.transient_by_subsystem[cosmosim::core::memorySubsystemIndex(cosmosim::core::MemorySubsystem::kPmMesh)] > 0);
}

void testRestartPayloadHashIgnoresTransientScratchButCoversPersistentTruth() {
  cosmosim::core::SimulationState state = makeParticleState();
  cosmosim::core::IntegratorState integrator_state;
  cosmosim::core::HierarchicalTimeBinScheduler scheduler(2);
  scheduler.reset(state.particles.size(), 0, 0);

  cosmosim::io::RestartWritePayload payload;
  payload.persistent_state.simulation_state = &state;
  payload.integrator_state = &integrator_state;
  payload.scheduler = &scheduler;
  payload.normalized_config_text = "stage6_final_acceptance = true\n";
  payload.normalized_config_hash_hex = cosmosim::core::stableConfigHashHex(payload.normalized_config_text);
  payload.provenance = cosmosim::core::makeProvenanceRecord(payload.normalized_config_hash_hex, "stage6f");
  payload.distributed_gravity_state.schema_version = 2;
  payload.distributed_gravity_state.world_size = 1;
  payload.distributed_gravity_state.pm_grid_nx = 4;
  payload.distributed_gravity_state.pm_grid_ny = 4;
  payload.distributed_gravity_state.pm_grid_nz = 4;
  payload.distributed_gravity_state.owning_rank_by_item.assign(state.particles.size(), 0);
  payload.distributed_gravity_state.pm_slab_begin_x_by_rank = {0};
  payload.distributed_gravity_state.pm_slab_end_x_by_rank = {4};

  const std::uint64_t base_hash = cosmosim::io::restartPayloadIntegrityHash(payload);

  cosmosim::core::TransientStepWorkspace workspace;
  workspace.particle_position_x_comoving.reserve(64);
  workspace.hydro_recon_gradient_x.reserve(32);
  static_cast<void>(workspace.scratch.allocateBytes(1024, alignof(double)));
  const std::uint64_t hash_after_transient_mutation = cosmosim::io::restartPayloadIntegrityHash(payload);
  assert(hash_after_transient_mutation == base_hash);

  state.particles.mass_code[2] += 0.25;
  const std::uint64_t hash_after_persistent_mutation = cosmosim::io::restartPayloadIntegrityHash(payload);
  assert(hash_after_persistent_mutation != base_hash);
}

}  // namespace

int main() {
  requireStage6CompileTimeContracts();
  testActiveKernelViewsRejectStaleGenerations();
  testMemoryReportsCoverStage6CategoriesWithoutCountingViewsAsOwners();
  testTreePmRuntimeMemoryReportAccountsLiveSolverBuffers();
  testDiagnosticsStateWrappersDelegateToNarrowViews();
  testRestartPayloadHashIgnoresTransientScratchButCoversPersistentTruth();
  return 0;
}
