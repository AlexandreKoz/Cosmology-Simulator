#include "cosmosim/workflows/hydro_amr_runtime.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <numeric>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "cosmosim/amr/amr_hydro_orchestrator.hpp"
#include "cosmosim/core/profiling.hpp"
#include "cosmosim/hydro/hydro_boundary_conditions.hpp"
#include "cosmosim/hydro/hydro_cartesian_patch.hpp"
#include "cosmosim/hydro/hydro_reconstruction.hpp"
#include "cosmosim/hydro/hydro_riemann.hpp"
#include "cosmosim/workflows/runtime_services.hpp"
#include "workflows/internal/runtime_stage_resource_access.hpp"
#include "workflows/internal/cartesian_gas_cell_layout.hpp"
#include "workflows/internal/gas_cell_ownership.hpp"
#include "workflows/internal/amr_migration_payload.hpp"
#include "workflows/internal/particle_ghost_runtime.hpp"

#if COSMOSIM_ENABLE_MPI
#include <mpi.h>
#endif

namespace cosmosim::workflows {
namespace {

constexpr double k_gamma_adiabatic = 5.0 / 3.0;
constexpr double k_pressure_floor = 1.0e-10;
constexpr double k_density_floor = 1.0e-10;
[[nodiscard]] internal::CartesianGasCellRowLayout requireCartesianGasCellRowLayout(
    const core::SimulationState& state,
    const core::SimulationConfig& config,
    std::string_view caller) {
  internal::CartesianGasCellLayoutBuildResult result =
      internal::buildCartesianGasCellRowLayout(state, config);
  if (!result.ok()) {
    throw std::runtime_error(
        std::string(caller) + ": fixed Cartesian hydro geometry rejected: " + result.diagnostic);
  }
  return std::move(result.layout);
}

struct HydroGhostConservedSnapshot {
  std::vector<std::uint32_t> cell_indices;
  std::vector<std::uint64_t> gas_cell_ids;
  std::vector<std::uint64_t> parent_particle_ids;
  std::vector<std::uint32_t> owner_ranks;
  std::vector<hydro::HydroConservedState> conserved_state;
};

struct HydroConservativeGhostSyncReport {
  std::size_t restored_ghost_cells = 0;
  double rejected_remote_delta_l1 = 0.0;
  std::vector<parallel::HydroConservativeFluxCorrectionRecord> correction_records;
};

struct HydroRemoteGhostBoundaryReport {
  std::size_t boundary_metadata_sent = 0;
  std::size_t boundary_metadata_received = 0;
  std::size_t requested_cells = 0;
  std::size_t received_cells = 0;
  std::size_t imported_ghosts = 0;
  std::size_t interface_faces = 0;
  std::size_t stale_or_invalid_payloads = 0;
  std::uint64_t request_bytes = 0;
  std::uint64_t payload_bytes = 0;
};

struct HydroBoundaryCellAdvertisement {
  std::uint64_t gas_cell_id = 0;
  int owner_rank = 0;
  std::uint64_t hydro_sync_epoch = 0;
  std::uint64_t decomposition_epoch = 0;
  double center_x_comoving = 0.0;
  double center_y_comoving = 0.0;
  double center_z_comoving = 0.0;
};

[[nodiscard]] bool finiteBoundaryAdvertisement(const HydroBoundaryCellAdvertisement& record) {
  return record.gas_cell_id != 0U &&
      record.owner_rank >= 0 &&
      std::isfinite(record.center_x_comoving) &&
      std::isfinite(record.center_y_comoving) &&
      std::isfinite(record.center_z_comoving);
}

[[nodiscard]] std::vector<HydroBoundaryCellAdvertisement> exchangeHydroBoundaryAdvertisements(
    const parallel::MpiContext& mpi_context,
    std::span<const HydroBoundaryCellAdvertisement> local_records) {
  for (const HydroBoundaryCellAdvertisement& record : local_records) {
    if (!finiteBoundaryAdvertisement(record)) {
      throw std::invalid_argument("hydro boundary advertisement is invalid");
    }
    if (record.owner_rank != mpi_context.worldRank()) {
      throw std::invalid_argument("hydro boundary advertisement owner rank does not match MPI context");
    }
  }
  if (!mpi_context.isEnabled()) {
    return std::vector<HydroBoundaryCellAdvertisement>(local_records.begin(), local_records.end());
  }
#if COSMOSIM_ENABLE_MPI
  const int world_size = mpi_context.worldSize();
  const std::uint64_t local_count = static_cast<std::uint64_t>(local_records.size());
  std::vector<std::uint64_t> counts64(static_cast<std::size_t>(world_size), 0U);
  MPI_Allgather(
      const_cast<std::uint64_t*>(&local_count),
      1,
      MPI_UINT64_T,
      counts64.data(),
      1,
      MPI_UINT64_T,
      MPI_COMM_WORLD);
  std::vector<int> recv_counts(static_cast<std::size_t>(world_size), 0);
  std::vector<int> recv_displs(static_cast<std::size_t>(world_size), 0);
  std::uint64_t total_records = 0;
  for (int rank = 0; rank < world_size; ++rank) {
    const std::uint64_t bytes = counts64[static_cast<std::size_t>(rank)] * sizeof(HydroBoundaryCellAdvertisement);
    if (bytes > static_cast<std::uint64_t>(std::numeric_limits<int>::max())) {
      throw std::overflow_error("hydro boundary advertisement exchange byte count exceeds MPI int limit");
    }
    recv_counts[static_cast<std::size_t>(rank)] = static_cast<int>(bytes);
    if (rank > 0) {
      recv_displs[static_cast<std::size_t>(rank)] =
          recv_displs[static_cast<std::size_t>(rank - 1)] + recv_counts[static_cast<std::size_t>(rank - 1)];
    }
    total_records += counts64[static_cast<std::size_t>(rank)];
  }
  std::vector<HydroBoundaryCellAdvertisement> result(static_cast<std::size_t>(total_records));
  MPI_Allgatherv(
      const_cast<HydroBoundaryCellAdvertisement*>(local_records.data()),
      static_cast<int>(local_records.size() * sizeof(HydroBoundaryCellAdvertisement)),
      MPI_BYTE,
      result.data(),
      recv_counts.data(),
      recv_displs.data(),
      MPI_BYTE,
      MPI_COMM_WORLD);
  for (const HydroBoundaryCellAdvertisement& record : result) {
    if (!finiteBoundaryAdvertisement(record) || record.owner_rank >= world_size) {
      throw std::runtime_error("hydro boundary advertisement exchange returned invalid metadata");
    }
  }
  return result;
#else
  throw std::runtime_error("hydro boundary advertisement exchange requires MPI support when MPI context is enabled");
#endif
}

[[nodiscard]] std::uint8_t hydroAxisCode(hydro::HydroFaceAxis axis) {
  switch (axis) {
    case hydro::HydroFaceAxis::kX:
      return 0U;
    case hydro::HydroFaceAxis::kY:
      return 1U;
    case hydro::HydroFaceAxis::kZ:
      return 2U;
  }
  return 0U;
}

[[nodiscard]] std::uint8_t hydroSideCode(hydro::HydroFaceSide side) {
  return side == hydro::HydroFaceSide::kLower ? 0U : 1U;
}

[[nodiscard]] std::uint64_t hydroInterfaceFaceKey(
    std::uint64_t local_gas_cell_id,
    std::uint64_t remote_gas_cell_id,
    hydro::HydroFaceAxis axis) {
  const std::uint64_t lo = std::min(local_gas_cell_id, remote_gas_cell_id);
  const std::uint64_t hi = std::max(local_gas_cell_id, remote_gas_cell_id);
  std::uint64_t key = 1469598103934665603ULL;
  auto mix = [&](std::uint64_t value) {
    key ^= value;
    key *= 1099511628211ULL;
  };
  mix(lo);
  mix(hi);
  mix(static_cast<std::uint64_t>(hydroAxisCode(axis) + 1U));
  return key == 0U ? 1U : key;
}

[[nodiscard]] bool nearlySameCoordinate(double lhs, double rhs, double tolerance) {
  return std::abs(lhs - rhs) <= tolerance;
}

[[nodiscard]] HydroGhostConservedSnapshot snapshotHydroGhostConservedCells(
    const core::SimulationState& state,
    const hydro::HydroConservedStateSoa& conserved,
    const hydro::HydroPatchGeometry& geometry,
    std::span<const std::uint32_t> geometry_row_by_dense_row,
    std::uint32_t world_rank) {
  HydroGhostConservedSnapshot snapshot;
  state.requireGasCellIdentityMapCoversDenseRows("snapshotHydroGhostConservedCells");
  const auto particle_row_by_id = internal::buildParticleRowById(state);
  snapshot.cell_indices.reserve(state.cells.size());
  snapshot.gas_cell_ids.reserve(state.cells.size());
  snapshot.parent_particle_ids.reserve(state.cells.size());
  snapshot.owner_ranks.reserve(state.cells.size());
  snapshot.conserved_state.reserve(state.cells.size());
  for (std::uint32_t cell_index = 0; cell_index < state.cells.size(); ++cell_index) {
    const std::uint32_t owner_rank = internal::gasCellOwnerRankForLocalRow(
        state,
        cell_index,
        particle_row_by_id,
        "snapshotHydroGhostConservedCells");
    if (owner_rank == world_rank) {
      continue;
    }
    const core::GasCellIdentityRecord& identity = internal::gasCellIdentityRecordForLocalRow(
        state,
        cell_index,
        "snapshotHydroGhostConservedCells");
    snapshot.cell_indices.push_back(cell_index);
    snapshot.gas_cell_ids.push_back(identity.gas_cell_id);
    snapshot.parent_particle_ids.push_back(identity.parent_particle_id.value_or(0U));
    snapshot.owner_ranks.push_back(owner_rank);
    if (cell_index >= geometry_row_by_dense_row.size()) {
      throw std::out_of_range("snapshotHydroGhostConservedCells: dense row is outside Cartesian geometry map");
    }
    snapshot.conserved_state.push_back(conserved.loadCell(geometry_row_by_dense_row[cell_index]));
  }
  for (const hydro::HydroGhostCell& ghost : geometry.ghost_cells) {
    if (ghost.boundary_kind != hydro::HydroBoundaryKind::kImportedMpi) {
      continue;
    }
    if (ghost.origin_rank < 0 || ghost.origin_rank == static_cast<int>(world_rank) ||
        ghost.origin_gas_cell_id == 0U || ghost.ghost_cell >= conserved.size()) {
      throw std::runtime_error("snapshotHydroGhostConservedCells: imported MPI ghost metadata is invalid");
    }
    snapshot.cell_indices.push_back(static_cast<std::uint32_t>(ghost.ghost_cell));
    snapshot.gas_cell_ids.push_back(ghost.origin_gas_cell_id);
    snapshot.parent_particle_ids.push_back(0U);
    snapshot.owner_ranks.push_back(static_cast<std::uint32_t>(ghost.origin_rank));
    snapshot.conserved_state.push_back(conserved.loadCell(ghost.ghost_cell));
  }
  return snapshot;
}

[[nodiscard]] HydroConservativeGhostSyncReport restoreHydroGhostConservedCells(
    hydro::HydroConservedStateSoa& conserved,
    const HydroGhostConservedSnapshot& snapshot,
    std::uint32_t world_rank) {
  if (snapshot.cell_indices.size() != snapshot.conserved_state.size() ||
      snapshot.cell_indices.size() != snapshot.gas_cell_ids.size() ||
      snapshot.cell_indices.size() != snapshot.parent_particle_ids.size() ||
      snapshot.cell_indices.size() != snapshot.owner_ranks.size()) {
    throw std::invalid_argument("hydro ghost snapshot is inconsistent");
  }
  HydroConservativeGhostSyncReport report;
  report.correction_records.reserve(snapshot.cell_indices.size());
  for (std::size_t i = 0; i < snapshot.cell_indices.size(); ++i) {
    const std::uint32_t cell_index = snapshot.cell_indices[i];
    if (cell_index >= conserved.size()) {
      throw std::out_of_range("hydro ghost snapshot cell index is outside conserved state");
    }
    const hydro::HydroConservedState before = snapshot.conserved_state[i];
    const hydro::HydroConservedState after = conserved.loadCell(cell_index);
    const double dm = after.mass_density_comoving - before.mass_density_comoving;
    const double dmx = after.momentum_density_x_comoving - before.momentum_density_x_comoving;
    const double dmy = after.momentum_density_y_comoving - before.momentum_density_y_comoving;
    const double dmz = after.momentum_density_z_comoving - before.momentum_density_z_comoving;
    const double de = after.total_energy_density_comoving - before.total_energy_density_comoving;
    report.rejected_remote_delta_l1 += std::abs(dm);
    report.rejected_remote_delta_l1 += std::abs(dmx);
    report.rejected_remote_delta_l1 += std::abs(dmy);
    report.rejected_remote_delta_l1 += std::abs(dmz);
    report.rejected_remote_delta_l1 += std::abs(de);
    if (dm != 0.0 || dmx != 0.0 || dmy != 0.0 || dmz != 0.0 || de != 0.0) {
      parallel::HydroConservativeFluxCorrectionRecord record;
      record.gas_cell_id = snapshot.gas_cell_ids[i];
      record.parent_particle_id = snapshot.parent_particle_ids[i];
      record.source_rank = static_cast<int>(world_rank);
      record.owner_rank = static_cast<int>(snapshot.owner_ranks[i]);
      record.delta_mass_density_comoving = dm;
      record.delta_momentum_density_x_comoving = dmx;
      record.delta_momentum_density_y_comoving = dmy;
      record.delta_momentum_density_z_comoving = dmz;
      record.delta_total_energy_density_comoving = de;
      parallel::validateHydroConservativeFluxCorrectionRecord(record);
      report.correction_records.push_back(record);
    }
    conserved.storeCell(cell_index, before);
    ++report.restored_ghost_cells;
  }
  return report;
}



class HydroAmrRuntimeImpl final : public HydroAmrRuntime {
 public:
  HydroAmrRuntimeImpl(
      const core::SimulationConfig& config,
      const core::ModePolicy& mode_policy,
      const GravityAccelerationProvider& gravity_callback,
      const RuntimeServices& services)
      : m_config(config),
        m_mode_policy(mode_policy),
        m_gravity_callback(gravity_callback),
        m_mpi_context(services.mpi_context),
        m_solver(k_gamma_adiabatic),
        m_reconstruction(hydro::HydroReconstructionPolicy{
            .limiter = hydro::HydroSlopeLimiter::kMonotonizedCentral,
            .dt_over_dx_code = 0.0,
            .rho_floor = k_density_floor,
            .pressure_floor = k_pressure_floor,
            .enable_muscl_hancock_predictor = true,
            .adiabatic_index = k_gamma_adiabatic,
        }) {}

  [[nodiscard]] const hydro::HydroProfileEvent& lastHydroProfile() const noexcept { return m_last_hydro_profile; }
  [[nodiscard]] const internal::SolverGhostRefreshReport&
  lastGhostRefreshReport() const noexcept {
    return m_last_ghost_refresh;
  }
  [[nodiscard]] const HydroRemoteGhostBoundaryReport& lastRemoteGhostBoundaryReport() const noexcept {
    return m_last_remote_ghost_boundary_report;
  }
  [[nodiscard]] const core::HydroCflDiagnostics& lastHydroCflDiagnostics() const noexcept {
    return m_last_hydro_cfl_diagnostics;
  }
  [[nodiscard]] std::uint64_t ghostExchangeBytesRecent() const noexcept override {
    return m_last_ghost_refresh.sent_bytes + m_last_ghost_refresh.received_bytes;
  }
  [[nodiscard]] std::size_t remoteImportedGhostCount() const noexcept override {
    return m_last_remote_ghost_boundary_report.imported_ghosts;
  }
  [[nodiscard]] std::size_t remoteInterfaceFaceCount() const noexcept override {
    return m_last_remote_ghost_boundary_report.interface_faces;
  }
  [[nodiscard]] std::size_t remoteStaleInvalidPayloadCount() const noexcept override {
    return m_last_remote_ghost_boundary_report.stale_or_invalid_payloads;
  }

  void execute(HydroAmrStageView& view) override {
    view.requireFresh();
    core::StepContext& context = internal::RuntimeStageAccess::hydroAmrContext(
        view,
        {{RuntimeResourceKey::kParticleVelocity, RuntimeResourceAccessMode::kReadWrite},
         {RuntimeResourceKey::kHydroConservedState, RuntimeResourceAccessMode::kReadWrite},
         {RuntimeResourceKey::kHydroPrimitiveState, RuntimeResourceAccessMode::kReadWrite},
         {RuntimeResourceKey::kAmrPatchState, RuntimeResourceAccessMode::kReadWrite},
         {RuntimeResourceKey::kGravityAcceleration, RuntimeResourceAccessMode::kRead},
         {RuntimeResourceKey::kMigrationOwnership, RuntimeResourceAccessMode::kRead},
         {RuntimeResourceKey::kIntegratorTruth, RuntimeResourceAccessMode::kRead}});
    if (context.stage != core::IntegrationStage::kHydroUpdate) {
      throw std::logic_error("hydro handler received an unregistered stage");
    }
    const auto hydro_cell_rank_count_after_refresh = [&]() {
      const std::uint64_t local_has_hydro_cells = context.state.cells.size() == 0 ? 0ULL : 1ULL;
      return m_mpi_context.allreduceSumUint64(local_has_hydro_cells);
    };
    m_current_world_rank = static_cast<std::uint32_t>(std::max(m_mpi_context.worldRank(), 0));
    m_last_ghost_refresh = internal::refreshParticleGhostsForSolver(
        context, m_mpi_context, "hydro.godunov", &m_ghost_cache_lifecycle);
    const std::uint64_t hydro_cell_rank_count = hydro_cell_rank_count_after_refresh();
    const std::uint64_t local_has_production_amr =
        amr::hasProductionAmrHydroCoverage(context.state) ? 1ULL : 0ULL;
    const std::uint64_t production_amr_rank_count =
        m_mpi_context.allreduceSumUint64(local_has_production_amr);
    const bool all_ranks_have_hydro_cells =
        !m_mpi_context.isEnabled() ||
        hydro_cell_rank_count == static_cast<std::uint64_t>(m_mpi_context.worldSize());
    m_last_hydro_profile = {};

    // Production AMR exchanges are collective.  Ranks with zero local cells
    // must still enter the control-plane and zero-record payload exchanges
    // whenever any peer owns a production patch.
    if (production_amr_rank_count > 0U) {
      if (context.state.cells.size() != 0U) {
        context.state.requireGasCellIdentityMapCoversDenseRows("hydro callback");
      }
      runProductionAmrHydroPath(context);
      return;
    }
    if (context.state.cells.size() == 0U) {
      return;
    }

    context.state.requireGasCellIdentityMapCoversDenseRows("hydro callback");
    const auto particle_row_by_id = internal::buildParticleRowById(context.state);
    std::vector<std::uint8_t> parent_mirror_updated(context.state.particles.size(), 0U);

    if (m_mpi_context.isEnabled()) {
      m_cached_geometry_key.reset();
    }
    rebuildGeometryIfNeeded(context.state, context.integrator_state.dt_time_code);
    verifyAcceptedHydroCfl(context, context.active_set.cell_indices);

    m_conserved.resize(m_geometry.totalCellStorageCount());
    for (std::size_t cell_index = 0; cell_index < context.state.cells.size(); ++cell_index) {
      const double rho = std::max(context.state.gas_cells.density_code[cell_index], k_density_floor);
      double pressure = context.state.gas_cells.pressure_code[cell_index];
      if (pressure <= 0.0) {
        const double internal_energy = std::max(context.state.gas_cells.internal_energy_code[cell_index], k_pressure_floor);
        pressure = std::max((k_gamma_adiabatic - 1.0) * rho * internal_energy, k_pressure_floor);
      }
      const hydro::HydroPrimitiveState primitive{
          .rho_comoving = rho,
          .vel_x_peculiar = context.state.gas_cells.velocity_x_peculiar[cell_index],
          .vel_y_peculiar = context.state.gas_cells.velocity_y_peculiar[cell_index],
          .vel_z_peculiar = context.state.gas_cells.velocity_z_peculiar[cell_index],
          .pressure_comoving = pressure,
      };
      if (cell_index >= m_geometry_row_by_dense_row.size()) {
        throw std::out_of_range("hydro callback dense row is outside Cartesian geometry map");
      }
      m_conserved.storeCell(
          m_geometry_row_by_dense_row[cell_index],
          hydro::HydroCoreSolver::conservedFromPrimitive(primitive, k_gamma_adiabatic));
    }
    m_last_remote_ghost_boundary_report = refreshRemoteHydroBoundaryGhosts(
        context,
        particle_row_by_id);
    hydro::fillHydroBoundaryGhostCells(m_conserved, m_geometry, k_gamma_adiabatic);
    const std::uint32_t world_rank = static_cast<std::uint32_t>(m_mpi_context.worldRank());
    const HydroGhostConservedSnapshot ghost_conserved_snapshot = snapshotHydroGhostConservedCells(
        context.state, m_conserved, m_geometry, m_geometry_row_by_dense_row, world_rank);
    const hydro::HydroActiveSetView active_view = buildActiveFaceView(context.active_set.cell_indices);

    // Hydro runs after drift/force refresh while IntegratorState still carries
    // the committed step-begin epoch. TimelineStep is authoritative for the
    // step-end force/source epoch and for H in inverse code time.
    const double hydro_scale_factor = context.cosmology_background != nullptr
        ? context.timeline_step.scale_factor_end
        : context.integrator_state.current_scale_factor;
    const double hubble_rate_code = context.cosmology_background != nullptr
        ? context.timeline_step.hubble_end_code
        : 0.0;

    hydro::HydroUpdateContext update{
        .dt_code = context.integrator_state.dt_time_code,
        .scale_factor = std::max(1.0e-12, hydro_scale_factor),
        .hubble_rate_code = hubble_rate_code,
    };

    std::vector<double> metallicity(context.state.cells.size(), 0.0);
    std::vector<double> temperature(context.state.cells.size(), 0.0);
    std::vector<double> hydrogen_number_density(context.state.cells.size(), 0.0);
    const auto& dense_accel_x = m_gravity_callback.cellAccelX();
    const auto& dense_accel_y = m_gravity_callback.cellAccelY();
    const auto& dense_accel_z = m_gravity_callback.cellAccelZ();
    if (dense_accel_x.size() != context.state.cells.size() ||
        dense_accel_y.size() != context.state.cells.size() ||
        dense_accel_z.size() != context.state.cells.size()) {
      throw std::runtime_error("hydro callback gravity-cell acceleration lanes do not cover dense gas-cell rows");
    }
    m_ordered_cell_accel_x.assign(context.state.cells.size(), 0.0);
    m_ordered_cell_accel_y.assign(context.state.cells.size(), 0.0);
    m_ordered_cell_accel_z.assign(context.state.cells.size(), 0.0);
    for (std::size_t geometry_row = 0; geometry_row < m_dense_row_by_geometry_row.size(); ++geometry_row) {
      const std::uint32_t dense_row = m_dense_row_by_geometry_row[geometry_row];
      m_ordered_cell_accel_x[geometry_row] = dense_accel_x[dense_row];
      m_ordered_cell_accel_y[geometry_row] = dense_accel_y[dense_row];
      m_ordered_cell_accel_z[geometry_row] = dense_accel_z[dense_row];
    }
    hydro::HydroSourceContext source_context{
        .update = update,
        .gravity_accel_x_peculiar = m_ordered_cell_accel_x,
        .gravity_accel_y_peculiar = m_ordered_cell_accel_y,
        .gravity_accel_z_peculiar = m_ordered_cell_accel_z,
        .hydrogen_number_density_cgs = hydrogen_number_density,
        .metallicity_mass_fraction = metallicity,
        .temperature_k = temperature,
        .redshift = std::max(0.0, (update.scale_factor > 0.0 ? (1.0 / update.scale_factor - 1.0) : 0.0)),
    };

    hydro::ComovingGravityExpansionSource gravity_source;
    std::array<const hydro::HydroSourceTerm*, 1> sources{&gravity_source};
    m_solver.advancePatchActiveSetWithScratch(
        m_conserved,
        m_geometry,
        active_view,
        update,
        m_reconstruction,
        m_riemann_solver,
        sources,
        source_context,
        m_scratch,
        &m_primitive_cache,
        &m_last_hydro_profile);
    if (context.profiler_session != nullptr) {
      const hydro::HydroConservationReport& conservation = m_last_hydro_profile.conservation;
      context.profiler_session->recordEvent(core::RuntimeEvent{
          .event_kind = "hydro.conservation",
          .severity = core::RuntimeEventSeverity::kInfo,
          .subsystem = "hydro.godunov",
          .step_index = context.integrator_state.step_index,
          .simulation_time_code = context.integrator_state.current_time_code,
          .scale_factor = context.integrator_state.current_scale_factor,
          .message = "computed volume-integrated hydro conservation totals for the local update stage",
          .payload = {{"cell_count", std::to_string(conservation.cell_count)},
                      {"before_mass", std::to_string(conservation.before.mass)},
                      {"after_mass", std::to_string(conservation.after.mass)},
                      {"flux_delta_mass", std::to_string(conservation.flux_delta.mass)},
                      {"source_delta_mass", std::to_string(conservation.source_delta.mass)},
                      {"floor_delta_mass", std::to_string(conservation.floor_delta.mass)},
                      {"residual_mass", std::to_string(conservation.residual.mass)},
                      {"residual_momentum_x", std::to_string(conservation.residual.momentum_x)},
                      {"residual_momentum_y", std::to_string(conservation.residual.momentum_y)},
                      {"residual_momentum_z", std::to_string(conservation.residual.momentum_z)},
                      {"residual_total_energy", std::to_string(conservation.residual.total_energy)},
                      {"residual_internal_energy", std::to_string(conservation.residual.internal_energy)},
                      {"flux_delta_total_energy", std::to_string(conservation.flux_delta.total_energy)},
                      {"source_delta_total_energy", std::to_string(conservation.source_delta.total_energy)},
                      {"floor_delta_total_energy", std::to_string(conservation.floor_delta.total_energy)},
                      {"internal_energy_floor_count", std::to_string(conservation.internal_energy_floor_count)},
                      {"mass_tolerance", "1e-10"},
                      {"momentum_tolerance", "1e-10"},
                      {"energy_tolerance", "1e-10"}}});
    }

    const HydroConservativeGhostSyncReport ghost_sync_report = restoreHydroGhostConservedCells(
        m_conserved, ghost_conserved_snapshot, world_rank);
    std::vector<parallel::HydroConservativeFluxCorrectionRecord> global_flux_corrections;
    if (all_ranks_have_hydro_cells) {
      global_flux_corrections =
          parallel::executeBlockingHydroConservativeFluxCorrectionExchange(
              m_mpi_context, ghost_sync_report.correction_records, context.integrator_state.step_index);
    } else if (!ghost_sync_report.correction_records.empty()) {
      throw std::runtime_error(
          "hydro conservative flux correction produced MPI correction records while at least one rank has no hydro cells");
    }
    std::size_t applied_flux_corrections = 0;
    double applied_flux_delta_l1 = 0.0;
    for (const parallel::HydroConservativeFluxCorrectionRecord& correction : global_flux_corrections) {
      if (correction.owner_rank != static_cast<int>(world_rank)) {
        continue;
      }
      const auto correction_row = context.state.rowForGasCellId(correction.gas_cell_id);
      if (!correction_row.has_value()) {
        throw std::runtime_error("hydro conservative flux correction targeted an unknown local gas_cell_id");
      }
      const std::uint32_t cell_index = *correction_row;
      const std::uint32_t owner_rank = internal::gasCellOwnerRankForLocalRow(
          context.state,
          cell_index,
          particle_row_by_id,
          "hydro conservative flux correction");
      if (owner_rank != world_rank) {
        throw std::runtime_error("hydro conservative flux correction targeted a non-authoritative local ghost cell");
      }
      if (cell_index >= m_geometry_row_by_dense_row.size()) {
        throw std::out_of_range("hydro flux correction dense row is outside Cartesian geometry map");
      }
      const std::uint32_t geometry_row = m_geometry_row_by_dense_row[cell_index];
      hydro::HydroConservedState owned = m_conserved.loadCell(geometry_row);
      owned.mass_density_comoving += correction.delta_mass_density_comoving;
      owned.momentum_density_x_comoving += correction.delta_momentum_density_x_comoving;
      owned.momentum_density_y_comoving += correction.delta_momentum_density_y_comoving;
      owned.momentum_density_z_comoving += correction.delta_momentum_density_z_comoving;
      owned.total_energy_density_comoving += correction.delta_total_energy_density_comoving;
      m_conserved.storeCell(geometry_row, owned);
      ++applied_flux_corrections;
      applied_flux_delta_l1 += std::abs(correction.delta_mass_density_comoving);
      applied_flux_delta_l1 += std::abs(correction.delta_momentum_density_x_comoving);
      applied_flux_delta_l1 += std::abs(correction.delta_momentum_density_y_comoving);
      applied_flux_delta_l1 += std::abs(correction.delta_momentum_density_z_comoving);
      applied_flux_delta_l1 += std::abs(correction.delta_total_energy_density_comoving);
    }
    if (context.profiler_session != nullptr) {
      context.profiler_session->recordEvent(core::RuntimeEvent{
          .event_kind = "hydro.conservative_ghost_sync",
          .severity = core::RuntimeEventSeverity::kInfo,
          .subsystem = "hydro.godunov",
          .step_index = context.integrator_state.step_index,
          .simulation_time_code = context.integrator_state.current_time_code,
          .scale_factor = context.integrator_state.current_scale_factor,
          .message = "restored imported hydro ghosts and applied conservative owner-side boundary corrections",
          .payload = {{"restored_ghost_cells", std::to_string(ghost_sync_report.restored_ghost_cells)},
                      {"local_correction_records", std::to_string(ghost_sync_report.correction_records.size())},
                      {"global_correction_records", std::to_string(global_flux_corrections.size())},
                      {"applied_flux_corrections", std::to_string(applied_flux_corrections)},
                      {"rejected_remote_delta_l1", std::to_string(ghost_sync_report.rejected_remote_delta_l1)},
                      {"applied_flux_delta_l1", std::to_string(applied_flux_delta_l1)},
                      {"imported_mpi_ghosts", std::to_string(m_last_remote_ghost_boundary_report.imported_ghosts)},
                      {"requested_remote_cells", std::to_string(m_last_remote_ghost_boundary_report.requested_cells)},
                      {"received_remote_cells", std::to_string(m_last_remote_ghost_boundary_report.received_cells)},
                      {"remote_interface_faces", std::to_string(m_last_remote_ghost_boundary_report.interface_faces)},
                      {"remote_stale_or_invalid_payloads", std::to_string(m_last_remote_ghost_boundary_report.stale_or_invalid_payloads)},
                      {"remote_request_bytes", std::to_string(m_last_remote_ghost_boundary_report.request_bytes)},
                      {"remote_payload_bytes", std::to_string(m_last_remote_ghost_boundary_report.payload_bytes)}}});
    }

    for (std::size_t cell_index = 0; cell_index < context.state.cells.size(); ++cell_index) {
      const std::uint32_t gas_cell_owner_rank = internal::gasCellOwnerRankForLocalRow(
          context.state,
          static_cast<std::uint32_t>(cell_index),
          particle_row_by_id,
          "hydro callback primitive store");
      if (gas_cell_owner_rank != world_rank) {
        continue;
      }
      if (cell_index >= m_geometry_row_by_dense_row.size()) {
        throw std::out_of_range("hydro primitive store dense row is outside Cartesian geometry map");
      }
      const std::uint32_t geometry_row = m_geometry_row_by_dense_row[cell_index];
      const hydro::HydroPrimitiveState primitive =
          hydro::HydroCoreSolver::primitiveFromConserved(m_conserved.loadCell(geometry_row), k_gamma_adiabatic);
      context.state.gas_cells.density_code[cell_index] = primitive.rho_comoving;
      context.state.gas_cells.pressure_code[cell_index] = primitive.pressure_comoving;
      context.state.gas_cells.internal_energy_code[cell_index] =
          primitive.pressure_comoving / ((k_gamma_adiabatic - 1.0) * std::max(primitive.rho_comoving, k_density_floor));
      context.state.cells.mass_code[cell_index] = primitive.rho_comoving * m_geometry.cell_volume_comoving;
      context.state.gas_cells.velocity_x_peculiar[cell_index] = primitive.vel_x_peculiar;
      context.state.gas_cells.velocity_y_peculiar[cell_index] = primitive.vel_y_peculiar;
      context.state.gas_cells.velocity_z_peculiar[cell_index] = primitive.vel_z_peculiar;
    }

    internal::synchronizeParentParticleCompatibilityMirrors(
        context.state, world_rank, "hydro callback parent compatibility mirror");
  }

 private:
  [[nodiscard]] static hydro::HydroBoundaryKind hydroBoundaryKindFromModePolicy(
      core::BoundaryCondition boundary_condition) {
    switch (boundary_condition) {
      case core::BoundaryCondition::kPeriodic:
        return hydro::HydroBoundaryKind::kPeriodic;
      case core::BoundaryCondition::kOpen:
        return hydro::HydroBoundaryKind::kOpen;
      case core::BoundaryCondition::kReflective:
        return hydro::HydroBoundaryKind::kReflective;
    }
    return hydro::HydroBoundaryKind::kOpen;
  }

  [[nodiscard]] std::array<double, 3> ghostCenterForBoundaryFace(
      const core::SimulationState& state,
      const hydro::HydroGhostCell& ghost) const {
    const std::uint32_t dense_row = m_dense_row_by_geometry_row.at(ghost.owner_real_cell);
    std::array<double, 3> center{
        state.cells.center_x_comoving[dense_row],
        state.cells.center_y_comoving[dense_row],
        state.cells.center_z_comoving[dense_row]};
    const double sign = ghost.side == hydro::HydroFaceSide::kLower ? -1.0 : 1.0;
    switch (ghost.axis) {
      case hydro::HydroFaceAxis::kX:
        center[0] += sign * m_geometry.cell_width_x_comoving;
        break;
      case hydro::HydroFaceAxis::kY:
        center[1] += sign * m_geometry.cell_width_y_comoving;
        break;
      case hydro::HydroFaceAxis::kZ:
        center[2] += sign * m_geometry.cell_width_z_comoving;
        break;
    }
    return center;
  }

  [[nodiscard]] bool isBoundaryGeometryRow(std::size_t geometry_row) const {
    const auto ijk = m_geometry.cellIjk(geometry_row);
    return ijk[0] == 0U || ijk[0] + 1U == m_geometry.nx ||
        ijk[1] == 0U || ijk[1] + 1U == m_geometry.ny ||
        ijk[2] == 0U || ijk[2] + 1U == m_geometry.nz;
  }

  [[nodiscard]] std::vector<HydroBoundaryCellAdvertisement> buildLocalHydroBoundaryAdvertisements(
      const core::SimulationState& state,
      const std::unordered_map<std::uint64_t, std::uint32_t>& particle_row_by_id,
      std::uint64_t hydro_sync_epoch,
      std::uint64_t decomposition_epoch) const {
    std::vector<HydroBoundaryCellAdvertisement> records;
    records.reserve(state.cells.size());
    const std::uint32_t world_rank = static_cast<std::uint32_t>(m_mpi_context.worldRank());
    for (std::uint32_t dense_row = 0; dense_row < state.cells.size(); ++dense_row) {
      if (internal::gasCellOwnerRankForLocalRow(state, dense_row, particle_row_by_id, "hydro boundary advertisement") !=
          world_rank) {
        continue;
      }
      if (dense_row >= m_geometry_row_by_dense_row.size()) {
        throw std::out_of_range("hydro boundary advertisement dense row is outside geometry map");
      }
      const std::size_t geometry_row = m_geometry_row_by_dense_row[dense_row];
      if (!isBoundaryGeometryRow(geometry_row)) {
        continue;
      }
      const core::GasCellIdentityRecord& identity = internal::gasCellIdentityRecordForLocalRow(
          state,
          dense_row,
          "hydro boundary advertisement");
      records.push_back(HydroBoundaryCellAdvertisement{
          .gas_cell_id = identity.gas_cell_id,
          .owner_rank = static_cast<int>(world_rank),
          .hydro_sync_epoch = hydro_sync_epoch,
          .decomposition_epoch = decomposition_epoch,
          .center_x_comoving = state.cells.center_x_comoving[dense_row],
          .center_y_comoving = state.cells.center_y_comoving[dense_row],
          .center_z_comoving = state.cells.center_z_comoving[dense_row]});
    }
    return records;
  }

  [[nodiscard]] std::optional<HydroBoundaryCellAdvertisement> findRemoteBoundaryAdvertisement(
      std::span<const HydroBoundaryCellAdvertisement> advertisements,
      const std::array<double, 3>& owner_center,
      hydro::HydroFaceAxis axis,
      hydro::HydroFaceSide side,
      std::uint64_t hydro_sync_epoch) const {
    const double tolerance = 1.0e-8 * std::max({
        1.0,
        std::abs(m_geometry.cell_width_x_comoving),
        std::abs(m_geometry.cell_width_y_comoving),
        std::abs(m_geometry.cell_width_z_comoving)});
    std::optional<HydroBoundaryCellAdvertisement> match;
    double best_normal_distance = std::numeric_limits<double>::infinity();
    for (const HydroBoundaryCellAdvertisement& record : advertisements) {
      if (record.owner_rank == m_mpi_context.worldRank()) {
        continue;
      }
      if (record.hydro_sync_epoch != hydro_sync_epoch) {
        continue;
      }
      double normal_distance = 0.0;
      bool transverse_match = false;
      switch (axis) {
        case hydro::HydroFaceAxis::kX:
          normal_distance = record.center_x_comoving - owner_center[0];
          transverse_match = nearlySameCoordinate(record.center_y_comoving, owner_center[1], tolerance) &&
              nearlySameCoordinate(record.center_z_comoving, owner_center[2], tolerance);
          break;
        case hydro::HydroFaceAxis::kY:
          normal_distance = record.center_y_comoving - owner_center[1];
          transverse_match = nearlySameCoordinate(record.center_x_comoving, owner_center[0], tolerance) &&
              nearlySameCoordinate(record.center_z_comoving, owner_center[2], tolerance);
          break;
        case hydro::HydroFaceAxis::kZ:
          normal_distance = record.center_z_comoving - owner_center[2];
          transverse_match = nearlySameCoordinate(record.center_x_comoving, owner_center[0], tolerance) &&
              nearlySameCoordinate(record.center_y_comoving, owner_center[1], tolerance);
          break;
      }
      const bool correct_side = side == hydro::HydroFaceSide::kUpper
          ? normal_distance > tolerance
          : normal_distance < -tolerance;
      if (!transverse_match || !correct_side) {
        continue;
      }
      const double abs_distance = std::abs(normal_distance);
      if (abs_distance + tolerance < best_normal_distance) {
        best_normal_distance = abs_distance;
        match = record;
      } else if (std::abs(abs_distance - best_normal_distance) <= tolerance) {
        throw std::runtime_error("hydro remote boundary discovery found duplicate remote cells for one face");
      }
    }
    return match;
  }

  [[nodiscard]] HydroRemoteGhostBoundaryReport refreshRemoteHydroBoundaryGhosts(
      core::StepContext& context,
      const std::unordered_map<std::uint64_t, std::uint32_t>& particle_row_by_id) {
    HydroRemoteGhostBoundaryReport report;
    if (!m_mpi_context.isEnabled() || m_mpi_context.worldSize() <= 1) {
      return report;
    }
    const std::uint64_t hydro_sync_epoch =
        context.integrator_state.step_index * core::integrationStageCount() +
        core::integrationStageIndex(context.stage) + 1U;
    const std::uint64_t decomposition_epoch = context.state.gasCellIdentityGeneration();
    const auto local_advertisements = buildLocalHydroBoundaryAdvertisements(
        context.state,
        particle_row_by_id,
        hydro_sync_epoch,
        decomposition_epoch);
    const auto advertisements = exchangeHydroBoundaryAdvertisements(m_mpi_context, local_advertisements);
    report.boundary_metadata_sent = local_advertisements.size();
    report.boundary_metadata_received = advertisements.size();

    std::vector<parallel::HydroGhostCellRequest> requests;
    std::unordered_map<std::uint64_t, std::size_t> ghost_slot_by_face_key;
    m_remote_non_authority_faces.clear();
    for (std::size_t face_index = 0; face_index < m_geometry.faces.size(); ++face_index) {
      hydro::HydroFace& face = m_geometry.faces[face_index];
      if (face.ghost_cell_slot != hydro::k_invalid_ghost_cell_slot ||
          face.neighbor_cell == hydro::k_invalid_cell_index ||
          face.owner_cell >= m_dense_row_by_geometry_row.size() ||
          face.neighbor_cell >= m_dense_row_by_geometry_row.size()) {
        continue;
      }
      const std::uint32_t owner_dense_row = m_dense_row_by_geometry_row[face.owner_cell];
      const std::uint32_t neighbor_dense_row = m_dense_row_by_geometry_row[face.neighbor_cell];
      const std::uint32_t owner_rank = internal::gasCellOwnerRankForLocalRow(
          context.state,
          owner_dense_row,
          particle_row_by_id,
          "hydro internal remote boundary discovery");
      const std::uint32_t neighbor_rank = internal::gasCellOwnerRankForLocalRow(
          context.state,
          neighbor_dense_row,
          particle_row_by_id,
          "hydro internal remote boundary discovery");
      if (owner_rank == neighbor_rank) {
        continue;
      }
      if (owner_rank != static_cast<std::uint32_t>(m_mpi_context.worldRank())) {
        m_remote_non_authority_faces.insert(face_index);
        continue;
      }
      const core::GasCellIdentityRecord& owner_identity = internal::gasCellIdentityRecordForLocalRow(
          context.state,
          owner_dense_row,
          "hydro internal remote boundary discovery");
      const core::GasCellIdentityRecord& remote_identity = internal::gasCellIdentityRecordForLocalRow(
          context.state,
          neighbor_dense_row,
          "hydro internal remote boundary discovery");
      const std::uint64_t face_key =
          hydroInterfaceFaceKey(owner_identity.gas_cell_id, remote_identity.gas_cell_id, face.axis);
      const std::size_t ghost_slot = m_geometry.ghost_cells.size();
      const std::size_t ghost_cell = m_geometry.cellCount() + ghost_slot;
      m_geometry.ghost_cells.push_back(hydro::HydroGhostCell{
          .owner_real_cell = face.owner_cell,
          .source_real_cell = face.owner_cell,
          .ghost_cell = ghost_cell,
          .ghost_slot = ghost_slot,
          .boundary_kind = hydro::HydroBoundaryKind::kImportedMpi,
          .axis = face.axis,
          .side = hydro::HydroFaceSide::kUpper,
          .mutation_rights = hydro::HydroGhostMutationRights::kReadOnlyImported,
          .origin_gas_cell_id = remote_identity.gas_cell_id,
          .origin_rank = static_cast<int>(neighbor_rank),
          .hydro_sync_epoch = hydro_sync_epoch,
          .decomposition_epoch = decomposition_epoch});
      face.neighbor_cell = ghost_cell;
      face.neighbor_plus_cell = hydro::k_invalid_cell_index;
      face.ghost_cell_slot = ghost_slot;
      if (!ghost_slot_by_face_key.emplace(face_key, ghost_slot).second) {
        throw std::runtime_error("hydro internal remote boundary discovery produced duplicate face keys");
      }
      requests.push_back(parallel::HydroGhostCellRequest{
          .descriptor = parallel::HydroGhostCellDescriptor{
              .gas_cell_id = remote_identity.gas_cell_id,
              .owner_rank = static_cast<int>(neighbor_rank),
              .consumer_rank = m_mpi_context.worldRank(),
              .hydro_sync_epoch = hydro_sync_epoch,
              .decomposition_epoch = decomposition_epoch,
              .boundary_state_only = true},
          .face_key = face_key,
          .axis = hydroAxisCode(face.axis),
          .side = 1U});
    }
    if (m_conserved.size() < m_geometry.totalCellStorageCount()) {
      m_conserved.resize(m_geometry.totalCellStorageCount());
    }
    for (hydro::HydroGhostCell& ghost : m_geometry.ghost_cells) {
      if (ghost.boundary_kind == hydro::HydroBoundaryKind::kImportedMpi &&
          ghost.hydro_sync_epoch != hydro_sync_epoch) {
        ghost.boundary_kind = hydroBoundaryKindFromModePolicy(m_mode_policy.hydro_boundary);
        ghost.mutation_rights = hydro::HydroGhostMutationRights::kWritablePhysicalBoundaryScratch;
        ghost.origin_gas_cell_id = 0U;
        ghost.origin_rank = -1;
        ghost.hydro_sync_epoch = 0U;
        ghost.decomposition_epoch = 0U;
      }
      const std::uint32_t owner_dense_row = m_dense_row_by_geometry_row.at(ghost.owner_real_cell);
      if (internal::gasCellOwnerRankForLocalRow(
              context.state,
              owner_dense_row,
              particle_row_by_id,
              "hydro remote boundary discovery") != static_cast<std::uint32_t>(m_mpi_context.worldRank())) {
        continue;
      }
      const std::array<double, 3> owner_center{
          context.state.cells.center_x_comoving[owner_dense_row],
          context.state.cells.center_y_comoving[owner_dense_row],
          context.state.cells.center_z_comoving[owner_dense_row]};
      const auto remote = findRemoteBoundaryAdvertisement(
          advertisements,
          owner_center,
          ghost.axis,
          ghost.side,
          hydro_sync_epoch);
      if (!remote.has_value()) {
        continue;
      }
      const core::GasCellIdentityRecord& local_identity = internal::gasCellIdentityRecordForLocalRow(
          context.state,
          owner_dense_row,
          "hydro remote boundary discovery");
      const std::uint64_t face_key = hydroInterfaceFaceKey(local_identity.gas_cell_id, remote->gas_cell_id, ghost.axis);
      if (!ghost_slot_by_face_key.emplace(face_key, ghost.ghost_slot).second) {
        throw std::runtime_error("hydro remote boundary discovery produced duplicate face keys");
      }
      ghost.boundary_kind = hydro::HydroBoundaryKind::kImportedMpi;
      ghost.mutation_rights = hydro::HydroGhostMutationRights::kReadOnlyImported;
      ghost.source_real_cell = ghost.owner_real_cell;
      ghost.origin_gas_cell_id = remote->gas_cell_id;
      ghost.origin_rank = remote->owner_rank;
      ghost.hydro_sync_epoch = hydro_sync_epoch;
      ghost.decomposition_epoch = decomposition_epoch;
      requests.push_back(parallel::HydroGhostCellRequest{
          .descriptor = parallel::HydroGhostCellDescriptor{
              .gas_cell_id = remote->gas_cell_id,
              .owner_rank = remote->owner_rank,
              .consumer_rank = m_mpi_context.worldRank(),
              .hydro_sync_epoch = hydro_sync_epoch,
              .decomposition_epoch = decomposition_epoch,
              .boundary_state_only = true},
          .face_key = face_key,
          .axis = hydroAxisCode(ghost.axis),
          .side = hydroSideCode(ghost.side)});
    }

    const auto global_requests = parallel::executeBlockingHydroGhostCellRequestExchange(
        m_mpi_context,
        requests,
        hydro_sync_epoch);
    std::vector<parallel::HydroGhostCellPayloadRecord> local_payloads;
    for (const parallel::HydroGhostCellRequest& request : global_requests) {
      if (request.descriptor.owner_rank != m_mpi_context.worldRank()) {
        continue;
      }
      const auto row = context.state.rowForGasCellId(request.descriptor.gas_cell_id);
      if (!row.has_value()) {
        throw std::runtime_error("hydro ghost request targeted an unknown owner gas_cell_id");
      }
      if (internal::gasCellOwnerRankForLocalRow(
              context.state,
              *row,
              particle_row_by_id,
              "hydro ghost payload response") != static_cast<std::uint32_t>(m_mpi_context.worldRank())) {
        throw std::runtime_error("hydro ghost request targeted a local non-authoritative gas cell");
      }
      if (*row >= m_geometry_row_by_dense_row.size()) {
        throw std::out_of_range("hydro ghost payload response dense row is outside geometry map");
      }
      const hydro::HydroConservedState state = m_conserved.loadCell(m_geometry_row_by_dense_row[*row]);
      local_payloads.push_back(parallel::HydroGhostCellPayloadRecord{
          .descriptor = request.descriptor,
          .face_key = request.face_key,
          .mass_density_comoving = state.mass_density_comoving,
          .momentum_density_x_comoving = state.momentum_density_x_comoving,
          .momentum_density_y_comoving = state.momentum_density_y_comoving,
          .momentum_density_z_comoving = state.momentum_density_z_comoving,
          .total_energy_density_comoving = state.total_energy_density_comoving});
    }
    const auto global_payloads = parallel::executeBlockingHydroGhostCellPayloadExchange(
        m_mpi_context,
        local_payloads,
        hydro_sync_epoch);

    std::unordered_set<std::uint64_t> received_face_keys;
    for (const parallel::HydroGhostCellPayloadRecord& payload : global_payloads) {
      if (payload.descriptor.consumer_rank != m_mpi_context.worldRank()) {
        continue;
      }
      const auto slot_it = ghost_slot_by_face_key.find(payload.face_key);
      if (slot_it == ghost_slot_by_face_key.end()) {
        ++report.stale_or_invalid_payloads;
        continue;
      }
      if (!received_face_keys.insert(payload.face_key).second) {
        throw std::runtime_error("hydro ghost payload exchange returned duplicate payload for one face");
      }
      hydro::HydroGhostCell& ghost = m_geometry.ghost_cells.at(slot_it->second);
      if (payload.descriptor.gas_cell_id != ghost.origin_gas_cell_id ||
          payload.descriptor.owner_rank != ghost.origin_rank ||
          payload.descriptor.hydro_sync_epoch != ghost.hydro_sync_epoch ||
          payload.descriptor.decomposition_epoch != ghost.decomposition_epoch) {
        ++report.stale_or_invalid_payloads;
        continue;
      }
      m_conserved.storeCell(ghost.ghost_cell, hydro::HydroConservedState{
          .mass_density_comoving = payload.mass_density_comoving,
          .momentum_density_x_comoving = payload.momentum_density_x_comoving,
          .momentum_density_y_comoving = payload.momentum_density_y_comoving,
          .momentum_density_z_comoving = payload.momentum_density_z_comoving,
          .total_energy_density_comoving = payload.total_energy_density_comoving});
    }
    if (received_face_keys.size() != requests.size()) {
      throw std::runtime_error("hydro remote boundary ghost refresh did not receive every requested ghost cell");
    }
    report.requested_cells = requests.size();
    report.received_cells = received_face_keys.size();
    report.imported_ghosts = received_face_keys.size();
    report.interface_faces = requests.size();
    report.request_bytes = static_cast<std::uint64_t>(requests.size() * sizeof(parallel::HydroGhostCellRequest));
    report.payload_bytes =
        static_cast<std::uint64_t>(local_payloads.size() * sizeof(parallel::HydroGhostCellPayloadRecord));
    return report;
  }

  void runProductionAmrHydroPath(core::StepContext& context) {
    // Hydro runs after drift/force refresh while IntegratorState still carries
    // the committed step-begin epoch. TimelineStep is authoritative for the
    // step-end force/source epoch and for H in inverse code time.
    const double hydro_scale_factor = context.cosmology_background != nullptr
        ? context.timeline_step.scale_factor_end
        : context.integrator_state.current_scale_factor;
    const double hubble_rate_code = context.cosmology_background != nullptr
        ? context.timeline_step.hubble_end_code
        : 0.0;

    const hydro::HydroUpdateContext update{
        .dt_code = context.integrator_state.dt_time_code,
        .scale_factor = std::max(1.0e-12, hydro_scale_factor),
        .hubble_rate_code = hubble_rate_code,
    };
    std::vector<double> metallicity(context.state.cells.size(), 0.0);
    std::vector<double> temperature(context.state.cells.size(), 0.0);
    std::vector<double> hydrogen_number_density(context.state.cells.size(), 0.0);
    const hydro::HydroSourceContext source_context{
        .update = update,
        .gravity_accel_x_peculiar = m_gravity_callback.cellAccelX(),
        .gravity_accel_y_peculiar = m_gravity_callback.cellAccelY(),
        .gravity_accel_z_peculiar = m_gravity_callback.cellAccelZ(),
        .hydrogen_number_density_cgs = hydrogen_number_density,
        .metallicity_mass_fraction = metallicity,
        .temperature_k = temperature,
        .redshift = std::max(0.0, (update.scale_factor > 0.0 ? (1.0 / update.scale_factor - 1.0) : 0.0)),
    };
    hydro::ComovingGravityExpansionSource gravity_source;
    std::array<const hydro::HydroSourceTerm*, 1> sources{&gravity_source};
    const amr::ProductionAmrHydroOptions amr_options{
        .physical_boundary_kind = hydroBoundaryKindFromModePolicy(m_mode_policy.hydro_boundary),
        .adiabatic_index = k_gamma_adiabatic,
        .density_floor = k_density_floor,
        .pressure_floor = k_pressure_floor,
        .state_time_code = context.integrator_state.current_time_code,
        .ghost_fill_time_code = context.integrator_state.current_time_code};
    amr::ProductionAmrHydroDiagnostics amr_diagnostics;
    std::size_t remote_patch_count = 0;
    std::size_t remote_flux_register_count = 0;
    std::size_t inbound_flux_register_count = 0;
    parallel::DirectedAmrExchangeDiagnostics directed_amr_diagnostics{};
    if (m_mpi_context.isEnabled() && m_mpi_context.worldSize() > 1) {
      const auto local_patch_records =
          internal::buildMigrationAmrPatchPayloadRecords(
              context.state, m_mpi_context.worldRank());
      const auto local_cell_records =
          internal::buildMigrationAmrPatchCellPayloadRecords(
              context.state, m_mpi_context.worldRank());
      const auto directed_amr_exchange = parallel::executeBlockingDirectedAmrPatchPayloadExchange(
          m_mpi_context,
          local_patch_records,
          local_cell_records,
          context.integrator_state.step_index);
      directed_amr_diagnostics = directed_amr_exchange.diagnostics;

      std::unordered_map<std::uint64_t, int> owner_rank_by_patch_id;
      owner_rank_by_patch_id.reserve(local_patch_records.size() + directed_amr_exchange.patch_payloads_received.size());
      for (const auto& record : local_patch_records) {
        const auto [it, inserted] = owner_rank_by_patch_id.emplace(record.patch_id, record.owner_rank);
        if (!inserted && it->second != record.owner_rank) {
          throw std::runtime_error("directed AMR hydro found conflicting local patch ownership metadata");
        }
      }

      std::unordered_map<std::uint64_t, parallel::AmrPatchPayloadRecord> patch_record_by_id;
      patch_record_by_id.reserve(directed_amr_exchange.patch_payloads_received.size());
      for (const auto& record : directed_amr_exchange.patch_payloads_received) {
        if (record.owner_rank == m_mpi_context.worldRank()) {
          throw std::runtime_error("directed AMR hydro imported a local patch as a remote ghost");
        }
        const auto [it, inserted] = patch_record_by_id.emplace(record.patch_id, record);
        if (!inserted) {
          throw std::runtime_error("directed AMR hydro found duplicate remote patch payload records");
        }
        owner_rank_by_patch_id.emplace(record.patch_id, record.owner_rank);
      }
      std::unordered_map<std::uint64_t, std::vector<parallel::AmrPatchCellPayloadRecord>> cells_by_patch_id;
      for (const auto& record : directed_amr_exchange.patch_cell_payloads_received) {
        cells_by_patch_id[record.patch_id].push_back(record);
      }
      std::vector<amr::DistributedAmrRemotePatch> remote_patches;
      remote_patches.reserve(directed_amr_exchange.patch_payloads_received.size());
      for (const auto& patch_record : directed_amr_exchange.patch_payloads_received) {
        auto cells_it = cells_by_patch_id.find(patch_record.patch_id);
        if (cells_it == cells_by_patch_id.end() ||
            cells_it->second.size() != patch_record.cell_count) {
          throw std::runtime_error("directed AMR hydro remote patch cell payload coverage mismatch");
        }
        std::sort(
            cells_it->second.begin(),
            cells_it->second.end(),
            [](const auto& lhs, const auto& rhs) { return lhs.local_cell_offset < rhs.local_cell_offset; });
        amr::DistributedAmrRemotePatch remote;
        remote.patch = amr::PatchDescriptor{
            .patch_id = patch_record.patch_id,
            .parent_patch_id = patch_record.parent_patch_id,
            .level = static_cast<std::uint8_t>(patch_record.level),
            .morton_key = patch_record.morton_key,
            .origin_comov = {patch_record.origin_x_comoving, patch_record.origin_y_comoving, patch_record.origin_z_comoving},
            .extent_comov = {patch_record.extent_x_comoving, patch_record.extent_y_comoving, patch_record.extent_z_comoving},
            .cell_dims = {patch_record.cell_dim_x, patch_record.cell_dim_y, patch_record.cell_dim_z}};
        remote.owner_rank = patch_record.owner_rank;
        remote.ghost_hydro_epoch = patch_record.decomposition_epoch;
        remote.expected_ghost_hydro_epoch = patch_record.decomposition_epoch;
        remote.gas_cell_ids.resize(patch_record.cell_count);
        remote.conserved_cells.resize(patch_record.cell_count);
        for (std::size_t i = 0; i < cells_it->second.size(); ++i) {
          const auto& cell_record = cells_it->second[i];
          if (cell_record.local_cell_offset != i || cell_record.owner_rank != patch_record.owner_rank) {
            throw std::runtime_error("directed AMR hydro remote patch cell payload has stale offset or owner metadata");
          }
          remote.gas_cell_ids[i] = cell_record.gas_cell_id;
          remote.conserved_cells[i] = hydro::HydroCoreSolver::conservedFromPrimitive(
              hydro::HydroPrimitiveState{
                  .rho_comoving = cell_record.density_code,
                  .vel_x_peculiar = cell_record.velocity_x_peculiar,
                  .vel_y_peculiar = cell_record.velocity_y_peculiar,
                  .vel_z_peculiar = cell_record.velocity_z_peculiar,
                  .pressure_comoving = cell_record.pressure_code},
              k_gamma_adiabatic);
        }
        remote_patches.push_back(std::move(remote));
      }
      remote_patch_count = remote_patches.size();
      std::vector<amr::FluxRegisterEntry> outbound_remote_entries;
      if (amr::hasProductionAmrHydroCoverage(context.state)) {
        amr_diagnostics = amr::advanceDistributedProductionAmrHydro(
            context.state,
            context.active_set.cell_indices,
            update,
            source_context,
            m_solver,
            m_riemann_solver,
            sources,
            amr_options,
            amr::DistributedAmrHydroExchange{
                .local_rank = m_mpi_context.worldRank(),
                .ghost_hydro_epoch = context.state.gasCellIdentityGeneration(),
                .expected_ghost_hydro_epoch = context.state.gasCellIdentityGeneration(),
                .remote_patches = remote_patches,
                .outbound_remote_flux_registers = &outbound_remote_entries});
      } else if (context.state.cells.size() != 0U || context.state.patches.size() != 0U) {
        throw std::runtime_error(
            "distributed AMR hydro rank has partial local patch coverage; only an explicitly empty rank may skip advancement");
      }
      std::vector<parallel::AmrFluxRegisterPayloadRecord> local_flux_payloads;
      local_flux_payloads.reserve(outbound_remote_entries.size());
      for (const amr::FluxRegisterEntry& entry : outbound_remote_entries) {
        const auto owner_it = owner_rank_by_patch_id.find(entry.coarse_patch_id);
        if (owner_it == owner_rank_by_patch_id.end()) {
          throw std::runtime_error("distributed AMR hydro remote reflux entry has unknown coarse patch owner");
        }
        parallel::AmrFluxRegisterPayloadRecord record;
        record.register_key = entry.register_key;
        record.coarse_patch_id = entry.coarse_patch_id;
        record.coarse_gas_cell_id = entry.coarse_gas_cell_id;
        record.coarse_cell_index = static_cast<std::uint64_t>(entry.coarse_cell_index);
        record.level = entry.level;
        record.axis = hydroAxisCode(entry.axis);
        record.orientation = hydroSideCode(entry.orientation);
        record.source_rank = m_mpi_context.worldRank();
        record.owner_rank = owner_it->second;
        record.gas_cell_identity_generation = context.state.gasCellIdentityGeneration();
        record.patch_geometry_generation = context.state.cellIndexGeneration();
        record.coarse_mass_flux_code = entry.coarse_face_flux_code.mass_code;
        record.coarse_momentum_x_flux_code = entry.coarse_face_flux_code.momentum_x_code;
        record.coarse_momentum_y_flux_code = entry.coarse_face_flux_code.momentum_y_code;
        record.coarse_momentum_z_flux_code = entry.coarse_face_flux_code.momentum_z_code;
        record.coarse_total_energy_flux_code = entry.coarse_face_flux_code.total_energy_code;
        record.fine_mass_flux_code = entry.fine_face_flux_code.mass_code;
        record.fine_momentum_x_flux_code = entry.fine_face_flux_code.momentum_x_code;
        record.fine_momentum_y_flux_code = entry.fine_face_flux_code.momentum_y_code;
        record.fine_momentum_z_flux_code = entry.fine_face_flux_code.momentum_z_code;
        record.fine_total_energy_flux_code = entry.fine_face_flux_code.total_energy_code;
        record.face_area_comov = entry.face_area_comov;
        record.coarse_area_comov = entry.coarse_area_comov;
        record.fine_area_comov = entry.fine_area_comov;
        record.dt_code = entry.dt_code;
        record.coarse_face_count = entry.coarse_face_count;
        record.fine_face_count = entry.fine_face_count;
        parallel::validateAmrFluxRegisterPayloadRecord(record);
        local_flux_payloads.push_back(record);
      }
      const auto global_flux_payloads = parallel::executeBlockingAmrFluxRegisterPayloadExchange(
          m_mpi_context, local_flux_payloads, context.integrator_state.step_index);
      directed_amr_diagnostics.outbound_reflux_count = static_cast<std::uint64_t>(local_flux_payloads.size());
      directed_amr_diagnostics.inbound_reflux_count = static_cast<std::uint64_t>(global_flux_payloads.size());
      directed_amr_diagnostics.directed_flux_records_sent = static_cast<std::uint64_t>(local_flux_payloads.size());
      directed_amr_diagnostics.directed_flux_records_received = static_cast<std::uint64_t>(global_flux_payloads.size());
      directed_amr_diagnostics.flux_payload_bytes =
          (static_cast<std::uint64_t>(local_flux_payloads.size()) + static_cast<std::uint64_t>(global_flux_payloads.size())) *
          sizeof(parallel::AmrFluxRegisterPayloadRecord);
      std::vector<amr::FluxRegisterEntry> inbound_entries;
      std::unordered_set<std::uint64_t> inbound_keys;
      for (const auto& payload : global_flux_payloads) {
        if (payload.owner_rank != m_mpi_context.worldRank()) {
          continue;
        }
        if (!inbound_keys.insert(payload.register_key).second) {
          throw std::runtime_error("distributed AMR hydro received duplicate reflux register key for local owner");
        }
        inbound_entries.push_back(amr::FluxRegisterEntry{
            .register_key = payload.register_key,
            .coarse_patch_id = payload.coarse_patch_id,
            .coarse_gas_cell_id = payload.coarse_gas_cell_id,
            .coarse_cell_index = static_cast<std::size_t>(payload.coarse_cell_index),
            .level = payload.level,
            .axis = payload.axis == 0U ? hydro::HydroFaceAxis::kX : (payload.axis == 1U ? hydro::HydroFaceAxis::kY : hydro::HydroFaceAxis::kZ),
            .orientation = payload.orientation == 0U ? hydro::HydroFaceSide::kLower : hydro::HydroFaceSide::kUpper,
            .coarse_face_flux_code = amr::ConservedState{
                .mass_code = payload.coarse_mass_flux_code,
                .momentum_x_code = payload.coarse_momentum_x_flux_code,
                .momentum_y_code = payload.coarse_momentum_y_flux_code,
                .momentum_z_code = payload.coarse_momentum_z_flux_code,
                .total_energy_code = payload.coarse_total_energy_flux_code},
            .fine_face_flux_code = amr::ConservedState{
                .mass_code = payload.fine_mass_flux_code,
                .momentum_x_code = payload.fine_momentum_x_flux_code,
                .momentum_y_code = payload.fine_momentum_y_flux_code,
                .momentum_z_code = payload.fine_momentum_z_flux_code,
                .total_energy_code = payload.fine_total_energy_flux_code},
            .face_area_comov = payload.face_area_comov,
            .coarse_area_comov = payload.coarse_area_comov,
            .fine_area_comov = payload.fine_area_comov,
            .dt_code = payload.dt_code,
            .coarse_face_count = payload.coarse_face_count,
            .fine_face_count = payload.fine_face_count});
      }
      inbound_flux_register_count = inbound_entries.size();
      remote_flux_register_count = local_flux_payloads.size();
      if (!inbound_entries.empty()) {
        const auto local_descriptors = amr::buildProductionAmrPatchDescriptors(context.state);
        const amr::RefluxDiagnostics remote_reflux = amr::applyFluxRegistersToSimulationState(
            context.state, inbound_entries, local_descriptors, k_gamma_adiabatic);
        amr_diagnostics.reflux.complete_register_count += remote_reflux.complete_register_count;
        amr_diagnostics.reflux.skipped_incomplete_register_count += remote_reflux.skipped_incomplete_register_count;
        amr_diagnostics.reflux.skipped_area_mismatch_count += remote_reflux.skipped_area_mismatch_count;
        amr_diagnostics.reflux.skipped_missing_target_count += remote_reflux.skipped_missing_target_count;
        amr_diagnostics.reflux.corrected_cells += remote_reflux.corrected_cells;
        amr_diagnostics.reflux.corrected_mass_code += remote_reflux.corrected_mass_code;
        amr_diagnostics.reflux.corrected_momentum_x_code += remote_reflux.corrected_momentum_x_code;
        amr_diagnostics.reflux.corrected_momentum_y_code += remote_reflux.corrected_momentum_y_code;
        amr_diagnostics.reflux.corrected_momentum_z_code += remote_reflux.corrected_momentum_z_code;
        amr_diagnostics.reflux.corrected_total_energy_code += remote_reflux.corrected_total_energy_code;
        amr_diagnostics.reflux.corrected_energy_code += remote_reflux.corrected_energy_code;
        amr_diagnostics.reflux.corrected_internal_energy_code += remote_reflux.corrected_internal_energy_code;
      }
    } else {
      amr_diagnostics = amr::advanceProductionAmrHydro(
          context.state,
          context.active_set.cell_indices,
          update,
          source_context,
          m_solver,
          m_riemann_solver,
          sources,
          amr_options);
    }
    internal::synchronizeParentParticleCompatibilityMirrors(
        context.state,
        m_current_world_rank,
        "production AMR hydro parent compatibility mirror");
    if (context.profiler_session != nullptr) {
      context.profiler_session->recordEvent(core::RuntimeEvent{
          .event_kind = "hydro.amr_production_stage",
          .severity = core::RuntimeEventSeverity::kInfo,
          .subsystem = "hydro.amr",
          .step_index = context.integrator_state.step_index,
          .simulation_time_code = context.integrator_state.current_time_code,
          .scale_factor = context.integrator_state.current_scale_factor,
          .message = "advanced production AMR hydro patches from SimulationState-authoritative gas rows",
          .payload = {{"patch_count", std::to_string(amr_diagnostics.patch_count)},
                      {"advanced_patch_count", std::to_string(amr_diagnostics.advanced_patch_count)},
                      {"active_cell_count", std::to_string(amr_diagnostics.active_cell_count)},
                      {"active_face_count", std::to_string(amr_diagnostics.active_face_count)},
                      {"flux_register_entry_count", std::to_string(amr_diagnostics.flux_register_entry_count)},
                      {"reflux_corrected_cells", std::to_string(amr_diagnostics.reflux.corrected_cells)},
                      {"reflux_corrected_mass", std::to_string(amr_diagnostics.reflux.corrected_mass_code)},
                      {"reflux_corrected_momentum_x", std::to_string(amr_diagnostics.reflux.corrected_momentum_x_code)},
                      {"reflux_corrected_momentum_y", std::to_string(amr_diagnostics.reflux.corrected_momentum_y_code)},
                      {"reflux_corrected_momentum_z", std::to_string(amr_diagnostics.reflux.corrected_momentum_z_code)},
                      {"reflux_corrected_total_energy", std::to_string(amr_diagnostics.reflux.corrected_total_energy_code)},
                      {"reflux_corrected_energy", std::to_string(amr_diagnostics.reflux.corrected_energy_code)},
                      {"reflux_corrected_internal_energy", std::to_string(amr_diagnostics.reflux.corrected_internal_energy_code)},
                      {"reflux_complete_register_count", std::to_string(amr_diagnostics.reflux.complete_register_count)},
                      {"reflux_skipped_incomplete_register_count", std::to_string(amr_diagnostics.reflux.skipped_incomplete_register_count)},
                      {"reflux_skipped_area_mismatch_count", std::to_string(amr_diagnostics.reflux.skipped_area_mismatch_count)},
                      {"reflux_skipped_missing_target_count", std::to_string(amr_diagnostics.reflux.skipped_missing_target_count)},
                      {"distributed_remote_patch_count", std::to_string(remote_patch_count)},
                      {"distributed_remote_flux_register_count", std::to_string(remote_flux_register_count)},
                      {"distributed_inbound_flux_register_count", std::to_string(inbound_flux_register_count)},
                      {"directed_amr_candidate_peer_count", std::to_string(directed_amr_diagnostics.candidate_peer_count)},
                      {"directed_amr_neighbor_peer_count", std::to_string(directed_amr_diagnostics.neighbor_peer_count)},
                      {"directed_amr_patch_descriptor_records_sent", std::to_string(directed_amr_diagnostics.directed_patch_descriptor_records_sent)},
                      {"directed_amr_patch_descriptor_records_received", std::to_string(directed_amr_diagnostics.directed_patch_descriptor_records_received)},
                      {"directed_amr_patch_cell_records_sent", std::to_string(directed_amr_diagnostics.directed_patch_cell_records_sent)},
                      {"directed_amr_patch_cell_records_received", std::to_string(directed_amr_diagnostics.directed_patch_cell_records_received)},
                      {"directed_amr_flux_records_sent", std::to_string(directed_amr_diagnostics.directed_flux_records_sent)},
                      {"directed_amr_flux_records_received", std::to_string(directed_amr_diagnostics.directed_flux_records_received)},
                      {"directed_amr_control_plane_bytes", std::to_string(directed_amr_diagnostics.control_plane_bytes)},
                      {"directed_amr_patch_descriptor_bytes", std::to_string(directed_amr_diagnostics.patch_descriptor_bytes)},
                      {"directed_amr_patch_cell_payload_bytes", std::to_string(directed_amr_diagnostics.patch_cell_payload_bytes)},
                      {"directed_amr_flux_payload_bytes", std::to_string(directed_amr_diagnostics.flux_payload_bytes)},
                      {"directed_amr_remote_interface_count", std::to_string(directed_amr_diagnostics.remote_interface_count)},
                      {"directed_amr_inbound_reflux_count", std::to_string(directed_amr_diagnostics.inbound_reflux_count)},
                      {"directed_amr_outbound_reflux_count", std::to_string(directed_amr_diagnostics.outbound_reflux_count)}}});
    }
  }

  struct HydroGeometryCacheKey {
    std::size_t cell_count = 0;
    std::size_t nx = 0;
    std::size_t ny = 0;
    std::size_t nz = 0;
    double origin_x_comoving = 0.0;
    double origin_y_comoving = 0.0;
    double origin_z_comoving = 0.0;
    double cell_width_x_comoving = 0.0;
    double cell_width_y_comoving = 0.0;
    double cell_width_z_comoving = 0.0;
    double dt_time_code = 0.0;
    core::BoundaryCondition hydro_boundary = core::BoundaryCondition::kPeriodic;
    std::uint64_t patch_signature = 0;
    std::uint64_t gas_cell_identity_generation = 0;
    std::uint64_t cell_index_generation = 0;
    std::uint64_t row_mapping_signature = 0;

    [[nodiscard]] bool operator==(const HydroGeometryCacheKey& rhs) const noexcept {
      return cell_count == rhs.cell_count &&
          nx == rhs.nx &&
          ny == rhs.ny &&
          nz == rhs.nz &&
          origin_x_comoving == rhs.origin_x_comoving &&
          origin_y_comoving == rhs.origin_y_comoving &&
          origin_z_comoving == rhs.origin_z_comoving &&
          cell_width_x_comoving == rhs.cell_width_x_comoving &&
          cell_width_y_comoving == rhs.cell_width_y_comoving &&
          cell_width_z_comoving == rhs.cell_width_z_comoving &&
          dt_time_code == rhs.dt_time_code &&
          hydro_boundary == rhs.hydro_boundary &&
          patch_signature == rhs.patch_signature &&
          gas_cell_identity_generation == rhs.gas_cell_identity_generation &&
          cell_index_generation == rhs.cell_index_generation &&
          row_mapping_signature == rhs.row_mapping_signature;
    }
  };

  [[nodiscard]] static std::uint64_t patchSignature(const core::PatchSoa& patches) {
    std::uint64_t signature = 1469598103934665603ULL;
    const auto mix = [&signature](std::uint64_t value) {
      signature ^= value;
      signature *= 1099511628211ULL;
    };
    mix(static_cast<std::uint64_t>(patches.size()));
    for (std::size_t patch = 0; patch < patches.size(); ++patch) {
      mix(patches.patch_id[patch]);
      mix(static_cast<std::uint64_t>(static_cast<std::int64_t>(patches.level[patch])));
      mix(patches.first_cell[patch]);
      mix(patches.cell_count[patch]);
      mix(patches.owning_rank[patch]);
    }
    return signature;
  }

  void rebuildGeometryIfNeeded(const core::SimulationState& state, double dt_time_code) {
    const std::size_t cell_count = state.cells.size();
    if (cell_count == 0) {
      m_cached_geometry_key.reset();
      m_geometry = {};
      return;
    }

    const internal::CartesianGasCellRowLayout layout = requireCartesianGasCellRowLayout(
        state, m_config, "reference hydro workflow");
    const hydro::HydroCartesianPatchSpec& spec = layout.spec;
    const HydroGeometryCacheKey key{
        .cell_count = cell_count,
        .nx = spec.nx,
        .ny = spec.ny,
        .nz = spec.nz,
        .origin_x_comoving = spec.origin_x_comoving,
        .origin_y_comoving = spec.origin_y_comoving,
        .origin_z_comoving = spec.origin_z_comoving,
        .cell_width_x_comoving = spec.cell_width_x_comoving,
        .cell_width_y_comoving = spec.cell_width_y_comoving,
        .cell_width_z_comoving = spec.cell_width_z_comoving,
        .dt_time_code = dt_time_code,
        .hydro_boundary = m_mode_policy.hydro_boundary,
        .patch_signature = patchSignature(state.patches),
        .gas_cell_identity_generation = state.gasCellIdentityGeneration(),
        .cell_index_generation = state.cellIndexGeneration(),
        .row_mapping_signature = layout.mapping_signature,
    };
    if (m_cached_geometry_key.has_value() && *m_cached_geometry_key == key) {
      return;
    }

    m_cached_geometry_key = key;
    m_dense_row_by_geometry_row = layout.dense_row_by_geometry_row;
    m_geometry_row_by_dense_row = layout.geometry_row_by_dense_row;
    m_geometry = hydro::makeCartesianPatchGeometry(spec);
    hydro::appendCartesianBoundaryGhostFaces(
        m_geometry,
        hydroBoundaryKindFromModePolicy(m_mode_policy.hydro_boundary));
    const double dx = std::max(1.0e-6, spec.cell_width_x_comoving);
    const double dy = std::max(1.0e-6, spec.cell_width_y_comoving);
    const double dz = std::max(1.0e-6, spec.cell_width_z_comoving);
    m_reconstruction = hydro::MusclHancockReconstruction(hydro::HydroReconstructionPolicy{
        .limiter = hydro::HydroSlopeLimiter::kMonotonizedCentral,
        .dt_over_dx_code = dt_time_code / dx,
        .dt_over_cell_width_code = {dt_time_code / dx, dt_time_code / dy, dt_time_code / dz},
        .rho_floor = k_density_floor,
        .pressure_floor = k_pressure_floor,
        .enable_muscl_hancock_predictor = true,
        .adiabatic_index = k_gamma_adiabatic,
    });
  }

  [[nodiscard]] hydro::HydroActiveSetView buildActiveFaceView(std::span<const std::uint32_t> active_cells) {
    m_active_cells.clear();
    m_active_cells.reserve(active_cells.size());
    for (const std::uint32_t dense_row : active_cells) {
      if (dense_row >= m_geometry_row_by_dense_row.size()) {
        throw std::out_of_range("active gas-cell row is outside Cartesian row map");
      }
      m_active_cells.push_back(m_geometry_row_by_dense_row[dense_row]);
    }
    std::unordered_set<std::size_t> active_cell_lookup(m_active_cells.begin(), m_active_cells.end());
    m_active_faces.clear();
    for (std::size_t face_index = 0; face_index < m_geometry.faces.size(); ++face_index) {
      const hydro::HydroFace& face = m_geometry.faces[face_index];
      if (m_remote_non_authority_faces.contains(face_index)) {
        continue;
      }
      if (face.ghost_cell_slot != hydro::k_invalid_ghost_cell_slot) {
        const hydro::HydroGhostCell& ghost = m_geometry.ghost_cells.at(face.ghost_cell_slot);
        if (ghost.boundary_kind == hydro::HydroBoundaryKind::kImportedMpi &&
            m_current_world_rank > static_cast<std::uint32_t>(ghost.origin_rank)) {
          continue;
        }
      }
      if (active_cell_lookup.contains(face.owner_cell) ||
          (face.neighbor_cell != hydro::k_invalid_cell_index && active_cell_lookup.contains(face.neighbor_cell))) {
        m_active_faces.push_back(face_index);
      }
    }
    return hydro::HydroActiveSetView{.active_cells = m_active_cells, .active_faces = m_active_faces};
  }

  [[nodiscard]] std::optional<std::pair<std::uint64_t, std::uint32_t>> patchIdentityForCellRow(
      const core::SimulationState& state,
      std::uint32_t cell_index) const {
    const auto* identity = state.gas_cell_identity.findByLocalRow(cell_index);
    if (identity == nullptr || identity->gas_cell_id == 0U) {
      throw std::runtime_error("hydro CFL diagnostics rejected missing gas-cell identity");
    }
    if (state.patches.size() == 0U) {
      return std::nullopt;
    }
    if (identity->owning_patch_id == 0U || cell_index >= m_geometry_row_by_dense_row.size()) {
      throw std::runtime_error("hydro CFL diagnostics rejected incomplete fixed-patch identity metadata");
    }
    const auto patch_it = std::find(
        state.patches.patch_id.begin(), state.patches.patch_id.end(), identity->owning_patch_id);
    if (patch_it == state.patches.patch_id.end()) {
      throw std::runtime_error("hydro CFL diagnostics rejected an identity patch absent from PatchSoa");
    }
    const std::size_t patch_index = static_cast<std::size_t>(std::distance(state.patches.patch_id.begin(), patch_it));
    if (state.cells.patch_index[cell_index] != patch_index) {
      throw std::runtime_error("hydro CFL diagnostics rejected stale dense patch-index metadata");
    }
    return std::pair<std::uint64_t, std::uint32_t>{
        identity->owning_patch_id,
        m_geometry_row_by_dense_row[cell_index]};
  }

  [[nodiscard]] core::HydroCflDiagnostics hydroCflDiagnosticsForCell(
      const core::SimulationState& state,
      std::uint32_t cell_index,
      double accepted_dt_time_code) const {
    const core::DirectionalCflTimeStepInput input{
        .cell_width_axis_code = {
            m_geometry.cell_width_x_comoving,
            m_geometry.cell_width_y_comoving,
            m_geometry.cell_width_z_comoving},
        .velocity_axis_code = {
            state.gas_cells.velocity_x_peculiar[cell_index],
            state.gas_cells.velocity_y_peculiar[cell_index],
            state.gas_cells.velocity_z_peculiar[cell_index]},
        .sound_speed_code = std::max(state.gas_cells.sound_speed_code[cell_index], 0.0),
    };
    const auto patch_identity = patchIdentityForCellRow(state, cell_index);
    return core::makeHydroCflDiagnostics(
        cell_index,
        input,
        0.4,
        accepted_dt_time_code,
        state.gas_cell_identity.findByLocalRow(cell_index) != nullptr
            ? state.gas_cell_identity.findByLocalRow(cell_index)->gas_cell_id
            : 0U,
        patch_identity.has_value() ? std::optional<std::uint64_t>(patch_identity->first) : std::nullopt,
        patch_identity.has_value() ? std::optional<std::uint32_t>(patch_identity->second) : std::nullopt);
  }

  void verifyAcceptedHydroCfl(
      core::StepContext& context,
      std::span<const std::uint32_t> active_cells) {
    if (context.integrator_state.dt_time_code <= 0.0) {
      throw std::invalid_argument("hydro CFL guard requires dt_time_code > 0");
    }
    const std::size_t cell_count = context.state.cells.size();
    if (cell_count == 0) {
      return;
    }

    std::vector<std::uint32_t> all_cells;
    std::span<const std::uint32_t> cells_to_check = active_cells;
    if (cells_to_check.empty()) {
      all_cells.resize(cell_count);
      std::iota(all_cells.begin(), all_cells.end(), 0U);
      cells_to_check = all_cells;
    }

    core::HydroCflDiagnostics worst;
    bool have_worst = false;
    for (const std::size_t cell : cells_to_check) {
      if (cell >= cell_count) {
        throw std::out_of_range("hydro CFL guard active cell index out of range");
      }
      const core::HydroCflDiagnostics diagnostics = hydroCflDiagnosticsForCell(
          context.state,
          static_cast<std::uint32_t>(cell),
          context.integrator_state.dt_time_code);
      if (!have_worst || diagnostics.safety_factor < worst.safety_factor) {
        worst = diagnostics;
        have_worst = true;
      }
    }
    if (!have_worst) {
      return;
    }
    m_last_hydro_cfl_diagnostics = worst;
    if (context.profiler_session != nullptr) {
      context.profiler_session->recordEvent(core::RuntimeEvent{
          .event_kind = "hydro.cfl_guard",
          .severity = core::RuntimeEventSeverity::kInfo,
          .subsystem = "hydro.godunov",
          .step_index = context.integrator_state.step_index,
          .simulation_time_code = context.integrator_state.current_time_code,
          .scale_factor = context.integrator_state.current_scale_factor,
          .message = "checked accepted hydro timestep against local directional CFL bound",
          .payload = {{"local_row", std::to_string(worst.local_row)},
                      {"gas_cell_id", worst.has_gas_cell_id ? std::to_string(worst.gas_cell_id) : "n/a"},
                      {"patch_id", worst.has_patch_id ? std::to_string(worst.patch_id) : "n/a"},
                      {"patch_row", worst.has_patch_row ? std::to_string(worst.patch_row) : "n/a"},
                      {"proposed_dt_time_code", std::to_string(worst.proposed_dt_time_code)},
                      {"accepted_dt_time_code", std::to_string(worst.accepted_dt_time_code)},
                      {"cfl_number", std::to_string(worst.cfl_number)},
                      {"safety_factor", std::to_string(worst.safety_factor)},
                      {"velocity_x_code", std::to_string(worst.velocity_axis_code[0])},
                      {"velocity_y_code", std::to_string(worst.velocity_axis_code[1])},
                      {"velocity_z_code", std::to_string(worst.velocity_axis_code[2])},
                      {"sound_speed_code", std::to_string(worst.sound_speed_code)}}});
    }
    core::assertHydroCflStable(worst);
  }

  const core::SimulationConfig& m_config;
  const core::ModePolicy& m_mode_policy;
  const GravityAccelerationProvider& m_gravity_callback;
  const parallel::MpiContext& m_mpi_context;
  std::uint32_t m_current_world_rank = 0;
  hydro::HydroCoreSolver m_solver;
  hydro::MusclHancockReconstruction m_reconstruction;
  hydro::HllcRiemannSolver m_riemann_solver;
  hydro::HydroConservedStateSoa m_conserved;
  hydro::HydroScratchBuffers m_scratch;
  hydro::HydroPrimitiveCacheSoa m_primitive_cache;
  hydro::HydroPatchGeometry m_geometry;
  // Dense CellSoa row <-> physical Cartesian solver-row maps.
  std::vector<std::uint32_t> m_dense_row_by_geometry_row;
  std::vector<std::uint32_t> m_geometry_row_by_dense_row;
  std::vector<double> m_ordered_cell_accel_x;
  std::vector<double> m_ordered_cell_accel_y;
  std::vector<double> m_ordered_cell_accel_z;
  hydro::HydroProfileEvent m_last_hydro_profile{};
  core::HydroCflDiagnostics m_last_hydro_cfl_diagnostics{};
  internal::SolverGhostRefreshReport m_last_ghost_refresh{};
  HydroRemoteGhostBoundaryReport m_last_remote_ghost_boundary_report{};
  parallel::GhostCacheLifecycle m_ghost_cache_lifecycle{};
  std::vector<std::size_t> m_active_cells;
  std::vector<std::size_t> m_active_faces;
  std::unordered_set<std::size_t> m_remote_non_authority_faces;
  std::optional<HydroGeometryCacheKey> m_cached_geometry_key;
};
}  // namespace

std::unique_ptr<HydroAmrRuntime> makeHydroAmrRuntime(
    const core::SimulationConfig& config,
    const core::ModePolicy& mode_policy,
    const GravityAccelerationProvider& gravity_acceleration,
    const RuntimeServices& services) {
  return std::make_unique<HydroAmrRuntimeImpl>(
      config, mode_policy, gravity_acceleration, services);
}

}  // namespace cosmosim::workflows
