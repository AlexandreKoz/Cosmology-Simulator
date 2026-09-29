#include "workflows/internal/particle_ghost_runtime.hpp"

#include <cstddef>
#include <cstdint>
#include <exception>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "cosmosim/core/profiling.hpp"
#include "cosmosim/core/simulation_state.hpp"

namespace cosmosim::workflows::internal {
namespace {

[[nodiscard]] parallel::GhostLayerEpoch makeRuntimeGhostLayerEpoch(const core::StepContext& context) {
  return parallel::GhostLayerEpoch{
      .decomposition_epoch = context.state.particleIndexGeneration(),
      .ghost_sync_epoch = context.integrator_state.step_index * core::integrationStageCount() +
          core::integrationStageIndex(context.stage) + 1U,
      .particle_index_generation = context.state.particleIndexGeneration(),
  };
}

[[nodiscard]] bool detectLocalParticleGhostDemand(
    const core::SimulationState& state,
    int world_rank,
    int world_size) {
  const std::size_t particle_count = state.particles.size();
  if (world_rank < 0 || world_size <= 0 || world_rank >= world_size) {
    throw std::invalid_argument("particle ghost demand scan received an invalid MPI rank topology");
  }
  if (state.particles.position_y_comoving.size() != particle_count ||
      state.particles.position_z_comoving.size() != particle_count ||
      state.particles.velocity_x_peculiar.size() != particle_count ||
      state.particles.velocity_y_peculiar.size() != particle_count ||
      state.particles.velocity_z_peculiar.size() != particle_count ||
      state.particles.mass_code.size() != particle_count) {
    throw std::invalid_argument(
        "particle ghost demand scan requires component-consistent gravity/kinematic particle lanes");
  }
  if ((!state.hasHomogeneousDmoMetadata() &&
       state.particle_sidecar.owning_rank.size() != particle_count) ||
      state.particle_sidecar.particle_id.size() != particle_count) {
    throw std::invalid_argument(
        "particle ghost demand scan requires ownership and particle-ID sidecars aligned with particle state");
  }

  bool local_has_ghost_demand = false;
  for (std::size_t particle_index = 0; particle_index < particle_count; ++particle_index) {
    const std::uint32_t owner_rank = state.particleOwningRank(particle_index);
    if (owner_rank >= static_cast<std::uint32_t>(world_size)) {
      throw std::invalid_argument("particle ghost demand scan found owning_rank outside the MPI world");
    }
    local_has_ghost_demand = local_has_ghost_demand ||
        owner_rank != static_cast<std::uint32_t>(world_rank);
  }
  return local_has_ghost_demand;
}

[[nodiscard]] std::vector<parallel::LocalGhostDescriptor> buildParticleGhostDescriptors(
    const core::SimulationState& state,
    int world_rank,
    const parallel::GhostLayerEpoch& epoch) {
  std::vector<parallel::LocalGhostDescriptor> descriptors;
  descriptors.reserve(state.particles.size());
  for (std::size_t particle_index = 0; particle_index < state.particles.size(); ++particle_index) {
    const int owner_rank = static_cast<int>(state.particleOwningRank(particle_index));
    descriptors.push_back(parallel::LocalGhostDescriptor{
        .residency = (owner_rank == world_rank) ? parallel::LocalIndexResidency::kOwned
                                                : parallel::LocalIndexResidency::kGhost,
        .owning_rank = owner_rank,
        .particle_id = state.particle_sidecar.particle_id[particle_index],
        .epoch = epoch,
    });
  }
  return descriptors;
}

[[nodiscard]] parallel::ReadOnlyGhostExchangeView buildParticleGhostPayloadView(
    const core::SimulationState& state,
    const parallel::GhostLayerEpoch& epoch) {
  parallel::ReadOnlyGhostExchangeView view{
      .epoch = epoch,
      .entity_id = state.particle_sidecar.particle_id,
      .position_x_comoving = state.particles.position_x_comoving,
      .position_y_comoving = state.particles.position_y_comoving,
      .position_z_comoving = state.particles.position_z_comoving,
      .mass_code = state.particles.mass_code,
      .density_code = {},
      .velocity_x_code = state.particles.velocity_x_peculiar,
      .velocity_y_code = state.particles.velocity_y_peculiar,
      .velocity_z_code = state.particles.velocity_z_peculiar,
      .pressure_code = {},
      .internal_energy_code = {},
  };
  parallel::requireFreshGhostExchangeView(view, epoch);
  if (view.position_x_comoving.size() != view.size() ||
      view.position_y_comoving.size() != view.size() ||
      view.position_z_comoving.size() != view.size() ||
      view.mass_code.size() != view.size() ||
      view.velocity_x_code.size() != view.size() ||
      view.velocity_y_code.size() != view.size() ||
      view.velocity_z_code.size() != view.size()) {
    throw std::runtime_error("particle ghost payload view is missing required gravity/kinematic lanes");
  }
  return view;
}

[[nodiscard]] parallel::GhostRefreshCommitReport commitParticleGhostPayloadToState(
    core::SimulationState& state,
    int world_rank,
    std::span<const parallel::LocalGhostDescriptor> descriptors,
    const parallel::GhostExchangePlan& plan,
    const parallel::BlockingGhostExchangeResult& result,
    const parallel::GhostLayerEpoch& epoch) {
  parallel::validateBlockingGhostExchangeContracts(plan, descriptors, world_rank, epoch);
  if (!result.received_ghosts.isConsistent() ||
      !result.received_ghosts.epoch.matches(epoch)) {
    throw std::invalid_argument("received particle ghost payload is inconsistent or stale");
  }
  if (result.received_ghosts.hasHydroPayload()) {
    throw std::invalid_argument(
        "generic particle ghost payload must not contain hydro lanes; use gas_cell_id keyed hydro ghost exchange");
  }

  const std::size_t received_count = result.received_ghosts.size();
  if (result.received_ghosts.position_x_comoving.size() != received_count ||
      result.received_ghosts.position_y_comoving.size() != received_count ||
      result.received_ghosts.position_z_comoving.size() != received_count ||
      result.received_ghosts.mass_code.size() != received_count ||
      result.received_ghosts.velocity_x_code.size() != received_count ||
      result.received_ghosts.velocity_y_code.size() != received_count ||
      result.received_ghosts.velocity_z_code.size() != received_count) {
    throw std::invalid_argument("received particle ghost payload is missing required gravity/kinematic lanes");
  }

  std::size_t expected_count = 0U;
  for (const auto& indices : plan.recv_local_indices_by_neighbor) {
    expected_count += indices.size();
  }
  if (received_count != expected_count) {
    throw std::invalid_argument("received particle ghost payload count does not match exchange plan");
  }

  parallel::GhostRefreshCommitReport report;
  std::size_t result_row = 0U;
  for (std::size_t slot = 0; slot < plan.recv_local_indices_by_neighbor.size(); ++slot) {
    const int peer_rank = plan.neighbor_ranks[slot];
    for (const std::uint32_t local_index : plan.recv_local_indices_by_neighbor[slot]) {
      if (local_index >= descriptors.size() || local_index >= state.particles.size()) {
        throw std::out_of_range("particle ghost commit target row is outside local state");
      }
      const auto& descriptor = descriptors[local_index];
      if (descriptor.residency != parallel::LocalIndexResidency::kGhost ||
          descriptor.owning_rank != peer_rank || descriptor.owning_rank == world_rank) {
        throw std::invalid_argument("particle ghost commit target does not match remote ownership plan");
      }
      if (!descriptor.epoch.matches(epoch)) {
        throw std::invalid_argument("particle ghost commit target descriptor is stale");
      }
      if (state.particle_sidecar.particle_id[local_index] != descriptor.particle_id ||
          result.received_ghosts.entity_id[result_row] != descriptor.particle_id) {
        throw std::runtime_error("particle ghost commit particle-ID authority mismatch");
      }

      state.particles.position_x_comoving[local_index] = result.received_ghosts.position_x_comoving[result_row];
      state.particles.position_y_comoving[local_index] = result.received_ghosts.position_y_comoving[result_row];
      state.particles.position_z_comoving[local_index] = result.received_ghosts.position_z_comoving[result_row];
      state.particles.mass_code[local_index] = result.received_ghosts.mass_code[result_row];
      state.particles.velocity_x_peculiar[local_index] = result.received_ghosts.velocity_x_code[result_row];
      state.particles.velocity_y_peculiar[local_index] = result.received_ghosts.velocity_y_code[result_row];
      state.particles.velocity_z_peculiar[local_index] = result.received_ghosts.velocity_z_code[result_row];
      ++result_row;
      ++report.updated_ghost_slots;
    }
  }
  report.committed_payload_bytes = static_cast<std::uint64_t>(report.updated_ghost_slots) *
      static_cast<std::uint64_t>(parallel::ghostRefreshPayloadRecordBytes());
  return report;
}

void recordParticleGhostRefreshEvent(
    const core::StepContext& context,
    const parallel::GhostLayerEpoch& epoch,
    std::string_view subsystem_name,
    const parallel::GhostCacheLifecycle* lifecycle,
    std::size_t neighbor_count,
    std::uint64_t sent_bytes,
    std::uint64_t received_bytes,
    std::size_t committed_slots,
    bool globally_empty) {
  if (context.profiler_session == nullptr) {
    return;
  }
  context.profiler_session->recordEvent(core::RuntimeEvent{
      .event_kind = "parallel.blocking_ghost_refresh",
      .severity = core::RuntimeEventSeverity::kInfo,
      .subsystem = std::string(subsystem_name),
      .step_index = context.integrator_state.step_index,
      .simulation_time_code = context.integrator_state.current_time_code,
      .scale_factor = context.integrator_state.current_scale_factor,
      .message = globally_empty
          ? "globally empty particle ghost refresh committed without population-scale staging"
          : "blocking particle ghost refresh completed before solver access",
      .payload = {
          {"neighbor_count", std::to_string(neighbor_count)},
          {"sent_bytes", std::to_string(sent_bytes)},
          {"received_bytes", std::to_string(received_bytes)},
          {"committed_ghost_slots", std::to_string(committed_slots)},
          {"global_ghost_demand", globally_empty ? "0" : "1"},
          {"ghost_sync_epoch", std::to_string(epoch.ghost_sync_epoch)},
          {"ghost_cache_refresh_count", lifecycle != nullptr ? std::to_string(lifecycle->refresh_count) : "0"},
          {"ghost_cache_invalidation_count", lifecycle != nullptr ? std::to_string(lifecycle->invalidation_count) : "0"},
      },
  });
}

}  // namespace

[[nodiscard]] SolverGhostRefreshReport refreshParticleGhostsForSolver(
    core::StepContext& context,
    const parallel::MpiContext& mpi_context,
    std::string_view subsystem_name,
    parallel::GhostCacheLifecycle* lifecycle) {
  const int world_rank = mpi_context.worldRank();
  const int world_size = mpi_context.worldSize();
  const parallel::GhostLayerEpoch epoch = makeRuntimeGhostLayerEpoch(context);
  if (lifecycle != nullptr) {
    parallel::invalidateGhostCache(*lifecycle, epoch);
  }

  bool local_has_ghost_demand = false;
  std::exception_ptr local_demand_failure;
  try {
    local_has_ghost_demand = detectLocalParticleGhostDemand(
        context.state, world_rank, world_size);
  } catch (...) {
    local_demand_failure = std::current_exception();
  }
  mpi_context.rethrowCollectivePreparationFailure(
      local_demand_failure,
      "particle ghost empty-demand discovery");
  const bool global_has_ghost_demand =
      mpi_context.allreduceMaxUint64(local_has_ghost_demand ? 1ULL : 0ULL) != 0ULL;

  if (!global_has_ghost_demand) {
    if (lifecycle != nullptr) {
      parallel::markGhostCacheCommitted(*lifecycle, epoch);
      parallel::requireValidGhostCache(*lifecycle, epoch, std::string(subsystem_name));
    }
    recordParticleGhostRefreshEvent(
        context, epoch, subsystem_name, lifecycle, 0U, 0U, 0U, 0U, true);
    return {};
  }

  std::vector<parallel::LocalGhostDescriptor> descriptors;
  parallel::ReadOnlyGhostExchangeView authoritative_view;
  std::exception_ptr local_materialization_failure;
  try {
    descriptors = buildParticleGhostDescriptors(context.state, world_rank, epoch);
    authoritative_view = buildParticleGhostPayloadView(context.state, epoch);
  } catch (...) {
    local_materialization_failure = std::current_exception();
  }
  mpi_context.rethrowCollectivePreparationFailure(
      local_materialization_failure,
      "particle ghost descriptor/view preparation");

  const auto exchange = parallel::executeBlockingGhostRefreshExchangeFromDescriptors(
      mpi_context, descriptors, authoritative_view, epoch);
  const auto commit_report = commitParticleGhostPayloadToState(
      context.state, world_rank, descriptors, exchange.plan, exchange.result, epoch);
  if (lifecycle != nullptr) {
    parallel::markGhostCacheCommitted(*lifecycle, epoch);
    parallel::requireValidGhostCache(*lifecycle, epoch, std::string(subsystem_name));
  }

  recordParticleGhostRefreshEvent(
      context,
      epoch,
      subsystem_name,
      lifecycle,
      exchange.plan.neighbor_ranks.size(),
      exchange.result.sent_bytes,
      exchange.result.received_bytes,
      commit_report.updated_ghost_slots,
      false);
  return SolverGhostRefreshReport{
      .sent_bytes = exchange.result.sent_bytes,
      .received_bytes = exchange.result.received_bytes,
      .committed_slots = commit_report.updated_ghost_slots,
  };
}

}  // namespace cosmosim::workflows::internal
