#include "workflows/internal/particle_ghost_runtime.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
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

[[nodiscard]] std::vector<parallel::LocalGhostDescriptor> buildParticleGhostDescriptors(
    const core::SimulationState& state,
    int world_rank,
    const parallel::GhostLayerEpoch& epoch) {
  std::vector<parallel::LocalGhostDescriptor> descriptors;
  descriptors.reserve(state.particles.size());
  for (std::size_t particle_index = 0; particle_index < state.particles.size(); ++particle_index) {
    const int owner_rank = static_cast<int>(state.particle_sidecar.owning_rank[particle_index]);
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

[[nodiscard]] parallel::GhostExchangeBufferSoA buildParticleGhostPayloadState(
    const core::SimulationState& state,
    const parallel::GhostLayerEpoch& epoch) {
  parallel::GhostExchangeBufferSoA payload;
  const std::size_t particle_count = state.particles.size();
  payload.epoch = epoch;
  payload.entity_id.assign(state.particle_sidecar.particle_id.begin(), state.particle_sidecar.particle_id.end());
  payload.position_x_comoving.assign(state.particles.position_x_comoving.begin(), state.particles.position_x_comoving.end());
  payload.position_y_comoving.assign(state.particles.position_y_comoving.begin(), state.particles.position_y_comoving.end());
  payload.position_z_comoving.assign(state.particles.position_z_comoving.begin(), state.particles.position_z_comoving.end());
  payload.mass_code.assign(state.particles.mass_code.begin(), state.particles.mass_code.end());
  payload.velocity_x_code.assign(state.particles.velocity_x_peculiar.begin(), state.particles.velocity_x_peculiar.end());
  payload.velocity_y_code.assign(state.particles.velocity_y_peculiar.begin(), state.particles.velocity_y_peculiar.end());
  payload.velocity_z_code.assign(state.particles.velocity_z_peculiar.begin(), state.particles.velocity_z_peculiar.end());
  payload.density_code.assign(particle_count, 0.0);
  payload.pressure_code.assign(particle_count, 0.0);
  payload.internal_energy_code.assign(particle_count, 0.0);

  for (std::size_t cell_index = 0; cell_index < state.gas_cells.size(); ++cell_index) {
    const std::uint64_t parent_id = state.gas_cells.parent_particle_id[cell_index];
    const auto it = std::find(state.particle_sidecar.particle_id.begin(), state.particle_sidecar.particle_id.end(), parent_id);
    if (it == state.particle_sidecar.particle_id.end()) {
      continue;
    }
    const auto particle_index = static_cast<std::size_t>(std::distance(state.particle_sidecar.particle_id.begin(), it));
    payload.density_code[particle_index] = state.gas_cells.density_code[cell_index];
    payload.pressure_code[particle_index] = state.gas_cells.pressure_code[cell_index];
    payload.internal_energy_code[particle_index] = state.gas_cells.internal_energy_code[cell_index];
  }
  if (!payload.isConsistent() || !payload.hasGravityPayload()) {
    throw std::runtime_error("particle ghost payload construction produced inconsistent gravity lanes");
  }
  return payload;
}

void applyCommittedParticleGhostPayload(
    core::SimulationState& state,
    int world_rank,
    const std::vector<parallel::LocalGhostDescriptor>& descriptors,
    const parallel::GhostExchangeBufferSoA& payload) {
  if (!payload.isConsistent() || !payload.hasGravityPayload() || !payload.hasHydroPayload()) {
    throw std::invalid_argument("committed particle ghost payload must contain gravity and hydro lanes");
  }
  if (descriptors.size() != state.particles.size() || payload.size() < descriptors.size()) {
    throw std::invalid_argument("committed particle ghost payload shape does not match SimulationState particles");
  }

  std::unordered_map<std::uint64_t, std::size_t> particle_row_by_id;
  particle_row_by_id.reserve(state.particles.size());
  for (std::size_t particle_index = 0; particle_index < state.particles.size(); ++particle_index) {
    particle_row_by_id.emplace(state.particle_sidecar.particle_id[particle_index], particle_index);
  }

  for (std::size_t particle_index = 0; particle_index < descriptors.size(); ++particle_index) {
    const auto& descriptor = descriptors[particle_index];
    if (descriptor.residency != parallel::LocalIndexResidency::kGhost) {
      continue;
    }
    if (descriptor.owning_rank == world_rank) {
      throw std::logic_error("local ghost descriptor cannot be owned by local rank during commit");
    }
    if (payload.entity_id[particle_index] != descriptor.particle_id) {
      throw std::runtime_error("committed ghost payload entity_id drifted from descriptor particle_id");
    }
    // These rows are imported, non-authoritative ghost copies. Updating them is
    // allowed only because the authoritative owner remains descriptor.owning_rank.
    state.particles.position_x_comoving[particle_index] = payload.position_x_comoving[particle_index];
    state.particles.position_y_comoving[particle_index] = payload.position_y_comoving[particle_index];
    state.particles.position_z_comoving[particle_index] = payload.position_z_comoving[particle_index];
    state.particles.mass_code[particle_index] = payload.mass_code[particle_index];
    state.particles.velocity_x_peculiar[particle_index] = payload.velocity_x_code[particle_index];
    state.particles.velocity_y_peculiar[particle_index] = payload.velocity_y_code[particle_index];
    state.particles.velocity_z_peculiar[particle_index] = payload.velocity_z_code[particle_index];
  }

  for (std::size_t cell_index = 0; cell_index < state.gas_cells.size(); ++cell_index) {
    const std::uint64_t parent_id = state.gas_cells.parent_particle_id[cell_index];
    const auto particle_it = particle_row_by_id.find(parent_id);
    if (particle_it == particle_row_by_id.end()) {
      continue;
    }
    const std::size_t particle_index = particle_it->second;
    if (descriptors[particle_index].residency != parallel::LocalIndexResidency::kGhost) {
      continue;
    }
    state.cells.center_x_comoving[cell_index] = payload.position_x_comoving[particle_index];
    state.cells.center_y_comoving[cell_index] = payload.position_y_comoving[particle_index];
    state.cells.center_z_comoving[cell_index] = payload.position_z_comoving[particle_index];
    state.cells.mass_code[cell_index] = payload.mass_code[particle_index];
    state.gas_cells.density_code[cell_index] = payload.density_code[particle_index];
    state.gas_cells.pressure_code[cell_index] = payload.pressure_code[particle_index];
    state.gas_cells.internal_energy_code[cell_index] = payload.internal_energy_code[particle_index];
  }
}

}  // namespace

[[nodiscard]] SolverGhostRefreshReport refreshParticleGhostsForSolver(
    core::StepContext& context,
    const parallel::MpiContext& mpi_context,
    std::string_view subsystem_name,
    parallel::GhostCacheLifecycle* lifecycle) {
  const int world_rank = mpi_context.worldRank();
  const parallel::GhostLayerEpoch epoch = makeRuntimeGhostLayerEpoch(context);
  if (lifecycle != nullptr) {
    parallel::invalidateGhostCache(*lifecycle, epoch);
  }
  std::vector<parallel::LocalGhostDescriptor> descriptors = buildParticleGhostDescriptors(
      context.state, world_rank, epoch);
  parallel::GhostExchangeBufferSoA ghost_storage = buildParticleGhostPayloadState(context.state, epoch);

  const auto exchange = parallel::executeBlockingGhostRefreshExchangeFromDescriptors(
      mpi_context, descriptors, ghost_storage, epoch);
  const auto commit_report = parallel::commitBlockingGhostRefreshResult(
      ghost_storage, descriptors, exchange.plan, exchange.result, epoch);
  if (lifecycle != nullptr) {
    parallel::markGhostCacheCommitted(*lifecycle, epoch);
    parallel::requireValidGhostCache(*lifecycle, epoch, std::string(subsystem_name));
  }
  applyCommittedParticleGhostPayload(context.state, world_rank, descriptors, ghost_storage);

  if (context.profiler_session != nullptr) {
    context.profiler_session->recordEvent(core::RuntimeEvent{
        .event_kind = "parallel.blocking_ghost_refresh",
        .severity = core::RuntimeEventSeverity::kInfo,
        .subsystem = std::string(subsystem_name),
        .step_index = context.integrator_state.step_index,
        .simulation_time_code = context.integrator_state.current_time_code,
        .scale_factor = context.integrator_state.current_scale_factor,
        .message = "blocking correctness-first particle ghost refresh completed before solver access",
        .payload = {
            {"neighbor_count", std::to_string(exchange.plan.neighbor_ranks.size())},
            {"sent_bytes", std::to_string(exchange.result.sent_bytes)},
            {"received_bytes", std::to_string(exchange.result.received_bytes)},
            {"committed_ghost_slots", std::to_string(commit_report.updated_ghost_slots)},
            {"ghost_sync_epoch", std::to_string(epoch.ghost_sync_epoch)},
            {"ghost_cache_refresh_count", lifecycle != nullptr ? std::to_string(lifecycle->refresh_count) : "0"},
            {"ghost_cache_invalidation_count", lifecycle != nullptr ? std::to_string(lifecycle->invalidation_count) : "0"},
        },
    });
  }
  return SolverGhostRefreshReport{
      .sent_bytes = exchange.result.sent_bytes,
      .received_bytes = exchange.result.received_bytes,
      .committed_slots = commit_report.updated_ghost_slots,
  };
}

}  // namespace cosmosim::workflows::internal
