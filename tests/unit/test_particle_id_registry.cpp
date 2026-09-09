#include <array>
#include <cassert>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include "cosmosim/core/memory_governor.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/parallel/distributed_memory.hpp"
#include "cosmosim/physics/star_formation.hpp"
#include "../../src/workflows/internal/particle_id_registry.hpp"

namespace {
void testShardedRegistryAndLeaseLifetime() {
  using namespace cosmosim;
  parallel::MpiContext mpi_context(false, 1, 0);
  core::MemoryGovernor governor(core::MemoryGovernorPolicy{.hard_limit_bytes = 1U << 20U});
  auto registry = workflows::internal::makeDistributedParticleIdRegistry(mpi_context, &governor);
  core::SimulationState state;
  state.resizeParticles(2U);
  constexpr std::uint64_t k_birth_key = 0x123456789ULL;
  const std::uint64_t collision_id = physics::starFormationParticleIdFromBirthKey(k_birth_key, 0U);
  state.particle_sidecar.particle_id[0] = collision_id;
  state.particle_sidecar.particle_id[1] = 987654321U;
  const std::array keys{k_birth_key};
  const auto first = registry->precommit(state, keys);
  assert(first.size() == 1U);
  assert(first[0] == physics::starFormationParticleIdFromBirthKey(k_birth_key, 1U));
  const auto live = governor.snapshot().committed_bytes;
  assert(live >= 3U * sizeof(std::uint64_t));
  registry->finishPrecommit();
  assert(governor.snapshot().committed_bytes < live);
  const auto second = registry->precommit(state, std::array<std::uint64_t,1>{k_birth_key + 1U});
  assert(second.size() == 1U);
  registry->finishPrecommit();
  assert(registry->precommit(state, {}).empty());
  registry->finishPrecommit();
  bool duplicate_rejected = false;
  try {
    (void)registry->precommit(state, std::array{k_birth_key + 2U, k_birth_key + 2U});
  } catch (const std::runtime_error&) { duplicate_rejected = true; }
  assert(duplicate_rejected);
  registry->finishPrecommit();
}

void testInvalidInitialPopulationAndTightAdmission() {
  using namespace cosmosim;
  parallel::MpiContext mpi_context(false, 1, 0);
  core::SimulationState state;
  state.resizeParticles(2U);
  state.particle_sidecar.particle_id[0] = 17U;
  state.particle_sidecar.particle_id[1] = 17U;
  core::MemoryGovernor governor(core::MemoryGovernorPolicy{.hard_limit_bytes = 1U << 20U});
  auto registry = workflows::internal::makeDistributedParticleIdRegistry(mpi_context, &governor);
  bool rejected = false;
  try { (void)registry->precommit(state, std::array<std::uint64_t,1>{19U}); }
  catch (const std::runtime_error&) { rejected = true; }
  assert(rejected);
  state.particle_sidecar.particle_id[1] = 18U;
  const auto recovered = registry->precommit(state, std::array<std::uint64_t,1>{19U});
  assert(recovered.size() == 1U);
  registry->finishPrecommit();

  core::MemoryGovernor tiny(core::MemoryGovernorPolicy{.hard_limit_bytes = 1U});
  auto limited = workflows::internal::makeDistributedParticleIdRegistry(mpi_context, &tiny);
  bool memory_rejected = false;
  try { (void)limited->precommit(state, std::array<std::uint64_t,1>{20U}); }
  catch (const core::MemoryAdmissionError&) { memory_rejected = true; }
  assert(memory_rejected);
  assert(tiny.snapshot().reserved_bytes == 0U);
  assert(tiny.snapshot().committed_bytes == 0U);
}
}

int main() {
  testShardedRegistryAndLeaseLifetime();
  testInvalidInitialPopulationAndTightAdmission();
}
