#include <cassert>
#include <cstdint>
#include <limits>
#include <string>
#include <stdexcept>
#include <vector>

#include "cosmosim/core/memory_accounting.hpp"
#include "cosmosim/core/memory_governor.hpp"
#include "cosmosim/core/retained_capacity_transaction.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/core/time_scheduler.hpp"
#include "../../src/workflows/internal/retained_population_growth.hpp"

namespace {
namespace core = cosmosim::core;

void testReplacementAndFailureReconciliation() {
  std::vector<std::uint64_t> first{7U, 11U};
  std::vector<std::uint64_t> second{13U};
  const auto owned = [&]() {
    return core::checkedMemoryBytesAdd(
        core::ownedCapacityBytesForContainer(first),
        core::ownedCapacityBytesForContainer(second), "test retained total");
  };
  const std::uint64_t original = owned();
  core::MemoryGovernor governor(core::MemoryGovernorPolicy{
      .hard_limit_bytes = original + 4U * sizeof(std::uint64_t) - 1U});
  governor.setBaselineOwnedBytes(original);
  core::RetainedCapacityTransaction rejected(owned);
  rejected.add(first, 4U);
  bool denied = false;
  try {
    rejected.execute(&governor, core::MemoryClass::kCanonicalPersistent, "test.growth");
  } catch (const core::MemoryAdmissionError&) { denied = true; }
  assert(denied);
  assert(first.size() == 2U && first[0] == 7U && first[1] == 11U);
  assert(governor.snapshot().baseline_owned_bytes == original);
  assert(governor.snapshot().committed_bytes == 0U);

  core::MemoryGovernor exact(core::MemoryGovernorPolicy{
      .hard_limit_bytes = original + 4U * sizeof(std::uint64_t)});
  exact.setBaselineOwnedBytes(original);
  core::RetainedCapacityTransaction growth(owned);
  growth.add(first, 4U);
  assert(growth.replacementBytes() == 4U * sizeof(std::uint64_t));
  growth.execute(&exact, core::MemoryClass::kCanonicalPersistent, "test.exact");
  assert(first.capacity() == 4U && first.size() == 2U);
  assert(first[0] == 7U && first[1] == 11U);
  assert(exact.snapshot().baseline_owned_bytes == owned());
  assert(exact.snapshot().committed_bytes == 0U);
  assert(exact.snapshot().peak_accounted_bytes <= exact.policy().hard_limit_bytes);

  // Multiple retained lanes transfer their physical backing exactly once.
  core::MemoryGovernor failure(core::MemoryGovernorPolicy{
      .hard_limit_bytes = owned() + 8U * sizeof(std::uint64_t)});
  failure.setBaselineOwnedBytes(owned());
  core::RetainedCapacityTransaction partial(owned);
  partial.add(first, 6U);
  partial.add(second, 2U);
  partial.execute(&failure, core::MemoryClass::kCanonicalPersistent, "test.partial");
  assert(first.capacity() == 6U && second.capacity() == 2U);
  assert(first[0] == 7U && second[0] == 13U);
  assert(failure.snapshot().baseline_owned_bytes == owned());

  core::RetainedCapacityTransaction reuse(owned);
  reuse.add(first, 6U);
  reuse.add(second, 2U);
  assert(reuse.replacementBytes() == 0U);
  reuse.execute(&failure, core::MemoryClass::kCanonicalPersistent, "test.reuse");
  assert(failure.snapshot().committed_bytes == 0U);

  bool overflow = false;
  try {
    core::RetainedCapacityTransaction huge(owned);
    huge.add(first, std::numeric_limits<std::size_t>::max());
  } catch (const std::exception&) { overflow = true; }
  assert(overflow);
}

struct ThrowingValue {
  std::uint64_t value = 0U;
  static inline bool throw_on_copy = false;
  ThrowingValue() = default;
  explicit ThrowingValue(std::uint64_t v) : value(v) {}
  ThrowingValue(const ThrowingValue& other) : value(other.value) {
    if (throw_on_copy) throw std::runtime_error("injected replacement copy failure");
  }
  ThrowingValue& operator=(const ThrowingValue&) = default;
  ThrowingValue(ThrowingValue&&) noexcept = default;
  ThrowingValue& operator=(ThrowingValue&&) noexcept = default;
};

void testPartialGrowthException() {
  std::vector<std::uint64_t> first{7U, 11U};
  std::vector<ThrowingValue> second;
  second.emplace_back(13U);
  const auto owned = [&]() {
    return core::checkedMemoryBytesAdd(
        core::ownedCapacityBytesForContainer(first),
        core::ownedCapacityBytesForContainer(second), "injected growth baseline");
  };
  const auto initial = owned();
  core::RetainedCapacityTransaction plan(owned);
  plan.add(first, 4U);
  plan.add(second, 3U);
  core::MemoryGovernor governor(core::MemoryGovernorPolicy{
      .hard_limit_bytes = initial + plan.replacementBytes()});
  governor.setBaselineOwnedBytes(initial);
  ThrowingValue::throw_on_copy = true;
  bool failed = false;
  try {
    plan.execute(&governor, core::MemoryClass::kCanonicalPersistent, "test.partial_failure");
  } catch (const std::runtime_error& error) {
    failed = std::string(error.what()) == "injected replacement copy failure";
  }
  ThrowingValue::throw_on_copy = false;
  assert(failed);
  assert(first.size() == 2U && first.capacity() == 4U);
  assert(first[0] == 7U && first[1] == 11U);
  assert(second.size() == 1U && second[0].value == 13U);
  assert(governor.snapshot().baseline_owned_bytes == owned());
  assert(governor.snapshot().committed_bytes == 0U);
  assert(governor.snapshot().peak_accounted_bytes <= governor.policy().hard_limit_bytes);
  bool reused_plan = false;
  try {
    plan.execute(&governor, core::MemoryClass::kCanonicalPersistent, "test.invalid_retry");
  } catch (const std::logic_error&) { reused_plan = true; }
  assert(reused_plan);
  // A later transaction can retry the failed lane without charging the first
  // lane's already-retained capacity a second time.
  core::RetainedCapacityTransaction retry(owned);
  retry.add(first, 4U);
  retry.add(second, 3U);
  assert(retry.replacementBytes() == 3U * sizeof(ThrowingValue));
  retry.execute(&governor, core::MemoryClass::kCanonicalPersistent, "test.retry");
  assert(second.size() == 1U && second.capacity() == 3U);
  assert(governor.snapshot().baseline_owned_bytes == owned());
  assert(governor.snapshot().committed_bytes == 0U);
}

void testPopulationAndSchedulerCapacity() {
  core::SimulationState state;
  state.resizeParticles(1U);
  state.particle_sidecar.particle_id[0] = 1U;
  state.particle_sidecar.species_tag[0] =
      static_cast<std::uint32_t>(core::ParticleSpecies::kDarkMatter);
  state.species.count_by_species[core::particleSpeciesIndex(core::ParticleSpecies::kDarkMatter)] = 1U;
  state.rebuildSpeciesIndex();
  core::HierarchicalTimeBinScheduler scheduler;
  scheduler.reset(1U, 0U, 0U);
  const auto owned = [&]() {
    return core::checkedMemoryBytesAdd(
        core::memoryReportBaselineOwnedBytes(core::collectSimulationMemoryReport(state)),
        scheduler.ownedCapacityBytes(), "test population baseline");
  };
  const std::uint64_t initial = owned();
  core::RetainedCapacityTransaction plan(owned);
  cosmosim::workflows::internal::planParticlePopulationGrowth(
      plan, state, 3U, 2U, core::ParticleSpecies::kStar);
  scheduler.planAppendCapacity(plan, 2U, 0U);
  const std::uint64_t replacement = plan.replacementBytes();
  assert(replacement > 0U);
  core::MemoryGovernor governor(core::MemoryGovernorPolicy{
      .hard_limit_bytes = core::checkedMemoryBytesAdd(initial, replacement,
                                                      "test exact population ceiling")});
  governor.setBaselineOwnedBytes(initial);
  plan.execute(&governor, core::MemoryClass::kCanonicalPersistent, "test.population");
  assert(state.particles.size() == 1U);
  assert(state.star_particles.size() == 0U);
  assert(scheduler.elementCount() == 1U);
  assert(governor.snapshot().baseline_owned_bytes == owned());
  assert(governor.snapshot().committed_bytes == 0U);
  state.resizeParticles(3U);
  state.star_particles.resize(2U);
  scheduler.appendElements(2U, 0U, 1U);
  assert(scheduler.elementCount() == 3U);
  assert(scheduler.hotMetadata().next_activation_tick[1] == 1U);
  assert(scheduler.hotMetadata().next_activation_tick[2] == 1U);
  assert(governor.snapshot().baseline_owned_bytes == owned());
  assert(governor.snapshot().peak_accounted_bytes <= governor.policy().hard_limit_bytes);
  core::RetainedCapacityTransaction repeat(owned);
  cosmosim::workflows::internal::planParticlePopulationGrowth(
      repeat, state, 3U, 2U, core::ParticleSpecies::kStar);
  assert(repeat.replacementBytes() == 0U);

  // DMO and star populations must not acquire BH, gas, or tracer optional lanes.
  assert(state.black_holes.size() == 0U);
  assert(state.gas_cells.size() == 0U);
  assert(state.tracers.size() == 0U);
}

}  // namespace
int main() {
  testReplacementAndFailureReconciliation();
  testPartialGrowthException();
  testPopulationAndSchedulerCapacity();
  return 0;
}
