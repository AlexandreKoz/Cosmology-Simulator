// Transactional (row-aligned) versus ordinary (ID-remapped) force-endpoint
// restoration (2026-10-27; AUTHORED, NOT RUN). Exercises the production
// helpers shared by GravityRuntimeImpl:
//   * internal::remapForceEndpointById  -- ordinary restart import path;
//   * internal::capture/restoreForceEndpointRowAligned -- rollback path.
// Checks values (bitwise), ordering, repeated restores, every refusal
// invariant with "nothing written" semantics, and that a restore performs
// no heap allocation at all (global operator new is counted in this test).
#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <new>
#include <numeric>
#include <random>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

#include "../support/test_require.hpp"
#include "cosmosim/io/restart_checkpoint.hpp"
#include "cosmosim/workflows/gravity_runtime.hpp"
#include "workflows/internal/gravity_force_endpoint.hpp"

namespace {
std::atomic<std::uint64_t> g_allocations{0U};
}  // namespace

// Counting replacements for the ordinary (unaligned) global allocation
// functions of this test executable.
void* operator new(std::size_t bytes) {
  g_allocations.fetch_add(1U, std::memory_order_relaxed);
  if (void* pointer = std::malloc(bytes == 0U ? 1U : bytes)) return pointer;
  throw std::bad_alloc();
}
void* operator new[](std::size_t bytes) {
  g_allocations.fetch_add(1U, std::memory_order_relaxed);
  if (void* pointer = std::malloc(bytes == 0U ? 1U : bytes)) return pointer;
  throw std::bad_alloc();
}
void operator delete(void* pointer) noexcept { std::free(pointer); }
void operator delete[](void* pointer) noexcept { std::free(pointer); }
void operator delete(void* pointer, std::size_t) noexcept { std::free(pointer); }
void operator delete[](void* pointer, std::size_t) noexcept { std::free(pointer); }

namespace {
using namespace cosmosim;
using workflows::TransactionalForceEndpoint;
namespace internal = workflows::internal;

constexpr std::size_t k_rows = 257U;
constexpr std::size_t k_lanes = TransactionalForceEndpoint::k_lane_count;

// Ten dense runtime-owned lanes (total xyz, tree xyz, PM xyz, nearest r2).
struct Lanes {
  std::vector<std::vector<double>> lane = std::vector<std::vector<double>>(k_lanes, std::vector<double>(k_rows, 0.0));
  internal::ForceEndpointLanes view() {
    return {lane[0], lane[1], lane[2], lane[3], lane[4], lane[5], lane[6], lane[7], lane[8], lane[9]};
  }
  internal::ConstForceEndpointLanes constView() const {
    return {lane[0], lane[1], lane[2], lane[3], lane[4], lane[5], lane[6], lane[7], lane[8], lane[9]};
  }
  void fill(double value) { for (auto& values : lane) std::fill(values.begin(), values.end(), value); }
  bool operator==(const Lanes& other) const {
    if (lane.size() != other.lane.size()) return false;
    for (std::size_t l = 0; l < lane.size(); ++l) {
      // Some negative-path fixtures deliberately truncate one lane; compare
      // extents before element access to avoid undefined behavior in the test.
      if (lane[l].size() != other.lane[l].size()) return false;
      for (std::size_t row = 0; row < lane[l].size(); ++row) {
        // Bitwise comparison (also distinguishes signed zero / NaN payloads).
        if (std::memcmp(&lane[l][row], &other.lane[l][row], sizeof(double)) != 0) return false;
      }
    }
    return true;
  }
};

// Value of lane `l` for stable particle ID `id` (finite; nearest r2 may be +inf).
double reference(std::size_t l, std::uint64_t id) {
  if (l == 9U && id % 7U == 0U) return std::numeric_limits<double>::infinity();
  return std::ldexp(static_cast<double>(id * 131U + l * 17U + 1U), -7) * (l % 2U == 0U ? 1.0 : -1.0);
}

// Checkpoint cache whose rows are a permutation of the state rows.
io::GravityForceCachePersistentState permutedCache(const std::vector<std::uint64_t>& state_ids, std::uint32_t seed) {
  std::vector<std::size_t> order(state_ids.size());
  std::iota(order.begin(), order.end(), 0U);
  std::mt19937 rng(seed);
  std::shuffle(order.begin(), order.end(), rng);
  io::GravityForceCachePersistentState cache;
  cache.valid = true;
  cache.hierarchical_endpoint_valid = true;
  cache.hierarchical_pm_field_version = 9U;
  std::vector<std::vector<double>*> lanes{
      &cache.particle_accel_x_comoving, &cache.particle_accel_y_comoving, &cache.particle_accel_z_comoving,
      &cache.hierarchical_tree_accel_x_comoving, &cache.hierarchical_tree_accel_y_comoving,
      &cache.hierarchical_tree_accel_z_comoving, &cache.hierarchical_pm_accel_x_comoving,
      &cache.hierarchical_pm_accel_y_comoving, &cache.hierarchical_pm_accel_z_comoving,
      &cache.hierarchical_nearest_pair_r2_comoving};
  for (const auto row : order) {
    const std::uint64_t id = state_ids[row];
    cache.particle_id.push_back(id);
    for (std::size_t l = 0; l < k_lanes; ++l) lanes[l]->push_back(reference(l, id));
  }
  return cache;
}

template <class F>
bool refusedWithoutWrite(F&& action, const Lanes& destination, const std::string& fragment) {
  const Lanes before = destination;
  bool refused = false;
  try {
    action();
  } catch (const std::exception& error) {
    refused = std::string(error.what()).find(fragment) != std::string::npos;
  }
  return refused && destination == before;
}

void testRowAlignedAgreesWithIdRemap() {
  std::vector<std::uint64_t> ids(k_rows);
  for (std::size_t row = 0; row < k_rows; ++row) ids[row] = 1000U + 3U * (k_rows - row);  // non-monotone
  const auto cache = permutedCache(ids, 17U);

  // Ordinary restart import: rows matched by ID even though the cache order differs.
  Lanes imported;
  internal::remapForceEndpointById(cache, ids, imported.view(), true);
  for (std::size_t row = 0; row < k_rows; ++row) {
    for (std::size_t l = 0; l < k_lanes; ++l) {
      const double expected = reference(l, ids[row]);
      CHUI_REQUIRE(std::memcmp(&imported.lane[l][row], &expected, sizeof(double)) == 0);
    }
  }
  // Without the hierarchical endpoint only the totals are written.
  Lanes totals_only;
  totals_only.fill(-5.0);
  internal::remapForceEndpointById(cache, ids, totals_only.view(), false);
  for (std::size_t row = 0; row < k_rows; ++row) {
    CHUI_REQUIRE(totals_only.lane[0][row] == reference(0U, ids[row]));
    CHUI_REQUIRE(totals_only.lane[5][row] == -5.0);
  }

  // Transactional capture from the imported (live) lanes, then restore into
  // lanes the rejected attempt overwrote: bitwise identical to the import.
  const internal::RowAlignedEndpointStamp stamp{.particle_count = k_rows, .particle_index_generation = 42U,
                                                .layout_certificate = 0xC0FFEEULL};
  TransactionalForceEndpoint endpoint;
  endpoint.reserveRows(k_rows);
  CHUI_REQUIRE(endpoint.minimumLaneCapacity() >= k_rows);
  CHUI_REQUIRE(endpoint.retainedCapacityBytes() >= k_lanes * k_rows * sizeof(double));
  const auto capacity_before_capture = endpoint.retainedCapacityBytes();
  internal::captureForceEndpointRowAligned(imported.constView(), stamp, 9U, endpoint);
  CHUI_REQUIRE(endpoint.valid && endpoint.particle_count == k_rows && endpoint.pm_field_version == 9U);
  CHUI_REQUIRE(endpoint.particle_index_generation == 42U && endpoint.layout_certificate == 0xC0FFEEULL);
  CHUI_REQUIRE(endpoint.retainedCapacityBytes() == capacity_before_capture);  // capacity reused

  Lanes live = imported;
  for (int rejected_cycle = 0; rejected_cycle < 3; ++rejected_cycle) {
    live.fill(1.0e300 * (rejected_cycle + 1));  // the rejected attempt's forces
    const auto allocations_before = g_allocations.load();
    internal::restoreForceEndpointRowAligned(endpoint, stamp, live.view());
    CHUI_REQUIRE(g_allocations.load() == allocations_before);  // no ID map, no allocation
    CHUI_REQUIRE(live == imported);                             // values and row order
  }
  // The endpoint itself is untouched by restores (repeatable rollback).
  CHUI_REQUIRE(endpoint.valid && endpoint.total_x_comoving == imported.lane[0]);

  endpoint.release();
  CHUI_REQUIRE(!endpoint.valid && endpoint.retainedCapacityBytes() == 0U);
}

void testRowAlignedRefusals() {
  std::vector<std::uint64_t> ids(k_rows);
  std::iota(ids.begin(), ids.end(), 1U);
  Lanes source;
  internal::remapForceEndpointById(permutedCache(ids, 5U), ids, source.view(), true);
  const internal::RowAlignedEndpointStamp stamp{.particle_count = k_rows, .particle_index_generation = 7U,
                                                .layout_certificate = 99U};
  TransactionalForceEndpoint endpoint;
  internal::captureForceEndpointRowAligned(source.constView(), stamp, 3U, endpoint);
  Lanes destination;
  destination.fill(2.5);

  auto changed = stamp;
  changed.particle_count = k_rows - 1U;
  CHUI_REQUIRE(refusedWithoutWrite([&] { internal::restoreForceEndpointRowAligned(endpoint, changed, destination.view()); },
                                   destination, "particle count changed"));
  changed = stamp;
  changed.particle_index_generation = 8U;  // reorder/migration/compaction happened
  CHUI_REQUIRE(refusedWithoutWrite([&] { internal::restoreForceEndpointRowAligned(endpoint, changed, destination.view()); },
                                   destination, "row layout changed"));
  changed = stamp;
  changed.layout_certificate = 100U;  // same count, different identity (row identity not proven)
  CHUI_REQUIRE(refusedWithoutWrite([&] { internal::restoreForceEndpointRowAligned(endpoint, changed, destination.view()); },
                                   destination, "certificate"));
  // A runtime lane of the wrong extent is refused before any lane is written.
  Lanes short_lane = destination;
  short_lane.lane[7].resize(k_rows - 1U);
  CHUI_REQUIRE(refusedWithoutWrite([&] { internal::restoreForceEndpointRowAligned(endpoint, stamp, short_lane.view()); },
                                   short_lane, "runtime endpoint lane extent"));
  // A truncated captured lane is refused.
  auto truncated = endpoint;
  truncated.pm_y_comoving.pop_back();
  CHUI_REQUIRE(refusedWithoutWrite([&] { internal::restoreForceEndpointRowAligned(truncated, stamp, destination.view()); },
                                   destination, "captured endpoint lane extent"));
  // No captured endpoint.
  TransactionalForceEndpoint empty;
  CHUI_REQUIRE(refusedWithoutWrite([&] { internal::restoreForceEndpointRowAligned(empty, stamp, destination.view()); },
                                   destination, "no captured endpoint"));
  // Validation alone never writes and accepts the matching stamp.
  internal::validateRowAlignedForceEndpoint(endpoint, stamp, destination.constView());
  CHUI_REQUIRE(destination.lane[0][0] == 2.5);

  // Capture refusals leave the endpoint invalid.
  TransactionalForceEndpoint refused_capture;
  Lanes ragged = source;
  ragged.lane[2].pop_back();
  bool capture_refused = false;
  try {
    internal::captureForceEndpointRowAligned(ragged.constView(), stamp, 3U, refused_capture);
  } catch (const std::logic_error&) {
    capture_refused = true;
  }
  CHUI_REQUIRE(capture_refused && !refused_capture.valid);
  capture_refused = false;
  try {
    internal::captureForceEndpointRowAligned(source.constView(), stamp, 0U, refused_capture);  // no PM field
  } catch (const std::logic_error&) {
    capture_refused = true;
  }
  CHUI_REQUIRE(capture_refused && !refused_capture.valid);
}

void testIdRemapRefusals() {
  std::vector<std::uint64_t> ids(k_rows);
  std::iota(ids.begin(), ids.end(), 1U);
  Lanes destination;
  destination.fill(4.0);
  auto missing = permutedCache(ids, 3U);
  missing.particle_id[11] = 999999U;  // a state ID absent from the checkpoint
  CHUI_REQUIRE(refusedWithoutWrite([&] { internal::remapForceEndpointById(missing, ids, destination.view(), true); },
                                   destination, "missing a persistent particle ID"));
  auto duplicate = permutedCache(ids, 3U);
  duplicate.particle_id[12] = duplicate.particle_id[13];
  CHUI_REQUIRE(refusedWithoutWrite([&] { internal::remapForceEndpointById(duplicate, ids, destination.view(), true); },
                                   destination, "zero or duplicate"));
  auto zero = permutedCache(ids, 3U);
  zero.particle_id[0] = 0U;
  CHUI_REQUIRE(refusedWithoutWrite([&] { internal::remapForceEndpointById(zero, ids, destination.view(), true); },
                                   destination, "zero or duplicate"));
  auto non_finite = permutedCache(ids, 3U);
  non_finite.particle_accel_y_comoving[4] = std::numeric_limits<double>::quiet_NaN();
  CHUI_REQUIRE(refusedWithoutWrite([&] { internal::remapForceEndpointById(non_finite, ids, destination.view(), true); },
                                   destination, "non-finite"));
  auto short_endpoint = permutedCache(ids, 3U);
  short_endpoint.hierarchical_pm_accel_z_comoving.pop_back();
  CHUI_REQUIRE(refusedWithoutWrite([&] { internal::remapForceEndpointById(short_endpoint, ids, destination.view(), true); },
                                   destination, "hierarchical endpoint extent"));
}

}  // namespace

int main() {
  try {
    testRowAlignedAgreesWithIdRemap();
    testRowAlignedRefusals();
    testIdRemapRefusals();
  } catch (const std::exception& error) {
    std::cerr << "gravity force endpoint: " << error.what() << '\n';
    return 1;
  }
  return 0;
}
