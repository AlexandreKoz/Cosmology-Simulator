#include "workflows/internal/gravity_force_endpoint.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace cosmosim::workflows {

std::uint64_t TransactionalForceEndpoint::retainedCapacityBytes() const noexcept {
  std::uint64_t bytes = 0U;
  for (const auto* lane : {&total_x_comoving, &total_y_comoving, &total_z_comoving,
                           &tree_x_comoving, &tree_y_comoving, &tree_z_comoving,
                           &pm_x_comoving, &pm_y_comoving, &pm_z_comoving,
                           &nearest_pair_r2_comoving}) {
    bytes += static_cast<std::uint64_t>(lane->capacity()) * sizeof(double);
  }
  return bytes;
}

std::size_t TransactionalForceEndpoint::minimumLaneCapacity() const noexcept {
  std::size_t minimum = total_x_comoving.capacity();
  for (const auto* lane : {&total_y_comoving, &total_z_comoving,
                           &tree_x_comoving, &tree_y_comoving, &tree_z_comoving,
                           &pm_x_comoving, &pm_y_comoving, &pm_z_comoving,
                           &nearest_pair_r2_comoving}) {
    minimum = std::min(minimum, lane->capacity());
  }
  return minimum;
}

void TransactionalForceEndpoint::reserveRows(std::size_t rows) {
  for (auto* lane : {&total_x_comoving, &total_y_comoving, &total_z_comoving,
                     &tree_x_comoving, &tree_y_comoving, &tree_z_comoving,
                     &pm_x_comoving, &pm_y_comoving, &pm_z_comoving,
                     &nearest_pair_r2_comoving}) {
    lane->reserve(rows);
  }
}

void TransactionalForceEndpoint::release() noexcept {
  for (auto* lane : {&total_x_comoving, &total_y_comoving, &total_z_comoving,
                     &tree_x_comoving, &tree_y_comoving, &tree_z_comoving,
                     &pm_x_comoving, &pm_y_comoving, &pm_z_comoving,
                     &nearest_pair_r2_comoving}) {
    std::vector<double>().swap(*lane);
  }
  particle_count = 0U;
  particle_index_generation = 0U;
  layout_certificate = 0U;
  pm_field_version = 0U;
  valid = false;
}

}  // namespace cosmosim::workflows

namespace cosmosim::workflows::internal {
namespace {

[[nodiscard]] std::array<std::span<double>, TransactionalForceEndpoint::k_lane_count> lanesOf(
    const ForceEndpointLanes& lanes) {
  return {lanes.total_x, lanes.total_y, lanes.total_z, lanes.tree_x, lanes.tree_y,
          lanes.tree_z, lanes.pm_x, lanes.pm_y, lanes.pm_z, lanes.nearest_pair_r2};
}

[[nodiscard]] std::array<std::span<const double>, TransactionalForceEndpoint::k_lane_count> lanesOf(
    const ConstForceEndpointLanes& lanes) {
  return {lanes.total_x, lanes.total_y, lanes.total_z, lanes.tree_x, lanes.tree_y,
          lanes.tree_z, lanes.pm_x, lanes.pm_y, lanes.pm_z, lanes.nearest_pair_r2};
}

[[nodiscard]] std::array<const std::vector<double>*, TransactionalForceEndpoint::k_lane_count> lanesOf(
    const TransactionalForceEndpoint& endpoint) {
  return {&endpoint.total_x_comoving, &endpoint.total_y_comoving, &endpoint.total_z_comoving,
          &endpoint.tree_x_comoving, &endpoint.tree_y_comoving, &endpoint.tree_z_comoving,
          &endpoint.pm_x_comoving, &endpoint.pm_y_comoving, &endpoint.pm_z_comoving,
          &endpoint.nearest_pair_r2_comoving};
}

[[nodiscard]] std::array<std::vector<double>*, TransactionalForceEndpoint::k_lane_count> lanesOf(
    TransactionalForceEndpoint& endpoint) {
  return {&endpoint.total_x_comoving, &endpoint.total_y_comoving, &endpoint.total_z_comoving,
          &endpoint.tree_x_comoving, &endpoint.tree_y_comoving, &endpoint.tree_z_comoving,
          &endpoint.pm_x_comoving, &endpoint.pm_y_comoving, &endpoint.pm_z_comoving,
          &endpoint.nearest_pair_r2_comoving};
}

[[noreturn]] void refuse(std::string_view detail) {
  throw std::logic_error("transactional force endpoint restore refused: " + std::string(detail));
}

}  // namespace

void captureForceEndpointRowAligned(
    const ConstForceEndpointLanes& source, const RowAlignedEndpointStamp& stamp,
    std::uint64_t pm_field_version, TransactionalForceEndpoint& endpoint) {
  endpoint.valid = false;
  const auto source_lanes = lanesOf(source);
  for (const auto lane : source_lanes) {
    if (static_cast<std::uint64_t>(lane.size()) != stamp.particle_count) {
      throw std::logic_error(
          "transactional force endpoint capture: a source lane does not cover every particle row");
    }
  }
  if (pm_field_version == 0U) {
    throw std::logic_error("transactional force endpoint capture: no committed PM field version");
  }
  const auto destination = lanesOf(endpoint);
  for (std::size_t lane = 0U; lane < destination.size(); ++lane) {
    destination[lane]->assign(source_lanes[lane].begin(), source_lanes[lane].end());
  }
  endpoint.particle_count = stamp.particle_count;
  endpoint.particle_index_generation = stamp.particle_index_generation;
  endpoint.layout_certificate = stamp.layout_certificate;
  endpoint.pm_field_version = pm_field_version;
  endpoint.valid = true;
}

void validateRowAlignedForceEndpoint(
    const TransactionalForceEndpoint& endpoint, const RowAlignedEndpointStamp& expected,
    const ConstForceEndpointLanes& destination) {
  if (!endpoint.valid) refuse("no captured endpoint");
  if (endpoint.particle_count != expected.particle_count) {
    refuse("particle count changed since the cycle opening (" + std::to_string(endpoint.particle_count) +
           " captured, " + std::to_string(expected.particle_count) + " now)");
  }
  if (endpoint.particle_index_generation != expected.particle_index_generation) {
    refuse("particle row layout changed since the cycle opening (particle-index generation " +
           std::to_string(endpoint.particle_index_generation) + " captured, " +
           std::to_string(expected.particle_index_generation) + " now)");
  }
  if (endpoint.layout_certificate != expected.layout_certificate) {
    refuse("particle identity/layout certificate differs from the cycle opening (row identity not proven)");
  }
  for (const auto* lane : lanesOf(endpoint)) {
    if (static_cast<std::uint64_t>(lane->size()) != expected.particle_count) {
      refuse("captured endpoint lane extent does not match the particle count");
    }
  }
  for (const auto lane : lanesOf(destination)) {
    if (static_cast<std::uint64_t>(lane.size()) != expected.particle_count) {
      refuse("runtime endpoint lane extent does not match the particle count");
    }
  }
}

void restoreForceEndpointRowAligned(
    const TransactionalForceEndpoint& endpoint, const RowAlignedEndpointStamp& expected,
    const ForceEndpointLanes& destination) {
  validateRowAlignedForceEndpoint(endpoint, expected, asConst(destination));
  const auto source = lanesOf(endpoint);
  const auto target = lanesOf(destination);
  for (std::size_t lane = 0U; lane < source.size(); ++lane) {
    std::copy(source[lane]->begin(), source[lane]->end(), target[lane].begin());
  }
}

void remapForceEndpointById(
    const io::GravityForceCachePersistentState& cache,
    std::span<const std::uint64_t> state_particle_ids,
    const ForceEndpointLanes& destination, bool include_hierarchical_endpoint) {
  const std::size_t n = state_particle_ids.size();
  const auto requireTotals = [&](const std::vector<double>& lane) {
    if (lane.size() != n) {
      throw std::runtime_error(
          "restart gravity force cache particle extent does not match authoritative SimulationState");
    }
    for (const double value : lane) {
      if (!std::isfinite(value)) {
        throw std::runtime_error("restart gravity force cache contains non-finite values");
      }
    }
  };
  requireTotals(cache.particle_accel_x_comoving);
  requireTotals(cache.particle_accel_y_comoving);
  requireTotals(cache.particle_accel_z_comoving);
  const std::array<const std::vector<double>*, 7> endpoint_lanes{
      &cache.hierarchical_tree_accel_x_comoving, &cache.hierarchical_tree_accel_y_comoving,
      &cache.hierarchical_tree_accel_z_comoving, &cache.hierarchical_pm_accel_x_comoving,
      &cache.hierarchical_pm_accel_y_comoving, &cache.hierarchical_pm_accel_z_comoving,
      &cache.hierarchical_nearest_pair_r2_comoving};
  if (include_hierarchical_endpoint) {
    for (const auto* lane : endpoint_lanes) {
      if (lane->size() != n) {
        throw std::runtime_error(
            "restart hierarchical endpoint extent does not match authoritative SimulationState");
      }
    }
  }
  const auto target = lanesOf(destination);
  const std::size_t written_lanes = include_hierarchical_endpoint ? target.size() : 3U;
  for (std::size_t lane = 0U; lane < written_lanes; ++lane) {
    if (target[lane].size() != n) {
      throw std::logic_error("restart gravity force cache destination lane does not cover every particle row");
    }
  }
  if (cache.particle_id.size() != n) {
    throw std::runtime_error(
        "restart gravity force cache particle identity extent does not match acceleration lanes");
  }
  std::unordered_map<std::uint64_t, std::size_t> index_by_id;
  index_by_id.reserve(n);
  for (std::size_t index = 0; index < n; ++index) {
    if (cache.particle_id[index] == 0U || !index_by_id.emplace(cache.particle_id[index], index).second) {
      throw std::runtime_error(
          "restart gravity force cache particle identity lane contains a zero or duplicate stable ID");
    }
  }
  // Resolve every row before writing so a missing ID leaves the destination untouched.
  std::vector<std::size_t> cached_row_of(n);
  for (std::size_t row = 0; row < n; ++row) {
    const auto it = index_by_id.find(state_particle_ids[row]);
    if (it == index_by_id.end()) {
      throw std::runtime_error("restart gravity force cache is missing a persistent particle ID");
    }
    cached_row_of[row] = it->second;
  }
  const std::array<const std::vector<double>*, 3> totals{
      &cache.particle_accel_x_comoving, &cache.particle_accel_y_comoving, &cache.particle_accel_z_comoving};
  for (std::size_t row = 0; row < n; ++row) {
    const std::size_t cached_row = cached_row_of[row];
    for (std::size_t lane = 0U; lane < 3U; ++lane) target[lane][row] = (*totals[lane])[cached_row];
    if (include_hierarchical_endpoint) {
      for (std::size_t lane = 0U; lane < endpoint_lanes.size(); ++lane) {
        target[3U + lane][row] = (*endpoint_lanes[lane])[cached_row];
      }
    }
  }
}

}  // namespace cosmosim::workflows::internal
