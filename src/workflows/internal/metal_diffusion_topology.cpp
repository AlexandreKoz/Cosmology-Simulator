#include "workflows/internal/metal_diffusion_topology.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

namespace cosmosim::workflows::internal {
namespace {

struct TopologyFace {
  physics::MetalDiffusionFace face;
  std::uint8_t axis = 0U;
};

struct PatchBoundaryFace {
  std::uint32_t cell = 0U;
  std::uint8_t axis = 0U;
  std::int8_t side = 0;
  double plane = 0.0;
  double transverse_u_min = 0.0;
  double transverse_u_max = 0.0;
  double transverse_v_min = 0.0;
  double transverse_v_max = 0.0;
  double normal_width = 0.0;
};

struct NeighborAccumulator {
  double velocity_area_sum[3]{};
  double distance_area_sum = 0.0;
  double area_sum = 0.0;
};

[[nodiscard]] bool patchHasOwnedLeafCells(
    const core::SimulationState& state,
    std::span<const std::uint8_t> owned_leaf_mask,
    std::uint32_t patch_index,
    std::uint32_t world_rank) {
  if (patch_index >= state.patches.size() ||
      state.patches.owning_rank[patch_index] != world_rank) {
    return false;
  }
  const std::uint32_t first = state.patches.first_cell[patch_index];
  const std::uint32_t count = state.patches.cell_count[patch_index];
  if (count == 0U || first >= owned_leaf_mask.size() ||
      static_cast<std::size_t>(first) + static_cast<std::size_t>(count) >
          owned_leaf_mask.size()) {
    return false;
  }
  return std::any_of(
      owned_leaf_mask.begin() + static_cast<std::ptrdiff_t>(first),
      owned_leaf_mask.begin() + static_cast<std::ptrdiff_t>(first + count),
      [](std::uint8_t value) { return value != 0U; });
}

}  // namespace

MetalDiffusionTopologyResult buildMetalDiffusionTopology(
    const core::SimulationState& state,
    std::span<const std::uint8_t> owned_leaf_mask,
    std::uint32_t world_rank,
    double scale_factor,
    core::BoundaryCondition hydro_boundary,
    std::span<physics::MetalDiffusionCell> diffusion_cells) {
  if (!std::isfinite(scale_factor) || !(scale_factor > 0.0)) {
    throw std::invalid_argument(
        "metal diffusion topology requires finite positive scale factor");
  }
  if (owned_leaf_mask.size() != state.cells.size() ||
      diffusion_cells.size() != state.cells.size()) {
    throw std::invalid_argument(
        "metal diffusion topology cell/mask extent mismatch");
  }

  MetalDiffusionTopologyResult result;
  std::vector<TopologyFace> topology_faces;
  std::vector<PatchBoundaryFace> patch_boundary_faces;
  double topology_scale = 1.0;
  std::array<double, 3> domain_min{
      std::numeric_limits<double>::infinity(),
      std::numeric_limits<double>::infinity(),
      std::numeric_limits<double>::infinity()};
  std::array<double, 3> domain_max{
      -std::numeric_limits<double>::infinity(),
      -std::numeric_limits<double>::infinity(),
      -std::numeric_limits<double>::infinity()};

  for (std::uint32_t patch_index = 0;
       patch_index < state.patches.size(); ++patch_index) {
    if (state.patches.owning_rank[patch_index] != world_rank) {
      continue;
    }
    const bool has_owned_leaf_cells = patchHasOwnedLeafCells(
        state, owned_leaf_mask, patch_index, world_rank);
    const std::uint32_t first = state.patches.first_cell[patch_index];
    const std::uint32_t nx = state.patches.cell_dim_x[patch_index];
    const std::uint32_t ny = state.patches.cell_dim_y[patch_index];
    const std::uint32_t nz = state.patches.cell_dim_z[patch_index];
    if (nx == 0U || ny == 0U || nz == 0U) {
      throw std::runtime_error(
          "owned leaf diffusion patch has zero Cartesian dimension");
    }
    const std::uint64_t expected_count =
        static_cast<std::uint64_t>(nx) * static_cast<std::uint64_t>(ny) *
        static_cast<std::uint64_t>(nz);
    if (expected_count != state.patches.cell_count[patch_index]) {
      throw std::runtime_error(
          "owned leaf diffusion patch cell_count disagrees with Cartesian dimensions");
    }
    if (static_cast<std::uint64_t>(first) + expected_count > state.cells.size()) {
      throw std::runtime_error(
          "owned leaf diffusion patch cell range exceeds dense gas storage");
    }

    const double dx = state.patches.extent_x_comoving[patch_index] /
        static_cast<double>(nx) * scale_factor;
    const double dy = state.patches.extent_y_comoving[patch_index] /
        static_cast<double>(ny) * scale_factor;
    const double dz = state.patches.extent_z_comoving[patch_index] /
        static_cast<double>(nz) * scale_factor;
    const double ox = state.patches.origin_x_comoving[patch_index] * scale_factor;
    const double oy = state.patches.origin_y_comoving[patch_index] * scale_factor;
    const double oz = state.patches.origin_z_comoving[patch_index] * scale_factor;
    if (!(dx > 0.0) || !(dy > 0.0) || !(dz > 0.0) ||
        !std::isfinite(dx) || !std::isfinite(dy) || !std::isfinite(dz) ||
        !std::isfinite(ox) || !std::isfinite(oy) || !std::isfinite(oz)) {
      throw std::runtime_error(
          "owned leaf diffusion patch has non-finite or non-positive geometry");
    }

    domain_min[0] = std::min(domain_min[0], ox);
    domain_min[1] = std::min(domain_min[1], oy);
    domain_min[2] = std::min(domain_min[2], oz);
    domain_max[0] = std::max(
        domain_max[0], ox + dx * static_cast<double>(nx));
    domain_max[1] = std::max(
        domain_max[1], oy + dy * static_cast<double>(ny));
    domain_max[2] = std::max(
        domain_max[2], oz + dz * static_cast<double>(nz));
    topology_scale = std::max({
        topology_scale, std::abs(ox), std::abs(oy), std::abs(oz),
        std::abs(dx * static_cast<double>(nx)),
        std::abs(dy * static_cast<double>(ny)),
        std::abs(dz * static_cast<double>(nz))});

    // Covered/coarse patches still define the physical domain envelope used
    // for periodic wrap detection, but only authoritative leaf cells may
    // contribute diffusion faces.
    if (!has_owned_leaf_cells) {
      continue;
    }

    const auto row = [first, nx, ny](
                         std::uint32_t i, std::uint32_t j,
                         std::uint32_t k) -> std::uint32_t {
      return first + i + nx * (j + ny * k);
    };
    const auto addTopologyFace = [&](std::uint32_t left,
                                     std::uint32_t right,
                                     std::uint8_t axis,
                                     double area,
                                     double distance,
                                     physics::MetalDiffusionBoundaryKind kind) {
      if (left >= owned_leaf_mask.size() || right >= owned_leaf_mask.size() ||
          owned_leaf_mask[left] == 0U || owned_leaf_mask[right] == 0U) {
        return;
      }
      const physics::MetalDiffusionFace face{
          .left_cell = left,
          .right_cell = right,
          .area_code = area,
          .center_distance_code = distance,
          .boundary_kind = kind};
      result.faces.push_back(face);
      topology_faces.push_back(TopologyFace{.face = face, .axis = axis});
    };
    const auto addBoundaryFace = [&](std::uint32_t cell,
                                     std::uint8_t axis,
                                     std::int8_t side,
                                     double plane,
                                     double u0,
                                     double u1,
                                     double v0,
                                     double v1,
                                     double normal_width) {
      if (cell >= owned_leaf_mask.size() || owned_leaf_mask[cell] == 0U) {
        return;
      }
      patch_boundary_faces.push_back(PatchBoundaryFace{
          .cell = cell,
          .axis = axis,
          .side = side,
          .plane = plane,
          .transverse_u_min = u0,
          .transverse_u_max = u1,
          .transverse_v_min = v0,
          .transverse_v_max = v1,
          .normal_width = normal_width});
    };

    for (std::uint32_t k = 0; k < nz; ++k) {
      for (std::uint32_t j = 0; j < ny; ++j) {
        for (std::uint32_t i = 0; i < nx; ++i) {
          const std::uint32_t left = row(i, j, k);
          if (i + 1U < nx) {
            addTopologyFace(
                left, row(i + 1U, j, k), 0U, dy * dz, dx,
                physics::MetalDiffusionBoundaryKind::kInternal);
            ++result.interior_patch_face_count;
          }
          if (j + 1U < ny) {
            addTopologyFace(
                left, row(i, j + 1U, k), 1U, dx * dz, dy,
                physics::MetalDiffusionBoundaryKind::kInternal);
            ++result.interior_patch_face_count;
          }
          if (k + 1U < nz) {
            addTopologyFace(
                left, row(i, j, k + 1U), 2U, dx * dy, dz,
                physics::MetalDiffusionBoundaryKind::kInternal);
            ++result.interior_patch_face_count;
          }

          const double x0 = ox + static_cast<double>(i) * dx;
          const double x1 = x0 + dx;
          const double y0 = oy + static_cast<double>(j) * dy;
          const double y1 = y0 + dy;
          const double z0 = oz + static_cast<double>(k) * dz;
          const double z1 = z0 + dz;
          if (i == 0U) {
            addBoundaryFace(left, 0U, -1, x0, y0, y1, z0, z1, dx);
          }
          if (i + 1U == nx) {
            addBoundaryFace(left, 0U, +1, x1, y0, y1, z0, z1, dx);
          }
          if (j == 0U) {
            addBoundaryFace(left, 1U, -1, y0, x0, x1, z0, z1, dy);
          }
          if (j + 1U == ny) {
            addBoundaryFace(left, 1U, +1, y1, x0, x1, z0, z1, dy);
          }
          if (k == 0U) {
            addBoundaryFace(left, 2U, -1, z0, x0, x1, y0, y1, dz);
          }
          if (k + 1U == nz) {
            addBoundaryFace(left, 2U, +1, z1, x0, x1, y0, y1, dz);
          }
        }
      }
    }
  }

  if (result.faces.empty() && patch_boundary_faces.empty()) {
    return result;
  }

  const double face_tolerance = 1.0e-10 * topology_scale;
  const double quantization = std::max(face_tolerance, 1.0e-14);
  std::map<std::pair<std::uint8_t, std::int64_t>, std::vector<std::size_t>>
      lower_faces_by_plane;
  for (std::size_t index = 0; index < patch_boundary_faces.size(); ++index) {
    const PatchBoundaryFace& candidate = patch_boundary_faces[index];
    if (candidate.side >= 0) {
      continue;
    }
    const auto bin = static_cast<std::int64_t>(
        std::llround(candidate.plane / quantization));
    lower_faces_by_plane[{candidate.axis, bin}].push_back(index);
  }

  std::set<std::tuple<std::uint32_t, std::uint32_t, std::uint8_t>>
      unique_cross_patch_faces;
  const auto addOverlappingInterface = [&](const PatchBoundaryFace& upper,
                                           const PatchBoundaryFace& lower,
                                           physics::MetalDiffusionBoundaryKind kind,
                                           std::uint64_t& counter) {
    if (upper.cell == lower.cell || upper.axis != lower.axis) {
      return;
    }
    const double overlap_u =
        std::min(upper.transverse_u_max, lower.transverse_u_max) -
        std::max(upper.transverse_u_min, lower.transverse_u_min);
    const double overlap_v =
        std::min(upper.transverse_v_max, lower.transverse_v_max) -
        std::max(upper.transverse_v_min, lower.transverse_v_min);
    if (!(overlap_u > face_tolerance) || !(overlap_v > face_tolerance)) {
      return;
    }
    const std::uint32_t lo = std::min(upper.cell, lower.cell);
    const std::uint32_t hi = std::max(upper.cell, lower.cell);
    if (!unique_cross_patch_faces.emplace(lo, hi, upper.axis).second) {
      return;
    }
    const physics::MetalDiffusionFace face{
        .left_cell = upper.cell,
        .right_cell = lower.cell,
        .area_code = overlap_u * overlap_v,
        .center_distance_code = 0.5 *
            (upper.normal_width + lower.normal_width),
        .boundary_kind = kind};
    result.faces.push_back(face);
    topology_faces.push_back(TopologyFace{.face = face, .axis = upper.axis});
    ++counter;
  };

  for (const PatchBoundaryFace& upper : patch_boundary_faces) {
    if (upper.side <= 0) {
      continue;
    }
    const auto upper_bin = static_cast<std::int64_t>(
        std::llround(upper.plane / quantization));
    for (std::int64_t delta_bin = -1; delta_bin <= 1; ++delta_bin) {
      const auto found = lower_faces_by_plane.find(
          {upper.axis, upper_bin + delta_bin});
      if (found == lower_faces_by_plane.end()) {
        continue;
      }
      for (const std::size_t lower_index : found->second) {
        const PatchBoundaryFace& lower = patch_boundary_faces[lower_index];
        if (std::abs(upper.plane - lower.plane) > face_tolerance) {
          continue;
        }
        addOverlappingInterface(
            upper, lower, physics::MetalDiffusionBoundaryKind::kInternal,
            result.cross_patch_face_count);
      }
    }
  }

  if (hydro_boundary == core::BoundaryCondition::kPeriodic) {
    for (std::uint8_t axis = 0U; axis < 3U; ++axis) {
      if (!std::isfinite(domain_min[axis]) ||
          !std::isfinite(domain_max[axis]) ||
          !(domain_max[axis] > domain_min[axis])) {
        throw std::runtime_error(
            "periodic metal diffusion requires finite positive domain extent");
      }
      for (const PatchBoundaryFace& upper : patch_boundary_faces) {
        if (upper.axis != axis || upper.side <= 0 ||
            std::abs(upper.plane - domain_max[axis]) > face_tolerance) {
          continue;
        }
        for (const PatchBoundaryFace& lower : patch_boundary_faces) {
          if (lower.axis != axis || lower.side >= 0 ||
              std::abs(lower.plane - domain_min[axis]) > face_tolerance) {
            continue;
          }
          addOverlappingInterface(
              upper, lower, physics::MetalDiffusionBoundaryKind::kPeriodic,
              result.periodic_wrap_face_count);
        }
      }
    }
  }

  std::vector<std::array<std::array<NeighborAccumulator, 2>, 3>>
      neighbor_accumulators(state.cells.size());
  const std::array<std::span<const double>, 3> velocity_fields{
      state.gas_cells.velocity_x_peculiar,
      state.gas_cells.velocity_y_peculiar,
      state.gas_cells.velocity_z_peculiar};
  for (const TopologyFace& topology_face : topology_faces) {
    const physics::MetalDiffusionFace& face = topology_face.face;
    const std::uint8_t axis = topology_face.axis;
    for (std::uint8_t side = 0U; side < 2U; ++side) {
      const std::uint32_t cell =
          side == 0U ? face.left_cell : face.right_cell;
      const std::uint32_t neighbor =
          side == 0U ? face.right_cell : face.left_cell;
      const std::size_t direction = side == 0U ? 1U : 0U;
      NeighborAccumulator& accumulator =
          neighbor_accumulators[cell][axis][direction];
      for (std::size_t component = 0; component < 3U; ++component) {
        accumulator.velocity_area_sum[component] +=
            face.area_code * velocity_fields[component][neighbor];
      }
      accumulator.distance_area_sum +=
          face.area_code * face.center_distance_code;
      accumulator.area_sum += face.area_code;
    }
  }

  for (std::uint32_t cell_index = 0;
       cell_index < state.cells.size(); ++cell_index) {
    if (owned_leaf_mask[cell_index] == 0U) {
      continue;
    }
    for (std::size_t axis = 0; axis < 3U; ++axis) {
      const NeighborAccumulator& lower =
          neighbor_accumulators[cell_index][axis][0];
      const NeighborAccumulator& upper =
          neighbor_accumulators[cell_index][axis][1];
      for (std::size_t component = 0; component < 3U; ++component) {
        const double self_velocity = velocity_fields[component][cell_index];
        double derivative = 0.0;
        if (lower.area_sum > 0.0 && upper.area_sum > 0.0) {
          const double lower_velocity =
              lower.velocity_area_sum[component] / lower.area_sum;
          const double upper_velocity =
              upper.velocity_area_sum[component] / upper.area_sum;
          const double lower_distance =
              lower.distance_area_sum / lower.area_sum;
          const double upper_distance =
              upper.distance_area_sum / upper.area_sum;
          derivative = (upper_velocity - lower_velocity) /
              std::max(lower_distance + upper_distance, 1.0e-30);
        } else if (upper.area_sum > 0.0) {
          derivative =
              (upper.velocity_area_sum[component] / upper.area_sum -
               self_velocity) /
              std::max(
                  upper.distance_area_sum / upper.area_sum, 1.0e-30);
        } else if (lower.area_sum > 0.0) {
          derivative =
              (self_velocity -
               lower.velocity_area_sum[component] / lower.area_sum) /
              std::max(
                  lower.distance_area_sum / lower.area_sum, 1.0e-30);
        }
        diffusion_cells[cell_index].velocity_gradient.grad[component][axis] =
            derivative;
      }
    }
  }

  return result;
}

}  // namespace cosmosim::workflows::internal
