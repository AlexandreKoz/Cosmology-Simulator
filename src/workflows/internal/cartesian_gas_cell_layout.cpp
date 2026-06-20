#include "workflows/internal/cartesian_gas_cell_layout.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <sstream>
#include <unordered_map>
#include <utility>
#include <vector>

namespace cosmosim::workflows::internal {
namespace {

constexpr double k_coordinate_relative_tolerance = 1.0e-10;
constexpr std::uint32_t k_invalid_row = std::numeric_limits<std::uint32_t>::max();

[[nodiscard]] bool nearlyEqual(double lhs, double rhs) {
  const double scale = std::max({1.0, std::abs(lhs), std::abs(rhs)});
  return std::abs(lhs - rhs) <= k_coordinate_relative_tolerance * scale;
}

[[nodiscard]] bool finitePositive(double value) {
  return std::isfinite(value) && value > 0.0;
}

[[nodiscard]] std::optional<std::vector<double>> sortedUniqueCoordinates(
    const cosmosim::core::AlignedVector<double>& values,
    std::string* diagnostic) {
  std::vector<double> sorted(values.begin(), values.end());
  if (std::any_of(sorted.begin(), sorted.end(), [](double value) { return !std::isfinite(value); })) {
    *diagnostic = "gas-cell center coordinates contain a non-finite value";
    return std::nullopt;
  }
  std::sort(sorted.begin(), sorted.end());
  std::vector<double> unique;
  unique.reserve(sorted.size());
  for (const double value : sorted) {
    if (unique.empty() || !nearlyEqual(unique.back(), value)) {
      unique.push_back(value);
    }
  }
  return unique;
}

[[nodiscard]] std::optional<std::size_t> uniqueCoordinateIndex(
    const std::vector<double>& coordinates,
    double value,
    std::string* diagnostic,
    const char* axis) {
  std::optional<std::size_t> match;
  for (std::size_t index = 0; index < coordinates.size(); ++index) {
    if (!nearlyEqual(coordinates[index], value)) {
      continue;
    }
    if (match.has_value()) {
      *diagnostic = std::string("ambiguous ") + axis + " coordinate match while mapping a gas-cell center";
      return std::nullopt;
    }
    match = index;
  }
  if (!match.has_value()) {
    *diagnostic = std::string("off-lattice ") + axis + " coordinate while mapping a gas-cell center";
  }
  return match;
}

[[nodiscard]] std::optional<double> uniformSpacing(
    const std::vector<double>& coordinates,
    double singleton_width,
    std::string* diagnostic,
    const char* axis) {
  if (coordinates.empty()) {
    *diagnostic = std::string("empty ") + axis + " coordinate axis";
    return std::nullopt;
  }
  if (coordinates.size() == 1U) {
    if (!finitePositive(singleton_width)) {
      *diagnostic = std::string("single-cell ") + axis + " axis has no positive physical width";
      return std::nullopt;
    }
    return singleton_width;
  }
  const double expected = coordinates[1] - coordinates[0];
  if (!finitePositive(expected)) {
    *diagnostic = std::string("non-increasing ") + axis + " coordinate axis";
    return std::nullopt;
  }
  for (std::size_t index = 2; index < coordinates.size(); ++index) {
    const double actual = coordinates[index] - coordinates[index - 1U];
    if (!nearlyEqual(actual, expected)) {
      *diagnostic = std::string("nonuniform ") + axis + " coordinate spacing in fixed Cartesian hydro layout";
      return std::nullopt;
    }
  }
  return expected;
}

[[nodiscard]] std::optional<std::size_t> patchIndexById(
    const cosmosim::core::SimulationState& state,
    std::uint64_t patch_id) {
  for (std::size_t patch_index = 0; patch_index < state.patches.size(); ++patch_index) {
    if (state.patches.patch_id[patch_index] == patch_id) {
      return patch_index;
    }
  }
  return std::nullopt;
}

[[nodiscard]] bool hasExplicitPatchGeometry(
    const cosmosim::core::SimulationState& state,
    std::size_t patch_index) {
  return patch_index < state.patches.size() &&
      state.patches.cell_dim_x[patch_index] > 0U &&
      state.patches.cell_dim_y[patch_index] > 0U &&
      state.patches.cell_dim_z[patch_index] > 0U &&
      finitePositive(state.patches.extent_x_comoving[patch_index]) &&
      finitePositive(state.patches.extent_y_comoving[patch_index]) &&
      finitePositive(state.patches.extent_z_comoving[patch_index]);
}

[[nodiscard]] bool hasAnyExplicitPatchGeometry(
    const cosmosim::core::SimulationState& state) {
  for (std::size_t patch_index = 0; patch_index < state.patches.size(); ++patch_index) {
    if (hasExplicitPatchGeometry(state, patch_index)) {
      return true;
    }
  }
  return false;
}

[[nodiscard]] std::optional<cosmosim::hydro::HydroCartesianPatchSpec> explicitSinglePatchSpec(
    const cosmosim::core::SimulationState& state,
    std::string* diagnostic) {
  if (state.patches.size() != 1U || !hasExplicitPatchGeometry(state, 0U)) {
    return std::nullopt;
  }
  const std::size_t nx = state.patches.cell_dim_x[0];
  const std::size_t ny = state.patches.cell_dim_y[0];
  const std::size_t nz = state.patches.cell_dim_z[0];
  if (nx > std::numeric_limits<std::size_t>::max() / ny || nx * ny > std::numeric_limits<std::size_t>::max() / nz ||
      nx * ny * nz != state.cells.size()) {
    *diagnostic = "explicit Cartesian patch dimensions do not cover exactly the dense gas-cell rows";
    return std::nullopt;
  }
  return cosmosim::hydro::HydroCartesianPatchSpec{
      .nx = nx,
      .ny = ny,
      .nz = nz,
      .origin_x_comoving = state.patches.origin_x_comoving[0],
      .origin_y_comoving = state.patches.origin_y_comoving[0],
      .origin_z_comoving = state.patches.origin_z_comoving[0],
      .cell_width_x_comoving = state.patches.extent_x_comoving[0] / static_cast<double>(nx),
      .cell_width_y_comoving = state.patches.extent_y_comoving[0] / static_cast<double>(ny),
      .cell_width_z_comoving = state.patches.extent_z_comoving[0] / static_cast<double>(nz),
  };
}

[[nodiscard]] std::optional<std::array<std::size_t, 3>> centerIndicesForSpec(
    const cosmosim::hydro::HydroCartesianPatchSpec& spec,
    double x,
    double y,
    double z,
    std::string* diagnostic) {
  const auto axis_index = [&](double coordinate, double origin, double width, std::size_t count, const char* axis)
      -> std::optional<std::size_t> {
    const double scaled = (coordinate - origin) / width - 0.5;
    const long long rounded = std::llround(scaled);
    if (rounded < 0 || rounded >= static_cast<long long>(count) ||
        !nearlyEqual(coordinate, origin + (static_cast<double>(rounded) + 0.5) * width)) {
      *diagnostic = std::string("explicit Cartesian patch geometry is inconsistent with a gas-cell ") + axis + " center";
      return std::nullopt;
    }
    return static_cast<std::size_t>(rounded);
  };
  const auto i = axis_index(x, spec.origin_x_comoving, spec.cell_width_x_comoving, spec.nx, "x");
  const auto j = axis_index(y, spec.origin_y_comoving, spec.cell_width_y_comoving, spec.ny, "y");
  const auto k = axis_index(z, spec.origin_z_comoving, spec.cell_width_z_comoving, spec.nz, "z");
  if (!i.has_value() || !j.has_value() || !k.has_value()) {
    return std::nullopt;
  }
  return std::array<std::size_t, 3>{*i, *j, *k};
}

void mixSignature(std::uint64_t* signature, std::uint64_t value) {
  *signature ^= value;
  *signature *= 1099511628211ULL;
}

}  // namespace

CartesianGasCellLayoutBuildResult buildCartesianGasCellRowLayout(
    const core::SimulationState& state,
    const core::SimulationConfig& config) {
  CartesianGasCellLayoutBuildResult result;
  const std::size_t cell_count = state.cells.size();
  if (cell_count == 0U) {
    result.diagnostic = "fixed Cartesian hydro layout requires at least one gas cell";
    return result;
  }
  if (state.cells.center_x_comoving.size() != cell_count ||
      state.cells.center_y_comoving.size() != cell_count ||
      state.cells.center_z_comoving.size() != cell_count) {
    result.diagnostic = "gas-cell center lanes do not cover every dense gas-cell row";
    return result;
  }
  try {
    state.requireGasCellIdentityMapCoversDenseRows("fixed Cartesian hydro layout construction");
  } catch (const std::exception& exception) {
    result.diagnostic = exception.what();
    return result;
  }

  std::string diagnostic;
  const auto explicit_spec = explicitSinglePatchSpec(state, &diagnostic);
  if (!diagnostic.empty()) {
    result.diagnostic = std::move(diagnostic);
    return result;
  }
  // Explicit multi-patch geometry belongs to the AMR production path. Never
  // ignore it and silently reinterpret the same rows as a single fixed patch.
  if (!explicit_spec.has_value() && hasAnyExplicitPatchGeometry(state)) {
    result.diagnostic =
        "fixed Cartesian hydro cannot ignore partial or multi-patch explicit geometry; "
        "require complete AMR coverage or one internally consistent fixed patch";
    return result;
  }

  hydro::HydroCartesianPatchSpec spec;
  std::vector<double> x_coordinates;
  std::vector<double> y_coordinates;
  std::vector<double> z_coordinates;
  if (explicit_spec.has_value()) {
    spec = *explicit_spec;
  } else {
    const auto x_axis = sortedUniqueCoordinates(state.cells.center_x_comoving, &diagnostic);
    const auto y_axis = sortedUniqueCoordinates(state.cells.center_y_comoving, &diagnostic);
    const auto z_axis = sortedUniqueCoordinates(state.cells.center_z_comoving, &diagnostic);
    if (!x_axis.has_value() || !y_axis.has_value() || !z_axis.has_value()) {
      result.diagnostic = std::move(diagnostic);
      return result;
    }
    x_coordinates = *x_axis;
    y_coordinates = *y_axis;
    z_coordinates = *z_axis;
    if (x_coordinates.size() > std::numeric_limits<std::size_t>::max() / y_coordinates.size() ||
        x_coordinates.size() * y_coordinates.size() > std::numeric_limits<std::size_t>::max() / z_coordinates.size() ||
        x_coordinates.size() * y_coordinates.size() * z_coordinates.size() != cell_count) {
      result.diagnostic = "gas-cell centers do not form a complete Cartesian product; missing or duplicate physical occupancy";
      return result;
    }
    const auto dx = uniformSpacing(x_coordinates, config.cosmology.box_size_x_mpc_comoving, &diagnostic, "x");
    const auto dy = uniformSpacing(y_coordinates, config.cosmology.box_size_y_mpc_comoving, &diagnostic, "y");
    const auto dz = uniformSpacing(z_coordinates, config.cosmology.box_size_z_mpc_comoving, &diagnostic, "z");
    if (!dx.has_value() || !dy.has_value() || !dz.has_value()) {
      result.diagnostic = std::move(diagnostic);
      return result;
    }
    spec = hydro::HydroCartesianPatchSpec{
        .nx = x_coordinates.size(), .ny = y_coordinates.size(), .nz = z_coordinates.size(),
        .origin_x_comoving = x_coordinates.front() - 0.5 * *dx,
        .origin_y_comoving = y_coordinates.front() - 0.5 * *dy,
        .origin_z_comoving = z_coordinates.front() - 0.5 * *dz,
        .cell_width_x_comoving = *dx,
        .cell_width_y_comoving = *dy,
        .cell_width_z_comoving = *dz};
  }

  hydro::HydroPatchGeometry geometry;
  try {
    geometry = hydro::makeCartesianPatchGeometry(spec);
  } catch (const std::exception& exception) {
    result.diagnostic = std::string("invalid fixed Cartesian patch specification: ") + exception.what();
    return result;
  }
  if (geometry.cellCount() != cell_count) {
    result.diagnostic = "fixed Cartesian patch cell count disagrees with dense gas-cell coverage";
    return result;
  }

  result.layout.spec = spec;
  result.layout.dense_row_by_geometry_row.assign(cell_count, k_invalid_row);
  result.layout.geometry_row_by_dense_row.assign(cell_count, k_invalid_row);
  for (std::uint32_t dense_row = 0; dense_row < cell_count; ++dense_row) {
    const auto* identity = state.gas_cell_identity.findByLocalRow(dense_row);
    if (identity == nullptr || identity->gas_cell_id == 0U) {
      result.diagnostic = "gas-cell identity map is incomplete while building fixed Cartesian hydro layout";
      return result;
    }
    if (state.patches.size() != 0U) {
      const auto patch_index = patchIndexById(state, identity->owning_patch_id);
      if (!patch_index.has_value()) {
        result.diagnostic = "gas-cell identity owns a patch absent from PatchSoa";
        return result;
      }
      if (state.cells.patch_index[dense_row] != *patch_index) {
        result.diagnostic = "gas-cell identity patch ownership disagrees with the dense cell patch-index mirror";
        return result;
      }
    }

    std::optional<std::array<std::size_t, 3>> indices;
    if (explicit_spec.has_value()) {
      indices = centerIndicesForSpec(
          spec,
          state.cells.center_x_comoving[dense_row],
          state.cells.center_y_comoving[dense_row],
          state.cells.center_z_comoving[dense_row],
          &diagnostic);
    } else {
      const auto i = uniqueCoordinateIndex(x_coordinates, state.cells.center_x_comoving[dense_row], &diagnostic, "x");
      const auto j = uniqueCoordinateIndex(y_coordinates, state.cells.center_y_comoving[dense_row], &diagnostic, "y");
      const auto k = uniqueCoordinateIndex(z_coordinates, state.cells.center_z_comoving[dense_row], &diagnostic, "z");
      if (i.has_value() && j.has_value() && k.has_value()) {
        indices = std::array<std::size_t, 3>{*i, *j, *k};
      }
    }
    if (!indices.has_value()) {
      result.diagnostic = diagnostic.empty() ? "failed to locate a gas-cell center on the Cartesian lattice" : diagnostic;
      return result;
    }
    const std::size_t geometry_row = geometry.linearCellIndex((*indices)[0], (*indices)[1], (*indices)[2]);
    if (result.layout.dense_row_by_geometry_row[geometry_row] != k_invalid_row) {
      result.diagnostic = "duplicate physical Cartesian occupancy while mapping dense gas-cell rows";
      return result;
    }
    result.layout.dense_row_by_geometry_row[geometry_row] = dense_row;
    result.layout.geometry_row_by_dense_row[dense_row] = static_cast<std::uint32_t>(geometry_row);
  }
  if (std::any_of(result.layout.dense_row_by_geometry_row.begin(), result.layout.dense_row_by_geometry_row.end(),
                  [](std::uint32_t row) { return row == k_invalid_row; })) {
    result.diagnostic = "missing physical Cartesian cell after mapping dense gas-cell rows";
    return result;
  }

  std::uint64_t signature = 1469598103934665603ULL;
  mixSignature(&signature, static_cast<std::uint64_t>(cell_count));
  for (const std::uint32_t row : result.layout.dense_row_by_geometry_row) {
    mixSignature(&signature, row);
  }
  result.layout.mapping_signature = signature;
  return result;
}

}  // namespace cosmosim::workflows::internal
