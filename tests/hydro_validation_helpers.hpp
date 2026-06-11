#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

#include "cosmosim/hydro/hydro_cartesian_patch.hpp"
#include "cosmosim/hydro/hydro_core_solver.hpp"

namespace cosmosim::tests::hydro_validation {

constexpr double k_pi = 3.141592653589793238462643383279502884;

struct CellCenter {
  double x_comoving = 0.0;
  double y_comoving = 0.0;
  double z_comoving = 0.0;
};

struct RadialBin {
  double radius_min_code = 0.0;
  double radius_max_code = 0.0;
  double mean_density = 0.0;
  double mean_pressure = 0.0;
  double mean_radial_velocity = 0.0;
  double density_variance = 0.0;
  std::size_t count = 0;
};

struct HydroStepConfig {
  double gamma = 5.0 / 3.0;
  double dt_code = 1.0e-4;
  std::size_t step_count = 1;
  double rho_floor = 1.0e-10;
  double pressure_floor = 1.0e-10;
  hydro::HydroSlopeLimiter limiter = hydro::HydroSlopeLimiter::kMonotonizedCentral;
};

inline void requireOrThrow(bool condition, const std::string& message) {
  if (!condition) {
    throw std::runtime_error(message);
  }
}

[[nodiscard]] inline hydro::HydroPatchGeometry makePeriodicCartesianPatch(
    std::size_t nx,
    std::size_t ny,
    std::size_t nz) {
  hydro::HydroPatchGeometry geometry = hydro::makeCartesianPatchGeometry(hydro::HydroCartesianPatchSpec{
      .nx = nx,
      .ny = ny,
      .nz = nz,
      .origin_x_comoving = 0.0,
      .origin_y_comoving = 0.0,
      .origin_z_comoving = 0.0,
      .cell_width_x_comoving = 1.0 / static_cast<double>(nx),
      .cell_width_y_comoving = 1.0 / static_cast<double>(ny),
      .cell_width_z_comoving = 1.0 / static_cast<double>(nz)});

  if (nx > 1U) {
    for (std::size_t k = 0; k < nz; ++k) {
      for (std::size_t j = 0; j < ny; ++j) {
        geometry.faces.push_back(hydro::HydroFace{
            .owner_cell = geometry.linearCellIndex(nx - 1U, j, k),
            .neighbor_cell = geometry.linearCellIndex(0U, j, k),
            .owner_minus_cell = geometry.linearCellIndex(nx - 2U, j, k),
            .neighbor_plus_cell = geometry.linearCellIndex(1U, j, k),
            .area_comoving = geometry.cell_width_y_comoving * geometry.cell_width_z_comoving,
            .normal_x = 1.0,
            .axis = hydro::HydroFaceAxis::kX});
      }
    }
  }
  if (ny > 1U) {
    for (std::size_t k = 0; k < nz; ++k) {
      for (std::size_t i = 0; i < nx; ++i) {
        geometry.faces.push_back(hydro::HydroFace{
            .owner_cell = geometry.linearCellIndex(i, ny - 1U, k),
            .neighbor_cell = geometry.linearCellIndex(i, 0U, k),
            .owner_minus_cell = geometry.linearCellIndex(i, ny - 2U, k),
            .neighbor_plus_cell = geometry.linearCellIndex(i, 1U, k),
            .area_comoving = geometry.cell_width_x_comoving * geometry.cell_width_z_comoving,
            .normal_y = 1.0,
            .axis = hydro::HydroFaceAxis::kY});
      }
    }
  }
  if (nz > 1U) {
    for (std::size_t j = 0; j < ny; ++j) {
      for (std::size_t i = 0; i < nx; ++i) {
        geometry.faces.push_back(hydro::HydroFace{
            .owner_cell = geometry.linearCellIndex(i, j, nz - 1U),
            .neighbor_cell = geometry.linearCellIndex(i, j, 0U),
            .owner_minus_cell = geometry.linearCellIndex(i, j, nz - 2U),
            .neighbor_plus_cell = geometry.linearCellIndex(i, j, 1U),
            .area_comoving = geometry.cell_width_x_comoving * geometry.cell_width_y_comoving,
            .normal_z = 1.0,
            .axis = hydro::HydroFaceAxis::kZ});
      }
    }
  }
  return geometry;
}

[[nodiscard]] inline std::vector<std::size_t> allRealCells(const hydro::HydroPatchGeometry& geometry) {
  std::vector<std::size_t> cells(geometry.cellCount());
  for (std::size_t i = 0; i < cells.size(); ++i) {
    cells[i] = i;
  }
  return cells;
}

[[nodiscard]] inline CellCenter cellCenter(
    const hydro::HydroPatchGeometry& geometry,
    std::size_t cell_index) {
  const std::array<std::size_t, 3> ijk = geometry.cellIjk(cell_index);
  return CellCenter{
      .x_comoving = geometry.origin_x_comoving +
          (static_cast<double>(ijk[0]) + 0.5) * geometry.cell_width_x_comoving,
      .y_comoving = geometry.origin_y_comoving +
          (static_cast<double>(ijk[1]) + 0.5) * geometry.cell_width_y_comoving,
      .z_comoving = geometry.origin_z_comoving +
          (static_cast<double>(ijk[2]) + 0.5) * geometry.cell_width_z_comoving};
}

inline void storePrimitive(
    hydro::HydroConservedStateSoa& conserved,
    std::size_t cell_index,
    const hydro::HydroPrimitiveState& primitive,
    double gamma) {
  conserved.storeCell(cell_index, hydro::HydroCoreSolver::conservedFromPrimitive(primitive, gamma));
}

inline void requireFinitePositiveState(
    const hydro::HydroConservedStateSoa& conserved,
    double gamma,
    const std::string& case_name) {
  for (std::size_t cell = 0; cell < conserved.size(); ++cell) {
    const hydro::HydroConservedState state = conserved.loadCell(cell);
    requireOrThrow(std::isfinite(state.mass_density_comoving), case_name + ": non-finite density");
    requireOrThrow(std::isfinite(state.momentum_density_x_comoving), case_name + ": non-finite momentum x");
    requireOrThrow(std::isfinite(state.momentum_density_y_comoving), case_name + ": non-finite momentum y");
    requireOrThrow(std::isfinite(state.momentum_density_z_comoving), case_name + ": non-finite momentum z");
    requireOrThrow(std::isfinite(state.total_energy_density_comoving), case_name + ": non-finite total energy");
    requireOrThrow(state.mass_density_comoving > 0.0, case_name + ": non-positive density");
    const hydro::HydroPrimitiveState primitive = hydro::HydroCoreSolver::primitiveFromConserved(state, gamma);
    const double kinetic_density = 0.5 *
        (state.momentum_density_x_comoving * state.momentum_density_x_comoving +
         state.momentum_density_y_comoving * state.momentum_density_y_comoving +
         state.momentum_density_z_comoving * state.momentum_density_z_comoving) /
        state.mass_density_comoving;
    const double internal_energy_density = state.total_energy_density_comoving - kinetic_density;
    requireOrThrow(std::isfinite(primitive.pressure_comoving), case_name + ": non-finite pressure");
    requireOrThrow(primitive.pressure_comoving > 0.0, case_name + ": non-positive pressure");
    requireOrThrow(internal_energy_density > 0.0, case_name + ": non-positive internal energy");
  }
}

[[nodiscard]] inline hydro::HydroConservationTotals conservedTotals(
    const hydro::HydroConservedStateSoa& conserved,
    const hydro::HydroPatchGeometry& geometry) {
  const std::vector<std::size_t> cells = allRealCells(geometry);
  return hydro::HydroCoreSolver::conservationTotals(conserved, geometry, cells);
}

inline void advanceHydroSteps(
    hydro::HydroConservedStateSoa& conserved,
    const hydro::HydroPatchGeometry& geometry,
    const HydroStepConfig& config,
    std::span<const hydro::HydroSourceTerm* const> source_terms = {},
    hydro::HydroSourceContext source_context = {}) {
  hydro::HydroUpdateContext update{
      .dt_code = config.dt_code,
      .scale_factor = 1.0,
      .hubble_rate_code = 0.0};
  source_context.update = update;
  hydro::HydroCoreSolver solver(config.gamma);
  hydro::MusclHancockReconstruction reconstruction(hydro::HydroReconstructionPolicy{
      .limiter = config.limiter,
      .dt_over_dx_code = config.dt_code / geometry.cell_width_x_comoving,
      .dt_over_cell_width_code = {
          config.dt_code / geometry.cell_width_x_comoving,
          config.dt_code / geometry.cell_width_y_comoving,
          config.dt_code / geometry.cell_width_z_comoving},
      .rho_floor = config.rho_floor,
      .pressure_floor = config.pressure_floor,
      .enable_muscl_hancock_predictor = true,
      .adiabatic_index = config.gamma});
  hydro::HllcRiemannSolver riemann;
  hydro::HydroScratchBuffers scratch;
  hydro::HydroPrimitiveCacheSoa primitive_cache(conserved.size());
  for (std::size_t step = 0; step < config.step_count; ++step) {
    solver.advancePatchWithScratch(
        conserved,
        geometry,
        update,
        reconstruction,
        riemann,
        source_terms,
        source_context,
        scratch,
        &primitive_cache,
        nullptr);
  }
}

[[nodiscard]] inline std::vector<RadialBin> radialBins(
    const hydro::HydroConservedStateSoa& conserved,
    const hydro::HydroPatchGeometry& geometry,
    double gamma,
    std::span<const double> edges) {
  requireOrThrow(edges.size() >= 2U, "radialBins requires at least two edges");
  std::vector<RadialBin> bins(edges.size() - 1U);
  std::vector<double> density_square_sum(bins.size(), 0.0);
  for (std::size_t i = 0; i < bins.size(); ++i) {
    bins[i].radius_min_code = edges[i];
    bins[i].radius_max_code = edges[i + 1U];
  }

  for (std::size_t cell = 0; cell < geometry.cellCount(); ++cell) {
    const CellCenter center = cellCenter(geometry, cell);
    const double dx = center.x_comoving - 0.5;
    const double dy = center.y_comoving - 0.5;
    const double dz = center.z_comoving - 0.5;
    const double radius = std::sqrt(dx * dx + dy * dy + dz * dz);
    auto it = std::upper_bound(edges.begin(), edges.end(), radius);
    if (it == edges.begin() || it == edges.end()) {
      continue;
    }
    const std::size_t bin_index = static_cast<std::size_t>(it - edges.begin() - 1);
    hydro::HydroPrimitiveState primitive =
        hydro::HydroCoreSolver::primitiveFromConserved(conserved.loadCell(cell), gamma);
    const double inv_radius = radius > 1.0e-14 ? 1.0 / radius : 0.0;
    const double vr = (primitive.vel_x_peculiar * dx + primitive.vel_y_peculiar * dy +
                       primitive.vel_z_peculiar * dz) *
        inv_radius;
    bins[bin_index].mean_density += primitive.rho_comoving;
    bins[bin_index].mean_pressure += primitive.pressure_comoving;
    bins[bin_index].mean_radial_velocity += vr;
    density_square_sum[bin_index] += primitive.rho_comoving * primitive.rho_comoving;
    ++bins[bin_index].count;
  }

  for (std::size_t i = 0; i < bins.size(); ++i) {
    if (bins[i].count == 0U) {
      continue;
    }
    const double inv_count = 1.0 / static_cast<double>(bins[i].count);
    bins[i].mean_density *= inv_count;
    bins[i].mean_pressure *= inv_count;
    bins[i].mean_radial_velocity *= inv_count;
    bins[i].density_variance =
        density_square_sum[i] * inv_count - bins[i].mean_density * bins[i].mean_density;
  }
  return bins;
}

}  // namespace cosmosim::tests::hydro_validation
