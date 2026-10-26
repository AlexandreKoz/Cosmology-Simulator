// In-place r2c/c2r PM storage regression (2026-10-26).
//
// The PM working set aliases the padded real view (rows of 2 (nz/2 + 1)
// doubles) onto the complex spectrum. These owner-executed checks pin the
// spectral reference operator through that layout for even, odd and
// non-cubic meshes, every Fourier axis, the uniform-density zero mode, linear
// normalization, isolated-boundary antisymmetry, memory ownership, and (under
// MPI) distributed-versus-serial forces on uneven and zero-width slabs.
//
// Tolerances: a single resolved mesh eigenmode is exact for the spectral
// operator (Ψ_k = -4πG δρ_k / k², A = -i k Ψ_k), so differences are FFT
// roundoff; 1e-9 relative is a deliberately loose roundoff bound, not a
// physics accuracy target.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/gravity/pm_solver.hpp"
#include "cosmosim/parallel/distributed_memory.hpp"

#if COSMOSIM_ENABLE_MPI
#include <mpi.h>
#endif

namespace {

using cosmosim::gravity::PmGridShape;
using cosmosim::gravity::PmGridStorage;
using cosmosim::gravity::PmSolveOptions;
using cosmosim::gravity::PmSolver;

constexpr double k_pi = 3.141592653589793238462643383279502884;
constexpr double k_roundoff_rel = 1.0e-9;

void requireOrThrow(bool condition, const std::string& message) {
  if (!condition) throw std::runtime_error(message);
}

[[nodiscard]] std::string shapeName(const PmGridShape& s) {
  return std::to_string(s.nx) + "x" + std::to_string(s.ny) + "x" + std::to_string(s.nz);
}

[[nodiscard]] PmSolveOptions boxOptions(double lx, double ly, double lz, double g) {
  PmSolveOptions options;
  options.box_size_x_mpc_comoving = lx;
  options.box_size_y_mpc_comoving = ly;
  options.box_size_z_mpc_comoving = lz;
  options.gravitational_constant_code = g;
  return options;
}

// δρ = A sin(k_m x_axis) with integer mode m along `axis`; cell centres.
double modeDensity(const PmGridShape& s, const std::array<double, 3>& l, int axis, int m, double amplitude,
                   std::size_t ix, std::size_t iy, std::size_t iz) {
  const std::array<std::size_t, 3> idx{ix, iy, iz};
  const std::array<std::size_t, 3> n{s.nx, s.ny, s.nz};
  const double x = (static_cast<double>(idx[axis]) + 0.5) / static_cast<double>(n[axis]) * l[axis];
  const double k = 2.0 * k_pi * static_cast<double>(m) / l[axis];
  return amplitude * std::sin(k * x);
}

// Analytic spectral response: A_axis = (4πG A / k) cos(k x), transverse 0.
void checkFourierModeResponse(const PmGridShape& s, const std::array<double, 3>& l, int axis, int m, double g) {
  const double amplitude = 0.25;
  PmGridStorage grid(s);
  PmSolver solver(s);
  const PmSolveOptions options = boxOptions(l[0], l[1], l[2], g);
  grid.ensureDensityStorage();
  for (std::size_t ix = 0; ix < s.nx; ++ix)
    for (std::size_t iy = 0; iy < s.ny; ++iy)
      for (std::size_t iz = 0; iz < s.nz; ++iz)
        grid.density()[grid.linearIndex(ix, iy, iz)] = modeDensity(s, l, axis, m, amplitude, ix, iy, iz);
  solver.solvePoissonPeriodic(grid, options, nullptr);
  const std::array<std::size_t, 3> n{s.nx, s.ny, s.nz};
  const double k = 2.0 * k_pi * static_cast<double>(m) / l[axis];
  const double force_amp = 4.0 * k_pi * g * amplitude / k;
  const double phi_amp = -4.0 * k_pi * g * amplitude / (k * k);
  double max_err = 0.0, max_transverse = 0.0, max_phi_err = 0.0;
  for (std::size_t ix = 0; ix < s.nx; ++ix)
    for (std::size_t iy = 0; iy < s.ny; ++iy)
      for (std::size_t iz = 0; iz < s.nz; ++iz) {
        const std::array<std::size_t, 3> idx{ix, iy, iz};
        const double x = (static_cast<double>(idx[axis]) + 0.5) / static_cast<double>(n[axis]) * l[axis];
        const std::size_t c = grid.linearIndex(ix, iy, iz);
        const std::array<double, 3> f{grid.force_x()[c], grid.force_y()[c], grid.force_z()[c]};
        max_err = std::max(max_err, std::abs(f[axis] - force_amp * std::cos(k * x)));
        for (int t = 0; t < 3; ++t) if (t != axis) max_transverse = std::max(max_transverse, std::abs(f[t]));
        max_phi_err = std::max(max_phi_err, std::abs(grid.potential()[c] - phi_amp * std::sin(k * x)));
      }
  const std::string label = "Fourier mode " + shapeName(s) + " axis=" + std::to_string(axis) +
      " m=" + std::to_string(m) + " backend=" + std::string(PmSolver::fftBackendName());
  requireOrThrow(max_err <= k_roundoff_rel * force_amp, label + ": force amplitude/phase error " + std::to_string(max_err));
  requireOrThrow(max_transverse <= k_roundoff_rel * force_amp, label + ": transverse leakage " + std::to_string(max_transverse));
  requireOrThrow(max_phi_err <= k_roundoff_rel * std::abs(phi_amp), label + ": potential error " + std::to_string(max_phi_err));
}

void testFourierModesAcrossLayouts() {
  // Even, odd and mixed-parity nz exercise padding of 2 (even) and 1 (odd)
  // doubles per row; non-cubic boxes exercise per-axis k.
  const std::array<PmGridShape, 4> shapes{PmGridShape{16, 8, 8}, PmGridShape{15, 9, 7},
                                          PmGridShape{12, 10, 6}, PmGridShape{8, 6, 9}};
  const std::array<double, 3> box{1.5, 0.75, 2.25};
  for (const PmGridShape& s : shapes) {
    for (int axis = 0; axis < 3; ++axis) {
      const std::array<std::size_t, 3> n{s.nx, s.ny, s.nz};
      for (const int m : {1, static_cast<int>((n[axis] - 1U) / 2U)}) {
        if (m < 1) continue;
        checkFourierModeResponse(s, box, axis, m, 1.0);
      }
    }
  }
}

void testUniformDensityGivesZeroForce() {
  for (const PmGridShape& s : {PmGridShape{8, 8, 8}, PmGridShape{9, 5, 7}}) {
    PmGridStorage grid(s);
    PmSolver solver(s);
    grid.ensureDensityStorage();
    std::fill(grid.density().begin(), grid.density().end(), 3.7);
    solver.solvePoissonPeriodic(grid, boxOptions(1.0, 1.0, 1.0, 1.0), nullptr);
    double max_force = 0.0;
    for (std::size_t c = 0; c < grid.localCellCount(); ++c) {
      max_force = std::max({max_force, std::abs(grid.force_x()[c]), std::abs(grid.force_y()[c]),
                            std::abs(grid.force_z()[c])});
    }
    requireOrThrow(max_force <= 1.0e-12, "uniform density must give zero periodic force on " + shapeName(s));
  }
}

void testLinearNormalization() {
  // Doubling G and the density amplitude must quadruple the force exactly
  // (linear operator, no hidden per-solve state).
  const PmGridShape s{12, 6, 10};
  const std::array<double, 3> box{1.0, 2.0, 0.5};
  std::vector<double> base;
  for (const double scale : {1.0, 2.0}) {
    PmGridStorage grid(s);
    PmSolver solver(s);
    grid.ensureDensityStorage();
    for (std::size_t ix = 0; ix < s.nx; ++ix)
      for (std::size_t iy = 0; iy < s.ny; ++iy)
        for (std::size_t iz = 0; iz < s.nz; ++iz)
          grid.density()[grid.linearIndex(ix, iy, iz)] =
              scale * (modeDensity(s, box, 0, 1, 0.3, ix, iy, iz) + modeDensity(s, box, 2, 2, 0.1, ix, iy, iz));
    solver.solvePoissonPeriodic(grid, boxOptions(box[0], box[1], box[2], scale), nullptr);
    if (base.empty()) {
      base.assign(grid.force_z().begin(), grid.force_z().end());
    } else {
      for (std::size_t c = 0; c < base.size(); ++c) {
        requireOrThrow(std::abs(grid.force_z()[c] - 4.0 * base[c]) <= 1.0e-12 * std::max(1.0, std::abs(4.0 * base[c])),
                       "PM force must scale linearly with G and density");
      }
    }
  }
}

void testRepeatedSolvesReuseInPlacePlan() {
  // The second solve runs on the cached plan whose real view still holds the
  // previous inverse transform; results must be identical.
  const PmGridShape s{10, 8, 7};
  const std::array<double, 3> box{1.0, 1.0, 1.0};
  PmGridStorage grid(s);
  PmSolver solver(s);
  std::vector<double> first;
  for (int pass = 0; pass < 2; ++pass) {
    grid.ensureDensityStorage();
    for (std::size_t ix = 0; ix < s.nx; ++ix)
      for (std::size_t iy = 0; iy < s.ny; ++iy)
        for (std::size_t iz = 0; iz < s.nz; ++iz)
          grid.density()[grid.linearIndex(ix, iy, iz)] = modeDensity(s, box, 1, 2, 0.2, ix, iy, iz);
    solver.solvePoissonPeriodic(grid, boxOptions(1.0, 1.0, 1.0, 1.0), nullptr);
    if (pass == 0) {
      first.assign(grid.force_y().begin(), grid.force_y().end());
    } else {
      requireOrThrow(std::equal(first.begin(), first.end(), grid.force_y().begin()),
                     "repeated solve on the cached in-place plan must be bitwise reproducible");
    }
  }
}

void testIsolatedBoundaryAntisymmetry() {
  // Isolated/open PM uses its own zero-padded workspace (not the in-place
  // periodic alias). A point mass at a cell centre gives antisymmetric x force.
  const PmGridShape s{16, 16, 16};
  PmGridStorage grid(s);
  PmSolver solver(s);
  PmSolveOptions options = boxOptions(1.0, 1.0, 1.0, 1.0);
  options.boundary_condition = cosmosim::gravity::PmBoundaryCondition::kIsolatedOpen;
  grid.ensureDensityStorage();
  std::fill(grid.density().begin(), grid.density().end(), 0.0);
  const std::size_t c = 8U;
  grid.density()[grid.linearIndex(c, c, c)] = 1.0;
  solver.solvePoissonIsolatedOpen(grid, options, nullptr);
  for (std::size_t d = 1; d < 6; ++d) {
    const double right = grid.force_x()[grid.linearIndex(c + d, c, c)];
    const double left = grid.force_x()[grid.linearIndex(c - d, c, c)];
    requireOrThrow(std::isfinite(right) && std::isfinite(left), "isolated PM force must be finite");
    requireOrThrow(right < 0.0 && left > 0.0, "isolated PM force must point toward the mass");
    requireOrThrow(std::abs(right + left) <= 1.0e-10 * std::abs(right),
                   "isolated PM point-mass force must be antisymmetric");
  }
}

void testInPlaceMemoryOwnership() {
  for (const PmGridShape& s : {PmGridShape{16, 8, 8}, PmGridShape{15, 9, 7}}) {
    const auto layout = cosmosim::parallel::makePmSlabLayout(s.nx, s.ny, s.nz, 1, 0);
    const auto estimate = cosmosim::gravity::estimatePmPlanResourcesMemory(
        s, layout, cosmosim::core::PmDecompositionMode::kSlab);
    requireOrThrow(estimate.real_view_aliases_complex, "in-place estimate must report aliasing");
    requireOrThrow(estimate.real_array_bytes == 0U, "in-place layout owns no separate real array");
    const std::uint64_t complex_cells = s.nx * s.ny * (s.nz / 2U + 1U);
    requireOrThrow(estimate.allocated_local_complex_cells == complex_cells,
                   "serial in-place plan allocates nx*ny*(nz/2+1) complex cells");
    // Padded real view covers nx*ny rows of 2(nz/2+1) doubles exactly.
    requireOrThrow(2U * estimate.allocated_local_complex_cells == s.nx * s.ny * 2U * (s.nz / 2U + 1U),
                   "padded real view must equal the complex allocation");
    const auto fd = cosmosim::gravity::estimatePmPlanResourcesMemory(
        s, layout, cosmosim::core::PmDecompositionMode::kSlab,
        cosmosim::gravity::PmForceGradientBackend::kFiniteDifference4);
    requireOrThrow(fd.finite_difference_halo_bytes == 0U, "serial FD gradient needs no halo planes");
  }
}

#if COSMOSIM_ENABLE_MPI && COSMOSIM_ENABLE_FFTW
// Distributed spectral forces must equal the serial reference on every owned
// cell, including uneven and zero-width slab ownership.
void checkDistributedMatchesSerial(const PmGridShape& s, cosmosim::gravity::PmForceGradientBackend backend) {
  int world_size = 1, world_rank = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  if (world_size < 2) return;
  const std::array<double, 3> box{1.25, 1.0, 0.75};
  PmSolveOptions options = boxOptions(box[0], box[1], box[2], 1.0);
  options.force_gradient_backend = backend;
  const auto density_at = [&](std::size_t ix, std::size_t iy, std::size_t iz) {
    return modeDensity(s, box, 0, 1, 0.2, ix, iy, iz) + modeDensity(s, box, 1, 1, 0.1, ix, iy, iz) +
        0.05 * std::cos(static_cast<double>(ix * 7U + iy * 3U + iz * 5U));
  };
  const auto layout = cosmosim::parallel::makePmSlabLayout(s.nx, s.ny, s.nz, world_size, world_rank);
  PmGridStorage grid(s, layout);
  PmSolver solver(s);
  grid.ensureDensityStorage();
  for (std::size_t lx = 0; lx < layout.local_nx(); ++lx) {
    const std::size_t ix = layout.globalXFromLocal(lx);
    for (std::size_t iy = 0; iy < s.ny; ++iy)
      for (std::size_t iz = 0; iz < s.nz; ++iz) grid.density()[grid.linearIndex(ix, iy, iz)] = density_at(ix, iy, iz);
  }
  solver.solvePoissonPeriodic(grid, options, nullptr);
  PmGridStorage reference(s);
  PmSolver reference_solver(s);
  reference.ensureDensityStorage();
  for (std::size_t ix = 0; ix < s.nx; ++ix)
    for (std::size_t iy = 0; iy < s.ny; ++iy)
      for (std::size_t iz = 0; iz < s.nz; ++iz) reference.density()[reference.linearIndex(ix, iy, iz)] = density_at(ix, iy, iz);
  reference_solver.solvePoissonPeriodic(reference, options, nullptr);
  double max_ref = 0.0, max_err = 0.0;
  for (std::size_t lx = 0; lx < layout.local_nx(); ++lx) {
    const std::size_t ix = layout.globalXFromLocal(lx);
    for (std::size_t iy = 0; iy < s.ny; ++iy)
      for (std::size_t iz = 0; iz < s.nz; ++iz) {
        const std::size_t l = grid.linearIndex(ix, iy, iz);
        const std::size_t r = reference.linearIndex(ix, iy, iz);
        for (const auto& [got, want] : {std::pair{grid.force_x()[l], reference.force_x()[r]},
                                        std::pair{grid.force_y()[l], reference.force_y()[r]},
                                        std::pair{grid.force_z()[l], reference.force_z()[r]}}) {
          max_ref = std::max(max_ref, std::abs(want));
          max_err = std::max(max_err, std::abs(got - want));
        }
      }
  }
  int local_ok = max_err <= k_roundoff_rel * std::max(max_ref, 1.0) ? 1 : 0;
  int all_ok = 0;
  MPI_Allreduce(&local_ok, &all_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  requireOrThrow(all_ok == 1, "distributed in-place PM differs from serial on " + shapeName(s) +
                                  " backend=" + std::to_string(static_cast<int>(backend)) +
                                  " (rank " + std::to_string(world_rank) + " max_err=" +
                                  std::to_string(max_err) + ")");
}

void testDistributedInPlaceLayouts() {
  int world_size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  const auto spectral = cosmosim::gravity::PmForceGradientBackend::kSpectral;
  const auto ws = static_cast<std::size_t>(world_size);
  checkDistributedMatchesSerial(PmGridShape{4U * ws, 6, 8}, spectral);            // even slabs, even nz
  checkDistributedMatchesSerial(PmGridShape{4U * ws + 1U, 5, 7}, spectral);       // uneven slabs, odd nz
  checkDistributedMatchesSerial(PmGridShape{ws + 1U, 6, 6}, spectral);            // one extra plane
  if (world_size >= 4) checkDistributedMatchesSerial(PmGridShape{3, 8, 8}, spectral);  // zero-width rank
  // FD gradient: every rank owns >= 2 planes; compared against its own
  // serial FD reference (not the spectral one — distinct operator).
  checkDistributedMatchesSerial(PmGridShape{4U * ws, 6, 8},
                                cosmosim::gravity::PmForceGradientBackend::kFiniteDifference4);
}
#endif

}  // namespace

int main() {
#if COSMOSIM_ENABLE_MPI
  MPI_Init(nullptr, nullptr);
#endif
  int failed = 0;
  try {
    // Serial-layout solves are executed by every rank: under MPI the PM API
    // still performs its collective entry checks, so rank-gated serial work
    // would itself diverge.
    testFourierModesAcrossLayouts();
    testUniformDensityGivesZeroForce();
    testLinearNormalization();
    testRepeatedSolvesReuseInPlacePlan();
    testIsolatedBoundaryAntisymmetry();
    testInPlaceMemoryOwnership();
#if COSMOSIM_ENABLE_MPI && COSMOSIM_ENABLE_FFTW
    testDistributedInPlaceLayouts();
#endif
  } catch (const std::exception& error) {
    std::fprintf(stderr, "in-place FFT layout test failed: %s\n", error.what());
    failed = 1;
  }
#if COSMOSIM_ENABLE_MPI
  int any_failed = 0;
  MPI_Allreduce(&failed, &any_failed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  failed = any_failed;
  MPI_Finalize();
#endif
  return failed == 0 ? 0 : 1;
}
