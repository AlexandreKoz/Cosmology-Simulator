// 20261012 regression repairs: NOT RUN — deferred to maintainer.
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstdint>
#include <iostream>
#include <limits>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/gravity/tree_pm_coupling.hpp"

#if defined(_OPENMP)
#include <omp.h>
#endif

#if COSMOSIM_ENABLE_MPI
#include <mpi.h>
#endif

namespace {

void requireOrThrow(bool condition, const std::string& message) {
  if (!condition) {
    throw std::runtime_error(message);
  }
}

[[nodiscard]] constexpr cosmosim::gravity::PmGridShape pmShapeForAvailableBackend(
    std::size_t fftw_n,
    std::size_t fallback_n) noexcept {
#if COSMOSIM_ENABLE_FFTW
  return {fftw_n, fftw_n, fftw_n};
#else
  return {fallback_n, fallback_n, fallback_n};
#endif
}

[[nodiscard]] constexpr double rcutCellsForAvailableBackend(
    double fftw_rcut_cells,
    double fallback_rcut_cells) noexcept {
  // Naive-DFT fixtures intentionally use tiny meshes. Keep their diagnostic
  // split inside the production minimum-image cutoff domain without raising
  // the O(N_mesh^2) fallback cost; FFTW fixtures retain their named profile.
#if COSMOSIM_ENABLE_FFTW
  return fftw_rcut_cells;
#else
  return fallback_rcut_cells;
#endif
}

[[nodiscard]] double minimumImageDelta(double delta, double box_size_comoving) {
  return delta - box_size_comoving * std::nearbyint(delta / box_size_comoving);
}

struct ForceField {
  std::vector<double> ax;
  std::vector<double> ay;
  std::vector<double> az;
};

enum class ForceEvaluation {
  kFullTreePm,
  kShortRangeOnly,
};

// Full forces refresh PM through the ordinary API. Split-only fixtures must
// explicitly select the no-refresh surface; cache-lifetime tests below call
// solveActiveSetWithPmCadence directly on their retained coordinator.
void solveForEvaluation(
    cosmosim::gravity::TreePmCoordinator& coordinator,
    std::span<const double> pos_x, std::span<const double> pos_y,
    std::span<const double> pos_z, std::span<const double> mass,
    const cosmosim::gravity::TreePmForceAccumulatorView& accumulator,
    const cosmosim::gravity::TreePmOptions& options, ForceEvaluation evaluation,
    cosmosim::gravity::TreePmProfileEvent* profile = nullptr,
    cosmosim::gravity::TreePmDiagnostics* diagnostics = nullptr,
    const cosmosim::gravity::TreeSofteningView& softening = {}) {
  requireOrThrow(options.short_range_only == (evaluation == ForceEvaluation::kShortRangeOnly),
      "fixture force surface disagrees with short_range_only");
  if (evaluation == ForceEvaluation::kShortRangeOnly) {
    coordinator.solveActiveSetWithPmCadence(pos_x, pos_y, pos_z, mass,
        accumulator, options, false, profile, diagnostics, softening);
  } else {
    coordinator.solveActiveSet(pos_x, pos_y, pos_z, mass,
        accumulator, options, profile, diagnostics, softening);
  }
}

struct AxisBoxLengths {
  double x = 1.0;
  double y = 1.0;
  double z = 1.0;
};

[[nodiscard]] ForceField computeDirectReference(
    std::span<const double> pos_x,
    std::span<const double> pos_y,
    std::span<const double> pos_z,
    std::span<const double> mass,
    double box_size_comoving,
    const cosmosim::gravity::TreeSofteningPolicy& softening,
    double gravitational_constant_code) {
  ForceField field{
      std::vector<double>(pos_x.size(), 0.0),
      std::vector<double>(pos_x.size(), 0.0),
      std::vector<double>(pos_x.size(), 0.0),
  };

  for (std::size_t target = 0; target < pos_x.size(); ++target) {
    for (std::size_t source = 0; source < pos_x.size(); ++source) {
      if (target == source) {
        continue;
      }
      const double dx = minimumImageDelta(pos_x[source] - pos_x[target], box_size_comoving);
      const double dy = minimumImageDelta(pos_y[source] - pos_y[target], box_size_comoving);
      const double dz = minimumImageDelta(pos_z[source] - pos_z[target], box_size_comoving);
      const double r2 = dx * dx + dy * dy + dz * dz;
      const double softened = cosmosim::gravity::softenedInvR3(r2, softening);
      const double factor = gravitational_constant_code * mass[source] * softened;
      field.ax[target] += factor * dx;
      field.ay[target] += factor * dy;
      field.az[target] += factor * dz;
    }
  }

  return field;
}

[[nodiscard]] double relativeL2Error(const ForceField& value, const ForceField& reference) {
  double ref_norm = 0.0;
  double err_norm = 0.0;
  for (std::size_t i = 0; i < value.ax.size(); ++i) {
    const double dax = value.ax[i] - reference.ax[i];
    const double day = value.ay[i] - reference.ay[i];
    const double daz = value.az[i] - reference.az[i];
    err_norm += dax * dax + day * day + daz * daz;
    ref_norm += reference.ax[i] * reference.ax[i] + reference.ay[i] * reference.ay[i] + reference.az[i] * reference.az[i];
  }
  return std::sqrt(err_norm / std::max(ref_norm, 1.0e-30));
}

[[nodiscard]] ForceField solveTreePm(
    std::span<const double> pos_x,
    std::span<const double> pos_y,
    std::span<const double> pos_z,
    std::span<const double> mass,
    const cosmosim::gravity::PmGridShape& pm_shape,
    cosmosim::gravity::TreePmOptions options,
    cosmosim::gravity::TreePmDiagnostics* diagnostics,
    cosmosim::gravity::TreePmProfileEvent* profile = nullptr,
    std::span<const double> previous_acceleration_magnitude_code = {},
    ForceEvaluation evaluation = ForceEvaluation::kFullTreePm) {
  std::vector<std::uint32_t> active(pos_x.size(), 0U);
  for (std::size_t i = 0; i < pos_x.size(); ++i) {
    active[i] = static_cast<std::uint32_t>(i);
  }

  ForceField field{
      std::vector<double>(pos_x.size(), 0.0),
      std::vector<double>(pos_x.size(), 0.0),
      std::vector<double>(pos_x.size(), 0.0),
  };

  cosmosim::gravity::TreePmForceAccumulatorView accumulator{
      active,
      field.ax,
      field.ay,
      field.az,
      previous_acceleration_magnitude_code,
  };

  cosmosim::gravity::TreePmCoordinator coordinator(pm_shape);
  solveForEvaluation(coordinator, pos_x, pos_y, pos_z, mass, accumulator,
      options, evaluation, profile, diagnostics);
  return field;
}

[[nodiscard]] ForceField solveTreePmActiveSubset(
    std::span<const double> pos_x,
    std::span<const double> pos_y,
    std::span<const double> pos_z,
    std::span<const double> mass,
    std::span<const std::uint32_t> active,
    const cosmosim::gravity::PmGridShape& pm_shape,
    cosmosim::gravity::TreePmOptions options,
    cosmosim::gravity::TreePmDiagnostics* diagnostics,
    cosmosim::gravity::TreePmProfileEvent* profile = nullptr,
    std::span<const double> previous_acceleration_magnitude_code = {}) {
  ForceField field{
      std::vector<double>(active.size(), 0.0),
      std::vector<double>(active.size(), 0.0),
      std::vector<double>(active.size(), 0.0),
  };
  cosmosim::gravity::TreePmForceAccumulatorView accumulator{
      active, field.ax, field.ay, field.az, previous_acceleration_magnitude_code};
  cosmosim::gravity::TreePmCoordinator coordinator(pm_shape);
  solveForEvaluation(coordinator, pos_x, pos_y, pos_z, mass, accumulator,
      options, ForceEvaluation::kFullTreePm, profile, diagnostics);
  return field;
}

void requireEquivalentForces(
    const ForceField& lhs,
    const ForceField& rhs,
    double tolerance,
    const std::string& context) {
  requireOrThrow(lhs.ax.size() == rhs.ax.size() && lhs.ay.size() == lhs.ax.size() &&
      lhs.az.size() == lhs.ax.size() && rhs.ay.size() == rhs.ax.size() && rhs.az.size() == rhs.ax.size(),
      context + ": force extent mismatch");
  requireOrThrow(std::isfinite(tolerance) && tolerance >= 0.0, context + ": invalid comparison tolerance");
  for (std::size_t i = 0; i < lhs.ax.size(); ++i) {
    const std::array lhs_components{lhs.ax[i], lhs.ay[i], lhs.az[i]};
    const std::array rhs_components{rhs.ax[i], rhs.ay[i], rhs.az[i]};
    for (std::size_t axis = 0; axis < lhs_components.size(); ++axis) {
      const double scale = std::max({1.0, std::abs(lhs_components[axis]), std::abs(rhs_components[axis])});
      if (!std::isfinite(lhs_components[axis]) || !std::isfinite(rhs_components[axis]) ||
          !(std::abs(lhs_components[axis] - rhs_components[axis]) <= tolerance * scale)) {
        std::ostringstream message;
        message << context << ": force mismatch at target=" << i << ", axis=" << axis
                << ", lhs=" << lhs_components[axis] << ", rhs=" << rhs_components[axis]
                << ", tolerance=" << tolerance * scale;
        throw std::runtime_error(message.str());
      }
    }
  }
}

void testForceComparatorRejectsNonfiniteComponents() {
  const ForceField reference{{1.0}, {2.0}, {3.0}};
  for (const double invalid : {std::numeric_limits<double>::quiet_NaN(),
                              std::numeric_limits<double>::infinity()}) {
    auto candidate = reference;
    candidate.ay[0] = invalid;
    bool rejected = false;
    try { requireEquivalentForces(candidate, reference, 1e-12, "adversarial comparator"); }
    catch (const std::runtime_error& error) {
      rejected = std::string(error.what()).find("adversarial comparator: force mismatch") != std::string::npos;
    }
    requireOrThrow(rejected, "force comparator accepted a nonfinite acceleration");
  }
}

void testPmPlanRetainedStorageMatchesAllocation() {
  int ranks = 1, rank = 0;
#if COSMOSIM_ENABLE_MPI
  MPI_Comm_size(MPI_COMM_WORLD, &ranks);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
#if !COSMOSIM_ENABLE_FFTW
  if (ranks > 1) return; // Distributed FFT storage needs FFTW-MPI.
#endif
  using namespace cosmosim::gravity;
  const PmGridShape shape{10, 7, 9};
  const auto layout = cosmosim::parallel::makePmSlabLayout(10, 7, 9, ranks, rank);
  for (const auto mode : {cosmosim::core::PmDecompositionMode::kSlab,
                          cosmosim::core::PmDecompositionMode::kPencil}) {
    PmGridStorage grid(shape, layout);
    PmSolver solver(shape);
    PmSolveOptions options;
    options.box_size_mpc_comoving = 1.0;
    options.decomposition_mode = mode;
    grid.ensureDensityStorage();
    solver.solvePoissonPeriodicForcesOnly(grid, options);
    const auto estimate = estimatePmPlanResourcesMemory(shape, layout, mode);
    requireOrThrow(estimate.used_backend_allocation_query == (ranks > 1),
        "active distributed PM estimate did not query backend allocation geometry");
    const auto allocated = estimate.allocated_local_complex_cells;
    // In-place layout: the padded real view aliases fourier; a spectral solve
    // owns fourier plus the retained potential_k copy.
    const auto arrays = 2U * allocated * sizeof(std::complex<double>);
    constexpr auto axes = 2U * (10U + 7U + 5U) * sizeof(double);
    cosmosim::core::MemoryReportBuilder builder;
    solver.appendMemoryReport(builder);
    const auto report = std::move(builder).finish();
    bool saw_arrays = false, saw_axes = false;
    for (const auto& entry : report.entries) {
      if (entry.label == "pm_solver.plan_cache_owned_arrays") {
        requireOrThrow(entry.current_size_bytes == arrays, "PM retained FFT arrays disagree with allocation geometry");
        saw_arrays = true;
      } else if (entry.label == "pm_solver.plan_spectral_axis_metadata") {
        requireOrThrow(entry.current_size_bytes == axes, "PM retained global spectral-axis metadata is miscounted");
        saw_axes = true;
      } else continue;
      requireOrThrow(entry.owned_capacity_bytes >= entry.current_size_bytes &&
          entry.high_water_bytes >= entry.owned_capacity_bytes && !entry.estimate_only,
          "PM live report lost retained physical capacity");
    }
    requireOrThrow(saw_arrays && saw_axes && estimate.total_owned_bytes == arrays + axes &&
        estimate.scalar_spectral_array_bytes == 0U, "PM plan owned-memory contract is incomplete");
    solver.shutdownBackendResources();
  }
}

// 20261025 optimization campaign: authored, NOT RUN — deferred to maintainer.
// Finite-difference PM gradient backend on one Fourier mode along each axis.
// For rho = 1 + A cos(k x) the spectral backend is exact and the fourth-order
// stencil multiplies the force by T(kh) = (8 sin(kh) - sin(2kh)) / (6 kh).
// The FD solve must reproduce T(kh) x spectral (slab halos included when
// MPI ranks > 1) and must never allocate the spectral phi_k copy.
void testPmFiniteDifferenceGradientBackend() {
  int ranks = 1, rank = 0;
#if COSMOSIM_ENABLE_MPI
  MPI_Comm_size(MPI_COMM_WORLD, &ranks);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
#if !COSMOSIM_ENABLE_FFTW
  if (ranks > 1) return;
#endif
  using namespace cosmosim::gravity;
  constexpr std::size_t n = 16U;
  if (static_cast<std::size_t>(ranks) * 2U > n) return;  // every rank needs two x planes
  const PmGridShape shape{n, n, n};
  const auto layout = cosmosim::parallel::makePmSlabLayout(n, n, n, ranks, rank);
  constexpr double pi = 3.14159265358979323846;
  for (int axis = 0; axis < 3; ++axis) {
    const auto fill = [&](PmGridStorage& grid) {
      grid.ensureDensityStorage();
      auto density = grid.density();
      for (std::size_t lx = 0; lx < layout.local_nx(); ++lx) {
        const std::size_t gx = layout.globalXFromLocal(lx);
        for (std::size_t iy = 0; iy < n; ++iy) {
          for (std::size_t iz = 0; iz < n; ++iz) {
            const std::size_t coordinate = axis == 0 ? gx : (axis == 1 ? iy : iz);
            density[grid.linearIndex(gx, iy, iz)] =
                1.0 + 0.1 * std::cos(2.0 * pi * static_cast<double>(coordinate) / static_cast<double>(n));
          }
        }
      }
    };
    PmSolveOptions options;
    options.box_size_mpc_comoving = 1.0;
    PmGridStorage spectral_grid(shape, layout);
    PmSolver spectral_solver(shape);
    fill(spectral_grid);
    spectral_solver.solvePoissonPeriodicForcesOnly(spectral_grid, options);
    options.force_gradient_backend = PmForceGradientBackend::kFiniteDifference4;
    PmGridStorage fd_grid(shape, layout);
    PmSolver fd_solver(shape);
    fill(fd_grid);
    fd_solver.solvePoissonPeriodicForcesOnly(fd_grid, options);
    const double kh = 2.0 * pi / static_cast<double>(n);
    const double transfer = (8.0 * std::sin(kh) - std::sin(2.0 * kh)) / (6.0 * kh);
    const std::array<std::span<const double>, 3> spectral{
        std::as_const(spectral_grid).force_x(), std::as_const(spectral_grid).force_y(),
        std::as_const(spectral_grid).force_z()};
    const std::array<std::span<const double>, 3> fd{
        std::as_const(fd_grid).force_x(), std::as_const(fd_grid).force_y(), std::as_const(fd_grid).force_z()};
    double scale = 0.0;
    for (const double value : spectral[static_cast<std::size_t>(axis)]) scale = std::max(scale, std::abs(value));
    requireOrThrow(scale > 0.0, "single-mode spectral PM force vanished");
    for (std::size_t component = 0; component < 3U; ++component) {
      for (std::size_t i = 0; i < spectral[component].size(); ++i) {
        const double expected = transfer * spectral[component][i];
        requireOrThrow(std::abs(fd[component][i] - expected) <= 1.0e-9 * scale,
            "finite-difference PM force differs from T(kh) x spectral force");
      }
    }
    cosmosim::core::MemoryReportBuilder builder;
    fd_solver.appendMemoryReport(builder);
    const auto report = std::move(builder).finish();
    const auto estimate = estimatePmPlanResourcesMemory(
        shape, layout, cosmosim::core::PmDecompositionMode::kSlab, PmForceGradientBackend::kFiniteDifference4);
    for (const auto& entry : report.entries) {
      if (entry.label != "pm_solver.plan_cache_owned_arrays") continue;
      requireOrThrow(entry.current_size_bytes ==
          estimate.real_array_bytes + estimate.complex_spectral_array_bytes + estimate.finite_difference_halo_bytes,
          "finite-difference PM plan retained a spectral phi_k copy or misreported halo planes");
    }
    fd_solver.shutdownBackendResources();
    spectral_solver.shutdownBackendResources();
  }
}

void testShortRangeSurfaceRejectsPmRequests() {
#if COSMOSIM_ENABLE_MPI
  int world_size = 1; MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  if (world_size != 1) return;
#endif
  using namespace cosmosim::gravity;
  const std::vector<double> x{0.25, 0.30}, mass(2, 1.0);
  const std::vector<std::uint32_t> active{0U, 1U};
  ForceField force{std::vector<double>(2), std::vector<double>(2), std::vector<double>(2)};
  TreePmForceAccumulatorView accumulator{active, force.ax, force.ay, force.az};
  TreePmOptions options;
  options.pm_options.box_size_mpc_comoving = 1.0;
  options.short_range_only = true;
  options.tree_options.softening.epsilon_comoving = 0.01;
  options.split_policy = makeTreePmSplitPolicyFromMeshSpacing(1.25, 3.0, 0.125);
  TreePmCoordinator coordinator({8, 8, 8});
  TreePmDiagnostics diagnostics;
  solveForEvaluation(coordinator, x, x, x, mass, accumulator, options,
      ForceEvaluation::kShortRangeOnly, nullptr, &diagnostics);
  requireOrThrow(diagnostics.pm_solve_count == 0U && diagnostics.pm_reuse_count == 0U,
      "valid short-range request touched PM");
  requireOrThrow(std::isfinite(force.ax[0]) && force.ax[0] > 0.0,
      "valid short-range request never evaluated the pair");
  const auto rejected = [&](const auto& action) {
    bool refused = false;
    try { action(); }
    catch (const std::invalid_argument& error) {
      const std::string message = error.what();
      refused = message.find("short-range split surface") != std::string::npos &&
          message.find("PM refresh/capture or zoom correction") != std::string::npos;
    }
    requireOrThrow(refused, "contradictory short-range/PM request was not rejected by its guard");
  };
  rejected([&] { coordinator.solveActiveSet(x, x, x, mass, accumulator, options); });
  rejected([&] { coordinator.solveActiveSetWithPmCadence(
      x, x, x, mass, accumulator, options, true); });
  const std::vector<std::uint8_t> high_res(2, 1U);
  options.zoom_focused_pm_shape = {8, 8, 8};
  options.source_is_high_res = high_res;
  options.active_is_high_res = high_res;
  options.zoom_region_radius_comoving = 0.1;
  options.enable_zoom_long_range_correction = true;
  rejected([&] { coordinator.solveActiveSetWithPmCadence(
      x, x, x, mass, accumulator, options, false); });
  options.enable_zoom_long_range_correction = false;
  ForceField captured{std::vector<double>(2), std::vector<double>(2), std::vector<double>(2)};
  accumulator.long_range_accel_x_comoving = captured.ax;
  accumulator.long_range_accel_y_comoving = captured.ay;
  accumulator.long_range_accel_z_comoving = captured.az;
  rejected([&] { coordinator.solveActiveSetWithPmCadence(
      x, x, x, mass, accumulator, options, false); });
}

void testCoincidentPairVectorIsZeroWithFiniteSoftening() {
#if COSMOSIM_ENABLE_MPI
  int world_size = 1; MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  if (world_size != 1) return; // The routing fixtures below cover distributed targets.
#endif
  using namespace cosmosim::gravity;
  const std::vector<double> x(2, 0.25), mass(2, 1.0);
  TreePmOptions options;
  options.pm_options.box_size_mpc_comoving = 1.0;
  options.pm_options.gravitational_constant_code = 1.0;
  options.tree_options.gravitational_constant_code = 1.0;
  options.tree_options.max_leaf_size = 8;
  options.split_policy = makeTreePmSplitPolicyFromMeshSpacing(1.25, 3.0, 0.125);
  options.short_range_only = true;
  for (const double epsilon : {0.0, 0.01, 1e-110}) {
    options.tree_options.softening.epsilon_comoving = epsilon;
    for (const bool lookup : {false, true}) {
      options.gaussian_pair_lookup_enabled = lookup;
      TreePmDiagnostics diagnostics;
      TreePmProfileEvent profile;
      const auto force = solveTreePm(x, x, x, mass, {8, 8, 8}, options,
          &diagnostics, &profile, {}, ForceEvaluation::kShortRangeOnly);
      requireOrThrow(diagnostics.pm_solve_count == 0U && diagnostics.pm_reuse_count == 0U,
          "short-range fixture refreshed or reused PM");
      requireOrThrow(profile.tree_profile.particle_particle_interactions == 2U,
          "coincident pairs must remain counted direct interactions");
      for (std::size_t i = 0; i < 2U; ++i) {
        requireOrThrow(force.ax[i] == 0.0 && force.ay[i] == 0.0 && force.az[i] == 0.0,
            "coincident short-range force must be the exact zero vector for both kernels");
      }
    }
  }
}

void testRecoveryDirectPairsCutoffThreadsAndMemory() {
#if COSMOSIM_ENABLE_MPI
  int world_size=1; MPI_Comm_size(MPI_COMM_WORLD,&world_size);
  if (world_size!=1) return;  // Separate distributed fixtures below own MPI routing.
#endif
  using namespace cosmosim::gravity;
  constexpr std::size_t n=96U;
  std::vector<double> x(n),y(n),z(n),mass(n,1.0),epsilon(n);
  std::vector<std::uint8_t> overrides(n,1U);
  std::vector<std::uint32_t> active(n);
  for (std::size_t i=0;i<n;++i) {
    active[i]=static_cast<std::uint32_t>(i);
    x[i]=std::fmod(0.173*i,2.0);y[i]=std::fmod(0.257*i,3.0);z[i]=std::fmod(0.317*i,4.0);
    epsilon[i]=0.002+0.0001*(i%7U);
    // Exercise general multiple-image wrapping and negative coordinates.
    if (i%2U) {x[i]-=8.0;y[i]+=9.0;z[i]-=20.0;}
  }
  TreePmOptions options;
  options.pm_options.box_size_x_mpc_comoving=2.0;options.pm_options.box_size_y_mpc_comoving=3.0;
  options.pm_options.box_size_z_mpc_comoving=4.0;
  options.pm_options.gravitational_constant_code=1.0;
  options.tree_options.gravitational_constant_code=1.0;
  options.tree_options.max_leaf_size=1;
  options.tree_options.multipole_order=TreeMultipoleOrder::kMonopole; // strict descends to pairs
  options.split_policy=makeTreePmSplitPolicyFromMeshSpacing(1.25,3.0,0.125);
  options.short_range_only=true;
  const double cutoff=options.split_policy.cutoff_radius_comoving;
  x[0]=0.0;y[0]=0.0;z[0]=0.0;
  x[1]=cutoff;y[1]=0.0;z[1]=0.0;
  x[2]=std::nextafter(cutoff,1.0);y[2]=0.0;z[2]=0.0;
  x[3]=std::nextafter(cutoff,0.0);y[3]=0.0;z[3]=0.0;
  x[4]=1.0;y[4]=0.0;z[4]=0.0; // exact half-box tie
  ForceField direct{std::vector<double>(n),std::vector<double>(n),std::vector<double>(n)};
  std::uint64_t expected_pairs=0;
  for (std::size_t i=0;i<n;++i) for (std::size_t j=0;j<n;++j) {
    if (i==j) continue;
    const double dx=minimumImageDelta(x[j]-x[i],2.0),dy=minimumImageDelta(y[j]-y[i],3.0),dz=minimumImageDelta(z[j]-z[i],4.0);
    const double r2=dx*dx+dy*dy+dz*dz;
    if (r2>cutoff*cutoff) continue;
    ++expected_pairs;
    const double f=treePmSoftenedShortRangeInvR3(r2,std::max(epsilon[i],epsilon[j]),options.split_policy.split_scale_comoving);
    direct.ax[i]+=f*dx;direct.ay[i]+=f*dy;direct.az[i]+=f*dz;
  }
  cosmosim::core::MemoryGovernor governor;
  {
    const PmGridShape shape{8,8,8};
    const auto layout=cosmosim::parallel::makePmSlabLayout(8,8,8,1,0);
    TreePmCoordinator coordinator(shape,layout,cosmosim::parallel::MpiContext(false,1,0),&governor);
    ForceField field{std::vector<double>(n),std::vector<double>(n),std::vector<double>(n)};
    const TreeSofteningView soften{.source_particle_epsilon_comoving=epsilon,
        .source_particle_epsilon_override_mask=overrides};
    TreePmForceAccumulatorView acc{active,field.ax,field.ay,field.az};
    TreePmDiagnostics analytic,lookup;
    solveForEvaluation(coordinator,x,y,z,mass,acc,options,ForceEvaluation::kShortRangeOnly,nullptr,&analytic,soften);
    requireOrThrow(analytic.local_traversal.direct_pair_evaluations==expected_pairs,"cutoff/identity changed pair set");
    requireOrThrow(relativeL2Error(field,direct)<1e-13,"reverse-octant traversal changed physical interactions");
    const auto baseline=field;
    options.gaussian_pair_lookup_enabled=true;options.full_mac_diagnostics=false;
    solveForEvaluation(coordinator,x,y,z,mass,acc,options,ForceEvaluation::kShortRangeOnly,nullptr,&lookup,soften);
    requireOrThrow(lookup.local_traversal.direct_pair_evaluations==expected_pairs,"lookup changed pair set");
    requireOrThrow(lookup.local_traversal.visited_nodes==analytic.local_traversal.visited_nodes,"fast MAC changed traversal decisions");
    requireOrThrow(relativeL2Error(field,baseline)<5e-9,"lookup exceeded force regression budget");
    bool accounted=false;
    for (const auto& entry:coordinator.memoryReport().entries) if (entry.label=="treepm.gaussian_pair_table") {
      requireOrThrow(entry.owned_capacity_bytes==sizeof(TreePmGaussianCoefficientTable) && entry.governed_commitment,
          "lookup retained storage is not governed physical capacity");accounted=true;
    }
    requireOrThrow(accounted,"lookup table absent from memory report");
#if defined(_OPENMP)
    const int old_threads=omp_get_max_threads();
    omp_set_num_threads(1);
    solveForEvaluation(coordinator,x,y,z,mass,acc,options,ForceEvaluation::kShortRangeOnly,nullptr,&lookup,soften);
    const auto one_thread=field;
    omp_set_num_threads(4);
    solveForEvaluation(coordinator,x,y,z,mass,acc,options,ForceEvaluation::kShortRangeOnly,nullptr,&lookup,soften);
    requireOrThrow(field.ax==one_thread.ax && field.ay==one_thread.ay && field.az==one_thread.az,
        "target force accumulation changed with worker count");
    omp_set_num_threads(old_threads);
#endif
    const auto retained=governor.snapshot().committed_bytes;
    solveForEvaluation(coordinator,x,y,z,mass,acc,options,ForceEvaluation::kShortRangeOnly,nullptr,nullptr,soften);
    requireOrThrow(governor.snapshot().committed_bytes==retained,"lookup table retained another allocation across epochs");
    coordinator.shutdownMpiResources();
  }
  requireOrThrow(governor.snapshot().committed_bytes==0U,"lookup commitment outlived coordinator");
}

void testPeriodicTreePmAgainstDirectReference() {
  constexpr std::size_t particle_count = 48;
  constexpr double box_size_comoving = 1.0;

  std::vector<double> pos_x(particle_count, 0.0);
  std::vector<double> pos_y(particle_count, 0.0);
  std::vector<double> pos_z(particle_count, 0.0);
  std::vector<double> mass(particle_count, 1.0);

  for (std::size_t i = 0; i < particle_count; ++i) {
    pos_x[i] = std::fmod((17.0 * static_cast<double>(i) + 3.0) * 0.017, box_size_comoving);
    pos_y[i] = std::fmod((29.0 * static_cast<double>(i) + 5.0) * 0.013, box_size_comoving);
    pos_z[i] = std::fmod((41.0 * static_cast<double>(i) + 7.0) * 0.011, box_size_comoving);
    mass[i] = 0.8 + 0.01 * static_cast<double>(i % 11U);
  }

  cosmosim::gravity::TreePmOptions options;
  options.pm_options.box_size_mpc_comoving = box_size_comoving;
  options.pm_options.scale_factor = 1.0;
  options.pm_options.gravitational_constant_code = 1.0;
  options.pm_options.enable_window_deconvolution = true;
  options.tree_options.opening_theta = 0.55;
  options.tree_options.max_leaf_size = 8;
  options.tree_options.gravitational_constant_code = 1.0;
  options.tree_options.softening.epsilon_comoving = 0.01;

  const cosmosim::gravity::PmGridShape pm_shape = pmShapeForAvailableBackend(32, 8);
  const double mesh_spacing = box_size_comoving / static_cast<double>(pm_shape.nx);
  options.split_policy = cosmosim::gravity::makeTreePmSplitPolicyFromMeshSpacing(
      1.25, rcutCellsForAvailableBackend(4.5, 3.9), mesh_spacing);

  cosmosim::gravity::TreePmDiagnostics diagnostics;
  const ForceField tree_pm_force = solveTreePm(pos_x, pos_y, pos_z, mass, pm_shape, options, &diagnostics);

  const ForceField reference = computeDirectReference(
      pos_x,
      pos_y,
      pos_z,
      mass,
      box_size_comoving,
      options.tree_options.softening,
      options.tree_options.gravitational_constant_code);

  const double rel_l2 = relativeL2Error(tree_pm_force, reference);
#if COSMOSIM_ENABLE_FFTW
  const double max_rel_l2 = 0.75;
#else
  const double max_rel_l2 = 1.8;
#endif

  std::ostringstream diag;
  diag << "TreePM periodic validation failed (direct reference is minimum-image periodic, not Ewald): build_flag.COSMOSIM_ENABLE_FFTW="
       << (COSMOSIM_ENABLE_FFTW ? "ON" : "OFF")
       << ", treePmSupportedByBuild()=" << (cosmosim::gravity::treePmSupportedByBuild() ? "true" : "false")
       << ", rel_l2=" << rel_l2 << " (required <= " << max_rel_l2 << ')'
       << ", pm_grid=" << pm_shape.nx
       << ", mesh_spacing=" << diagnostics.mesh_spacing_comoving
       << ", asmth_cells=" << diagnostics.asmth_cells
       << ", rcut_cells=" << diagnostics.rcut_cells
       << ", split_scale_comoving=" << diagnostics.split_scale_comoving
       << ", cutoff_radius_comoving=" << diagnostics.cutoff_radius_comoving
       << ", residual_pruned_nodes=" << diagnostics.residual_pruned_nodes
       << ", residual_pair_skips_cutoff=" << diagnostics.residual_pair_skips_cutoff
       << ", composition_error_at_split=" << diagnostics.composition_error_at_split
       << ", max_relative_composition_error=" << diagnostics.max_relative_composition_error;

  requireOrThrow(rel_l2 < max_rel_l2, diag.str());
  requireOrThrow(diagnostics.composition_error_at_split < 1.0e-12, diag.str());
  requireOrThrow(diagnostics.max_relative_composition_error < 1.0e-12, diag.str());
  requireOrThrow(
      diagnostics.residual_pruned_nodes > 0 || diagnostics.residual_pair_skips_cutoff > 0,
      diag.str());
  requireOrThrow(diagnostics.residual_pair_skips_cutoff > 0, diag.str());
}

void testPeriodicTreeGeometryIsSeamSafeAndTranslationInvariant() {
  constexpr std::uint8_t x_seam = 1U << 0U;
  constexpr std::uint8_t y_seam = 1U << 1U;
  constexpr std::uint8_t z_seam = 1U << 2U;
  const AxisBoxLengths box{.x = 1.0, .y = 1.5, .z = 2.0};
  const cosmosim::gravity::PmGridShape pm_shape{8U, 8U, 8U};

  enum class TargetPlacement {
    kAllClusterParticles,
    kInteriorTarget,
    kSeamTargetWithInteriorSources,
    kSeamTargetWithSeamSources,
  };
  struct Case {
    std::string label;
    std::uint8_t source_seams = 0U;
    TargetPlacement target_placement = TargetPlacement::kAllClusterParticles;
    cosmosim::gravity::TreeMultipoleOrder multipole_order = cosmosim::gravity::TreeMultipoleOrder::kMonopole;
    cosmosim::gravity::TreeOpeningCriterion opening_criterion =
        cosmosim::gravity::TreeOpeningCriterion::kBarnesHutGeometric;
    bool direct_leaf = false;
    bool require_compact_root = false;
    bool require_internal_acceptance = false;
  };

  std::vector<Case> cases;
  for (const auto order : {cosmosim::gravity::TreeMultipoleOrder::kMonopole,
                           cosmosim::gravity::TreeMultipoleOrder::kQuadrupole}) {
    for (const auto criterion : {cosmosim::gravity::TreeOpeningCriterion::kBarnesHutGeometric,
                                 cosmosim::gravity::TreeOpeningCriterion::kBarnesHutComDistance}) {
      for (const bool direct_leaf : {true, false}) {
        cases.push_back(Case{
            .label = "x_seam_full_numerics_matrix",
            .source_seams = x_seam,
            .target_placement = TargetPlacement::kAllClusterParticles,
            .multipole_order = order,
            .opening_criterion = criterion,
            .direct_leaf = direct_leaf,
            .require_compact_root = true,
        });
      }
    }
  }
  cases.push_back(Case{.label = "y_seam", .source_seams = y_seam, .require_compact_root = true});
  cases.push_back(Case{
      .label = "z_seam_quadrupole",
      .source_seams = z_seam,
      .multipole_order = cosmosim::gravity::TreeMultipoleOrder::kQuadrupole,
      .require_compact_root = true});
  cases.push_back(Case{
      .label = "xy_edge_com_mac",
      .source_seams = static_cast<std::uint8_t>(x_seam | y_seam),
      .opening_criterion = cosmosim::gravity::TreeOpeningCriterion::kBarnesHutComDistance,
      .require_compact_root = true});
  cases.push_back(Case{
      .label = "xyz_corner_quadrupole",
      .source_seams = static_cast<std::uint8_t>(x_seam | y_seam | z_seam),
      .multipole_order = cosmosim::gravity::TreeMultipoleOrder::kQuadrupole,
      .opening_criterion = cosmosim::gravity::TreeOpeningCriterion::kBarnesHutComDistance,
      .require_compact_root = true});
  cases.push_back(Case{
      .label = "interior_target_seam_sources",
      .source_seams = static_cast<std::uint8_t>(x_seam | y_seam),
      .target_placement = TargetPlacement::kInteriorTarget,
      .multipole_order = cosmosim::gravity::TreeMultipoleOrder::kQuadrupole,
      .opening_criterion = cosmosim::gravity::TreeOpeningCriterion::kBarnesHutComDistance,
      .require_internal_acceptance = true});
  cases.push_back(Case{
      .label = "seam_target_interior_sources",
      .source_seams = 0U,
      .target_placement = TargetPlacement::kSeamTargetWithInteriorSources,
      .multipole_order = cosmosim::gravity::TreeMultipoleOrder::kQuadrupole});
  cases.push_back(Case{
      .label = "seam_target_seam_sources",
      .source_seams = static_cast<std::uint8_t>(x_seam | y_seam | z_seam),
      .target_placement = TargetPlacement::kSeamTargetWithSeamSources,
      .opening_criterion = cosmosim::gravity::TreeOpeningCriterion::kBarnesHutComDistance});

  for (const Case& test_case : cases) {
    constexpr std::size_t source_count = 12U;
    std::vector<double> pos_x(source_count, 0.0);
    std::vector<double> pos_y(source_count, 0.0);
    std::vector<double> pos_z(source_count, 0.0);
    std::vector<double> mass(source_count, 1.0);
    const auto coordinate = [](std::size_t i, double length, bool crosses_seam) {
      const double offset = (0.006 + 0.0015 * static_cast<double>(i % 4U)) * length;
      if (crosses_seam) {
        return (i % 2U) == 0U ? offset : length - offset;
      }
      return (0.5 + 0.004 * static_cast<double>(static_cast<int>(i % 5U) - 2)) * length;
    };
    for (std::size_t i = 0; i < source_count; ++i) {
      pos_x[i] = coordinate(i, box.x, (test_case.source_seams & x_seam) != 0U);
      pos_y[i] = coordinate(i + 1U, box.y, (test_case.source_seams & y_seam) != 0U);
      pos_z[i] = coordinate(i + 2U, box.z, (test_case.source_seams & z_seam) != 0U);
      mass[i] = 0.7 + 0.07 * static_cast<double>(i % 6U);
    }

    std::vector<std::uint32_t> active;
    if (test_case.target_placement == TargetPlacement::kAllClusterParticles) {
      active.resize(source_count);
      for (std::size_t i = 0; i < source_count; ++i) {
        active[i] = static_cast<std::uint32_t>(i);
      }
    } else {
      const bool target_on_seam =
          test_case.target_placement == TargetPlacement::kSeamTargetWithInteriorSources ||
          test_case.target_placement == TargetPlacement::kSeamTargetWithSeamSources;
      pos_x.push_back(target_on_seam ? 0.998 * box.x : 0.47 * box.x);
      pos_y.push_back(target_on_seam ? 0.002 * box.y : 0.44 * box.y);
      pos_z.push_back(target_on_seam ? 0.999 * box.z : 0.53 * box.z);
      mass.push_back(0.0);
      active.push_back(static_cast<std::uint32_t>(source_count));
    }

    std::vector<double> lifted_x = pos_x;
    std::vector<double> lifted_y = pos_y;
    std::vector<double> lifted_z = pos_z;
    for (std::size_t i = 0; i < pos_x.size(); ++i) {
      lifted_x[i] += static_cast<double>(static_cast<int>(i % 3U) - 1) * box.x;
      lifted_y[i] += static_cast<double>(static_cast<int>((i + 1U) % 3U) - 1) * box.y;
      lifted_z[i] += static_cast<double>(static_cast<int>((i + 2U) % 3U) - 1) * box.z;
    }

    cosmosim::gravity::TreePmOptions options;
    options.pm_options.box_size_x_mpc_comoving = box.x;
    options.pm_options.box_size_y_mpc_comoving = box.y;
    options.pm_options.box_size_z_mpc_comoving = box.z;
    options.pm_options.scale_factor = 1.0;
    options.pm_options.gravitational_constant_code = 1.0;
    options.tree_options.opening_theta = 0.8;
    options.tree_options.opening_criterion = test_case.opening_criterion;
    options.tree_options.multipole_order = test_case.multipole_order;
    options.tree_options.max_leaf_size = test_case.direct_leaf ? 64U : 1U;
    options.tree_options.gravitational_constant_code = 1.0;
    options.tree_options.softening.epsilon_comoving = 2.0e-3;
    options.split_policy = cosmosim::gravity::makeTreePmSplitPolicyFromMeshSpacing(
        1.25, 3.9, box.x / static_cast<double>(pm_shape.nx));

    cosmosim::gravity::TreePmDiagnostics wrapped_diagnostics;
    cosmosim::gravity::TreePmDiagnostics lifted_diagnostics;
    cosmosim::gravity::TreePmProfileEvent wrapped_profile;
    const ForceField wrapped_force = solveTreePmActiveSubset(
        pos_x, pos_y, pos_z, mass, active, pm_shape, options, &wrapped_diagnostics, &wrapped_profile);
    const ForceField lifted_force = solveTreePmActiveSubset(
        lifted_x, lifted_y, lifted_z, mass, active, pm_shape, options, &lifted_diagnostics);
    requireEquivalentForces(wrapped_force, lifted_force, 5.0e-11, test_case.label);
    auto accelerated_options = options;
    accelerated_options.radial_coefficients_enabled = true;
    accelerated_options.full_mac_diagnostics = false;
    cosmosim::gravity::TreePmDiagnostics accelerated_diagnostics;
    const auto accelerated_force = solveTreePmActiveSubset(pos_x, pos_y, pos_z, mass, active,
        pm_shape, accelerated_options, &accelerated_diagnostics);
    requireEquivalentForces(accelerated_force, wrapped_force, 5.0e-11,
        test_case.label + ": shared screened derivatives and fast descent");
    requireOrThrow(accelerated_diagnostics.local_traversal.accepted_internal_multipoles ==
            wrapped_diagnostics.local_traversal.accepted_internal_multipoles &&
        accelerated_diagnostics.local_traversal.direct_pair_evaluations ==
            wrapped_diagnostics.local_traversal.direct_pair_evaluations,
        test_case.label + ": fast diagnostics changed accepted work");

    const double geometry_scale = std::max({1.0,
                                            wrapped_diagnostics.tree_root_half_size_comoving,
                                            lifted_diagnostics.tree_root_half_size_comoving});
    requireOrThrow(
        std::abs(wrapped_diagnostics.tree_root_half_size_comoving -
                 lifted_diagnostics.tree_root_half_size_comoving) <= 5.0e-12 * geometry_scale,
        test_case.label + ": periodic translation changed root extent");
    if (test_case.require_compact_root) {
      requireOrThrow(
          wrapped_diagnostics.tree_root_half_size_comoving < 0.03 * box.z,
          test_case.label + ": seam cluster produced a noncompact tree root");
    }
    if (test_case.require_internal_acceptance) {
      requireOrThrow(
          wrapped_profile.tree_profile.particle_particle_interactions < source_count,
          test_case.label + ": distant compact seam cluster was not accepted as an internal node");
    }
  }
}

void testRelativeForceErrorMacUsesOwnerHistoryAcrossPeriodicSeam() {
  const AxisBoxLengths box{.x = 1.0, .y = 1.4, .z = 1.8};
  const cosmosim::gravity::PmGridShape pm_shape{8U, 8U, 8U};
  std::vector<double> pos_x{0.006, 0.994, 0.009, 0.991, 0.012, 0.988, 0.015, 0.985, 0.35};
  std::vector<double> pos_y{0.69, 0.71, 0.68, 0.72, 0.695, 0.705, 0.685, 0.715, 0.70};
  std::vector<double> pos_z{0.89, 0.91, 0.88, 0.92, 0.895, 0.905, 0.885, 0.915, 0.90};
  std::vector<double> mass(pos_x.size(), 1.0);
  mass.back() = 0.0;
  const std::vector<std::uint32_t> active{static_cast<std::uint32_t>(pos_x.size() - 1U)};

  cosmosim::gravity::TreePmOptions relative_options;
  relative_options.pm_options.box_size_x_mpc_comoving = box.x;
  relative_options.pm_options.box_size_y_mpc_comoving = box.y;
  relative_options.pm_options.box_size_z_mpc_comoving = box.z;
  relative_options.pm_options.gravitational_constant_code = 1.0;
  relative_options.tree_options.opening_criterion = cosmosim::gravity::TreeOpeningCriterion::kRelativeForceError;
  relative_options.tree_options.opening_theta = 0.03;
  relative_options.tree_options.multipole_order = cosmosim::gravity::TreeMultipoleOrder::kQuadrupole;
  relative_options.tree_options.max_leaf_size = 1U;
  relative_options.tree_options.relative_force_tolerance = 0.005;
  relative_options.tree_options.relative_force_acceleration_floor_code = 1.0e-12;
  relative_options.tree_options.softening.epsilon_comoving = 1.0e-3;
  relative_options.split_policy = cosmosim::gravity::makeTreePmSplitPolicyFromMeshSpacing(
      1.25, 3.9, box.x / static_cast<double>(pm_shape.nx));

  cosmosim::gravity::TreePmOptions com_options = relative_options;
  com_options.tree_options.opening_criterion = cosmosim::gravity::TreeOpeningCriterion::kBarnesHutComDistance;
  cosmosim::gravity::TreePmProfileEvent missing_profile;
  cosmosim::gravity::TreePmProfileEvent com_profile;
  const ForceField missing_history = solveTreePmActiveSubset(
      pos_x, pos_y, pos_z, mass, active, pm_shape, relative_options, nullptr, &missing_profile);
  const ForceField com_fallback = solveTreePmActiveSubset(
      pos_x, pos_y, pos_z, mass, active, pm_shape, com_options, nullptr, &com_profile);
  requireEquivalentForces(missing_history, com_fallback, 1.0e-13, "relative MAC missing-history fallback");
  requireOrThrow(
      missing_profile.tree_profile.visited_nodes == com_profile.tree_profile.visited_nodes,
      "relative MAC missing history must reproduce COM-distance tree work");

  const std::vector<double> previous_acceleration{1.0e8};
  cosmosim::gravity::TreePmProfileEvent history_profile;
  const ForceField with_history = solveTreePmActiveSubset(
      pos_x,
      pos_y,
      pos_z,
      mass,
      active,
      pm_shape,
      relative_options,
      nullptr,
      &history_profile,
      previous_acceleration);
  requireOrThrow(
      history_profile.tree_profile.particle_particle_interactions <
          missing_profile.tree_profile.particle_particle_interactions,
      "relative MAC owner history should permit less tree work for a distant compact node");

  std::vector<double> lifted_x = pos_x;
  std::vector<double> lifted_y = pos_y;
  std::vector<double> lifted_z = pos_z;
  for (std::size_t i = 0; i < pos_x.size(); ++i) {
    lifted_x[i] += (i % 2U == 0U ? box.x : -box.x);
    lifted_y[i] += (i % 3U == 0U ? box.y : -box.y);
    lifted_z[i] += (i % 2U == 0U ? -box.z : box.z);
  }
  const ForceField lifted_history = solveTreePmActiveSubset(
      lifted_x,
      lifted_y,
      lifted_z,
      mass,
      active,
      pm_shape,
      relative_options,
      nullptr,
      nullptr,
      previous_acceleration);
  requireEquivalentForces(with_history, lifted_history, 5.0e-11, "relative MAC periodic translation");
}

void testIndependentTargetAndEmptySourceContractSingleRank() {
  const cosmosim::gravity::PmGridShape pm_shape{8U, 8U, 8U};
  cosmosim::gravity::TreePmOptions options;
  options.pm_options.box_size_mpc_comoving = 1.0;
  options.pm_options.gravitational_constant_code = 1.0;
  options.tree_options.softening.epsilon_comoving = 1.0e-3;
  options.split_policy = cosmosim::gravity::makeTreePmSplitPolicyFromMeshSpacing(1.25, 3.9, 0.125);

  const std::vector<double> empty;
  const std::vector<std::uint32_t> target_identity{std::numeric_limits<std::uint32_t>::max()};
  const std::vector<double> target_x{0.01};
  const std::vector<double> target_y{0.25};
  const std::vector<double> target_z{0.75};
  ForceField empty_source_force{{0.0}, {0.0}, {0.0}};
  const cosmosim::gravity::TreePmForceAccumulatorView empty_source_target{
      .active_particle_index = target_identity,
      .accel_x_comoving = empty_source_force.ax,
      .accel_y_comoving = empty_source_force.ay,
      .accel_z_comoving = empty_source_force.az,
      .target_pos_x_comoving = target_x,
      .target_pos_y_comoving = target_y,
      .target_pos_z_comoving = target_z,
  };
  cosmosim::gravity::TreePmCoordinator empty_source_coordinator(pm_shape);
  empty_source_coordinator.solveActiveSet(
      empty, empty, empty, empty, empty_source_target, options, nullptr, nullptr);
  requireOrThrow(
      empty_source_force.ax[0] == 0.0 && empty_source_force.ay[0] == 0.0 && empty_source_force.az[0] == 0.0,
      "an independent target with no sources must receive exactly zero force");

  const std::vector<std::uint32_t> no_targets;
  ForceField all_empty_force;
  const cosmosim::gravity::TreePmForceAccumulatorView all_empty{
      .active_particle_index = no_targets,
      .accel_x_comoving = all_empty_force.ax,
      .accel_y_comoving = all_empty_force.ay,
      .accel_z_comoving = all_empty_force.az,
  };
  cosmosim::gravity::TreePmCoordinator all_empty_coordinator(pm_shape);
  all_empty_coordinator.solveActiveSet(empty, empty, empty, empty, all_empty, options, nullptr, nullptr);
}

void testResidualCutoffPrunesTreePath() {
  constexpr double box_size_comoving = 1.0;
  std::vector<double> pos_x = {0.10, 0.18, 0.62, 0.88};
  std::vector<double> pos_y = {0.10, 0.12, 0.60, 0.90};
  std::vector<double> pos_z = {0.10, 0.09, 0.58, 0.86};
  std::vector<double> mass = {1.0, 1.0, 1.0, 1.0};

  const cosmosim::gravity::PmGridShape pm_shape = pmShapeForAvailableBackend(32, 8);
  const double mesh_spacing = box_size_comoving / static_cast<double>(pm_shape.nx);

  cosmosim::gravity::TreePmOptions options;
  options.pm_options.box_size_mpc_comoving = box_size_comoving;
  options.pm_options.scale_factor = 1.0;
  options.pm_options.gravitational_constant_code = 1.0;
  options.tree_options.opening_theta = 0.55;
  options.tree_options.max_leaf_size = 2;
  options.tree_options.gravitational_constant_code = 1.0;
  options.tree_options.softening.epsilon_comoving = 1.0e-3;
  options.split_policy = cosmosim::gravity::makeTreePmSplitPolicyFromMeshSpacing(1.0e-3, 3.0, mesh_spacing);

  cosmosim::gravity::TreePmDiagnostics diagnostics;
  [[maybe_unused]] const ForceField field = solveTreePm(pos_x, pos_y, pos_z, mass, pm_shape, options, &diagnostics);

  std::ostringstream msg;
  msg << "Residual cutoff pruning counters invalid: pruned_nodes=" << diagnostics.residual_pruned_nodes
      << ", pair_skips=" << diagnostics.residual_pair_skips_cutoff
      << ", pair_evaluations=" << diagnostics.residual_pair_evaluations;

  requireOrThrow(diagnostics.residual_pruned_nodes > 0, msg.str());
  requireOrThrow(diagnostics.residual_pair_skips_cutoff > 0, msg.str());
  requireOrThrow(diagnostics.residual_pair_evaluations > 0, msg.str());
}

void testPmOnlyAndTreePmConsistency() {
  constexpr std::size_t particle_count = 32;
  constexpr double box_size_comoving = 1.0;

  std::vector<double> pos_x(particle_count, 0.0);
  std::vector<double> pos_y(particle_count, 0.0);
  std::vector<double> pos_z(particle_count, 0.0);
  std::vector<double> mass(particle_count, 1.0);

  for (std::size_t i = 0; i < particle_count; ++i) {
    pos_x[i] = std::fmod((13.0 * static_cast<double>(i) + 1.0) * 0.023, box_size_comoving);
    pos_y[i] = std::fmod((19.0 * static_cast<double>(i) + 2.0) * 0.019, box_size_comoving);
    pos_z[i] = std::fmod((23.0 * static_cast<double>(i) + 3.0) * 0.017, box_size_comoving);
    mass[i] = 0.9 + 0.05 * static_cast<double>(i % 5U);
  }
  // Guarantee one source inside r_cut on every backend: the FFTW 32^3 mesh
  // gives r_cut = 4.5/32 = 0.14, below every separation of the quasi-random
  // lattice above, so the residual path otherwise evaluates no direct pair.
  pos_x[1] = std::fmod(pos_x[0] + 0.03, box_size_comoving);
  pos_y[1] = std::fmod(pos_y[0] + 0.02, box_size_comoving);
  pos_z[1] = std::fmod(pos_z[0] + 0.01, box_size_comoving);

  const cosmosim::gravity::PmGridShape pm_shape = pmShapeForAvailableBackend(32, 8);
  const double mesh_spacing = box_size_comoving / static_cast<double>(pm_shape.nx);

  cosmosim::gravity::TreePmOptions base_options;
  base_options.pm_options.box_size_mpc_comoving = box_size_comoving;
  base_options.pm_options.scale_factor = 1.0;
  base_options.pm_options.gravitational_constant_code = 1.0;
  base_options.pm_options.enable_window_deconvolution = true;
  base_options.tree_options.opening_theta = 0.55;
  base_options.tree_options.max_leaf_size = 8;
  base_options.tree_options.gravitational_constant_code = 1.0;
  base_options.tree_options.softening.epsilon_comoving = 0.01;

  const ForceField direct = computeDirectReference(
      pos_x,
      pos_y,
      pos_z,
      mass,
      box_size_comoving,
      base_options.tree_options.softening,
      base_options.tree_options.gravitational_constant_code);

  cosmosim::gravity::TreePmOptions pm_only_options = base_options;
  pm_only_options.split_policy = cosmosim::gravity::makeTreePmSplitPolicyFromMeshSpacing(
      0.2, rcutCellsForAvailableBackend(4.5, 3.9), mesh_spacing);
  cosmosim::gravity::TreePmDiagnostics pm_only_diagnostics;
  const ForceField pm_only = solveTreePm(pos_x, pos_y, pos_z, mass, pm_shape, pm_only_options, &pm_only_diagnostics);

  cosmosim::gravity::TreePmOptions split_options = base_options;
  split_options.split_policy = cosmosim::gravity::makeTreePmSplitPolicyFromMeshSpacing(
      1.25, rcutCellsForAvailableBackend(4.5, 3.9), mesh_spacing);
  cosmosim::gravity::TreePmDiagnostics split_diagnostics;
  const ForceField split = solveTreePm(pos_x, pos_y, pos_z, mass, pm_shape, split_options, &split_diagnostics);

  const double pm_only_rel = relativeL2Error(pm_only, direct);
  const double split_rel = relativeL2Error(split, direct);

  std::ostringstream message;
  message << "TreePM consistency failure: split_rel=" << split_rel
          << ", pm_only_rel=" << pm_only_rel << ", split_pruned_nodes=" << split_diagnostics.residual_pruned_nodes
          << ", split_pair_skips=" << split_diagnostics.residual_pair_skips_cutoff;

  requireOrThrow(split_rel < 0.90, message.str());
  requireOrThrow(pm_only_rel < 2.2, message.str());
  // This integration fixture uses a minimum-image direct sum, not the
  // repository's Ewald-periodic reference. PM contains periodic image forces,
  // so being closer to this deliberately incomplete reference is not a
  // monotonic measure of TreePM split complementarity. Keep the broad seam /
  // split regression gates here; periodic force accuracy and complementarity
  // are certified independently by test_tree_pm_ewald_accuracy.
  requireOrThrow(split_diagnostics.residual_pair_evaluations > 0U, message.str());
}

void testCosmologicalTreePmNormalizationAndSplitComplementarity() {
  constexpr double box_size_comoving = 1.0;
  constexpr double non_unit_scale_factor = 0.5;
  const std::vector<double> pos_x = {
      0.030, 0.052, 0.105, 0.188, 0.327, 0.469,
      0.531, 0.674, 0.812, 0.895, 0.948, 0.970,
  };
  const std::vector<double> pos_y = {
      0.120, 0.137, 0.490, 0.715, 0.238, 0.881,
      0.862, 0.401, 0.643, 0.292, 0.108, 0.091,
  };
  const std::vector<double> pos_z = {
      0.760, 0.744, 0.205, 0.349, 0.917, 0.558,
      0.579, 0.083, 0.426, 0.681, 0.816, 0.834,
  };
  const std::vector<double> mass = {
      0.80, 1.15, 0.95, 1.30, 0.72, 1.08,
      0.91, 1.22, 0.86, 1.11, 0.77, 1.04,
  };

  const cosmosim::gravity::PmGridShape pm_shape = pmShapeForAvailableBackend(32, 8);
  const double mesh_spacing = box_size_comoving / static_cast<double>(pm_shape.nx);
  cosmosim::gravity::TreePmOptions unit_options;
  unit_options.pm_options.box_size_mpc_comoving = box_size_comoving;
  unit_options.pm_options.scale_factor = 1.0;
  unit_options.pm_options.gravitational_constant_code = 1.0;
  unit_options.pm_options.enable_window_deconvolution = true;
  unit_options.tree_options.opening_theta = 0.55;
  unit_options.tree_options.max_leaf_size = 2U;
  unit_options.tree_options.gravitational_constant_code = 1.0;
  unit_options.tree_options.softening.epsilon_comoving = 2.0e-3;
  unit_options.split_policy =
      cosmosim::gravity::makeTreePmSplitPolicyFromMeshSpacing(
          1.25, rcutCellsForAvailableBackend(6.25, 3.9), mesh_spacing);

  cosmosim::gravity::TreePmDiagnostics unit_diagnostics;
  const ForceField unit_force =
      solveTreePm(pos_x, pos_y, pos_z, mass, pm_shape, unit_options, &unit_diagnostics);

  cosmosim::gravity::TreePmOptions scaled_options = unit_options;
  scaled_options.pm_options.scale_factor = non_unit_scale_factor;
  cosmosim::gravity::TreePmDiagnostics scaled_diagnostics;
  const ForceField scaled_force =
      solveTreePm(pos_x, pos_y, pos_z, mass, pm_shape, scaled_options, &scaled_diagnostics);

  const double expected_force_scale = 1.0;
  ForceField expected_scaled_force = unit_force;
  for (std::size_t i = 0; i < expected_scaled_force.ax.size(); ++i) {
    expected_scaled_force.ax[i] *= expected_force_scale;
    expected_scaled_force.ay[i] *= expected_force_scale;
    expected_scaled_force.az[i] *= expected_force_scale;
  }
  const double total_scaling_rel_l2 = relativeL2Error(scaled_force, expected_scaled_force);
  const auto component_scaling_error = [expected_force_scale](double scaled, double unit) {
    const double expected = expected_force_scale * unit;
    return std::abs(scaled - expected) / std::max(std::abs(expected), 1.0e-30);
  };
  const auto norm_scaling_error = [expected_force_scale](double scaled, double unit) {
    return std::abs(scaled - expected_force_scale * unit) /
        std::max(expected_force_scale * unit, 1.0e-30);
  };
  const double pm_norm_scaling_error = norm_scaling_error(
      scaled_diagnostics.force_l2_pm_global, unit_diagnostics.force_l2_pm_global);
  const double tree_norm_scaling_error = norm_scaling_error(
      scaled_diagnostics.force_l2_tree_short_range,
      unit_diagnostics.force_l2_tree_short_range);
  const double total_norm_scaling_error = norm_scaling_error(
      scaled_diagnostics.force_l2_total, unit_diagnostics.force_l2_total);
  double max_component_scaling_error = 0.0;
  for (std::size_t i = 0; i < unit_force.ax.size(); ++i) {
    for (const auto [scaled, unit] : std::array{
             std::pair{scaled_force.ax[i], unit_force.ax[i]},
             std::pair{scaled_force.ay[i], unit_force.ay[i]},
             std::pair{scaled_force.az[i], unit_force.az[i]}}) {
      if (std::abs(expected_force_scale * unit) > 1.0e-12 * unit_diagnostics.force_l2_total) {
        max_component_scaling_error =
            std::max(max_component_scaling_error, component_scaling_error(scaled, unit));
      }
    }
  }

  const double pm_fraction = unit_diagnostics.force_l2_pm_global /
      (unit_diagnostics.force_l2_pm_global + unit_diagnostics.force_l2_tree_short_range);
  const double tree_fraction = unit_diagnostics.force_l2_tree_short_range /
      (unit_diagnostics.force_l2_pm_global + unit_diagnostics.force_l2_tree_short_range);
  std::ostringstream message;
  message << "cosmological TreePM normalization/split regression failed: a="
          << non_unit_scale_factor
          << ", expected_force_scale=" << expected_force_scale
          << ", total_scaling_rel_l2=" << total_scaling_rel_l2
          << ", max_component_scaling_error=" << max_component_scaling_error
          << ", pm_norm_scaling_error=" << pm_norm_scaling_error
          << ", tree_norm_scaling_error=" << tree_norm_scaling_error
          << ", total_norm_scaling_error=" << total_norm_scaling_error
          << ", unit_pm_l2=" << unit_diagnostics.force_l2_pm_global
          << ", unit_tree_l2=" << unit_diagnostics.force_l2_tree_short_range
          << ", pm_fraction=" << pm_fraction
          << ", tree_fraction=" << tree_fraction
          << ", unit_composition_error=" << unit_diagnostics.max_relative_composition_error
          << ", scaled_composition_error=" << scaled_diagnostics.max_relative_composition_error;

  requireOrThrow(unit_diagnostics.force_l2_total > 0.0, message.str());
  requireOrThrow(pm_fraction > 1.0e-2 && tree_fraction > 1.0e-2, message.str());
  requireOrThrow(total_scaling_rel_l2 <= 2.0e-12, message.str());
  requireOrThrow(max_component_scaling_error <= 2.0e-11, message.str());
  requireOrThrow(pm_norm_scaling_error <= 2.0e-12, message.str());
  requireOrThrow(tree_norm_scaling_error <= 2.0e-12, message.str());
  requireOrThrow(total_norm_scaling_error <= 2.0e-12, message.str());
  requireOrThrow(unit_diagnostics.max_relative_composition_error <= 1.0e-12, message.str());
  requireOrThrow(scaled_diagnostics.max_relative_composition_error <= 1.0e-12, message.str());

  bool emit_metrics = true;
#if COSMOSIM_ENABLE_MPI
  int world_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  emit_metrics = world_rank == 0;
#endif
  if (emit_metrics) {
    std::cout << "TREEPM_COSMOLOGICAL_NORMALIZATION_PASS"
              << " a=" << non_unit_scale_factor
              << " expected_force_scale=" << expected_force_scale
              << " total_scaling_rel_l2=" << total_scaling_rel_l2
              << " max_component_scaling_error=" << max_component_scaling_error
              << " pm_norm_scaling_error=" << pm_norm_scaling_error
              << " tree_norm_scaling_error=" << tree_norm_scaling_error
              << " total_norm_scaling_error=" << total_norm_scaling_error
              << " pm_fraction=" << pm_fraction
              << " tree_fraction=" << tree_fraction
              << '\n';
  }
}

void testTreePmRejectsMismatchedGravitationalConstants() {
  const std::vector<double> pos_x = {0.20, 0.35};
  const std::vector<double> pos_y = {0.30, 0.42};
  const std::vector<double> pos_z = {0.40, 0.53};
  const std::vector<double> mass = {1.0, 1.0};
  const cosmosim::gravity::PmGridShape pm_shape{8U, 8U, 8U};
  cosmosim::gravity::TreePmOptions options;
  options.pm_options.box_size_mpc_comoving = 1.0;
  options.pm_options.gravitational_constant_code = 1.0;
  options.tree_options.gravitational_constant_code = 1.25;
  options.split_policy =
      cosmosim::gravity::makeTreePmSplitPolicyFromMeshSpacing(1.25, 3.9, 0.125);

  bool rejected = false;
  try {
    [[maybe_unused]] const ForceField field =
        solveTreePm(pos_x, pos_y, pos_z, mass, pm_shape, options, nullptr);
  } catch (const std::invalid_argument& error) {
    rejected = std::string(error.what()).find("gravitational_constant_code values must match") !=
        std::string::npos;
  }
  requireOrThrow(
      rejected,
      "TreePM must reject mismatched PM/tree gravitational constants before composing the force split");
}

void testPeriodicTreePmRejectsAmbiguousMinimumImageCutoff() {
  const std::vector<double> pos_x = {0.20, 0.35};
  const std::vector<double> pos_y = {0.30, 0.42};
  const std::vector<double> pos_z = {0.40, 0.53};
  const std::vector<double> mass = {1.0, 1.0};
  const cosmosim::gravity::PmGridShape pm_shape{8U, 8U, 8U};
  cosmosim::gravity::TreePmOptions options;
  options.pm_options.box_size_mpc_comoving = 1.0;
  options.pm_options.gravitational_constant_code = 1.0;
  options.tree_options.gravitational_constant_code = 1.0;
  options.split_policy =
      cosmosim::gravity::makeTreePmSplitPolicyFromMeshSpacing(1.25, 4.0, 0.125);

  bool rejected = false;
  try {
    [[maybe_unused]] const ForceField field =
        solveTreePm(pos_x, pos_y, pos_z, mass, pm_shape, options, nullptr);
  } catch (const std::invalid_argument& error) {
    rejected = std::string(error.what()).find("< half the shortest box axis") !=
        std::string::npos;
  }
  requireOrThrow(
      rejected,
      "periodic TreePM must reject r_cut=L_min/2 before choosing one of two equidistant images");
}

void testActiveSubsetMatchesFullSolveSingleRank() {
  constexpr double box_size_comoving = 1.0;
  const cosmosim::gravity::PmGridShape pm_shape = pmShapeForAvailableBackend(16, 8);
  const double mesh_spacing = box_size_comoving / static_cast<double>(pm_shape.nx);

  const std::vector<double> pos_x = {0.05, 0.12, 0.19, 0.28, 0.34, 0.43, 0.55, 0.61, 0.74, 0.88};
  const std::vector<double> pos_y = {0.08, 0.16, 0.24, 0.31, 0.37, 0.46, 0.52, 0.67, 0.79, 0.91};
  const std::vector<double> pos_z = {0.03, 0.11, 0.21, 0.29, 0.41, 0.49, 0.58, 0.69, 0.82, 0.94};
  std::vector<double> mass(pos_x.size(), 1.0);
  for (std::size_t i = 0; i < mass.size(); ++i) {
    mass[i] = 0.95 + 0.015 * static_cast<double>(i % 5U);
  }

  cosmosim::gravity::TreePmOptions options;
  options.pm_options.box_size_mpc_comoving = box_size_comoving;
  options.pm_options.scale_factor = 1.0;
  options.pm_options.gravitational_constant_code = 1.0;
  options.pm_options.enable_window_deconvolution = true;
  options.tree_options.opening_theta = 0.6;
  options.tree_options.max_leaf_size = 4;
  options.tree_options.gravitational_constant_code = 1.0;
  options.tree_options.softening.epsilon_comoving = 1.0e-3;
  options.split_policy = cosmosim::gravity::makeTreePmSplitPolicyFromMeshSpacing(
      1.5, rcutCellsForAvailableBackend(4.0, 3.9), mesh_spacing);

  const ForceField full_field = solveTreePm(pos_x, pos_y, pos_z, mass, pm_shape, options, nullptr);

  std::vector<std::uint32_t> active_subset = {1U, 4U, 7U};
  ForceField subset_field{
      std::vector<double>(active_subset.size(), 0.0),
      std::vector<double>(active_subset.size(), 0.0),
      std::vector<double>(active_subset.size(), 0.0),
  };
  cosmosim::gravity::TreePmForceAccumulatorView subset_accumulator{
      active_subset,
      subset_field.ax,
      subset_field.ay,
      subset_field.az,
  };
  cosmosim::gravity::TreePmCoordinator subset_coordinator(pm_shape);
  subset_coordinator.solveActiveSet(pos_x, pos_y, pos_z, mass, subset_accumulator, options, nullptr, nullptr);

  for (std::size_t i = 0; i < active_subset.size(); ++i) {
    const std::uint32_t global_index = active_subset[i];
    requireOrThrow(std::abs(subset_field.ax[i] - full_field.ax[global_index]) <= 1.0e-12, "active-subset ax mismatch");
    requireOrThrow(std::abs(subset_field.ay[i] - full_field.ay[global_index]) <= 1.0e-12, "active-subset ay mismatch");
    requireOrThrow(std::abs(subset_field.az[i] - full_field.az[global_index]) <= 1.0e-12, "active-subset az mismatch");
  }
}

void testLongRangeCadenceCacheIsFailClosedAndOwnershipCompatible() {
#if COSMOSIM_ENABLE_MPI
  int world_size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  if (world_size != 1) {
    return;
  }
#endif
  const cosmosim::gravity::PmGridShape shape =
      pmShapeForAvailableBackend(8U, 8U);
  const std::vector<double> pos_x{0.15, 0.62};
  const std::vector<double> pos_y{0.28, 0.71};
  const std::vector<double> pos_z{0.41, 0.84};
  const std::vector<double> mass{1.0, 1.5};
  const std::vector<std::uint32_t> active{0U, 1U};
  ForceField force{
      std::vector<double>(active.size(), 0.0),
      std::vector<double>(active.size(), 0.0),
      std::vector<double>(active.size(), 0.0)};
  const cosmosim::gravity::TreePmForceAccumulatorView accumulator{
      .active_particle_index = active,
      .accel_x_comoving = force.ax,
      .accel_y_comoving = force.ay,
      .accel_z_comoving = force.az,
  };
  cosmosim::gravity::TreePmOptions options;
  options.pm_options.box_size_mpc_comoving = 1.0;
  options.pm_options.box_size_x_mpc_comoving = 1.0;
  options.pm_options.box_size_y_mpc_comoving = 1.0;
  options.pm_options.box_size_z_mpc_comoving = 1.0;
  options.pm_options.scale_factor = 0.5;
  options.pm_options.gravitational_constant_code = 1.0;
  options.tree_options.gravitational_constant_code = 1.0;
  options.tree_options.softening.epsilon_comoving = 1.0e-3;
  options.split_policy = cosmosim::gravity::makeTreePmSplitPolicyFromMeshSpacing(
      1.25, 3.9, 1.0 / static_cast<double>(shape.nx));
  options.decomposition_epoch = cosmosim::gravity::DecompositionEpoch{3U};
  options.force_epoch = cosmosim::gravity::ForceEvaluationEpoch{.sequence = 7U, .scale_factor = 1.0};

  cosmosim::gravity::TreePmCoordinator coordinator(shape);
  cosmosim::gravity::TreePmDiagnostics diagnostics;
  coordinator.solveActiveSetWithPmCadence(
      pos_x, pos_y, pos_z, mass, accumulator, options, true, nullptr,
      &diagnostics);
  requireOrThrow(diagnostics.pm_solve_count == 1U,
                 "authorized cadence refresh did not solve PM");

  force.ax.assign(active.size(), 0.0);
  force.ay.assign(active.size(), 0.0);
  force.az.assign(active.size(), 0.0);
  options.decomposition_epoch = cosmosim::gravity::DecompositionEpoch{4U};
  coordinator.solveActiveSetWithPmCadence(
      pos_x, pos_y, pos_z, mass, accumulator, options, false, nullptr,
      &diagnostics);
  requireOrThrow(diagnostics.pm_reuse_count == 1U,
                 "particle-ownership epoch incorrectly invalidated fixed-slab PM cache");

  options.force_epoch = cosmosim::gravity::ForceEvaluationEpoch{.sequence = 8U, .scale_factor = 1.0};
  bool incompatible_reuse_threw = false;
  try {
    coordinator.solveActiveSetWithPmCadence(
        pos_x, pos_y, pos_z, mass, accumulator, options, false, nullptr,
        &diagnostics);
  } catch (const std::runtime_error& ex) {
    incompatible_reuse_threw =
        std::string(ex.what()).find("without a compatible PM field") !=
        std::string::npos;
  }
  requireOrThrow(incompatible_reuse_threw,
                 "incompatible PM cache reuse did not fail closed");
}

#if COSMOSIM_ENABLE_MPI
void testRankLocalSerialTreePmInsideMpiWorldSkipsDistributedArena() {
  int world_size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  if (world_size < 2 || world_size > 8) {
    return;
  }

  const cosmosim::gravity::PmGridShape pm_shape =
      pmShapeForAvailableBackend(8U, 8U);
  const std::vector<double> pos_x{0.17, 0.63};
  const std::vector<double> pos_y{0.29, 0.74};
  const std::vector<double> pos_z{0.41, 0.86};
  const std::vector<double> mass{1.0, 1.25};
  const std::vector<std::uint32_t> active{0U, 1U};
  ForceField force{
      std::vector<double>(active.size(), 0.0),
      std::vector<double>(active.size(), 0.0),
      std::vector<double>(active.size(), 0.0)};
  const cosmosim::gravity::TreePmForceAccumulatorView accumulator{
      .active_particle_index = active,
      .accel_x_comoving = force.ax,
      .accel_y_comoving = force.ay,
      .accel_z_comoving = force.az,
  };

  cosmosim::gravity::TreePmOptions options;
  options.pm_options.box_size_mpc_comoving = 1.0;
  options.pm_options.scale_factor = 1.0;
  options.pm_options.gravitational_constant_code = 1.0;
  options.tree_options.gravitational_constant_code = 1.0;
  options.tree_options.softening.epsilon_comoving = 1.0e-3;
  options.split_policy = cosmosim::gravity::makeTreePmSplitPolicyFromMeshSpacing(
      1.25, rcutCellsForAvailableBackend(3.0, 3.0),
      1.0 / static_cast<double>(pm_shape.nx));

  cosmosim::gravity::TreePmCoordinator coordinator(pm_shape);
  requireOrThrow(
      coordinator.slabLayout().world_size == 1 && coordinator.ownsFullPmDomain(),
      "rank-local-serial TreePM regression requires one-rank full-domain PM ownership");

  cosmosim::gravity::TreePmDiagnostics diagnostics;
  coordinator.solveActiveSet(
      pos_x, pos_y, pos_z, mass, accumulator, options, nullptr, &diagnostics);

  requireOrThrow(
      diagnostics.pm_halo_value_count == 0U,
      "rank-local-serial TreePM unexpectedly entered PM halo exchange");
  requireOrThrow(
      diagnostics.residual_remote_request_packets == 0U &&
          diagnostics.residual_remote_response_packets == 0U,
      "rank-local-serial TreePM unexpectedly entered distributed residual exchange");

  const auto memory_report = coordinator.memoryReport();
  const auto physical_arena = std::find_if(
      memory_report.entries.begin(), memory_report.entries.end(),
      [](const auto& entry) {
        return entry.label == "treepm.communication_arena";
      });
  requireOrThrow(
      physical_arena == memory_report.entries.end(),
      "rank-local-serial TreePM configured the distributed gravity communication arena");
}

void testDivergentTreePmLayoutFailsCoordinately() {
  int world_size = 1;
  int world_rank = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  if (world_size < 2 || world_size > 8) {
    return;
  }

  const cosmosim::gravity::PmGridShape pm_shape{8U, 8U, 8U};
  cosmosim::gravity::TreePmOptions options;
  options.pm_options.box_size_mpc_comoving = 1.0;
  options.pm_options.gravitational_constant_code = 1.0;
  options.tree_options.gravitational_constant_code = 1.0;
  options.tree_options.softening.epsilon_comoving = 1.0e-3;
  options.split_policy = cosmosim::gravity::makeTreePmSplitPolicyFromMeshSpacing(
      1.25, 3.9, 1.0 / static_cast<double>(pm_shape.nx));

  const std::vector<double> empty;
  const std::vector<std::uint32_t> no_targets;
  ForceField no_force;
  const cosmosim::gravity::TreePmForceAccumulatorView accumulator{
      .active_particle_index = no_targets,
      .accel_x_comoving = no_force.ax,
      .accel_y_comoving = no_force.ay,
      .accel_z_comoving = no_force.az,
  };

  bool rejected = false;
  std::string diagnostic;
  try {
    if (world_rank == 0) {
      // Deliberately advertise a rank-local serial layout while every peer
      // advertises the communicator-wide slab layout. All ranks must vote the
      // mismatch before any rank enters a layout-selected PM collective.
      cosmosim::gravity::TreePmCoordinator coordinator(pm_shape);
      coordinator.solveActiveSet(
          empty, empty, empty, empty, accumulator, options, nullptr, nullptr);
    } else {
      const auto layout = cosmosim::parallel::makePmSlabLayout(
          pm_shape.nx, pm_shape.ny, pm_shape.nz, world_size, world_rank);
      cosmosim::gravity::TreePmCoordinator coordinator(pm_shape, layout);
      coordinator.solveActiveSet(
          empty, empty, empty, empty, accumulator, options, nullptr, nullptr);
    }
  } catch (const std::exception& ex) {
    rejected = true;
    diagnostic = ex.what();
  }

  const std::uint64_t local_rejection = rejected ? 1U : 0U;
  std::uint64_t rejection_count = 0U;
  MPI_Allreduce(
      &local_rejection,
      &rejection_count,
      1,
      MPI_UINT64_T,
      MPI_SUM,
      MPI_COMM_WORLD);
  requireOrThrow(
      rejection_count == static_cast<std::uint64_t>(world_size),
      "divergent TreePM layouts were not rejected on every rank");
  requireOrThrow(
      diagnostic.find("world metadata") != std::string::npos ||
          diagnostic.find("disagree") != std::string::npos ||
          diagnostic.find("peer rank rejected") != std::string::npos,
      "divergent TreePM layout rejection omitted the coordinated preflight diagnostic");

  const cosmosim::gravity::PmGridShape rank_local_shape =
      world_rank == 0 ? cosmosim::gravity::PmGridShape{8U, 8U, 8U}
                      : cosmosim::gravity::PmGridShape{10U, 8U, 8U};
  const auto rank_local_layout = cosmosim::parallel::makePmSlabLayout(
      rank_local_shape.nx,
      rank_local_shape.ny,
      rank_local_shape.nz,
      world_size,
      world_rank);
  bool shape_rejected = false;
  std::string shape_diagnostic;
  try {
    cosmosim::gravity::TreePmCoordinator coordinator(
        rank_local_shape, rank_local_layout);
    coordinator.solveActiveSet(
        empty, empty, empty, empty, accumulator, options, nullptr, nullptr);
  } catch (const std::exception& ex) {
    shape_rejected = true;
    shape_diagnostic = ex.what();
  }
  const std::uint64_t local_shape_rejection = shape_rejected ? 1U : 0U;
  std::uint64_t shape_rejection_count = 0U;
  MPI_Allreduce(
      &local_shape_rejection,
      &shape_rejection_count,
      1,
      MPI_UINT64_T,
      MPI_SUM,
      MPI_COMM_WORLD);
  requireOrThrow(
      shape_rejection_count == static_cast<std::uint64_t>(world_size),
      "divergent TreePM mesh shapes were not rejected on every rank");
  requireOrThrow(
      shape_diagnostic.find("mesh shape") != std::string::npos ||
          shape_diagnostic.find("diverged") != std::string::npos,
      "divergent TreePM mesh shape rejection omitted the protocol diagnostic");
}

void testDistributedShortRangeExportImportMatchesSingleRankReference() {
  int world_size = 1;
  int world_rank = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  if (world_size != 2) {
    return;
  }

  constexpr double box_size_comoving = 1.0;
  const cosmosim::gravity::PmGridShape pm_shape = pmShapeForAvailableBackend(16, 8);
  const double mesh_spacing = box_size_comoving / static_cast<double>(pm_shape.nx);

  const std::vector<double> global_x = {
      0.08, 0.14, 0.22, 0.30, 0.41, 0.47, 0.492, 0.498,
      0.502, 0.508, 0.53, 0.61, 0.72, 0.81, 0.90, 0.96,
  };
  const std::vector<double> global_y = {
      0.11, 0.19, 0.27, 0.33, 0.39, 0.46, 0.49, 0.52,
      0.55, 0.58, 0.63, 0.67, 0.71, 0.76, 0.84, 0.93,
  };
  const std::vector<double> global_z = {
      0.07, 0.13, 0.20, 0.28, 0.35, 0.44, 0.495, 0.499,
      0.501, 0.505, 0.57, 0.64, 0.70, 0.79, 0.87, 0.95,
  };
  std::vector<double> global_mass(global_x.size(), 1.0);
  for (std::size_t i = 0; i < global_mass.size(); ++i) {
    global_mass[i] = 0.9 + 0.02 * static_cast<double>(i % 7U);
  }

  std::vector<double> local_x;
  std::vector<double> local_y;
  std::vector<double> local_z;
  std::vector<double> local_mass;
  std::vector<std::size_t> local_to_global;
  for (std::size_t i = 0; i < global_x.size(); ++i) {
    const int owner_rank = (global_x[i] < 0.5) ? 0 : 1;
    if (owner_rank == world_rank) {
      local_x.push_back(global_x[i]);
      local_y.push_back(global_y[i]);
      local_z.push_back(global_z[i]);
      local_mass.push_back(global_mass[i]);
      local_to_global.push_back(i);
    }
  }

  std::vector<std::uint32_t> local_active(local_x.size(), 0U);
  for (std::size_t i = 0; i < local_active.size(); ++i) {
    local_active[i] = static_cast<std::uint32_t>(i);
  }
  ForceField local_field{
      std::vector<double>(local_x.size(), 0.0),
      std::vector<double>(local_x.size(), 0.0),
      std::vector<double>(local_x.size(), 0.0),
  };
  cosmosim::gravity::TreePmForceAccumulatorView local_accumulator{
      local_active,
      local_field.ax,
      local_field.ay,
      local_field.az,
  };

  cosmosim::gravity::TreePmOptions options;
  options.pm_options.box_size_mpc_comoving = box_size_comoving;
  options.pm_options.scale_factor = 1.0;
  options.pm_options.gravitational_constant_code = 1.0;
  options.tree_options.opening_theta = 0.55;
  options.tree_options.max_leaf_size = 4;
  options.tree_options.gravitational_constant_code = 1.0;
  options.tree_options.softening.epsilon_comoving = 1.0e-3;
  options.split_policy = cosmosim::gravity::makeTreePmSplitPolicyFromMeshSpacing(
      2.0, rcutCellsForAvailableBackend(4.0, 3.9), mesh_spacing);
  // Keep the exchange deliberately tiny (two packets), but respect the
  // current 96-byte short-range request wire format.
  options.tree_exchange_batch_bytes =
      2U * cosmosim::gravity::kTreePmShortRangeRequestWireBytes;

  options.full_residual_lookup_enabled=true; options.leaf_simd_enabled=true;
  options.radial_coefficients_enabled=true; options.full_mac_diagnostics=false;
  const auto local_layout =
      cosmosim::parallel::makePmSlabLayout(pm_shape.nx, pm_shape.ny, pm_shape.nz, world_size, world_rank);
  cosmosim::gravity::TreePmCoordinator distributed(pm_shape, local_layout);
  cosmosim::gravity::TreePmDiagnostics distributed_diag;
  distributed.solveActiveSet(local_x, local_y, local_z, local_mass, local_accumulator, options, nullptr, &distributed_diag);

  std::vector<std::uint32_t> ref_active(global_x.size(), 0U);
  for (std::size_t i = 0; i < ref_active.size(); ++i) {
    ref_active[i] = static_cast<std::uint32_t>(i);
  }
  ForceField reference_field{
      std::vector<double>(global_x.size(), 0.0),
      std::vector<double>(global_x.size(), 0.0),
      std::vector<double>(global_x.size(), 0.0),
  };
  cosmosim::gravity::TreePmForceAccumulatorView ref_accumulator{
      ref_active,
      reference_field.ax,
      reference_field.ay,
      reference_field.az,
  };
  cosmosim::gravity::TreePmCoordinator single_rank(pm_shape);
  cosmosim::gravity::TreePmDiagnostics reference_diag;
  auto reference_options=options; reference_options.full_residual_lookup_enabled=false;
  reference_options.leaf_simd_enabled=false; reference_options.radial_coefficients_enabled=false;
  reference_options.full_mac_diagnostics=true;
  single_rank.solveActiveSet(global_x, global_y, global_z, global_mass, ref_accumulator, reference_options, nullptr, &reference_diag);

  double local_diff2 = 0.0;
  double local_ref2 = 0.0;
  for (std::size_t i = 0; i < local_x.size(); ++i) {
    const std::size_t gi = local_to_global[i];
    const double dx = local_field.ax[i] - reference_field.ax[gi];
    const double dy = local_field.ay[i] - reference_field.ay[gi];
    const double dz = local_field.az[i] - reference_field.az[gi];
    local_diff2 += dx * dx + dy * dy + dz * dz;
    local_ref2 += reference_field.ax[gi] * reference_field.ax[gi] +
        reference_field.ay[gi] * reference_field.ay[gi] +
        reference_field.az[gi] * reference_field.az[gi];
  }

  double global_diff2 = 0.0;
  double global_ref2 = 0.0;
  MPI_Allreduce(&local_diff2, &global_diff2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&local_ref2, &global_ref2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  const double rel_l2 = std::sqrt(global_diff2 / std::max(global_ref2, 1.0e-30));

  std::ostringstream msg;
  msg << "Distributed TreePM short-range export/import drift: rel_l2=" << rel_l2
      << ", world_rank=" << world_rank
      << ", residual_pruned_nodes=" << distributed_diag.residual_pruned_nodes
      << ", residual_pair_skips_cutoff=" << distributed_diag.residual_pair_skips_cutoff
      << ", reference_pair_skips_cutoff=" << reference_diag.residual_pair_skips_cutoff;
  requireOrThrow(rel_l2 <= 1.0e-9, msg.str());
  requireOrThrow(distributed_diag.residual_pair_skips_cutoff > 0, msg.str());
}

void testDistributedActiveSubsetMatchesSingleRankReference() {
  int world_size = 1;
  int world_rank = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  if (world_size != 2) {
    return;
  }

  constexpr double box_size_comoving = 1.0;
  const cosmosim::gravity::PmGridShape pm_shape = pmShapeForAvailableBackend(16, 8);
  const double mesh_spacing = box_size_comoving / static_cast<double>(pm_shape.nx);

  const std::vector<double> global_x = {
      0.06, 0.13, 0.22, 0.29, 0.38, 0.46, 0.495, 0.499,
      0.501, 0.507, 0.56, 0.62, 0.71, 0.80, 0.89, 0.97,
  };
  const std::vector<double> global_y = {
      0.09, 0.18, 0.26, 0.32, 0.40, 0.48, 0.50, 0.53,
      0.55, 0.59, 0.64, 0.68, 0.73, 0.78, 0.85, 0.92,
  };
  const std::vector<double> global_z = {
      0.04, 0.14, 0.23, 0.30, 0.36, 0.45, 0.494, 0.498,
      0.502, 0.509, 0.58, 0.65, 0.74, 0.82, 0.88, 0.96,
  };
  std::vector<double> global_mass(global_x.size(), 1.0);
  for (std::size_t i = 0; i < global_mass.size(); ++i) {
    global_mass[i] = 0.92 + 0.02 * static_cast<double>(i % 6U);
  }

  std::vector<double> local_x;
  std::vector<double> local_y;
  std::vector<double> local_z;
  std::vector<double> local_mass;
  std::vector<std::size_t> local_to_global;
  for (std::size_t i = 0; i < global_x.size(); ++i) {
    const int owner_rank = (global_x[i] < 0.5) ? 0 : 1;
    if (owner_rank == world_rank) {
      local_x.push_back(global_x[i]);
      local_y.push_back(global_y[i]);
      local_z.push_back(global_z[i]);
      local_mass.push_back(global_mass[i]);
      local_to_global.push_back(i);
    }
  }

  std::vector<std::uint32_t> local_active;
  std::vector<std::size_t> local_active_to_global;
  for (std::size_t li = 0; li < local_to_global.size(); ++li) {
    if ((local_to_global[li] % 2U) == 0U) {
      local_active.push_back(static_cast<std::uint32_t>(li));
      local_active_to_global.push_back(local_to_global[li]);
    }
  }

  ForceField local_field{
      std::vector<double>(local_active.size(), 0.0),
      std::vector<double>(local_active.size(), 0.0),
      std::vector<double>(local_active.size(), 0.0),
  };
  cosmosim::gravity::TreePmForceAccumulatorView local_accumulator{
      local_active,
      local_field.ax,
      local_field.ay,
      local_field.az,
  };

  cosmosim::gravity::TreePmOptions options;
  options.pm_options.box_size_mpc_comoving = box_size_comoving;
  options.pm_options.scale_factor = 1.0;
  options.pm_options.gravitational_constant_code = 1.0;
  options.tree_options.opening_theta = 0.55;
  options.tree_options.max_leaf_size = 4;
  options.tree_options.gravitational_constant_code = 1.0;
  options.tree_options.softening.epsilon_comoving = 1.0e-3;
  options.split_policy = cosmosim::gravity::makeTreePmSplitPolicyFromMeshSpacing(
      2.0, rcutCellsForAvailableBackend(4.0, 3.9), mesh_spacing);
  // Keep the exchange deliberately tiny (two packets), but respect the
  // current 96-byte short-range request wire format.
  options.tree_exchange_batch_bytes =
      2U * cosmosim::gravity::kTreePmShortRangeRequestWireBytes;

  const auto local_layout =
      cosmosim::parallel::makePmSlabLayout(pm_shape.nx, pm_shape.ny, pm_shape.nz, world_size, world_rank);
  cosmosim::gravity::TreePmCoordinator distributed(pm_shape, local_layout);
  distributed.solveActiveSet(local_x, local_y, local_z, local_mass, local_accumulator, options, nullptr, nullptr);

  std::vector<std::uint32_t> ref_active(global_x.size(), 0U);
  for (std::size_t i = 0; i < ref_active.size(); ++i) {
    ref_active[i] = static_cast<std::uint32_t>(i);
  }
  ForceField reference_field{
      std::vector<double>(global_x.size(), 0.0),
      std::vector<double>(global_x.size(), 0.0),
      std::vector<double>(global_x.size(), 0.0),
  };
  cosmosim::gravity::TreePmForceAccumulatorView ref_accumulator{
      ref_active,
      reference_field.ax,
      reference_field.ay,
      reference_field.az,
  };
  cosmosim::gravity::TreePmCoordinator single_rank(pm_shape);
  single_rank.solveActiveSet(global_x, global_y, global_z, global_mass, ref_accumulator, options, nullptr, nullptr);

  double local_diff2 = 0.0;
  double local_ref2 = 0.0;
  for (std::size_t i = 0; i < local_active_to_global.size(); ++i) {
    const std::size_t gi = local_active_to_global[i];
    const double dx = local_field.ax[i] - reference_field.ax[gi];
    const double dy = local_field.ay[i] - reference_field.ay[gi];
    const double dz = local_field.az[i] - reference_field.az[gi];
    local_diff2 += dx * dx + dy * dy + dz * dz;
    local_ref2 += reference_field.ax[gi] * reference_field.ax[gi] +
        reference_field.ay[gi] * reference_field.ay[gi] +
        reference_field.az[gi] * reference_field.az[gi];
  }
  double global_diff2 = 0.0;
  double global_ref2 = 0.0;
  MPI_Allreduce(&local_diff2, &global_diff2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&local_ref2, &global_ref2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  const double rel_l2 = std::sqrt(global_diff2 / std::max(global_ref2, 1.0e-30));
  requireOrThrow(rel_l2 <= 1.0e-9, "distributed active-subset vs single-rank mismatch");
}

void testDistributedEmptyRanksIndependentTargetsAndMigration() {
  int world_size = 1;
  int world_rank = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  if (world_size < 2 || world_size > 8) {
    return;
  }

  const cosmosim::gravity::PmGridShape pm_shape{8U, 8U, 8U};
  cosmosim::gravity::TreePmOptions options;
  options.pm_options.box_size_mpc_comoving = 1.0;
  options.pm_options.gravitational_constant_code = 1.0;
  options.pm_options.assignment_scheme = cosmosim::gravity::PmAssignmentScheme::kTsc;
  options.pm_options.enable_window_deconvolution = true;
  options.tree_options.opening_criterion = cosmosim::gravity::TreeOpeningCriterion::kRelativeForceError;
  options.tree_options.relative_force_tolerance = 0.005;
  options.tree_options.relative_force_acceleration_floor_code = 1.0e-12;
  options.tree_options.max_leaf_size = 1U;
  options.tree_options.softening.epsilon_comoving = 1.0e-3;
  options.split_policy = cosmosim::gravity::makeTreePmSplitPolicyFromMeshSpacing(1.25, 3.9, 0.125);
  options.decomposition_epoch = cosmosim::gravity::DecompositionEpoch{7U};
  options.source_generation = cosmosim::gravity::GravitySourceGeneration{11U};
  options.pm_field_version = cosmosim::gravity::PmFieldVersion{11U};
  options.force_epoch = cosmosim::gravity::ForceEvaluationEpoch{.sequence = 11U, .scale_factor = 1.0};
  options.tree_exchange_batch_bytes = 256U;

  const std::size_t global_source_count = world_size == 2 ? 1U : 2U;
  std::vector<double> global_x(global_source_count, 0.0);
  std::vector<double> global_y(global_source_count, 0.0);
  std::vector<double> global_z(global_source_count, 0.0);
  std::vector<double> global_mass(global_source_count, 0.0);
  for (std::size_t i = 0; i < global_source_count; ++i) {
    global_x[i] = 0.992 - 0.014 * static_cast<double>(i);
    global_y[i] = 0.30 + 0.025 * static_cast<double>(i);
    global_z[i] = 0.70 - 0.018 * static_cast<double>(i);
    global_mass[i] = 0.9 + 0.2 * static_cast<double>(i);
  }

  std::vector<double> local_x;
  std::vector<double> local_y;
  std::vector<double> local_z;
  std::vector<double> local_mass;
  if (world_rank == world_size - 1) {
    local_x = global_x;
    local_y = global_y;
    local_z = global_z;
    local_mass = global_mass;
  }
  const std::vector<std::uint32_t> independent_target = world_rank == 0
      ? std::vector<std::uint32_t>{std::numeric_limits<std::uint32_t>::max()}
      : std::vector<std::uint32_t>{};
  const std::vector<double> target_x = world_rank == 0 ? std::vector<double>{0.008} : std::vector<double>{};
  const std::vector<double> target_y = world_rank == 0 ? std::vector<double>{0.31} : std::vector<double>{};
  const std::vector<double> target_z = world_rank == 0 ? std::vector<double>{0.69} : std::vector<double>{};
  const std::vector<double> previous_acceleration = world_rank == 0
      ? std::vector<double>{1.0e6}
      : std::vector<double>{};
  ForceField distributed_target_force{
      std::vector<double>(independent_target.size(), 0.0),
      std::vector<double>(independent_target.size(), 0.0),
      std::vector<double>(independent_target.size(), 0.0),
  };
  const cosmosim::gravity::TreePmForceAccumulatorView target_accumulator{
      .active_particle_index = independent_target,
      .accel_x_comoving = distributed_target_force.ax,
      .accel_y_comoving = distributed_target_force.ay,
      .accel_z_comoving = distributed_target_force.az,
      .previous_acceleration_magnitude_code = previous_acceleration,
      .target_pos_x_comoving = target_x,
      .target_pos_y_comoving = target_y,
      .target_pos_z_comoving = target_z,
  };
  const auto layout =
      cosmosim::parallel::makePmSlabLayout(pm_shape.nx, pm_shape.ny, pm_shape.nz, world_size, world_rank);
  cosmosim::gravity::TreePmCoordinator distributed(pm_shape, layout);
  cosmosim::gravity::TreePmDiagnostics distributed_diagnostics;
  distributed.solveActiveSet(
      local_x,
      local_y,
      local_z,
      local_mass,
      target_accumulator,
      options,
      nullptr,
      &distributed_diagnostics);

  const std::vector<std::uint32_t> reference_target{std::numeric_limits<std::uint32_t>::max()};
  const std::vector<double> reference_target_x{0.008};
  const std::vector<double> reference_target_y{0.31};
  const std::vector<double> reference_target_z{0.69};
  const std::vector<double> reference_previous{1.0e6};
  ForceField reference_force{{0.0}, {0.0}, {0.0}};
  const cosmosim::gravity::TreePmForceAccumulatorView reference_accumulator{
      .active_particle_index = reference_target,
      .accel_x_comoving = reference_force.ax,
      .accel_y_comoving = reference_force.ay,
      .accel_z_comoving = reference_force.az,
      .previous_acceleration_magnitude_code = reference_previous,
      .target_pos_x_comoving = reference_target_x,
      .target_pos_y_comoving = reference_target_y,
      .target_pos_z_comoving = reference_target_z,
  };
  cosmosim::gravity::TreePmCoordinator reference(pm_shape);
  reference.solveActiveSet(
      global_x, global_y, global_z, global_mass, reference_accumulator, options, nullptr, nullptr);

  double local_max_error = 0.0;
  if (world_rank == 0) {
    local_max_error = std::max({
        std::abs(distributed_target_force.ax[0] - reference_force.ax[0]),
        std::abs(distributed_target_force.ay[0] - reference_force.ay[0]),
        std::abs(distributed_target_force.az[0] - reference_force.az[0]),
    });
    requireOrThrow(
        distributed_diagnostics.residual_remote_request_packets > 0U,
        "empty-source rank with an independent target did not issue a remote tree request");
  } else {
    requireOrThrow(independent_target.empty(), "non-target rank unexpectedly owns an active target");
  }
  double global_max_error = 0.0;
  MPI_Allreduce(&local_max_error, &global_max_error, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  requireOrThrow(global_max_error <= 2.0e-9, "empty-rank independent target differs from single-rank TreePM");

  // Exercise an all-empty collective round and distinguish explicit zero-source
  // descriptors from missing rank participation.
  const std::vector<double> empty;
  const std::vector<std::uint32_t> no_targets;
  ForceField no_force;
  const cosmosim::gravity::TreePmForceAccumulatorView empty_accumulator{
      .active_particle_index = no_targets,
      .accel_x_comoving = no_force.ax,
      .accel_y_comoving = no_force.ay,
      .accel_z_comoving = no_force.az,
  };
  cosmosim::gravity::TreePmCoordinator all_empty(pm_shape, layout);
  options.source_generation = cosmosim::gravity::GravitySourceGeneration{12U};
  options.pm_field_version = cosmosim::gravity::PmFieldVersion{12U};
  options.force_epoch = cosmosim::gravity::ForceEvaluationEpoch{.sequence = 12U, .scale_factor = 1.0};
  all_empty.solveActiveSet(empty, empty, empty, empty, empty_accumulator, options, nullptr, nullptr);

  // Reuse one coordinator across a decomposition transition in which every
  // nonzero rank becomes source-empty after the first force epoch.
  cosmosim::gravity::TreePmCoordinator migration(pm_shape, layout);
  std::vector<double> before_x{(static_cast<double>(world_rank) + 0.5) / static_cast<double>(world_size)};
  std::vector<double> before_y{0.2 + 0.1 * static_cast<double>(world_rank)};
  std::vector<double> before_z{0.4};
  std::vector<double> before_mass{1.0};
  options.decomposition_epoch = cosmosim::gravity::DecompositionEpoch{20U};
  options.source_generation = cosmosim::gravity::GravitySourceGeneration{20U};
  options.pm_field_version = cosmosim::gravity::PmFieldVersion{20U};
  options.force_epoch = cosmosim::gravity::ForceEvaluationEpoch{.sequence = 20U, .scale_factor = 1.0};
  migration.solveActiveSet(
      before_x, before_y, before_z, before_mass, empty_accumulator, options, nullptr, nullptr);
  std::vector<double> after_x;
  std::vector<double> after_y;
  std::vector<double> after_z;
  std::vector<double> after_mass;
  if (world_rank == 0) {
    for (int rank = 0; rank < world_size; ++rank) {
      after_x.push_back((static_cast<double>(rank) + 0.5) / static_cast<double>(world_size));
      after_y.push_back(0.2 + 0.1 * static_cast<double>(rank));
      after_z.push_back(0.4);
      after_mass.push_back(1.0);
    }
  }
  options.decomposition_epoch = cosmosim::gravity::DecompositionEpoch{21U};
  options.source_generation = cosmosim::gravity::GravitySourceGeneration{21U};
  options.pm_field_version = cosmosim::gravity::PmFieldVersion{21U};
  options.force_epoch = cosmosim::gravity::ForceEvaluationEpoch{.sequence = 21U, .scale_factor = 1.0};
  migration.solveActiveSet(
      after_x, after_y, after_z, after_mass, empty_accumulator, options, nullptr, nullptr);
}
#endif


void testZoomMembershipPeriodicWrapAffectsCorrection() {
  constexpr double box_size_comoving = 1.0;
  std::vector<double> pos_x = {0.98, 0.02, 0.50};
  std::vector<double> pos_y = {0.50, 0.50, 0.50};
  std::vector<double> pos_z = {0.50, 0.50, 0.50};
  std::vector<double> mass = {1.0, 1.0, 2.0};
  std::vector<std::uint8_t> source_is_high_res = {1, 1, 0};
  std::vector<std::uint8_t> active_is_high_res = source_is_high_res;

  cosmosim::gravity::TreePmOptions opts;
  opts.pm_options.box_size_mpc_comoving = box_size_comoving;
  opts.pm_options.box_size_x_mpc_comoving = box_size_comoving;
  opts.pm_options.box_size_y_mpc_comoving = box_size_comoving;
  opts.pm_options.box_size_z_mpc_comoving = box_size_comoving;
  opts.pm_options.scale_factor = 1.0;
  opts.pm_options.gravitational_constant_code = 1.0;
  opts.tree_options.opening_theta = 0.55;
  opts.tree_options.max_leaf_size = 4;
  opts.tree_options.gravitational_constant_code = 1.0;
  opts.tree_options.softening.epsilon_comoving = 1.0e-3;
  const cosmosim::gravity::PmGridShape coarse_shape = pmShapeForAvailableBackend(16, 8);
  const double mesh_spacing = box_size_comoving / static_cast<double>(coarse_shape.nx);
  opts.split_policy = cosmosim::gravity::makeTreePmSplitPolicyFromMeshSpacing(
      1.25, rcutCellsForAvailableBackend(4.5, 3.9), mesh_spacing);
  opts.enable_zoom_long_range_correction = true;
  opts.zoom_focused_pm_shape = pmShapeForAvailableBackend(24, 8);
  opts.source_is_high_res = source_is_high_res;
  opts.active_is_high_res = active_is_high_res;
  opts.zoom_region_center_x_comoving = 0.0;
  opts.zoom_region_center_y_comoving = 0.5;
  opts.zoom_region_center_z_comoving = 0.5;
  opts.zoom_region_radius_comoving = 0.06;
  opts.zoom_contamination_radius_comoving = 0.1;
  cosmosim::gravity::TreePmDiagnostics diag;
  const ForceField force = solveTreePm(pos_x, pos_y, pos_z, mass, coarse_shape, opts, &diag);
  requireOrThrow(diag.force_l2_pm_zoom_correction > 0.0, "wrapped zoom pair should produce focused correction");
  requireOrThrow(std::abs(force.ax[0]) > 0.0 || std::abs(force.ax[1]) > 0.0, "wrapped high-res targets should feel non-zero force");
}

void testZoomFocusedPmCorrectionAndContaminationDiagnostics() {
  constexpr double box_size_comoving = 1.0;
  std::vector<double> pos_x = {0.49, 0.52, 0.51, 0.20, 0.80};
  std::vector<double> pos_y = {0.50, 0.50, 0.53, 0.20, 0.80};
  std::vector<double> pos_z = {0.50, 0.50, 0.47, 0.20, 0.80};
  std::vector<double> mass = {1.0, 1.1, 0.9, 2.0, 2.5};
  std::vector<std::uint8_t> source_is_high_res = {1, 1, 1, 0, 0};
  std::vector<std::uint8_t> active_is_high_res = source_is_high_res;

  cosmosim::gravity::TreePmOptions no_zoom;
  no_zoom.pm_options.box_size_mpc_comoving = box_size_comoving;
  no_zoom.pm_options.scale_factor = 1.0;
  no_zoom.pm_options.gravitational_constant_code = 1.0;
  no_zoom.tree_options.opening_theta = 0.55;
  no_zoom.tree_options.max_leaf_size = 4;
  no_zoom.tree_options.gravitational_constant_code = 1.0;
  no_zoom.tree_options.softening.epsilon_comoving = 1.0e-3;
  const cosmosim::gravity::PmGridShape coarse_shape = pmShapeForAvailableBackend(16, 8);
  const double mesh_spacing = box_size_comoving / static_cast<double>(coarse_shape.nx);
  no_zoom.split_policy = cosmosim::gravity::makeTreePmSplitPolicyFromMeshSpacing(
      1.25, rcutCellsForAvailableBackend(4.5, 3.9), mesh_spacing);
  cosmosim::gravity::TreePmDiagnostics no_zoom_diag;
  const ForceField no_zoom_force =
      solveTreePm(pos_x, pos_y, pos_z, mass, coarse_shape, no_zoom, &no_zoom_diag);

  cosmosim::gravity::TreePmOptions with_zoom = no_zoom;
  with_zoom.enable_zoom_long_range_correction = true;
  with_zoom.zoom_focused_pm_shape = pmShapeForAvailableBackend(32, 8);
  with_zoom.source_is_high_res = source_is_high_res;
  with_zoom.active_is_high_res = active_is_high_res;
  with_zoom.zoom_region_center_x_comoving = 0.5;
  with_zoom.zoom_region_center_y_comoving = 0.5;
  with_zoom.zoom_region_center_z_comoving = 0.5;
  with_zoom.zoom_region_radius_comoving = 0.08;
  with_zoom.zoom_contamination_radius_comoving = 0.55;
  cosmosim::gravity::TreePmDiagnostics with_zoom_diag;
  const ForceField with_zoom_force =
      solveTreePm(pos_x, pos_y, pos_z, mass, coarse_shape, with_zoom, &with_zoom_diag);

  const double highres_delta = std::abs(with_zoom_force.ax[0] - no_zoom_force.ax[0]) +
      std::abs(with_zoom_force.ay[1] - no_zoom_force.ay[1]) +
      std::abs(with_zoom_force.az[2] - no_zoom_force.az[2]);
  const double lowres_delta = std::abs(with_zoom_force.ax[3] - no_zoom_force.ax[3]) +
      std::abs(with_zoom_force.ay[4] - no_zoom_force.ay[4]);
  requireOrThrow(highres_delta > 1.0e-6, "zoom PM correction should change high-res target forces");
  requireOrThrow(lowres_delta < 1.0e-12, "zoom PM correction must not modify low-res target forces");
  requireOrThrow(with_zoom_diag.zoom_high_res_source_count == 3, "unexpected zoom high-res source count");
  requireOrThrow(with_zoom_diag.zoom_low_res_source_count == 2, "unexpected zoom low-res source count");
  requireOrThrow(with_zoom_diag.zoom_low_res_contamination_count >= 1, "expected at least one low-res contaminant");
  requireOrThrow(with_zoom_diag.force_l2_pm_zoom_correction > 0.0, "zoom PM correction norm should be positive");
  requireOrThrow(std::isfinite(with_zoom_diag.force_l2_pm_global), "pm-global force norm must be finite");
  requireOrThrow(std::isfinite(with_zoom_diag.force_l2_tree_short_range), "tree-short force norm must be finite");
  requireOrThrow(std::isfinite(with_zoom_diag.force_l2_total), "total force norm must be finite");
}

void testBoundedLeafAccelerationReuseAndFallback() {
#if COSMOSIM_ENABLE_MPI
  int world_size=1; MPI_Comm_size(MPI_COMM_WORLD,&world_size);
  if (world_size != 1) return;
#endif
  using namespace cosmosim::gravity;
  constexpr std::size_t n=35U; // two vector blocks and a partial tail
  std::vector<double> x(n),y(n),z(n),mass(n,1.0);
  std::vector<TreeLocalIndex> active(n);
  for (std::size_t i=0;i<n;++i) {
    active[i]=static_cast<TreeLocalIndex>(i);
    x[i]=std::fmod(.173*i,1.0); y[i]=std::fmod(.257*i,2.0); z[i]=std::fmod(.317*i,3.0);
  }
  x[0]=-0.0; // canonical signed zero must not corrupt periodic radix ordering
  x[1]=x[0];y[1]=y[0];z[1]=z[0]; // distinct identities, coincident geometry
  TreePmOptions options;
  options.pm_options.box_size_x_mpc_comoving=1.0;
  options.pm_options.box_size_y_mpc_comoving=2.0;
  options.pm_options.box_size_z_mpc_comoving=3.0;
  options.tree_options.softening.epsilon_comoving=.01;
  options.tree_options.max_leaf_size=64U;
  options.tree_options.multipole_order=TreeMultipoleOrder::kMonopole;
  options.split_policy=makeTreePmSplitPolicyFromMeshSpacing(1.25,3.0,.125);
  options.short_range_only=true;
  options.source_generation=GravitySourceGeneration{1U}; options.source_layout_generation=1U;
  TreePmDiagnostics reference_d;
  const auto reference=solveTreePm(x,y,z,mass,{8,8,8},options,&reference_d,nullptr,{},ForceEvaluation::kShortRangeOnly);
  for (const bool table : {false,true}) {
    for (const bool simd : {false,true}) {
      options.full_residual_lookup_enabled=table; options.leaf_simd_enabled=simd;
      TreePmDiagnostics d;
      const auto value=solveTreePm(x,y,z,mass,{8,8,8},options,&d,nullptr,{},ForceEvaluation::kShortRangeOnly);
      requireEquivalentForces(value,reference,1e-10,"bounded homogeneous leaf kernel");
      requireOrThrow(d.local_traversal.direct_pair_evaluations == reference_d.local_traversal.direct_pair_evaluations,
          "leaf masks changed cutoff/self interaction set");
      if (table) requireOrThrow(d.local_traversal.residual_lookup_targets == n,"homogeneous table never reached production");
    }
  }
  options.full_residual_lookup_enabled=true; options.leaf_simd_enabled=true;
  cosmosim::core::MemoryGovernor limited({.hard_limit_bytes=128U*1024U});
  {
    TreePmCoordinator solver({8,8,8},cosmosim::parallel::makePmSlabLayout(8,8,8,1,0),
        cosmosim::parallel::MpiContext(false,1,0),&limited);
    ForceField value{std::vector<double>(n),std::vector<double>(n),std::vector<double>(n)};
    TreePmForceAccumulatorView acc{active,value.ax,value.ay,value.az};
    TreePmDiagnostics d;TreePmProfileEvent profile;
    profile = {};
    solveForEvaluation(solver,x,y,z,mass,acc,options,ForceEvaluation::kShortRangeOnly,&profile,&d);
    requireEquivalentForces(value,reference,1e-12,"governed table fallback");
    requireOrThrow(d.local_traversal.residual_lookup_targets==0U,"denied table still reported active");
    requireOrThrow(limited.snapshot().rejection_count>0U,"memory fallback was not exercised");
    for (const auto& entry : solver.memoryReport().entries) {
      if (entry.label=="treepm.worker_leaf_scratch") {
        requireOrThrow(entry.governed_commitment,"leaf scratch physical commitment missing");
        requireOrThrow(entry.owned_capacity_bytes%sizeof(TreePmLeafScratch)==0U,"leaf scratch capacity charge wrong");
      }
    }
    options.identical_source_tree_reuse_enabled=true;
    profile = {};
    solveForEvaluation(solver,x,y,z,mass,acc,options,ForceEvaluation::kShortRangeOnly,&profile,&d);
    requireOrThrow(profile.tree_profile.rebuild_reason==5U && profile.tree_profile.build_count==0U,
        "identical sources did not reuse installed tree");
    options.identical_source_tree_reuse_enabled=false; options.topology_refit_enabled=true;
    options.source_generation=GravitySourceGeneration{2U};
    profile = {};
    solveForEvaluation(solver,x,y,z,mass,acc,options,ForceEvaluation::kShortRangeOnly,&profile,&d);
    requireOrThrow(profile.tree_profile.refit_attempted==1U && profile.tree_profile.refit_accepted==1U,
        "generation-advanced legal leaf membership did not refit");
    requireEquivalentForces(value,reference,1e-12,"refit must refresh source geometry");
    options.source_layout_generation=2U;
    profile = {};
    solveForEvaluation(solver,x,y,z,mass,acc,options,ForceEvaluation::kShortRangeOnly,&profile,&d);
    requireOrThrow(profile.tree_profile.build_count==1U && profile.tree_profile.rebuild_reason==8U,
        "changed dense source layout did not invalidate tree reuse/refit");
  }
  requireOrThrow(limited.snapshot().committed_bytes==0U && limited.snapshot().reserved_bytes==0U,
      "optional leaf/table ownership outlived coordinator");
  cosmosim::core::MemoryGovernor roomy({.hard_limit_bytes=2U*1024U*1024U});
  {
    TreePmCoordinator solver({8,8,8},cosmosim::parallel::makePmSlabLayout(8,8,8,1,0),
        cosmosim::parallel::MpiContext(false,1,0),&roomy);
    ForceField value{std::vector<double>(n),std::vector<double>(n),std::vector<double>(n)};
    TreePmForceAccumulatorView acc{active,value.ax,value.ay,value.az};
    solveForEvaluation(solver,x,y,z,mass,acc,options,ForceEvaluation::kShortRangeOnly);
    requireOrThrow(roomy.snapshot().committed_bytes>=sizeof(TreePmResidualCoefficientTable),
        "table admission did not charge fixed physical storage");
    solver.trimOptionalAccelerationCaches(2U*1024U*1024U);
    for (const auto& entry : solver.memoryReport().entries)
      requireOrThrow(entry.label!="treepm.full_residual_table" && entry.label!="treepm.worker_leaf_scratch",
          "mandatory phase could not evict optional caches");
    requireEquivalentForces(value,reference,1e-10,"governed admitted accelerator");
  }
  requireOrThrow(roomy.snapshot().committed_bytes==0U,"cache eviction retained a commitment");
}

// 20261025 optimization campaign: authored, NOT RUN — deferred to maintainer.
// Group walk with per-target source-leaf culling and decoupled target groups.
// With a vanishing maximum opening angle no node passes the MAC, so both the
// group walk and the per-target walk reduce to the exact Gaussian-table direct
// residual within r_cut. Clusters straddle the periodic seams and a corner, so
// culling must use minimum-image leaf geometry. Forces agree to roundoff
// (summation order differs), the in-cutoff pair sets are identical, and the
// culling must reject some leaves.
void testGroupWalkLeafCullingMatchesPerTargetWalk() {
#if COSMOSIM_ENABLE_MPI
  int world_size=1; MPI_Comm_size(MPI_COMM_WORLD,&world_size);
  if (world_size != 1) return;
#endif
  using namespace cosmosim::gravity;
  constexpr std::size_t clusters=4U, per_cluster=96U, n=clusters*per_cluster;
  const std::array<std::array<double,3>,clusters> centers{{
      {0.01,0.5,0.5}, {0.995,0.99,0.005}, {0.5,0.5,0.5}, {0.3,0.7,0.2}}};
  std::vector<double> x(n),y(n),z(n),mass(n,1.0);
  std::vector<TreeLocalIndex> active(n);
  std::uint64_t state=0x9e3779b97f4a7c15ULL;
  const auto uniform=[&state]() {
    state=state*6364136223846793005ULL+1442695040888963407ULL;
    return static_cast<double>(state>>11U)*(1.0/9007199254740992.0);
  };
  const auto wrap=[](double v) { v-=std::floor(v); return v>=1.0 ? 0.0 : v; };
  for (std::size_t i=0;i<n;++i) {
    const auto& c=centers[i/per_cluster];
    active[i]=static_cast<TreeLocalIndex>(i);
    // Alternate a dense core (groups fill to the requested size) with a
    // halo wider than r_cut (so some listed leaves are culled per target).
    const double spread=(i%2U==0U) ? 0.02 : 0.3;
    x[i]=wrap(c[0]+spread*(uniform()-0.5));
    y[i]=wrap(c[1]+spread*(uniform()-0.5));
    z[i]=wrap(c[2]+spread*(uniform()-0.5));
  }
  TreePmOptions options;
  options.pm_options.box_size_x_mpc_comoving=1.0;
  options.pm_options.box_size_y_mpc_comoving=1.0;
  options.pm_options.box_size_z_mpc_comoving=1.0;
  options.tree_options.softening.epsilon_comoving=.002;
  options.tree_options.max_leaf_size=8U;
  options.tree_options.relative_force_tolerance=1.0e-3;
  options.split_policy=makeTreePmSplitPolicyFromMeshSpacing(1.25,6.25,.03125);
  options.short_range_only=true;
  options.acceptance_policy=TreePmAcceptancePolicy::kAdaptiveRelative;
  options.adaptive_error_estimator=TreePmAdaptiveErrorEstimator::kRelativeMonopole;
  options.adaptive_maximum_opening_angle=1.0e-6;
  options.gaussian_pair_lookup_enabled=true;
  options.source_generation=GravitySourceGeneration{1U}; options.source_layout_generation=1U;
  const std::vector<double> previous(n,1.0);
  const PmGridShape shape{32,32,32};
  TreePmDiagnostics reference_d;
  const auto reference=solveTreePm(x,y,z,mass,shape,options,&reference_d,nullptr,previous,ForceEvaluation::kShortRangeOnly);
  requireOrThrow(reference_d.group_walk_status==static_cast<std::uint64_t>(TreePmGroupWalkStatus::kNotRequested),
      "per-target reference reported a group walk");
  options.group_walk_enabled=true;
  for (const std::size_t group_size : {std::size_t{0U}, std::size_t{32U}}) {
    options.group_walk_target_group_size=group_size;
    TreePmDiagnostics d;
    const auto grouped=solveTreePm(x,y,z,mass,shape,options,&d,nullptr,previous,ForceEvaluation::kShortRangeOnly);
    const auto& c=d.local_traversal;
    requireOrThrow(d.group_walk_requested==1U &&
        d.group_walk_status==static_cast<std::uint64_t>(TreePmGroupWalkStatus::kActive), "group walk not active");
    requireOrThrow(c.group_walk_grouped_targets==n && c.group_walk_fallback_targets==0U,
        "eligible targets fell back to the per-target walk");
    requireOrThrow(c.group_walk_culled_leaves>0U && c.group_walk_culled_leaves<c.group_walk_candidate_leaves,
        "per-target leaf culling never rejected (or rejected every) source leaf");
    requireOrThrow(c.group_walk_candidate_particles>=c.direct_pair_evaluations,
        "culling dropped a source the pair kernel evaluated");
    requireOrThrow(c.direct_pair_evaluations==reference_d.local_traversal.direct_pair_evaluations,
        "leaf culling changed the in-cutoff direct pair set");
    if (group_size!=0U) {
      requireOrThrow(c.group_walk_target_group_size_max<=group_size,"target group exceeded the configured size");
      requireOrThrow(c.group_walk_target_group_size_max>options.tree_options.max_leaf_size,
          "target grouping is still bounded by the source leaf size");
    }
    requireEquivalentForces(grouped,reference,1.0e-9,"culled group walk versus per-target walk");
  }
  // Missing acceleration history forces a counted per-target fallback with
  // unchanged forces.
  options.group_walk_target_group_size=32U;
  TreePmDiagnostics fallback_d;
  const auto fallback=solveTreePm(x,y,z,mass,shape,options,&fallback_d,nullptr,{},ForceEvaluation::kShortRangeOnly);
  requireOrThrow(fallback_d.local_traversal.group_walk_fallback_history==n &&
      fallback_d.local_traversal.group_walk_fallback_targets==n &&
      fallback_d.local_traversal.group_walk_groups==0U, "history fallback was not attributed");
  TreePmDiagnostics geometric_d;
  options.group_walk_enabled=false;
  const auto geometric=solveTreePm(x,y,z,mass,shape,options,&geometric_d,nullptr,{},ForceEvaluation::kShortRangeOnly);
  requireEquivalentForces(fallback,geometric,0.0,"history fallback must be the per-target walk");
}

}  // namespace

int main() {
#if COSMOSIM_ENABLE_MPI
  MPI_Init(nullptr, nullptr);
#endif
  try {
  testPmPlanRetainedStorageMatchesAllocation();
  testPmFiniteDifferenceGradientBackend();
  testShortRangeSurfaceRejectsPmRequests();
  testBoundedLeafAccelerationReuseAndFallback();
  testGroupWalkLeafCullingMatchesPerTargetWalk();
  testForceComparatorRejectsNonfiniteComponents();
  testCoincidentPairVectorIsZeroWithFiniteSoftening();
  testRecoveryDirectPairsCutoffThreadsAndMemory();
  testPeriodicTreePmAgainstDirectReference();
  testPeriodicTreeGeometryIsSeamSafeAndTranslationInvariant();
  testRelativeForceErrorMacUsesOwnerHistoryAcrossPeriodicSeam();
  testIndependentTargetAndEmptySourceContractSingleRank();
  testResidualCutoffPrunesTreePath();
  testPmOnlyAndTreePmConsistency();
  testCosmologicalTreePmNormalizationAndSplitComplementarity();
  testTreePmRejectsMismatchedGravitationalConstants();
  testPeriodicTreePmRejectsAmbiguousMinimumImageCutoff();
  testZoomMembershipPeriodicWrapAffectsCorrection();
  testZoomFocusedPmCorrectionAndContaminationDiagnostics();
  testActiveSubsetMatchesFullSolveSingleRank();
  testLongRangeCadenceCacheIsFailClosedAndOwnershipCompatible();
#if COSMOSIM_ENABLE_MPI
  testRankLocalSerialTreePmInsideMpiWorldSkipsDistributedArena();
  testDivergentTreePmLayoutFailsCoordinately();
  testDistributedShortRangeExportImportMatchesSingleRankReference();
  testDistributedActiveSubsetMatchesSingleRankReference();
  testDistributedEmptyRanksIndependentTargetsAndMigration();
#endif
  } catch (const std::exception& error) {
    std::cerr << "periodic TreePM regression failed: " << error.what() << '\n';
#if COSMOSIM_ENABLE_MPI
    MPI_Abort(MPI_COMM_WORLD, 1);
#endif
    return 1;
  }
#if COSMOSIM_ENABLE_MPI
  MPI_Finalize();
#endif
  return 0;
}
