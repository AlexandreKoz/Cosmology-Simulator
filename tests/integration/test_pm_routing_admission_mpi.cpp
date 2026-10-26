// Distributed PM routing admission and workload-plan validity (2026-10-26).
//
// Owner-executed MPI regression for the rank-consistent planning transaction
// of PmSolver::planDensityRouting / planForceInterpolationRouting and the
// re-measured plan check in assignDensity / interpolateForces. Every case is
// written so that a rank-divergent collective shows up as a hang; a watchdog
// thread aborts the job with a named case after k_watchdog_seconds so the
// failure is attributable even without a ctest TIMEOUT.
//
// Intended launches: 2 and 4 ranks (see CMakeLists.txt). Single-rank runs
// exercise only the serial validation paths.

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <mutex>
#include <span>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/gravity/pm_solver.hpp"
#include "cosmosim/parallel/distributed_memory.hpp"

#if COSMOSIM_ENABLE_MPI
#include <mpi.h>
#endif

namespace {

#if COSMOSIM_ENABLE_MPI

using cosmosim::gravity::PmGridShape;
using cosmosim::gravity::PmGridStorage;
using cosmosim::gravity::PmSolveOptions;
using cosmosim::gravity::PmSolver;

constexpr int k_watchdog_seconds = 120;
// MPI_Abort from the watchdog thread is only legal under MPI_THREAD_MULTIPLE;
// otherwise the watchdog calls std::abort, which mpiexec propagates.
bool g_mpi_thread_multiple = false;

void requireOrThrow(bool condition, const std::string& message) {
  if (!condition) throw std::runtime_error(message);
}

// Aborts the whole job if a case does not finish (deadlock detector).
class CaseWatchdog {
 public:
  explicit CaseWatchdog(std::string name) : m_name(std::move(name)), m_thread([this]() { watch(); }) {}
  ~CaseWatchdog() {
    {
      std::lock_guard<std::mutex> lock(m_mutex);
      m_done = true;
    }
    m_cv.notify_all();
    m_thread.join();
  }
  CaseWatchdog(const CaseWatchdog&) = delete;
  CaseWatchdog& operator=(const CaseWatchdog&) = delete;

 private:
  void watch() {
    std::unique_lock<std::mutex> lock(m_mutex);
    if (!m_cv.wait_for(lock, std::chrono::seconds(k_watchdog_seconds), [this]() { return m_done; })) {
      std::fprintf(stderr, "PM routing admission case '%s' exceeded %d s: probable rank-divergent collective\n",
                   m_name.c_str(), k_watchdog_seconds);
      std::fflush(stderr);
      if (g_mpi_thread_multiple) MPI_Abort(MPI_COMM_WORLD, 3);
      std::abort();
    }
  }
  std::string m_name;
  std::mutex m_mutex;
  std::condition_variable m_cv;
  bool m_done = false;
  std::thread m_thread;
};

struct World {
  int size = 1;
  int rank = 0;
};

[[nodiscard]] World world() {
  World w;
  MPI_Comm_size(MPI_COMM_WORLD, &w.size);
  MPI_Comm_rank(MPI_COMM_WORLD, &w.rank);
  return w;
}

// Every rank must observe the same outcome: all throw, or none.
// Returns the local exception message ("" when none).
template <typename Call>
std::string requireCollectiveRejection(Call&& call, const std::string& message) {
  std::string what;
  int local_rejected = 0;
  try {
    call();
  } catch (const std::exception& error) {
    local_rejected = 1;
    what = error.what();
  }
  int all_rejected = 0;
  int any_rejected = 0;
  MPI_Allreduce(&local_rejected, &all_rejected, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&local_rejected, &any_rejected, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  requireOrThrow(all_rejected == 1 && any_rejected == 1, message);
  return what;
}

template <typename Call>
void requireCollectiveSuccess(Call&& call, const std::string& message) {
  int local_failed = 0;
  std::string what;
  try {
    call();
  } catch (const std::exception& error) {
    local_failed = 1;
    what = error.what();
  }
  int any_failed = 0;
  MPI_Allreduce(&local_failed, &any_failed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  requireOrThrow(any_failed == 0, message + (what.empty() ? std::string() : ": " + what));
}

struct Particles {
  std::vector<double> x, y, z, mass;
  [[nodiscard]] std::size_t size() const noexcept { return x.size(); }
  void push(double px, double py, double pz, double m) {
    x.push_back(px); y.push_back(py); z.push_back(pz); mass.push_back(m);
  }
};

// Deterministic global population identical on every rank.
[[nodiscard]] Particles globalPopulation(std::size_t count, double box, std::uint64_t seed) {
  Particles p;
  std::uint64_t state = seed;
  const auto next = [&state]() {
    state = state * 6364136223846793005ULL + 1442695040888963407ULL;
    return static_cast<double>(state >> 11U) * (1.0 / 9007199254740992.0);
  };
  for (std::size_t i = 0; i < count; ++i) {
    p.push(next() * box, next() * box, next() * box, 0.5 + next());
  }
  return p;
}

enum class Distribution { kRoundRobin, kAllOnRankZero, kNoneOnLastRank };

[[nodiscard]] Particles localShare(const Particles& all, const World& w, Distribution distribution) {
  Particles local;
  for (std::size_t i = 0; i < all.size(); ++i) {
    bool mine = false;
    switch (distribution) {
      case Distribution::kRoundRobin: mine = static_cast<int>(i % static_cast<std::size_t>(w.size)) == w.rank; break;
      case Distribution::kAllOnRankZero: mine = w.rank == 0; break;
      case Distribution::kNoneOnLastRank:
        mine = w.rank != w.size - 1 &&
            static_cast<int>(i % static_cast<std::size_t>(std::max(w.size - 1, 1))) == w.rank;
        break;
    }
    if (mine) local.push(all.x[i], all.y[i], all.z[i], all.mass[i]);
  }
  return local;
}

[[nodiscard]] PmSolveOptions periodicOptions(double box) {
  PmSolveOptions options;
  options.box_size_mpc_comoving = box;
  return options;
}

// Distributed deposition (planned) must match the serial reference and the
// measured reservations must bound the observed wire high-water marks.
void checkPlannedDensityMatchesReference(
    const PmGridShape& shape, const Particles& all, Distribution distribution, const std::string& label) {
  const World w = world();
  const double box = 1.0;
  const PmSolveOptions options = periodicOptions(box);
  const auto layout = cosmosim::parallel::makePmSlabLayout(shape.nx, shape.ny, shape.nz, w.size, w.rank);
  PmGridStorage grid(shape, layout);
  PmSolver solver(shape);
  const Particles local = localShare(all, w, distribution);
  PmSolver::PmRoutingWorkloadPlan plan{};
  cosmosim::gravity::PmProfileEvent profile{};
  requireCollectiveSuccess([&]() {
    plan = solver.planDensityRouting(grid, local.x, local.y, local.z, local.mass, options);
    solver.assignDensity(grid, local.x, local.y, local.z, local.mass, options, &profile);
  }, label + ": planned distributed deposition failed");
  if (w.size > 1 && !grid.ownsFullDomain()) {
    requireOrThrow(plan.distributed, label + ": distributed plan expected");
    requireOrThrow(profile.routed_send_buffer_high_water_bytes <= plan.send_reservation_bytes,
                   label + ": send high-water exceeds the planned reservation");
    requireOrThrow(profile.routed_receive_buffer_high_water_bytes <= plan.receive_reservation_bytes,
                   label + ": receive high-water exceeds the planned reservation");
    requireOrThrow(plan.send_reservation_bytes + plan.receive_reservation_bytes <= plan.configured_ceiling_bytes,
                   label + ": workload reservation exceeds the configured ceiling");
  }

  PmGridStorage reference_grid(shape);
  PmSolver reference_solver(shape);
  reference_solver.assignDensity(reference_grid, all.x, all.y, all.z, all.mass, options, nullptr);
  for (std::size_t local_ix = 0; local_ix < layout.local_nx(); ++local_ix) {
    const std::size_t ix = layout.globalXFromLocal(local_ix);
    for (std::size_t iy = 0; iy < shape.ny; ++iy) {
      for (std::size_t iz = 0; iz < shape.nz; ++iz) {
        const double got = grid.density()[grid.linearIndex(ix, iy, iz)];
        const double want = reference_grid.density()[reference_grid.linearIndex(ix, iy, iz)];
        requireOrThrow(std::abs(got - want) <= 1.0e-12 * std::max(1.0, std::abs(want)),
                       label + ": distributed density differs from serial reference");
      }
    }
  }
}

// Planned view-based interpolation must match the serial reference.
void checkPlannedInterpolationMatchesReference(
    const PmGridShape& shape, const Particles& all, Distribution distribution, const std::string& label) {
  const World w = world();
  const PmSolveOptions options = periodicOptions(1.0);
  const auto layout = cosmosim::parallel::makePmSlabLayout(shape.nx, shape.ny, shape.nz, w.size, w.rank);
  PmGridStorage grid(shape, layout);
  PmGridStorage reference_grid(shape);
  PmSolver solver(shape);
  PmSolver reference_solver(shape);
  for (std::size_t ix = 0; ix < shape.nx; ++ix) {
    for (std::size_t iy = 0; iy < shape.ny; ++iy) {
      for (std::size_t iz = 0; iz < shape.nz; ++iz) {
        const double value = 1.0 + static_cast<double>(ix + 10U * iy + 100U * iz);
        const std::size_t r = reference_grid.linearIndex(ix, iy, iz);
        reference_grid.force_x()[r] = value;
        reference_grid.force_y()[r] = -2.0 * value;
        reference_grid.force_z()[r] = 0.5 * value;
        if (layout.ownsGlobalX(ix)) {
          const std::size_t l = grid.linearIndex(ix, iy, iz);
          grid.force_x()[l] = value;
          grid.force_y()[l] = -2.0 * value;
          grid.force_z()[l] = 0.5 * value;
        }
      }
    }
  }
  const Particles local = localShare(all, w, distribution);
  std::vector<std::uint32_t> active(local.size());
  for (std::size_t i = 0; i < active.size(); ++i) active[i] = static_cast<std::uint32_t>(i);
  std::vector<double> ax(local.size()), ay(local.size()), az(local.size());
  const PmSolver::PmForceTargetView view{
      .active_particle_index = active,
      .pos_x_comoving = local.x,
      .pos_y_comoving = local.y,
      .pos_z_comoving = local.z,
      .accel_x_comoving = ax,
      .accel_y_comoving = ay,
      .accel_z_comoving = az,
  };
  requireCollectiveSuccess([&]() {
    (void)solver.planForceInterpolationRouting(grid, view, options);
    solver.interpolateForces(grid, view, options, nullptr);
  }, label + ": planned distributed interpolation failed");
  std::vector<double> rx(local.size()), ry(local.size()), rz(local.size());
  reference_solver.interpolateForces(reference_grid, local.x, local.y, local.z, rx, ry, rz, options, nullptr);
  for (std::size_t i = 0; i < local.size(); ++i) {
    requireOrThrow(std::abs(ax[i] - rx[i]) <= 1.0e-12 * std::max(1.0, std::abs(rx[i])) &&
                   std::abs(ay[i] - ry[i]) <= 1.0e-12 * std::max(1.0, std::abs(ry[i])) &&
                   std::abs(az[i] - rz[i]) <= 1.0e-12 * std::max(1.0, std::abs(rz[i])),
                   label + ": distributed interpolation differs from serial reference");
  }
}

void testInvalidCoordinateExtentsOnOneRank() {
  const CaseWatchdog watchdog("invalid_coordinate_extents_on_one_rank");
  const World w = world();
  const PmGridShape shape{8, 4, 4};
  const auto layout = cosmosim::parallel::makePmSlabLayout(shape.nx, shape.ny, shape.nz, w.size, w.rank);
  PmGridStorage grid(shape, layout);
  PmSolver solver(shape);
  const PmSolveOptions options = periodicOptions(1.0);
  const std::vector<double> two{0.1, 0.6};
  const std::vector<double> one{0.1};
  const int bad_rank = w.size - 1;
  const std::string what = requireCollectiveRejection([&]() {
    (void)solver.planDensityRouting(grid, two, w.rank == bad_rank ? std::span<const double>(one) : two, two, two,
                                    options);
  }, "invalid coordinate extents on one rank were not rejected on every rank");
  if (w.size > 1) {
    if (w.rank == bad_rank) {
      requireOrThrow(what.find("source extents must match") != std::string::npos,
                     "failing rank must keep its original diagnostic");
    } else {
      requireOrThrow(what.find("peer rank " + std::to_string(bad_rank)) != std::string::npos,
                     "peer ranks must name the failing rank");
    }
  }
  // The failed transaction must not leave a plan behind: an unplanned routing
  // call now proceeds collectively.
  requireCollectiveSuccess([&]() { solver.assignDensity(grid, two, two, two, two, options, nullptr); },
                           "routing after a failed plan transaction");
}

void testInvalidBoxOnOneRank() {
  const CaseWatchdog watchdog("invalid_box_on_one_rank");
  const World w = world();
  const PmGridShape shape{8, 4, 4};
  const auto layout = cosmosim::parallel::makePmSlabLayout(shape.nx, shape.ny, shape.nz, w.size, w.rank);
  PmGridStorage grid(shape, layout);
  PmSolver solver(shape);
  PmSolveOptions options = periodicOptions(1.0);
  if (w.rank == 0) options.box_size_mpc_comoving = -1.0;
  const std::vector<double> one{0.25};
  (void)requireCollectiveRejection([&]() {
    (void)solver.planDensityRouting(grid, one, one, one, one, options);
  }, "a rank-local invalid box length was not rejected on every rank");
}

void testInvalidSourceIndexOnOneRank() {
  const CaseWatchdog watchdog("invalid_source_index_on_one_rank");
  const World w = world();
  const PmGridShape shape{8, 4, 4};
  const auto layout = cosmosim::parallel::makePmSlabLayout(shape.nx, shape.ny, shape.nz, w.size, w.rank);
  PmGridStorage grid(shape, layout);
  PmSolver solver(shape);
  const PmSolveOptions options = periodicOptions(1.0);
  const std::vector<double> coordinates{0.1, 0.4, 0.9};
  std::vector<double> accel(3, 0.0);
  const std::vector<std::uint32_t> active{0U, 1U};
  const std::vector<cosmosim::gravity::TreeLocalIndex> good_rows{0U, 2U};
  const std::vector<cosmosim::gravity::TreeLocalIndex> bad_rows{0U, 7U};
  const bool bad = w.rank == w.size - 1;
  const PmSolver::PmForceTargetView view{
      .active_particle_index = active,
      .pos_x_comoving = coordinates,
      .pos_y_comoving = coordinates,
      .pos_z_comoving = coordinates,
      .accel_x_comoving = accel,
      .accel_y_comoving = accel,
      .accel_z_comoving = accel,
      .coordinate_layout = PmSolver::PmForceCoordinateLayout::kIndexedSource,
      .output_layout = PmSolver::PmForceOutputLayout::kIndexedGlobal,
      .coordinate_source_index = bad ? std::span<const cosmosim::gravity::TreeLocalIndex>(bad_rows)
                                     : std::span<const cosmosim::gravity::TreeLocalIndex>(good_rows),
  };
  const std::string what = requireCollectiveRejection([&]() {
    (void)solver.planForceInterpolationRouting(grid, view, options);
  }, "an out-of-range source index on one rank was not rejected on every rank");
  if (bad) {
    requireOrThrow(what.find("coordinate_source_index") != std::string::npos,
                   "failing rank must report the out-of-range source index");
  }
}

void testZeroParticleRank() {
  const CaseWatchdog watchdog("zero_particle_rank");
  const Particles all = globalPopulation(97U, 1.0, 11U);
  checkPlannedDensityMatchesReference(PmGridShape{8, 4, 4}, all, Distribution::kNoneOnLastRank,
                                      "zero-particle rank (density)");
  checkPlannedInterpolationMatchesReference(PmGridShape{8, 4, 4}, all, Distribution::kNoneOnLastRank,
                                            "zero-particle rank (interpolation)");
}

void testUnevenSlabOwnership() {
  const CaseWatchdog watchdog("uneven_slab_ownership");
  const World w = world();
  const Particles all = globalPopulation(131U, 1.0, 23U);
  // nx = P + 1: one rank owns two planes; nx = 3 with four ranks leaves a
  // zero-width slab rank.
  const std::size_t nx = static_cast<std::size_t>(w.size) + 1U;
  checkPlannedDensityMatchesReference(PmGridShape{nx, 4, 4}, all, Distribution::kRoundRobin, "uneven slabs");
  checkPlannedInterpolationMatchesReference(PmGridShape{nx, 4, 4}, all, Distribution::kRoundRobin,
                                            "uneven slabs");
  if (w.size >= 4) {
    checkPlannedDensityMatchesReference(PmGridShape{3, 4, 4}, all, Distribution::kRoundRobin,
                                        "zero-width slab rank");
    checkPlannedInterpolationMatchesReference(PmGridShape{3, 4, 4}, all, Distribution::kRoundRobin,
                                              "zero-width slab rank");
  }
}

void testSkewedDistributionOnOneSlab() {
  const CaseWatchdog watchdog("skewed_distribution_on_one_slab");
  // Every particle sits in the last x slab and is owned by rank 0, so rank 0
  // sends everything to one peer while the others send nothing.
  Particles all = globalPopulation(257U, 1.0, 37U);
  for (double& x : all.x) x = 0.97 + 0.02 * (x - std::floor(x));
  checkPlannedDensityMatchesReference(PmGridShape{16, 4, 4}, all, Distribution::kAllOnRankZero,
                                      "skewed (density)");
  checkPlannedInterpolationMatchesReference(PmGridShape{16, 4, 4}, all, Distribution::kAllOnRankZero,
                                            "skewed (interpolation)");
}

void testPlanRejectedAfterCoordinateMutation() {
  const CaseWatchdog watchdog("plan_after_coordinate_mutation");
  const World w = world();
  if (w.size < 2) return;
  if (w.size > 4) return;  // needs >= 4 owned x planes per rank
  const PmGridShape shape{16, 4, 4};
  const auto layout = cosmosim::parallel::makePmSlabLayout(shape.nx, shape.ny, shape.nz, w.size, w.rank);
  PmGridStorage grid(shape, layout);
  PmSolver solver(shape);
  const PmSolveOptions options = periodicOptions(1.0);
  // Plan with every local particle two planes inside this rank's own slab:
  // its CIC stencil stays local under either cell- or node-centred
  // conventions, so the planned remote workload is exactly zero. Then move
  // them all into a peer's slab interior before depositing.
  Particles local;
  const double own_center = (static_cast<double>(layout.owned_x.begin_x) + 2.0) / static_cast<double>(shape.nx);
  for (int i = 0; i < 16; ++i) local.push(own_center, 0.3, 0.3, 1.0);
  requireCollectiveSuccess([&]() {
    (void)solver.planDensityRouting(grid, local.x, local.y, local.z, local.mass, options);
  }, "planning before mutation");
  const auto peer = cosmosim::parallel::makePmSlabLayout(shape.nx, shape.ny, shape.nz, w.size, (w.rank + 1) % w.size);
  const double peer_center = (static_cast<double>(peer.owned_x.begin_x) + 2.0) / static_cast<double>(shape.nx);
  for (double& x : local.x) x = peer_center;
  const std::string what = requireCollectiveRejection([&]() {
    solver.assignDensity(grid, local.x, local.y, local.z, local.mass, options, nullptr);
  }, "a plan invalidated by coordinate mutation was consumed");
  requireOrThrow(what.find("does not bound the current workload") != std::string::npos,
                 "stale-plan rejection must name the plan violation");
  // The plan was withdrawn: the same deposition now succeeds unplanned.
  requireCollectiveSuccess([&]() {
    solver.assignDensity(grid, local.x, local.y, local.z, local.mass, options, nullptr);
  }, "unplanned deposition after a rejected stale plan");
}

void testPlanRejectedAfterRoutingOptionChange() {
  const CaseWatchdog watchdog("plan_after_routing_option_change");
  const World w = world();
  if (w.size < 2) return;
  const PmGridShape shape{8, 4, 4};
  const auto layout = cosmosim::parallel::makePmSlabLayout(shape.nx, shape.ny, shape.nz, w.size, w.rank);
  PmGridStorage grid(shape, layout);
  PmSolver solver(shape);
  const Particles local = localShare(globalPopulation(64U, 1.0, 5U), w, Distribution::kRoundRobin);
  PmSolveOptions options = periodicOptions(1.0);
  requireCollectiveSuccess([&]() {
    (void)solver.planDensityRouting(grid, local.x, local.y, local.z, local.mass, options);
  }, "planning before option change");
  options.routing_exchange_batch_bytes = 4096U;
  (void)requireCollectiveRejection([&]() {
    solver.assignDensity(grid, local.x, local.y, local.z, local.mass, options, nullptr);
  }, "a plan made under different routing options was consumed");
  // Assignment scheme change is also part of the key.
  options = periodicOptions(1.0);
  requireCollectiveSuccess([&]() {
    (void)solver.planDensityRouting(grid, local.x, local.y, local.z, local.mass, options);
  }, "planning before scheme change");
  options.assignment_scheme = cosmosim::gravity::PmAssignmentScheme::kTsc;
  (void)requireCollectiveRejection([&]() {
    solver.assignDensity(grid, local.x, local.y, local.z, local.mass, options, nullptr);
  }, "a plan made under a different assignment scheme was consumed");
}

void testPlanRejectedAfterDecompositionChange() {
  const CaseWatchdog watchdog("plan_after_decomposition_change");
  const World w = world();
  if (w.size < 2) return;
  const PmGridShape shape{8, 4, 4};
  const auto layout = cosmosim::parallel::makePmSlabLayout(shape.nx, shape.ny, shape.nz, w.size, w.rank);
  PmGridStorage grid(shape, layout);
  PmSolver solver(shape);
  const Particles local = localShare(globalPopulation(64U, 1.0, 9U), w, Distribution::kRoundRobin);
  const PmSolveOptions options = periodicOptions(1.0);
  requireCollectiveSuccess([&]() {
    (void)solver.planDensityRouting(grid, local.x, local.y, local.z, local.mass, options);
  }, "planning before decomposition change");
  // Routing against a different grid decomposition (one rank now owns the
  // whole mesh, the others a distributed slab) must be rejected collectively
  // — by layout validation or by the plan key — and never hang.
  const auto other_layout = w.rank == 0
      ? cosmosim::parallel::makePmSlabLayout(shape.nx, shape.ny, shape.nz, 1, 0)
      : layout;
  PmGridStorage other_grid(shape, other_layout);
  (void)requireCollectiveRejection([&]() {
    solver.assignDensity(other_grid, local.x, local.y, local.z, local.mass, options, nullptr);
  }, "a plan made under a different decomposition was consumed");
}

void testPlanPresenceMustBeCollective() {
  const CaseWatchdog watchdog("plan_presence_must_be_collective");
  const World w = world();
  if (w.size < 2) return;
  const PmGridShape shape{8, 4, 4};
  const auto layout = cosmosim::parallel::makePmSlabLayout(shape.nx, shape.ny, shape.nz, w.size, w.rank);
  PmGridStorage grid(shape, layout);
  PmSolver solver(shape);
  PmSolver planning_solver(shape);
  const Particles local = localShare(globalPopulation(32U, 1.0, 3U), w, Distribution::kRoundRobin);
  const PmSolveOptions options = periodicOptions(1.0);
  // Collective planning on a helper solver keeps collectives aligned; only
  // rank 0's routing solver then holds a plan (via a second collective plan
  // whose result is discarded on the other ranks by using planning_solver).
  requireCollectiveSuccess([&]() {
    if (w.rank == 0) {
      (void)solver.planDensityRouting(grid, local.x, local.y, local.z, local.mass, options);
    } else {
      (void)planning_solver.planDensityRouting(grid, local.x, local.y, local.z, local.mass, options);
    }
  }, "split planning");
  const std::string what = requireCollectiveRejection([&]() {
    solver.assignDensity(grid, local.x, local.y, local.z, local.mass, options, nullptr);
  }, "rank-divergent plan presence was not rejected collectively");
  requireOrThrow(what.find("presence differs across ranks") != std::string::npos,
                 "plan presence divergence must be named");
}

#endif

}  // namespace

int main() {
#if COSMOSIM_ENABLE_MPI
  int provided = 0;
  MPI_Init_thread(nullptr, nullptr, MPI_THREAD_MULTIPLE, &provided);
  g_mpi_thread_multiple = provided >= MPI_THREAD_MULTIPLE;
  int failed = 0;
  try {
    testInvalidCoordinateExtentsOnOneRank();
    testInvalidBoxOnOneRank();
    testInvalidSourceIndexOnOneRank();
    testZeroParticleRank();
    testUnevenSlabOwnership();
    testSkewedDistributionOnOneSlab();
    testPlanRejectedAfterCoordinateMutation();
    testPlanRejectedAfterRoutingOptionChange();
    testPlanRejectedAfterDecompositionChange();
    testPlanPresenceMustBeCollective();
  } catch (const std::exception& error) {
    std::fprintf(stderr, "PM routing admission test failed: %s\n", error.what());
    failed = 1;
  }
  int any_failed = 0;
  MPI_Allreduce(&failed, &any_failed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Finalize();
  return any_failed == 0 ? 0 : 1;
#else
  return 0;
#endif
}
