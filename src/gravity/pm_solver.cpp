#include "cosmosim/gravity/pm_solver.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <numeric>
#include <optional>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
#include <limits>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/core/device_buffer.hpp"
#include "cosmosim/gravity/tree_pm_split_kernel.hpp"
#if COSMOSIM_ENABLE_CUDA
#include <cuda_runtime.h>
#include "cosmosim/gravity/pm_cuda_kernels.hpp"
#endif

#if COSMOSIM_ENABLE_FFTW
#include <fftw3.h>
#if COSMOSIM_ENABLE_MPI
#include <fftw3-mpi.h>
#include <mpi.h>
#endif
#endif

namespace cosmosim::gravity {
namespace {

constexpr double k_pi = 3.141592653589793238462643383279502884;

[[nodiscard]] std::size_t wrapIndex(std::ptrdiff_t i, std::size_t n) {
  const std::ptrdiff_t n_signed = static_cast<std::ptrdiff_t>(n);
  std::ptrdiff_t wrapped = i % n_signed;
  if (wrapped < 0) {
    wrapped += n_signed;
  }
  return static_cast<std::size_t>(wrapped);
}

[[nodiscard]] double wrapPosition(double x, double box_size) {
  const double wrapped = std::fmod(x, box_size);
  return wrapped < 0.0 ? wrapped + box_size : wrapped;
}

[[nodiscard]] bool positionInsideOpenDomain(double x, double box_size) {
  return x >= 0.0 && x < box_size;
}

[[nodiscard]] double sinc(double x) {
  if (std::abs(x) < 1.0e-12) {
    return 1.0;
  }
  return std::sin(x) / x;
}

struct PmAxisStencil1d {
  std::array<std::ptrdiff_t, 3> offsets{};
  std::array<double, 3> weights{};
  std::size_t count = 0;
};

struct PmDensityContributionRecord {
  std::uint32_t global_ix = 0;
  std::uint32_t global_iy = 0;
  std::uint32_t global_iz = 0;
  double mass_contribution = 0.0;
};

struct PmInterpolationRequestRecord {
  std::uint32_t particle_index = 0;
  std::uint32_t global_ix = 0;
  std::uint32_t global_iy = 0;
  std::uint32_t global_iz = 0;
  double weight = 0.0;
};

struct PmForceContributionRecord {
  std::uint32_t particle_index = 0;
  double accel_x = 0.0;
  double accel_y = 0.0;
  double accel_z = 0.0;
};

struct PmPotentialContributionRecord {
  std::uint32_t particle_index = 0;
  double potential = 0.0;
};

[[nodiscard]] PmAxisStencil1d makeAxisStencil(double grid_position, PmAssignmentScheme scheme) {
  PmAxisStencil1d stencil{};
  if (scheme == PmAssignmentScheme::kCic) {
    const double base = std::floor(grid_position);
    const std::ptrdiff_t i0 = static_cast<std::ptrdiff_t>(base);
    const double t = grid_position - base;
    stencil.offsets = {i0, i0 + 1, 0};
    stencil.weights = {1.0 - t, t, 0.0};
    stencil.count = 2;
    return stencil;
  }

  const double center = std::floor(grid_position + 0.5);
  const std::ptrdiff_t ic = static_cast<std::ptrdiff_t>(center);
  const double delta = grid_position - center;
  const double w_m1 = 0.5 * std::pow(0.5 - delta, 2.0);
  const double w_0 = 0.75 - delta * delta;
  const double w_p1 = 0.5 * std::pow(0.5 + delta, 2.0);
  stencil.offsets = {ic - 1, ic, ic + 1};
  stencil.weights = {w_m1, w_0, w_p1};
  stencil.count = 3;
  return stencil;
}

[[nodiscard]] int assignmentWindowExponent(PmAssignmentScheme scheme) {
  switch (scheme) {
    case PmAssignmentScheme::kCic:
      return 2;
    case PmAssignmentScheme::kTsc:
      return 3;
  }
  throw std::invalid_argument("Unknown PM assignment scheme in assignmentWindowExponent");
}

[[nodiscard]] std::uint64_t bytesForGridSweep(std::size_t cell_count) {
  return static_cast<std::uint64_t>(cell_count * sizeof(double));
}

[[nodiscard]] std::uint64_t bytesForParticles(std::size_t particle_count) {
  return static_cast<std::uint64_t>(particle_count * sizeof(double) * 4U);
}

struct BoxLengths {
  double lx = 0.0;
  double ly = 0.0;
  double lz = 0.0;
};

[[nodiscard]] BoxLengths effectiveBoxLengths(const PmSolveOptions& options) {
  const double scalar = options.box_size_mpc_comoving;
  const double lx = options.box_size_x_mpc_comoving > 0.0 ? options.box_size_x_mpc_comoving : scalar;
  const double ly = options.box_size_y_mpc_comoving > 0.0 ? options.box_size_y_mpc_comoving : scalar;
  const double lz = options.box_size_z_mpc_comoving > 0.0 ? options.box_size_z_mpc_comoving : scalar;
  return BoxLengths{.lx = lx, .ly = ly, .lz = lz};
}

void validateOptions(const PmGridShape& shape, const PmSolveOptions& options) {
  if (!shape.isValid()) {
    throw std::invalid_argument("PM grid shape must be non-zero in all dimensions");
  }
  const BoxLengths lengths = effectiveBoxLengths(options);
  if (lengths.lx <= 0.0 || lengths.ly <= 0.0 || lengths.lz <= 0.0) {
    throw std::invalid_argument("PM solve requires positive box lengths on all axes");
  }
  if (options.scale_factor <= 0.0) {
    throw std::invalid_argument("PM solve requires scale_factor > 0");
  }
  if (options.gravitational_constant_code <= 0.0) {
    throw std::invalid_argument("PM solve requires gravitational_constant_code > 0");
  }
  if (options.execution_policy == core::ExecutionPolicy::kCuda && options.data_residency == PmDataResidencyPolicy::kHostOnly) {
    throw std::invalid_argument(
        "execution_policy=cuda requires data_residency=kPreferDevice for explicit host/device ownership");
  }
  if (options.execution_policy == core::ExecutionPolicy::kCuda && options.assignment_scheme != PmAssignmentScheme::kCic) {
    throw std::invalid_argument(
        "execution_policy=cuda currently supports only assignment_scheme=cic in this build");
  }
  if (options.execution_policy == core::ExecutionPolicy::kCuda &&
      (std::abs(lengths.lx - lengths.ly) > 1.0e-12 || std::abs(lengths.lx - lengths.lz) > 1.0e-12)) {
    throw std::invalid_argument(
        "execution_policy=cuda currently requires cubic PM box lengths in this build");
  }
  if (options.boundary_condition == PmBoundaryCondition::kIsolatedOpen &&
      options.enable_window_deconvolution) {
    throw std::invalid_argument("isolated PM currently does not support window deconvolution");
  }
}

void validateSingleRankFullDomainGridContract(const PmGridStorage& grid, std::string_view callsite) {
  if (!grid.ownsFullDomain()) {
    throw std::invalid_argument(
        std::string(callsite) +
        " requires full-domain PM slab ownership in single-rank mode; use a valid distributed slab layout for multi-rank PM execution");
  }
}

}  // namespace

class PmSolver::Impl {
 public:
  struct PlanKey {
    int world_size = 1;
    int world_rank = 0;
    std::size_t owned_begin_x = 0;
    std::size_t owned_end_x = 0;
    core::PmDecompositionMode decomposition_mode = core::PmDecompositionMode::kSlab;

    [[nodiscard]] bool operator==(const PlanKey& other) const noexcept {
      return world_size == other.world_size && world_rank == other.world_rank && owned_begin_x == other.owned_begin_x &&
          owned_end_x == other.owned_end_x && decomposition_mode == other.decomposition_mode;
    }
  };

  struct PlanKeyHasher {
    [[nodiscard]] std::size_t operator()(const PlanKey& key) const noexcept {
      std::size_t seed = static_cast<std::size_t>(key.world_size * 1315423911U + key.world_rank);
      seed ^= key.owned_begin_x + 0x9e3779b97f4a7c15ULL + (seed << 6U) + (seed >> 2U);
      seed ^= key.owned_end_x + 0x9e3779b97f4a7c15ULL + (seed << 6U) + (seed >> 2U);
      seed ^= static_cast<std::size_t>(key.decomposition_mode) + 0x9e3779b97f4a7c15ULL + (seed << 6U) + (seed >> 2U);
      return seed;
    }
  };

  struct PlanResources {
    parallel::PmSlabLayout layout{};
    std::vector<double> real;
    std::vector<std::complex<double>> fourier;
    std::vector<std::complex<double>> potential_k;
    std::vector<std::complex<double>> working_k;
#if COSMOSIM_ENABLE_FFTW
    fftw_plan forward_plan = nullptr;
    fftw_plan inverse_plan = nullptr;
#endif
    bool is_distributed = false;
    bool spectral_transposed = false;
    std::size_t transposed_local_ny = 0;
    std::size_t transposed_begin_y = 0;
  };

  struct DensityExchangeBuffers {
    std::vector<std::vector<PmDensityContributionRecord>> send_records_by_rank;
    std::vector<PmDensityContributionRecord> send_flat_records;
    std::vector<PmDensityContributionRecord> recv_flat_records;
    std::vector<int> send_counts;
    std::vector<int> send_displs;
    std::vector<int> recv_counts;
    std::vector<int> recv_displs;
    std::vector<int> send_counts_bytes;
    std::vector<int> send_displs_bytes;
    std::vector<int> recv_counts_bytes;
    std::vector<int> recv_displs_bytes;
  };

  template <typename ContributionRecord>
  struct InterpolationExchangeBuffers {
    std::vector<std::vector<PmInterpolationRequestRecord>> send_requests_by_rank;
    std::vector<PmInterpolationRequestRecord> send_requests_flat;
    std::vector<PmInterpolationRequestRecord> recv_requests_flat;
    std::vector<std::vector<ContributionRecord>> send_contribs_by_rank;
    std::vector<ContributionRecord> send_contribs_flat;
    std::vector<ContributionRecord> recv_contribs_flat;
    std::vector<int> send_counts;
    std::vector<int> send_displs;
    std::vector<int> recv_counts;
    std::vector<int> recv_displs;
    std::vector<int> send_counts_bytes;
    std::vector<int> send_displs_bytes;
    std::vector<int> recv_counts_bytes;
    std::vector<int> recv_displs_bytes;
    std::vector<int> send_contrib_counts;
    std::vector<int> send_contrib_displs;
    std::vector<int> recv_contrib_counts;
    std::vector<int> recv_contrib_displs;
    std::vector<int> send_contrib_counts_bytes;
    std::vector<int> send_contrib_displs_bytes;
    std::vector<int> recv_contrib_counts_bytes;
    std::vector<int> recv_contrib_displs_bytes;
  };

  explicit Impl(PmGridShape shape) : m_shape(shape) {}

  ~Impl() {
#if COSMOSIM_ENABLE_FFTW
    for (auto& [_, plan] : m_plan_cache) {
      if (plan.forward_plan != nullptr) {
        fftw_destroy_plan(plan.forward_plan);
      }
      if (plan.inverse_plan != nullptr) {
        fftw_destroy_plan(plan.inverse_plan);
      }
    }
    if (m_isolated_workspace.forward_plan != nullptr) {
      fftw_destroy_plan(m_isolated_workspace.forward_plan);
    }
    if (m_isolated_workspace.inverse_plan != nullptr) {
      fftw_destroy_plan(m_isolated_workspace.inverse_plan);
    }
#endif
  }

  struct IsolatedOpenWorkspace {
    std::size_t nx = 0;
    std::size_t ny = 0;
    std::size_t nz = 0;
    std::vector<std::complex<double>> rho_k;
    std::vector<std::complex<double>> kernel_k;
    std::vector<std::complex<double>> scratch;
    std::vector<double> potential_real;
    bool kernel_ready = false;
    double dx = 0.0;
    double dy = 0.0;
    double dz = 0.0;
    double split_scale = -1.0;
#if COSMOSIM_ENABLE_FFTW
    fftw_plan forward_plan = nullptr;
    fftw_plan inverse_plan = nullptr;
#endif
  };

  void ensureIsolatedWorkspace(std::size_t nx, std::size_t ny, std::size_t nz) {
    const bool shape_changed = nx != m_isolated_workspace.nx || ny != m_isolated_workspace.ny || nz != m_isolated_workspace.nz;
    if (!shape_changed) {
      return;
    }
#if COSMOSIM_ENABLE_FFTW
    if (m_isolated_workspace.forward_plan != nullptr) {
      fftw_destroy_plan(m_isolated_workspace.forward_plan);
      m_isolated_workspace.forward_plan = nullptr;
    }
    if (m_isolated_workspace.inverse_plan != nullptr) {
      fftw_destroy_plan(m_isolated_workspace.inverse_plan);
      m_isolated_workspace.inverse_plan = nullptr;
    }
#endif
    m_isolated_workspace.nx = nx;
    m_isolated_workspace.ny = ny;
    m_isolated_workspace.nz = nz;
    const std::size_t n_tot = nx * ny * nz;
    m_isolated_workspace.rho_k.assign(n_tot, {0.0, 0.0});
    m_isolated_workspace.kernel_k.assign(n_tot, {0.0, 0.0});
    m_isolated_workspace.scratch.assign(n_tot, {0.0, 0.0});
    m_isolated_workspace.potential_real.assign(n_tot, 0.0);
    m_isolated_workspace.kernel_ready = false;
    m_isolated_workspace.dx = 0.0;
    m_isolated_workspace.dy = 0.0;
    m_isolated_workspace.dz = 0.0;
    m_isolated_workspace.split_scale = -1.0;
#if COSMOSIM_ENABLE_FFTW
    m_isolated_workspace.forward_plan = fftw_plan_dft_3d(
        static_cast<int>(nx),
        static_cast<int>(ny),
        static_cast<int>(nz),
        reinterpret_cast<fftw_complex*>(m_isolated_workspace.scratch.data()),
        reinterpret_cast<fftw_complex*>(m_isolated_workspace.scratch.data()),
        FFTW_FORWARD,
        FFTW_ESTIMATE);
    m_isolated_workspace.inverse_plan = fftw_plan_dft_3d(
        static_cast<int>(nx),
        static_cast<int>(ny),
        static_cast<int>(nz),
        reinterpret_cast<fftw_complex*>(m_isolated_workspace.scratch.data()),
        reinterpret_cast<fftw_complex*>(m_isolated_workspace.scratch.data()),
        FFTW_BACKWARD,
        FFTW_ESTIMATE);
    if (m_isolated_workspace.forward_plan == nullptr || m_isolated_workspace.inverse_plan == nullptr) {
      throw std::runtime_error("failed to create isolated PM FFTW plans");
    }
#endif
  }

  void isolatedForward(std::span<std::complex<double>> field) {
    std::copy(field.begin(), field.end(), m_isolated_workspace.scratch.begin());
#if COSMOSIM_ENABLE_FFTW
    fftw_execute(m_isolated_workspace.forward_plan);
    std::copy(m_isolated_workspace.scratch.begin(), m_isolated_workspace.scratch.end(), field.begin());
#else
    naiveComplexDft(field, /*forward=*/true);
#endif
  }

  void isolatedInverse(std::span<std::complex<double>> field) {
    std::copy(field.begin(), field.end(), m_isolated_workspace.scratch.begin());
#if COSMOSIM_ENABLE_FFTW
    fftw_execute(m_isolated_workspace.inverse_plan);
    std::copy(m_isolated_workspace.scratch.begin(), m_isolated_workspace.scratch.end(), field.begin());
#else
    naiveComplexDft(field, /*forward=*/false);
#endif
  }

  [[nodiscard]] IsolatedOpenWorkspace& isolatedWorkspace() { return m_isolated_workspace; }

  [[nodiscard]] PlanResources& planForLayout(const parallel::PmSlabLayout& layout, core::PmDecompositionMode decomposition_mode) {
    const PlanKey key{
        .world_size = layout.world_size,
        .world_rank = layout.world_rank,
        .owned_begin_x = layout.owned_x.begin_x,
        .owned_end_x = layout.owned_x.end_x,
        .decomposition_mode = decomposition_mode,
    };
    auto it = m_plan_cache.find(key);
    if (it != m_plan_cache.end()) {
      m_active_key = key;
      return it->second;
    }

    PlanResources plan{};
    plan.layout = layout;
    plan.real.assign(layout.localCellCount(), 0.0);
    const std::size_t nz_complex = m_shape.nz / 2U + 1U;
    const std::size_t expected_local_complex_size = layout.local_nx() * m_shape.ny * nz_complex;
    std::size_t allocated_local_complex_size = expected_local_complex_size;

#if COSMOSIM_ENABLE_FFTW
#if COSMOSIM_ENABLE_MPI
    if (layout.world_size > 1) {
      int mpi_world_size = 1;
      int mpi_world_rank = 0;
      MPI_Comm_size(MPI_COMM_WORLD, &mpi_world_size);
      MPI_Comm_rank(MPI_COMM_WORLD, &mpi_world_rank);
      if (mpi_world_size != layout.world_size || mpi_world_rank != layout.world_rank) {
        throw std::invalid_argument("PM slab layout world metadata must match MPI_COMM_WORLD for distributed FFT path");
      }
      fftw_mpi_init();
      ptrdiff_t backend_local_nx = 0;
      ptrdiff_t backend_begin_x = 0;
      ptrdiff_t backend_alloc_local = 0;
      ptrdiff_t backend_local_ny = 0;
      ptrdiff_t backend_begin_y = 0;
      if (decomposition_mode == core::PmDecompositionMode::kPencil) {
        backend_alloc_local = fftw_mpi_local_size_3d_transposed(
            static_cast<ptrdiff_t>(m_shape.nx),
            static_cast<ptrdiff_t>(m_shape.ny),
            static_cast<ptrdiff_t>(nz_complex),
            MPI_COMM_WORLD,
            &backend_local_nx,
            &backend_begin_x,
            &backend_local_ny,
            &backend_begin_y);
      } else {
        backend_alloc_local = fftw_mpi_local_size_3d(
            static_cast<ptrdiff_t>(m_shape.nx),
            static_cast<ptrdiff_t>(m_shape.ny),
            static_cast<ptrdiff_t>(nz_complex),
            MPI_COMM_WORLD,
            &backend_local_nx,
            &backend_begin_x);
      }
      if (backend_local_nx != static_cast<ptrdiff_t>(layout.local_nx()) ||
          backend_begin_x != static_cast<ptrdiff_t>(layout.owned_x.begin_x)) {
        throw std::invalid_argument(
            "PM slab layout is incompatible with FFTW MPI ownership for this communicator; the configured slab partition does not match backend local_nx/local_0_start");
      }
      if (backend_alloc_local <= 0) {
        throw std::runtime_error("FFTW MPI reported non-positive local allocation size for distributed PM plan");
      }
      allocated_local_complex_size = static_cast<std::size_t>(backend_alloc_local);
      plan.is_distributed = true;
      if (decomposition_mode == core::PmDecompositionMode::kPencil) {
        plan.spectral_transposed = true;
        plan.transposed_local_ny = static_cast<std::size_t>(backend_local_ny);
        plan.transposed_begin_y = static_cast<std::size_t>(backend_begin_y);
      }
    }
#endif
    plan.fourier.assign(allocated_local_complex_size, std::complex<double>(0.0, 0.0));
    plan.potential_k.assign(allocated_local_complex_size, std::complex<double>(0.0, 0.0));
    plan.working_k.assign(allocated_local_complex_size, std::complex<double>(0.0, 0.0));
#if COSMOSIM_ENABLE_MPI
    if (layout.world_size > 1) {
      if (allocated_local_complex_size < expected_local_complex_size) {
        throw std::runtime_error("FFTW MPI local allocation is smaller than the expected slab-local Fourier extent");
      }
      plan.forward_plan = fftw_mpi_plan_dft_r2c_3d(
          static_cast<ptrdiff_t>(m_shape.nx),
          static_cast<ptrdiff_t>(m_shape.ny),
          static_cast<ptrdiff_t>(m_shape.nz),
          plan.real.data(),
          reinterpret_cast<fftw_complex*>(plan.fourier.data()),
          MPI_COMM_WORLD,
          decomposition_mode == core::PmDecompositionMode::kPencil
              ? (FFTW_MEASURE | FFTW_MPI_TRANSPOSED_OUT)
              : FFTW_MEASURE);
      plan.inverse_plan = fftw_mpi_plan_dft_c2r_3d(
          static_cast<ptrdiff_t>(m_shape.nx),
          static_cast<ptrdiff_t>(m_shape.ny),
          static_cast<ptrdiff_t>(m_shape.nz),
          reinterpret_cast<fftw_complex*>(plan.fourier.data()),
          plan.real.data(),
          MPI_COMM_WORLD,
          decomposition_mode == core::PmDecompositionMode::kPencil
              ? (FFTW_MEASURE | FFTW_MPI_TRANSPOSED_IN)
              : FFTW_MEASURE);
    } else
#endif
    {
      plan.fourier.assign(expected_local_complex_size, std::complex<double>(0.0, 0.0));
      plan.potential_k.assign(expected_local_complex_size, std::complex<double>(0.0, 0.0));
      plan.working_k.assign(expected_local_complex_size, std::complex<double>(0.0, 0.0));
      plan.forward_plan = fftw_plan_dft_r2c_3d(
          static_cast<int>(m_shape.nx),
          static_cast<int>(m_shape.ny),
          static_cast<int>(m_shape.nz),
          plan.real.data(),
          reinterpret_cast<fftw_complex*>(plan.fourier.data()),
          FFTW_MEASURE);
      plan.inverse_plan = fftw_plan_dft_c2r_3d(
          static_cast<int>(m_shape.nx),
          static_cast<int>(m_shape.ny),
          static_cast<int>(m_shape.nz),
          reinterpret_cast<fftw_complex*>(plan.fourier.data()),
          plan.real.data(),
          FFTW_MEASURE);
    }
    if (plan.forward_plan == nullptr || plan.inverse_plan == nullptr) {
      throw std::runtime_error("Failed to create FFTW plans for PM solver");
    }
#else
    if (!layout.ownsFullDomain()) {
      throw std::invalid_argument(
          "PM solver naive DFT fallback requires full-domain slab ownership; distributed PM requires COSMOSIM_ENABLE_FFTW=ON and COSMOSIM_ENABLE_MPI=ON");
    }
    plan.fourier.assign(expected_local_complex_size, std::complex<double>(0.0, 0.0));
    plan.potential_k.assign(expected_local_complex_size, std::complex<double>(0.0, 0.0));
    plan.working_k.assign(expected_local_complex_size, std::complex<double>(0.0, 0.0));
#endif

    auto [insert_it, inserted] = m_plan_cache.emplace(key, std::move(plan));
    (void)inserted;
    ++m_plan_build_count;
    m_active_key = key;
    return insert_it->second;
  }

  [[nodiscard]] std::span<double> realGrid() { return activePlan().real; }
  [[nodiscard]] std::span<std::complex<double>> fourierGrid() { return activePlan().fourier; }
  [[nodiscard]] std::span<std::complex<double>> potentialScratch() { return activePlan().potential_k; }
  [[nodiscard]] std::span<std::complex<double>> workingScratch() { return activePlan().working_k; }

  double forwardFft() {
    const auto start = std::chrono::steady_clock::now();
#if COSMOSIM_ENABLE_FFTW
    fftw_execute(activePlan().forward_plan);
#else
    naiveForwardDft();
#endif
    const auto stop = std::chrono::steady_clock::now();
    return std::chrono::duration<double, std::milli>(stop - start).count();
  }

  double inverseFft() {
    const auto start = std::chrono::steady_clock::now();
#if COSMOSIM_ENABLE_FFTW
    fftw_execute(activePlan().inverse_plan);
#else
    naiveInverseDft();
#endif
    const auto stop = std::chrono::steady_clock::now();
    return std::chrono::duration<double, std::milli>(stop - start).count();
  }

  [[nodiscard]] std::size_t planCount() const noexcept { return m_plan_cache.size(); }
  [[nodiscard]] std::size_t planBuildCount() const noexcept { return m_plan_build_count; }

  [[nodiscard]] DensityExchangeBuffers& densityExchangeBuffersForLayout(const parallel::PmSlabLayout& layout) {
    if (m_density_exchange.world_size != layout.world_size || m_density_exchange.world_rank != layout.world_rank) {
      m_density_exchange.world_size = layout.world_size;
      m_density_exchange.world_rank = layout.world_rank;
      m_density_exchange.buffers.send_records_by_rank.assign(static_cast<std::size_t>(layout.world_size), {});
      m_density_exchange.buffers.send_counts.assign(static_cast<std::size_t>(layout.world_size), 0);
      m_density_exchange.buffers.send_displs.assign(static_cast<std::size_t>(layout.world_size), 0);
      m_density_exchange.buffers.recv_counts.assign(static_cast<std::size_t>(layout.world_size), 0);
      m_density_exchange.buffers.recv_displs.assign(static_cast<std::size_t>(layout.world_size), 0);
      m_density_exchange.buffers.send_counts_bytes.assign(static_cast<std::size_t>(layout.world_size), 0);
      m_density_exchange.buffers.send_displs_bytes.assign(static_cast<std::size_t>(layout.world_size), 0);
      m_density_exchange.buffers.recv_counts_bytes.assign(static_cast<std::size_t>(layout.world_size), 0);
      m_density_exchange.buffers.recv_displs_bytes.assign(static_cast<std::size_t>(layout.world_size), 0);
    }
    return m_density_exchange.buffers;
  }

  [[nodiscard]] InterpolationExchangeBuffers<PmForceContributionRecord>&
  forceInterpolationExchangeBuffersForLayout(const parallel::PmSlabLayout& layout) {
    if (m_force_exchange.world_size != layout.world_size || m_force_exchange.world_rank != layout.world_rank) {
      m_force_exchange.world_size = layout.world_size;
      m_force_exchange.world_rank = layout.world_rank;
      resetInterpolationExchangeForWorld(layout.world_size, m_force_exchange.buffers);
    }
    return m_force_exchange.buffers;
  }

  [[nodiscard]] InterpolationExchangeBuffers<PmPotentialContributionRecord>&
  potentialInterpolationExchangeBuffersForLayout(const parallel::PmSlabLayout& layout) {
    if (m_potential_exchange.world_size != layout.world_size ||
        m_potential_exchange.world_rank != layout.world_rank) {
      m_potential_exchange.world_size = layout.world_size;
      m_potential_exchange.world_rank = layout.world_rank;
      resetInterpolationExchangeForWorld(layout.world_size, m_potential_exchange.buffers);
    }
    return m_potential_exchange.buffers;
  }

 private:
  template <typename ContributionRecord>
  static void resetInterpolationExchangeForWorld(
      int world_size,
      InterpolationExchangeBuffers<ContributionRecord>& buffers) {
    buffers.send_requests_by_rank.assign(static_cast<std::size_t>(world_size), {});
    buffers.send_contribs_by_rank.assign(static_cast<std::size_t>(world_size), {});
    buffers.send_counts.assign(static_cast<std::size_t>(world_size), 0);
    buffers.send_displs.assign(static_cast<std::size_t>(world_size), 0);
    buffers.recv_counts.assign(static_cast<std::size_t>(world_size), 0);
    buffers.recv_displs.assign(static_cast<std::size_t>(world_size), 0);
    buffers.send_counts_bytes.assign(static_cast<std::size_t>(world_size), 0);
    buffers.send_displs_bytes.assign(static_cast<std::size_t>(world_size), 0);
    buffers.recv_counts_bytes.assign(static_cast<std::size_t>(world_size), 0);
    buffers.recv_displs_bytes.assign(static_cast<std::size_t>(world_size), 0);
    buffers.send_contrib_counts.assign(static_cast<std::size_t>(world_size), 0);
    buffers.send_contrib_displs.assign(static_cast<std::size_t>(world_size), 0);
    buffers.recv_contrib_counts.assign(static_cast<std::size_t>(world_size), 0);
    buffers.recv_contrib_displs.assign(static_cast<std::size_t>(world_size), 0);
    buffers.send_contrib_counts_bytes.assign(static_cast<std::size_t>(world_size), 0);
    buffers.send_contrib_displs_bytes.assign(static_cast<std::size_t>(world_size), 0);
    buffers.recv_contrib_counts_bytes.assign(static_cast<std::size_t>(world_size), 0);
    buffers.recv_contrib_displs_bytes.assign(static_cast<std::size_t>(world_size), 0);
  }

  [[nodiscard]] PlanResources& activePlan() {
    if (!m_active_key.has_value()) {
      throw std::logic_error("PM solver plan has not been initialized for the active slab layout");
    }
    return m_plan_cache.at(*m_active_key);
  }

#if !COSMOSIM_ENABLE_FFTW
  void naiveComplexDft(std::span<std::complex<double>> field, bool forward) {
    const std::size_t nx = m_isolated_workspace.nx;
    const std::size_t ny = m_isolated_workspace.ny;
    const std::size_t nz = m_isolated_workspace.nz;
    const std::size_t total = nx * ny * nz;
    std::vector<std::complex<double>> out(total, {0.0, 0.0});
    const double sign = forward ? -1.0 : 1.0;
    for (std::size_t kx = 0; kx < nx; ++kx) {
      for (std::size_t ky = 0; ky < ny; ++ky) {
        for (std::size_t kz = 0; kz < nz; ++kz) {
          std::complex<double> acc(0.0, 0.0);
          for (std::size_t x = 0; x < nx; ++x) {
            for (std::size_t y = 0; y < ny; ++y) {
              for (std::size_t z = 0; z < nz; ++z) {
                const double phase = sign * 2.0 * k_pi *
                    (static_cast<double>(kx * x) / static_cast<double>(nx) +
                     static_cast<double>(ky * y) / static_cast<double>(ny) +
                     static_cast<double>(kz * z) / static_cast<double>(nz));
                const std::complex<double> euler(std::cos(phase), std::sin(phase));
                acc += field[(x * ny + y) * nz + z] * euler;
              }
            }
          }
          out[(kx * ny + ky) * nz + kz] = acc;
        }
      }
    }
    if (!forward) {
      const double inv_total = 1.0 / static_cast<double>(total);
      for (auto& v : out) {
        v *= inv_total;
      }
    }
    std::copy(out.begin(), out.end(), field.begin());
  }

  void naiveForwardDft() {
    if (!activePlan().layout.ownsFullDomain()) {
      throw std::logic_error(
          "PM solver naiveForwardDft requires full-domain ownership; distributed PM is unavailable without FFTW/MPI");
    }
    auto& real = activePlan().real;
    auto& fourier = activePlan().fourier;
    const std::size_t nz_complex = m_shape.nz / 2U + 1U;
    for (std::size_t ix = 0; ix < m_shape.nx; ++ix) {
      for (std::size_t iy = 0; iy < m_shape.ny; ++iy) {
        for (std::size_t iz = 0; iz < nz_complex; ++iz) {
          std::complex<double> acc(0.0, 0.0);
          for (std::size_t x = 0; x < m_shape.nx; ++x) {
            for (std::size_t y = 0; y < m_shape.ny; ++y) {
              for (std::size_t z = 0; z < m_shape.nz; ++z) {
                const double phase = -2.0 * k_pi *
                    (static_cast<double>(ix * x) / static_cast<double>(m_shape.nx) +
                     static_cast<double>(iy * y) / static_cast<double>(m_shape.ny) +
                     static_cast<double>(iz * z) / static_cast<double>(m_shape.nz));
                const std::complex<double> euler(std::cos(phase), std::sin(phase));
                const std::size_t rindex = (x * m_shape.ny + y) * m_shape.nz + z;
                acc += real[rindex] * euler;
              }
            }
          }
          const std::size_t cindex = (ix * m_shape.ny + iy) * nz_complex + iz;
          fourier[cindex] = acc;
        }
      }
    }
  }

  void naiveInverseDft() {
    if (!activePlan().layout.ownsFullDomain()) {
      throw std::logic_error(
          "PM solver naiveInverseDft requires full-domain ownership; distributed PM is unavailable without FFTW/MPI");
    }
    auto& real = activePlan().real;
    auto& fourier = activePlan().fourier;
    const std::size_t total = m_shape.cellCount();
    std::fill(real.begin(), real.end(), 0.0);
    const std::size_t nz_complex = m_shape.nz / 2U + 1U;

    for (std::size_t x = 0; x < m_shape.nx; ++x) {
      for (std::size_t y = 0; y < m_shape.ny; ++y) {
        for (std::size_t z = 0; z < m_shape.nz; ++z) {
          std::complex<double> acc(0.0, 0.0);
          for (std::size_t ix = 0; ix < m_shape.nx; ++ix) {
            for (std::size_t iy = 0; iy < m_shape.ny; ++iy) {
              for (std::size_t iz = 0; iz < nz_complex; ++iz) {
                const double phase = 2.0 * k_pi *
                    (static_cast<double>(ix * x) / static_cast<double>(m_shape.nx) +
                     static_cast<double>(iy * y) / static_cast<double>(m_shape.ny) +
                     static_cast<double>(iz * z) / static_cast<double>(m_shape.nz));
                const std::complex<double> euler(std::cos(phase), std::sin(phase));
                const std::size_t cindex = (ix * m_shape.ny + iy) * nz_complex + iz;
                if (iz == 0 || (m_shape.nz % 2U == 0U && iz == m_shape.nz / 2U)) {
                  acc += fourier[cindex] * euler;
                } else {
                  acc += fourier[cindex] * euler + std::conj(fourier[cindex]) * std::conj(euler);
                }
              }
            }
          }
          real[(x * m_shape.ny + y) * m_shape.nz + z] = acc.real() / static_cast<double>(total);
        }
      }
    }
  }
#endif

  PmGridShape m_shape;
  std::unordered_map<PlanKey, PlanResources, PlanKeyHasher> m_plan_cache;
  std::optional<PlanKey> m_active_key;
  std::size_t m_plan_build_count = 0;
  struct {
    int world_size = 1;
    int world_rank = 0;
    DensityExchangeBuffers buffers;
  } m_density_exchange{};
  struct {
    int world_size = 1;
    int world_rank = 0;
    InterpolationExchangeBuffers<PmForceContributionRecord> buffers;
  } m_force_exchange{};
  struct {
    int world_size = 1;
    int world_rank = 0;
    InterpolationExchangeBuffers<PmPotentialContributionRecord> buffers;
  } m_potential_exchange{};
  IsolatedOpenWorkspace m_isolated_workspace{};
};

std::size_t PmGridShape::cellCount() const {
  return nx * ny * nz;
}

bool PmGridShape::isValid() const {
  return nx > 0 && ny > 0 && nz > 0;
}

void PmProfiler::reset() {
  m_totals = {};
}

void PmProfiler::append(const PmProfileEvent& event) {
  m_totals.bytes_moved += event.bytes_moved;
  m_totals.assign_ms += event.assign_ms;
  m_totals.fft_forward_ms += event.fft_forward_ms;
  m_totals.poisson_ms += event.poisson_ms;
  m_totals.gradient_ms += event.gradient_ms;
  m_totals.fft_inverse_ms += event.fft_inverse_ms;
  m_totals.fft_transpose_ms += event.fft_transpose_ms;
  m_totals.interpolate_ms += event.interpolate_ms;
  m_totals.transfer_h2d_ms += event.transfer_h2d_ms;
  m_totals.transfer_d2h_ms += event.transfer_d2h_ms;
  m_totals.device_kernel_ms += event.device_kernel_ms;
  m_totals.fft_transpose_bytes += event.fft_transpose_bytes;
}

const PmProfileEvent& PmProfiler::totals() const {
  return m_totals;
}

PmGridStorage::PmGridStorage(PmGridShape shape)
    : PmGridStorage(
          shape,
          parallel::makePmSlabLayout(shape.nx, shape.ny, shape.nz, /*world_size=*/1, /*world_rank=*/0)) {}

PmGridStorage::PmGridStorage(PmGridShape shape, parallel::PmSlabLayout layout)
    : m_shape(shape),
      m_layout(std::move(layout)),
      m_density(m_layout.localCellCount(), 0.0),
      m_potential(m_layout.localCellCount(), 0.0),
      m_force_x(m_layout.localCellCount(), 0.0),
      m_force_y(m_layout.localCellCount(), 0.0),
      m_force_z(m_layout.localCellCount(), 0.0) {
  if (!m_shape.isValid()) {
    throw std::invalid_argument("PM grid shape must be valid");
  }
  if (!m_layout.isValid()) {
    throw std::invalid_argument("PM slab layout must be valid");
  }
  if (m_layout.global_nx != m_shape.nx || m_layout.global_ny != m_shape.ny || m_layout.global_nz != m_shape.nz) {
    throw std::invalid_argument("PM slab layout global shape must match PM grid shape");
  }
}

const PmGridShape& PmGridStorage::shape() const {
  return m_shape;
}

const parallel::PmSlabLayout& PmGridStorage::slabLayout() const {
  return m_layout;
}

bool PmGridStorage::ownsFullDomain() const noexcept {
  return m_layout.ownsFullDomain();
}

std::size_t PmGridStorage::localCellCount() const noexcept {
  return m_layout.localCellCount();
}

std::span<double> PmGridStorage::density() {
  return m_density;
}

std::span<const double> PmGridStorage::density() const {
  return m_density;
}

std::span<double> PmGridStorage::potential() {
  return m_potential;
}

std::span<const double> PmGridStorage::potential() const {
  return m_potential;
}

std::span<double> PmGridStorage::force_x() {
  return m_force_x;
}

std::span<const double> PmGridStorage::force_x() const {
  return m_force_x;
}

std::span<double> PmGridStorage::force_y() {
  return m_force_y;
}

std::span<const double> PmGridStorage::force_y() const {
  return m_force_y;
}

std::span<double> PmGridStorage::force_z() {
  return m_force_z;
}

std::span<const double> PmGridStorage::force_z() const {
  return m_force_z;
}

std::size_t PmGridStorage::linearIndex(std::size_t ix, std::size_t iy, std::size_t iz) const {
  return m_layout.localLinearIndex(ix, iy, iz);
}

void PmGridStorage::clear() {
  std::fill(m_density.begin(), m_density.end(), 0.0);
  std::fill(m_potential.begin(), m_potential.end(), 0.0);
  std::fill(m_force_x.begin(), m_force_x.end(), 0.0);
  std::fill(m_force_y.begin(), m_force_y.end(), 0.0);
  std::fill(m_force_z.begin(), m_force_z.end(), 0.0);
}

PmSolver::PmSolver(PmGridShape shape) : m_shape(shape), m_impl(std::make_unique<Impl>(shape)) {
  if (!shape.isValid()) {
    throw std::invalid_argument("PM solver requires valid shape");
  }
}

PmSolver::~PmSolver() = default;
PmSolver::PmSolver(PmSolver&&) noexcept = default;
PmSolver& PmSolver::operator=(PmSolver&&) noexcept = default;

const PmGridShape& PmSolver::shape() const {
  return m_shape;
}

void PmSolver::assignDensity(
    PmGridStorage& grid,
    std::span<const double> pos_x,
    std::span<const double> pos_y,
    std::span<const double> pos_z,
    std::span<const double> mass,
    const PmSolveOptions& options,
    PmProfileEvent* profile) const {
  validateOptions(m_shape, options);
  if (grid.shape().cellCount() != m_shape.cellCount()) {
    throw std::invalid_argument("PM solver/grid shape mismatch in assignDensity");
  }
  if (pos_x.size() != pos_y.size() || pos_x.size() != pos_z.size() || pos_x.size() != mass.size()) {
    throw std::invalid_argument("Particle coordinate/mass spans must match in assignDensity");
  }
  if (!grid.slabLayout().isValid()) {
    throw std::invalid_argument("PmSolver::assignDensity requires a valid PM slab layout");
  }
  const bool distributed_slabs = grid.slabLayout().world_size > 1;
#if !COSMOSIM_ENABLE_MPI
  if (distributed_slabs) {
    throw std::invalid_argument(
        "PmSolver::assignDensity distributed slabs require COSMOSIM_ENABLE_MPI=ON");
  }
#else
  if (distributed_slabs) {
    int mpi_world_size = 1;
    int mpi_world_rank = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_world_rank);
    if (mpi_world_size != grid.slabLayout().world_size || mpi_world_rank != grid.slabLayout().world_rank) {
      throw std::invalid_argument(
          "PmSolver::assignDensity slab layout world metadata must match MPI_COMM_WORLD");
    }
  } else {
    validateSingleRankFullDomainGridContract(grid, "PmSolver::assignDensity");
  }
#endif

  const auto start = std::chrono::steady_clock::now();
  std::fill(grid.density().begin(), grid.density().end(), 0.0);

  const BoxLengths lengths = effectiveBoxLengths(options);
  const double inv_dx = static_cast<double>(m_shape.nx) / lengths.lx;
  const double inv_dy = static_cast<double>(m_shape.ny) / lengths.ly;
  const double inv_dz = static_cast<double>(m_shape.nz) / lengths.lz;

  const auto accumulate_owned = [&](const PmDensityContributionRecord& record) {
    if (record.global_ix >= m_shape.nx || record.global_iy >= m_shape.ny || record.global_iz >= m_shape.nz) {
      throw std::invalid_argument("PmSolver::assignDensity received out-of-range PM contribution record");
    }
    if (!grid.slabLayout().ownsGlobalX(record.global_ix)) {
      throw std::invalid_argument("PmSolver::assignDensity received contribution for non-owned PM slab x-index");
    }
    grid.density()[grid.linearIndex(record.global_ix, record.global_iy, record.global_iz)] += record.mass_contribution;
  };

  if (!distributed_slabs) {
    for (std::size_t p = 0; p < mass.size(); ++p) {
      const double x = wrapPosition(pos_x[p], lengths.lx) * inv_dx;
      const double y = wrapPosition(pos_y[p], lengths.ly) * inv_dy;
      const double z = wrapPosition(pos_z[p], lengths.lz) * inv_dz;

      const PmAxisStencil1d sx = makeAxisStencil(x, options.assignment_scheme);
      const PmAxisStencil1d sy = makeAxisStencil(y, options.assignment_scheme);
      const PmAxisStencil1d sz = makeAxisStencil(z, options.assignment_scheme);

      for (std::size_t dx = 0; dx < sx.count; ++dx) {
        const std::size_t ix = wrapIndex(sx.offsets[dx], m_shape.nx);
        for (std::size_t dy = 0; dy < sy.count; ++dy) {
          const std::size_t iy = wrapIndex(sy.offsets[dy], m_shape.ny);
          for (std::size_t dz = 0; dz < sz.count; ++dz) {
            const std::size_t iz = wrapIndex(sz.offsets[dz], m_shape.nz);
            const double weight = sx.weights[dx] * sy.weights[dy] * sz.weights[dz];
            accumulate_owned(PmDensityContributionRecord{
                .global_ix = static_cast<std::uint32_t>(ix),
                .global_iy = static_cast<std::uint32_t>(iy),
                .global_iz = static_cast<std::uint32_t>(iz),
                .mass_contribution = mass[p] * weight,
            });
          }
        }
      }
    }
  } else {
#if COSMOSIM_ENABLE_MPI
    auto& exchange = m_impl->densityExchangeBuffersForLayout(grid.slabLayout());
    for (auto& per_rank : exchange.send_records_by_rank) {
      per_rank.clear();
    }
    exchange.send_flat_records.clear();
    exchange.recv_flat_records.clear();
    std::fill(exchange.send_counts.begin(), exchange.send_counts.end(), 0);
    std::fill(exchange.send_displs.begin(), exchange.send_displs.end(), 0);
    std::fill(exchange.recv_counts.begin(), exchange.recv_counts.end(), 0);
    std::fill(exchange.recv_displs.begin(), exchange.recv_displs.end(), 0);
    std::fill(exchange.send_counts_bytes.begin(), exchange.send_counts_bytes.end(), 0);
    std::fill(exchange.send_displs_bytes.begin(), exchange.send_displs_bytes.end(), 0);
    std::fill(exchange.recv_counts_bytes.begin(), exchange.recv_counts_bytes.end(), 0);
    std::fill(exchange.recv_displs_bytes.begin(), exchange.recv_displs_bytes.end(), 0);

    for (std::size_t p = 0; p < mass.size(); ++p) {
      const double x = wrapPosition(pos_x[p], lengths.lx) * inv_dx;
      const double y = wrapPosition(pos_y[p], lengths.ly) * inv_dy;
      const double z = wrapPosition(pos_z[p], lengths.lz) * inv_dz;
      const PmAxisStencil1d sx = makeAxisStencil(x, options.assignment_scheme);
      const PmAxisStencil1d sy = makeAxisStencil(y, options.assignment_scheme);
      const PmAxisStencil1d sz = makeAxisStencil(z, options.assignment_scheme);

      for (std::size_t dx = 0; dx < sx.count; ++dx) {
        const std::size_t ix = wrapIndex(sx.offsets[dx], m_shape.nx);
        const int destination_rank = parallel::pmOwnerRankForGlobalX(m_shape.nx, grid.slabLayout().world_size, ix);
        auto& batch = exchange.send_records_by_rank[static_cast<std::size_t>(destination_rank)];
        for (std::size_t dy = 0; dy < sy.count; ++dy) {
          const std::size_t iy = wrapIndex(sy.offsets[dy], m_shape.ny);
          for (std::size_t dz = 0; dz < sz.count; ++dz) {
            const std::size_t iz = wrapIndex(sz.offsets[dz], m_shape.nz);
            const double weight = sx.weights[dx] * sy.weights[dy] * sz.weights[dz];
            batch.push_back(PmDensityContributionRecord{
                .global_ix = static_cast<std::uint32_t>(ix),
                .global_iy = static_cast<std::uint32_t>(iy),
                .global_iz = static_cast<std::uint32_t>(iz),
                .mass_contribution = mass[p] * weight,
            });
          }
        }
      }
    }

    std::size_t total_send_records = 0;
    for (int rank = 0; rank < grid.slabLayout().world_size; ++rank) {
      const std::size_t count = exchange.send_records_by_rank[static_cast<std::size_t>(rank)].size();
      if (count > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw std::invalid_argument("PmSolver::assignDensity send contribution count exceeds MPI int limit");
      }
      exchange.send_counts[static_cast<std::size_t>(rank)] = static_cast<int>(count);
      exchange.send_displs[static_cast<std::size_t>(rank)] = static_cast<int>(total_send_records);
      total_send_records += count;
    }
    exchange.send_flat_records.reserve(total_send_records);
    for (int rank = 0; rank < grid.slabLayout().world_size; ++rank) {
      const auto& records = exchange.send_records_by_rank[static_cast<std::size_t>(rank)];
      exchange.send_flat_records.insert(exchange.send_flat_records.end(), records.begin(), records.end());
    }

    MPI_Alltoall(
        exchange.send_counts.data(),
        1,
        MPI_INT,
        exchange.recv_counts.data(),
        1,
        MPI_INT,
        MPI_COMM_WORLD);

    std::size_t total_recv_records = 0;
    for (int rank = 0; rank < grid.slabLayout().world_size; ++rank) {
      exchange.recv_displs[static_cast<std::size_t>(rank)] = static_cast<int>(total_recv_records);
      total_recv_records += static_cast<std::size_t>(exchange.recv_counts[static_cast<std::size_t>(rank)]);
    }
    exchange.recv_flat_records.resize(total_recv_records);

    const int record_bytes = static_cast<int>(sizeof(PmDensityContributionRecord));
    for (int rank = 0; rank < grid.slabLayout().world_size; ++rank) {
      const auto idx = static_cast<std::size_t>(rank);
      exchange.send_counts_bytes[idx] = exchange.send_counts[idx] * record_bytes;
      exchange.send_displs_bytes[idx] = exchange.send_displs[idx] * record_bytes;
      exchange.recv_counts_bytes[idx] = exchange.recv_counts[idx] * record_bytes;
      exchange.recv_displs_bytes[idx] = exchange.recv_displs[idx] * record_bytes;
    }

    MPI_Alltoallv(
        reinterpret_cast<const std::uint8_t*>(exchange.send_flat_records.data()),
        exchange.send_counts_bytes.data(),
        exchange.send_displs_bytes.data(),
        MPI_BYTE,
        reinterpret_cast<std::uint8_t*>(exchange.recv_flat_records.data()),
        exchange.recv_counts_bytes.data(),
        exchange.recv_displs_bytes.data(),
        MPI_BYTE,
        MPI_COMM_WORLD);

    for (const PmDensityContributionRecord& record : exchange.recv_flat_records) {
      accumulate_owned(record);
    }
#endif
  }

  const double cell_volume =
      (lengths.lx * lengths.ly * lengths.lz) / static_cast<double>(m_shape.cellCount());
  for (double& density_cell : grid.density()) {
    density_cell /= cell_volume;
  }

  const auto stop = std::chrono::steady_clock::now();
  if (profile != nullptr) {
    profile->assign_ms += std::chrono::duration<double, std::milli>(stop - start).count();
    profile->bytes_moved += bytesForGridSweep(m_shape.cellCount());
    profile->bytes_moved += bytesForParticles(pos_x.size());
  }
}

void PmSolver::solvePoissonPeriodic(PmGridStorage& grid, const PmSolveOptions& options, PmProfileEvent* profile) {
  validateOptions(m_shape, options);
  if (grid.shape().nx != m_shape.nx || grid.shape().ny != m_shape.ny || grid.shape().nz != m_shape.nz) {
    throw std::invalid_argument("PM solver/grid shape mismatch in solvePoissonPeriodic");
  }
  if (!grid.slabLayout().isValid()) {
    throw std::invalid_argument("PmSolver::solvePoissonPeriodic requires valid PM slab layout");
  }
  if (grid.slabLayout().world_size > 1) {
#if !(COSMOSIM_ENABLE_FFTW && COSMOSIM_ENABLE_MPI)
    throw std::invalid_argument(
        "PmSolver::solvePoissonPeriodic distributed slabs require COSMOSIM_ENABLE_FFTW=ON and COSMOSIM_ENABLE_MPI=ON");
#endif
  } else if (!grid.ownsFullDomain()) {
    throw std::invalid_argument("PmSolver::solvePoissonPeriodic single-rank path requires full-domain slab ownership");
  }

  auto& plan = m_impl->planForLayout(grid.slabLayout(), options.decomposition_mode);
  auto real = m_impl->realGrid();
  auto fourier = m_impl->fourierGrid();
  auto potential_k = m_impl->potentialScratch();
  auto working_k = m_impl->workingScratch();

  std::copy(grid.density().begin(), grid.density().end(), real.begin());

  double local_density_sum = std::accumulate(real.begin(), real.end(), 0.0);
  double global_density_sum = local_density_sum;
#if COSMOSIM_ENABLE_MPI
  if (grid.slabLayout().world_size > 1) {
    MPI_Allreduce(&local_density_sum, &global_density_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  }
#endif
  const double mean_density = global_density_sum / static_cast<double>(m_shape.cellCount());
  for (double& value : real) {
    value -= mean_density;
  }

  const double forward_fft_ms = m_impl->forwardFft();
  if (profile != nullptr) {
    profile->fft_forward_ms += forward_fft_ms;
    if (plan.spectral_transposed) {
      profile->fft_transpose_ms += forward_fft_ms;
    }
  }

  const auto poisson_start = std::chrono::steady_clock::now();
  const std::size_t nz_complex = m_shape.nz / 2U + 1U;
  const double prefactor = -4.0 * k_pi * options.gravitational_constant_code * options.scale_factor * options.scale_factor;
  const BoxLengths lengths = effectiveBoxLengths(options);
  const double dkx = 2.0 * k_pi / lengths.lx;
  const double dky = 2.0 * k_pi / lengths.ly;
  const double dkz = 2.0 * k_pi / lengths.lz;

  const std::size_t global_x_begin = plan.layout.owned_x.begin_x;
  auto apply_poisson_stencil = [&](std::size_t ix, std::size_t iy, std::size_t base_index) {
    const std::ptrdiff_t nx_mode = ix <= m_shape.nx / 2U ? static_cast<std::ptrdiff_t>(ix)
                                                          : static_cast<std::ptrdiff_t>(ix) - static_cast<std::ptrdiff_t>(m_shape.nx);
    const std::ptrdiff_t ny_mode = iy <= m_shape.ny / 2U ? static_cast<std::ptrdiff_t>(iy)
                                                          : static_cast<std::ptrdiff_t>(iy) - static_cast<std::ptrdiff_t>(m_shape.ny);
    const double kx = dkx * static_cast<double>(nx_mode);
    const double ky = dky * static_cast<double>(ny_mode);
    for (std::size_t iz = 0; iz < nz_complex; ++iz) {
        const double kz = dkz * static_cast<double>(iz);
        const std::size_t index = base_index + iz;
        const double k2 = kx * kx + ky * ky + kz * kz;

        // Discrete mode mapping for r2c layout:
        //   kx(ix) = 2π/L * (ix <= Nx/2 ? ix : ix - Nx)
        //   ky(iy) = 2π/L * (iy <= Ny/2 ? iy : iy - Ny)
        //   kz(iz) = 2π/L * iz, iz in [0, Nz/2]
        // The negative kz branch is represented by Hermitian conjugates and is
        // not explicitly stored in the reduced-complex array.
        if (k2 == 0.0) {
          // Zero-mode pinning policy in periodic boxes: potential mean is fixed
          // to zero, so the solver enforces phi_{k=0}=0 and a_{k=0}=0.
          fourier[index] = {0.0, 0.0};
          continue;
        }

        double window_correction = 1.0;
        if (options.enable_window_deconvolution) {
          const int window_exponent = assignmentWindowExponent(options.assignment_scheme);
          const double wx = std::pow(
              sinc(0.5 * kx * lengths.lx / static_cast<double>(m_shape.nx)),
              static_cast<double>(window_exponent));
          const double wy = std::pow(
              sinc(0.5 * ky * lengths.ly / static_cast<double>(m_shape.ny)),
              static_cast<double>(window_exponent));
          const double wz = std::pow(
              sinc(0.5 * kz * lengths.lz / static_cast<double>(m_shape.nz)),
              static_cast<double>(window_exponent));
          const double transfer_window = wx * wy * wz;
          window_correction = 1.0 / std::max(transfer_window * transfer_window, 1.0e-12);
        }

        double split_filter = 1.0;
        if (options.tree_pm_split_scale_comoving > 0.0) {
          const double wave_number_comoving = std::sqrt(k2);
          split_filter = treePmGaussianFourierLongRangeFilter(wave_number_comoving, options.tree_pm_split_scale_comoving);
        }

        fourier[index] *= prefactor * window_correction * split_filter / k2;
    }
  };
  if (plan.spectral_transposed) {
    for (std::size_t local_iy = 0; local_iy < plan.transposed_local_ny; ++local_iy) {
      const std::size_t iy = plan.transposed_begin_y + local_iy;
      for (std::size_t ix = 0; ix < m_shape.nx; ++ix) {
        const std::size_t base = (local_iy * m_shape.nx + ix) * nz_complex;
        apply_poisson_stencil(ix, iy, base);
      }
    }
  } else {
    for (std::size_t local_ix = 0; local_ix < plan.layout.local_nx(); ++local_ix) {
      const std::size_t ix = global_x_begin + local_ix;
      for (std::size_t iy = 0; iy < m_shape.ny; ++iy) {
        const std::size_t base = (local_ix * m_shape.ny + iy) * nz_complex;
        apply_poisson_stencil(ix, iy, base);
      }
    }
  }
  const auto poisson_stop = std::chrono::steady_clock::now();

  if (profile != nullptr) {
    profile->poisson_ms += std::chrono::duration<double, std::milli>(poisson_stop - poisson_start).count();
  }

  std::copy(fourier.begin(), fourier.end(), potential_k.begin());

  auto inverse_into = [this, profile, &plan](std::span<const std::complex<double>> src, std::span<double> dst) {
    auto fourier_dst = m_impl->fourierGrid();
    std::copy(src.begin(), src.end(), fourier_dst.begin());
    const double fft_time = m_impl->inverseFft();
    auto real_values = m_impl->realGrid();
#if COSMOSIM_ENABLE_FFTW
    const double normalization = 1.0 / static_cast<double>(m_shape.cellCount());
    for (std::size_t i = 0; i < dst.size(); ++i) {
      dst[i] = real_values[i] * normalization;
    }
#else
    std::copy(real_values.begin(), real_values.end(), dst.begin());
#endif
    if (profile != nullptr) {
      profile->fft_inverse_ms += fft_time;
      if (plan.spectral_transposed) {
        profile->fft_transpose_ms += fft_time;
      }
    }
  };

  const auto grad_start = std::chrono::steady_clock::now();
  inverse_into(potential_k, grid.potential());

  auto fill_gradient_k = [&](std::span<std::complex<double>> dst, int axis) {
    if (plan.spectral_transposed) {
      for (std::size_t local_iy = 0; local_iy < plan.transposed_local_ny; ++local_iy) {
        const std::size_t iy = plan.transposed_begin_y + local_iy;
        const std::ptrdiff_t ny_mode = iy <= m_shape.ny / 2U ? static_cast<std::ptrdiff_t>(iy)
                                                              : static_cast<std::ptrdiff_t>(iy) - static_cast<std::ptrdiff_t>(m_shape.ny);
        const double ky = dky * static_cast<double>(ny_mode);
        for (std::size_t ix = 0; ix < m_shape.nx; ++ix) {
          const std::ptrdiff_t nx_mode = ix <= m_shape.nx / 2U ? static_cast<std::ptrdiff_t>(ix)
                                                                : static_cast<std::ptrdiff_t>(ix) - static_cast<std::ptrdiff_t>(m_shape.nx);
          const double kx = dkx * static_cast<double>(nx_mode);
          for (std::size_t iz = 0; iz < nz_complex; ++iz) {
            const double kz = dkz * static_cast<double>(iz);
            const std::size_t index = (local_iy * m_shape.nx + ix) * nz_complex + iz;
            const double component_k = axis == 0 ? kx : (axis == 1 ? ky : kz);
            dst[index] = std::complex<double>(0.0, -component_k) * potential_k[index];
          }
        }
      }
    } else {
      for (std::size_t local_ix = 0; local_ix < plan.layout.local_nx(); ++local_ix) {
        const std::size_t ix = global_x_begin + local_ix;
        const std::ptrdiff_t nx_mode = ix <= m_shape.nx / 2U ? static_cast<std::ptrdiff_t>(ix)
                                                              : static_cast<std::ptrdiff_t>(ix) - static_cast<std::ptrdiff_t>(m_shape.nx);
        const double kx = dkx * static_cast<double>(nx_mode);
        for (std::size_t iy = 0; iy < m_shape.ny; ++iy) {
          const std::ptrdiff_t ny_mode = iy <= m_shape.ny / 2U ? static_cast<std::ptrdiff_t>(iy)
                                                                : static_cast<std::ptrdiff_t>(iy) - static_cast<std::ptrdiff_t>(m_shape.ny);
          const double ky = dky * static_cast<double>(ny_mode);
          for (std::size_t iz = 0; iz < nz_complex; ++iz) {
            const double kz = dkz * static_cast<double>(iz);
            const std::size_t index = (local_ix * m_shape.ny + iy) * nz_complex + iz;
            const double component_k = axis == 0 ? kx : (axis == 1 ? ky : kz);
            dst[index] = std::complex<double>(0.0, -component_k) * potential_k[index];
          }
        }
      }
    }
  };

  fill_gradient_k(working_k, 0);
  inverse_into(working_k, grid.force_x());
  fill_gradient_k(working_k, 1);
  inverse_into(working_k, grid.force_y());
  fill_gradient_k(working_k, 2);
  inverse_into(working_k, grid.force_z());
  const auto grad_stop = std::chrono::steady_clock::now();
  if (profile != nullptr) {
    profile->gradient_ms += std::chrono::duration<double, std::milli>(grad_stop - grad_start).count();
  }

  if (profile != nullptr) {
    profile->bytes_moved += bytesForGridSweep(m_shape.cellCount()) * 6U;
    if (plan.spectral_transposed) {
      profile->fft_transpose_bytes +=
          static_cast<std::uint64_t>(sizeof(std::complex<double>)) * static_cast<std::uint64_t>(fourier.size()) * 8ULL;
    }
  }
}

void PmSolver::solvePoissonIsolatedOpen(PmGridStorage& grid, const PmSolveOptions& options, PmProfileEvent* profile) {
  validateOptions(m_shape, options);
  if (grid.shape().nx != m_shape.nx || grid.shape().ny != m_shape.ny || grid.shape().nz != m_shape.nz) {
    throw std::invalid_argument("PM solver/grid shape mismatch in solvePoissonIsolatedOpen");
  }
  if (!grid.slabLayout().isValid() || !grid.ownsFullDomain() || grid.slabLayout().world_size != 1) {
    throw std::invalid_argument("PmSolver::solvePoissonIsolatedOpen currently requires single-rank full-domain ownership");
  }

  const BoxLengths lengths = effectiveBoxLengths(options);
  const double dx = lengths.lx / static_cast<double>(m_shape.nx);
  const double dy = lengths.ly / static_cast<double>(m_shape.ny);
  const double dz = lengths.lz / static_cast<double>(m_shape.nz);
  if (m_shape.nx < 3 || m_shape.ny < 3 || m_shape.nz < 3) {
    throw std::invalid_argument("PmSolver::solvePoissonIsolatedOpen requires nx,ny,nz >= 3 for one-sided boundary gradients");
  }

  const auto poisson_start = std::chrono::steady_clock::now();
  const std::size_t pad_nx = 2U * m_shape.nx;
  const std::size_t pad_ny = 2U * m_shape.ny;
  const std::size_t pad_nz = 2U * m_shape.nz;
  const std::size_t pad_total = pad_nx * pad_ny * pad_nz;
  m_impl->ensureIsolatedWorkspace(pad_nx, pad_ny, pad_nz);
  auto& ws = m_impl->isolatedWorkspace();

  std::fill(ws.rho_k.begin(), ws.rho_k.end(), std::complex<double>(0.0, 0.0));
  for (std::size_t ix = 0; ix < m_shape.nx; ++ix) {
    for (std::size_t iy = 0; iy < m_shape.ny; ++iy) {
      for (std::size_t iz = 0; iz < m_shape.nz; ++iz) {
        const std::size_t src = grid.linearIndex(ix, iy, iz);
        const std::size_t dst = (ix * pad_ny + iy) * pad_nz + iz;
        ws.rho_k[dst] = {grid.density()[src], 0.0};
      }
    }
  }
  m_impl->isolatedForward(ws.rho_k);

  if (!ws.kernel_ready || ws.dx != dx || ws.dy != dy || ws.dz != dz ||
      ws.split_scale != options.tree_pm_split_scale_comoving) {
    ws.dx = dx;
    ws.dy = dy;
    ws.dz = dz;
    ws.split_scale = options.tree_pm_split_scale_comoving;
    std::fill(ws.kernel_k.begin(), ws.kernel_k.end(), std::complex<double>(0.0, 0.0));
    for (std::size_t ix = 0; ix < pad_nx; ++ix) {
      const std::ptrdiff_t sx = ix <= m_shape.nx
          ? static_cast<std::ptrdiff_t>(ix)
          : static_cast<std::ptrdiff_t>(ix) - static_cast<std::ptrdiff_t>(pad_nx);
      const double rx = static_cast<double>(sx) * dx;
      for (std::size_t iy = 0; iy < pad_ny; ++iy) {
        const std::ptrdiff_t sy = iy <= m_shape.ny
            ? static_cast<std::ptrdiff_t>(iy)
            : static_cast<std::ptrdiff_t>(iy) - static_cast<std::ptrdiff_t>(pad_ny);
        const double ry = static_cast<double>(sy) * dy;
        for (std::size_t iz = 0; iz < pad_nz; ++iz) {
          const std::ptrdiff_t sz = iz <= m_shape.nz
              ? static_cast<std::ptrdiff_t>(iz)
              : static_cast<std::ptrdiff_t>(iz) - static_cast<std::ptrdiff_t>(pad_nz);
          const double rz = static_cast<double>(sz) * dz;
          const double r = std::sqrt(rx * rx + ry * ry + rz * rz);
          double kernel = 0.0;  // self term policy: zero self-contribution at r=0.
          if (r > 0.0) {
            kernel = -1.0 / r;
            if (options.tree_pm_split_scale_comoving > 0.0) {
              kernel = -std::erf(0.5 * r / options.tree_pm_split_scale_comoving) / r;
            }
          }
          ws.kernel_k[(ix * pad_ny + iy) * pad_nz + iz] = {kernel, 0.0};
        }
      }
    }
    m_impl->isolatedForward(ws.kernel_k);
    ws.kernel_ready = true;
  }

  for (std::size_t i = 0; i < pad_total; ++i) {
    ws.rho_k[i] *= ws.kernel_k[i];
  }
  m_impl->isolatedInverse(ws.rho_k);

  const double cell_volume = dx * dy * dz;
  const double prefactor = options.gravitational_constant_code * options.scale_factor * options.scale_factor * cell_volume;
  for (std::size_t ix = 0; ix < m_shape.nx; ++ix) {
    for (std::size_t iy = 0; iy < m_shape.ny; ++iy) {
      for (std::size_t iz = 0; iz < m_shape.nz; ++iz) {
        const std::size_t physical_index = grid.linearIndex(ix, iy, iz);
        const std::size_t padded_index = (ix * pad_ny + iy) * pad_nz + iz;
        double value = ws.rho_k[padded_index].real();
#if COSMOSIM_ENABLE_FFTW
        value /= static_cast<double>(pad_total);
#endif
        grid.potential()[physical_index] = prefactor * value;
      }
    }
  }

  for (std::size_t ix = 0; ix < m_shape.nx; ++ix) {
    for (std::size_t iy = 0; iy < m_shape.ny; ++iy) {
      for (std::size_t iz = 0; iz < m_shape.nz; ++iz) {
        const std::size_t center = grid.linearIndex(ix, iy, iz);
        const auto grad_x = [&]() {
          if (ix == 0) {
            return (-3.0 * grid.potential()[grid.linearIndex(0, iy, iz)] +
                    4.0 * grid.potential()[grid.linearIndex(1, iy, iz)] -
                    grid.potential()[grid.linearIndex(2, iy, iz)]) / (2.0 * dx);
          }
          if (ix + 1U == m_shape.nx) {
            return (3.0 * grid.potential()[grid.linearIndex(m_shape.nx - 1U, iy, iz)] -
                    4.0 * grid.potential()[grid.linearIndex(m_shape.nx - 2U, iy, iz)] +
                    grid.potential()[grid.linearIndex(m_shape.nx - 3U, iy, iz)]) / (2.0 * dx);
          }
          return (grid.potential()[grid.linearIndex(ix + 1U, iy, iz)] -
                  grid.potential()[grid.linearIndex(ix - 1U, iy, iz)]) / (2.0 * dx);
        }();
        const auto grad_y = [&]() {
          if (iy == 0) {
            return (-3.0 * grid.potential()[grid.linearIndex(ix, 0, iz)] +
                    4.0 * grid.potential()[grid.linearIndex(ix, 1, iz)] -
                    grid.potential()[grid.linearIndex(ix, 2, iz)]) / (2.0 * dy);
          }
          if (iy + 1U == m_shape.ny) {
            return (3.0 * grid.potential()[grid.linearIndex(ix, m_shape.ny - 1U, iz)] -
                    4.0 * grid.potential()[grid.linearIndex(ix, m_shape.ny - 2U, iz)] +
                    grid.potential()[grid.linearIndex(ix, m_shape.ny - 3U, iz)]) / (2.0 * dy);
          }
          return (grid.potential()[grid.linearIndex(ix, iy + 1U, iz)] -
                  grid.potential()[grid.linearIndex(ix, iy - 1U, iz)]) / (2.0 * dy);
        }();
        const auto grad_z = [&]() {
          if (iz == 0) {
            return (-3.0 * grid.potential()[grid.linearIndex(ix, iy, 0)] +
                    4.0 * grid.potential()[grid.linearIndex(ix, iy, 1)] -
                    grid.potential()[grid.linearIndex(ix, iy, 2)]) / (2.0 * dz);
          }
          if (iz + 1U == m_shape.nz) {
            return (3.0 * grid.potential()[grid.linearIndex(ix, iy, m_shape.nz - 1U)] -
                    4.0 * grid.potential()[grid.linearIndex(ix, iy, m_shape.nz - 2U)] +
                    grid.potential()[grid.linearIndex(ix, iy, m_shape.nz - 3U)]) / (2.0 * dz);
          }
          return (grid.potential()[grid.linearIndex(ix, iy, iz + 1U)] -
                  grid.potential()[grid.linearIndex(ix, iy, iz - 1U)]) / (2.0 * dz);
        }();
        grid.force_x()[center] = -grad_x;
        grid.force_y()[center] = -grad_y;
        grid.force_z()[center] = -grad_z;
      }
    }
  }

  if (profile != nullptr) {
    const auto stop = std::chrono::steady_clock::now();
    profile->poisson_ms += std::chrono::duration<double, std::milli>(stop - poisson_start).count();
    profile->gradient_ms += 0.0;
  }
}

void PmSolver::interpolateForces(
    const PmGridStorage& grid,
    std::span<const double> pos_x,
    std::span<const double> pos_y,
    std::span<const double> pos_z,
    std::span<double> accel_x,
    std::span<double> accel_y,
    std::span<double> accel_z,
    const PmSolveOptions& options,
    PmProfileEvent* profile) const {
  validateOptions(m_shape, options);
  if (pos_x.size() != pos_y.size() || pos_x.size() != pos_z.size() || pos_x.size() != accel_x.size() ||
      pos_x.size() != accel_y.size() || pos_x.size() != accel_z.size()) {
    throw std::invalid_argument("Particle coordinate/acceleration spans must match in interpolateForces");
  }
  if (!grid.slabLayout().isValid()) {
    throw std::invalid_argument("PmSolver::interpolateForces requires a valid PM slab layout");
  }

  const bool distributed_slabs = grid.slabLayout().world_size > 1;
#if !COSMOSIM_ENABLE_MPI
  if (distributed_slabs) {
    throw std::invalid_argument("PmSolver::interpolateForces distributed slabs require COSMOSIM_ENABLE_MPI=ON");
  }
  validateSingleRankFullDomainGridContract(grid, "PmSolver::interpolateForces");
#else
  if (distributed_slabs) {
    int mpi_world_size = 1;
    int mpi_world_rank = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_world_rank);
    if (mpi_world_size != grid.slabLayout().world_size || mpi_world_rank != grid.slabLayout().world_rank) {
      throw std::invalid_argument("PmSolver::interpolateForces slab layout world metadata must match MPI_COMM_WORLD");
    }
  } else {
    validateSingleRankFullDomainGridContract(grid, "PmSolver::interpolateForces");
  }
#endif

  const auto start = std::chrono::steady_clock::now();

  const BoxLengths lengths = effectiveBoxLengths(options);
  const double inv_dx = static_cast<double>(m_shape.nx) / lengths.lx;
  const double inv_dy = static_cast<double>(m_shape.ny) / lengths.ly;
  const double inv_dz = static_cast<double>(m_shape.nz) / lengths.lz;

  if (!distributed_slabs) {
    for (std::size_t p = 0; p < pos_x.size(); ++p) {
      const double x = wrapPosition(pos_x[p], lengths.lx) * inv_dx;
      const double y = wrapPosition(pos_y[p], lengths.ly) * inv_dy;
      const double z = wrapPosition(pos_z[p], lengths.lz) * inv_dz;

      const PmAxisStencil1d sx = makeAxisStencil(x, options.assignment_scheme);
      const PmAxisStencil1d sy = makeAxisStencil(y, options.assignment_scheme);
      const PmAxisStencil1d sz = makeAxisStencil(z, options.assignment_scheme);

      double gx = 0.0;
      double gy = 0.0;
      double gz = 0.0;

      for (std::size_t dx = 0; dx < sx.count; ++dx) {
        const std::size_t ix = wrapIndex(sx.offsets[dx], m_shape.nx);
        for (std::size_t dy = 0; dy < sy.count; ++dy) {
          const std::size_t iy = wrapIndex(sy.offsets[dy], m_shape.ny);
          for (std::size_t dz = 0; dz < sz.count; ++dz) {
            const std::size_t iz = wrapIndex(sz.offsets[dz], m_shape.nz);
            const double weight = sx.weights[dx] * sy.weights[dy] * sz.weights[dz];
            const std::size_t index = grid.linearIndex(ix, iy, iz);
            gx += weight * grid.force_x()[index];
            gy += weight * grid.force_y()[index];
            gz += weight * grid.force_z()[index];
          }
        }
      }

      accel_x[p] = gx;
      accel_y[p] = gy;
      accel_z[p] = gz;
    }
  } else {
#if COSMOSIM_ENABLE_MPI
    const int world_size = grid.slabLayout().world_size;
    auto& exchange = m_impl->forceInterpolationExchangeBuffersForLayout(grid.slabLayout());
    for (auto& per_rank : exchange.send_requests_by_rank) {
      per_rank.clear();
    }
    for (auto& per_rank : exchange.send_contribs_by_rank) {
      per_rank.clear();
    }
    exchange.send_requests_flat.clear();
    exchange.recv_requests_flat.clear();
    exchange.send_contribs_flat.clear();
    exchange.recv_contribs_flat.clear();
    std::fill(exchange.send_counts.begin(), exchange.send_counts.end(), 0);
    std::fill(exchange.send_displs.begin(), exchange.send_displs.end(), 0);
    std::fill(exchange.recv_counts.begin(), exchange.recv_counts.end(), 0);
    std::fill(exchange.recv_displs.begin(), exchange.recv_displs.end(), 0);
    std::fill(exchange.send_counts_bytes.begin(), exchange.send_counts_bytes.end(), 0);
    std::fill(exchange.send_displs_bytes.begin(), exchange.send_displs_bytes.end(), 0);
    std::fill(exchange.recv_counts_bytes.begin(), exchange.recv_counts_bytes.end(), 0);
    std::fill(exchange.recv_displs_bytes.begin(), exchange.recv_displs_bytes.end(), 0);
    std::fill(exchange.send_contrib_counts.begin(), exchange.send_contrib_counts.end(), 0);
    std::fill(exchange.send_contrib_displs.begin(), exchange.send_contrib_displs.end(), 0);
    std::fill(exchange.recv_contrib_counts.begin(), exchange.recv_contrib_counts.end(), 0);
    std::fill(exchange.recv_contrib_displs.begin(), exchange.recv_contrib_displs.end(), 0);
    std::fill(exchange.send_contrib_counts_bytes.begin(), exchange.send_contrib_counts_bytes.end(), 0);
    std::fill(exchange.send_contrib_displs_bytes.begin(), exchange.send_contrib_displs_bytes.end(), 0);
    std::fill(exchange.recv_contrib_counts_bytes.begin(), exchange.recv_contrib_counts_bytes.end(), 0);
    std::fill(exchange.recv_contrib_displs_bytes.begin(), exchange.recv_contrib_displs_bytes.end(), 0);

    for (std::size_t p = 0; p < pos_x.size(); ++p) {
      const double x = wrapPosition(pos_x[p], lengths.lx) * inv_dx;
      const double y = wrapPosition(pos_y[p], lengths.ly) * inv_dy;
      const double z = wrapPosition(pos_z[p], lengths.lz) * inv_dz;
      const PmAxisStencil1d sx = makeAxisStencil(x, options.assignment_scheme);
      const PmAxisStencil1d sy = makeAxisStencil(y, options.assignment_scheme);
      const PmAxisStencil1d sz = makeAxisStencil(z, options.assignment_scheme);
      for (std::size_t dx = 0; dx < sx.count; ++dx) {
        const std::size_t ix = wrapIndex(sx.offsets[dx], m_shape.nx);
        const int destination_rank = parallel::pmOwnerRankForGlobalX(m_shape.nx, world_size, ix);
        auto& batch = exchange.send_requests_by_rank[static_cast<std::size_t>(destination_rank)];
        for (std::size_t dy = 0; dy < sy.count; ++dy) {
          const std::size_t iy = wrapIndex(sy.offsets[dy], m_shape.ny);
          for (std::size_t dz = 0; dz < sz.count; ++dz) {
            const std::size_t iz = wrapIndex(sz.offsets[dz], m_shape.nz);
            const double weight = sx.weights[dx] * sy.weights[dy] * sz.weights[dz];
            batch.push_back(PmInterpolationRequestRecord{
                .particle_index = static_cast<std::uint32_t>(p),
                .global_ix = static_cast<std::uint32_t>(ix),
                .global_iy = static_cast<std::uint32_t>(iy),
                .global_iz = static_cast<std::uint32_t>(iz),
                .weight = weight,
            });
          }
        }
      }
    }

    std::size_t total_send = 0;
    for (int rank = 0; rank < world_size; ++rank) {
      const std::size_t count = exchange.send_requests_by_rank[static_cast<std::size_t>(rank)].size();
      if (count > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw std::invalid_argument("PmSolver::interpolateForces request count exceeds MPI int limit");
      }
      exchange.send_counts[static_cast<std::size_t>(rank)] = static_cast<int>(count);
      exchange.send_displs[static_cast<std::size_t>(rank)] = static_cast<int>(total_send);
      total_send += count;
    }

    exchange.send_requests_flat.reserve(total_send);
    for (int rank = 0; rank < world_size; ++rank) {
      const auto& source = exchange.send_requests_by_rank[static_cast<std::size_t>(rank)];
      exchange.send_requests_flat.insert(exchange.send_requests_flat.end(), source.begin(), source.end());
    }

    MPI_Alltoall(
        exchange.send_counts.data(), 1, MPI_INT, exchange.recv_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);

    std::size_t total_recv = 0;
    for (int rank = 0; rank < world_size; ++rank) {
      exchange.recv_displs[static_cast<std::size_t>(rank)] = static_cast<int>(total_recv);
      total_recv += static_cast<std::size_t>(exchange.recv_counts[static_cast<std::size_t>(rank)]);
    }
    exchange.recv_requests_flat.resize(total_recv);

    const int request_bytes = static_cast<int>(sizeof(PmInterpolationRequestRecord));
    for (int rank = 0; rank < world_size; ++rank) {
      const auto r = static_cast<std::size_t>(rank);
      exchange.send_counts_bytes[r] = exchange.send_counts[r] * request_bytes;
      exchange.send_displs_bytes[r] = exchange.send_displs[r] * request_bytes;
      exchange.recv_counts_bytes[r] = exchange.recv_counts[r] * request_bytes;
      exchange.recv_displs_bytes[r] = exchange.recv_displs[r] * request_bytes;
    }
    MPI_Alltoallv(
        reinterpret_cast<const std::uint8_t*>(exchange.send_requests_flat.data()),
        exchange.send_counts_bytes.data(),
        exchange.send_displs_bytes.data(),
        MPI_BYTE,
        reinterpret_cast<std::uint8_t*>(exchange.recv_requests_flat.data()),
        exchange.recv_counts_bytes.data(),
        exchange.recv_displs_bytes.data(),
        MPI_BYTE,
        MPI_COMM_WORLD);

    for (int source_rank = 0; source_rank < world_size; ++source_rank) {
      auto& batch = exchange.send_contribs_by_rank[static_cast<std::size_t>(source_rank)];
      const int begin = exchange.recv_displs[static_cast<std::size_t>(source_rank)];
      const int count = exchange.recv_counts[static_cast<std::size_t>(source_rank)];
      for (int i = 0; i < count; ++i) {
        const auto& request = exchange.recv_requests_flat[static_cast<std::size_t>(begin + i)];
        if (request.particle_index >= pos_x.size()) {
          throw std::invalid_argument("PmSolver::interpolateForces request particle index out of range");
        }
        if (request.global_ix >= m_shape.nx || request.global_iy >= m_shape.ny || request.global_iz >= m_shape.nz) {
          throw std::invalid_argument("PmSolver::interpolateForces request PM index out of range");
        }
        if (!grid.slabLayout().ownsGlobalX(request.global_ix)) {
          throw std::invalid_argument("PmSolver::interpolateForces request x-index is not owned by receiving slab");
        }
        const std::size_t index = grid.linearIndex(request.global_ix, request.global_iy, request.global_iz);
        batch.push_back(PmForceContributionRecord{
            .particle_index = request.particle_index,
            .accel_x = request.weight * grid.force_x()[index],
            .accel_y = request.weight * grid.force_y()[index],
            .accel_z = request.weight * grid.force_z()[index],
        });
      }
    }

    std::size_t total_send_contribs = 0;
    for (int rank = 0; rank < world_size; ++rank) {
      const std::size_t count = exchange.send_contribs_by_rank[static_cast<std::size_t>(rank)].size();
      if (count > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw std::invalid_argument("PmSolver::interpolateForces contribution count exceeds MPI int limit");
      }
      exchange.send_contrib_counts[static_cast<std::size_t>(rank)] = static_cast<int>(count);
      exchange.send_contrib_displs[static_cast<std::size_t>(rank)] = static_cast<int>(total_send_contribs);
      total_send_contribs += count;
    }

    exchange.send_contribs_flat.reserve(total_send_contribs);
    for (int rank = 0; rank < world_size; ++rank) {
      const auto& source = exchange.send_contribs_by_rank[static_cast<std::size_t>(rank)];
      exchange.send_contribs_flat.insert(exchange.send_contribs_flat.end(), source.begin(), source.end());
    }

    MPI_Alltoall(
        exchange.send_contrib_counts.data(),
        1,
        MPI_INT,
        exchange.recv_contrib_counts.data(),
        1,
        MPI_INT,
        MPI_COMM_WORLD);

    std::size_t total_recv_contribs = 0;
    for (int rank = 0; rank < world_size; ++rank) {
      exchange.recv_contrib_displs[static_cast<std::size_t>(rank)] = static_cast<int>(total_recv_contribs);
      total_recv_contribs += static_cast<std::size_t>(exchange.recv_contrib_counts[static_cast<std::size_t>(rank)]);
    }
    exchange.recv_contribs_flat.resize(total_recv_contribs);

    const int contrib_bytes = static_cast<int>(sizeof(PmForceContributionRecord));
    for (int rank = 0; rank < world_size; ++rank) {
      const auto r = static_cast<std::size_t>(rank);
      exchange.send_contrib_counts_bytes[r] = exchange.send_contrib_counts[r] * contrib_bytes;
      exchange.send_contrib_displs_bytes[r] = exchange.send_contrib_displs[r] * contrib_bytes;
      exchange.recv_contrib_counts_bytes[r] = exchange.recv_contrib_counts[r] * contrib_bytes;
      exchange.recv_contrib_displs_bytes[r] = exchange.recv_contrib_displs[r] * contrib_bytes;
    }
    MPI_Alltoallv(
        reinterpret_cast<const std::uint8_t*>(exchange.send_contribs_flat.data()),
        exchange.send_contrib_counts_bytes.data(),
        exchange.send_contrib_displs_bytes.data(),
        MPI_BYTE,
        reinterpret_cast<std::uint8_t*>(exchange.recv_contribs_flat.data()),
        exchange.recv_contrib_counts_bytes.data(),
        exchange.recv_contrib_displs_bytes.data(),
        MPI_BYTE,
        MPI_COMM_WORLD);

    std::fill(accel_x.begin(), accel_x.end(), 0.0);
    std::fill(accel_y.begin(), accel_y.end(), 0.0);
    std::fill(accel_z.begin(), accel_z.end(), 0.0);
    for (const auto& contribution : exchange.recv_contribs_flat) {
      if (contribution.particle_index >= pos_x.size()) {
        throw std::invalid_argument("PmSolver::interpolateForces response particle index out of range");
      }
      const std::size_t p = static_cast<std::size_t>(contribution.particle_index);
      accel_x[p] += contribution.accel_x;
      accel_y[p] += contribution.accel_y;
      accel_z[p] += contribution.accel_z;
    }
#endif
  }

  const auto stop = std::chrono::steady_clock::now();
  if (profile != nullptr) {
    profile->interpolate_ms += std::chrono::duration<double, std::milli>(stop - start).count();
    profile->bytes_moved += bytesForParticles(pos_x.size());
  }
}

void PmSolver::interpolatePotential(
    const PmGridStorage& grid,
    std::span<const double> pos_x,
    std::span<const double> pos_y,
    std::span<const double> pos_z,
    std::span<double> potential,
    const PmSolveOptions& options,
    PmProfileEvent* profile) const {
  validateOptions(m_shape, options);
  if (pos_x.size() != pos_y.size() || pos_x.size() != pos_z.size() || pos_x.size() != potential.size()) {
    throw std::invalid_argument("Particle coordinate/potential spans must match in interpolatePotential");
  }
  if (!grid.slabLayout().isValid()) {
    throw std::invalid_argument("PmSolver::interpolatePotential requires a valid PM slab layout");
  }

  const bool distributed_slabs = grid.slabLayout().world_size > 1;
#if !COSMOSIM_ENABLE_MPI
  if (distributed_slabs) {
    throw std::invalid_argument("PmSolver::interpolatePotential distributed slabs require COSMOSIM_ENABLE_MPI=ON");
  }
  validateSingleRankFullDomainGridContract(grid, "PmSolver::interpolatePotential");
#else
  if (distributed_slabs) {
    int mpi_world_size = 1;
    int mpi_world_rank = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_world_rank);
    if (mpi_world_size != grid.slabLayout().world_size || mpi_world_rank != grid.slabLayout().world_rank) {
      throw std::invalid_argument(
          "PmSolver::interpolatePotential slab layout world metadata must match MPI_COMM_WORLD");
    }
  } else {
    validateSingleRankFullDomainGridContract(grid, "PmSolver::interpolatePotential");
  }
#endif

  const auto start = std::chrono::steady_clock::now();
  const BoxLengths lengths = effectiveBoxLengths(options);
  const double inv_dx = static_cast<double>(m_shape.nx) / lengths.lx;
  const double inv_dy = static_cast<double>(m_shape.ny) / lengths.ly;
  const double inv_dz = static_cast<double>(m_shape.nz) / lengths.lz;

  if (!distributed_slabs) {
    for (std::size_t p = 0; p < pos_x.size(); ++p) {
      const double x = wrapPosition(pos_x[p], lengths.lx) * inv_dx;
      const double y = wrapPosition(pos_y[p], lengths.ly) * inv_dy;
      const double z = wrapPosition(pos_z[p], lengths.lz) * inv_dz;

      const PmAxisStencil1d sx = makeAxisStencil(x, options.assignment_scheme);
      const PmAxisStencil1d sy = makeAxisStencil(y, options.assignment_scheme);
      const PmAxisStencil1d sz = makeAxisStencil(z, options.assignment_scheme);

      double phi = 0.0;
      for (std::size_t dx = 0; dx < sx.count; ++dx) {
        const std::size_t ix = wrapIndex(sx.offsets[dx], m_shape.nx);
        for (std::size_t dy = 0; dy < sy.count; ++dy) {
          const std::size_t iy = wrapIndex(sy.offsets[dy], m_shape.ny);
          for (std::size_t dz = 0; dz < sz.count; ++dz) {
            const std::size_t iz = wrapIndex(sz.offsets[dz], m_shape.nz);
            const double weight = sx.weights[dx] * sy.weights[dy] * sz.weights[dz];
            phi += weight * grid.potential()[grid.linearIndex(ix, iy, iz)];
          }
        }
      }
      potential[p] = phi;
    }
  } else {
#if COSMOSIM_ENABLE_MPI
    const int world_size = grid.slabLayout().world_size;
    auto& exchange = m_impl->potentialInterpolationExchangeBuffersForLayout(grid.slabLayout());
    for (auto& per_rank : exchange.send_requests_by_rank) {
      per_rank.clear();
    }
    for (auto& per_rank : exchange.send_contribs_by_rank) {
      per_rank.clear();
    }
    exchange.send_requests_flat.clear();
    exchange.recv_requests_flat.clear();
    exchange.send_contribs_flat.clear();
    exchange.recv_contribs_flat.clear();
    std::fill(exchange.send_counts.begin(), exchange.send_counts.end(), 0);
    std::fill(exchange.send_displs.begin(), exchange.send_displs.end(), 0);
    std::fill(exchange.recv_counts.begin(), exchange.recv_counts.end(), 0);
    std::fill(exchange.recv_displs.begin(), exchange.recv_displs.end(), 0);
    std::fill(exchange.send_counts_bytes.begin(), exchange.send_counts_bytes.end(), 0);
    std::fill(exchange.send_displs_bytes.begin(), exchange.send_displs_bytes.end(), 0);
    std::fill(exchange.recv_counts_bytes.begin(), exchange.recv_counts_bytes.end(), 0);
    std::fill(exchange.recv_displs_bytes.begin(), exchange.recv_displs_bytes.end(), 0);
    std::fill(exchange.send_contrib_counts.begin(), exchange.send_contrib_counts.end(), 0);
    std::fill(exchange.send_contrib_displs.begin(), exchange.send_contrib_displs.end(), 0);
    std::fill(exchange.recv_contrib_counts.begin(), exchange.recv_contrib_counts.end(), 0);
    std::fill(exchange.recv_contrib_displs.begin(), exchange.recv_contrib_displs.end(), 0);
    std::fill(exchange.send_contrib_counts_bytes.begin(), exchange.send_contrib_counts_bytes.end(), 0);
    std::fill(exchange.send_contrib_displs_bytes.begin(), exchange.send_contrib_displs_bytes.end(), 0);
    std::fill(exchange.recv_contrib_counts_bytes.begin(), exchange.recv_contrib_counts_bytes.end(), 0);
    std::fill(exchange.recv_contrib_displs_bytes.begin(), exchange.recv_contrib_displs_bytes.end(), 0);

    for (std::size_t p = 0; p < pos_x.size(); ++p) {
      const double x = wrapPosition(pos_x[p], lengths.lx) * inv_dx;
      const double y = wrapPosition(pos_y[p], lengths.ly) * inv_dy;
      const double z = wrapPosition(pos_z[p], lengths.lz) * inv_dz;
      const PmAxisStencil1d sx = makeAxisStencil(x, options.assignment_scheme);
      const PmAxisStencil1d sy = makeAxisStencil(y, options.assignment_scheme);
      const PmAxisStencil1d sz = makeAxisStencil(z, options.assignment_scheme);
      for (std::size_t dx = 0; dx < sx.count; ++dx) {
        const std::size_t ix = wrapIndex(sx.offsets[dx], m_shape.nx);
        const int destination_rank = parallel::pmOwnerRankForGlobalX(m_shape.nx, world_size, ix);
        auto& batch = exchange.send_requests_by_rank[static_cast<std::size_t>(destination_rank)];
        for (std::size_t dy = 0; dy < sy.count; ++dy) {
          const std::size_t iy = wrapIndex(sy.offsets[dy], m_shape.ny);
          for (std::size_t dz = 0; dz < sz.count; ++dz) {
            const std::size_t iz = wrapIndex(sz.offsets[dz], m_shape.nz);
            const double weight = sx.weights[dx] * sy.weights[dy] * sz.weights[dz];
            batch.push_back(PmInterpolationRequestRecord{
                .particle_index = static_cast<std::uint32_t>(p),
                .global_ix = static_cast<std::uint32_t>(ix),
                .global_iy = static_cast<std::uint32_t>(iy),
                .global_iz = static_cast<std::uint32_t>(iz),
                .weight = weight,
            });
          }
        }
      }
    }

    std::size_t total_send = 0;
    for (int rank = 0; rank < world_size; ++rank) {
      const std::size_t count = exchange.send_requests_by_rank[static_cast<std::size_t>(rank)].size();
      if (count > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw std::invalid_argument("PmSolver::interpolatePotential request count exceeds MPI int limit");
      }
      exchange.send_counts[static_cast<std::size_t>(rank)] = static_cast<int>(count);
      exchange.send_displs[static_cast<std::size_t>(rank)] = static_cast<int>(total_send);
      total_send += count;
    }

    exchange.send_requests_flat.reserve(total_send);
    for (int rank = 0; rank < world_size; ++rank) {
      const auto& source = exchange.send_requests_by_rank[static_cast<std::size_t>(rank)];
      exchange.send_requests_flat.insert(exchange.send_requests_flat.end(), source.begin(), source.end());
    }

    MPI_Alltoall(
        exchange.send_counts.data(), 1, MPI_INT, exchange.recv_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);

    std::size_t total_recv = 0;
    for (int rank = 0; rank < world_size; ++rank) {
      exchange.recv_displs[static_cast<std::size_t>(rank)] = static_cast<int>(total_recv);
      total_recv += static_cast<std::size_t>(exchange.recv_counts[static_cast<std::size_t>(rank)]);
    }
    exchange.recv_requests_flat.resize(total_recv);

    const int request_bytes = static_cast<int>(sizeof(PmInterpolationRequestRecord));
    for (int rank = 0; rank < world_size; ++rank) {
      const auto r = static_cast<std::size_t>(rank);
      exchange.send_counts_bytes[r] = exchange.send_counts[r] * request_bytes;
      exchange.send_displs_bytes[r] = exchange.send_displs[r] * request_bytes;
      exchange.recv_counts_bytes[r] = exchange.recv_counts[r] * request_bytes;
      exchange.recv_displs_bytes[r] = exchange.recv_displs[r] * request_bytes;
    }
    MPI_Alltoallv(
        reinterpret_cast<const std::uint8_t*>(exchange.send_requests_flat.data()),
        exchange.send_counts_bytes.data(),
        exchange.send_displs_bytes.data(),
        MPI_BYTE,
        reinterpret_cast<std::uint8_t*>(exchange.recv_requests_flat.data()),
        exchange.recv_counts_bytes.data(),
        exchange.recv_displs_bytes.data(),
        MPI_BYTE,
        MPI_COMM_WORLD);

    for (int source_rank = 0; source_rank < world_size; ++source_rank) {
      auto& batch = exchange.send_contribs_by_rank[static_cast<std::size_t>(source_rank)];
      const int begin = exchange.recv_displs[static_cast<std::size_t>(source_rank)];
      const int count = exchange.recv_counts[static_cast<std::size_t>(source_rank)];
      for (int i = 0; i < count; ++i) {
        const auto& request = exchange.recv_requests_flat[static_cast<std::size_t>(begin + i)];
        if (request.particle_index >= pos_x.size()) {
          throw std::invalid_argument("PmSolver::interpolatePotential request particle index out of range");
        }
        if (request.global_ix >= m_shape.nx || request.global_iy >= m_shape.ny || request.global_iz >= m_shape.nz) {
          throw std::invalid_argument("PmSolver::interpolatePotential request PM index out of range");
        }
        if (!grid.slabLayout().ownsGlobalX(request.global_ix)) {
          throw std::invalid_argument("PmSolver::interpolatePotential request x-index is not owned by receiving slab");
        }
        const std::size_t index = grid.linearIndex(request.global_ix, request.global_iy, request.global_iz);
        batch.push_back(PmPotentialContributionRecord{
            .particle_index = request.particle_index,
            .potential = request.weight * grid.potential()[index],
        });
      }
    }

    std::size_t total_send_contribs = 0;
    for (int rank = 0; rank < world_size; ++rank) {
      const std::size_t count = exchange.send_contribs_by_rank[static_cast<std::size_t>(rank)].size();
      if (count > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw std::invalid_argument("PmSolver::interpolatePotential contribution count exceeds MPI int limit");
      }
      exchange.send_contrib_counts[static_cast<std::size_t>(rank)] = static_cast<int>(count);
      exchange.send_contrib_displs[static_cast<std::size_t>(rank)] = static_cast<int>(total_send_contribs);
      total_send_contribs += count;
    }

    exchange.send_contribs_flat.reserve(total_send_contribs);
    for (int rank = 0; rank < world_size; ++rank) {
      const auto& source = exchange.send_contribs_by_rank[static_cast<std::size_t>(rank)];
      exchange.send_contribs_flat.insert(exchange.send_contribs_flat.end(), source.begin(), source.end());
    }

    MPI_Alltoall(
        exchange.send_contrib_counts.data(),
        1,
        MPI_INT,
        exchange.recv_contrib_counts.data(),
        1,
        MPI_INT,
        MPI_COMM_WORLD);

    std::size_t total_recv_contribs = 0;
    for (int rank = 0; rank < world_size; ++rank) {
      exchange.recv_contrib_displs[static_cast<std::size_t>(rank)] = static_cast<int>(total_recv_contribs);
      total_recv_contribs += static_cast<std::size_t>(exchange.recv_contrib_counts[static_cast<std::size_t>(rank)]);
    }
    exchange.recv_contribs_flat.resize(total_recv_contribs);

    const int contrib_bytes = static_cast<int>(sizeof(PmPotentialContributionRecord));
    for (int rank = 0; rank < world_size; ++rank) {
      const auto r = static_cast<std::size_t>(rank);
      exchange.send_contrib_counts_bytes[r] = exchange.send_contrib_counts[r] * contrib_bytes;
      exchange.send_contrib_displs_bytes[r] = exchange.send_contrib_displs[r] * contrib_bytes;
      exchange.recv_contrib_counts_bytes[r] = exchange.recv_contrib_counts[r] * contrib_bytes;
      exchange.recv_contrib_displs_bytes[r] = exchange.recv_contrib_displs[r] * contrib_bytes;
    }
    MPI_Alltoallv(
        reinterpret_cast<const std::uint8_t*>(exchange.send_contribs_flat.data()),
        exchange.send_contrib_counts_bytes.data(),
        exchange.send_contrib_displs_bytes.data(),
        MPI_BYTE,
        reinterpret_cast<std::uint8_t*>(exchange.recv_contribs_flat.data()),
        exchange.recv_contrib_counts_bytes.data(),
        exchange.recv_contrib_displs_bytes.data(),
        MPI_BYTE,
        MPI_COMM_WORLD);

    std::fill(potential.begin(), potential.end(), 0.0);
    for (const auto& contribution : exchange.recv_contribs_flat) {
      if (contribution.particle_index >= pos_x.size()) {
        throw std::invalid_argument("PmSolver::interpolatePotential response particle index out of range");
      }
      potential[static_cast<std::size_t>(contribution.particle_index)] += contribution.potential;
    }
#endif
  }

  const auto stop = std::chrono::steady_clock::now();
  if (profile != nullptr) {
    profile->interpolate_ms += std::chrono::duration<double, std::milli>(stop - start).count();
    profile->bytes_moved += bytesForParticles(pos_x.size());
  }
}

void PmSolver::solveForParticles(
    PmGridStorage& grid,
    std::span<const double> pos_x,
    std::span<const double> pos_y,
    std::span<const double> pos_z,
    std::span<const double> mass,
    std::span<double> accel_x,
    std::span<double> accel_y,
    std::span<double> accel_z,
    const PmSolveOptions& options,
    PmProfileEvent* profile) {
  validateOptions(m_shape, options);

  if (options.execution_policy == core::ExecutionPolicy::kCuda) {
#if COSMOSIM_ENABLE_CUDA
    cudaStream_t stream = nullptr;
    if (cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking) != cudaSuccess) {
      throw std::runtime_error("Failed to create CUDA stream for PM solve");
    }

    try {
      const auto copy_h2d_start = std::chrono::steady_clock::now();
      core::DeviceBufferDouble pos_x_device(pos_x.size());
      core::DeviceBufferDouble pos_y_device(pos_y.size());
      core::DeviceBufferDouble pos_z_device(pos_z.size());
      core::DeviceBufferDouble mass_device(mass.size());
      core::DeviceBufferDouble density_device(m_shape.cellCount());

      core::DeviceBufferDouble force_x_device(m_shape.cellCount());
      core::DeviceBufferDouble force_y_device(m_shape.cellCount());
      core::DeviceBufferDouble force_z_device(m_shape.cellCount());
      core::DeviceBufferDouble accel_x_device(accel_x.size());
      core::DeviceBufferDouble accel_y_device(accel_y.size());
      core::DeviceBufferDouble accel_z_device(accel_z.size());

      pos_x_device.copyFromHost(pos_x, stream);
      pos_y_device.copyFromHost(pos_y, stream);
      pos_z_device.copyFromHost(pos_z, stream);
      mass_device.copyFromHost(mass, stream);
      if (cudaStreamSynchronize(stream) != cudaSuccess) {
        throw std::runtime_error("Failed while synchronizing H2D particle copy");
      }
      const auto copy_h2d_stop = std::chrono::steady_clock::now();
      if (profile != nullptr) {
        profile->transfer_h2d_ms += std::chrono::duration<double, std::milli>(copy_h2d_stop - copy_h2d_start).count();
      }

      const auto kernel_assign_start = std::chrono::steady_clock::now();
      pmCudaAssignDensityCic(
          PmCudaAssignLaunch{pos_x.size(), m_shape.nx, m_shape.ny, m_shape.nz, options.box_size_mpc_comoving},
          pos_x_device.data(),
          pos_y_device.data(),
          pos_z_device.data(),
          mass_device.data(),
          density_device.data(),
          stream);
      if (cudaStreamSynchronize(stream) != cudaSuccess) {
        throw std::runtime_error("Failed while synchronizing PM assignment kernel");
      }
      const auto kernel_assign_stop = std::chrono::steady_clock::now();
      if (profile != nullptr) {
        profile->device_kernel_ms +=
            std::chrono::duration<double, std::milli>(kernel_assign_stop - kernel_assign_start).count();
      }

      const auto copy_density_start = std::chrono::steady_clock::now();
      density_device.copyToHost(grid.density(), stream);
      if (cudaStreamSynchronize(stream) != cudaSuccess) {
        throw std::runtime_error("Failed while synchronizing D2H density copy");
      }
      const auto copy_density_stop = std::chrono::steady_clock::now();
      if (profile != nullptr) {
        profile->transfer_d2h_ms += std::chrono::duration<double, std::milli>(copy_density_stop - copy_density_start).count();
      }

      const BoxLengths lengths = effectiveBoxLengths(options);
      const double cell_volume =
          (lengths.lx * lengths.ly * lengths.lz) / static_cast<double>(m_shape.cellCount());
      for (double& density_cell : grid.density()) {
        density_cell /= cell_volume;
      }

      solvePoissonPeriodic(grid, options, profile);

      const auto copy_forces_start = std::chrono::steady_clock::now();
      force_x_device.copyFromHost(grid.force_x(), stream);
      force_y_device.copyFromHost(grid.force_y(), stream);
      force_z_device.copyFromHost(grid.force_z(), stream);
      if (cudaStreamSynchronize(stream) != cudaSuccess) {
        throw std::runtime_error("Failed while synchronizing H2D force copy");
      }
      const auto copy_forces_stop = std::chrono::steady_clock::now();
      if (profile != nullptr) {
        profile->transfer_h2d_ms += std::chrono::duration<double, std::milli>(copy_forces_stop - copy_forces_start).count();
      }

      const auto kernel_interp_start = std::chrono::steady_clock::now();
      pmCudaInterpolateForcesCic(
          PmCudaInterpLaunch{pos_x.size(), m_shape.nx, m_shape.ny, m_shape.nz, options.box_size_mpc_comoving},
          pos_x_device.data(),
          pos_y_device.data(),
          pos_z_device.data(),
          force_x_device.data(),
          force_y_device.data(),
          force_z_device.data(),
          accel_x_device.data(),
          accel_y_device.data(),
          accel_z_device.data(),
          stream);
      if (cudaStreamSynchronize(stream) != cudaSuccess) {
        throw std::runtime_error("Failed while synchronizing PM interpolation kernel");
      }
      const auto kernel_interp_stop = std::chrono::steady_clock::now();
      if (profile != nullptr) {
        profile->device_kernel_ms +=
            std::chrono::duration<double, std::milli>(kernel_interp_stop - kernel_interp_start).count();
      }

      const auto copy_accel_start = std::chrono::steady_clock::now();
      accel_x_device.copyToHost(accel_x, stream);
      accel_y_device.copyToHost(accel_y, stream);
      accel_z_device.copyToHost(accel_z, stream);
      if (cudaStreamSynchronize(stream) != cudaSuccess) {
        throw std::runtime_error("Failed while synchronizing D2H acceleration copy");
      }
      const auto copy_accel_stop = std::chrono::steady_clock::now();
      if (profile != nullptr) {
        profile->transfer_d2h_ms += std::chrono::duration<double, std::milli>(copy_accel_stop - copy_accel_start).count();
      }
    } catch (...) {
      (void)cudaStreamDestroy(stream);
      throw;
    }

    (void)cudaStreamDestroy(stream);
    if (profile != nullptr) {
      profile->bytes_moved += bytesForParticles(pos_x.size()) * 2U;
      profile->bytes_moved += bytesForGridSweep(m_shape.cellCount()) * 4U;
    }
    return;
#else
    throw std::runtime_error("PM solve requested execution_policy=cuda, but this build has COSMOSIM_ENABLE_CUDA=OFF");
#endif
  }

  assignDensity(grid, pos_x, pos_y, pos_z, mass, options, profile);
  if (options.boundary_condition == PmBoundaryCondition::kPeriodic) {
    solvePoissonPeriodic(grid, options, profile);
  } else {
    solvePoissonIsolatedOpen(grid, options, profile);
  }
  interpolateForces(grid, pos_x, pos_y, pos_z, accel_x, accel_y, accel_z, options, profile);
}

bool PmSolver::fftBackendAvailable() {
#if COSMOSIM_ENABLE_FFTW
  return true;
#else
  return false;
#endif
}

bool PmSolver::cudaBackendAvailable() {
#if COSMOSIM_ENABLE_CUDA
  int device_count = 0;
  return cudaGetDeviceCount(&device_count) == cudaSuccess && device_count > 0;
#else
  return false;
#endif
}

std::string PmSolver::fftBackendName() {
#if COSMOSIM_ENABLE_FFTW
  return "fftw";
#else
  return "naive_dft";
#endif
}

std::size_t PmSolver::cachedPlanCount() const {
  return m_impl->planCount();
}

std::size_t PmSolver::planBuildCount() const {
  return m_impl->planBuildCount();
}

bool treePmSupportedByBuild() {
  return true;
}

void requireTreePmSupportOrThrow(core::GravitySolver gravity_solver) {
  static_cast<void>(gravity_solver);
}

}  // namespace cosmosim::gravity
