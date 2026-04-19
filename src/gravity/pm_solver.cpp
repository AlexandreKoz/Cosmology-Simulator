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
#include <string>
#include <utility>
#include <vector>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/core/device_buffer.hpp"
#include "cosmosim/gravity/tree_pm_split_kernel.hpp"
#if COSMOSIM_ENABLE_CUDA
#include <cuda_runtime.h>
#include "cosmosim/gravity/pm_cuda_kernels.hpp"
#endif

#if COSMOSIM_ENABLE_FFTW
#include <fftw3.h>
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

void validateOptions(const PmGridShape& shape, const PmSolveOptions& options) {
  if (!shape.isValid()) {
    throw std::invalid_argument("PM grid shape must be non-zero in all dimensions");
  }
  if (options.box_size_mpc_comoving <= 0.0) {
    throw std::invalid_argument("PM solve requires box_size_mpc_comoving > 0");
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
}

void validateSingleRankFullDomainGridContract(const PmGridStorage& grid, std::string_view callsite) {
  if (!grid.ownsFullDomain()) {
    throw std::invalid_argument(
        std::string(callsite) +
        " currently requires a full-domain PM slab on this rank; distributed FFT/deposit/gather is not implemented in this phase");
  }
}

}  // namespace

class PmSolver::Impl {
 public:
  explicit Impl(PmGridShape shape)
      : m_shape(shape),
        m_real(shape.cellCount(), 0.0),
        m_fourier(shape.nx * shape.ny * (shape.nz / 2U + 1U), std::complex<double>(0.0, 0.0)),
        m_potential_k(m_shape.nx * m_shape.ny * (m_shape.nz / 2U + 1U), std::complex<double>(0.0, 0.0)),
        m_working_k(m_shape.nx * m_shape.ny * (m_shape.nz / 2U + 1U), std::complex<double>(0.0, 0.0)) {
#if COSMOSIM_ENABLE_FFTW
    m_forward_plan = fftw_plan_dft_r2c_3d(
        static_cast<int>(m_shape.nx),
        static_cast<int>(m_shape.ny),
        static_cast<int>(m_shape.nz),
        m_real.data(),
        reinterpret_cast<fftw_complex*>(m_fourier.data()),
        FFTW_MEASURE);
    m_inverse_plan = fftw_plan_dft_c2r_3d(
        static_cast<int>(m_shape.nx),
        static_cast<int>(m_shape.ny),
        static_cast<int>(m_shape.nz),
        reinterpret_cast<fftw_complex*>(m_fourier.data()),
        m_real.data(),
        FFTW_MEASURE);
    if (m_forward_plan == nullptr || m_inverse_plan == nullptr) {
      throw std::runtime_error("Failed to create FFTW plans for PM solver");
    }
#endif
  }

  ~Impl() {
#if COSMOSIM_ENABLE_FFTW
    if (m_forward_plan != nullptr) {
      fftw_destroy_plan(m_forward_plan);
      m_forward_plan = nullptr;
    }
    if (m_inverse_plan != nullptr) {
      fftw_destroy_plan(m_inverse_plan);
      m_inverse_plan = nullptr;
    }
#endif
  }

  [[nodiscard]] std::span<double> realGrid() { return m_real; }
  [[nodiscard]] std::span<std::complex<double>> fourierGrid() { return m_fourier; }
  [[nodiscard]] std::span<std::complex<double>> potentialScratch() { return m_potential_k; }
  [[nodiscard]] std::span<std::complex<double>> workingScratch() { return m_working_k; }

  double forwardFft() {
    const auto start = std::chrono::steady_clock::now();
#if COSMOSIM_ENABLE_FFTW
    fftw_execute(m_forward_plan);
#else
    naiveForwardDft();
#endif
    const auto stop = std::chrono::steady_clock::now();
    return std::chrono::duration<double, std::milli>(stop - start).count();
  }

  double inverseFft() {
    const auto start = std::chrono::steady_clock::now();
#if COSMOSIM_ENABLE_FFTW
    fftw_execute(m_inverse_plan);
#else
    naiveInverseDft();
#endif
    const auto stop = std::chrono::steady_clock::now();
    return std::chrono::duration<double, std::milli>(stop - start).count();
  }

 private:
#if !COSMOSIM_ENABLE_FFTW
  void naiveForwardDft() {
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
                acc += m_real[rindex] * euler;
              }
            }
          }
          const std::size_t cindex = (ix * m_shape.ny + iy) * nz_complex + iz;
          m_fourier[cindex] = acc;
        }
      }
    }
  }

  void naiveInverseDft() {
    const std::size_t total = m_shape.cellCount();
    std::fill(m_real.begin(), m_real.end(), 0.0);
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
                  acc += m_fourier[cindex] * euler;
                } else {
                  acc += m_fourier[cindex] * euler + std::conj(m_fourier[cindex]) * std::conj(euler);
                }
              }
            }
          }
          m_real[(x * m_shape.ny + y) * m_shape.nz + z] = acc.real() / static_cast<double>(total);
        }
      }
    }
  }
#endif

  PmGridShape m_shape;
  std::vector<double> m_real;
  std::vector<std::complex<double>> m_fourier;
  std::vector<std::complex<double>> m_potential_k;
  std::vector<std::complex<double>> m_working_k;
#if COSMOSIM_ENABLE_FFTW
  fftw_plan m_forward_plan = nullptr;
  fftw_plan m_inverse_plan = nullptr;
#endif
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
  m_totals.interpolate_ms += event.interpolate_ms;
  m_totals.transfer_h2d_ms += event.transfer_h2d_ms;
  m_totals.transfer_d2h_ms += event.transfer_d2h_ms;
  m_totals.device_kernel_ms += event.device_kernel_ms;
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
  validateSingleRankFullDomainGridContract(grid, "PmSolver::assignDensity");
  if (pos_x.size() != pos_y.size() || pos_x.size() != pos_z.size() || pos_x.size() != mass.size()) {
    throw std::invalid_argument("Particle coordinate/mass spans must match in assignDensity");
  }

  const auto start = std::chrono::steady_clock::now();
  std::fill(grid.density().begin(), grid.density().end(), 0.0);

  const double inv_dx = static_cast<double>(m_shape.nx) / options.box_size_mpc_comoving;
  const double inv_dy = static_cast<double>(m_shape.ny) / options.box_size_mpc_comoving;
  const double inv_dz = static_cast<double>(m_shape.nz) / options.box_size_mpc_comoving;

  for (std::size_t p = 0; p < mass.size(); ++p) {
    const double x = wrapPosition(pos_x[p], options.box_size_mpc_comoving) * inv_dx;
    const double y = wrapPosition(pos_y[p], options.box_size_mpc_comoving) * inv_dy;
    const double z = wrapPosition(pos_z[p], options.box_size_mpc_comoving) * inv_dz;

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
          grid.density()[grid.linearIndex(ix, iy, iz)] += mass[p] * weight;
        }
      }
    }
  }

  const double cell_volume = std::pow(options.box_size_mpc_comoving, 3.0) /
      static_cast<double>(m_shape.cellCount());
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
  if (grid.shape().cellCount() != m_shape.cellCount()) {
    throw std::invalid_argument("PM solver/grid shape mismatch in solvePoissonPeriodic");
  }
  validateSingleRankFullDomainGridContract(grid, "PmSolver::solvePoissonPeriodic");

  auto real = m_impl->realGrid();
  std::copy(grid.density().begin(), grid.density().end(), real.begin());

  const double mean_density = std::accumulate(real.begin(), real.end(), 0.0) /
      static_cast<double>(m_shape.cellCount());
  for (double& value : real) {
    value -= mean_density;
  }

  if (profile != nullptr) {
    profile->fft_forward_ms += m_impl->forwardFft();
  } else {
    (void)m_impl->forwardFft();
  }

  const auto poisson_start = std::chrono::steady_clock::now();
  auto fourier = m_impl->fourierGrid();
  const std::size_t nz_complex = m_shape.nz / 2U + 1U;
  const double prefactor = -4.0 * k_pi * options.gravitational_constant_code * options.scale_factor * options.scale_factor;
  const double dkx = 2.0 * k_pi / options.box_size_mpc_comoving;
  const double dky = 2.0 * k_pi / options.box_size_mpc_comoving;
  const double dkz = 2.0 * k_pi / options.box_size_mpc_comoving;

  for (std::size_t ix = 0; ix < m_shape.nx; ++ix) {
    const std::ptrdiff_t nx_mode = ix <= m_shape.nx / 2U ? static_cast<std::ptrdiff_t>(ix)
                                                          : static_cast<std::ptrdiff_t>(ix) - static_cast<std::ptrdiff_t>(m_shape.nx);
    const double kx = dkx * static_cast<double>(nx_mode);
    for (std::size_t iy = 0; iy < m_shape.ny; ++iy) {
      const std::ptrdiff_t ny_mode = iy <= m_shape.ny / 2U ? static_cast<std::ptrdiff_t>(iy)
                                                            : static_cast<std::ptrdiff_t>(iy) - static_cast<std::ptrdiff_t>(m_shape.ny);
      const double ky = dky * static_cast<double>(ny_mode);
      for (std::size_t iz = 0; iz < nz_complex; ++iz) {
        const double kz = dkz * static_cast<double>(iz);
        const std::size_t index = (ix * m_shape.ny + iy) * nz_complex + iz;
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
              sinc(0.5 * kx * options.box_size_mpc_comoving / static_cast<double>(m_shape.nx)),
              static_cast<double>(window_exponent));
          const double wy = std::pow(
              sinc(0.5 * ky * options.box_size_mpc_comoving / static_cast<double>(m_shape.ny)),
              static_cast<double>(window_exponent));
          const double wz = std::pow(
              sinc(0.5 * kz * options.box_size_mpc_comoving / static_cast<double>(m_shape.nz)),
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
    }
  }
  const auto poisson_stop = std::chrono::steady_clock::now();

  if (profile != nullptr) {
    profile->poisson_ms += std::chrono::duration<double, std::milli>(poisson_stop - poisson_start).count();
  }

  auto potential_k = m_impl->potentialScratch();
  std::copy(fourier.begin(), fourier.end(), potential_k.begin());
  auto working_k = m_impl->workingScratch();

  auto inverse_into = [this, profile](std::span<const std::complex<double>> src, std::span<double> dst) {
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
    }
  };

  const auto grad_start = std::chrono::steady_clock::now();
  inverse_into(potential_k, grid.potential());

  auto fill_gradient_k = [&](std::span<std::complex<double>> dst, int axis) {
    for (std::size_t ix = 0; ix < m_shape.nx; ++ix) {
      const std::ptrdiff_t nx_mode = ix <= m_shape.nx / 2U ? static_cast<std::ptrdiff_t>(ix)
                                                            : static_cast<std::ptrdiff_t>(ix) - static_cast<std::ptrdiff_t>(m_shape.nx);
      const double kx = dkx * static_cast<double>(nx_mode);
      for (std::size_t iy = 0; iy < m_shape.ny; ++iy) {
        const std::ptrdiff_t ny_mode = iy <= m_shape.ny / 2U ? static_cast<std::ptrdiff_t>(iy)
                                                              : static_cast<std::ptrdiff_t>(iy) - static_cast<std::ptrdiff_t>(m_shape.ny);
        const double ky = dky * static_cast<double>(ny_mode);
        for (std::size_t iz = 0; iz < nz_complex; ++iz) {
          const double kz = dkz * static_cast<double>(iz);
          const std::size_t index = (ix * m_shape.ny + iy) * nz_complex + iz;
          const double component_k = axis == 0 ? kx : (axis == 1 ? ky : kz);
          dst[index] = std::complex<double>(0.0, -component_k) * potential_k[index];
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
  validateSingleRankFullDomainGridContract(grid, "PmSolver::interpolateForces");

  const auto start = std::chrono::steady_clock::now();

  const double inv_dx = static_cast<double>(m_shape.nx) / options.box_size_mpc_comoving;
  const double inv_dy = static_cast<double>(m_shape.ny) / options.box_size_mpc_comoving;
  const double inv_dz = static_cast<double>(m_shape.nz) / options.box_size_mpc_comoving;

  for (std::size_t p = 0; p < pos_x.size(); ++p) {
    const double x = wrapPosition(pos_x[p], options.box_size_mpc_comoving) * inv_dx;
    const double y = wrapPosition(pos_y[p], options.box_size_mpc_comoving) * inv_dy;
    const double z = wrapPosition(pos_z[p], options.box_size_mpc_comoving) * inv_dz;

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
  validateSingleRankFullDomainGridContract(grid, "PmSolver::interpolatePotential");

  const auto start = std::chrono::steady_clock::now();
  const double inv_dx = static_cast<double>(m_shape.nx) / options.box_size_mpc_comoving;
  const double inv_dy = static_cast<double>(m_shape.ny) / options.box_size_mpc_comoving;
  const double inv_dz = static_cast<double>(m_shape.nz) / options.box_size_mpc_comoving;

  for (std::size_t p = 0; p < pos_x.size(); ++p) {
    const double x = wrapPosition(pos_x[p], options.box_size_mpc_comoving) * inv_dx;
    const double y = wrapPosition(pos_y[p], options.box_size_mpc_comoving) * inv_dy;
    const double z = wrapPosition(pos_z[p], options.box_size_mpc_comoving) * inv_dz;

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

      const double cell_volume = std::pow(options.box_size_mpc_comoving, 3.0) / static_cast<double>(m_shape.cellCount());
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
  solvePoissonPeriodic(grid, options, profile);
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

bool treePmSupportedByBuild() {
  return true;
}

void requireTreePmSupportOrThrow(core::GravitySolver gravity_solver) {
  static_cast<void>(gravity_solver);
}

}  // namespace cosmosim::gravity
