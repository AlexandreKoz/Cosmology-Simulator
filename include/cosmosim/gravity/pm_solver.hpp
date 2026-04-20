#pragma once

#include <complex>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include "cosmosim/core/config.hpp"
#include "cosmosim/core/execution_policy.hpp"
#include "cosmosim/parallel/distributed_memory.hpp"

namespace cosmosim::gravity {

enum class PmAssignmentScheme {
  kCic,
  kTsc,
};

enum class PmBoundaryCondition {
  kPeriodic,
  kIsolatedOpen,
};

struct PmGridShape {
  std::size_t nx = 0;
  std::size_t ny = 0;
  std::size_t nz = 0;

  [[nodiscard]] std::size_t cellCount() const;
  [[nodiscard]] bool isValid() const;
};

enum class PmDataResidencyPolicy {
  kHostOnly,
  kPreferDevice,
};

struct PmSolveOptions {
  double box_size_x_mpc_comoving = 0.0;
  double box_size_y_mpc_comoving = 0.0;
  double box_size_z_mpc_comoving = 0.0;
  // Legacy scalar compatibility lane; if axis-aware lengths are zero, this value is used for all axes.
  double box_size_mpc_comoving = 0.0;
  double scale_factor = 1.0;
  double gravitational_constant_code = 1.0;
  PmAssignmentScheme assignment_scheme = PmAssignmentScheme::kCic;
  bool enable_window_deconvolution = false;
  core::ExecutionPolicy execution_policy = core::ExecutionPolicy::kHostSerial;
  PmDataResidencyPolicy data_residency = PmDataResidencyPolicy::kHostOnly;
  core::PmDecompositionMode decomposition_mode = core::PmDecompositionMode::kSlab;
  PmBoundaryCondition boundary_condition = PmBoundaryCondition::kPeriodic;
  // Optional TreePM long-range Gaussian split scale. <=0 disables filtering.
  double tree_pm_split_scale_comoving = 0.0;
};

struct PmProfileEvent {
  std::uint64_t bytes_moved = 0;
  double assign_ms = 0.0;
  double fft_forward_ms = 0.0;
  double poisson_ms = 0.0;
  double gradient_ms = 0.0;
  double fft_inverse_ms = 0.0;
  double fft_transpose_ms = 0.0;
  double interpolate_ms = 0.0;
  double transfer_h2d_ms = 0.0;
  double transfer_d2h_ms = 0.0;
  double device_kernel_ms = 0.0;
  std::uint64_t fft_transpose_bytes = 0;
};

class PmProfiler {
 public:
  void reset();
  void append(const PmProfileEvent& event);
  [[nodiscard]] const PmProfileEvent& totals() const;

 private:
  PmProfileEvent m_totals{};
};

class PmGridStorage {
 public:
  explicit PmGridStorage(PmGridShape shape);
  PmGridStorage(PmGridShape shape, parallel::PmSlabLayout layout);

  [[nodiscard]] const PmGridShape& shape() const;
  [[nodiscard]] const parallel::PmSlabLayout& slabLayout() const;
  [[nodiscard]] bool ownsFullDomain() const noexcept;
  [[nodiscard]] std::size_t localCellCount() const noexcept;

  [[nodiscard]] std::span<double> density();
  [[nodiscard]] std::span<const double> density() const;

  [[nodiscard]] std::span<double> potential();
  [[nodiscard]] std::span<const double> potential() const;

  [[nodiscard]] std::span<double> force_x();
  [[nodiscard]] std::span<const double> force_x() const;

  [[nodiscard]] std::span<double> force_y();
  [[nodiscard]] std::span<const double> force_y() const;

  [[nodiscard]] std::span<double> force_z();
  [[nodiscard]] std::span<const double> force_z() const;

  [[nodiscard]] std::size_t linearIndex(std::size_t ix, std::size_t iy, std::size_t iz) const;
  void clear();

 private:
  PmGridShape m_shape;
  parallel::PmSlabLayout m_layout;
  std::vector<double> m_density;
  std::vector<double> m_potential;
  std::vector<double> m_force_x;
  std::vector<double> m_force_y;
  std::vector<double> m_force_z;
};

class PmSolver {
 public:
  explicit PmSolver(PmGridShape shape);
  ~PmSolver();
  PmSolver(PmSolver&&) noexcept;
  PmSolver& operator=(PmSolver&&) noexcept;
  PmSolver(const PmSolver&) = delete;
  PmSolver& operator=(const PmSolver&) = delete;

  [[nodiscard]] const PmGridShape& shape() const;

  void assignDensity(
      PmGridStorage& grid,
      std::span<const double> pos_x,
      std::span<const double> pos_y,
      std::span<const double> pos_z,
      std::span<const double> mass,
      const PmSolveOptions& options,
      PmProfileEvent* profile = nullptr) const;

  // Periodic comoving Poisson solve contract on the mesh:
  //   ∇²φ(x) = 4π G a² [ρ(x) - ρ̄]
  // Fourier convention (for k != 0):
  //   φ_k = -4π G a² δρ_k / k², with δρ = ρ - ρ̄
  //   a_i(k) = -i k_i φ_k
  // Periodic zero mode policy:
  //   φ_{k=0} = 0 and therefore a_{k=0} = 0.
  //
  // After return, grid.potential() and grid.force_{x,y,z}() are populated and
  // available for direct mesh inspection and interpolation.
  void solvePoissonPeriodic(PmGridStorage& grid, const PmSolveOptions& options, PmProfileEvent* profile = nullptr);
  void solvePoissonIsolatedOpen(PmGridStorage& grid, const PmSolveOptions& options, PmProfileEvent* profile = nullptr);

  void interpolateForces(
      const PmGridStorage& grid,
      std::span<const double> pos_x,
      std::span<const double> pos_y,
      std::span<const double> pos_z,
      std::span<double> accel_x,
      std::span<double> accel_y,
      std::span<double> accel_z,
      const PmSolveOptions& options,
      PmProfileEvent* profile = nullptr) const;
  // For slab-distributed meshes, interpolation uses reverse communication:
  // particle owners send weighted stencil requests to slab owners and receive
  // weighted contributions back, then accumulate locally in particle order.

  // Assignment-scheme transpose gather of mesh potential values to particles.
  // This uses the same geometric stencil/convention as interpolateForces.
  void interpolatePotential(
      const PmGridStorage& grid,
      std::span<const double> pos_x,
      std::span<const double> pos_y,
      std::span<const double> pos_z,
      std::span<double> potential,
      const PmSolveOptions& options,
      PmProfileEvent* profile = nullptr) const;

  void solveForParticles(
      PmGridStorage& grid,
      std::span<const double> pos_x,
      std::span<const double> pos_y,
      std::span<const double> pos_z,
      std::span<const double> mass,
      std::span<double> accel_x,
      std::span<double> accel_y,
      std::span<double> accel_z,
      const PmSolveOptions& options,
      PmProfileEvent* profile = nullptr);

  [[nodiscard]] static bool fftBackendAvailable();
  [[nodiscard]] static bool cudaBackendAvailable();
  [[nodiscard]] static std::string fftBackendName();
  [[nodiscard]] std::size_t cachedPlanCount() const;
  [[nodiscard]] std::size_t planBuildCount() const;

 private:
  class Impl;
  PmGridShape m_shape;
  std::unique_ptr<Impl> m_impl;
};

[[nodiscard]] bool treePmSupportedByBuild();
void requireTreePmSupportOrThrow(core::GravitySolver gravity_solver);

}  // namespace cosmosim::gravity
