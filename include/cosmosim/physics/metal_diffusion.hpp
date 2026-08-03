#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include "cosmosim/core/config.hpp"

namespace cosmosim::physics {

// The diffusion operator advances authoritative metal masses while gas masses
// remain fixed. Geometry and all diffusivities use one consistent physical or
// code-unit system supplied by the caller; therefore kappa has L^2/T and face
// distances, areas, volumes, and dt must use the matching L and T units.
struct MetalDiffusionConfig {
  bool enabled = false;
  core::MetalDiffusionModel model = core::MetalDiffusionModel::kNone;
  core::MetalDiffusionTimeIntegrator time_integrator =
      core::MetalDiffusionTimeIntegrator::kExplicitSubcycling;
  double smagorinsky_coefficient = 0.05;
  double parabolic_cfl = 0.4;
  std::uint32_t max_subcycles = 128;
  std::uint32_t max_rkl_stages = 64;
  double diffusivity_floor_code = 0.0;
  double diffusivity_ceiling_code = 1.0e30;
};

struct MetalDiffusionVelocityGradient {
  // grad[velocity component][coordinate derivative], in inverse time.
  std::array<std::array<double, 3>, 3> grad{};
};

struct MetalDiffusionCell {
  std::uint64_t gas_cell_id = 0;
  double gas_mass_code = 0.0;
  double metal_mass_code = 0.0;
  double density_code = 0.0;
  double volume_code = 0.0;
  double filter_length_code = 0.0;
  MetalDiffusionVelocityGradient velocity_gradient;
  bool is_owned_leaf = true;
};

enum class MetalDiffusionBoundaryKind : std::uint8_t {
  kInternal = 0,
  kPeriodic = 1,
  kReflective = 2,
  kOpen = 3,
};

struct MetalDiffusionFace {
  std::uint32_t left_cell = 0;
  std::uint32_t right_cell = 0;
  double area_code = 0.0;
  double center_distance_code = 0.0;
  MetalDiffusionBoundaryKind boundary_kind = MetalDiffusionBoundaryKind::kInternal;
};

struct MetalDiffusionStepReport {
  std::uint32_t subcycles = 0;
  std::uint32_t rkl_stages = 0;
  std::uint64_t faces_evaluated = 0;
  std::uint64_t limited_faces = 0;
  double stable_dt_code = 0.0;
  double metal_mass_before_code = 0.0;
  double metal_mass_after_code = 0.0;
  double open_boundary_loss_code = 0.0;
  double conservation_residual_code = 0.0;
  double minimum_metallicity = 0.0;
  double maximum_metallicity = 0.0;
};

[[nodiscard]] double traceFreeStrainMagnitude(
    const MetalDiffusionVelocityGradient& gradient) noexcept;

[[nodiscard]] double smagorinskyMetalDiffusivityCode(
    const MetalDiffusionConfig& config,
    const MetalDiffusionCell& cell) noexcept;

class MetalDiffusionModel {
 public:
  explicit MetalDiffusionModel(MetalDiffusionConfig config);

  [[nodiscard]] const MetalDiffusionConfig& config() const noexcept;

  [[nodiscard]] double stableExplicitTimestepCode(
      std::span<const MetalDiffusionCell> cells,
      std::span<const MetalDiffusionFace> faces) const;

  // Mutates only MetalDiffusionCell::metal_mass_code. Pairwise face updates are
  // equal and opposite. A symmetric donor/receiver limiter enforces
  // 0 <= M_Z <= M_gas without non-conservative clipping.
  [[nodiscard]] MetalDiffusionStepReport advance(
      std::span<MetalDiffusionCell> cells,
      std::span<const MetalDiffusionFace> faces,
      double dt_code) const;

 private:
  [[nodiscard]] std::vector<double> evaluateDerivative(
      std::span<const MetalDiffusionCell> cells,
      std::span<const MetalDiffusionFace> faces,
      std::span<const double> metal_mass_code,
      std::uint64_t* faces_evaluated) const;

  void conservativeLimitedEuler(
      std::span<const MetalDiffusionCell> cells,
      std::span<const MetalDiffusionFace> faces,
      std::span<double> metal_mass_code,
      double dt_code,
      std::uint64_t* faces_evaluated,
      std::uint64_t* limited_faces) const;

  void explicitSspRk2(
      std::span<const MetalDiffusionCell> cells,
      std::span<const MetalDiffusionFace> faces,
      std::span<double> metal_mass_code,
      double dt_code,
      std::uint64_t* faces_evaluated,
      std::uint64_t* limited_faces) const;

  void rkl2(
      std::span<const MetalDiffusionCell> cells,
      std::span<const MetalDiffusionFace> faces,
      std::span<double> metal_mass_code,
      double dt_code,
      std::uint32_t stages,
      std::uint64_t* faces_evaluated) const;

  MetalDiffusionConfig m_config;
};

[[nodiscard]] MetalDiffusionConfig makeMetalDiffusionConfig(
    const core::PhysicsConfig& physics_config);

}  // namespace cosmosim::physics
