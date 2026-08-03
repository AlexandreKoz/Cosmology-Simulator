#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <vector>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/physics/metal_diffusion.hpp"

int main() {
  constexpr std::size_t kCellCount = 1U << 17;
  constexpr std::size_t kIterations = 8;
  const double dx = 1.0 / static_cast<double>(kCellCount);
  std::vector<cosmosim::physics::MetalDiffusionCell> cells(kCellCount);
  std::vector<cosmosim::physics::MetalDiffusionFace> faces(kCellCount);
  for (std::size_t i = 0; i < kCellCount; ++i) {
    auto& cell = cells[i];
    cell.gas_cell_id = i + 1U;
    cell.gas_mass_code = dx;
    cell.metal_mass_code = dx * (0.02 + 0.01 * std::sin(6.283185307179586 * i / kCellCount));
    cell.density_code = 1.0;
    cell.volume_code = dx;
    cell.filter_length_code = 1.0;
    cell.velocity_gradient.grad[0][0] = 0.5;
    cell.velocity_gradient.grad[1][1] = -0.5;
    faces[i] = {
        .left_cell = static_cast<std::uint32_t>(i),
        .right_cell = static_cast<std::uint32_t>((i + 1U) % kCellCount),
        .area_code = 1.0,
        .center_distance_code = dx,
        .boundary_kind = i + 1U == kCellCount
            ? cosmosim::physics::MetalDiffusionBoundaryKind::kPeriodic
            : cosmosim::physics::MetalDiffusionBoundaryKind::kInternal,
    };
  }
  cosmosim::physics::MetalDiffusionConfig config;
  config.enabled = true;
  config.model = cosmosim::core::MetalDiffusionModel::kSmagorinsky;
  config.time_integrator = cosmosim::core::MetalDiffusionTimeIntegrator::kExplicitSubcycling;
  config.smagorinsky_coefficient = 1.0e-8;
  config.max_subcycles = 16;
  cosmosim::physics::MetalDiffusionModel model(config);
  const double dt = 0.25 * model.stableExplicitTimestepCode(cells, faces);

  const auto start = std::chrono::steady_clock::now();
  std::uint64_t evaluated = 0;
  for (std::size_t iteration = 0; iteration < kIterations; ++iteration) {
    const auto report = model.advance(cells, faces, dt);
    evaluated += report.faces_evaluated;
  }
  const double seconds = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - start).count();
  std::cout << "bench_metal_diffusion"
            << " build_type=" << COSMOSIM_BUILD_TYPE
            << " cells=" << kCellCount
            << " faces=" << kCellCount
            << " iterations=" << kIterations
            << " seconds=" << seconds
            << " faces_per_s=" << evaluated / std::max(seconds, 1.0e-12)
            << " persistent_scalar_bytes_per_cell=" << sizeof(double)
            << " scratch_bytes_per_cell_approx=" << 5U * sizeof(double)
            << '\n';
}
