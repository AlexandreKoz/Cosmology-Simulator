#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "cosmosim/analysis/diagnostics.hpp"
#include "validation_tolerance.hpp"

namespace {

void requireOrThrow(bool condition, const std::string& message) {
  if (!condition) {
    throw std::runtime_error(message);
  }
}

cosmosim::core::SimulationConfig makeConfig() {
  auto config = cosmosim::core::makeUnvalidatedSimulationConfigForTests();
  config.output.run_name = "validation_power_spectrum_mesh_consistency";
  config.cosmology.box_size_mpc_comoving = 1.0;
  return config;
}

cosmosim::core::SimulationState makeDeterministicState() {
  cosmosim::core::SimulationState state;
  constexpr std::size_t PARTICLE_COUNT = 256;
  state.resizeParticles(PARTICLE_COUNT);
  state.resizeCells(PARTICLE_COUNT);
  state.species.count_by_species = {PARTICLE_COUNT, 0, 0, 0, 0};

  for (std::size_t i = 0; i < PARTICLE_COUNT; ++i) {
    const double phase = static_cast<double>(i + 1U);
    const double x = std::fmod((37.0 * phase + 3.0) * 0.013 + 0.03 * std::sin(0.17 * phase), 1.0);
    const double y = std::fmod((53.0 * phase + 5.0) * 0.011 + 0.02 * std::cos(0.11 * phase), 1.0);
    const double z = std::fmod((61.0 * phase + 7.0) * 0.009 + 0.02 * std::sin(0.07 * phase + 0.3), 1.0);
    state.particle_sidecar.particle_id[i] = 1000U + static_cast<std::uint64_t>(i);
    state.particle_sidecar.species_tag[i] =
        static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kDarkMatter);
    state.particles.mass_code[i] = 0.9 + 0.03 * static_cast<double>(i % 7U);
    state.particles.position_x_comoving[i] = (x < 0.0) ? x + 1.0 : x;
    state.particles.position_y_comoving[i] = (y < 0.0) ? y + 1.0 : y;
    state.particles.position_z_comoving[i] = (z < 0.0) ? z + 1.0 : z;
    state.cells.center_x_comoving[i] = state.particles.position_x_comoving[i];
    state.cells.center_y_comoving[i] = state.particles.position_y_comoving[i];
    state.cells.center_z_comoving[i] = state.particles.position_z_comoving[i];
    state.cells.mass_code[i] = 1.0;
    state.gas_cells.density_code[i] = 1.0;
  }

  state.rebuildSpeciesIndex();
  return state;
}

[[nodiscard]] bool nearlyEqual(double lhs, double rhs, double relative_tolerance = 1.0e-12) {
  const double scale = std::max({1.0, std::abs(lhs), std::abs(rhs)});
  return std::abs(lhs - rhs) <= relative_tolerance * scale;
}

struct CommonKComparison {
  double relative_l2 = 0.0;
  std::size_t compared_bin_count = 0;
};

[[nodiscard]] CommonKComparison compareCommonPhysicalKBins(
    const cosmosim::analysis::PowerSpectrumEstimate& coarse,
    const cosmosim::analysis::PowerSpectrumEstimate& fine,
    double max_shared_k_code) {
  requireOrThrow(
      nearlyEqual(coarse.box_size_code, fine.box_size_code),
      "power-spectrum comparison requires identical box-size units");
  requireOrThrow(
      nearlyEqual(coarse.k_fundamental_code, fine.k_fundamental_code),
      "power-spectrum comparison requires identical fundamental k");
  requireOrThrow(
      coarse.options.mass_assignment == fine.options.mass_assignment &&
          coarse.options.window_correction == fine.options.window_correction &&
          coarse.options.shot_noise_policy == fine.options.shot_noise_policy,
      "power-spectrum comparison requires identical assignment/window/noise conventions");

  double difference_norm2 = 0.0;
  double reference_norm2 = 0.0;
  std::size_t compared_bin_count = 0;
  for (const auto& coarse_bin : coarse.bins) {
    if (coarse_bin.k_upper_code > max_shared_k_code) {
      break;
    }
    if (coarse_bin.bin_index >= fine.bins.size()) {
      throw std::runtime_error("fine power spectrum lacks a corresponding common-k bin");
    }
    const auto& fine_bin = fine.bins[coarse_bin.bin_index];

    // Array position is never accepted as proof of physical equivalence.  The
    // bin edges, measured k center, and mode count must all describe the same
    // Fourier shell before the powers are subtracted.
    requireOrThrow(
        nearlyEqual(coarse_bin.k_lower_code, fine_bin.k_lower_code) &&
            nearlyEqual(coarse_bin.k_upper_code, fine_bin.k_upper_code),
        "power-spectrum mesh comparison encountered incompatible physical k-bin edges");
    if (coarse_bin.empty || fine_bin.empty) {
      requireOrThrow(
          coarse_bin.empty == fine_bin.empty,
          "power-spectrum common-k bin is empty on only one mesh");
      continue;
    }
    requireOrThrow(
        nearlyEqual(coarse_bin.k_center_code, fine_bin.k_center_code),
        "power-spectrum mesh comparison encountered mismatched physical k centers");
    requireOrThrow(
        coarse_bin.mode_count == fine_bin.mode_count,
        "power-spectrum common-k bin contains a different Fourier-mode set");

    const double difference =
        coarse_bin.power_code_volume - fine_bin.power_code_volume;
    difference_norm2 += difference * difference;
    reference_norm2 += fine_bin.power_code_volume * fine_bin.power_code_volume;
    ++compared_bin_count;
  }

  requireOrThrow(
      compared_bin_count >= 2U,
      "power-spectrum mesh comparison needs at least two populated shared-k bins");
  requireOrThrow(reference_norm2 > 0.0, "power-spectrum shared-k reference norm is zero");
  return CommonKComparison{
      .relative_l2 = std::sqrt(difference_norm2 / reference_norm2),
      .compared_bin_count = compared_bin_count,
  };
}

void testPowerSpectrumMeshConsistency(
    const cosmosim::validation::ValidationToleranceTable& tolerances) {
  const cosmosim::core::SimulationState state = makeDeterministicState();
  cosmosim::analysis::DiagnosticsEngine diagnostics(makeConfig());

  // mesh_n/bin_count is scaled together so both estimates have the same fixed
  // physical bin width.  Only the conservative low-k domain below 75% of the
  // coarse axis Nyquist is compared; high-k modes resolved only by the fine
  // mesh are intentionally excluded.
  const cosmosim::analysis::PowerSpectrumEstimateOptions coarse_options{
      .mesh_n = 8,
      .bin_count = 12,
      .mass_assignment = cosmosim::analysis::PowerSpectrumMassAssignment::kCloudInCell,
      .window_correction =
          cosmosim::analysis::PowerSpectrumWindowCorrection::kDeconvolveAssignmentWindow,
      .shot_noise_policy =
          cosmosim::analysis::PowerSpectrumShotNoisePolicy::kReportWithoutSubtraction,
  };
  const cosmosim::analysis::PowerSpectrumEstimateOptions fine_options{
      .mesh_n = 16,
      .bin_count = 24,
      .mass_assignment = coarse_options.mass_assignment,
      .window_correction = coarse_options.window_correction,
      .shot_noise_policy = coarse_options.shot_noise_policy,
  };

  const auto coarse = diagnostics.computePowerSpectrumEstimate(state, coarse_options);
  const auto fine = diagnostics.computePowerSpectrumEstimate(state, fine_options);
  requireOrThrow(!coarse.bins.empty(), "coarse power spectrum must not be empty");
  requireOrThrow(!fine.bins.empty(), "fine power spectrum must not be empty");
  requireOrThrow(
      nearlyEqual(coarse.bins.front().k_upper_code, fine.bins.front().k_upper_code),
      "power-spectrum test setup failed to construct common physical k-bin widths");

  const double max_shared_k_code = 0.75 * coarse.k_axis_nyquist_code;
  const CommonKComparison comparison =
      compareCommonPhysicalKBins(coarse, fine, max_shared_k_code);
  requireOrThrow(std::isfinite(comparison.relative_l2), "power-spectrum relative L2 must be finite");

  std::ostringstream message;
  message << "power-spectrum common-k mesh consistency exceeded tolerance: rel_l2="
          << comparison.relative_l2
          << ", compared_bins=" << comparison.compared_bin_count
          << ", max_shared_k_code=" << max_shared_k_code;
  requireOrThrow(
      comparison.relative_l2 <= tolerances.require(
          "analysis_power_spectrum_mesh_consistency.max_relative_l2_error"),
      message.str());

  // Negative control: equal array indices with incompatible physical k must be
  // rejected before any power subtraction can occur.
  auto mismatched_fine = fine;
  std::size_t altered_bin = 0U;
  while (altered_bin < mismatched_fine.bins.size() && mismatched_fine.bins[altered_bin].empty) {
    ++altered_bin;
  }
  requireOrThrow(altered_bin < mismatched_fine.bins.size(), "negative-control spectrum has no populated bin");
  mismatched_fine.bins[altered_bin].k_lower_code += 0.125 * fine.k_fundamental_code;
  bool mismatch_rejected = false;
  try {
    static_cast<void>(compareCommonPhysicalKBins(coarse, mismatched_fine, max_shared_k_code));
  } catch (const std::runtime_error&) {
    mismatch_rejected = true;
  }
  requireOrThrow(
      mismatch_rejected,
      "power-spectrum negative control failed: incompatible k coordinates were compared by index");
}

}  // namespace

int main() {
  const auto tolerances = cosmosim::validation::ValidationToleranceTable::loadFromFile(
      std::filesystem::path(COSMOSIM_SOURCE_DIR) /
      "validation/reference/validation_tolerances_v1.txt");
  testPowerSpectrumMeshConsistency(tolerances);
  return 0;
}
