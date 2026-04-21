#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

#include "cosmosim/analysis/diagnostics.hpp"

namespace {

void requireOrThrow(bool condition, const std::string& message) {
  if (!condition) {
    throw std::runtime_error(message);
  }
}

cosmosim::core::SimulationConfig makeConfig(int mesh_n) {
  cosmosim::core::SimulationConfig config;
  config.output.run_name = "validation_phase3_power_spectrum";
  config.output.output_directory = "validation_outputs";
  config.cosmology.box_size_mpc_comoving = 1.0;
  config.analysis.power_spectrum_mesh_n = mesh_n;
  config.analysis.power_spectrum_bin_count = 4;
  return config;
}

cosmosim::core::SimulationState makeDeterministicState() {
  cosmosim::core::SimulationState state;
  constexpr std::size_t count = 256;
  state.resizeParticles(count);
  state.resizeCells(count);
  state.species.count_by_species = {count, 0, 0, 0, 0};

  for (std::size_t i = 0; i < count; ++i) {
    const double phase = static_cast<double>(i + 1U);
    const double x = std::fmod((37.0 * phase + 3.0) * 0.013 + 0.03 * std::sin(0.17 * phase), 1.0);
    const double y = std::fmod((53.0 * phase + 5.0) * 0.011 + 0.02 * std::cos(0.11 * phase), 1.0);
    const double z = std::fmod((61.0 * phase + 7.0) * 0.009 + 0.02 * std::sin(0.07 * phase + 0.3), 1.0);
    state.particle_sidecar.particle_id[i] = 1000U + static_cast<std::uint64_t>(i);
    state.particle_sidecar.species_tag[i] = static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kDarkMatter);
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

double relativeL2(
    std::span<const cosmosim::analysis::PowerSpectrumBin> coarse,
    std::span<const cosmosim::analysis::PowerSpectrumBin> fine) {
  const std::size_t n = std::min(coarse.size(), fine.size());
  double diff2 = 0.0;
  double ref2 = 0.0;
  for (std::size_t i = 0; i < n; ++i) {
    const double dp = coarse[i].power_code_volume - fine[i].power_code_volume;
    diff2 += dp * dp;
    ref2 += fine[i].power_code_volume * fine[i].power_code_volume;
  }
  return std::sqrt(diff2 / std::max(ref2, 1.0e-30));
}

void writeSpectrumJson(
    const std::filesystem::path& path,
    std::span<const cosmosim::analysis::PowerSpectrumBin> bins) {
  std::filesystem::create_directories(path.parent_path());
  std::ofstream out(path);
  out << "{\n  \"power_spectrum\": [\n";
  for (std::size_t i = 0; i < bins.size(); ++i) {
    const auto& bin = bins[i];
    out << "    {\"k_code\": " << bin.k_center_code
        << ", \"p_code_volume\": " << bin.power_code_volume
        << ", \"mode_count\": " << bin.mode_count << "}";
    if (i + 1U < bins.size()) {
      out << ',';
    }
    out << "\n";
  }
  out << "  ]\n}\n";
}

}  // namespace

int main() {
  const cosmosim::core::SimulationState state = makeDeterministicState();

  cosmosim::analysis::DiagnosticsEngine coarse_engine(makeConfig(8));
  cosmosim::analysis::DiagnosticsEngine fine_engine(makeConfig(12));
  const auto coarse = coarse_engine.computePowerSpectrum(state, 8, 4);
  const auto fine = fine_engine.computePowerSpectrum(state, 12, 4);
  requireOrThrow(!coarse.empty(), "coarse power spectrum must not be empty");
  requireOrThrow(!fine.empty(), "fine power spectrum must not be empty");

  const double rel_l2 = relativeL2(coarse, fine);
  requireOrThrow(std::isfinite(rel_l2), "power spectrum relative L2 must be finite");
  {
    std::ostringstream msg;
    msg << "deterministic power spectrum mesh-consistency exceeded tolerance: rel_l2=" << rel_l2;
    requireOrThrow(rel_l2 <= 0.45, msg.str());
  }

  const std::filesystem::path root = COSMOSIM_SOURCE_DIR;
  writeSpectrumJson(root / "validation" / "reference" / "phase3" / "power_spectrum_low.json", coarse);
  writeSpectrumJson(root / "validation" / "reference" / "phase3" / "power_spectrum_high.json", fine);
  return 0;
}
