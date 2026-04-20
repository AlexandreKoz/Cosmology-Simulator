#include "cosmosim/analysis/diagnostics.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <numbers>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace cosmosim::analysis {
namespace {

constexpr double k_two_pi = 2.0 * std::numbers::pi_v<double>;

[[nodiscard]] std::size_t flatten(std::size_t ix, std::size_t iy, std::size_t iz, std::size_t n) {
  return (ix * n + iy) * n + iz;
}

[[nodiscard]] std::size_t flatten2(std::size_t ix, std::size_t iy, std::size_t n) { return ix * n + iy; }

[[nodiscard]] bool finite3(double x, double y, double z) {
  return std::isfinite(x) && std::isfinite(y) && std::isfinite(z);
}

[[nodiscard]] double magnitude(const std::array<double, 3>& v) {
  return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

[[nodiscard]] const char* diagnosticClassLabel(DiagnosticClass cls) {
  switch (cls) {
    case DiagnosticClass::kRunHealth:
      return "run_health";
    case DiagnosticClass::kScienceLight:
      return "science_light";
    case DiagnosticClass::kScienceHeavy:
      return "science_heavy";
  }
  throw std::logic_error("unhandled DiagnosticClass enum value during serialization");
}

[[nodiscard]] const char* diagnosticTierLabel(DiagnosticTier tier) {
  switch (tier) {
    case DiagnosticTier::kInfrastructureHealth:
      return "infrastructure_health";
    case DiagnosticTier::kValidatedScience:
      return "validated_science";
    case DiagnosticTier::kReferenceScience:
      return "reference_science";
  }
  throw std::logic_error("unhandled DiagnosticTier enum value during serialization");
}

[[nodiscard]] const char* diagnosticMaturityLabel(DiagnosticMaturity maturity) {
  switch (maturity) {
    case DiagnosticMaturity::kProduction:
      return "production";
    case DiagnosticMaturity::kValidated:
      return "validated";
    case DiagnosticMaturity::kProvisional:
      return "provisional";
  }
  throw std::logic_error("unhandled DiagnosticMaturity enum value during serialization");
}

[[nodiscard]] const char* diagnosticScalabilityLabel(DiagnosticScalability scalability) {
  switch (scalability) {
    case DiagnosticScalability::kCheap:
      return "cheap";
    case DiagnosticScalability::kModerate:
      return "moderate";
    case DiagnosticScalability::kHeavyReference:
      return "heavy_reference";
  }
  throw std::logic_error("unhandled DiagnosticScalability enum value during serialization");
}

[[nodiscard]] const char* executionPolicyLabel(core::AnalysisConfig::DiagnosticsExecutionPolicy policy) {
  switch (policy) {
    case core::AnalysisConfig::DiagnosticsExecutionPolicy::kRunHealthOnly:
      return "run_health_only";
    case core::AnalysisConfig::DiagnosticsExecutionPolicy::kRunHealthAndLightScience:
      return "run_health_and_light_science";
    case core::AnalysisConfig::DiagnosticsExecutionPolicy::kAllIncludingProvisional:
      return "all_including_provisional";
  }
  throw std::logic_error("unhandled DiagnosticsExecutionPolicy enum value during serialization");
}

}  // namespace

DiagnosticsEngine::DiagnosticsEngine(core::SimulationConfig config) : m_config(std::move(config)) {}

RunHealthCounters DiagnosticsEngine::computeRunHealth(const core::SimulationState& state) const {
  RunHealthCounters counters;
  counters.particle_count = static_cast<std::uint64_t>(state.particles.size());
  counters.cell_count = static_cast<std::uint64_t>(state.cells.size());
  counters.star_count = static_cast<std::uint64_t>(state.star_particles.size());
  counters.ownership_invariants_ok = state.validateOwnershipInvariants();
  counters.unique_particle_ids_ok = state.validateUniqueParticleIds();
  counters.gravity_softening_sidecar_size_ok =
      state.particle_sidecar.gravity_softening_comoving.empty() ||
      state.particle_sidecar.gravity_softening_comoving.size() == state.particles.size();

  for (std::size_t i = 0; i < state.particles.size(); ++i) {
    if (!finite3(
            state.particles.position_x_comoving[i],
            state.particles.position_y_comoving[i],
            state.particles.position_z_comoving[i]) ||
        !finite3(
            state.particles.velocity_x_peculiar[i],
            state.particles.velocity_y_peculiar[i],
            state.particles.velocity_z_peculiar[i]) ||
        !std::isfinite(state.particles.mass_code[i])) {
      ++counters.non_finite_particles;
    }
    if (state.particles.mass_code[i] <= 0.0) {
      ++counters.non_positive_particle_mass;
    }
  }

  for (std::size_t i = 0; i < state.cells.size(); ++i) {
    if (!finite3(state.cells.center_x_comoving[i], state.cells.center_y_comoving[i], state.cells.center_z_comoving[i]) ||
        !std::isfinite(state.cells.mass_code[i]) || !std::isfinite(state.gas_cells.density_code[i])) {
      ++counters.non_finite_cells;
    }
  }
  for (std::size_t i = 0; i < state.particle_sidecar.gravity_softening_comoving.size(); ++i) {
    const double softening = state.particle_sidecar.gravity_softening_comoving[i];
    if (!std::isfinite(softening) || softening < 0.0) {
      ++counters.non_finite_gravity_softening;
    }
  }

  return counters;
}

std::vector<PowerSpectrumBin> DiagnosticsEngine::computePowerSpectrum(
    const core::SimulationState& state,
    std::size_t mesh_n,
    std::size_t bin_count) const {
  if (mesh_n == 0 || bin_count == 0) {
    throw std::invalid_argument("power spectrum requires mesh_n > 0 and bin_count > 0");
  }

  const double box_size_mpc_comov = m_config.cosmology.box_size_mpc_comoving;
  const double cell_size = box_size_mpc_comov / static_cast<double>(mesh_n);
  const double cell_volume = cell_size * cell_size * cell_size;
  const std::size_t mesh_count = mesh_n * mesh_n * mesh_n;

  std::vector<double> mass_density(mesh_count, 0.0);
  double total_mass = 0.0;
  for (std::size_t p = 0; p < state.particles.size(); ++p) {
    const double wrapped_x = std::fmod(std::fmod(state.particles.position_x_comoving[p], box_size_mpc_comov) + box_size_mpc_comov,
                                       box_size_mpc_comov);
    const double wrapped_y = std::fmod(std::fmod(state.particles.position_y_comoving[p], box_size_mpc_comov) + box_size_mpc_comov,
                                       box_size_mpc_comov);
    const double wrapped_z = std::fmod(std::fmod(state.particles.position_z_comoving[p], box_size_mpc_comov) + box_size_mpc_comov,
                                       box_size_mpc_comov);

    const std::size_t ix = std::min(mesh_n - 1, static_cast<std::size_t>(wrapped_x / cell_size));
    const std::size_t iy = std::min(mesh_n - 1, static_cast<std::size_t>(wrapped_y / cell_size));
    const std::size_t iz = std::min(mesh_n - 1, static_cast<std::size_t>(wrapped_z / cell_size));

    const double mass = state.particles.mass_code[p];
    mass_density[flatten(ix, iy, iz, mesh_n)] += mass / cell_volume;
    total_mass += mass;
  }

  const double mean_density = (mesh_count > 0) ? total_mass / (box_size_mpc_comov * box_size_mpc_comov * box_size_mpc_comov)
                                                : 0.0;
  if (mean_density <= 0.0) {
    return {};
  }

  std::vector<double> delta(mesh_count, 0.0);
  for (std::size_t i = 0; i < mesh_count; ++i) {
    delta[i] = mass_density[i] / mean_density - 1.0;
  }

  const std::size_t nyquist_mode = mesh_n / 2;
  const double k_fundamental = k_two_pi / box_size_mpc_comov;
  const double k_max = std::sqrt(3.0) * static_cast<double>(nyquist_mode) * k_fundamental;
  const double bin_width = k_max / static_cast<double>(bin_count);

  std::vector<double> power_sum(bin_count, 0.0);
  std::vector<double> k_sum(bin_count, 0.0);
  std::vector<std::uint64_t> mode_count(bin_count, 0);

  for (std::size_t nx = 0; nx < mesh_n; ++nx) {
    const int kx_int = (nx <= nyquist_mode) ? static_cast<int>(nx) : static_cast<int>(nx) - static_cast<int>(mesh_n);
    for (std::size_t ny = 0; ny < mesh_n; ++ny) {
      const int ky_int = (ny <= nyquist_mode) ? static_cast<int>(ny) : static_cast<int>(ny) - static_cast<int>(mesh_n);
      for (std::size_t nz = 0; nz < mesh_n; ++nz) {
        const int kz_int = (nz <= nyquist_mode) ? static_cast<int>(nz) : static_cast<int>(nz) - static_cast<int>(mesh_n);

        if (kx_int == 0 && ky_int == 0 && kz_int == 0) {
          continue;
        }

        const double k_mag =
            std::sqrt(static_cast<double>(kx_int * kx_int + ky_int * ky_int + kz_int * kz_int)) * k_fundamental;
        const std::size_t bin_index = std::min(bin_count - 1, static_cast<std::size_t>(k_mag / bin_width));

        std::complex<double> delta_k{0.0, 0.0};
        for (std::size_t ix = 0; ix < mesh_n; ++ix) {
          for (std::size_t iy = 0; iy < mesh_n; ++iy) {
            for (std::size_t iz = 0; iz < mesh_n; ++iz) {
              const double phase =
                  -k_two_pi * (static_cast<double>(kx_int * static_cast<int>(ix) + ky_int * static_cast<int>(iy) +
                                                   kz_int * static_cast<int>(iz)) /
                               static_cast<double>(mesh_n));
              delta_k += delta[flatten(ix, iy, iz, mesh_n)] * std::complex<double>(std::cos(phase), std::sin(phase));
            }
          }
        }

        delta_k /= static_cast<double>(mesh_count);
        const double power = box_size_mpc_comov * box_size_mpc_comov * box_size_mpc_comov * std::norm(delta_k);
        power_sum[bin_index] += power;
        k_sum[bin_index] += k_mag;
        ++mode_count[bin_index];
      }
    }
  }

  std::vector<PowerSpectrumBin> bins;
  bins.reserve(bin_count);
  for (std::size_t i = 0; i < bin_count; ++i) {
    if (mode_count[i] == 0) {
      continue;
    }
    bins.push_back(PowerSpectrumBin{
        .k_center_code = k_sum[i] / static_cast<double>(mode_count[i]),
        .power_code_volume = power_sum[i] / static_cast<double>(mode_count[i]),
        .mode_count = mode_count[i],
    });
  }
  return bins;
}

std::vector<StarFormationHistoryBin> DiagnosticsEngine::computeStarFormationHistory(
    const core::SimulationState& state,
    std::size_t bin_count) const {
  if (bin_count == 0) {
    throw std::invalid_argument("sf history requires bin_count > 0");
  }

  std::vector<double> mass_sum(bin_count, 0.0);
  for (std::size_t i = 0; i < state.star_particles.size(); ++i) {
    const double a = std::clamp(state.star_particles.formation_scale_factor[i], 0.0, 1.0);
    const std::size_t bin = std::min(bin_count - 1, static_cast<std::size_t>(a * static_cast<double>(bin_count)));
    mass_sum[bin] += state.star_particles.birth_mass_code[i];
  }

  std::vector<StarFormationHistoryBin> result;
  result.reserve(bin_count);
  for (std::size_t i = 0; i < bin_count; ++i) {
    result.push_back(StarFormationHistoryBin{
        .scale_factor_center = (static_cast<double>(i) + 0.5) / static_cast<double>(bin_count),
        .formed_mass_code = mass_sum[i],
    });
  }
  return result;
}

AngularMomentumBudget DiagnosticsEngine::computeAngularMomentumBudget(const core::SimulationState& state) const {
  AngularMomentumBudget budget;
  for (std::size_t i = 0; i < state.particles.size(); ++i) {
    const double mass = state.particles.mass_code[i];
    const std::array<double, 3> r = {
        state.particles.position_x_comoving[i],
        state.particles.position_y_comoving[i],
        state.particles.position_z_comoving[i]};
    const std::array<double, 3> v = {
        state.particles.velocity_x_peculiar[i],
        state.particles.velocity_y_peculiar[i],
        state.particles.velocity_z_peculiar[i]};

    const std::array<double, 3> l = {
        mass * (r[1] * v[2] - r[2] * v[1]),
        mass * (r[2] * v[0] - r[0] * v[2]),
        mass * (r[0] * v[1] - r[1] * v[0])};

    for (std::size_t j = 0; j < 3; ++j) {
      budget.total_l_code[j] += l[j];
    }

    const auto species = static_cast<core::ParticleSpecies>(state.particle_sidecar.species_tag[i]);
    std::array<double, 3>* target = nullptr;
    switch (species) {
      case core::ParticleSpecies::kDarkMatter:
        target = &budget.dark_matter_l_code;
        break;
      case core::ParticleSpecies::kGas:
        target = &budget.gas_l_code;
        break;
      case core::ParticleSpecies::kStar:
        target = &budget.star_l_code;
        break;
      case core::ParticleSpecies::kBlackHole:
        target = &budget.black_hole_l_code;
        break;
      case core::ParticleSpecies::kTracer:
        break;
    }
    if (target != nullptr) {
      for (std::size_t j = 0; j < 3; ++j) {
        (*target)[j] += l[j];
      }
    }
  }
  return budget;
}

std::vector<double> DiagnosticsEngine::computeGasXySliceDensity(
    const core::SimulationState& state,
    std::size_t grid_n) const {
  const double box = m_config.cosmology.box_size_mpc_comoving;
  const double dz_half = box / static_cast<double>(grid_n) * 0.5;
  const double z_mid = box * 0.5;
  const double cell_size = box / static_cast<double>(grid_n);

  std::vector<double> slice(grid_n * grid_n, 0.0);
  std::vector<std::uint32_t> count(grid_n * grid_n, 0);
  for (std::size_t i = 0; i < state.cells.size(); ++i) {
    const double z = state.cells.center_z_comoving[i];
    if (std::abs(z - z_mid) > dz_half) {
      continue;
    }

    const std::size_t ix = std::min(grid_n - 1, static_cast<std::size_t>(std::fmod(std::fmod(state.cells.center_x_comoving[i], box) + box, box) / cell_size));
    const std::size_t iy = std::min(grid_n - 1, static_cast<std::size_t>(std::fmod(std::fmod(state.cells.center_y_comoving[i], box) + box, box) / cell_size));
    const std::size_t idx = flatten2(ix, iy, grid_n);
    slice[idx] += state.gas_cells.density_code[i];
    ++count[idx];
  }

  for (std::size_t i = 0; i < slice.size(); ++i) {
    if (count[i] > 0) {
      slice[i] /= static_cast<double>(count[i]);
    }
  }
  return slice;
}

std::vector<double> DiagnosticsEngine::computeGasXyProjectionDensity(
    const core::SimulationState& state,
    std::size_t grid_n) const {
  const double box = m_config.cosmology.box_size_mpc_comoving;
  const double cell_size = box / static_cast<double>(grid_n);

  std::vector<double> projection(grid_n * grid_n, 0.0);
  for (std::size_t i = 0; i < state.cells.size(); ++i) {
    const std::size_t ix = std::min(grid_n - 1, static_cast<std::size_t>(std::fmod(std::fmod(state.cells.center_x_comoving[i], box) + box, box) / cell_size));
    const std::size_t iy = std::min(grid_n - 1, static_cast<std::size_t>(std::fmod(std::fmod(state.cells.center_y_comoving[i], box) + box, box) / cell_size));
    projection[flatten2(ix, iy, grid_n)] += state.gas_cells.density_code[i];
  }

  return projection;
}

DiagnosticsBundle DiagnosticsEngine::generateBundle(
    const core::SimulationState& state,
    std::uint64_t step_index,
    double scale_factor,
    DiagnosticClass diagnostic_class) const {
  DiagnosticsBundle bundle;
  bundle.step_index = step_index;
  bundle.scale_factor = scale_factor;
  bundle.diagnostic_class = diagnostic_class;
  bundle.diagnostics_execution_policy = m_config.analysis.diagnostics_execution_policy;
  bundle.health = computeRunHealth(state);
  bundle.records.push_back(DiagnosticRecord{
      .name = "run_health_counters",
      .tier = DiagnosticTier::kInfrastructureHealth,
      .maturity = DiagnosticMaturity::kProduction,
      .scalability = DiagnosticScalability::kCheap,
      .executed = true,
      .policy_note = "always_on_for_integrity_checks",
  });
  bundle.records.push_back(DiagnosticRecord{
      .name = "gravity_health_summary",
      .tier = DiagnosticTier::kInfrastructureHealth,
      .maturity = DiagnosticMaturity::kProduction,
      .scalability = DiagnosticScalability::kCheap,
      .executed = true,
      .policy_note = "always_on_gravity_state_sanity",
  });

  if (diagnostic_class == DiagnosticClass::kScienceLight || diagnostic_class == DiagnosticClass::kScienceHeavy) {
    bundle.star_formation_history =
        computeStarFormationHistory(state, static_cast<std::size_t>(m_config.analysis.sf_history_bin_count));
    bundle.records.push_back(DiagnosticRecord{
        .name = "star_formation_history",
        .tier = DiagnosticTier::kValidatedScience,
        .maturity = DiagnosticMaturity::kValidated,
        .scalability = DiagnosticScalability::kModerate,
        .executed = true,
        .policy_note = "validated_lightweight_science",
    });
    bundle.angular_momentum = computeAngularMomentumBudget(state);
    bundle.records.push_back(DiagnosticRecord{
        .name = "angular_momentum_budget",
        .tier = DiagnosticTier::kValidatedScience,
        .maturity = DiagnosticMaturity::kValidated,
        .scalability = DiagnosticScalability::kModerate,
        .executed = true,
        .policy_note = "validated_lightweight_science",
    });
    bundle.quicklook_grid_n = static_cast<std::size_t>(m_config.analysis.quicklook_grid_n);
    bundle.xy_slice_density_code = computeGasXySliceDensity(state, bundle.quicklook_grid_n);
    bundle.xy_projection_density_code = computeGasXyProjectionDensity(state, bundle.quicklook_grid_n);
    bundle.records.push_back(DiagnosticRecord{
        .name = "gas_xy_slice_density",
        .tier = DiagnosticTier::kValidatedScience,
        .maturity = DiagnosticMaturity::kValidated,
        .scalability = DiagnosticScalability::kModerate,
        .executed = true,
        .policy_note = "validated_lightweight_science",
    });
    bundle.records.push_back(DiagnosticRecord{
        .name = "gas_xy_projection_density",
        .tier = DiagnosticTier::kValidatedScience,
        .maturity = DiagnosticMaturity::kValidated,
        .scalability = DiagnosticScalability::kModerate,
        .executed = true,
        .policy_note = "validated_lightweight_science",
    });
  }

  if (diagnostic_class == DiagnosticClass::kScienceHeavy) {
    const bool heavy_allowed =
        (m_config.analysis.diagnostics_execution_policy ==
         core::AnalysisConfig::DiagnosticsExecutionPolicy::kAllIncludingProvisional);
    if (heavy_allowed) {
      bundle.power_spectrum = computePowerSpectrum(
          state,
          static_cast<std::size_t>(m_config.analysis.power_spectrum_mesh_n),
          static_cast<std::size_t>(m_config.analysis.power_spectrum_bin_count));
    }
    bundle.records.push_back(DiagnosticRecord{
        .name = "power_spectrum",
        .tier = DiagnosticTier::kReferenceScience,
        .maturity = DiagnosticMaturity::kProvisional,
        .scalability = DiagnosticScalability::kHeavyReference,
        .executed = heavy_allowed,
        .policy_note = heavy_allowed ? "reference_only_non_default" : "blocked_by_execution_policy",
    });
  }

  return bundle;
}

void DiagnosticsEngine::writeBundle(const DiagnosticsBundle& bundle) const {
  const std::filesystem::path output_path = bundlePath(bundle);
  std::filesystem::create_directories(output_path.parent_path());

  std::ofstream out(output_path);
  if (!out) {
    throw std::runtime_error("failed to open diagnostics bundle path: " + output_path.string());
  }
  out << std::setprecision(17);
  out << "{\n";
  out << "  \"schema_version\": 1,\n";
  out << "  \"run_name\": \"" << m_config.output.run_name << "\",\n";
  out << "  \"diagnostic_class\": \"" << diagnosticClassLabel(bundle.diagnostic_class) << "\",\n";
  out << "  \"diagnostics_execution_policy\": \"" << executionPolicyLabel(bundle.diagnostics_execution_policy)
      << "\",\n";
  out << "  \"step_index\": " << bundle.step_index << ",\n";
  out << "  \"scale_factor\": " << bundle.scale_factor << ",\n";
  out << "  \"units\": {\"frame\": \"comoving\", \"mass\": \"code\", \"length\": \"code\"},\n";
  out << "  \"health\": {\"particle_count\": " << bundle.health.particle_count
      << ", \"cell_count\": " << bundle.health.cell_count << ", \"star_count\": " << bundle.health.star_count
      << ", \"ownership_invariants_ok\": " << (bundle.health.ownership_invariants_ok ? "true" : "false")
      << ", \"unique_particle_ids_ok\": " << (bundle.health.unique_particle_ids_ok ? "true" : "false")
      << ", \"gravity_softening_sidecar_size_ok\": "
      << (bundle.health.gravity_softening_sidecar_size_ok ? "true" : "false")
      << ", \"non_finite_particles\": " << bundle.health.non_finite_particles
      << ", \"non_finite_cells\": " << bundle.health.non_finite_cells
      << ", \"non_finite_gravity_softening\": " << bundle.health.non_finite_gravity_softening
      << ", \"non_positive_particle_mass\": " << bundle.health.non_positive_particle_mass << "},\n";
  out << "  \"diagnostic_records\": [";
  for (std::size_t i = 0; i < bundle.records.size(); ++i) {
    const auto& record = bundle.records[i];
    out << (i == 0 ? "" : ", ") << "{\"name\": \"" << record.name << "\", \"tier\": \""
        << diagnosticTierLabel(record.tier) << "\", \"maturity\": \"" << diagnosticMaturityLabel(record.maturity)
        << "\", \"scalability\": \"" << diagnosticScalabilityLabel(record.scalability)
        << "\", \"executed\": " << (record.executed ? "true" : "false") << ", \"policy_note\": \""
        << record.policy_note << "\"}";
  }
  out << "],\n";

  out << "  \"angular_momentum_budget\": {\"total_norm\": " << magnitude(bundle.angular_momentum.total_l_code)
      << ", \"gas_norm\": " << magnitude(bundle.angular_momentum.gas_l_code)
      << ", \"star_norm\": " << magnitude(bundle.angular_momentum.star_l_code)
      << ", \"dark_matter_norm\": " << magnitude(bundle.angular_momentum.dark_matter_l_code)
      << ", \"black_hole_norm\": " << magnitude(bundle.angular_momentum.black_hole_l_code) << "},\n";

  out << "  \"power_spectrum\": [";
  for (std::size_t i = 0; i < bundle.power_spectrum.size(); ++i) {
    const auto& b = bundle.power_spectrum[i];
    out << (i == 0 ? "" : ", ") << "{\"k_code\": " << b.k_center_code << ", \"p_code_volume\": "
        << b.power_code_volume << ", \"mode_count\": " << b.mode_count << "}";
  }
  out << "],\n";

  out << "  \"sf_history\": [";
  for (std::size_t i = 0; i < bundle.star_formation_history.size(); ++i) {
    const auto& b = bundle.star_formation_history[i];
    out << (i == 0 ? "" : ", ")
        << "{\"scale_factor\": " << b.scale_factor_center << ", \"formed_mass_code\": " << b.formed_mass_code
        << "}";
  }
  out << "]\n";
  out << "}\n";

  if (!bundle.xy_projection_density_code.empty()) {
    const std::filesystem::path quicklook = quicklookPath(bundle);
    std::ofstream csv(quicklook);
    if (!csv) {
      throw std::runtime_error("failed to open quicklook path: " + quicklook.string());
    }
    const std::size_t n = bundle.quicklook_grid_n;
    for (std::size_t ix = 0; ix < n; ++ix) {
      for (std::size_t iy = 0; iy < n; ++iy) {
        if (iy > 0) {
          csv << ',';
        }
        csv << bundle.xy_projection_density_code[flatten2(ix, iy, n)];
      }
      csv << '\n';
    }
  }
}

void DiagnosticsEngine::enforceRetentionPolicy() const {
  const std::filesystem::path output_dir = diagnosticsOutputDirectory();
  if (!std::filesystem::exists(output_dir)) {
    return;
  }

  std::vector<std::filesystem::path> json_files;
  for (const auto& entry : std::filesystem::directory_iterator(output_dir)) {
    if (entry.is_regular_file() && entry.path().extension() == ".json") {
      json_files.push_back(entry.path());
    }
  }

  std::sort(json_files.begin(), json_files.end());
  const std::size_t keep = static_cast<std::size_t>(m_config.analysis.retention_bundle_count);
  if (json_files.size() <= keep) {
    return;
  }

  const std::size_t remove_count = json_files.size() - keep;
  for (std::size_t i = 0; i < remove_count; ++i) {
    const std::filesystem::path to_remove = json_files[i];
    std::filesystem::remove(to_remove);

    std::filesystem::path quicklook = to_remove;
    quicklook.replace_extension(".csv");
    std::filesystem::remove(quicklook);
  }
}

std::filesystem::path DiagnosticsEngine::diagnosticsOutputDirectory() const {
  return std::filesystem::path(m_config.output.output_directory) / m_config.output.run_name / "diagnostics";
}

std::filesystem::path DiagnosticsEngine::bundlePath(const DiagnosticsBundle& bundle) const {
  std::ostringstream filename;
  filename << m_config.analysis.diagnostics_stem << "_" << diagnosticClassLabel(bundle.diagnostic_class) << "_step_"
           << std::setw(8) << std::setfill('0') << bundle.step_index << ".json";
  return diagnosticsOutputDirectory() / filename.str();
}

std::filesystem::path DiagnosticsEngine::quicklookPath(const DiagnosticsBundle& bundle) const {
  std::filesystem::path json = bundlePath(bundle);
  json.replace_extension(".csv");
  return json;
}

DiagnosticsCallback::DiagnosticsCallback(core::SimulationConfig config)
    : m_engine(config), m_config(std::move(config)) {}

std::string_view DiagnosticsCallback::callbackName() const { return "analysis_diagnostics"; }

void DiagnosticsCallback::onStage(core::StepContext& context) {
  if (context.stage != core::IntegrationStage::kAnalysisHooks || !m_config.analysis.enable_diagnostics) {
    return;
  }

  const std::uint64_t step = context.integrator_state.step_index;

  if (step % static_cast<std::uint64_t>(m_config.analysis.run_health_interval_steps) == 0) {
    runDiagnostics(context, DiagnosticClass::kRunHealth);
  }
  if (step % static_cast<std::uint64_t>(m_config.analysis.science_light_interval_steps) == 0) {
    if (m_config.analysis.diagnostics_execution_policy !=
        core::AnalysisConfig::DiagnosticsExecutionPolicy::kRunHealthOnly) {
      runDiagnostics(context, DiagnosticClass::kScienceLight);
    }
  }
  if (step % static_cast<std::uint64_t>(m_config.analysis.science_heavy_interval_steps) == 0) {
    if (m_config.analysis.diagnostics_execution_policy ==
        core::AnalysisConfig::DiagnosticsExecutionPolicy::kAllIncludingProvisional) {
      runDiagnostics(context, DiagnosticClass::kScienceHeavy);
    }
  }

  m_engine.enforceRetentionPolicy();
}

const DiagnosticsTiming& DiagnosticsCallback::timing() const noexcept { return m_timing; }

void DiagnosticsCallback::runDiagnostics(core::StepContext& context, DiagnosticClass diagnostic_class) {
  const auto begin = std::chrono::steady_clock::now();
  const DiagnosticsBundle bundle = m_engine.generateBundle(
      context.state,
      context.integrator_state.step_index,
      context.integrator_state.current_scale_factor,
      diagnostic_class);
  m_engine.writeBundle(bundle);
  const auto end = std::chrono::steady_clock::now();
  const double elapsed_ms =
      static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()) / 1000.0;

  if (diagnostic_class == DiagnosticClass::kRunHealth) {
    m_timing.cumulative_run_health_ms += elapsed_ms;
    ++m_timing.run_health_calls;
  } else if (diagnostic_class == DiagnosticClass::kScienceLight) {
    m_timing.cumulative_light_ms += elapsed_ms;
    ++m_timing.light_calls;
  } else {
    m_timing.cumulative_heavy_ms += elapsed_ms;
    ++m_timing.heavy_calls;
  }
}

}  // namespace cosmosim::analysis
