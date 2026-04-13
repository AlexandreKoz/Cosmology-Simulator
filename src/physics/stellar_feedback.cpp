#include "cosmosim/physics/stellar_feedback.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

namespace cosmosim::physics {
namespace {

constexpr double k_mass_floor = 1.0e-20;
constexpr double k_energy_floor = 1.0e-30;
constexpr double k_distance_floor = 1.0e-12;
constexpr double k_u01_norm = 1.0 / static_cast<double>(std::numeric_limits<std::uint64_t>::max());

struct WeightedCellDistance {
  std::uint32_t cell_index = 0;
  double distance2 = 0.0;
  double dx = 0.0;
  double dy = 0.0;
  double dz = 0.0;
};

[[nodiscard]] std::uint64_t splitmix64(std::uint64_t x) {
  x += 0x9e3779b97f4a7c15ull;
  x = (x ^ (x >> 30u)) * 0xbf58476d1ce4e5b9ull;
  x = (x ^ (x >> 27u)) * 0x94d049bb133111ebull;
  return x ^ (x >> 31u);
}

[[nodiscard]] double uniform01(std::uint64_t seed, std::uint32_t star_index, std::uint64_t step_seed) {
  std::uint64_t mixed = seed;
  mixed ^= splitmix64(step_seed + 0x9e3779b97f4a7c15ull);
  mixed ^= splitmix64(static_cast<std::uint64_t>(star_index) + 0xbf58476d1ce4e5b9ull);
  return static_cast<double>(splitmix64(mixed)) * k_u01_norm;
}

}  // namespace

void StellarFeedbackModuleState::ensureStarCapacity(std::size_t star_count) {
  if (last_returned_mass_cumulative_code.size() < star_count) {
    last_returned_mass_cumulative_code.resize(star_count, 0.0);
    carry_mass_code.resize(star_count, 0.0);
    carry_metals_code.resize(star_count, 0.0);
    carry_thermal_energy_erg.resize(star_count, 0.0);
    carry_kinetic_energy_erg.resize(star_count, 0.0);
    carry_momentum_code.resize(star_count, 0.0);
  }
}

StellarFeedbackModel::StellarFeedbackModel(StellarFeedbackConfig config) : m_config(std::move(config)) {
  if (m_config.epsilon_thermal < 0.0 || m_config.epsilon_kinetic < 0.0 || m_config.epsilon_momentum < 0.0) {
    throw std::invalid_argument("StellarFeedbackModel: efficiencies must be >= 0");
  }
  if (m_config.sn_energy_erg_per_mass_code <= 0.0 || m_config.momentum_code_per_mass_code < 0.0) {
    throw std::invalid_argument("StellarFeedbackModel: budget scales must be physically positive");
  }
  if (m_config.stochastic_event_probability <= 0.0 || m_config.stochastic_event_probability > 1.0) {
    throw std::invalid_argument("StellarFeedbackModel: stochastic_event_probability must be in (0,1]");
  }
  if (m_config.neighbor_count == 0) {
    throw std::invalid_argument("StellarFeedbackModel: neighbor_count must be > 0");
  }
}

const StellarFeedbackConfig& StellarFeedbackModel::config() const noexcept {
  return m_config;
}

StellarFeedbackBudget StellarFeedbackModel::computeBudget(
    double source_mass_code,
    double returned_mass_code,
    double returned_metals_code) const {
  StellarFeedbackBudget budget;
  budget.source_mass_code = std::max(source_mass_code, 0.0);
  budget.returned_mass_code = std::max(returned_mass_code, 0.0);
  budget.returned_metals_code = std::max(returned_metals_code, 0.0);

  budget.total_energy_erg = budget.source_mass_code * m_config.sn_energy_erg_per_mass_code;
  budget.thermal_energy_erg = budget.total_energy_erg * m_config.epsilon_thermal;
  budget.kinetic_energy_erg = budget.total_energy_erg * m_config.epsilon_kinetic;
  budget.momentum_budget_code = budget.source_mass_code * m_config.momentum_code_per_mass_code * m_config.epsilon_momentum;

  if (m_config.mode == StellarFeedbackMode::kThermal) {
    budget.kinetic_energy_erg = 0.0;
    budget.momentum_budget_code = 0.0;
  } else if (m_config.mode == StellarFeedbackMode::kKinetic) {
    budget.thermal_energy_erg = 0.0;
    budget.momentum_budget_code = 0.0;
  } else if (m_config.mode == StellarFeedbackMode::kMomentum) {
    budget.thermal_energy_erg = 0.0;
    budget.kinetic_energy_erg = 0.0;
  }

  return budget;
}

std::vector<StellarFeedbackTarget> StellarFeedbackModel::selectTargets(
    const core::SimulationState& state,
    std::uint32_t particle_index) const {
  if (particle_index >= state.particles.size() || state.cells.size() == 0) {
    return {};
  }

  const double px = state.particles.position_x_comoving[particle_index];
  const double py = state.particles.position_y_comoving[particle_index];
  const double pz = state.particles.position_z_comoving[particle_index];

  std::vector<WeightedCellDistance> distances;
  distances.reserve(state.cells.size());
  for (std::uint32_t cell_index = 0; cell_index < state.cells.size(); ++cell_index) {
    const double dx = state.cells.center_x_comoving[cell_index] - px;
    const double dy = state.cells.center_y_comoving[cell_index] - py;
    const double dz = state.cells.center_z_comoving[cell_index] - pz;
    const double d2 = dx * dx + dy * dy + dz * dz;
    distances.push_back(WeightedCellDistance{.cell_index = cell_index, .distance2 = d2, .dx = dx, .dy = dy, .dz = dz});
  }

  const std::size_t keep_count = std::min<std::size_t>(m_config.neighbor_count, distances.size());
  std::partial_sort(
      distances.begin(),
      distances.begin() + static_cast<std::ptrdiff_t>(keep_count),
      distances.end(),
      [](const WeightedCellDistance& lhs, const WeightedCellDistance& rhs) { return lhs.distance2 < rhs.distance2; });

  double weight_sum = 0.0;
  for (std::size_t i = 0; i < keep_count; ++i) {
    weight_sum += 1.0 / std::sqrt(std::max(distances[i].distance2, k_distance_floor));
  }
  if (weight_sum <= 0.0) {
    return {};
  }

  std::vector<StellarFeedbackTarget> targets;
  targets.reserve(keep_count);
  for (std::size_t i = 0; i < keep_count; ++i) {
    const double inv_r = 1.0 / std::sqrt(std::max(distances[i].distance2, k_distance_floor));
    targets.push_back(StellarFeedbackTarget{
        .cell_index = distances[i].cell_index,
        .weight = inv_r / weight_sum,
        .radial_dx_comoving = distances[i].dx,
        .radial_dy_comoving = distances[i].dy,
        .radial_dz_comoving = distances[i].dz,
    });
  }
  return targets;
}

StellarFeedbackStepReport StellarFeedbackModel::apply(
    core::SimulationState& state,
    StellarFeedbackModuleState& module_state,
    std::span<const std::uint32_t> active_star_indices,
    std::span<const double> returned_mass_delta_code,
    std::span<const double> returned_metals_delta_code,
    double dt_code) const {
  StellarFeedbackStepReport report;
  if (!m_config.enabled || dt_code <= 0.0 || state.star_particles.size() == 0) {
    state.sidecars.upsert(buildMetadataSidecar(report));
    return report;
  }

  module_state.ensureStarCapacity(state.star_particles.size());
  const std::uint64_t step_seed = state.metadata.step_index;

  for (const std::uint32_t star_index : active_star_indices) {
    ++report.counters.scanned_stars;
    if (star_index >= state.star_particles.size()) {
      continue;
    }

    const std::uint32_t particle_index = state.star_particles.particle_index[star_index];
    if (particle_index >= state.particles.size()) {
      continue;
    }

    const double returned_mass =
        (star_index < returned_mass_delta_code.size()) ? std::max(returned_mass_delta_code[star_index], 0.0) : 0.0;
    const double returned_metals =
        (star_index < returned_metals_delta_code.size()) ? std::max(returned_metals_delta_code[star_index], 0.0) : 0.0;

    const double source_mass = m_config.use_returned_mass_budget
                                   ? returned_mass
                                   : std::max(state.star_particles.birth_mass_code[star_index] * dt_code, 0.0);

    StellarFeedbackStarReport star_report;
    star_report.star_index = star_index;
    star_report.particle_index = particle_index;
    star_report.budget = computeBudget(source_mass, returned_mass, returned_metals);

    star_report.budget.returned_mass_code += module_state.carry_mass_code[star_index];
    star_report.budget.returned_metals_code += module_state.carry_metals_code[star_index];
    star_report.budget.thermal_energy_erg += module_state.carry_thermal_energy_erg[star_index];
    star_report.budget.kinetic_energy_erg += module_state.carry_kinetic_energy_erg[star_index];
    star_report.budget.momentum_budget_code += module_state.carry_momentum_code[star_index];
    module_state.carry_mass_code[star_index] = 0.0;
    module_state.carry_metals_code[star_index] = 0.0;
    module_state.carry_thermal_energy_erg[star_index] = 0.0;
    module_state.carry_kinetic_energy_erg[star_index] = 0.0;
    module_state.carry_momentum_code[star_index] = 0.0;

    if (star_report.budget.returned_mass_code <= k_mass_floor &&
        star_report.budget.thermal_energy_erg <= k_energy_floor &&
        star_report.budget.kinetic_energy_erg <= k_energy_floor &&
        star_report.budget.momentum_budget_code <= k_mass_floor) {
      continue;
    }

    star_report.stochastic_event_fired = true;
    if (m_config.variant == StellarFeedbackVariant::kStochastic) {
      star_report.stochastic_event_fired = stochasticEventFires(star_index, step_seed);
    }

    std::vector<StellarFeedbackTarget> targets = selectTargets(state, particle_index);
    star_report.target_count = targets.size();
    report.counters.target_cells_visited += targets.size();

    if (!star_report.stochastic_event_fired || targets.empty()) {
      star_report.unresolved_mass_code = star_report.budget.returned_mass_code;
      star_report.unresolved_metals_code = star_report.budget.returned_metals_code;
      star_report.unresolved_thermal_energy_erg = star_report.budget.thermal_energy_erg;
      star_report.unresolved_kinetic_energy_erg = star_report.budget.kinetic_energy_erg;
      star_report.unresolved_momentum_code = star_report.budget.momentum_budget_code;
    } else {
      for (const StellarFeedbackTarget& target : targets) {
        const std::uint32_t cell_index = target.cell_index;
        const double weight = target.weight;

        const double mass_add = star_report.budget.returned_mass_code * weight;
        const double metals_add = star_report.budget.returned_metals_code * weight;
        const double thermal_add = star_report.budget.thermal_energy_erg * weight;
        const double kinetic_add = star_report.budget.kinetic_energy_erg * weight;
        const double momentum_add = star_report.budget.momentum_budget_code * weight;

        state.cells.mass_code[cell_index] += mass_add;
        state.gas_cells.density_code[cell_index] += mass_add;

        if (m_config.variant == StellarFeedbackVariant::kDelayedCooling) {
          star_report.delayed_cooling_applied = true;
          star_report.unresolved_thermal_energy_erg += thermal_add;
        } else {
          state.gas_cells.internal_energy_code[cell_index] += thermal_add;
          star_report.deposited_thermal_energy_erg += thermal_add;
        }

        state.gas_cells.internal_energy_code[cell_index] += kinetic_add;
        star_report.deposited_kinetic_energy_erg += kinetic_add;

        // Momentum storage is intentionally sidecar-only for now since gas velocity fields are not in CellSoa.
        star_report.unresolved_momentum_code += momentum_add;
        star_report.deposited_mass_code += mass_add;
        star_report.deposited_metals_code += metals_add;
      }
    }

    module_state.carry_mass_code[star_index] += star_report.unresolved_mass_code;
    module_state.carry_metals_code[star_index] += star_report.unresolved_metals_code;
    module_state.carry_thermal_energy_erg[star_index] += star_report.unresolved_thermal_energy_erg;
    module_state.carry_kinetic_energy_erg[star_index] += star_report.unresolved_kinetic_energy_erg;
    module_state.carry_momentum_code[star_index] += star_report.unresolved_momentum_code;

    ++report.counters.feedback_stars;
    report.counters.source_mass_code += star_report.budget.source_mass_code;
    report.counters.deposited_mass_code += star_report.deposited_mass_code;
    report.counters.deposited_metals_code += star_report.deposited_metals_code;
    report.counters.deposited_thermal_energy_erg += star_report.deposited_thermal_energy_erg;
    report.counters.deposited_kinetic_energy_erg += star_report.deposited_kinetic_energy_erg;
    report.counters.deposited_momentum_code += star_report.deposited_momentum_code;
    report.counters.unresolved_mass_code += star_report.unresolved_mass_code;
    report.counters.unresolved_metals_code += star_report.unresolved_metals_code;
    report.counters.unresolved_thermal_energy_erg += star_report.unresolved_thermal_energy_erg;
    report.counters.unresolved_kinetic_energy_erg += star_report.unresolved_kinetic_energy_erg;
    report.counters.unresolved_momentum_code += star_report.unresolved_momentum_code;
    report.star_reports.push_back(star_report);
  }

  state.sidecars.upsert(buildMetadataSidecar(report));
  return report;
}

core::ModuleSidecarBlock StellarFeedbackModel::buildMetadataSidecar(const StellarFeedbackStepReport& report) const {
  std::ostringstream stream;
  stream << "module=stellar_feedback\n";
  stream << "schema_version=" << m_config.metadata_schema_version << "\n";
  stream << "mode=" << modeToString(m_config.mode) << "\n";
  stream << "variant=" << variantToString(m_config.variant) << "\n";
  stream << "use_returned_mass_budget=" << (m_config.use_returned_mass_budget ? "true" : "false") << "\n";
  stream << "epsilon_thermal=" << m_config.epsilon_thermal << "\n";
  stream << "epsilon_kinetic=" << m_config.epsilon_kinetic << "\n";
  stream << "epsilon_momentum=" << m_config.epsilon_momentum << "\n";
  stream << "sn_energy_erg_per_mass_code=" << m_config.sn_energy_erg_per_mass_code << "\n";
  stream << "momentum_code_per_mass_code=" << m_config.momentum_code_per_mass_code << "\n";
  stream << "neighbor_count=" << m_config.neighbor_count << "\n";
  stream << "delayed_cooling_time_code=" << m_config.delayed_cooling_time_code << "\n";
  stream << "stochastic_event_probability=" << m_config.stochastic_event_probability << "\n";
  stream << "scanned_stars=" << report.counters.scanned_stars << "\n";
  stream << "feedback_stars=" << report.counters.feedback_stars << "\n";
  stream << "target_cells_visited=" << report.counters.target_cells_visited << "\n";
  stream << "source_mass_code=" << report.counters.source_mass_code << "\n";
  stream << "deposited_mass_code=" << report.counters.deposited_mass_code << "\n";
  stream << "deposited_metals_code=" << report.counters.deposited_metals_code << "\n";
  stream << "deposited_thermal_energy_erg=" << report.counters.deposited_thermal_energy_erg << "\n";
  stream << "deposited_kinetic_energy_erg=" << report.counters.deposited_kinetic_energy_erg << "\n";
  stream << "unresolved_mass_code=" << report.counters.unresolved_mass_code << "\n";
  stream << "unresolved_metals_code=" << report.counters.unresolved_metals_code << "\n";
  stream << "unresolved_thermal_energy_erg=" << report.counters.unresolved_thermal_energy_erg << "\n";
  stream << "unresolved_kinetic_energy_erg=" << report.counters.unresolved_kinetic_energy_erg << "\n";
  stream << "unresolved_momentum_code=" << report.counters.unresolved_momentum_code << "\n";

  const std::string text = stream.str();
  core::ModuleSidecarBlock block;
  block.module_name = "stellar_feedback";
  block.schema_version = m_config.metadata_schema_version;
  block.payload.resize(text.size());
  for (std::size_t i = 0; i < text.size(); ++i) {
    block.payload[i] = static_cast<std::byte>(text[i]);
  }
  return block;
}

std::string StellarFeedbackModel::modeToString(StellarFeedbackMode mode) {
  switch (mode) {
    case StellarFeedbackMode::kThermal:
      return "thermal";
    case StellarFeedbackMode::kKinetic:
      return "kinetic";
    case StellarFeedbackMode::kMomentum:
      return "momentum";
    case StellarFeedbackMode::kThermalKineticMomentum:
      return "thermal_kinetic_momentum";
  }
  return "thermal_kinetic_momentum";
}

std::string StellarFeedbackModel::variantToString(StellarFeedbackVariant variant) {
  switch (variant) {
    case StellarFeedbackVariant::kNone:
      return "none";
    case StellarFeedbackVariant::kDelayedCooling:
      return "delayed_cooling";
    case StellarFeedbackVariant::kStochastic:
      return "stochastic";
  }
  return "none";
}

bool StellarFeedbackModel::stochasticEventFires(std::uint32_t star_index, std::uint64_t step_seed) const {
  const double event_u01 = uniform01(m_config.random_seed, star_index, step_seed);
  return event_u01 < m_config.stochastic_event_probability;
}

StellarFeedbackConfig makeStellarFeedbackConfig(const core::PhysicsConfig& physics_config) {
  StellarFeedbackConfig config;
  config.enabled = physics_config.enable_feedback;
  config.use_returned_mass_budget = physics_config.fb_use_returned_mass_budget;
  config.epsilon_thermal = physics_config.fb_epsilon_thermal;
  config.epsilon_kinetic = physics_config.fb_epsilon_kinetic;
  config.epsilon_momentum = physics_config.fb_epsilon_momentum;
  config.sn_energy_erg_per_mass_code = physics_config.fb_sn_energy_erg_per_mass_code;
  config.momentum_code_per_mass_code = physics_config.fb_momentum_code_per_mass_code;
  config.neighbor_count = physics_config.fb_neighbor_count;
  config.delayed_cooling_time_code = physics_config.fb_delayed_cooling_time_code;
  config.stochastic_event_probability = physics_config.fb_stochastic_event_probability;
  config.random_seed = physics_config.fb_random_seed;

  switch (physics_config.fb_mode) {
    case core::FeedbackMode::kThermal:
      config.mode = StellarFeedbackMode::kThermal;
      break;
    case core::FeedbackMode::kKinetic:
      config.mode = StellarFeedbackMode::kKinetic;
      break;
    case core::FeedbackMode::kMomentum:
      config.mode = StellarFeedbackMode::kMomentum;
      break;
    case core::FeedbackMode::kThermalKineticMomentum:
      config.mode = StellarFeedbackMode::kThermalKineticMomentum;
      break;
  }

  switch (physics_config.fb_variant) {
    case core::FeedbackVariant::kDelayedCooling:
      config.variant = StellarFeedbackVariant::kDelayedCooling;
      break;
    case core::FeedbackVariant::kStochastic:
      config.variant = StellarFeedbackVariant::kStochastic;
      break;
    case core::FeedbackVariant::kNone:
      config.variant = StellarFeedbackVariant::kNone;
      break;
  }

  return config;
}

}  // namespace cosmosim::physics
