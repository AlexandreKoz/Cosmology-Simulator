#include "cosmosim/workflows/time_coordinator.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <numeric>
#include <optional>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

#include "cosmosim/amr/amr_hydro_orchestrator.hpp"
#include "cosmosim/core/constants.hpp"
#include "cosmosim/core/profiling.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/core/units.hpp"
#include "cosmosim/io/restart_checkpoint.hpp"
#include "cosmosim/workflows/gravity_runtime.hpp"
#include "cosmosim/workflows/hydro_amr_runtime.hpp"
#include "cosmosim/workflows/runtime_module_registry.hpp"
#include "cosmosim/workflows/runtime_resources.hpp"
#include "cosmosim/workflows/runtime_services.hpp"
#include "workflows/internal/cartesian_gas_cell_layout.hpp"
#include "cosmosim/workflows/migration_balance_runtime.hpp"
#include "cosmosim/workflows/output_restart_runtime.hpp"

namespace cosmosim::workflows {
namespace {

constexpr RuntimeEpochField k_particle_stage_epochs =
    RuntimeEpochField::kParticleIndex |
    RuntimeEpochField::kSchedulerTick |
    RuntimeEpochField::kStepIndex;
constexpr RuntimeEpochField k_state_stage_epochs =
    RuntimeEpochField::kParticleIndex |
    RuntimeEpochField::kCellIndex |
    RuntimeEpochField::kGasIdentity |
    RuntimeEpochField::kSchedulerTick |
    RuntimeEpochField::kStepIndex;

[[nodiscard]] std::string formatRuntimeDouble(double value) {
  std::ostringstream stream;
  stream << std::scientific
         << std::setprecision(std::numeric_limits<double>::max_digits10)
         << value;
  return stream.str();
}

void ensureSchedulerCoversState(
    std::size_t required_size,
    core::HierarchicalTimeBinScheduler& scheduler,
    std::string_view scheduler_name) {
  if (required_size > static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max())) {
    throw std::overflow_error(std::string(scheduler_name) + " element count exceeds uint32 range");
  }
  const std::uint32_t required_count = static_cast<std::uint32_t>(required_size);
  if (required_count <= scheduler.elementCount()) {
    return;
  }
  const std::uint8_t new_element_bin = scheduler.maxBin() > 0 ? 1U : 0U;
  const std::uint64_t bin_period = 1ULL << new_element_bin;
  if (scheduler.currentTick() > std::numeric_limits<std::uint64_t>::max() - bin_period) {
    throw std::overflow_error(std::string(scheduler_name) + " next activation tick overflows uint64");
  }
  const std::uint64_t first_activation_tick =
      ((scheduler.currentTick() / bin_period) + 1ULL) * bin_period;
  scheduler.appendElements(
      required_count - scheduler.elementCount(), new_element_bin, first_activation_tick);
}

void initializeSchedulerBins(
    const core::SimulationState& state,
    core::HierarchicalTimeBinScheduler& particle_scheduler,
    core::HierarchicalTimeBinScheduler& gas_cell_scheduler) {
  if (state.particles.size() > static_cast<std::size_t>(
          std::numeric_limits<std::uint32_t>::max()) ||
      state.cells.size() > static_cast<std::size_t>(
          std::numeric_limits<std::uint32_t>::max())) {
    throw std::overflow_error(
        "reference workflow scheduler element count exceeds uint32 range");
  }
  const std::uint32_t particle_count =
      static_cast<std::uint32_t>(state.particles.size());
  const std::uint32_t cell_count =
      static_cast<std::uint32_t>(state.cells.size());
  const std::uint8_t particle_default_bin =
      particle_scheduler.maxBin() > 0 ? 1U : 0U;
  const std::uint8_t gas_default_bin =
      gas_cell_scheduler.maxBin() > 0 ? 1U : 0U;
  particle_scheduler.reset(particle_count, particle_default_bin, 0U);
  gas_cell_scheduler.reset(cell_count, gas_default_bin, 0U);
  if (cell_count > 0U) {
    state.requireGasCellIdentityMapCoversDenseRows(
        "initialize independent gas-cell scheduler");
    for (std::uint32_t cell_index = 0; cell_index < cell_count; ++cell_index) {
      gas_cell_scheduler.setElementBin(
          cell_index, 0U, gas_cell_scheduler.currentTick());
    }
  }
}

void ensureSchedulersCoverState(
    const core::SimulationState& state,
    core::HierarchicalTimeBinScheduler& particle_scheduler,
    core::HierarchicalTimeBinScheduler& gas_cell_scheduler) {
  ensureSchedulerCoversState(state.particles.size(), particle_scheduler, "particle scheduler");
  ensureSchedulerCoversState(state.cells.size(), gas_cell_scheduler, "gas-cell scheduler");
}

void syncTimeBinsFromSchedulers(
    const core::HierarchicalTimeBinScheduler& particle_scheduler,
    const core::HierarchicalTimeBinScheduler& gas_cell_scheduler,
    core::SimulationState& state) {
  core::syncTimeBinMirrorsFromScheduler(
      particle_scheduler, state, core::TimeBinMirrorDomain::kParticles);
  core::syncGasCellTimeBinMirrorsFromGasCellScheduler(gas_cell_scheduler, state);
}

[[nodiscard]] double newtonGCodeFromUnits(const core::UnitSystem& units) {
  return core::newtonGravitationalConstantCode(units);
}

[[nodiscard]] gravity::TreeSofteningSpeciesPolicy speciesSofteningByTag(
    const core::SimulationConfig& config) {
  gravity::TreeSofteningSpeciesPolicy values{};
  values.enabled = config.numerics.gravity_softening_gas_kpc_comoving > 0.0 ||
      config.numerics.gravity_softening_dark_matter_kpc_comoving > 0.0 ||
      config.numerics.gravity_softening_star_kpc_comoving > 0.0 ||
      config.numerics.gravity_softening_black_hole_kpc_comoving > 0.0 ||
      config.numerics.gravity_softening_tracer_kpc_comoving > 0.0;
  values.epsilon_comoving_by_species.fill(config.numerics.gravity_softening_kpc_comoving * 1.0e-3);
  values.epsilon_comoving_by_species[static_cast<std::size_t>(core::ParticleSpecies::kDarkMatter)] =
      (config.numerics.gravity_softening_dark_matter_kpc_comoving > 0.0)
      ? (config.numerics.gravity_softening_dark_matter_kpc_comoving * 1.0e-3)
      : values.epsilon_comoving_by_species[static_cast<std::size_t>(core::ParticleSpecies::kDarkMatter)];
  values.epsilon_comoving_by_species[static_cast<std::size_t>(core::ParticleSpecies::kGas)] =
      (config.numerics.gravity_softening_gas_kpc_comoving > 0.0)
      ? (config.numerics.gravity_softening_gas_kpc_comoving * 1.0e-3)
      : values.epsilon_comoving_by_species[static_cast<std::size_t>(core::ParticleSpecies::kGas)];
  values.epsilon_comoving_by_species[static_cast<std::size_t>(core::ParticleSpecies::kStar)] =
      (config.numerics.gravity_softening_star_kpc_comoving > 0.0)
      ? (config.numerics.gravity_softening_star_kpc_comoving * 1.0e-3)
      : values.epsilon_comoving_by_species[static_cast<std::size_t>(core::ParticleSpecies::kStar)];
  values.epsilon_comoving_by_species[static_cast<std::size_t>(core::ParticleSpecies::kBlackHole)] =
      (config.numerics.gravity_softening_black_hole_kpc_comoving > 0.0)
      ? (config.numerics.gravity_softening_black_hole_kpc_comoving * 1.0e-3)
      : values.epsilon_comoving_by_species[static_cast<std::size_t>(core::ParticleSpecies::kBlackHole)];
  values.epsilon_comoving_by_species[static_cast<std::size_t>(core::ParticleSpecies::kTracer)] =
      (config.numerics.gravity_softening_tracer_kpc_comoving > 0.0)
      ? (config.numerics.gravity_softening_tracer_kpc_comoving * 1.0e-3)
      : values.epsilon_comoving_by_species[static_cast<std::size_t>(core::ParticleSpecies::kTracer)];
  return values;
}

struct AdaptiveTimeStepCriteriaStorage {
  std::vector<std::uint32_t> gas_particle_index_by_cell;
  std::vector<std::uint64_t> gas_cell_id_by_cell;
  std::vector<std::uint64_t> patch_id_by_cell;
  std::vector<std::uint32_t> patch_row_by_cell;
  std::vector<double> cell_width_x_code;
  std::vector<double> cell_width_y_code;
  std::vector<double> cell_width_z_code;
};

struct LocalGasCellCflMetadata {
  std::vector<double> cell_width_x_code;
  std::vector<double> cell_width_y_code;
  std::vector<double> cell_width_z_code;
  std::vector<std::uint64_t> patch_id_by_cell;
  std::vector<std::uint32_t> patch_row_by_cell;
};

[[nodiscard]] std::unordered_map<std::uint64_t, std::uint32_t> buildParticleRowById(
    const core::SimulationState& state) {
  std::unordered_map<std::uint64_t, std::uint32_t> row_by_id;
  row_by_id.reserve(state.particles.size());
  for (std::uint32_t row = 0; row < state.particles.size(); ++row) {
    row_by_id.emplace(state.particle_sidecar.particle_id[row], row);
  }
  return row_by_id;
}

[[nodiscard]] std::optional<std::uint32_t> parentParticleRowForGasCellRow(
    const core::SimulationState& state,
    std::uint32_t cell_row,
    const std::unordered_map<std::uint64_t, std::uint32_t>& particle_row_by_id,
    std::string_view caller) {
  const auto* record = state.gas_cell_identity.findByLocalRow(cell_row);
  if (record == nullptr) {
    throw std::runtime_error(std::string(caller) + ": gas-cell identity map is missing a dense local row");
  }
  if (!record->parent_particle_id.has_value()) {
    return std::nullopt;
  }
  const auto parent_it = particle_row_by_id.find(*record->parent_particle_id);
  if (parent_it == particle_row_by_id.end()) {
    return std::nullopt;
  }
  return parent_it->second;
}

[[nodiscard]] const core::GasCellIdentityRecord& gasCellIdentityRecordForLocalRow(
    const core::SimulationState& state,
    std::uint32_t cell_row,
    std::string_view caller) {
  const auto* record = state.gas_cell_identity.findByLocalRow(cell_row);
  if (record == nullptr) {
    throw std::runtime_error(std::string(caller) + ": gas-cell identity map is missing a dense local row");
  }
  return *record;
}

[[nodiscard]] std::uint32_t gasCellOwnerRankForLocalRow(
    const core::SimulationState& state,
    std::uint32_t cell_row,
    const std::unordered_map<std::uint64_t, std::uint32_t>& particle_row_by_id,
    std::string_view caller) {
  if (cell_row >= state.cells.size()) {
    throw std::out_of_range(std::string(caller) + ": gas-cell row is outside CellSoa");
  }
  const std::uint32_t patch_index = state.cells.patch_index[cell_row];
  if (patch_index < state.patches.size()) {
    if (patch_index >= state.patches.owning_rank.size()) {
      throw std::runtime_error(std::string(caller) + ": patch owning-rank lane is shorter than patch table");
    }
    return state.patches.owning_rank[patch_index];
  }

  const auto& identity = gasCellIdentityRecordForLocalRow(state, cell_row, caller);
  if (!identity.parent_particle_id.has_value()) {
    throw std::runtime_error(std::string(caller) + ": parentless gas cell has no valid patch owner");
  }
  const auto parent_it = particle_row_by_id.find(*identity.parent_particle_id);
  if (parent_it == particle_row_by_id.end()) {
    throw std::runtime_error(std::string(caller) + ": gas-cell parent_particle_id is not a local particle");
  }
  return state.particle_sidecar.owning_rank[parent_it->second];
}

// Particle lanes are compatibility mirrors for parent-associated gas only.  The
// gas-cell state remains authoritative; this helper derives a deterministic
// aggregate by stable gas_cell_id so dense-row reorders and split/merge layouts
// cannot select an arbitrary child as the parent mirror.
void synchronizeParentParticleCompatibilityMirrors(
    core::SimulationState& state,
    std::uint32_t world_rank,
    std::string_view caller) {
  state.requireGasCellIdentityMapCoversDenseRows(caller);
  const auto particle_row_by_id = buildParticleRowById(state);
  struct ParentMirrorAccumulator {
    double mass_code = 0.0;
    double momentum_x_code = 0.0;
    double momentum_y_code = 0.0;
    double momentum_z_code = 0.0;
  };
  std::unordered_map<std::uint32_t, std::vector<std::uint32_t>> rows_by_parent;
  for (std::uint32_t cell_row = 0; cell_row < state.cells.size(); ++cell_row) {
    if (gasCellOwnerRankForLocalRow(state, cell_row, particle_row_by_id, caller) != world_rank) {
      continue;
    }
    const auto parent_row = parentParticleRowForGasCellRow(state, cell_row, particle_row_by_id, caller);
    if (parent_row.has_value()) {
      rows_by_parent[*parent_row].push_back(cell_row);
    }
  }
  for (auto& [parent_row, rows] : rows_by_parent) {
    std::sort(rows.begin(), rows.end(), [&](std::uint32_t lhs, std::uint32_t rhs) {
      return gasCellIdentityRecordForLocalRow(state, lhs, caller).gas_cell_id <
          gasCellIdentityRecordForLocalRow(state, rhs, caller).gas_cell_id;
    });
    ParentMirrorAccumulator aggregate;
    for (const std::uint32_t cell_row : rows) {
      const double mass_code = state.cells.mass_code[cell_row];
      aggregate.mass_code += mass_code;
      aggregate.momentum_x_code += mass_code * state.gas_cells.velocity_x_peculiar[cell_row];
      aggregate.momentum_y_code += mass_code * state.gas_cells.velocity_y_peculiar[cell_row];
      aggregate.momentum_z_code += mass_code * state.gas_cells.velocity_z_peculiar[cell_row];
    }
    if (aggregate.mass_code > 0.0) {
      state.particles.mass_code[parent_row] = aggregate.mass_code;
      state.particles.velocity_x_peculiar[parent_row] = aggregate.momentum_x_code / aggregate.mass_code;
      state.particles.velocity_y_peculiar[parent_row] = aggregate.momentum_y_code / aggregate.mass_code;
      state.particles.velocity_z_peculiar[parent_row] = aggregate.momentum_z_code / aggregate.mass_code;
    }
  }
}

[[nodiscard]] internal::CartesianGasCellRowLayout requireCartesianGasCellRowLayout(
    const core::SimulationState& state,
    const core::SimulationConfig& config,
    std::string_view caller) {
  internal::CartesianGasCellLayoutBuildResult result =
      internal::buildCartesianGasCellRowLayout(state, config);
  if (!result.ok()) {
    throw std::runtime_error(
        std::string(caller) + ": fixed Cartesian hydro geometry rejected: " + result.diagnostic);
  }
  return std::move(result.layout);
}

[[nodiscard]] LocalGasCellCflMetadata buildLocalGasCellCflMetadata(
    const core::SimulationState& state,
    const core::SimulationConfig& config) {
  LocalGasCellCflMetadata metadata;
  const std::size_t cell_count = state.cells.size();
  metadata.cell_width_x_code.resize(cell_count);
  metadata.cell_width_y_code.resize(cell_count);
  metadata.cell_width_z_code.resize(cell_count);
  metadata.patch_id_by_cell.assign(cell_count, 0U);
  metadata.patch_row_by_cell.assign(cell_count, std::numeric_limits<std::uint32_t>::max());
  if (cell_count == 0U) {
    return metadata;
  }
  state.requireGasCellIdentityMapCoversDenseRows("hydro CFL metadata construction");

  if (amr::hasProductionAmrHydroCoverage(state)) {
    const std::vector<amr::PatchDescriptor> descriptors = amr::buildProductionAmrPatchDescriptors(state);
    for (const amr::PatchDescriptor& descriptor : descriptors) {
      const amr::AmrHydroPatchGeometry patch_geometry = amr::buildAmrHydroPatchGeometry(state, descriptor);
      for (const amr::AmrHydroCellDescriptor& cell : patch_geometry.real_cells) {
        if (cell.local_cell_row >= cell_count || cell.patch_local_cell == hydro::k_invalid_cell_index) {
          throw std::runtime_error("AMR CFL metadata received an invalid physical cell descriptor");
        }
        const auto* identity = state.gas_cell_identity.findByLocalRow(cell.local_cell_row);
        if (identity == nullptr || identity->gas_cell_id != cell.gas_cell_id ||
            identity->owning_patch_id != descriptor.patch_id) {
          throw std::runtime_error("AMR CFL metadata rejected stale gas-cell identity or patch ownership");
        }
        const std::uint32_t row = cell.local_cell_row;
        metadata.cell_width_x_code[row] = patch_geometry.geometry.cell_width_x_comoving;
        metadata.cell_width_y_code[row] = patch_geometry.geometry.cell_width_y_comoving;
        metadata.cell_width_z_code[row] = patch_geometry.geometry.cell_width_z_comoving;
        metadata.patch_id_by_cell[row] = descriptor.patch_id;
        metadata.patch_row_by_cell[row] = static_cast<std::uint32_t>(cell.patch_local_cell);
      }
    }
    for (std::uint32_t row = 0; row < cell_count; ++row) {
      if (!std::isfinite(metadata.cell_width_x_code[row]) || metadata.cell_width_x_code[row] <= 0.0 ||
          !std::isfinite(metadata.cell_width_y_code[row]) || metadata.cell_width_y_code[row] <= 0.0 ||
          !std::isfinite(metadata.cell_width_z_code[row]) || metadata.cell_width_z_code[row] <= 0.0 ||
          metadata.patch_id_by_cell[row] == 0U ||
          metadata.patch_row_by_cell[row] == std::numeric_limits<std::uint32_t>::max()) {
        throw std::runtime_error("AMR CFL metadata did not cover every authoritative gas-cell row");
      }
    }
    return metadata;
  }

  const internal::CartesianGasCellRowLayout layout = requireCartesianGasCellRowLayout(
      state, config, "hydro CFL metadata construction");
  for (std::uint32_t row = 0; row < cell_count; ++row) {
    metadata.cell_width_x_code[row] = layout.spec.cell_width_x_comoving;
    metadata.cell_width_y_code[row] = layout.spec.cell_width_y_comoving;
    metadata.cell_width_z_code[row] = layout.spec.cell_width_z_comoving;
    const auto* identity = state.gas_cell_identity.findByLocalRow(row);
    if (identity == nullptr || identity->gas_cell_id == 0U) {
      throw std::runtime_error("fixed-patch CFL metadata rejected incomplete gas-cell identity coverage");
    }
    if (state.patches.size() != 0U) {
      if (identity->owning_patch_id == 0U) {
        throw std::runtime_error("fixed-patch CFL metadata rejected a gas cell without explicit patch ownership");
      }
      bool patch_found = false;
      for (std::size_t patch = 0; patch < state.patches.size(); ++patch) {
        if (state.patches.patch_id[patch] == identity->owning_patch_id) {
          if (state.cells.patch_index[row] != patch) {
            throw std::runtime_error("fixed-patch CFL metadata rejected stale dense patch-index mirror");
          }
          patch_found = true;
          break;
        }
      }
      if (!patch_found) {
        throw std::runtime_error("fixed-patch CFL metadata rejected an identity patch absent from PatchSoa");
      }
      metadata.patch_id_by_cell[row] = identity->owning_patch_id;
      metadata.patch_row_by_cell[row] = layout.geometry_row_by_dense_row.at(row);
    }
  }
  return metadata;
}

[[nodiscard]] core::AdaptiveTimeStepCriteriaView buildAdaptiveTimeStepCriteriaView(
    const core::SimulationState& state,
    const core::SimulationConfig& config,
    std::span<const double> particle_accel_x,
    std::span<const double> particle_accel_y,
    std::span<const double> particle_accel_z,
    std::span<const double> cell_accel_x,
    std::span<const double> cell_accel_y,
    std::span<const double> cell_accel_z,
    AdaptiveTimeStepCriteriaStorage& storage) {
  state.requireGasCellIdentityMapCoversDenseRows("adaptive time-bin view construction");
  const auto particle_row_by_id = buildParticleRowById(state);
  constexpr std::uint32_t k_no_parent_particle = std::numeric_limits<std::uint32_t>::max();
  storage.gas_particle_index_by_cell.clear();
  storage.gas_particle_index_by_cell.reserve(state.cells.size());
  storage.gas_cell_id_by_cell.assign(state.gas_cells.gas_cell_id.begin(), state.gas_cells.gas_cell_id.end());
  LocalGasCellCflMetadata cfl_metadata = buildLocalGasCellCflMetadata(state, config);
  for (std::uint32_t cell_index = 0; cell_index < state.cells.size(); ++cell_index) {
    const auto parent_row = parentParticleRowForGasCellRow(
        state, cell_index, particle_row_by_id, "adaptive time-bin view construction");
    storage.gas_particle_index_by_cell.push_back(parent_row.value_or(k_no_parent_particle));
  }
  storage.patch_id_by_cell = std::move(cfl_metadata.patch_id_by_cell);
  storage.patch_row_by_cell = std::move(cfl_metadata.patch_row_by_cell);
  storage.cell_width_x_code = std::move(cfl_metadata.cell_width_x_code);
  storage.cell_width_y_code = std::move(cfl_metadata.cell_width_y_code);
  storage.cell_width_z_code = std::move(cfl_metadata.cell_width_z_code);
  return core::AdaptiveTimeStepCriteriaView{
      .particles = core::TimeStepParticleCriteriaView{
          .velocity_x_peculiar = state.particles.velocity_x_peculiar,
          .velocity_y_peculiar = state.particles.velocity_y_peculiar,
          .velocity_z_peculiar = state.particles.velocity_z_peculiar,
          .species_tag = state.particle_sidecar.species_tag,
          .gravity_softening_comoving = state.particle_sidecar.gravity_softening_comoving,
          .accel_x_comoving = particle_accel_x,
          .accel_y_comoving = particle_accel_y,
          .accel_z_comoving = particle_accel_z,
          .black_hole_particle_index = state.black_holes.particle_index,
          .black_hole_subgrid_mass_code = state.black_holes.subgrid_mass_code,
          .black_hole_accretion_rate_code = state.black_holes.accretion_rate_code,
      },
      .gas_cells = core::TimeStepGasCellCriteriaView{
          .gas_particle_index_by_cell = storage.gas_particle_index_by_cell,
          .gas_cell_id_by_cell = storage.gas_cell_id_by_cell,
          .patch_id_by_cell = storage.patch_id_by_cell,
          .patch_row_by_cell = storage.patch_row_by_cell,
          .cell_width_x_code = storage.cell_width_x_code,
          .cell_width_y_code = storage.cell_width_y_code,
          .cell_width_z_code = storage.cell_width_z_code,
          .cell_mass_code = state.cells.mass_code,
          .velocity_x_peculiar = state.gas_cells.velocity_x_peculiar,
          .velocity_y_peculiar = state.gas_cells.velocity_y_peculiar,
          .velocity_z_peculiar = state.gas_cells.velocity_z_peculiar,
          .density_code = state.gas_cells.density_code,
          .temperature_code = state.gas_cells.temperature_code,
          .sound_speed_code = state.gas_cells.sound_speed_code,
          .accel_x_comoving = cell_accel_x,
          .accel_y_comoving = cell_accel_y,
          .accel_z_comoving = cell_accel_z,
      },
  };
}

enum class SchedulerElementFamily {
  kParticles,
  kGasCells,
};

void updateAdaptiveTimeBinsFromView(
    const core::AdaptiveTimeStepCriteriaView& view,
    core::HierarchicalTimeBinScheduler& scheduler,
    const core::IntegratorState& integrator_state,
    const core::SimulationConfig& config,
    SchedulerElementFamily element_family,
    std::span<const std::uint32_t> requested_elements,
    bool update_all_elements) {
  if (integrator_state.dt_time_code <= 0.0) {
    throw std::invalid_argument("adaptive time-bin update requires dt_time_code > 0");
  }
  const core::TimeStepLimits limits{
      .min_dt_time_code = integrator_state.dt_time_code,
      .max_dt_time_code = integrator_state.dt_time_code * static_cast<double>(1ULL << scheduler.maxBin()),
      .max_bin = scheduler.maxBin(),
  };
  const std::size_t particle_count = view.particles.velocity_x_peculiar.size();
  const std::size_t cell_count = view.gas_cells.cell_mass_code.size();
  if (view.particles.velocity_y_peculiar.size() != particle_count ||
      view.particles.velocity_z_peculiar.size() != particle_count ||
      view.particles.species_tag.size() != particle_count) {
    throw std::invalid_argument("particle timestep criteria view has mismatched extents");
  }
  if (view.gas_cells.density_code.size() != cell_count ||
      view.gas_cells.velocity_x_peculiar.size() != cell_count ||
      view.gas_cells.velocity_y_peculiar.size() != cell_count ||
      view.gas_cells.velocity_z_peculiar.size() != cell_count ||
      view.gas_cells.temperature_code.size() != cell_count ||
      view.gas_cells.sound_speed_code.size() != cell_count ||
      view.gas_cells.gas_particle_index_by_cell.size() != cell_count ||
      view.gas_cells.gas_cell_id_by_cell.size() != cell_count ||
      view.gas_cells.patch_id_by_cell.size() != cell_count ||
      view.gas_cells.patch_row_by_cell.size() != cell_count ||
      view.gas_cells.cell_width_x_code.size() != cell_count ||
      view.gas_cells.cell_width_y_code.size() != cell_count ||
      view.gas_cells.cell_width_z_code.size() != cell_count) {
    throw std::invalid_argument("gas-cell timestep criteria view has mismatched extents");
  }
  const auto species_softening = speciesSofteningByTag(config);
  const double global_softening = config.numerics.gravity_softening_kpc_comoving * 1.0e-3;
  const core::UnitSystem runtime_units = core::makeUnitSystem(
      config.units.length_unit,
      config.units.mass_unit,
      config.units.velocity_unit);
  const double newton_g_code = newtonGCodeFromUnits(runtime_units);
  const double gravity_scale_factor = std::max(integrator_state.current_scale_factor, 1.0e-12);
  std::optional<double> cosmology_dt;
  if (integrator_state.current_scale_factor > 0.0 && config.cosmology.hubble_param > 0.0) {
    core::CosmologyBackgroundConfig background_config;
    background_config.hubble_param = config.cosmology.hubble_param;
    background_config.omega_matter = config.cosmology.omega_matter;
    background_config.omega_lambda = config.cosmology.omega_lambda;
    const core::LambdaCdmBackground background(background_config);
    cosmology_dt = core::computeCosmologyExpansionTimeStep(
        background,
        integrator_state.current_scale_factor,
        config.numerics.cosmology_max_delta_ln_a,
        config.numerics.cosmology_max_hubble_time_fraction,
        integrator_state.time_si_per_code);
  }
  const auto star_formation_dt_for_cell = [&](std::uint32_t cell_index) -> std::optional<double> {
    if (!config.physics.enable_star_formation || cell_index >= cell_count) {
      return std::nullopt;
    }
    const double gas_mass = view.gas_cells.cell_mass_code[cell_index];
    const double rho = view.gas_cells.density_code[cell_index];
    const double temperature = view.gas_cells.temperature_code[cell_index];
    if (gas_mass <= 0.0 || rho < config.physics.sf_density_threshold_code ||
        temperature > config.physics.sf_temperature_threshold_k || config.physics.sf_epsilon_ff <= 0.0) {
      return std::nullopt;
    }
    const double t_ff_code = std::sqrt(3.0 * core::constants::k_pi /
        (32.0 * newton_g_code * std::max(rho, 1.0e-30)));
    return std::max(1.0e-12, config.numerics.source_max_fractional_change * t_ff_code / config.physics.sf_epsilon_ff);
  };
  const auto black_hole_dt_for_particle = [&](std::uint32_t particle_index) -> std::optional<double> {
    if (!config.physics.enable_black_hole_agn || particle_index >= view.particles.species_tag.size() ||
        view.particles.species_tag[particle_index] != static_cast<std::uint32_t>(core::ParticleSpecies::kBlackHole)) {
      return std::nullopt;
    }
    const std::size_t bh_count = view.particles.black_hole_particle_index.size();
    if (view.particles.black_hole_subgrid_mass_code.size() != bh_count ||
        view.particles.black_hole_accretion_rate_code.size() != bh_count) {
      throw std::invalid_argument("black-hole timestep criteria view has mismatched extents");
    }
    for (std::size_t bh_index = 0; bh_index < bh_count; ++bh_index) {
      if (view.particles.black_hole_particle_index[bh_index] != particle_index) {
        continue;
      }
      const double mass = std::max(view.particles.black_hole_subgrid_mass_code[bh_index], 1.0e-30);
      const double mdot = std::max(view.particles.black_hole_accretion_rate_code[bh_index], 0.0);
      if (mdot <= 0.0) {
        return std::nullopt;
      }
      return std::max(1.0e-12, config.numerics.source_max_fractional_change * mass / mdot);
    }
    return std::nullopt;
  };
  const std::size_t expected_element_count =
      element_family == SchedulerElementFamily::kParticles ? particle_count : cell_count;
  if (scheduler.elementCount() != expected_element_count) {
    throw std::invalid_argument(
        "adaptive time-bin scheduler extent does not match its declared element family");
  }
  const auto for_each_requested_element = [&](auto&& callback) {
    if (update_all_elements) {
      for (std::uint32_t element = 0; element < expected_element_count;
           ++element) {
        callback(element);
      }
      return;
    }
    for (const std::uint32_t element : requested_elements) {
      if (element >= expected_element_count) {
        throw std::out_of_range(
            "active timestep-criteria element is outside scheduler extent");
      }
      callback(element);
    }
  };
  if (element_family == SchedulerElementFamily::kGasCells) {
    for_each_requested_element([&](std::uint32_t cell_index) {
      const std::uint32_t gas_index = view.gas_cells.gas_particle_index_by_cell[cell_index];
      if (gas_index != std::numeric_limits<std::uint32_t>::max() && gas_index >= particle_count) {
        throw std::out_of_range("gas-particle index in timestep criteria view is out of range");
      }
      const double vx = view.gas_cells.velocity_x_peculiar[cell_index];
      const double vy = view.gas_cells.velocity_y_peculiar[cell_index];
      const double vz = view.gas_cells.velocity_z_peculiar[cell_index];
      const core::DirectionalCflTimeStepInput hydro_cfl_input{
          .cell_width_axis_code = {
              view.gas_cells.cell_width_x_code[cell_index],
              view.gas_cells.cell_width_y_code[cell_index],
              view.gas_cells.cell_width_z_code[cell_index]},
          .velocity_axis_code = {vx, vy, vz},
          .sound_speed_code = std::max(view.gas_cells.sound_speed_code[cell_index], 0.0),
      };
      const double cfl_dt =
          core::computeDirectionalCflTimeStep(hydro_cfl_input, 0.4);
      const double ax = (cell_index < view.gas_cells.accel_x_comoving.size()) ? view.gas_cells.accel_x_comoving[cell_index] : 0.0;
      const double ay = (cell_index < view.gas_cells.accel_y_comoving.size()) ? view.gas_cells.accel_y_comoving[cell_index] : 0.0;
      const double az = (cell_index < view.gas_cells.accel_z_comoving.size()) ? view.gas_cells.accel_z_comoving[cell_index] : 0.0;
      const double amag = std::sqrt(ax * ax + ay * ay + az * az);
      const double eps = gas_index != std::numeric_limits<std::uint32_t>::max() &&
              !view.particles.gravity_softening_comoving.empty()
          ? view.particles.gravity_softening_comoving[gas_index]
          : species_softening.epsilon_comoving_by_species[static_cast<std::size_t>(core::ParticleSpecies::kGas)];
      const double gravity_dt = core::computeComovingGravityTimeStep(
          {.softening_length_comoving_code = std::max(eps, 1.0e-12),
           .scale_free_acceleration_magnitude_code = amag,
           .scale_factor = gravity_scale_factor},
          0.2);
      scheduler.submitCandidateTimeStep(
          cell_index, cfl_dt, limits, core::TimeStepCandidateSource::kHydroCfl, "gas_cell_hydro_cfl");
      scheduler.submitCandidateTimeStep(
          cell_index, gravity_dt, limits, core::TimeStepCandidateSource::kGravityAcceleration, "gas_cell_gravity_acceleration");
      if (cosmology_dt.has_value()) {
        scheduler.submitCandidateTimeStep(
            cell_index, *cosmology_dt, limits, core::TimeStepCandidateSource::kCosmologyExpansion, "gas_cell_cosmology_expansion");
      }
      if (const auto source_dt = star_formation_dt_for_cell(cell_index); source_dt.has_value()) {
        scheduler.submitCandidateTimeStep(
            cell_index, *source_dt, limits, core::TimeStepCandidateSource::kSourceTerm, "gas_cell_star_formation_source");
      }
    });
    return;
  }

  for_each_requested_element([&](std::uint32_t particle_index) {
    const double ax = (particle_index < view.particles.accel_x_comoving.size()) ? view.particles.accel_x_comoving[particle_index] : 0.0;
    const double ay = (particle_index < view.particles.accel_y_comoving.size()) ? view.particles.accel_y_comoving[particle_index] : 0.0;
    const double az = (particle_index < view.particles.accel_z_comoving.size()) ? view.particles.accel_z_comoving[particle_index] : 0.0;
    const double amag = std::sqrt(ax * ax + ay * ay + az * az);
    const double eps = !view.particles.gravity_softening_comoving.empty()
        ? view.particles.gravity_softening_comoving[particle_index]
        : ((static_cast<std::size_t>(view.particles.species_tag[particle_index]) < species_softening.epsilon_comoving_by_species.size())
              ? species_softening.epsilon_comoving_by_species[static_cast<std::size_t>(view.particles.species_tag[particle_index])]
              : global_softening);
    const double gravity_dt = core::computeComovingGravityTimeStep(
        {.softening_length_comoving_code = std::max(eps, 1.0e-12),
         .scale_free_acceleration_magnitude_code = amag,
         .scale_factor = gravity_scale_factor},
        0.2);
    scheduler.submitCandidateTimeStep(
        particle_index,
        gravity_dt,
        limits,
        core::TimeStepCandidateSource::kGravityAcceleration,
        "particle_gravity_acceleration");
    if (cosmology_dt.has_value()) {
      scheduler.submitCandidateTimeStep(
          particle_index,
          *cosmology_dt,
          limits,
          core::TimeStepCandidateSource::kCosmologyExpansion,
          "particle_cosmology_expansion");
    }
    if (const auto source_dt = black_hole_dt_for_particle(particle_index); source_dt.has_value()) {
      scheduler.submitCandidateTimeStep(
          particle_index,
          *source_dt,
          limits,
          core::TimeStepCandidateSource::kSourceTerm,
          "particle_black_hole_source");
    }
  });
}

void updateAdaptiveTimeBinFamilies(
    core::SimulationState& state,
    core::HierarchicalTimeBinScheduler& particle_scheduler,
    core::HierarchicalTimeBinScheduler& gas_cell_scheduler,
    const core::IntegratorState& integrator_state,
    const core::SimulationConfig& config,
    std::span<const double> particle_accel_x,
    std::span<const double> particle_accel_y,
    std::span<const double> particle_accel_z,
    std::span<const double> cell_accel_x,
    std::span<const double> cell_accel_y,
    std::span<const double> cell_accel_z,
    std::span<const std::uint32_t> active_particle_indices,
    std::span<const std::uint32_t> active_cell_indices,
    bool update_all_elements) {
  AdaptiveTimeStepCriteriaStorage storage;
  const core::AdaptiveTimeStepCriteriaView view = buildAdaptiveTimeStepCriteriaView(
      state,
      config,
      particle_accel_x,
      particle_accel_y,
      particle_accel_z,
      cell_accel_x,
      cell_accel_y,
      cell_accel_z,
      storage);
  updateAdaptiveTimeBinsFromView(
      view,
      particle_scheduler,
      integrator_state,
      config,
      SchedulerElementFamily::kParticles,
      active_particle_indices,
      update_all_elements);
  updateAdaptiveTimeBinsFromView(
      view,
      gas_cell_scheduler,
      integrator_state,
      config,
      SchedulerElementFamily::kGasCells,
      active_cell_indices,
      update_all_elements);
}



}  // namespace

RungZeroTimeState::RungZeroTimeState(std::uint8_t max_bin)
    : m_particle_scheduler(max_bin),
      m_gas_cell_scheduler(max_bin) {}

core::HierarchicalTimeBinScheduler&
RungZeroTimeState::particleScheduler() noexcept {
  return m_particle_scheduler;
}

const core::HierarchicalTimeBinScheduler&
RungZeroTimeState::particleScheduler() const noexcept {
  return m_particle_scheduler;
}

core::HierarchicalTimeBinScheduler&
RungZeroTimeState::gasCellScheduler() noexcept {
  return m_gas_cell_scheduler;
}

const core::HierarchicalTimeBinScheduler&
RungZeroTimeState::gasCellScheduler() const noexcept {
  return m_gas_cell_scheduler;
}

core::IntegratorState& RungZeroTimeState::integratorState() noexcept {
  return m_integrator_state;
}

const core::IntegratorState& RungZeroTimeState::integratorState() const noexcept {
  return m_integrator_state;
}

internal::PendingOutputBoundary& RungZeroTimeState::pendingOutput() noexcept {
  return m_pending_output;
}

const internal::PendingOutputBoundary&
RungZeroTimeState::pendingOutput() const noexcept {
  return m_pending_output;
}

RungZeroTimeState initializeRungZeroTimeState(
    const core::SimulationConfig& config,
    const ReferenceWorkflowOptions& options,
    core::SimulationState& state,
    const core::UnitSystem& units,
    const core::LambdaCdmBackground* cosmology_background,
    const io::RestartReadResult* restart_state) {
  if (config.numerics.hierarchical_max_rung != 0) {
    throw std::logic_error(
        "ReferenceWorkflow production KDK requires hierarchical_max_rung=0");
  }
  const std::uint8_t max_bin = static_cast<std::uint8_t>(
      std::max(0, std::min(config.numerics.hierarchical_max_rung, 12)));
  RungZeroTimeState time_state(max_bin);
  if (restart_state != nullptr) {
    time_state.m_particle_scheduler.importPersistentState(
        restart_state->scheduler_state);
    time_state.m_gas_cell_scheduler.importPersistentState(
        restart_state->gas_cell_scheduler_state);
    time_state.m_integrator_state = restart_state->integrator_state;
    if (restart_state->diagnostics.restart_schema_version <
        io::restartSchema().version) {
      time_state.m_integrator_state.pm_refresh_enabled = true;
    }
    if (time_state.m_particle_scheduler.elementCount() !=
        state.particles.size()) {
      throw std::runtime_error(
          "ReferenceWorkflow restart payload particle scheduler coverage does not match SimulationState");
    }
    if (time_state.m_gas_cell_scheduler.elementCount() != state.cells.size()) {
      throw std::runtime_error(
          "ReferenceWorkflow restart payload gas-cell scheduler coverage does not match SimulationState");
    }
    state.requireGasCellIdentityMapCoversDenseRows(
        "ReferenceWorkflow restart resume");
  } else {
    initializeSchedulerBins(
        state,
        time_state.m_particle_scheduler,
        time_state.m_gas_cell_scheduler);
    if (!options.initial_particle_scheduler_identity_records.empty()) {
      core::rebuildSchedulerFromParticleIdentityRecords(
          time_state.m_particle_scheduler,
          options.initial_particle_scheduler_identity_records,
          state.particle_sidecar.particle_id);
    }
    core::IntegratorState& integrator_state = time_state.m_integrator_state;
    integrator_state.step_index = options.step_index;
    integrator_state.current_time_code = config.numerics.t_code_begin;
    integrator_state.time_si_per_code = units.timeSiPerCode();
    integrator_state.current_scale_factor = cosmology_background != nullptr
        ? config.numerics.a_begin
        : 1.0;
    integrator_state.current_redshift = integrator_state.current_scale_factor > 0.0
        ? 1.0 / integrator_state.current_scale_factor - 1.0
        : 0.0;
    integrator_state.current_hubble_rate_code = cosmology_background != nullptr
        ? cosmology_background->hubbleSi(
              integrator_state.current_scale_factor) *
              integrator_state.time_si_per_code
        : 0.0;
    integrator_state.dt_time_code = options.dt_time_code > 0.0
        ? options.dt_time_code
        : std::max(
              1.0e-6,
              (config.numerics.t_code_end - config.numerics.t_code_begin) /
                  static_cast<double>(
                      std::max(config.numerics.max_global_steps, 1)));
    integrator_state.time_bins.hierarchical_enabled = true;
    integrator_state.time_bins.max_bin =
        time_state.m_particle_scheduler.maxBin();
    integrator_state.pm_refresh_enabled = true;
    integrator_state.pm_sync_state.reset(static_cast<std::uint64_t>(
        std::max(config.numerics.treepm_update_cadence_steps, 1)));
    for (std::size_t particle_index = 0;
         particle_index < state.particles.size(); ++particle_index) {
      state.particle_sidecar.last_drift_time_code[particle_index] =
          integrator_state.current_time_code;
      state.particle_sidecar.last_drift_scale_factor[particle_index] =
          integrator_state.current_scale_factor;
    }
  }
  syncTimeBinsFromSchedulers(
      time_state.m_particle_scheduler,
      time_state.m_gas_cell_scheduler,
      state);
  const io::OutputCadencePersistentState* restored_output =
      restart_state != nullptr ? &restart_state->output_cadence_state : nullptr;
  time_state.m_pending_output = internal::initializeOutputCadence(
      config,
      options,
      time_state.m_integrator_state,
      restored_output);
  return time_state;
}

TimeCoordinator::TimeCoordinator(
    const RuntimeServices& services,
    RungZeroTimeState& time_state,
    GravityRuntime& gravity,
    HydroAmrRuntime& hydro_amr,
    RuntimeExecutionPlan execution_plan,
    const internal::MigrationBalanceRuntime& migration_balance) noexcept
    : m_services(services),
      m_time_state(time_state),
      m_gravity(gravity),
      m_hydro_amr(hydro_amr),
      m_execution_plan(std::move(execution_plan)),
      m_migration_balance(migration_balance) {}

void TimeCoordinator::runRungZeroSegment(
    const core::SimulationConfig& config,
    const ReferenceWorkflowOptions& options,
    core::SimulationState& state,
    const core::LambdaCdmBackground* cosmology_background,
    std::span<const std::uint64_t> expected_global_particle_ids,
    ReferenceWorkflowReport& report,
    core::ProfilerSession& profiler,
    const core::ModePolicy& mode_policy,
    bool restoring_from_restart) {
  core::HierarchicalTimeBinScheduler& particle_scheduler =
      m_time_state.m_particle_scheduler;
  core::HierarchicalTimeBinScheduler& gas_cell_scheduler =
      m_time_state.m_gas_cell_scheduler;
  core::IntegratorState& integrator_state = m_time_state.m_integrator_state;
  internal::PendingOutputBoundary& pending_output = m_time_state.m_pending_output;
  if (!restoring_from_restart) {
    updateAdaptiveTimeBins(
        state,
        particle_scheduler,
        gas_cell_scheduler,
        integrator_state,
        config,
        {},
        {},
        {},
        {},
        {},
        {},
        {},
        {},
        true);
  }
  ensureSchedulersCoverState(state, particle_scheduler, gas_cell_scheduler);

  const std::uint64_t run_start_step_index = integrator_state.step_index;
  const std::uint64_t configured_segment_steps = options.max_steps_override > 0
      ? options.max_steps_override
      : static_cast<std::uint64_t>(std::max(config.numerics.max_global_steps, 0));
  const std::uint64_t target_step_index = integrator_state.step_index + configured_segment_steps;
  core::TransientStepWorkspace workspace;
  while (integrator_state.step_index < target_step_index &&
         integrator_state.current_time_code < config.numerics.t_code_end) {
    const double remaining_time_code =
        config.numerics.t_code_end - integrator_state.current_time_code;
    if (!std::isfinite(remaining_time_code) || remaining_time_code <= 0.0) {
      throw std::runtime_error(
          "ReferenceWorkflow computed an invalid remaining endpoint interval");
    }
    double ordered_step_limit_time_code = config.numerics.t_code_end;
    bool limited_by_output_event = false;
    if (pending_output.snapshot_interval_time_code > 0.0 &&
        pending_output.next_snapshot_time_code < ordered_step_limit_time_code &&
        pending_output.next_snapshot_time_code > integrator_state.current_time_code) {
      ordered_step_limit_time_code = pending_output.next_snapshot_time_code;
      limited_by_output_event = true;
    }
    const double ordered_remaining_time_code =
        ordered_step_limit_time_code - integrator_state.current_time_code;
    if (!std::isfinite(ordered_remaining_time_code) || ordered_remaining_time_code <= 0.0) {
      throw std::runtime_error("ReferenceWorkflow computed an invalid ordered timeline interval");
    }
    if (integrator_state.dt_time_code > ordered_remaining_time_code) {
      const double unclipped_dt_time_code = integrator_state.dt_time_code;
      integrator_state.dt_time_code = ordered_remaining_time_code;
      if (limited_by_output_event) {
        pending_output.restart_resume_dt_time_code = unclipped_dt_time_code;
      }
      profiler.recordEvent(core::RuntimeEvent{
          .event_kind = limited_by_output_event ? "time.output_event_clip" : "time.endpoint_clip",
          .severity = core::RuntimeEventSeverity::kInfo,
          .subsystem = "core.time",
          .step_index = integrator_state.step_index,
          .simulation_time_code = integrator_state.current_time_code,
          .scale_factor = integrator_state.current_scale_factor,
          .message = limited_by_output_event
              ? "integration interval clipped to an ordered output event"
              : "integration interval clipped to the configured endpoint",
          .payload = {{"unclipped_dt_time_code", formatRuntimeDouble(unclipped_dt_time_code)},
                      {"clipped_dt_time_code", formatRuntimeDouble(integrator_state.dt_time_code)},
                      {limited_by_output_event ? "output_event_time_code" : "endpoint_time_code",
                       formatRuntimeDouble(ordered_step_limit_time_code)}},
      });
    }

    const std::span<const std::uint32_t> active_particles =
        particle_scheduler.beginSubstep();
    const std::span<const std::uint32_t> active_cells =
        gas_cell_scheduler.beginSubstep();
    if (particle_scheduler.currentTick() != gas_cell_scheduler.currentTick()) {
      throw std::runtime_error("particle and gas-cell schedulers lost their shared integer timeline");
    }
    const bool local_has_active_work = !active_particles.empty() || !active_cells.empty();
    const std::uint64_t active_work_rank_count =
        m_services.mpi_context.allreduceSumUint64(local_has_active_work ? 1ULL : 0ULL);
    if (active_work_rank_count == 0ULL) {
      particle_scheduler.endSubstep();
      gas_cell_scheduler.endSubstep();
      const bool particle_decomposition_changed = m_migration_balance.rebalance(
          state,
          particle_scheduler,
          gas_cell_scheduler,
          parallel::DecompositionRuntimeMeasurements{},
          active_particles,
          expected_global_particle_ids,
          integrator_state.step_index);
      if (particle_decomposition_changed) {
        m_gravity.commitParticleDecompositionChange();
      }
      ensureSchedulersCoverState(state, particle_scheduler, gas_cell_scheduler);
      syncTimeBinsFromSchedulers(particle_scheduler, gas_cell_scheduler, state);
      continue;
    }

    workspace.clear();
    profiler.counters().addCount("workflow_workspace_reuses", 1U);
    profiler.counters().setCount("scheduler_active_index_copy_bytes", 0U);
    internal::latchOutputRequestForCompletedStep(
        config,
        options,
        integrator_state.step_index + 1U,
        integrator_state.current_time_code + integrator_state.dt_time_code,
        pending_output);
    const double resume_dt_after_step = pending_output.restart_resume_dt_time_code;
    const core::StepBoundaryKind requested_boundary =
        internal::requestedBoundaryForPendingOutput(pending_output);
    const std::uint64_t global_active_particle_count =
        m_services.mpi_context.allreduceSumUint64(
            static_cast<std::uint64_t>(active_particles.size()));
    const std::uint64_t global_particle_count =
        m_services.mpi_context.allreduceSumUint64(
            static_cast<std::uint64_t>(state.particles.size()));
    const std::uint64_t global_active_cell_count =
        m_services.mpi_context.allreduceSumUint64(
            static_cast<std::uint64_t>(active_cells.size()));
    const std::uint64_t global_cell_count =
        m_services.mpi_context.allreduceSumUint64(
            static_cast<std::uint64_t>(state.cells.size()));
    core::ActiveSetDescriptor active_set = core::makeSchedulerActiveSetDescriptor(
        particle_scheduler, state, active_particles, active_cells);
    active_set.has_global_synchronization_metadata = true;
    active_set.globally_complete_active_set =
        global_active_particle_count == global_particle_count &&
        global_active_cell_count == global_cell_count;
    executeSingleStep(
        state,
        integrator_state,
        active_set,
        cosmology_background,
        &workspace,
        &mode_policy,
        &profiler,
        particle_scheduler.currentTick(),
        requested_boundary);
    if (resume_dt_after_step > 0.0) {
      integrator_state.dt_time_code = resume_dt_after_step;
    }

    const std::array runtime_reports{
        core::collectSimulationMemoryReport(state, &workspace),
        m_gravity.memoryReport()};
    profiler.setMemoryReport(core::mergeMemoryReports(runtime_reports));
    state.metadata.step_index = integrator_state.step_index;
    state.metadata.scale_factor = integrator_state.current_scale_factor;
    ensureSchedulersCoverState(state, particle_scheduler, gas_cell_scheduler);
    updateAdaptiveTimeBins(
        state,
        particle_scheduler,
        gas_cell_scheduler,
        integrator_state,
        config,
        m_gravity.particleAccelX(),
        m_gravity.particleAccelY(),
        m_gravity.particleAccelZ(),
        m_gravity.cellAccelX(),
        m_gravity.cellAccelY(),
        m_gravity.cellAccelZ(),
        active_particles,
        active_cells,
        false);
    profiler.counters().addCount(
        "timestep_particle_criteria_evaluations",
        static_cast<std::uint64_t>(active_particles.size()));
    profiler.counters().addCount(
        "timestep_gas_cell_criteria_evaluations",
        static_cast<std::uint64_t>(active_cells.size()));
    ensureSchedulersCoverState(state, particle_scheduler, gas_cell_scheduler);
    particle_scheduler.endSubstep();
    gas_cell_scheduler.endSubstep();

    parallel::DecompositionRuntimeMeasurements rebalance_measurements =
        m_gravity.lastRuntimeDecompositionMeasurements();
    const hydro::HydroProfileEvent& hydro_profile = m_hydro_amr.lastHydroProfile();
    const std::uint64_t hydro_ghost_exchange_bytes =
        m_hydro_amr.ghostExchangeBytesRecent();
    rebalance_measurements.hydro_face_fluxes_recent = hydro_profile.face_count;
    rebalance_measurements.hydro_wall_ms_recent = hydro_profile.total_ms;
    rebalance_measurements.ghost_exchange_bytes_recent += hydro_ghost_exchange_bytes;
    rebalance_measurements.has_measurements = rebalance_measurements.has_measurements ||
        hydro_profile.face_count > 0 || hydro_profile.total_ms > 0.0 ||
        hydro_ghost_exchange_bytes > 0;
    const bool particle_decomposition_changed = m_migration_balance.rebalance(
        state,
        particle_scheduler,
        gas_cell_scheduler,
        rebalance_measurements,
        active_particles,
        expected_global_particle_ids,
        integrator_state.step_index);
    if (particle_decomposition_changed) {
      m_gravity.commitParticleDecompositionChange();
    }
    ensureSchedulersCoverState(state, particle_scheduler, gas_cell_scheduler);
    syncTimeBinsFromSchedulers(particle_scheduler, gas_cell_scheduler, state);
    if (integrator_state.last_completed_restart_safe &&
        integrator_state.last_completed_boundary_kind !=
            core::StepBoundaryKind::kLocalActiveBinStep) {
      executeOutputBoundary(
          state,
          integrator_state,
          &profiler,
          requested_boundary);
    }
  }

  report.completed_steps = integrator_state.step_index - run_start_step_index;
  report.final_time_code = integrator_state.current_time_code;
  report.final_scale_factor = integrator_state.current_scale_factor;
  report.final_hydro_cfl_diagnostics = m_hydro_amr.lastHydroCflDiagnostics();
  report.final_hydro_imported_mpi_ghosts =
      static_cast<std::uint64_t>(m_hydro_amr.remoteImportedGhostCount());
  report.final_hydro_remote_interface_faces =
      static_cast<std::uint64_t>(m_hydro_amr.remoteInterfaceFaceCount());
  report.final_hydro_remote_stale_payloads =
      static_cast<std::uint64_t>(m_hydro_amr.remoteStaleInvalidPayloadCount());
}

void TimeCoordinator::executeSingleStep(
    core::SimulationState& state,
    core::IntegratorState& integrator_state,
    core::ActiveSetDescriptor active_set,
    const core::LambdaCdmBackground* cosmology_background,
    core::TransientStepWorkspace* workspace,
    const core::ModePolicy* mode_policy,
    core::ProfilerSession* profiler_session,
    std::optional<std::uint64_t> expected_scheduler_tick,
    core::StepBoundaryKind requested_boundary_kind) {
  m_lifecycle.executeSingleStepWithDispatcher(
      state,
      integrator_state,
      active_set,
      [this](core::StepContext& context, bool require_output_safe_boundary) {
        dispatchStage(context, require_output_safe_boundary);
      },
      cosmology_background,
      workspace,
      mode_policy,
      profiler_session,
      expected_scheduler_tick,
      requested_boundary_kind);
}

void TimeCoordinator::executeOutputBoundary(
    core::SimulationState& state,
    core::IntegratorState& integrator_state,
    core::ProfilerSession* profiler_session,
    core::StepBoundaryKind requested_boundary_kind) {
  m_lifecycle.executeOutputBoundaryWithDispatcher(
      state,
      integrator_state,
      [this](core::StepContext& context, bool require_output_safe_boundary) {
        dispatchStage(context, require_output_safe_boundary);
      },
      profiler_session,
      requested_boundary_kind);
}

void TimeCoordinator::updateAdaptiveTimeBins(
    core::SimulationState& state,
    core::HierarchicalTimeBinScheduler& particle_scheduler,
    core::HierarchicalTimeBinScheduler& gas_cell_scheduler,
    const core::IntegratorState& integrator_state,
    const core::SimulationConfig& config,
    std::span<const double> particle_accel_x,
    std::span<const double> particle_accel_y,
    std::span<const double> particle_accel_z,
    std::span<const double> cell_accel_x,
    std::span<const double> cell_accel_y,
    std::span<const double> cell_accel_z,
    std::span<const std::uint32_t> active_particle_indices,
    std::span<const std::uint32_t> active_cell_indices,
    bool update_all_elements) {
  updateAdaptiveTimeBinFamilies(
      state,
      particle_scheduler,
      gas_cell_scheduler,
      integrator_state,
      config,
      particle_accel_x,
      particle_accel_y,
      particle_accel_z,
      cell_accel_x,
      cell_accel_y,
      cell_accel_z,
      active_particle_indices,
      active_cell_indices,
      update_all_elements);
}

void TimeCoordinator::dispatchStage(
    core::StepContext& context,
    bool require_output_safe_boundary) {
  (void)require_output_safe_boundary;
  const std::string stage_name =
      "stage." + std::string(core::integrationStageName(context.stage));
  COSMOSIM_PROFILE_SCOPE(context.profiler_session, stage_name);
  if (context.profiler_session != nullptr) {
    context.profiler_session->counters().addCount(
        stage_name + ".invocations", 1U);
  }

  SimulationRuntimeEpochSource epoch_source(
      context.state, m_time_state.m_particle_scheduler, context.integrator_state);
  AnalysisStageView audit_view(
      RuntimeResourceLease(epoch_source, k_state_stage_epochs), context);
  m_execution_plan.executeAuditStage(context.stage, audit_view);

  switch (context.stage) {
    case core::IntegrationStage::kGravityKickPre:
    case core::IntegrationStage::kForceRefresh:
    case core::IntegrationStage::kGravityKickPost: {
      GravityStageView view(
          RuntimeResourceLease(epoch_source, k_state_stage_epochs), context);
      m_execution_plan.executeStage(context.stage, view);
      break;
    }
    case core::IntegrationStage::kDrift: {
      DriftParticleStageView view(
          RuntimeResourceLease(epoch_source, k_particle_stage_epochs), context);
      m_execution_plan.executeStage(context.stage, view);
      break;
    }
    case core::IntegrationStage::kHydroUpdate: {
      HydroAmrStageView view(
          RuntimeResourceLease(epoch_source, k_state_stage_epochs), context);
      m_execution_plan.executeStage(context.stage, view);
      break;
    }
    case core::IntegrationStage::kSourceTerms: {
      SourceMutationStageView view(
          RuntimeResourceLease(epoch_source, k_state_stage_epochs), context);
      m_execution_plan.executeStage(context.stage, view);
      break;
    }
    case core::IntegrationStage::kAnalysisHooks: {
      AnalysisStageView view(
          RuntimeResourceLease(epoch_source, k_state_stage_epochs), context);
      m_execution_plan.executeStage(context.stage, view);
      break;
    }
    case core::IntegrationStage::kOutputCheck: {
      OutputRestartStageView view(
          RuntimeResourceLease(epoch_source, k_state_stage_epochs), context);
      m_execution_plan.executeStage(context.stage, view);
      break;
    }
  }
}

}  // namespace cosmosim::workflows
