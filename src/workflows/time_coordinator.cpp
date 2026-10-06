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
#include "cosmosim/core/memory_governor.hpp"
#include "cosmosim/core/retained_capacity_transaction.hpp"
#include "cosmosim/core/profiling.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/core/units.hpp"
#include "cosmosim/io/restart_checkpoint.hpp"
#include "cosmosim/physics/effective_multiphase_ism.hpp"
#include "cosmosim/physics/metal_diffusion.hpp"
#include "cosmosim/workflows/gravity_runtime.hpp"
#include "cosmosim/workflows/hydro_amr_runtime.hpp"
#include "cosmosim/workflows/source_runtime.hpp"
#include "cosmosim/workflows/runtime_module_registry.hpp"
#include "cosmosim/workflows/runtime_resources.hpp"
#include "cosmosim/workflows/runtime_services.hpp"
#include "cosmosim/workflows/runtime_console_reporter.hpp"
#include "workflows/internal/cartesian_gas_cell_layout.hpp"
#include "workflows/internal/gas_cell_ownership.hpp"
#include "workflows/internal/runtime_stage_resource_access.hpp"
#include "workflows/internal/star_formation_geometry.hpp"
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
      0U;
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
  std::vector<double> velocity_divergence_code;
  std::vector<double> metal_diffusion_dt_code;
};

struct LocalGasCellCflMetadata {
  std::vector<double> cell_width_x_code;
  std::vector<double> cell_width_y_code;
  std::vector<double> cell_width_z_code;
  std::vector<std::uint64_t> patch_id_by_cell;
  std::vector<std::uint32_t> patch_row_by_cell;
};

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
    double scale_factor,
    AdaptiveTimeStepCriteriaStorage& storage) {
  state.requireGasCellIdentityMapCoversDenseRows("adaptive time-bin view construction");
  const auto particle_row_by_id = internal::buildParticleRowById(state);
  constexpr std::uint32_t k_no_parent_particle = std::numeric_limits<std::uint32_t>::max();
  storage.gas_particle_index_by_cell.clear();
  storage.gas_particle_index_by_cell.reserve(state.cells.size());
  storage.gas_cell_id_by_cell.assign(state.gas_cells.gas_cell_id.begin(), state.gas_cells.gas_cell_id.end());
  LocalGasCellCflMetadata cfl_metadata = buildLocalGasCellCflMetadata(state, config);
  for (std::uint32_t cell_index = 0; cell_index < state.cells.size(); ++cell_index) {
    const auto parent_row = internal::parentParticleRowForGasCellRow(
        state, cell_index, particle_row_by_id, "adaptive time-bin view construction");
    storage.gas_particle_index_by_cell.push_back(parent_row.value_or(k_no_parent_particle));
  }
  storage.patch_id_by_cell = std::move(cfl_metadata.patch_id_by_cell);
  storage.patch_row_by_cell = std::move(cfl_metadata.patch_row_by_cell);
  storage.cell_width_x_code = std::move(cfl_metadata.cell_width_x_code);
  storage.cell_width_y_code = std::move(cfl_metadata.cell_width_y_code);
  storage.cell_width_z_code = std::move(cfl_metadata.cell_width_z_code);
  storage.velocity_divergence_code.assign(
      state.cells.size(), std::numeric_limits<double>::quiet_NaN());
  storage.metal_diffusion_dt_code.assign(
      state.cells.size(), std::numeric_limits<double>::infinity());
  const double length_to_physical = config.units.coordinate_frame ==
          core::CoordinateFrame::kComoving
      ? std::max(scale_factor, 1.0e-12)
      : 1.0;
  for (std::uint32_t cell_index = 0; cell_index < state.cells.size(); ++cell_index) {
    const internal::StarFormationPatchCellGeometry geometry =
        internal::starFormationPatchCellGeometry(state, cell_index);
    if (!geometry.valid) {
      continue;
    }
    const std::array<std::span<const double>, 3> velocity_fields{
        state.gas_cells.velocity_x_peculiar,
        state.gas_cells.velocity_y_peculiar,
        state.gas_cells.velocity_z_peculiar};
    const std::array<double, 3> spacing_phys_code{
        geometry.dx_stored * length_to_physical,
        geometry.dy_stored * length_to_physical,
        geometry.dz_stored * length_to_physical};
    physics::MetalDiffusionVelocityGradient velocity_gradient;
    bool gradient_valid = true;
    for (std::size_t component = 0; component < velocity_fields.size(); ++component) {
      for (std::size_t axis = 0; axis < spacing_phys_code.size(); ++axis) {
        velocity_gradient.grad[component][axis] =
            internal::starFormationDerivativeAtCell(
                velocity_fields[component], geometry, static_cast<int>(axis),
                spacing_phys_code[axis], cell_index);
        gradient_valid = gradient_valid &&
            std::isfinite(velocity_gradient.grad[component][axis]);
      }
    }
    if (gradient_valid) {
      storage.velocity_divergence_code[cell_index] =
          velocity_gradient.grad[0][0] + velocity_gradient.grad[1][1] +
          velocity_gradient.grad[2][2];
      if (config.physics.enable_metal_diffusion &&
          config.physics.metal_diffusion_model ==
              core::MetalDiffusionModel::kSmagorinsky) {
        const double filter_length = std::cbrt(
            spacing_phys_code[0] * spacing_phys_code[1] * spacing_phys_code[2]);
        const double strain = physics::traceFreeStrainMagnitude(velocity_gradient);
        const double kappa = std::clamp(
            config.physics.metal_diffusion_coefficient * filter_length *
                filter_length * strain,
            config.physics.metal_diffusion_coefficient_floor_code,
            config.physics.metal_diffusion_coefficient_ceiling_code);
        if (kappa > 0.0) {
          const double minimum_spacing_squared = std::min({
              spacing_phys_code[0] * spacing_phys_code[0],
              spacing_phys_code[1] * spacing_phys_code[1],
              spacing_phys_code[2] * spacing_phys_code[2]});
          storage.metal_diffusion_dt_code[cell_index] =
              config.physics.metal_diffusion_cfl * minimum_spacing_squared /
              (6.0 * kappa);
        }
      }
    }
  }
  return core::AdaptiveTimeStepCriteriaView{
      .particles = core::TimeStepParticleCriteriaView{
          .velocity_x_peculiar = state.particles.velocity_x_peculiar,
          .velocity_y_peculiar = state.particles.velocity_y_peculiar,
          .velocity_z_peculiar = state.particles.velocity_z_peculiar,
          .species_tag = state.particle_sidecar.species_tag,
          .homogeneous_dmo_species = state.hasHomogeneousDmoMetadata(),
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
          .velocity_divergence_code = storage.velocity_divergence_code,
          .metal_diffusion_dt_code = storage.metal_diffusion_dt_code,
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

[[nodiscard]] double updateAdaptiveTimeBinsFromView(
    const core::AdaptiveTimeStepCriteriaView& view,
    core::HierarchicalTimeBinScheduler& scheduler,
    const core::IntegratorState& integrator_state,
    const core::SimulationConfig& config,
    const physics::EffectiveMultiphaseEosTable* effective_eos_table,
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
      (!view.particles.homogeneous_dmo_species && view.particles.species_tag.size() != particle_count) ||
      (view.particles.homogeneous_dmo_species && !view.particles.species_tag.empty())) {
    throw std::invalid_argument("particle timestep criteria view has mismatched extents");
  }
  if (view.gas_cells.density_code.size() != cell_count ||
      view.gas_cells.velocity_x_peculiar.size() != cell_count ||
      view.gas_cells.velocity_y_peculiar.size() != cell_count ||
      view.gas_cells.velocity_z_peculiar.size() != cell_count ||
      view.gas_cells.temperature_code.size() != cell_count ||
      view.gas_cells.sound_speed_code.size() != cell_count ||
      view.gas_cells.velocity_divergence_code.size() != cell_count ||
      view.gas_cells.metal_diffusion_dt_code.size() != cell_count ||
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
  const core::ModePolicy mode_policy = core::buildModePolicy(config.mode);
  if (mode_policy.cosmological_comoving_frame &&
      config.units.coordinate_frame == core::CoordinateFrame::kComoving &&
      integrator_state.current_scale_factor > 0.0 &&
      config.cosmology.hubble_param > 0.0) {
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
    if (!config.physics.enable_star_formation || cell_index >= cell_count ||
        (config.physics.star_formation_model !=
             core::StarFormationModelKind::kEffectiveMultiphaseTngLike &&
         config.physics.sf_epsilon_ff <= 0.0)) {
      return std::nullopt;
    }
    const double gas_mass = view.gas_cells.cell_mass_code[cell_index];
    const double stored_density = view.gas_cells.density_code[cell_index];
    const double temperature = view.gas_cells.temperature_code[cell_index];
    if (!(gas_mass > 0.0) || !(stored_density > 0.0) ||
        !std::isfinite(gas_mass) || !std::isfinite(stored_density) || !std::isfinite(temperature)) {
      return std::nullopt;
    }
    if (config.physics.star_formation_model ==
        core::StarFormationModelKind::kLegacySchmidtThreshold) {
      const double scale_factor = std::max(integrator_state.current_scale_factor, 1.0e-12);
      const double physical_density = config.units.coordinate_frame == core::CoordinateFrame::kComoving
          ? stored_density / (scale_factor * scale_factor * scale_factor)
          : stored_density;
      if (physical_density < config.physics.sf_density_threshold_code ||
          temperature > config.physics.sf_temperature_threshold_k) {
        return std::nullopt;
      }
      const double t_ff_code = std::sqrt(3.0 * core::constants::k_pi /
          (32.0 * newton_g_code * std::max(physical_density, 1.0e-30)));
      return std::max(
          1.0e-12,
          config.numerics.source_max_fractional_change * t_ff_code /
              config.physics.sf_epsilon_ff);
    }

    const double scale_factor = std::max(integrator_state.current_scale_factor, 1.0e-12);
    const double physical_density = config.units.coordinate_frame == core::CoordinateFrame::kComoving
        ? stored_density / (scale_factor * scale_factor * scale_factor)
        : stored_density;
    if (config.physics.star_formation_model ==
        core::StarFormationModelKind::kEffectiveMultiphaseTngLike) {
      if (effective_eos_table == nullptr) return std::nullopt;
      const auto equilibrium = effective_eos_table->lookup(physical_density);
      if (!equilibrium.above_threshold || !equilibrium.valid ||
          equilibrium.entry.cold_mass_fraction <= 0.0) {
        return std::nullopt;
      }
      const double long_lived_factor =
          config.physics.sf_effective_birth_mass_convention ==
              core::EffectiveIsmBirthMassConvention::kLongLivedMass
          ? (1.0 - config.physics.sf_effective_massive_star_fraction)
          : 1.0;
      const double rate_per_mass = long_lived_factor *
          equilibrium.entry.cold_mass_fraction /
          std::max(equilibrium.entry.star_formation_timescale_code, 1.0e-30);
      const double maximum_fraction = std::clamp(
          std::min(config.physics.sf_max_fractional_mass_conversion,
                   config.numerics.source_max_fractional_change),
          1.0e-12, 1.0 - 1.0e-12);
      return std::max(1.0e-12, -std::log1p(-maximum_fraction) / rate_per_mass);
    }
    if (config.physics.sf_temperature_safety_ceiling_k > 0.0 &&
        temperature > config.physics.sf_temperature_safety_ceiling_k) {
      return std::nullopt;
    }
    const double t_ff_code = std::sqrt(3.0 * core::constants::k_pi /
        (32.0 * newton_g_code * std::max(physical_density, 1.0e-30)));
    double collapse_time_code = t_ff_code;
    if (config.physics.sf_collapse_timescale ==
        core::StarFormationCollapseTimescale::kMinimumFreeFallOrCompression) {
      const double divergence_code = view.gas_cells.velocity_divergence_code[cell_index];
      if (std::isfinite(divergence_code) && divergence_code < 0.0) {
        collapse_time_code = std::min(t_ff_code, -1.0 / divergence_code);
      }
    }
    const double maximum_fraction = std::clamp(
        std::min(
            config.physics.sf_max_fractional_mass_conversion,
            config.numerics.source_max_fractional_change),
        1.0e-12,
        1.0 - 1.0e-12);
    const double dt_limit = -collapse_time_code * std::log1p(-maximum_fraction) /
        config.physics.sf_epsilon_ff;
    return std::max(1.0e-12, dt_limit);
  };
  const std::size_t bh_count = view.particles.black_hole_particle_index.size();
  if (view.particles.black_hole_subgrid_mass_code.size() != bh_count ||
      view.particles.black_hole_accretion_rate_code.size() != bh_count) {
    throw std::invalid_argument("black-hole timestep criteria view has mismatched extents");
  }
  std::vector<std::size_t> bh_row_by_particle;
  if (config.physics.enable_black_hole_agn && bh_count != 0U) {
    bh_row_by_particle.assign(
        particle_count, std::numeric_limits<std::size_t>::max());
    for (std::size_t bh_index = 0; bh_index < bh_count; ++bh_index) {
      const std::uint32_t particle_index =
          view.particles.black_hole_particle_index[bh_index];
      if (particle_index >= particle_count) {
        throw std::invalid_argument(
            "black-hole sidecar references a particle row outside the particle extent");
      }
      if (bh_row_by_particle[particle_index] !=
          std::numeric_limits<std::size_t>::max()) {
        throw std::invalid_argument(
            "black-hole sidecar contains duplicate particle-row ownership");
      }
      bh_row_by_particle[particle_index] = bh_index;
    }
  }
  const auto black_hole_dt_for_particle = [&](std::uint32_t particle_index) -> std::optional<double> {
    if (!config.physics.enable_black_hole_agn || view.particles.homogeneous_dmo_species ||
        particle_index >= view.particles.species_tag.size() ||
        view.particles.species_tag[particle_index] != static_cast<std::uint32_t>(core::ParticleSpecies::kBlackHole)) {
      return std::nullopt;
    }
    if (bh_row_by_particle.empty()) {
      return std::nullopt;
    }
    const std::size_t bh_index = bh_row_by_particle[particle_index];
    if (bh_index == std::numeric_limits<std::size_t>::max()) {
      return std::nullopt;
    }
    const double mass = std::max(
        view.particles.black_hole_subgrid_mass_code[bh_index], 1.0e-30);
    const double mdot = std::max(
        view.particles.black_hole_accretion_rate_code[bh_index], 0.0);
    if (mdot <= 0.0) {
      return std::nullopt;
    }
    return std::max(
        1.0e-12, config.numerics.source_max_fractional_change * mass / mdot);
  };
  double local_min_dt_time_code = std::numeric_limits<double>::infinity();
  const auto consider_dt = [&](double value) {
    if (std::isinf(value) && value > 0.0) {
      return;  // Explicit no-limit sentinel; another criterion or endpoint must bound the step.
    }
    if (!std::isfinite(value) || value <= 0.0) {
      throw std::runtime_error(
          "adaptive timestep criterion produced an invalid value; expected finite positive or +infinity");
    }
    local_min_dt_time_code = std::min(local_min_dt_time_code, value);
  };
  if (cosmology_dt.has_value()) consider_dt(*cosmology_dt);
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
          .coordinate_frame = config.units.coordinate_frame,
          .scale_factor = gravity_scale_factor,
      };
      const double cfl_dt =
          core::computeDirectionalCflTimeStep(hydro_cfl_input, 0.4);
      consider_dt(cfl_dt);
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
      consider_dt(gravity_dt);
      scheduler.submitCandidateTimeStep(
          cell_index, cfl_dt, limits, core::TimeStepCandidateSource::kHydroCfl, "gas_cell_hydro_cfl");
      scheduler.submitCandidateTimeStep(
          cell_index, gravity_dt, limits, core::TimeStepCandidateSource::kGravityAcceleration, "gas_cell_gravity_acceleration");
      if (cosmology_dt.has_value()) {
        scheduler.submitCandidateTimeStep(
            cell_index, *cosmology_dt, limits, core::TimeStepCandidateSource::kCosmologyExpansion, "gas_cell_cosmology_expansion");
      }
      if (const auto source_dt = star_formation_dt_for_cell(cell_index); source_dt.has_value()) {
        consider_dt(*source_dt);
        scheduler.submitCandidateTimeStep(
            cell_index, *source_dt, limits, core::TimeStepCandidateSource::kSourceTerm, "gas_cell_star_formation_source");
      }
      const double diffusion_dt = view.gas_cells.metal_diffusion_dt_code[cell_index];
      if (config.physics.enable_metal_diffusion && std::isfinite(diffusion_dt) &&
          diffusion_dt > 0.0) {
        consider_dt(diffusion_dt);
        scheduler.submitCandidateTimeStep(
            cell_index, diffusion_dt, limits, core::TimeStepCandidateSource::kSourceTerm,
            "gas_cell_metal_diffusion_parabolic");
      }
    });
    return local_min_dt_time_code;
  }

  for_each_requested_element([&](std::uint32_t particle_index) {
    const double ax = (particle_index < view.particles.accel_x_comoving.size()) ? view.particles.accel_x_comoving[particle_index] : 0.0;
    const double ay = (particle_index < view.particles.accel_y_comoving.size()) ? view.particles.accel_y_comoving[particle_index] : 0.0;
    const double az = (particle_index < view.particles.accel_z_comoving.size()) ? view.particles.accel_z_comoving[particle_index] : 0.0;
    const double amag = std::sqrt(ax * ax + ay * ay + az * az);
    const std::size_t particle_species = view.particles.homogeneous_dmo_species
        ? static_cast<std::size_t>(core::ParticleSpecies::kDarkMatter)
        : static_cast<std::size_t>(view.particles.species_tag[particle_index]);
    const double eps = !view.particles.gravity_softening_comoving.empty()
        ? view.particles.gravity_softening_comoving[particle_index]
        : ((particle_species < species_softening.epsilon_comoving_by_species.size())
              ? species_softening.epsilon_comoving_by_species[particle_species]
              : global_softening);
    const double gravity_dt = core::computeComovingGravityTimeStep(
        {.softening_length_comoving_code = std::max(eps, 1.0e-12),
         .scale_free_acceleration_magnitude_code = amag,
         .scale_factor = gravity_scale_factor},
        0.2);
    consider_dt(gravity_dt);
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
      consider_dt(*source_dt);
      scheduler.submitCandidateTimeStep(
          particle_index,
          *source_dt,
          limits,
          core::TimeStepCandidateSource::kSourceTerm,
          "particle_black_hole_source");
    }
  });
  return local_min_dt_time_code;
}

[[nodiscard]] double updateAdaptiveTimeBinFamilies(
    core::SimulationState& state,
    core::HierarchicalTimeBinScheduler& particle_scheduler,
    core::HierarchicalTimeBinScheduler& gas_cell_scheduler,
    const core::IntegratorState& integrator_state,
    const core::SimulationConfig& config,
    const physics::EffectiveMultiphaseEosTable* effective_eos_table,
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
      integrator_state.current_scale_factor,
      storage);
  const double particle_min_dt = updateAdaptiveTimeBinsFromView(
      view,
      particle_scheduler,
      integrator_state,
      config,
      effective_eos_table,
      SchedulerElementFamily::kParticles,
      active_particle_indices,
      update_all_elements);
  const double gas_min_dt = updateAdaptiveTimeBinsFromView(
      view,
      gas_cell_scheduler,
      integrator_state,
      config,
      effective_eos_table,
      SchedulerElementFamily::kGasCells,
      active_cell_indices,
      update_all_elements);
  return std::min(particle_min_dt, gas_min_dt);
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
  if (config.numerics.hierarchical_max_rung > 0 &&
      (state.cells.size() != 0U || state.star_particles.size() != 0U || state.black_holes.size() != 0U || state.tracers.size() != 0U)) {
    throw std::logic_error("hierarchical KDK currently requires a zero-gas DMO population");
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
        : (config.numerics.t_code_end - config.numerics.t_code_begin);
    integrator_state.time_bins.hierarchical_enabled = true;
    integrator_state.time_bins.max_bin =
        time_state.m_particle_scheduler.maxBin();
    integrator_state.pm_refresh_enabled = true;
    integrator_state.pm_sync_state.reset(static_cast<std::uint64_t>(
        std::max(config.numerics.treepm_update_cadence_steps, 1)));
    state.updateAllParticleDriftEpoch(integrator_state.current_time_code,
        integrator_state.current_scale_factor);
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
    SourceRuntime& source,
    RuntimeExecutionPlan execution_plan,
    const internal::MigrationBalanceRuntime& migration_balance,
    std::shared_ptr<const physics::EffectiveMultiphaseEosTable> effective_eos_table) noexcept
    : m_services(services),
      m_time_state(time_state),
      m_gravity(gravity),
      m_hydro_amr(hydro_amr),
      m_source(source),
      m_execution_plan(std::move(execution_plan)),
      m_migration_balance(migration_balance),
      m_effective_eos_table(std::move(effective_eos_table)) {}

void TimeCoordinator::runRungZeroSegment(
    const core::SimulationConfig& config,
    const ReferenceWorkflowOptions& options,
    core::SimulationState& state,
    const core::LambdaCdmBackground* cosmology_background,
    std::vector<std::uint64_t>& expected_global_particle_ids,
    ReferenceWorkflowReport& report,
    core::ProfilerSession& profiler,
    const core::ModePolicy& mode_policy,
    bool restoring_from_restart) {
  if (config.numerics.hierarchical_max_rung > 0) {
    runHierarchicalSegment(config, options, state, cosmology_background,
        expected_global_particle_ids, report, profiler, mode_policy);
    return;
  }
  core::HierarchicalTimeBinScheduler& particle_scheduler =
      m_time_state.m_particle_scheduler;
  core::HierarchicalTimeBinScheduler& gas_cell_scheduler =
      m_time_state.m_gas_cell_scheduler;
  core::IntegratorState& integrator_state = m_time_state.m_integrator_state;
  internal::PendingOutputBoundary& pending_output = m_time_state.m_pending_output;
  if (!restoring_from_restart) {
    static_cast<void>(updateAdaptiveTimeBins(
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
        true));
  }
  ensureSchedulersCoverState(state, particle_scheduler, gas_cell_scheduler);

  const auto install_authoritative_domain_geometry = [&]() {
    std::exception_ptr install_failure;
    try {
      core::MemoryReservation geometry_reservation;
      {
        const auto leaves = m_migration_balance.authoritativeTopDomainLeaves(
            state, m_gravity.decompositionEpoch(), &geometry_reservation);
        m_gravity.installAuthoritativeTopDomainLeaves(
            leaves, state.gravitySourceGeneration());
      }
      if (geometry_reservation.valid()) {
        const std::array geometry_reports{
            core::collectSimulationMemoryReport(state),
            core::collectSchedulerMemoryReport(
                particle_scheduler, gas_cell_scheduler),
            m_gravity.memoryReport(),
            m_hydro_amr.memoryReport(),
            m_source.memoryReport(),
            m_services.profiler.retainedEventMemoryReport()};
        geometry_reservation.reconcileBaselineOwnedAndRelease(
            core::memoryReportBaselineOwnedBytes(
                core::mergeMemoryReports(geometry_reports)));
      }
    } catch (...) {
      install_failure = std::current_exception();
    }
    FailureCoordinator(m_services).rethrowCollectiveFailure(
        install_failure, "authoritative top-domain installation");
  };
  // Domain geometry is a decomposition-owned contract. Install at the segment
  // boundary (including restart) and again only after ownership decomposition
  // changes; each install replaces both the stable seed leaf set and the
  // current published set. commitParticleDecompositionChange() explicitly
  // invalidates geometry freshness first, so stale routing geometry is never
  // consumable between commit and reinstall. Between those events
  // GravityRuntime refits published leaf bounds from the unchanged seed set
  // against the current source generation before each force solve; TreePM
  // validates freshness (generation equality) and the existing conservative
  // coverage check before routing.
  install_authoritative_domain_geometry();

  const std::uint64_t run_start_step_index = integrator_state.step_index;
  const std::uint64_t configured_segment_steps = options.max_steps_override > 0
      ? options.max_steps_override
      : static_cast<std::uint64_t>(std::max(config.numerics.max_global_steps, 0));
  const std::uint64_t target_step_index = integrator_state.step_index + configured_segment_steps;
  core::TransientStepWorkspace workspace(m_services.memory_governor);
  while (integrator_state.step_index < target_step_index &&
         integrator_state.current_time_code < config.numerics.t_code_end) {
    const RuntimeConsoleReporter::Clock::time_point console_step_begin =
        RuntimeConsoleReporter::Clock::now();
    const std::uint64_t pm_refresh_count_before = m_gravity.longRangeRefreshCount();
    const std::uint64_t pm_reuse_count_before = m_gravity.longRangeReuseCount();
    // Rung zero has one physical timestep authority. Re-evaluate all local
    // criteria before constructing KDK stage times, then agree the minimum
    // collectively so every rank advances the identical interval.
    double local_physical_dt = std::numeric_limits<double>::infinity();
    std::exception_ptr local_timestep_failure;
    try {
      local_physical_dt = updateAdaptiveTimeBins(
          state, particle_scheduler, gas_cell_scheduler, integrator_state, config,
          m_gravity.particleAccelX(), m_gravity.particleAccelY(), m_gravity.particleAccelZ(),
          m_gravity.cellAccelX(), m_gravity.cellAccelY(), m_gravity.cellAccelZ(),
          {}, {}, true);
      const bool explicit_no_limit = std::isinf(local_physical_dt) && local_physical_dt > 0.0;
      if ((!std::isfinite(local_physical_dt) && !explicit_no_limit) ||
          local_physical_dt <= 0.0) {
        throw std::runtime_error(
            "local global-timestep candidate is invalid; expected finite positive or +infinity");
      }
    } catch (...) {
      local_timestep_failure = std::current_exception();
    }
    m_services.mpi_context.rethrowCollectivePreparationFailure(
        local_timestep_failure, "global timestep candidate selection");
    double accepted_dt = m_services.mpi_context.allreduceMinDouble(local_physical_dt);
    if (options.dt_time_code > 0.0) accepted_dt = std::min(accepted_dt, options.dt_time_code);
    const double step_begin_scale_factor = integrator_state.current_scale_factor;
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
    if (accepted_dt > ordered_remaining_time_code) {
      const double unclipped_dt_time_code = accepted_dt;
      accepted_dt = ordered_remaining_time_code;
      if (limited_by_output_event && std::isfinite(unclipped_dt_time_code)) {
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
                      {"clipped_dt_time_code", formatRuntimeDouble(accepted_dt)},
                      {limited_by_output_event ? "output_event_time_code" : "endpoint_time_code",
                       formatRuntimeDouble(ordered_step_limit_time_code)}},
      });
    }
    if (!std::isfinite(accepted_dt) || accepted_dt <= 0.0) {
      throw std::runtime_error("ReferenceWorkflow global timestep selection produced an invalid dt");
    }
    integrator_state.dt_time_code = accepted_dt;
    const double accepted_dt_time_code_for_console = accepted_dt;

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
        install_authoritative_domain_geometry();
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
    const bool console_snapshot_requested = pending_output.snapshot_due;
    const bool console_restart_requested = pending_output.checkpoint_due;
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
    m_newly_created_particle_ids.clear();
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
    if (cosmology_background != nullptr) {
      const double a_next = integrator_state.current_scale_factor;
      if (!std::isfinite(a_next) || a_next <= 0.0 || a_next + 1.0e-14 < step_begin_scale_factor) {
        throw std::runtime_error("cosmological step violated finite positive monotonic scale-factor invariant");
      }
      const double delta_ln_a = std::log(a_next / step_begin_scale_factor);
      const double allowed_delta = config.numerics.cosmology_max_delta_ln_a * (1.0 + 1.0e-10);
      if (!std::isfinite(delta_ln_a) || delta_ln_a > allowed_delta) {
        throw std::runtime_error("cosmological step exceeded numerics.cosmology_max_delta_ln_a");
      }
      if ((config.numerics.integrator_time_variable == core::IntegratorTimeVariable::kScaleFactor ||
           config.numerics.integrator_time_variable == core::IntegratorTimeVariable::kLogScaleFactor) &&
          a_next > config.numerics.a_end + 1.0e-12 * std::max(1.0, config.numerics.a_end)) {
        throw std::runtime_error("cosmological step crossed the authoritative scale-factor endpoint");
      }
    }
    expected_global_particle_ids.insert(
        expected_global_particle_ids.end(),
        m_newly_created_particle_ids.begin(),
        m_newly_created_particle_ids.end());
    if (resume_dt_after_step > 0.0) {
      integrator_state.dt_time_code = resume_dt_after_step;
    }

    const std::array runtime_reports{
        core::collectSimulationMemoryReport(state, &workspace),
        core::collectSchedulerMemoryReport(
            particle_scheduler, gas_cell_scheduler),
        m_gravity.memoryReport(),
        m_hydro_amr.memoryReport(),
        m_source.memoryReport(),
        profiler.retainedEventMemoryReport()};
    core::MemoryReport merged_runtime_memory_report =
        core::mergeMemoryReports(runtime_reports);
    if (m_services.memory_governor != nullptr) {
      const std::uint64_t governor_baseline_bytes =
          core::memoryReportBaselineOwnedBytes(merged_runtime_memory_report);
      m_services.memory_governor->setBaselineOwnedBytes(governor_baseline_bytes);
      core::attachMemoryGovernorSnapshot(
          merged_runtime_memory_report, *m_services.memory_governor);
    }
    core::attachProcessMemoryObservation(
        merged_runtime_memory_report, core::observeProcessMemory());
    attachDistributedMemoryTelemetry(
        merged_runtime_memory_report, m_services);
    profiler.setMemoryReport(std::move(merged_runtime_memory_report));
    state.metadata.step_index = integrator_state.step_index;
    state.metadata.scale_factor = integrator_state.current_scale_factor;
    ensureSchedulersCoverState(state, particle_scheduler, gas_cell_scheduler);
    static_cast<void>(updateAdaptiveTimeBins(
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
        false));
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
      install_authoritative_domain_geometry();
      if (m_services.console_reporter != nullptr) {
        m_services.console_reporter->emitDecomposition(
            integrator_state.step_index, m_gravity.decompositionEpoch());
      }
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

    if (m_services.console_reporter != nullptr) {
      RuntimeConsoleStepStatus console_status;
      console_status.step_index = integrator_state.step_index;
      console_status.t_code = integrator_state.current_time_code;
      console_status.dt_time_code = accepted_dt_time_code_for_console;
      if (mode_policy.cosmological_comoving_frame) {
        console_status.a_scale = integrator_state.current_scale_factor;
        if (std::isfinite(integrator_state.current_scale_factor) &&
            integrator_state.current_scale_factor > 0.0) {
          console_status.redshift =
              1.0 / integrator_state.current_scale_factor - 1.0;
        }
      }
      console_status.active_particle_count = global_active_particle_count;
      console_status.total_particle_count = global_particle_count;
      console_status.active_cell_count = global_active_cell_count;
      console_status.total_cell_count = global_cell_count;
      console_status.wall_step_seconds = std::chrono::duration<double>(
          RuntimeConsoleReporter::Clock::now() - console_step_begin).count();
      const bool pm_refreshed =
          m_gravity.longRangeRefreshCount() > pm_refresh_count_before;
      const bool pm_reused =
          m_gravity.longRangeReuseCount() > pm_reuse_count_before;
      console_status.pm_activity = pm_refreshed
          ? "refresh"
          : (pm_reused ? "reuse" : "none");
      const bool snapshot_flushed =
          console_snapshot_requested && !pending_output.snapshot_due;
      const bool restart_flushed =
          console_restart_requested && !pending_output.checkpoint_due;
      if (snapshot_flushed && restart_flushed) {
        console_status.output_activity = "snapshot+restart";
      } else if (snapshot_flushed) {
        console_status.output_activity = "snapshot";
      } else if (restart_flushed) {
        console_status.output_activity = "restart";
      } else {
        console_status.output_activity = "none";
      }
      if (const core::MemoryReport* memory_report = profiler.memoryReport();
          memory_report != nullptr && memory_report->governor_snapshot.has_value()) {
        console_status.memory_headroom_bytes =
            memory_report->governor_snapshot->headroom_bytes;
        console_status.memory_pressure = std::string(core::memoryPressureLabel(
            memory_report->governor_snapshot->pressure));
      }
      m_services.console_reporter->emitStep(console_status);
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

void TimeCoordinator::runHierarchicalSegment(
    const core::SimulationConfig& config, const ReferenceWorkflowOptions& options,
    core::SimulationState& state, const core::LambdaCdmBackground* background,
    std::vector<std::uint64_t>& expected_ids, ReferenceWorkflowReport& report,
    core::ProfilerSession& profiler, const core::ModePolicy& mode_policy) {
  auto& particles = m_time_state.m_particle_scheduler;
  auto& cells = m_time_state.m_gas_cell_scheduler;
  auto& integrator = m_time_state.m_integrator_state;
  auto& pending = m_time_state.m_pending_output;
  const auto& mpi = m_services.mpi_context;
  std::exception_ptr eligibility_failure;
  try {
    if (state.cells.size() != 0U || state.star_particles.size() != 0U ||
        state.black_holes.size() != 0U || state.tracers.size() != 0U ||
        !state.hasHomogeneousDmoMetadata() || particles.maxBin() == 0U ||
        particles.maxBin() != cells.maxBin()) {
      throw std::logic_error("hierarchical workflow requires compact collisionless DMO on one scheduler hierarchy");
    }
  } catch (...) { eligibility_failure = std::current_exception(); }
  FailureCoordinator(m_services).rethrowCollectiveFailure(eligibility_failure, "hierarchical DMO eligibility");
  core::TransientStepWorkspace workspace(m_services.memory_governor);
  const auto owned_runtime_bytes = [&]() {
    const std::array reports{core::collectSimulationMemoryReport(state, &workspace),
        core::collectSchedulerMemoryReport(particles, cells), m_gravity.memoryReport(),
        m_hydro_amr.memoryReport(), m_source.memoryReport(), profiler.retainedEventMemoryReport()};
    return core::memoryReportBaselineOwnedBytes(core::mergeMemoryReports(reports));
  };
  const auto reconcile_runtime = [&]() {
    core::MemoryReport memory;
    std::exception_ptr memory_failure;
    try {
      const std::array reports{core::collectSimulationMemoryReport(state, &workspace),
          core::collectSchedulerMemoryReport(particles, cells), m_gravity.memoryReport(),
          m_hydro_amr.memoryReport(), m_source.memoryReport(), profiler.retainedEventMemoryReport()};
      memory = core::mergeMemoryReports(reports);
      if (m_services.memory_governor != nullptr) {
        m_services.memory_governor->setBaselineOwnedBytes(core::memoryReportBaselineOwnedBytes(memory));
        core::attachMemoryGovernorSnapshot(memory, *m_services.memory_governor);
      }
      core::attachProcessMemoryObservation(memory, core::observeProcessMemory());
    } catch (...) { memory_failure = std::current_exception(); }
    FailureCoordinator(m_services).rethrowCollectiveFailure(memory_failure, "hierarchical local memory report preparation");
    attachDistributedMemoryTelemetry(memory, m_services);
    profiler.setMemoryReport(std::move(memory));
  };
  const auto install_geometry = [&]() {
    std::exception_ptr failure;
    try {
      core::MemoryReservation reservation;
      {
        const auto leaves = m_migration_balance.authoritativeTopDomainLeaves(state, m_gravity.decompositionEpoch(), &reservation);
        m_gravity.installAuthoritativeTopDomainLeaves(leaves, state.gravitySourceGeneration());
      }
      if (reservation.valid()) reservation.reconcileBaselineOwnedAndRelease(owned_runtime_bytes());
    } catch (...) { failure = std::current_exception(); }
    FailureCoordinator(m_services).rethrowCollectiveFailure(failure, "hierarchical domain geometry installation");
  };
  install_geometry();
  const std::uint64_t requested_steps = options.max_steps_override > 0U
      ? options.max_steps_override : static_cast<std::uint64_t>(config.numerics.max_global_steps);
  if (integrator.step_index > std::numeric_limits<std::uint64_t>::max() - requested_steps) {
    throw std::overflow_error("hierarchical coarse segment step count overflow");
  }
  const std::uint64_t run_start_step = integrator.step_index;
  const std::uint64_t final_step = integrator.step_index + requested_steps;
  while (integrator.step_index < final_step && integrator.current_time_code < config.numerics.t_code_end) {
    const auto console_begin = RuntimeConsoleReporter::Clock::now();
    const double block_a_begin = integrator.current_scale_factor;
    std::exception_ptr preparation_failure;
    try {
      workspace.prepareGravityParticleIndexScratch(state.particles.size());
      for (std::size_t row = 0U; row < state.particles.size(); ++row) {
        workspace.gravity_particle_index_scratch[row] = core::checkedIntegralNarrow<std::uint32_t>(row, "hierarchical identity row");
      }
      core::RetainedCapacityTransaction mirror_plan(owned_runtime_bytes);
      mirror_plan.add(state.particles.time_bin, state.particles.size());
      mirror_plan.execute(m_services.memory_governor, core::MemoryClass::kPhaseResident, "time.hierarchical_bin_mirror");
      state.particles.time_bin.resize(state.particles.size(), 0U);
    } catch (...) { preparation_failure = std::current_exception(); }
    FailureCoordinator(m_services).rethrowCollectiveFailure(preparation_failure, "hierarchical workspace admission");
    const auto all_rows = std::span<const std::uint32_t>(workspace.gravity_particle_index_scratch);
    // Mandatory fresh synchronized endpoint bootstrap in BOTH uninterrupted
    // and restarted runs. Existing total force history remains the MAC scale;
    // derived split caches/PM meshes need no new restart payload.
    core::StepContext sync{.state = state, .integrator_state = integrator,
        .active_set = {.particle_indices = all_rows, .cell_indices = {}},
        .active_gravity_particles = {}, .has_active_gravity_particles = false,
        .workspace = &workspace, .cosmology_background = background,
        .mode_policy = &mode_policy, .profiler_session = &profiler,
        .timeline_step = {.time_begin_code = integrator.current_time_code,
            .time_end_code = integrator.current_time_code,
            .scale_factor_begin = integrator.current_scale_factor,
            .scale_factor_end = integrator.current_scale_factor},
        .boundary = {.kind = core::StepBoundaryKind::kGlobalSynchronizationPoint},
        .pm_refresh_directive = {.force_refresh_surface = true, .cadence_opportunity_allowed = true,
            .sync_event_requested = true, .force_evaluation_scale_factor = integrator.current_scale_factor,
            .reason = core::PmRefreshDirective::Reason::kScheduledForceRefreshStage},
        .stage = core::IntegrationStage::kForceRefresh};
    sync.hierarchical_kdk = {.enabled = true, .include_long_range_force = true,
        .force_only_synchronization = true, .coarse_source_generation = state.gravitySourceGeneration(),
        .evaluation_tick = particles.currentTick()};
    dispatchStage(sync, false);
    const auto& d = sync.pm_refresh_directive;
    if (!d.has_sync_event || !d.solver_executed || !d.refresh_long_range_field) {
      throw std::logic_error("hierarchical synchronized bootstrap did not commit a fresh split force");
    }
    const core::PmSyncEvent event{d.gravity_kick_opportunity, d.refresh_long_range_field,
        d.field_version, d.last_refresh_opportunity, d.field_built_step_index, d.field_built_scale_factor};
    integrator.pm_sync_state.commitKickOpportunity(event);
    integrator.pm_sync_state.commitRefresh(event);
    integrator.pm_long_range_field_valid = true;
    integrator.pm_source_generation = state.gravitySourceGeneration();
    const auto ax = m_gravity.particleAccelX();
    const auto ay = m_gravity.particleAccelY();
    const auto az = m_gravity.particleAccelZ();
    const auto softening = speciesSofteningByTag(config);
    const double eps = softening.epsilon_comoving_by_species[static_cast<std::size_t>(core::ParticleSpecies::kDarkMatter)];
    const auto criterion = [&](std::uint32_t row) {
      return core::computeComovingGravityTimeStep({.softening_length_comoving_code = std::max(eps, 1.0e-12),
          .scale_free_acceleration_magnitude_code = std::sqrt(ax[row] * ax[row] + ay[row] * ay[row] + az[row] * az[row]),
          .scale_factor = integrator.current_scale_factor}, 0.2);
    };
    double local_min = std::numeric_limits<double>::infinity();
    preparation_failure = {};
    try {
      if (ax.size() != state.particles.size() || ay.size() != ax.size() || az.size() != ax.size()) {
        throw std::logic_error("hierarchical timestep criteria lack synchronized total force rows");
      }
      for (const auto row : all_rows) local_min = std::min(local_min, criterion(row));
    } catch (...) { preparation_failure = std::current_exception(); }
    FailureCoordinator(m_services).rethrowCollectiveFailure(preparation_failure, "hierarchical timestep criteria");
    const auto periods = particles.binPeriodTicks(particles.maxBin());
    double coarse_dt = mpi.allreduceMinDouble(local_min) * static_cast<double>(periods);
    if (background != nullptr) coarse_dt = std::min(coarse_dt, core::computeCosmologyExpansionTimeStep(
        *background, integrator.current_scale_factor, config.numerics.cosmology_max_delta_ln_a,
        config.numerics.cosmology_max_hubble_time_fraction, integrator.time_si_per_code));
    if (options.dt_time_code > 0.0) coarse_dt = std::min(coarse_dt, options.dt_time_code);
    double boundary_time = config.numerics.t_code_end;
    if (pending.snapshot_interval_time_code > 0.0 && pending.next_snapshot_time_code > integrator.current_time_code) {
      boundary_time = std::min(boundary_time, pending.next_snapshot_time_code);
    }
    coarse_dt = std::min(coarse_dt, boundary_time - integrator.current_time_code);
    const double quantum = coarse_dt / static_cast<double>(periods);
    if (!std::isfinite(quantum) || quantum <= 0.0 || integrator.current_time_code + quantum <= integrator.current_time_code) {
      throw std::runtime_error("hierarchical quantum is not finite positive representable");
    }
    integrator.dt_time_code = quantum;
    const core::TimeStepLimits limits{.min_dt_time_code = quantum,
        .max_dt_time_code = coarse_dt, .max_bin = particles.maxBin()};
    preparation_failure = {};
    try {
      for (const auto row : all_rows) particles.submitCandidateTimeStep(row, criterion(row), limits,
          core::TimeStepCandidateSource::kGravityAcceleration, "hierarchical_synchronized_gravity");
      core::RetainedCapacityTransaction bins_plan(owned_runtime_bytes);
      particles.planSynchronizedCandidateCapacity(bins_plan);
      bins_plan.execute(m_services.memory_governor, core::MemoryClass::kPhaseResident, "time.hierarchical_bin_membership");
      particles.commitSynchronizedCandidates();
      syncTimeBinsFromSchedulers(particles, cells, state);
    } catch (...) { preparation_failure = std::current_exception(); }
    FailureCoordinator(m_services).rethrowCollectiveFailure(preparation_failure, "hierarchical synchronized rung assignment");
    preparation_failure = {};
    try {
      internal::latchOutputRequestForCompletedStep(config, options, integrator.step_index + 1U,
        integrator.current_time_code + coarse_dt, pending);
      profiler.recordEvent(core::RuntimeEvent{.event_kind = "time.hierarchical_block", .severity = core::RuntimeEventSeverity::kInfo,
        .subsystem = "core.time", .step_index = integrator.step_index,
        .simulation_time_code = integrator.current_time_code, .scale_factor = integrator.current_scale_factor,
        .message = "source-implemented synchronized power-of-two DMO KDK block; qualification pending",
        .payload = {{"quantum_time_code", formatRuntimeDouble(quantum)},
            {"coarse_interval_time_code", formatRuntimeDouble(coarse_dt)}, {"fine_ticks", std::to_string(periods)},
            {"pm_policy", "coarse_endpoint_half_kicks"}, {"rung_assignment", "fixed_within_block"}}});
    } catch (...) { preparation_failure = std::current_exception(); }
    FailureCoordinator(m_services).rethrowCollectiveFailure(preparation_failure, "hierarchical output and event preparation");
    core::MemoryReservation timeline_reservation;
    preparation_failure = {};
    try {
      if (m_services.memory_governor != nullptr) {
        timeline_reservation = m_services.memory_governor->reserve(core::MemoryClass::kPhaseResident,
            core::checkedSizeMultiply(core::k_hierarchical_timeline_capacity, sizeof(double), "hierarchical timeline bytes"),
            "time.hierarchical_timeline");
        timeline_reservation.commit();
      }
    } catch (...) { preparation_failure = std::current_exception(); }
    FailureCoordinator(m_services).rethrowCollectiveFailure(preparation_failure, "hierarchical timeline admission");
    parallel::DecompositionRuntimeMeasurements block_work{};
    m_lifecycle.executeHierarchicalBlockWithDispatcher(state, integrator, particles, cells,
        [&](core::StepContext& context, bool safe_output) {
          if (context.stage == core::IntegrationStage::kAnalysisHooks && context.hierarchical_kdk.synchronization_end) {
            preparation_failure = {};
            try {
              if (background != nullptr &&
                  (std::log(integrator.current_scale_factor / block_a_begin) >
                      config.numerics.cosmology_max_delta_ln_a * (1.0 + 1.0e-10) ||
                   ((config.numerics.integrator_time_variable == core::IntegratorTimeVariable::kScaleFactor ||
                     config.numerics.integrator_time_variable == core::IntegratorTimeVariable::kLogScaleFactor) &&
                    integrator.current_scale_factor > config.numerics.a_end + 1.0e-12 * std::max(1.0, config.numerics.a_end)))) {
                throw std::runtime_error("hierarchical block crossed the configured cosmological interval or endpoint");
              }
            } catch (...) { preparation_failure = std::current_exception(); }
            FailureCoordinator(m_services).rethrowCollectiveFailure(preparation_failure, "hierarchical cosmological endpoint");
            std::exception_ptr boundary_failure;
            try { reconcile_runtime(); } catch (...) { boundary_failure = std::current_exception(); }
            FailureCoordinator(m_services).rethrowCollectiveFailure(boundary_failure, "hierarchical pre-migration memory reconciliation");
            const bool migrated = m_migration_balance.rebalance(state, particles, cells,
                block_work, context.active_set.particle_indices,
                expected_ids, integrator.step_index);
            if (migrated) {
              m_gravity.commitParticleDecompositionChange();
              install_geometry();
              // Ownership commit invalidates the former target-row span. Analysis
              // and output consume coherent state, never the departed active view.
              context.active_set = {};
            }
            boundary_failure = {};
            try {
              syncTimeBinsFromSchedulers(particles, cells, state);
              reconcile_runtime();
            } catch (...) { boundary_failure = std::current_exception(); }
            FailureCoordinator(m_services).rethrowCollectiveFailure(boundary_failure, "hierarchical post-migration mirrors");
          }
          dispatchStage(context, safe_output);
          if (context.stage == core::IntegrationStage::kForceRefresh) {
            std::exception_ptr work_failure;
            try {
              const auto current = m_gravity.lastRuntimeDecompositionMeasurements();
              const auto add = [](std::uint64_t a, std::uint64_t b) {
                return core::checkedMemoryBytesAdd(a, b, "hierarchical block work counter");
              };
              block_work.tree_pair_evaluations_recent = add(block_work.tree_pair_evaluations_recent, current.tree_pair_evaluations_recent);
              block_work.incoming_tree_pair_evaluations_recent = add(block_work.incoming_tree_pair_evaluations_recent, current.incoming_tree_pair_evaluations_recent);
              block_work.tree_remote_request_bytes_recent = add(block_work.tree_remote_request_bytes_recent, current.tree_remote_request_bytes_recent);
              block_work.pm_mesh_cells_touched_recent = add(block_work.pm_mesh_cells_touched_recent, current.pm_mesh_cells_touched_recent);
              block_work.pm_fft_transpose_bytes_recent = add(block_work.pm_fft_transpose_bytes_recent, current.pm_fft_transpose_bytes_recent);
              block_work.ghost_exchange_bytes_recent = add(block_work.ghost_exchange_bytes_recent, current.ghost_exchange_bytes_recent);
              block_work.tree_wall_ms_recent += current.tree_wall_ms_recent;
              block_work.pm_wall_ms_recent += current.pm_wall_ms_recent;
              block_work.has_measurements = true;
              block_work.has_spatial_tree_work = current.has_spatial_tree_work;
              block_work.spatial_tree_work_per_target = current.spatial_tree_work_per_target;
            } catch (...) { work_failure = std::current_exception(); }
            FailureCoordinator(m_services).rethrowCollectiveFailure(work_failure, "hierarchical block work aggregation");
          }
        }, background, workspace, &mode_policy, &profiler,
        [&](std::exception_ptr failure, std::string_view phase) {
          FailureCoordinator(m_services).rethrowCollectiveFailure(failure, phase);
        });
    if (timeline_reservation.valid()) timeline_reservation.release();
    preparation_failure = {};
    try { reconcile_runtime(); } catch (...) { preparation_failure = std::current_exception(); }
    FailureCoordinator(m_services).rethrowCollectiveFailure(preparation_failure, "hierarchical block memory report");
    profiler.counters().addCount("time.hierarchical_coarse_blocks");
    profiler.counters().addCount("time.hierarchical_fine_drifts", periods);
    const auto global_particles = mpi.allreduceSumUint64(static_cast<std::uint64_t>(state.particles.size()));
    if (m_services.console_reporter != nullptr) {
      RuntimeConsoleStepStatus status;
      status.step_index = integrator.step_index;
      status.t_code = integrator.current_time_code;
      status.dt_time_code = coarse_dt;
      status.a_scale = integrator.current_scale_factor;
      status.redshift = integrator.current_redshift;
      status.active_particle_count = global_particles;
      status.total_particle_count = global_particles;
      status.wall_step_seconds = std::chrono::duration<double>(RuntimeConsoleReporter::Clock::now() - console_begin).count();
      status.pm_activity = "coarse_endpoints";
      status.output_activity = "coarse_sync";
      m_services.console_reporter->emitStep(status);
    }
  }
  report.completed_steps = integrator.step_index - run_start_step;
  report.final_time_code = integrator.current_time_code;
  report.final_scale_factor = integrator.current_scale_factor;
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

double TimeCoordinator::updateAdaptiveTimeBins(
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
  if (config.physics.enable_star_formation &&
      config.physics.star_formation_model ==
          core::StarFormationModelKind::kEffectiveMultiphaseTngLike &&
      m_effective_eos_table == nullptr) {
    const core::UnitSystem runtime_units = core::makeUnitSystem(
        config.units.length_unit, config.units.mass_unit, config.units.velocity_unit);
    m_effective_eos_table = std::make_shared<physics::EffectiveMultiphaseEosTable>(
        physics::makeEffectiveMultiphaseEosConfig(config.physics),
        runtime_units,
        physics::makeEffectiveIsmReferenceCoolingProvider(config.physics));
  }
  return updateAdaptiveTimeBinFamilies(
      state,
      particle_scheduler,
      gas_cell_scheduler,
      integrator_state,
      config,
      m_effective_eos_table.get(),
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

void TimeCoordinator::coordinatePmRefreshDirective(core::StepContext& context) {
  if (context.stage != core::IntegrationStage::kGravityKickPre &&
      context.stage != core::IntegrationStage::kForceRefresh &&
      context.stage != core::IntegrationStage::kGravityKickPost) {
    return;
  }

  if (context.hierarchical_kdk.enabled) {
    auto& d = context.pm_refresh_directive;
    const auto& mpi = m_services.mpi_context;
    const auto votes = mpi.allreduceSumUint64(d.sync_event_requested ? 1ULL : 0ULL);
    if (votes != 0U && votes != static_cast<std::uint64_t>(mpi.worldSize())) {
      throw std::logic_error("hierarchical PM synchronization directive diverged across ranks");
    }
    if (d.sync_event_requested) {
      std::exception_ptr directive_failure;
      try {
        if (!d.cadence_opportunity_allowed || !context.boundary.pm_refresh_allowed ||
            !context.hierarchical_kdk.include_long_range_force) {
          throw std::logic_error("hierarchical PM refresh requested outside a synchronized split endpoint");
        }
      } catch (...) { directive_failure = std::current_exception(); }
      FailureCoordinator(m_services).rethrowCollectiveFailure(directive_failure, "hierarchical PM endpoint eligibility");
      const bool all_valid = mpi.allreduceSumUint64(context.integrator_state.pm_long_range_field_valid ? 1ULL : 0ULL) ==
          static_cast<std::uint64_t>(mpi.worldSize());
      directive_failure = {};
      try {
        core::materializePmRefreshDirective(context, all_valid);
        if (!d.refresh_long_range_field) throw std::logic_error("hierarchical PM endpoint cadence must refresh");
      } catch (...) { directive_failure = std::current_exception(); }
      FailureCoordinator(m_services).rethrowCollectiveFailure(directive_failure, "hierarchical PM endpoint materialization");
    }
    return;
  }

  auto& directive = context.pm_refresh_directive;
  const auto& mpi_context = m_services.mpi_context;
  const std::uint64_t world_size = static_cast<std::uint64_t>(
      std::max(mpi_context.worldSize(), 1));

  const auto materialize_single_rank = [&]() {
    if (!directive.sync_event_requested) {
      return;
    }
    if (!directive.cadence_opportunity_allowed) {
      if (directive.reason == core::PmRefreshDirective::Reason::kSourceMutationForceRefresh) {
        throw std::runtime_error(
            "gravity source state changed after the scheduled PM refresh on a boundary where long-range repair is illegal");
      }
      throw std::runtime_error(
          "PM cadence event requested outside an integrator-authorized refresh opportunity");
    }
    core::materializePmRefreshDirective(
        context, context.integrator_state.pm_long_range_field_valid);
  };

  if (!mpi_context.isEnabled() || world_size <= 1U) {
    materialize_single_rank();
    return;
  }

  // All ranks reach this fixed-size reduction before any branch-dependent
  // TreePM work. Local predicates are evidence only; the values below define
  // the one distributed refresh decision for this stage.
  std::array<std::uint64_t, 5> votes{
      context.integrator_state.pm_refresh_enabled ? 1ULL : 0ULL,
      context.boundary.pm_refresh_allowed ? 1ULL : 0ULL,
      context.boundary.local_substep ? 1ULL : 0ULL,
      context.integrator_state.pm_long_range_field_valid ? 1ULL : 0ULL,
      context.integrator_state.pm_source_generation != context.state.gravitySourceGeneration() ? 1ULL : 0ULL,
  };
  mpi_context.allreduceSumUint64sInPlace(votes);

  const bool pm_enablement_consistent = votes[0] == 0U || votes[0] == world_size;
  const bool global_pm_enabled = votes[0] == world_size;
  const bool global_refresh_allowed = votes[1] == world_size;
  const bool any_local_substep = votes[2] != 0U;
  const bool global_field_valid = votes[3] == world_size;
  const bool field_validity_mixed = votes[3] != 0U && votes[3] != world_size;
  const bool global_source_changed = votes[4] != 0U;

  std::exception_ptr coordination_failure;
  try {
    if (!pm_enablement_consistent) {
      throw std::runtime_error(
          "TreePM PM-refresh enablement diverged across ranks before cadence coordination");
    }
    if (!global_pm_enabled && !global_field_valid) {
      throw std::runtime_error(
          "TreePM long-range PM field is invalid while PM refresh is disabled");
    }

    // Rebuild the directive from global semantics rather than preserving a
    // rank-local candidate. The event itself remains uncommitted until the
    // gravity callback succeeds and the core integrator accepts it.
    directive.has_sync_event = false;
    directive.refresh_long_range_field = false;
    directive.solver_executed = false;
    directive.sync_stage = core::PmSyncStage::kNone;
    directive.gravity_kick_opportunity = 0U;
    directive.field_version = 0U;
    directive.last_refresh_opportunity = 0U;
    directive.field_built_step_index = 0U;
    directive.field_built_scale_factor = 1.0;

    switch (context.stage) {
      case core::IntegrationStage::kGravityKickPre: {
        const bool global_initial_bootstrap_needed =
            global_pm_enabled && !global_field_valid;
        const bool global_initial_bootstrap_allowed =
            !any_local_substep ||
            (context.integrator_state.step_index == 0U &&
             global_initial_bootstrap_needed);
        directive.initial_cache_bootstrap_allowed =
            global_initial_bootstrap_allowed;
        directive.sync_event_requested = global_initial_bootstrap_needed;
        directive.cadence_opportunity_allowed =
            global_initial_bootstrap_needed && global_initial_bootstrap_allowed;
        if (global_initial_bootstrap_needed) {
          if (!global_initial_bootstrap_allowed) {
            throw std::runtime_error(
                "TreePM distributed initial PM bootstrap is not legal at this integration boundary");
          }
          directive.reason =
              core::PmRefreshDirective::Reason::kInitialForceBootstrap;
          directive.force_evaluation_scale_factor =
              context.timeline_step.scale_factor_begin;
        } else {
          directive.reason = core::PmRefreshDirective::Reason::kNone;
        }
        break;
      }
      case core::IntegrationStage::kForceRefresh: {
        directive.force_refresh_surface = true;
        directive.requires_predicted_inactive_sources = any_local_substep;
        directive.reason =
            core::PmRefreshDirective::Reason::kScheduledForceRefreshStage;
        directive.force_evaluation_scale_factor =
            context.timeline_step.scale_factor_end;
        directive.cadence_opportunity_allowed =
            global_pm_enabled && global_refresh_allowed;
        directive.sync_event_requested = directive.cadence_opportunity_allowed;
        if (global_pm_enabled && !global_field_valid &&
            !global_refresh_allowed) {
          throw std::runtime_error(
              "TreePM distributed PM field is invalid but no rank-global refresh boundary is available");
        }
        break;
      }
      case core::IntegrationStage::kGravityKickPost: {
        directive.force_refresh_surface = global_source_changed;
        directive.cadence_opportunity_allowed = false;
        directive.sync_event_requested = false;
        directive.reason = core::PmRefreshDirective::Reason::kNone;
        if (global_source_changed) {
          if (!global_pm_enabled) {
            throw std::runtime_error(
                "TreePM gravity sources changed while PM refresh is disabled");
          }
          if (!global_refresh_allowed) {
            throw std::runtime_error(
                "gravity source state changed after the scheduled PM refresh on a boundary where long-range repair is illegal");
          }
          directive.reason =
              core::PmRefreshDirective::Reason::kSourceMutationForceRefresh;
          directive.force_evaluation_scale_factor =
              context.timeline_step.scale_factor_end;
          directive.cadence_opportunity_allowed = true;
          directive.sync_event_requested = true;
        } else if (global_pm_enabled && !global_field_valid) {
          if (field_validity_mixed) {
            throw std::runtime_error(
                "TreePM long-range PM validity diverged across ranks without an authorized refresh request");
          }
          throw std::runtime_error(
              "TreePM post-kick reached an invalid long-range PM field without a source-mutation refresh request");
        }
        break;
      }
      default:
        break;
    }

    if (directive.sync_event_requested) {
      core::materializePmRefreshDirective(context, global_field_valid);
    }
  } catch (...) {
    coordination_failure = std::current_exception();
  }

  FailureCoordinator(m_services).rethrowCollectiveFailure(
      coordination_failure, "TreePM PM-refresh rank-global directive coordination");
}

void TimeCoordinator::dispatchStage(
    core::StepContext& context,
    bool require_output_safe_boundary) {
  (void)require_output_safe_boundary;
  context.particle_scheduler = &m_time_state.m_particle_scheduler;
  context.gas_cell_scheduler = &m_time_state.m_gas_cell_scheduler;
  context.newly_created_particle_ids = &m_newly_created_particle_ids;
  const std::string stage_name =
      "stage." + std::string(core::integrationStageName(context.stage));
  COSMOSIM_PROFILE_SCOPE(context.profiler_session, stage_name);
  if (context.profiler_session != nullptr) {
    context.profiler_session->counters().addCount(
        stage_name + ".invocations", 1U);
  }

  coordinatePmRefreshDirective(context);

  SimulationRuntimeEpochSource epoch_source(
      context.state, m_time_state.m_particle_scheduler, context.integrator_state);
  internal::RuntimeStageResourceBundle stage_resources(context);
  AnalysisStageView audit_view(
      RuntimeResourceLease(epoch_source, k_state_stage_epochs), stage_resources);
  if (!context.hierarchical_kdk.force_only_synchronization) {
    m_execution_plan.executeAuditStage(context.stage, audit_view);
  }
  if (context.hierarchical_kdk.enabled &&
      (context.stage == core::IntegrationStage::kHydroUpdate ||
       context.stage == core::IntegrationStage::kSourceTerms ||
       ((context.stage == core::IntegrationStage::kAnalysisHooks || context.stage == core::IntegrationStage::kOutputCheck) &&
        !context.hierarchical_kdk.synchronization_end))) return;

  switch (context.stage) {
    case core::IntegrationStage::kGravityKickPre:
    case core::IntegrationStage::kForceRefresh:
    case core::IntegrationStage::kGravityKickPost: {
      GravityStageView view(
          RuntimeResourceLease(epoch_source, k_state_stage_epochs), stage_resources);
      m_execution_plan.executeStage(context.stage, view);
      break;
    }
    case core::IntegrationStage::kDrift: {
      DriftParticleStageView view(
          RuntimeResourceLease(epoch_source, k_particle_stage_epochs), stage_resources);
      m_execution_plan.executeStage(context.stage, view);
      break;
    }
    case core::IntegrationStage::kHydroUpdate: {
      HydroAmrStageView view(
          RuntimeResourceLease(epoch_source, k_state_stage_epochs), stage_resources);
      m_execution_plan.executeStage(context.stage, view);
      break;
    }
    case core::IntegrationStage::kSourceTerms: {
      SourceMutationStageView view(
          RuntimeResourceLease(epoch_source, k_state_stage_epochs), stage_resources);
      m_execution_plan.executeStage(context.stage, view);
      break;
    }
    case core::IntegrationStage::kAnalysisHooks: {
      AnalysisStageView view(
          RuntimeResourceLease(epoch_source, k_state_stage_epochs), stage_resources);
      m_execution_plan.executeStage(context.stage, view);
      break;
    }
    case core::IntegrationStage::kOutputCheck: {
      OutputRestartStageView view(
          RuntimeResourceLease(epoch_source, k_state_stage_epochs), stage_resources);
      m_execution_plan.executeStage(context.stage, view);
      break;
    }
  }
}

}  // namespace cosmosim::workflows
