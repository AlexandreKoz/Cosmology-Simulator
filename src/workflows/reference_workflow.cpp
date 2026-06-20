#include "cosmosim/workflows/reference_workflow.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <iterator>
#include <numeric>
#include <optional>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <cstring>
#include <type_traits>
#include <utility>
#include <vector>

#include "cosmosim/analysis/diagnostics.hpp"
#include "workflows/internal/cartesian_gas_cell_layout.hpp"
#include "cosmosim/amr/amr_hydro_orchestrator.hpp"
#include "cosmosim/core/build_config.hpp"
#include "cosmosim/core/cosmology.hpp"
#include "cosmosim/core/constants.hpp"
#include "cosmosim/core/cuda_runtime.hpp"
#include "cosmosim/core/profiling.hpp"
#include "cosmosim/core/simulation_mode.hpp"
#include "cosmosim/core/time_integration.hpp"
#include "cosmosim/core/units.hpp"
#include "cosmosim/gravity/tree_pm_coupling.hpp"
#include "cosmosim/hydro/hydro_boundary_conditions.hpp"
#include "cosmosim/hydro/hydro_cartesian_patch.hpp"
#include "cosmosim/hydro/hydro_core_solver.hpp"
#include "cosmosim/hydro/hydro_reconstruction.hpp"
#include "cosmosim/hydro/hydro_riemann.hpp"
#include "cosmosim/io/ic_reader.hpp"
#include "cosmosim/io/restart_checkpoint.hpp"
#include "cosmosim/io/snapshot_hdf5.hpp"
#include "cosmosim/physics/black_hole_agn.hpp"
#include "cosmosim/physics/star_formation.hpp"

#if COSMOSIM_ENABLE_HDF5
#include <hdf5.h>
#endif

#if COSMOSIM_ENABLE_MPI
#include <mpi.h>
#endif

namespace cosmosim::workflows {
namespace {

constexpr double k_gamma_adiabatic = 5.0 / 3.0;
constexpr double k_pressure_floor = 1.0e-10;
constexpr double k_density_floor = 1.0e-10;
constexpr std::size_t k_default_generated_particle_axis = 6;

constexpr std::array<core::IntegrationStage, core::integrationStageCount()> k_all_integration_stages = {
    core::IntegrationStage::kGravityKickPre,
    core::IntegrationStage::kDrift,
    core::IntegrationStage::kForceRefresh,
    core::IntegrationStage::kHydroUpdate,
    core::IntegrationStage::kSourceTerms,
    core::IntegrationStage::kGravityKickPost,
    core::IntegrationStage::kAnalysisHooks,
    core::IntegrationStage::kOutputCheck,
};
constexpr std::array<core::IntegrationStage, 1> k_drift_stage = {core::IntegrationStage::kDrift};
constexpr std::array<core::IntegrationStage, 3> k_gravity_stages = {
    core::IntegrationStage::kGravityKickPre,
    core::IntegrationStage::kForceRefresh,
    core::IntegrationStage::kGravityKickPost,
};
constexpr std::array<core::IntegrationStage, 1> k_hydro_stage = {core::IntegrationStage::kHydroUpdate};

[[nodiscard]] gravity::PmAssignmentScheme toPmAssignmentScheme(
    core::TreePmAssignmentScheme assignment_scheme) {
  switch (assignment_scheme) {
    case core::TreePmAssignmentScheme::kCic:
      return gravity::PmAssignmentScheme::kCic;
    case core::TreePmAssignmentScheme::kTsc:
      return gravity::PmAssignmentScheme::kTsc;
  }
  throw std::runtime_error("unhandled TreePm assignment scheme enum value");
}

[[nodiscard]] std::string treePmAssignmentSchemeName(core::TreePmAssignmentScheme assignment_scheme) {
  switch (assignment_scheme) {
    case core::TreePmAssignmentScheme::kCic:
      return "cic";
    case core::TreePmAssignmentScheme::kTsc:
      return "tsc";
  }
  throw std::runtime_error("unhandled TreePm assignment scheme enum value");
}



gravity::TreePmCoordinator makeRuntimeAwareTreePmCoordinator(
    const core::SimulationConfig& config,
    const gravity::PmGridShape& pm_grid_shape) {
  parallel::MpiContext mpi_context;
  mpi_context.validateExpectedWorldSizeOrThrow(config.parallel.mpi_ranks_expected);
  const auto layout = parallel::makePmSlabLayout(
      pm_grid_shape.nx,
      pm_grid_shape.ny,
      pm_grid_shape.nz,
      mpi_context.worldSize(),
      mpi_context.worldRank());
  return gravity::TreePmCoordinator(pm_grid_shape, layout, mpi_context);
}


[[nodiscard]] double newtonGCodeFromUnits(const core::UnitSystem& units) {
  return core::constants::k_newton_g_si * units.mass_si_per_code * units.timeSiPerCode() * units.timeSiPerCode() /
      (units.length_si_per_code * units.length_si_per_code * units.length_si_per_code);
}

[[nodiscard]] physics::StarFormationConfig makeRuntimeStarFormationConfig(
    const core::PhysicsConfig& physics_config,
    const core::UnitSystem& units) {
  physics::StarFormationConfig config = physics::makeStarFormationConfig(physics_config);
  config.newton_g_code = newtonGCodeFromUnits(units);
  return config;
}

[[nodiscard]] physics::BlackHoleAgnConfig makeRuntimeBlackHoleAgnConfig(
    const core::PhysicsConfig& physics_config,
    const core::UnitSystem& units) {
  physics::BlackHoleAgnConfig config = physics::makeBlackHoleAgnConfig(physics_config);
  config.proton_mass_code = config.proton_mass_si / units.mass_si_per_code;
  config.thomson_cross_section_code = config.thomson_cross_section_si /
      (units.length_si_per_code * units.length_si_per_code);
  config.newton_g_code = config.newton_g_si * units.mass_si_per_code * units.timeSiPerCode() * units.timeSiPerCode() /
      (units.length_si_per_code * units.length_si_per_code * units.length_si_per_code);
  config.speed_of_light_code = config.speed_of_light_si / units.velocity_si_per_code;
  return config;
}

[[nodiscard]] io::StochasticPersistentState buildStochasticPersistentState(
    const core::SimulationConfig& config,
    const core::IntegratorState& integrator_state,
    std::uint32_t rank_local_seed_offset) {
  io::StochasticPersistentState stochastic_state;
  if (config.physics.enable_star_formation && config.physics.sf_stochastic_spawning) {
    stochastic_state.modules.push_back(io::StochasticModulePersistentState{
        .module_name = "star_formation",
        .schema_version = 1,
        .rng_policy = "stateless_splitmix64(seed,step_index,cell_index,rank_local_seed_offset)",
        .random_seed = config.physics.sf_random_seed,
        .rank_local_seed_offset = rank_local_seed_offset,
        .last_committed_step_index = integrator_state.step_index,
        .deterministic_from_serialized_inputs = true,
    });
  }
  if (config.physics.enable_feedback && config.physics.fb_variant == core::FeedbackVariant::kStochastic) {
    stochastic_state.modules.push_back(io::StochasticModulePersistentState{
        .module_name = "stellar_feedback",
        .schema_version = 1,
        .rng_policy = "stateless_splitmix64(seed,step_index,star_index)",
        .random_seed = config.physics.fb_random_seed,
        .rank_local_seed_offset = rank_local_seed_offset,
        .last_committed_step_index = integrator_state.step_index,
        .deterministic_from_serialized_inputs = true,
    });
  }
  return stochastic_state;
}

[[nodiscard]] bool stochasticStatesEquivalent(
    io::StochasticPersistentState lhs,
    io::StochasticPersistentState rhs) {
  const auto less_by_name = [](
      const io::StochasticModulePersistentState& a,
      const io::StochasticModulePersistentState& b) {
    return a.module_name < b.module_name;
  };
  std::sort(lhs.modules.begin(), lhs.modules.end(), less_by_name);
  std::sort(rhs.modules.begin(), rhs.modules.end(), less_by_name);
  if (lhs.modules.size() != rhs.modules.size()) {
    return false;
  }
  for (std::size_t i = 0; i < lhs.modules.size(); ++i) {
    const auto& a = lhs.modules[i];
    const auto& b = rhs.modules[i];
    if (a.module_name != b.module_name || a.schema_version != b.schema_version ||
        a.rng_policy != b.rng_policy || a.random_seed != b.random_seed ||
        a.rank_local_seed_offset != b.rank_local_seed_offset ||
        a.last_committed_step_index != b.last_committed_step_index ||
        a.deterministic_from_serialized_inputs != b.deterministic_from_serialized_inputs) {
      return false;
    }
  }
  return true;
}

[[nodiscard]] bool gravityForceCachesEquivalent(
    const io::GravityForceCachePersistentState& lhs,
    const io::GravityForceCachePersistentState& rhs) {
  return lhs.valid == rhs.valid &&
      lhs.particle_id == rhs.particle_id &&
      lhs.gas_cell_id == rhs.gas_cell_id &&
      lhs.particle_accel_x_comoving == rhs.particle_accel_x_comoving &&
      lhs.particle_accel_y_comoving == rhs.particle_accel_y_comoving &&
      lhs.particle_accel_z_comoving == rhs.particle_accel_z_comoving &&
      lhs.cell_accel_x_comoving == rhs.cell_accel_x_comoving &&
      lhs.cell_accel_y_comoving == rhs.cell_accel_y_comoving &&
      lhs.cell_accel_z_comoving == rhs.cell_accel_z_comoving;
}

[[nodiscard]] std::string pmDecompositionModeName(core::PmDecompositionMode mode) {
  switch (mode) {
    case core::PmDecompositionMode::kSlab:
      return "slab";
    case core::PmDecompositionMode::kPencil:
      return "pencil";
  }
  throw std::runtime_error("unhandled PM decomposition mode enum value");
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

[[nodiscard]] bool hasSpeciesSpecificSoftening(const core::SimulationConfig& config) {
  return config.numerics.gravity_softening_gas_kpc_comoving > 0.0 ||
      config.numerics.gravity_softening_dark_matter_kpc_comoving > 0.0 ||
      config.numerics.gravity_softening_star_kpc_comoving > 0.0 ||
      config.numerics.gravity_softening_black_hole_kpc_comoving > 0.0 ||
      config.numerics.gravity_softening_tracer_kpc_comoving > 0.0;
}

[[nodiscard]] std::string describeSofteningPolicy(
    const core::SimulationConfig& config,
    const core::SimulationState* state = nullptr) {
  const bool per_species = hasSpeciesSpecificSoftening(config);
  const bool per_particle = state != nullptr && std::any_of(
      state->particle_sidecar.has_gravity_softening_override.begin(),
      state->particle_sidecar.has_gravity_softening_override.end(),
      [](std::uint8_t flag) { return flag != 0U; });
  if (per_particle && per_species) {
    return "comoving_species_plus_particle_override";
  }
  if (per_particle) {
    return "comoving_particle_override";
  }
  if (per_species) {
    return "comoving_species";
  }
  return "comoving_fixed";
}

void maybeInitializeParticleSofteningFromSpeciesPolicy(
    core::SimulationState& state,
    const core::SimulationConfig& config) {
  if (!state.particle_sidecar.gravity_softening_comoving.empty()) {
    return;
  }
  if (!hasSpeciesSpecificSoftening(config)) {
    return;
  }
  const auto by_species = speciesSofteningByTag(config);
  state.particle_sidecar.gravity_softening_comoving.resize(state.particles.size(), 0.0);
  for (std::size_t i = 0; i < state.particles.size(); ++i) {
    const std::size_t species_tag = static_cast<std::size_t>(state.particle_sidecar.species_tag[i]);
    if (species_tag < by_species.epsilon_comoving_by_species.size()) {
      state.particle_sidecar.gravity_softening_comoving[i] = by_species.epsilon_comoving_by_species[species_tag];
    } else {
      state.particle_sidecar.gravity_softening_comoving[i] = config.numerics.gravity_softening_kpc_comoving * 1.0e-3;
    }
  }
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
    SchedulerElementFamily element_family) {
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
  if (element_family == SchedulerElementFamily::kGasCells) {
    for (std::uint32_t cell_index = 0; cell_index < cell_count; ++cell_index) {
      core::TimeStepCriteriaRegistry registry;
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
      registry.registerCflHook([=](std::uint32_t) {
        return core::computeDirectionalCflTimeStep(hydro_cfl_input, 0.4);
      });
      const double ax = (cell_index < view.gas_cells.accel_x_comoving.size()) ? view.gas_cells.accel_x_comoving[cell_index] : 0.0;
      const double ay = (cell_index < view.gas_cells.accel_y_comoving.size()) ? view.gas_cells.accel_y_comoving[cell_index] : 0.0;
      const double az = (cell_index < view.gas_cells.accel_z_comoving.size()) ? view.gas_cells.accel_z_comoving[cell_index] : 0.0;
      const double amag = std::sqrt(ax * ax + ay * ay + az * az);
      const double eps = gas_index != std::numeric_limits<std::uint32_t>::max() &&
              !view.particles.gravity_softening_comoving.empty()
          ? view.particles.gravity_softening_comoving[gas_index]
          : species_softening.epsilon_comoving_by_species[static_cast<std::size_t>(core::ParticleSpecies::kGas)];
      registry.registerGravityHook([=](std::uint32_t) {
        return core::computeGravityTimeStep({.softening_length_code = std::max(eps, 1.0e-12), .acceleration_magnitude_code = amag}, 0.2);
      });
      const double cfl_dt = registry.hooks().cfl_hook(cell_index);
      const double gravity_dt = registry.hooks().gravity_hook(cell_index);
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
    }
    return;
  }

  for (std::uint32_t particle_index = 0; particle_index < particle_count; ++particle_index) {
    core::TimeStepCriteriaRegistry registry;
    const double ax = (particle_index < view.particles.accel_x_comoving.size()) ? view.particles.accel_x_comoving[particle_index] : 0.0;
    const double ay = (particle_index < view.particles.accel_y_comoving.size()) ? view.particles.accel_y_comoving[particle_index] : 0.0;
    const double az = (particle_index < view.particles.accel_z_comoving.size()) ? view.particles.accel_z_comoving[particle_index] : 0.0;
    const double amag = std::sqrt(ax * ax + ay * ay + az * az);
    const double eps = !view.particles.gravity_softening_comoving.empty()
        ? view.particles.gravity_softening_comoving[particle_index]
        : ((static_cast<std::size_t>(view.particles.species_tag[particle_index]) < species_softening.epsilon_comoving_by_species.size())
              ? species_softening.epsilon_comoving_by_species[static_cast<std::size_t>(view.particles.species_tag[particle_index])]
              : global_softening);
    registry.registerGravityHook([=](std::uint32_t) {
      return core::computeGravityTimeStep({.softening_length_code = std::max(eps, 1.0e-12), .acceleration_magnitude_code = amag}, 0.2);
    });
    scheduler.submitCandidateTimeStep(
        particle_index,
        registry.hooks().gravity_hook(particle_index),
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
  }
}

void updateAdaptiveTimeBins(
    core::SimulationState& state,
    core::HierarchicalTimeBinScheduler& scheduler,
    const core::IntegratorState& integrator_state,
    const core::SimulationConfig& config,
    std::span<const double> particle_accel_x,
    std::span<const double> particle_accel_y,
    std::span<const double> particle_accel_z,
    std::span<const double> cell_accel_x,
    std::span<const double> cell_accel_y,
    std::span<const double> cell_accel_z) {
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
      view, scheduler, integrator_state, config, SchedulerElementFamily::kParticles);
}

void updateGasCellAdaptiveTimeBins(
    core::SimulationState& state,
    core::HierarchicalTimeBinScheduler& gas_cell_scheduler,
    const core::IntegratorState& integrator_state,
    const core::SimulationConfig& config,
    std::span<const double> particle_accel_x,
    std::span<const double> particle_accel_y,
    std::span<const double> particle_accel_z,
    std::span<const double> cell_accel_x,
    std::span<const double> cell_accel_y,
    std::span<const double> cell_accel_z) {
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
      view, gas_cell_scheduler, integrator_state, config, SchedulerElementFamily::kGasCells);
}



struct GasCellMigrationRecord {
  std::uint64_t particle_id = 0;
  double center_x_comoving = 0.0;
  double center_y_comoving = 0.0;
  double center_z_comoving = 0.0;
  double mass_code = 0.0;
  std::uint8_t time_bin = 0;
  std::uint32_t patch_index = 0;
  double velocity_x_peculiar = 0.0;
  double velocity_y_peculiar = 0.0;
  double velocity_z_peculiar = 0.0;
  double density_code = 0.0;
  double pressure_code = 0.0;
  double internal_energy_code = 0.0;
  double temperature_code = 0.0;
  double sound_speed_code = 0.0;
};

static_assert(std::is_trivially_copyable_v<GasCellMigrationRecord>);

[[nodiscard]] std::string formatIndexedRankedFileStem(
    std::string_view stem,
    std::uint64_t index,
    int world_size,
    int world_rank) {
  std::ostringstream out;
  out << stem << '_' << std::setw(3) << std::setfill('0') << index;
  if (world_size > 1) {
    out << "_rank" << std::setw(3) << std::setfill('0') << world_rank;
  }
  out << ".hdf5";
  return out.str();
}

[[nodiscard]] std::unordered_map<std::uint64_t, GasCellMigrationRecord> collectLocalGasCellRecords(
    const core::SimulationState& state,
    std::span<const std::uint32_t> gas_particle_indices) {
  core::legacyRequireParticleBoundGasCellContract(state, "collectLocalGasCellRecords legacy import path");
  std::unordered_map<std::uint64_t, GasCellMigrationRecord> records;
  const auto gas_globals = state.particle_species_index.globalIndices(core::ParticleSpecies::kGas);
  std::vector<std::uint8_t> keep_mask(state.particles.size(), 0U);
  for (const std::uint32_t particle_index : gas_particle_indices) {
    if (particle_index >= state.particles.size()) {
      throw std::out_of_range("gas particle index out of range while collecting gas migration records");
    }
    keep_mask[particle_index] = 1U;
  }
  for (std::size_t cell_index = 0; cell_index < gas_globals.size(); ++cell_index) {
    const std::uint32_t particle_index = gas_globals[cell_index];
    if (keep_mask[particle_index] == 0U) {
      continue;
    }
    GasCellMigrationRecord record;
    record.particle_id = state.particle_sidecar.particle_id[particle_index];
    record.center_x_comoving = state.cells.center_x_comoving[cell_index];
    record.center_y_comoving = state.cells.center_y_comoving[cell_index];
    record.center_z_comoving = state.cells.center_z_comoving[cell_index];
    record.mass_code = state.cells.mass_code[cell_index];
    record.time_bin = state.cells.time_bin[cell_index];
    record.patch_index = state.cells.patch_index[cell_index];
    record.velocity_x_peculiar = state.gas_cells.velocity_x_peculiar[cell_index];
    record.velocity_y_peculiar = state.gas_cells.velocity_y_peculiar[cell_index];
    record.velocity_z_peculiar = state.gas_cells.velocity_z_peculiar[cell_index];
    record.density_code = state.gas_cells.density_code[cell_index];
    record.pressure_code = state.gas_cells.pressure_code[cell_index];
    record.internal_energy_code = state.gas_cells.internal_energy_code[cell_index];
    record.temperature_code = state.gas_cells.temperature_code[cell_index];
    record.sound_speed_code = state.gas_cells.sound_speed_code[cell_index];
    records.emplace(record.particle_id, record);
  }
  return records;
}

[[nodiscard]] std::vector<std::uint64_t> gasParticleIdByOldCellIndex(const core::SimulationState& state) {
  core::legacyRequireParticleBoundGasCellContract(state, "gasParticleIdByOldCellIndex legacy import path");
  std::vector<std::uint64_t> cell_particle_ids(state.cells.size(), 0ULL);
  for (std::size_t cell_index = 0; cell_index < state.cells.size(); ++cell_index) {
    cell_particle_ids[cell_index] =
        core::parentParticleIdForGasCellRow(state, static_cast<std::uint32_t>(cell_index)).value();
  }
  return cell_particle_ids;
}

void rebuildLocalGasStateFromParticleIds(
    core::SimulationState& state,
    const std::unordered_map<std::uint64_t, GasCellMigrationRecord>& gas_records_by_particle_id,
    std::span<const std::uint64_t> old_cell_particle_id) {
  const auto gas_globals = state.particle_species_index.globalIndices(core::ParticleSpecies::kGas);
  core::CellSoa rebuilt_cells;
  core::GasCellSidecar rebuilt_gas;
  rebuilt_cells.resize(gas_globals.size());
  rebuilt_gas.resize(gas_globals.size());

  std::unordered_map<std::uint64_t, std::uint32_t> new_cell_index_by_particle_id;
  new_cell_index_by_particle_id.reserve(gas_globals.size());
  for (std::size_t cell_index = 0; cell_index < gas_globals.size(); ++cell_index) {
    const std::uint64_t particle_id = state.particle_sidecar.particle_id[gas_globals[cell_index]];
    const auto found = gas_records_by_particle_id.find(particle_id);
    if (found == gas_records_by_particle_id.end()) {
      throw std::runtime_error("missing gas-cell migration record for local gas particle after ownership compaction");
    }
    const GasCellMigrationRecord& record = found->second;
    rebuilt_gas.gas_cell_id[cell_index] = particle_id;
    rebuilt_gas.parent_particle_id[cell_index] = particle_id;
    rebuilt_cells.center_x_comoving[cell_index] = record.center_x_comoving;
    rebuilt_cells.center_y_comoving[cell_index] = record.center_y_comoving;
    rebuilt_cells.center_z_comoving[cell_index] = record.center_z_comoving;
    rebuilt_cells.mass_code[cell_index] = record.mass_code;
    rebuilt_cells.time_bin[cell_index] = record.time_bin;
    rebuilt_cells.patch_index[cell_index] = record.patch_index;
    rebuilt_gas.velocity_x_peculiar[cell_index] = record.velocity_x_peculiar;
    rebuilt_gas.velocity_y_peculiar[cell_index] = record.velocity_y_peculiar;
    rebuilt_gas.velocity_z_peculiar[cell_index] = record.velocity_z_peculiar;
    rebuilt_gas.density_code[cell_index] = record.density_code;
    rebuilt_gas.pressure_code[cell_index] = record.pressure_code;
    rebuilt_gas.internal_energy_code[cell_index] = record.internal_energy_code;
    rebuilt_gas.temperature_code[cell_index] = record.temperature_code;
    rebuilt_gas.sound_speed_code[cell_index] = record.sound_speed_code;
    new_cell_index_by_particle_id.emplace(particle_id, static_cast<std::uint32_t>(cell_index));
  }

  state.cells = std::move(rebuilt_cells);
  state.gas_cells = std::move(rebuilt_gas);
  state.bumpCellIndexGeneration();

  const auto remap_host_cell = [&](std::uint32_t old_cell_index) {
    if (old_cell_index >= old_cell_particle_id.size()) {
      throw std::runtime_error(
          "rebuildLocalGasStateFromParticleIds: sidecar host_cell_index does not refer to an old gas cell");
    }
    const std::uint64_t host_particle_id = old_cell_particle_id[old_cell_index];
    const auto found = new_cell_index_by_particle_id.find(host_particle_id);
    if (found == new_cell_index_by_particle_id.end()) {
      throw std::runtime_error(
          "rebuildLocalGasStateFromParticleIds: sidecar host gas cell was removed during ownership compaction");
    }
    return found->second;
  };

  for (std::size_t row = 0; row < state.black_holes.size(); ++row) {
    state.black_holes.host_cell_index[row] = remap_host_cell(state.black_holes.host_cell_index[row]);
  }
  for (std::size_t row = 0; row < state.tracers.size(); ++row) {
    state.tracers.host_cell_index[row] = remap_host_cell(state.tracers.host_cell_index[row]);
  }
}

void compactStateToCurrentOwner(
    core::SimulationState& state,
    int world_rank) {
  if (world_rank < 0) {
    throw std::invalid_argument("compactStateToCurrentOwner requires non-negative world rank");
  }
  std::vector<std::uint32_t> kept_gas_particle_indices;
  kept_gas_particle_indices.reserve(state.cells.size());
  std::vector<std::uint32_t> outbound_indices;
  outbound_indices.reserve(state.particles.size());
  for (std::size_t particle_index = 0; particle_index < state.particles.size(); ++particle_index) {
    const bool owned_here = state.particle_sidecar.owning_rank[particle_index] == static_cast<std::uint32_t>(world_rank);
    if (owned_here) {
      if (state.particle_sidecar.species_tag[particle_index] == static_cast<std::uint32_t>(core::ParticleSpecies::kGas)) {
        kept_gas_particle_indices.push_back(static_cast<std::uint32_t>(particle_index));
      }
    } else {
      outbound_indices.push_back(static_cast<std::uint32_t>(particle_index));
    }
  }
  const auto gas_records = collectLocalGasCellRecords(state, kept_gas_particle_indices);
  const auto old_cell_particle_id = gasParticleIdByOldCellIndex(state);
  core::ParticleMigrationCommit commit;
  commit.world_rank = world_rank;
  commit.outbound_local_indices = std::move(outbound_indices);
  state.commitParticleMigration(commit);
  rebuildLocalGasStateFromParticleIds(state, gas_records, old_cell_particle_id);
}


[[nodiscard]] parallel::DecompositionConfig makeWorkflowDecompositionConfig(
    const core::SimulationConfig& config,
    int world_size) {
  parallel::DecompositionConfig decomposition_config;
  decomposition_config.world_size = world_size;
  decomposition_config.domain_x_min_comov = 0.0;
  decomposition_config.domain_x_max_comov = config.cosmology.box_size_x_mpc_comoving;
  decomposition_config.domain_y_min_comov = 0.0;
  decomposition_config.domain_y_max_comov = config.cosmology.box_size_y_mpc_comoving;
  decomposition_config.domain_z_min_comov = 0.0;
  decomposition_config.domain_z_max_comov = config.cosmology.box_size_z_mpc_comoving;
  decomposition_config.owned_particle_weight = 0.0;
  decomposition_config.active_target_weight = 0.0;
  decomposition_config.remote_tree_interaction_weight = 0.0;
  decomposition_config.work_weight = 0.0;
  decomposition_config.memory_weight = 0.0;
  decomposition_config.component_weights = parallel::DecompositionWeightCoefficients{
      .particle_count = config.parallel.decomposition_particle_count_weight,
      .gas_cell = config.parallel.decomposition_gas_cell_weight,
      .tree_interaction = config.parallel.decomposition_tree_interaction_weight,
      .pm_mesh = config.parallel.decomposition_pm_mesh_weight,
      .amr_patch = config.parallel.decomposition_amr_patch_weight,
      .active_fraction = config.parallel.decomposition_active_fraction_weight,
      .memory_pressure = config.parallel.decomposition_memory_pressure_weight,
      .gpu_occupancy = config.parallel.decomposition_gpu_occupancy_weight,
      .generic_work = config.parallel.decomposition_generic_work_weight,
  };
  return decomposition_config;
}

[[nodiscard]] parallel::DecompositionFeedbackCoefficients makeWorkflowFeedbackCoefficients(
    const core::SimulationConfig& config) {
  return parallel::DecompositionFeedbackCoefficients{
      .measured_tree_pair = config.parallel.decomposition_measured_tree_pair_weight,
      .measured_pm_cell = config.parallel.decomposition_measured_pm_cell_weight,
      .measured_amr_cell = config.parallel.decomposition_measured_amr_cell_weight,
      .measured_hydro_face = config.parallel.decomposition_measured_hydro_face_weight,
      .measured_wall_ms = config.parallel.decomposition_measured_wall_ms_weight,
  };
}

[[nodiscard]] std::uint64_t estimateParticleMemoryBytesForDecomposition(
    const core::SimulationState& state,
    std::uint32_t species_tag) {
  std::uint64_t bytes = sizeof(double) * 7U + sizeof(std::uint64_t) * 2U + sizeof(std::uint32_t) * 3U;
  if (!state.particle_sidecar.gravity_softening_comoving.empty()) {
    bytes += sizeof(double);
  }
  if (!state.particle_sidecar.has_gravity_softening_override.empty()) {
    bytes += sizeof(std::uint8_t);
  }
  if (species_tag == static_cast<std::uint32_t>(core::ParticleSpecies::kGas)) {
    bytes += sizeof(double) * 8U + sizeof(std::uint64_t) * 2U + sizeof(std::uint32_t) * 2U;
  } else if (species_tag == static_cast<std::uint32_t>(core::ParticleSpecies::kStar)) {
    bytes += sizeof(std::uint32_t) + sizeof(double) * 13U;
  } else if (species_tag == static_cast<std::uint32_t>(core::ParticleSpecies::kBlackHole)) {
    bytes += sizeof(std::uint32_t) * 2U + sizeof(double) * 8U;
  } else if (species_tag == static_cast<std::uint32_t>(core::ParticleSpecies::kTracer)) {
    bytes += sizeof(std::uint64_t) * 2U + sizeof(std::uint32_t) * 2U + sizeof(double) * 3U;
  }
  for (const core::ModuleSidecarBlock* block_ptr : state.sidecars.blocksSortedByName()) {
    const core::ModuleSidecarBlock& block = *block_ptr;
    if (!block.particle_indexed || block.row_stride_bytes == 0U) {
      continue;
    }
    const bool species_mask_requires_row = (block.required_species_mask & (1U << species_tag)) != 0U ||
        (block.requirement.kind == core::ModuleSidecarRequirementKind::kSpeciesMask &&
         (block.requirement.species_mask & (1U << species_tag)) != 0U);
    const bool predicate_may_require_row =
        (block.requirement.kind == core::ModuleSidecarRequirementKind::kGasDensityAtLeast &&
         species_tag == static_cast<std::uint32_t>(core::ParticleSpecies::kGas)) ||
        (block.requirement.kind == core::ModuleSidecarRequirementKind::kBlackHoleAccretionAtLeast &&
         species_tag == static_cast<std::uint32_t>(core::ParticleSpecies::kBlackHole)) ||
        block.requirement.kind == core::ModuleSidecarRequirementKind::kParticleFlagMask;
    if (species_mask_requires_row || predicate_may_require_row) {
      bytes += block.row_stride_bytes;
    }
  }
  return bytes;
}

[[nodiscard]] std::vector<parallel::DecompositionItem> buildRuntimeDecompositionItems(
    const core::SimulationState& state,
    const core::SimulationConfig& config,
    int world_rank,
    std::span<const std::uint32_t> active_particle_indices) {
  std::vector<std::uint8_t> active_mask(state.particles.size(), 0U);
  for (const std::uint32_t pidx : active_particle_indices) {
    if (pidx < active_mask.size()) {
      active_mask[pidx] = 1U;
    }
  }
  std::vector<std::uint32_t> patch_cell_count(state.patches.size(), 0U);
  for (std::size_t cell_index = 0; cell_index < state.cells.size(); ++cell_index) {
    const std::uint32_t patch_index = state.cells.patch_index[cell_index];
    if (patch_index < patch_cell_count.size()) {
      ++patch_cell_count[patch_index];
    }
  }

  std::vector<parallel::DecompositionItem> items;
  items.reserve(state.particles.size() + state.patches.size());
  for (std::size_t particle_index = 0; particle_index < state.particles.size(); ++particle_index) {
    const std::uint32_t species_tag = state.particle_sidecar.species_tag[particle_index];
    parallel::DecompositionItem item;
    item.entity_id = state.particle_sidecar.particle_id[particle_index];
    item.kind = parallel::DecompositionEntityKind::kParticle;
    item.current_owner_rank = state.particle_sidecar.owning_rank.empty()
        ? world_rank
        : static_cast<int>(state.particle_sidecar.owning_rank[particle_index]);
    item.x_comov = state.particles.position_x_comoving[particle_index];
    item.y_comov = state.particles.position_y_comoving[particle_index];
    item.z_comov = state.particles.position_z_comoving[particle_index];
    item.active_target_count_recent = active_mask[particle_index] != 0U ? 1U : 0U;
    item.remote_tree_interactions_recent = 1U;
    item.memory_bytes = estimateParticleMemoryBytesForDecomposition(state, species_tag);
    double amr_patch_cost = 0.0;
    if (species_tag == static_cast<std::uint32_t>(core::ParticleSpecies::kGas) && !state.cells.patch_index.empty()) {
      std::vector<std::uint32_t> seen_patch_indices;
      for (const std::uint32_t cell_row : state.gas_cell_identity.rowsForParentParticleId(item.entity_id)) {
        if (cell_row >= state.cells.patch_index.size()) {
          continue;
        }
        const std::uint32_t patch_index = state.cells.patch_index[cell_row];
        if (patch_index >= patch_cell_count.size() ||
            std::find(seen_patch_indices.begin(), seen_patch_indices.end(), patch_index) != seen_patch_indices.end()) {
          continue;
        }
        seen_patch_indices.push_back(patch_index);
        amr_patch_cost += static_cast<double>(patch_cell_count[patch_index]) *
            (1.0 + static_cast<double>(std::max<std::int32_t>(state.patches.level[patch_index], 0)));
      }
    }
    item.work_components = parallel::DecompositionWorkComponents{
        .particle_count_cost = 1.0,
        .gas_cell_cost = species_tag == static_cast<std::uint32_t>(core::ParticleSpecies::kGas) ? 1.0 : 0.0,
        .tree_interaction_cost = 1.0,
        .pm_mesh_cost = 1.0,
        .amr_patch_cost = amr_patch_cost,
        .active_fraction_cost = static_cast<double>(item.active_target_count_recent),
        .memory_pressure_cost = static_cast<double>(item.memory_bytes),
        .gpu_occupancy_cost = 0.0,
        .generic_work_cost = 1.0,
        .has_explicit_components = true,
    };
    items.push_back(item);
  }

  for (std::size_t patch_index = 0; patch_index < state.patches.size(); ++patch_index) {
    if (state.patches.cell_count[patch_index] == 0U) {
      continue;
    }
    const std::uint32_t first_cell = state.patches.first_cell[patch_index];
    if (first_cell >= state.cells.size()) {
      continue;
    }
    parallel::DecompositionItem item;
    item.entity_id = state.patches.patch_id[patch_index];
    item.kind = parallel::DecompositionEntityKind::kAmrPatch;
    item.current_owner_rank = state.patches.owning_rank.empty()
        ? world_rank
        : static_cast<int>(state.patches.owning_rank[patch_index]);
    item.x_comov = state.cells.center_x_comoving[first_cell];
    item.y_comov = state.cells.center_y_comoving[first_cell];
    item.z_comov = state.cells.center_z_comoving[first_cell];
    item.memory_bytes = static_cast<std::uint64_t>(state.patches.cell_count[patch_index]) *
        static_cast<std::uint64_t>(sizeof(double) * 8U + sizeof(std::uint32_t) * 2U);
    item.work_components = parallel::DecompositionWorkComponents{
        .amr_patch_cost = static_cast<double>(state.patches.cell_count[patch_index]) *
            (1.0 + static_cast<double>(std::max<std::int32_t>(state.patches.level[patch_index], 0))),
        .memory_pressure_cost = static_cast<double>(item.memory_bytes),
        .generic_work_cost = static_cast<double>(state.patches.cell_count[patch_index]),
        .has_explicit_components = true,
    };
    items.push_back(item);
  }
  parallel::applyRuntimeDecompositionFeedback(items, parallel::DecompositionRuntimeMeasurements{}, makeWorkflowFeedbackCoefficients(config));
  return items;
}


void syncTimeBinsFromScheduler(
    const core::HierarchicalTimeBinScheduler& scheduler,
    core::SimulationState& state);

namespace migration_wire {

void appendBytes(std::vector<std::uint8_t>& out, const void* data, std::size_t bytes) {
  const auto* first = static_cast<const std::uint8_t*>(data);
  out.insert(out.end(), first, first + bytes);
}

template <typename T>
void appendPod(std::vector<std::uint8_t>& out, const T& value) {
  static_assert(std::is_trivially_copyable_v<T>);
  appendBytes(out, &value, sizeof(T));
}

template <typename T>
T readPod(std::span<const std::uint8_t> bytes, std::size_t& offset, std::string_view label) {
  static_assert(std::is_trivially_copyable_v<T>);
  if (offset + sizeof(T) > bytes.size()) {
    throw std::runtime_error("particle migration wire packet truncated while reading " + std::string(label));
  }
  T value{};
  std::memcpy(&value, bytes.data() + offset, sizeof(T));
  offset += sizeof(T);
  return value;
}

void appendString(std::vector<std::uint8_t>& out, std::string_view value) {
  const std::uint64_t size = static_cast<std::uint64_t>(value.size());
  appendPod(out, size);
  appendBytes(out, value.data(), value.size());
}

std::string readString(std::span<const std::uint8_t> bytes, std::size_t& offset, std::string_view label) {
  const std::uint64_t size = readPod<std::uint64_t>(bytes, offset, label);
  if (size > static_cast<std::uint64_t>(bytes.size() - offset)) {
    throw std::runtime_error("particle migration wire packet truncated while reading string " + std::string(label));
  }
  std::string value(reinterpret_cast<const char*>(bytes.data() + offset), static_cast<std::size_t>(size));
  offset += static_cast<std::size_t>(size);
  return value;
}

void appendBytePayload(std::vector<std::uint8_t>& out, std::span<const std::byte> payload) {
  const std::uint64_t size = static_cast<std::uint64_t>(payload.size());
  appendPod(out, size);
  if (!payload.empty()) {
    appendBytes(out, payload.data(), payload.size());
  }
}

std::vector<std::byte> readBytePayload(std::span<const std::uint8_t> bytes, std::size_t& offset, std::string_view label) {
  const std::uint64_t size = readPod<std::uint64_t>(bytes, offset, label);
  if (size > static_cast<std::uint64_t>(bytes.size() - offset)) {
    throw std::runtime_error("particle migration wire packet truncated while reading payload " + std::string(label));
  }
  std::vector<std::byte> payload(static_cast<std::size_t>(size));
  if (!payload.empty()) {
    std::memcpy(payload.data(), bytes.data() + offset, payload.size());
  }
  offset += static_cast<std::size_t>(size);
  return payload;
}

void appendModulePayload(std::vector<std::uint8_t>& out, const core::ModuleParticleSidecarPayload& payload) {
  appendString(out, payload.module_name);
  appendPod(out, payload.schema_version);
  appendPod(out, payload.row_stride_bytes);
  appendPod(out, payload.required_species_mask);
  appendPod(out, static_cast<std::uint32_t>(payload.requirement.kind));
  appendPod(out, payload.requirement.species_mask);
  appendPod(out, payload.requirement.particle_flags_mask);
  appendPod(out, payload.requirement.threshold_code);
  appendBytePayload(out, std::span<const std::byte>(payload.payload.data(), payload.payload.size()));
}

core::ModuleParticleSidecarPayload readModulePayload(std::span<const std::uint8_t> bytes, std::size_t& offset) {
  core::ModuleParticleSidecarPayload payload;
  payload.module_name = readString(bytes, offset, "module_name");
  payload.schema_version = readPod<std::uint32_t>(bytes, offset, "module_schema_version");
  payload.row_stride_bytes = readPod<std::uint32_t>(bytes, offset, "module_row_stride_bytes");
  payload.required_species_mask = readPod<std::uint32_t>(bytes, offset, "module_required_species_mask");
  payload.requirement.kind = static_cast<core::ModuleSidecarRequirementKind>(
      readPod<std::uint32_t>(bytes, offset, "module_requirement_kind"));
  payload.requirement.species_mask = readPod<std::uint32_t>(bytes, offset, "module_requirement_species_mask");
  payload.requirement.particle_flags_mask = readPod<std::uint32_t>(bytes, offset, "module_requirement_particle_flags_mask");
  payload.requirement.threshold_code = readPod<double>(bytes, offset, "module_requirement_threshold_code");
  payload.payload = readBytePayload(bytes, offset, "module_payload");
  return payload;
}

void appendMigrationRecord(std::vector<std::uint8_t>& out, const core::ParticleMigrationRecord& record) {
  appendPod(out, record.particle_id);
  appendPod(out, record.sfc_key);
  appendPod(out, record.species_tag);
  appendPod(out, record.particle_flags);
  appendPod(out, record.owning_rank);
  appendPod(out, record.position_x_comoving);
  appendPod(out, record.position_y_comoving);
  appendPod(out, record.position_z_comoving);
  appendPod(out, record.velocity_x_peculiar);
  appendPod(out, record.velocity_y_peculiar);
  appendPod(out, record.velocity_z_peculiar);
  appendPod(out, record.mass_code);
  appendPod(out, record.time_bin);
  appendPod(out, static_cast<std::uint8_t>(record.has_scheduler_fields ? 1U : 0U));
  appendPod(out, record.scheduler_fields);
  appendPod(out, record.last_drift_time_code);
  appendPod(out, record.last_drift_scale_factor);
  appendPod(out, static_cast<std::uint8_t>(record.has_gravity_softening_value ? 1U : 0U));
  appendPod(out, static_cast<std::uint8_t>(record.has_gravity_softening_override ? 1U : 0U));
  appendPod(out, record.gravity_softening_comoving);
  appendPod(out, static_cast<std::uint8_t>(record.has_gas_cell_fields ? 1U : 0U));
  appendPod(out, record.gas_cell_fields);
  appendPod(out, static_cast<std::uint8_t>(record.has_star_fields ? 1U : 0U));
  appendPod(out, record.star_fields);
  appendPod(out, static_cast<std::uint8_t>(record.has_black_hole_fields ? 1U : 0U));
  appendPod(out, record.black_hole_fields);
  appendPod(out, static_cast<std::uint8_t>(record.has_tracer_fields ? 1U : 0U));
  appendPod(out, record.tracer_fields);
  const std::uint64_t module_count = static_cast<std::uint64_t>(record.module_sidecar_payloads.size());
  appendPod(out, module_count);
  for (const auto& payload : record.module_sidecar_payloads) {
    appendModulePayload(out, payload);
  }
}

core::ParticleMigrationRecord readMigrationRecord(std::span<const std::uint8_t> bytes, std::size_t& offset) {
  core::ParticleMigrationRecord record;
  record.particle_id = readPod<std::uint64_t>(bytes, offset, "particle_id");
  record.sfc_key = readPod<std::uint64_t>(bytes, offset, "sfc_key");
  record.species_tag = readPod<std::uint32_t>(bytes, offset, "species_tag");
  record.particle_flags = readPod<std::uint32_t>(bytes, offset, "particle_flags");
  record.owning_rank = readPod<std::uint32_t>(bytes, offset, "owning_rank");
  record.position_x_comoving = readPod<double>(bytes, offset, "position_x_comoving");
  record.position_y_comoving = readPod<double>(bytes, offset, "position_y_comoving");
  record.position_z_comoving = readPod<double>(bytes, offset, "position_z_comoving");
  record.velocity_x_peculiar = readPod<double>(bytes, offset, "velocity_x_peculiar");
  record.velocity_y_peculiar = readPod<double>(bytes, offset, "velocity_y_peculiar");
  record.velocity_z_peculiar = readPod<double>(bytes, offset, "velocity_z_peculiar");
  record.mass_code = readPod<double>(bytes, offset, "mass_code");
  record.time_bin = readPod<std::uint8_t>(bytes, offset, "time_bin");
  record.has_scheduler_fields = readPod<std::uint8_t>(bytes, offset, "has_scheduler_fields") != 0U;
  record.scheduler_fields = readPod<core::SchedulerMigrationFields>(bytes, offset, "scheduler_fields");
  record.last_drift_time_code = readPod<double>(bytes, offset, "last_drift_time_code");
  record.last_drift_scale_factor = readPod<double>(bytes, offset, "last_drift_scale_factor");
  record.has_gravity_softening_value = readPod<std::uint8_t>(bytes, offset, "has_gravity_softening_value") != 0U;
  record.has_gravity_softening_override = readPod<std::uint8_t>(bytes, offset, "has_gravity_softening_override") != 0U;
  record.gravity_softening_comoving = readPod<double>(bytes, offset, "gravity_softening_comoving");
  record.has_gas_cell_fields = readPod<std::uint8_t>(bytes, offset, "has_gas_cell_fields") != 0U;
  record.gas_cell_fields = readPod<core::GasCellMigrationFields>(bytes, offset, "gas_cell_fields");
  record.has_star_fields = readPod<std::uint8_t>(bytes, offset, "has_star_fields") != 0U;
  record.star_fields = readPod<core::StarParticleMigrationFields>(bytes, offset, "star_fields");
  record.has_black_hole_fields = readPod<std::uint8_t>(bytes, offset, "has_black_hole_fields") != 0U;
  record.black_hole_fields = readPod<core::BlackHoleParticleMigrationFields>(bytes, offset, "black_hole_fields");
  record.has_tracer_fields = readPod<std::uint8_t>(bytes, offset, "has_tracer_fields") != 0U;
  record.tracer_fields = readPod<core::TracerParticleMigrationFields>(bytes, offset, "tracer_fields");
  const std::uint64_t module_count = readPod<std::uint64_t>(bytes, offset, "module_count");
  if (module_count > 1'000'000ULL) {
    throw std::runtime_error("particle migration wire packet has unreasonable module payload count");
  }
  record.module_sidecar_payloads.reserve(static_cast<std::size_t>(module_count));
  for (std::uint64_t i = 0; i < module_count; ++i) {
    record.module_sidecar_payloads.push_back(readModulePayload(bytes, offset));
  }
  return record;
}

std::vector<std::uint8_t> serializeMigrationRecords(std::span<const core::ParticleMigrationRecord> records) {
  std::vector<std::uint8_t> bytes;
  appendPod(bytes, static_cast<std::uint64_t>(records.size()));
  for (const core::ParticleMigrationRecord& record : records) {
    appendMigrationRecord(bytes, record);
  }
  return bytes;
}

std::vector<core::ParticleMigrationRecord> deserializeMigrationRecords(std::span<const std::uint8_t> bytes) {
  std::size_t offset = 0;
  const std::uint64_t count = readPod<std::uint64_t>(bytes, offset, "migration_record_count");
  if (count > 100'000'000ULL) {
    throw std::runtime_error("particle migration wire packet has unreasonable record count");
  }
  std::vector<core::ParticleMigrationRecord> records;
  records.reserve(static_cast<std::size_t>(count));
  for (std::uint64_t i = 0; i < count; ++i) {
    records.push_back(readMigrationRecord(bytes, offset));
  }
  if (offset != bytes.size()) {
    throw std::runtime_error("particle migration wire packet has trailing bytes");
  }
  return records;
}

}  // namespace migration_wire

std::vector<core::ParticleMigrationRecord> exchangeRuntimeParticleMigrationRecords(
    const parallel::MpiContext& mpi_context,
    const std::vector<std::vector<core::ParticleMigrationRecord>>& records_by_rank) {
  if (records_by_rank.size() != static_cast<std::size_t>(mpi_context.worldSize())) {
    throw std::invalid_argument("particle migration exchange records_by_rank size must match world size");
  }
  if (mpi_context.worldSize() == 1) {
    return records_by_rank.empty() ? std::vector<core::ParticleMigrationRecord>{} : records_by_rank[0];
  }
  if (!mpi_context.isEnabled()) {
    throw std::runtime_error("runtime particle migration execution requires MPI for world_size > 1");
  }
#if COSMOSIM_ENABLE_MPI
  std::vector<std::vector<std::uint8_t>> send_payloads(records_by_rank.size());
  std::vector<int> send_counts(records_by_rank.size(), 0);
  for (std::size_t rank = 0; rank < records_by_rank.size(); ++rank) {
    send_payloads[rank] = migration_wire::serializeMigrationRecords(records_by_rank[rank]);
    if (send_payloads[rank].size() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
      throw std::overflow_error("particle migration payload exceeds MPI int byte-count limit");
    }
    send_counts[rank] = static_cast<int>(send_payloads[rank].size());
  }
  std::vector<int> recv_counts(send_counts.size(), 0);
  MPI_Alltoall(send_counts.data(), 1, MPI_INT, recv_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);

  std::vector<int> send_offsets(send_counts.size(), 0);
  std::vector<int> recv_offsets(recv_counts.size(), 0);
  int total_send = 0;
  int total_recv = 0;
  for (std::size_t rank = 0; rank < send_counts.size(); ++rank) {
    send_offsets[rank] = total_send;
    recv_offsets[rank] = total_recv;
    total_send += send_counts[rank];
    total_recv += recv_counts[rank];
  }
  std::vector<std::uint8_t> send_bytes(static_cast<std::size_t>(total_send));
  for (std::size_t rank = 0; rank < send_payloads.size(); ++rank) {
    if (!send_payloads[rank].empty()) {
      std::memcpy(send_bytes.data() + send_offsets[rank], send_payloads[rank].data(), send_payloads[rank].size());
    }
  }
  std::vector<std::uint8_t> recv_bytes(static_cast<std::size_t>(total_recv));
  MPI_Alltoallv(
      send_bytes.data(),
      send_counts.data(),
      send_offsets.data(),
      MPI_BYTE,
      recv_bytes.data(),
      recv_counts.data(),
      recv_offsets.data(),
      MPI_BYTE,
      MPI_COMM_WORLD);

  std::vector<core::ParticleMigrationRecord> inbound;
  for (std::size_t rank = 0; rank < recv_counts.size(); ++rank) {
    const auto offset = static_cast<std::size_t>(recv_offsets[rank]);
    const auto count = static_cast<std::size_t>(recv_counts[rank]);
    auto decoded = migration_wire::deserializeMigrationRecords(
        std::span<const std::uint8_t>(recv_bytes.data() + offset, count));
    inbound.insert(inbound.end(), std::make_move_iterator(decoded.begin()), std::make_move_iterator(decoded.end()));
  }
  return inbound;
#else
  throw std::runtime_error("runtime particle migration execution requires an MPI-enabled build");
#endif
}

void requireGlobalOwnedParticlePartitionIdentity(
    const core::SimulationState& state,
    const parallel::MpiContext& mpi_context,
    std::span<const std::uint64_t> expected_global_particle_ids,
    std::string_view caller) {
  std::vector<std::uint64_t> local_owned_ids;
  local_owned_ids.reserve(state.particles.size());
  const std::uint32_t world_rank = static_cast<std::uint32_t>(mpi_context.worldRank());
  for (std::size_t particle_index = 0; particle_index < state.particles.size(); ++particle_index) {
    if (state.particle_sidecar.owning_rank[particle_index] == world_rank) {
      local_owned_ids.push_back(state.particle_sidecar.particle_id[particle_index]);
    }
  }
  const parallel::ExactOwnershipPartitionReport report =
      parallel::validateExactGlobalOwnershipPartition(mpi_context, local_owned_ids, expected_global_particle_ids);
  if (!report.valid()) {
    throw std::runtime_error(std::string(caller) +
        ": exact distributed authoritative ownership table has duplicate=" +
        std::to_string(report.duplicate_particle_ids.size()) +
        ", missing=" + std::to_string(report.missing_expected_particle_ids.size()) +
        ", extra=" + std::to_string(report.extra_particle_ids.size()) + " particle IDs");
  }
}

void recordRuntimeRebalanceDecision(
    core::ProfilerSession* profiler,
    const parallel::RuntimeRebalancePlan& rebalance,
    std::uint64_t step_index) {
  if (profiler == nullptr) {
    return;
  }
  profiler->recordEvent(core::RuntimeEvent{
      .event_kind = "parallel.decomposition.runtime_rebalance",
      .severity = rebalance.should_rebalance ? core::RuntimeEventSeverity::kWarning : core::RuntimeEventSeverity::kInfo,
      .subsystem = "parallel.domain_decomposition",
      .step_index = step_index,
      .message = rebalance.should_rebalance
          ? "runtime work-weighted decomposition produced migration intent"
          : "runtime work-weighted decomposition remained below rebalance threshold",
      .payload = {{"reason", rebalance.reason},
                  {"should_rebalance", rebalance.should_rebalance ? "true" : "false"},
                  {"particle_migration_count", std::to_string(rebalance.particle_migrations.size())},
                  {"amr_patch_ownership_update_count", std::to_string(rebalance.amr_patch_ownership_updates.size())},
                  {"migrated_load_fraction", std::to_string(rebalance.migrated_load_fraction)},
                  {"current_weighted_imbalance_ratio", std::to_string(rebalance.current_metrics.weighted_imbalance_ratio)},
                  {"target_weighted_imbalance_ratio", std::to_string(rebalance.target_decomposition.metrics.weighted_imbalance_ratio)},
                  {"current_memory_imbalance_ratio", std::to_string(rebalance.current_metrics.memory_imbalance_ratio)},
                  {"target_memory_imbalance_ratio", std::to_string(rebalance.target_decomposition.metrics.memory_imbalance_ratio)}}});
}

[[nodiscard]] std::vector<parallel::AmrPatchPayloadRecord> buildLocalAmrPatchPayloadRecords(
    const core::SimulationState& state,
    int world_rank) {
  std::vector<parallel::AmrPatchPayloadRecord> records;
  if (world_rank < 0 || state.patches.size() == 0) {
    return records;
  }
  records.reserve(state.patches.size());
  for (std::size_t patch_index = 0; patch_index < state.patches.size(); ++patch_index) {
    if (state.patches.owning_rank[patch_index] != static_cast<std::uint32_t>(world_rank)) {
      continue;
    }
    const std::uint32_t first_cell = state.patches.first_cell[patch_index];
    const std::uint32_t cell_count = state.patches.cell_count[patch_index];
    if (cell_count == 0U) {
      continue;
    }
    if (static_cast<std::uint64_t>(first_cell) + static_cast<std::uint64_t>(cell_count) > state.cells.size()) {
      throw std::runtime_error("AMR patch payload build found a patch range outside CellSoa");
    }
    parallel::AmrPatchPayloadRecord record;
    record.patch_id = state.patches.patch_id[patch_index];
    record.owner_rank = world_rank;
    record.level = static_cast<std::uint32_t>(std::max<std::int32_t>(state.patches.level[patch_index], 0));
    record.first_cell = first_cell;
    record.cell_count = cell_count;
    for (std::uint32_t offset = 0; offset < cell_count; ++offset) {
      const std::uint32_t cell_index = first_cell + offset;
      record.cell_mass_sum_code += state.cells.mass_code[cell_index];
      if (cell_index < state.gas_cells.internal_energy_code.size()) {
        record.gas_internal_energy_sum_code += state.gas_cells.internal_energy_code[cell_index];
      }
    }
    parallel::validateAmrPatchPayloadRecord(record);
    records.push_back(record);
  }
  return records;
}


[[nodiscard]] std::vector<parallel::AmrPatchCellPayloadRecord> buildLocalAmrPatchCellPayloadRecords(
    const core::SimulationState& state,
    int world_rank) {
  std::vector<parallel::AmrPatchCellPayloadRecord> records;
  if (world_rank < 0 || state.patches.size() == 0) {
    return records;
  }
  std::size_t total_owned_cells = 0;
  for (std::size_t patch_index = 0; patch_index < state.patches.size(); ++patch_index) {
    if (state.patches.owning_rank[patch_index] == static_cast<std::uint32_t>(world_rank)) {
      total_owned_cells += state.patches.cell_count[patch_index];
    }
  }
  records.reserve(total_owned_cells);
  for (std::size_t patch_index = 0; patch_index < state.patches.size(); ++patch_index) {
    if (state.patches.owning_rank[patch_index] != static_cast<std::uint32_t>(world_rank)) {
      continue;
    }
    const std::uint32_t first_cell = state.patches.first_cell[patch_index];
    const std::uint32_t cell_count = state.patches.cell_count[patch_index];
    if (static_cast<std::uint64_t>(first_cell) + static_cast<std::uint64_t>(cell_count) > state.cells.size()) {
      throw std::runtime_error("AMR patch cell payload build found a patch range outside CellSoa");
    }
    for (std::uint32_t offset = 0; offset < cell_count; ++offset) {
      const std::uint32_t cell_index = first_cell + offset;
      parallel::AmrPatchCellPayloadRecord record;
      record.patch_id = state.patches.patch_id[patch_index];
      record.owner_rank = world_rank;
      record.local_cell_offset = offset;
      record.patch_index = static_cast<std::uint32_t>(patch_index);
      record.center_x_comoving = state.cells.center_x_comoving[cell_index];
      record.center_y_comoving = state.cells.center_y_comoving[cell_index];
      record.center_z_comoving = state.cells.center_z_comoving[cell_index];
      record.mass_code = state.cells.mass_code[cell_index];
      record.time_bin = state.cells.time_bin[cell_index];
      const core::GasCellIdentityRecord& identity = gasCellIdentityRecordForLocalRow(
          state,
          cell_index,
          "buildLocalAmrPatchCellPayloadRecords");
      record.gas_cell_id = identity.gas_cell_id;
      record.parent_particle_id = identity.parent_particle_id.value_or(0U);
      record.density_code = cell_index < state.gas_cells.density_code.size() ? state.gas_cells.density_code[cell_index] : 0.0;
      record.pressure_code = cell_index < state.gas_cells.pressure_code.size() ? state.gas_cells.pressure_code[cell_index] : 0.0;
      record.internal_energy_code = cell_index < state.gas_cells.internal_energy_code.size()
          ? state.gas_cells.internal_energy_code[cell_index]
          : 0.0;
      record.temperature_code = cell_index < state.gas_cells.temperature_code.size()
          ? state.gas_cells.temperature_code[cell_index]
          : 0.0;
      record.sound_speed_code = cell_index < state.gas_cells.sound_speed_code.size()
          ? state.gas_cells.sound_speed_code[cell_index]
          : 0.0;
      parallel::validateAmrPatchCellPayloadRecord(record);
      records.push_back(record);
    }
  }
  return records;
}

void exchangeAndValidateAmrPatchPayloads(
    core::SimulationState& state,
    const parallel::MpiContext& mpi_context,
    int world_rank,
    std::uint64_t step_index,
    core::ProfilerSession* profiler) {
  const std::vector<parallel::AmrPatchPayloadRecord> local_records = buildLocalAmrPatchPayloadRecords(state, world_rank);
  const std::vector<parallel::AmrPatchCellPayloadRecord> local_cell_records =
      buildLocalAmrPatchCellPayloadRecords(state, world_rank);
  const std::vector<parallel::AmrPatchPayloadRecord> global_records = parallel::executeBlockingAmrPatchPayloadExchange(
      mpi_context,
      local_records,
      step_index);
  const std::vector<parallel::AmrPatchCellPayloadRecord> global_cell_records =
      parallel::executeBlockingAmrPatchCellPayloadExchange(mpi_context, local_cell_records, step_index);

  std::unordered_map<std::uint64_t, int> owner_by_patch_id;
  std::unordered_map<std::uint64_t, std::uint32_t> expected_cell_count_by_patch_id;
  std::unordered_map<std::uint64_t, std::uint32_t> observed_cell_count_by_patch_id;
  owner_by_patch_id.reserve(global_records.size());
  expected_cell_count_by_patch_id.reserve(global_records.size());
  observed_cell_count_by_patch_id.reserve(global_records.size());
  for (const parallel::AmrPatchPayloadRecord& record : global_records) {
    const auto [it, inserted] = owner_by_patch_id.emplace(record.patch_id, record.owner_rank);
    if (!inserted && it->second != record.owner_rank) {
      throw std::runtime_error("AMR patch payload exchange detected duplicate authoritative patch ownership");
    }
    expected_cell_count_by_patch_id.emplace(record.patch_id, record.cell_count);
  }
  for (const parallel::AmrPatchCellPayloadRecord& record : global_cell_records) {
    const auto owner_it = owner_by_patch_id.find(record.patch_id);
    if (owner_it == owner_by_patch_id.end()) {
      throw std::runtime_error("AMR patch cell payload exchange found a cell for an unknown patch");
    }
    if (owner_it->second != record.owner_rank) {
      throw std::runtime_error("AMR patch cell payload owner does not match authoritative patch owner");
    }
    ++observed_cell_count_by_patch_id[record.patch_id];
  }
  for (const auto& [patch_id, expected_count] : expected_cell_count_by_patch_id) {
    const std::uint32_t observed_count = observed_cell_count_by_patch_id[patch_id];
    if (observed_count != expected_count) {
      throw std::runtime_error("AMR patch cell payload exchange coverage mismatch for patch_id=" +
                               std::to_string(patch_id));
    }
  }


  std::unordered_map<std::uint64_t, std::size_t> local_patch_index_by_id;
  local_patch_index_by_id.reserve(state.patches.size());
  for (std::size_t patch_index = 0; patch_index < state.patches.size(); ++patch_index) {
    if (state.patches.owning_rank[patch_index] == static_cast<std::uint32_t>(world_rank)) {
      local_patch_index_by_id.emplace(state.patches.patch_id[patch_index], patch_index);
    }
  }
  std::size_t applied_owned_patch_cells = 0;
  for (const parallel::AmrPatchCellPayloadRecord& record : global_cell_records) {
    if (record.owner_rank != world_rank) {
      continue;
    }
    const auto patch_it = local_patch_index_by_id.find(record.patch_id);
    if (patch_it == local_patch_index_by_id.end()) {
      throw std::runtime_error("AMR patch cell payload for this rank has no local authoritative patch");
    }
    const std::size_t patch_index = patch_it->second;
    if (record.local_cell_offset >= state.patches.cell_count[patch_index]) {
      throw std::runtime_error("AMR patch cell payload offset exceeds local patch cell count");
    }
    const std::uint32_t cell_index = state.patches.first_cell[patch_index] + record.local_cell_offset;
    if (cell_index >= state.cells.size()) {
      throw std::runtime_error("AMR patch cell payload target cell is outside local CellSoa");
    }
    state.cells.center_x_comoving[cell_index] = record.center_x_comoving;
    state.cells.center_y_comoving[cell_index] = record.center_y_comoving;
    state.cells.center_z_comoving[cell_index] = record.center_z_comoving;
    state.cells.mass_code[cell_index] = record.mass_code;
    state.cells.time_bin[cell_index] = record.time_bin;
    if (cell_index < state.gas_cells.density_code.size()) {
      state.gas_cells.gas_cell_id[cell_index] = record.gas_cell_id;
      state.gas_cells.parent_particle_id[cell_index] = record.parent_particle_id;
      state.gas_cells.density_code[cell_index] = record.density_code;
      state.gas_cells.pressure_code[cell_index] = record.pressure_code;
      state.gas_cells.internal_energy_code[cell_index] = record.internal_energy_code;
      state.gas_cells.temperature_code[cell_index] = record.temperature_code;
      state.gas_cells.sound_speed_code[cell_index] = record.sound_speed_code;
    }
    ++applied_owned_patch_cells;
  }
  if (profiler != nullptr) {
    profiler->recordEvent(core::RuntimeEvent{
        .event_kind = "amr.patch_payload_exchange",
        .severity = core::RuntimeEventSeverity::kInfo,
        .subsystem = "amr.patch_exchange",
        .step_index = step_index,
        .message = "blocking AMR patch payload exchange validated authoritative patch ownership and cell summaries",
        .payload = {{"local_patch_payloads", std::to_string(local_records.size())},
                    {"global_patch_payloads", std::to_string(global_records.size())},
                    {"global_patch_cell_payloads", std::to_string(global_cell_records.size())},
                    {"applied_owned_patch_cells", std::to_string(applied_owned_patch_cells)},
                    {"payload_bytes", std::to_string(global_records.size() * sizeof(parallel::AmrPatchPayloadRecord) +
                                                       global_cell_records.size() * sizeof(parallel::AmrPatchCellPayloadRecord))}}});
  }
}

void applyMeasuredRuntimeRebalancePlan(
    core::SimulationState& state,
    core::HierarchicalTimeBinScheduler& scheduler,
    const core::SimulationConfig& config,
    const parallel::MpiContext& mpi_context,
    int world_rank,
    const parallel::DecompositionRuntimeMeasurements& measurements,
    std::span<const std::uint32_t> active_particle_indices,
    std::span<const std::uint64_t> expected_global_particle_ids,
    core::ProfilerSession* profiler,
    std::uint64_t step_index) {
  if (!config.parallel.decomposition_runtime_rebalance_enabled || mpi_context.worldSize() <= 1) {
    return;
  }
  if (world_rank < 0 || world_rank >= mpi_context.worldSize()) {
    throw std::invalid_argument("runtime rebalance world_rank is outside MPI world");
  }
  auto local_items = buildRuntimeDecompositionItems(state, config, world_rank, active_particle_indices);
  parallel::applyRuntimeDecompositionFeedback(local_items, measurements, makeWorkflowFeedbackCoefficients(config));
  const std::vector<parallel::DecompositionItem> global_items =
      parallel::gatherDecompositionItemsAcrossRanks(mpi_context, local_items);
  parallel::RuntimeRebalanceConfig rebalance_config{
      .world_size = mpi_context.worldSize(),
      .imbalance_trigger_ratio = config.parallel.decomposition_rebalance_imbalance_trigger,
      .memory_trigger_ratio = config.parallel.decomposition_rebalance_memory_trigger,
      .max_migrated_load_fraction = config.parallel.decomposition_rebalance_max_migrated_load_fraction,
      .allow_particle_migration = true,
      .allow_amr_patch_reassignment = true,
  };
  const auto rebalance = parallel::buildRuntimeRebalancePlan(
      global_items,
      makeWorkflowDecompositionConfig(config, mpi_context.worldSize()),
      rebalance_config);
  if (!rebalance.should_rebalance) {
    exchangeAndValidateAmrPatchPayloads(state, mpi_context, world_rank, step_index, profiler);
    recordRuntimeRebalanceDecision(profiler, rebalance, step_index);
    return;
  }

  std::unordered_map<std::uint64_t, std::uint32_t> local_index_by_particle_id;
  local_index_by_particle_id.reserve(state.particles.size());
  for (std::size_t particle_index = 0; particle_index < state.particles.size(); ++particle_index) {
    local_index_by_particle_id.emplace(
        state.particle_sidecar.particle_id[particle_index],
        static_cast<std::uint32_t>(particle_index));
  }

  std::unordered_map<std::uint32_t, int> outbound_target_by_local_index;
  outbound_target_by_local_index.reserve(rebalance.particle_migrations.size());
  const auto add_particle_migration = [&](std::uint32_t local_index, int new_owner_rank, std::string_view source_label) {
    if (local_index >= state.particles.size()) {
      throw std::out_of_range("runtime rebalance attempted to migrate a particle index outside SimulationState");
    }
    if (new_owner_rank < 0 || new_owner_rank >= mpi_context.worldSize()) {
      throw std::invalid_argument("runtime rebalance produced particle migration target outside MPI world");
    }
    if (new_owner_rank == world_rank) {
      return;
    }
    if (state.particle_sidecar.owning_rank[local_index] != static_cast<std::uint32_t>(world_rank)) {
      throw std::runtime_error("runtime rebalance attempted to migrate a non-authoritative local particle via " + std::string(source_label));
    }
    const auto [it, inserted] = outbound_target_by_local_index.emplace(local_index, new_owner_rank);
    if (!inserted && it->second != new_owner_rank) {
      throw std::runtime_error("runtime rebalance produced conflicting destinations for one particle");
    }
  };

  for (const auto& intent : rebalance.particle_migrations) {
    if (intent.old_owner_rank != world_rank || intent.new_owner_rank == world_rank) {
      continue;
    }
    const auto found = local_index_by_particle_id.find(intent.particle_id);
    if (found == local_index_by_particle_id.end()) {
      throw std::runtime_error("runtime rebalance migration intent references a particle ID missing on the owning rank");
    }
    add_particle_migration(found->second, intent.new_owner_rank, "particle_sfc_intent");
  }

  for (const auto& update : rebalance.amr_patch_ownership_updates) {
    if (update.new_owner_rank < 0 || update.new_owner_rank >= mpi_context.worldSize()) {
      throw std::runtime_error("runtime rebalance produced AMR patch owner outside MPI world");
    }
    for (std::size_t patch_index = 0; patch_index < state.patches.size(); ++patch_index) {
      if (state.patches.patch_id[patch_index] != update.patch_id) {
        continue;
      }
      if (state.patches.owning_rank[patch_index] == static_cast<std::uint32_t>(world_rank) &&
          update.new_owner_rank != world_rank) {
        const std::uint32_t first_cell = state.patches.first_cell[patch_index];
        const std::uint32_t cell_count = state.patches.cell_count[patch_index];
        if (static_cast<std::uint64_t>(first_cell) + static_cast<std::uint64_t>(cell_count) > state.cells.size()) {
          throw std::runtime_error("runtime rebalance AMR patch cell range is outside CellSoa");
        }
        for (std::uint32_t cell_offset = 0; cell_offset < cell_count; ++cell_offset) {
          const std::uint32_t cell_index = first_cell + cell_offset;
          const auto gas_particle_index = parentParticleRowForGasCellRow(
              state,
              cell_index,
              local_index_by_particle_id,
              "runtime rebalance AMR patch ownership update");
          if (gas_particle_index.has_value() &&
              state.particle_sidecar.owning_rank[*gas_particle_index] == static_cast<std::uint32_t>(world_rank)) {
            add_particle_migration(*gas_particle_index, update.new_owner_rank, "amr_patch_ownership_update");
          }
        }
      }
      state.patches.owning_rank[patch_index] = static_cast<std::uint32_t>(update.new_owner_rank);
    }
  }

  std::vector<std::uint32_t> outbound_local_indices;
  outbound_local_indices.reserve(outbound_target_by_local_index.size());
  for (const auto& [local_index, target_rank] : outbound_target_by_local_index) {
    (void)target_rank;
    outbound_local_indices.push_back(local_index);
  }
  std::sort(outbound_local_indices.begin(), outbound_local_indices.end());

  std::vector<std::uint32_t> preserved_indices;
  preserved_indices.reserve(state.particles.size() - outbound_local_indices.size());
  for (std::uint32_t particle_index = 0; particle_index < state.particles.size(); ++particle_index) {
    if (!std::binary_search(outbound_local_indices.begin(), outbound_local_indices.end(), particle_index)) {
      preserved_indices.push_back(particle_index);
    }
  }
  std::vector<core::ParticleMigrationRecord> scheduler_records =
      state.packParticleMigrationRecords(preserved_indices, scheduler);

  std::vector<std::vector<core::ParticleMigrationRecord>> outbound_records_by_rank(
      static_cast<std::size_t>(mpi_context.worldSize()));
  for (const std::uint32_t local_index : outbound_local_indices) {
    const int target_rank = outbound_target_by_local_index.at(local_index);
    auto records = state.packParticleMigrationRecords(std::span<const std::uint32_t>(&local_index, 1), scheduler);
    if (records.size() != 1U) {
      throw std::runtime_error("runtime rebalance particle migration packing returned an unexpected record count");
    }
    records[0].owning_rank = static_cast<std::uint32_t>(target_rank);
    outbound_records_by_rank[static_cast<std::size_t>(target_rank)].push_back(std::move(records[0]));
  }

  std::vector<core::ParticleMigrationRecord> inbound_records =
      exchangeRuntimeParticleMigrationRecords(mpi_context, outbound_records_by_rank);
  for (const core::ParticleMigrationRecord& record : inbound_records) {
    if (record.owning_rank != static_cast<std::uint32_t>(world_rank)) {
      throw std::runtime_error("runtime particle migration exchange delivered a record to the wrong destination rank");
    }
  }

  if (!outbound_local_indices.empty() || !inbound_records.empty()) {
    scheduler_records.insert(scheduler_records.end(), inbound_records.begin(), inbound_records.end());
    core::ParticleMigrationCommit commit;
    commit.world_rank = world_rank;
    commit.outbound_local_indices = outbound_local_indices;
    commit.inbound_records = inbound_records;
    state.commitParticleMigration(commit);

    std::vector<std::uint64_t> destination_particle_ids(state.particle_sidecar.particle_id.begin(),
                                                        state.particle_sidecar.particle_id.end());
    core::rebuildSchedulerFromParticleMigrationRecords(scheduler, scheduler_records, destination_particle_ids);
    syncTimeBinsFromScheduler(scheduler, state);
    requireGlobalOwnedParticlePartitionIdentity(
        state, mpi_context, expected_global_particle_ids, "runtime rebalance particle migration commit");
  }

  exchangeAndValidateAmrPatchPayloads(state, mpi_context, world_rank, step_index, profiler);

  recordRuntimeRebalanceDecision(profiler, rebalance, step_index);
  if (profiler != nullptr && (!outbound_local_indices.empty() || !inbound_records.empty())) {
    profiler->recordEvent(core::RuntimeEvent{
        .event_kind = "parallel.decomposition.runtime_migration_commit",
        .severity = core::RuntimeEventSeverity::kWarning,
        .subsystem = "parallel.domain_decomposition",
        .step_index = step_index,
        .message = "runtime rebalance committed authoritative particle migration at a safe scheduler boundary",
        .payload = {
            {"outbound_particle_count", std::to_string(outbound_local_indices.size())},
            {"inbound_particle_count", std::to_string(inbound_records.size())},
            {"current_weighted_imbalance_ratio", std::to_string(rebalance.current_metrics.weighted_imbalance_ratio)},
            {"target_weighted_imbalance_ratio", std::to_string(rebalance.target_decomposition.metrics.weighted_imbalance_ratio)},
            {"current_memory_imbalance_ratio", std::to_string(rebalance.current_metrics.memory_imbalance_ratio)},
            {"target_memory_imbalance_ratio", std::to_string(rebalance.target_decomposition.metrics.memory_imbalance_ratio)},
        },
    });
  }
}

void applyInitialGravityAwareDecomposition(
    core::SimulationState& state,
    const core::SimulationConfig& config,
    int world_size,
    int world_rank,
    core::ProfilerSession* profiler) {
  if (world_size <= 1) {
    return;
  }
  constexpr std::size_t k_density_grid = 16;
  const std::size_t grid_cells = k_density_grid * k_density_grid * k_density_grid;
  std::vector<std::uint32_t> occupancy(grid_cells, 0U);
  std::vector<std::uint32_t> active_occupancy(grid_cells, 0U);
  std::vector<std::uint32_t> gas_occupancy(grid_cells, 0U);
  std::vector<std::uint32_t> pm_x_occupancy(static_cast<std::size_t>(std::max(config.numerics.treepm_pm_grid_nx, 1)), 0U);
  const auto wrap = [](double x, double box) {
    if (box <= 0.0) {
      return x;
    }
    double wrapped = std::fmod(x, box);
    if (wrapped < 0.0) {
      wrapped += box;
    }
    return wrapped;
  };
  const auto density_cell_index = [&](double x, double y, double z) {
    const double box_x = config.cosmology.box_size_x_mpc_comoving;
    const double box_y = config.cosmology.box_size_y_mpc_comoving;
    const double box_z = config.cosmology.box_size_z_mpc_comoving;
    const std::size_t ix = (box_x > 0.0)
        ? std::min<std::size_t>(
              k_density_grid - 1U,
              static_cast<std::size_t>((wrap(x, box_x) / box_x) * static_cast<double>(k_density_grid)))
        : 0U;
    const std::size_t iy = (box_y > 0.0)
        ? std::min<std::size_t>(
              k_density_grid - 1U,
              static_cast<std::size_t>((wrap(y, box_y) / box_y) * static_cast<double>(k_density_grid)))
        : 0U;
    const std::size_t iz = (box_z > 0.0)
        ? std::min<std::size_t>(
              k_density_grid - 1U,
              static_cast<std::size_t>((wrap(z, box_z) / box_z) * static_cast<double>(k_density_grid)))
        : 0U;
    return (ix * k_density_grid + iy) * k_density_grid + iz;
  };
  const auto pm_x_index = [&](double x) {
    const double box_x = config.cosmology.box_size_x_mpc_comoving;
    const std::size_t nx = pm_x_occupancy.size();
    if (box_x <= 0.0 || nx == 0) {
      return std::size_t{0};
    }
    const double scaled = (wrap(x, box_x) / box_x) * static_cast<double>(nx);
    return std::min<std::size_t>(nx - 1U, static_cast<std::size_t>(scaled));
  };

  for (std::size_t particle_index = 0; particle_index < state.particles.size(); ++particle_index) {
    const std::size_t cell = density_cell_index(
        state.particles.position_x_comoving[particle_index],
        state.particles.position_y_comoving[particle_index],
        state.particles.position_z_comoving[particle_index]);
    ++occupancy[cell];
    if (!state.particles.time_bin.empty() && state.particles.time_bin[particle_index] == 0U) {
      ++active_occupancy[cell];
    }
    if (state.particle_sidecar.species_tag[particle_index] == static_cast<std::uint32_t>(core::ParticleSpecies::kGas)) {
      ++gas_occupancy[cell];
    }
    ++pm_x_occupancy[pm_x_index(state.particles.position_x_comoving[particle_index])];
  }

  std::vector<std::uint32_t> patch_cell_count(state.patches.size(), 0U);
  for (std::size_t cell_index = 0; cell_index < state.cells.size(); ++cell_index) {
    const std::uint32_t patch_index = state.cells.patch_index[cell_index];
    if (patch_index < patch_cell_count.size()) {
      ++patch_cell_count[patch_index];
    }
  }

  const auto particle_memory_bytes = [&](std::uint32_t species_tag) {
    std::uint64_t bytes = sizeof(double) * 7U + sizeof(std::uint64_t) * 2U + sizeof(std::uint32_t) * 3U;
    if (!state.particle_sidecar.gravity_softening_comoving.empty()) {
      bytes += sizeof(double);
    }
    if (!state.particle_sidecar.has_gravity_softening_override.empty()) {
      bytes += sizeof(std::uint8_t);
    }
    if (species_tag == static_cast<std::uint32_t>(core::ParticleSpecies::kGas)) {
      bytes += sizeof(double) * 8U + sizeof(std::uint64_t) * 2U + sizeof(std::uint32_t);
    } else if (species_tag == static_cast<std::uint32_t>(core::ParticleSpecies::kStar)) {
      bytes += sizeof(std::uint32_t) + sizeof(double) * 13U;
    } else if (species_tag == static_cast<std::uint32_t>(core::ParticleSpecies::kBlackHole)) {
      bytes += sizeof(std::uint32_t) * 2U + sizeof(double) * 8U;
    } else if (species_tag == static_cast<std::uint32_t>(core::ParticleSpecies::kTracer)) {
      bytes += sizeof(std::uint64_t) * 2U + sizeof(std::uint32_t) * 2U + sizeof(double) * 3U;
    }
    return bytes;
  };

  std::vector<parallel::DecompositionItem> items;
  items.reserve(state.particles.size() + state.patches.size());
  constexpr std::uint32_t k_invalid_patch_index = std::numeric_limits<std::uint32_t>::max();
  std::vector<std::uint32_t> patch_index_by_item;
  patch_index_by_item.reserve(state.particles.size() + state.patches.size());
  for (std::size_t particle_index = 0; particle_index < state.particles.size(); ++particle_index) {
    parallel::DecompositionItem item;
    item.entity_id = state.particle_sidecar.particle_id[particle_index];
    item.kind = parallel::DecompositionEntityKind::kParticle;
    item.current_owner_rank = static_cast<int>(state.particle_sidecar.owning_rank[particle_index]);
    item.x_comov = state.particles.position_x_comoving[particle_index];
    item.y_comov = state.particles.position_y_comoving[particle_index];
    item.z_comov = state.particles.position_z_comoving[particle_index];
    const std::size_t density_cell = density_cell_index(item.x_comov, item.y_comov, item.z_comov);
    const std::uint32_t local_density = occupancy[density_cell];
    const std::uint32_t local_active = active_occupancy[density_cell];
    const std::uint32_t local_gas = gas_occupancy[density_cell];
    const std::uint32_t pm_load = pm_x_occupancy[pm_x_index(item.x_comov)];
    const std::uint32_t species_tag = state.particle_sidecar.species_tag[particle_index];
    double amr_patch_cost = 0.0;
    if (species_tag == static_cast<std::uint32_t>(core::ParticleSpecies::kGas) && !state.cells.patch_index.empty()) {
      std::vector<std::uint32_t> seen_patch_indices;
      for (const std::uint32_t cell_row : state.gas_cell_identity.rowsForParentParticleId(item.entity_id)) {
        if (cell_row >= state.cells.patch_index.size()) {
          continue;
        }
        const std::uint32_t patch_index = state.cells.patch_index[cell_row];
        if (patch_index >= patch_cell_count.size() ||
            std::find(seen_patch_indices.begin(), seen_patch_indices.end(), patch_index) != seen_patch_indices.end()) {
          continue;
        }
        seen_patch_indices.push_back(patch_index);
        amr_patch_cost += static_cast<double>(patch_cell_count[patch_index]);
      }
    }
    const double local_density_d = static_cast<double>(std::max<std::uint32_t>(local_density, 1U));
    item.active_target_count_recent = local_active;
    item.remote_tree_interactions_recent = static_cast<std::uint64_t>(
        std::llround(local_density_d * std::log2(local_density_d + 1.0)));
    item.work_units = 1.0;
    item.memory_bytes = particle_memory_bytes(species_tag);
    item.work_components = parallel::DecompositionWorkComponents{
        .particle_count_cost = 1.0,
        .gas_cell_cost = (species_tag == static_cast<std::uint32_t>(core::ParticleSpecies::kGas))
            ? (1.0 + static_cast<double>(local_gas))
            : 0.0,
        .tree_interaction_cost = static_cast<double>(item.remote_tree_interactions_recent),
        .pm_mesh_cost = static_cast<double>(pm_load),
        .amr_patch_cost = amr_patch_cost,
        .active_fraction_cost = static_cast<double>(local_active),
        .memory_pressure_cost = static_cast<double>(item.memory_bytes),
        .gpu_occupancy_cost = 0.0,
        .generic_work_cost = 1.0 + std::sqrt(local_density_d),
        .has_explicit_components = true,
    };
    items.push_back(item);
    patch_index_by_item.push_back(k_invalid_patch_index);
  }

  if (state.patches.owning_rank.size() != state.patches.size()) {
    state.patches.owning_rank.assign(state.patches.size(), static_cast<std::uint32_t>(std::max(world_rank, 0)));
  }
  for (std::size_t patch_index = 0; patch_index < state.patches.size(); ++patch_index) {
    if (state.patches.cell_count[patch_index] == 0U) {
      continue;
    }
    const std::uint32_t first_cell = state.patches.first_cell[patch_index];
    if (first_cell >= state.cells.size()) {
      continue;
    }
    parallel::DecompositionItem patch_item;
    patch_item.entity_id = state.patches.patch_id[patch_index];
    patch_item.kind = parallel::DecompositionEntityKind::kAmrPatch;
    patch_item.current_owner_rank = static_cast<int>(state.patches.owning_rank[patch_index]);
    patch_item.x_comov = state.cells.center_x_comoving[first_cell];
    patch_item.y_comov = state.cells.center_y_comoving[first_cell];
    patch_item.z_comov = state.cells.center_z_comoving[first_cell];
    patch_item.memory_bytes = static_cast<std::uint64_t>(state.patches.cell_count[patch_index]) *
        static_cast<std::uint64_t>(sizeof(double) * 8U + sizeof(std::uint32_t));
    patch_item.work_components = parallel::DecompositionWorkComponents{
        .amr_patch_cost = static_cast<double>(state.patches.cell_count[patch_index]) *
            (1.0 + static_cast<double>(std::max(state.patches.level[patch_index], 0))),
        .memory_pressure_cost = static_cast<double>(patch_item.memory_bytes),
        .generic_work_cost = static_cast<double>(state.patches.cell_count[patch_index]),
        .has_explicit_components = true,
    };
    items.push_back(patch_item);
    patch_index_by_item.push_back(static_cast<std::uint32_t>(patch_index));
  }

  parallel::DecompositionConfig decomposition_config;
  decomposition_config.world_size = world_size;
  decomposition_config.domain_x_min_comov = 0.0;
  decomposition_config.domain_x_max_comov = config.cosmology.box_size_x_mpc_comoving;
  decomposition_config.domain_y_min_comov = 0.0;
  decomposition_config.domain_y_max_comov = config.cosmology.box_size_y_mpc_comoving;
  decomposition_config.domain_z_min_comov = 0.0;
  decomposition_config.domain_z_max_comov = config.cosmology.box_size_z_mpc_comoving;
  decomposition_config.owned_particle_weight = 0.0;
  decomposition_config.active_target_weight = 0.0;
  decomposition_config.remote_tree_interaction_weight = 0.0;
  decomposition_config.work_weight = 0.0;
  decomposition_config.memory_weight = 0.0;
  decomposition_config.component_weights = parallel::DecompositionWeightCoefficients{
      .particle_count = config.parallel.decomposition_particle_count_weight,
      .gas_cell = config.parallel.decomposition_gas_cell_weight,
      .tree_interaction = config.parallel.decomposition_tree_interaction_weight,
      .pm_mesh = config.parallel.decomposition_pm_mesh_weight,
      .amr_patch = config.parallel.decomposition_amr_patch_weight,
      .active_fraction = config.parallel.decomposition_active_fraction_weight,
      .memory_pressure = config.parallel.decomposition_memory_pressure_weight,
      .gpu_occupancy = config.parallel.decomposition_gpu_occupancy_weight,
      .generic_work = config.parallel.decomposition_generic_work_weight,
  };
  const auto plan = parallel::buildMortonSfcDecomposition(items, decomposition_config);
  for (std::size_t item_index = 0; item_index < state.particles.size(); ++item_index) {
    state.particle_sidecar.owning_rank[item_index] = static_cast<std::uint32_t>(plan.owning_rank_by_item[item_index]);
  }
  for (std::size_t item_index = 0; item_index < patch_index_by_item.size(); ++item_index) {
    const std::uint32_t patch_index = patch_index_by_item[item_index];
    if (patch_index != k_invalid_patch_index && patch_index < state.patches.owning_rank.size()) {
      state.patches.owning_rank[patch_index] = static_cast<std::uint32_t>(plan.owning_rank_by_item[item_index]);
    }
  }
  parallel::recordDistributedProfiling(profiler, plan.metrics, 0, 0);
  if (profiler != nullptr) {
    profiler->recordEvent(core::RuntimeEvent{
        .event_kind = "parallel.decomposition.work_weighted",
        .severity = core::RuntimeEventSeverity::kInfo,
        .subsystem = "parallel.domain_decomposition",
        .message = "initial domain decomposition used explicit work-weight components",
        .payload = {{"world_size", std::to_string(world_size)},
                    {"item_count", std::to_string(items.size())},
                    {"weighted_imbalance_ratio", std::to_string(plan.metrics.weighted_imbalance_ratio)},
                    {"memory_imbalance_ratio", std::to_string(plan.metrics.memory_imbalance_ratio)}}});
  }
  compactStateToCurrentOwner(state, world_rank);
}



[[nodiscard]] double wrapPeriodicPosition(double position_comoving, double box_size_mpc_comoving) {
  if (box_size_mpc_comoving <= 0.0) {
    return position_comoving;
  }
  double wrapped = std::fmod(position_comoving, box_size_mpc_comoving);
  if (wrapped < 0.0) {
    wrapped += box_size_mpc_comoving;
  }
  if (wrapped >= box_size_mpc_comoving) {
    wrapped -= box_size_mpc_comoving;
  }
  return wrapped;
}

void seedParticleOwnershipFromPmSlabs(
    core::SimulationState& state,
    std::size_t pm_grid_nx,
    int world_size,
    double box_size_x_mpc_comoving) {
  if (pm_grid_nx == 0) {
    throw std::invalid_argument("seedParticleOwnershipFromPmSlabs requires pm_grid_nx > 0");
  }
  if (world_size <= 0) {
    throw std::invalid_argument("seedParticleOwnershipFromPmSlabs requires world_size > 0");
  }
  for (std::size_t i = 0; i < state.particles.size(); ++i) {
    const double wrapped_x = wrapPeriodicPosition(state.particles.position_x_comoving[i], box_size_x_mpc_comoving);
    std::size_t global_x = 0;
    if (box_size_x_mpc_comoving > 0.0) {
      const double scaled = (wrapped_x / box_size_x_mpc_comoving) * static_cast<double>(pm_grid_nx);
      global_x = std::min(pm_grid_nx - 1U, static_cast<std::size_t>(scaled));
    }
    state.particle_sidecar.owning_rank[i] = static_cast<std::uint32_t>(
        parallel::pmOwnerRankForGlobalX(pm_grid_nx, world_size, global_x));
  }
}


[[nodiscard]] parallel::GhostLayerEpoch makeRuntimeGhostLayerEpoch(const core::StepContext& context) {
  return parallel::GhostLayerEpoch{
      .decomposition_epoch = context.state.particleIndexGeneration(),
      .ghost_sync_epoch = context.integrator_state.step_index * core::integrationStageCount() +
          core::integrationStageIndex(context.stage) + 1U,
      .particle_index_generation = context.state.particleIndexGeneration(),
  };
}

[[nodiscard]] std::vector<parallel::LocalGhostDescriptor> buildParticleGhostDescriptors(
    const core::SimulationState& state,
    int world_rank,
    const parallel::GhostLayerEpoch& epoch) {
  std::vector<parallel::LocalGhostDescriptor> descriptors;
  descriptors.reserve(state.particles.size());
  for (std::size_t particle_index = 0; particle_index < state.particles.size(); ++particle_index) {
    const int owner_rank = static_cast<int>(state.particle_sidecar.owning_rank[particle_index]);
    descriptors.push_back(parallel::LocalGhostDescriptor{
        .residency = (owner_rank == world_rank) ? parallel::LocalIndexResidency::kOwned
                                                : parallel::LocalIndexResidency::kGhost,
        .owning_rank = owner_rank,
        .particle_id = state.particle_sidecar.particle_id[particle_index],
        .epoch = epoch,
    });
  }
  return descriptors;
}

[[nodiscard]] parallel::GhostExchangeBufferSoA buildParticleGhostPayloadState(
    const core::SimulationState& state,
    const parallel::GhostLayerEpoch& epoch) {
  parallel::GhostExchangeBufferSoA payload;
  const std::size_t particle_count = state.particles.size();
  payload.epoch = epoch;
  payload.entity_id.assign(state.particle_sidecar.particle_id.begin(), state.particle_sidecar.particle_id.end());
  payload.position_x_comoving.assign(state.particles.position_x_comoving.begin(), state.particles.position_x_comoving.end());
  payload.position_y_comoving.assign(state.particles.position_y_comoving.begin(), state.particles.position_y_comoving.end());
  payload.position_z_comoving.assign(state.particles.position_z_comoving.begin(), state.particles.position_z_comoving.end());
  payload.mass_code.assign(state.particles.mass_code.begin(), state.particles.mass_code.end());
  payload.velocity_x_code.assign(state.particles.velocity_x_peculiar.begin(), state.particles.velocity_x_peculiar.end());
  payload.velocity_y_code.assign(state.particles.velocity_y_peculiar.begin(), state.particles.velocity_y_peculiar.end());
  payload.velocity_z_code.assign(state.particles.velocity_z_peculiar.begin(), state.particles.velocity_z_peculiar.end());
  payload.density_code.assign(particle_count, 0.0);
  payload.pressure_code.assign(particle_count, 0.0);
  payload.internal_energy_code.assign(particle_count, 0.0);

  for (std::size_t cell_index = 0; cell_index < state.gas_cells.size(); ++cell_index) {
    const std::uint64_t parent_id = state.gas_cells.parent_particle_id[cell_index];
    const auto it = std::find(state.particle_sidecar.particle_id.begin(), state.particle_sidecar.particle_id.end(), parent_id);
    if (it == state.particle_sidecar.particle_id.end()) {
      continue;
    }
    const auto particle_index = static_cast<std::size_t>(std::distance(state.particle_sidecar.particle_id.begin(), it));
    payload.density_code[particle_index] = state.gas_cells.density_code[cell_index];
    payload.pressure_code[particle_index] = state.gas_cells.pressure_code[cell_index];
    payload.internal_energy_code[particle_index] = state.gas_cells.internal_energy_code[cell_index];
  }
  if (!payload.isConsistent() || !payload.hasGravityPayload()) {
    throw std::runtime_error("particle ghost payload construction produced inconsistent gravity lanes");
  }
  return payload;
}

void applyCommittedParticleGhostPayload(
    core::SimulationState& state,
    int world_rank,
    const std::vector<parallel::LocalGhostDescriptor>& descriptors,
    const parallel::GhostExchangeBufferSoA& payload) {
  if (!payload.isConsistent() || !payload.hasGravityPayload() || !payload.hasHydroPayload()) {
    throw std::invalid_argument("committed particle ghost payload must contain gravity and hydro lanes");
  }
  if (descriptors.size() != state.particles.size() || payload.size() < descriptors.size()) {
    throw std::invalid_argument("committed particle ghost payload shape does not match SimulationState particles");
  }

  std::unordered_map<std::uint64_t, std::size_t> particle_row_by_id;
  particle_row_by_id.reserve(state.particles.size());
  for (std::size_t particle_index = 0; particle_index < state.particles.size(); ++particle_index) {
    particle_row_by_id.emplace(state.particle_sidecar.particle_id[particle_index], particle_index);
  }

  for (std::size_t particle_index = 0; particle_index < descriptors.size(); ++particle_index) {
    const auto& descriptor = descriptors[particle_index];
    if (descriptor.residency != parallel::LocalIndexResidency::kGhost) {
      continue;
    }
    if (descriptor.owning_rank == world_rank) {
      throw std::logic_error("local ghost descriptor cannot be owned by local rank during commit");
    }
    if (payload.entity_id[particle_index] != descriptor.particle_id) {
      throw std::runtime_error("committed ghost payload entity_id drifted from descriptor particle_id");
    }
    // These rows are imported, non-authoritative ghost copies. Updating them is
    // allowed only because the authoritative owner remains descriptor.owning_rank.
    state.particles.position_x_comoving[particle_index] = payload.position_x_comoving[particle_index];
    state.particles.position_y_comoving[particle_index] = payload.position_y_comoving[particle_index];
    state.particles.position_z_comoving[particle_index] = payload.position_z_comoving[particle_index];
    state.particles.mass_code[particle_index] = payload.mass_code[particle_index];
    state.particles.velocity_x_peculiar[particle_index] = payload.velocity_x_code[particle_index];
    state.particles.velocity_y_peculiar[particle_index] = payload.velocity_y_code[particle_index];
    state.particles.velocity_z_peculiar[particle_index] = payload.velocity_z_code[particle_index];
  }

  for (std::size_t cell_index = 0; cell_index < state.gas_cells.size(); ++cell_index) {
    const std::uint64_t parent_id = state.gas_cells.parent_particle_id[cell_index];
    const auto particle_it = particle_row_by_id.find(parent_id);
    if (particle_it == particle_row_by_id.end()) {
      continue;
    }
    const std::size_t particle_index = particle_it->second;
    if (descriptors[particle_index].residency != parallel::LocalIndexResidency::kGhost) {
      continue;
    }
    state.cells.center_x_comoving[cell_index] = payload.position_x_comoving[particle_index];
    state.cells.center_y_comoving[cell_index] = payload.position_y_comoving[particle_index];
    state.cells.center_z_comoving[cell_index] = payload.position_z_comoving[particle_index];
    state.cells.mass_code[cell_index] = payload.mass_code[particle_index];
    state.gas_cells.density_code[cell_index] = payload.density_code[particle_index];
    state.gas_cells.pressure_code[cell_index] = payload.pressure_code[particle_index];
    state.gas_cells.internal_energy_code[cell_index] = payload.internal_energy_code[particle_index];
  }
}

struct SolverGhostRefreshReport {
  std::uint64_t sent_bytes = 0;
  std::uint64_t received_bytes = 0;
  std::size_t committed_slots = 0;
};

struct HydroGhostConservedSnapshot {
  std::vector<std::uint32_t> cell_indices;
  std::vector<std::uint64_t> gas_cell_ids;
  std::vector<std::uint64_t> parent_particle_ids;
  std::vector<std::uint32_t> owner_ranks;
  std::vector<hydro::HydroConservedState> conserved_state;
};

struct HydroConservativeGhostSyncReport {
  std::size_t restored_ghost_cells = 0;
  double rejected_remote_delta_l1 = 0.0;
  std::vector<parallel::HydroConservativeFluxCorrectionRecord> correction_records;
};

[[nodiscard]] HydroGhostConservedSnapshot snapshotHydroGhostConservedCells(
    const core::SimulationState& state,
    const hydro::HydroConservedStateSoa& conserved,
    std::span<const std::uint32_t> geometry_row_by_dense_row,
    std::uint32_t world_rank) {
  HydroGhostConservedSnapshot snapshot;
  state.requireGasCellIdentityMapCoversDenseRows("snapshotHydroGhostConservedCells");
  const auto particle_row_by_id = buildParticleRowById(state);
  snapshot.cell_indices.reserve(state.cells.size());
  snapshot.gas_cell_ids.reserve(state.cells.size());
  snapshot.parent_particle_ids.reserve(state.cells.size());
  snapshot.owner_ranks.reserve(state.cells.size());
  snapshot.conserved_state.reserve(state.cells.size());
  for (std::uint32_t cell_index = 0; cell_index < state.cells.size(); ++cell_index) {
    const std::uint32_t owner_rank = gasCellOwnerRankForLocalRow(
        state,
        cell_index,
        particle_row_by_id,
        "snapshotHydroGhostConservedCells");
    if (owner_rank == world_rank) {
      continue;
    }
    const core::GasCellIdentityRecord& identity = gasCellIdentityRecordForLocalRow(
        state,
        cell_index,
        "snapshotHydroGhostConservedCells");
    snapshot.cell_indices.push_back(cell_index);
    snapshot.gas_cell_ids.push_back(identity.gas_cell_id);
    snapshot.parent_particle_ids.push_back(identity.parent_particle_id.value_or(0U));
    snapshot.owner_ranks.push_back(owner_rank);
    if (cell_index >= geometry_row_by_dense_row.size()) {
      throw std::out_of_range("snapshotHydroGhostConservedCells: dense row is outside Cartesian geometry map");
    }
    snapshot.conserved_state.push_back(conserved.loadCell(geometry_row_by_dense_row[cell_index]));
  }
  return snapshot;
}

[[nodiscard]] HydroConservativeGhostSyncReport restoreHydroGhostConservedCells(
    hydro::HydroConservedStateSoa& conserved,
    const HydroGhostConservedSnapshot& snapshot,
    std::uint32_t world_rank) {
  if (snapshot.cell_indices.size() != snapshot.conserved_state.size() ||
      snapshot.cell_indices.size() != snapshot.gas_cell_ids.size() ||
      snapshot.cell_indices.size() != snapshot.parent_particle_ids.size() ||
      snapshot.cell_indices.size() != snapshot.owner_ranks.size()) {
    throw std::invalid_argument("hydro ghost snapshot is inconsistent");
  }
  HydroConservativeGhostSyncReport report;
  report.correction_records.reserve(snapshot.cell_indices.size());
  for (std::size_t i = 0; i < snapshot.cell_indices.size(); ++i) {
    const std::uint32_t cell_index = snapshot.cell_indices[i];
    if (cell_index >= conserved.size()) {
      throw std::out_of_range("hydro ghost snapshot cell index is outside conserved state");
    }
    const hydro::HydroConservedState before = snapshot.conserved_state[i];
    const hydro::HydroConservedState after = conserved.loadCell(cell_index);
    const double dm = after.mass_density_comoving - before.mass_density_comoving;
    const double dmx = after.momentum_density_x_comoving - before.momentum_density_x_comoving;
    const double dmy = after.momentum_density_y_comoving - before.momentum_density_y_comoving;
    const double dmz = after.momentum_density_z_comoving - before.momentum_density_z_comoving;
    const double de = after.total_energy_density_comoving - before.total_energy_density_comoving;
    report.rejected_remote_delta_l1 += std::abs(dm);
    report.rejected_remote_delta_l1 += std::abs(dmx);
    report.rejected_remote_delta_l1 += std::abs(dmy);
    report.rejected_remote_delta_l1 += std::abs(dmz);
    report.rejected_remote_delta_l1 += std::abs(de);
    if (dm != 0.0 || dmx != 0.0 || dmy != 0.0 || dmz != 0.0 || de != 0.0) {
      parallel::HydroConservativeFluxCorrectionRecord record;
      record.gas_cell_id = snapshot.gas_cell_ids[i];
      record.parent_particle_id = snapshot.parent_particle_ids[i];
      record.source_rank = static_cast<int>(world_rank);
      record.owner_rank = static_cast<int>(snapshot.owner_ranks[i]);
      record.delta_mass_density_comoving = dm;
      record.delta_momentum_density_x_comoving = dmx;
      record.delta_momentum_density_y_comoving = dmy;
      record.delta_momentum_density_z_comoving = dmz;
      record.delta_total_energy_density_comoving = de;
      parallel::validateHydroConservativeFluxCorrectionRecord(record);
      report.correction_records.push_back(record);
    }
    conserved.storeCell(cell_index, before);
    ++report.restored_ghost_cells;
  }
  return report;
}

[[nodiscard]] SolverGhostRefreshReport refreshParticleGhostsForSolver(
    core::StepContext& context,
    const parallel::MpiContext& mpi_context,
    std::string_view subsystem_name,
    parallel::GhostCacheLifecycle* lifecycle = nullptr) {
  const int world_rank = mpi_context.worldRank();
  const parallel::GhostLayerEpoch epoch = makeRuntimeGhostLayerEpoch(context);
  if (lifecycle != nullptr) {
    parallel::invalidateGhostCache(*lifecycle, epoch);
  }
  std::vector<parallel::LocalGhostDescriptor> descriptors = buildParticleGhostDescriptors(
      context.state, world_rank, epoch);
  parallel::GhostExchangeBufferSoA ghost_storage = buildParticleGhostPayloadState(context.state, epoch);

  const auto exchange = parallel::executeBlockingGhostRefreshExchangeFromDescriptors(
      mpi_context, descriptors, ghost_storage, epoch);
  const auto commit_report = parallel::commitBlockingGhostRefreshResult(
      ghost_storage, descriptors, exchange.plan, exchange.result, epoch);
  if (lifecycle != nullptr) {
    parallel::markGhostCacheCommitted(*lifecycle, epoch);
    parallel::requireValidGhostCache(*lifecycle, epoch, std::string(subsystem_name));
  }
  applyCommittedParticleGhostPayload(context.state, world_rank, descriptors, ghost_storage);

  if (context.profiler_session != nullptr) {
    context.profiler_session->recordEvent(core::RuntimeEvent{
        .event_kind = "parallel.blocking_ghost_refresh",
        .severity = core::RuntimeEventSeverity::kInfo,
        .subsystem = std::string(subsystem_name),
        .step_index = context.integrator_state.step_index,
        .simulation_time_code = context.integrator_state.current_time_code,
        .scale_factor = context.integrator_state.current_scale_factor,
        .message = "blocking correctness-first particle ghost refresh completed before solver access",
        .payload = {
            {"neighbor_count", std::to_string(exchange.plan.neighbor_ranks.size())},
            {"sent_bytes", std::to_string(exchange.result.sent_bytes)},
            {"received_bytes", std::to_string(exchange.result.received_bytes)},
            {"committed_ghost_slots", std::to_string(commit_report.updated_ghost_slots)},
            {"ghost_sync_epoch", std::to_string(epoch.ghost_sync_epoch)},
            {"ghost_cache_refresh_count", lifecycle != nullptr ? std::to_string(lifecycle->refresh_count) : "0"},
            {"ghost_cache_invalidation_count", lifecycle != nullptr ? std::to_string(lifecycle->invalidation_count) : "0"},
        },
    });
  }
  return SolverGhostRefreshReport{
      .sent_bytes = exchange.result.sent_bytes,
      .received_bytes = exchange.result.received_bytes,
      .committed_slots = commit_report.updated_ghost_slots,
  };
}

[[nodiscard]] core::ProvenanceRecord makeGravityAwareProvenanceRecord(
    const core::FrozenConfig& frozen_config,
    const core::SimulationConfig& config) {
  core::ProvenanceRecord record = core::makeProvenanceRecord(
      frozen_config.provenance.config_hash_hex, "unknown");
  record.config_schema_name = "cosmosim_config";
  record.config_schema_version = std::to_string(frozen_config.config.schema_version);
  record.raw_input_config = frozen_config.raw_text;
  record.normalized_config = frozen_config.normalized_text;
  record.derived_runtime_state =
      core::serializeDerivedRuntimeConfig(core::deriveRuntimeConfig(frozen_config));
  const double dx = config.cosmology.box_size_x_mpc_comoving /
      static_cast<double>(config.numerics.treepm_pm_grid_nx);
  const double dy = config.cosmology.box_size_y_mpc_comoving /
      static_cast<double>(config.numerics.treepm_pm_grid_ny);
  const double dz = config.cosmology.box_size_z_mpc_comoving /
      static_cast<double>(config.numerics.treepm_pm_grid_nz);
  const double mesh_spacing_mpc_comoving = std::cbrt(dx * dy * dz);
  const gravity::TreePmSplitPolicy split_policy = gravity::makeTreePmSplitPolicyFromMeshSpacing(
      config.numerics.treepm_asmth_cells,
      config.numerics.treepm_rcut_cells,
      mesh_spacing_mpc_comoving);
  record.gravity_treepm_pm_grid = config.numerics.treepm_pm_grid_nx;
  record.gravity_treepm_pm_grid_nx = config.numerics.treepm_pm_grid_nx;
  record.gravity_treepm_pm_grid_ny = config.numerics.treepm_pm_grid_ny;
  record.gravity_treepm_pm_grid_nz = config.numerics.treepm_pm_grid_nz;
  record.gravity_treepm_assignment_scheme =
      treePmAssignmentSchemeName(config.numerics.treepm_assignment_scheme);
  record.gravity_treepm_window_deconvolution =
      config.numerics.treepm_enable_window_deconvolution;
  record.gravity_treepm_asmth_cells = config.numerics.treepm_asmth_cells;
  record.gravity_treepm_rcut_cells = config.numerics.treepm_rcut_cells;
  record.gravity_treepm_mesh_spacing_mpc_comoving = mesh_spacing_mpc_comoving;
  record.gravity_treepm_mesh_spacing_x_mpc_comoving = dx;
  record.gravity_treepm_mesh_spacing_y_mpc_comoving = dy;
  record.gravity_treepm_mesh_spacing_z_mpc_comoving = dz;
  record.gravity_treepm_split_scale_mpc_comoving = split_policy.split_scale_comoving;
  record.gravity_treepm_cutoff_radius_mpc_comoving = split_policy.cutoff_radius_comoving;
  record.gravity_treepm_update_cadence_steps = config.numerics.treepm_update_cadence_steps;
  record.gravity_treepm_pm_decomposition_mode =
      pmDecompositionModeName(config.numerics.treepm_pm_decomposition_mode);
  record.gravity_treepm_tree_exchange_batch_bytes =
      config.numerics.treepm_tree_exchange_batch_bytes;
  record.gravity_softening_policy = describeSofteningPolicy(config);
  record.gravity_softening_kernel = "plummer";
  record.gravity_softening_epsilon_kpc_comoving = config.numerics.gravity_softening_kpc_comoving;
  record.gravity_pm_fft_backend = gravity::PmSolver::fftBackendName();
  switch (config.mode.zoom_long_range_strategy) {
    case core::ZoomLongRangeStrategy::kDisabled:
      record.zoom_long_range_strategy = "disabled";
      break;
    case core::ZoomLongRangeStrategy::kGlobalCoarsePlusFocusedHighResCorrection:
      record.zoom_long_range_strategy = "global_coarse_plus_focused_highres_correction";
      break;
  }
  record.zoom_region_center_x_mpc_comoving = config.mode.zoom_region_center_x_mpc_comoving;
  record.zoom_region_center_y_mpc_comoving = config.mode.zoom_region_center_y_mpc_comoving;
  record.zoom_region_center_z_mpc_comoving = config.mode.zoom_region_center_z_mpc_comoving;
  record.zoom_region_radius_mpc_comoving = config.mode.zoom_region_radius_mpc_comoving;
  record.zoom_focused_pm_grid = std::to_string(config.mode.zoom_focused_pm_grid_nx) + "x" +
      std::to_string(config.mode.zoom_focused_pm_grid_ny) + "x" +
      std::to_string(config.mode.zoom_focused_pm_grid_nz);
  record.zoom_contamination_radius_mpc_comoving = config.mode.zoom_contamination_radius_mpc_comoving;
  return record;
}

[[nodiscard]] std::string pmSlabSignature(const parallel::DistributedRestartState& distributed_state) {
  std::ostringstream stream;
  for (std::size_t rank = 0; rank < distributed_state.pm_slab_begin_x_by_rank.size(); ++rank) {
    if (rank > 0) {
      stream << ';';
    }
    stream << rank << ':' << distributed_state.pm_slab_begin_x_by_rank[rank] << '-' <<
        distributed_state.pm_slab_end_x_by_rank[rank];
  }
  return stream.str();
}

[[nodiscard]] std::filesystem::path computeRunDirectory(
    const core::SimulationConfig& config,
    const std::filesystem::path* output_root_override) {
  const std::filesystem::path output_root =
      (output_root_override != nullptr && !output_root_override->empty())
      ? *output_root_override
      : std::filesystem::path(config.output.output_directory);
  return output_root / config.output.run_name;
}

[[nodiscard]] std::filesystem::path resolveConfigRelativePath(
    const core::FrozenConfig& frozen_config,
    const std::filesystem::path& candidate) {
  if (candidate.is_absolute()) {
    return candidate;
  }

  const std::filesystem::path source_path(frozen_config.provenance.source_name);
  if (!source_path.empty() && source_path.has_parent_path()) {
    return source_path.parent_path() / candidate;
  }
  return candidate;
}

[[nodiscard]] bool repeatedCanonicalOrder(const std::vector<std::string>& observed) {
  const auto canonical = core::StageScheduler::kickDriftKickOrder();
  if (observed.empty() || observed.size() % canonical.size() != 0) {
    return false;
  }

  for (std::size_t chunk = 0; chunk < observed.size(); chunk += canonical.size()) {
    for (std::size_t i = 0; i < canonical.size(); ++i) {
      if (observed[chunk + i] != core::integrationStageName(canonical[i])) {
        return false;
      }
    }
  }
  return true;
}

[[nodiscard]] std::string formatIndexedFileStem(std::string_view stem, std::uint64_t index) {
  std::ostringstream out;
  out << stem << '_' << std::setw(3) << std::setfill('0') << index << ".hdf5";
  return out.str();
}

void fnv1aMix(std::uint64_t& hash, const void* data, std::size_t size) {
  constexpr std::uint64_t k_fnv_offset_basis = 1469598103934665603ULL;
  constexpr std::uint64_t k_fnv_prime = 1099511628211ULL;
  if (hash == 0) {
    hash = k_fnv_offset_basis;
  }
  const auto* bytes = static_cast<const unsigned char*>(data);
  for (std::size_t i = 0; i < size; ++i) {
    hash ^= static_cast<std::uint64_t>(bytes[i]);
    hash *= k_fnv_prime;
  }
}

void fnv1aMix(std::uint64_t& hash, double value) { fnv1aMix(hash, &value, sizeof(value)); }

void fnv1aMix(std::uint64_t& hash, std::uint64_t value) { fnv1aMix(hash, &value, sizeof(value)); }

[[nodiscard]] std::uint64_t computeStateDigest(const core::SimulationState& state, const core::IntegratorState& integrator_state) {
  std::uint64_t hash = 0;
  fnv1aMix(hash, static_cast<std::uint64_t>(state.particles.size()));
  fnv1aMix(hash, static_cast<std::uint64_t>(state.cells.size()));
  for (const double value : state.particles.position_x_comoving) {
    fnv1aMix(hash, value);
  }
  for (const double value : state.particles.position_y_comoving) {
    fnv1aMix(hash, value);
  }
  for (const double value : state.particles.position_z_comoving) {
    fnv1aMix(hash, value);
  }
  for (const double value : state.particles.velocity_x_peculiar) {
    fnv1aMix(hash, value);
  }
  for (const double value : state.particles.velocity_y_peculiar) {
    fnv1aMix(hash, value);
  }
  for (const double value : state.particles.velocity_z_peculiar) {
    fnv1aMix(hash, value);
  }
  for (const double value : state.particles.mass_code) {
    fnv1aMix(hash, value);
  }
  for (const double value : state.gas_cells.density_code) {
    fnv1aMix(hash, value);
  }
  for (const double value : state.gas_cells.internal_energy_code) {
    fnv1aMix(hash, value);
  }
  fnv1aMix(hash, integrator_state.current_time_code);
  fnv1aMix(hash, integrator_state.current_scale_factor);
  fnv1aMix(hash, integrator_state.step_index);
  return hash;
}

[[nodiscard]] std::uint64_t computeParticleIdSum(const core::SimulationState& state) {
  std::uint64_t sum = 0;
  for (const std::uint64_t particle_id : state.particle_sidecar.particle_id) {
    sum += particle_id;
  }
  return sum;
}

[[nodiscard]] std::uint64_t computeParticleIdSquareSum(const core::SimulationState& state) {
  std::uint64_t sum = 0;
  for (const std::uint64_t particle_id : state.particle_sidecar.particle_id) {
    sum += particle_id * particle_id;
  }
  return sum;
}

[[nodiscard]] std::uint64_t computeParticleIdXor(const core::SimulationState& state) {
  std::uint64_t value = 0;
  for (const std::uint64_t particle_id : state.particle_sidecar.particle_id) {
    value ^= particle_id;
  }
  return value;
}

[[nodiscard]] bool moduleSidecarsEqualForRestart(
    const core::ModuleSidecarRegistry& lhs,
    const core::ModuleSidecarRegistry& rhs) {
  const auto lhs_blocks = lhs.blocksSortedByName();
  const auto rhs_blocks = rhs.blocksSortedByName();
  if (lhs_blocks.size() != rhs_blocks.size()) {
    return false;
  }
  for (std::size_t i = 0; i < lhs_blocks.size(); ++i) {
    if (lhs_blocks[i]->module_name != rhs_blocks[i]->module_name ||
        lhs_blocks[i]->schema_version != rhs_blocks[i]->schema_version ||
        lhs_blocks[i]->payload != rhs_blocks[i]->payload) {
      return false;
    }
  }
  return true;
}

[[nodiscard]] bool restartRuntimeStateExactlyEquivalent(
    const core::SimulationState& restored,
    const core::SimulationState& reference) {
  return restored.metadata.serialize() == reference.metadata.serialize() &&
      restored.particles.position_x_comoving == reference.particles.position_x_comoving &&
      restored.particles.position_y_comoving == reference.particles.position_y_comoving &&
      restored.particles.position_z_comoving == reference.particles.position_z_comoving &&
      restored.particles.velocity_x_peculiar == reference.particles.velocity_x_peculiar &&
      restored.particles.velocity_y_peculiar == reference.particles.velocity_y_peculiar &&
      restored.particles.velocity_z_peculiar == reference.particles.velocity_z_peculiar &&
      restored.particles.mass_code == reference.particles.mass_code &&
      restored.particles.time_bin == reference.particles.time_bin &&
      restored.particle_sidecar.particle_id == reference.particle_sidecar.particle_id &&
      restored.particle_sidecar.sfc_key == reference.particle_sidecar.sfc_key &&
      restored.particle_sidecar.species_tag == reference.particle_sidecar.species_tag &&
      restored.particle_sidecar.particle_flags == reference.particle_sidecar.particle_flags &&
      restored.particle_sidecar.owning_rank == reference.particle_sidecar.owning_rank &&
      restored.particle_sidecar.gravity_softening_comoving == reference.particle_sidecar.gravity_softening_comoving &&
      restored.particle_sidecar.has_gravity_softening_override == reference.particle_sidecar.has_gravity_softening_override &&
      restored.cells.center_x_comoving == reference.cells.center_x_comoving &&
      restored.cells.center_y_comoving == reference.cells.center_y_comoving &&
      restored.cells.center_z_comoving == reference.cells.center_z_comoving &&
      restored.cells.mass_code == reference.cells.mass_code &&
      restored.cells.time_bin == reference.cells.time_bin &&
      restored.cells.patch_index == reference.cells.patch_index &&
      restored.gas_cells.gas_cell_id == reference.gas_cells.gas_cell_id &&
      restored.gas_cells.parent_particle_id == reference.gas_cells.parent_particle_id &&
      restored.gas_cells.density_code == reference.gas_cells.density_code &&
      restored.gas_cells.pressure_code == reference.gas_cells.pressure_code &&
      restored.gas_cells.internal_energy_code == reference.gas_cells.internal_energy_code &&
      restored.gas_cells.temperature_code == reference.gas_cells.temperature_code &&
      restored.gas_cells.sound_speed_code == reference.gas_cells.sound_speed_code &&
      restored.patches.patch_id == reference.patches.patch_id &&
      restored.patches.level == reference.patches.level &&
      restored.patches.first_cell == reference.patches.first_cell &&
      restored.patches.cell_count == reference.patches.cell_count &&
      restored.patches.owning_rank == reference.patches.owning_rank &&
      restored.star_particles.particle_index == reference.star_particles.particle_index &&
      restored.star_particles.formation_scale_factor == reference.star_particles.formation_scale_factor &&
      restored.star_particles.birth_mass_code == reference.star_particles.birth_mass_code &&
      restored.star_particles.metallicity_mass_fraction == reference.star_particles.metallicity_mass_fraction &&
      restored.star_particles.stellar_age_years_last == reference.star_particles.stellar_age_years_last &&
      restored.star_particles.stellar_returned_mass_cumulative_code == reference.star_particles.stellar_returned_mass_cumulative_code &&
      restored.star_particles.stellar_returned_metals_cumulative_code == reference.star_particles.stellar_returned_metals_cumulative_code &&
      restored.star_particles.stellar_feedback_energy_cumulative_erg == reference.star_particles.stellar_feedback_energy_cumulative_erg &&
      restored.star_particles.stellar_returned_mass_channel_cumulative_code == reference.star_particles.stellar_returned_mass_channel_cumulative_code &&
      restored.star_particles.stellar_returned_metals_channel_cumulative_code == reference.star_particles.stellar_returned_metals_channel_cumulative_code &&
      restored.star_particles.stellar_feedback_energy_channel_cumulative_erg == reference.star_particles.stellar_feedback_energy_channel_cumulative_erg &&
      restored.black_holes.particle_index == reference.black_holes.particle_index &&
      restored.black_holes.host_cell_index == reference.black_holes.host_cell_index &&
      restored.black_holes.subgrid_mass_code == reference.black_holes.subgrid_mass_code &&
      restored.black_holes.accretion_rate_code == reference.black_holes.accretion_rate_code &&
      restored.black_holes.feedback_energy_code == reference.black_holes.feedback_energy_code &&
      restored.black_holes.eddington_ratio == reference.black_holes.eddington_ratio &&
      restored.black_holes.cumulative_accreted_mass_code == reference.black_holes.cumulative_accreted_mass_code &&
      restored.black_holes.cumulative_feedback_energy_code == reference.black_holes.cumulative_feedback_energy_code &&
      restored.black_holes.duty_cycle_active_time_code == reference.black_holes.duty_cycle_active_time_code &&
      restored.black_holes.duty_cycle_total_time_code == reference.black_holes.duty_cycle_total_time_code &&
      restored.tracers.particle_index == reference.tracers.particle_index &&
      restored.tracers.parent_particle_id == reference.tracers.parent_particle_id &&
      restored.tracers.injection_step == reference.tracers.injection_step &&
      restored.tracers.host_cell_index == reference.tracers.host_cell_index &&
      restored.tracers.mass_fraction_of_host == reference.tracers.mass_fraction_of_host &&
      restored.tracers.last_host_mass_code == reference.tracers.last_host_mass_code &&
      restored.tracers.cumulative_exchanged_mass_code == reference.tracers.cumulative_exchanged_mass_code &&
      restored.species.count_by_species == reference.species.count_by_species &&
      moduleSidecarsEqualForRestart(restored.sidecars, reference.sidecars);
}

void ensureRunDirectory(const std::filesystem::path& run_directory) {
  std::filesystem::create_directories(run_directory);
}

void flushCommonArtifacts(
    const core::FrozenConfig& frozen_config,
    core::ProfilerSession& profiler,
    ReferenceWorkflowReport& report) {
  ensureRunDirectory(report.run_directory);

  report.normalized_config_snapshot_path = report.run_directory / "normalized_config.param.txt";
  if (!std::filesystem::exists(report.normalized_config_snapshot_path)) {
    core::writeNormalizedConfigSnapshot(frozen_config, report.run_directory);
  }
  report.normalized_config_snapshot_written = true;

  report.profiler_json_path = report.run_directory / "profile.json";
  report.profiler_csv_path = report.run_directory / "profile.csv";
  report.operational_report_json_path = report.run_directory / "operational_events.json";
  core::writeProfilerReportJson(profiler, report.profiler_json_path);
  core::writeProfilerReportCsv(profiler, report.profiler_csv_path);
  core::writeOperationalReportJson(
      profiler,
      report.operational_report_json_path,
      frozen_config.config.output.run_name,
      frozen_config.provenance.config_hash_hex);
}

[[nodiscard]] io::IcReadResult loadInitialConditions(const core::FrozenConfig& frozen_config) {
  const core::SimulationConfig& config = frozen_config.config;
  if (config.mode.ic_file == "generated") {
    return io::convertGeneratedIsolatedIcToState(config, k_default_generated_particle_axis);
  }

  const std::filesystem::path ic_path =
      resolveConfigRelativePath(frozen_config, std::filesystem::path(config.mode.ic_file));
  return io::readGadgetArepoHdf5Ic(ic_path, config);
}

[[nodiscard]] bool hasHdf5Extension(const std::filesystem::path& path) {
  const std::string ext = path.extension().string();
  return ext == ".h5" || ext == ".hdf5";
}

[[nodiscard]] std::unordered_set<std::uint64_t> loadZoomParticleIdsFromText(const std::filesystem::path& path) {
  std::ifstream in(path);
  if (!in) {
    throw std::runtime_error("failed to open zoom region file: " + path.string());
  }
  std::unordered_set<std::uint64_t> ids;
  std::string line;
  while (std::getline(in, line)) {
    const auto comment = line.find('#');
    if (comment != std::string::npos) {
      line.erase(comment);
    }
    std::istringstream stream(line);
    std::uint64_t id = 0;
    while (stream >> id) {
      ids.insert(id);
    }
  }
  return ids;
}

#if COSMOSIM_ENABLE_HDF5
[[nodiscard]] std::unordered_set<std::uint64_t> loadZoomParticleIdsFromHdf5(const std::filesystem::path& path) {
  hid_t file = H5Fopen(path.string().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file < 0) {
    throw std::runtime_error("failed to open zoom region HDF5 file: " + path.string());
  }
  auto close_file = [&]() { H5Fclose(file); };
  const char* candidates[] = {"/Region/ParticleIDs", "/ParticleIDs", "Region/ParticleIDs", "ParticleIDs"};
  hid_t dataset = -1;
  for (const char* candidate : candidates) {
    if (H5Lexists(file, candidate, H5P_DEFAULT) > 0) {
      dataset = H5Dopen2(file, candidate, H5P_DEFAULT);
      if (dataset >= 0) {
        break;
      }
    }
  }
  if (dataset < 0) {
    close_file();
    throw std::runtime_error("zoom region HDF5 file lacks ParticleIDs dataset: " + path.string());
  }
  hid_t space = H5Dget_space(dataset);
  if (space < 0) {
    H5Dclose(dataset); close_file();
    throw std::runtime_error("failed to query zoom region HDF5 dataspace");
  }
  hsize_t dims[1] = {0};
  if (H5Sget_simple_extent_ndims(space) != 1 || H5Sget_simple_extent_dims(space, dims, nullptr) != 1) {
    H5Sclose(space); H5Dclose(dataset); close_file();
    throw std::runtime_error("zoom region ParticleIDs dataset must be 1D");
  }
  std::vector<std::uint64_t> values(static_cast<std::size_t>(dims[0]), 0ULL);
  const hid_t dtype = H5Dget_type(dataset);
  const H5T_class_t cls = H5Tget_class(dtype);
  if (cls != H5T_INTEGER) {
    H5Tclose(dtype); H5Sclose(space); H5Dclose(dataset); close_file();
    throw std::runtime_error("zoom region ParticleIDs dataset must be integer typed");
  }
  if (!values.empty() && H5Dread(dataset, H5T_NATIVE_UINT64, H5S_ALL, H5S_ALL, H5P_DEFAULT, values.data()) < 0) {
    H5Tclose(dtype); H5Sclose(space); H5Dclose(dataset); close_file();
    throw std::runtime_error("failed to read zoom region ParticleIDs dataset");
  }
  H5Tclose(dtype); H5Sclose(space); H5Dclose(dataset); close_file();
  return std::unordered_set<std::uint64_t>(values.begin(), values.end());
}
#endif


void materializeRootHydroPatchIfMissing(
    core::SimulationState& state,
    const core::SimulationConfig& config) {
  if (state.cells.size() == 0U || state.patches.size() != 0U) {
    return;
  }
  state.patches.resize(1U);
  state.patches.patch_id[0] = 1ULL;
  state.patches.level[0] = 0;
  state.patches.first_cell[0] = 0U;
  state.patches.cell_count[0] = static_cast<std::uint32_t>(state.cells.size());
  state.patches.owning_rank[0] = 0U;
  if (state.cells.patch_index.size() != state.cells.size()) {
    state.cells.patch_index.resize(state.cells.size());
  }
  std::fill(state.cells.patch_index.begin(), state.cells.patch_index.end(), 0U);

  if (state.gas_cells.size() == state.cells.size()) {
    if (state.gas_cell_identity.empty()) {
      // Import-only bridge: bootstrap truth once from persisted sidecar lanes.
      state.refreshGasCellIdentityMapFromSidecarLanes();
    } else {
      std::vector<core::GasCellIdentityRecord> records(
          state.gas_cell_identity.records().begin(), state.gas_cell_identity.records().end());
      for (core::GasCellIdentityRecord& record : records) {
        record.owning_patch_id = state.patches.patch_id[0];
      }
      state.replaceGasCellIdentityRecords(std::move(records));
    }
  }

  // Materialize explicit geometry only after identity coverage exists and only
  // when the shared strict builder can prove the Cartesian topology.  Do not
  // manufacture a near-cubic shape from the number of storage rows.
  const internal::CartesianGasCellLayoutBuildResult layout =
      internal::buildCartesianGasCellRowLayout(state, config);
  if (layout.ok()) {
    state.patches.origin_x_comoving[0] = layout.layout.spec.origin_x_comoving;
    state.patches.origin_y_comoving[0] = layout.layout.spec.origin_y_comoving;
    state.patches.origin_z_comoving[0] = layout.layout.spec.origin_z_comoving;
    state.patches.extent_x_comoving[0] =
        layout.layout.spec.cell_width_x_comoving * static_cast<double>(layout.layout.spec.nx);
    state.patches.extent_y_comoving[0] =
        layout.layout.spec.cell_width_y_comoving * static_cast<double>(layout.layout.spec.ny);
    state.patches.extent_z_comoving[0] =
        layout.layout.spec.cell_width_z_comoving * static_cast<double>(layout.layout.spec.nz);
    state.patches.cell_dim_x[0] = static_cast<std::uint16_t>(layout.layout.spec.nx);
    state.patches.cell_dim_y[0] = static_cast<std::uint16_t>(layout.layout.spec.ny);
    state.patches.cell_dim_z[0] = static_cast<std::uint16_t>(layout.layout.spec.nz);
  }
  state.bumpCellIndexGeneration();
}

void finalizeStateMetadata(const core::FrozenConfig& frozen_config, core::SimulationState& state) {
  state.metadata.run_name = frozen_config.config.output.run_name;
  state.metadata.normalized_config_hash = frozen_config.provenance.config_hash;
  state.metadata.normalized_config_hash_hex = frozen_config.provenance.config_hash_hex;
  state.metadata.snapshot_stem = frozen_config.config.output.output_stem;
  state.metadata.restart_stem = frozen_config.config.output.restart_stem;
  state.metadata.scale_factor = (state.metadata.scale_factor > 0.0)
      ? state.metadata.scale_factor
      : std::max(1.0, frozen_config.config.numerics.t_code_begin);
  state.rebuildSpeciesIndex();
  materializeRootHydroPatchIfMissing(state, frozen_config.config);
}

void initializeSchedulerBins(
    const core::SimulationState& state,
    core::HierarchicalTimeBinScheduler& particle_scheduler,
    core::HierarchicalTimeBinScheduler& gas_cell_scheduler) {
  if (state.particles.size() > static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max()) ||
      state.cells.size() > static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max())) {
    throw std::overflow_error("reference workflow scheduler element count exceeds uint32 range");
  }
  const std::uint32_t particle_count = static_cast<std::uint32_t>(state.particles.size());
  const std::uint32_t cell_count = static_cast<std::uint32_t>(state.cells.size());
  if (particle_count == 0U && cell_count == 0U) {
    throw std::runtime_error("reference workflow requires non-empty initial conditions");
  }

  const std::uint8_t particle_default_bin = particle_scheduler.maxBin() > 0 ? 1U : 0U;
  const std::uint8_t gas_default_bin = gas_cell_scheduler.maxBin() > 0 ? 1U : 0U;
  particle_scheduler.reset(particle_count, particle_default_bin, 0U);
  gas_cell_scheduler.reset(cell_count, gas_default_bin, 0U);
  if (cell_count > 0U) {
    state.requireGasCellIdentityMapCoversDenseRows("initialize independent gas-cell scheduler");
    for (std::uint32_t cell_index = 0; cell_index < cell_count; ++cell_index) {
      gas_cell_scheduler.setElementBin(cell_index, 0U, gas_cell_scheduler.currentTick());
    }
  }
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
  scheduler.appendElements(
      required_count - scheduler.elementCount(), new_element_bin, scheduler.currentTick() + 1U);
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

// Particle migration/rebalance paths own only particle scheduler state.  Gas-cell
// mirrors are synchronized by the caller that owns the gas-cell scheduler.
void syncTimeBinsFromScheduler(
    const core::HierarchicalTimeBinScheduler& scheduler,
    core::SimulationState& state) {
  core::syncTimeBinMirrorsFromScheduler(
      scheduler, state, core::TimeBinMirrorDomain::kParticles);
}


class StageAuditCallback final : public core::IntegrationCallback {
 public:
  explicit StageAuditCallback(std::vector<std::string>* stage_sequence)
      : m_stage_sequence(stage_sequence) {}

  [[nodiscard]] std::string_view callbackName() const override { return "stage_audit"; }
  [[nodiscard]] std::span<const core::IntegrationStage> integrationStages() const override {
    return k_all_integration_stages;
  }
  [[nodiscard]] std::span<const core::StageContract> stageContracts() const override { return m_contracts; }

  void onStage(core::StepContext& context) override {
    m_stage_sequence->push_back(std::string(core::integrationStageName(context.stage)));
  }

 private:
  static constexpr std::array<core::StageContract, 8> m_contracts{{
      {.stage = core::IntegrationStage::kGravityKickPre, .required_inputs = core::StageDataDomain::kDiagnostics, .mutated_state = core::StageDataDomain::kDiagnostics, .produced_outputs = core::StageDataDomain::kDiagnostics, .allowed_side_effects = core::StageDataDomain::kDiagnostics, .sync_requirements = core::StageSyncRequirement::kNone, .active_set_family = core::StageActiveSetFamily::kNone, .restart_safety = core::StageSafety::kSafe, .output_safety = core::StageSafety::kSafe, .owner = core::StageSubsystem::kAnalysis},
      {.stage = core::IntegrationStage::kDrift, .required_inputs = core::StageDataDomain::kDiagnostics, .mutated_state = core::StageDataDomain::kDiagnostics, .produced_outputs = core::StageDataDomain::kDiagnostics, .allowed_side_effects = core::StageDataDomain::kDiagnostics, .sync_requirements = core::StageSyncRequirement::kNone, .active_set_family = core::StageActiveSetFamily::kNone, .restart_safety = core::StageSafety::kSafe, .output_safety = core::StageSafety::kSafe, .owner = core::StageSubsystem::kAnalysis},
      {.stage = core::IntegrationStage::kForceRefresh, .required_inputs = core::StageDataDomain::kDiagnostics, .mutated_state = core::StageDataDomain::kDiagnostics, .produced_outputs = core::StageDataDomain::kDiagnostics, .allowed_side_effects = core::StageDataDomain::kDiagnostics, .sync_requirements = core::StageSyncRequirement::kNone, .active_set_family = core::StageActiveSetFamily::kNone, .restart_safety = core::StageSafety::kSafe, .output_safety = core::StageSafety::kSafe, .owner = core::StageSubsystem::kAnalysis},
      {.stage = core::IntegrationStage::kHydroUpdate, .required_inputs = core::StageDataDomain::kDiagnostics, .mutated_state = core::StageDataDomain::kDiagnostics, .produced_outputs = core::StageDataDomain::kDiagnostics, .allowed_side_effects = core::StageDataDomain::kDiagnostics, .sync_requirements = core::StageSyncRequirement::kNone, .active_set_family = core::StageActiveSetFamily::kNone, .restart_safety = core::StageSafety::kSafe, .output_safety = core::StageSafety::kSafe, .owner = core::StageSubsystem::kAnalysis},
      {.stage = core::IntegrationStage::kSourceTerms, .required_inputs = core::StageDataDomain::kDiagnostics, .mutated_state = core::StageDataDomain::kDiagnostics, .produced_outputs = core::StageDataDomain::kDiagnostics, .allowed_side_effects = core::StageDataDomain::kDiagnostics, .sync_requirements = core::StageSyncRequirement::kNone, .active_set_family = core::StageActiveSetFamily::kNone, .restart_safety = core::StageSafety::kSafe, .output_safety = core::StageSafety::kSafe, .owner = core::StageSubsystem::kAnalysis},
      {.stage = core::IntegrationStage::kGravityKickPost, .required_inputs = core::StageDataDomain::kDiagnostics, .mutated_state = core::StageDataDomain::kDiagnostics, .produced_outputs = core::StageDataDomain::kDiagnostics, .allowed_side_effects = core::StageDataDomain::kDiagnostics, .sync_requirements = core::StageSyncRequirement::kNone, .active_set_family = core::StageActiveSetFamily::kNone, .restart_safety = core::StageSafety::kSafe, .output_safety = core::StageSafety::kSafe, .owner = core::StageSubsystem::kAnalysis},
      {.stage = core::IntegrationStage::kAnalysisHooks, .required_inputs = core::StageDataDomain::kDiagnostics, .mutated_state = core::StageDataDomain::kDiagnostics, .produced_outputs = core::StageDataDomain::kDiagnostics, .allowed_side_effects = core::StageDataDomain::kDiagnostics, .sync_requirements = core::StageSyncRequirement::kLocalOnly, .active_set_family = core::StageActiveSetFamily::kNone, .restart_safety = core::StageSafety::kSafe, .output_safety = core::StageSafety::kSafe, .owner = core::StageSubsystem::kAnalysis},
      {.stage = core::IntegrationStage::kOutputCheck, .required_inputs = core::StageDataDomain::kOutputState, .mutated_state = core::StageDataDomain::kDiagnostics, .produced_outputs = core::StageDataDomain::kDiagnostics, .allowed_side_effects = core::StageDataDomain::kDiagnostics, .sync_requirements = core::StageSyncRequirement::kGlobal, .active_set_family = core::StageActiveSetFamily::kOutputState, .restart_safety = core::StageSafety::kSafe, .output_safety = core::StageSafety::kSafe, .owner = core::StageSubsystem::kAnalysis},
  }};
  std::vector<std::string>* m_stage_sequence = nullptr;
};

class DriftCallback final : public core::IntegrationCallback {
 public:
  [[nodiscard]] std::string_view callbackName() const override { return "drift"; }
  [[nodiscard]] std::span<const core::IntegrationStage> integrationStages() const override { return k_drift_stage; }
  [[nodiscard]] std::span<const core::StageContract> stageContracts() const override { return m_contracts; }

  void onStage(core::StepContext& context) override {
    if (context.stage != core::IntegrationStage::kDrift) {
      throw std::logic_error("drift handler received an unregistered stage");
    }

    const double drift_factor = context.timeline_step.drift_factor_code;
    parallel::MpiContext mpi_context;
    const std::uint32_t world_rank = static_cast<std::uint32_t>(mpi_context.worldRank());
    for (const std::uint32_t particle_index : context.active_set.particle_indices) {
      if (particle_index >= context.state.particles.size()) {
        throw std::out_of_range("drift callback active particle index out of range");
      }
      if (particle_index >= context.state.particle_sidecar.owning_rank.size()) {
        throw std::out_of_range("drift callback owning-rank sidecar index out of range");
      }
      if (context.state.particle_sidecar.owning_rank[particle_index] != world_rank) {
        continue;
      }
      context.state.particles.position_x_comoving[particle_index] +=
          context.state.particles.velocity_x_peculiar[particle_index] * drift_factor;
      context.state.particles.position_y_comoving[particle_index] +=
          context.state.particles.velocity_y_peculiar[particle_index] * drift_factor;
      context.state.particles.position_z_comoving[particle_index] +=
          context.state.particles.velocity_z_peculiar[particle_index] * drift_factor;
    }

    context.state.requireGasCellIdentityMapCoversDenseRows("drift callback");
    // Hydro cell centers are Eulerian geometry in both supported production
    // paths: fixed Cartesian H1 patches derive topology from these centers and
    // AMR patches derive it from explicit patch geometry.  A dense gas-cell row
    // and an optional parent particle are storage/compatibility relations, not
    // permission to advect physical cell centers during particle drift.  Moving
    // mesh geometry is not implemented by this callback; attempting to mirror
    // parent-particle positions here corrupts fixed Cartesian topology and makes
    // a valid row-order-independent patch appear nonuniform.

  }

 private:
  static constexpr std::array<core::StageContract, 1> m_contracts{{
      {.stage = core::IntegrationStage::kDrift,
       .required_inputs = core::StageDataDomain::kParticles,
       .mutated_state = core::StageDataDomain::kParticles | core::StageDataDomain::kGasCells,
       .produced_outputs = core::StageDataDomain::kParticles | core::StageDataDomain::kGasCells,
       .allowed_side_effects = core::StageDataDomain::kNone,
       .sync_requirements = core::StageSyncRequirement::kLocalOnly,
       .active_set_family = core::StageActiveSetFamily::kActiveParticles,
       .restart_safety = core::StageSafety::kUnsafe,
       .output_safety = core::StageSafety::kUnsafe,
       .owner = core::StageSubsystem::kCore},
  }};
};

class GravityStageCallback final : public core::IntegrationCallback {
 public:
  struct GravityHealthSummary {
    std::uint64_t cheap_checks_executed = 0;
    std::uint64_t heavy_checks_executed = 0;
    std::uint64_t warning_count = 0;
    std::uint64_t fatal_count = 0;
    std::uint64_t pm_field_non_finite_count = 0;
    std::uint64_t force_non_finite_count = 0;
    std::uint64_t force_abnormal_count = 0;
    std::uint64_t illegal_sync_state_count = 0;
    std::uint64_t zoom_sanity_failure_count = 0;
    std::uint64_t decomposition_sanity_failure_count = 0;
  };

  enum class PmSyncSurface : std::uint8_t {
    kKickPre = 0,
    kForceRefresh = 1,
    kKickPost = 2,
  };

  struct PmLongRangeKickDecision {
    PmSyncSurface sync_surface = PmSyncSurface::kKickPre;
    std::uint64_t gravity_kick_opportunity = 0;
    bool refresh_long_range_field = false;
    std::uint64_t field_version = 0;
    std::uint64_t last_refresh_opportunity = 0;
    std::uint64_t field_built_step_index = 0;
    double field_built_scale_factor = 1.0;
  };

  GravityStageCallback(
      const core::SimulationConfig& config,
      const core::ModePolicy& mode_policy,
      std::filesystem::path zoom_region_path = {})
      : m_config(config),
        m_mode_policy(mode_policy),
        m_pm_update_cadence_steps(static_cast<std::uint64_t>(config.numerics.treepm_update_cadence_steps)),
        m_pm_grid_shape(gravity::PmGridShape{
            static_cast<std::size_t>(config.numerics.treepm_pm_grid_nx),
            static_cast<std::size_t>(config.numerics.treepm_pm_grid_ny),
            static_cast<std::size_t>(config.numerics.treepm_pm_grid_nz)}),
        m_mesh_spacing_x_mpc_comoving(
            config.cosmology.box_size_x_mpc_comoving / static_cast<double>(m_pm_grid_shape.nx)),
        m_mesh_spacing_y_mpc_comoving(
            config.cosmology.box_size_y_mpc_comoving / static_cast<double>(m_pm_grid_shape.ny)),
        m_mesh_spacing_z_mpc_comoving(
            config.cosmology.box_size_z_mpc_comoving / static_cast<double>(m_pm_grid_shape.nz)),
        m_tree_pm_coordinator(makeRuntimeAwareTreePmCoordinator(config, m_pm_grid_shape)),
        m_zoom_region_path(std::move(zoom_region_path)) {
    if (!m_pm_grid_shape.isValid()) {
      throw std::runtime_error("numerics.treepm_pm_grid_n{xyz} must all be > 0");
    }
    m_tree_pm_options.pm_options.box_size_x_mpc_comoving = config.cosmology.box_size_x_mpc_comoving;
    m_tree_pm_options.pm_options.box_size_y_mpc_comoving = config.cosmology.box_size_y_mpc_comoving;
    m_tree_pm_options.pm_options.box_size_z_mpc_comoving = config.cosmology.box_size_z_mpc_comoving;
    m_tree_pm_options.pm_options.box_size_mpc_comoving = config.cosmology.box_size_x_mpc_comoving;
    m_tree_pm_options.pm_options.scale_factor = 1.0;
    m_tree_pm_options.pm_options.gravitational_constant_code = 1.0;
    m_tree_pm_options.pm_options.assignment_scheme =
        toPmAssignmentScheme(config.numerics.treepm_assignment_scheme);
    m_tree_pm_options.pm_options.enable_window_deconvolution =
        config.numerics.treepm_enable_window_deconvolution;
    m_tree_pm_options.pm_options.decomposition_mode = config.numerics.treepm_pm_decomposition_mode;
    m_tree_pm_options.pm_options.boundary_condition = mode_policy.gravity_boundary ==
            core::GravityBoundaryModel::kPeriodicPoisson
        ? gravity::PmBoundaryCondition::kPeriodic
        : gravity::PmBoundaryCondition::kIsolatedOpen;
    parallel::MpiContext mpi_context;
    const core::CudaRuntimeInfo cuda_runtime = core::queryCudaRuntime();
    m_runtime_topology = parallel::buildDistributedExecutionTopology(
        m_pm_grid_shape.nx,
        m_pm_grid_shape.ny,
        m_pm_grid_shape.nz,
        mpi_context,
        config.parallel.mpi_ranks_expected,
        config.parallel.gpu_devices,
        cuda_runtime.runtime_available,
        cuda_runtime.visible_device_count,
        pmDecompositionModeName(config.numerics.treepm_pm_decomposition_mode));
    if (m_runtime_topology.usesCuda()) {
      core::setCudaDeviceOrThrow(m_runtime_topology.device_assignment.assigned_device_index);
      m_tree_pm_options.pm_options.execution_policy = core::ExecutionPolicy::kCuda;
      m_tree_pm_options.pm_options.data_residency = gravity::PmDataResidencyPolicy::kPreferDevice;
    }

    m_tree_pm_options.tree_options.opening_theta = 0.7;
    m_tree_pm_options.tree_options.gravitational_constant_code = 1.0;
    m_tree_pm_options.tree_options.softening.kernel = gravity::TreeSofteningKernel::kPlummer;
    m_tree_pm_options.tree_options.softening.epsilon_comoving =
        config.numerics.gravity_softening_kpc_comoving * 1.0e-3;
    m_tree_pm_species_softening = speciesSofteningByTag(config);
    m_tree_pm_options.split_policy = gravity::makeTreePmSplitPolicyFromMeshSpacing(
        config.numerics.treepm_asmth_cells,
        config.numerics.treepm_rcut_cells,
        std::cbrt(
            m_mesh_spacing_x_mpc_comoving * m_mesh_spacing_y_mpc_comoving *
            m_mesh_spacing_z_mpc_comoving));
    m_tree_pm_options.enable_zoom_long_range_correction =
        config.mode.zoom_long_range_strategy ==
        core::ZoomLongRangeStrategy::kGlobalCoarsePlusFocusedHighResCorrection;
    m_tree_pm_options.zoom_focused_pm_shape = gravity::PmGridShape{
        static_cast<std::size_t>(std::max(config.mode.zoom_focused_pm_grid_nx, 0)),
        static_cast<std::size_t>(std::max(config.mode.zoom_focused_pm_grid_ny, 0)),
        static_cast<std::size_t>(std::max(config.mode.zoom_focused_pm_grid_nz, 0))};
    m_tree_pm_options.zoom_region_center_x_comoving = config.mode.zoom_region_center_x_mpc_comoving;
    m_tree_pm_options.zoom_region_center_y_comoving = config.mode.zoom_region_center_y_mpc_comoving;
    m_tree_pm_options.zoom_region_center_z_comoving = config.mode.zoom_region_center_z_mpc_comoving;
    m_tree_pm_options.zoom_region_radius_comoving = config.mode.zoom_region_radius_mpc_comoving;
    m_tree_pm_options.zoom_contamination_radius_comoving =
        config.mode.zoom_contamination_radius_mpc_comoving > 0.0
        ? config.mode.zoom_contamination_radius_mpc_comoving
        : config.mode.zoom_region_radius_mpc_comoving;
    m_tree_pm_options.tree_exchange_batch_bytes = config.numerics.treepm_tree_exchange_batch_bytes;
    m_pm_assignment_scheme = treePmAssignmentSchemeName(config.numerics.treepm_assignment_scheme);
    m_pm_backend = gravity::PmSolver::fftBackendName();
    if (m_runtime_topology.usesCuda()) {
      m_pm_backend += "+cuda_cic";
    }
  }

  [[nodiscard]] std::string_view callbackName() const override { return "gravity"; }
  [[nodiscard]] std::span<const core::IntegrationStage> integrationStages() const override { return k_gravity_stages; }
  [[nodiscard]] std::span<const core::StageContract> stageContracts() const override { return m_contracts; }
  [[nodiscard]] std::size_t pmGridSize() const noexcept { return m_pm_grid_shape.nx; }
  [[nodiscard]] const gravity::PmGridShape& pmGridShape() const noexcept { return m_pm_grid_shape; }
  [[nodiscard]] std::uint64_t longRangeRefreshCount() const noexcept { return m_long_range_refresh_count; }
  [[nodiscard]] std::uint64_t longRangeReuseCount() const noexcept { return m_long_range_reuse_count; }
  [[nodiscard]] int pmCadenceSteps() const noexcept { return static_cast<int>(m_pm_update_cadence_steps); }
  [[nodiscard]] std::span<const ReferenceWorkflowReport::TreePmCadenceRecord> cadenceRecords() const noexcept {
    return m_cadence_records;
  }
  [[nodiscard]] const parallel::DistributedExecutionTopology& runtimeTopology() const noexcept { return m_runtime_topology; }
  [[nodiscard]] core::MemoryReport memoryReport() const { return m_tree_pm_coordinator.memoryReport(); }
  [[nodiscard]] const parallel::DecompositionRuntimeMeasurements& lastRuntimeDecompositionMeasurements() const noexcept {
    return m_last_decomposition_measurements;
  }

  [[nodiscard]] std::span<const double> cellAccelX() const noexcept { return m_cell_accel_x; }
  [[nodiscard]] std::span<const double> cellAccelY() const noexcept { return m_cell_accel_y; }
  [[nodiscard]] std::span<const double> cellAccelZ() const noexcept { return m_cell_accel_z; }
  [[nodiscard]] std::span<const double> particleAccelX() const noexcept { return m_particle_accel_x; }
  [[nodiscard]] std::span<const double> particleAccelY() const noexcept { return m_particle_accel_y; }
  [[nodiscard]] std::span<const double> particleAccelZ() const noexcept { return m_particle_accel_z; }

  [[nodiscard]] io::GravityForceCachePersistentState exportRestartForceCache(
      const core::SimulationState& state) const {
    io::GravityForceCachePersistentState cache;
    cache.valid = m_force_cache_valid;
    if (!cache.valid) {
      return cache;
    }
    if (m_particle_accel_x.size() != state.particles.size() ||
        m_particle_accel_y.size() != state.particles.size() ||
        m_particle_accel_z.size() != state.particles.size() ||
        m_cell_accel_x.size() != state.cells.size() ||
        m_cell_accel_y.size() != state.cells.size() ||
        m_cell_accel_z.size() != state.cells.size()) {
      throw std::runtime_error(
          "gravity force cache cannot be checkpointed because its dense lanes do not match SimulationState");
    }
    cache.particle_id.assign(state.particle_sidecar.particle_id.begin(), state.particle_sidecar.particle_id.end());
    cache.gas_cell_id.assign(state.gas_cells.gas_cell_id.begin(), state.gas_cells.gas_cell_id.end());
    cache.particle_accel_x_comoving = m_particle_accel_x;
    cache.particle_accel_y_comoving = m_particle_accel_y;
    cache.particle_accel_z_comoving = m_particle_accel_z;
    cache.cell_accel_x_comoving = m_cell_accel_x;
    cache.cell_accel_y_comoving = m_cell_accel_y;
    cache.cell_accel_z_comoving = m_cell_accel_z;
    return cache;
  }

  void importRestartForceCache(
      const io::GravityForceCachePersistentState& cache,
      const core::SimulationState& state) {
    if (!cache.valid) {
      m_particle_accel_x.clear();
      m_particle_accel_y.clear();
      m_particle_accel_z.clear();
      m_cell_accel_x.clear();
      m_cell_accel_y.clear();
      m_cell_accel_z.clear();
      m_force_cache_valid = false;
      return;
    }
    const auto requireTriplet = [](const std::vector<double>& x,
                                   const std::vector<double>& y,
                                   const std::vector<double>& z,
                                   std::size_t expected,
                                   std::string_view label) {
      if (x.size() != expected || y.size() != expected || z.size() != expected) {
        throw std::runtime_error(
            "restart gravity force cache " + std::string(label) +
            " extent does not match authoritative SimulationState");
      }
      for (const double value : x) {
        if (!std::isfinite(value)) {
          throw std::runtime_error("restart gravity force cache contains non-finite values");
        }
      }
      for (const double value : y) {
        if (!std::isfinite(value)) {
          throw std::runtime_error("restart gravity force cache contains non-finite values");
        }
      }
      for (const double value : z) {
        if (!std::isfinite(value)) {
          throw std::runtime_error("restart gravity force cache contains non-finite values");
        }
      }
    };
    const auto buildIndexById = [](std::span<const std::uint64_t> ids,
                                   std::size_t expected,
                                   std::string_view label) {
      if (ids.size() != expected) {
        throw std::runtime_error("restart gravity force cache " + std::string(label) +
                                 " identity extent does not match acceleration lanes");
      }
      std::unordered_map<std::uint64_t, std::size_t> index_by_id;
      index_by_id.reserve(ids.size());
      for (std::size_t index = 0; index < ids.size(); ++index) {
        if (ids[index] == 0U || !index_by_id.emplace(ids[index], index).second) {
          throw std::runtime_error("restart gravity force cache " + std::string(label) +
                                   " identity lane contains a zero or duplicate stable ID");
        }
      }
      return index_by_id;
    };
    requireTriplet(
        cache.particle_accel_x_comoving, cache.particle_accel_y_comoving,
        cache.particle_accel_z_comoving, state.particles.size(), "particle");
    requireTriplet(
        cache.cell_accel_x_comoving, cache.cell_accel_y_comoving,
        cache.cell_accel_z_comoving, state.cells.size(), "gas-cell");
    const auto particle_index_by_id = buildIndexById(
        cache.particle_id, state.particles.size(), "particle");
    const auto gas_index_by_id = buildIndexById(
        cache.gas_cell_id, state.cells.size(), "gas-cell");

    m_particle_accel_x.resize(state.particles.size());
    m_particle_accel_y.resize(state.particles.size());
    m_particle_accel_z.resize(state.particles.size());
    for (std::size_t row = 0; row < state.particles.size(); ++row) {
      const auto it = particle_index_by_id.find(state.particle_sidecar.particle_id[row]);
      if (it == particle_index_by_id.end()) {
        throw std::runtime_error("restart gravity force cache is missing a persistent particle ID");
      }
      const std::size_t cached_row = it->second;
      m_particle_accel_x[row] = cache.particle_accel_x_comoving[cached_row];
      m_particle_accel_y[row] = cache.particle_accel_y_comoving[cached_row];
      m_particle_accel_z[row] = cache.particle_accel_z_comoving[cached_row];
    }
    m_cell_accel_x.resize(state.cells.size());
    m_cell_accel_y.resize(state.cells.size());
    m_cell_accel_z.resize(state.cells.size());
    for (std::size_t row = 0; row < state.cells.size(); ++row) {
      const auto it = gas_index_by_id.find(state.gas_cells.gas_cell_id[row]);
      if (it == gas_index_by_id.end()) {
        throw std::runtime_error("restart gravity force cache is missing a persistent gas-cell ID");
      }
      const std::size_t cached_row = it->second;
      m_cell_accel_x[row] = cache.cell_accel_x_comoving[cached_row];
      m_cell_accel_y[row] = cache.cell_accel_y_comoving[cached_row];
      m_cell_accel_z[row] = cache.cell_accel_z_comoving[cached_row];
    }
    m_force_cache_valid = true;
  }

  [[nodiscard]] static PmSyncSurface toPmSyncSurface(core::IntegrationStage stage) {
    if (stage == core::IntegrationStage::kGravityKickPre) {
      return PmSyncSurface::kKickPre;
    }
    if (stage == core::IntegrationStage::kForceRefresh) {
      return PmSyncSurface::kForceRefresh;
    }
    if (stage == core::IntegrationStage::kGravityKickPost) {
      return PmSyncSurface::kKickPost;
    }
    throw std::invalid_argument("TreePM sync surface is only defined on legal gravity force-refresh/kick stages");
  }

  [[nodiscard]] static std::string_view pmSyncSurfaceName(PmSyncSurface surface) {
    switch (surface) {
      case PmSyncSurface::kKickPre:
        return "kick_pre";
      case PmSyncSurface::kForceRefresh:
        return "force_refresh";
      case PmSyncSurface::kKickPost:
        return "kick_post";
    }
    return "unknown";
  }
  [[nodiscard]] static std::string_view pmRefreshReasonName(core::PmRefreshDirective::Reason reason) {
    switch (reason) {
      case core::PmRefreshDirective::Reason::kNone:
        return "none";
      case core::PmRefreshDirective::Reason::kInitialForceBootstrap:
        return "initial_force_bootstrap";
      case core::PmRefreshDirective::Reason::kScheduledForceRefreshStage:
        return "scheduled_force_refresh_stage";
    }
    return "unknown";
  }

  void onStage(core::StepContext& context) override {
    const bool is_kick_stage = context.stage == core::IntegrationStage::kGravityKickPre ||
        context.stage == core::IntegrationStage::kGravityKickPost;
    const bool needs_initial_force_bootstrap = context.stage == core::IntegrationStage::kGravityKickPre && !m_force_cache_valid;
    if (needs_initial_force_bootstrap && !context.pm_refresh_directive.initial_cache_bootstrap_allowed) {
      throw std::runtime_error(
          "TreePM initial force bootstrap requested outside an integrator-authorized global boundary");
    }
    const bool is_force_refresh_stage = context.pm_refresh_directive.force_refresh_surface || needs_initial_force_bootstrap;
    if (!is_kick_stage && !is_force_refresh_stage) {
      throw std::logic_error("gravity handler received an unregistered stage");
    }
    if (context.pm_refresh_directive.has_sync_event &&
        context.pm_refresh_directive.refresh_long_range_field &&
        !context.boundary.pm_refresh_allowed) {
      throw std::runtime_error("TreePM long-range PM refresh reached an illegal integration boundary");
    }
    if (context.stage == core::IntegrationStage::kForceRefresh && !context.pm_refresh_directive.force_refresh_surface) {
      throw std::runtime_error("TreePM force-refresh stage lacks an integrator-issued PM refresh directive");
    }
    if (is_kick_stage && !is_force_refresh_stage) {
      applyCachedKick(context);
      return;
    }

    parallel::MpiContext mpi_context;
    const std::uint64_t world_size = static_cast<std::uint64_t>(mpi_context.worldSize());
    if (m_runtime_topology.world_size != mpi_context.worldSize()) {
      throw std::runtime_error("reference workflow runtime topology world_size drifted from MPI context");
    }

    const auto requireKickConsensus = [&](std::uint64_t local_value, std::string_view name) {
      const std::uint64_t global_sum = mpi_context.allreduceSumUint64(local_value);
      if (global_sum != local_value * world_size) {
        throw std::runtime_error(
            "TreePM cadence rank-consensus failure for " + std::string(name) +
            ": local=" + std::to_string(local_value) + ", reduced_sum=" + std::to_string(global_sum) +
            ", world_size=" + std::to_string(world_size));
      }
    };

    const std::size_t particle_count = context.state.particles.size();

    const SolverGhostRefreshReport gravity_ghost_refresh = refreshParticleGhostsForSolver(
        context, mpi_context, "gravity.treepm", &m_ghost_cache_lifecycle);

    rebuildOwnedParticleCompactView(
        context,
        context.active_set.particle_indices,
        context.pm_refresh_directive.requires_predicted_inactive_sources);
    m_active_accel_x.assign(m_local_active_indices.size(), 0.0);
    m_active_accel_y.assign(m_local_active_indices.size(), 0.0);
    m_active_accel_z.assign(m_local_active_indices.size(), 0.0);
    m_active_is_high_res.assign(m_local_active_indices.size(), 0U);
    m_active_slot_by_particle.assign(particle_count, -1);
    for (std::size_t i = 0; i < m_local_active_global_indices.size(); ++i) {
      m_active_slot_by_particle[m_local_active_global_indices[i]] = static_cast<int>(i);
    }
    ensureZoomMembershipLoaded(context.state);
    const double box_size_x = m_config.cosmology.box_size_x_mpc_comoving;
    const double box_size_y = m_config.cosmology.box_size_y_mpc_comoving;
    const double box_size_z = m_config.cosmology.box_size_z_mpc_comoving;
    m_source_is_high_res.assign(m_local_source_x.size(), 0U);
    for (std::size_t i = 0; i < m_local_source_x.size(); ++i) {
      if (!m_zoom_high_res_particle_ids.empty()) {
        const std::uint32_t global_index = m_local_to_global[i];
        const std::uint64_t particle_id = context.state.particle_sidecar.particle_id[global_index];
        m_source_is_high_res[i] = m_zoom_high_res_particle_ids.contains(particle_id) ? 1U : 0U;
        continue;
      }
      const double dx = m_local_source_x[i] - m_tree_pm_options.zoom_region_center_x_comoving;
      const double dy = m_local_source_y[i] - m_tree_pm_options.zoom_region_center_y_comoving;
      const double dz = m_local_source_z[i] - m_tree_pm_options.zoom_region_center_z_comoving;
      const double wrapped_dx = dx - box_size_x * std::nearbyint(dx / box_size_x);
      const double wrapped_dy = dy - box_size_y * std::nearbyint(dy / box_size_y);
      const double wrapped_dz = dz - box_size_z * std::nearbyint(dz / box_size_z);
      const double r = std::sqrt(wrapped_dx * wrapped_dx + wrapped_dy * wrapped_dy + wrapped_dz * wrapped_dz);
      m_source_is_high_res[i] = (r <= m_tree_pm_options.zoom_region_radius_comoving) ? 1U : 0U;
    }
    for (std::size_t i = 0; i < m_local_active_indices.size(); ++i) {
      m_active_is_high_res[i] = m_source_is_high_res[m_local_active_indices[i]];
    }
    m_tree_pm_options.source_is_high_res = m_source_is_high_res;
    m_tree_pm_options.active_is_high_res = m_active_is_high_res;

    gravity::TreePmForceAccumulatorView accumulator{
        .active_particle_index = m_local_active_indices,
        .accel_x_comoving = m_active_accel_x,
        .accel_y_comoving = m_active_accel_y,
        .accel_z_comoving = m_active_accel_z,
    };

    m_tree_pm_options.pm_options.scale_factor = std::max(1.0e-12, context.integrator_state.current_scale_factor);
    if (!m_mode_policy.cosmological_comoving_frame) {
      m_tree_pm_options.pm_options.scale_factor = 1.0;
    }

    // The integrator owns PM cadence state.  At legal global PM-refresh
    // boundaries it issues a concrete sync event; local active-bin force
    // evaluations may recompute short-range forces only and must reuse the
    // already-valid long-range field without advancing PM cadence truth.
    PmLongRangeKickDecision decision{
        .sync_surface = toPmSyncSurface(context.stage),
        .gravity_kick_opportunity = context.integrator_state.pm_sync_state.gravityKickOpportunity(),
        .refresh_long_range_field = false,
        .field_version = context.integrator_state.pm_sync_state.fieldVersion(),
        .last_refresh_opportunity = context.integrator_state.pm_sync_state.lastRefreshOpportunity(),
        .field_built_step_index = context.integrator_state.pm_sync_state.lastRefreshStepIndex(),
        .field_built_scale_factor = context.integrator_state.pm_sync_state.lastRefreshScaleFactor(),
    };

    if (context.pm_refresh_directive.has_sync_event) {
      if (!context.pm_refresh_directive.cadence_opportunity_allowed) {
        throw std::runtime_error(
            "TreePM received a PM sync event outside an integrator-authorized refresh opportunity");
      }
      const core::PmSyncEvent sync_event{
          .gravity_kick_opportunity = context.pm_refresh_directive.gravity_kick_opportunity,
          .refresh_long_range_field = context.pm_refresh_directive.refresh_long_range_field,
          .field_version = context.pm_refresh_directive.field_version,
          .last_refresh_opportunity = context.pm_refresh_directive.last_refresh_opportunity,
          .field_built_step_index = context.pm_refresh_directive.field_built_step_index,
          .field_built_scale_factor = context.pm_refresh_directive.field_built_scale_factor,
      };
      requireKickConsensus(sync_event.gravity_kick_opportunity, "gravity_kick_opportunity");
      const std::uint64_t refresh_vote = sync_event.refresh_long_range_field ? 1ULL : 0ULL;
      const std::uint64_t refresh_vote_sum = mpi_context.allreduceSumUint64(refresh_vote);
      if (refresh_vote_sum != 0ULL && refresh_vote_sum != world_size) {
        throw std::runtime_error(
            "TreePM long-range cadence decision diverged across ranks: refresh_vote_sum=" +
            std::to_string(refresh_vote_sum) + ", world_size=" + std::to_string(world_size));
      }
      const bool refresh_long_range = (refresh_vote_sum == world_size);
      if (refresh_long_range != sync_event.refresh_long_range_field) {
        throw std::runtime_error(
            "TreePM long-range cadence local decision drifted from rank consensus");
      }

      decision.gravity_kick_opportunity = sync_event.gravity_kick_opportunity;
      decision.refresh_long_range_field = refresh_long_range;
      decision.field_version = sync_event.field_version;
      decision.last_refresh_opportunity = sync_event.last_refresh_opportunity;
      decision.field_built_step_index = sync_event.field_built_step_index;
      decision.field_built_scale_factor = sync_event.field_built_scale_factor;
    } else {
      if (!context.integrator_state.pm_long_range_field_valid || decision.field_version == 0U) {
        throw std::runtime_error(
            "TreePM local force refresh attempted before a valid long-range PM field existed");
      }
      if (context.pm_refresh_directive.cadence_opportunity_allowed) {
        throw std::runtime_error(
            "TreePM cadence opportunity was declared without an accompanying integrator-owned PM sync event");
      }
    }

    requireKickConsensus(decision.field_version, "long_range_field_version");
    requireKickConsensus(decision.last_refresh_opportunity, "last_long_range_refresh_opportunity");

    if (m_particle_accel_x.size() != particle_count) {
      m_particle_accel_x.assign(particle_count, 0.0);
      m_particle_accel_y.assign(particle_count, 0.0);
      m_particle_accel_z.assign(particle_count, 0.0);
    }
    if (m_cell_accel_x.size() != context.state.cells.size()) {
      m_cell_accel_x.assign(context.state.cells.size(), 0.0);
      m_cell_accel_y.assign(context.state.cells.size(), 0.0);
      m_cell_accel_z.assign(context.state.cells.size(), 0.0);
    }
    const gravity::TreeSofteningView softening_view{
        .source_species_tag = std::span<const std::uint32_t>(m_local_source_species_tag.data(), m_local_source_species_tag.size()),
        .source_particle_epsilon_comoving = std::span<const double>(
            m_local_source_softening_comoving.empty() ? nullptr : m_local_source_softening_comoving.data(),
            m_local_source_softening_comoving.size()),
        .source_particle_epsilon_override_mask = std::span<const std::uint8_t>(
            m_local_source_softening_override_mask.empty() ? nullptr : m_local_source_softening_override_mask.data(),
            m_local_source_softening_override_mask.size()),
        .target_species_tag = std::span<const std::uint32_t>(m_active_target_species_tag.data(), m_active_target_species_tag.size()),
        .target_particle_epsilon_comoving = std::span<const double>(
            m_active_target_softening_comoving.empty() ? nullptr : m_active_target_softening_comoving.data(),
            m_active_target_softening_comoving.size()),
        .target_particle_epsilon_override_mask = std::span<const std::uint8_t>(
            m_active_target_softening_override_mask.empty() ? nullptr : m_active_target_softening_override_mask.data(),
            m_active_target_softening_override_mask.size()),
        .species_policy = m_tree_pm_species_softening,
    };
    gravity::TreePmProfileEvent tree_pm_profile;
    m_tree_pm_coordinator.solveActiveSetWithPmCadence(
        m_local_source_x,
        m_local_source_y,
        m_local_source_z,
        m_local_source_mass,
        accumulator,
        m_tree_pm_options,
        decision.refresh_long_range_field,
        &tree_pm_profile,
        &m_last_tree_pm_diagnostics,
        softening_view);

    m_last_decomposition_measurements = parallel::DecompositionRuntimeMeasurements{
        .tree_pair_evaluations_recent = m_last_tree_pm_diagnostics.residual_pair_evaluations +
            tree_pm_profile.tree_profile.particle_particle_interactions,
        .tree_remote_request_bytes_recent = m_last_tree_pm_diagnostics.residual_remote_request_bytes +
            m_last_tree_pm_diagnostics.residual_remote_response_bytes,
        .pm_mesh_cells_touched_recent = static_cast<std::uint64_t>(m_tree_pm_coordinator.slabLayout().localCellCount()),
        .pm_fft_transpose_bytes_recent = tree_pm_profile.pm_profile.fft_transpose_bytes,
        .amr_patch_cells_updated_recent = static_cast<std::uint64_t>(context.state.cells.size()),
        .hydro_face_fluxes_recent = 0,
        .ghost_exchange_bytes_recent = gravity_ghost_refresh.sent_bytes + gravity_ghost_refresh.received_bytes,
        .tree_wall_ms_recent = tree_pm_profile.tree_profile.build_ms + tree_pm_profile.tree_profile.multipole_ms +
            tree_pm_profile.tree_profile.traversal_ms + tree_pm_profile.tree_short_range_ms,
        .pm_wall_ms_recent = tree_pm_profile.pm_profile.assign_ms + tree_pm_profile.pm_profile.fft_forward_ms +
            tree_pm_profile.pm_profile.poisson_ms + tree_pm_profile.pm_profile.gradient_ms +
            tree_pm_profile.pm_profile.fft_inverse_ms + tree_pm_profile.pm_profile.fft_transpose_ms +
            tree_pm_profile.pm_profile.interpolate_ms,
        .gpu_kernel_ms_recent = tree_pm_profile.pm_profile.device_kernel_ms,
        .accelerator_occupancy_fraction_recent = (tree_pm_profile.pm_profile.device_kernel_ms > 0.0) ? 1.0 : 0.0,
        .has_measurements = true,
    };

    const bool allow_heavy_reference_checks =
        m_config.analysis.diagnostics_execution_policy ==
        core::AnalysisConfig::DiagnosticsExecutionPolicy::kAllIncludingProvisional;
    GravityHealthSummary gravity_health = validateGravityHealth(
        context,
        decision,
        allow_heavy_reference_checks);

    if (decision.refresh_long_range_field) {
      m_has_long_range_field = true;
      ++m_long_range_refresh_count;
    } else {
      ++m_long_range_reuse_count;
    }
    context.pm_refresh_directive.solver_executed = true;
    m_last_committed_field_version = decision.field_version;

    const std::uint64_t inactive_particles_skipped = static_cast<std::uint64_t>(
        context.state.particles.size() - m_local_active_global_indices.size());
    m_cadence_records.push_back(ReferenceWorkflowReport::TreePmCadenceRecord{
        .step_index = context.integrator_state.step_index,
        .stage_name = std::string(core::integrationStageName(context.stage)),
        .pm_sync_surface = std::string(pmSyncSurfaceName(decision.sync_surface)),
        .gravity_kick_opportunity = decision.gravity_kick_opportunity,
        .field_version = decision.field_version,
        .last_refresh_opportunity = decision.last_refresh_opportunity,
        .field_built_step_index = decision.field_built_step_index,
        .field_built_scale_factor = decision.field_built_scale_factor,
        .field_age_in_kick_opportunities =
            decision.gravity_kick_opportunity - decision.last_refresh_opportunity,
        .active_particles_kicked = static_cast<std::uint64_t>(m_local_active_global_indices.size()),
        .inactive_particles_skipped = inactive_particles_skipped,
        .refreshed_long_range_field = decision.refresh_long_range_field,
    });
    if (context.profiler_session != nullptr) {
      context.profiler_session->recordEvent(core::RuntimeEvent{
          .event_kind = "gravity.pm_long_range_field",
          .severity = core::RuntimeEventSeverity::kInfo,
          .subsystem = "gravity.treepm",
          .step_index = context.integrator_state.step_index,
          .simulation_time_code = context.integrator_state.current_time_code,
          .scale_factor = context.integrator_state.current_scale_factor,
          .message = decision.refresh_long_range_field
              ? "PM long-range field refreshed for gravity kick"
              : "PM long-range field reused for gravity kick",
          .payload = {{"stage", std::string(core::integrationStageName(context.stage))},
                      {"pm_sync_surface", std::string(pmSyncSurfaceName(decision.sync_surface))},
                      {"gravity_kick_opportunity", std::to_string(context.integrator_state.pm_sync_state.gravityKickOpportunity())},
                      {"field_version", std::to_string(context.integrator_state.pm_sync_state.fieldVersion())},
                      {"field_built_step_index", std::to_string(context.integrator_state.pm_sync_state.lastRefreshStepIndex())},
                      {"field_built_scale_factor", std::to_string(context.integrator_state.pm_sync_state.lastRefreshScaleFactor())},
                      {"pm_update_cadence_steps", std::to_string(context.integrator_state.pm_sync_state.cadenceSteps())},
                      {"pm_grid", std::to_string(m_pm_grid_shape.nx) + "x" + std::to_string(m_pm_grid_shape.ny) +
                              "x" + std::to_string(m_pm_grid_shape.nz)},
                      {"pm_assignment_scheme", m_pm_assignment_scheme},
                      {"pm_window_deconvolution", m_config.numerics.treepm_enable_window_deconvolution ? "true" : "false"},
                      {"asmth_cells", std::to_string(m_config.numerics.treepm_asmth_cells)},
                      {"rcut_cells", std::to_string(m_config.numerics.treepm_rcut_cells)},
                      {"mesh_spacing_x_mpc_comoving", std::to_string(m_mesh_spacing_x_mpc_comoving)},
                      {"mesh_spacing_y_mpc_comoving", std::to_string(m_mesh_spacing_y_mpc_comoving)},
                      {"mesh_spacing_z_mpc_comoving", std::to_string(m_mesh_spacing_z_mpc_comoving)},
                      {"split_scale_mpc_comoving", std::to_string(m_tree_pm_options.split_policy.split_scale_comoving)},
                      {"cutoff_radius_mpc_comoving", std::to_string(m_tree_pm_options.split_policy.cutoff_radius_comoving)},
                      {"softening_policy", describeSofteningPolicy(m_config, &context.state)},
                      {"softening_kernel", "plummer"},
                      {"softening_epsilon_kpc_comoving", std::to_string(m_config.numerics.gravity_softening_kpc_comoving)},
                      {"pm_fft_backend", m_pm_backend},
                      {"active_particles_kicked", std::to_string(m_local_active_global_indices.size())},
                      {"inactive_particles_skipped", std::to_string(inactive_particles_skipped)},
                      {"ghost_refresh_sent_bytes", std::to_string(gravity_ghost_refresh.sent_bytes)},
                      {"ghost_refresh_received_bytes", std::to_string(gravity_ghost_refresh.received_bytes)},
                      {"ghost_refresh_committed_slots", std::to_string(gravity_ghost_refresh.committed_slots)},
                      {"predicted_inactive_source_particles", std::to_string(m_source_predicted_inactive_count)},
                      {"predicted_inactive_sources_required",
                          context.pm_refresh_directive.requires_predicted_inactive_sources ? "true" : "false"},
                      {"pm_refresh_reason", std::string(pmRefreshReasonName(context.pm_refresh_directive.reason))},
                      {"pm_refresh_force_eval_scale_factor", std::to_string(context.pm_refresh_directive.force_evaluation_scale_factor)},
                      {"refreshed_long_range_field", decision.refresh_long_range_field ? "true" : "false"}},
      });
      context.profiler_session->recordEvent(core::RuntimeEvent{
          .event_kind = "gravity.zoom_force_diagnostics",
          .severity = core::RuntimeEventSeverity::kInfo,
          .subsystem = "gravity.treepm",
          .step_index = context.integrator_state.step_index,
          .simulation_time_code = context.integrator_state.current_time_code,
          .scale_factor = context.integrator_state.current_scale_factor,
          .message = "zoom force decomposition and contamination diagnostics",
          .payload = {
              {"force_l2_pm_global", std::to_string(m_last_tree_pm_diagnostics.force_l2_pm_global)},
              {"force_l2_pm_zoom_correction", std::to_string(m_last_tree_pm_diagnostics.force_l2_pm_zoom_correction)},
              {"force_l2_tree_short_range", std::to_string(m_last_tree_pm_diagnostics.force_l2_tree_short_range)},
              {"force_l2_total", std::to_string(m_last_tree_pm_diagnostics.force_l2_total)},
              {"zoom_high_res_source_count", std::to_string(m_last_tree_pm_diagnostics.zoom_high_res_source_count)},
              {"zoom_low_res_source_count", std::to_string(m_last_tree_pm_diagnostics.zoom_low_res_source_count)},
              {"zoom_low_res_contamination_count", std::to_string(m_last_tree_pm_diagnostics.zoom_low_res_contamination_count)},
              {"zoom_low_res_contamination_mass_code", std::to_string(m_last_tree_pm_diagnostics.zoom_low_res_contamination_mass_code)},
              {"zoom_membership_source", m_zoom_membership_source},
          },
      });
      context.profiler_session->recordEvent(core::RuntimeEvent{
          .event_kind = "gravity.health_summary",
          .severity = (gravity_health.fatal_count > 0U)
              ? core::RuntimeEventSeverity::kFatal
              : (gravity_health.warning_count > 0U)
                    ? core::RuntimeEventSeverity::kWarning
                    : core::RuntimeEventSeverity::kInfo,
          .subsystem = "gravity.treepm",
          .step_index = context.integrator_state.step_index,
          .simulation_time_code = context.integrator_state.current_time_code,
          .scale_factor = context.integrator_state.current_scale_factor,
          .message = "gravity health summary across PM field, force, sync, zoom, and decomposition checks",
          .payload = {
              {"cheap_checks_executed", std::to_string(gravity_health.cheap_checks_executed)},
              {"heavy_checks_executed", std::to_string(gravity_health.heavy_checks_executed)},
              {"warning_count", std::to_string(gravity_health.warning_count)},
              {"fatal_count", std::to_string(gravity_health.fatal_count)},
              {"pm_field_non_finite_count", std::to_string(gravity_health.pm_field_non_finite_count)},
              {"force_non_finite_count", std::to_string(gravity_health.force_non_finite_count)},
              {"force_abnormal_count", std::to_string(gravity_health.force_abnormal_count)},
              {"illegal_sync_state_count", std::to_string(gravity_health.illegal_sync_state_count)},
              {"zoom_sanity_failure_count", std::to_string(gravity_health.zoom_sanity_failure_count)},
              {"decomposition_sanity_failure_count", std::to_string(gravity_health.decomposition_sanity_failure_count)},
              {"heavy_reference_checks_opt_in", allow_heavy_reference_checks ? "true" : "false"},
          },
      });
    }

    m_force_cache_valid = true;
    if (is_kick_stage) {
      applyActiveKickFromFreshForce(context);
    }

    for (std::size_t active_slot = 0; active_slot < m_local_active_global_indices.size(); ++active_slot) {
      const std::uint32_t particle_index = m_local_active_global_indices[active_slot];
      m_particle_accel_x[particle_index] = m_active_accel_x[active_slot];
      m_particle_accel_y[particle_index] = m_active_accel_y[active_slot];
      m_particle_accel_z[particle_index] = m_active_accel_z[active_slot];
    }

    context.state.requireGasCellIdentityMapCoversDenseRows("gravity callback gas-cell acceleration sync");
    const auto particle_row_by_id = buildParticleRowById(context.state);
    for (std::size_t cell_index = 0; cell_index < context.state.cells.size(); ++cell_index) {
      const std::optional<std::uint32_t> gas_particle_index = parentParticleRowForGasCellRow(
          context.state,
          static_cast<std::uint32_t>(cell_index),
          particle_row_by_id,
          "gravity callback gas-cell acceleration sync");
      if (!gas_particle_index.has_value()) {
        m_cell_accel_x[cell_index] = 0.0;
        m_cell_accel_y[cell_index] = 0.0;
        m_cell_accel_z[cell_index] = 0.0;
        continue;
      }
      const int active_slot = m_active_slot_by_particle[*gas_particle_index];
      if (active_slot < 0) {
        continue;
      }
      m_cell_accel_x[cell_index] = m_active_accel_x[static_cast<std::size_t>(active_slot)];
      m_cell_accel_y[cell_index] = m_active_accel_y[static_cast<std::size_t>(active_slot)];
      m_cell_accel_z[cell_index] = m_active_accel_z[static_cast<std::size_t>(active_slot)];
    }
  }

 private:
  static constexpr std::array<core::StageContract, 3> m_contracts{{
      {.stage = core::IntegrationStage::kGravityKickPre,
       .required_inputs = core::StageDataDomain::kParticles | core::StageDataDomain::kPmField,
       .mutated_state = core::StageDataDomain::kParticles | core::StageDataDomain::kDiagnostics,
       .produced_outputs = core::StageDataDomain::kParticles | core::StageDataDomain::kDiagnostics,
       .allowed_side_effects = core::StageDataDomain::kDiagnostics,
       .sync_requirements = core::StageSyncRequirement::kGlobal,
       .active_set_family = core::StageActiveSetFamily::kActiveParticles,
       .restart_safety = core::StageSafety::kUnsafe,
       .output_safety = core::StageSafety::kUnsafe,
       .owner = core::StageSubsystem::kGravity},
      {.stage = core::IntegrationStage::kForceRefresh,
       .required_inputs = core::StageDataDomain::kParticles | core::StageDataDomain::kGhostCells | core::StageDataDomain::kPmField | core::StageDataDomain::kTreeState,
       .mutated_state = core::StageDataDomain::kPmField | core::StageDataDomain::kTreeState | core::StageDataDomain::kDiagnostics,
       .produced_outputs = core::StageDataDomain::kPmField | core::StageDataDomain::kTreeState | core::StageDataDomain::kDiagnostics,
       .allowed_side_effects = core::StageDataDomain::kDiagnostics,
       .sync_requirements = core::StageSyncRequirement::kForceEvaluation,
       .active_set_family = core::StageActiveSetFamily::kActiveParticles,
       .restart_safety = core::StageSafety::kUnsafe,
       .output_safety = core::StageSafety::kUnsafe,
       .owner = core::StageSubsystem::kGravity},
      {.stage = core::IntegrationStage::kGravityKickPost,
       .required_inputs = core::StageDataDomain::kParticles | core::StageDataDomain::kPmField,
       .mutated_state = core::StageDataDomain::kParticles | core::StageDataDomain::kDiagnostics,
       .produced_outputs = core::StageDataDomain::kParticles | core::StageDataDomain::kDiagnostics,
       .allowed_side_effects = core::StageDataDomain::kDiagnostics,
       .sync_requirements = core::StageSyncRequirement::kLocalOnly,
       .active_set_family = core::StageActiveSetFamily::kActiveParticles,
       .restart_safety = core::StageSafety::kUnsafe,
       .output_safety = core::StageSafety::kUnsafe,
       .owner = core::StageSubsystem::kGravity},
  }};

  [[nodiscard]] static double kickFactorForStage(const core::StepContext& context) {
    if (context.stage == core::IntegrationStage::kGravityKickPre) {
      return context.timeline_step.first_kick_factor_code;
    }
    if (context.stage == core::IntegrationStage::kGravityKickPost) {
      return context.timeline_step.second_kick_factor_code;
    }
    throw std::invalid_argument("gravity kick factor requested outside KDK kick stage");
  }

  [[nodiscard]] static double hubbleDragFactorForStage(const core::StepContext& context) {
    if (context.stage == core::IntegrationStage::kGravityKickPre) {
      return context.timeline_step.first_hubble_drag_factor;
    }
    if (context.stage == core::IntegrationStage::kGravityKickPost) {
      return context.timeline_step.second_hubble_drag_factor;
    }
    throw std::invalid_argument("Hubble drag factor requested outside KDK kick stage");
  }

  static void applyPeculiarVelocityKick(
      core::SimulationState& state,
      std::uint32_t particle_index,
      double accel_x_code,
      double accel_y_code,
      double accel_z_code,
      double kick_factor_code,
      double hubble_drag_factor) {
    state.particles.velocity_x_peculiar[particle_index] =
        hubble_drag_factor * state.particles.velocity_x_peculiar[particle_index] + kick_factor_code * accel_x_code;
    state.particles.velocity_y_peculiar[particle_index] =
        hubble_drag_factor * state.particles.velocity_y_peculiar[particle_index] + kick_factor_code * accel_y_code;
    state.particles.velocity_z_peculiar[particle_index] =
        hubble_drag_factor * state.particles.velocity_z_peculiar[particle_index] + kick_factor_code * accel_z_code;
  }

  void applyActiveKickFromFreshForce(core::StepContext& context) {
    const double kick_factor = kickFactorForStage(context);
    const double hubble_drag = hubbleDragFactorForStage(context);
    for (std::size_t active_slot = 0; active_slot < m_local_active_global_indices.size(); ++active_slot) {
      const std::uint32_t particle_index = m_local_active_global_indices[active_slot];
      applyPeculiarVelocityKick(
          context.state,
          particle_index,
          m_active_accel_x[active_slot],
          m_active_accel_y[active_slot],
          m_active_accel_z[active_slot],
          kick_factor,
          hubble_drag);
    }
  }

  void applyCachedKick(core::StepContext& context) {
    if (!m_force_cache_valid) {
      throw std::runtime_error("gravity kick requested before a coherent force-refresh boundary");
    }
    const std::size_t particle_count = context.state.particles.size();
    if (m_particle_accel_x.size() != particle_count ||
        m_particle_accel_y.size() != particle_count ||
        m_particle_accel_z.size() != particle_count) {
      m_particle_accel_x.resize(particle_count, 0.0);
      m_particle_accel_y.resize(particle_count, 0.0);
      m_particle_accel_z.resize(particle_count, 0.0);
    }
    rebuildOwnedParticleCompactView(context, context.active_set.particle_indices, false);
    const double kick_factor = kickFactorForStage(context);
    const double hubble_drag = hubbleDragFactorForStage(context);
    for (const std::uint32_t particle_index : m_local_active_global_indices) {
      applyPeculiarVelocityKick(
          context.state,
          particle_index,
          m_particle_accel_x[particle_index],
          m_particle_accel_y[particle_index],
          m_particle_accel_z[particle_index],
          kick_factor,
          hubble_drag);
    }
  }

  void emitGravityEvent(
      core::StepContext& context,
      core::RuntimeEventSeverity severity,
      std::string message,
      std::unordered_map<std::string, std::string> payload = {}) const {
    if (context.profiler_session == nullptr) {
      return;
    }
    context.profiler_session->recordEvent(core::RuntimeEvent{
        .event_kind = "gravity.health_check",
        .severity = severity,
        .subsystem = "gravity.treepm",
        .step_index = context.integrator_state.step_index,
        .simulation_time_code = context.integrator_state.current_time_code,
        .scale_factor = context.integrator_state.current_scale_factor,
        .message = std::move(message),
        .payload = std::move(payload),
    });
  }

  [[nodiscard]] GravityHealthSummary validateGravityHealth(
      core::StepContext& context,
      const PmLongRangeKickDecision& decision,
      bool allow_heavy_reference_checks) {
    GravityHealthSummary summary;
    auto failFatal = [&](std::string message, std::unordered_map<std::string, std::string> payload = {}) {
      ++summary.fatal_count;
      emitGravityEvent(context, core::RuntimeEventSeverity::kFatal, message, payload);
      throw std::runtime_error("fatal gravity-state check failed: " + message);
    };
    auto addWarning = [&](std::string message, std::unordered_map<std::string, std::string> payload = {}) {
      ++summary.warning_count;
      emitGravityEvent(context, core::RuntimeEventSeverity::kWarning, message, payload);
    };

    ++summary.cheap_checks_executed;
    if (m_local_active_indices.size() != m_local_active_global_indices.size()) {
      ++summary.decomposition_sanity_failure_count;
      failFatal("local active index vectors diverged after compact-view rebuild");
    }
    if (m_local_active_indices.size() > m_local_source_x.size()) {
      ++summary.decomposition_sanity_failure_count;
      failFatal("active local subset exceeds local source population");
    }

    ++summary.cheap_checks_executed;
    if (!std::isfinite(m_last_tree_pm_diagnostics.force_l2_pm_global) ||
        !std::isfinite(m_last_tree_pm_diagnostics.force_l2_pm_zoom_correction) ||
        !std::isfinite(m_last_tree_pm_diagnostics.force_l2_tree_short_range) ||
        !std::isfinite(m_last_tree_pm_diagnostics.force_l2_total)) {
      ++summary.pm_field_non_finite_count;
      failFatal("PM/tree force norm diagnostics contain NaN or Inf");
    }

    ++summary.cheap_checks_executed;
    for (std::size_t i = 0; i < m_active_accel_x.size(); ++i) {
      if (!std::isfinite(m_active_accel_x[i]) ||
          !std::isfinite(m_active_accel_y[i]) ||
          !std::isfinite(m_active_accel_z[i])) {
        ++summary.force_non_finite_count;
        failFatal(
            "active-set acceleration contains NaN or Inf",
            {{"active_slot", std::to_string(i)}});
      }
    }

    ++summary.cheap_checks_executed;
    if (context.integrator_state.pm_sync_state.lastRefreshOpportunity() > context.integrator_state.pm_sync_state.gravityKickOpportunity() ||
        context.integrator_state.pm_sync_state.fieldVersion() < m_last_committed_field_version) {
      ++summary.illegal_sync_state_count;
      failFatal("gravity sync-state regressed across kick opportunities");
    }
    if (decision.refresh_long_range_field) {
      if (decision.field_version != context.integrator_state.pm_sync_state.fieldVersion() + 1U) {
        ++summary.illegal_sync_state_count;
        failFatal("refresh decision must increment field_version by exactly one");
      }
    } else if (decision.field_version != context.integrator_state.pm_sync_state.fieldVersion()) {
      ++summary.illegal_sync_state_count;
      failFatal("reuse decision must not mutate field_version");
    }

    ++summary.cheap_checks_executed;
    if (m_mode_policy.zoom_region_expected &&
        m_tree_pm_options.zoom_region_radius_comoving <= 0.0) {
      ++summary.zoom_sanity_failure_count;
      failFatal("zoom mode requires a strictly positive zoom region radius");
    }
    if (!m_mode_policy.zoom_region_expected &&
        m_tree_pm_options.enable_zoom_long_range_correction &&
        m_last_tree_pm_diagnostics.zoom_high_res_source_count == 0U &&
        m_last_tree_pm_diagnostics.zoom_low_res_source_count == 0U) {
      ++summary.zoom_sanity_failure_count;
      addWarning("zoom long-range correction enabled without detected zoom population");
    }

    if (allow_heavy_reference_checks) {
      ++summary.heavy_checks_executed;
      const double total_force = m_last_tree_pm_diagnostics.force_l2_total;
      const double pm_component = m_last_tree_pm_diagnostics.force_l2_pm_global;
      const double tree_component = m_last_tree_pm_diagnostics.force_l2_tree_short_range;
      if (total_force > 0.0 && std::isfinite(total_force)) {
        const double component_sum = pm_component + tree_component;
        const double relative_gap = std::abs(component_sum - total_force) / total_force;
        if (relative_gap > 1.0) {
          ++summary.force_abnormal_count;
          addWarning(
              "force decomposition appears abnormal under heavy reference check",
              {{"relative_gap", std::to_string(relative_gap)}});
        }
      }
      ++summary.heavy_checks_executed;
      if (m_mode_policy.zoom_region_expected &&
          m_last_tree_pm_diagnostics.zoom_low_res_contamination_count > 0U) {
        ++summary.force_abnormal_count;
        addWarning(
            "zoom contamination detected in low-resolution source particles",
            {{"zoom_low_res_contamination_count",
              std::to_string(m_last_tree_pm_diagnostics.zoom_low_res_contamination_count)}});
      }
    }

    return summary;
  }

  void ensureZoomMembershipLoaded(const core::SimulationState& state) {
    if (m_zoom_membership_loaded) {
      return;
    }
    m_zoom_membership_loaded = true;
    if (!m_mode_policy.zoom_region_expected || m_config.mode.zoom_region_file.empty()) {
      m_zoom_membership_source = m_mode_policy.zoom_region_expected ? "configured_spherical_region" : "disabled";
      return;
    }
    const std::filesystem::path path = m_zoom_region_path.empty()
        ? std::filesystem::path(m_config.mode.zoom_region_file)
        : m_zoom_region_path;
    if (hasHdf5Extension(path)) {
#if COSMOSIM_ENABLE_HDF5
      m_zoom_high_res_particle_ids = loadZoomParticleIdsFromHdf5(path);
      m_zoom_membership_source = "particle_id_file_hdf5";
#else
      throw std::runtime_error("zoom_region_file requires HDF5 support in this build");
#endif
    } else {
      m_zoom_high_res_particle_ids = loadZoomParticleIdsFromText(path);
      m_zoom_membership_source = "particle_id_file_text";
    }
    for (const std::uint64_t id : state.particle_sidecar.particle_id) {
      (void)id;
    }
  }

  void rebuildOwnedParticleCompactView(
      const core::StepContext& context,
      std::span<const std::uint32_t> active_particles,
      bool predict_inactive_sources = false) {
    const core::SimulationState& state = context.state;
    const std::uint32_t world_rank = static_cast<std::uint32_t>(m_runtime_topology.world_rank);
    const std::size_t particle_count = state.particles.size();
    m_owned_local_index_by_global.assign(particle_count, -1);
    m_local_source_x.clear();
    m_local_source_y.clear();
    m_local_source_z.clear();
    m_local_source_mass.clear();
    m_local_source_species_tag.clear();
    m_local_source_softening_comoving.clear();
    m_local_source_softening_override_mask.clear();
    m_local_to_global.clear();
    m_source_predicted_inactive_count = 0;

    std::vector<std::uint8_t> active_mask;
    if (predict_inactive_sources) {
      if (state.particle_sidecar.last_drift_time_code.size() != particle_count ||
          state.particle_sidecar.last_drift_scale_factor.size() != particle_count) {
        throw std::runtime_error("PM source prediction requires per-particle drift epoch sidecars");
      }
      active_mask.assign(particle_count, 0U);
      for (const std::uint32_t global_index : active_particles) {
        if (global_index >= particle_count) {
          throw std::out_of_range("gravity callback active particle index out of range");
        }
        active_mask[global_index] = 1U;
      }
    }

    const bool periodic_sources =
        m_tree_pm_options.pm_options.boundary_condition == gravity::PmBoundaryCondition::kPeriodic;
    const double box_size_x = m_config.cosmology.box_size_x_mpc_comoving;
    const double box_size_y = m_config.cosmology.box_size_y_mpc_comoving;
    const double box_size_z = m_config.cosmology.box_size_z_mpc_comoving;
    const bool has_softening_values = !state.particle_sidecar.gravity_softening_comoving.empty();
    const bool has_softening_masks = !state.particle_sidecar.has_gravity_softening_override.empty();
    for (std::size_t global_index = 0; global_index < particle_count; ++global_index) {
      const bool is_owned_source = state.particle_sidecar.owning_rank[global_index] == world_rank;
      const bool is_imported_ghost_source = !is_owned_source &&
          state.particle_sidecar.owning_rank[global_index] <
              static_cast<std::uint32_t>(std::max(m_runtime_topology.world_size, 1));
      if (!is_owned_source && !is_imported_ghost_source) {
        continue;
      }
      double source_x = state.particles.position_x_comoving[global_index];
      double source_y = state.particles.position_y_comoving[global_index];
      double source_z = state.particles.position_z_comoving[global_index];
      if (predict_inactive_sources && active_mask[global_index] == 0U) {
        const double source_time = state.particle_sidecar.last_drift_time_code[global_index];
        const double source_scale = state.particle_sidecar.last_drift_scale_factor[global_index];
        if (source_time > context.timeline_step.time_end_code + 1.0e-12 || source_scale <= 0.0) {
          throw std::runtime_error("inactive PM source has an invalid or future drift epoch");
        }
        double inactive_drift_factor = context.timeline_step.time_end_code - source_time;
        if (context.cosmology_background != nullptr) {
          inactive_drift_factor = core::computeComovingDriftFactor(
              *context.cosmology_background,
              source_scale,
              context.timeline_step.scale_factor_end,
              64) / context.timeline_step.time_si_per_code;
        }
        if (!std::isfinite(inactive_drift_factor) || inactive_drift_factor < -1.0e-14) {
          throw std::runtime_error("inactive PM source prediction produced an invalid drift factor");
        }
        source_x += state.particles.velocity_x_peculiar[global_index] * inactive_drift_factor;
        source_y += state.particles.velocity_y_peculiar[global_index] * inactive_drift_factor;
        source_z += state.particles.velocity_z_peculiar[global_index] * inactive_drift_factor;
        ++m_source_predicted_inactive_count;
      }
      if (periodic_sources) {
        source_x = wrapPeriodicPosition(source_x, box_size_x);
        source_y = wrapPeriodicPosition(source_y, box_size_y);
        source_z = wrapPeriodicPosition(source_z, box_size_z);
      }
      if (is_owned_source) {
        m_owned_local_index_by_global[global_index] = static_cast<int>(m_local_to_global.size());
      }
      m_local_to_global.push_back(static_cast<std::uint32_t>(global_index));
      m_local_source_x.push_back(source_x);
      m_local_source_y.push_back(source_y);
      m_local_source_z.push_back(source_z);
      m_local_source_mass.push_back(state.particles.mass_code[global_index]);
      m_local_source_species_tag.push_back(state.particle_sidecar.species_tag[global_index]);
      if (has_softening_values) {
        m_local_source_softening_comoving.push_back(state.particle_sidecar.gravity_softening_comoving[global_index]);
      }
      if (has_softening_masks) {
        m_local_source_softening_override_mask.push_back(
            state.particle_sidecar.hasGravitySofteningOverride(global_index) ? 1U : 0U);
      }
    }

    m_local_active_indices.clear();
    m_local_active_global_indices.clear();
    m_active_target_species_tag.clear();
    m_active_target_softening_comoving.clear();
    m_active_target_softening_override_mask.clear();
    for (const std::uint32_t global_index : active_particles) {
      if (global_index >= particle_count) {
        throw std::out_of_range("gravity callback active particle index out of range");
      }
      const int local_index = m_owned_local_index_by_global[global_index];
      if (local_index < 0) {
        continue;
      }
      const auto local_u32 = static_cast<std::uint32_t>(local_index);
      m_local_active_indices.push_back(local_u32);
      m_local_active_global_indices.push_back(global_index);
      m_active_target_species_tag.push_back(m_local_source_species_tag[local_u32]);
      if (has_softening_values) {
        m_active_target_softening_comoving.push_back(m_local_source_softening_comoving[local_u32]);
      }
      if (has_softening_masks) {
        m_active_target_softening_override_mask.push_back(m_local_source_softening_override_mask[local_u32]);
      }
    }
  }

  const core::SimulationConfig& m_config;
  const core::ModePolicy& m_mode_policy;
  parallel::DistributedExecutionTopology m_runtime_topology{};
  std::uint64_t m_pm_update_cadence_steps = 1;
  gravity::PmGridShape m_pm_grid_shape{};
  double m_mesh_spacing_x_mpc_comoving = 0.0;
  double m_mesh_spacing_y_mpc_comoving = 0.0;
  double m_mesh_spacing_z_mpc_comoving = 0.0;
  std::string m_pm_assignment_scheme = "unknown";
  std::string m_pm_backend = "unknown";
  gravity::TreePmCoordinator m_tree_pm_coordinator;
  gravity::TreePmOptions m_tree_pm_options;
  gravity::TreeSofteningSpeciesPolicy m_tree_pm_species_softening{};
  std::vector<double> m_active_accel_x;
  std::vector<double> m_active_accel_y;
  std::vector<double> m_active_accel_z;
  std::vector<double> m_particle_accel_x;
  std::vector<double> m_particle_accel_y;
  std::vector<double> m_particle_accel_z;
  std::vector<double> m_cell_accel_x;
  std::vector<double> m_cell_accel_y;
  std::vector<double> m_cell_accel_z;
  std::vector<int> m_active_slot_by_particle;
  std::vector<int> m_owned_local_index_by_global;
  std::vector<double> m_local_source_x;
  std::vector<double> m_local_source_y;
  std::vector<double> m_local_source_z;
  std::vector<double> m_local_source_mass;
  std::vector<std::uint32_t> m_local_source_species_tag;
  std::vector<double> m_local_source_softening_comoving;
  std::vector<std::uint8_t> m_local_source_softening_override_mask;
  std::vector<std::uint32_t> m_active_target_species_tag;
  std::vector<double> m_active_target_softening_comoving;
  std::vector<std::uint8_t> m_active_target_softening_override_mask;
  std::vector<std::uint8_t> m_source_is_high_res;
  std::vector<std::uint8_t> m_active_is_high_res;
  std::vector<std::uint32_t> m_local_to_global;
  std::filesystem::path m_zoom_region_path;
  bool m_zoom_membership_loaded = false;
  std::unordered_set<std::uint64_t> m_zoom_high_res_particle_ids;
  std::string m_zoom_membership_source = "disabled";
  std::vector<std::uint32_t> m_local_active_indices;
  std::vector<std::uint32_t> m_local_active_global_indices;
  gravity::TreePmDiagnostics m_last_tree_pm_diagnostics{};
  parallel::DecompositionRuntimeMeasurements m_last_decomposition_measurements{};
  bool m_has_long_range_field = false;
  std::uint64_t m_last_committed_field_version = 0;
  std::uint64_t m_long_range_refresh_count = 0;
  std::uint64_t m_long_range_reuse_count = 0;
  bool m_force_cache_valid = false;
  parallel::GhostCacheLifecycle m_ghost_cache_lifecycle{};
  std::uint64_t m_source_predicted_inactive_count = 0;
  std::vector<ReferenceWorkflowReport::TreePmCadenceRecord> m_cadence_records;
};

class HydroStageCallback final : public core::IntegrationCallback {
 public:
  HydroStageCallback(const core::SimulationConfig& config, const core::ModePolicy& mode_policy, const GravityStageCallback& gravity_callback)
      : m_config(config),
        m_mode_policy(mode_policy),
        m_gravity_callback(gravity_callback),
        m_solver(k_gamma_adiabatic),
        m_reconstruction(hydro::HydroReconstructionPolicy{
            .limiter = hydro::HydroSlopeLimiter::kMonotonizedCentral,
            .dt_over_dx_code = 0.0,
            .rho_floor = k_density_floor,
            .pressure_floor = k_pressure_floor,
            .enable_muscl_hancock_predictor = true,
            .adiabatic_index = k_gamma_adiabatic,
        }) {}

  [[nodiscard]] std::string_view callbackName() const override { return "hydro"; }
  [[nodiscard]] std::span<const core::IntegrationStage> integrationStages() const override { return k_hydro_stage; }
  [[nodiscard]] std::span<const core::StageContract> stageContracts() const override { return m_contracts; }
  [[nodiscard]] const hydro::HydroProfileEvent& lastHydroProfile() const noexcept { return m_last_hydro_profile; }
  [[nodiscard]] const SolverGhostRefreshReport& lastGhostRefreshReport() const noexcept { return m_last_ghost_refresh; }
  [[nodiscard]] const core::HydroCflDiagnostics& lastHydroCflDiagnostics() const noexcept {
    return m_last_hydro_cfl_diagnostics;
  }

  void onStage(core::StepContext& context) override {
    if (context.stage != core::IntegrationStage::kHydroUpdate) {
      throw std::logic_error("hydro handler received an unregistered stage");
    }
    if (context.state.cells.size() == 0) {
      return;
    }

    parallel::MpiContext mpi_context;
    m_current_world_rank = static_cast<std::uint32_t>(std::max(mpi_context.worldRank(), 0));
    m_last_ghost_refresh = refreshParticleGhostsForSolver(
        context, mpi_context, "hydro.godunov", &m_ghost_cache_lifecycle);
    m_last_hydro_profile = {};

    context.state.requireGasCellIdentityMapCoversDenseRows("hydro callback");
    if (amr::hasProductionAmrHydroCoverage(context.state)) {
      runProductionAmrHydroPath(context);
      return;
    }
    const auto particle_row_by_id = buildParticleRowById(context.state);
    std::vector<std::uint8_t> parent_mirror_updated(context.state.particles.size(), 0U);

    rebuildGeometryIfNeeded(context.state, context.integrator_state.dt_time_code);
    const hydro::HydroActiveSetView active_view = buildActiveFaceView(context.active_set.cell_indices);
    verifyAcceptedHydroCfl(context, context.active_set.cell_indices);

    m_conserved.resize(m_geometry.totalCellStorageCount());
    for (std::size_t cell_index = 0; cell_index < context.state.cells.size(); ++cell_index) {
      const double rho = std::max(context.state.gas_cells.density_code[cell_index], k_density_floor);
      double pressure = context.state.gas_cells.pressure_code[cell_index];
      if (pressure <= 0.0) {
        const double internal_energy = std::max(context.state.gas_cells.internal_energy_code[cell_index], k_pressure_floor);
        pressure = std::max((k_gamma_adiabatic - 1.0) * rho * internal_energy, k_pressure_floor);
      }
      const hydro::HydroPrimitiveState primitive{
          .rho_comoving = rho,
          .vel_x_peculiar = context.state.gas_cells.velocity_x_peculiar[cell_index],
          .vel_y_peculiar = context.state.gas_cells.velocity_y_peculiar[cell_index],
          .vel_z_peculiar = context.state.gas_cells.velocity_z_peculiar[cell_index],
          .pressure_comoving = pressure,
      };
      if (cell_index >= m_geometry_row_by_dense_row.size()) {
        throw std::out_of_range("hydro callback dense row is outside Cartesian geometry map");
      }
      m_conserved.storeCell(
          m_geometry_row_by_dense_row[cell_index],
          hydro::HydroCoreSolver::conservedFromPrimitive(primitive, k_gamma_adiabatic));
    }
    hydro::fillHydroBoundaryGhostCells(m_conserved, m_geometry, k_gamma_adiabatic);
    const std::uint32_t world_rank = static_cast<std::uint32_t>(mpi_context.worldRank());
    const HydroGhostConservedSnapshot ghost_conserved_snapshot = snapshotHydroGhostConservedCells(
        context.state, m_conserved, m_geometry_row_by_dense_row, world_rank);

    double hubble_rate_code = 0.0;
    if (context.cosmology_background != nullptr && context.integrator_state.current_scale_factor > 0.0) {
      hubble_rate_code =
          core::computeScaleFactorRate(*context.cosmology_background, context.integrator_state.current_scale_factor) /
          context.integrator_state.current_scale_factor;
    }

    hydro::HydroUpdateContext update{
        .dt_code = context.integrator_state.dt_time_code,
        .scale_factor = std::max(1.0e-12, context.integrator_state.current_scale_factor),
        .hubble_rate_code = hubble_rate_code,
    };

    std::vector<double> metallicity(context.state.cells.size(), 0.0);
    std::vector<double> temperature(context.state.cells.size(), 0.0);
    std::vector<double> hydrogen_number_density(context.state.cells.size(), 0.0);
    const auto& dense_accel_x = m_gravity_callback.cellAccelX();
    const auto& dense_accel_y = m_gravity_callback.cellAccelY();
    const auto& dense_accel_z = m_gravity_callback.cellAccelZ();
    if (dense_accel_x.size() != context.state.cells.size() ||
        dense_accel_y.size() != context.state.cells.size() ||
        dense_accel_z.size() != context.state.cells.size()) {
      throw std::runtime_error("hydro callback gravity-cell acceleration lanes do not cover dense gas-cell rows");
    }
    m_ordered_cell_accel_x.assign(context.state.cells.size(), 0.0);
    m_ordered_cell_accel_y.assign(context.state.cells.size(), 0.0);
    m_ordered_cell_accel_z.assign(context.state.cells.size(), 0.0);
    for (std::size_t geometry_row = 0; geometry_row < m_dense_row_by_geometry_row.size(); ++geometry_row) {
      const std::uint32_t dense_row = m_dense_row_by_geometry_row[geometry_row];
      m_ordered_cell_accel_x[geometry_row] = dense_accel_x[dense_row];
      m_ordered_cell_accel_y[geometry_row] = dense_accel_y[dense_row];
      m_ordered_cell_accel_z[geometry_row] = dense_accel_z[dense_row];
    }
    hydro::HydroSourceContext source_context{
        .update = update,
        .gravity_accel_x_peculiar = m_ordered_cell_accel_x,
        .gravity_accel_y_peculiar = m_ordered_cell_accel_y,
        .gravity_accel_z_peculiar = m_ordered_cell_accel_z,
        .hydrogen_number_density_cgs = hydrogen_number_density,
        .metallicity_mass_fraction = metallicity,
        .temperature_k = temperature,
        .redshift = std::max(0.0, (update.scale_factor > 0.0 ? (1.0 / update.scale_factor - 1.0) : 0.0)),
    };

    hydro::ComovingGravityExpansionSource gravity_source;
    std::array<const hydro::HydroSourceTerm*, 1> sources{&gravity_source};
    m_solver.advancePatchActiveSetWithScratch(
        m_conserved,
        m_geometry,
        active_view,
        update,
        m_reconstruction,
        m_riemann_solver,
        sources,
        source_context,
        m_scratch,
        &m_primitive_cache,
        &m_last_hydro_profile);
    if (context.profiler_session != nullptr) {
      const hydro::HydroConservationReport& conservation = m_last_hydro_profile.conservation;
      context.profiler_session->recordEvent(core::RuntimeEvent{
          .event_kind = "hydro.conservation",
          .severity = core::RuntimeEventSeverity::kInfo,
          .subsystem = "hydro.godunov",
          .step_index = context.integrator_state.step_index,
          .simulation_time_code = context.integrator_state.current_time_code,
          .scale_factor = context.integrator_state.current_scale_factor,
          .message = "computed volume-integrated hydro conservation totals for the local update stage",
          .payload = {{"cell_count", std::to_string(conservation.cell_count)},
                      {"before_mass", std::to_string(conservation.before.mass)},
                      {"after_mass", std::to_string(conservation.after.mass)},
                      {"flux_delta_mass", std::to_string(conservation.flux_delta.mass)},
                      {"source_delta_mass", std::to_string(conservation.source_delta.mass)},
                      {"floor_delta_mass", std::to_string(conservation.floor_delta.mass)},
                      {"residual_mass", std::to_string(conservation.residual.mass)},
                      {"residual_momentum_x", std::to_string(conservation.residual.momentum_x)},
                      {"residual_momentum_y", std::to_string(conservation.residual.momentum_y)},
                      {"residual_momentum_z", std::to_string(conservation.residual.momentum_z)},
                      {"residual_total_energy", std::to_string(conservation.residual.total_energy)},
                      {"residual_internal_energy", std::to_string(conservation.residual.internal_energy)},
                      {"flux_delta_total_energy", std::to_string(conservation.flux_delta.total_energy)},
                      {"source_delta_total_energy", std::to_string(conservation.source_delta.total_energy)},
                      {"floor_delta_total_energy", std::to_string(conservation.floor_delta.total_energy)},
                      {"internal_energy_floor_count", std::to_string(conservation.internal_energy_floor_count)},
                      {"mass_tolerance", "1e-10"},
                      {"momentum_tolerance", "1e-10"},
                      {"energy_tolerance", "1e-10"}}});
    }

    const HydroConservativeGhostSyncReport ghost_sync_report = restoreHydroGhostConservedCells(
        m_conserved, ghost_conserved_snapshot, world_rank);
    const std::vector<parallel::HydroConservativeFluxCorrectionRecord> global_flux_corrections =
        parallel::executeBlockingHydroConservativeFluxCorrectionExchange(
            mpi_context, ghost_sync_report.correction_records, context.integrator_state.step_index);
    std::size_t applied_flux_corrections = 0;
    double applied_flux_delta_l1 = 0.0;
    for (const parallel::HydroConservativeFluxCorrectionRecord& correction : global_flux_corrections) {
      if (correction.owner_rank != static_cast<int>(world_rank)) {
        continue;
      }
      const auto correction_row = context.state.rowForGasCellId(correction.gas_cell_id);
      if (!correction_row.has_value()) {
        throw std::runtime_error("hydro conservative flux correction targeted an unknown local gas_cell_id");
      }
      const std::uint32_t cell_index = *correction_row;
      const std::uint32_t owner_rank = gasCellOwnerRankForLocalRow(
          context.state,
          cell_index,
          particle_row_by_id,
          "hydro conservative flux correction");
      if (owner_rank != world_rank) {
        throw std::runtime_error("hydro conservative flux correction targeted a non-authoritative local ghost cell");
      }
      if (cell_index >= m_geometry_row_by_dense_row.size()) {
        throw std::out_of_range("hydro flux correction dense row is outside Cartesian geometry map");
      }
      const std::uint32_t geometry_row = m_geometry_row_by_dense_row[cell_index];
      hydro::HydroConservedState owned = m_conserved.loadCell(geometry_row);
      owned.mass_density_comoving += correction.delta_mass_density_comoving;
      owned.momentum_density_x_comoving += correction.delta_momentum_density_x_comoving;
      owned.momentum_density_y_comoving += correction.delta_momentum_density_y_comoving;
      owned.momentum_density_z_comoving += correction.delta_momentum_density_z_comoving;
      owned.total_energy_density_comoving += correction.delta_total_energy_density_comoving;
      m_conserved.storeCell(geometry_row, owned);
      ++applied_flux_corrections;
      applied_flux_delta_l1 += std::abs(correction.delta_mass_density_comoving);
      applied_flux_delta_l1 += std::abs(correction.delta_momentum_density_x_comoving);
      applied_flux_delta_l1 += std::abs(correction.delta_momentum_density_y_comoving);
      applied_flux_delta_l1 += std::abs(correction.delta_momentum_density_z_comoving);
      applied_flux_delta_l1 += std::abs(correction.delta_total_energy_density_comoving);
    }
    if (context.profiler_session != nullptr) {
      context.profiler_session->recordEvent(core::RuntimeEvent{
          .event_kind = "hydro.conservative_ghost_sync",
          .severity = core::RuntimeEventSeverity::kInfo,
          .subsystem = "hydro.godunov",
          .step_index = context.integrator_state.step_index,
          .simulation_time_code = context.integrator_state.current_time_code,
          .scale_factor = context.integrator_state.current_scale_factor,
          .message = "restored imported hydro ghosts and applied conservative owner-side boundary corrections",
          .payload = {{"restored_ghost_cells", std::to_string(ghost_sync_report.restored_ghost_cells)},
                      {"local_correction_records", std::to_string(ghost_sync_report.correction_records.size())},
                      {"global_correction_records", std::to_string(global_flux_corrections.size())},
                      {"applied_flux_corrections", std::to_string(applied_flux_corrections)},
                      {"rejected_remote_delta_l1", std::to_string(ghost_sync_report.rejected_remote_delta_l1)},
                      {"applied_flux_delta_l1", std::to_string(applied_flux_delta_l1)}}});
    }

    for (std::size_t cell_index = 0; cell_index < context.state.cells.size(); ++cell_index) {
      const std::uint32_t gas_cell_owner_rank = gasCellOwnerRankForLocalRow(
          context.state,
          static_cast<std::uint32_t>(cell_index),
          particle_row_by_id,
          "hydro callback primitive store");
      if (gas_cell_owner_rank != world_rank) {
        continue;
      }
      if (cell_index >= m_geometry_row_by_dense_row.size()) {
        throw std::out_of_range("hydro primitive store dense row is outside Cartesian geometry map");
      }
      const std::uint32_t geometry_row = m_geometry_row_by_dense_row[cell_index];
      const hydro::HydroPrimitiveState primitive =
          hydro::HydroCoreSolver::primitiveFromConserved(m_conserved.loadCell(geometry_row), k_gamma_adiabatic);
      context.state.gas_cells.density_code[cell_index] = primitive.rho_comoving;
      context.state.gas_cells.pressure_code[cell_index] = primitive.pressure_comoving;
      context.state.gas_cells.internal_energy_code[cell_index] =
          primitive.pressure_comoving / ((k_gamma_adiabatic - 1.0) * std::max(primitive.rho_comoving, k_density_floor));
      context.state.cells.mass_code[cell_index] = primitive.rho_comoving * m_geometry.cell_volume_comoving;
      context.state.gas_cells.velocity_x_peculiar[cell_index] = primitive.vel_x_peculiar;
      context.state.gas_cells.velocity_y_peculiar[cell_index] = primitive.vel_y_peculiar;
      context.state.gas_cells.velocity_z_peculiar[cell_index] = primitive.vel_z_peculiar;
    }

    synchronizeParentParticleCompatibilityMirrors(
        context.state, world_rank, "hydro callback parent compatibility mirror");
  }

 private:
  static constexpr std::array<core::StageContract, 1> m_contracts{{
      {.stage = core::IntegrationStage::kHydroUpdate,
       .required_inputs = core::StageDataDomain::kGasCells | core::StageDataDomain::kParticles | core::StageDataDomain::kGhostCells,
       .mutated_state = core::StageDataDomain::kGasCells | core::StageDataDomain::kParticles,
       .produced_outputs = core::StageDataDomain::kGasCells | core::StageDataDomain::kParticles,
       .allowed_side_effects = core::StageDataDomain::kNone,
       .sync_requirements = core::StageSyncRequirement::kLocalOnly,
       .active_set_family = core::StageActiveSetFamily::kGasCells,
       .restart_safety = core::StageSafety::kUnsafe,
       .output_safety = core::StageSafety::kUnsafe,
       .owner = core::StageSubsystem::kHydro},
  }};

  [[nodiscard]] static hydro::HydroBoundaryKind hydroBoundaryKindFromModePolicy(
      core::BoundaryCondition boundary_condition) {
    switch (boundary_condition) {
      case core::BoundaryCondition::kPeriodic:
        return hydro::HydroBoundaryKind::kPeriodic;
      case core::BoundaryCondition::kOpen:
        return hydro::HydroBoundaryKind::kOpen;
      case core::BoundaryCondition::kReflective:
        return hydro::HydroBoundaryKind::kReflective;
    }
    return hydro::HydroBoundaryKind::kOpen;
  }

  void runProductionAmrHydroPath(core::StepContext& context) {
    double hubble_rate_code = 0.0;
    if (context.cosmology_background != nullptr && context.integrator_state.current_scale_factor > 0.0) {
      hubble_rate_code =
          core::computeScaleFactorRate(*context.cosmology_background, context.integrator_state.current_scale_factor) /
          context.integrator_state.current_scale_factor;
    }

    const hydro::HydroUpdateContext update{
        .dt_code = context.integrator_state.dt_time_code,
        .scale_factor = std::max(1.0e-12, context.integrator_state.current_scale_factor),
        .hubble_rate_code = hubble_rate_code,
    };
    std::vector<double> metallicity(context.state.cells.size(), 0.0);
    std::vector<double> temperature(context.state.cells.size(), 0.0);
    std::vector<double> hydrogen_number_density(context.state.cells.size(), 0.0);
    const hydro::HydroSourceContext source_context{
        .update = update,
        .gravity_accel_x_peculiar = m_gravity_callback.cellAccelX(),
        .gravity_accel_y_peculiar = m_gravity_callback.cellAccelY(),
        .gravity_accel_z_peculiar = m_gravity_callback.cellAccelZ(),
        .hydrogen_number_density_cgs = hydrogen_number_density,
        .metallicity_mass_fraction = metallicity,
        .temperature_k = temperature,
        .redshift = std::max(0.0, (update.scale_factor > 0.0 ? (1.0 / update.scale_factor - 1.0) : 0.0)),
    };
    hydro::ComovingGravityExpansionSource gravity_source;
    std::array<const hydro::HydroSourceTerm*, 1> sources{&gravity_source};
    const amr::ProductionAmrHydroDiagnostics amr_diagnostics = amr::advanceProductionAmrHydro(
        context.state,
        context.active_set.cell_indices,
        update,
        source_context,
        m_solver,
        m_riemann_solver,
        sources,
        amr::ProductionAmrHydroOptions{
            .physical_boundary_kind = hydroBoundaryKindFromModePolicy(m_mode_policy.hydro_boundary),
            .adiabatic_index = k_gamma_adiabatic,
            .density_floor = k_density_floor,
            .pressure_floor = k_pressure_floor});
    synchronizeParentParticleCompatibilityMirrors(
        context.state,
        m_current_world_rank,
        "production AMR hydro parent compatibility mirror");
    if (context.profiler_session != nullptr) {
      context.profiler_session->recordEvent(core::RuntimeEvent{
          .event_kind = "hydro.amr_production_stage",
          .severity = core::RuntimeEventSeverity::kInfo,
          .subsystem = "hydro.amr",
          .step_index = context.integrator_state.step_index,
          .simulation_time_code = context.integrator_state.current_time_code,
          .scale_factor = context.integrator_state.current_scale_factor,
          .message = "advanced production AMR hydro patches from SimulationState-authoritative gas rows",
          .payload = {{"patch_count", std::to_string(amr_diagnostics.patch_count)},
                      {"advanced_patch_count", std::to_string(amr_diagnostics.advanced_patch_count)},
                      {"active_cell_count", std::to_string(amr_diagnostics.active_cell_count)},
                      {"active_face_count", std::to_string(amr_diagnostics.active_face_count)},
                      {"flux_register_entry_count", std::to_string(amr_diagnostics.flux_register_entry_count)},
                      {"reflux_corrected_cells", std::to_string(amr_diagnostics.reflux.corrected_cells)},
                      {"reflux_corrected_mass", std::to_string(amr_diagnostics.reflux.corrected_mass_code)},
                      {"reflux_corrected_momentum_x", std::to_string(amr_diagnostics.reflux.corrected_momentum_x_code)},
                      {"reflux_corrected_momentum_y", std::to_string(amr_diagnostics.reflux.corrected_momentum_y_code)},
                      {"reflux_corrected_momentum_z", std::to_string(amr_diagnostics.reflux.corrected_momentum_z_code)},
                      {"reflux_corrected_total_energy", std::to_string(amr_diagnostics.reflux.corrected_total_energy_code)},
                      {"reflux_corrected_energy", std::to_string(amr_diagnostics.reflux.corrected_energy_code)},
                      {"reflux_corrected_internal_energy", std::to_string(amr_diagnostics.reflux.corrected_internal_energy_code)},
                      {"reflux_complete_register_count", std::to_string(amr_diagnostics.reflux.complete_register_count)},
                      {"reflux_skipped_incomplete_register_count", std::to_string(amr_diagnostics.reflux.skipped_incomplete_register_count)},
                      {"reflux_skipped_area_mismatch_count", std::to_string(amr_diagnostics.reflux.skipped_area_mismatch_count)},
                      {"reflux_skipped_missing_target_count", std::to_string(amr_diagnostics.reflux.skipped_missing_target_count)}}});
    }
  }

  struct HydroGeometryCacheKey {
    std::size_t cell_count = 0;
    std::size_t nx = 0;
    std::size_t ny = 0;
    std::size_t nz = 0;
    double origin_x_comoving = 0.0;
    double origin_y_comoving = 0.0;
    double origin_z_comoving = 0.0;
    double cell_width_x_comoving = 0.0;
    double cell_width_y_comoving = 0.0;
    double cell_width_z_comoving = 0.0;
    double dt_time_code = 0.0;
    core::BoundaryCondition hydro_boundary = core::BoundaryCondition::kPeriodic;
    std::uint64_t patch_signature = 0;
    std::uint64_t gas_cell_identity_generation = 0;
    std::uint64_t cell_index_generation = 0;
    std::uint64_t row_mapping_signature = 0;

    [[nodiscard]] bool operator==(const HydroGeometryCacheKey& rhs) const noexcept {
      return cell_count == rhs.cell_count &&
          nx == rhs.nx &&
          ny == rhs.ny &&
          nz == rhs.nz &&
          origin_x_comoving == rhs.origin_x_comoving &&
          origin_y_comoving == rhs.origin_y_comoving &&
          origin_z_comoving == rhs.origin_z_comoving &&
          cell_width_x_comoving == rhs.cell_width_x_comoving &&
          cell_width_y_comoving == rhs.cell_width_y_comoving &&
          cell_width_z_comoving == rhs.cell_width_z_comoving &&
          dt_time_code == rhs.dt_time_code &&
          hydro_boundary == rhs.hydro_boundary &&
          patch_signature == rhs.patch_signature &&
          gas_cell_identity_generation == rhs.gas_cell_identity_generation &&
          cell_index_generation == rhs.cell_index_generation &&
          row_mapping_signature == rhs.row_mapping_signature;
    }
  };

  [[nodiscard]] static std::uint64_t patchSignature(const core::PatchSoa& patches) {
    std::uint64_t signature = 1469598103934665603ULL;
    const auto mix = [&signature](std::uint64_t value) {
      signature ^= value;
      signature *= 1099511628211ULL;
    };
    mix(static_cast<std::uint64_t>(patches.size()));
    for (std::size_t patch = 0; patch < patches.size(); ++patch) {
      mix(patches.patch_id[patch]);
      mix(static_cast<std::uint64_t>(static_cast<std::int64_t>(patches.level[patch])));
      mix(patches.first_cell[patch]);
      mix(patches.cell_count[patch]);
      mix(patches.owning_rank[patch]);
    }
    return signature;
  }

  void rebuildGeometryIfNeeded(const core::SimulationState& state, double dt_time_code) {
    const std::size_t cell_count = state.cells.size();
    if (cell_count == 0) {
      m_cached_geometry_key.reset();
      m_geometry = {};
      return;
    }

    const internal::CartesianGasCellRowLayout layout = requireCartesianGasCellRowLayout(
        state, m_config, "reference hydro workflow");
    const hydro::HydroCartesianPatchSpec& spec = layout.spec;
    const HydroGeometryCacheKey key{
        .cell_count = cell_count,
        .nx = spec.nx,
        .ny = spec.ny,
        .nz = spec.nz,
        .origin_x_comoving = spec.origin_x_comoving,
        .origin_y_comoving = spec.origin_y_comoving,
        .origin_z_comoving = spec.origin_z_comoving,
        .cell_width_x_comoving = spec.cell_width_x_comoving,
        .cell_width_y_comoving = spec.cell_width_y_comoving,
        .cell_width_z_comoving = spec.cell_width_z_comoving,
        .dt_time_code = dt_time_code,
        .hydro_boundary = m_mode_policy.hydro_boundary,
        .patch_signature = patchSignature(state.patches),
        .gas_cell_identity_generation = state.gasCellIdentityGeneration(),
        .cell_index_generation = state.cellIndexGeneration(),
        .row_mapping_signature = layout.mapping_signature,
    };
    if (m_cached_geometry_key.has_value() && *m_cached_geometry_key == key) {
      return;
    }

    m_cached_geometry_key = key;
    m_dense_row_by_geometry_row = layout.dense_row_by_geometry_row;
    m_geometry_row_by_dense_row = layout.geometry_row_by_dense_row;
    m_geometry = hydro::makeCartesianPatchGeometry(spec);
    hydro::appendCartesianBoundaryGhostFaces(
        m_geometry,
        hydroBoundaryKindFromModePolicy(m_mode_policy.hydro_boundary));
    const double dx = std::max(1.0e-6, spec.cell_width_x_comoving);
    const double dy = std::max(1.0e-6, spec.cell_width_y_comoving);
    const double dz = std::max(1.0e-6, spec.cell_width_z_comoving);
    m_reconstruction = hydro::MusclHancockReconstruction(hydro::HydroReconstructionPolicy{
        .limiter = hydro::HydroSlopeLimiter::kMonotonizedCentral,
        .dt_over_dx_code = dt_time_code / dx,
        .dt_over_cell_width_code = {dt_time_code / dx, dt_time_code / dy, dt_time_code / dz},
        .rho_floor = k_density_floor,
        .pressure_floor = k_pressure_floor,
        .enable_muscl_hancock_predictor = true,
        .adiabatic_index = k_gamma_adiabatic,
    });
  }

  [[nodiscard]] hydro::HydroActiveSetView buildActiveFaceView(std::span<const std::uint32_t> active_cells) {
    m_active_cells.clear();
    m_active_cells.reserve(active_cells.size());
    for (const std::uint32_t dense_row : active_cells) {
      if (dense_row >= m_geometry_row_by_dense_row.size()) {
        throw std::out_of_range("active gas-cell row is outside Cartesian row map");
      }
      m_active_cells.push_back(m_geometry_row_by_dense_row[dense_row]);
    }
    std::unordered_set<std::size_t> active_cell_lookup(m_active_cells.begin(), m_active_cells.end());
    m_active_faces.clear();
    for (std::size_t face_index = 0; face_index < m_geometry.faces.size(); ++face_index) {
      const hydro::HydroFace& face = m_geometry.faces[face_index];
      if (active_cell_lookup.contains(face.owner_cell) ||
          (face.neighbor_cell != hydro::k_invalid_cell_index && active_cell_lookup.contains(face.neighbor_cell))) {
        m_active_faces.push_back(face_index);
      }
    }
    return hydro::HydroActiveSetView{.active_cells = m_active_cells, .active_faces = m_active_faces};
  }

  [[nodiscard]] std::optional<std::pair<std::uint64_t, std::uint32_t>> patchIdentityForCellRow(
      const core::SimulationState& state,
      std::uint32_t cell_index) const {
    const auto* identity = state.gas_cell_identity.findByLocalRow(cell_index);
    if (identity == nullptr || identity->gas_cell_id == 0U) {
      throw std::runtime_error("hydro CFL diagnostics rejected missing gas-cell identity");
    }
    if (state.patches.size() == 0U) {
      return std::nullopt;
    }
    if (identity->owning_patch_id == 0U || cell_index >= m_geometry_row_by_dense_row.size()) {
      throw std::runtime_error("hydro CFL diagnostics rejected incomplete fixed-patch identity metadata");
    }
    const auto patch_it = std::find(
        state.patches.patch_id.begin(), state.patches.patch_id.end(), identity->owning_patch_id);
    if (patch_it == state.patches.patch_id.end()) {
      throw std::runtime_error("hydro CFL diagnostics rejected an identity patch absent from PatchSoa");
    }
    const std::size_t patch_index = static_cast<std::size_t>(std::distance(state.patches.patch_id.begin(), patch_it));
    if (state.cells.patch_index[cell_index] != patch_index) {
      throw std::runtime_error("hydro CFL diagnostics rejected stale dense patch-index metadata");
    }
    return std::pair<std::uint64_t, std::uint32_t>{
        identity->owning_patch_id,
        m_geometry_row_by_dense_row[cell_index]};
  }

  [[nodiscard]] core::HydroCflDiagnostics hydroCflDiagnosticsForCell(
      const core::SimulationState& state,
      std::uint32_t cell_index,
      double accepted_dt_time_code) const {
    const core::DirectionalCflTimeStepInput input{
        .cell_width_axis_code = {
            m_geometry.cell_width_x_comoving,
            m_geometry.cell_width_y_comoving,
            m_geometry.cell_width_z_comoving},
        .velocity_axis_code = {
            state.gas_cells.velocity_x_peculiar[cell_index],
            state.gas_cells.velocity_y_peculiar[cell_index],
            state.gas_cells.velocity_z_peculiar[cell_index]},
        .sound_speed_code = std::max(state.gas_cells.sound_speed_code[cell_index], 0.0),
    };
    const auto patch_identity = patchIdentityForCellRow(state, cell_index);
    return core::makeHydroCflDiagnostics(
        cell_index,
        input,
        0.4,
        accepted_dt_time_code,
        state.gas_cell_identity.findByLocalRow(cell_index) != nullptr
            ? state.gas_cell_identity.findByLocalRow(cell_index)->gas_cell_id
            : 0U,
        patch_identity.has_value() ? std::optional<std::uint64_t>(patch_identity->first) : std::nullopt,
        patch_identity.has_value() ? std::optional<std::uint32_t>(patch_identity->second) : std::nullopt);
  }

  void verifyAcceptedHydroCfl(
      core::StepContext& context,
      std::span<const std::uint32_t> active_cells) {
    if (context.integrator_state.dt_time_code <= 0.0) {
      throw std::invalid_argument("hydro CFL guard requires dt_time_code > 0");
    }
    const std::size_t cell_count = context.state.cells.size();
    if (cell_count == 0) {
      return;
    }

    std::vector<std::uint32_t> all_cells;
    std::span<const std::uint32_t> cells_to_check = active_cells;
    if (cells_to_check.empty()) {
      all_cells.resize(cell_count);
      std::iota(all_cells.begin(), all_cells.end(), 0U);
      cells_to_check = all_cells;
    }

    core::HydroCflDiagnostics worst;
    bool have_worst = false;
    for (const std::size_t cell : cells_to_check) {
      if (cell >= cell_count) {
        throw std::out_of_range("hydro CFL guard active cell index out of range");
      }
      const core::HydroCflDiagnostics diagnostics = hydroCflDiagnosticsForCell(
          context.state,
          static_cast<std::uint32_t>(cell),
          context.integrator_state.dt_time_code);
      if (!have_worst || diagnostics.safety_factor < worst.safety_factor) {
        worst = diagnostics;
        have_worst = true;
      }
    }
    if (!have_worst) {
      return;
    }
    m_last_hydro_cfl_diagnostics = worst;
    if (context.profiler_session != nullptr) {
      context.profiler_session->recordEvent(core::RuntimeEvent{
          .event_kind = "hydro.cfl_guard",
          .severity = core::RuntimeEventSeverity::kInfo,
          .subsystem = "hydro.godunov",
          .step_index = context.integrator_state.step_index,
          .simulation_time_code = context.integrator_state.current_time_code,
          .scale_factor = context.integrator_state.current_scale_factor,
          .message = "checked accepted hydro timestep against local directional CFL bound",
          .payload = {{"local_row", std::to_string(worst.local_row)},
                      {"gas_cell_id", worst.has_gas_cell_id ? std::to_string(worst.gas_cell_id) : "n/a"},
                      {"patch_id", worst.has_patch_id ? std::to_string(worst.patch_id) : "n/a"},
                      {"patch_row", worst.has_patch_row ? std::to_string(worst.patch_row) : "n/a"},
                      {"proposed_dt_time_code", std::to_string(worst.proposed_dt_time_code)},
                      {"accepted_dt_time_code", std::to_string(worst.accepted_dt_time_code)},
                      {"cfl_number", std::to_string(worst.cfl_number)},
                      {"safety_factor", std::to_string(worst.safety_factor)},
                      {"velocity_x_code", std::to_string(worst.velocity_axis_code[0])},
                      {"velocity_y_code", std::to_string(worst.velocity_axis_code[1])},
                      {"velocity_z_code", std::to_string(worst.velocity_axis_code[2])},
                      {"sound_speed_code", std::to_string(worst.sound_speed_code)}}});
    }
    core::assertHydroCflStable(worst);
  }

  const core::SimulationConfig& m_config;
  const core::ModePolicy& m_mode_policy;
  const GravityStageCallback& m_gravity_callback;
  std::uint32_t m_current_world_rank = 0;
  hydro::HydroCoreSolver m_solver;
  hydro::MusclHancockReconstruction m_reconstruction;
  hydro::HllcRiemannSolver m_riemann_solver;
  hydro::HydroConservedStateSoa m_conserved;
  hydro::HydroScratchBuffers m_scratch;
  hydro::HydroPrimitiveCacheSoa m_primitive_cache;
  hydro::HydroPatchGeometry m_geometry;
  // Dense CellSoa row <-> physical Cartesian solver-row maps.
  std::vector<std::uint32_t> m_dense_row_by_geometry_row;
  std::vector<std::uint32_t> m_geometry_row_by_dense_row;
  std::vector<double> m_ordered_cell_accel_x;
  std::vector<double> m_ordered_cell_accel_y;
  std::vector<double> m_ordered_cell_accel_z;
  hydro::HydroProfileEvent m_last_hydro_profile{};
  core::HydroCflDiagnostics m_last_hydro_cfl_diagnostics{};
  SolverGhostRefreshReport m_last_ghost_refresh{};
  parallel::GhostCacheLifecycle m_ghost_cache_lifecycle{};
  std::vector<std::size_t> m_active_cells;
  std::vector<std::size_t> m_active_faces;
  std::optional<HydroGeometryCacheKey> m_cached_geometry_key;
};

bool maybeWriteOutputs(
    const core::FrozenConfig& frozen_config,
    const core::SimulationConfig& config,
    const core::SimulationState& state,
    const core::IntegratorState& integrator_state,
    const core::HierarchicalTimeBinScheduler& scheduler,
    const core::HierarchicalTimeBinScheduler& gas_cell_scheduler,
    const GravityStageCallback& gravity_callback,
    ReferenceWorkflowReport& report,
    core::ProfilerSession& profiler,
    bool write_outputs_enabled,
    bool snapshot_due,
    bool checkpoint_due) {
  if (!write_outputs_enabled || (!snapshot_due && !checkpoint_due)) {
    return false;
  }

#if !COSMOSIM_ENABLE_HDF5
  throw std::runtime_error(
      "runtime outputs requested, but this build lacks HDF5 support. Reconfigure with COSMOSIM_ENABLE_HDF5=ON.");
#else
  if (integrator_state.step_index == 0) {
    return false;
  }
  if (checkpoint_due) {
    core::assertCanWriteCheckpointAtBoundary(integrator_state, scheduler.currentTick());
    if (gas_cell_scheduler.currentTick() != scheduler.currentTick()) {
      throw std::runtime_error("output checkpoint requires particle and gas-cell schedulers on one integer tick");
    }
  }
  if (snapshot_due) {
    core::assertCanWriteSnapshotAtBoundary(integrator_state);
  }
  if ((snapshot_due || checkpoint_due) && !core::isOutputSafeBoundary(integrator_state.last_completed_boundary_kind)) {
    return false;
  }

  bool output_flushed = false;
  if (snapshot_due) {
    io::SnapshotWritePayload snapshot_payload;
    snapshot_payload.state = &state;
    snapshot_payload.config = &config;
    snapshot_payload.normalized_config_text = frozen_config.normalized_text;
    snapshot_payload.provenance =
        makeGravityAwareProvenanceRecord(frozen_config, config);
    report.snapshot_path = report.run_directory / formatIndexedRankedFileStem(config.output.output_stem, integrator_state.step_index, gravity_callback.runtimeTopology().world_size, gravity_callback.runtimeTopology().world_rank);
    io::writeGadgetArepoSnapshotHdf5(report.snapshot_path, snapshot_payload);
    report.snapshot_roundtrip_executed = true;
    const io::SnapshotReadResult snapshot_read = io::readGadgetArepoSnapshotHdf5(report.snapshot_path, config);
    report.snapshot_roundtrip_ok = snapshot_read.state.particles.size() == state.particles.size();
    profiler.recordEvent(core::RuntimeEvent{
        .event_kind = "snapshot.write.complete",
        .severity = report.snapshot_roundtrip_ok ? core::RuntimeEventSeverity::kInfo : core::RuntimeEventSeverity::kWarning,
        .subsystem = "io.snapshot",
        .step_index = integrator_state.step_index,
        .simulation_time_code = integrator_state.current_time_code,
        .scale_factor = integrator_state.current_scale_factor,
        .message = "snapshot output written and verified",
        .payload = {{"path", report.snapshot_path.string()}},
    });
    output_flushed = true;
  }

  if (checkpoint_due) {
    io::RestartWritePayload restart_payload;
    restart_payload.persistent_state.simulation_state = &state;
    restart_payload.integrator_state = &integrator_state;
    restart_payload.scheduler = &scheduler;
    restart_payload.gas_cell_scheduler = &gas_cell_scheduler;
    const io::GravityForceCachePersistentState gravity_force_cache =
        gravity_callback.exportRestartForceCache(state);
    restart_payload.gravity_force_cache = &gravity_force_cache;
    restart_payload.provenance =
        makeGravityAwareProvenanceRecord(frozen_config, config);
    restart_payload.normalized_config_text = frozen_config.normalized_text;
    restart_payload.normalized_config_hash_hex = frozen_config.provenance.config_hash_hex;
    restart_payload.distributed_gravity_state.schema_version = 2;
    restart_payload.distributed_gravity_state.decomposition_epoch = integrator_state.step_index;
    restart_payload.distributed_gravity_state.world_size = gravity_callback.runtimeTopology().world_size;
    restart_payload.distributed_gravity_state.pm_grid_nx = gravity_callback.pmGridShape().nx;
    restart_payload.distributed_gravity_state.pm_grid_ny = gravity_callback.pmGridShape().ny;
    restart_payload.distributed_gravity_state.pm_grid_nz = gravity_callback.pmGridShape().nz;
    restart_payload.distributed_gravity_state.pm_decomposition_mode =
        pmDecompositionModeName(config.numerics.treepm_pm_decomposition_mode);
    restart_payload.distributed_gravity_state.gravity_kick_opportunity =
        integrator_state.pm_sync_state.gravityKickOpportunity();
    restart_payload.distributed_gravity_state.pm_update_cadence_steps =
        integrator_state.pm_sync_state.cadenceSteps();
    restart_payload.distributed_gravity_state.long_range_field_version =
        integrator_state.pm_sync_state.fieldVersion();
    restart_payload.distributed_gravity_state.last_long_range_refresh_opportunity =
        integrator_state.pm_sync_state.lastRefreshOpportunity();
    restart_payload.distributed_gravity_state.long_range_field_built_step_index =
        integrator_state.pm_sync_state.lastRefreshStepIndex();
    restart_payload.distributed_gravity_state.long_range_field_built_scale_factor =
        integrator_state.pm_sync_state.lastRefreshScaleFactor();
    restart_payload.distributed_gravity_state.long_range_restart_policy = "deterministic_rebuild";
    restart_payload.output_cadence_state.output_enabled = write_outputs_enabled;
    restart_payload.output_cadence_state.write_restarts = config.output.write_restarts;
    restart_payload.output_cadence_state.snapshot_due = false;
    restart_payload.output_cadence_state.checkpoint_due = false;
    restart_payload.output_cadence_state.last_completed_step_index = integrator_state.step_index;
    restart_payload.output_cadence_state.snapshot_interval_steps =
        static_cast<std::uint64_t>(std::max(config.output.snapshot_interval_steps, 0));
    if (config.output.snapshot_interval_steps > 0) {
      const std::uint64_t interval = static_cast<std::uint64_t>(config.output.snapshot_interval_steps);
      restart_payload.output_cadence_state.next_snapshot_step_index =
          ((integrator_state.step_index / interval) + 1U) * interval;
    }
    restart_payload.output_cadence_state.snapshot_stem = config.output.output_stem;
    restart_payload.output_cadence_state.restart_stem = config.output.restart_stem;
    restart_payload.stochastic_state = buildStochasticPersistentState(
        config,
        integrator_state,
        static_cast<std::uint32_t>(std::max(gravity_callback.runtimeTopology().world_rank, 0)));
    restart_payload.distributed_gravity_state.owning_rank_by_item.reserve(state.particle_sidecar.owning_rank.size());
    for (const std::uint32_t owner : state.particle_sidecar.owning_rank) {
      restart_payload.distributed_gravity_state.owning_rank_by_item.push_back(static_cast<int>(owner));
    }
    const std::size_t world_size = static_cast<std::size_t>(gravity_callback.runtimeTopology().world_size);
    restart_payload.distributed_gravity_state.pm_slab_begin_x_by_rank.resize(world_size, 0);
    restart_payload.distributed_gravity_state.pm_slab_end_x_by_rank.resize(world_size, 0);
    for (std::size_t rank = 0; rank < world_size; ++rank) {
      const parallel::PmSlabRange range = parallel::pmOwnedXRangeForRank(
          gravity_callback.pmGridShape().nx,
          static_cast<int>(world_size),
          static_cast<int>(rank));
      restart_payload.distributed_gravity_state.pm_slab_begin_x_by_rank[rank] = range.begin_x;
      restart_payload.distributed_gravity_state.pm_slab_end_x_by_rank[rank] = range.end_x;
    }
    restart_payload.provenance.gravity_treepm_decomposition_epoch =
        restart_payload.distributed_gravity_state.decomposition_epoch;
    restart_payload.provenance.gravity_treepm_restart_world_size =
        restart_payload.distributed_gravity_state.world_size;
    restart_payload.provenance.gravity_treepm_restart_pm_grid =
        std::to_string(restart_payload.distributed_gravity_state.pm_grid_nx) + "x" +
        std::to_string(restart_payload.distributed_gravity_state.pm_grid_ny) + "x" +
        std::to_string(restart_payload.distributed_gravity_state.pm_grid_nz);
    restart_payload.provenance.gravity_treepm_restart_slab_signature =
        pmSlabSignature(restart_payload.distributed_gravity_state);
    restart_payload.provenance.gravity_treepm_restart_kick_opportunity =
        restart_payload.distributed_gravity_state.gravity_kick_opportunity;
    restart_payload.provenance.gravity_treepm_restart_field_version =
        restart_payload.distributed_gravity_state.long_range_field_version;
    restart_payload.provenance.gravity_treepm_long_range_restart_policy =
        restart_payload.distributed_gravity_state.long_range_restart_policy;

    report.restart_path = report.run_directory / formatIndexedRankedFileStem(config.output.restart_stem, integrator_state.step_index, gravity_callback.runtimeTopology().world_size, gravity_callback.runtimeTopology().world_rank);
    io::writeRestartCheckpointHdf5(report.restart_path, restart_payload);
    report.restart_roundtrip_executed = true;
    const io::RestartReadResult restart_read = io::readRestartCheckpointHdf5(report.restart_path);
    const auto compatibility = parallel::evaluateDistributedRestartCompatibility(
        restart_read.distributed_gravity_state,
        gravity_callback.runtimeTopology());
    const bool restart_rank_qualified_name =
        gravity_callback.runtimeTopology().world_size == 1 ||
        report.restart_path.filename().string().find("_rank") != std::string::npos;
    report.restart_roundtrip_ok = restartRuntimeStateExactlyEquivalent(restart_read.state, state) &&
        restart_read.integrator_state.pm_refresh_enabled == integrator_state.pm_refresh_enabled &&
        gravityForceCachesEquivalent(restart_read.gravity_force_cache, gravity_force_cache) &&
        restart_read.scheduler_state.current_tick == scheduler.currentTick() &&
        restart_read.distributed_gravity_state.owning_rank_by_item.size() == state.particle_sidecar.owning_rank.size() &&
        restart_read.distributed_gravity_state.owning_rank_by_item == restart_payload.distributed_gravity_state.owning_rank_by_item &&
        restart_read.distributed_gravity_state.pm_slab_begin_x_by_rank == restart_payload.distributed_gravity_state.pm_slab_begin_x_by_rank &&
        restart_read.distributed_gravity_state.pm_slab_end_x_by_rank == restart_payload.distributed_gravity_state.pm_slab_end_x_by_rank &&
        restart_read.distributed_gravity_state.gravity_kick_opportunity == restart_payload.distributed_gravity_state.gravity_kick_opportunity &&
        restart_read.distributed_gravity_state.pm_update_cadence_steps == restart_payload.distributed_gravity_state.pm_update_cadence_steps &&
        restart_read.distributed_gravity_state.long_range_field_version == restart_payload.distributed_gravity_state.long_range_field_version &&
        restart_read.distributed_gravity_state.last_long_range_refresh_opportunity ==
            restart_payload.distributed_gravity_state.last_long_range_refresh_opportunity &&
        restart_read.distributed_gravity_state.long_range_field_built_step_index ==
            restart_payload.distributed_gravity_state.long_range_field_built_step_index &&
        std::abs(
            restart_read.distributed_gravity_state.long_range_field_built_scale_factor -
            restart_payload.distributed_gravity_state.long_range_field_built_scale_factor) <= 1.0e-12 &&
        restart_read.distributed_gravity_state.long_range_restart_policy ==
            restart_payload.distributed_gravity_state.long_range_restart_policy &&
        restart_read.output_cadence_state.output_enabled == restart_payload.output_cadence_state.output_enabled &&
        restart_read.output_cadence_state.write_restarts == restart_payload.output_cadence_state.write_restarts &&
        restart_read.output_cadence_state.next_snapshot_step_index ==
            restart_payload.output_cadence_state.next_snapshot_step_index &&
        restart_read.output_cadence_state.snapshot_stem == restart_payload.output_cadence_state.snapshot_stem &&
        restart_read.output_cadence_state.restart_stem == restart_payload.output_cadence_state.restart_stem &&
        stochasticStatesEquivalent(restart_read.stochastic_state, restart_payload.stochastic_state) &&
        restart_rank_qualified_name &&
        compatibility.compatible();
    profiler.recordEvent(core::RuntimeEvent{
        .event_kind = "restart.write.complete",
        .severity = report.restart_roundtrip_ok ? core::RuntimeEventSeverity::kInfo : core::RuntimeEventSeverity::kWarning,
        .subsystem = "io.restart",
        .step_index = integrator_state.step_index,
        .simulation_time_code = integrator_state.current_time_code,
        .scale_factor = integrator_state.current_scale_factor,
        .message = "restart checkpoint written and verified",
        .payload = {{"path", report.restart_path.string()},
                    {"restart_schema", io::restartSchema().name},
                    {"restart_schema_version", std::to_string(io::restartSchema().version)},
                    {"boundary_kind", std::string(core::stepBoundaryKindName(integrator_state.last_completed_boundary_kind))},
                    {"scheduler_current_tick", std::to_string(scheduler.currentTick())},
                    {"scheduler_max_bin", std::to_string(scheduler.maxBin())},
                    {"pm_cadence_steps", std::to_string(integrator_state.pm_sync_state.cadenceSteps())},
                    {"pm_gravity_kick_opportunity", std::to_string(integrator_state.pm_sync_state.gravityKickOpportunity())},
                    {"pm_field_version", std::to_string(integrator_state.pm_sync_state.fieldVersion())},
                    {"pm_long_range_field_valid", integrator_state.pm_long_range_field_valid ? "true" : "false"},
                    {"output_next_snapshot_step_index", std::to_string(restart_payload.output_cadence_state.next_snapshot_step_index)},
                    {"stochastic_module_count", std::to_string(restart_payload.stochastic_state.modules.size())}},
    });
    profiler.recordEvent(core::RuntimeEvent{
        .event_kind = "restart.read.complete",
        .severity = report.restart_roundtrip_ok ? core::RuntimeEventSeverity::kInfo : core::RuntimeEventSeverity::kWarning,
        .subsystem = "io.restart",
        .step_index = restart_read.integrator_state.step_index,
        .simulation_time_code = restart_read.integrator_state.current_time_code,
        .scale_factor = restart_read.integrator_state.current_scale_factor,
        .message = "restart checkpoint read and validated",
        .payload = {{"path", report.restart_path.string()},
                    {"restart_schema", restart_read.diagnostics.restart_schema_name},
                    {"restart_schema_version", std::to_string(restart_read.diagnostics.restart_schema_version)},
                    {"boundary_kind", restart_read.diagnostics.last_completed_boundary_kind},
                    {"restart_safe", restart_read.diagnostics.restart_safe ? "true" : "false"},
                    {"scheduler_current_tick", std::to_string(restart_read.diagnostics.scheduler_current_tick)},
                    {"scheduler_max_bin", std::to_string(restart_read.diagnostics.scheduler_max_bin)},
                    {"scheduler_element_count", std::to_string(restart_read.diagnostics.scheduler_element_count)},
                    {"scheduler_active_count", std::to_string(restart_read.diagnostics.scheduler_active_count)},
                    {"scheduler_pending_transition_count", std::to_string(restart_read.diagnostics.scheduler_pending_transition_count)},
                    {"pm_cadence_steps", std::to_string(restart_read.diagnostics.pm_cadence_steps)},
                    {"pm_gravity_kick_opportunity", std::to_string(restart_read.diagnostics.pm_gravity_kick_opportunity)},
                    {"pm_field_version", std::to_string(restart_read.diagnostics.pm_field_version)},
                    {"pm_long_range_field_valid", restart_read.diagnostics.pm_long_range_field_valid ? "true" : "false"},
                    {"output_next_snapshot_step_index", std::to_string(restart_read.diagnostics.output_next_snapshot_step_index)},
                    {"stochastic_module_count", std::to_string(restart_read.diagnostics.stochastic_module_count)},
                    {"payload_hash_hex", restart_read.payload_hash_hex}},
    });
    output_flushed = true;
  }
  return output_flushed;
#endif
}

}  // namespace

ReferenceWorkflowRunner::ReferenceWorkflowRunner(core::FrozenConfig frozen_config)
    : m_frozen_config(std::move(frozen_config)) {}

const core::FrozenConfig& ReferenceWorkflowRunner::frozenConfig() const noexcept {
  return m_frozen_config;
}

ReferenceWorkflowReport ReferenceWorkflowRunner::run(const ReferenceWorkflowOptions& options) const {
  return runImpl(nullptr, options);
}

ReferenceWorkflowReport ReferenceWorkflowRunner::run(
    const std::filesystem::path& output_root_override,
    const ReferenceWorkflowOptions& options) const {
  return runImpl(&output_root_override, options);
}

struct PendingOutputBoundary {
  bool snapshot_due = false;
  bool checkpoint_due = false;
};

void latchOutputRequestForCompletedStep(
    const core::SimulationConfig& config,
    const ReferenceWorkflowOptions& options,
    std::uint64_t next_step_index,
    PendingOutputBoundary& pending) {
  if (!options.write_outputs || config.output.snapshot_interval_steps <= 0) {
    return;
  }
  if ((next_step_index % static_cast<std::uint64_t>(config.output.snapshot_interval_steps)) != 0U) {
    return;
  }
  pending.snapshot_due = true;
  pending.checkpoint_due = pending.checkpoint_due || config.output.write_restarts;
}

[[nodiscard]] core::StepBoundaryKind requestedBoundaryForPendingOutput(const PendingOutputBoundary& pending) {
  if (pending.checkpoint_due) {
    return core::StepBoundaryKind::kCheckpointPoint;
  }
  if (pending.snapshot_due) {
    return core::StepBoundaryKind::kSnapshotPoint;
  }
  return core::StepBoundaryKind::kGlobalSynchronizationPoint;
}

class OutputBoundaryCallback final : public core::IntegrationCallback {
 public:
  OutputBoundaryCallback(
      const core::FrozenConfig& frozen_config,
      const core::SimulationConfig& config,
      const core::HierarchicalTimeBinScheduler& scheduler,
      const core::HierarchicalTimeBinScheduler& gas_cell_scheduler,
      const GravityStageCallback& gravity_callback,
      ReferenceWorkflowReport& report,
      core::ProfilerSession& profiler,
      PendingOutputBoundary& pending_output,
      bool write_outputs_enabled)
      : m_frozen_config(frozen_config),
        m_config(config),
        m_scheduler(scheduler),
        m_gas_cell_scheduler(gas_cell_scheduler),
        m_gravity_callback(gravity_callback),
        m_report(report),
        m_profiler(profiler),
        m_pending_output(pending_output),
        m_write_outputs_enabled(write_outputs_enabled) {}

  [[nodiscard]] std::string_view callbackName() const override { return "output_boundary"; }
  [[nodiscard]] std::span<const core::IntegrationStage> integrationStages() const override {
    static constexpr std::array stages{core::IntegrationStage::kOutputCheck};
    return stages;
  }
  [[nodiscard]] std::span<const core::StageContract> stageContracts() const override { return m_contracts; }

  void onStage(core::StepContext& context) override {
    if (context.stage != core::IntegrationStage::kOutputCheck) {
      throw std::logic_error("output boundary handler received an unregistered stage");
    }
    if (!m_pending_output.snapshot_due && !m_pending_output.checkpoint_due) {
      return;
    }
    if (!context.boundary.output_safe || !context.boundary.restart_safe) {
      if (m_pending_output.checkpoint_due) {
        core::assertCanWriteCheckpointAtBoundary(context.integrator_state, m_scheduler.currentTick());
      }
      return;
    }
    const bool output_flushed = maybeWriteOutputs(
        m_frozen_config,
        m_config,
        context.state,
        context.integrator_state,
        m_scheduler,
        m_gas_cell_scheduler,
        m_gravity_callback,
        m_report,
        m_profiler,
        m_write_outputs_enabled,
        m_pending_output.snapshot_due,
        m_pending_output.checkpoint_due);
    if (output_flushed) {
      m_pending_output = {};
    }
  }

 private:
  static constexpr std::array<core::StageContract, 1> m_contracts{{
      {.stage = core::IntegrationStage::kOutputCheck,
       .required_inputs = core::StageDataDomain::kParticles | core::StageDataDomain::kGasCells | core::StageDataDomain::kOutputState | core::StageDataDomain::kRestartState,
       .mutated_state = core::StageDataDomain::kOutputState | core::StageDataDomain::kRestartState | core::StageDataDomain::kDiagnostics,
       .produced_outputs = core::StageDataDomain::kOutputState | core::StageDataDomain::kRestartState | core::StageDataDomain::kDiagnostics,
       .allowed_side_effects = core::StageDataDomain::kOutputState | core::StageDataDomain::kRestartState | core::StageDataDomain::kDiagnostics,
       .sync_requirements = core::StageSyncRequirement::kGlobal,
       .active_set_family = core::StageActiveSetFamily::kOutputState,
       .restart_safety = core::StageSafety::kSafe,
       .output_safety = core::StageSafety::kSafe,
       .owner = core::StageSubsystem::kOutput},
  }};

  const core::FrozenConfig& m_frozen_config;
  const core::SimulationConfig& m_config;
  const core::HierarchicalTimeBinScheduler& m_scheduler;
  const core::HierarchicalTimeBinScheduler& m_gas_cell_scheduler;
  const GravityStageCallback& m_gravity_callback;
  ReferenceWorkflowReport& m_report;
  core::ProfilerSession& m_profiler;
  PendingOutputBoundary& m_pending_output;
  bool m_write_outputs_enabled = false;
};

ReferenceWorkflowReport ReferenceWorkflowRunner::runImpl(
    const std::filesystem::path* output_root_override,
    const ReferenceWorkflowOptions& options) const {
  const core::SimulationConfig& config = m_frozen_config.config;

  ReferenceWorkflowReport report;
  report.run_directory = computeRunDirectory(config, output_root_override);
  report.config_compatible = true;
  report.schema_compatible =
      config.schema_version == 1 && io::gadgetArepoSchemaMap().schema_version >= 2 &&
      io::isRestartSchemaCompatible(io::restartSchema().version);

  core::ProfilerSession profiler(true);
  profiler.recordEvent(core::RuntimeEvent{
      .event_kind = "config.freeze",
      .severity = core::RuntimeEventSeverity::kInfo,
      .subsystem = "core.config",
      .step_index = options.step_index,
      .simulation_time_code = config.numerics.t_code_begin,
      .scale_factor = 1.0,
      .message = "frozen configuration accepted for runtime workflow",
      .payload = {{"schema_version", std::to_string(config.schema_version)},
                  {"config_hash_hex", m_frozen_config.provenance.config_hash_hex},
                  {"source_name", m_frozen_config.provenance.source_name}},
  });

  try {
    ensureRunDirectory(report.run_directory);
    core::writeNormalizedConfigSnapshot(m_frozen_config, report.run_directory);
    report.normalized_config_snapshot_path = report.run_directory / "normalized_config.param.txt";
    report.normalized_config_snapshot_written = true;

    if (!report.schema_compatible) {
      throw std::runtime_error("runtime workflow schema compatibility validation failed");
    }

#if !COSMOSIM_ENABLE_HDF5
    if (options.write_outputs) {
      throw std::runtime_error(
          "this build lacks HDF5 support, so the config-driven runtime cannot emit snapshots/restarts; reconfigure with COSMOSIM_ENABLE_HDF5=ON or use a no-output internal test path");
    }
#endif

    const core::ModePolicy mode_policy = core::buildModePolicy(config.mode);
    core::validateModePolicy(config, mode_policy);

    if (options.initial_state_override != nullptr && options.restart_state_override != nullptr) {
      throw std::invalid_argument(
          "ReferenceWorkflowOptions cannot provide both initial_state_override and restart_state_override");
    }

    io::IcReadResult ic_result;
    if (options.restart_state_override != nullptr) {
      // Restart payloads are authoritative for production continuation.  Do not
      // recreate time-bin or integrator state from configuration defaults.
      ic_result.state = options.restart_state_override->state;
      ic_result.report.defaulted_fields.push_back("restart_state_override=checkpoint_payload");
    } else if (options.initial_state_override != nullptr) {
      // A test/benchmark override still traverses every production workflow
      // boundary after import.  It exists so identity-decoupled states that
      // cannot be represented by legacy GADGET IC particle rows are testable
      // without adding a second execution path.
      ic_result.state = *options.initial_state_override;
      ic_result.report.defaulted_fields.push_back("initial_state_override=caller_supplied");
    } else {
      ic_result = loadInitialConditions(m_frozen_config);
    }
    core::SimulationState state = std::move(ic_result.state);
    profiler.setMemoryReport(core::collectSimulationMemoryReport(state));
    const core::MemoryReport* startup_memory_report = profiler.memoryReport();
    if (startup_memory_report != nullptr) {
      profiler.recordEvent(core::RuntimeEvent{
          .event_kind = "memory.startup_snapshot",
          .severity = core::RuntimeEventSeverity::kInfo,
          .subsystem = "core.memory",
          .step_index = options.step_index,
          .simulation_time_code = config.numerics.t_code_begin,
          .scale_factor = 1.0,
          .message = "startup memory accounting snapshot from owned host buffers",
          .payload = {{"persistent_total_bytes", std::to_string(startup_memory_report->totals.persistent_total_bytes)},
                      {"transient_total_bytes", std::to_string(startup_memory_report->totals.transient_total_bytes)},
                      {"unknown_external_allocations", "true"}},
      });
    }
    finalizeStateMetadata(m_frozen_config, state);
    maybeInitializeParticleSofteningFromSpeciesPolicy(state, config);
    const auto expected_global_identity = parallel::summarizeLocalOwnedParticleIds(state.particle_sidecar.particle_id);
    const std::vector<std::uint64_t> expected_global_particle_ids(
        state.particle_sidecar.particle_id.begin(), state.particle_sidecar.particle_id.end());
    parallel::MpiContext mpi_context;
    report.world_size = mpi_context.worldSize();
    report.world_rank = mpi_context.worldRank();
    seedParticleOwnershipFromPmSlabs(
        state,
        static_cast<std::size_t>(config.numerics.treepm_pm_grid_nx),
        mpi_context.worldSize(),
        config.cosmology.box_size_x_mpc_comoving);
    applyInitialGravityAwareDecomposition(state, config, mpi_context.worldSize(), mpi_context.worldRank(), &profiler);

    report.local_particle_count = static_cast<std::uint64_t>(state.particles.size());
    report.local_cell_count = static_cast<std::uint64_t>(state.cells.size());
    report.local_particle_id_sum = computeParticleIdSum(state);
    const std::uint64_t local_particle_id_square_sum = computeParticleIdSquareSum(state);
    report.local_particle_id_xor = computeParticleIdXor(state);
    const auto local_identity = parallel::summarizeLocalOwnedParticleIds(state.particle_sidecar.particle_id);
    report.local_particle_ids_unique = local_identity.local_particle_ids_unique;
    report.global_particle_count = mpi_context.allreduceSumUint64(report.local_particle_count);
    report.global_cell_count = mpi_context.allreduceSumUint64(report.local_cell_count);
    report.global_particle_id_sum = mpi_context.allreduceSumUint64(report.local_particle_id_sum);
    const std::uint64_t global_particle_id_square_sum =
        mpi_context.allreduceSumUint64(local_particle_id_square_sum);
    report.global_particle_id_xor = mpi_context.allreduceXorUint64(report.local_particle_id_xor);
    const bool all_ranks_have_unique_local_ids =
        mpi_context.allreduceSumUint64(report.local_particle_ids_unique ? 1ULL : 0ULL) ==
        static_cast<std::uint64_t>(mpi_context.worldSize());
    const parallel::LocalOwnershipIdentitySummary reduced_identity{
        .local_owned_count = report.global_particle_count,
        .local_particle_id_sum = report.global_particle_id_sum,
        .local_particle_id_square_sum = global_particle_id_square_sum,
        .local_particle_id_xor = report.global_particle_id_xor,
        .local_particle_ids_unique = all_ranks_have_unique_local_ids,
    };
    report.global_particle_partition_identity_match = parallel::partitionIdentityMatchesGeneratedSet(
        reduced_identity,
        expected_global_identity.local_owned_count,
        expected_global_identity.local_particle_id_sum,
        expected_global_identity.local_particle_id_square_sum,
        expected_global_identity.local_particle_id_xor);
    if (!report.global_particle_partition_identity_match) {
      throw std::runtime_error(
          "distributed ownership invariant failed after initial decomposition: duplicate, missing, or extra authoritative particle IDs");
    }
    profiler.recordEvent(core::RuntimeEvent{
        .event_kind = "parallel.ownership.partition_summary",
        .severity = core::RuntimeEventSeverity::kInfo,
        .subsystem = "parallel.ownership",
        .step_index = options.step_index,
        .simulation_time_code = config.numerics.t_code_begin,
        .scale_factor = 1.0,
        .message = "initial ownership-partition summary recorded",
        .payload = {{"world_size", std::to_string(report.world_size)},
                    {"world_rank", std::to_string(report.world_rank)},
                    {"local_particle_count", std::to_string(report.local_particle_count)},
                    {"global_particle_count", std::to_string(report.global_particle_count)},
                    {"local_cell_count", std::to_string(report.local_cell_count)},
                    {"global_cell_count", std::to_string(report.global_cell_count)},
                    {"local_particle_id_sum", std::to_string(report.local_particle_id_sum)},
                    {"global_particle_id_sum", std::to_string(report.global_particle_id_sum)},
                    {"local_particle_id_xor", std::to_string(report.local_particle_id_xor)},
                    {"global_particle_id_xor", std::to_string(report.global_particle_id_xor)},
                    {"local_particle_ids_unique", report.local_particle_ids_unique ? "true" : "false"},
                    {"global_particle_partition_identity_match",
                     report.global_particle_partition_identity_match ? "true" : "false"}},
    });

    core::CosmologyBackgroundConfig background_config;
    background_config.hubble_param = config.cosmology.hubble_param;
    background_config.omega_matter = config.cosmology.omega_matter;
    background_config.omega_lambda = config.cosmology.omega_lambda;
    const std::optional<core::LambdaCdmBackground> background = mode_policy.cosmological_comoving_frame
        ? std::optional<core::LambdaCdmBackground>(core::LambdaCdmBackground(background_config))
        : std::nullopt;

    const bool restoring_from_restart = options.restart_state_override != nullptr;
    core::HierarchicalTimeBinScheduler scheduler(
        static_cast<std::uint8_t>(std::max(0, std::min(config.numerics.hierarchical_max_rung, 12))));
    core::HierarchicalTimeBinScheduler gas_cell_scheduler(
        static_cast<std::uint8_t>(std::max(0, std::min(config.numerics.hierarchical_max_rung, 12))));

    const core::UnitSystem runtime_units = core::makeUnitSystem(
        config.units.length_unit,
        config.units.mass_unit,
        config.units.velocity_unit);

    core::IntegratorState integrator_state;
    if (restoring_from_restart) {
      const io::RestartReadResult& restart = *options.restart_state_override;
      scheduler.importPersistentState(restart.scheduler_state);
      gas_cell_scheduler.importPersistentState(restart.gas_cell_scheduler_state);
      integrator_state = restart.integrator_state;
      // v14-v19 checkpoints predate serialization of this frozen workflow
      // policy bit.  Preserve their established compatibility behavior while
      // v20+ resumes use the checkpoint-authoritative value verbatim.
      if (restart.diagnostics.restart_schema_version < io::restartSchema().version) {
        integrator_state.pm_refresh_enabled = true;
      }
      if (scheduler.elementCount() != state.particles.size()) {
        throw std::runtime_error(
            "ReferenceWorkflow restart payload particle scheduler coverage does not match SimulationState");
      }
      if (gas_cell_scheduler.elementCount() != state.cells.size()) {
        throw std::runtime_error(
            "ReferenceWorkflow restart payload gas-cell scheduler coverage does not match SimulationState");
      }
      state.requireGasCellIdentityMapCoversDenseRows("ReferenceWorkflow restart resume");
    } else {
      initializeSchedulerBins(state, scheduler, gas_cell_scheduler);
      integrator_state.step_index = options.step_index;
      integrator_state.current_time_code = config.numerics.t_code_begin;
      integrator_state.time_si_per_code = runtime_units.timeSiPerCode();
      integrator_state.current_scale_factor = background.has_value()
          ? config.numerics.a_begin
          : 1.0;
      integrator_state.current_redshift = (integrator_state.current_scale_factor > 0.0)
          ? (1.0 / integrator_state.current_scale_factor - 1.0)
          : 0.0;
      integrator_state.current_hubble_rate_code = background.has_value()
          ? background->hubbleSi(integrator_state.current_scale_factor) * integrator_state.time_si_per_code
          : 0.0;
      integrator_state.dt_time_code = options.dt_time_code > 0.0
          ? options.dt_time_code
          : std::max(
                1.0e-6,
                (config.numerics.t_code_end - config.numerics.t_code_begin) /
                    static_cast<double>(std::max(config.numerics.max_global_steps, 1)));
      integrator_state.time_bins.hierarchical_enabled = true;
      integrator_state.time_bins.max_bin = scheduler.maxBin();
      integrator_state.pm_refresh_enabled = true;
      integrator_state.pm_sync_state.reset(static_cast<std::uint64_t>(std::max(config.numerics.treepm_update_cadence_steps, 1)));
      for (std::size_t particle_index = 0; particle_index < state.particles.size(); ++particle_index) {
        state.particle_sidecar.last_drift_time_code[particle_index] = integrator_state.current_time_code;
        state.particle_sidecar.last_drift_scale_factor[particle_index] = integrator_state.current_scale_factor;
      }
    }
    syncTimeBinsFromSchedulers(scheduler, gas_cell_scheduler, state);
    PendingOutputBoundary pending_output;

    core::StepOrchestrator orchestrator;
    StageAuditCallback stage_audit(&report.stage_sequence);
    DriftCallback drift_callback;
    const std::filesystem::path zoom_region_path =
        config.mode.zoom_region_file.empty()
        ? std::filesystem::path{}
        : resolveConfigRelativePath(m_frozen_config, std::filesystem::path(config.mode.zoom_region_file));
    GravityStageCallback gravity_callback(config, mode_policy, zoom_region_path);
    if (restoring_from_restart) {
      const io::RestartReadResult& restart = *options.restart_state_override;
      if (restart.gravity_force_cache.valid) {
        gravity_callback.importRestartForceCache(restart.gravity_force_cache, state);
      } else {
        // A pre-v20 checkpoint has no cache that can safely service the next
        // KDK pre-kick.  Do not pretend its persisted PM validity bit implies
        // a live force cache in this process; force the legal bootstrap path.
        integrator_state.pm_long_range_field_valid = false;
      }
    }
    report.treepm_pm_grid = gravity_callback.pmGridSize();
    report.treepm_pm_grid_nx = gravity_callback.pmGridShape().nx;
    report.treepm_pm_grid_ny = gravity_callback.pmGridShape().ny;
    report.treepm_pm_grid_nz = gravity_callback.pmGridShape().nz;
    report.treepm_pm_grid_shape = std::to_string(report.treepm_pm_grid_nx) + "x" +
        std::to_string(report.treepm_pm_grid_ny) + "x" + std::to_string(report.treepm_pm_grid_nz);
    report.treepm_update_cadence_steps = static_cast<int>(integrator_state.pm_sync_state.cadenceSteps());
    profiler.recordEvent(core::RuntimeEvent{
        .event_kind = "gravity.treepm_setup",
        .severity = core::RuntimeEventSeverity::kInfo,
        .subsystem = "gravity.treepm",
        .step_index = integrator_state.step_index,
        .simulation_time_code = integrator_state.current_time_code,
        .scale_factor = integrator_state.current_scale_factor,
        .message = "TreePM runtime configuration initialized",
        .payload = {
            {"pm_grid", std::to_string(config.numerics.treepm_pm_grid_nx) + "x" +
                    std::to_string(config.numerics.treepm_pm_grid_ny) + "x" +
                    std::to_string(config.numerics.treepm_pm_grid_nz)},
            {"pm_assignment_scheme", treePmAssignmentSchemeName(config.numerics.treepm_assignment_scheme)},
            {"pm_window_deconvolution", config.numerics.treepm_enable_window_deconvolution ? "true" : "false"},
            {"asmth_cells", std::to_string(config.numerics.treepm_asmth_cells)},
            {"rcut_cells", std::to_string(config.numerics.treepm_rcut_cells)},
            {"mesh_spacing_x_mpc_comoving", std::to_string(config.cosmology.box_size_x_mpc_comoving / static_cast<double>(config.numerics.treepm_pm_grid_nx))},
            {"mesh_spacing_y_mpc_comoving", std::to_string(config.cosmology.box_size_y_mpc_comoving / static_cast<double>(config.numerics.treepm_pm_grid_ny))},
            {"mesh_spacing_z_mpc_comoving", std::to_string(config.cosmology.box_size_z_mpc_comoving / static_cast<double>(config.numerics.treepm_pm_grid_nz))},
            {"split_scale_mpc_comoving", std::to_string(config.numerics.treepm_asmth_cells *
                std::cbrt((config.cosmology.box_size_x_mpc_comoving / static_cast<double>(config.numerics.treepm_pm_grid_nx)) *
                          (config.cosmology.box_size_y_mpc_comoving / static_cast<double>(config.numerics.treepm_pm_grid_ny)) *
                          (config.cosmology.box_size_z_mpc_comoving / static_cast<double>(config.numerics.treepm_pm_grid_nz))))},
            {"cutoff_radius_mpc_comoving", std::to_string(config.numerics.treepm_rcut_cells *
                std::cbrt((config.cosmology.box_size_x_mpc_comoving / static_cast<double>(config.numerics.treepm_pm_grid_nx)) *
                          (config.cosmology.box_size_y_mpc_comoving / static_cast<double>(config.numerics.treepm_pm_grid_ny)) *
                          (config.cosmology.box_size_z_mpc_comoving / static_cast<double>(config.numerics.treepm_pm_grid_nz))))},
            {"pm_update_cadence_steps", std::to_string(config.numerics.treepm_update_cadence_steps)},
            {"zoom_long_range_strategy",
                config.mode.zoom_long_range_strategy == core::ZoomLongRangeStrategy::kDisabled
                    ? "disabled"
                    : "global_coarse_plus_focused_highres_correction"},
            {"zoom_region_center_x_mpc_comoving", std::to_string(config.mode.zoom_region_center_x_mpc_comoving)},
            {"zoom_region_center_y_mpc_comoving", std::to_string(config.mode.zoom_region_center_y_mpc_comoving)},
            {"zoom_region_center_z_mpc_comoving", std::to_string(config.mode.zoom_region_center_z_mpc_comoving)},
            {"zoom_region_radius_mpc_comoving", std::to_string(config.mode.zoom_region_radius_mpc_comoving)},
            {"zoom_focused_pm_grid", std::to_string(config.mode.zoom_focused_pm_grid_nx) + "x" +
                std::to_string(config.mode.zoom_focused_pm_grid_ny) + "x" +
                std::to_string(config.mode.zoom_focused_pm_grid_nz)},
            {"zoom_contamination_radius_mpc_comoving", std::to_string(config.mode.zoom_contamination_radius_mpc_comoving)},
            {"softening_policy", describeSofteningPolicy(config, &state)},
            {"softening_kernel", "plummer"},
            {"softening_epsilon_kpc_comoving", std::to_string(config.numerics.gravity_softening_kpc_comoving)},
            {"pm_fft_backend", gravity::PmSolver::fftBackendName()},
        },
    });
    {
      const std::array startup_reports{
          core::collectSimulationMemoryReport(state),
          gravity_callback.memoryReport()};
      profiler.setMemoryReport(core::mergeMemoryReports(startup_reports));
      const core::MemoryReport* runtime_memory_report = profiler.memoryReport();
      if (runtime_memory_report != nullptr) {
        profiler.recordEvent(core::RuntimeEvent{
            .event_kind = "memory.runtime_startup_snapshot",
            .severity = core::RuntimeEventSeverity::kInfo,
            .subsystem = "core.memory",
            .step_index = integrator_state.step_index,
            .simulation_time_code = integrator_state.current_time_code,
            .scale_factor = integrator_state.current_scale_factor,
            .message = "startup memory accounting snapshot including gravity solver workspaces",
            .payload = {{"persistent_total_bytes", std::to_string(runtime_memory_report->totals.persistent_total_bytes)},
                        {"transient_total_bytes", std::to_string(runtime_memory_report->totals.transient_total_bytes)},
                        {"unknown_total_bytes", std::to_string(runtime_memory_report->totals.unknown_total_bytes)}},
        });
      }
    }
    HydroStageCallback hydro_callback(config, mode_policy, gravity_callback);
    analysis::DiagnosticsCallback diagnostics_callback(config);
    physics::StarFormationCallback star_formation_callback(
        physics::StarFormationModel(makeRuntimeStarFormationConfig(config.physics, runtime_units)),
        static_cast<std::uint32_t>(std::max(mpi_context.worldRank(), 0)));
    physics::BlackHoleAgnCallback bh_callback(
        physics::BlackHoleAgnModel(makeRuntimeBlackHoleAgnConfig(config.physics, runtime_units)));
    OutputBoundaryCallback output_boundary_callback(
        m_frozen_config,
        config,
        scheduler,
        gas_cell_scheduler,
        gravity_callback,
        report,
        profiler,
        pending_output,
        options.write_outputs);

    orchestrator.registerCallback(stage_audit);
    orchestrator.registerCallback(drift_callback);
    orchestrator.registerCallback(gravity_callback);
    orchestrator.registerCallback(hydro_callback);
    orchestrator.registerCallback(star_formation_callback);
    orchestrator.registerCallback(bh_callback);
    orchestrator.registerCallback(diagnostics_callback);
    orchestrator.registerCallback(output_boundary_callback);

    if (!restoring_from_restart) {
      updateAdaptiveTimeBins(
          state,
          scheduler,
          integrator_state,
          config,
          {},
          {},
          {},
          {},
          {},
          {});
      updateGasCellAdaptiveTimeBins(
          state,
          gas_cell_scheduler,
          integrator_state,
          config,
          {},
          {},
          {},
          {},
          {},
          {});
    }
    ensureSchedulersCoverState(state, scheduler, gas_cell_scheduler);

    const std::uint64_t run_start_step_index = integrator_state.step_index;
    const std::uint64_t configured_segment_steps = options.max_steps_override > 0
        ? options.max_steps_override
        : static_cast<std::uint64_t>(std::max(config.numerics.max_global_steps, 0));
    const std::uint64_t target_step_index = integrator_state.step_index + configured_segment_steps;
    while (integrator_state.step_index < target_step_index &&
           integrator_state.current_time_code < config.numerics.t_code_end) {
      const std::span<const std::uint32_t> active_particle_scheduler_elements = scheduler.beginSubstep();
      const std::span<const std::uint32_t> active_gas_cell_scheduler_elements = gas_cell_scheduler.beginSubstep();
      if (scheduler.currentTick() != gas_cell_scheduler.currentTick()) {
        throw std::runtime_error("particle and gas-cell schedulers lost their shared integer timeline");
      }
      std::vector<std::uint32_t> active_particles(
          active_particle_scheduler_elements.begin(), active_particle_scheduler_elements.end());
      std::vector<std::uint32_t> active_cells(
          active_gas_cell_scheduler_elements.begin(), active_gas_cell_scheduler_elements.end());
      if (active_particles.empty() && active_cells.empty()) {
        scheduler.endSubstep();
        gas_cell_scheduler.endSubstep();
        continue;
      }

      core::TransientStepWorkspace workspace;
      latchOutputRequestForCompletedStep(
          config,
          options,
          integrator_state.step_index + 1U,
          pending_output);
      const core::StepBoundaryKind requested_boundary = requestedBoundaryForPendingOutput(pending_output);
      orchestrator.executeSchedulerSubstep(
          state,
          integrator_state,
          scheduler,
          active_particles,
          active_cells,
          background.has_value() ? &background.value() : nullptr,
          &workspace,
          &mode_policy,
          &profiler,
          requested_boundary);

      {
        const std::array runtime_reports{
            core::collectSimulationMemoryReport(state, &workspace),
            gravity_callback.memoryReport()};
        profiler.setMemoryReport(core::mergeMemoryReports(runtime_reports));
      }
      state.metadata.step_index = integrator_state.step_index;
      state.metadata.scale_factor = integrator_state.current_scale_factor;
      updateAdaptiveTimeBins(
          state,
          scheduler,
          integrator_state,
          config,
          gravity_callback.particleAccelX(),
          gravity_callback.particleAccelY(),
          gravity_callback.particleAccelZ(),
          gravity_callback.cellAccelX(),
          gravity_callback.cellAccelY(),
          gravity_callback.cellAccelZ());
      updateGasCellAdaptiveTimeBins(
          state,
          gas_cell_scheduler,
          integrator_state,
          config,
          gravity_callback.particleAccelX(),
          gravity_callback.particleAccelY(),
          gravity_callback.particleAccelZ(),
          gravity_callback.cellAccelX(),
          gravity_callback.cellAccelY(),
          gravity_callback.cellAccelZ());
      ensureSchedulersCoverState(state, scheduler, gas_cell_scheduler);
      scheduler.endSubstep();
      gas_cell_scheduler.endSubstep();

      parallel::DecompositionRuntimeMeasurements rebalance_measurements =
          gravity_callback.lastRuntimeDecompositionMeasurements();
      const hydro::HydroProfileEvent& hydro_profile = hydro_callback.lastHydroProfile();
      const SolverGhostRefreshReport& hydro_ghost_report = hydro_callback.lastGhostRefreshReport();
      rebalance_measurements.hydro_face_fluxes_recent = hydro_profile.face_count;
      rebalance_measurements.hydro_wall_ms_recent = hydro_profile.total_ms;
      rebalance_measurements.ghost_exchange_bytes_recent +=
          hydro_ghost_report.sent_bytes + hydro_ghost_report.received_bytes;
      rebalance_measurements.has_measurements = rebalance_measurements.has_measurements ||
          hydro_profile.face_count > 0 || hydro_profile.total_ms > 0.0 ||
          hydro_ghost_report.sent_bytes > 0 || hydro_ghost_report.received_bytes > 0;
      applyMeasuredRuntimeRebalancePlan(
          state,
          scheduler,
          config,
          mpi_context,
          mpi_context.worldRank(),
          rebalance_measurements,
          active_particles,
          expected_global_particle_ids,
          &profiler,
          integrator_state.step_index);
      ensureSchedulersCoverState(state, scheduler, gas_cell_scheduler);
      syncTimeBinsFromSchedulers(scheduler, gas_cell_scheduler, state);
      orchestrator.executeOutputBoundary(
          state,
          integrator_state,
          &profiler,
          requested_boundary);
    }

    report.completed_steps = integrator_state.step_index - run_start_step_index;
    report.final_hydro_cfl_diagnostics = hydro_callback.lastHydroCflDiagnostics();
    report.final_state_digest = computeStateDigest(state, integrator_state);
    report.local_particle_count = static_cast<std::uint64_t>(state.particles.size());
    report.local_cell_count = static_cast<std::uint64_t>(state.cells.size());
    report.local_particle_id_sum = computeParticleIdSum(state);
    const std::uint64_t final_local_particle_id_square_sum = computeParticleIdSquareSum(state);
    report.local_particle_id_xor = computeParticleIdXor(state);
    const auto final_local_identity = parallel::summarizeLocalOwnedParticleIds(state.particle_sidecar.particle_id);
    report.local_particle_ids_unique = final_local_identity.local_particle_ids_unique;
    report.global_particle_count = mpi_context.allreduceSumUint64(report.local_particle_count);
    report.global_cell_count = mpi_context.allreduceSumUint64(report.local_cell_count);
    report.global_particle_id_sum = mpi_context.allreduceSumUint64(report.local_particle_id_sum);
    const std::uint64_t final_global_particle_id_square_sum =
        mpi_context.allreduceSumUint64(final_local_particle_id_square_sum);
    report.global_particle_id_xor = mpi_context.allreduceXorUint64(report.local_particle_id_xor);
    const bool final_all_ranks_have_unique_local_ids =
        mpi_context.allreduceSumUint64(report.local_particle_ids_unique ? 1ULL : 0ULL) ==
        static_cast<std::uint64_t>(mpi_context.worldSize());
    const parallel::LocalOwnershipIdentitySummary final_reduced_identity{
        .local_owned_count = report.global_particle_count,
        .local_particle_id_sum = report.global_particle_id_sum,
        .local_particle_id_square_sum = final_global_particle_id_square_sum,
        .local_particle_id_xor = report.global_particle_id_xor,
        .local_particle_ids_unique = final_all_ranks_have_unique_local_ids,
    };
    report.global_particle_partition_identity_match = parallel::partitionIdentityMatchesGeneratedSet(
        final_reduced_identity,
        expected_global_identity.local_owned_count,
        expected_global_identity.local_particle_id_sum,
        expected_global_identity.local_particle_id_square_sum,
        expected_global_identity.local_particle_id_xor);
    if (!report.global_particle_partition_identity_match) {
      throw std::runtime_error(
          "distributed ownership invariant failed at run completion: duplicate, missing, or extra authoritative particle IDs");
    }
    report.treepm_long_range_refresh_count = gravity_callback.longRangeRefreshCount();
    report.treepm_long_range_reuse_count = gravity_callback.longRangeReuseCount();
    report.treepm_cadence_records.assign(
        gravity_callback.cadenceRecords().begin(),
        gravity_callback.cadenceRecords().end());
    report.canonical_stage_order = repeatedCanonicalOrder(report.stage_sequence);
    profiler.recordEvent(core::RuntimeEvent{
        .event_kind = "run.complete",
        .severity = core::RuntimeEventSeverity::kInfo,
        .subsystem = "workflows.reference",
        .step_index = integrator_state.step_index,
        .simulation_time_code = integrator_state.current_time_code,
        .scale_factor = integrator_state.current_scale_factor,
        .message = "runtime workflow completed",
        .payload = {{"completed_steps", std::to_string(report.completed_steps)},
                    {"run_directory", report.run_directory.string()},
                    {"treepm_update_cadence_steps", std::to_string(report.treepm_update_cadence_steps)},
                    {"treepm_long_range_refresh_count", std::to_string(report.treepm_long_range_refresh_count)},
                    {"treepm_long_range_reuse_count", std::to_string(report.treepm_long_range_reuse_count)},
                    {"local_particle_count", std::to_string(report.local_particle_count)},
                    {"global_particle_count", std::to_string(report.global_particle_count)},
                    {"local_particle_id_sum", std::to_string(report.local_particle_id_sum)},
                    {"global_particle_id_sum", std::to_string(report.global_particle_id_sum)},
                    {"local_particle_id_xor", std::to_string(report.local_particle_id_xor)},
                    {"global_particle_id_xor", std::to_string(report.global_particle_id_xor)},
                    {"local_particle_ids_unique", report.local_particle_ids_unique ? "true" : "false"},
                    {"global_particle_partition_identity_match",
                     report.global_particle_partition_identity_match ? "true" : "false"},
                    {"final_state_digest", std::to_string(report.final_state_digest)}},
    });
    flushCommonArtifacts(m_frozen_config, profiler, report);
    return report;
  } catch (const std::exception& ex) {
    profiler.recordEvent(core::RuntimeEvent{
        .event_kind = "run.failure",
        .severity = core::RuntimeEventSeverity::kFatal,
        .subsystem = "workflows.reference",
        .step_index = std::nullopt,
        .simulation_time_code = std::nullopt,
        .scale_factor = std::nullopt,
        .message = "runtime workflow failed",
        .payload = {{"error", ex.what()}, {"run_directory", report.run_directory.string()}},
    });
    try {
      flushCommonArtifacts(m_frozen_config, profiler, report);
    } catch (...) {
    }
    throw;
  }
}

}  // namespace cosmosim::workflows
