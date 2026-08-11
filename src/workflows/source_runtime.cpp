#include "cosmosim/workflows/source_runtime.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <numeric>
#include <span>
#include <stdexcept>
#include <unordered_set>
#include <vector>

#include "cosmosim/physics/black_hole_agn.hpp"
#include "cosmosim/physics/effective_multiphase_ism.hpp"
#include "cosmosim/parallel/distributed_memory.hpp"
#include "cosmosim/physics/star_formation.hpp"
#include "cosmosim/physics/stellar_evolution.hpp"
#include "cosmosim/physics/stellar_feedback.hpp"
#include "cosmosim/physics/metal_diffusion.hpp"
#include "workflows/internal/gas_cell_ownership.hpp"
#include "workflows/internal/runtime_stage_resource_access.hpp"
#include "workflows/internal/star_formation_geometry.hpp"

#if COSMOSIM_ENABLE_MPI
#include <mpi.h>
#endif

namespace cosmosim::workflows {
namespace {

[[nodiscard]] double newtonGCodeFromUnits(const core::UnitSystem& units) {
  return core::newtonGravitationalConstantCode(units);
}

[[nodiscard]] physics::StarFormationConfig makeRuntimeStarFormationConfig(
    const core::SimulationConfig& config,
    const core::UnitSystem& units) {
  physics::StarFormationConfig result = physics::makeStarFormationConfig(config.physics);
  result.newton_g_code = newtonGCodeFromUnits(units);
  result.density_is_comoving = config.units.coordinate_frame == core::CoordinateFrame::kComoving;
  result.geometry_is_comoving = config.units.coordinate_frame == core::CoordinateFrame::kComoving;
  return result;
}

[[nodiscard]] physics::StellarFeedbackConfig makeRuntimeStellarFeedbackConfig(
    const core::SimulationConfig& config,
    const core::UnitSystem& units) {
  physics::StellarFeedbackConfig result =
      physics::makeStellarFeedbackConfig(config.physics);
  // Mass and metal return is part of stellar evolution even when energetic
  // feedback is disabled. Keep the deposition transaction enabled and zero only
  // the energy/momentum coupling in that case.
  result.enabled = config.physics.enable_stellar_evolution;
  if (!config.physics.enable_feedback) {
    result.epsilon_thermal = 0.0;
    result.epsilon_kinetic = 0.0;
    result.epsilon_momentum = 0.0;
  }
  constexpr double k_joule_per_erg = 1.0e-7;
  result.total_energy_code_per_erg = k_joule_per_erg /
      (units.mass_si_per_code * units.velocity_si_per_code *
       units.velocity_si_per_code);
  return result;
}

[[nodiscard]] std::shared_ptr<const physics::EffectiveMultiphaseEosTable>
makeRuntimeEffectiveEosTable(
    const core::SimulationConfig& config,
    const core::UnitSystem& units) {
  if (!config.physics.enable_star_formation ||
      config.physics.star_formation_model !=
          core::StarFormationModelKind::kEffectiveMultiphaseTngLike) {
    return nullptr;
  }
  return std::make_shared<const physics::EffectiveMultiphaseEosTable>(
      physics::makeEffectiveMultiphaseEosConfig(config.physics),
      units,
      physics::makeEffectiveIsmReferenceCoolingProvider(config.physics));
}

[[nodiscard]] std::vector<std::uint64_t> allGatherUint64(
    const parallel::MpiContext& mpi_context,
    std::span<const std::uint64_t> local_values) {
  if (!mpi_context.isEnabled()) {
    return std::vector<std::uint64_t>(local_values.begin(), local_values.end());
  }
#if COSMOSIM_ENABLE_MPI
  const int world_size = mpi_context.worldSize();
  const std::uint64_t local_count = static_cast<std::uint64_t>(local_values.size());
  std::vector<std::uint64_t> counts(static_cast<std::size_t>(world_size), 0U);
  MPI_Allgather(
      const_cast<std::uint64_t*>(&local_count), 1, MPI_UINT64_T,
      counts.data(), 1, MPI_UINT64_T, MPI_COMM_WORLD);
  std::vector<int> recv_counts(static_cast<std::size_t>(world_size), 0);
  std::vector<int> displacements(static_cast<std::size_t>(world_size), 0);
  std::uint64_t total = 0U;
  for (int rank = 0; rank < world_size; ++rank) {
    if (counts[static_cast<std::size_t>(rank)] >
        static_cast<std::uint64_t>(std::numeric_limits<int>::max())) {
      throw std::overflow_error("particle-ID precommit count exceeds MPI int range");
    }
    recv_counts[static_cast<std::size_t>(rank)] =
        static_cast<int>(counts[static_cast<std::size_t>(rank)]);
    if (rank > 0) {
      displacements[static_cast<std::size_t>(rank)] =
          displacements[static_cast<std::size_t>(rank - 1)] +
          recv_counts[static_cast<std::size_t>(rank - 1)];
    }
    total += counts[static_cast<std::size_t>(rank)];
  }
  std::vector<std::uint64_t> gathered(static_cast<std::size_t>(total), 0U);
  MPI_Allgatherv(
      const_cast<std::uint64_t*>(local_values.data()),
      static_cast<int>(local_values.size()), MPI_UINT64_T,
      gathered.data(), recv_counts.data(), displacements.data(),
      MPI_UINT64_T, MPI_COMM_WORLD);
  return gathered;
#else
  throw std::runtime_error("distributed particle-ID precommit requires an MPI build");
#endif
}

class DistributedParticleIdRegistry final : public physics::ParticleIdPrecommit {
 public:
  explicit DistributedParticleIdRegistry(const parallel::MpiContext& mpi_context)
      : m_mpi_context(mpi_context) {}

  [[nodiscard]] std::vector<std::uint64_t> precommit(
      const core::SimulationState& state,
      std::span<const std::uint64_t> birth_keys) override {
    if (!m_initialized) {
      m_occupied = allGatherUint64(m_mpi_context, state.particle_sidecar.particle_id);
      std::sort(m_occupied.begin(), m_occupied.end());
      if ((!m_occupied.empty() && m_occupied.front() == 0U) ||
          std::adjacent_find(m_occupied.begin(), m_occupied.end()) != m_occupied.end()) {
        throw std::runtime_error(
            "ParticleIdRegistry: zero or duplicate existing ID during distributed initialization");
      }
      m_initialized = true;
    }

    std::vector<std::uint64_t> global_birth_keys = allGatherUint64(m_mpi_context, birth_keys);
    std::sort(global_birth_keys.begin(), global_birth_keys.end());
    if ((!global_birth_keys.empty() && global_birth_keys.front() == 0U) ||
        std::adjacent_find(global_birth_keys.begin(), global_birth_keys.end()) !=
            global_birth_keys.end()) {
      throw std::runtime_error(
          "ParticleIdRegistry: duplicate immutable birth key across owner ranks before mutation");
    }

    const std::vector<std::uint64_t> global_ids =
        physics::precommitStarParticleIdsExact(m_occupied, global_birth_keys);
    const std::size_t occupied_before = m_occupied.size();
    m_occupied.insert(m_occupied.end(), global_ids.begin(), global_ids.end());
    std::inplace_merge(
        m_occupied.begin(), m_occupied.begin() + static_cast<std::ptrdiff_t>(occupied_before),
        m_occupied.end());
    if (std::adjacent_find(m_occupied.begin(), m_occupied.end()) != m_occupied.end()) {
      throw std::runtime_error(
          "ParticleIdRegistry: exact distributed precommit produced a duplicate ID");
    }

    std::vector<std::uint64_t> local_ids;
    local_ids.reserve(birth_keys.size());
    for (const std::uint64_t birth_key : birth_keys) {
      const auto it = std::lower_bound(
          global_birth_keys.begin(), global_birth_keys.end(), birth_key);
      if (it == global_birth_keys.end() || *it != birth_key) {
        throw std::runtime_error(
            "ParticleIdRegistry: local birth key missing from global precommit");
      }
      const std::size_t index = static_cast<std::size_t>(it - global_birth_keys.begin());
      local_ids.push_back(global_ids[index]);
    }
    return local_ids;
  }

 private:
  const parallel::MpiContext& m_mpi_context;
  bool m_initialized = false;
  std::vector<std::uint64_t> m_occupied;
};

[[nodiscard]] physics::BlackHoleAgnConfig makeRuntimeBlackHoleAgnConfig(
    const core::PhysicsConfig& physics_config,
    const core::UnitSystem& units) {
  physics::BlackHoleAgnConfig config = physics::makeBlackHoleAgnConfig(physics_config);
  config.proton_mass_code = config.proton_mass_si / units.mass_si_per_code;
  config.thomson_cross_section_code = config.thomson_cross_section_si /
      (units.length_si_per_code * units.length_si_per_code);
  config.newton_g_code = config.newton_g_si * units.mass_si_per_code *
      units.timeSiPerCode() * units.timeSiPerCode() /
      (units.length_si_per_code * units.length_si_per_code * units.length_si_per_code);
  config.speed_of_light_code = config.speed_of_light_si / units.velocity_si_per_code;
  return config;
}

using PatchCellGeometry = internal::StarFormationPatchCellGeometry;
using internal::starFormationDerivativeAtCell;
using internal::starFormationPatchCellGeometry;
using internal::starFormationPatchIsLeaf;

class SourceRuntimeImpl final : public SourceRuntime {
 public:
  SourceRuntimeImpl(
      const core::SimulationConfig& config,
      const core::UnitSystem& units,
      std::uint32_t world_rank,
      const parallel::MpiContext& mpi_context)
      : m_units(units),
        m_effective_eos_table(makeRuntimeEffectiveEosTable(config, units)),
        m_star_formation(makeRuntimeStarFormationConfig(config, units), m_effective_eos_table),
        m_stellar_evolution(
            physics::makeStellarEvolutionConfig(config.physics),
            physics::loadStellarEvolutionTable(config.physics)),
        m_stellar_feedback(makeRuntimeStellarFeedbackConfig(config, units)),
        m_metal_diffusion(physics::makeMetalDiffusionConfig(config.physics)),
        m_black_hole(makeRuntimeBlackHoleAgnConfig(config.physics, units)),
        m_particle_id_registry(mpi_context),
        m_world_rank(world_rank),
        m_coordinate_frame(config.units.coordinate_frame),
        m_is_cosmological(
            config.mode.mode == core::SimulationMode::kCosmoCube ||
            config.mode.mode == core::SimulationMode::kZoomIn) {
    const double hubble_si = config.cosmology.hubble_param *
        core::constants::k_hubble_100_km_s_mpc_si;
    const double hubble_code = hubble_si * units.timeSiPerCode();
    const double g_code = newtonGCodeFromUnits(units);
    m_mean_baryon_density0_code = 3.0 * hubble_code * hubble_code /
        (8.0 * core::constants::k_pi * g_code) * config.cosmology.omega_baryon;
  }

  void execute(SourceMutationStageView& view) override {
    view.requireFresh();
    core::StepContext& context = internal::RuntimeStageAccess::sourceContext(
        view,
        {{RuntimeResourceKey::kSourceMutationState, RuntimeResourceAccessMode::kReadWrite},
         {RuntimeResourceKey::kParticlePosition, RuntimeResourceAccessMode::kReadWrite},
         {RuntimeResourceKey::kParticleVelocity, RuntimeResourceAccessMode::kReadWrite},
         {RuntimeResourceKey::kParticleIdentity, RuntimeResourceAccessMode::kReadWrite},
         {RuntimeResourceKey::kParticleSpeciesIndex, RuntimeResourceAccessMode::kReadWrite},
         {RuntimeResourceKey::kHydroConservedState, RuntimeResourceAccessMode::kReadWrite},
         {RuntimeResourceKey::kHydroPrimitiveState, RuntimeResourceAccessMode::kRead},
         {RuntimeResourceKey::kAmrPatchState, RuntimeResourceAccessMode::kRead},
         {RuntimeResourceKey::kEffectiveIsmThermodynamics, RuntimeResourceAccessMode::kRead},
         {RuntimeResourceKey::kMigrationOwnership, RuntimeResourceAccessMode::kReadWrite},
         {RuntimeResourceKey::kIntegratorTruth, RuntimeResourceAccessMode::kRead}});
    if (context.stage != core::IntegrationStage::kSourceTerms) {
      throw std::logic_error("source runtime received a non-source stage");
    }

    // SourceTerms executes after drift/hydro on the step-end state while the
    // committed IntegratorState remains at the step beginning until commitStep().
    // Make the physical evaluation epoch explicit so every source model sees
    // the state and cosmological conversion from the same instant.
    const double source_evaluation_scale_factor = m_is_cosmological
        ? context.timeline_step.scale_factor_end
        : 1.0;
    if (!std::isfinite(source_evaluation_scale_factor) ||
        source_evaluation_scale_factor <= 0.0) {
      throw std::runtime_error(
          "source runtime requires a finite positive source-stage scale factor");
    }

    const std::size_t cell_count = context.state.cells.size();
    if (cell_count > 0U && m_star_formation.config().enabled) {
      std::span<const std::uint32_t> active_cells = context.active_set.cell_indices;
      if (!context.active_set.cells_are_subset && active_cells.empty()) {
        m_full_cell_indices.resize(cell_count);
        std::iota(m_full_cell_indices.begin(), m_full_cell_indices.end(), 0U);
        active_cells = m_full_cell_indices;
      }
      buildStarFormationInputs(context, active_cells, source_evaluation_scale_factor);
      const std::size_t particle_count_before_birth = context.state.particles.size();
      const physics::StarFormationStepReport report = m_star_formation.applyFromInputs(
          context.state,
          m_star_formation_inputs,
          context.integrator_state.dt_time_code,
          source_evaluation_scale_factor,
          context.integrator_state.step_index,
          &m_particle_id_registry);
      if (report.counters.spawned_particles > 0U) {
        const std::size_t particle_count_after_birth = context.state.particles.size();
        if (particle_count_after_birth < particle_count_before_birth ||
            particle_count_after_birth - particle_count_before_birth !=
                static_cast<std::size_t>(report.counters.spawned_particles)) {
          throw std::runtime_error(
              "source runtime star-formation report disagrees with appended particle rows");
        }
        // Gas cells are the hydro and gas-gravity mass authority; generic gas
        // particles are compatibility mirrors only. Refresh legacy mirrors once
        // after the batch so I/O/adapter consumers observe the post-birth mass.
        internal::synchronizeParentParticleCompatibilityMirrors(
            context.state, m_world_rank, "SourceRuntime star-formation batch");
        if (context.newly_created_particle_ids != nullptr) {
          context.newly_created_particle_ids->insert(
              context.newly_created_particle_ids->end(),
              context.state.particle_sidecar.particle_id.begin() +
                  static_cast<std::ptrdiff_t>(particle_count_before_birth),
              context.state.particle_sidecar.particle_id.end());
        }
      }
      if (report.counters.spawned_particles > 0U && context.particle_scheduler != nullptr) {
        const std::uint64_t current_tick = context.particle_scheduler->currentTick();
        if (current_tick == std::numeric_limits<std::uint64_t>::max()) {
          throw std::overflow_error("source runtime newborn activation tick overflows uint64");
        }
        // Newborn stars have no retroactive kick. They join bin zero at the next
        // legal scheduler tick and become visible to the next force boundary.
        context.particle_scheduler->appendElements(
            static_cast<std::uint32_t>(report.counters.spawned_particles),
            0U,
            current_tick + 1U);
      }
    }

    executeStellarEvolutionAndEnrichment(context);
    executeMetalDiffusion(context, source_evaluation_scale_factor);

    const physics::BlackHoleAgnStepReport black_hole_report = m_black_hole.apply(
        context.state,
        m_seed_candidates,
        context.integrator_state.dt_time_code,
        source_evaluation_scale_factor,
        m_coordinate_frame == core::CoordinateFrame::kComoving,
        context.integrator_state.step_index,
        &m_particle_id_registry);
    if (black_hole_report.counters.gas_mass_removed_code > 0.0) {
      // Generic gas particles are compatibility mirrors only. Keep them coherent
      // for I/O/legacy consumers after the authoritative gas->BH transaction.
      internal::synchronizeParentParticleCompatibilityMirrors(
          context.state, m_world_rank, "SourceRuntime BH accretion/seeding batch");
    }
    if (!black_hole_report.seeded_particle_ids.empty()) {
      if (context.newly_created_particle_ids != nullptr) {
        context.newly_created_particle_ids->insert(
            context.newly_created_particle_ids->end(),
            black_hole_report.seeded_particle_ids.begin(),
            black_hole_report.seeded_particle_ids.end());
      }
      if (context.particle_scheduler != nullptr) {
        const std::uint64_t current_tick = context.particle_scheduler->currentTick();
        if (current_tick == std::numeric_limits<std::uint64_t>::max()) {
          throw std::overflow_error("source runtime BH activation tick overflows uint64");
        }
        context.particle_scheduler->appendElements(
            static_cast<std::uint32_t>(black_hole_report.seeded_particle_ids.size()),
            0U,
            current_tick + 1U);
      }
    }
  }

 private:
  void buildOwnedLeafCellMetadata(const core::StepContext& context) {
    const core::SimulationState& state = context.state;
    const std::size_t cell_count = state.cells.size();
    m_cell_volume_code.assign(cell_count, 0.0);
    m_owned_leaf_mask.assign(cell_count, 0U);
    m_feedback_candidate_cells.clear();
    m_feedback_candidate_cells.reserve(cell_count);
    std::unordered_set<std::uint64_t> patch_ids_with_children;
    patch_ids_with_children.reserve(state.patches.size());
    for (std::size_t patch_index = 0; patch_index < state.patches.size(); ++patch_index) {
      const std::uint64_t parent_id = state.patches.parent_patch_id[patch_index];
      if (parent_id != 0U) {
        patch_ids_with_children.insert(parent_id);
      }
    }
    for (std::uint32_t cell_index = 0; cell_index < cell_count; ++cell_index) {
      const PatchCellGeometry geometry = starFormationPatchCellGeometry(state, cell_index);
      bool owned_leaf = true;
      double volume = state.gas_cells.density_code[cell_index] > 0.0
          ? state.cells.mass_code[cell_index] /
                state.gas_cells.density_code[cell_index]
          : 0.0;
      if (geometry.valid) {
        volume = geometry.dx_stored * geometry.dy_stored * geometry.dz_stored;
        owned_leaf = state.patches.owning_rank[geometry.patch_index] == m_world_rank &&
            starFormationPatchIsLeaf(
                state, geometry.patch_index, patch_ids_with_children);
      }
      m_cell_volume_code[cell_index] = volume;
      m_owned_leaf_mask[cell_index] = owned_leaf ? 1U : 0U;
      if (owned_leaf) {
        m_feedback_candidate_cells.push_back(cell_index);
      }
    }
  }

  void buildActiveStarRows(const core::StepContext& context) {
    const core::SimulationState& state = context.state;
    m_active_star_indices.clear();
    if (!context.active_set.particles_are_subset ||
        context.active_set.particle_indices.empty()) {
      m_active_star_indices.resize(state.star_particles.size());
      std::iota(m_active_star_indices.begin(), m_active_star_indices.end(), 0U);
      return;
    }
    std::unordered_set<std::uint32_t> active_particles(
        context.active_set.particle_indices.begin(),
        context.active_set.particle_indices.end());
    for (std::uint32_t star_index = 0;
         star_index < state.star_particles.size(); ++star_index) {
      if (active_particles.contains(
              state.star_particles.particle_index[star_index])) {
        m_active_star_indices.push_back(star_index);
      }
    }
  }

  [[nodiscard]] double elapsedStellarEvolutionYears(
      const core::StepContext& context) const {
    constexpr double k_seconds_per_year = 31557600.0;
    double elapsed_si = context.timeline_step.dt_time_si;
    if (m_is_cosmological && context.cosmology_background != nullptr &&
        context.timeline_step.scale_factor_begin > 0.0 &&
        context.timeline_step.scale_factor_end >=
            context.timeline_step.scale_factor_begin) {
      elapsed_si = context.cosmology_background->cosmicTimeIntervalSi(
          context.timeline_step.scale_factor_begin,
          context.timeline_step.scale_factor_end);
    } else if (!(elapsed_si > 0.0)) {
      elapsed_si = context.integrator_state.dt_time_code * m_units.timeSiPerCode();
    }
    return std::max(elapsed_si / k_seconds_per_year, 0.0);
  }

  void executeStellarEvolutionAndEnrichment(core::StepContext& context) {
    if (!m_stellar_evolution.config().enabled ||
        context.state.star_particles.size() == 0U) {
      return;
    }
    buildActiveStarRows(context);
    if (m_active_star_indices.empty()) {
      return;
    }
    buildOwnedLeafCellMetadata(context);
    const double elapsed_years = elapsedStellarEvolutionYears(context);
    const physics::StellarEvolutionStepReport evolution_report =
        m_stellar_evolution.evaluateElapsedYears(
            context.state, m_active_star_indices, elapsed_years);

    const std::size_t star_count = context.state.star_particles.size();
    m_returned_mass_delta_code.assign(star_count, 0.0);
    m_returned_metals_delta_code.assign(star_count, 0.0);
    m_feedback_energy_delta_erg.assign(star_count, 0.0);
    for (const physics::StellarEvolutionStarBudget& budget :
         evolution_report.budgets) {
      m_returned_mass_delta_code[budget.star_index] =
          budget.interval.returned_mass_code;
      m_returned_metals_delta_code[budget.star_index] =
          budget.interval.returned_metals_code;
      m_feedback_energy_delta_erg[budget.star_index] =
          budget.interval.feedback_energy_erg;
    }

    const physics::StellarFeedbackGeometryView geometry_view{
        .particle_position_x_comoving =
            context.state.particles.position_x_comoving,
        .particle_position_y_comoving =
            context.state.particles.position_y_comoving,
        .particle_position_z_comoving =
            context.state.particles.position_z_comoving,
        .cell_center_x_comoving = context.state.cells.center_x_comoving,
        .cell_center_y_comoving = context.state.cells.center_y_comoving,
        .cell_center_z_comoving = context.state.cells.center_z_comoving,
        .gas_cell_id = context.state.gas_cells.gas_cell_id,
        .is_owned_leaf = m_owned_leaf_mask,
        .candidate_cell_indices = m_feedback_candidate_cells,
    };
    physics::StellarFeedbackDepositionView deposition_view{
        .cell_mass_code = context.state.cells.mass_code,
        .gas_density_code = context.state.gas_cells.density_code,
        .gas_internal_energy_code =
            context.state.gas_cells.internal_energy_code,
        .gas_metal_mass_code = context.state.gas_cells.metal_mass_code,
        .cell_volume_code = m_cell_volume_code,
    };
    (void)m_stellar_feedback.applyWithViews(
        context.state, m_stellar_feedback_state, geometry_view,
        deposition_view, m_active_star_indices, m_returned_mass_delta_code,
        m_returned_metals_delta_code, context.integrator_state.dt_time_code,
        m_feedback_energy_delta_erg);

    // The returned budget is now either deposited or durably attached to its
    // source star. Only after that transaction succeeds may stellar mass and
    // cumulative SSP bookkeeping advance.
    m_stellar_evolution.commitBudgets(context.state, evolution_report);
    internal::synchronizeParentParticleCompatibilityMirrors(
        context.state, m_world_rank,
        "SourceRuntime stellar-evolution enrichment batch");
  }

  void executeMetalDiffusion(
      core::StepContext& context,
      double source_evaluation_scale_factor) {
    if (!m_metal_diffusion.config().enabled || context.state.cells.size() == 0U) {
      return;
    }
    buildOwnedLeafCellMetadata(context);
    const double scale_factor = m_coordinate_frame == core::CoordinateFrame::kComoving
        ? source_evaluation_scale_factor
        : 1.0;
    m_diffusion_cells.clear();
    m_diffusion_cells.resize(context.state.cells.size());
    for (std::uint32_t cell_index = 0;
         cell_index < context.state.cells.size(); ++cell_index) {
      const PatchCellGeometry geometry =
          starFormationPatchCellGeometry(context.state, cell_index);
      const double volume_stored = m_cell_volume_code[cell_index];
      const double volume_phys = volume_stored * scale_factor * scale_factor *
          scale_factor;
      physics::MetalDiffusionCell cell;
      cell.gas_cell_id = context.state.gas_cells.gas_cell_id[cell_index];
      cell.gas_mass_code = context.state.cells.mass_code[cell_index];
      cell.metal_mass_code =
          context.state.gas_cells.metal_mass_code[cell_index];
      cell.volume_code = volume_phys;
      cell.density_code = volume_phys > 0.0
          ? cell.gas_mass_code / volume_phys : 0.0;
      cell.filter_length_code = volume_phys > 0.0
          ? std::cbrt(volume_phys) : 0.0;
      cell.is_owned_leaf = m_owned_leaf_mask[cell_index] != 0U;
      if (geometry.valid && cell.is_owned_leaf) {
        const std::array<std::span<const double>, 3> velocity_fields{
            context.state.gas_cells.velocity_x_peculiar,
            context.state.gas_cells.velocity_y_peculiar,
            context.state.gas_cells.velocity_z_peculiar};
        const std::array<double, 3> spacing{
            geometry.dx_stored * scale_factor,
            geometry.dy_stored * scale_factor,
            geometry.dz_stored * scale_factor};
        for (int component = 0; component < 3; ++component) {
          for (int axis = 0; axis < 3; ++axis) {
            cell.velocity_gradient.grad[component][axis] =
                starFormationDerivativeAtCell(
                    velocity_fields[component], geometry, axis,
                    spacing[axis], cell_index);
          }
        }
      }
      m_diffusion_cells[cell_index] = cell;
    }

    m_diffusion_faces.clear();
    for (std::uint32_t patch_index = 0;
         patch_index < context.state.patches.size(); ++patch_index) {
      if (context.state.patches.owning_rank[patch_index] != m_world_rank) {
        continue;
      }
      const std::uint32_t first = context.state.patches.first_cell[patch_index];
      const std::uint32_t nx = context.state.patches.cell_dim_x[patch_index];
      const std::uint32_t ny = context.state.patches.cell_dim_y[patch_index];
      const std::uint32_t nz = context.state.patches.cell_dim_z[patch_index];
      if (nx == 0U || ny == 0U || nz == 0U) {
        continue;
      }
      const double dx = context.state.patches.extent_x_comoving[patch_index] /
          static_cast<double>(nx) * scale_factor;
      const double dy = context.state.patches.extent_y_comoving[patch_index] /
          static_cast<double>(ny) * scale_factor;
      const double dz = context.state.patches.extent_z_comoving[patch_index] /
          static_cast<double>(nz) * scale_factor;
      const auto row = [first, nx, ny](std::uint32_t i, std::uint32_t j,
                                      std::uint32_t k) {
        return first + i + nx * (j + ny * k);
      };
      for (std::uint32_t k = 0; k < nz; ++k) {
        for (std::uint32_t j = 0; j < ny; ++j) {
          for (std::uint32_t i = 0; i < nx; ++i) {
            const std::uint32_t left = row(i, j, k);
            if (i + 1U < nx) {
              m_diffusion_faces.push_back({
                  .left_cell = left, .right_cell = row(i + 1U, j, k),
                  .area_code = dy * dz, .center_distance_code = dx});
            }
            if (j + 1U < ny) {
              m_diffusion_faces.push_back({
                  .left_cell = left, .right_cell = row(i, j + 1U, k),
                  .area_code = dx * dz, .center_distance_code = dy});
            }
            if (k + 1U < nz) {
              m_diffusion_faces.push_back({
                  .left_cell = left, .right_cell = row(i, j, k + 1U),
                  .area_code = dx * dy, .center_distance_code = dz});
            }
          }
        }
      }
    }
    if (m_diffusion_faces.empty()) {
      return;
    }
    const physics::MetalDiffusionStepReport report = m_metal_diffusion.advance(
        m_diffusion_cells, m_diffusion_faces,
        context.integrator_state.dt_time_code);
    for (std::uint32_t cell_index = 0;
         cell_index < context.state.cells.size(); ++cell_index) {
      if (m_diffusion_cells[cell_index].is_owned_leaf) {
        context.state.gas_cells.metal_mass_code[cell_index] =
            m_diffusion_cells[cell_index].metal_mass_code;
      }
    }
    std::ostringstream metadata;
    metadata << "module=metal_diffusion\n";
    metadata << "model=" << core::metalDiffusionModelToString(
        m_metal_diffusion.config().model) << "\n";
    metadata << "time_integrator=" <<
        core::metalDiffusionTimeIntegratorToString(
            m_metal_diffusion.config().time_integrator) << "\n";
    metadata << "stable_dt_code=" << report.stable_dt_code << "\n";
    metadata << "subcycles=" << report.subcycles << "\n";
    metadata << "rkl_stages=" << report.rkl_stages << "\n";
    metadata << "conservation_residual_code=" <<
        report.conservation_residual_code << "\n";
    const std::string text = metadata.str();
    core::ModuleSidecarBlock sidecar;
    sidecar.module_name = "metal_diffusion";
    sidecar.schema_version = 1U;
    sidecar.payload.resize(text.size());
    for (std::size_t i = 0; i < text.size(); ++i) {
      sidecar.payload[i] = static_cast<std::byte>(text[i]);
    }
    context.state.sidecars.upsert(std::move(sidecar));
  }

  void buildStarFormationInputs(
      const core::StepContext& context,
      std::span<const std::uint32_t> active_cells,
      double source_evaluation_scale_factor) {
    const core::SimulationState& state = context.state;
    m_star_formation_inputs.clear();
    m_star_formation_inputs.reserve(active_cells.size());

    std::unordered_set<std::uint64_t> patch_ids_with_children;
    patch_ids_with_children.reserve(state.patches.size());
    for (std::size_t patch_index = 0; patch_index < state.patches.size(); ++patch_index) {
      const std::uint64_t parent_id = state.patches.parent_patch_id[patch_index];
      if (parent_id != 0U) {
        patch_ids_with_children.insert(parent_id);
      }
    }

    const double scale_factor = source_evaluation_scale_factor;
    const double length_to_physical =
        m_coordinate_frame == core::CoordinateFrame::kComoving ? scale_factor : 1.0;

    for (const std::uint32_t cell_index : active_cells) {
      physics::StarFormationCellInput input;
      input.cell_index = cell_index;
      input.owning_rank = m_world_rank;
      input.is_active = true;
      if (cell_index >= state.cells.size() || cell_index >= state.gas_cells.size()) {
        input.is_active = false;
        m_star_formation_inputs.push_back(input);
        continue;
      }

      input.gas_cell_id = state.gas_cells.gas_cell_id[cell_index];
      if (input.gas_cell_id == 0U) {
        input.gas_cell_id = state.gasCellIdForLocalRow(cell_index).value_or(0U);
      }
      input.gas_mass_code = state.cells.mass_code[cell_index];
      input.gas_density_code = state.gas_cells.density_code[cell_index];
      input.gas_temperature_k = state.gas_cells.temperature_code[cell_index];
      input.gas_sound_speed_code = state.gas_cells.sound_speed_code[cell_index];
      input.gas_specific_internal_energy_code =
          state.gas_cells.internal_energy_code[cell_index];
      input.is_cosmological = m_is_cosmological;
      if (m_is_cosmological && m_mean_baryon_density0_code > 0.0 && scale_factor > 0.0) {
        const double density_phys = m_coordinate_frame == core::CoordinateFrame::kComoving
            ? input.gas_density_code / (scale_factor * scale_factor * scale_factor)
            : input.gas_density_code;
        const double mean_density_phys = m_mean_baryon_density0_code /
            (scale_factor * scale_factor * scale_factor);
        input.baryon_overdensity = density_phys / std::max(mean_density_phys, 1.0e-30);
      }
      input.velocity_x_peculiar = state.gas_cells.velocity_x_peculiar[cell_index];
      input.velocity_y_peculiar = state.gas_cells.velocity_y_peculiar[cell_index];
      input.velocity_z_peculiar = state.gas_cells.velocity_z_peculiar[cell_index];
      input.gas_metal_mass_code = state.gas_cells.metal_mass_code[cell_index];
      input.metallicity_mass_fraction = input.gas_mass_code > 0.0
          ? input.gas_metal_mass_code / input.gas_mass_code
          : 0.0;
      input.center_x_comoving = state.cells.center_x_comoving[cell_index];
      input.center_y_comoving = state.cells.center_y_comoving[cell_index];
      input.center_z_comoving = state.cells.center_z_comoving[cell_index];

      const PatchCellGeometry geometry = starFormationPatchCellGeometry(state, cell_index);
      if (geometry.valid) {
        input.cell_volume_code = geometry.dx_stored * geometry.dy_stored * geometry.dz_stored;
        input.is_owned = state.patches.owning_rank[geometry.patch_index] == m_world_rank;
        input.is_leaf = starFormationPatchIsLeaf(state, geometry.patch_index, patch_ids_with_children);
        input.is_ghost = !input.is_owned;

        const double dx_phys = geometry.dx_stored * length_to_physical;
        const double dy_phys = geometry.dy_stored * length_to_physical;
        const double dz_phys = geometry.dz_stored * length_to_physical;
        const std::array<std::span<const double>, 3> velocity_fields{
            state.gas_cells.velocity_x_peculiar,
            state.gas_cells.velocity_y_peculiar,
            state.gas_cells.velocity_z_peculiar};
        std::array<std::array<double, 3>, 3> gradient{};
        bool gradient_valid = true;
        for (int component = 0; component < 3; ++component) {
          gradient[component][0] = starFormationDerivativeAtCell(
              velocity_fields[component], geometry, 0, dx_phys, cell_index);
          gradient[component][1] = starFormationDerivativeAtCell(
              velocity_fields[component], geometry, 1, dy_phys, cell_index);
          gradient[component][2] = starFormationDerivativeAtCell(
              velocity_fields[component], geometry, 2, dz_phys, cell_index);
          for (int axis = 0; axis < 3; ++axis) {
            gradient_valid = gradient_valid && std::isfinite(gradient[component][axis]);
          }
        }
        if (gradient_valid) {
          input.velocity_divergence_code =
              gradient[0][0] + gradient[1][1] + gradient[2][2];
          input.velocity_gradient_frobenius_sq_code = 0.0;
          for (const auto& row : gradient) {
            for (const double value : row) {
              input.velocity_gradient_frobenius_sq_code += value * value;
            }
          }
        } else {
          input.velocity_divergence_code = std::numeric_limits<double>::quiet_NaN();
          input.velocity_gradient_frobenius_sq_code =
              std::numeric_limits<double>::quiet_NaN();
        }
      } else {
        input.cell_volume_code = input.gas_density_code > 0.0
            ? input.gas_mass_code / input.gas_density_code
            : 0.0;
        input.is_owned = true;
        input.is_leaf = true;
        input.is_ghost = false;
        // Adaptive production requires a real patch gradient. Legacy mode keeps the
        // historical converging-flow value of zero for non-AMR compatibility states.
        if (m_star_formation.config().model ==
            core::StarFormationModelKind::kAdaptiveBoundJeans) {
          input.velocity_divergence_code = std::numeric_limits<double>::quiet_NaN();
          input.velocity_gradient_frobenius_sq_code =
              std::numeric_limits<double>::quiet_NaN();
        } else {
          input.velocity_divergence_code = 0.0;
          input.velocity_gradient_frobenius_sq_code = 0.0;
        }
      }
      m_star_formation_inputs.push_back(input);
    }
  }

  core::UnitSystem m_units;
  std::shared_ptr<const physics::EffectiveMultiphaseEosTable> m_effective_eos_table;
  physics::StarFormationModel m_star_formation;
  physics::StellarEvolutionBookkeeper m_stellar_evolution;
  physics::StellarFeedbackModel m_stellar_feedback;
  physics::MetalDiffusionModel m_metal_diffusion;
  physics::StellarFeedbackModuleState m_stellar_feedback_state;
  physics::BlackHoleAgnModel m_black_hole;
  DistributedParticleIdRegistry m_particle_id_registry;
  std::uint32_t m_world_rank = 0;
  core::CoordinateFrame m_coordinate_frame = core::CoordinateFrame::kComoving;
  bool m_is_cosmological = false;
  double m_mean_baryon_density0_code = 0.0;
  std::vector<std::uint32_t> m_full_cell_indices;
  std::vector<physics::StarFormationCellInput> m_star_formation_inputs;
  std::vector<std::uint32_t> m_active_star_indices;
  std::vector<double> m_returned_mass_delta_code;
  std::vector<double> m_returned_metals_delta_code;
  std::vector<double> m_feedback_energy_delta_erg;
  std::vector<double> m_cell_volume_code;
  std::vector<std::uint8_t> m_owned_leaf_mask;
  std::vector<std::uint32_t> m_feedback_candidate_cells;
  std::vector<physics::MetalDiffusionCell> m_diffusion_cells;
  std::vector<physics::MetalDiffusionFace> m_diffusion_faces;
  std::vector<physics::BlackHoleSeedCandidate> m_seed_candidates;
};

}  // namespace

std::unique_ptr<SourceRuntime> makeSourceRuntime(
    const core::SimulationConfig& config,
    const core::UnitSystem& units,
    std::uint32_t world_rank,
    const parallel::MpiContext& mpi_context) {
  return std::make_unique<SourceRuntimeImpl>(config, units, world_rank, mpi_context);
}

}  // namespace cosmosim::workflows
