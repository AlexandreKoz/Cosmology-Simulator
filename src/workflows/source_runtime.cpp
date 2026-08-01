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
#include "cosmosim/physics/star_formation.hpp"
#include "workflows/internal/gas_cell_ownership.hpp"
#include "workflows/internal/runtime_stage_resource_access.hpp"
#include "workflows/internal/star_formation_geometry.hpp"

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
      std::uint32_t world_rank)
      : m_star_formation(makeRuntimeStarFormationConfig(config, units)),
        m_black_hole(makeRuntimeBlackHoleAgnConfig(config.physics, units)),
        m_world_rank(world_rank),
        m_coordinate_frame(config.units.coordinate_frame) {}

  void execute(SourceMutationStageView& view) override {
    view.requireFresh();
    core::StepContext& context = internal::RuntimeStageAccess::sourceContext(
        view,
        {{RuntimeResourceKey::kSourceMutationState, RuntimeResourceAccessMode::kReadWrite},
         {RuntimeResourceKey::kParticlePosition, RuntimeResourceAccessMode::kReadWrite},
         {RuntimeResourceKey::kParticleVelocity, RuntimeResourceAccessMode::kReadWrite},
         {RuntimeResourceKey::kHydroConservedState, RuntimeResourceAccessMode::kReadWrite},
         {RuntimeResourceKey::kMigrationOwnership, RuntimeResourceAccessMode::kReadWrite},
         {RuntimeResourceKey::kIntegratorTruth, RuntimeResourceAccessMode::kRead}});
    if (context.stage != core::IntegrationStage::kSourceTerms) {
      throw std::logic_error("source runtime received a non-source stage");
    }

    const std::size_t cell_count = context.state.cells.size();
    if (cell_count > 0U && m_star_formation.config().enabled) {
      std::span<const std::uint32_t> active_cells = context.active_set.cell_indices;
      if (!context.active_set.cells_are_subset && active_cells.empty()) {
        m_full_cell_indices.resize(cell_count);
        std::iota(m_full_cell_indices.begin(), m_full_cell_indices.end(), 0U);
        active_cells = m_full_cell_indices;
      }
      buildStarFormationInputs(context, active_cells);
      const std::size_t particle_count_before_birth = context.state.particles.size();
      const physics::StarFormationStepReport report = m_star_formation.applyFromInputs(
          context.state,
          m_star_formation_inputs,
          context.integrator_state.dt_time_code,
          context.integrator_state.current_scale_factor,
          context.integrator_state.step_index);
      if (report.counters.spawned_particles > 0U) {
        const std::size_t particle_count_after_birth = context.state.particles.size();
        if (particle_count_after_birth < particle_count_before_birth ||
            particle_count_after_birth - particle_count_before_birth !=
                static_cast<std::size_t>(report.counters.spawned_particles)) {
          throw std::runtime_error(
              "source runtime star-formation report disagrees with appended particle rows");
        }
        // Gas cells are the hydro mass authority, while generic gas particles are
        // compatibility/gravity mirrors. Refresh the affected aggregate authority
        // once after the batch so the next gravity boundary cannot see both the
        // pre-birth gas mass and the newborn stellar mass.
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

    (void)m_black_hole.apply(
        context.state,
        m_seed_candidates,
        context.integrator_state.dt_time_code,
        context.integrator_state.step_index);
  }

 private:
  void buildStarFormationInputs(
      const core::StepContext& context,
      std::span<const std::uint32_t> active_cells) {
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

    const double scale_factor = context.integrator_state.current_scale_factor;
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

  physics::StarFormationModel m_star_formation;
  physics::BlackHoleAgnModel m_black_hole;
  std::uint32_t m_world_rank = 0;
  core::CoordinateFrame m_coordinate_frame = core::CoordinateFrame::kComoving;
  std::vector<std::uint32_t> m_full_cell_indices;
  std::vector<physics::StarFormationCellInput> m_star_formation_inputs;
  std::vector<physics::BlackHoleSeedCandidate> m_seed_candidates;
};

}  // namespace

std::unique_ptr<SourceRuntime> makeSourceRuntime(
    const core::SimulationConfig& config,
    const core::UnitSystem& units,
    std::uint32_t world_rank) {
  return std::make_unique<SourceRuntimeImpl>(config, units, world_rank);
}

}  // namespace cosmosim::workflows
