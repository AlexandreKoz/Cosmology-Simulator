#include "cosmosim/workflows/source_runtime.hpp"

#include <cstdint>
#include <memory>
#include <numeric>
#include <span>
#include <stdexcept>
#include <vector>

#include "cosmosim/physics/black_hole_agn.hpp"
#include "cosmosim/physics/star_formation.hpp"

namespace cosmosim::workflows {
namespace {

[[nodiscard]] double newtonGCodeFromUnits(const core::UnitSystem& units) {
  return core::newtonGravitationalConstantCode(units);
}

[[nodiscard]] physics::StarFormationConfig makeRuntimeStarFormationConfig(
    const core::PhysicsConfig& physics_config,
    const core::UnitSystem& units) {
  physics::StarFormationConfig config =
      physics::makeStarFormationConfig(physics_config);
  config.newton_g_code = newtonGCodeFromUnits(units);
  return config;
}

[[nodiscard]] physics::BlackHoleAgnConfig makeRuntimeBlackHoleAgnConfig(
    const core::PhysicsConfig& physics_config,
    const core::UnitSystem& units) {
  physics::BlackHoleAgnConfig config =
      physics::makeBlackHoleAgnConfig(physics_config);
  config.proton_mass_code = config.proton_mass_si / units.mass_si_per_code;
  config.thomson_cross_section_code = config.thomson_cross_section_si /
      (units.length_si_per_code * units.length_si_per_code);
  config.newton_g_code = config.newton_g_si * units.mass_si_per_code *
      units.timeSiPerCode() * units.timeSiPerCode() /
      (units.length_si_per_code * units.length_si_per_code *
       units.length_si_per_code);
  config.speed_of_light_code =
      config.speed_of_light_si / units.velocity_si_per_code;
  return config;
}

class SourceRuntimeImpl final : public SourceRuntime {
 public:
  SourceRuntimeImpl(
      const core::SimulationConfig& config,
      const core::UnitSystem& units,
      std::uint32_t world_rank)
      : m_star_formation(
            makeRuntimeStarFormationConfig(config.physics, units)),
        m_black_hole(makeRuntimeBlackHoleAgnConfig(config.physics, units)),
        m_rank_local_seed_offset(world_rank) {}

  void execute(SourceMutationStageView& view) override {
    view.requireFresh();
    core::StepContext& context = stageContext(view);
    if (context.stage != core::IntegrationStage::kSourceTerms) {
      throw std::logic_error("source runtime received a non-source stage");
    }
    const std::size_t cell_count = context.state.cells.size();
    if (cell_count > 0U) {
      if (m_velocity_divergence_code.size() < cell_count) {
        m_velocity_divergence_code.resize(cell_count, 0.0);
      }
      if (m_metallicity_mass_fraction.size() < cell_count) {
        m_metallicity_mass_fraction.resize(cell_count, 0.0);
      }
      std::span<const std::uint32_t> active_cells =
          context.active_set.cell_indices;
      if (!context.active_set.cells_are_subset && active_cells.empty()) {
        m_full_cell_indices.resize(cell_count);
        std::iota(m_full_cell_indices.begin(), m_full_cell_indices.end(), 0U);
        active_cells = m_full_cell_indices;
      }
      (void)m_star_formation.apply(
          context.state,
          active_cells,
          m_velocity_divergence_code,
          m_metallicity_mass_fraction,
          context.integrator_state.dt_time_code,
          context.integrator_state.current_scale_factor,
          context.integrator_state.step_index,
          m_rank_local_seed_offset);
    }
    (void)m_black_hole.apply(
        context.state,
        m_seed_candidates,
        context.integrator_state.dt_time_code,
        context.integrator_state.step_index);
  }

 private:
  physics::StarFormationModel m_star_formation;
  physics::BlackHoleAgnModel m_black_hole;
  std::uint32_t m_rank_local_seed_offset = 0;
  std::vector<std::uint32_t> m_full_cell_indices;
  std::vector<double> m_velocity_divergence_code;
  std::vector<double> m_metallicity_mass_fraction;
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
