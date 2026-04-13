#include <array>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include "cosmosim/amr/amr_framework.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/physics/star_formation.hpp"
#include "validation_tolerance.hpp"

namespace {

void requireOrThrow(bool condition, const std::string& message) {
  if (!condition) {
    throw std::runtime_error(message);
  }
}

void seedRefinementHotspot(cosmosim::amr::AmrPatch& patch) {
  auto metrics = patch.metricsView();
  for (std::size_t i = 0; i < metrics.size(); ++i) {
    metrics[i].density_code = (i == 0) ? 16.0 : 0.05;
    metrics[i].sound_speed_code = 1.0;
    metrics[i].gradient_indicator = (i == 0) ? 2.0 : 0.1;
    metrics[i].particle_count = (i == 0) ? 16 : 0;
  }
}

void testAmrRefluxRegression(const cosmosim::validation::ValidationToleranceTable& tolerances) {
  cosmosim::amr::PatchHierarchy hierarchy;
  cosmosim::amr::PatchDescriptor root;
  root.cell_dims = {4, 4, 4};
  const std::uint64_t root_id = hierarchy.createRootPatch(root);

  auto* root_patch = hierarchy.findPatch(root_id);
  requireOrThrow(root_patch != nullptr, "amr_reflux: root patch missing");

  auto root_conserved = root_patch->conservedView();
  for (auto& cell : root_conserved) {
    cell.mass_code = 1.0;
    cell.total_energy_code = 2.0;
  }

  seedRefinementHotspot(*root_patch);

  cosmosim::amr::RefinementCriteria criteria;
  criteria.mass_threshold_code = 0.5;
  criteria.gradient_threshold = 1.0;
  criteria.particle_threshold = 8;
  criteria.jeans_resolution_cells = 0.5;

  cosmosim::amr::RefineDerefineManager manager(criteria);
  const auto regrid_diag = manager.regrid(hierarchy);
  requireOrThrow(regrid_diag.refined_patch_count == 1, "amr_reflux: expected one refined patch");

  std::array<cosmosim::amr::FluxRegisterEntry, 1> entries{{
      {.coarse_patch_id = root_id,
       .coarse_cell_index = 0,
       .coarse_face_flux_code = {.mass_code = 2.0, .total_energy_code = 4.0},
       .fine_face_flux_code = {.mass_code = 2.5, .total_energy_code = 4.8},
       .face_area_comov = 0.25,
       .dt_code = 0.1},
  }};

  const auto reflux_diag = cosmosim::amr::RefluxSynchronizer::apply(hierarchy, entries);
  requireOrThrow(reflux_diag.corrected_cells == 1, "amr_reflux: expected one corrected cell");

  root_patch = hierarchy.findPatch(root_id);
  requireOrThrow(root_patch != nullptr, "amr_reflux: root patch unavailable after reflux");

  const double expected_mass = 1.0 - ((2.5 - 2.0) * 0.25 * 0.1 / root_patch->cellVolumeComov());
  const double mass_error = std::abs(root_patch->conservedView()[0].mass_code - expected_mass);
  requireOrThrow(
      mass_error <= tolerances.require("amr_reflux.mass_correction_abs"),
      "amr_reflux: coarse-cell mass correction mismatch");
}

void testStarFormationMassBudgetRegression(const cosmosim::validation::ValidationToleranceTable& tolerances) {
  cosmosim::core::SimulationState state;
  state.resizeCells(8);
  for (std::size_t i = 0; i < state.cells.size(); ++i) {
    state.cells.mass_code[i] = 2.0;
    state.gas_cells.density_code[i] = 20.0;
    state.gas_cells.temperature_code[i] = 8.0e3;
  }

  cosmosim::physics::StarFormationConfig config;
  config.epsilon_ff = 0.02;
  config.stochastic_spawning = false;
  cosmosim::physics::StarFormationModel model(config);

  std::vector<std::uint32_t> active(state.cells.size());
  for (std::size_t i = 0; i < active.size(); ++i) {
    active[i] = static_cast<std::uint32_t>(i);
  }
  std::vector<double> div_v(state.cells.size(), -0.5);
  std::vector<double> metallicity(state.cells.size(), 0.01);

  double gas_before = 0.0;
  for (double mass : state.cells.mass_code) {
    gas_before += mass;
  }

  const auto report = model.apply(state, active, div_v, metallicity, 3.0e8, 0.9, 4242, 0);
  requireOrThrow(report.counters.spawned_particles == state.cells.size(), "star_formation: all cells should spawn");

  double gas_after = 0.0;
  for (double mass : state.cells.mass_code) {
    gas_after += mass;
  }

  double star_mass = 0.0;
  for (double mass : state.particles.mass_code) {
    star_mass += mass;
  }

  const double mass_budget_error = std::abs((gas_after + star_mass) - gas_before);
  requireOrThrow(
      mass_budget_error <= tolerances.require("star_formation.mass_budget_abs"),
      "star_formation: gas+stars mass budget mismatch");
}

}  // namespace

int main() {
  const auto tolerances = cosmosim::validation::ValidationToleranceTable::loadFromFile(
      std::string(COSMOSIM_SOURCE_DIR) + "/validation/reference/validation_tolerances_v1.txt");

  testAmrRefluxRegression(tolerances);
  testStarFormationMassBudgetRegression(tolerances);
  return 0;
}
