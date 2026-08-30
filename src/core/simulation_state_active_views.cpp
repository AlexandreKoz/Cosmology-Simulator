#include "cosmosim/core/simulation_state.hpp"

#include <algorithm>
#include <limits>
#include <unordered_set>

namespace cosmosim::core {

void ActiveIndexSet::clear() {
  particle_indices.clear();
  cell_indices.clear();
}

std::size_t ParticleActiveView::size() const noexcept { return particle_id.size(); }
bool ParticleActiveView::isConsistent() const noexcept {
  const std::size_t expected = particle_id.size();
  return species_tag.size() == expected && position_x_comoving.size() == expected &&
      position_y_comoving.size() == expected && position_z_comoving.size() == expected &&
      velocity_x_peculiar.size() == expected && velocity_y_peculiar.size() == expected &&
      velocity_z_peculiar.size() == expected && mass_code.size() == expected;
}
std::size_t CellActiveView::size() const noexcept { return center_x_comoving.size(); }
bool CellActiveView::isConsistent() const noexcept {
  const std::size_t expected = center_x_comoving.size();
  return center_y_comoving.size() == expected && center_z_comoving.size() == expected &&
      mass_code.size() == expected && patch_index.size() == expected &&
      density_code.size() == expected && pressure_code.size() == expected;
}
std::size_t GravityParticleKernelView::size() const noexcept { return particle_index.size(); }
bool GravityParticleKernelView::isConsistent() const noexcept {
  const std::size_t expected = particle_index.size();
  return position_x_comoving.size() == expected && position_y_comoving.size() == expected &&
      position_z_comoving.size() == expected && velocity_x_peculiar.size() == expected &&
      velocity_y_peculiar.size() == expected && velocity_z_peculiar.size() == expected &&
      mass_code.size() == expected;
}
std::size_t HydroCellKernelView::size() const noexcept { return cell_index.size(); }
bool HydroCellKernelView::isConsistent() const noexcept {
  const std::size_t expected = cell_index.size();
  return gas_cell_id.size() == expected && local_cell_row.size() == expected &&
      center_x_comoving.size() == expected && center_y_comoving.size() == expected &&
      center_z_comoving.size() == expected && mass_code.size() == expected &&
      density_code.size() == expected && pressure_code.size() == expected;
}

void TransientStepWorkspace::clear() {
  particle_id.clear();
  particle_species_tag.clear();
  particle_position_x_comoving.clear();
  particle_position_y_comoving.clear();
  particle_position_z_comoving.clear();
  particle_velocity_x_peculiar.clear();
  particle_velocity_y_peculiar.clear();
  particle_velocity_z_peculiar.clear();
  particle_mass_code.clear();
  gravity_particle_index.clear();

  hydro_cell_index.clear();
  hydro_gas_cell_id.clear();
  hydro_local_cell_row.clear();
  hydro_cell_center_x_comoving.clear();
  hydro_cell_center_y_comoving.clear();
  hydro_cell_center_z_comoving.clear();
  hydro_cell_mass_code.clear();
  hydro_cell_density_code.clear();
  hydro_cell_pressure_code.clear();

  cell_center_x_comoving.clear();
  cell_center_y_comoving.clear();
  cell_center_z_comoving.clear();
  cell_mass_code.clear();
  cell_patch_index.clear();
  cell_density_code.clear();
  cell_pressure_code.clear();

  hydro_recon_gradient_x.clear();
  hydro_recon_gradient_y.clear();
  hydro_recon_gradient_z.clear();

  gravity_particle_index_scratch = {};
  scratch.reset();
}

void TransientStepWorkspace::prepareGravityParticleIndexScratch(
    std::size_t particle_count) {
  if (particle_count == 0U) {
    gravity_particle_index_scratch = {};
    return;
  }
  auto* index_data = scratch.allocateArray<std::uint32_t>(particle_count);
  gravity_particle_index_scratch =
      std::span<std::uint32_t>(index_data, particle_count);
}

ParticleActiveView buildParticleActiveView(
    const SimulationState& state,
    std::span<const std::uint32_t> active_particle_indices,
    TransientStepWorkspace& workspace) {
  workspace.particle_id.resize(active_particle_indices.size());
  workspace.particle_species_tag.resize(active_particle_indices.size());
  workspace.particle_position_x_comoving.resize(active_particle_indices.size());
  workspace.particle_position_y_comoving.resize(active_particle_indices.size());
  workspace.particle_position_z_comoving.resize(active_particle_indices.size());
  workspace.particle_velocity_x_peculiar.resize(active_particle_indices.size());
  workspace.particle_velocity_y_peculiar.resize(active_particle_indices.size());
  workspace.particle_velocity_z_peculiar.resize(active_particle_indices.size());
  workspace.particle_mass_code.resize(active_particle_indices.size());

  for (const std::uint32_t source : active_particle_indices) {
    if (source >= state.particles.size()) {
      throw std::out_of_range("buildParticleActiveView: particle index out of range");
    }
  }

  gatherSpan<std::uint64_t>(state.particle_sidecar.particle_id, active_particle_indices, workspace.particle_id);
  gatherSpan<std::uint32_t>(state.particle_sidecar.species_tag, active_particle_indices, workspace.particle_species_tag);
  gatherSpan<double>(state.particles.position_x_comoving, active_particle_indices, workspace.particle_position_x_comoving);
  gatherSpan<double>(state.particles.position_y_comoving, active_particle_indices, workspace.particle_position_y_comoving);
  gatherSpan<double>(state.particles.position_z_comoving, active_particle_indices, workspace.particle_position_z_comoving);
  gatherSpan<double>(state.particles.velocity_x_peculiar, active_particle_indices, workspace.particle_velocity_x_peculiar);
  gatherSpan<double>(state.particles.velocity_y_peculiar, active_particle_indices, workspace.particle_velocity_y_peculiar);
  gatherSpan<double>(state.particles.velocity_z_peculiar, active_particle_indices, workspace.particle_velocity_z_peculiar);
  gatherSpan<double>(state.particles.mass_code, active_particle_indices, workspace.particle_mass_code);

  return ParticleActiveView{
      .particle_id = workspace.particle_id,
      .species_tag = workspace.particle_species_tag,
      .position_x_comoving = workspace.particle_position_x_comoving,
      .position_y_comoving = workspace.particle_position_y_comoving,
      .position_z_comoving = workspace.particle_position_z_comoving,
      .velocity_x_peculiar = workspace.particle_velocity_x_peculiar,
      .velocity_y_peculiar = workspace.particle_velocity_y_peculiar,
      .velocity_z_peculiar = workspace.particle_velocity_z_peculiar,
      .mass_code = workspace.particle_mass_code,
      .source_particle_index_generation = state.particleIndexGeneration(),
  };
}

CellActiveView buildCellActiveView(
    const SimulationState& state,
    std::span<const std::uint32_t> active_cell_indices,
    TransientStepWorkspace& workspace) {
  workspace.cell_center_x_comoving.resize(active_cell_indices.size());
  workspace.cell_center_y_comoving.resize(active_cell_indices.size());
  workspace.cell_center_z_comoving.resize(active_cell_indices.size());
  workspace.cell_mass_code.resize(active_cell_indices.size());
  workspace.cell_patch_index.resize(active_cell_indices.size());
  workspace.cell_density_code.resize(active_cell_indices.size());
  workspace.cell_pressure_code.resize(active_cell_indices.size());

  for (const std::uint32_t source : active_cell_indices) {
    if (source >= state.cells.size()) {
      throw std::out_of_range("buildCellActiveView: cell index out of range");
    }
  }

  gatherSpan<double>(state.cells.center_x_comoving, active_cell_indices, workspace.cell_center_x_comoving);
  gatherSpan<double>(state.cells.center_y_comoving, active_cell_indices, workspace.cell_center_y_comoving);
  gatherSpan<double>(state.cells.center_z_comoving, active_cell_indices, workspace.cell_center_z_comoving);
  gatherSpan<double>(state.cells.mass_code, active_cell_indices, workspace.cell_mass_code);
  gatherSpan<std::uint32_t>(state.cells.patch_index, active_cell_indices, workspace.cell_patch_index);
  gatherSpan<double>(state.gas_cells.density_code, active_cell_indices, workspace.cell_density_code);
  gatherSpan<double>(state.gas_cells.pressure_code, active_cell_indices, workspace.cell_pressure_code);

  return CellActiveView{
      .center_x_comoving = workspace.cell_center_x_comoving,
      .center_y_comoving = workspace.cell_center_y_comoving,
      .center_z_comoving = workspace.cell_center_z_comoving,
      .mass_code = workspace.cell_mass_code,
      .patch_index = workspace.cell_patch_index,
      .density_code = workspace.cell_density_code,
      .pressure_code = workspace.cell_pressure_code,
      .source_cell_index_generation = state.cellIndexGeneration(),
  };
}

GravityParticleKernelView buildGravityParticleKernelView(
    const SimulationState& state,
    std::span<const std::uint32_t> active_particle_indices,
    TransientStepWorkspace& workspace) {
  workspace.gravity_particle_index.resize(active_particle_indices.size());
  workspace.particle_position_x_comoving.resize(active_particle_indices.size());
  workspace.particle_position_y_comoving.resize(active_particle_indices.size());
  workspace.particle_position_z_comoving.resize(active_particle_indices.size());
  workspace.particle_velocity_x_peculiar.resize(active_particle_indices.size());
  workspace.particle_velocity_y_peculiar.resize(active_particle_indices.size());
  workspace.particle_velocity_z_peculiar.resize(active_particle_indices.size());
  workspace.particle_mass_code.resize(active_particle_indices.size());

  for (std::size_t i = 0; i < active_particle_indices.size(); ++i) {
    const auto source = active_particle_indices[i];
    if (source >= state.particles.size()) {
      throw std::out_of_range("buildGravityParticleKernelView: particle index out of range");
    }
    workspace.gravity_particle_index[i] = source;
  }

  gatherSpan<double>(
      state.particles.position_x_comoving,
      active_particle_indices,
      workspace.particle_position_x_comoving);
  gatherSpan<double>(
      state.particles.position_y_comoving,
      active_particle_indices,
      workspace.particle_position_y_comoving);
  gatherSpan<double>(
      state.particles.position_z_comoving,
      active_particle_indices,
      workspace.particle_position_z_comoving);
  gatherSpan<double>(
      state.particles.velocity_x_peculiar,
      active_particle_indices,
      workspace.particle_velocity_x_peculiar);
  gatherSpan<double>(
      state.particles.velocity_y_peculiar,
      active_particle_indices,
      workspace.particle_velocity_y_peculiar);
  gatherSpan<double>(
      state.particles.velocity_z_peculiar,
      active_particle_indices,
      workspace.particle_velocity_z_peculiar);
  gatherSpan<double>(state.particles.mass_code, active_particle_indices, workspace.particle_mass_code);
  return GravityParticleKernelView{
      .particle_index = workspace.gravity_particle_index,
      .position_x_comoving = workspace.particle_position_x_comoving,
      .position_y_comoving = workspace.particle_position_y_comoving,
      .position_z_comoving = workspace.particle_position_z_comoving,
      .velocity_x_peculiar = workspace.particle_velocity_x_peculiar,
      .velocity_y_peculiar = workspace.particle_velocity_y_peculiar,
      .velocity_z_peculiar = workspace.particle_velocity_z_peculiar,
      .mass_code = workspace.particle_mass_code,
      .source_particle_index_generation = state.particleIndexGeneration(),
  };
}

GravityParticleKernelView buildGravityParticleKernelViewAllParticles(
    const SimulationState& state,
    TransientStepWorkspace& workspace) {
  workspace.gravity_particle_index.resize(state.particles.size());
  workspace.particle_position_x_comoving.resize(state.particles.size());
  workspace.particle_position_y_comoving.resize(state.particles.size());
  workspace.particle_position_z_comoving.resize(state.particles.size());
  workspace.particle_velocity_x_peculiar.resize(state.particles.size());
  workspace.particle_velocity_y_peculiar.resize(state.particles.size());
  workspace.particle_velocity_z_peculiar.resize(state.particles.size());
  workspace.particle_mass_code.resize(state.particles.size());

  for (std::size_t i = 0; i < state.particles.size(); ++i) {
    workspace.gravity_particle_index[i] = checkedLocalParticleRow(i, "gravity all-particle view row");
  }
  std::copy(state.particles.position_x_comoving.begin(), state.particles.position_x_comoving.end(), workspace.particle_position_x_comoving.begin());
  std::copy(state.particles.position_y_comoving.begin(), state.particles.position_y_comoving.end(), workspace.particle_position_y_comoving.begin());
  std::copy(state.particles.position_z_comoving.begin(), state.particles.position_z_comoving.end(), workspace.particle_position_z_comoving.begin());
  std::copy(state.particles.velocity_x_peculiar.begin(), state.particles.velocity_x_peculiar.end(), workspace.particle_velocity_x_peculiar.begin());
  std::copy(state.particles.velocity_y_peculiar.begin(), state.particles.velocity_y_peculiar.end(), workspace.particle_velocity_y_peculiar.begin());
  std::copy(state.particles.velocity_z_peculiar.begin(), state.particles.velocity_z_peculiar.end(), workspace.particle_velocity_z_peculiar.begin());
  std::copy(state.particles.mass_code.begin(), state.particles.mass_code.end(), workspace.particle_mass_code.begin());
  return GravityParticleKernelView{
      .particle_index = workspace.gravity_particle_index,
      .position_x_comoving = workspace.particle_position_x_comoving,
      .position_y_comoving = workspace.particle_position_y_comoving,
      .position_z_comoving = workspace.particle_position_z_comoving,
      .velocity_x_peculiar = workspace.particle_velocity_x_peculiar,
      .velocity_y_peculiar = workspace.particle_velocity_y_peculiar,
      .velocity_z_peculiar = workspace.particle_velocity_z_peculiar,
      .mass_code = workspace.particle_mass_code,
      .source_particle_index_generation = state.particleIndexGeneration(),
  };
}

GravityParticleKernelView buildGravityParticleKernelViewAllParticlesDirect(
    SimulationState& state,
    TransientStepWorkspace& workspace) {
  (void)checkedLocalCount(
      state.particles.size(),
      kMaxLocalParticleCount,
      "particle",
      "buildGravityParticleKernelViewAllParticlesDirect");
  std::span<std::uint32_t> particle_index;
  if (workspace.gravity_particle_index_scratch.size() >= state.particles.size()) {
    particle_index = workspace.gravity_particle_index_scratch.first(
        state.particles.size());
  } else {
    // Compatibility path for standalone/core callers and for a within-step
    // population increase beyond the collectively admitted scratch extent.
    // Production scratch growth itself is coordinated by GravityRuntime.
    workspace.gravity_particle_index.resize(state.particles.size());
    particle_index = workspace.gravity_particle_index;
  }
  for (std::size_t i = 0; i < state.particles.size(); ++i) {
    particle_index[i] = checkedLocalParticleRow(i, "gravity all-particle view row");
  }
  return GravityParticleKernelView{
      .particle_index = particle_index,
      .position_x_comoving = state.particles.position_x_comoving,
      .position_y_comoving = state.particles.position_y_comoving,
      .position_z_comoving = state.particles.position_z_comoving,
      .velocity_x_peculiar = state.particles.velocity_x_peculiar,
      .velocity_y_peculiar = state.particles.velocity_y_peculiar,
      .velocity_z_peculiar = state.particles.velocity_z_peculiar,
      .mass_code = state.particles.mass_code,
      .source_particle_index_generation = state.particleIndexGeneration(),
  };
}

void scatterGravityParticleKernelView(const GravityParticleKernelView& view, SimulationState& state) {
  if (view.source_particle_index_generation != state.particleIndexGeneration()) {
    throw std::runtime_error("scatterGravityParticleKernelView: stale particle view generation");
  }
  if (!view.isConsistent()) {
    throw std::invalid_argument("scatterGravityParticleKernelView: view lanes have inconsistent sizes");
  }
  std::unordered_set<std::uint32_t> destinations;
  destinations.reserve(view.size());
  for (const std::uint32_t destination : view.particle_index) {
    if (destination >= state.particles.size()) {
      throw std::out_of_range("scatterGravityParticleKernelView: stale particle index");
    }
    if (!destinations.insert(destination).second) {
      throw std::invalid_argument("scatterGravityParticleKernelView: duplicate destination particle index");
    }
  }
  for (std::size_t i = 0; i < view.size(); ++i) {
    const auto destination = view.particle_index[i];
    state.particles.position_x_comoving[destination] = view.position_x_comoving[i];
    state.particles.position_y_comoving[destination] = view.position_y_comoving[i];
    state.particles.position_z_comoving[destination] = view.position_z_comoving[i];
    state.particles.velocity_x_peculiar[destination] = view.velocity_x_peculiar[i];
    state.particles.velocity_y_peculiar[destination] = view.velocity_y_peculiar[i];
    state.particles.velocity_z_peculiar[destination] = view.velocity_z_peculiar[i];
    state.particles.mass_code[destination] = view.mass_code[i];
  }
}

HydroCellKernelView buildHydroCellKernelView(
    const SimulationState& state,
    std::span<const std::uint32_t> active_cell_indices,
    TransientStepWorkspace& workspace) {
  if (!active_cell_indices.empty()) {
    state.requireGasCellIdentityMapCoversDenseRows("buildHydroCellKernelView");
  }
  workspace.hydro_cell_index.resize(active_cell_indices.size());
  workspace.hydro_gas_cell_id.resize(active_cell_indices.size());
  workspace.hydro_local_cell_row.resize(active_cell_indices.size());
  workspace.hydro_cell_center_x_comoving.resize(active_cell_indices.size());
  workspace.hydro_cell_center_y_comoving.resize(active_cell_indices.size());
  workspace.hydro_cell_center_z_comoving.resize(active_cell_indices.size());
  workspace.hydro_cell_mass_code.resize(active_cell_indices.size());
  workspace.hydro_cell_density_code.resize(active_cell_indices.size());
  workspace.hydro_cell_pressure_code.resize(active_cell_indices.size());

  for (std::size_t i = 0; i < active_cell_indices.size(); ++i) {
    const auto source = active_cell_indices[i];
    if (source >= state.cells.size()) {
      throw std::out_of_range("buildHydroCellKernelView: cell index out of range");
    }
    const auto gas_cell_id = state.gasCellIdForLocalRow(source);
    if (!gas_cell_id.has_value()) {
      throw std::runtime_error("buildHydroCellKernelView: active cell row is absent from GasCellIdentityMap");
    }
    workspace.hydro_cell_index[i] = source;
    workspace.hydro_gas_cell_id[i] = *gas_cell_id;
    workspace.hydro_local_cell_row[i] = source;
  }

  gatherSpan<double>(
      state.cells.center_x_comoving,
      active_cell_indices,
      workspace.hydro_cell_center_x_comoving);
  gatherSpan<double>(
      state.cells.center_y_comoving,
      active_cell_indices,
      workspace.hydro_cell_center_y_comoving);
  gatherSpan<double>(
      state.cells.center_z_comoving,
      active_cell_indices,
      workspace.hydro_cell_center_z_comoving);
  gatherSpan<double>(state.cells.mass_code, active_cell_indices, workspace.hydro_cell_mass_code);
  gatherSpan<double>(state.gas_cells.density_code, active_cell_indices, workspace.hydro_cell_density_code);
  gatherSpan<double>(state.gas_cells.pressure_code, active_cell_indices, workspace.hydro_cell_pressure_code);
  return HydroCellKernelView{
      .cell_index = workspace.hydro_cell_index,
      .gas_cell_id = workspace.hydro_gas_cell_id,
      .local_cell_row = workspace.hydro_local_cell_row,
      .center_x_comoving = workspace.hydro_cell_center_x_comoving,
      .center_y_comoving = workspace.hydro_cell_center_y_comoving,
      .center_z_comoving = workspace.hydro_cell_center_z_comoving,
      .mass_code = workspace.hydro_cell_mass_code,
      .density_code = workspace.hydro_cell_density_code,
      .pressure_code = workspace.hydro_cell_pressure_code,
      .source_cell_index_generation = state.cellIndexGeneration(),
      .source_gas_cell_identity_generation = state.gasCellIdentityGeneration(),
  };
}

void scatterHydroCellKernelView(const HydroCellKernelView& view, SimulationState& state) {
  if (view.source_cell_index_generation != state.cellIndexGeneration()) {
    throw std::runtime_error("scatterHydroCellKernelView: stale cell view generation");
  }
  if (view.source_gas_cell_identity_generation != state.gasCellIdentityGeneration()) {
    throw std::runtime_error("scatterHydroCellKernelView: stale gas-cell identity map generation");
  }
  if (!view.isConsistent()) {
    throw std::invalid_argument("scatterHydroCellKernelView: view lanes have inconsistent sizes");
  }
  if (view.size() > 0) {
    state.requireGasCellIdentityMapCoversDenseRows("scatterHydroCellKernelView");
  }

  std::vector<std::uint32_t> destinations(view.size());
  std::unordered_set<std::uint32_t> unique_destinations;
  unique_destinations.reserve(view.size());
  for (std::size_t i = 0; i < view.size(); ++i) {
    const auto destination = state.rowForGasCellId(view.gas_cell_id[i]);
    if (!destination.has_value() || *destination >= state.cells.size()) {
      throw std::out_of_range("scatterHydroCellKernelView: gas_cell_id does not resolve to a valid local cell row");
    }
    if (view.local_cell_row[i] != view.cell_index[i]) {
      throw std::invalid_argument("scatterHydroCellKernelView: local_cell_row/index mirrors disagree");
    }
    if (!unique_destinations.insert(*destination).second) {
      throw std::invalid_argument("scatterHydroCellKernelView: duplicate destination gas cell");
    }
    destinations[i] = *destination;
  }
  for (std::size_t i = 0; i < view.size(); ++i) {
    const std::uint32_t destination = destinations[i];
    state.cells.center_x_comoving[destination] = view.center_x_comoving[i];
    state.cells.center_y_comoving[destination] = view.center_y_comoving[i];
    state.cells.center_z_comoving[destination] = view.center_z_comoving[i];
    state.cells.mass_code[destination] = view.mass_code[i];
    state.gas_cells.density_code[destination] = view.density_code[i];
    state.gas_cells.pressure_code[destination] = view.pressure_code[i];
  }
}


}  // namespace cosmosim::core
