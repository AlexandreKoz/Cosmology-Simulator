#include "cosmosim/core/simulation_state.hpp"

#include <unordered_map>

namespace cosmosim::core {
namespace {
std::unordered_map<const void*, std::uint64_t> g_particle_view_generation;
std::unordered_map<const void*, std::uint64_t> g_cell_view_generation;
}  // namespace

void ActiveIndexSet::clear() {
  particle_indices.clear();
  cell_indices.clear();
}

std::size_t ParticleActiveView::size() const noexcept { return particle_id.size(); }
std::size_t CellActiveView::size() const noexcept { return center_x_comoving.size(); }
std::size_t GravityParticleKernelView::size() const noexcept { return particle_index.size(); }
std::size_t HydroCellKernelView::size() const noexcept { return cell_index.size(); }

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

  scratch.reset();
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
  g_particle_view_generation[static_cast<const void*>(workspace.gravity_particle_index.data())] =
      state.particleIndexGeneration();

  return GravityParticleKernelView{
      .particle_index = workspace.gravity_particle_index,
      .position_x_comoving = workspace.particle_position_x_comoving,
      .position_y_comoving = workspace.particle_position_y_comoving,
      .position_z_comoving = workspace.particle_position_z_comoving,
      .velocity_x_peculiar = workspace.particle_velocity_x_peculiar,
      .velocity_y_peculiar = workspace.particle_velocity_y_peculiar,
      .velocity_z_peculiar = workspace.particle_velocity_z_peculiar,
      .mass_code = workspace.particle_mass_code,
  };
}

void scatterGravityParticleKernelView(const GravityParticleKernelView& view, SimulationState& state) {
  const auto* view_key = static_cast<const void*>(view.particle_index.data());
  const auto it = g_particle_view_generation.find(view_key);
  if (it != g_particle_view_generation.end() && it->second != state.particleIndexGeneration()) {
    throw std::runtime_error("scatterGravityParticleKernelView: stale particle view generation");
  }
  for (std::size_t i = 0; i < view.size(); ++i) {
    const auto destination = view.particle_index[i];
    if (destination >= state.particles.size()) {
      throw std::out_of_range("scatterGravityParticleKernelView: stale particle index");
    }
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
  workspace.hydro_cell_index.resize(active_cell_indices.size());
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
    workspace.hydro_cell_index[i] = source;
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
  g_cell_view_generation[static_cast<const void*>(workspace.hydro_cell_index.data())] =
      state.cellIndexGeneration();

  return HydroCellKernelView{
      .cell_index = workspace.hydro_cell_index,
      .center_x_comoving = workspace.hydro_cell_center_x_comoving,
      .center_y_comoving = workspace.hydro_cell_center_y_comoving,
      .center_z_comoving = workspace.hydro_cell_center_z_comoving,
      .mass_code = workspace.hydro_cell_mass_code,
      .density_code = workspace.hydro_cell_density_code,
      .pressure_code = workspace.hydro_cell_pressure_code,
  };
}

void scatterHydroCellKernelView(const HydroCellKernelView& view, SimulationState& state) {
  const auto* view_key = static_cast<const void*>(view.cell_index.data());
  const auto it = g_cell_view_generation.find(view_key);
  if (it != g_cell_view_generation.end() && it->second != state.cellIndexGeneration()) {
    throw std::runtime_error("scatterHydroCellKernelView: stale cell view generation");
  }
  for (std::size_t i = 0; i < view.size(); ++i) {
    const auto destination = view.cell_index[i];
    if (destination >= state.cells.size()) {
      throw std::out_of_range("scatterHydroCellKernelView: stale cell index");
    }
    state.cells.center_x_comoving[destination] = view.center_x_comoving[i];
    state.cells.center_y_comoving[destination] = view.center_y_comoving[i];
    state.cells.center_z_comoving[destination] = view.center_z_comoving[i];
    state.cells.mass_code[destination] = view.mass_code[i];
    state.gas_cells.density_code[destination] = view.density_code[i];
    state.gas_cells.pressure_code[destination] = view.pressure_code[i];
  }
}

}  // namespace cosmosim::core
