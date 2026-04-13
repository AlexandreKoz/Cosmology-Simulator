#include "cosmosim/physics/tracer_support.hpp"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>

#include "cosmosim/core/build_config.hpp"

namespace cosmosim::physics {

TracerModel::TracerModel(TracerConfig config) : m_config(config) {}

const TracerConfig& TracerModel::config() const noexcept { return m_config; }

void TracerModel::inject(cosmosim::core::SimulationState& state, const TracerInjectionRequest& request) const {
  if (!m_config.enabled) {
    throw std::runtime_error("TracerModel.inject: tracers are disabled in config");
  }
  if (request.tracer_particle_index >= state.particles.size()) {
    throw std::out_of_range("TracerModel.inject: tracer_particle_index out of range");
  }
  if (request.host_cell_index >= state.cells.size()) {
    throw std::out_of_range("TracerModel.inject: host_cell_index out of range");
  }
  if (request.injected_mass_code < 0.0) {
    throw std::invalid_argument("TracerModel.inject: injected_mass_code must be >= 0");
  }
  const std::uint32_t tracer_species = static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kTracer);
  if (state.particle_sidecar.species_tag[request.tracer_particle_index] != tracer_species) {
    throw std::invalid_argument("TracerModel.inject: tracer particle must use species tag tracer");
  }

  const double host_mass_code = state.cells.mass_code[request.host_cell_index];
  if (host_mass_code <= 0.0) {
    throw std::invalid_argument("TracerModel.inject: host mass must be > 0");
  }

  const std::size_t local_index = state.tracers.size();
  state.tracers.resize(local_index + 1);
  state.tracers.particle_index[local_index] = request.tracer_particle_index;
  state.tracers.host_cell_index[local_index] = request.host_cell_index;
  state.tracers.parent_particle_id[local_index] = request.parent_particle_id;
  state.tracers.injection_step[local_index] = request.injection_step;
  state.tracers.mass_fraction_of_host[local_index] = request.injected_mass_code / host_mass_code;
  state.tracers.last_host_mass_code[local_index] = host_mass_code;
  state.tracers.cumulative_exchanged_mass_code[local_index] = 0.0;

  if (m_config.track_mass) {
    state.particles.mass_code[request.tracer_particle_index] = request.injected_mass_code;
  }

  std::ostringstream sidecar_text;
  sidecar_text << "tracer_count=" << state.tracers.size() << '\n';
  sidecar_text << "track_mass=" << (m_config.track_mass ? 1 : 0) << '\n';
  sidecar_text << "min_host_mass_code=" << m_config.min_host_mass_code << '\n';

  cosmosim::core::ModuleSidecarBlock sidecar;
  sidecar.module_name = "tracer_support";
  const std::string serialized = sidecar_text.str();
  sidecar.payload.resize(serialized.size());
  for (std::size_t i = 0; i < serialized.size(); ++i) {
    sidecar.payload[i] = static_cast<std::byte>(serialized[i]);
  }
  state.sidecars.upsert(std::move(sidecar));
}

TracerUpdateCounters TracerModel::updateMassFromHostCells(
    cosmosim::core::SimulationState& state,
    std::span<const std::uint32_t> active_cell_indices) const {
  TracerUpdateCounters counters;
  if (!m_config.enabled || !m_config.track_mass) {
    return counters;
  }

  std::vector<std::uint8_t> active_cell_mask;
  if (!active_cell_indices.empty()) {
    active_cell_mask.assign(state.cells.size(), 0U);
    for (const std::uint32_t cell_index : active_cell_indices) {
      if (cell_index < active_cell_mask.size()) {
        active_cell_mask[cell_index] = 1U;
      }
    }
  }

  for (std::size_t tracer_row = 0; tracer_row < state.tracers.size(); ++tracer_row) {
    const std::uint32_t tracer_particle_index = state.tracers.particle_index[tracer_row];
    const std::uint32_t host_cell_index = state.tracers.host_cell_index[tracer_row];

    if (host_cell_index >= state.cells.size() || tracer_particle_index >= state.particles.size()) {
      ++counters.skipped_invalid_host;
      continue;
    }
    if (!active_cell_mask.empty() && active_cell_mask[host_cell_index] == 0U) {
      ++counters.skipped_inactive_host;
      continue;
    }

    const double host_mass_code = state.cells.mass_code[host_cell_index];
    if (host_mass_code <= m_config.min_host_mass_code) {
      ++counters.skipped_low_host_mass;
      continue;
    }

    const double new_tracer_mass_code = state.tracers.mass_fraction_of_host[tracer_row] * host_mass_code;
    const double old_tracer_mass_code = state.particles.mass_code[tracer_particle_index];

    state.particles.mass_code[tracer_particle_index] = new_tracer_mass_code;
    state.tracers.cumulative_exchanged_mass_code[tracer_row] +=
        std::abs(new_tracer_mass_code - old_tracer_mass_code);
    state.tracers.last_host_mass_code[tracer_row] = host_mass_code;

    ++counters.updated_tracers;
    counters.cumulative_absolute_mass_delta_code +=
        std::abs(new_tracer_mass_code - old_tracer_mass_code);
  }

  std::ostringstream sidecar_text;
  sidecar_text << "updated_tracers=" << counters.updated_tracers << '\n';
  sidecar_text << "skipped_inactive_host=" << counters.skipped_inactive_host << '\n';
  sidecar_text << "skipped_low_host_mass=" << counters.skipped_low_host_mass << '\n';
  sidecar_text << "skipped_invalid_host=" << counters.skipped_invalid_host << '\n';
  sidecar_text << "cumulative_absolute_mass_delta_code=" << counters.cumulative_absolute_mass_delta_code << '\n';

  cosmosim::core::ModuleSidecarBlock sidecar;
  sidecar.module_name = "tracer_support";
  const std::string serialized = sidecar_text.str();
  sidecar.payload.resize(serialized.size());
  for (std::size_t i = 0; i < serialized.size(); ++i) {
    sidecar.payload[i] = static_cast<std::byte>(serialized[i]);
  }
  state.sidecars.upsert(std::move(sidecar));

  return counters;
}

TracerCallback::TracerCallback(TracerModel model) : m_model(std::move(model)) {}

std::string_view TracerCallback::callbackName() const { return "tracer_support"; }

void TracerCallback::onStage(cosmosim::core::StepContext& context) {
#if !COSMOSIM_ENABLE_TRACERS
  static_cast<void>(context);
  m_last_counters = {};
#else
  if (context.stage != cosmosim::core::IntegrationStage::kSourceTerms) {
    return;
  }

  if (!m_model.config().enabled) {
    m_last_counters = {};
    return;
  }

  if (context.active_set.cells_are_subset) {
    m_last_counters = m_model.updateMassFromHostCells(context.state, context.active_set.cell_indices);
    return;
  }

  m_all_cells_cache.resize(context.state.cells.size());
  for (std::size_t i = 0; i < m_all_cells_cache.size(); ++i) {
    m_all_cells_cache[i] = static_cast<std::uint32_t>(i);
  }
  m_last_counters = m_model.updateMassFromHostCells(context.state, m_all_cells_cache);
#endif
}

const TracerUpdateCounters& TracerCallback::lastUpdateCounters() const noexcept { return m_last_counters; }

}  // namespace cosmosim::physics
