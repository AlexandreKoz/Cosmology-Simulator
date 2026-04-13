#pragma once

#include <cstdint>
#include <span>
#include <string_view>
#include <vector>

#include "cosmosim/core/time_integration.hpp"

namespace cosmosim::physics {

struct TracerConfig {
  bool enabled = false;
  bool track_mass = true;
  double min_host_mass_code = 0.0;
};

struct TracerInjectionRequest {
  std::uint32_t tracer_particle_index = 0;
  std::uint32_t host_cell_index = 0;
  std::uint64_t parent_particle_id = 0;
  std::uint64_t injection_step = 0;
  double injected_mass_code = 0.0;
};

struct TracerUpdateCounters {
  std::uint64_t updated_tracers = 0;
  std::uint64_t skipped_inactive_host = 0;
  std::uint64_t skipped_low_host_mass = 0;
  std::uint64_t skipped_invalid_host = 0;
  double cumulative_absolute_mass_delta_code = 0.0;
};

class TracerModel {
 public:
  explicit TracerModel(TracerConfig config = {});

  [[nodiscard]] const TracerConfig& config() const noexcept;
  void inject(cosmosim::core::SimulationState& state, const TracerInjectionRequest& request) const;

  [[nodiscard]] TracerUpdateCounters updateMassFromHostCells(
      cosmosim::core::SimulationState& state,
      std::span<const std::uint32_t> active_cell_indices) const;

 private:
  TracerConfig m_config;
};

class TracerCallback final : public cosmosim::core::IntegrationCallback {
 public:
  explicit TracerCallback(TracerModel model);

  [[nodiscard]] std::string_view callbackName() const override;
  void onStage(cosmosim::core::StepContext& context) override;

  [[nodiscard]] const TracerUpdateCounters& lastUpdateCounters() const noexcept;

 private:
  TracerModel m_model;
  TracerUpdateCounters m_last_counters;
  std::vector<std::uint32_t> m_all_cells_cache;
};

}  // namespace cosmosim::physics
