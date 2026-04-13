#include <cassert>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>
#include <vector>

#include "cosmosim/core/profiling.hpp"
#include "cosmosim/core/time_integration.hpp"

namespace {

class CountingCallback final : public cosmosim::core::IntegrationCallback {
 public:
  std::string_view callbackName() const override { return "counting_callback"; }

  void onStage(cosmosim::core::StepContext& context) override {
    ++stage_calls;
    if (context.stage == cosmosim::core::IntegrationStage::kHydroUpdate) {
      ++hydro_stage_calls;
    }
  }

  std::uint64_t stage_calls = 0;
  std::uint64_t hydro_stage_calls = 0;
};

std::string readFile(const std::filesystem::path& path) {
  std::ifstream in(path);
  assert(in.good());
  return std::string(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
}

void testMiniRunProfileReportGeneration() {
  cosmosim::core::SimulationState state;
  state.particles.resize(8);
  state.cells.resize(4);

  std::vector<std::uint32_t> active_particles = {0, 1, 2, 3};
  std::vector<std::uint32_t> active_cells = {0, 1};

  cosmosim::core::ActiveSetDescriptor active_set{
      .particle_indices = active_particles,
      .cell_indices = active_cells,
      .particles_are_subset = true,
      .cells_are_subset = true,
  };

  cosmosim::core::IntegratorState integrator_state;
  integrator_state.dt_time_code = 0.01;

  CountingCallback callback;
  cosmosim::core::StepOrchestrator orchestrator;
  orchestrator.registerCallback(callback);

  cosmosim::core::ProfilerSession profiler(true);
  orchestrator.executeSingleStep(state, integrator_state, active_set, nullptr, nullptr, nullptr, &profiler);

  assert(callback.stage_calls == cosmosim::core::StageScheduler::kickDriftKickOrder().size());
  assert(callback.hydro_stage_calls == 1);

  assert(profiler.counters().count("step_invocations") == 1);
  assert(profiler.counters().count("active_particles") == 4);
  assert(profiler.counters().count("active_cells") == 2);
  assert(profiler.counters().count("stage.hydro_update.invocations") == 1);

  const auto temp_dir = std::filesystem::temp_directory_path();
  const auto json_path = temp_dir / "cosmosim_profile_integration.json";
  const auto csv_path = temp_dir / "cosmosim_profile_integration.csv";

  cosmosim::core::writeProfilerReportJson(profiler, json_path);
  cosmosim::core::writeProfilerReportCsv(profiler, csv_path);

  const std::string json_text = readFile(json_path);
  const std::string csv_text = readFile(csv_path);

  assert(json_text.find("\"schema_version\": 1") != std::string::npos);
  assert(json_text.find("\"stage.hydro_update.invocations\": 1") != std::string::npos);
  assert(csv_text.find("root/step_orchestrator.execute_single_step") != std::string::npos);

  std::filesystem::remove(json_path);
  std::filesystem::remove(csv_path);
}

}  // namespace

int main() {
  testMiniRunProfileReportGeneration();
  return 0;
}
