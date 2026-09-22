#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

#include <mpi.h>

#include "cosmosim/analysis/diagnostics.hpp"
#include "cosmosim/core/config.hpp"
#include "cosmosim/core/profiling.hpp"
#include "cosmosim/parallel/distributed_memory.hpp"
#include "cosmosim/workflows/runtime_services.hpp"

#include "../support/mpi_test_workspace.hpp"

namespace {

bool nearlyEqual(double lhs, double rhs, double tol = 1.0e-12) {
  return std::abs(lhs - rhs) <= tol * std::max({1.0, std::abs(lhs), std::abs(rhs)});
}

std::string readText(const std::filesystem::path& path) {
  std::ifstream in(path);
  std::ostringstream text;
  text << in.rdbuf();
  return text.str();
}

}  // namespace

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  cosmosim::parallel::MpiContext mpi;
  const int rank = mpi.worldRank();
  const int size = mpi.worldSize();
  if (size < 2 || size > 4) {
    MPI_Finalize();
    return 1;
  }

  const std::uint64_t active = rank == size - 1
      ? 0U
      : static_cast<std::uint64_t>(rank + 1);

  cosmosim::analysis::DiagnosticsBundle bundle;
  bundle.step_index = 7;
  bundle.scale_factor = 0.5;
  bundle.diagnostic_class = cosmosim::analysis::DiagnosticClass::kScienceLight;
  bundle.health.particle_count = active;
  bundle.health.cell_count = 2U * active;
  bundle.health.star_count = rank == 0 ? 1U : 0U;
  bundle.health.non_finite_particles = rank == 1 ? 2U : 0U;
  bundle.health.non_finite_cells = rank == 0 ? 1U : 0U;
  bundle.health.non_positive_particle_mass = rank == 1 ? 1U : 0U;
  bundle.health.ownership_invariants_ok = rank != 1;
  bundle.health.unique_particle_ids_ok = true;
  bundle.health.gravity_softening_sidecar_size_ok = rank != 0;

  const double contribution = static_cast<double>(active);
  for (std::size_t axis = 0; axis < 3; ++axis) {
    bundle.angular_momentum.total_l_code[axis] = contribution * static_cast<double>(axis + 1U);
    bundle.angular_momentum.gas_l_code[axis] = contribution * 0.5;
    bundle.angular_momentum.star_l_code[axis] = rank == 0 ? static_cast<double>(axis + 1U) : 0.0;
    bundle.angular_momentum.dark_matter_l_code[axis] = contribution * 2.0;
    bundle.angular_momentum.black_hole_l_code[axis] = 0.0;
  }

  bundle.star_formation_history = {
      {.scale_factor_center = 0.25, .formed_mass_code = contribution},
      {.scale_factor_center = 0.75, .formed_mass_code = 2.0 * contribution},
  };
  bundle.quicklook_grid_n = 2;
  bundle.xy_slice_density_code.assign(4, static_cast<double>(rank + 1));
  bundle.xy_slice_sample_count.assign(4, active);
  bundle.xy_projection_density_code.assign(4, 3.0 * contribution);
  bundle.power_spectrum = {{.k_center_code = 1.0, .power_code_volume = contribution, .mode_count = active}};
  bundle.records.push_back({
      .name = "power_spectrum",
      .tier = cosmosim::analysis::DiagnosticTier::kProvisionalScience,
      .maturity = cosmosim::analysis::DiagnosticMaturity::kProvisional,
      .scalability = cosmosim::analysis::DiagnosticScalability::kScalableFft,
      .executed = true,
      .policy_note = "local_test_value",
  });

  bundle.memory_report.totals.persistent_total_bytes =
      100U * static_cast<std::uint64_t>(rank + 1);

  cosmosim::analysis::reduceDiagnosticsBundleAcrossRanks(bundle, mpi);
  cosmosim::core::ProfilerSession profiler(false);
  cosmosim::workflows::RuntimeServices services{
      .mpi_context = mpi,
      .profiler = profiler,
      .memory_governor = nullptr,
      .console_reporter = nullptr,
      .deterministic_execution = true,
  };
  cosmosim::workflows::attachDistributedMemoryTelemetry(bundle.memory_report, services);

  const std::uint64_t expected_active_sum =
      static_cast<std::uint64_t>(size - 1) * static_cast<std::uint64_t>(size) / 2U;
  assert(bundle.globally_reduced);
  assert(bundle.contributing_rank_count == size);
  assert(bundle.health.particle_count == expected_active_sum);
  assert(bundle.health.cell_count == 2U * expected_active_sum);
  assert(bundle.health.star_count == 1U);
  assert(bundle.health.non_finite_particles == 2U);
  assert(bundle.health.non_finite_cells == 1U);
  assert(bundle.health.non_positive_particle_mass == 1U);
  assert(!bundle.health.ownership_invariants_ok);
  assert(bundle.health.ownership_invariants_failed_ranks == 1U);
  assert(bundle.health.unique_particle_ids_ok);
  assert(bundle.health.unique_particle_ids_failed_ranks == 0U);
  assert(!bundle.health.gravity_softening_sidecar_size_ok);
  assert(bundle.health.gravity_softening_sidecar_size_failed_ranks == 1U);

  for (std::size_t axis = 0; axis < 3; ++axis) {
    assert(nearlyEqual(bundle.angular_momentum.total_l_code[axis],
                       static_cast<double>(expected_active_sum) * static_cast<double>(axis + 1U)));
    assert(nearlyEqual(bundle.angular_momentum.dark_matter_l_code[axis],
                       2.0 * static_cast<double>(expected_active_sum)));
  }
  assert(nearlyEqual(bundle.star_formation_history[0].formed_mass_code,
                     static_cast<double>(expected_active_sum)));
  assert(nearlyEqual(bundle.star_formation_history[1].formed_mass_code,
                     2.0 * static_cast<double>(expected_active_sum)));

  double weighted_density_sum = 0.0;
  double sample_sum = 0.0;
  for (int r = 0; r < size - 1; ++r) {
    const double samples = static_cast<double>(r + 1);
    weighted_density_sum += samples * static_cast<double>(r + 1);
    sample_sum += samples;
  }
  const double expected_slice = weighted_density_sum / sample_sum;
  for (double value : bundle.xy_slice_density_code) {
    assert(nearlyEqual(value, expected_slice));
  }
  for (double value : bundle.xy_projection_density_code) {
    assert(nearlyEqual(value, 3.0 * static_cast<double>(expected_active_sum)));
  }

  assert(bundle.power_spectrum.empty());
  assert(!bundle.records.front().executed);
  assert(bundle.records.front().policy_note == "unsupported_under_mpi_requires_global_density_fft");

  const std::uint64_t expected_memory_sum =
      100U * static_cast<std::uint64_t>(size) * static_cast<std::uint64_t>(size + 1) / 2U;
  assert(bundle.memory_report.distributed.valid);
  assert(bundle.memory_report.distributed.global_sum_owned_bytes == expected_memory_sum);
  assert(bundle.memory_report.distributed.rank_max_owned_bytes == 100U * static_cast<std::uint64_t>(size));

  auto workspace = cosmosim::test_support::createMpiSharedWorkspace(
      "cosmosim_diagnostics_mpi_global_" + std::to_string(size));
  auto config = cosmosim::core::makeUnvalidatedSimulationConfigForTests();
  config.output.output_directory = workspace.root().string();
  config.output.run_name = "diagnostics_mpi_global";
  config.analysis.diagnostics_stem = "diagnostics";
  cosmosim::analysis::DiagnosticsEngine engine(config);
  if (mpi.isRoot()) {
    engine.writeBundle(bundle);
  }
  MPI_Barrier(MPI_COMM_WORLD);

  const auto diagnostic_dir = workspace.root() / "diagnostics_mpi_global" / "diagnostics";
  const auto json_path = diagnostic_dir / "diagnostics_science_light_step_00000007.json";
  assert(std::filesystem::exists(json_path));
  assert(!std::filesystem::exists(json_path.string() + ".part"));
  const auto sfr_path = workspace.root() / "diagnostics_mpi_global" / "sfr_history.csv";
  assert(std::filesystem::exists(sfr_path));
  assert(!std::filesystem::exists(sfr_path.string() + ".part"));
  if (mpi.isRoot()) {
    const std::string json = readText(json_path);
    assert(json.find("\"scope\": \"global\"") != std::string::npos);
    assert(json.find("\"contributing_rank_count\": " + std::to_string(size)) != std::string::npos);
    assert(json.find("unsupported_under_mpi_requires_global_density_fft") != std::string::npos);
  }

  MPI_Finalize();
  return 0;
}
