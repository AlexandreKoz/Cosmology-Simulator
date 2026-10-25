#include <cassert>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "cosmosim/core/config.hpp"

#include "../support/test_temp_workspace.hpp"

namespace {

void checkExample(const std::filesystem::path& path, cosmosim::core::SimulationMode expected_mode) {
  const auto frozen = cosmosim::core::loadFrozenConfigFromFile(path);
  assert(frozen.config.mode.mode == expected_mode);
  assert(!frozen.config.output.run_name.empty());
  assert(!frozen.provenance.config_hash_hex.empty());
  assert(frozen.config.numerics.treepm_pm_grid > 0);
  assert(frozen.config.numerics.treepm_pm_grid_nx > 0);
  assert(frozen.config.numerics.treepm_pm_grid_ny > 0);
  assert(frozen.config.numerics.treepm_pm_grid_nz > 0);
  assert(frozen.config.numerics.treepm_asmth_cells > 0.0);
  assert(frozen.config.numerics.treepm_rcut_cells > 0.0);
  assert(frozen.config.numerics.treepm_update_cadence_steps >= 1);
  assert(frozen.normalized_text.find("treepm_pm_grid_nx = ") != std::string::npos);
  assert(frozen.normalized_text.find("treepm_pm_grid = ") == std::string::npos);

  const auto workspace =
      cosmosim::test_support::TestTempWorkspace::createProcessLocal(
          "cosmosim_config_test_" + frozen.config.output.run_name);
  const std::filesystem::path run_dir = workspace.root();
  cosmosim::core::writeNormalizedConfigSnapshot(frozen, run_dir);
  assert(std::filesystem::exists(run_dir / "normalized_config.param.txt"));
}

void checkEveryStarFormingProfileSelectsModel(const std::filesystem::path& config_root) {
  for (const auto& entry : std::filesystem::recursive_directory_iterator(config_root)) {
    if (!entry.is_regular_file() || entry.path().extension() != ".txt") continue;
    std::ifstream stream(entry.path());
    std::ostringstream content;
    content << stream.rdbuf();
    const std::string text = content.str();
    if (text.find("enable_star_formation = true") != std::string::npos) {
      assert(text.find("star_formation_model = ") != std::string::npos);
    }
  }
}

}  // namespace

int main() {
  const std::filesystem::path source_dir = COSMOSIM_SOURCE_DIR;
  checkEveryStarFormingProfileSelectsModel(source_dir / "configs");
  checkExample(source_dir / "configs/cosmo_cube.param.txt", cosmosim::core::SimulationMode::kCosmoCube);
#ifdef _OPENMP
  checkExample(
      source_dir / "configs/chui_firstlight_32_smoke.param.txt",
      cosmosim::core::SimulationMode::kCosmoCube);
  const auto firstlight = cosmosim::core::loadFrozenConfigFromFile(
      source_dir / "configs/chui_firstlight_32_smoke.param.txt");
#else
  {
    // A non-OpenMP build must fail closed on the example's four-thread team
    // request, exactly as core config documents.
    bool thread_team_rejected = false;
    try {
      (void)cosmosim::core::loadFrozenConfigFromFile(
          source_dir / "configs/chui_firstlight_32_smoke.param.txt");
    } catch (const cosmosim::core::ConfigError& error) {
      thread_team_rejected =
          std::string(error.what()).find(
              "parallel.omp_threads > 1 requires a binary built with OpenMP support") !=
          std::string::npos;
    }
    assert(thread_team_rejected);
  }
  // Reload with a serial team size so every remaining firstlight invariant is
  // still exercised on a non-OpenMP build.
  const auto firstlight = cosmosim::core::loadFrozenConfigFromString(
      [&] {
        std::ifstream firstlight_stream(
            source_dir / "configs/chui_firstlight_32_smoke.param.txt");
        std::ostringstream firstlight_text;
        firstlight_text << firstlight_stream.rdbuf();
        std::string serial_team_text = firstlight_text.str();
        const std::size_t team_position = serial_team_text.find("omp_threads = 4");
        assert(team_position != std::string::npos);
        serial_team_text.replace(team_position, std::string("omp_threads = 4").size(),
                                 "omp_threads = 1");
        return serial_team_text;
      }(),
      "chui_firstlight_32_smoke");
#endif
  const auto firstlight_reparsed = cosmosim::core::loadFrozenConfigFromString(
      firstlight.normalized_text, "chui_firstlight_32_smoke_normalized");
  assert(firstlight_reparsed.normalized_text == firstlight.normalized_text);
  assert(firstlight_reparsed.provenance.config_hash_hex ==
         firstlight.provenance.config_hash_hex);
  assert(firstlight.config.mode.ic_bridge_source_length_unit_to_si ==
         3.0856775814913673e22);
  assert(firstlight.config.cosmology.box_size_x_mpc_comoving ==
         14.761890702961235);
  assert(firstlight.config.cosmology.box_size_x_mpc_comoving *
             firstlight.config.cosmology.hubble_param ==
         10.0);
  assert(firstlight.config.numerics.a_begin == 0.04);
  assert(firstlight.config.numerics.a_end == 0.041);
  checkExample(source_dir / "configs/zoom_in.param.txt", cosmosim::core::SimulationMode::kZoomIn);
  assert(cosmosim::core::loadFrozenConfigFromFile(source_dir / "configs/zoom_in.param.txt").config.physics.star_formation_model == cosmosim::core::StarFormationModelKind::kAdaptiveBoundJeans);
  assert(cosmosim::core::loadFrozenConfigFromFile(source_dir / "configs/cosmo_cube.param.txt").config.physics.star_formation_model == cosmosim::core::StarFormationModelKind::kEffectiveMultiphaseTngLike);
  checkExample(source_dir / "configs/isolated_galaxy.param.txt", cosmosim::core::SimulationMode::kIsolatedGalaxy);
  checkExample(source_dir / "configs/isolated_cluster.param.txt", cosmosim::core::SimulationMode::kIsolatedCluster);
  checkExample(
      source_dir / "configs/release/release_smoke_cosmo_cube.param.txt",
      cosmosim::core::SimulationMode::kCosmoCube);
  checkExample(
      source_dir / "configs/release/release_smoke_zoom_in.param.txt",
      cosmosim::core::SimulationMode::kZoomIn);
  checkExample(
      source_dir / "configs/release/release_smoke_isolated_galaxy.param.txt",
      cosmosim::core::SimulationMode::kIsolatedGalaxy);
  {
    // Measured fast TreePM continuation deck (docs/tree_pm_coupling.md).
    const auto fast = cosmosim::core::loadFrozenConfigFromFile(
        source_dir / "configs/treepm_evolved_snapshot_fast.param.txt");
    const auto& n = fast.config.numerics;
    assert(n.treepm_pm_grid_nx == 128 && n.treepm_pm_grid_ny == 128 && n.treepm_pm_grid_nz == 128);
    assert(n.treepm_adaptive_acceptance_enabled);
    assert(n.treepm_adaptive_error_estimator == cosmosim::core::TreePmAdaptiveErrorEstimator::kRelativeMonopole);
    assert(n.treepm_cutoff_node_policy == cosmosim::core::TreePmCutoffNodePolicy::kComDistance);
    assert(n.treepm_tree_relative_force_tolerance == 0.0025);
    assert(n.treepm_adaptive_maximum_opening_angle == 1.0);
    assert(n.hierarchical_max_rung == 0);
    assert(n.treepm_group_walk_enabled);
    assert(!n.treepm_leaf_simd_enabled && !n.treepm_full_residual_lookup_enabled);
    (void)n;
  }
  {
    // Unqualified DMO performance profile (docs/tree_pm_coupling.md). Every
    // performance choice is explicit in the deck, never a changed default.
    const auto perf = cosmosim::core::loadFrozenConfigFromFile(
        source_dir / "configs/performance/dmo_128_pm256_performance_unqualified.param.txt");
    const auto& n = perf.config.numerics;
    assert(n.treepm_pm_grid_nx == 256 && n.treepm_pm_grid_ny == 256 && n.treepm_pm_grid_nz == 256);
    assert(n.treepm_group_walk_enabled && n.treepm_group_walk_target_size == 32);
    assert(n.treepm_tree_max_leaf_size == 8);
    assert(n.treepm_tree_relative_force_tolerance == 0.0025);
    assert(n.gravity_timestep_eta == 0.2);
    assert(!perf.config.parallel.decomposition_spatial_work_enabled);
    assert(perf.config.analysis.run_health_interval_steps == 8);
    assert(perf.config.output.restart_interval_steps == 1024);
    assert(perf.config.output.snapshot_verification ==
           cosmosim::core::SnapshotVerificationPolicy::kStructural);
    (void)n;
  }
  return 0;
}
