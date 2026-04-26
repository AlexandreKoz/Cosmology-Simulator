#include <cassert>
#include <filesystem>
#include <string>
#include <vector>

#include "cosmosim/core/config.hpp"

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

  const std::filesystem::path run_dir = std::filesystem::temp_directory_path() / "cosmosim_config_test" /
                                        frozen.config.output.run_name;
  cosmosim::core::writeNormalizedConfigSnapshot(frozen, run_dir);
  assert(std::filesystem::exists(run_dir / "normalized_config.param.txt"));
}

}  // namespace

int main() {
  const std::filesystem::path source_dir = COSMOSIM_SOURCE_DIR;
  checkExample(source_dir / "configs/cosmo_cube.param.txt", cosmosim::core::SimulationMode::kCosmoCube);
  checkExample(source_dir / "configs/zoom_in.param.txt", cosmosim::core::SimulationMode::kZoomIn);
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
  return 0;
}
