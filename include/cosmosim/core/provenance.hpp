#pragma once

#include <cstdint>
#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

namespace cosmosim::core {

// Canonical provenance payload written alongside run outputs.
struct ProvenanceRecord {
  std::string schema_version = "provenance_v3";
  std::string git_sha = "unknown";
  std::string compiler_id = "unknown";
  std::string compiler_version = "unknown";
  std::string build_preset = "manual";
  std::string enabled_features;
  std::string config_hash_hex;
  std::string timestamp_utc;
  std::string hardware_summary;
  int author_rank = 0;
  int gravity_treepm_pm_grid = 0;
  std::string gravity_treepm_assignment_scheme = "unknown";
  bool gravity_treepm_window_deconvolution = false;
  double gravity_treepm_asmth_cells = 0.0;
  double gravity_treepm_rcut_cells = 0.0;
  double gravity_treepm_mesh_spacing_mpc_comoving = 0.0;
  double gravity_treepm_split_scale_mpc_comoving = 0.0;
  double gravity_treepm_cutoff_radius_mpc_comoving = 0.0;
  int gravity_treepm_update_cadence_steps = 0;
  std::string gravity_treepm_pm_decomposition_mode = "slab";
  std::uint64_t gravity_treepm_tree_exchange_batch_bytes = 0;
  std::string gravity_softening_policy = "unspecified";
  std::string gravity_softening_kernel = "unknown";
  double gravity_softening_epsilon_kpc_comoving = 0.0;
  std::string gravity_pm_fft_backend = "unknown";
  std::uint64_t gravity_treepm_decomposition_epoch = 0;
  int gravity_treepm_restart_world_size = 1;
  std::string gravity_treepm_restart_pm_grid = "0x0x0";
  std::string gravity_treepm_restart_slab_signature = "";
  std::uint64_t gravity_treepm_restart_kick_opportunity = 0;
  std::uint64_t gravity_treepm_restart_field_version = 0;
  std::string gravity_treepm_long_range_restart_policy = "deterministic_rebuild";
  std::string zoom_long_range_strategy = "disabled";
  double zoom_region_center_x_mpc_comoving = 0.0;
  double zoom_region_center_y_mpc_comoving = 0.0;
  double zoom_region_center_z_mpc_comoving = 0.0;
  double zoom_region_radius_mpc_comoving = 0.0;
  std::string zoom_focused_pm_grid = "0x0x0";
  double zoom_contamination_radius_mpc_comoving = 0.0;
};

[[nodiscard]] std::string collectCompilerId();
[[nodiscard]] std::string collectCompilerVersion();
[[nodiscard]] std::string collectHardwareSummary();
[[nodiscard]] std::string utcTimestampNowIso8601();

// Deterministic hashing for normalized configuration text.
[[nodiscard]] std::uint64_t stableConfigHash(const std::string& normalized_config_text);
[[nodiscard]] std::string stableConfigHashHex(const std::string& normalized_config_text);

[[nodiscard]] ProvenanceRecord makeProvenanceRecord(
    const std::string& config_hash_hex,
    const std::string& git_sha,
    int rank = 0);

void writeProvenanceRecord(
    const ProvenanceRecord& record,
    const std::filesystem::path& run_directory,
    const std::string& file_name = "provenance.meta.txt");

[[nodiscard]] ProvenanceRecord readProvenanceRecord(
    const std::filesystem::path& run_directory,
    const std::string& file_name = "provenance.meta.txt");

[[nodiscard]] std::string serializeProvenanceRecord(const ProvenanceRecord& record);
[[nodiscard]] ProvenanceRecord deserializeProvenanceRecord(std::string_view text);

}  // namespace cosmosim::core
