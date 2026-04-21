#include <cassert>
#include <filesystem>
#include <string>

#include "cosmosim/core/provenance.hpp"

int main() {
  const auto run_directory = std::filesystem::temp_directory_path() / "cosmosim_provenance_roundtrip";
  std::filesystem::remove_all(run_directory);

  cosmosim::core::ProvenanceRecord in;
  in.git_sha = "abc123def";
  in.compiler_id = "test_compiler";
  in.compiler_version = "1.2.3";
  in.build_preset = "test_preset";
  in.enabled_features = "mpi=0,hdf5=0,fftw=0,cuda=0,python=0";
  in.config_hash_hex = "0011223344556677";
  in.timestamp_utc = "2026-04-05T00:00:00Z";
  in.hardware_summary = "logical_threads=8";
  in.author_rank = 0;
  in.gravity_treepm_pm_grid = 64;
  in.gravity_treepm_pm_grid_nx = 64;
  in.gravity_treepm_pm_grid_ny = 48;
  in.gravity_treepm_pm_grid_nz = 32;
  in.gravity_treepm_assignment_scheme = "tsc";
  in.gravity_treepm_window_deconvolution = true;
  in.gravity_treepm_asmth_cells = 1.25;
  in.gravity_treepm_rcut_cells = 4.5;
  in.gravity_treepm_mesh_spacing_mpc_comoving = 0.5;
  in.gravity_treepm_mesh_spacing_x_mpc_comoving = 0.5;
  in.gravity_treepm_mesh_spacing_y_mpc_comoving = 0.625;
  in.gravity_treepm_mesh_spacing_z_mpc_comoving = 0.75;
  in.gravity_treepm_split_scale_mpc_comoving = 0.625;
  in.gravity_treepm_cutoff_radius_mpc_comoving = 2.25;
  in.gravity_treepm_update_cadence_steps = 2;
  in.gravity_treepm_pm_decomposition_mode = "slab";
  in.gravity_treepm_tree_exchange_batch_bytes = 4194304ULL;
  in.gravity_softening_policy = "comoving_fixed";
  in.gravity_softening_kernel = "plummer";
  in.gravity_softening_epsilon_kpc_comoving = 1.0;
  in.gravity_pm_fft_backend = "fftw3";

  cosmosim::core::writeProvenanceRecord(in, run_directory);

  const auto out = cosmosim::core::readProvenanceRecord(run_directory);
  assert(out.schema_version == in.schema_version);
  assert(out.git_sha == in.git_sha);
  assert(out.compiler_id == in.compiler_id);
  assert(out.compiler_version == in.compiler_version);
  assert(out.build_preset == in.build_preset);
  assert(out.enabled_features == in.enabled_features);
  assert(out.config_hash_hex == in.config_hash_hex);
  assert(out.timestamp_utc == in.timestamp_utc);
  assert(out.hardware_summary == in.hardware_summary);
  assert(out.author_rank == in.author_rank);
  assert(out.gravity_treepm_pm_grid == in.gravity_treepm_pm_grid);
  assert(out.gravity_treepm_pm_grid_nx == in.gravity_treepm_pm_grid_nx);
  assert(out.gravity_treepm_pm_grid_ny == in.gravity_treepm_pm_grid_ny);
  assert(out.gravity_treepm_pm_grid_nz == in.gravity_treepm_pm_grid_nz);
  assert(out.gravity_treepm_assignment_scheme == in.gravity_treepm_assignment_scheme);
  assert(out.gravity_treepm_window_deconvolution == in.gravity_treepm_window_deconvolution);
  assert(out.gravity_treepm_asmth_cells == in.gravity_treepm_asmth_cells);
  assert(out.gravity_treepm_rcut_cells == in.gravity_treepm_rcut_cells);
  assert(out.gravity_treepm_mesh_spacing_mpc_comoving == in.gravity_treepm_mesh_spacing_mpc_comoving);
  assert(out.gravity_treepm_mesh_spacing_x_mpc_comoving == in.gravity_treepm_mesh_spacing_x_mpc_comoving);
  assert(out.gravity_treepm_mesh_spacing_y_mpc_comoving == in.gravity_treepm_mesh_spacing_y_mpc_comoving);
  assert(out.gravity_treepm_mesh_spacing_z_mpc_comoving == in.gravity_treepm_mesh_spacing_z_mpc_comoving);
  assert(out.gravity_treepm_split_scale_mpc_comoving == in.gravity_treepm_split_scale_mpc_comoving);
  assert(out.gravity_treepm_cutoff_radius_mpc_comoving == in.gravity_treepm_cutoff_radius_mpc_comoving);
  assert(out.gravity_treepm_update_cadence_steps == in.gravity_treepm_update_cadence_steps);
  assert(out.gravity_treepm_pm_decomposition_mode == in.gravity_treepm_pm_decomposition_mode);
  assert(out.gravity_treepm_tree_exchange_batch_bytes == in.gravity_treepm_tree_exchange_batch_bytes);
  assert(out.gravity_softening_policy == in.gravity_softening_policy);
  assert(out.gravity_softening_kernel == in.gravity_softening_kernel);
  assert(out.gravity_softening_epsilon_kpc_comoving == in.gravity_softening_epsilon_kpc_comoving);
  assert(out.gravity_pm_fft_backend == in.gravity_pm_fft_backend);

  std::filesystem::remove_all(run_directory);
  return 0;
}
