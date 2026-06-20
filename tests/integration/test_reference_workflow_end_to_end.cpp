#include <cassert>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <array>
#include <cmath>
#include <cstdint>
#include <unordered_map>

#include "cosmosim/io/restart_checkpoint.hpp"

#include "cosmosim/cosmosim.hpp"
#include "cosmosim/core/build_config.hpp"

[[nodiscard]] std::string readFile(const std::filesystem::path& path) {
  std::ifstream in(path);
  assert(in.good());
  return std::string(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
}

[[nodiscard]] std::uint64_t xorRangeOneToN(std::uint64_t n) {
  switch (n & 3ULL) {
    case 0ULL: return n;
    case 1ULL: return 1ULL;
    case 2ULL: return n + 1ULL;
    default: return 0ULL;
  }
}

int main() {
  auto buildConfigText = [](int cadence_steps, const std::string& run_name, std::string_view assignment_scheme) {
    std::stringstream stream;
    stream << "schema_version = 1\n\n";
    stream << "[mode]\n";
    stream << "mode = zoom_in\n";
    stream << "ic_file = generated\n";
    stream << "zoom_high_res_region = false\n\n";
    stream << "[numerics]\n";
    stream << "time_begin_code = 0.01\n";
    stream << "time_end_code = 0.0102\n";
    stream << "max_global_steps = 2\n";
    stream << "hierarchical_max_rung = 1\n\n";
    stream << "treepm_pm_grid = 9\n";
    stream << "treepm_asmth_cells = 1.75\n";
    stream << "treepm_rcut_cells = 6.0\n";
    stream << "treepm_assignment_scheme = " << assignment_scheme << '\n';
    stream << "treepm_update_cadence_steps = " << cadence_steps << "\n\n";
    stream << "[output]\n";
    stream << "run_name = " << run_name << '\n';
    stream << "output_directory = integration_outputs\n";
    stream << "output_stem = snapshot\n";
    stream << "restart_stem = restart\n";
    return stream.str();
  };

  auto makeDecoupledGasCellState = [](const std::array<std::uint32_t, 3>& row_to_physical) {
    cosmosim::core::SimulationState state;
    state.resizeParticles(2);
    state.resizeCells(3);
    state.resizePatches(1);

    state.particle_sidecar.particle_id = {5001U, 6001U};
    state.particle_sidecar.sfc_key = {100U, 200U};
    state.particle_sidecar.species_tag = {
        static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kGas),
        static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kDarkMatter)};
    state.particle_sidecar.owning_rank = {0U, 0U};
    state.particle_sidecar.particle_flags = {0U, 0U};
    state.particles.position_x_comoving = {0.15, 0.45};
    state.particles.position_y_comoving = {0.05, 0.05};
    state.particles.position_z_comoving = {0.05, 0.05};
    state.particles.velocity_x_peculiar = {0.0, 0.0};
    state.particles.velocity_y_peculiar = {0.0, 0.0};
    state.particles.velocity_z_peculiar = {0.0, 0.0};
    state.particles.mass_code = {3.0, 10.0};
    state.particles.time_bin = {0U, 0U};
    state.species.count_by_species.fill(0U);
    state.species.count_by_species[static_cast<std::size_t>(cosmosim::core::ParticleSpecies::kGas)] = 1U;
    state.species.count_by_species[static_cast<std::size_t>(cosmosim::core::ParticleSpecies::kDarkMatter)] = 1U;
    state.rebuildSpeciesIndex();

    state.patches.patch_id[0] = 9001U;
    state.patches.level[0] = 0U;
    state.patches.first_cell[0] = 0U;
    state.patches.cell_count[0] = 3U;
    state.patches.owning_rank[0] = 0U;
    state.patches.origin_x_comoving[0] = 0.0;
    state.patches.origin_y_comoving[0] = 0.0;
    state.patches.origin_z_comoving[0] = 0.0;
    state.patches.extent_x_comoving[0] = 0.3;
    state.patches.extent_y_comoving[0] = 0.1;
    state.patches.extent_z_comoving[0] = 0.1;
    state.patches.cell_dim_x[0] = 3U;
    state.patches.cell_dim_y[0] = 1U;
    state.patches.cell_dim_z[0] = 1U;

    struct CellSeed {
      std::uint64_t gas_cell_id;
      std::optional<std::uint64_t> parent_particle_id;
      double density;
      double pressure;
      double internal_energy;
      double velocity_x;
    };
    const std::array<CellSeed, 3> physical_cells{{
        {8001U, std::nullopt, 1.0, 1.2, 1.8, -0.1},
        {8002U, 5001U, 1.1, 1.3, 1.9, 0.0},
        {8003U, 5001U, 0.9, 1.4, 2.0, 0.1},
    }};

    std::vector<cosmosim::core::GasCellIdentityRecord> records;
    records.reserve(physical_cells.size());
    for (std::uint32_t row = 0; row < row_to_physical.size(); ++row) {
      const std::uint32_t physical = row_to_physical[row];
      assert(physical < physical_cells.size());
      const CellSeed& seed = physical_cells[physical];
      state.cells.center_x_comoving[row] = 0.05 + 0.1 * static_cast<double>(physical);
      state.cells.center_y_comoving[row] = 0.05;
      state.cells.center_z_comoving[row] = 0.05;
      state.cells.mass_code[row] = 1.0;
      state.cells.time_bin[row] = 0U;
      state.cells.patch_index[row] = 0U;
      state.gas_cells.density_code[row] = seed.density;
      state.gas_cells.pressure_code[row] = seed.pressure;
      state.gas_cells.internal_energy_code[row] = seed.internal_energy;
      state.gas_cells.temperature_code[row] = 1.0;
      state.gas_cells.sound_speed_code[row] = 1.0;
      state.gas_cells.velocity_x_peculiar[row] = seed.velocity_x;
      state.gas_cells.velocity_y_peculiar[row] = 0.0;
      state.gas_cells.velocity_z_peculiar[row] = 0.0;
      records.push_back(cosmosim::core::GasCellIdentityRecord{
          .gas_cell_id = seed.gas_cell_id,
          .parent_particle_id = seed.parent_particle_id,
          .owning_patch_id = 9001U,
          .local_cell_row = row,
      });
    }
    state.replaceGasCellIdentityRecords(std::move(records));
    assert(state.gasCellIdentityMapMatchesSidecarLanes());
    assert(state.validateOwnershipInvariants());
    return state;
  };

  auto gasStateByStableId = [](const cosmosim::core::SimulationState& state) {
    std::unordered_map<std::uint64_t, std::array<double, 7>> values;
    state.requireGasCellIdentityMapCoversDenseRows("reference workflow H2 stable-id comparison");
    for (std::uint32_t row = 0; row < state.cells.size(); ++row) {
      const auto* record = state.gas_cell_identity.findByLocalRow(row);
      assert(record != nullptr);
      values.emplace(record->gas_cell_id, std::array<double, 7>{
          state.cells.mass_code[row],
          state.cells.center_x_comoving[row],
          state.gas_cells.density_code[row],
          state.gas_cells.pressure_code[row],
          state.gas_cells.internal_energy_code[row],
          state.gas_cells.velocity_x_peculiar[row],
          static_cast<double>(state.cells.time_bin[row])});
    }
    return values;
  };

  std::stringstream stream;
  stream << buildConfigText(1, "reference_integration_test", "cic");

  const cosmosim::core::FrozenConfig frozen =
      cosmosim::core::loadFrozenConfigFromString(stream.str(), "test_reference_workflow_end_to_end");

  cosmosim::workflows::ReferenceWorkflowRunner runner(frozen);
  const std::filesystem::path output_dir =
      std::filesystem::temp_directory_path() / "cosmosim_reference_workflow_test";
  const cosmosim::workflows::ReferenceWorkflowReport report =
      runner.run(output_dir, cosmosim::workflows::ReferenceWorkflowOptions{.write_outputs = false});

  const std::filesystem::path expected_run_dir = output_dir / "reference_integration_test";

  assert(report.config_compatible);
  assert(report.schema_compatible);
  assert(report.canonical_stage_order);
  assert(report.stage_sequence.size() == cosmosim::core::StageScheduler::kickDriftKickOrder().size() * 2U);
  assert(report.stage_sequence.front() == "gravity_kick_pre");
  assert(report.stage_sequence.back() == "output_check");
  assert(report.completed_steps == 2);
  assert(report.world_size == 1);
  assert(report.world_rank == 0);
  assert(report.local_particle_count == 42);
  assert(report.global_particle_count == 42);
  assert(report.local_cell_count == 6);
  assert(report.global_cell_count == 6);
  assert(report.local_particle_id_sum == (42ULL * 43ULL) / 2ULL);
  assert(report.global_particle_id_sum == (42ULL * 43ULL) / 2ULL);
  assert(report.local_particle_id_xor == xorRangeOneToN(42));
  assert(report.global_particle_id_xor == xorRangeOneToN(42));
  assert(report.treepm_pm_grid == 9);
  assert(report.treepm_pm_grid_nx == 9);
  assert(report.treepm_pm_grid_ny == 9);
  assert(report.treepm_pm_grid_nz == 9);
  assert(report.treepm_pm_grid_shape == "9x9x9");
  assert(report.treepm_update_cadence_steps == 1);
  assert(report.treepm_long_range_refresh_count == 3);
  assert(report.treepm_long_range_reuse_count == 0);
  assert(report.treepm_cadence_records.size() == 3);
  for (const auto& record : report.treepm_cadence_records) {
    assert(record.refreshed_long_range_field);
    assert(record.field_built_step_index == record.step_index);
  }
  assert(report.run_directory == expected_run_dir);
  assert(report.normalized_config_snapshot_written);
  assert(std::filesystem::exists(report.normalized_config_snapshot_path));
  assert(std::filesystem::exists(report.profiler_json_path));
  assert(std::filesystem::exists(report.profiler_csv_path));
  assert(std::filesystem::exists(report.operational_report_json_path));

  std::ifstream op_in(report.operational_report_json_path);
  const std::string op_text((std::istreambuf_iterator<char>(op_in)), std::istreambuf_iterator<char>());
  op_in.close();
  assert(op_text.find("\"event_kind\": \"config.freeze\"") != std::string::npos);
  assert(op_text.find("\"provenance_config_hash_hex\"") != std::string::npos);
  assert(op_text.find("\"status\": \"ok\"") != std::string::npos);
  assert(op_text.find("\"event_kind\": \"gravity.treepm_setup\"") != std::string::npos);
  assert(op_text.find("\"event_kind\": \"gravity.pm_long_range_field\"") != std::string::npos);
  assert(op_text.find("\"pm_grid\": \"9x9x9\"") != std::string::npos);
  assert(op_text.find("\"pm_assignment_scheme\": \"cic\"") != std::string::npos);
  assert(op_text.find("\"softening_kernel\": \"plummer\"") != std::string::npos);
  assert(op_text.find("\"pm_fft_backend\"") != std::string::npos);
  assert(op_text.find("\"pm_update_cadence_steps\": \"1\"") != std::string::npos);
  assert(op_text.find("\"event_kind\": \"gravity.health_summary\"") != std::string::npos);
  assert(op_text.find("\"heavy_reference_checks_opt_in\": \"false\"") != std::string::npos);

  std::stringstream tsc_stream;
  tsc_stream << buildConfigText(1, "reference_integration_test_tsc", "tsc");

  const cosmosim::core::FrozenConfig tsc_frozen =
      cosmosim::core::loadFrozenConfigFromString(tsc_stream.str(), "test_reference_workflow_tsc");
  cosmosim::workflows::ReferenceWorkflowRunner tsc_runner(tsc_frozen);
  const cosmosim::workflows::ReferenceWorkflowReport tsc_report =
      tsc_runner.run(output_dir, cosmosim::workflows::ReferenceWorkflowOptions{.write_outputs = false});
  assert(tsc_report.completed_steps == 2);
  assert(tsc_report.treepm_pm_grid == 9);
  assert(tsc_report.treepm_pm_grid_shape == "9x9x9");
  assert(tsc_report.treepm_long_range_refresh_count == 3);
  assert(tsc_report.treepm_long_range_reuse_count == 0);

  const std::string cadence_two_config =
      buildConfigText(2, "reference_integration_test_cadence_two", "cic");
  const cosmosim::core::FrozenConfig cadence_two_frozen =
      cosmosim::core::loadFrozenConfigFromString(cadence_two_config, "test_reference_workflow_cadence_two");
  cosmosim::workflows::ReferenceWorkflowRunner cadence_two_runner(cadence_two_frozen);
  const cosmosim::workflows::ReferenceWorkflowReport cadence_two_report_a =
      cadence_two_runner.run(output_dir, cosmosim::workflows::ReferenceWorkflowOptions{.write_outputs = false});
  const cosmosim::workflows::ReferenceWorkflowReport cadence_two_report_b =
      cadence_two_runner.run(output_dir, cosmosim::workflows::ReferenceWorkflowOptions{.write_outputs = false});

  assert(cadence_two_report_a.completed_steps == 2);
  assert(cadence_two_report_a.treepm_update_cadence_steps == 2);
  assert(cadence_two_report_a.treepm_long_range_refresh_count == 2);
  assert(cadence_two_report_a.treepm_long_range_reuse_count == 1);
  assert(cadence_two_report_a.treepm_cadence_records.size() == 3);
  const std::vector<std::uint64_t> expected_field_versions{1, 1, 2};
  const std::vector<std::uint64_t> expected_field_built_steps{0, 0, 1};
  const std::vector<bool> expected_refresh_flags{true, false, true};
  const std::vector<std::string> expected_stage_names{
      "gravity_kick_pre", "force_refresh", "force_refresh"};
  for (std::size_t i = 0; i < cadence_two_report_a.treepm_cadence_records.size(); ++i) {
    const auto& record = cadence_two_report_a.treepm_cadence_records[i];
    assert(record.stage_name == expected_stage_names[i]);
    assert(record.gravity_kick_opportunity == i + 1U);
    assert(record.field_version == expected_field_versions[i]);
    assert(record.field_built_step_index == expected_field_built_steps[i]);
    assert(record.refreshed_long_range_field == expected_refresh_flags[i]);
  }
  assert(cadence_two_report_a.final_state_digest == cadence_two_report_b.final_state_digest);
  assert(cadence_two_report_a.treepm_cadence_records.size() == cadence_two_report_b.treepm_cadence_records.size());
  for (std::size_t i = 0; i < cadence_two_report_a.treepm_cadence_records.size(); ++i) {
    const auto& lhs = cadence_two_report_a.treepm_cadence_records[i];
    const auto& rhs = cadence_two_report_b.treepm_cadence_records[i];
    assert(lhs.step_index == rhs.step_index);
    assert(lhs.stage_name == rhs.stage_name);
    assert(lhs.gravity_kick_opportunity == rhs.gravity_kick_opportunity);
    assert(lhs.field_version == rhs.field_version);
    assert(lhs.field_built_step_index == rhs.field_built_step_index);
    assert(lhs.refreshed_long_range_field == rhs.refreshed_long_range_field);
  }

  const std::string restart_diagnostics_config =
      buildConfigText(1, "reference_integration_test_restart_diagnostics", "cic") +
      "snapshot_interval_steps = 1\n"
      "write_restarts = true\n";
  const cosmosim::core::FrozenConfig restart_diagnostics_frozen =
      cosmosim::core::loadFrozenConfigFromString(
          restart_diagnostics_config,
          "test_reference_workflow_restart_diagnostics");
  cosmosim::workflows::ReferenceWorkflowRunner restart_diagnostics_runner(restart_diagnostics_frozen);
#if COSMOSIM_ENABLE_HDF5
  const cosmosim::workflows::ReferenceWorkflowReport restart_diagnostics_report =
      restart_diagnostics_runner.run(output_dir, cosmosim::workflows::ReferenceWorkflowOptions{.write_outputs = true});
  assert(restart_diagnostics_report.restart_roundtrip_executed);
  assert(restart_diagnostics_report.restart_roundtrip_ok);
  assert(std::filesystem::exists(restart_diagnostics_report.restart_path));
  const std::string restart_diag_events_text = readFile(restart_diagnostics_report.operational_report_json_path);
  assert(restart_diag_events_text.find("\"event_kind\": \"restart.write.complete\"") != std::string::npos);
  assert(restart_diag_events_text.find("\"event_kind\": \"restart.read.complete\"") != std::string::npos);
  assert(restart_diag_events_text.find("\"payload_hash_hex\"") != std::string::npos);
  assert(restart_diag_events_text.find("\"scheduler_pending_transition_count\"") != std::string::npos);
#else
  bool trapped_hdf5_disabled = false;
  try {
    (void)restart_diagnostics_runner.run(
        output_dir,
        cosmosim::workflows::ReferenceWorkflowOptions{.write_outputs = true});
  } catch (const std::runtime_error& error) {
    trapped_hdf5_disabled =
        std::string(error.what()).find("lacks HDF5 support") != std::string::npos;
  }
  assert(trapped_hdf5_disabled);
#endif

#if COSMOSIM_ENABLE_HDF5
  // This drives the real ReferenceWorkflow pipeline with a parentless cell,
  // two cells sharing one optional parent particle, and deliberately shuffled
  // dense rows.  The only comparison key is stable gas_cell_id.
  const std::string h2_identity_config =
      buildConfigText(1, "reference_h2_decoupled_gas_cells", "cic") +
      "snapshot_interval_steps = 1\n"
      "write_restarts = true\n";
  const cosmosim::core::FrozenConfig h2_identity_frozen =
      cosmosim::core::loadFrozenConfigFromString(
          h2_identity_config, "test_reference_workflow_h2_decoupled_gas_cells");

  const cosmosim::core::SimulationState canonical_h2_state =
      makeDecoupledGasCellState({0U, 1U, 2U});
  const cosmosim::core::SimulationState shuffled_h2_state =
      makeDecoupledGasCellState({2U, 0U, 1U});

  cosmosim::workflows::ReferenceWorkflowRunner h2_identity_runner(h2_identity_frozen);
  const auto canonical_h2_report = h2_identity_runner.run(
      output_dir,
      cosmosim::workflows::ReferenceWorkflowOptions{
          .write_outputs = true,
          .initial_state_override = &canonical_h2_state});

  std::string shuffled_h2_config = h2_identity_config;
  const std::string canonical_run_name = "reference_h2_decoupled_gas_cells";
  const std::string shuffled_run_name = "reference_h2_decoupled_gas_cells_shuffled";
  const std::size_t run_name_position = shuffled_h2_config.find(canonical_run_name);
  assert(run_name_position != std::string::npos);
  shuffled_h2_config.replace(run_name_position, canonical_run_name.size(), shuffled_run_name);
  const cosmosim::core::FrozenConfig shuffled_h2_frozen =
      cosmosim::core::loadFrozenConfigFromString(
          shuffled_h2_config, "test_reference_workflow_h2_decoupled_gas_cells_shuffled");
  cosmosim::workflows::ReferenceWorkflowRunner shuffled_h2_runner(shuffled_h2_frozen);
  const auto shuffled_h2_report = shuffled_h2_runner.run(
      output_dir,
      cosmosim::workflows::ReferenceWorkflowOptions{
          .write_outputs = true,
          .initial_state_override = &shuffled_h2_state});

  assert(canonical_h2_report.completed_steps == 2U);
  assert(shuffled_h2_report.completed_steps == 2U);
  assert(canonical_h2_report.restart_roundtrip_executed && canonical_h2_report.restart_roundtrip_ok);
  assert(shuffled_h2_report.restart_roundtrip_executed && shuffled_h2_report.restart_roundtrip_ok);

  const auto canonical_h2_restart = cosmosim::io::readRestartCheckpointHdf5(canonical_h2_report.restart_path);
  const auto shuffled_h2_restart = cosmosim::io::readRestartCheckpointHdf5(shuffled_h2_report.restart_path);
  const auto canonical_by_id = gasStateByStableId(canonical_h2_restart.state);
  const auto shuffled_by_id = gasStateByStableId(shuffled_h2_restart.state);
  assert(canonical_by_id.size() == 3U);
  assert(canonical_by_id == shuffled_by_id);
  const auto* parentless = canonical_h2_restart.state.gas_cell_identity.findByGasCellId(8001U);
  assert(parentless != nullptr && !parentless->parent_particle_id.has_value());
  assert(canonical_h2_restart.state.gas_cell_identity.rowsForParentParticleId(5001U).size() == 2U);
  assert(canonical_h2_restart.gas_cell_scheduler_ids.size() == 3U);
#endif

  auto invalidGravityConfig = []() {
    std::stringstream cfg;
    cfg << "schema_version = 1\n\n";
    cfg << "[mode]\n";
    cfg << "mode = zoom_in\n";
    cfg << "ic_file = generated\n";
    cfg << "zoom_high_res_region = false\n\n";
    cfg << "[numerics]\n";
    cfg << "time_begin_code = 0.01\n";
    cfg << "time_end_code = 0.0101\n";
    cfg << "max_global_steps = 1\n";
    cfg << "treepm_pm_grid = 9\n";
    cfg << "treepm_asmth_cells = nan\n";
    cfg << "treepm_rcut_cells = 6.0\n";
    cfg << "treepm_update_cadence_steps = 1\n\n";
    cfg << "[output]\n";
    cfg << "run_name = reference_integration_test_invalid_gravity\n";
    cfg << "output_directory = integration_outputs\n";
    cfg << "output_stem = snapshot\n";
    cfg << "restart_stem = restart\n";
    return cfg.str();
  };

  bool trapped_invalid_gravity_state = false;
  try {
    const cosmosim::core::FrozenConfig bad_frozen =
        cosmosim::core::loadFrozenConfigFromString(invalidGravityConfig(), "test_reference_workflow_invalid_gravity");
    cosmosim::workflows::ReferenceWorkflowRunner bad_runner(bad_frozen);
    (void)bad_runner.run(output_dir, cosmosim::workflows::ReferenceWorkflowOptions{.write_outputs = false});
  } catch (const cosmosim::core::ConfigError&) {
    trapped_invalid_gravity_state = true;
  } catch (const std::runtime_error& error) {
    const std::string message = error.what();
    trapped_invalid_gravity_state =
        message.find("fatal gravity-state check failed") != std::string::npos ||
        message.find("runtime workflow schema compatibility validation failed") != std::string::npos;
  }
  assert(trapped_invalid_gravity_state);


  const std::filesystem::path zoom_id_path = output_dir / "zoom_ids.txt";
  {
    std::ofstream zoom_out(zoom_id_path);
    zoom_out << "1\n2\n3\n";
  }
  std::stringstream zoom_stream;
  zoom_stream << "schema_version = 1\n\n";
  zoom_stream << "[mode]\n";
  zoom_stream << "mode = zoom_in\n";
  zoom_stream << "ic_file = generated\n";
  zoom_stream << "zoom_high_res_region = true\n";
  zoom_stream << "zoom_region_file = " << zoom_id_path.string() << "\n";
  zoom_stream << "zoom_region_radius = 0.05\n";
  zoom_stream << "zoom_long_range_strategy = global_coarse_plus_focused_highres_correction\n";
  zoom_stream << "zoom_focused_pm_grid_nx = 9\n";
  zoom_stream << "zoom_focused_pm_grid_ny = 9\n";
  zoom_stream << "zoom_focused_pm_grid_nz = 9\n\n";
  zoom_stream << "[numerics]\n";
  zoom_stream << "time_begin_code = 0.01\n";
  zoom_stream << "time_end_code = 0.0101\n";
  zoom_stream << "max_global_steps = 1\n";
  zoom_stream << "hierarchical_max_rung = 1\n";
  zoom_stream << "treepm_pm_grid = 9\n";
  zoom_stream << "treepm_asmth_cells = 1.25\n";
  zoom_stream << "treepm_rcut_cells = 4.5\n";
  zoom_stream << "treepm_assignment_scheme = cic\n";
  zoom_stream << "treepm_update_cadence_steps = 1\n\n";
  zoom_stream << "[output]\n";
  zoom_stream << "run_name = reference_integration_test_zoom_file\n";
  zoom_stream << "output_directory = integration_outputs\n";
  zoom_stream << "output_stem = snapshot\n";
  zoom_stream << "restart_stem = restart\n";

  const cosmosim::core::FrozenConfig zoom_frozen =
      cosmosim::core::loadFrozenConfigFromString(zoom_stream.str(), "test_reference_workflow_zoom_file");
  cosmosim::workflows::ReferenceWorkflowRunner zoom_runner(zoom_frozen);
  const cosmosim::workflows::ReferenceWorkflowReport zoom_report =
      zoom_runner.run(output_dir, cosmosim::workflows::ReferenceWorkflowOptions{.write_outputs = false});
  std::ifstream zoom_op_in(zoom_report.operational_report_json_path);
  const std::string zoom_op_text((std::istreambuf_iterator<char>(zoom_op_in)), std::istreambuf_iterator<char>());
  zoom_op_in.close();
  assert(zoom_op_text.find("\"zoom_membership_source\": \"particle_id_file_text\"") != std::string::npos);

  std::error_code cleanup_error;
  std::filesystem::remove_all(output_dir, cleanup_error);
  return 0;
}
