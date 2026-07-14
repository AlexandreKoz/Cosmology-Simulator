#pragma once

#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <optional>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "cosmosim/cosmosim.hpp"
#include "cosmosim/io/restart_checkpoint.hpp"

namespace cosmosim::test::workflow_hydro_fixture {

constexpr double k_gamma = 5.0 / 3.0;
constexpr std::uint64_t k_patch_id = 424242U;
constexpr std::size_t k_nx = 2U;
constexpr std::size_t k_ny = 2U;
constexpr std::size_t k_nz = 2U;
constexpr std::size_t k_cell_count = k_nx * k_ny * k_nz;

[[nodiscard]] inline std::uint32_t linearRow(
    const std::uint32_t i,
    const std::uint32_t j,
    const std::uint32_t k) {
  return i + static_cast<std::uint32_t>(k_nx) *
      (j + static_cast<std::uint32_t>(k_ny) * k);
}

[[nodiscard]] inline std::string makeConfigText(const std::string& run_name) {
  std::ostringstream stream;
  stream << "schema_version = 1\n\n";
  stream << "[mode]\n";
  stream << "mode = cosmo_cube\n";
  stream << "ic_file = generated\n\n";
  stream << "[numerics]\n";
  stream << "time_begin_code = 0.01\n";
  stream << "time_end_code = 0.0104\n";
  stream << "max_global_steps = 2\n";
  stream << "hierarchical_max_rung = 0\n";
  stream << "treepm_pm_grid = 9\n";
  stream << "treepm_asmth_cells = 1.25\n";
  // Keep the low-cost CI mesh inside the periodic one-minimum-image cutoff.
  stream << "treepm_rcut_cells = 3.9\n";
  stream << "treepm_assignment_scheme = cic\n";
  stream << "treepm_update_cadence_steps = 1\n\n";
  stream << "[output]\n";
  stream << "run_name = " << run_name << '\n';
  stream << "output_directory = integration_outputs\n";
  stream << "output_stem = snapshot\n";
  stream << "restart_stem = restart\n";
  stream << "snapshot_interval_steps = 1\n";
  stream << "write_restarts = true\n";
  return stream.str();
}

[[nodiscard]] inline core::FrozenConfig makeFrozenConfig(const std::string& run_name) {
  return core::loadFrozenConfigFromString(
      makeConfigText(run_name),
      "reference_workflow_hydro_h1_fixture_" + run_name);
}

// new_dense_row -> physical Cartesian row. Patch geometry is deliberately left
// legacy/unmaterialized: production fixed-patch workflow hydro must construct it
// from the real centers, never from dense row order or cell_count factorization.
[[nodiscard]] inline core::SimulationState makeState(
    const std::array<std::uint32_t, k_cell_count>& dense_to_physical) {
  core::SimulationState state;
  state.resizeParticles(k_cell_count);
  state.resizeCells(k_cell_count);
  state.resizePatches(1);

  state.species.count_by_species.fill(0U);
  state.species.count_by_species[static_cast<std::size_t>(core::ParticleSpecies::kGas)] = k_cell_count;
  for (std::uint32_t physical = 0; physical < k_cell_count; ++physical) {
    const std::uint32_t i = physical % static_cast<std::uint32_t>(k_nx);
    const std::uint32_t j = (physical / static_cast<std::uint32_t>(k_nx)) % static_cast<std::uint32_t>(k_ny);
    const std::uint32_t k = physical / static_cast<std::uint32_t>(k_nx * k_ny);
    state.particle_sidecar.particle_id[physical] = 10001U + physical;
    state.particle_sidecar.sfc_key[physical] = 17U + physical * 29U;
    state.particle_sidecar.species_tag[physical] = static_cast<std::uint32_t>(core::ParticleSpecies::kGas);
    state.particle_sidecar.owning_rank[physical] = 0U;
    state.particles.position_x_comoving[physical] = 0.25 + 0.5 * static_cast<double>(i);
    state.particles.position_y_comoving[physical] = 0.25 + 0.5 * static_cast<double>(j);
    state.particles.position_z_comoving[physical] = 0.25 + 0.5 * static_cast<double>(k);
    state.particles.velocity_x_peculiar[physical] = 0.0;
    state.particles.velocity_y_peculiar[physical] = 0.0;
    state.particles.velocity_z_peculiar[physical] = 0.0;
    state.particles.mass_code[physical] = 1.0;
    state.particles.time_bin[physical] = 0U;
  }
  state.rebuildSpeciesIndex();

  state.patches.patch_id[0] = k_patch_id;
  state.patches.level[0] = 0;
  state.patches.first_cell[0] = 0U;
  state.patches.cell_count[0] = static_cast<std::uint32_t>(k_cell_count);
  state.patches.owning_rank[0] = 0U;
  // Zero dimensions/extents intentionally mark fixed-patch compatibility state.
  // The H1 workflow must derive strict Cartesian geometry from centers.

  std::vector<core::GasCellIdentityRecord> records;
  records.reserve(k_cell_count);
  constexpr double volume = 0.125;
  for (std::uint32_t dense_row = 0; dense_row < k_cell_count; ++dense_row) {
    const std::uint32_t physical = dense_to_physical[dense_row];
    assert(physical < k_cell_count);
    const std::uint32_t i = physical % static_cast<std::uint32_t>(k_nx);
    const std::uint32_t j = (physical / static_cast<std::uint32_t>(k_nx)) % static_cast<std::uint32_t>(k_ny);
    const std::uint32_t k = physical / static_cast<std::uint32_t>(k_nx * k_ny);
    const double rho = 0.9 + 0.07 * static_cast<double>(physical);
    const double pressure = 0.8 + 0.03 * static_cast<double>(i) +
        0.05 * static_cast<double>(j) + 0.07 * static_cast<double>(k);
    const double vx = i == 0U ? -0.13 : 0.19;
    const double vy = j == 0U ? 0.09 : -0.06;
    const double vz = k == 0U ? -0.04 : 0.08;

    state.cells.center_x_comoving[dense_row] = 0.25 + 0.5 * static_cast<double>(i);
    state.cells.center_y_comoving[dense_row] = 0.25 + 0.5 * static_cast<double>(j);
    state.cells.center_z_comoving[dense_row] = 0.25 + 0.5 * static_cast<double>(k);
    state.cells.mass_code[dense_row] = rho * volume;
    state.cells.time_bin[dense_row] = 0U;
    state.cells.patch_index[dense_row] = 0U;
    state.gas_cells.density_code[dense_row] = rho;
    state.gas_cells.pressure_code[dense_row] = pressure;
    state.gas_cells.internal_energy_code[dense_row] = pressure / ((k_gamma - 1.0) * rho);
    state.gas_cells.temperature_code[dense_row] = pressure / rho;
    state.gas_cells.sound_speed_code[dense_row] = std::sqrt(k_gamma * pressure / rho);
    state.gas_cells.velocity_x_peculiar[dense_row] = vx;
    state.gas_cells.velocity_y_peculiar[dense_row] = vy;
    state.gas_cells.velocity_z_peculiar[dense_row] = vz;

    const std::uint64_t cell_id = 70001U + physical;
    std::optional<std::uint64_t> parent;
    if (physical == 1U || physical == 2U) {
      parent = 10002U;  // Intentional shared parent.
    } else if (physical != 0U) {
      parent = 10001U + physical;
    }
    records.push_back(core::GasCellIdentityRecord{
        .gas_cell_id = cell_id,
        .parent_particle_id = parent,
        .owning_patch_id = k_patch_id,
        .local_cell_row = dense_row,
    });
  }
  state.replaceGasCellIdentityRecords(std::move(records));
  assert(state.gasCellIdentityMapMatchesSidecarLanes());
  assert(state.validateOwnershipInvariants());
  return state;
}

struct CellSnapshot {
  std::array<double, 11> values{};
  std::uint64_t patch_id = 0;
  std::uint32_t dense_row = 0;
};

[[nodiscard]] inline std::unordered_map<std::uint64_t, CellSnapshot> byStableId(
    const core::SimulationState& state) {
  state.requireGasCellIdentityMapCoversDenseRows("reference workflow hydro fixture");
  std::unordered_map<std::uint64_t, CellSnapshot> result;
  for (std::uint32_t row = 0; row < state.cells.size(); ++row) {
    const core::GasCellIdentityRecord* record = state.gas_cell_identity.findByLocalRow(row);
    assert(record != nullptr);
    result.emplace(record->gas_cell_id, CellSnapshot{
        .values = {state.cells.center_x_comoving[row], state.cells.center_y_comoving[row],
                   state.cells.center_z_comoving[row], state.cells.mass_code[row],
                   state.gas_cells.density_code[row], state.gas_cells.pressure_code[row],
                   state.gas_cells.internal_energy_code[row], state.gas_cells.temperature_code[row],
                   state.gas_cells.velocity_x_peculiar[row], state.gas_cells.velocity_y_peculiar[row],
                   state.gas_cells.velocity_z_peculiar[row]},
        .patch_id = record->owning_patch_id,
        .dense_row = row,
    });
  }
  return result;
}

inline void assertNear(const double lhs, const double rhs, const double tolerance = 2.0e-11) {
  assert(std::abs(lhs - rhs) <= tolerance * std::max({1.0, std::abs(lhs), std::abs(rhs)}));
}

inline void assertEquivalentByStableId(
    const core::SimulationState& lhs,
    const core::SimulationState& rhs) {
  const auto lhs_by_id = byStableId(lhs);
  const auto rhs_by_id = byStableId(rhs);
  assert(lhs_by_id.size() == rhs_by_id.size());
  for (const auto& [gas_cell_id, lhs_cell] : lhs_by_id) {
    const auto rhs_it = rhs_by_id.find(gas_cell_id);
    assert(rhs_it != rhs_by_id.end());
    assert(lhs_cell.patch_id == rhs_it->second.patch_id);
    for (std::size_t component = 0; component < lhs_cell.values.size(); ++component) {
      assertNear(lhs_cell.values[component], rhs_it->second.values[component]);
    }
    const auto lhs_row = lhs.rowForGasCellId(gas_cell_id);
    const auto rhs_row = rhs.rowForGasCellId(gas_cell_id);
    assert(lhs_row.has_value() && rhs_row.has_value());
    assert(lhs.cells.time_bin[*lhs_row] == rhs.cells.time_bin[*rhs_row]);
  }
}

// Reorder dense gas-cell storage at a safe workflow boundary.  The physical
// Cartesian centers and stable gas_cell_id ownership remain unchanged; the
// gas scheduler is remapped by gas_cell_id rather than dense local row.
inline void reorderGasRowsAndScheduler(
    io::RestartReadResult& restart,
    const std::array<std::uint32_t, k_cell_count>& new_row_to_old_row) {
  core::SimulationState& state = restart.state;
  assert(state.cells.size() == k_cell_count);
  state.requireGasCellIdentityMapCoversDenseRows("workflow hydro fixture reorder");

  core::HierarchicalTimeBinScheduler old_scheduler(restart.gas_cell_scheduler_state.max_bin);
  old_scheduler.importPersistentState(restart.gas_cell_scheduler_state);
  std::array<std::uint32_t, k_cell_count> all_rows{};
  for (std::uint32_t row = 0; row < k_cell_count; ++row) {
    all_rows[row] = row;
  }
  const std::vector<core::TimeBinSchedulerIdentityRecord> scheduler_records =
      core::exportGasCellSchedulerIdentityRecords(old_scheduler, state, all_rows);

  const auto old_cells = state.cells;
  const auto old_gas_cells = state.gas_cells;
  const core::GasCellIdentityMap old_identity = state.gas_cell_identity;
  for (std::uint32_t new_row = 0; new_row < k_cell_count; ++new_row) {
    const std::uint32_t old_row = new_row_to_old_row[new_row];
    assert(old_row < k_cell_count);
    state.cells.center_x_comoving[new_row] = old_cells.center_x_comoving[old_row];
    state.cells.center_y_comoving[new_row] = old_cells.center_y_comoving[old_row];
    state.cells.center_z_comoving[new_row] = old_cells.center_z_comoving[old_row];
    state.cells.mass_code[new_row] = old_cells.mass_code[old_row];
    state.cells.time_bin[new_row] = old_cells.time_bin[old_row];
    state.cells.patch_index[new_row] = old_cells.patch_index[old_row];
    state.gas_cells.gas_cell_id[new_row] = old_gas_cells.gas_cell_id[old_row];
    state.gas_cells.parent_particle_id[new_row] = old_gas_cells.parent_particle_id[old_row];
    state.gas_cells.velocity_x_peculiar[new_row] = old_gas_cells.velocity_x_peculiar[old_row];
    state.gas_cells.velocity_y_peculiar[new_row] = old_gas_cells.velocity_y_peculiar[old_row];
    state.gas_cells.velocity_z_peculiar[new_row] = old_gas_cells.velocity_z_peculiar[old_row];
    state.gas_cells.density_code[new_row] = old_gas_cells.density_code[old_row];
    state.gas_cells.pressure_code[new_row] = old_gas_cells.pressure_code[old_row];
    state.gas_cells.internal_energy_code[new_row] = old_gas_cells.internal_energy_code[old_row];
    state.gas_cells.temperature_code[new_row] = old_gas_cells.temperature_code[old_row];
    state.gas_cells.sound_speed_code[new_row] = old_gas_cells.sound_speed_code[old_row];
  }

  std::vector<core::GasCellIdentityRecord> records;
  records.reserve(k_cell_count);
  for (std::uint32_t new_row = 0; new_row < k_cell_count; ++new_row) {
    const core::GasCellIdentityRecord* old_record =
        old_identity.findByLocalRow(new_row_to_old_row[new_row]);
    assert(old_record != nullptr);
    core::GasCellIdentityRecord record = *old_record;
    record.local_cell_row = new_row;
    records.push_back(record);
  }
  state.replaceGasCellIdentityRecords(std::move(records));
  state.requireGasCellIdentityMapCoversDenseRows("workflow hydro fixture reordered state");

  core::HierarchicalTimeBinScheduler remapped_scheduler(restart.gas_cell_scheduler_state.max_bin);
  // Preserve the shared timeline before replacing dense-row scheduler lanes.
  remapped_scheduler.importPersistentState(restart.gas_cell_scheduler_state);
  core::rebuildSchedulerFromGasCellIdentityRecords(
      remapped_scheduler, scheduler_records, state);
  core::syncGasCellTimeBinMirrorsFromGasCellScheduler(remapped_scheduler, state);
  restart.gas_cell_scheduler_state = remapped_scheduler.exportPersistentState();
  restart.gas_cell_scheduler_ids.resize(state.cells.size());
  for (std::uint32_t row = 0; row < state.cells.size(); ++row) {
    restart.gas_cell_scheduler_ids[row] = *state.gas_cell_identity.gasCellIdForLocalRow(row);
  }
}

struct GlobalConservation {
  double mass = 0.0;
  double momentum_x = 0.0;
  double momentum_y = 0.0;
  double momentum_z = 0.0;
  double total_energy = 0.0;
};

[[nodiscard]] inline GlobalConservation globalConservation(const core::SimulationState& state) {
  GlobalConservation totals;
  for (std::uint32_t row = 0; row < state.cells.size(); ++row) {
    const double mass = state.cells.mass_code[row];
    const double vx = state.gas_cells.velocity_x_peculiar[row];
    const double vy = state.gas_cells.velocity_y_peculiar[row];
    const double vz = state.gas_cells.velocity_z_peculiar[row];
    totals.mass += mass;
    totals.momentum_x += mass * vx;
    totals.momentum_y += mass * vy;
    totals.momentum_z += mass * vz;
    totals.total_energy += mass * state.gas_cells.internal_energy_code[row] +
        0.5 * mass * (vx * vx + vy * vy + vz * vz);
  }
  return totals;
}

inline void assertConservationEquivalent(const GlobalConservation& lhs, const GlobalConservation& rhs) {
  assertNear(lhs.mass, rhs.mass);
  assertNear(lhs.momentum_x, rhs.momentum_x);
  assertNear(lhs.momentum_y, rhs.momentum_y);
  assertNear(lhs.momentum_z, rhs.momentum_z);
  assertNear(lhs.total_energy, rhs.total_energy);
}

inline void assertEquivalentCfl(
    const core::HydroCflDiagnostics& lhs,
    const core::HydroCflDiagnostics& rhs) {
  assert(lhs.has_gas_cell_id && rhs.has_gas_cell_id);
  assert(lhs.has_patch_id && rhs.has_patch_id);
  assert(lhs.has_patch_row && rhs.has_patch_row);
  assert(lhs.gas_cell_id == rhs.gas_cell_id);
  assert(lhs.patch_id == rhs.patch_id);
  assert(lhs.patch_row == rhs.patch_row);
  assert(lhs.limiting_axis == rhs.limiting_axis);
  assertNear(lhs.proposed_dt_time_code, rhs.proposed_dt_time_code);
  assertNear(lhs.accepted_dt_time_code, rhs.accepted_dt_time_code);
  for (std::size_t axis = 0; axis < 3; ++axis) {
    assertNear(lhs.cell_width_axis_code[axis], rhs.cell_width_axis_code[axis]);
  }
}

inline void assertEquivalentGasSchedulerById(
    const io::RestartReadResult& lhs,
    const io::RestartReadResult& rhs) {
  assert(lhs.gas_cell_scheduler_ids.size() == lhs.gas_cell_scheduler_state.bin_index.size());
  assert(rhs.gas_cell_scheduler_ids.size() == rhs.gas_cell_scheduler_state.bin_index.size());
  std::unordered_map<std::uint64_t, std::size_t> lhs_rows;
  for (std::size_t row = 0; row < lhs.gas_cell_scheduler_ids.size(); ++row) {
    assert(lhs_rows.emplace(lhs.gas_cell_scheduler_ids[row], row).second);
  }
  for (std::size_t rhs_row = 0; rhs_row < rhs.gas_cell_scheduler_ids.size(); ++rhs_row) {
    const auto lhs_it = lhs_rows.find(rhs.gas_cell_scheduler_ids[rhs_row]);
    assert(lhs_it != lhs_rows.end());
    const std::size_t lhs_row = lhs_it->second;
    assert(lhs.gas_cell_scheduler_state.bin_index[lhs_row] == rhs.gas_cell_scheduler_state.bin_index[rhs_row]);
    assert(lhs.gas_cell_scheduler_state.next_activation_tick[lhs_row] == rhs.gas_cell_scheduler_state.next_activation_tick[rhs_row]);
    assert(lhs.gas_cell_scheduler_state.active_flag[lhs_row] == rhs.gas_cell_scheduler_state.active_flag[rhs_row]);
    assert(lhs.gas_cell_scheduler_state.pending_bin_index[lhs_row] == rhs.gas_cell_scheduler_state.pending_bin_index[rhs_row]);
  }
}

inline void cleanup(const std::filesystem::path& root) {
  std::error_code error;
  std::filesystem::remove_all(root, error);
}

}  // namespace cosmosim::test::workflow_hydro_fixture
