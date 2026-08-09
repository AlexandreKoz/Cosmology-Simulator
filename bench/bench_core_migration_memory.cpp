#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include "cosmosim/core/memory_accounting.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "core_benchmark_process_memory.hpp"

namespace {

using cosmosim::core::GasCellIdentityRecord;
using cosmosim::core::GasCellMigrationCommit;
using cosmosim::core::GasCellMigrationRecord;
using cosmosim::core::SimulationState;

struct Scenario {
  std::string name;
  std::size_t cell_count;
  double outbound_fraction;
  double inbound_fraction;
};

Scenario selectScenario(const std::string& name, std::size_t cell_count) {
  if (name == "mostly_local") {
    return {name, cell_count, 0.02, 0.02};
  }
  if (name == "moderate") {
    return {name, cell_count, 0.25, 0.25};
  }
  if (name == "heavy") {
    return {name, cell_count, 0.75, 0.75};
  }
  throw std::invalid_argument("scenario must be mostly_local, moderate, or heavy");
}

void seedState(SimulationState& state, std::size_t count) {
  state.resizeCells(count);
  std::vector<GasCellIdentityRecord> identities;
  identities.reserve(count);
  for (std::size_t i = 0; i < count; ++i) {
    const std::uint32_t row = cosmosim::core::checkedLocalCellRow(i, "migration benchmark seed row");
    const std::uint64_t id = 1000000U + static_cast<std::uint64_t>(i);
    identities.push_back(GasCellIdentityRecord{
        .gas_cell_id = id,
        .parent_particle_id = std::nullopt,
        .owning_patch_id = 0U,
        .local_cell_row = row,
    });
    state.gas_cells.gas_cell_id[row] = id;
    state.gas_cells.parent_particle_id[row] = 0U;
    state.cells.center_x_comoving[row] = static_cast<double>(i % 1024U);
    state.cells.center_y_comoving[row] = static_cast<double>((i / 1024U) % 1024U);
    state.cells.center_z_comoving[row] = static_cast<double>(i / (1024U * 1024U));
    state.cells.mass_code[row] = 1.0;
    state.gas_cells.density_code[row] = 1.0;
    state.gas_cells.pressure_code[row] = 1.0;
    state.gas_cells.internal_energy_code[row] = 1.0;
    state.gas_cells.temperature_code[row] = 1.0;
    state.gas_cells.sound_speed_code[row] = 1.0;
  }
  state.gas_cell_identity.assign(std::move(identities));
}

GasCellMigrationRecord makeInbound(std::uint64_t id) {
  GasCellMigrationRecord record;
  record.owning_rank = 0U;
  record.fields.gas_cell_id = id;
  record.fields.has_parent_particle = 0U;
  record.fields.parent_particle_id = 0U;
  record.fields.owning_patch_id = 0U;
  record.fields.cell_mass_code = 1.0;
  record.fields.density_code = 1.0;
  record.fields.pressure_code = 1.0;
  record.fields.internal_energy_code = 1.0;
  record.fields.temperature_code = 1.0;
  record.fields.sound_speed_code = 1.0;
  return record;
}

std::uint64_t knownStateBytes(const SimulationState& state) {
  const auto report = cosmosim::core::collectSimulationMemoryReport(state, nullptr);
  return report.totals.persistent_total_bytes + report.totals.transient_total_bytes;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const std::size_t cell_count = argc > 2 ? static_cast<std::size_t>(std::stoull(argv[2])) : 50000U;
    if (cell_count == 0U) {
      throw std::invalid_argument("cell_count must be positive");
    }
    const Scenario scenario = selectScenario(argc > 1 ? argv[1] : "moderate", cell_count);
    const auto process_baseline = cosmosim::bench::sampleProcessMemory();

    SimulationState state;
    seedState(state, scenario.cell_count);
    const std::uint64_t steady_state_report_bytes = knownStateBytes(state);
    const auto steady_memory = cosmosim::bench::sampleProcessMemory();

    const std::size_t outbound_count = static_cast<std::size_t>(
        static_cast<double>(scenario.cell_count) * scenario.outbound_fraction);
    const std::size_t inbound_count = static_cast<std::size_t>(
        static_cast<double>(scenario.cell_count) * scenario.inbound_fraction);

    GasCellMigrationCommit commit;
    commit.world_rank = 0;
    commit.outbound_local_cell_rows.reserve(outbound_count);
    for (std::size_t i = 0; i < outbound_count; ++i) {
      commit.outbound_local_cell_rows.push_back(
          cosmosim::core::checkedLocalCellRow(i, "migration benchmark outbound row"));
    }
    commit.inbound_records.reserve(inbound_count);
    for (std::size_t i = 0; i < inbound_count; ++i) {
      commit.inbound_records.push_back(makeInbound(9000000000ULL + static_cast<std::uint64_t>(i)));
    }

    const auto precommit_memory = cosmosim::bench::sampleProcessMemory();
    const std::uint64_t exchange_payload_bytes =
        static_cast<std::uint64_t>(commit.inbound_records.capacity()) * sizeof(GasCellMigrationRecord) +
        static_cast<std::uint64_t>(commit.outbound_local_cell_rows.capacity()) * sizeof(std::uint32_t);

    const auto commit_begin = std::chrono::steady_clock::now();
    state.commitGasCellMigration(commit);
    const auto commit_end = std::chrono::steady_clock::now();
    const double commit_ms =
        std::chrono::duration<double, std::milli>(commit_end - commit_begin).count();
    const auto postcommit_memory = cosmosim::bench::sampleProcessMemory();
    const std::uint64_t final_report_bytes = knownStateBytes(state);

    std::cout << "benchmark=core_migration_memory"
              << " scenario=" << scenario.name
              << " initial_cells=" << scenario.cell_count
              << " outbound_cells=" << outbound_count
              << " inbound_cells=" << inbound_count
              << " final_cells=" << state.cells.size()
              << " steady_state_report_bytes=" << steady_state_report_bytes
              << " exchange_payload_capacity_bytes=" << exchange_payload_bytes
              << " final_report_bytes=" << final_report_bytes
              << " commit_ms=" << commit_ms;

    if (process_baseline.current_rss_bytes) {
      std::cout << " process_baseline_rss_bytes=" << *process_baseline.current_rss_bytes;
    }
    if (steady_memory.current_rss_bytes) {
      std::cout << " steady_rss_bytes=" << *steady_memory.current_rss_bytes;
    }
    if (precommit_memory.current_rss_bytes) {
      std::cout << " precommit_rss_bytes=" << *precommit_memory.current_rss_bytes;
    }
    if (postcommit_memory.current_rss_bytes) {
      std::cout << " postcommit_rss_bytes=" << *postcommit_memory.current_rss_bytes;
    }
    if (postcommit_memory.peak_rss_bytes) {
      std::cout << " peak_rss_bytes=" << *postcommit_memory.peak_rss_bytes;
      if (steady_memory.current_rss_bytes && *steady_memory.current_rss_bytes > 0U) {
        std::cout << " peak_to_steady_rss_ratio="
                  << static_cast<double>(*postcommit_memory.peak_rss_bytes) /
                         static_cast<double>(*steady_memory.current_rss_bytes);
      }
    } else {
      std::cout << " peak_rss=unavailable";
    }
    std::cout << '\n';
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "bench_core_migration_memory failed: " << error.what() << '\n';
    return 1;
  }
}
