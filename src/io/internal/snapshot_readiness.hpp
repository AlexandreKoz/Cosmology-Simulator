#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <unordered_map>

#include "cosmosim/io/snapshot_hdf5.hpp"

namespace cosmosim::io::internal {

inline void updateSnapshotReadiness(
    const core::SimulationState& state,
    SnapshotIoReport* report) {
  if (report == nullptr) {
    return;
  }
  report->analysis_ready = true;
  report->evolution_readiness_reasons.clear();

  const auto reject = [&](std::string reason) {
    report->evolution_readiness_reasons.push_back(std::move(reason));
  };

  if (!report->unavailable_fields.empty()) {
    reject("one or more scientific fields are explicitly unavailable");
  }
  if (!state.validatePersistentParticleIds()) {
    reject("persistent particle IDs are not all nonzero and unique");
  }
  if (!state.validateOwnershipInvariants()) {
    reject("SimulationState ownership/sidecar invariants are not satisfied");
  }

  std::unordered_map<std::uint64_t, std::size_t> particle_row_by_id;
  particle_row_by_id.reserve(state.particle_sidecar.particle_id.size());
  for (std::size_t row = 0; row < state.particle_sidecar.particle_id.size(); ++row) {
    particle_row_by_id.emplace(state.particle_sidecar.particle_id[row], row);
  }
  std::unordered_map<std::uint64_t, std::size_t> patch_row_by_id;
  patch_row_by_id.reserve(state.patches.patch_id.size());
  for (std::size_t row = 0; row < state.patches.patch_id.size(); ++row) {
    if (state.patches.patch_id[row] != 0U) {
      patch_row_by_id.emplace(state.patches.patch_id[row], row);
    }
  }

  bool gas_owner_failure = false;
  for (const core::GasCellIdentityRecord& record : state.gas_cell_identity.records()) {
    if (record.gas_cell_id == 0U || record.local_cell_row >= state.cells.size()) {
      gas_owner_failure = true;
      break;
    }

    bool owner_resolves = false;
    if (record.owning_patch_id != 0U) {
      const auto patch_it = patch_row_by_id.find(record.owning_patch_id);
      if (patch_it != patch_row_by_id.end()) {
        const std::size_t patch_row = patch_it->second;
        owner_resolves = state.cells.patch_index[record.local_cell_row] == patch_row;
      }
    }
    if (!owner_resolves && record.parent_particle_id.has_value()) {
      owner_resolves =
          *record.parent_particle_id != 0U &&
          particle_row_by_id.contains(*record.parent_particle_id);
    }
    if (!owner_resolves) {
      gas_owner_failure = true;
      break;
    }
  }
  if (gas_owner_failure) {
    reject(
        "one or more gas cells have neither a resolvable authoritative patch "
        "owner nor a resolvable local parent particle");
  }

  report->evolution_ready = report->evolution_readiness_reasons.empty();
}

}  // namespace cosmosim::io::internal
