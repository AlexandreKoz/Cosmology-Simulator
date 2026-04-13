#include <cassert>
#include <stdexcept>

#include "cosmosim/core/provenance.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/core/time_integration.hpp"
#include "cosmosim/io/restart_checkpoint.hpp"

int main() {
  const auto& schema = cosmosim::io::restartSchema();
  assert(schema.name == "cosmosim_restart_v2");
  assert(schema.version == 2);
  assert(cosmosim::io::isRestartSchemaCompatible(2));
  assert(!cosmosim::io::isRestartSchemaCompatible(1));
  const auto& checklist = cosmosim::io::exactRestartCompletenessChecklist();
  assert(!checklist.empty());
  assert(checklist.front() == "simulation_state_lanes_and_metadata");
  assert(checklist.back() == "payload_integrity_hash_and_hex");

  bool missing_fields_threw = false;
  try {
    cosmosim::io::RestartWritePayload incomplete;
    (void)cosmosim::io::restartPayloadIntegrityHash(incomplete);
  } catch (const std::invalid_argument&) {
    missing_fields_threw = true;
  }
  assert(missing_fields_threw);

  cosmosim::core::SimulationState state;
  state.resizeParticles(1);
  state.resizeCells(0);
  state.resizePatches(0);
  state.particle_sidecar.species_tag[0] =
      static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kDarkMatter);
  state.species.count_by_species[0] = 1;
  state.rebuildSpeciesIndex();
  assert(state.validateOwnershipInvariants());

  cosmosim::core::IntegratorState integrator_state;
  cosmosim::core::HierarchicalTimeBinScheduler scheduler(1);
  scheduler.reset(1, 0, 0);

  cosmosim::io::RestartWritePayload payload;
  payload.state = &state;
  payload.integrator_state = &integrator_state;
  payload.scheduler = &scheduler;
  payload.normalized_config_text = "mode = toy\n";
  payload.normalized_config_hash_hex = "a";
  payload.provenance = cosmosim::core::makeProvenanceRecord("a", "deadbeef");

  const std::uint64_t hash_before = cosmosim::io::restartPayloadIntegrityHash(payload);
  integrator_state.time_bins.active_bin = 1;
  const std::uint64_t hash_after = cosmosim::io::restartPayloadIntegrityHash(payload);
  assert(hash_before != hash_after);

  bool missing_metadata_threw = false;
  try {
    cosmosim::io::RestartWritePayload invalid_metadata = payload;
    invalid_metadata.normalized_config_hash_hex.clear();
    (void)cosmosim::io::restartPayloadIntegrityHash(invalid_metadata);
  } catch (const std::invalid_argument&) {
    missing_metadata_threw = true;
  }
  assert(missing_metadata_threw);

  return 0;
}
