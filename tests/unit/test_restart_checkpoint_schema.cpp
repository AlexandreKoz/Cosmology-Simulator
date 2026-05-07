#include <cassert>
#include <stdexcept>
#include <string_view>

#include "cosmosim/core/provenance.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/core/time_integration.hpp"
#include "cosmosim/io/restart_checkpoint.hpp"

int main() {
  const auto& schema = cosmosim::io::restartSchema();
  assert(schema.name == "cosmosim_restart_v6");
  assert(schema.version == 6);
  assert(cosmosim::io::isRestartSchemaCompatible(6));
  assert(!cosmosim::io::isRestartSchemaCompatible(5));
  const auto& checklist = cosmosim::io::exactRestartCompletenessChecklist();
  assert(!checklist.empty());
  assert(checklist.front() == "simulation_state_lanes_and_metadata");
  bool saw_softening = false;
  bool saw_gas_identity = false;
  bool saw_species_sidecars = false;
  for (const std::string_view item : checklist) {
    saw_softening = saw_softening || item == "particle_identity_and_softening_override_lanes";
    saw_gas_identity = saw_gas_identity || item == "gas_cell_identity_lanes";
    saw_species_sidecars = saw_species_sidecars || item == "species_specific_sidecars";
  }
  assert(saw_softening);
  assert(saw_gas_identity);
  assert(saw_species_sidecars);
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
  payload.normalized_config_hash_hex = cosmosim::core::stableConfigHashHex(payload.normalized_config_text);
  payload.provenance = cosmosim::core::makeProvenanceRecord(payload.normalized_config_hash_hex, "deadbeef");
  payload.distributed_gravity_state.schema_version = 2;
  payload.distributed_gravity_state.world_size = 1;
  payload.distributed_gravity_state.pm_grid_nx = 4;
  payload.distributed_gravity_state.pm_grid_ny = 4;
  payload.distributed_gravity_state.pm_grid_nz = 4;
  payload.distributed_gravity_state.owning_rank_by_item = {0};
  payload.distributed_gravity_state.pm_slab_begin_x_by_rank = {0};
  payload.distributed_gravity_state.pm_slab_end_x_by_rank = {4};

  const std::uint64_t hash_before = cosmosim::io::restartPayloadIntegrityHash(payload);
  integrator_state.time_bins.active_bin = 1;
  const std::uint64_t hash_after = cosmosim::io::restartPayloadIntegrityHash(payload);
  assert(hash_before != hash_after);

  integrator_state.time_bins.active_bin = 0;
  const std::uint64_t hash_restored = cosmosim::io::restartPayloadIntegrityHash(payload);
  state.particle_sidecar.setGravitySofteningOverride(0, 0.125);
  assert(cosmosim::io::restartPayloadIntegrityHash(payload) != hash_restored);
  const std::uint64_t hash_with_softening = cosmosim::io::restartPayloadIntegrityHash(payload);
  state.sidecars.upsert(cosmosim::core::ModuleSidecarBlock{
      .module_name = "hash_probe",
      .schema_version = 1,
      .payload = {std::byte{0xAA}}});
  assert(cosmosim::io::restartPayloadIntegrityHash(payload) != hash_with_softening);
  const std::uint64_t hash_with_sidecar = cosmosim::io::restartPayloadIntegrityHash(payload);
  payload.distributed_gravity_state.decomposition_epoch = 1;
  assert(cosmosim::io::restartPayloadIntegrityHash(payload) != hash_with_sidecar);
  payload.distributed_gravity_state.decomposition_epoch = 0;

  bool missing_metadata_threw = false;
  try {
    cosmosim::io::RestartWritePayload invalid_metadata = payload;
    invalid_metadata.normalized_config_hash_hex.clear();
    (void)cosmosim::io::restartPayloadIntegrityHash(invalid_metadata);
  } catch (const std::invalid_argument&) {
    missing_metadata_threw = true;
  }
  assert(missing_metadata_threw);

  bool mismatched_metadata_threw = false;
  try {
    cosmosim::io::RestartWritePayload invalid_metadata = payload;
    invalid_metadata.normalized_config_hash_hex = "0000000000000000";
    invalid_metadata.provenance =
        cosmosim::core::makeProvenanceRecord(invalid_metadata.normalized_config_hash_hex, "deadbeef");
    (void)cosmosim::io::restartPayloadIntegrityHash(invalid_metadata);
  } catch (const std::invalid_argument&) {
    mismatched_metadata_threw = true;
  }
  assert(mismatched_metadata_threw);

  bool text_hash_mismatch_threw = false;
  try {
    cosmosim::io::RestartWritePayload invalid_metadata = payload;
    invalid_metadata.normalized_config_text = "mode = altered\n";
    invalid_metadata.normalized_config_hash_hex = payload.normalized_config_hash_hex;
    invalid_metadata.provenance = payload.provenance;
    (void)cosmosim::io::restartPayloadIntegrityHash(invalid_metadata);
  } catch (const std::invalid_argument&) {
    text_hash_mismatch_threw = true;
  }
  assert(text_hash_mismatch_threw);

  return 0;
}
