#include <cassert>
#include <stdexcept>
#include <string_view>
#include <type_traits>

#include "cosmosim/core/provenance.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/core/time_integration.hpp"
#include "cosmosim/io/restart_checkpoint.hpp"

template <typename T>
concept ExposesParticleSoaStorageMember = requires(T payload) {
  payload.particle_soa_storage;
};

template <typename T>
concept ExposesTransientWorkspaceMember = requires(T payload) {
  payload.transient_workspace;
};

int main() {
  static_assert(
      std::is_same_v<decltype(cosmosim::io::RestartWritePayload{}.persistent_state.simulation_state), const cosmosim::core::SimulationState*>,
      "Restart payload must serialize canonical SimulationState ownership root");
  static_assert(
      !ExposesParticleSoaStorageMember<cosmosim::io::RestartWritePayload>,
      "Restart payload must not expose ParticleSoaStorage as persistent restart truth");
  static_assert(
      !ExposesTransientWorkspaceMember<cosmosim::io::RestartWritePayload>,
      "Restart payload must not expose transient workspace state");

  const auto& schema = cosmosim::io::restartSchema();
  assert(schema.name == "cosmosim_restart_v14");
  assert(schema.version == 14);
  assert(cosmosim::io::isRestartSchemaCompatible(14));
  assert(!cosmosim::io::isRestartSchemaCompatible(13));
  assert(!cosmosim::io::isRestartSchemaCompatible(11));
  const auto& checklist = cosmosim::io::exactRestartCompletenessChecklist();
  assert(!checklist.empty());
  assert(checklist.front() == "simulation_state_lanes_and_metadata");
  bool saw_softening = false;
  bool saw_gas_identity = false;
  bool saw_species_sidecars = false;
  bool saw_output_cadence = false;
  bool saw_stochastic_state = false;
  bool saw_restart_diagnostics = false;
  for (const std::string_view item : checklist) {
    saw_softening = saw_softening || item == "particle_identity_softening_and_drift_epoch_lanes";
    saw_gas_identity = saw_gas_identity || item == "gas_cell_identity_lanes";
    saw_species_sidecars = saw_species_sidecars || item == "species_specific_sidecars";
    saw_output_cadence = saw_output_cadence || item == "output_cadence_persistent_state";
    saw_stochastic_state = saw_stochastic_state || item == "stochastic_module_persistent_state";
    saw_restart_diagnostics = saw_restart_diagnostics || item == "restart_diagnostics_summary";
  }
  assert(saw_softening);
  assert(saw_gas_identity);
  assert(saw_species_sidecars);
  assert(saw_output_cadence);
  assert(saw_stochastic_state);
  assert(saw_restart_diagnostics);
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
  payload.persistent_state.simulation_state = &state;
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
  payload.output_cadence_state.output_enabled = true;
  payload.output_cadence_state.write_restarts = true;
  payload.output_cadence_state.last_completed_step_index = integrator_state.step_index;
  payload.output_cadence_state.snapshot_interval_steps = 4;
  payload.output_cadence_state.next_snapshot_step_index = 4;
  payload.output_cadence_state.snapshot_stem = "snapshot";
  payload.output_cadence_state.restart_stem = "restart";
  assert(cosmosim::io::restartPayloadIntegrityHash(payload) != hash_restored);
  payload.output_cadence_state = {};
  payload.stochastic_state.modules.push_back(cosmosim::io::StochasticModulePersistentState{
      .module_name = "star_formation",
      .schema_version = 1,
      .rng_policy = "stateless_splitmix64(seed,step_index,cell_index,rank_local_seed_offset)",
      .random_seed = 123456789ull,
      .rank_local_seed_offset = 2,
      .last_committed_step_index = integrator_state.step_index,
      .deterministic_from_serialized_inputs = true,
  });
  assert(cosmosim::io::restartPayloadIntegrityHash(payload) != hash_restored);
  payload.stochastic_state = {};
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

  bool unsafe_half_step_payload_threw = false;
  try {
    cosmosim::core::IntegratorState half_step_integrator_state = integrator_state;
    half_step_integrator_state.inside_kdk_step = true;
    half_step_integrator_state.current_boundary_kind = cosmosim::core::StepBoundaryKind::kCheckpointPoint;
    half_step_integrator_state.last_completed_boundary_kind = cosmosim::core::StepBoundaryKind::kCheckpointPoint;
    half_step_integrator_state.last_completed_restart_safe = true;
    cosmosim::io::RestartWritePayload unsafe_payload = payload;
    unsafe_payload.integrator_state = &half_step_integrator_state;
    (void)cosmosim::io::restartPayloadIntegrityHash(unsafe_payload);
  } catch (const std::runtime_error& ex) {
    const std::string message = ex.what();
    unsafe_half_step_payload_threw =
        message.find("inside_kdk_step=true") != std::string::npos &&
        message.find("half-step restart is not represented") != std::string::npos;
  }
  assert(unsafe_half_step_payload_threw);

  bool missing_metadata_threw = false;
  try {
    cosmosim::io::RestartWritePayload invalid_metadata = payload;
    invalid_metadata.normalized_config_hash_hex.clear();
    (void)cosmosim::io::restartPayloadIntegrityHash(invalid_metadata);
  } catch (const std::invalid_argument&) {
    missing_metadata_threw = true;
  }
  assert(missing_metadata_threw);



  bool provenance_hash_message_threw = false;
  try {
    cosmosim::io::RestartWritePayload invalid_metadata = payload;
    invalid_metadata.provenance.normalized_config_hash_hex.clear();
    invalid_metadata.provenance.config_hash_hex.clear();
    (void)cosmosim::io::restartPayloadIntegrityHash(invalid_metadata);
  } catch (const std::invalid_argument& ex) {
    const std::string message = ex.what();
    provenance_hash_message_threw =
        message.find("provenance.normalized_config_hash_hex or provenance.config_hash_hex") != std::string::npos;
  }
  assert(provenance_hash_message_threw);

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

  {
    cosmosim::core::SimulationState gas_state;
    gas_state.resizeParticles(4);
    gas_state.resizeCells(2);
    gas_state.resizePatches(0);
    gas_state.species.count_by_species = {2, 2, 0, 0, 0};
    for (std::size_t i = 0; i < gas_state.particles.size(); ++i) {
      gas_state.particle_sidecar.particle_id[i] = 900 + i;
      gas_state.particle_sidecar.owning_rank[i] = 0;
      gas_state.particle_sidecar.species_tag[i] =
          i < 2 ? static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kDarkMatter)
                : static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kGas);
    }
    gas_state.rebuildSpeciesIndex();
    gas_state.refreshGasCellIdentityFromParticleOrder();

    cosmosim::core::HierarchicalTimeBinScheduler gas_scheduler(2);
    gas_scheduler.reset(4, 0, 0);
    gas_scheduler.setElementBin(2, 1, 0);
    gas_scheduler.setElementBin(3, 2, 0);
    cosmosim::core::syncTimeBinMirrorsFromScheduler(
        gas_scheduler,
        gas_state,
        cosmosim::core::TimeBinMirrorDomain::kParticlesAndCells);

    cosmosim::io::RestartWritePayload gas_payload = payload;
    gas_payload.persistent_state.simulation_state = &gas_state;
    gas_payload.scheduler = &gas_scheduler;
    gas_payload.distributed_gravity_state.owning_rank_by_item = {0, 0, 0, 0};
    assert(cosmosim::io::restartPayloadIntegrityHash(gas_payload) != 0);

    gas_state.cells.time_bin[0] = 0;
    bool stale_cell_mirror_threw = false;
    try {
      (void)cosmosim::io::restartPayloadIntegrityHash(gas_payload);
    } catch (const std::invalid_argument&) {
      stale_cell_mirror_threw = true;
    }
    assert(stale_cell_mirror_threw);
  }

  return 0;
}
