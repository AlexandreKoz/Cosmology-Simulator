#include "workflows/internal/output_verification.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <string>
#include <unordered_map>
#include <utility>

#include "cosmosim/core/config.hpp"
#include "cosmosim/core/provenance.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/io/io_contract.hpp"
#include "cosmosim/io/snapshot_hdf5.hpp"

namespace cosmosim::workflows::internal {

[[nodiscard]] bool snapshotScalarEquivalent(double lhs, double rhs) {
  if (!std::isfinite(lhs) || !std::isfinite(rhs)) {
    return false;
  }
  const double scale = std::max({1.0, std::abs(lhs), std::abs(rhs)});
  return std::abs(lhs - rhs) <= 32.0 * std::numeric_limits<double>::epsilon() * scale;
}

[[nodiscard]] SnapshotRoundtripVerification verifySnapshotRoundtrip(
    const io::SnapshotReadResult& restored,
    const core::SimulationState& expected,
    const core::SimulationConfig& config,
    std::string_view expected_normalized_config,
    const core::ProvenanceRecord& expected_provenance) {
  const auto fail = [](std::string detail) {
    return SnapshotRoundtripVerification{.ok = false, .detail = std::move(detail)};
  };
  if (restored.report.schema_name != io::gadgetArepoSchemaMap().schema_name ||
      restored.report.schema_version != io::gadgetArepoSchemaMap().schema_version) {
    return fail("snapshot schema identity mismatch");
  }
  if (restored.report.file_kind != io::sharedIoContractNames().science_snapshot_file_kind ||
      restored.report.restart_compatible) {
    return fail("snapshot file-kind contract mismatch");
  }
  if (restored.normalized_config_text != expected_normalized_config ||
      restored.provenance.config_hash_hex != expected_provenance.config_hash_hex ||
      restored.provenance.normalized_config_hash_hex !=
          expected_provenance.normalized_config_hash_hex ||
      restored.provenance.schema_version != expected_provenance.schema_version) {
    return fail("snapshot config/provenance identity mismatch");
  }
  if (!snapshotScalarEquivalent(restored.report.header_time, expected.metadata.scale_factor) ||
      !snapshotScalarEquivalent(
          restored.report.header_redshift,
          expected.metadata.scale_factor > 0.0
              ? 1.0 / expected.metadata.scale_factor - 1.0
              : 0.0) ||
      !snapshotScalarEquivalent(
          restored.report.header_box_size_x,
          config.cosmology.box_size_x_mpc_comoving) ||
      !snapshotScalarEquivalent(
          restored.report.header_box_size_y,
          config.cosmology.box_size_y_mpc_comoving) ||
      !snapshotScalarEquivalent(
          restored.report.header_box_size_z,
          config.cosmology.box_size_z_mpc_comoving) ||
      !snapshotScalarEquivalent(
          restored.report.header_omega_matter,
          config.cosmology.omega_matter) ||
      !snapshotScalarEquivalent(
          restored.report.header_omega_lambda,
          config.cosmology.omega_lambda) ||
      !snapshotScalarEquivalent(
          restored.report.header_omega_baryon,
          config.cosmology.omega_baryon) ||
      !snapshotScalarEquivalent(
          restored.report.header_hubble_param,
          config.cosmology.hubble_param)) {
    return fail("snapshot Header cosmology/time metadata mismatch");
  }
  if (restored.state.particles.size() != expected.particles.size() ||
      !restored.state.validateUniqueParticleIds() ||
      !expected.validateUniqueParticleIds()) {
    return fail("snapshot particle count or stable-ID uniqueness mismatch");
  }

  std::unordered_map<std::uint64_t, std::size_t> restored_index_by_id;
  restored_index_by_id.reserve(restored.state.particles.size());
  for (std::size_t i = 0; i < restored.state.particles.size(); ++i) {
    restored_index_by_id.emplace(restored.state.particle_sidecar.particle_id[i], i);
  }
  const bool expected_has_softening =
      expected.particle_sidecar.gravity_softening_comoving.size() == expected.particles.size();
  const bool expected_has_softening_mask =
      expected.particle_sidecar.has_gravity_softening_override.size() == expected.particles.size();
  if (expected_has_softening !=
          (restored.state.particle_sidecar.gravity_softening_comoving.size() ==
           restored.state.particles.size()) ||
      expected_has_softening_mask !=
          (restored.state.particle_sidecar.has_gravity_softening_override.size() ==
           restored.state.particles.size())) {
    return fail("snapshot gravity-softening sidecar presence mismatch");
  }
  for (std::size_t expected_i = 0; expected_i < expected.particles.size(); ++expected_i) {
    const std::uint64_t particle_id = expected.particle_sidecar.particle_id[expected_i];
    const auto restored_it = restored_index_by_id.find(particle_id);
    if (restored_it == restored_index_by_id.end()) {
      return fail("snapshot is missing particle_id=" + std::to_string(particle_id));
    }
    const std::size_t restored_i = restored_it->second;
    if (restored.state.particle_sidecar.species_tag[restored_i] !=
        expected.particle_sidecar.species_tag[expected_i]) {
      return fail("snapshot species mismatch for particle_id=" + std::to_string(particle_id));
    }
    const std::array expected_values{
        expected.particles.position_x_comoving[expected_i],
        expected.particles.position_y_comoving[expected_i],
        expected.particles.position_z_comoving[expected_i],
        expected.particles.velocity_x_peculiar[expected_i],
        expected.particles.velocity_y_peculiar[expected_i],
        expected.particles.velocity_z_peculiar[expected_i],
        expected.particles.mass_code[expected_i]};
    const std::array restored_values{
        restored.state.particles.position_x_comoving[restored_i],
        restored.state.particles.position_y_comoving[restored_i],
        restored.state.particles.position_z_comoving[restored_i],
        restored.state.particles.velocity_x_peculiar[restored_i],
        restored.state.particles.velocity_y_peculiar[restored_i],
        restored.state.particles.velocity_z_peculiar[restored_i],
        restored.state.particles.mass_code[restored_i]};
    for (std::size_t component = 0; component < expected_values.size(); ++component) {
      if (!snapshotScalarEquivalent(expected_values[component], restored_values[component])) {
        return fail(
            "snapshot phase-space/mass mismatch for particle_id=" +
            std::to_string(particle_id));
      }
    }
    if (expected_has_softening && !snapshotScalarEquivalent(
            expected.particle_sidecar.gravity_softening_comoving[expected_i],
            restored.state.particle_sidecar.gravity_softening_comoving[restored_i])) {
      return fail("snapshot softening mismatch for particle_id=" + std::to_string(particle_id));
    }
    if (expected_has_softening_mask &&
        expected.particle_sidecar.has_gravity_softening_override[expected_i] !=
            restored.state.particle_sidecar.has_gravity_softening_override[restored_i]) {
      return fail("snapshot softening mask mismatch for particle_id=" + std::to_string(particle_id));
    }
  }

  if (restored.state.tracers.size() != expected.tracers.size()) {
    return fail("snapshot tracer count mismatch");
  }
  std::unordered_map<std::uint64_t, std::size_t> restored_tracer_row_by_particle_id;
  for (std::size_t row = 0; row < restored.state.tracers.size(); ++row) {
    const std::uint32_t particle_index = restored.state.tracers.particle_index[row];
    if (particle_index >= restored.state.particles.size()) {
      return fail("snapshot tracer particle index is out of range");
    }
    restored_tracer_row_by_particle_id.emplace(
        restored.state.particle_sidecar.particle_id[particle_index], row);
  }
  for (std::size_t expected_row = 0; expected_row < expected.tracers.size(); ++expected_row) {
    const std::uint32_t particle_index = expected.tracers.particle_index[expected_row];
    if (particle_index >= expected.particles.size()) {
      return fail("source tracer particle index is out of range");
    }
    const std::uint64_t particle_id = expected.particle_sidecar.particle_id[particle_index];
    const auto restored_row_it = restored_tracer_row_by_particle_id.find(particle_id);
    if (restored_row_it == restored_tracer_row_by_particle_id.end()) {
      return fail("snapshot tracer identity mismatch");
    }
    const std::size_t restored_row = restored_row_it->second;
    if (restored.state.tracers.parent_particle_id[restored_row] !=
            expected.tracers.parent_particle_id[expected_row] ||
        restored.state.tracers.injection_step[restored_row] !=
            expected.tracers.injection_step[expected_row] ||
        restored.state.tracers.host_cell_index[restored_row] !=
            expected.tracers.host_cell_index[expected_row] ||
        !snapshotScalarEquivalent(
            restored.state.tracers.mass_fraction_of_host[restored_row],
            expected.tracers.mass_fraction_of_host[expected_row]) ||
        !snapshotScalarEquivalent(
            restored.state.tracers.cumulative_exchanged_mass_code[restored_row],
            expected.tracers.cumulative_exchanged_mass_code[expected_row])) {
      return fail("snapshot tracer fields mismatch for particle_id=" + std::to_string(particle_id));
    }
  }
  return SnapshotRoundtripVerification{.ok = true, .detail = "scientific fields and metadata match"};
}

}  // namespace cosmosim::workflows::internal
