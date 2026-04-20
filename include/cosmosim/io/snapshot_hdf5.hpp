#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

#include "cosmosim/core/config.hpp"
#include "cosmosim/core/provenance.hpp"
#include "cosmosim/core/simulation_state.hpp"

namespace cosmosim::io {

struct GadgetArepoFieldAliases {
  std::string_view canonical_name;
  std::vector<std::string_view> read_aliases;
};

struct GadgetArepoSchemaMap {
  std::string_view schema_name = "gadget_arepo_v3";
  std::uint32_t schema_version = 3;

  std::string_view header_group = "/Header";
  std::string_view config_group = "/Config";
  std::string_view provenance_group = "/Provenance";
  std::string_view config_normalized_attribute = "normalized";

  std::array<std::string_view, 6> part_type_group = {
      "/PartType0", "/PartType1", "/PartType2", "/PartType3", "/PartType4", "/PartType5"};

  GadgetArepoFieldAliases coordinates{"Coordinates", {"Coordinates", "Position", "POS"}};
  GadgetArepoFieldAliases velocities{"Velocities", {"Velocities", "Velocity", "VEL"}};
  GadgetArepoFieldAliases masses{"Masses", {"Masses", "Mass"}};
  GadgetArepoFieldAliases particle_ids{"ParticleIDs", {"ParticleIDs", "ParticleID", "ID"}};
};

[[nodiscard]] const GadgetArepoSchemaMap& gadgetArepoSchemaMap();

struct SnapshotIoPolicy {
  bool enable_compression = false;
  int compression_level = 4;
  std::size_t chunk_particle_count = 1u << 15;
  bool write_particle_type_alias_groups = true;
};

struct SnapshotWritePayload {
  const core::SimulationState* state = nullptr;
  const core::SimulationConfig* config = nullptr;
  std::string normalized_config_text;
  core::ProvenanceRecord provenance;
  std::string git_sha = "unknown";
};

struct SnapshotReadOptions {
  bool require_ids = true;
  bool require_velocities = true;
  bool allow_mass_table_fallback = true;
};

struct SnapshotIoReport {
  std::string schema_name = "gadget_arepo_v3";
  std::uint32_t schema_version = 3;
  bool compression_enabled = false;
  int compression_level = 0;
  std::size_t chunk_particle_count = 0;
  std::vector<std::string> present_aliases;
  std::vector<std::string> defaulted_fields;
};

struct SnapshotReadResult {
  core::SimulationState state;
  core::ProvenanceRecord provenance;
  std::string normalized_config_text;
  SnapshotIoReport report;
};

void writeGadgetArepoSnapshotHdf5(
    const std::filesystem::path& output_path,
    const SnapshotWritePayload& payload,
    const SnapshotIoPolicy& policy = {});

[[nodiscard]] SnapshotReadResult readGadgetArepoSnapshotHdf5(
    const std::filesystem::path& input_path,
    const core::SimulationConfig& config,
    const SnapshotReadOptions& options = {});

}  // namespace cosmosim::io
