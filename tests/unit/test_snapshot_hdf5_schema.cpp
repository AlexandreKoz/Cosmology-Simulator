#include <cassert>

#include "cosmosim/io/io_contract.hpp"
#include "cosmosim/io/snapshot_hdf5.hpp"

namespace {

void testCanonicalSchemaNames() {
  const auto& schema = cosmosim::io::gadgetArepoSchemaMap();
  assert(schema.header_group == "/Header");
  assert(schema.config_group == "/Config");
  assert(schema.provenance_group == "/Provenance");
  assert(schema.config_normalized_attribute == "normalized");
  assert(schema.schema_name == "gadget_arepo_v2");
  assert(schema.schema_version == 2);
  assert(schema.part_type_group[0] == "/PartType0");
  assert(schema.coordinates.canonical_name == "Coordinates");
  assert(schema.velocities.canonical_name == "Velocities");
  assert(schema.masses.canonical_name == "Masses");
  assert(schema.particle_ids.canonical_name == "ParticleIDs");

  bool has_particle_id_alias = false;
  for (std::string_view alias : schema.particle_ids.read_aliases) {
    if (alias == "ParticleID") {
      has_particle_id_alias = true;
    }
  }
  assert(has_particle_id_alias);

  const auto& shared = cosmosim::io::sharedIoContractNames();
  assert(shared.normalized_config_text_dataset == "normalized_config_text");
  assert(shared.normalized_config_hash_hex_attribute == "normalized_config_hash_hex");
  assert(shared.provenance_record_dataset == "provenance_record");
}

}  // namespace

int main() {
  testCanonicalSchemaNames();
  return 0;
}
