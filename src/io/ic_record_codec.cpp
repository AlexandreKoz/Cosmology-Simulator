#include "io/internal/ic_record_codec.hpp"

#include <stdexcept>

#include "io/internal/ic_byte_codec.hpp"

namespace cosmosim::io::internal {
void serializeIcRecord(
    const IcParticleRecord& record,
    std::vector<std::uint8_t>& output) {
  const std::size_t begin = output.size();
  appendLe64(output, record.id);
  appendLe32(output, record.species);
  appendDouble(output, record.x);
  appendDouble(output, record.y);
  appendDouble(output, record.z);
  appendDouble(output, record.vx);
  appendDouble(output, record.vy);
  appendDouble(output, record.vz);
  appendDouble(output, record.mass);
  appendDouble(output, record.gas_density);
  appendDouble(output, record.gas_internal_energy);
  appendDouble(output, record.gas_metallicity);
  appendDouble(output, record.star_formation);
  appendDouble(output, record.star_birth_mass);
  appendDouble(output, record.star_metallicity);
  appendDouble(output, record.bh_mass);
  appendDouble(output, record.bh_mdot);
  appendLe64(output, record.tracer_parent);
  appendLe64(output, record.tracer_injection);
  appendLe32(output, record.tracer_host);
  appendDouble(output, record.tracer_fraction);
  appendDouble(output, record.tracer_last_host_mass);
  appendDouble(output, record.tracer_exchanged_mass);
  if (output.size() - begin != kIcWireRecordBytes) {
    throw std::logic_error("IC wire record byte contract drifted");
  }
}

IcParticleRecord deserializeIcRecord(
    std::span<const std::uint8_t> bytes,
    std::size_t& offset) {
  if (offset + kIcWireRecordBytes > bytes.size()) {
    throw std::runtime_error("truncated IC wire record");
  }
  IcParticleRecord record;
  record.id = readLe64(bytes, offset);
  record.species = readLe32(bytes, offset);
  record.x = readDouble(bytes, offset);
  record.y = readDouble(bytes, offset);
  record.z = readDouble(bytes, offset);
  record.vx = readDouble(bytes, offset);
  record.vy = readDouble(bytes, offset);
  record.vz = readDouble(bytes, offset);
  record.mass = readDouble(bytes, offset);
  record.gas_density = readDouble(bytes, offset);
  record.gas_internal_energy = readDouble(bytes, offset);
  record.gas_metallicity = readDouble(bytes, offset);
  record.star_formation = readDouble(bytes, offset);
  record.star_birth_mass = readDouble(bytes, offset);
  record.star_metallicity = readDouble(bytes, offset);
  record.bh_mass = readDouble(bytes, offset);
  record.bh_mdot = readDouble(bytes, offset);
  record.tracer_parent = readLe64(bytes, offset);
  record.tracer_injection = readLe64(bytes, offset);
  record.tracer_host = readLe32(bytes, offset);
  record.tracer_fraction = readDouble(bytes, offset);
  record.tracer_last_host_mass = readDouble(bytes, offset);
  record.tracer_exchanged_mass = readDouble(bytes, offset);
  return record;
}

}  // namespace cosmosim::io::internal
