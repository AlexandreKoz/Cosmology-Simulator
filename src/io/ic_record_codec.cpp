#include "io/internal/ic_record_codec.hpp"

#include <limits>
#include <stdexcept>

#include "cosmosim/core/particle_species.hpp"
#include "io/internal/ic_byte_codec.hpp"

namespace cosmosim::io::internal {
namespace {

[[nodiscard]] std::size_t recordBytesForSpecies(std::uint32_t species) {
  using cosmosim::core::ParticleSpecies;
  if (species == static_cast<std::uint32_t>(ParticleSpecies::kDarkMatter)) return IC_WIRE_DM_BYTES;
  if (species == static_cast<std::uint32_t>(ParticleSpecies::kGas)) return IC_WIRE_GAS_BYTES;
  if (species == static_cast<std::uint32_t>(ParticleSpecies::kStar)) return IC_WIRE_STAR_BYTES;
  if (species == static_cast<std::uint32_t>(ParticleSpecies::kBlackHole)) return IC_WIRE_BH_BYTES;
  if (species == static_cast<std::uint32_t>(ParticleSpecies::kTracer)) return IC_WIRE_TRACER_BYTES;
  throw std::runtime_error("IC wire record has invalid species tag");
}

[[nodiscard]] std::size_t peekRecordBytes(
    std::span<const std::uint8_t> bytes,
    std::size_t offset,
    std::uint32_t* species_out = nullptr) {
  if (offset > bytes.size() || bytes.size() - offset < 16U) {
    throw std::runtime_error("truncated IC wire v2 header");
  }
  std::size_t cursor = offset;
  const std::uint32_t magic = readLe32(bytes, cursor);
  const std::uint32_t version = readLe32(bytes, cursor);
  const std::uint32_t encoded_bytes = readLe32(bytes, cursor);
  const std::uint32_t species = readLe32(bytes, cursor);
  if (magic != IC_WIRE_MAGIC) throw std::runtime_error("IC wire magic mismatch");
  if (version != IC_WIRE_VERSION) throw std::runtime_error("unsupported IC wire version");
  const std::size_t expected = recordBytesForSpecies(species);
  if (encoded_bytes != expected) {
    throw std::runtime_error("IC wire record byte count does not match species schema");
  }
  if (encoded_bytes > bytes.size() - offset) {
    throw std::runtime_error("truncated IC wire record");
  }
  if (species_out != nullptr) *species_out = species;
  return encoded_bytes;
}

}  // namespace

std::size_t serializedIcRecordBytes(std::uint32_t species) {
  return recordBytesForSpecies(species);
}

void serializeIcRecord(
    const IcParticleRecord& record,
    std::vector<std::uint8_t>& output) {
  const std::size_t record_bytes = recordBytesForSpecies(record.species);
  const std::size_t begin = output.size();
  appendLe32(output, IC_WIRE_MAGIC);
  appendLe32(output, IC_WIRE_VERSION);
  appendLe32(output, static_cast<std::uint32_t>(record_bytes));
  appendLe32(output, record.species);
  appendLe64(output, record.id);
  appendDouble(output, record.x);
  appendDouble(output, record.y);
  appendDouble(output, record.z);
  appendDouble(output, record.vx);
  appendDouble(output, record.vy);
  appendDouble(output, record.vz);
  appendDouble(output, record.mass);

  using cosmosim::core::ParticleSpecies;
  if (record.species == static_cast<std::uint32_t>(ParticleSpecies::kGas)) {
    appendDouble(output, record.gas_density);
    appendDouble(output, record.gas_internal_energy);
    appendDouble(output, record.gas_metallicity);
  } else if (record.species == static_cast<std::uint32_t>(ParticleSpecies::kStar)) {
    appendDouble(output, record.star_formation);
    appendDouble(output, record.star_birth_mass);
    appendDouble(output, record.star_metallicity);
  } else if (record.species == static_cast<std::uint32_t>(ParticleSpecies::kBlackHole)) {
    appendDouble(output, record.bh_mass);
    appendDouble(output, record.bh_mdot);
  } else if (record.species == static_cast<std::uint32_t>(ParticleSpecies::kTracer)) {
    appendLe64(output, record.tracer_parent);
    appendLe64(output, record.tracer_injection);
    appendLe32(output, record.tracer_host);
    appendDouble(output, record.tracer_fraction);
    appendDouble(output, record.tracer_last_host_mass);
    appendDouble(output, record.tracer_exchanged_mass);
  }
  if (output.size() - begin != record_bytes) {
    throw std::logic_error("IC wire v2 record byte contract drifted");
  }
}

IcParticleRecord deserializeIcRecord(
    std::span<const std::uint8_t> bytes,
    std::size_t& offset) {
  std::uint32_t species = 0U;
  const std::size_t record_bytes = peekRecordBytes(bytes, offset, &species);
  const std::size_t record_end = offset + record_bytes;
  const std::uint32_t magic = readLe32(bytes, offset);
  const std::uint32_t version = readLe32(bytes, offset);
  const std::uint32_t encoded_bytes = readLe32(bytes, offset);
  static_cast<void>(magic);
  static_cast<void>(version);
  static_cast<void>(encoded_bytes);

  IcParticleRecord record;
  record.species = readLe32(bytes, offset);
  record.id = readLe64(bytes, offset);
  record.x = readDouble(bytes, offset);
  record.y = readDouble(bytes, offset);
  record.z = readDouble(bytes, offset);
  record.vx = readDouble(bytes, offset);
  record.vy = readDouble(bytes, offset);
  record.vz = readDouble(bytes, offset);
  record.mass = readDouble(bytes, offset);

  using cosmosim::core::ParticleSpecies;
  if (species == static_cast<std::uint32_t>(ParticleSpecies::kGas)) {
    record.gas_density = readDouble(bytes, offset);
    record.gas_internal_energy = readDouble(bytes, offset);
    record.gas_metallicity = readDouble(bytes, offset);
  } else if (species == static_cast<std::uint32_t>(ParticleSpecies::kStar)) {
    record.star_formation = readDouble(bytes, offset);
    record.star_birth_mass = readDouble(bytes, offset);
    record.star_metallicity = readDouble(bytes, offset);
  } else if (species == static_cast<std::uint32_t>(ParticleSpecies::kBlackHole)) {
    record.bh_mass = readDouble(bytes, offset);
    record.bh_mdot = readDouble(bytes, offset);
  } else if (species == static_cast<std::uint32_t>(ParticleSpecies::kTracer)) {
    record.tracer_parent = readLe64(bytes, offset);
    record.tracer_injection = readLe64(bytes, offset);
    record.tracer_host = readLe32(bytes, offset);
    record.tracer_fraction = readDouble(bytes, offset);
    record.tracer_last_host_mass = readDouble(bytes, offset);
    record.tracer_exchanged_mass = readDouble(bytes, offset);
  }
  if (offset != record_end) {
    throw std::runtime_error("IC wire v2 decoder did not consume declared record bytes");
  }
  return record;
}

std::size_t validateIcWireBuffer(std::span<const std::uint8_t> bytes) {
  std::size_t offset = 0U;
  std::size_t count = 0U;
  while (offset < bytes.size()) {
    const std::size_t record_bytes = peekRecordBytes(bytes, offset);
    offset += record_bytes;
    if (count == std::numeric_limits<std::size_t>::max()) {
      throw std::overflow_error("IC wire record count overflow");
    }
    ++count;
  }
  return count;
}

std::pair<std::size_t, std::size_t> lastIcWireRecordSpan(
    std::span<const std::uint8_t> bytes) {
  std::size_t offset = 0U;
  std::size_t last_begin = bytes.size();
  std::size_t last_size = 0U;
  while (offset < bytes.size()) {
    last_begin = offset;
    last_size = peekRecordBytes(bytes, offset);
    offset += last_size;
  }
  return {last_begin, last_size};
}

}  // namespace cosmosim::io::internal
