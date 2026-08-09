#pragma once

#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <utility>
#include <vector>

namespace cosmosim::io::internal {

// Wire v2 is a compact, endian-stable, self-framing record contract.  Every
// record carries its byte length and species so dark-matter traffic does not
// pay for gas/star/BH/tracer sidecars.
inline constexpr std::uint32_t IC_WIRE_MAGIC = 0x32434943U;  // "CIC2" little-endian on wire
inline constexpr std::uint32_t IC_WIRE_VERSION = 2U;
inline constexpr std::size_t IC_WIRE_COMMON_BYTES = 80U;
inline constexpr std::size_t IC_WIRE_DM_BYTES = IC_WIRE_COMMON_BYTES;
inline constexpr std::size_t IC_WIRE_GAS_BYTES = 104U;
inline constexpr std::size_t IC_WIRE_STAR_BYTES = 104U;
inline constexpr std::size_t IC_WIRE_BH_BYTES = 96U;
inline constexpr std::size_t IC_WIRE_TRACER_BYTES = 124U;
inline constexpr std::size_t IC_WIRE_MAX_RECORD_BYTES = IC_WIRE_TRACER_BYTES;

struct IcParticleRecord {
  std::uint64_t id = 0;
  std::uint32_t species = 0;
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;
  double vx = 0.0;
  double vy = 0.0;
  double vz = 0.0;
  double mass = 0.0;
  double gas_density = 0.0;
  double gas_internal_energy = 0.0;
  double gas_metallicity = 0.0;
  double star_formation = 0.0;
  double star_birth_mass = 0.0;
  double star_metallicity = 0.0;
  double bh_mass = 0.0;
  double bh_mdot = 0.0;
  std::uint64_t tracer_parent = 0;
  std::uint64_t tracer_injection = 0;
  std::uint32_t tracer_host = std::numeric_limits<std::uint32_t>::max();
  double tracer_fraction = 0.0;
  double tracer_last_host_mass = 0.0;
  double tracer_exchanged_mass = 0.0;
};

[[nodiscard]] std::size_t serializedIcRecordBytes(std::uint32_t species);

void serializeIcRecord(
    const IcParticleRecord& record,
    std::vector<std::uint8_t>& output);

[[nodiscard]] IcParticleRecord deserializeIcRecord(
    std::span<const std::uint8_t> bytes,
    std::size_t& offset);

// Validate all framing without constructing records; returns record count.
[[nodiscard]] std::size_t validateIcWireBuffer(std::span<const std::uint8_t> bytes);

// Returns [begin, size] for the final complete record, or {size,0} when empty.
[[nodiscard]] std::pair<std::size_t, std::size_t> lastIcWireRecordSpan(
    std::span<const std::uint8_t> bytes);

}  // namespace cosmosim::io::internal
