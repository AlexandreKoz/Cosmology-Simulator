#pragma once

#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <vector>

namespace cosmosim::io::internal {

inline constexpr std::size_t kIcWireRecordBytes = 168U;

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

void serializeIcRecord(
    const IcParticleRecord& record,
    std::vector<std::uint8_t>& output);

[[nodiscard]] IcParticleRecord deserializeIcRecord(
    std::span<const std::uint8_t> bytes,
    std::size_t& offset);

}  // namespace cosmosim::io::internal
