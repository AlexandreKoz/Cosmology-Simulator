#include <cassert>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include "cosmosim/core/particle_species.hpp"
#include "io/internal/ic_record_codec.hpp"

namespace {

cosmosim::io::internal::IcParticleRecord makeRecord(std::uint32_t species) {
  cosmosim::io::internal::IcParticleRecord input;
  input.id = 0xfedcba9876543210ULL;
  input.species = species;
  input.x = 1.25;
  input.y = -2.5;
  input.z = 3.75;
  input.vx = 4.5;
  input.vy = -5.5;
  input.vz = 6.5;
  input.mass = 7.5;
  input.gas_density = 8.5;
  input.gas_internal_energy = 9.5;
  input.gas_metallicity = 0.015;
  input.star_formation = 0.25;
  input.star_birth_mass = 10.5;
  input.star_metallicity = 0.02;
  input.bh_mass = 11.5;
  input.bh_mdot = 0.125;
  input.tracer_parent = 42U;
  input.tracer_injection = 99U;
  input.tracer_host = 7U;
  input.tracer_fraction = 0.5;
  input.tracer_last_host_mass = 12.5;
  input.tracer_exchanged_mass = -0.75;
  return input;
}

void checkRoundtrip(std::uint32_t species) {
  const auto input = makeRecord(species);
  std::vector<std::uint8_t> bytes;
  cosmosim::io::internal::serializeIcRecord(input, bytes);
  assert(bytes.size() == cosmosim::io::internal::serializedIcRecordBytes(species));
  assert(bytes.size() <= cosmosim::io::internal::IC_WIRE_MAX_RECORD_BYTES);
  assert(cosmosim::io::internal::validateIcWireBuffer(bytes) == 1U);
  std::size_t offset = 0U;
  const auto output = cosmosim::io::internal::deserializeIcRecord(bytes, offset);
  assert(offset == bytes.size());
  assert(output.id == input.id && output.species == input.species);
  assert(output.x == input.x && output.y == input.y && output.z == input.z);
  assert(output.vx == input.vx && output.vy == input.vy && output.vz == input.vz);
  assert(output.mass == input.mass);

  using cosmosim::core::ParticleSpecies;
  if (species == static_cast<std::uint32_t>(ParticleSpecies::kGas)) {
    assert(output.gas_density == input.gas_density);
    assert(output.gas_internal_energy == input.gas_internal_energy);
    assert(output.gas_metallicity == input.gas_metallicity);
  } else if (species == static_cast<std::uint32_t>(ParticleSpecies::kStar)) {
    assert(output.star_formation == input.star_formation);
    assert(output.star_birth_mass == input.star_birth_mass);
    assert(output.star_metallicity == input.star_metallicity);
  } else if (species == static_cast<std::uint32_t>(ParticleSpecies::kBlackHole)) {
    assert(output.bh_mass == input.bh_mass && output.bh_mdot == input.bh_mdot);
  } else if (species == static_cast<std::uint32_t>(ParticleSpecies::kTracer)) {
    assert(output.tracer_parent == input.tracer_parent);
    assert(output.tracer_injection == input.tracer_injection);
    assert(output.tracer_host == input.tracer_host);
    assert(output.tracer_fraction == input.tracer_fraction);
    assert(output.tracer_last_host_mass == input.tracer_last_host_mass);
    assert(output.tracer_exchanged_mass == input.tracer_exchanged_mass);
  }

  bytes.pop_back();
  offset = 0U;
  bool rejected = false;
  try {
    (void)cosmosim::io::internal::deserializeIcRecord(bytes, offset);
  } catch (const std::runtime_error&) {
    rejected = true;
  }
  assert(rejected);
}

}  // namespace

int main() {
  for (std::uint32_t species = 0U;
       species < static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kCount);
       ++species) {
    checkRoundtrip(species);
  }
  assert(cosmosim::io::internal::serializedIcRecordBytes(
             static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kDarkMatter)) < 176U);

  bool invalid_rejected = false;
  try {
    std::vector<std::uint8_t> bytes;
    cosmosim::io::internal::serializeIcRecord(
        makeRecord(static_cast<std::uint32_t>(cosmosim::core::ParticleSpecies::kCount)), bytes);
  } catch (const std::runtime_error&) {
    invalid_rejected = true;
  }
  assert(invalid_rejected);
  return 0;
}
