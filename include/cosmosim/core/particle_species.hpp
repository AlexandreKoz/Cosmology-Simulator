#pragma once

#include <cstddef>
#include <cstdint>
#include <stdexcept>

namespace cosmosim::core {

// Canonical species tags used by persistent state, gravity softening, I/O, and accounting.
enum class ParticleSpecies : std::uint8_t {
  kDarkMatter = 0,
  kGas = 1,
  kStar = 2,
  kBlackHole = 3,
  kTracer = 4,
  kCount = 5,
};

inline constexpr std::size_t k_particle_species_count =
    static_cast<std::size_t>(ParticleSpecies::kCount);

[[nodiscard]] constexpr std::size_t particleSpeciesIndex(ParticleSpecies species) {
  const auto index = static_cast<std::size_t>(species);
  if (index >= k_particle_species_count) {
    throw std::out_of_range("ParticleSpecies value is outside the canonical species range");
  }
  return index;
}

[[nodiscard]] constexpr bool isValidParticleSpeciesTag(std::uint32_t species_tag) noexcept {
  return species_tag < k_particle_species_count;
}

}  // namespace cosmosim::core
