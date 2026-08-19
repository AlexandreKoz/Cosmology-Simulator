#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>

#include "cosmosim/core/particle_species.hpp"
#include "cosmosim/gravity/tree_softening.hpp"

namespace cosmosim::test_support {

// Test-only Plummer reference.  This deliberately spells out the documented
// force law rather than calling production softenedInvR3().
[[nodiscard]] inline double plummerInvR3Reference(
    double squared_distance_code,
    double epsilon_comoving) {
  if (!std::isfinite(squared_distance_code) || squared_distance_code < 0.0 ||
      !std::isfinite(epsilon_comoving) || epsilon_comoving < 0.0) {
    throw std::invalid_argument("independent Plummer reference received invalid input");
  }
  const double softened_r2 = squared_distance_code + epsilon_comoving * epsilon_comoving;
  if (softened_r2 <= 0.0) {
    return 0.0;
  }
  return 1.0 / (softened_r2 * std::sqrt(softened_r2));
}

[[nodiscard]] inline double plummerPotentialReference(
    double squared_distance_code,
    double epsilon_comoving) {
  if (!std::isfinite(squared_distance_code) || squared_distance_code < 0.0 ||
      !std::isfinite(epsilon_comoving) || epsilon_comoving < 0.0) {
    throw std::invalid_argument("independent Plummer potential reference received invalid input");
  }
  const double softened_r2 = squared_distance_code + epsilon_comoving * epsilon_comoving;
  if (softened_r2 <= 0.0) {
    return 0.0;
  }
  return 1.0 / std::sqrt(softened_r2);
}

[[nodiscard]] inline double pairSofteningMaxReference(
    double source_epsilon_comoving,
    double target_epsilon_comoving) {
  if (!std::isfinite(source_epsilon_comoving) || source_epsilon_comoving < 0.0 ||
      !std::isfinite(target_epsilon_comoving) || target_epsilon_comoving < 0.0) {
    throw std::invalid_argument("independent pair-softening reference received invalid input");
  }
  return std::max(source_epsilon_comoving, target_epsilon_comoving);
}

[[nodiscard]] inline double sourceSofteningReference(
    std::size_t source_index,
    const cosmosim::gravity::TreeSofteningPolicy& fallback,
    const cosmosim::gravity::TreeSofteningView& view) {
  if (!view.source_particle_epsilon_comoving.empty()) {
    if (source_index >= view.source_particle_epsilon_comoving.size()) {
      throw std::out_of_range("reference source softening index out of range");
    }
    if (!view.source_particle_epsilon_override_mask.empty()) {
      if (source_index >= view.source_particle_epsilon_override_mask.size()) {
        throw std::invalid_argument("reference source override-mask size mismatch");
      }
      if (view.source_particle_epsilon_override_mask[source_index] != 0U) {
        return view.source_particle_epsilon_comoving[source_index];
      }
    }
  }
  if (view.species_policy.enabled && !view.source_species_tag.empty()) {
    if (source_index >= view.source_species_tag.size()) {
      throw std::out_of_range("reference source species index out of range");
    }
    const std::uint32_t species = view.source_species_tag[source_index];
    if (!cosmosim::core::isValidParticleSpeciesTag(species)) {
      throw std::invalid_argument("reference source species tag out of range");
    }
    return view.species_policy.epsilon_comoving_by_species[species];
  }
  return fallback.epsilon_comoving;
}

[[nodiscard]] inline double targetSofteningReference(
    std::size_t target_active_slot,
    std::size_t target_source_index,
    const cosmosim::gravity::TreeSofteningPolicy& fallback,
    const cosmosim::gravity::TreeSofteningView& view) {
  if (!view.target_particle_epsilon_comoving.empty()) {
    if (target_active_slot >= view.target_particle_epsilon_comoving.size()) {
      throw std::out_of_range("reference target softening index out of range");
    }
    if (!view.target_particle_epsilon_override_mask.empty()) {
      if (target_active_slot >= view.target_particle_epsilon_override_mask.size()) {
        throw std::invalid_argument("reference target override-mask size mismatch");
      }
      if (view.target_particle_epsilon_override_mask[target_active_slot] != 0U) {
        return view.target_particle_epsilon_comoving[target_active_slot];
      }
    }
  } else if (!view.source_particle_epsilon_comoving.empty()) {
    if (target_source_index >= view.source_particle_epsilon_comoving.size()) {
      throw std::out_of_range("reference target source-index softening out of range");
    }
    if (!view.source_particle_epsilon_override_mask.empty()) {
      if (target_source_index >= view.source_particle_epsilon_override_mask.size()) {
        throw std::invalid_argument("reference target source override-mask mismatch");
      }
      if (view.source_particle_epsilon_override_mask[target_source_index] != 0U) {
        return view.source_particle_epsilon_comoving[target_source_index];
      }
    }
  }

  if (view.species_policy.enabled) {
    if (!view.target_species_tag.empty()) {
      if (target_active_slot >= view.target_species_tag.size()) {
        throw std::out_of_range("reference target species index out of range");
      }
      const std::uint32_t species = view.target_species_tag[target_active_slot];
      if (!cosmosim::core::isValidParticleSpeciesTag(species)) {
        throw std::invalid_argument("reference target species tag out of range");
      }
      return view.species_policy.epsilon_comoving_by_species[species];
    }
    if (!view.source_species_tag.empty()) {
      if (target_source_index >= view.source_species_tag.size()) {
        throw std::out_of_range("reference target source species index out of range");
      }
      const std::uint32_t species = view.source_species_tag[target_source_index];
      if (!cosmosim::core::isValidParticleSpeciesTag(species)) {
        throw std::invalid_argument("reference target source species tag out of range");
      }
      return view.species_policy.epsilon_comoving_by_species[species];
    }
  }
  return fallback.epsilon_comoving;
}

}  // namespace cosmosim::test_support
