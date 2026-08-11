#pragma once

#include <array>
#include <cmath>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <string>

#include "cosmosim/core/particle_species.hpp"

namespace cosmosim::gravity {

// Conservative softening modes for short-range gravity force evaluation.
enum class TreeSofteningKernel {
  kPlummer,
};

struct TreeSofteningPolicy {
  TreeSofteningKernel kernel = TreeSofteningKernel::kPlummer;
  double epsilon_comoving = 0.0;
};

struct TreeSofteningSpeciesPolicy {
  static constexpr std::size_t kSpeciesCapacity = core::k_particle_species_count;
  std::array<double, kSpeciesCapacity> epsilon_comoving_by_species{};
  bool enabled = false;
};

struct TreeSofteningView {
  // Source spans are indexed by the tree/source particle row.
  std::span<const std::uint32_t> source_species_tag{};
  std::span<const double> source_particle_epsilon_comoving{};
  std::span<const std::uint8_t> source_particle_epsilon_override_mask{};

  // Target spans are indexed by active-set slot. When these spans are absent,
  // target resolution may fall back to source-indexed spans using the target particle row
  // supplied by TreeGravitySolver::evaluateActiveSet. This keeps compact active views
  // cheap while preserving species/default/override semantics for target particles.
  std::span<const std::uint32_t> target_species_tag{};
  std::span<const double> target_particle_epsilon_comoving{};
  std::span<const std::uint8_t> target_particle_epsilon_override_mask{};
  TreeSofteningSpeciesPolicy species_policy{};
};

inline void validateTreeSofteningPolicy(const TreeSofteningPolicy& policy) {
  if (policy.kernel != TreeSofteningKernel::kPlummer) {
    throw std::invalid_argument("Unsupported tree softening kernel");
  }
  if (!std::isfinite(policy.epsilon_comoving) || policy.epsilon_comoving < 0.0) {
    throw std::invalid_argument("Tree softening epsilon_comoving must be finite and non-negative");
  }
}

inline double validatedSofteningEpsilon(double epsilon_comoving, const char* context) {
  if (!std::isfinite(epsilon_comoving) || epsilon_comoving < 0.0) {
    throw std::invalid_argument(std::string(context) + " must be finite and non-negative");
  }
  return epsilon_comoving;
}

[[nodiscard]] inline double resolveSourceSofteningEpsilon(
    std::size_t source_index,
    const TreeSofteningPolicy& fallback,
    const TreeSofteningView& view) {
  if (!view.source_particle_epsilon_comoving.empty()) {
    if (source_index >= view.source_particle_epsilon_comoving.size()) {
      throw std::out_of_range("source softening index out of range");
    }
    if (!view.source_particle_epsilon_override_mask.empty()) {
      if (source_index >= view.source_particle_epsilon_override_mask.size()) {
        throw std::invalid_argument("source softening override mask has incompatible size");
      }
      if (view.source_particle_epsilon_override_mask[source_index] != 0U) {
        return validatedSofteningEpsilon(
            view.source_particle_epsilon_comoving[source_index], "source softening override");
      }
    }
  }
  if (view.species_policy.enabled && !view.source_species_tag.empty()) {
    if (source_index >= view.source_species_tag.size()) {
      throw std::out_of_range("source species index out of range");
    }
    const std::uint32_t species_tag = view.source_species_tag[source_index];
    if (!core::isValidParticleSpeciesTag(species_tag)) {
      throw std::invalid_argument("source species tag is outside the canonical species range");
    }
    return validatedSofteningEpsilon(
        view.species_policy.epsilon_comoving_by_species[species_tag],
        "source species softening");
  }
  validateTreeSofteningPolicy(fallback);
  return fallback.epsilon_comoving;
}

[[nodiscard]] inline double resolveTargetSofteningEpsilon(
    std::size_t target_active_slot,
    std::size_t target_source_index,
    const TreeSofteningPolicy& fallback,
    const TreeSofteningView& view) {
  if (!view.target_particle_epsilon_comoving.empty()) {
    if (target_active_slot >= view.target_particle_epsilon_comoving.size()) {
      throw std::out_of_range("target softening index out of range");
    }
    if (!view.target_particle_epsilon_override_mask.empty()) {
      if (target_active_slot >= view.target_particle_epsilon_override_mask.size()) {
        throw std::invalid_argument("target softening override mask has incompatible size");
      }
      if (view.target_particle_epsilon_override_mask[target_active_slot] != 0U) {
        return validatedSofteningEpsilon(
            view.target_particle_epsilon_comoving[target_active_slot], "target softening override");
      }
    }
  } else if (!view.source_particle_epsilon_comoving.empty()) {
    if (target_source_index >= view.source_particle_epsilon_comoving.size()) {
      throw std::out_of_range("target source-index softening index out of range");
    }
    if (!view.source_particle_epsilon_override_mask.empty()) {
      if (target_source_index >= view.source_particle_epsilon_override_mask.size()) {
        throw std::invalid_argument("target source-index softening override mask has incompatible size");
      }
      if (view.source_particle_epsilon_override_mask[target_source_index] != 0U) {
        return validatedSofteningEpsilon(
            view.source_particle_epsilon_comoving[target_source_index], "target source-index softening override");
      }
    }
  }

  if (view.species_policy.enabled) {
    if (!view.target_species_tag.empty()) {
      if (target_active_slot >= view.target_species_tag.size()) {
        throw std::out_of_range("target species index out of range");
      }
      const std::uint32_t species_tag = view.target_species_tag[target_active_slot];
      if (!core::isValidParticleSpeciesTag(species_tag)) {
        throw std::invalid_argument("target species tag is outside the canonical species range");
      }
      return validatedSofteningEpsilon(
          view.species_policy.epsilon_comoving_by_species[species_tag],
          "target species softening");
    } else if (!view.source_species_tag.empty()) {
      if (target_source_index >= view.source_species_tag.size()) {
        throw std::out_of_range("target source-index species index out of range");
      }
      const std::uint32_t species_tag = view.source_species_tag[target_source_index];
      if (!core::isValidParticleSpeciesTag(species_tag)) {
        throw std::invalid_argument("target source-index species tag is outside the canonical species range");
      }
      return validatedSofteningEpsilon(
          view.species_policy.epsilon_comoving_by_species[species_tag],
          "target source-index species softening");
    }
  }
  validateTreeSofteningPolicy(fallback);
  return fallback.epsilon_comoving;
}

[[nodiscard]] inline double resolveTargetSofteningEpsilon(
    std::size_t target_active_slot,
    const TreeSofteningPolicy& fallback,
    const TreeSofteningView& view) {
  return resolveTargetSofteningEpsilon(target_active_slot, target_active_slot, fallback, view);
}

[[nodiscard]] inline double combineSofteningPairEpsilonUnchecked(
    double epsilon_source_comoving,
    double epsilon_target_comoving) noexcept {
  return std::max(epsilon_source_comoving, epsilon_target_comoving);
}

[[nodiscard]] inline double combineSofteningPairEpsilon(
    double epsilon_source_comoving,
    double epsilon_target_comoving) {
  validatedSofteningEpsilon(epsilon_source_comoving, "source pair softening");
  validatedSofteningEpsilon(epsilon_target_comoving, "target pair softening");
  return combineSofteningPairEpsilonUnchecked(epsilon_source_comoving, epsilon_target_comoving);
}

[[nodiscard]] inline double softenedInvR3Unchecked(
    double squared_distance,
    double epsilon_comoving) noexcept {
  const double epsilon2 = epsilon_comoving * epsilon_comoving;
  const double softened_r2 = squared_distance + epsilon2;
  if (softened_r2 <= 0.0) {
    return 0.0;
  }
  const double denominator = softened_r2 * std::sqrt(softened_r2);
  return 1.0 / denominator;
}

[[nodiscard]] inline double softenedInvR3(double squared_distance, double epsilon_comoving) {
  if (!std::isfinite(squared_distance) || squared_distance < 0.0) {
    throw std::invalid_argument("squared_distance must be finite and non-negative");
  }
  validatedSofteningEpsilon(epsilon_comoving, "softening epsilon_comoving");
  return softenedInvR3Unchecked(squared_distance, epsilon_comoving);
}

[[nodiscard]] inline double softenedInvR3(double squared_distance, const TreeSofteningPolicy& policy) {
  validateTreeSofteningPolicy(policy);
  return softenedInvR3(squared_distance, policy.epsilon_comoving);
}

}  // namespace cosmosim::gravity
