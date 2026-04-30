#pragma once

#include <array>
#include <cmath>
#include <cstdint>
#include <span>
#include <stdexcept>

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
  static constexpr std::size_t kSpeciesCapacity = 5;
  std::array<double, kSpeciesCapacity> epsilon_comoving_by_species{};
  bool enabled = false;
};

struct TreeSofteningView {
  std::span<const std::uint32_t> source_species_tag{};
  std::span<const double> source_particle_epsilon_comoving{};
  std::span<const std::uint8_t> source_particle_epsilon_override_mask{};
  std::span<const double> target_particle_epsilon_comoving{};
  std::span<const std::uint8_t> target_particle_epsilon_override_mask{};
  TreeSofteningSpeciesPolicy species_policy{};
};

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
        return view.source_particle_epsilon_comoving[source_index];
      }
    }
  }
  if (view.species_policy.enabled && !view.source_species_tag.empty()) {
    if (source_index >= view.source_species_tag.size()) {
      throw std::out_of_range("source species index out of range");
    }
    const std::size_t species = static_cast<std::size_t>(view.source_species_tag[source_index]);
    if (species < view.species_policy.epsilon_comoving_by_species.size()) {
      return view.species_policy.epsilon_comoving_by_species[species];
    }
  }
  return fallback.epsilon_comoving;
}

[[nodiscard]] inline double resolveTargetSofteningEpsilon(
    std::size_t target_active_slot,
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
        return view.target_particle_epsilon_comoving[target_active_slot];
      }
    }
  }
  return fallback.epsilon_comoving;
}

[[nodiscard]] inline double combineSofteningPairEpsilon(double epsilon_source_comoving, double epsilon_target_comoving) {
  // Conservative pair law: the interaction uses the larger of the two softenings.
  // This preserves symmetry and avoids over-hardening mixed-resolution/species pairs.
  return std::max(epsilon_source_comoving, epsilon_target_comoving);
}

[[nodiscard]] inline double softenedInvR3(double squared_distance, double epsilon_comoving) {
  const double epsilon2 = epsilon_comoving * epsilon_comoving;
  const double denominator = std::pow(squared_distance + epsilon2, 1.5);
  if (denominator <= 0.0) {
    return 0.0;
  }
  return 1.0 / denominator;
}

[[nodiscard]] inline double softenedInvR3(double squared_distance, const TreeSofteningPolicy& policy) {
  return softenedInvR3(squared_distance, policy.epsilon_comoving);
}

}  // namespace cosmosim::gravity
