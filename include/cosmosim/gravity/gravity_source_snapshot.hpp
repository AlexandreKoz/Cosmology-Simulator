#pragma once

#include <cstdint>
#include <span>
#include <stdexcept>

#include "cosmosim/gravity/gravity_state_identity.hpp"

namespace cosmosim::gravity {

enum class GravitySourceRepresentation : std::uint8_t {
  kMaterializedGeneric = 0,
  kBorrowedHomogeneousDmo = 1,
};

// Immutable view of the authoritative gravity source state represented at one
// physical force epoch. Coordinates may be predicted staging owned by the
// workflow; canonical SimulationState need not be mutated to make inactive
// elements appear current.
struct GravitySourceSnapshot {
  std::span<const double> pos_x_comoving;
  std::span<const double> pos_y_comoving;
  std::span<const double> pos_z_comoving;
  std::span<const double> mass_code;
  std::span<const std::uint32_t> species_tag;
  GravitySourceGeneration generation{};
  ForceEvaluationEpoch evaluation_epoch{};
  std::uint64_t particle_index_generation = 0U;
  GravitySourceRepresentation representation =
      GravitySourceRepresentation::kMaterializedGeneric;
  bool contains_predicted_coordinates = false;

  [[nodiscard]] std::size_t size() const noexcept { return mass_code.size(); }
  void validate() const {
    if (pos_x_comoving.size() != mass_code.size() ||
        pos_y_comoving.size() != mass_code.size() ||
        pos_z_comoving.size() != mass_code.size() ||
        (!species_tag.empty() && species_tag.size() != mass_code.size())) {
      throw std::invalid_argument("GravitySourceSnapshot source lanes have mismatched extents");
    }
    if (!generation.valid() || !evaluation_epoch.valid()) {
      throw std::invalid_argument("GravitySourceSnapshot requires valid source generation and force epoch");
    }
    if (representation == GravitySourceRepresentation::kBorrowedHomogeneousDmo) {
      if (!species_tag.empty()) {
        throw std::invalid_argument(
            "borrowed homogeneous DMO GravitySourceSnapshot requires implicit species");
      }
      if (contains_predicted_coordinates) {
        throw std::invalid_argument(
            "borrowed canonical GravitySourceSnapshot cannot contain workflow-predicted coordinates");
      }
    }
  }
};

}  // namespace cosmosim::gravity
