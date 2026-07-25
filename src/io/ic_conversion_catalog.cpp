#include "io/internal/ic_conversion_catalog.hpp"

#include <cmath>
#include <stdexcept>

namespace cosmosim::io::internal {

IcCoordinateFrame mapBridgeFrame(
    core::InitialConditionCoordinateFrame frame) {
  switch (frame) {
    case core::InitialConditionCoordinateFrame::kComoving:
      return IcCoordinateFrame::kComoving;
    case core::InitialConditionCoordinateFrame::kPhysical:
      return IcCoordinateFrame::kPhysical;
    case core::InitialConditionCoordinateFrame::kUnspecified:
      break;
  }
  throw std::invalid_argument("IC bridge coordinate frame is unspecified");
}

IcVelocityConvention mapBridgeVelocityConvention(
    core::InitialConditionVelocityConvention convention) {
  switch (convention) {
    case core::InitialConditionVelocityConvention::kPhysicalPeculiar:
      return IcVelocityConvention::kPhysicalPeculiar;
    case core::InitialConditionVelocityConvention::kSqrtAScaledPeculiar:
      return IcVelocityConvention::kSqrtAScaledPeculiar;
    case core::InitialConditionVelocityConvention::kComovingCoordinateRate:
      return IcVelocityConvention::kComovingCoordinateRate;
    case core::InitialConditionVelocityConvention::kUnspecified:
      break;
  }
  throw std::invalid_argument("IC bridge velocity convention is unspecified");
}

IcFieldConversionContract fieldConversionContract(
    std::string_view canonical_path,
    IcFieldSemantics semantics,
    const IcSourceConvention& convention) {
  IcFieldConversionContract contract;
  const bool physical_frame =
      convention.frame == IcCoordinateFrame::kPhysical;
  if (canonical_path.ends_with("/Density")) {
    contract.base_unit_to_si =
        convention.source_units.mass_si_per_code /
        std::pow(convention.source_units.length_si_per_code, 3.0);
    contract.hubble_exponent = convention.mass_hubble_exponent -
        3.0 * convention.length_hubble_exponent;
    contract.scale_factor_exponent =
        convention.mass_scale_factor_exponent -
        3.0 * convention.length_scale_factor_exponent;
    contract.length_power = -3;
    contract.mass_power = 1;
    contract.frame_scale_factor_exponent = physical_frame ? 3.0 : 0.0;
    contract.source_unit = "source_mass/source_length^3";
    contract.target_unit = "runtime_mass/runtime_length^3";
    return contract;
  }
  if (canonical_path.ends_with("/BH_Mdot")) {
    contract.base_unit_to_si =
        convention.source_units.mass_si_per_code /
        convention.source_units.timeSiPerCode();
    contract.hubble_exponent = convention.mass_hubble_exponent +
        convention.velocity_hubble_exponent -
        convention.length_hubble_exponent;
    contract.scale_factor_exponent =
        convention.mass_scale_factor_exponent +
        convention.velocity_scale_factor_exponent -
        convention.length_scale_factor_exponent;
    contract.mass_power = 1;
    contract.time_power = -1;
    contract.frame_scale_factor_exponent = 0.0;
    contract.velocity_convention_power = 0U;
    contract.source_unit = "source_mass/source_time";
    contract.target_unit = "runtime_mass/runtime_time";
    return contract;
  }

  switch (semantics) {
    case IcFieldSemantics::kCoordinate:
      contract.base_unit_to_si = convention.source_units.length_si_per_code;
      contract.hubble_exponent = convention.length_hubble_exponent;
      contract.scale_factor_exponent =
          convention.length_scale_factor_exponent;
      contract.length_power = 1;
      contract.frame_scale_factor_exponent = physical_frame ? -1.0 : 0.0;
      contract.source_unit = "source_length";
      contract.target_unit = "runtime_comoving_length";
      break;
    case IcFieldSemantics::kVelocity:
      contract.base_unit_to_si = convention.source_units.velocity_si_per_code;
      contract.hubble_exponent = convention.velocity_hubble_exponent;
      contract.scale_factor_exponent =
          convention.velocity_scale_factor_exponent;
      contract.length_power = 1;
      contract.time_power = -1;
      contract.velocity_convention_power = 1U;
      contract.source_unit = "source_velocity";
      contract.target_unit = "runtime_peculiar_velocity";
      break;
    case IcFieldSemantics::kExtensive:
      contract.base_unit_to_si = convention.source_units.mass_si_per_code;
      contract.hubble_exponent = convention.mass_hubble_exponent;
      contract.scale_factor_exponent = convention.mass_scale_factor_exponent;
      contract.mass_power = 1;
      contract.source_unit = "source_mass";
      contract.target_unit = "runtime_mass";
      break;
    case IcFieldSemantics::kSpecific:
      contract.base_unit_to_si =
          convention.source_units.velocity_si_per_code *
          convention.source_units.velocity_si_per_code;
      contract.hubble_exponent =
          2.0 * convention.velocity_hubble_exponent;
      contract.scale_factor_exponent =
          2.0 * convention.velocity_scale_factor_exponent;
      contract.length_power = 2;
      contract.time_power = -2;
      contract.velocity_convention_power = 0U;
      contract.source_unit = "source_velocity^2";
      contract.target_unit = "runtime_specific_energy";
      break;
    case IcFieldSemantics::kIdentifier:
    case IcFieldSemantics::kIntensive:
      break;
  }
  return contract;
}

}  // namespace cosmosim::io::internal
