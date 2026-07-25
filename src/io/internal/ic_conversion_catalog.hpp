#pragma once

#include <cstdint>
#include <string>
#include <string_view>

#include "cosmosim/core/config.hpp"
#include "cosmosim/core/units.hpp"
#include "cosmosim/io/ic_reader.hpp"

namespace cosmosim::io::internal {

struct IcSourceConvention {
  core::UnitSystem source_units;
  IcCoordinateFrame frame = IcCoordinateFrame::kComoving;
  IcVelocityConvention velocity = IcVelocityConvention::kPhysicalPeculiar;
  double length_hubble_exponent = 0.0;
  double length_scale_factor_exponent = 0.0;
  double mass_hubble_exponent = 0.0;
  double mass_scale_factor_exponent = 0.0;
  double velocity_hubble_exponent = 0.0;
  double velocity_scale_factor_exponent = 0.0;
};

struct IcFieldConversionContract {
  double base_unit_to_si = 1.0;
  double hubble_exponent = 0.0;
  double scale_factor_exponent = 0.0;
  std::int8_t length_power = 0;
  std::int8_t mass_power = 0;
  std::int8_t time_power = 0;
  double frame_scale_factor_exponent = 0.0;
  std::uint8_t velocity_convention_power = 0;
  std::string source_unit = "dimensionless";
  std::string target_unit = "dimensionless";
};

[[nodiscard]] IcCoordinateFrame mapBridgeFrame(
    core::InitialConditionCoordinateFrame frame);
[[nodiscard]] IcVelocityConvention mapBridgeVelocityConvention(
    core::InitialConditionVelocityConvention convention);
[[nodiscard]] IcFieldConversionContract fieldConversionContract(
    std::string_view canonical_path,
    IcFieldSemantics semantics,
    const IcSourceConvention& convention);

}  // namespace cosmosim::io::internal
