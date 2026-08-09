#include "io/internal/snapshot_conversion.hpp"

#include <cmath>
#include <stdexcept>

namespace cosmosim::io::internal {
namespace {

constexpr double k_solar_mass_si = 1.98847e30;
constexpr double k_julian_year_si = 31557600.0;

[[nodiscard]] bool usesArepoCosmologicalScaling(SnapshotDialect dialect) {
  return dialect == SnapshotDialect::kArepoFormat3 ||
      dialect == SnapshotDialect::kGadget4Hdf5;
}

void requireValidCosmology(const SnapshotConversionContext& context) {
  if (!(context.scale_factor > 0.0) || !std::isfinite(context.scale_factor)) {
    throw std::invalid_argument("snapshot conversion requires finite positive scale factor");
  }
  if (!(context.hubble_param > 0.0) || !std::isfinite(context.hubble_param)) {
    throw std::invalid_argument("snapshot conversion requires finite positive HubbleParam");
  }
}

}  // namespace

const char* snapshotDialectLabel(SnapshotDialect dialect) {
  switch (dialect) {
    case SnapshotDialect::kAuto: return "auto";
    case SnapshotDialect::kChuiNative: return "chui_native";
    case SnapshotDialect::kArepoFormat3: return "arepo_format3";
    case SnapshotDialect::kGadget4Hdf5: return "gadget4_hdf5";
  }
  throw std::logic_error("unhandled snapshot dialect");
}

SnapshotDialect resolveSnapshotWriteDialect(
    SnapshotDialect requested,
    const core::SimulationConfig& config) {
  if (requested != SnapshotDialect::kAuto) {
    if ((requested == SnapshotDialect::kArepoFormat3 ||
         requested == SnapshotDialect::kGadget4Hdf5) &&
        config.units.coordinate_frame != core::CoordinateFrame::kComoving) {
      throw std::invalid_argument(
          "AREPO/GADGET cosmological snapshot export requires the CHUI comoving coordinate frame; "
          "use chui_native for isolated physical-coordinate runs");
    }
    return requested;
  }
  return config.units.coordinate_frame == core::CoordinateFrame::kComoving
      ? SnapshotDialect::kArepoFormat3
      : SnapshotDialect::kChuiNative;
}

SnapshotConversionContext makeSnapshotConversionContext(
    SnapshotDialect dialect,
    const core::SimulationConfig& config,
    double scale_factor) {
  SnapshotConversionContext context;
  context.dialect = dialect;
  context.units = core::makeUnitSystem(
      config.units.length_unit, config.units.mass_unit, config.units.velocity_unit);
  context.scale_factor = scale_factor;
  context.hubble_param = config.cosmology.hubble_param;
  context.units.validate();
  if (usesArepoCosmologicalScaling(dialect)) {
    requireValidCosmology(context);
  }
  return context;
}

double SnapshotConversionContext::positionToStored(double value_code) const {
  return usesArepoCosmologicalScaling(dialect) ? value_code * hubble_param : value_code;
}

double SnapshotConversionContext::positionFromStored(double value_stored) const {
  return usesArepoCosmologicalScaling(dialect) ? value_stored / hubble_param : value_stored;
}

double SnapshotConversionContext::velocityToStored(double peculiar_velocity_code) const {
  return usesArepoCosmologicalScaling(dialect)
      ? peculiar_velocity_code / std::sqrt(scale_factor)
      : peculiar_velocity_code;
}

double SnapshotConversionContext::velocityFromStored(double value_stored) const {
  return usesArepoCosmologicalScaling(dialect)
      ? value_stored * std::sqrt(scale_factor)
      : value_stored;
}

double SnapshotConversionContext::massToStored(double mass_code) const {
  return usesArepoCosmologicalScaling(dialect) ? mass_code * hubble_param : mass_code;
}

double SnapshotConversionContext::massFromStored(double value_stored) const {
  return usesArepoCosmologicalScaling(dialect) ? value_stored / hubble_param : value_stored;
}

double SnapshotConversionContext::densityComovingToStored(double density_code) const {
  if (!usesArepoCosmologicalScaling(dialect)) {
    return density_code;
  }
  return density_code / (hubble_param * hubble_param);
}

double SnapshotConversionContext::densityComovingFromStored(double value_stored) const {
  if (!usesArepoCosmologicalScaling(dialect)) {
    return value_stored;
  }
  return value_stored * hubble_param * hubble_param;
}

double SnapshotConversionContext::pressureComovingToStored(double pressure_code) const {
  if (!usesArepoCosmologicalScaling(dialect)) {
    return pressure_code;
  }
  return pressure_code / (hubble_param * hubble_param);
}

double SnapshotConversionContext::pressureComovingFromStored(double value_stored) const {
  if (!usesArepoCosmologicalScaling(dialect)) {
    return value_stored;
  }
  return value_stored * hubble_param * hubble_param;
}

double SnapshotConversionContext::internalEnergyToStored(double internal_energy_code) const {
  return internal_energy_code;
}

double SnapshotConversionContext::internalEnergyFromStored(double value_stored) const {
  return value_stored;
}

double SnapshotConversionContext::softeningComovingToStored(double softening_code) const {
  return positionToStored(softening_code);
}

double SnapshotConversionContext::softeningComovingFromStored(double value_stored) const {
  return positionFromStored(value_stored);
}

double SnapshotConversionContext::starFormationRateCodeToStored(
    double mass_per_time_code) const {
  const double mass_per_second_si =
      mass_per_time_code * units.mass_si_per_code / units.timeSiPerCode();
  return mass_per_second_si * k_julian_year_si / k_solar_mass_si;
}

double SnapshotConversionContext::starFormationRateStoredToCode(
    double mass_per_year_msun) const {
  const double mass_per_second_si =
      mass_per_year_msun * k_solar_mass_si / k_julian_year_si;
  return mass_per_second_si * units.timeSiPerCode() / units.mass_si_per_code;
}

double SnapshotConversionContext::boxSizeMpcToStored(double box_size_mpc_comoving) const {
  const core::UnitSystem mpc = core::makeUnitSystem("mpc", "msun", "km_s");
  const double box_code =
      box_size_mpc_comoving * mpc.length_si_per_code / units.length_si_per_code;
  return positionToStored(box_code);
}

double SnapshotConversionContext::boxSizeStoredToMpc(double box_size_stored) const {
  const core::UnitSystem mpc = core::makeUnitSystem("mpc", "msun", "km_s");
  const double box_code = positionFromStored(box_size_stored);
  return box_code * units.length_si_per_code / mpc.length_si_per_code;
}

}  // namespace cosmosim::io::internal
