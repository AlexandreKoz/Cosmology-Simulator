#pragma once

#include <string>

#include "cosmosim/core/config.hpp"
#include "cosmosim/core/units.hpp"
#include "cosmosim/io/snapshot_hdf5.hpp"

namespace cosmosim::io::internal {

struct SnapshotConversionContext {
  SnapshotDialect dialect = SnapshotDialect::kChuiNative;
  core::UnitSystem units;
  double scale_factor = 1.0;
  double hubble_param = 1.0;

  [[nodiscard]] double positionToStored(double value_code) const;
  [[nodiscard]] double positionFromStored(double value_stored) const;
  [[nodiscard]] double velocityToStored(double peculiar_velocity_code) const;
  [[nodiscard]] double velocityFromStored(double value_stored) const;
  [[nodiscard]] double massToStored(double mass_code) const;
  [[nodiscard]] double massFromStored(double value_stored) const;
  [[nodiscard]] double densityComovingToStored(double density_code) const;
  [[nodiscard]] double densityComovingFromStored(double value_stored) const;
  [[nodiscard]] double pressureComovingToStored(double pressure_code) const;
  [[nodiscard]] double pressureComovingFromStored(double value_stored) const;
  [[nodiscard]] double internalEnergyToStored(double internal_energy_code) const;
  [[nodiscard]] double internalEnergyFromStored(double value_stored) const;
  [[nodiscard]] double softeningComovingToStored(double softening_code) const;
  [[nodiscard]] double softeningComovingFromStored(double value_stored) const;
  [[nodiscard]] double starFormationRateCodeToStored(double mass_per_time_code) const;
  [[nodiscard]] double starFormationRateStoredToCode(double mass_per_year_msun) const;
  [[nodiscard]] double boxSizeMpcToStored(double box_size_mpc_comoving) const;
  [[nodiscard]] double boxSizeStoredToMpc(double box_size_stored) const;
};

[[nodiscard]] SnapshotDialect resolveSnapshotWriteDialect(
    SnapshotDialect requested,
    const core::SimulationConfig& config);

[[nodiscard]] SnapshotConversionContext makeSnapshotConversionContext(
    SnapshotDialect dialect,
    const core::SimulationConfig& config,
    double scale_factor);

[[nodiscard]] const char* snapshotDialectLabel(SnapshotDialect dialect);

}  // namespace cosmosim::io::internal
