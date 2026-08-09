#pragma once

#include <array>
#include <cstddef>
#include <utility>

namespace cosmosim::io::internal {

enum class SnapshotFieldStorage {
  kFloat64,
  kUint64,
  kUint32,
  kUint8,
};

struct SnapshotFieldContract {
  const char* name = "";
  SnapshotFieldStorage storage = SnapshotFieldStorage::kFloat64;
  bool require_finite = false;
  bool metallicity_fraction = false;
};

inline constexpr std::array<SnapshotFieldContract, 14> kRequiredGasSnapshotFields{{
    {"InternalEnergy", SnapshotFieldStorage::kFloat64, true, false},
    {"Density", SnapshotFieldStorage::kFloat64, true, false},
    {"Metallicity", SnapshotFieldStorage::kFloat64, true, true},
    {"CHUI_TemperatureCode", SnapshotFieldStorage::kFloat64, true, false},
    {"CHUI_SoundSpeedCode", SnapshotFieldStorage::kFloat64, true, false},
    {"StarFormationRate", SnapshotFieldStorage::kFloat64, true, false},
    {"ColdCloudMassFraction", SnapshotFieldStorage::kFloat64, true, false},
    {"EffectivePressure", SnapshotFieldStorage::kFloat64, true, false},
    {"EffectiveInternalEnergy", SnapshotFieldStorage::kFloat64, true, false},
    {"IsOnEffectiveEos", SnapshotFieldStorage::kUint8, false, false},
    {"GasCellIDs", SnapshotFieldStorage::kUint64, false, false},
    {"ParentParticleIDs", SnapshotFieldStorage::kUint64, false, false},
    {"HasParentParticle", SnapshotFieldStorage::kUint8, false, false},
    {"OwningPatchIDs", SnapshotFieldStorage::kUint64, false, false},
}};

inline constexpr std::array<SnapshotFieldContract, 15> kRequiredStarSnapshotFields{{
    {"Metallicity", SnapshotFieldStorage::kFloat64, true, true},
    {"StellarFormationTime", SnapshotFieldStorage::kFloat64, true, false},
    {"BirthMass", SnapshotFieldStorage::kFloat64, true, false},
    {"StarFormationBirthKey", SnapshotFieldStorage::kUint64, false, false},
    {"ParentGasCellID", SnapshotFieldStorage::kUint64, false, false},
    {"BirthIntegrationTick", SnapshotFieldStorage::kUint64, false, false},
    {"BirthOrdinal", SnapshotFieldStorage::kUint32, false, false},
    {"CHUI_StellarAgeYearsLast", SnapshotFieldStorage::kFloat64, true, false},
    {"CHUI_StellarReturnedMassCumulative", SnapshotFieldStorage::kFloat64, true, false},
    {"CHUI_StellarReturnedMetalsCumulative", SnapshotFieldStorage::kFloat64, true, false},
    {"CHUI_StellarNewlySynthesizedMetalsCumulative", SnapshotFieldStorage::kFloat64, true, false},
    {"CHUI_StellarFeedbackEnergyCumulativeErg", SnapshotFieldStorage::kFloat64, true, false},
    {"CHUI_StellarDepositedMassCumulative", SnapshotFieldStorage::kFloat64, true, false},
    {"CHUI_StellarDepositedMetalsCumulative", SnapshotFieldStorage::kFloat64, true, false},
    {"CHUI_StellarDepositedFeedbackEnergyCumulativeErg", SnapshotFieldStorage::kFloat64, true, false},
}};

inline constexpr std::array<SnapshotFieldContract, 6> kRequiredTracerSnapshotFields{{
    {"TracerParentParticleID", SnapshotFieldStorage::kUint64, false, false},
    {"TracerInjectionStep", SnapshotFieldStorage::kUint64, false, false},
    {"TracerHostCellIndex", SnapshotFieldStorage::kUint32, false, false},
    {"TracerMassFractionOfHost", SnapshotFieldStorage::kFloat64, true, false},
    {"TracerLastHostMassCode", SnapshotFieldStorage::kFloat64, true, false},
    {"TracerCumulativeExchangedMassCode", SnapshotFieldStorage::kFloat64, true, false},
}};

inline constexpr std::array<SnapshotFieldContract, 9> kRequiredBlackHoleSnapshotFields{{
    {"CHUI_BHSubgridMass", SnapshotFieldStorage::kFloat64, true, false},
    {"CHUI_BHAccretionRateMsunPerYr", SnapshotFieldStorage::kFloat64, true, false},
    {"CHUI_BHFeedbackEnergyCode", SnapshotFieldStorage::kFloat64, true, false},
    {"CHUI_BHEddingtonRatio", SnapshotFieldStorage::kFloat64, true, false},
    {"CHUI_BHCumulativeAccretedMass", SnapshotFieldStorage::kFloat64, true, false},
    {"CHUI_BHCumulativeFeedbackEnergyCode", SnapshotFieldStorage::kFloat64, true, false},
    {"CHUI_BHDutyCycleActiveTimeCode", SnapshotFieldStorage::kFloat64, true, false},
    {"CHUI_BHDutyCycleTotalTimeCode", SnapshotFieldStorage::kFloat64, true, false},
    {"CHUI_BHHostCellIndex", SnapshotFieldStorage::kUint32, false, false},
}};

template <typename Function>
void forEachRequiredChuiSnapshotField(std::size_t part_type, Function&& function) {
  const auto visit = [&](const auto& fields) {
    for (const auto& field : fields) {
      function(field);
    }
  };
  switch (part_type) {
    case 0U: visit(kRequiredGasSnapshotFields); break;
    case 3U: visit(kRequiredTracerSnapshotFields); break;
    case 4U: visit(kRequiredStarSnapshotFields); break;
    case 5U: visit(kRequiredBlackHoleSnapshotFields); break;
    default: break;
  }
}

}  // namespace cosmosim::io::internal
