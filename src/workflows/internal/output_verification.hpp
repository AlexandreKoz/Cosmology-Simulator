#pragma once

#include <string>
#include <string_view>

namespace cosmosim::core {
struct ProvenanceRecord;
struct SimulationConfig;
struct SimulationState;
}

namespace cosmosim::io {
struct SnapshotReadResult;
}

namespace cosmosim::workflows::internal {

struct SnapshotRoundtripVerification {
  bool ok = false;
  std::string detail;
};

[[nodiscard]] SnapshotRoundtripVerification verifySnapshotRoundtrip(
    const io::SnapshotReadResult& restored,
    const core::SimulationState& expected,
    const core::SimulationConfig& config,
    std::string_view expected_normalized_config,
    const core::ProvenanceRecord& expected_provenance);

[[nodiscard]] bool restartRuntimeStateExactlyEquivalent(
    const core::SimulationState& restored,
    const core::SimulationState& reference);

}  // namespace cosmosim::workflows::internal
