#pragma once

#include <filesystem>

#include "cosmosim/io/snapshot_hdf5.hpp"

namespace cosmosim::io::internal {

void writeSnapshotMemberIntegritySidecar(
    const std::filesystem::path& member_path,
    const SnapshotSetMemberInfo& member,
    bool durable_publication);

}  // namespace cosmosim::io::internal
