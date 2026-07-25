#pragma once

#include <filesystem>
#include <functional>
#include <span>

#include "cosmosim/core/config.hpp"
#include "cosmosim/io/ic_reader.hpp"
#include "io/internal/ic_record_codec.hpp"

namespace cosmosim::io::internal {

using IcManifestReadyCallback =
    std::function<void(const IcManifest& manifest)>;
using IcRecordBatchCallback =
    std::function<void(std::span<const IcParticleRecord> records)>;

[[nodiscard]] IcImportReport streamGadgetArepoHdf5Ic(
    const std::filesystem::path& ic_path,
    const core::SimulationConfig& config,
    const IcImportOptions& options,
    const IcManifestReadyCallback& on_manifest_ready,
    const IcRecordBatchCallback& on_record_batch);

}  // namespace cosmosim::io::internal
