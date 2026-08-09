#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <span>
#include <vector>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/io/ic_reader.hpp"
#include "cosmosim/parallel/distributed_memory.hpp"
#include "io/internal/ic_record_codec.hpp"

namespace cosmosim::io::distributed_audit_internal {

#if COSMOSIM_ENABLE_HDF5 && COSMOSIM_ENABLE_MPI

[[nodiscard]] int ownerForX(double x, double box_size, int world_size);

[[nodiscard]] std::uint64_t nestedByteCapacity(
    const std::vector<std::vector<std::uint8_t>>& buckets);

void validateChunkCoverage(
    const parallel::MpiContext& mpi_context,
    std::size_t file_index,
    std::size_t type_index,
    std::size_t start,
    std::size_t count,
    int reader_rank);

[[nodiscard]] std::uint64_t exactDistributedChunkReconciliation(
    const parallel::MpiContext& mpi_context,
    std::span<const internal::IcParticleRecord> source_records,
    std::span<const internal::IcParticleRecord> final_records,
    IcImportCounters& counters);

[[nodiscard]] std::uint64_t exactDistributedIdAudit(
    const parallel::MpiContext& mpi_context,
    std::span<const std::uint64_t> local_ids,
    std::size_t batch_count,
    const std::filesystem::path& scratch_root,
    IcImportCounters& counters);

void validateDistributedTotals(
    const parallel::MpiContext& mpi_context,
    const core::SimulationState& state,
    const IcManifest& manifest,
    const std::array<double, 5>& local_source_mass);

#endif

}  // namespace cosmosim::io::distributed_audit_internal
