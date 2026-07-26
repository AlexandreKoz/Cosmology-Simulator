#include "io/internal/ic_distributed_audit.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

#include "io/internal/ic_byte_codec.hpp"
#include "io/internal/ic_failure_protocol.hpp"
#include "io/internal/ic_file_set_common.hpp"
#include "io/internal/ic_mpi_collectives.hpp"

#if COSMOSIM_ENABLE_MPI
#include <mpi.h>
#endif

namespace cosmosim::io::distributed_audit_internal {

#if COSMOSIM_ENABLE_HDF5 && COSMOSIM_ENABLE_MPI

using file_set_internal::ParticleRecord;
using file_set_internal::checkedCounterAdd;
using file_set_internal::kParticleTypeCount;
using file_set_internal::speciesTag;
using failure_protocol_internal::alltoallBytes;
using failure_protocol_internal::injectIcTestFault;
using failure_protocol_internal::runCollectivePhase;
using failure_protocol_internal::runCollectivePhaseVoid;
using mpi_collective_internal::mpiAllreduce;

[[nodiscard]] int ownerForX(double x, double box_size, int world_size) {
  if (!std::isfinite(x) || !(box_size > 0.0) || x < 0.0 || x > box_size * (1.0 + 1.0e-12)) throw std::runtime_error("converted IC coordinate is outside the periodic box");
  if (x >= box_size) x = 0.0;
  const double fraction=x/box_size;const int owner=std::min(world_size-1,static_cast<int>(fraction*world_size));return std::max(0,owner);
}

[[nodiscard]] std::uint64_t checkedCapacityBytes(
    std::size_t capacity,
    std::size_t element_size,
    std::string_view label) {
  if (element_size != 0U &&
      capacity > std::numeric_limits<std::uint64_t>::max() / element_size) {
    throw std::overflow_error(std::string(label) + " capacity overflow");
  }
  return static_cast<std::uint64_t>(capacity) * element_size;
}

[[nodiscard]] std::uint64_t nestedByteCapacity(
    const std::vector<std::vector<std::uint8_t>>& buckets) {
  std::uint64_t bytes = checkedCapacityBytes(
      buckets.capacity(), sizeof(std::vector<std::uint8_t>),
      "IC nested bucket table");
  for (const auto& bucket : buckets) {
    if (bytes > std::numeric_limits<std::uint64_t>::max() -
            bucket.capacity()) {
      throw std::overflow_error("IC nested byte capacity overflow");
    }
    bytes += static_cast<std::uint64_t>(bucket.capacity());
  }
  return bytes;
}

void validateChunkCoverage(
    const parallel::MpiContext& mpi_context,
    std::size_t file_index,
    std::size_t type_index,
    std::size_t start,
    std::size_t count,
    int reader_rank) {
  const int local_completed =
      mpi_context.worldRank() == reader_rank ? 1 : 0;
  int global_completed = 0;
  mpiAllreduce(
      &local_completed, &global_completed, 1, MPI_INT, MPI_SUM,
      MPI_COMM_WORLD);
  if (global_completed != 1) {
    throw std::runtime_error(
        "distributed IC chunk coverage is not exactly one reader");
  }
  const std::uint64_t token =
      (static_cast<std::uint64_t>(file_index) * 0x9e3779b97f4a7c15ULL) ^
      (static_cast<std::uint64_t>(type_index) * 0xbf58476d1ce4e5b9ULL) ^
      (static_cast<std::uint64_t>(start) * 0x94d049bb133111ebULL) ^
      static_cast<std::uint64_t>(count);
  std::uint64_t minimum = 0U;
  std::uint64_t maximum = 0U;
  mpiAllreduce(&token, &minimum, 1, MPI_UINT64_T, MPI_MIN, MPI_COMM_WORLD);
  mpiAllreduce(&token, &maximum, 1, MPI_UINT64_T, MPI_MAX, MPI_COMM_WORLD);
  if (minimum != maximum) {
    throw std::runtime_error(
        "distributed IC ranks disagree on the active chunk token");
  }
}

[[nodiscard]] std::uint64_t exactDistributedChunkReconciliation(
    const parallel::MpiContext& mpi_context,
    std::span<const ParticleRecord> source_records,
    std::span<const ParticleRecord> final_records,
    IcImportCounters& counters) {
  std::vector<std::vector<std::uint8_t>> buckets =
      runCollectivePhase<std::vector<std::vector<std::uint8_t>>>(
          mpi_context, "IC source-to-final audit bucket preparation", [&]() {
            injectIcTestFault(mpi_context, "source_final_reconciliation");
            std::vector<std::vector<std::uint8_t>> prepared(
                static_cast<std::size_t>(mpi_context.worldSize()));
            const auto append_entry = [&](std::uint64_t id,
                                          std::uint32_t source_count,
                                          std::uint32_t final_count) {
              const int validator = static_cast<int>(
                  (id ^ (id >> 33U) ^ (id << 11U)) %
                  static_cast<std::uint64_t>(mpi_context.worldSize()));
              auto& bucket = prepared[static_cast<std::size_t>(validator)];
              internal::appendLe64(bucket, id);
              internal::appendLe32(bucket, source_count);
              internal::appendLe32(bucket, final_count);
            };
            for (const ParticleRecord& record : source_records) {
              append_entry(record.id, 1U, 0U);
            }
            for (const ParticleRecord& record : final_records) {
              append_entry(record.id, 0U, 1U);
            }
            return prepared;
          });

  std::uint64_t sent = 0U;
  std::uint64_t received_bytes = 0U;
  std::uint64_t exchange_peak = 0U;
  const std::vector<std::uint8_t> received = alltoallBytes(
      mpi_context, buckets, sent, received_bytes, &exchange_peak);
  runCollectivePhaseVoid(
      mpi_context, "IC source-to-final audit accounting", [&]() {
        injectIcTestFault(mpi_context, "source_final_accounting");
        if (mpi_context.isRoot()) {
          checkedCounterAdd(
              counters.exact_audit_exchange_count, 1U,
              "exact_audit_exchange_count");
        }
        checkedCounterAdd(counters.bytes_sent, sent, "bytes_sent");
        checkedCounterAdd(
            counters.bytes_received, received_bytes, "bytes_received");
      });

  struct BalanceEntry {
    std::uint64_t id = 0U;
    std::uint32_t source_count = 0U;
    std::uint32_t final_count = 0U;
  };
  std::vector<BalanceEntry> entries =
      runCollectivePhase<std::vector<BalanceEntry>>(
          mpi_context, "IC source-to-final audit decode", [&]() {
            if (received.size() % 16U != 0U) {
              throw std::runtime_error(
                  "IC source-to-final audit wire size is corrupt");
            }
            std::vector<BalanceEntry> decoded;
            decoded.reserve(received.size() / 16U);
            std::size_t offset = 0U;
            while (offset < received.size()) {
              decoded.push_back(BalanceEntry{
                  .id = internal::readLe64(received, offset),
                  .source_count = internal::readLe32(received, offset),
                  .final_count = internal::readLe32(received, offset)});
            }
            return decoded;
          });

  const int local_bad = runCollectivePhase<int>(
      mpi_context, "IC source-to-final audit reduction", [&]() {
        std::sort(
            entries.begin(), entries.end(),
            [](const BalanceEntry& lhs, const BalanceEntry& rhs) {
              return lhs.id < rhs.id;
            });
        std::size_t begin = 0U;
        while (begin < entries.size()) {
          std::size_t end = begin + 1U;
          std::uint64_t source_count = entries[begin].source_count;
          std::uint64_t final_count = entries[begin].final_count;
          while (end < entries.size() &&
                 entries[end].id == entries[begin].id) {
            source_count += entries[end].source_count;
            final_count += entries[end].final_count;
            ++end;
          }
          if (source_count != 1U || final_count != 1U) {
            return 1;
          }
          begin = end;
        }
        return 0;
      });
  int any_bad = 0;
  mpiAllreduce(&local_bad, &any_bad, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  if (any_bad != 0) {
    throw std::runtime_error(
        "distributed IC source and final ID multisets do not reconcile "
        "exactly");
  }
  return runCollectivePhase<std::uint64_t>(
      mpi_context, "IC source-to-final audit capacity accounting", [&]() {
        injectIcTestFault(mpi_context, "exact_audit_capacity_accounting");
        std::uint64_t exchange_storage = nestedByteCapacity(buckets);
        checkedCounterAdd(
            exchange_storage, exchange_peak,
            "source-to-final audit exchange storage");
        std::uint64_t decoded_storage =
            static_cast<std::uint64_t>(received.capacity());
        checkedCounterAdd(
            decoded_storage,
            checkedCapacityBytes(
                entries.capacity(), sizeof(BalanceEntry),
                "source-to-final audit decoded storage"),
            "source-to-final audit decoded storage");
        return std::max(exchange_storage, decoded_storage);
      });
}

[[nodiscard]] std::uint64_t exactDistributedIdAudit(
    const parallel::MpiContext& mpi_context,
    std::span<const std::uint64_t> local_ids,
    std::size_t batch_count,
    IcImportCounters& counters) {
  const std::size_t validated_batch_count = runCollectivePhase<std::size_t>(
      mpi_context, "IC global ID audit configuration", [&]() {
        if (batch_count == 0U) {
          throw std::invalid_argument(
              "IC ID audit batch count must be positive");
        }
        return batch_count;
      });

  struct TemporaryAuditDirectory {
    std::filesystem::path path;
    ~TemporaryAuditDirectory() {
      std::error_code error;
      std::filesystem::remove_all(path, error);
    }
  };

  const std::filesystem::path audit_directory =
      runCollectivePhase<std::filesystem::path>(
          mpi_context, "IC global ID audit temporary storage", [&]() {
            const auto nonce = std::chrono::high_resolution_clock::now()
                                   .time_since_epoch()
                                   .count();
            std::filesystem::path path =
                std::filesystem::temp_directory_path() /
                ("cosmosim_ic_id_audit_" +
                 std::to_string(mpi_context.worldRank()) + "_" +
                 std::to_string(nonce));
            if (!std::filesystem::create_directories(path)) {
              throw std::runtime_error(
                  "failed to create temporary distributed IC ID-audit "
                  "directory");
            }
            return path;
          });
  TemporaryAuditDirectory audit_cleanup{audit_directory};

  const std::uint64_t local_rounds =
      runCollectivePhase<std::uint64_t>(
          mpi_context, "IC global ID audit round accounting", [&]() {
            if (local_ids.empty()) {
              return std::uint64_t{0};
            }
            return 1U + static_cast<std::uint64_t>(
                (local_ids.size() - 1U) / validated_batch_count);
          });
  std::uint64_t rounds = 0U;
  mpiAllreduce(
      &local_rounds, &rounds, 1, MPI_UINT64_T, MPI_MAX, MPI_COMM_WORLD);
  counters.distributed_id_audit_round_count = rounds;

  std::vector<std::optional<std::filesystem::path>> sorted_run_levels;
  std::uint64_t next_run_index = 0U;
  std::uint64_t peak_bytes = 0U;

  const auto pathStorageBytes = [&]() {
    std::uint64_t bytes = checkedCapacityBytes(
        sorted_run_levels.capacity(),
        sizeof(std::optional<std::filesystem::path>),
        "ID-audit path table");
    for (const auto& path : sorted_run_levels) {
      if (path.has_value()) {
        if (path->native().capacity() >
            std::numeric_limits<std::uint64_t>::max() /
                sizeof(std::filesystem::path::value_type)) {
          throw std::overflow_error("ID-audit path storage overflow");
        }
        checkedCounterAdd(
            bytes,
            static_cast<std::uint64_t>(path->native().capacity()) *
                sizeof(std::filesystem::path::value_type),
            "ID-audit path storage");
      }
    }
    return bytes;
  };

  const auto writeSortedRun = [&](std::span<const std::uint64_t> ids) {
    const std::filesystem::path path =
        audit_directory / ("run_" + std::to_string(next_run_index++) + ".bin");
    std::ofstream output;
    output.rdbuf()->pubsetbuf(nullptr, 0);
    output.open(path, std::ios::binary | std::ios::trunc);
    if (!output) {
      throw std::runtime_error(
          "failed to create distributed IC ID-audit run");
    }
    for (const std::uint64_t id : ids) {
      output.write(
          reinterpret_cast<const char*>(&id),
          static_cast<std::streamsize>(sizeof(id)));
    }
    output.close();
    if (!output) {
      throw std::runtime_error(
          "failed to write distributed IC ID-audit run");
    }
    return path;
  };

  const auto mergeSortedRuns = [&](const std::filesystem::path& lhs_path,
                                   const std::filesystem::path& rhs_path) {
    const std::filesystem::path output_path =
        audit_directory / ("run_" + std::to_string(next_run_index++) + ".bin");
    std::ifstream lhs;
    std::ifstream rhs;
    std::ofstream output;
    lhs.rdbuf()->pubsetbuf(nullptr, 0);
    rhs.rdbuf()->pubsetbuf(nullptr, 0);
    output.rdbuf()->pubsetbuf(nullptr, 0);
    lhs.open(lhs_path, std::ios::binary);
    rhs.open(rhs_path, std::ios::binary);
    output.open(output_path, std::ios::binary | std::ios::trunc);
    if (!lhs || !rhs || !output) {
      throw std::runtime_error(
          "failed to open distributed IC ID-audit merge run");
    }

    std::uint64_t lhs_value = 0U;
    std::uint64_t rhs_value = 0U;
    bool has_lhs = static_cast<bool>(lhs.read(
        reinterpret_cast<char*>(&lhs_value),
        static_cast<std::streamsize>(sizeof(lhs_value))));
    bool has_rhs = static_cast<bool>(rhs.read(
        reinterpret_cast<char*>(&rhs_value),
        static_cast<std::streamsize>(sizeof(rhs_value))));
    std::optional<std::uint64_t> previous;
    bool duplicate = false;
    while (has_lhs || has_rhs) {
      const bool take_lhs = has_lhs && (!has_rhs || lhs_value <= rhs_value);
      const std::uint64_t value = take_lhs ? lhs_value : rhs_value;
      if (previous.has_value() && *previous == value) {
        duplicate = true;
        break;
      }
      previous = value;
      output.write(
          reinterpret_cast<const char*>(&value),
          static_cast<std::streamsize>(sizeof(value)));
      if (!output) {
        throw std::runtime_error(
            "failed to write distributed IC ID-audit merge run");
      }
      if (take_lhs) {
        has_lhs = static_cast<bool>(lhs.read(
            reinterpret_cast<char*>(&lhs_value),
            static_cast<std::streamsize>(sizeof(lhs_value))));
      } else {
        has_rhs = static_cast<bool>(rhs.read(
            reinterpret_cast<char*>(&rhs_value),
            static_cast<std::streamsize>(sizeof(rhs_value))));
      }
    }
    if (!lhs.eof() && lhs.fail()) {
      throw std::runtime_error(
          "failed to read distributed IC ID-audit left run");
    }
    if (!rhs.eof() && rhs.fail()) {
      throw std::runtime_error(
          "failed to read distributed IC ID-audit right run");
    }
    output.close();
    if (!duplicate && !output) {
      throw std::runtime_error(
          "failed to finalize distributed IC ID-audit merge run");
    }
    if (!duplicate) {
      std::error_code error;
      std::filesystem::remove(lhs_path, error);
      error.clear();
      std::filesystem::remove(rhs_path, error);
    }
    return std::pair{output_path, duplicate};
  };

  for (std::uint64_t round = 0U; round < rounds; ++round) {
    std::vector<std::vector<std::uint8_t>> buckets =
        runCollectivePhase<std::vector<std::vector<std::uint8_t>>>(
            mpi_context, "IC global ID audit bucket preparation", [&]() {
              const std::size_t begin = static_cast<std::size_t>(
                  std::min<std::uint64_t>(
                      round * validated_batch_count, local_ids.size()));
              const std::size_t end = std::min(
                  local_ids.size(), begin + validated_batch_count);
              std::vector<std::vector<std::uint8_t>> prepared(
                  static_cast<std::size_t>(mpi_context.worldSize()));
              for (std::size_t index = begin; index < end; ++index) {
                const std::uint64_t id = local_ids[index];
                const int validator = static_cast<int>(
                    (id ^ (id >> 33U) ^ (id << 11U)) %
                    static_cast<std::uint64_t>(mpi_context.worldSize()));
                internal::appendLe64(
                    prepared[static_cast<std::size_t>(validator)], id);
              }
              return prepared;
            });
    std::uint64_t sent = 0U;
    std::uint64_t received_bytes = 0U;
    std::uint64_t exchange_peak = 0U;
    const std::vector<std::uint8_t> received = alltoallBytes(
        mpi_context, buckets, sent, received_bytes, &exchange_peak);
    runCollectivePhaseVoid(
        mpi_context, "IC global ID audit exchange accounting", [&]() {
          injectIcTestFault(mpi_context, "global_id_audit_accounting");
          if (mpi_context.isRoot()) {
            checkedCounterAdd(
                counters.exact_audit_exchange_count, 1U,
                "exact_audit_exchange_count");
          }
          checkedCounterAdd(counters.bytes_sent, sent, "bytes_sent");
          checkedCounterAdd(
              counters.bytes_received, received_bytes, "bytes_received");
        });

    std::uint64_t round_vector_capacity = 0U;
    const int local_duplicate = runCollectivePhase<int>(
        mpi_context, "IC global ID audit sorted-run update", [&]() {
          if (received.size() % 8U != 0U) {
            throw std::runtime_error(
                "IC ID validation wire size is corrupt");
          }
          std::vector<std::uint64_t> round_ids;
          round_ids.reserve(received.size() / 8U);
          std::size_t offset = 0U;
          while (offset < received.size()) {
            round_ids.push_back(internal::readLe64(received, offset));
          }
          round_vector_capacity = checkedCapacityBytes(
              round_ids.capacity(), sizeof(std::uint64_t),
              "ID-audit round vector");
          std::sort(round_ids.begin(), round_ids.end());
          if (std::adjacent_find(round_ids.begin(), round_ids.end()) !=
              round_ids.end()) {
            return 1;
          }
          if (round_ids.empty()) {
            return 0;
          }

          std::filesystem::path current = writeSortedRun(round_ids);
          std::size_t level = 0U;
          for (;;) {
            if (sorted_run_levels.size() <= level) {
              sorted_run_levels.resize(level + 1U);
            }
            if (!sorted_run_levels[level].has_value()) {
              sorted_run_levels[level] = std::move(current);
              break;
            }
            const auto [merged_path, duplicate] = mergeSortedRuns(
                *sorted_run_levels[level], current);
            sorted_run_levels[level].reset();
            if (duplicate) {
              return 1;
            }
            current = merged_path;
            ++level;
          }
          return 0;
        });
    int any_duplicate = 0;
    mpiAllreduce(
        &local_duplicate, &any_duplicate, 1, MPI_INT, MPI_MAX,
        MPI_COMM_WORLD);
    if (any_duplicate != 0) {
      throw std::runtime_error(
          "distributed IC import contains duplicate particle IDs across "
          "files or ranks");
    }
    runCollectivePhaseVoid(
        mpi_context, "IC global ID audit capacity accounting", [&]() {
          injectIcTestFault(
              mpi_context, "duplicate_audit_capacity_accounting");
          std::uint64_t round_peak = nestedByteCapacity(buckets);
          checkedCounterAdd(
              round_peak, exchange_peak, "global ID audit exchange storage");
          checkedCounterAdd(
              round_peak, round_vector_capacity,
              "global ID audit round-vector storage");
          injectIcTestFault(mpi_context, "audit_path_storage_accounting");
          checkedCounterAdd(
              round_peak, pathStorageBytes(),
              "global ID audit path storage");
          injectIcTestFault(mpi_context, "audit_peak_staging_accounting");
          peak_bytes = std::max(peak_bytes, round_peak);
        });
  }

  const int local_cross_level_duplicate = runCollectivePhase<int>(
      mpi_context, "IC global ID audit final external merge", [&]() {
        std::optional<std::filesystem::path> merged;
        for (auto& level : sorted_run_levels) {
          if (!level.has_value()) {
            continue;
          }
          if (!merged.has_value()) {
            merged = std::move(*level);
            level.reset();
            continue;
          }
          const auto [merged_path, duplicate] =
              mergeSortedRuns(*merged, *level);
          level.reset();
          if (duplicate) {
            return 1;
          }
          merged = merged_path;
        }
        return 0;
      });
  int any_cross_level_duplicate = 0;
  mpiAllreduce(
      &local_cross_level_duplicate, &any_cross_level_duplicate, 1, MPI_INT,
      MPI_MAX, MPI_COMM_WORLD);
  if (any_cross_level_duplicate != 0) {
    throw std::runtime_error(
        "distributed IC import contains duplicate particle IDs across files "
        "or ranks");
  }
  return runCollectivePhase<std::uint64_t>(
      mpi_context, "IC global ID audit final capacity accounting", [&]() {
        injectIcTestFault(mpi_context, "audit_path_storage_accounting");
        return std::max(peak_bytes, pathStorageBytes());
      });
}

void validateDistributedTotals(
    const parallel::MpiContext& mpi_context,
    const core::SimulationState& state,
    const IcManifest& manifest,
    const std::array<double, 5>& local_source_mass) {
  struct LocalTotals {
    std::array<std::uint64_t, 5> counts{};
    std::array<double, 5> masses{};
  };
  const LocalTotals local = runCollectivePhase<LocalTotals>(
      mpi_context, "IC distributed local totals", [&]() {
        LocalTotals totals;
        std::array<double, 5> compensation{};
        for (std::size_t index = 0; index < state.particles.size(); ++index) {
          const std::uint32_t species =
              state.particle_sidecar.species_tag[index];
          if (species >= totals.counts.size()) {
            throw std::runtime_error(
                "invalid species tag after distributed IC routing");
          }
          ++totals.counts[species];
          const double value = state.particles.mass_code[index];
          if (!std::isfinite(value) || !(value > 0.0)) {
            throw std::runtime_error(
                "distributed IC final mass must be finite and positive");
          }
          const double corrected = value - compensation[species];
          const double updated = totals.masses[species] + corrected;
          compensation[species] =
              (updated - totals.masses[species]) - corrected;
          totals.masses[species] = updated;
        }
        return totals;
      });

  std::array<std::uint64_t, 5> global_counts{};
  std::array<double, 5> global_final_mass{};
  std::array<double, 5> global_source_mass{};
  mpiAllreduce(
      local.counts.data(), global_counts.data(), 5, MPI_UINT64_T, MPI_SUM,
      MPI_COMM_WORLD);
  mpiAllreduce(
      local.masses.data(), global_final_mass.data(), 5, MPI_DOUBLE, MPI_SUM,
      MPI_COMM_WORLD);
  mpiAllreduce(
      local_source_mass.data(), global_source_mass.data(), 5, MPI_DOUBLE,
      MPI_SUM, MPI_COMM_WORLD);

  const std::array<std::uint64_t, 5> expected_counts =
      runCollectivePhase<std::array<std::uint64_t, 5>>(
          mpi_context, "IC distributed expected totals", [&]() {
            std::array<std::uint64_t, 5> expected{};
            for (std::size_t type = 0; type < kParticleTypeCount; ++type) {
              if (manifest.num_part_total[type] == 0U) {
                continue;
              }
              const std::uint32_t species =
                  speciesTag(manifest.species_policy[type]);
              if (expected[species] >
                  std::numeric_limits<std::uint64_t>::max() -
                      manifest.num_part_total[type]) {
                throw std::overflow_error(
                    "distributed IC expected species count overflow");
              }
              expected[species] += manifest.num_part_total[type];
            }
            return expected;
          });
  if (global_counts != expected_counts) {
    throw std::runtime_error(
        "distributed IC global species counts do not match the manifest");
  }
  for (std::size_t species = 0; species < expected_counts.size(); ++species) {
    const double tolerance = 1.0e-12 * std::max(
        {1.0, std::abs(global_source_mass[species]),
         std::abs(global_final_mass[species])});
    if (std::abs(
            global_source_mass[species] - global_final_mass[species]) >
        tolerance) {
      throw std::runtime_error(
          "distributed IC global species mass total changed during routing");
    }
  }
}

#endif

}  // namespace cosmosim::io::distributed_audit_internal
