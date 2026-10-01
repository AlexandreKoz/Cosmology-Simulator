#include "cosmosim/parallel/distributed_memory.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <iomanip>
#include <limits>
#include <memory>
#include <numeric>
#include <optional>
#include <sstream>
#include <streambuf>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/core/memory_governor.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "parallel/internal/distributed_memory_support.hpp"

#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
#include <mpi.h>
#endif

namespace cosmosim::parallel {
namespace {

[[nodiscard]] double absoluteValue(double value) {
  return value < 0.0 ? -value : value;
}

[[nodiscard]] double stableRelativeError(double measured, double reference, double absolute_error) {
  (void)measured;
  const double denom = std::max(absoluteValue(reference), std::numeric_limits<double>::min());
  return absolute_error / denom;
}

[[nodiscard]] std::string rankConfigValueString(std::uint64_t value) {
  return std::to_string(value);
}

[[nodiscard]] std::string rankConfigValueString(int value) {
  return std::to_string(value);
}

[[nodiscard]] std::string rankConfigValueString(bool value) {
  return value ? "true" : "false";
}

void appendRankConfigMismatch(
    RankConfigConsensus* consensus,
    RankConfigMismatchProperty property,
    int baseline_rank,
    int rank,
    std::string baseline_value,
    std::string rank_value) {
  consensus->mismatches.push_back(RankConfigMismatch{
      .property = property,
      .baseline_rank = baseline_rank,
      .rank = rank,
      .baseline_value = std::move(baseline_value),
      .rank_value = std::move(rank_value),
  });
}

}  // namespace

LocalOwnershipIdentitySummary summarizeLocalOwnedParticleIds(
    std::span<const std::uint64_t> local_particle_ids,
    bool local_particle_ids_unique) {
  LocalOwnershipIdentitySummary summary;
  summary.local_owned_count = static_cast<std::uint64_t>(local_particle_ids.size());
  summary.local_particle_ids_unique = local_particle_ids_unique;
  for (const std::uint64_t particle_id : local_particle_ids) {
    summary.local_particle_id_sum += particle_id;
    summary.local_particle_id_square_sum += particle_id * particle_id;
    summary.local_particle_id_xor ^= particle_id;
  }
  return summary;
}

LocalOwnershipIdentitySummary summarizeLocalOwnedParticleIds(
    std::span<const std::uint64_t> local_particle_ids) {
  core::OwnershipValidationWorkspace scratch;
  const core::ParticleIdValidationResult validation =
      core::validateParticleIdsExact(local_particle_ids, scratch);
  return summarizeLocalOwnedParticleIds(local_particle_ids, validation.unique);
}

ExactOwnershipPartitionReport validateExactGlobalOwnershipPartition(
    const MpiContext& mpi_context,
    std::span<const std::uint64_t> local_owned_particle_ids,
    std::span<const std::uint64_t> expected_local_reference_particle_ids,
    core::MemoryGovernor* memory_governor) {
  ExactOwnershipPartitionReport report;
  report.global_owned_count = mpi_context.allreduceSumUint64(
      static_cast<std::uint64_t>(local_owned_particle_ids.size()));

  constexpr std::uint64_t k_workspace = k_exact_ownership_validation_workspace_limit_bytes;
  constexpr std::uint64_t k_comparison_storage_bytes = 24ULL * 1024ULL * 1024ULL;
  constexpr std::uint64_t k_send_stage_bytes = 15ULL * 1024ULL * 1024ULL;
  constexpr std::uint64_t k_receive_stage_bytes = 15ULL * 1024ULL * 1024ULL;
  constexpr std::uint64_t k_metadata_bytes = 2ULL * 1024ULL * 1024ULL;
  constexpr std::uint64_t k_reserved_slack_bytes =
      k_workspace - k_comparison_storage_bytes - k_send_stage_bytes -
      k_receive_stage_bytes - k_metadata_bytes;
  static_assert(k_reserved_slack_bytes == 8ULL * 1024ULL * 1024ULL);

  core::MemoryReservation validation_reservation;
  std::exception_ptr local_admission_failure;
  try {
    if (memory_governor != nullptr) {
      validation_reservation = memory_governor->reserve(
          core::MemoryClass::kDiagnostic, k_workspace,
          "parallel.exact_ownership_validation");
      validation_reservation.commit();
    }
  } catch (...) {
    local_admission_failure = std::current_exception();
  }
  mpi_context.rethrowCollectivePreparationFailure(
      local_admission_failure, "exact ownership validation workspace admission");

  const int world_size = mpi_context.worldSize();
  const int world_rank = mpi_context.worldRank();
  if (world_size <= 0 || world_rank < 0 || world_rank >= world_size) {
    throw std::invalid_argument("exact ownership validation requires a valid MPI topology");
  }
  const std::size_t rank_count = static_cast<std::size_t>(world_size);
  const std::uint64_t metadata_need = core::checkedMemoryBytesAdd(
      static_cast<std::uint64_t>(rank_count) * 64ULL, 4096ULL,
      "exact ownership metadata allowance");
  if (metadata_need > k_metadata_bytes) {
    throw std::length_error("exact ownership rank metadata exceeds the fixed workspace allowance");
  }

  const std::size_t comparison_capacity_ids = static_cast<std::size_t>(
      k_comparison_storage_bytes / sizeof(std::uint64_t));
  const std::size_t stage_capacity_ids = static_cast<std::size_t>(
      std::min(k_send_stage_bytes, k_receive_stage_bytes) / sizeof(std::uint64_t));
  const std::size_t per_peer_chunk_ids = stage_capacity_ids / rank_count;
  if (mpi_context.isEnabled() &&
      (per_peer_chunk_ids == 0U ||
       per_peer_chunk_ids > static_cast<std::size_t>(std::numeric_limits<int>::max()))) {
    throw std::length_error(
        "exact ownership MPI topology leaves no bounded per-peer staging capacity");
  }

  std::unique_ptr<std::uint64_t[]> comparison_storage;
  std::unique_ptr<std::uint64_t[]> send_stage;
  std::unique_ptr<std::uint64_t[]> receive_stage;
  std::unique_ptr<std::uint64_t[]> local_current_count_by_owner;
  std::unique_ptr<std::uint64_t[]> local_expected_count_by_owner;
  std::unique_ptr<std::uint64_t[]> current_count_by_owner;
  std::unique_ptr<std::uint64_t[]> expected_count_by_owner;
  std::unique_ptr<std::uint64_t[]> ordinal_by_owner;
  std::unique_ptr<std::uint64_t[]> written_by_owner;
  std::unique_ptr<int[]> send_counts;
  std::unique_ptr<int[]> receive_counts;
  std::unique_ptr<int[]> send_displacements;
  std::unique_ptr<int[]> receive_displacements;
  std::exception_ptr local_workspace_failure;
  try {
    comparison_storage = std::make_unique<std::uint64_t[]>(comparison_capacity_ids);
    if (mpi_context.isEnabled()) {
      send_stage = std::make_unique<std::uint64_t[]>(stage_capacity_ids);
      receive_stage = std::make_unique<std::uint64_t[]>(stage_capacity_ids);
    }
    local_current_count_by_owner = std::make_unique<std::uint64_t[]>(rank_count);
    local_expected_count_by_owner = std::make_unique<std::uint64_t[]>(rank_count);
    current_count_by_owner = std::make_unique<std::uint64_t[]>(rank_count);
    expected_count_by_owner = std::make_unique<std::uint64_t[]>(rank_count);
    ordinal_by_owner = std::make_unique<std::uint64_t[]>(rank_count);
    written_by_owner = std::make_unique<std::uint64_t[]>(rank_count);
    send_counts = std::make_unique<int[]>(rank_count);
    receive_counts = std::make_unique<int[]>(rank_count);
    send_displacements = std::make_unique<int[]>(rank_count);
    receive_displacements = std::make_unique<int[]>(rank_count);
  } catch (...) {
    local_workspace_failure = std::current_exception();
  }
  mpi_context.rethrowCollectivePreparationFailure(
      local_workspace_failure, "exact ownership validation workspace allocation");

  const std::span<std::uint64_t> local_current_counts(
      local_current_count_by_owner.get(), rank_count);
  const std::span<std::uint64_t> local_expected_counts(
      local_expected_count_by_owner.get(), rank_count);
  const std::span<std::uint64_t> current_counts(
      current_count_by_owner.get(), rank_count);
  const std::span<std::uint64_t> expected_counts(
      expected_count_by_owner.get(), rank_count);
  const std::span<std::uint64_t> ordinals(ordinal_by_owner.get(), rank_count);
  const std::span<std::uint64_t> written(written_by_owner.get(), rank_count);

  const auto mix_id = [](std::uint64_t value) noexcept {
    value += 0x9e3779b97f4a7c15ULL;
    value = (value ^ (value >> 30U)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27U)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31U);
  };
  const auto append_sample = [](std::vector<std::uint64_t>& samples, std::uint64_t id) {
    if (samples.size() < ExactOwnershipPartitionReport::k_max_diagnostic_samples) {
      samples.push_back(id);
    }
  };

  struct BucketPrefix {
    std::uint8_t bits = 0U;
    std::uint64_t value = 0U;
  };
  std::array<BucketPrefix, 65U> pending{};
  std::size_t pending_size = 1U;
  std::uint64_t local_duplicate_count = 0U;
  std::uint64_t local_missing_count = 0U;
  std::uint64_t local_extra_count = 0U;

  const auto matches_prefix = [](std::uint64_t id, const BucketPrefix& prefix) noexcept {
    if (prefix.bits == 0U) return true;
    if (prefix.bits == 64U) return id == prefix.value;
    return (id & ((1ULL << prefix.bits) - 1ULL)) == prefix.value;
  };

  while (pending_size != 0U) {
    const BucketPrefix prefix = pending[--pending_size];
    std::fill(local_current_counts.begin(), local_current_counts.end(), 0U);
    std::fill(local_expected_counts.begin(), local_expected_counts.end(), 0U);
    for (const std::uint64_t id : local_owned_particle_ids) {
      if (!matches_prefix(id, prefix)) continue;
      ++local_current_counts[static_cast<std::size_t>(
          mix_id(id) % static_cast<std::uint64_t>(world_size))];
    }
    for (const std::uint64_t id : expected_local_reference_particle_ids) {
      if (!matches_prefix(id, prefix)) continue;
      ++local_expected_counts[static_cast<std::size_t>(
          mix_id(id) % static_cast<std::uint64_t>(world_size))];
    }
    std::copy(local_current_counts.begin(), local_current_counts.end(), current_counts.begin());
    std::copy(local_expected_counts.begin(), local_expected_counts.end(), expected_counts.begin());
    mpi_context.allreduceSumUint64sInPlace(current_counts);
    mpi_context.allreduceSumUint64sInPlace(expected_counts);

    std::uint64_t max_owner_combined_count = 0U;
    for (std::size_t owner = 0; owner < rank_count; ++owner) {
      max_owner_combined_count = std::max(
          max_owner_combined_count,
          core::checkedMemoryBytesAdd(
              current_counts[owner], expected_counts[owner],
              "exact ownership bucket owner item count"));
    }
    if (max_owner_combined_count > comparison_capacity_ids) {
      if (prefix.bits == 64U) {
        const std::uint64_t current_count = current_counts[static_cast<std::size_t>(world_rank)];
        const std::uint64_t expected_count = expected_counts[static_cast<std::size_t>(world_rank)];
        if (current_count > 1U) {
          local_duplicate_count += current_count - 1U;
          append_sample(report.duplicate_particle_ids, prefix.value);
        }
        if (current_count == 0U && expected_count != 0U) {
          ++local_missing_count;
          append_sample(report.missing_expected_particle_ids, prefix.value);
        } else if (current_count != 0U && expected_count == 0U) {
          ++local_extra_count;
          append_sample(report.extra_particle_ids, prefix.value);
        }
        continue;
      }
      if (pending_size + 2U > pending.size()) {
        throw std::logic_error("exact ownership prefix stack overflow");
      }
      const std::uint8_t child_bits = static_cast<std::uint8_t>(prefix.bits + 1U);
      pending[pending_size++] = BucketPrefix{
          child_bits, prefix.value | (1ULL << prefix.bits)};
      pending[pending_size++] = BucketPrefix{child_bits, prefix.value};
      continue;
    }

    const std::size_t current_size = core::checkedIntegralNarrow<std::size_t>(
        current_counts[static_cast<std::size_t>(world_rank)],
        "exact ownership current comparison count");
    const std::size_t expected_size = core::checkedIntegralNarrow<std::size_t>(
        expected_counts[static_cast<std::size_t>(world_rank)],
        "exact ownership expected comparison count");
    if (current_size > comparison_capacity_ids ||
        expected_size > comparison_capacity_ids - current_size) {
      throw std::logic_error(
          "exact ownership prefix planning exceeded fixed comparison storage");
    }
    std::span<std::uint64_t> current(
        comparison_storage.get(), current_size);
    std::span<std::uint64_t> expected(
        comparison_storage.get() + current_size, expected_size);

    const auto collect_partition = [&](std::span<const std::uint64_t> ids,
                                       std::span<const std::uint64_t> local_count_by_owner,
                                       std::span<std::uint64_t> partition) {
      std::size_t partition_write = 0U;
      if (!mpi_context.isEnabled()) {
        for (const std::uint64_t id : ids) {
          if (matches_prefix(id, prefix)) {
            if (partition_write >= partition.size()) {
              throw std::logic_error(
                  "exact ownership serial collection exceeded counted partition");
            }
            partition[partition_write++] = id;
          }
        }
        if (partition_write != partition.size()) {
          throw std::logic_error("exact ownership serial count/collection mismatch");
        }
        return;
      }
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
      std::uint64_t local_round_count = 0U;
      for (const std::uint64_t count : local_count_by_owner) {
        local_round_count = std::max(
            local_round_count,
            (count + static_cast<std::uint64_t>(per_peer_chunk_ids) - 1U) /
                static_cast<std::uint64_t>(per_peer_chunk_ids));
      }
      const std::uint64_t global_round_count =
          mpi_context.allreduceMaxUint64(local_round_count);

      for (std::uint64_t round = 0U; round < global_round_count; ++round) {
        std::exception_ptr local_round_plan_failure;
        int send_total = 0;
        const std::uint64_t round_begin =
            round * static_cast<std::uint64_t>(per_peer_chunk_ids);
        try {
          for (std::size_t peer = 0; peer < rank_count; ++peer) {
            const std::uint64_t remaining = local_count_by_owner[peer] > round_begin
                ? local_count_by_owner[peer] - round_begin : 0U;
            const std::uint64_t count = std::min<std::uint64_t>(
                remaining, static_cast<std::uint64_t>(per_peer_chunk_ids));
            send_counts[peer] = core::checkedIntegralNarrow<int>(
                count, "exact ownership round send count");
            send_displacements[peer] = send_total;
            send_total = core::checkedIntegralNarrow<int>(
                static_cast<std::uint64_t>(send_total) + count,
                "exact ownership round send displacement");
          }
          if (static_cast<std::size_t>(send_total) > stage_capacity_ids) {
            throw std::logic_error(
                "exact ownership send plan exceeded fixed stage storage");
          }
        } catch (...) {
          local_round_plan_failure = std::current_exception();
        }
        mpi_context.rethrowCollectivePreparationFailure(
            local_round_plan_failure, "exact ownership streaming send plan");

        if (MPI_Alltoall(
                send_counts.get(), 1, MPI_INT,
                receive_counts.get(), 1, MPI_INT,
                MPI_COMM_WORLD) != MPI_SUCCESS) {
          throw std::runtime_error("exact ownership streaming count exchange failed");
        }

        std::exception_ptr local_payload_prep_failure;
        int receive_total = 0;
        try {
          for (std::size_t peer = 0; peer < rank_count; ++peer) {
            receive_displacements[peer] = receive_total;
            if (receive_counts[peer] < 0) {
              throw std::runtime_error(
                  "exact ownership streaming received a negative count");
            }
            receive_total = core::checkedIntegralNarrow<int>(
                static_cast<std::uint64_t>(receive_total) +
                    static_cast<std::uint64_t>(receive_counts[peer]),
                "exact ownership round receive displacement");
          }
          if (static_cast<std::size_t>(receive_total) > stage_capacity_ids) {
            throw std::logic_error(
                "exact ownership receive plan exceeded fixed stage storage");
          }

          std::fill(ordinals.begin(), ordinals.end(), 0U);
          std::fill(written.begin(), written.end(), 0U);
          for (const std::uint64_t id : ids) {
            if (!matches_prefix(id, prefix)) continue;
            const std::size_t owner = static_cast<std::size_t>(
                mix_id(id) % static_cast<std::uint64_t>(world_size));
            const std::uint64_t ordinal = ordinals[owner]++;
            if (ordinal < round_begin ||
                ordinal >= round_begin +
                    static_cast<std::uint64_t>(send_counts[owner])) {
              continue;
            }
            const std::size_t slot =
                static_cast<std::size_t>(send_displacements[owner]) +
                static_cast<std::size_t>(written[owner]++);
            if (slot >= static_cast<std::size_t>(send_total)) {
              throw std::logic_error(
                  "exact ownership send staging exceeded planned displacement");
            }
            send_stage[slot] = id;
          }
          for (std::size_t peer = 0; peer < rank_count; ++peer) {
            if (written[peer] != static_cast<std::uint64_t>(send_counts[peer])) {
              throw std::logic_error(
                  "exact ownership send staging count disagrees with prefix count");
            }
          }
        } catch (...) {
          local_payload_prep_failure = std::current_exception();
        }
        mpi_context.rethrowCollectivePreparationFailure(
            local_payload_prep_failure,
            "exact ownership streaming payload preparation");

        if (MPI_Alltoallv(
                send_total == 0 ? nullptr : send_stage.get(),
                send_counts.get(), send_displacements.get(), MPI_UINT64_T,
                receive_total == 0 ? nullptr : receive_stage.get(),
                receive_counts.get(), receive_displacements.get(), MPI_UINT64_T,
                MPI_COMM_WORLD) != MPI_SUCCESS) {
          throw std::runtime_error("exact ownership streaming payload exchange failed");
        }

        std::exception_ptr local_consume_failure;
        try {
          if (partition_write + static_cast<std::size_t>(receive_total) >
              partition.size()) {
            throw std::logic_error(
                "exact ownership streaming receive exceeded counted destination size");
          }
          std::copy_n(
              receive_stage.get(), static_cast<std::size_t>(receive_total),
              partition.begin() + static_cast<std::ptrdiff_t>(partition_write));
          partition_write += static_cast<std::size_t>(receive_total);
        } catch (...) {
          local_consume_failure = std::current_exception();
        }
        mpi_context.rethrowCollectivePreparationFailure(
            local_consume_failure, "exact ownership streaming payload consume");
      }

      std::exception_ptr local_count_failure;
      if (partition_write != partition.size()) {
        local_count_failure = std::make_exception_ptr(std::logic_error(
            "exact ownership streaming receive count mismatch"));
      }
      mpi_context.rethrowCollectivePreparationFailure(
          local_count_failure, "exact ownership streaming final count");
#else
      (void)local_count_by_owner;
      (void)partition;
      throw std::runtime_error(
          "exact ownership validation requires MPI support when MPI context is enabled");
#endif
    };

    collect_partition(
        local_owned_particle_ids, local_current_counts, current);
    collect_partition(
        expected_local_reference_particle_ids, local_expected_counts, expected);

    std::sort(current.begin(), current.end());
    std::sort(expected.begin(), expected.end());
    for (std::size_t i = 0U; i < current.size();) {
      std::size_t next = i + 1U;
      while (next < current.size() && current[next] == current[i]) ++next;
      const std::uint64_t multiplicity = static_cast<std::uint64_t>(next - i);
      if (multiplicity > 1U) {
        local_duplicate_count += multiplicity - 1U;
        append_sample(report.duplicate_particle_ids, current[i]);
      }
      i = next;
    }

    std::size_t i = 0U;
    std::size_t j = 0U;
    while (i < current.size() || j < expected.size()) {
      if (j == expected.size() ||
          (i < current.size() && current[i] < expected[j])) {
        const std::uint64_t id = current[i];
        ++local_extra_count;
        append_sample(report.extra_particle_ids, id);
        while (i < current.size() && current[i] == id) ++i;
      } else if (i == current.size() || expected[j] < current[i]) {
        const std::uint64_t id = expected[j];
        ++local_missing_count;
        append_sample(report.missing_expected_particle_ids, id);
        while (j < expected.size() && expected[j] == id) ++j;
      } else {
        const std::uint64_t id = current[i];
        while (i < current.size() && current[i] == id) ++i;
        while (j < expected.size() && expected[j] == id) ++j;
      }
    }
  }

  report.duplicate_count = mpi_context.allreduceSumUint64(local_duplicate_count);
  report.missing_count = mpi_context.allreduceSumUint64(local_missing_count);
  report.extra_count = mpi_context.allreduceSumUint64(local_extra_count);
  report.local_particle_ids_unique = report.duplicate_count == 0U;
  report.globally_unique = report.duplicate_count == 0U;
  report.matches_expected_ids = report.missing_count == 0U && report.extra_count == 0U;
  return report;
}

bool partitionIdentityMatchesGeneratedSet(
    const LocalOwnershipIdentitySummary& reduced_global_summary,
    std::uint64_t expected_global_count,
    std::uint64_t expected_particle_id_sum,
    std::uint64_t expected_particle_id_square_sum,
    std::uint64_t expected_particle_id_xor) {
  return reduced_global_summary.local_particle_ids_unique &&
      reduced_global_summary.local_owned_count == expected_global_count &&
      reduced_global_summary.local_particle_id_sum == expected_particle_id_sum &&
      reduced_global_summary.local_particle_id_square_sum == expected_particle_id_square_sum &&
      reduced_global_summary.local_particle_id_xor == expected_particle_id_xor;
}

bool partitionIdentityMatchesGeneratedSet(
    const LocalOwnershipIdentitySummary& reduced_global_summary,
    std::uint64_t expected_global_count,
    std::uint64_t expected_particle_id_sum,
    std::uint64_t expected_particle_id_xor) {
  return partitionIdentityMatchesGeneratedSet(
      reduced_global_summary,
      expected_global_count,
      expected_particle_id_sum,
      reduced_global_summary.local_particle_id_square_sum,
      expected_particle_id_xor);
}

double deterministicRankOrderedSum(std::span<const double> per_rank_values) {
  long double sum = 0.0;
  for (const double value : per_rank_values) {
    sum += static_cast<long double>(value);
  }
  return static_cast<double>(sum);
}

ReductionAgreement compareReductionAgreement(
    std::span<const double> per_rank_values,
    double measured_sum) {
  const double deterministic_baseline_sum = deterministicRankOrderedSum(per_rank_values);
  const double absolute_error = absoluteValue(measured_sum - deterministic_baseline_sum);
  return ReductionAgreement{
      .deterministic_baseline_sum = deterministic_baseline_sum,
      .measured_sum = measured_sum,
      .absolute_error = absolute_error,
      .relative_error = stableRelativeError(measured_sum, deterministic_baseline_sum, absolute_error),
  };
}

bool satisfiesReductionAgreement(
    const ReductionAgreement& agreement,
    const ReductionAgreementPolicy& policy) {
  if (policy.absolute_tolerance < 0.0 || policy.relative_tolerance < 0.0) {
    throw std::invalid_argument("reduction agreement tolerances must be non-negative");
  }
  const bool absolute_ok = agreement.absolute_error <= policy.absolute_tolerance;
  const bool relative_ok = agreement.relative_error <= policy.relative_tolerance;
  switch (policy.mode) {
    case ReductionAgreementMode::kAbsoluteOnly:
      return absolute_ok;
    case ReductionAgreementMode::kRelativeOnly:
      return relative_ok;
    case ReductionAgreementMode::kAbsoluteAndRelative:
      return absolute_ok && relative_ok;
    case ReductionAgreementMode::kAbsoluteOrRelative:
      return absolute_ok || relative_ok;
  }
  throw std::invalid_argument("unknown reduction agreement mode");
}

bool RankConfigConsensus::allConsistent() const noexcept {
  return normalized_config_hash_match && mpi_ranks_expected_match && deterministic_reduction_match;
}

RankConfigConsensus evaluateRankConfigConsensus(std::span<const RankConfigDigest> digests) {
  RankConfigConsensus consensus;
  if (digests.empty()) {
    return consensus;
  }

  const RankConfigDigest baseline = digests.front();
  for (const RankConfigDigest& digest : digests) {
    bool rank_matches = true;
    if (digest.normalized_config_hash != baseline.normalized_config_hash) {
      consensus.normalized_config_hash_match = false;
      rank_matches = false;
      appendRankConfigMismatch(
          &consensus,
          RankConfigMismatchProperty::kNormalizedConfigHash,
          baseline.world_rank,
          digest.world_rank,
          rankConfigValueString(baseline.normalized_config_hash),
          rankConfigValueString(digest.normalized_config_hash));
    }
    if (digest.mpi_ranks_expected != baseline.mpi_ranks_expected) {
      consensus.mpi_ranks_expected_match = false;
      rank_matches = false;
      appendRankConfigMismatch(
          &consensus,
          RankConfigMismatchProperty::kMpiRanksExpected,
          baseline.world_rank,
          digest.world_rank,
          rankConfigValueString(baseline.mpi_ranks_expected),
          rankConfigValueString(digest.mpi_ranks_expected));
    }
    if (digest.deterministic_reduction != baseline.deterministic_reduction) {
      consensus.deterministic_reduction_match = false;
      rank_matches = false;
      appendRankConfigMismatch(
          &consensus,
          RankConfigMismatchProperty::kDeterministicReduction,
          baseline.world_rank,
          digest.world_rank,
          rankConfigValueString(baseline.deterministic_reduction),
          rankConfigValueString(digest.deterministic_reduction));
    }
    if (!rank_matches) {
      consensus.mismatched_ranks.push_back(digest.world_rank);
    }
  }
  return consensus;
}


}  // namespace cosmosim::parallel
