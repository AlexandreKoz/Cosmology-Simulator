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
#include <numeric>
#include <optional>
#include <sstream>
#include <streambuf>
#include <unordered_map>
#include <unordered_set>
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
LocalOwnershipIdentitySummary summarizeLocalOwnedParticleIds(std::span<const std::uint64_t> local_particle_ids) {
  LocalOwnershipIdentitySummary summary;
  summary.local_owned_count = static_cast<std::uint64_t>(local_particle_ids.size());
  std::unordered_set<std::uint64_t> seen;
  seen.reserve(local_particle_ids.size());
  for (const std::uint64_t particle_id : local_particle_ids) {
    summary.local_particle_id_sum += particle_id;
    summary.local_particle_id_square_sum += particle_id * particle_id;
    summary.local_particle_id_xor ^= particle_id;
    if (!seen.insert(particle_id).second) {
      summary.local_particle_ids_unique = false;
    }
  }
  return summary;
}

ExactOwnershipPartitionReport validateExactGlobalOwnershipPartition(
    const MpiContext& mpi_context,
    std::span<const std::uint64_t> local_owned_particle_ids,
    std::span<const std::uint64_t> expected_local_reference_particle_ids) {
  ExactOwnershipPartitionReport report;
  report.global_owned_count = mpi_context.allreduceSumUint64(
      static_cast<std::uint64_t>(local_owned_particle_ids.size()));

  const auto sorted_unique = [](std::span<const std::uint64_t> values) {
    std::vector<std::uint64_t> sorted(values.begin(), values.end());
    std::sort(sorted.begin(), sorted.end());
    const bool unique =
        std::adjacent_find(sorted.begin(), sorted.end()) == sorted.end();
    return std::pair{std::move(sorted), unique};
  };

  if (!mpi_context.isEnabled()) {
    auto [current, current_unique] = sorted_unique(local_owned_particle_ids);
    auto [expected, expected_unique] =
        sorted_unique(expected_local_reference_particle_ids);
    report.local_particle_ids_unique = current_unique;
    report.globally_unique = current_unique;
    for (auto it = current.begin(); it != current.end();) {
      const auto range = std::equal_range(it, current.end(), *it);
      if (std::distance(range.first, range.second) > 1) {
        report.duplicate_particle_ids.push_back(*it);
      }
      it = range.second;
    }
    current.erase(std::unique(current.begin(), current.end()), current.end());
    expected.erase(std::unique(expected.begin(), expected.end()), expected.end());
    std::set_difference(
        expected.begin(), expected.end(), current.begin(), current.end(),
        std::back_inserter(report.missing_expected_particle_ids));
    std::set_difference(
        current.begin(), current.end(), expected.begin(), expected.end(),
        std::back_inserter(report.extra_particle_ids));
    report.matches_expected_ids = expected_unique &&
        report.missing_expected_particle_ids.empty() &&
        report.extra_particle_ids.empty();
    return report;
  }

#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  const int world_size = mpi_context.worldSize();
  const auto exchange_by_hash = [&](std::span<const std::uint64_t> ids) {
    std::vector<std::vector<std::uint8_t>> send_payloads;
    std::exception_ptr local_preparation_failure;
    try {
      send_payloads.resize(static_cast<std::size_t>(world_size));
      for (const std::uint64_t id : ids) {
        const std::uint64_t mixed = id ^ (id >> 33U) ^ (id << 11U);
        const std::size_t owner = static_cast<std::size_t>(
            mixed % static_cast<std::uint64_t>(world_size));
        auto& payload = send_payloads[owner];
        const std::size_t old_size = payload.size();
        const std::size_t new_size = core::checkedSizeAdd(
            old_size, sizeof(id), "exact ownership hash bucket byte growth");
        payload.resize(new_size);
        std::memcpy(payload.data() + old_size, &id, sizeof(id));
      }
    } catch (...) {
      local_preparation_failure = std::current_exception();
    }
    mpi_context.rethrowCollectivePreparationFailure(
        local_preparation_failure,
        "exact ownership hash local payload preparation");

    std::vector<std::vector<std::uint8_t>> recv_payloads =
        exchangeBoundedAlltoallBytes(mpi_context, send_payloads);

    std::vector<std::uint64_t> recv;
    std::exception_ptr local_decode_failure;
    try {
      std::size_t total_ids = 0U;
      for (const auto& payload : recv_payloads) {
        if (payload.size() % sizeof(std::uint64_t) != 0U) {
          throw std::runtime_error(
              "exact ownership hash exchange returned partial uint64 record bytes");
        }
        total_ids = core::checkedSizeAdd(
            total_ids, payload.size() / sizeof(std::uint64_t),
            "exact ownership hash receive record total");
      }
      recv.resize(total_ids);
      std::size_t destination = 0U;
      for (const auto& payload : recv_payloads) {
        const std::size_t count = payload.size() / sizeof(std::uint64_t);
        if (!payload.empty()) {
          std::memcpy(
              recv.data() + destination, payload.data(), payload.size());
        }
        destination += count;
      }
    } catch (...) {
      local_decode_failure = std::current_exception();
    }
    mpi_context.rethrowCollectivePreparationFailure(
        local_decode_failure,
        "exact ownership hash receive reassembly");
    return recv;
  };

  std::vector<std::uint64_t> current_partition =
      exchange_by_hash(local_owned_particle_ids);
  std::vector<std::uint64_t> expected_partition =
      exchange_by_hash(expected_local_reference_particle_ids);
  std::sort(current_partition.begin(), current_partition.end());
  std::sort(expected_partition.begin(), expected_partition.end());

  for (auto it = current_partition.begin(); it != current_partition.end();) {
    const auto range = std::equal_range(it, current_partition.end(), *it);
    if (std::distance(range.first, range.second) > 1) {
      report.duplicate_particle_ids.push_back(*it);
    }
    it = range.second;
  }
  const bool local_unique = report.duplicate_particle_ids.empty();
  current_partition.erase(
      std::unique(current_partition.begin(), current_partition.end()),
      current_partition.end());
  expected_partition.erase(
      std::unique(expected_partition.begin(), expected_partition.end()),
      expected_partition.end());
  std::set_difference(
      expected_partition.begin(), expected_partition.end(),
      current_partition.begin(), current_partition.end(),
      std::back_inserter(report.missing_expected_particle_ids));
  std::set_difference(
      current_partition.begin(), current_partition.end(),
      expected_partition.begin(), expected_partition.end(),
      std::back_inserter(report.extra_particle_ids));

  const int local_duplicate = local_unique ? 0 : 1;
  const int local_mismatch =
      (report.missing_expected_particle_ids.empty() &&
       report.extra_particle_ids.empty()) ? 0 : 1;
  int any_duplicate = 0;
  int any_mismatch = 0;
  MPI_Allreduce(
      &local_duplicate, &any_duplicate, 1, MPI_INT, MPI_MAX,
      MPI_COMM_WORLD);
  MPI_Allreduce(
      &local_mismatch, &any_mismatch, 1, MPI_INT, MPI_MAX,
      MPI_COMM_WORLD);
  report.local_particle_ids_unique = any_duplicate == 0;
  report.globally_unique = any_duplicate == 0;
  report.matches_expected_ids = any_mismatch == 0;
  return report;
#else
  throw std::runtime_error(
      "exact global ownership validation requires MPI support when MPI context is enabled");
#endif
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
