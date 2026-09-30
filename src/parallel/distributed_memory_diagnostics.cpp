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

  constexpr std::uint64_t k_workspace = k_exact_ownership_validation_workspace_limit_bytes;
  constexpr std::uint64_t k_id_bytes = sizeof(std::uint64_t);
  constexpr std::uint64_t k_min_ids_per_bucket = 1024U;
  const std::uint64_t max_combined_ids = std::max<std::uint64_t>(
      k_min_ids_per_bucket, k_workspace / (4U * k_id_bytes));

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
  std::vector<BucketPrefix> pending{{}};
  std::uint64_t local_duplicate_count = 0U;
  std::uint64_t local_missing_count = 0U;
  std::uint64_t local_extra_count = 0U;

  const int world_size = mpi_context.worldSize();
  const auto matches_prefix = [&](std::uint64_t id, const BucketPrefix& prefix) {
    if (prefix.bits == 0U) return true;
    // The low hash bits are used by the hash-owner mapping. Partition on an
    // independent high-bit radix so each pass bounds every owner's partition
    // instead of accidentally selecting a single owner.
    const std::uint64_t mixed = mix_id(id) >> 32U;
    const std::uint64_t mask = (1ULL << prefix.bits) - 1ULL;
    return (mixed & mask) == prefix.value;
  };

  while (!pending.empty()) {
    const BucketPrefix prefix = pending.back();
    pending.pop_back();
    std::vector<std::uint64_t> current_count_by_owner(static_cast<std::size_t>(world_size), 0U);
    std::vector<std::uint64_t> expected_count_by_owner(static_cast<std::size_t>(world_size), 0U);
    for (const auto id : local_owned_particle_ids) {
      if (!matches_prefix(id, prefix)) continue;
      const std::size_t owner = static_cast<std::size_t>(mix_id(id) % static_cast<std::uint64_t>(world_size));
      ++current_count_by_owner[owner];
    }
    for (const auto id : expected_local_reference_particle_ids) {
      if (!matches_prefix(id, prefix)) continue;
      const std::size_t owner = static_cast<std::size_t>(mix_id(id) % static_cast<std::uint64_t>(world_size));
      ++expected_count_by_owner[owner];
    }
    mpi_context.allreduceSumUint64sInPlace(current_count_by_owner);
    mpi_context.allreduceSumUint64sInPlace(expected_count_by_owner);
    std::uint64_t max_owner_combined_count = 0U;
    for (std::size_t owner = 0; owner < current_count_by_owner.size(); ++owner) {
      max_owner_combined_count = std::max(
          max_owner_combined_count,
          core::checkedMemoryBytesAdd(current_count_by_owner[owner], expected_count_by_owner[owner],
                                      "exact ownership bucket owner count"));
    }
    if (max_owner_combined_count > max_combined_ids) {
      if (prefix.bits >= 31U) {
        throw std::runtime_error("exact ownership validator could not refine a bucket below the hard workspace bound");
      }
      const std::uint8_t child_bits = static_cast<std::uint8_t>(prefix.bits + 1U);
      pending.push_back(BucketPrefix{child_bits, prefix.value | (1ULL << prefix.bits)});
      pending.push_back(BucketPrefix{child_bits, prefix.value});
      continue;
    }

    auto collect_partition = [&](std::span<const std::uint64_t> ids) {
      std::vector<std::uint64_t> partition;
      if (!mpi_context.isEnabled()) {
        partition.reserve(core::checkedIntegralNarrow<std::size_t>(
            static_cast<std::uint64_t>(ids.size()), "exact ownership serial reserve"));
        for (const auto id : ids) if (matches_prefix(id, prefix)) partition.push_back(id);
        return partition;
      }
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
      std::vector<std::vector<std::uint8_t>> send_payloads(static_cast<std::size_t>(world_size));
      for (const auto id : ids) {
        if (!matches_prefix(id, prefix)) continue;
        const std::size_t owner = static_cast<std::size_t>(mix_id(id) % static_cast<std::uint64_t>(world_size));
        auto& payload = send_payloads[owner];
        const std::size_t old_size = payload.size();
        payload.resize(core::checkedSizeAdd(old_size, sizeof(id), "exact ownership bounded bucket payload"));
        std::memcpy(payload.data() + old_size, &id, sizeof(id));
      }
      const auto recv_payloads = exchangeBoundedAlltoallBytes(mpi_context, send_payloads);
      std::size_t count = 0U;
      for (const auto& payload : recv_payloads) {
        if ((payload.size() % sizeof(std::uint64_t)) != 0U) {
          throw std::runtime_error("exact ownership bounded bucket exchange returned partial ID bytes");
        }
        count = core::checkedSizeAdd(count, payload.size() / sizeof(std::uint64_t),
                                     "exact ownership bounded receive count");
      }
      const std::uint64_t bytes = core::checkedIntegralNarrow<std::uint64_t>(
          core::checkedSizeMultiply(count, sizeof(std::uint64_t), "exact ownership bounded receive bytes"),
          "exact ownership bounded receive byte width");
      if (bytes > k_workspace) {
        throw std::runtime_error("exact ownership bounded receive exceeded the hard workspace contract");
      }
      partition.resize(count);
      std::size_t offset = 0U;
      for (const auto& payload : recv_payloads) {
        if (!payload.empty()) std::memcpy(partition.data() + offset, payload.data(), payload.size());
        offset += payload.size() / sizeof(std::uint64_t);
      }
      return partition;
#else
      throw std::runtime_error("exact ownership validation requires MPI support when MPI context is enabled");
#endif
    };

    std::vector<std::uint64_t> current = collect_partition(local_owned_particle_ids);
    std::vector<std::uint64_t> expected = collect_partition(expected_local_reference_particle_ids);
    const std::uint64_t local_workspace_bytes = core::checkedMemoryBytesAdd(
        core::checkedIntegralNarrow<std::uint64_t>(current.capacity() * sizeof(std::uint64_t), "exact ownership current capacity"),
        core::checkedIntegralNarrow<std::uint64_t>(expected.capacity() * sizeof(std::uint64_t), "exact ownership expected capacity"),
        "exact ownership combined bucket capacity");
    if (local_workspace_bytes > k_workspace) {
      if (prefix.bits >= 63U) {
        throw std::runtime_error("exact ownership validator bucket capacity exceeded the hard workspace bound");
      }
      const std::uint8_t child_bits = static_cast<std::uint8_t>(prefix.bits + 1U);
      pending.push_back(BucketPrefix{child_bits, prefix.value | (1ULL << prefix.bits)});
      pending.push_back(BucketPrefix{child_bits, prefix.value});
      continue;
    }
    std::sort(current.begin(), current.end());
    std::sort(expected.begin(), expected.end());
    for (auto it = current.begin(); it != current.end();) {
      const auto range = std::equal_range(it, current.end(), *it);
      const auto multiplicity = static_cast<std::uint64_t>(std::distance(range.first, range.second));
      if (multiplicity > 1U) {
        local_duplicate_count += multiplicity - 1U;
        append_sample(report.duplicate_particle_ids, *it);
      }
      it = range.second;
    }
    current.erase(std::unique(current.begin(), current.end()), current.end());
    expected.erase(std::unique(expected.begin(), expected.end()), expected.end());
    std::size_t i = 0U;
    std::size_t j = 0U;
    while (i < current.size() || j < expected.size()) {
      if (j == expected.size() || (i < current.size() && current[i] < expected[j])) {
        ++local_extra_count;
        append_sample(report.extra_particle_ids, current[i++]);
      } else if (i == current.size() || expected[j] < current[i]) {
        ++local_missing_count;
        append_sample(report.missing_expected_particle_ids, expected[j++]);
      } else {
        ++i;
        ++j;
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
