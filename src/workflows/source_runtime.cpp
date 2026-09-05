#include "cosmosim/workflows/source_runtime.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <exception>
#include <limits>
#include <memory>
#include <numeric>
#include <span>
#include <stdexcept>
#include <unordered_set>
#include <vector>

#include "cosmosim/core/checked_arithmetic.hpp"
#include "cosmosim/core/memory_governor.hpp"
#include "cosmosim/core/memory_accounting.hpp"
#include "cosmosim/physics/black_hole_agn.hpp"
#include "cosmosim/physics/effective_multiphase_ism.hpp"
#include "cosmosim/parallel/distributed_memory.hpp"
#include "cosmosim/physics/star_formation.hpp"
#include "cosmosim/physics/stellar_evolution.hpp"
#include "cosmosim/physics/stellar_feedback.hpp"
#include "cosmosim/physics/metal_diffusion.hpp"
#include "cosmosim/workflows/runtime_services.hpp"
#include "workflows/internal/gas_cell_ownership.hpp"
#include "workflows/internal/metal_diffusion_topology.hpp"
#include "workflows/internal/runtime_stage_resource_access.hpp"
#include "workflows/internal/star_formation_geometry.hpp"

#if COSMOSIM_ENABLE_MPI
#include <mpi.h>
#endif

namespace cosmosim::workflows {
namespace {

[[nodiscard]] double newtonGCodeFromUnits(const core::UnitSystem& units) {
  return core::newtonGravitationalConstantCode(units);
}

[[nodiscard]] physics::StarFormationConfig makeRuntimeStarFormationConfig(
    const core::SimulationConfig& config,
    const core::UnitSystem& units) {
  physics::StarFormationConfig result = physics::makeStarFormationConfig(config.physics);
  result.newton_g_code = newtonGCodeFromUnits(units);
  result.density_is_comoving = config.units.coordinate_frame == core::CoordinateFrame::kComoving;
  result.geometry_is_comoving = config.units.coordinate_frame == core::CoordinateFrame::kComoving;
  return result;
}

[[nodiscard]] physics::StellarFeedbackConfig makeRuntimeStellarFeedbackConfig(
    const core::SimulationConfig& config,
    const core::UnitSystem& units) {
  physics::StellarFeedbackConfig result =
      physics::makeStellarFeedbackConfig(config.physics);
  // Mass and metal return is part of stellar evolution even when energetic
  // feedback is disabled. Keep the deposition transaction enabled and zero only
  // the energy/momentum coupling in that case.
  result.enabled = config.physics.enable_stellar_evolution;
  if (!config.physics.enable_feedback) {
    result.epsilon_thermal = 0.0;
    result.epsilon_kinetic = 0.0;
    result.epsilon_momentum = 0.0;
  }
  constexpr double k_joule_per_erg = 1.0e-7;
  result.total_energy_code_per_erg = k_joule_per_erg /
      (units.mass_si_per_code * units.velocity_si_per_code *
       units.velocity_si_per_code);
  return result;
}

[[nodiscard]] std::shared_ptr<const physics::EffectiveMultiphaseEosTable>
makeRuntimeEffectiveEosTable(
    const core::SimulationConfig& config,
    const core::UnitSystem& units) {
  if (!config.physics.enable_star_formation ||
      config.physics.star_formation_model !=
          core::StarFormationModelKind::kEffectiveMultiphaseTngLike) {
    return nullptr;
  }
  return std::make_shared<const physics::EffectiveMultiphaseEosTable>(
      physics::makeEffectiveMultiphaseEosConfig(config.physics),
      units,
      physics::makeEffectiveIsmReferenceCoolingProvider(config.physics));
}

struct ShardedUint64Exchange {
  std::vector<std::uint64_t> values;
  std::vector<std::size_t> recv_counts;
  std::vector<std::size_t> recv_displacements;
};

[[nodiscard]] ShardedUint64Exchange exchangeShardedUint64Records(
    const parallel::MpiContext& mpi_context,
    const std::vector<std::vector<std::uint64_t>>& values_by_rank,
    std::size_t record_width) {
  const int world_size = mpi_context.worldSize();
  std::exception_ptr local_failure;
  try {
    if (record_width == 0U) {
      throw std::invalid_argument("sharded uint64 exchange requires nonzero record width");
    }
    if (world_size <= 0 || values_by_rank.size() != static_cast<std::size_t>(world_size)) {
      throw std::invalid_argument("sharded uint64 exchange rank extent mismatch");
    }
    for (const auto& values : values_by_rank) {
      if (values.size() % record_width != 0U) {
        throw std::invalid_argument("sharded uint64 exchange received a partial record");
      }
    }
  } catch (...) {
    local_failure = std::current_exception();
  }
  mpi_context.rethrowCollectivePreparationFailure(
      local_failure, "sharded uint64 local framing validation");

  if (!mpi_context.isEnabled()) {
    return ShardedUint64Exchange{
        .values = values_by_rank.front(),
        .recv_counts = {values_by_rank.front().size()},
        .recv_displacements = {0U},
    };
  }
#if COSMOSIM_ENABLE_MPI
  std::vector<std::uint64_t> send_counts64;
  std::vector<std::uint64_t> recv_counts64;
  local_failure = nullptr;
  try {
    send_counts64.resize(static_cast<std::size_t>(world_size), 0U);
    recv_counts64.resize(static_cast<std::size_t>(world_size), 0U);
    for (std::size_t rank = 0; rank < values_by_rank.size(); ++rank) {
      send_counts64[rank] = core::checkedIntegralNarrow<std::uint64_t>(
          values_by_rank[rank].size(), "sharded uint64 logical send count");
    }
  } catch (...) {
    local_failure = std::current_exception();
  }
  mpi_context.rethrowCollectivePreparationFailure(
      local_failure, "sharded uint64 count-buffer preparation");

  if (MPI_Alltoall(
          send_counts64.data(), 1, MPI_UINT64_T,
          recv_counts64.data(), 1, MPI_UINT64_T,
          MPI_COMM_WORLD) != MPI_SUCCESS) {
    throw std::runtime_error("sharded uint64 MPI_Alltoall count exchange failed");
  }

  parallel::BoundedMpiTransferPlan send_plan;
  parallel::BoundedMpiTransferPlan recv_plan;
  std::vector<std::uint64_t> recv_values;
  std::vector<std::uint64_t> send_round_buffer;
  std::vector<std::uint64_t> recv_round_buffer;
  parallel::BoundedMpiRoundLayout zero_round;
  local_failure = nullptr;
  try {
    std::vector<std::size_t> send_counts(send_counts64.size(), 0U);
    std::vector<std::size_t> recv_counts(recv_counts64.size(), 0U);
    for (std::size_t rank = 0; rank < send_counts.size(); ++rank) {
      send_counts[rank] = core::checkedIntegralNarrow<std::size_t>(
          send_counts64[rank], "sharded uint64 send count");
      recv_counts[rank] = core::checkedIntegralNarrow<std::size_t>(
          recv_counts64[rank], "sharded uint64 receive count");
      if (recv_counts[rank] % record_width != 0U) {
        throw std::runtime_error("sharded uint64 receive count violates record framing");
      }
    }
    const std::size_t round_bytes = parallel::mpiTransportRoundLimitBytes();
    const std::size_t round_elements = round_bytes / sizeof(std::uint64_t);
    if (round_elements == 0U) {
      throw std::runtime_error("sharded uint64 MPI transport round is smaller than one uint64 element");
    }
    send_plan = parallel::planBoundedMpiTransferRounds(
        send_counts,
        static_cast<std::size_t>(std::numeric_limits<int>::max()),
        round_elements);
    recv_plan = parallel::planBoundedMpiTransferRounds(
        recv_counts,
        static_cast<std::size_t>(std::numeric_limits<int>::max()),
        round_elements);
    recv_values.resize(recv_plan.logical_total_count, 0U);

    std::size_t maximum_send_round = 0U;
    for (const auto& round : send_plan.rounds) {
      maximum_send_round = std::max(maximum_send_round, round.round_count);
    }
    std::size_t maximum_recv_round = 0U;
    for (const auto& round : recv_plan.rounds) {
      maximum_recv_round = std::max(maximum_recv_round, round.round_count);
    }
    send_round_buffer.resize(maximum_send_round, 0U);
    recv_round_buffer.resize(maximum_recv_round, 0U);
    zero_round.counts.resize(static_cast<std::size_t>(world_size), 0);
    zero_round.displacements.resize(static_cast<std::size_t>(world_size), 0);
    zero_round.logical_offsets.resize(static_cast<std::size_t>(world_size), 0U);
  } catch (...) {
    local_failure = std::current_exception();
  }
  mpi_context.rethrowCollectivePreparationFailure(
      local_failure, "sharded uint64 post-count payload preparation");

  const std::uint64_t local_round_count = static_cast<std::uint64_t>(
      std::max(send_plan.rounds.size(), recv_plan.rounds.size()));
  const std::uint64_t global_round_count =
      mpi_context.allreduceMaxUint64(local_round_count);
  for (std::uint64_t round_index = 0U; round_index < global_round_count;
       ++round_index) {
    const auto& send_round = round_index < send_plan.rounds.size()
        ? send_plan.rounds[static_cast<std::size_t>(round_index)]
        : zero_round;
    const auto& recv_round = round_index < recv_plan.rounds.size()
        ? recv_plan.rounds[static_cast<std::size_t>(round_index)]
        : zero_round;
    for (std::size_t peer = 0; peer < send_round.counts.size(); ++peer) {
      const std::size_t count = static_cast<std::size_t>(send_round.counts[peer]);
      if (count == 0U) {
        continue;
      }
      const std::size_t copy_bytes = count * sizeof(std::uint64_t);
      std::memcpy(
          send_round_buffer.data() +
              static_cast<std::size_t>(send_round.displacements[peer]),
          values_by_rank[peer].data() + send_round.logical_offsets[peer],
          copy_bytes);
    }
    if (MPI_Alltoallv(
            send_round_buffer.empty() ? nullptr : send_round_buffer.data(),
            send_round.counts.data(), send_round.displacements.data(), MPI_UINT64_T,
            recv_round_buffer.empty() ? nullptr : recv_round_buffer.data(),
            recv_round.counts.data(), recv_round.displacements.data(), MPI_UINT64_T,
            MPI_COMM_WORLD) != MPI_SUCCESS) {
      throw std::runtime_error("sharded uint64 bounded MPI_Alltoallv payload exchange failed");
    }
    for (std::size_t peer = 0; peer < recv_round.counts.size(); ++peer) {
      const std::size_t count = static_cast<std::size_t>(recv_round.counts[peer]);
      if (count == 0U) {
        continue;
      }
      // The bounded plans already proved these prefixes and round extents.
      const std::size_t destination =
          recv_plan.logical_displacements[peer] +
          recv_round.logical_offsets[peer];
      const std::size_t copy_bytes = count * sizeof(std::uint64_t);
      std::memcpy(
          recv_values.data() + destination,
          recv_round_buffer.data() +
              static_cast<std::size_t>(recv_round.displacements[peer]),
          copy_bytes);
    }
  }
  return ShardedUint64Exchange{
      .values = std::move(recv_values),
      .recv_counts = std::move(recv_plan.logical_counts),
      .recv_displacements = std::move(recv_plan.logical_displacements),
  };
#else
  throw std::runtime_error("distributed sharded particle-ID registry requires an MPI build");
#endif
}

[[nodiscard]] std::size_t shardForUint64(std::uint64_t value, int world_size) {
  if (world_size <= 0) {
    throw std::invalid_argument("particle-ID sharding requires positive world size");
  }
  // Multiplicative mixing keeps sequential imported IDs from concentrating on
  // a subset of shards while remaining deterministic for any rank count.
  value ^= value >> 33U;
  value *= 0xff51afd7ed558ccdULL;
  value ^= value >> 33U;
  value *= 0xc4ceb9fe1a85ec53ULL;
  value ^= value >> 33U;
  return static_cast<std::size_t>(value % static_cast<std::uint64_t>(world_size));
}

class DistributedParticleIdRegistry final : public physics::ParticleIdPrecommit {
 public:
  explicit DistributedParticleIdRegistry(const parallel::MpiContext& mpi_context)
      : m_mpi_context(mpi_context) {}

  [[nodiscard]] std::vector<std::uint64_t> precommit(
      const core::SimulationState& state,
      std::span<const std::uint64_t> birth_keys) override {
    if (!m_initialized) {
      initializeOccupiedShard(state.particle_sidecar.particle_id);
    }
    validateBirthKeysDistributed(birth_keys);

    struct LocalCandidate {
      std::uint64_t particle_id = 0U;
      std::uint64_t birth_key = 0U;
      std::uint32_t collision_ordinal = 0U;
      bool resolved = false;
    };
    std::vector<LocalCandidate> candidates;
    candidates.reserve(birth_keys.size());
    for (const std::uint64_t birth_key : birth_keys) {
      candidates.push_back(LocalCandidate{
          .particle_id = physics::starFormationParticleIdFromBirthKey(birth_key, 0U),
          .birth_key = birth_key,
      });
    }

    std::unordered_set<std::uint64_t> batch_reserved_on_shard;
    const int world_size = m_mpi_context.worldSize();
    for (std::uint32_t pass = 0U; pass < 1024U; ++pass) {
      std::vector<std::vector<std::uint64_t>> requests(
          static_cast<std::size_t>(world_size));
      for (std::size_t local_index = 0; local_index < candidates.size(); ++local_index) {
        const LocalCandidate& candidate = candidates[local_index];
        if (candidate.resolved) {
          continue;
        }
        auto& target = requests[shardForUint64(candidate.particle_id, world_size)];
        target.push_back(candidate.particle_id);
        target.push_back(candidate.birth_key);
        target.push_back(static_cast<std::uint64_t>(local_index));
        target.push_back(static_cast<std::uint64_t>(candidate.collision_ordinal));
      }
      const ShardedUint64Exchange received =
          exchangeShardedUint64Records(m_mpi_context, requests, 4U);

      struct ReceivedCandidate {
        std::uint64_t particle_id = 0U;
        std::uint64_t birth_key = 0U;
        std::uint64_t origin_index = 0U;
        std::uint32_t ordinal = 0U;
        int origin_rank = 0;
      };
      std::vector<ReceivedCandidate> shard_candidates;
      shard_candidates.reserve(received.values.size() / 4U);
      for (int origin_rank = 0; origin_rank < world_size; ++origin_rank) {
        const std::size_t begin = received.recv_displacements[static_cast<std::size_t>(origin_rank)];
        const std::size_t count = received.recv_counts[static_cast<std::size_t>(origin_rank)];
        for (std::size_t offset = 0; offset < count; offset += 4U) {
          const std::size_t index = begin + offset;
          const std::uint64_t ordinal64 = received.values[index + 3U];
          if (ordinal64 > std::numeric_limits<std::uint32_t>::max()) {
            throw std::runtime_error("particle-ID collision ordinal wire value overflow");
          }
          shard_candidates.push_back(ReceivedCandidate{
              .particle_id = received.values[index],
              .birth_key = received.values[index + 1U],
              .origin_index = received.values[index + 2U],
              .ordinal = static_cast<std::uint32_t>(ordinal64),
              .origin_rank = origin_rank,
          });
        }
      }
      std::sort(shard_candidates.begin(), shard_candidates.end(),
                [](const ReceivedCandidate& lhs, const ReceivedCandidate& rhs) {
                  if (lhs.particle_id != rhs.particle_id) return lhs.particle_id < rhs.particle_id;
                  if (lhs.birth_key != rhs.birth_key) return lhs.birth_key < rhs.birth_key;
                  if (lhs.origin_rank != rhs.origin_rank) return lhs.origin_rank < rhs.origin_rank;
                  return lhs.origin_index < rhs.origin_index;
                });

      // decision record: local candidate index, 0=resolved, 1=rehash.
      std::vector<std::vector<std::uint64_t>> decisions(
          static_cast<std::size_t>(world_size));
      std::size_t begin = 0U;
      while (begin < shard_candidates.size()) {
        std::size_t end_group = begin + 1U;
        while (end_group < shard_candidates.size() &&
               shard_candidates[end_group].particle_id ==
                   shard_candidates[begin].particle_id) {
          ++end_group;
        }
        const std::uint64_t particle_id = shard_candidates[begin].particle_id;
        const bool collides_reserved = m_occupied_shard.contains(particle_id) ||
            batch_reserved_on_shard.contains(particle_id);
        const std::size_t first_rehash = collides_reserved ? begin : begin + 1U;
        if (!collides_reserved) {
          const ReceivedCandidate& winner = shard_candidates[begin];
          auto& result = decisions[static_cast<std::size_t>(winner.origin_rank)];
          result.push_back(winner.origin_index);
          result.push_back(0U);
          batch_reserved_on_shard.insert(particle_id);
        }
        for (std::size_t index = first_rehash; index < end_group; ++index) {
          const ReceivedCandidate& collision = shard_candidates[index];
          auto& result = decisions[static_cast<std::size_t>(collision.origin_rank)];
          result.push_back(collision.origin_index);
          result.push_back(1U);
        }
        begin = end_group;
      }

      const ShardedUint64Exchange returned =
          exchangeShardedUint64Records(m_mpi_context, decisions, 2U);
      for (std::size_t index = 0; index < returned.values.size(); index += 2U) {
        const std::uint64_t local_index64 = returned.values[index];
        const std::uint64_t action = returned.values[index + 1U];
        if (local_index64 >= candidates.size() || action > 1U) {
          throw std::runtime_error("particle-ID shard returned an invalid collision decision");
        }
        LocalCandidate& candidate = candidates[static_cast<std::size_t>(local_index64)];
        if (action == 0U) {
          candidate.resolved = true;
          continue;
        }
        if (candidate.collision_ordinal == std::numeric_limits<std::uint32_t>::max()) {
          throw std::runtime_error("ParticleIdRegistry: deterministic collision ordinal overflow");
        }
        ++candidate.collision_ordinal;
        candidate.particle_id = physics::starFormationParticleIdFromBirthKey(
            candidate.birth_key, candidate.collision_ordinal);
      }

      std::uint64_t local_unresolved = 0U;
      for (const LocalCandidate& candidate : candidates) {
        local_unresolved += candidate.resolved ? 0U : 1U;
      }
      if (m_mpi_context.allreduceSumUint64(local_unresolved) == 0U) {
        m_occupied_shard.insert(
            batch_reserved_on_shard.begin(), batch_reserved_on_shard.end());
        std::vector<std::uint64_t> result;
        result.reserve(candidates.size());
        for (const LocalCandidate& candidate : candidates) {
          result.push_back(candidate.particle_id);
        }
        return result;
      }
    }
    throw std::runtime_error(
        "ParticleIdRegistry: deterministic sharded collision resolution exhausted");
  }

 private:
  void initializeOccupiedShard(std::span<const std::uint64_t> local_ids) {
    const int world_size = m_mpi_context.worldSize();
    std::vector<std::vector<std::uint64_t>> ids_by_rank(
        static_cast<std::size_t>(world_size));
    for (const std::uint64_t id : local_ids) {
      ids_by_rank[shardForUint64(id, world_size)].push_back(id);
    }
    const ShardedUint64Exchange received =
        exchangeShardedUint64Records(m_mpi_context, ids_by_rank, 1U);
    std::vector<std::uint64_t> shard_ids = received.values;
    std::sort(shard_ids.begin(), shard_ids.end());
    const bool invalid_local =
        (!shard_ids.empty() && shard_ids.front() == 0U) ||
        std::adjacent_find(shard_ids.begin(), shard_ids.end()) != shard_ids.end();
    if (m_mpi_context.allreduceSumUint64(invalid_local ? 1ULL : 0ULL) != 0U) {
      throw std::runtime_error(
          "ParticleIdRegistry: zero or duplicate existing ID during sharded initialization");
    }
    m_occupied_shard.insert(shard_ids.begin(), shard_ids.end());
    m_initialized = true;
  }

  void validateBirthKeysDistributed(std::span<const std::uint64_t> birth_keys) const {
    const int world_size = m_mpi_context.worldSize();
    std::vector<std::vector<std::uint64_t>> keys_by_rank(
        static_cast<std::size_t>(world_size));
    for (const std::uint64_t key : birth_keys) {
      keys_by_rank[shardForUint64(key, world_size)].push_back(key);
    }
    const ShardedUint64Exchange received =
        exchangeShardedUint64Records(m_mpi_context, keys_by_rank, 1U);
    std::vector<std::uint64_t> shard_keys = received.values;
    std::sort(shard_keys.begin(), shard_keys.end());
    const bool invalid_local =
        (!shard_keys.empty() && shard_keys.front() == 0U) ||
        std::adjacent_find(shard_keys.begin(), shard_keys.end()) != shard_keys.end();
    if (m_mpi_context.allreduceSumUint64(invalid_local ? 1ULL : 0ULL) != 0U) {
      throw std::runtime_error(
          "ParticleIdRegistry: zero or duplicate immutable birth key across owner ranks");
    }
  }

  const parallel::MpiContext& m_mpi_context;
  bool m_initialized = false;
  // Only IDs whose hash shard belongs to this rank are retained. Aggregate
  // memory remains O(N), not O(N * ranks), while collision decisions stay
  // deterministic because each candidate ID has one authoritative shard.
  std::unordered_set<std::uint64_t> m_occupied_shard;
};

[[nodiscard]] physics::BlackHoleAgnConfig makeRuntimeBlackHoleAgnConfig(
    const core::PhysicsConfig& physics_config,
    const core::UnitSystem& units) {
  physics::BlackHoleAgnConfig config = physics::makeBlackHoleAgnConfig(physics_config);
  config.proton_mass_code = config.proton_mass_si / units.mass_si_per_code;
  config.thomson_cross_section_code = config.thomson_cross_section_si /
      (units.length_si_per_code * units.length_si_per_code);
  config.newton_g_code = config.newton_g_si * units.mass_si_per_code *
      units.timeSiPerCode() * units.timeSiPerCode() /
      (units.length_si_per_code * units.length_si_per_code * units.length_si_per_code);
  config.speed_of_light_code = config.speed_of_light_si / units.velocity_si_per_code;
  return config;
}

using PatchCellGeometry = internal::StarFormationPatchCellGeometry;
using internal::starFormationDerivativeAtCell;
using internal::starFormationPatchCellGeometry;
using internal::starFormationPatchIsLeaf;

class FeedbackCellVolumeProvider final
    : public physics::StellarFeedbackCellVolumeProvider {
 public:
  explicit FeedbackCellVolumeProvider(const core::SimulationState& state)
      : m_state(state) {}

  [[nodiscard]] double cellVolumeCode(std::uint32_t cell_index) const override {
    if (cell_index >= m_state.cells.size() || cell_index >= m_state.gas_cells.size()) {
      return 0.0;
    }
    const PatchCellGeometry geometry =
        starFormationPatchCellGeometry(m_state, cell_index);
    if (geometry.valid) {
      return geometry.dx_stored * geometry.dy_stored * geometry.dz_stored;
    }
    const double density = m_state.gas_cells.density_code[cell_index];
    return density > 0.0 ? m_state.cells.mass_code[cell_index] / density : 0.0;
  }

 private:
  const core::SimulationState& m_state;
};

class SourceRuntimeImpl final : public SourceRuntime {
 public:
  SourceRuntimeImpl(
      const core::SimulationConfig& config,
      const core::ModePolicy& mode_policy,
      const core::UnitSystem& units,
      std::uint32_t world_rank,
      const parallel::MpiContext& mpi_context,
      std::shared_ptr<const physics::EffectiveMultiphaseEosTable> effective_eos_table,
      const RuntimeServices* runtime_services)
      : m_units(units),
        m_effective_eos_table(
            effective_eos_table != nullptr
                ? std::move(effective_eos_table)
                : makeRuntimeEffectiveEosTable(config, units)),
        m_star_formation(makeRuntimeStarFormationConfig(config, units), m_effective_eos_table),
        m_stellar_evolution(
            physics::makeStellarEvolutionConfig(config.physics),
            physics::loadStellarEvolutionTable(config.physics)),
        m_stellar_feedback(makeRuntimeStellarFeedbackConfig(config, units)),
        m_metal_diffusion(physics::makeMetalDiffusionConfig(config.physics)),
        m_black_hole(makeRuntimeBlackHoleAgnConfig(config.physics, units)),
        m_particle_id_registry(mpi_context),
        m_mpi_context(mpi_context),
        m_runtime_services(runtime_services),
        m_memory_governor(
            runtime_services != nullptr ? runtime_services->memory_governor : nullptr),
        m_world_rank(world_rank),
        m_coordinate_frame(config.units.coordinate_frame),
        m_hydro_boundary(mode_policy.hydro_boundary),
        m_bh_enabled(config.physics.enable_black_hole_agn),
        m_bh_seeding_requested(config.physics.bh_enable_seeding),
        m_is_cosmological(
            config.mode.mode == core::SimulationMode::kCosmoCube ||
            config.mode.mode == core::SimulationMode::kZoomIn) {
    if (m_bh_seeding_requested) {
      if (!m_bh_enabled) {
        throw std::invalid_argument(
            "physics.bh_enable_seeding=true requires physics.enable_black_hole_agn=true");
      }
      throw std::runtime_error(
          "BH seeding requested, but ReferenceWorkflow has no authoritative halo/candidate provider; "
          "disable bh_enable_seeding or supply a production candidate provider before claiming seeding support");
    }

    const double hubble_si = config.cosmology.hubble_param *
        core::constants::k_hubble_100_km_s_mpc_si;
    const double hubble_code = hubble_si * units.timeSiPerCode();
    const double g_code = newtonGCodeFromUnits(units);
    m_mean_baryon_density0_code = 3.0 * hubble_code * hubble_code /
        (8.0 * core::constants::k_pi * g_code) * config.cosmology.omega_baryon;
  }

  [[nodiscard]] core::MemoryReport memoryReport() const override {
    core::MemoryReportBuilder builder;
    const auto add_container = [&builder](
        core::MemorySubsystem subsystem,
        core::MemoryClass memory_class,
        std::string_view label,
        const auto& container,
        bool governed) {
      builder.addEntry(core::MemoryEntry{
          .subsystem = subsystem,
          .lifetime = core::MemoryLifetime::kTransient,
          .memory_class = memory_class,
          .label = std::string(label),
          .current_size_bytes = core::currentSizeBytesForContainer(container),
          .owned_capacity_bytes = core::ownedCapacityBytesForContainer(container),
          .high_water_bytes = core::ownedCapacityBytesForContainer(container),
          .estimate_only = false,
          .governed_commitment = governed,
      });
    };
    add_container(
        core::MemorySubsystem::kScratch,
        core::MemoryClass::kPhaseResident,
        "sources.star_formation.contiguous_cell_batch",
        m_contiguous_cell_batch,
        false);
    add_container(
        core::MemorySubsystem::kScratch,
        core::MemoryClass::kPhaseResident,
        "sources.star_formation.input_batch",
        m_star_formation_inputs,
        m_star_formation_input_reservation.committed());
    add_container(
        core::MemorySubsystem::kActiveSets,
        core::MemoryClass::kPhaseResident,
        "sources.stellar_evolution.active_star_rows",
        m_active_star_indices,
        false);
    add_container(
        core::MemorySubsystem::kScratch,
        core::MemoryClass::kPhaseResident,
        "sources.stellar_evolution.contiguous_star_batch",
        m_contiguous_star_batch,
        false);
    add_container(
        core::MemorySubsystem::kScratch,
        core::MemoryClass::kPhaseResident,
        "sources.stellar_feedback.event_batch",
        m_feedback_events,
        m_feedback_event_reservation.committed());
    builder.addEntry(core::MemoryEntry{
        .subsystem = core::MemorySubsystem::kScratch,
        .lifetime = core::MemoryLifetime::kTransient,
        .memory_class = core::MemoryClass::kPhaseResident,
        .label = "sources.stellar_feedback.spatial_index",
        .current_size_bytes = m_feedback_spatial_index.ownedCapacityBytes(),
        .owned_capacity_bytes = m_feedback_spatial_index.ownedCapacityBytes(),
        .high_water_bytes = m_feedback_spatial_index.highWaterBytes(),
        .estimate_only = false,
        .governed_commitment = m_feedback_index_reservation.committed(),
    });
    const bool diffusion_governed = m_diffusion_phase_reservation.committed();
    add_container(
        core::MemorySubsystem::kScratch,
        core::MemoryClass::kPhaseResident,
        "sources.metal_diffusion.owned_leaf_mask",
        m_owned_leaf_mask,
        diffusion_governed);
    add_container(
        core::MemorySubsystem::kScratch,
        core::MemoryClass::kPhaseResident,
        "sources.metal_diffusion.rho_kappa",
        m_diffusion_rho_kappa_code,
        diffusion_governed);
    add_container(
        core::MemorySubsystem::kScratch,
        core::MemoryClass::kPhaseResident,
        "sources.metal_diffusion.faces",
        m_diffusion_faces,
        false);
    builder.addEntry(core::MemoryEntry{
        .subsystem = core::MemorySubsystem::kScratch,
        .lifetime = core::MemoryLifetime::kTransient,
        .memory_class = core::MemoryClass::kPhaseResident,
        .label = "sources.metal_diffusion.workspace",
        .current_size_bytes = m_diffusion_workspace.ownedCapacityBytes(),
        .owned_capacity_bytes = m_diffusion_workspace.ownedCapacityBytes(),
        .high_water_bytes = m_diffusion_workspace.highWaterBytes(),
        .estimate_only = false,
        .governed_commitment = diffusion_governed,
    });
    builder.addEntry(core::MemoryEntry{
        .subsystem = core::MemorySubsystem::kScratch,
        .lifetime = core::MemoryLifetime::kTransient,
        .memory_class = core::MemoryClass::kScratchArena,
        .label = "sources.metal_diffusion.topology_construction_scratch",
        .owned_capacity_bytes = 0U,
        .high_water_bytes = m_diffusion_topology_scratch_high_water_bytes,
        .estimate_only = true,
        .governed_commitment = false,
        .uncertainty_note =
            "measured container-capacity high-water for released topology construction scratch",
    });
    return std::move(builder).finish();
  }

  void execute(SourceMutationStageView& view) override {
    view.requireFresh();
    core::StepContext& context = internal::RuntimeStageAccess::sourceContext(
        view,
        {{RuntimeResourceKey::kSourceMutationState, RuntimeResourceAccessMode::kReadWrite},
         {RuntimeResourceKey::kParticlePosition, RuntimeResourceAccessMode::kReadWrite},
         {RuntimeResourceKey::kParticleVelocity, RuntimeResourceAccessMode::kReadWrite},
         {RuntimeResourceKey::kParticleIdentity, RuntimeResourceAccessMode::kReadWrite},
         {RuntimeResourceKey::kParticleSpeciesIndex, RuntimeResourceAccessMode::kReadWrite},
         {RuntimeResourceKey::kHydroConservedState, RuntimeResourceAccessMode::kReadWrite},
         {RuntimeResourceKey::kHydroPrimitiveState, RuntimeResourceAccessMode::kRead},
         {RuntimeResourceKey::kAmrPatchState, RuntimeResourceAccessMode::kRead},
         {RuntimeResourceKey::kEffectiveIsmThermodynamics, RuntimeResourceAccessMode::kRead},
         {RuntimeResourceKey::kMigrationOwnership, RuntimeResourceAccessMode::kReadWrite},
         {RuntimeResourceKey::kIntegratorTruth, RuntimeResourceAccessMode::kRead}});
    if (context.stage != core::IntegrationStage::kSourceTerms) {
      throw std::logic_error("source runtime received a non-source stage");
    }

    // SourceTerms executes after drift/hydro on the step-end state while the
    // committed IntegratorState remains at the step beginning until commitStep().
    // Make the physical evaluation epoch explicit so every source model sees
    // the state and cosmological conversion from the same instant.
    const double source_evaluation_scale_factor = m_is_cosmological
        ? context.timeline_step.scale_factor_end
        : 1.0;
    if (!std::isfinite(source_evaluation_scale_factor) ||
        source_evaluation_scale_factor <= 0.0) {
      throw std::runtime_error(
          "source runtime requires a finite positive source-stage scale factor");
    }

    const std::size_t cell_count = context.state.cells.size();
    const bool local_has_mutable_gravity_sources =
        cell_count > 0U || context.state.star_particles.size() > 0U ||
        context.state.black_holes.size() > 0U;
    const std::uint64_t gravity_source_mutation_rank_count =
        m_mpi_context.allreduceSumUint64(local_has_mutable_gravity_sources ? 1ULL : 0ULL);
    if (m_star_formation.config().enabled) {
      constexpr std::size_t k_star_formation_cell_batch_max = 4096U;
      const bool all_cells_active =
          !context.active_set.cells_are_subset &&
          context.active_set.cell_indices.empty();
      const std::span<const std::uint32_t> active_cells =
          context.active_set.cell_indices;
      const std::size_t local_active_count = all_cells_active
          ? cell_count : active_cells.size();
      const std::uint64_t local_batch_count = local_active_count == 0U
          ? 0U
          : 1U + static_cast<std::uint64_t>(
              (local_active_count - 1U) / k_star_formation_cell_batch_max);
      const std::uint64_t global_batch_count =
          m_mpi_context.allreduceMaxUint64(local_batch_count);
      const std::size_t required_input_capacity = std::min<std::size_t>(
          local_active_count, k_star_formation_cell_batch_max);

      core::MemoryReservation replacement_input_reservation;
      std::exception_ptr reservation_failure;
      try {
        if (m_memory_governor != nullptr &&
            required_input_capacity > m_star_formation_inputs.capacity()) {
          const std::size_t input_bytes = core::checkedSizeMultiply(
              required_input_capacity,
              sizeof(physics::StarFormationCellInput),
              "star-formation bounded input batch");
          replacement_input_reservation = m_memory_governor->reserve(
              core::MemoryClass::kPhaseResident,
              static_cast<std::uint64_t>(input_bytes),
              "sources.star_formation.input_batch");
        }
      } catch (...) {
        reservation_failure = std::current_exception();
      }
      if (m_runtime_services != nullptr) {
        FailureCoordinator(*m_runtime_services).rethrowCollectiveFailure(
            reservation_failure,
            "star-formation memory preflight");
      } else if (reservation_failure != nullptr) {
        std::rethrow_exception(reservation_failure);
      }
      if (replacement_input_reservation.pending()) {
        replacement_input_reservation.commit();
      }
      if (required_input_capacity > m_star_formation_inputs.capacity()) {
        m_star_formation_inputs.reserve(required_input_capacity);
        if (m_memory_governor != nullptr) {
          m_star_formation_input_reservation.release();
          m_star_formation_input_reservation =
              std::move(replacement_input_reservation);
        }
      }
      m_contiguous_cell_batch.reserve(required_input_capacity);

      std::unordered_set<std::uint64_t> patch_ids_with_children;
      patch_ids_with_children.reserve(context.state.patches.size());
      for (std::size_t patch_index = 0U;
           patch_index < context.state.patches.size(); ++patch_index) {
        const std::uint64_t parent_id =
            context.state.patches.parent_patch_id[patch_index];
        if (parent_id != 0U) {
          patch_ids_with_children.insert(parent_id);
        }
      }

      const std::size_t particle_count_before_birth =
          context.state.particles.size();
      std::uint64_t local_spawned_particles = 0U;
      for (std::uint64_t batch_round = 0U;
           batch_round < global_batch_count; ++batch_round) {
        const std::size_t batch_begin = batch_round < local_batch_count
            ? core::checkedSizeMultiply(
                core::checkedIntegralNarrow<std::size_t>(
                    batch_round,
                    "star-formation batch round"),
                k_star_formation_cell_batch_max,
                "star-formation batch begin")
            : local_active_count;
        const std::size_t batch_size = batch_begin < local_active_count
            ? std::min<std::size_t>(
                k_star_formation_cell_batch_max,
                local_active_count - batch_begin)
            : 0U;
        std::span<const std::uint32_t> cell_batch;
        if (batch_size == 0U) {
          m_contiguous_cell_batch.clear();
          cell_batch = {};
        } else if (all_cells_active) {
          m_contiguous_cell_batch.resize(batch_size);
          for (std::size_t i = 0U; i < batch_size; ++i) {
            m_contiguous_cell_batch[i] =
                core::checkedIntegralNarrow<std::uint32_t>(
                    batch_begin + i,
                    "star-formation contiguous active gas row");
          }
          cell_batch = m_contiguous_cell_batch;
        } else {
          cell_batch = active_cells.subspan(batch_begin, batch_size);
        }
        buildStarFormationInputs(
            context,
            cell_batch,
            source_evaluation_scale_factor,
            patch_ids_with_children);
        const physics::StarFormationStepReport batch_report =
            m_star_formation.applyFromInputs(
                context.state,
                m_star_formation_inputs,
                context.integrator_state.dt_time_code,
                source_evaluation_scale_factor,
                context.integrator_state.step_index,
                &m_particle_id_registry);
        if (batch_report.counters.spawned_particles >
            std::numeric_limits<std::uint64_t>::max() - local_spawned_particles) {
          throw std::overflow_error(
              "source runtime star-formation batch spawn count overflows uint64");
        }
        local_spawned_particles += batch_report.counters.spawned_particles;
      }

      if (local_spawned_particles > 0U) {
        const std::size_t particle_count_after_birth =
            context.state.particles.size();
        if (particle_count_after_birth < particle_count_before_birth ||
            particle_count_after_birth - particle_count_before_birth !=
                static_cast<std::size_t>(local_spawned_particles)) {
          throw std::runtime_error(
              "source runtime star-formation report disagrees with appended particle rows");
        }
        internal::synchronizeParentParticleCompatibilityMirrors(
            context.state, m_world_rank, "SourceRuntime star-formation batch");
        if (context.newly_created_particle_ids != nullptr) {
          context.newly_created_particle_ids->insert(
              context.newly_created_particle_ids->end(),
              context.state.particle_sidecar.particle_id.begin() +
                  static_cast<std::ptrdiff_t>(particle_count_before_birth),
              context.state.particle_sidecar.particle_id.end());
        }
      }
      if (local_spawned_particles > 0U &&
          context.particle_scheduler != nullptr) {
        const std::uint64_t current_tick =
            context.particle_scheduler->currentTick();
        if (current_tick == std::numeric_limits<std::uint64_t>::max()) {
          throw std::overflow_error(
              "source runtime newborn activation tick overflows uint64");
        }
        context.particle_scheduler->appendElements(
            core::checkedIntegralNarrow<std::uint32_t>(
                local_spawned_particles,
                "source runtime newborn scheduler append count"),
            0U,
            current_tick + 1U);
      }
    }

    executeStellarEvolutionAndEnrichment(context);
    executeMetalDiffusion(context, source_evaluation_scale_factor);

    physics::BlackHoleAgnStepReport black_hole_report{};
    if (m_bh_enabled) {
      black_hole_report = m_black_hole.apply(
          context.state,
          {},
          context.integrator_state.dt_time_code,
          source_evaluation_scale_factor,
          m_coordinate_frame == core::CoordinateFrame::kComoving,
          context.integrator_state.step_index,
          &m_particle_id_registry);
    }
    if (black_hole_report.counters.gas_mass_removed_code > 0.0) {
      // Generic gas particles are compatibility mirrors only. Keep them coherent
      // for I/O/legacy consumers after the authoritative gas->BH transaction.
      internal::synchronizeParentParticleCompatibilityMirrors(
          context.state, m_world_rank, "SourceRuntime BH accretion/seeding batch");
    }
    if (!black_hole_report.seeded_particle_ids.empty()) {
      if (context.newly_created_particle_ids != nullptr) {
        context.newly_created_particle_ids->insert(
            context.newly_created_particle_ids->end(),
            black_hole_report.seeded_particle_ids.begin(),
            black_hole_report.seeded_particle_ids.end());
      }
      if (context.particle_scheduler != nullptr) {
        const std::uint64_t current_tick = context.particle_scheduler->currentTick();
        if (current_tick == std::numeric_limits<std::uint64_t>::max()) {
          throw std::overflow_error("source runtime BH activation tick overflows uint64");
        }
        context.particle_scheduler->appendElements(
            static_cast<std::uint32_t>(black_hole_report.seeded_particle_ids.size()),
            0U,
            current_tick + 1U);
      }
    }
    if (gravity_source_mutation_rank_count > 0U) {
      // Source models can change gas/star/BH gravitating mass and membership.
      // Advance once per globally ordered source stage on every rank so PM/tree
      // validity follows authoritative state rather than PM refresh sequencing.
      context.state.bumpGravitySourceGeneration();
    }
  }

 private:
  void buildOwnedLeafCellMetadata(const core::StepContext& context) {
    const core::SimulationState& state = context.state;
    const std::size_t cell_count = state.cells.size();
    m_owned_leaf_mask.assign(cell_count, 0U);
    std::unordered_set<std::uint64_t> patch_ids_with_children;
    patch_ids_with_children.reserve(state.patches.size());
    for (std::size_t patch_index = 0U; patch_index < state.patches.size(); ++patch_index) {
      const std::uint64_t parent_id = state.patches.parent_patch_id[patch_index];
      if (parent_id != 0U) {
        patch_ids_with_children.insert(parent_id);
      }
    }
    for (std::uint32_t cell_index = 0U; cell_index < cell_count; ++cell_index) {
      const PatchCellGeometry geometry = starFormationPatchCellGeometry(state, cell_index);
      bool owned_leaf = true;
      if (geometry.valid) {
        owned_leaf = state.patches.owning_rank[geometry.patch_index] == m_world_rank &&
            starFormationPatchIsLeaf(
                state, geometry.patch_index, patch_ids_with_children);
      }
      m_owned_leaf_mask[cell_index] = owned_leaf ? 1U : 0U;
    }
  }

  void buildActiveStarRows(const core::StepContext& context) {
    const core::SimulationState& state = context.state;
    m_active_star_indices.clear();
    if (!context.active_set.particles_are_subset ||
        context.active_set.particle_indices.empty()) {
      m_active_star_indices.resize(state.star_particles.size());
      std::iota(m_active_star_indices.begin(), m_active_star_indices.end(), 0U);
      return;
    }
    std::unordered_set<std::uint32_t> active_particles(
        context.active_set.particle_indices.begin(),
        context.active_set.particle_indices.end());
    for (std::uint32_t star_index = 0;
         star_index < state.star_particles.size(); ++star_index) {
      if (active_particles.contains(
              state.star_particles.particle_index[star_index])) {
        m_active_star_indices.push_back(star_index);
      }
    }
  }

  [[nodiscard]] double elapsedStellarEvolutionYears(
      const core::StepContext& context) const {
    constexpr double k_seconds_per_year = 31557600.0;
    double elapsed_si = context.timeline_step.dt_time_si;
    if (m_is_cosmological && context.cosmology_background != nullptr &&
        context.timeline_step.scale_factor_begin > 0.0 &&
        context.timeline_step.scale_factor_end >=
            context.timeline_step.scale_factor_begin) {
      elapsed_si = context.cosmology_background->cosmicTimeIntervalSi(
          context.timeline_step.scale_factor_begin,
          context.timeline_step.scale_factor_end);
    } else if (!(elapsed_si > 0.0)) {
      elapsed_si = context.integrator_state.dt_time_code * m_units.timeSiPerCode();
    }
    return std::max(elapsed_si / k_seconds_per_year, 0.0);
  }

  void executeStellarEvolutionAndEnrichment(core::StepContext& context) {
    if (!m_stellar_evolution.config().enabled) {
      return;
    }
    const bool all_stars_active =
        !context.active_set.particles_are_subset ||
        context.active_set.particle_indices.empty();
    if (all_stars_active) {
      m_active_star_indices.clear();
    } else {
      buildActiveStarRows(context);
    }
    const std::size_t active_star_count = all_stars_active
        ? context.state.star_particles.size()
        : m_active_star_indices.size();
    constexpr std::size_t k_feedback_event_batch_max = 4096U;
    const std::size_t required_event_capacity = std::min<std::size_t>(
        active_star_count, k_feedback_event_batch_max);
    const std::size_t feedback_index_bytes_size = active_star_count == 0U
        ? 0U
        : core::checkedSizeMultiply(
              context.state.cells.size(),
              sizeof(std::uint32_t),
              "stellar-feedback spatial-index staging");
    const std::size_t event_bytes_size =
        required_event_capacity > m_feedback_events.capacity()
        ? core::checkedSizeMultiply(
              required_event_capacity,
              sizeof(physics::StellarFeedbackEvent),
              "stellar-feedback event batch")
        : 0U;
    core::MemoryReservation feedback_index_rebuild_reservation;
    core::MemoryReservation replacement_event_reservation;
    std::exception_ptr reservation_failure;
    try {
      if (m_memory_governor != nullptr && feedback_index_bytes_size != 0U) {
        feedback_index_rebuild_reservation = m_memory_governor->reserve(
            core::MemoryClass::kPhaseResident,
            static_cast<std::uint64_t>(feedback_index_bytes_size),
            "sources.stellar_feedback.spatial_index_rebuild");
      }
      if (m_memory_governor != nullptr && event_bytes_size != 0U) {
        replacement_event_reservation = m_memory_governor->reserve(
            core::MemoryClass::kPhaseResident,
            static_cast<std::uint64_t>(event_bytes_size),
            "sources.stellar_feedback.event_batch");
      }
    } catch (...) {
      reservation_failure = std::current_exception();
    }
    if (m_runtime_services != nullptr) {
      FailureCoordinator(*m_runtime_services).rethrowCollectiveFailure(
          reservation_failure,
          "stellar-feedback memory preflight");
    } else if (reservation_failure != nullptr) {
      std::rethrow_exception(reservation_failure);
    }
    if (active_star_count == 0U) {
      return;
    }
    if (feedback_index_rebuild_reservation.pending()) {
      feedback_index_rebuild_reservation.commit();
    }
    if (replacement_event_reservation.pending()) {
      replacement_event_reservation.commit();
    }
    const double elapsed_years = elapsedStellarEvolutionYears(context);

    std::unordered_set<std::uint64_t> patch_ids_with_children;
    patch_ids_with_children.reserve(context.state.patches.size());
    for (std::size_t patch_index = 0;
         patch_index < context.state.patches.size(); ++patch_index) {
      const std::uint64_t parent_id =
          context.state.patches.parent_patch_id[patch_index];
      if (parent_id != 0U) {
        patch_ids_with_children.insert(parent_id);
      }
    }
    std::vector<std::uint32_t> owned_leaf_cells;
    owned_leaf_cells.reserve(context.state.cells.size());
    for (std::uint32_t cell_index = 0;
         cell_index < context.state.cells.size(); ++cell_index) {
      const PatchCellGeometry geometry =
          starFormationPatchCellGeometry(context.state, cell_index);
      bool owned_leaf = true;
      if (geometry.valid) {
        owned_leaf =
            context.state.patches.owning_rank[geometry.patch_index] == m_world_rank &&
            starFormationPatchIsLeaf(
                context.state, geometry.patch_index, patch_ids_with_children);
      }
      if (owned_leaf) {
        owned_leaf_cells.push_back(cell_index);
      }
    }

    const physics::StellarFeedbackGeometryView geometry_view{
        .particle_position_x_comoving =
            context.state.particles.position_x_comoving,
        .particle_position_y_comoving =
            context.state.particles.position_y_comoving,
        .particle_position_z_comoving =
            context.state.particles.position_z_comoving,
        .cell_center_x_comoving = context.state.cells.center_x_comoving,
        .cell_center_y_comoving = context.state.cells.center_y_comoving,
        .cell_center_z_comoving = context.state.cells.center_z_comoving,
        .gas_cell_id = context.state.gas_cells.gas_cell_id,
    };
    m_feedback_spatial_index.rebuildFromCellIndices(
        geometry_view, std::move(owned_leaf_cells));
    if (m_memory_governor != nullptr) {
      // The moved candidate allocation has become the retained index and the
      // old index buffer has been released by vector move-assignment. Transfer
      // the committed rebuild reservation only after that ownership change.
      m_feedback_index_reservation.release();
      m_feedback_index_reservation = std::move(feedback_index_rebuild_reservation);
    }
    FeedbackCellVolumeProvider volume_provider(context.state);
    physics::StellarFeedbackDepositionView deposition_view{
        .cell_mass_code = context.state.cells.mass_code,
        .gas_density_code = context.state.gas_cells.density_code,
        .gas_internal_energy_code =
            context.state.gas_cells.internal_energy_code,
        .gas_metal_mass_code = context.state.gas_cells.metal_mass_code,
        .cell_volume_provider = &volume_provider,
    };
    m_feedback_events.clear();
    if (required_event_capacity > m_feedback_events.capacity()) {
      m_feedback_events.reserve(required_event_capacity);
      if (m_memory_governor != nullptr) {
        m_feedback_event_reservation.release();
        m_feedback_event_reservation = std::move(replacement_event_reservation);
      }
    }
    m_contiguous_star_batch.reserve(required_event_capacity);
    for (std::size_t batch_begin = 0U;
         batch_begin < active_star_count;
         batch_begin += k_feedback_event_batch_max) {
      const std::size_t batch_size = std::min<std::size_t>(
          k_feedback_event_batch_max, active_star_count - batch_begin);
      std::span<const std::uint32_t> star_batch;
      if (all_stars_active) {
        m_contiguous_star_batch.resize(batch_size);
        for (std::size_t i = 0U; i < batch_size; ++i) {
          m_contiguous_star_batch[i] = core::checkedIntegralNarrow<std::uint32_t>(
              batch_begin + i,
              "stellar-feedback contiguous active star row");
        }
        star_batch = m_contiguous_star_batch;
      } else {
        star_batch = std::span<const std::uint32_t>(m_active_star_indices)
            .subspan(batch_begin, batch_size);
      }
      const physics::StellarEvolutionStepReport evolution_report =
          m_stellar_evolution.evaluateElapsedYears(
              context.state, star_batch, elapsed_years);
      m_feedback_events.clear();
      for (const physics::StellarEvolutionStarBudget& budget :
           evolution_report.budgets) {
        m_feedback_events.push_back(physics::StellarFeedbackEvent{
            .star_index = budget.star_index,
            .returned_mass_code = budget.interval.returned_mass_code,
            .returned_metals_code = budget.interval.returned_metals_code,
            .feedback_energy_erg = budget.interval.feedback_energy_erg,
        });
      }
      if (!m_feedback_events.empty()) {
        (void)m_stellar_feedback.applyEventsWithViews(
            context.state,
            m_stellar_feedback_state,
            geometry_view,
            &m_feedback_spatial_index,
            deposition_view,
            m_feedback_events,
            context.integrator_state.dt_time_code);
      }

      // Each bounded event batch is deposited or durably carried before its
      // matching stellar-evolution ledger advances. No processed event history
      // survives beyond this batch.
      m_stellar_evolution.commitBudgets(context.state, evolution_report);
    }
    internal::synchronizeParentParticleCompatibilityMirrors(
        context.state, m_world_rank,
        "SourceRuntime stellar-evolution enrichment batch");
  }

  void executeMetalDiffusion(
      core::StepContext& context,
      double source_evaluation_scale_factor) {
    if (!m_metal_diffusion.config().enabled || context.state.cells.size() == 0U) {
      return;
    }
    if (m_mpi_context.isEnabled() && m_mpi_context.worldSize() > 1) {
      throw std::runtime_error(
          "metal diffusion is not permitted on multiple MPI ranks until the directed gas-cell interface exchange "
          "and equal-and-opposite remote flux commit are available; refusing a decomposition-dependent result");
    }

    const std::size_t cell_count = context.state.cells.size();
    const std::uint64_t workspace_bytes =
        m_metal_diffusion.requiredWorkspaceBytes(cell_count);
    const std::size_t mask_bytes = core::checkedSizeMultiply(
        cell_count, sizeof(std::uint8_t), "metal-diffusion owned-leaf mask");
    const std::size_t rho_kappa_bytes = core::checkedSizeMultiply(
        cell_count, sizeof(double), "metal-diffusion rho-kappa field");
    std::uint64_t phase_bytes = core::checkedMemoryBytesAdd(
        workspace_bytes,
        static_cast<std::uint64_t>(mask_bytes),
        "metal-diffusion phase workspace plus ownership mask");
    phase_bytes = core::checkedMemoryBytesAdd(
        phase_bytes,
        static_cast<std::uint64_t>(rho_kappa_bytes),
        "metal-diffusion phase plus rho-kappa field");

    core::MemoryReservation replacement_phase_reservation;
    std::exception_ptr reservation_failure;
    try {
      if (m_memory_governor != nullptr &&
          (m_owned_leaf_mask.capacity() < cell_count ||
           m_diffusion_rho_kappa_code.capacity() < cell_count ||
           m_diffusion_workspace.ownedCapacityBytes() < workspace_bytes)) {
        replacement_phase_reservation = m_memory_governor->reserve(
            core::MemoryClass::kPhaseResident,
            phase_bytes,
            "sources.metal_diffusion.cell_phase");
      }
    } catch (...) {
      reservation_failure = std::current_exception();
    }
    if (m_runtime_services != nullptr) {
      FailureCoordinator(*m_runtime_services).rethrowCollectiveFailure(
          reservation_failure,
          "metal-diffusion memory preflight");
    } else if (reservation_failure != nullptr) {
      std::rethrow_exception(reservation_failure);
    }
    if (replacement_phase_reservation.pending()) {
      replacement_phase_reservation.commit();
    }

    buildOwnedLeafCellMetadata(context);
    const double scale_factor = m_coordinate_frame == core::CoordinateFrame::kComoving
        ? source_evaluation_scale_factor
        : 1.0;
    internal::MetalDiffusionTopologyResult topology =
        internal::buildMetalDiffusionTopology(
            context.state,
            m_owned_leaf_mask,
            m_world_rank,
            scale_factor,
            m_hydro_boundary);
    m_diffusion_topology_scratch_high_water_bytes = std::max(
        m_diffusion_topology_scratch_high_water_bytes,
        topology.construction_scratch_high_water_bytes);
    m_diffusion_faces = std::move(topology.faces);
    m_diffusion_rho_kappa_code = std::move(topology.strain_magnitude_code);
    if (m_diffusion_rho_kappa_code.size() != cell_count) {
      throw std::runtime_error(
          "metal diffusion topology returned the wrong strain-field extent");
    }

    for (std::uint32_t cell_index = 0U; cell_index < cell_count; ++cell_index) {
      if (m_owned_leaf_mask[cell_index] == 0U) {
        m_diffusion_rho_kappa_code[cell_index] = 0.0;
        continue;
      }
      const PatchCellGeometry geometry =
          starFormationPatchCellGeometry(context.state, cell_index);
      double volume_stored = context.state.gas_cells.density_code[cell_index] > 0.0
          ? context.state.cells.mass_code[cell_index] /
                context.state.gas_cells.density_code[cell_index]
          : 0.0;
      if (geometry.valid) {
        volume_stored = geometry.dx_stored * geometry.dy_stored *
            geometry.dz_stored;
      }
      const double volume_phys = volume_stored * scale_factor * scale_factor *
          scale_factor;
      const double gas_mass = context.state.cells.mass_code[cell_index];
      const double density = volume_phys > 0.0 ? gas_mass / volume_phys : 0.0;
      const double filter_length = volume_phys > 0.0
          ? std::cbrt(volume_phys) : 0.0;
      m_diffusion_rho_kappa_code[cell_index] =
          physics::smagorinskyRhoKappaCode(
              m_metal_diffusion.config(),
              density,
              filter_length,
              m_diffusion_rho_kappa_code[cell_index]);
    }

    if (m_diffusion_faces.empty()) {
      if (m_memory_governor != nullptr &&
          replacement_phase_reservation.committed()) {
        m_diffusion_phase_reservation.release();
        m_diffusion_phase_reservation = std::move(replacement_phase_reservation);
      }
      return;
    }
    const physics::MetalDiffusionStepReport report =
        m_metal_diffusion.advanceFromView(
            physics::MetalDiffusionFieldView{
                .gas_mass_code = context.state.cells.mass_code,
                .metal_mass_code = context.state.gas_cells.metal_mass_code,
                .rho_kappa_code = m_diffusion_rho_kappa_code,
            },
            m_diffusion_faces,
            context.integrator_state.dt_time_code,
            m_diffusion_workspace);
    if (m_memory_governor != nullptr &&
        replacement_phase_reservation.committed()) {
      m_diffusion_phase_reservation.release();
      m_diffusion_phase_reservation = std::move(replacement_phase_reservation);
    }

    std::ostringstream metadata;
    metadata << "module=metal_diffusion\n";
    metadata << "model=" << core::metalDiffusionModelToString(
        m_metal_diffusion.config().model) << "\n";
    metadata << "time_integrator=" <<
        core::metalDiffusionTimeIntegratorToString(
            m_metal_diffusion.config().time_integrator) << "\n";
    metadata << "stable_dt_code=" << report.stable_dt_code << "\n";
    metadata << "subcycles=" << report.subcycles << "\n";
    metadata << "rkl_stages=" << report.rkl_stages << "\n";
    metadata << "conservation_residual_code=" <<
        report.conservation_residual_code << "\n";
    const std::string text = metadata.str();
    core::ModuleSidecarBlock sidecar;
    sidecar.module_name = "metal_diffusion";
    sidecar.schema_version = 1U;
    sidecar.payload.resize(text.size());
    for (std::size_t i = 0U; i < text.size(); ++i) {
      sidecar.payload[i] = static_cast<std::byte>(text[i]);
    }
    context.state.sidecars.upsert(std::move(sidecar));
  }

  void buildStarFormationInputs(
      const core::StepContext& context,
      std::span<const std::uint32_t> active_cells,
      double source_evaluation_scale_factor,
      const std::unordered_set<std::uint64_t>& patch_ids_with_children) {
    const core::SimulationState& state = context.state;
    m_star_formation_inputs.clear();

    const double scale_factor = source_evaluation_scale_factor;
    const double length_to_physical =
        m_coordinate_frame == core::CoordinateFrame::kComoving ? scale_factor : 1.0;

    for (const std::uint32_t cell_index : active_cells) {
      physics::StarFormationCellInput input;
      input.cell_index = cell_index;
      input.owning_rank = m_world_rank;
      input.is_active = true;
      if (cell_index >= state.cells.size() || cell_index >= state.gas_cells.size()) {
        input.is_active = false;
        m_star_formation_inputs.push_back(input);
        continue;
      }

      input.gas_cell_id = state.gas_cells.gas_cell_id[cell_index];
      if (input.gas_cell_id == 0U) {
        input.gas_cell_id = state.gasCellIdForLocalRow(cell_index).value_or(0U);
      }
      input.gas_mass_code = state.cells.mass_code[cell_index];
      input.gas_density_code = state.gas_cells.density_code[cell_index];
      input.gas_temperature_k = state.gas_cells.temperature_code[cell_index];
      input.gas_sound_speed_code = state.gas_cells.sound_speed_code[cell_index];
      input.gas_specific_internal_energy_code =
          state.gas_cells.internal_energy_code[cell_index];
      input.is_cosmological = m_is_cosmological;
      if (m_is_cosmological && m_mean_baryon_density0_code > 0.0 && scale_factor > 0.0) {
        const double density_phys = m_coordinate_frame == core::CoordinateFrame::kComoving
            ? input.gas_density_code / (scale_factor * scale_factor * scale_factor)
            : input.gas_density_code;
        const double mean_density_phys = m_mean_baryon_density0_code /
            (scale_factor * scale_factor * scale_factor);
        input.baryon_overdensity = density_phys / std::max(mean_density_phys, 1.0e-30);
      }
      input.velocity_x_peculiar = state.gas_cells.velocity_x_peculiar[cell_index];
      input.velocity_y_peculiar = state.gas_cells.velocity_y_peculiar[cell_index];
      input.velocity_z_peculiar = state.gas_cells.velocity_z_peculiar[cell_index];
      input.gas_metal_mass_code = state.gas_cells.metal_mass_code[cell_index];
      input.metallicity_mass_fraction = input.gas_mass_code > 0.0
          ? input.gas_metal_mass_code / input.gas_mass_code
          : 0.0;
      input.center_x_comoving = state.cells.center_x_comoving[cell_index];
      input.center_y_comoving = state.cells.center_y_comoving[cell_index];
      input.center_z_comoving = state.cells.center_z_comoving[cell_index];

      const PatchCellGeometry geometry = starFormationPatchCellGeometry(state, cell_index);
      if (geometry.valid) {
        input.cell_volume_code = geometry.dx_stored * geometry.dy_stored * geometry.dz_stored;
        input.is_owned = state.patches.owning_rank[geometry.patch_index] == m_world_rank;
        input.is_leaf = starFormationPatchIsLeaf(state, geometry.patch_index, patch_ids_with_children);
        input.is_ghost = !input.is_owned;

        const double dx_phys = geometry.dx_stored * length_to_physical;
        const double dy_phys = geometry.dy_stored * length_to_physical;
        const double dz_phys = geometry.dz_stored * length_to_physical;
        const std::array<std::span<const double>, 3> velocity_fields{
            state.gas_cells.velocity_x_peculiar,
            state.gas_cells.velocity_y_peculiar,
            state.gas_cells.velocity_z_peculiar};
        std::array<std::array<double, 3>, 3> gradient{};
        bool gradient_valid = true;
        for (std::size_t component = 0; component < velocity_fields.size(); ++component) {
          gradient[component][0] = starFormationDerivativeAtCell(
              velocity_fields[component], geometry, 0, dx_phys, cell_index);
          gradient[component][1] = starFormationDerivativeAtCell(
              velocity_fields[component], geometry, 1, dy_phys, cell_index);
          gradient[component][2] = starFormationDerivativeAtCell(
              velocity_fields[component], geometry, 2, dz_phys, cell_index);
          for (std::size_t axis = 0; axis < gradient[component].size(); ++axis) {
            gradient_valid = gradient_valid && std::isfinite(gradient[component][axis]);
          }
        }
        if (gradient_valid) {
          input.velocity_divergence_code =
              gradient[0][0] + gradient[1][1] + gradient[2][2];
          input.velocity_gradient_frobenius_sq_code = 0.0;
          for (const auto& row : gradient) {
            for (const double value : row) {
              input.velocity_gradient_frobenius_sq_code += value * value;
            }
          }
        } else {
          input.velocity_divergence_code = std::numeric_limits<double>::quiet_NaN();
          input.velocity_gradient_frobenius_sq_code =
              std::numeric_limits<double>::quiet_NaN();
        }
      } else {
        input.cell_volume_code = input.gas_density_code > 0.0
            ? input.gas_mass_code / input.gas_density_code
            : 0.0;
        input.is_owned = true;
        input.is_leaf = true;
        input.is_ghost = false;
        // Adaptive production requires a real patch gradient. Legacy mode keeps the
        // historical converging-flow value of zero for non-AMR compatibility states.
        if (m_star_formation.config().model ==
            core::StarFormationModelKind::kAdaptiveBoundJeans) {
          input.velocity_divergence_code = std::numeric_limits<double>::quiet_NaN();
          input.velocity_gradient_frobenius_sq_code =
              std::numeric_limits<double>::quiet_NaN();
        } else {
          input.velocity_divergence_code = 0.0;
          input.velocity_gradient_frobenius_sq_code = 0.0;
        }
      }
      m_star_formation_inputs.push_back(input);
    }
  }

  core::UnitSystem m_units;
  std::shared_ptr<const physics::EffectiveMultiphaseEosTable> m_effective_eos_table;
  physics::StarFormationModel m_star_formation;
  physics::StellarEvolutionBookkeeper m_stellar_evolution;
  physics::StellarFeedbackModel m_stellar_feedback;
  physics::MetalDiffusionModel m_metal_diffusion;
  physics::StellarFeedbackModuleState m_stellar_feedback_state;
  physics::BlackHoleAgnModel m_black_hole;
  DistributedParticleIdRegistry m_particle_id_registry;
  const parallel::MpiContext& m_mpi_context;
  const RuntimeServices* m_runtime_services = nullptr;
  core::MemoryGovernor* m_memory_governor = nullptr;
  core::MemoryReservation m_feedback_index_reservation;
  core::MemoryReservation m_feedback_event_reservation;
  core::MemoryReservation m_star_formation_input_reservation;
  core::MemoryReservation m_diffusion_phase_reservation;
  std::uint32_t m_world_rank = 0;
  core::CoordinateFrame m_coordinate_frame = core::CoordinateFrame::kComoving;
  core::BoundaryCondition m_hydro_boundary = core::BoundaryCondition::kOpen;
  bool m_bh_enabled = false;
  bool m_bh_seeding_requested = false;
  bool m_is_cosmological = false;
  double m_mean_baryon_density0_code = 0.0;
  std::vector<std::uint32_t> m_contiguous_cell_batch;
  std::vector<physics::StarFormationCellInput> m_star_formation_inputs;
  std::vector<std::uint32_t> m_active_star_indices;
  std::vector<std::uint32_t> m_contiguous_star_batch;
  std::vector<physics::StellarFeedbackEvent> m_feedback_events;
  physics::StellarFeedbackSpatialIndex m_feedback_spatial_index;
  std::vector<std::uint8_t> m_owned_leaf_mask;
  std::vector<double> m_diffusion_rho_kappa_code;
  std::vector<physics::MetalDiffusionFace> m_diffusion_faces;
  physics::MetalDiffusionWorkspace m_diffusion_workspace;
  std::uint64_t m_diffusion_topology_scratch_high_water_bytes = 0U;
};

}  // namespace

std::unique_ptr<SourceRuntime> makeSourceRuntime(
    const core::SimulationConfig& config,
    const core::ModePolicy& mode_policy,
    const core::UnitSystem& units,
    std::uint32_t world_rank,
    const parallel::MpiContext& mpi_context,
    std::shared_ptr<const physics::EffectiveMultiphaseEosTable> effective_eos_table,
    const RuntimeServices* runtime_services) {
  return std::make_unique<SourceRuntimeImpl>(
      config,
      mode_policy,
      units,
      world_rank,
      mpi_context,
      std::move(effective_eos_table),
      runtime_services);
}

}  // namespace cosmosim::workflows
