#include "cosmosim/workflows/source_runtime.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <numeric>
#include <span>
#include <stdexcept>
#include <unordered_set>
#include <vector>

#include "cosmosim/physics/black_hole_agn.hpp"
#include "cosmosim/physics/effective_multiphase_ism.hpp"
#include "cosmosim/parallel/distributed_memory.hpp"
#include "cosmosim/physics/star_formation.hpp"
#include "cosmosim/physics/stellar_evolution.hpp"
#include "cosmosim/physics/stellar_feedback.hpp"
#include "cosmosim/physics/metal_diffusion.hpp"
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
  std::vector<int> recv_counts;
  std::vector<int> recv_displacements;
};

[[nodiscard]] ShardedUint64Exchange exchangeShardedUint64Records(
    const parallel::MpiContext& mpi_context,
    const std::vector<std::vector<std::uint64_t>>& values_by_rank,
    std::size_t record_width) {
  if (record_width == 0U) {
    throw std::invalid_argument("sharded uint64 exchange requires nonzero record width");
  }
  const int world_size = mpi_context.worldSize();
  if (world_size <= 0 || values_by_rank.size() != static_cast<std::size_t>(world_size)) {
    throw std::invalid_argument("sharded uint64 exchange rank extent mismatch");
  }
  for (const auto& values : values_by_rank) {
    if (values.size() % record_width != 0U) {
      throw std::invalid_argument("sharded uint64 exchange received a partial record");
    }
  }

  if (!mpi_context.isEnabled()) {
    return ShardedUint64Exchange{
        .values = values_by_rank.front(),
        .recv_counts = {static_cast<int>(values_by_rank.front().size())},
        .recv_displacements = {0},
    };
  }
#if COSMOSIM_ENABLE_MPI
  std::vector<int> send_counts(static_cast<std::size_t>(world_size), 0);
  std::vector<int> send_displacements(static_cast<std::size_t>(world_size), 0);
  std::size_t total_send = 0U;
  for (int rank = 0; rank < world_size; ++rank) {
    const std::size_t count = values_by_rank[static_cast<std::size_t>(rank)].size();
    if (count > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
      throw std::overflow_error("sharded uint64 send count exceeds MPI int range");
    }
    send_counts[static_cast<std::size_t>(rank)] = static_cast<int>(count);
    if (rank > 0) {
      send_displacements[static_cast<std::size_t>(rank)] =
          send_displacements[static_cast<std::size_t>(rank - 1)] +
          send_counts[static_cast<std::size_t>(rank - 1)];
    }
    total_send += count;
  }
  if (total_send > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    throw std::overflow_error("sharded uint64 aggregate send count exceeds MPI int range");
  }
  std::vector<std::uint64_t> send_values;
  send_values.reserve(total_send);
  for (const auto& values : values_by_rank) {
    send_values.insert(send_values.end(), values.begin(), values.end());
  }

  std::vector<int> recv_counts(static_cast<std::size_t>(world_size), 0);
  const int counts_status = MPI_Alltoall(
      send_counts.data(), 1, MPI_INT,
      recv_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);
  if (counts_status != MPI_SUCCESS) {
    throw std::runtime_error("sharded uint64 MPI_Alltoall count exchange failed");
  }
  std::vector<int> recv_displacements(static_cast<std::size_t>(world_size), 0);
  std::size_t total_recv = 0U;
  for (int rank = 0; rank < world_size; ++rank) {
    const int count = recv_counts[static_cast<std::size_t>(rank)];
    if (count < 0 || static_cast<std::size_t>(count) % record_width != 0U) {
      throw std::runtime_error("sharded uint64 receive count violates record framing");
    }
    if (rank > 0) {
      recv_displacements[static_cast<std::size_t>(rank)] =
          recv_displacements[static_cast<std::size_t>(rank - 1)] +
          recv_counts[static_cast<std::size_t>(rank - 1)];
    }
    total_recv += static_cast<std::size_t>(count);
  }
  if (total_recv > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    throw std::overflow_error("sharded uint64 aggregate receive count exceeds MPI int range");
  }
  std::vector<std::uint64_t> recv_values(total_recv, 0U);
  const int payload_status = MPI_Alltoallv(
      send_values.data(), send_counts.data(), send_displacements.data(), MPI_UINT64_T,
      recv_values.data(), recv_counts.data(), recv_displacements.data(), MPI_UINT64_T,
      MPI_COMM_WORLD);
  if (payload_status != MPI_SUCCESS) {
    throw std::runtime_error("sharded uint64 MPI_Alltoallv payload exchange failed");
  }
  return ShardedUint64Exchange{
      .values = std::move(recv_values),
      .recv_counts = std::move(recv_counts),
      .recv_displacements = std::move(recv_displacements),
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
        const int begin = received.recv_displacements[static_cast<std::size_t>(origin_rank)];
        const int count = received.recv_counts[static_cast<std::size_t>(origin_rank)];
        for (int offset = 0; offset < count; offset += 4) {
          const std::size_t index = static_cast<std::size_t>(begin + offset);
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

class SourceRuntimeImpl final : public SourceRuntime {
 public:
  SourceRuntimeImpl(
      const core::SimulationConfig& config,
      const core::ModePolicy& mode_policy,
      const core::UnitSystem& units,
      std::uint32_t world_rank,
      const parallel::MpiContext& mpi_context)
      : m_units(units),
        m_effective_eos_table(makeRuntimeEffectiveEosTable(config, units)),
        m_star_formation(makeRuntimeStarFormationConfig(config, units), m_effective_eos_table),
        m_stellar_evolution(
            physics::makeStellarEvolutionConfig(config.physics),
            physics::loadStellarEvolutionTable(config.physics)),
        m_stellar_feedback(makeRuntimeStellarFeedbackConfig(config, units)),
        m_metal_diffusion(physics::makeMetalDiffusionConfig(config.physics)),
        m_black_hole(makeRuntimeBlackHoleAgnConfig(config.physics, units)),
        m_particle_id_registry(mpi_context),
        m_mpi_context(mpi_context),
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
    if (cell_count > 0U && m_star_formation.config().enabled) {
      std::span<const std::uint32_t> active_cells = context.active_set.cell_indices;
      if (!context.active_set.cells_are_subset && active_cells.empty()) {
        m_full_cell_indices.resize(cell_count);
        std::iota(m_full_cell_indices.begin(), m_full_cell_indices.end(), 0U);
        active_cells = m_full_cell_indices;
      }
      buildStarFormationInputs(context, active_cells, source_evaluation_scale_factor);
      const std::size_t particle_count_before_birth = context.state.particles.size();
      const physics::StarFormationStepReport report = m_star_formation.applyFromInputs(
          context.state,
          m_star_formation_inputs,
          context.integrator_state.dt_time_code,
          source_evaluation_scale_factor,
          context.integrator_state.step_index,
          &m_particle_id_registry);
      if (report.counters.spawned_particles > 0U) {
        const std::size_t particle_count_after_birth = context.state.particles.size();
        if (particle_count_after_birth < particle_count_before_birth ||
            particle_count_after_birth - particle_count_before_birth !=
                static_cast<std::size_t>(report.counters.spawned_particles)) {
          throw std::runtime_error(
              "source runtime star-formation report disagrees with appended particle rows");
        }
        // Gas cells are the hydro and gas-gravity mass authority; generic gas
        // particles are compatibility mirrors only. Refresh legacy mirrors once
        // after the batch so I/O/adapter consumers observe the post-birth mass.
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
      if (report.counters.spawned_particles > 0U && context.particle_scheduler != nullptr) {
        const std::uint64_t current_tick = context.particle_scheduler->currentTick();
        if (current_tick == std::numeric_limits<std::uint64_t>::max()) {
          throw std::overflow_error("source runtime newborn activation tick overflows uint64");
        }
        // Newborn stars have no retroactive kick. They join bin zero at the next
        // legal scheduler tick and become visible to the next force boundary.
        context.particle_scheduler->appendElements(
            static_cast<std::uint32_t>(report.counters.spawned_particles),
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
          m_seed_candidates,
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
    m_cell_volume_code.assign(cell_count, 0.0);
    m_owned_leaf_mask.assign(cell_count, 0U);
    m_feedback_candidate_cells.clear();
    m_feedback_candidate_cells.reserve(cell_count);
    std::unordered_set<std::uint64_t> patch_ids_with_children;
    patch_ids_with_children.reserve(state.patches.size());
    for (std::size_t patch_index = 0; patch_index < state.patches.size(); ++patch_index) {
      const std::uint64_t parent_id = state.patches.parent_patch_id[patch_index];
      if (parent_id != 0U) {
        patch_ids_with_children.insert(parent_id);
      }
    }
    for (std::uint32_t cell_index = 0; cell_index < cell_count; ++cell_index) {
      const PatchCellGeometry geometry = starFormationPatchCellGeometry(state, cell_index);
      bool owned_leaf = true;
      double volume = state.gas_cells.density_code[cell_index] > 0.0
          ? state.cells.mass_code[cell_index] /
                state.gas_cells.density_code[cell_index]
          : 0.0;
      if (geometry.valid) {
        volume = geometry.dx_stored * geometry.dy_stored * geometry.dz_stored;
        owned_leaf = state.patches.owning_rank[geometry.patch_index] == m_world_rank &&
            starFormationPatchIsLeaf(
                state, geometry.patch_index, patch_ids_with_children);
      }
      m_cell_volume_code[cell_index] = volume;
      m_owned_leaf_mask[cell_index] = owned_leaf ? 1U : 0U;
      if (owned_leaf) {
        m_feedback_candidate_cells.push_back(cell_index);
      }
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
    if (!m_stellar_evolution.config().enabled ||
        context.state.star_particles.size() == 0U) {
      return;
    }
    buildActiveStarRows(context);
    if (m_active_star_indices.empty()) {
      return;
    }
    buildOwnedLeafCellMetadata(context);
    const double elapsed_years = elapsedStellarEvolutionYears(context);
    const physics::StellarEvolutionStepReport evolution_report =
        m_stellar_evolution.evaluateElapsedYears(
            context.state, m_active_star_indices, elapsed_years);

    const std::size_t star_count = context.state.star_particles.size();
    m_returned_mass_delta_code.assign(star_count, 0.0);
    m_returned_metals_delta_code.assign(star_count, 0.0);
    m_feedback_energy_delta_erg.assign(star_count, 0.0);
    for (const physics::StellarEvolutionStarBudget& budget :
         evolution_report.budgets) {
      m_returned_mass_delta_code[budget.star_index] =
          budget.interval.returned_mass_code;
      m_returned_metals_delta_code[budget.star_index] =
          budget.interval.returned_metals_code;
      m_feedback_energy_delta_erg[budget.star_index] =
          budget.interval.feedback_energy_erg;
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
        .is_owned_leaf = m_owned_leaf_mask,
        .candidate_cell_indices = m_feedback_candidate_cells,
    };
    physics::StellarFeedbackDepositionView deposition_view{
        .cell_mass_code = context.state.cells.mass_code,
        .gas_density_code = context.state.gas_cells.density_code,
        .gas_internal_energy_code =
            context.state.gas_cells.internal_energy_code,
        .gas_metal_mass_code = context.state.gas_cells.metal_mass_code,
        .cell_volume_code = m_cell_volume_code,
    };
    (void)m_stellar_feedback.applyWithViews(
        context.state, m_stellar_feedback_state, geometry_view,
        deposition_view, m_active_star_indices, m_returned_mass_delta_code,
        m_returned_metals_delta_code, context.integrator_state.dt_time_code,
        m_feedback_energy_delta_erg);

    // The returned budget is now either deposited or durably attached to its
    // source star. Only after that transaction succeeds may stellar mass and
    // cumulative SSP bookkeeping advance.
    m_stellar_evolution.commitBudgets(context.state, evolution_report);
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
    buildOwnedLeafCellMetadata(context);
    const double scale_factor = m_coordinate_frame == core::CoordinateFrame::kComoving
        ? source_evaluation_scale_factor
        : 1.0;
    m_diffusion_cells.clear();
    m_diffusion_cells.resize(context.state.cells.size());
    for (std::uint32_t cell_index = 0;
         cell_index < context.state.cells.size(); ++cell_index) {
      const PatchCellGeometry geometry =
          starFormationPatchCellGeometry(context.state, cell_index);
      const double volume_stored = m_cell_volume_code[cell_index];
      const double volume_phys = volume_stored * scale_factor * scale_factor *
          scale_factor;
      physics::MetalDiffusionCell cell;
      cell.gas_cell_id = context.state.gas_cells.gas_cell_id[cell_index];
      cell.gas_mass_code = context.state.cells.mass_code[cell_index];
      cell.metal_mass_code =
          context.state.gas_cells.metal_mass_code[cell_index];
      cell.volume_code = volume_phys;
      cell.density_code = volume_phys > 0.0
          ? cell.gas_mass_code / volume_phys : 0.0;
      cell.filter_length_code = volume_phys > 0.0
          ? std::cbrt(volume_phys) : 0.0;
      cell.is_owned_leaf = m_owned_leaf_mask[cell_index] != 0U;
      if (geometry.valid && cell.is_owned_leaf) {
        const std::array<std::span<const double>, 3> velocity_fields{
            context.state.gas_cells.velocity_x_peculiar,
            context.state.gas_cells.velocity_y_peculiar,
            context.state.gas_cells.velocity_z_peculiar};
        const std::array<double, 3> spacing{
            geometry.dx_stored * scale_factor,
            geometry.dy_stored * scale_factor,
            geometry.dz_stored * scale_factor};
        for (std::size_t component = 0; component < velocity_fields.size(); ++component) {
          for (std::size_t axis = 0; axis < spacing.size(); ++axis) {
            cell.velocity_gradient.grad[component][axis] =
                starFormationDerivativeAtCell(
                    velocity_fields[component], geometry, static_cast<int>(axis),
                    spacing[axis], cell_index);
          }
        }
      }
      m_diffusion_cells[cell_index] = cell;
    }

    const internal::MetalDiffusionTopologyResult topology =
        internal::buildMetalDiffusionTopology(
            context.state,
            m_owned_leaf_mask,
            m_world_rank,
            scale_factor,
            m_hydro_boundary,
            m_diffusion_cells);
    m_diffusion_faces = topology.faces;
    if (m_diffusion_faces.empty()) {
      return;
    }
    const physics::MetalDiffusionStepReport report = m_metal_diffusion.advance(
        m_diffusion_cells, m_diffusion_faces,
        context.integrator_state.dt_time_code);
    for (std::uint32_t cell_index = 0;
         cell_index < context.state.cells.size(); ++cell_index) {
      if (m_diffusion_cells[cell_index].is_owned_leaf) {
        context.state.gas_cells.metal_mass_code[cell_index] =
            m_diffusion_cells[cell_index].metal_mass_code;
      }
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
    for (std::size_t i = 0; i < text.size(); ++i) {
      sidecar.payload[i] = static_cast<std::byte>(text[i]);
    }
    context.state.sidecars.upsert(std::move(sidecar));
  }

  void buildStarFormationInputs(
      const core::StepContext& context,
      std::span<const std::uint32_t> active_cells,
      double source_evaluation_scale_factor) {
    const core::SimulationState& state = context.state;
    m_star_formation_inputs.clear();
    m_star_formation_inputs.reserve(active_cells.size());

    std::unordered_set<std::uint64_t> patch_ids_with_children;
    patch_ids_with_children.reserve(state.patches.size());
    for (std::size_t patch_index = 0; patch_index < state.patches.size(); ++patch_index) {
      const std::uint64_t parent_id = state.patches.parent_patch_id[patch_index];
      if (parent_id != 0U) {
        patch_ids_with_children.insert(parent_id);
      }
    }

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
  std::uint32_t m_world_rank = 0;
  core::CoordinateFrame m_coordinate_frame = core::CoordinateFrame::kComoving;
  core::BoundaryCondition m_hydro_boundary = core::BoundaryCondition::kOpen;
  bool m_bh_enabled = false;
  bool m_bh_seeding_requested = false;
  bool m_is_cosmological = false;
  double m_mean_baryon_density0_code = 0.0;
  std::vector<std::uint32_t> m_full_cell_indices;
  std::vector<physics::StarFormationCellInput> m_star_formation_inputs;
  std::vector<std::uint32_t> m_active_star_indices;
  std::vector<double> m_returned_mass_delta_code;
  std::vector<double> m_returned_metals_delta_code;
  std::vector<double> m_feedback_energy_delta_erg;
  std::vector<double> m_cell_volume_code;
  std::vector<std::uint8_t> m_owned_leaf_mask;
  std::vector<std::uint32_t> m_feedback_candidate_cells;
  std::vector<physics::MetalDiffusionCell> m_diffusion_cells;
  std::vector<physics::MetalDiffusionFace> m_diffusion_faces;
  std::vector<physics::BlackHoleSeedCandidate> m_seed_candidates;
};

}  // namespace

std::unique_ptr<SourceRuntime> makeSourceRuntime(
    const core::SimulationConfig& config,
    const core::ModePolicy& mode_policy,
    const core::UnitSystem& units,
    std::uint32_t world_rank,
    const parallel::MpiContext& mpi_context) {
  return std::make_unique<SourceRuntimeImpl>(
      config, mode_policy, units, world_rank, mpi_context);
}

}  // namespace cosmosim::workflows
