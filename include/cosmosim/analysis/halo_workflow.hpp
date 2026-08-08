#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <limits>
#include <string>
#include <string_view>
#include <span>
#include <vector>

#include "cosmosim/core/config.hpp"
#include "cosmosim/core/simulation_state.hpp"

namespace cosmosim::parallel {
class MpiContext;
}

namespace cosmosim::analysis {

constexpr std::uint64_t k_unbound_halo_id = std::numeric_limits<std::uint64_t>::max();

enum class HaloWorkflowMode : std::uint8_t {
  kPostProcess = 0,
  kOnTheFly = 1,
};

struct HaloCatalogEntry {
  std::uint64_t halo_id = 0;
  std::uint64_t snapshot_step_index = 0;
  double snapshot_scale_factor = 1.0;
  std::uint64_t particle_count = 0;
  double total_mass_code = 0.0;
  std::array<double, 3> center_of_mass_comov{};
  std::array<double, 3> bulk_velocity_peculiar{};
  std::uint64_t min_particle_id = 0;
};

struct SubhaloCandidateEntry {
  std::uint64_t subhalo_id = 0;
  std::uint64_t host_halo_id = 0;
  std::uint32_t rank_in_host = 0;
  std::uint64_t particle_count = 0;
  double bound_mass_code = 0.0;
};

struct MergerTreeNodePlan {
  std::uint64_t tree_node_id = 0;
  std::uint64_t halo_id = 0;
  std::uint64_t snapshot_step_index = 0;
  double snapshot_scale_factor = 1.0;
  std::uint64_t descendant_tree_node_id = 0;
  bool has_descendant = false;
};

struct HaloCatalog {
  std::uint32_t schema_version = 1;
  std::uint64_t snapshot_step_index = 0;
  double snapshot_scale_factor = 1.0;
  std::uint64_t normalized_config_hash = 0;
  std::string run_name;
  std::string halo_finder = "fof";
  std::vector<HaloCatalogEntry> halos;
  std::vector<std::uint64_t> halo_id_by_particle;
  std::vector<SubhaloCandidateEntry> subhalo_candidates;
};

struct MergerTreePlan {
  std::uint32_t schema_version = 1;
  std::string run_name;
  std::uint64_t normalized_config_hash = 0;
  std::vector<MergerTreeNodePlan> nodes;
};

enum class FofNeighborSearch : std::uint8_t {
  kSpatialHash = 0,
  kAllPairsReference = 1,
};

struct FofConfig {
  double linking_length_factor_mean_interparticle = 0.2;
  std::uint64_t min_group_size = 16;
  bool include_gas = true;
  bool include_stars = true;
  bool include_black_holes = true;
  FofNeighborSearch neighbor_search = FofNeighborSearch::kSpatialHash;
};

struct FofProfilingCounters {
  std::uint64_t candidate_particle_count = 0;
  std::uint64_t pair_checks = 0;
  std::uint64_t pair_links = 0;
  std::uint64_t occupied_spatial_cells = 0;
  FofNeighborSearch neighbor_search = FofNeighborSearch::kSpatialHash;
  double linking_length_comov = 0.0;
};

struct HaloParticleView {
  std::span<const double> position_x_comoving;
  std::span<const double> position_y_comoving;
  std::span<const double> position_z_comoving;
  std::span<const double> velocity_x_peculiar;
  std::span<const double> velocity_y_peculiar;
  std::span<const double> velocity_z_peculiar;
  std::span<const double> mass_code;
  std::span<const std::uint32_t> species_tag;
  std::span<const std::uint64_t> particle_id;
  std::uint64_t normalized_config_hash = 0;
};

[[nodiscard]] HaloParticleView buildHaloParticleView(const core::SimulationState& state);

struct DistributedHaloCatalogResult {
  // Full global catalog exists on root_rank. Other ranks receive local stable-ID
  // membership labels without replicating the full global particle payload.
  HaloCatalog root_catalog;
  std::vector<std::uint64_t> local_halo_id_by_particle;
  int root_rank = 0;
  bool has_global_catalog = false;
};

struct HaloWorkflowReport {
  HaloCatalog catalog;
  MergerTreePlan tree_plan;
  // In distributed mode this contains membership labels for the caller's local
  // particle ordering on every rank; the full catalog/tree are root-owned.
  std::vector<std::uint64_t> local_halo_id_by_particle;
  bool has_global_catalog = true;
  FofProfilingCounters profiling;
  std::filesystem::path halo_catalog_path;
  std::filesystem::path merger_tree_plan_path;
};

class HaloCatalogSchema {
 public:
  [[nodiscard]] static std::uint32_t schemaVersion() noexcept;
  [[nodiscard]] static std::string_view catalogFormatName() noexcept;
  [[nodiscard]] static std::string_view mergerTreePlanFormatName() noexcept;
  [[nodiscard]] static std::vector<std::string_view> haloFields();
  [[nodiscard]] static std::vector<std::string_view> subhaloFields();
  [[nodiscard]] static std::vector<std::string_view> mergerTreeFields();
};

// Production FOF uses a periodic spatial hash/cell-linked neighbor search. The
// all-pairs implementation remains available explicitly as a small-N numerical oracle.
class FofHaloFinder {
 public:
  explicit FofHaloFinder(FofConfig config);

  [[nodiscard]] HaloCatalog buildCatalogFromView(
      const HaloParticleView& particles,
      const core::SimulationConfig& config,
      std::uint64_t snapshot_step_index,
      double snapshot_scale_factor,
      FofProfilingCounters* profiling = nullptr) const;

  [[nodiscard]] HaloCatalog buildCatalogReferenceAllPairsFromView(
      const HaloParticleView& particles,
      const core::SimulationConfig& config,
      std::uint64_t snapshot_step_index,
      double snapshot_scale_factor,
      FofProfilingCounters* profiling = nullptr) const;

  [[nodiscard]] HaloCatalog buildCatalog(
      const core::SimulationState& state,
      const core::SimulationConfig& config,
      std::uint64_t snapshot_step_index,
      double snapshot_scale_factor,
      FofProfilingCounters* profiling = nullptr) const;

  // Correctness-first MPI ownership path. Local particles are packed to root,
  // grouped once with the same production spatial hash, and stable membership
  // labels are broadcast back by particle ID. This closes rank-boundary halo
  // splitting without returning to all-pairs candidate search.
  [[nodiscard]] DistributedHaloCatalogResult buildDistributedCatalogFromView(
      const HaloParticleView& local_particles,
      const core::SimulationConfig& config,
      const parallel::MpiContext& mpi_context,
      std::uint64_t snapshot_step_index,
      double snapshot_scale_factor,
      int root_rank = 0,
      FofProfilingCounters* root_profiling = nullptr) const;

 private:
  [[nodiscard]] bool includeSpecies(core::ParticleSpecies species) const noexcept;

  FofConfig m_config;
};

class MergerTreePlanner {
 public:
  [[nodiscard]] MergerTreePlan buildPlan(
      const HaloCatalog& current,
      const HaloCatalog* previous) const;
};

class HaloWorkflowPlanner {
 public:
  explicit HaloWorkflowPlanner(core::SimulationConfig config);

  [[nodiscard]] static std::string_view ownershipBoundary();

  [[nodiscard]] HaloWorkflowMode workflowMode() const;

  [[nodiscard]] HaloWorkflowReport runSnapshotWorkflow(
      const core::SimulationState& state,
      std::uint64_t snapshot_step_index,
      double snapshot_scale_factor,
      const HaloCatalog* previous_catalog = nullptr) const;

  [[nodiscard]] HaloWorkflowReport runSnapshotWorkflowDistributed(
      const core::SimulationState& local_state,
      const parallel::MpiContext& mpi_context,
      std::uint64_t snapshot_step_index,
      double snapshot_scale_factor,
      const HaloCatalog* previous_root_catalog = nullptr,
      int root_rank = 0) const;

  void writeHaloCatalog(const HaloCatalog& catalog, const std::filesystem::path& path) const;
  void writeMergerTreePlan(const MergerTreePlan& tree_plan, const std::filesystem::path& path) const;

 private:
  [[nodiscard]] std::filesystem::path analysisDirectory() const;
  [[nodiscard]] std::filesystem::path haloCatalogPath(std::uint64_t snapshot_step_index) const;
  [[nodiscard]] std::filesystem::path mergerTreePlanPath(std::uint64_t snapshot_step_index) const;

  [[nodiscard]] static FofConfig fofConfigFromSimulationConfig(const core::SimulationConfig& config);

  core::SimulationConfig m_config;
  FofHaloFinder m_fof;
  MergerTreePlanner m_tree_planner;
};

}  // namespace cosmosim::analysis
