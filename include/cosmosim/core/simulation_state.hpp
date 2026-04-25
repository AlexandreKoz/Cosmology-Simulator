#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

#include "cosmosim/core/soa_storage.hpp"

namespace cosmosim::core {

// Canonical species tags used in sidecar accounting and invariant checks.
enum class ParticleSpecies : std::uint8_t {
  kDarkMatter = 0,
  kGas = 1,
  kStar = 2,
  kBlackHole = 3,
  kTracer = 4,
};

struct ParticleSoa {
  // Shared gravity-hot particle fields; species-specific cold data must live in sidecars.
  AlignedVector<double> position_x_comoving;
  AlignedVector<double> position_y_comoving;
  AlignedVector<double> position_z_comoving;
  AlignedVector<double> velocity_x_peculiar;
  AlignedVector<double> velocity_y_peculiar;
  AlignedVector<double> velocity_z_peculiar;
  AlignedVector<double> mass_code;
  AlignedVector<std::uint8_t> time_bin;

  void resize(std::size_t count);
  [[nodiscard]] std::size_t size() const noexcept;
  [[nodiscard]] bool isConsistent() const noexcept;
};

struct ParticleSidecar {
  // Shared metadata sidecar: IDs, species ownership, and rank ownership.
  AlignedVector<std::uint64_t> particle_id;
  // Space-filling-curve key used for locality-preserving reorder/grouping.
  AlignedVector<std::uint64_t> sfc_key;
  AlignedVector<std::uint32_t> species_tag;
  AlignedVector<std::uint32_t> particle_flags;
  AlignedVector<std::uint32_t> owning_rank;
  // Optional per-particle gravity softening sidecar (comoving code units).
  // Empty means "use species/global policy only".
  AlignedVector<double> gravity_softening_comoving;

  void resize(std::size_t count);
  [[nodiscard]] std::size_t size() const noexcept;
  [[nodiscard]] bool isConsistent() const noexcept;
};

struct CellSoa {
  // Gravity-facing gas-cell skeleton kept separate from hydro thermodynamics.
  AlignedVector<double> center_x_comoving;
  AlignedVector<double> center_y_comoving;
  AlignedVector<double> center_z_comoving;
  AlignedVector<double> mass_code;
  AlignedVector<std::uint8_t> time_bin;
  AlignedVector<std::uint32_t> patch_index;

  void resize(std::size_t count);
  [[nodiscard]] std::size_t size() const noexcept;
  [[nodiscard]] bool isConsistent() const noexcept;
};

struct GasCellSidecar {
  // Gas-cell thermodynamic and reconstruction state excluded from gravity-hot paths.
  AlignedVector<double> density_code;
  AlignedVector<double> pressure_code;
  AlignedVector<double> internal_energy_code;
  AlignedVector<double> temperature_code;
  AlignedVector<double> sound_speed_code;
  AlignedVector<double> recon_gradient_x;
  AlignedVector<double> recon_gradient_y;
  AlignedVector<double> recon_gradient_z;

  void resize(std::size_t count);
  [[nodiscard]] std::size_t size() const noexcept;
  [[nodiscard]] bool isConsistent() const noexcept;
};

struct StarParticleSidecar {
  // Stellar-formation metadata decoupled from common particle skeleton.
  AlignedVector<std::uint32_t> particle_index;
  AlignedVector<double> formation_scale_factor;
  AlignedVector<double> birth_mass_code;
  AlignedVector<double> metallicity_mass_fraction;
  // Explicit stellar evolution bookkeeping lanes for auditable enrichment.
  AlignedVector<double> stellar_age_years_last;
  AlignedVector<double> stellar_returned_mass_cumulative_code;
  AlignedVector<double> stellar_returned_metals_cumulative_code;
  AlignedVector<double> stellar_feedback_energy_cumulative_erg;
  std::array<AlignedVector<double>, 3> stellar_returned_mass_channel_cumulative_code;
  std::array<AlignedVector<double>, 3> stellar_returned_metals_channel_cumulative_code;
  std::array<AlignedVector<double>, 3> stellar_feedback_energy_channel_cumulative_erg;

  void resize(std::size_t count);
  [[nodiscard]] std::size_t size() const noexcept;
  [[nodiscard]] bool isConsistent() const noexcept;
};

struct BlackHoleParticleSidecar {
  // Black-hole subgrid metadata sidecar.
  AlignedVector<std::uint32_t> particle_index;
  AlignedVector<std::uint32_t> host_cell_index;
  AlignedVector<double> subgrid_mass_code;
  AlignedVector<double> accretion_rate_code;
  AlignedVector<double> feedback_energy_code;
  AlignedVector<double> eddington_ratio;
  AlignedVector<double> cumulative_accreted_mass_code;
  AlignedVector<double> cumulative_feedback_energy_code;
  AlignedVector<double> duty_cycle_active_time_code;
  AlignedVector<double> duty_cycle_total_time_code;

  void resize(std::size_t count);
  [[nodiscard]] std::size_t size() const noexcept;
  [[nodiscard]] bool isConsistent() const noexcept;
};

struct TracerParticleSidecar {
  // Tracer attachment metadata sidecar.
  AlignedVector<std::uint32_t> particle_index;
  AlignedVector<std::uint64_t> parent_particle_id;
  AlignedVector<std::uint64_t> injection_step;
  AlignedVector<std::uint32_t> host_cell_index;
  AlignedVector<double> mass_fraction_of_host;
  AlignedVector<double> last_host_mass_code;
  AlignedVector<double> cumulative_exchanged_mass_code;

  void resize(std::size_t count);
  [[nodiscard]] std::size_t size() const noexcept;
  [[nodiscard]] bool isConsistent() const noexcept;
};

struct PatchSoa {
  // AMR patch descriptors and contiguous cell ranges [first_cell, first_cell + cell_count).
  AlignedVector<std::uint64_t> patch_id;
  AlignedVector<std::int32_t> level;
  AlignedVector<std::uint32_t> first_cell;
  AlignedVector<std::uint32_t> cell_count;

  void resize(std::size_t count);
  [[nodiscard]] std::size_t size() const noexcept;
  [[nodiscard]] bool isConsistent() const noexcept;
};

struct SpeciesContainer {
  // Explicit species counts; used as an auditable ownership ledger.
  std::array<std::uint64_t, 5> count_by_species{};

  [[nodiscard]] std::uint64_t totalCount() const noexcept;
  [[nodiscard]] bool isConsistentWith(const ParticleSidecar& sidecar) const noexcept;
};

struct ParticleSpeciesIndex {
  // Explicit species-local to global particle index mapping.
  std::array<AlignedVector<std::uint32_t>, 5> global_index_by_species;
  AlignedVector<std::uint32_t> local_index_by_global;

  void rebuild(const ParticleSidecar& sidecar);
  [[nodiscard]] std::size_t count(ParticleSpecies species) const noexcept;
  [[nodiscard]] std::span<const std::uint32_t> globalIndices(ParticleSpecies species) const noexcept;
  [[nodiscard]] std::uint32_t localIndex(std::uint32_t global_index) const;
  [[nodiscard]] std::uint32_t globalIndex(ParticleSpecies species, std::uint32_t local_index) const;
};

struct StateMetadata {
  // Schema/provenance fields that must remain stable across restart/snapshot workflows.
  std::uint32_t schema_version = 2;
  std::string run_name = "cosmosim_run";
  std::uint64_t normalized_config_hash = 0;
  std::string normalized_config_hash_hex;
  std::uint64_t step_index = 0;
  double scale_factor = 1.0;
  std::string snapshot_stem = "snapshot";
  std::string restart_stem = "restart";

  [[nodiscard]] std::string serialize() const;
  [[nodiscard]] static StateMetadata deserialize(std::string_view text);
};

struct ModuleSidecarBlock {
  // Opaque module payload with an independent schema version.
  std::string module_name;
  std::uint32_t schema_version = 1;
  std::vector<std::byte> payload;
};

class ModuleSidecarRegistry {
 public:
  // Insert or replace sidecar payload for a module.
  void upsert(ModuleSidecarBlock block);
  [[nodiscard]] const ModuleSidecarBlock* find(std::string_view module_name) const;
  [[nodiscard]] std::size_t size() const noexcept;
  [[nodiscard]] std::vector<const ModuleSidecarBlock*> blocksSortedByName() const;

 private:
  std::unordered_map<std::string, ModuleSidecarBlock> m_sidecars;
};

struct ParticleTransferPacket {
  // Species-specific transfer packet for MPI or host-device staging.
  ParticleSpecies species = ParticleSpecies::kDarkMatter;
  AlignedVector<std::uint64_t> particle_id;
  AlignedVector<double> position_x_comoving;
  AlignedVector<double> position_y_comoving;
  AlignedVector<double> position_z_comoving;
  AlignedVector<double> velocity_x_peculiar;
  AlignedVector<double> velocity_y_peculiar;
  AlignedVector<double> velocity_z_peculiar;
  AlignedVector<double> mass_code;
  AlignedVector<std::uint8_t> time_bin;
  AlignedVector<std::uint32_t> owning_rank;
};

struct StarParticleMigrationFields {
  double formation_scale_factor = 0.0;
  double birth_mass_code = 0.0;
  double metallicity_mass_fraction = 0.0;
  double stellar_age_years_last = 0.0;
  double stellar_returned_mass_cumulative_code = 0.0;
  double stellar_returned_metals_cumulative_code = 0.0;
  double stellar_feedback_energy_cumulative_erg = 0.0;
  std::array<double, 3> stellar_returned_mass_channel_cumulative_code{};
  std::array<double, 3> stellar_returned_metals_channel_cumulative_code{};
  std::array<double, 3> stellar_feedback_energy_channel_cumulative_erg{};
};

struct BlackHoleParticleMigrationFields {
  std::uint32_t host_cell_index = 0;
  double subgrid_mass_code = 0.0;
  double accretion_rate_code = 0.0;
  double feedback_energy_code = 0.0;
  double eddington_ratio = 0.0;
  double cumulative_accreted_mass_code = 0.0;
  double cumulative_feedback_energy_code = 0.0;
  double duty_cycle_active_time_code = 0.0;
  double duty_cycle_total_time_code = 0.0;
};

struct TracerParticleMigrationFields {
  std::uint64_t parent_particle_id = 0;
  std::uint64_t injection_step = 0;
  std::uint32_t host_cell_index = 0;
  double mass_fraction_of_host = 0.0;
  double last_host_mass_code = 0.0;
  double cumulative_exchanged_mass_code = 0.0;
};

struct ParticleMigrationRecord {
  std::uint64_t particle_id = 0;
  std::uint64_t sfc_key = 0;
  std::uint32_t species_tag = 0;
  std::uint32_t particle_flags = 0;
  std::uint32_t owning_rank = 0;
  double position_x_comoving = 0.0;
  double position_y_comoving = 0.0;
  double position_z_comoving = 0.0;
  double velocity_x_peculiar = 0.0;
  double velocity_y_peculiar = 0.0;
  double velocity_z_peculiar = 0.0;
  double mass_code = 0.0;
  std::uint8_t time_bin = 0;
  bool has_star_fields = false;
  StarParticleMigrationFields star_fields{};
  bool has_black_hole_fields = false;
  BlackHoleParticleMigrationFields black_hole_fields{};
  bool has_tracer_fields = false;
  TracerParticleMigrationFields tracer_fields{};
};

struct ParticleMigrationCommit {
  int world_rank = 0;
  std::vector<std::uint32_t> outbound_local_indices;
  std::vector<ParticleMigrationRecord> inbound_records;
  std::vector<std::uint32_t> stale_local_ghost_indices;
};

class SimulationState {
 public:
  // Single ownership root for persistent run state.
  ParticleSoa particles;
  ParticleSidecar particle_sidecar;
  CellSoa cells;
  GasCellSidecar gas_cells;
  PatchSoa patches;
  SpeciesContainer species;
  ParticleSpeciesIndex particle_species_index;
  StarParticleSidecar star_particles;
  BlackHoleParticleSidecar black_holes;
  TracerParticleSidecar tracers;
  StateMetadata metadata;
  ModuleSidecarRegistry sidecars;

  void resizeParticles(std::size_t count);
  void resizeCells(std::size_t count);
  void resizePatches(std::size_t count);

  [[nodiscard]] bool validateOwnershipInvariants() const;
  [[nodiscard]] bool validateUniqueParticleIds() const;
  void rebuildSpeciesIndex();

  [[nodiscard]] ParticleTransferPacket packSpeciesTransferPacket(ParticleSpecies species_tag) const;
  [[nodiscard]] std::vector<ParticleMigrationRecord> packParticleMigrationRecords(
      std::span<const std::uint32_t> local_indices) const;
  void commitParticleMigration(const ParticleMigrationCommit& commit);
  [[nodiscard]] std::uint64_t particleIndexGeneration() const noexcept;
  [[nodiscard]] std::uint64_t cellIndexGeneration() const noexcept;
  void bumpParticleIndexGeneration() noexcept;
  void bumpCellIndexGeneration() noexcept;

 private:
  std::uint64_t m_particle_index_generation = 0;
  std::uint64_t m_cell_index_generation = 0;
};

struct ActiveIndexSet {
  // Per-step compact active index lists assembled by scheduler/rung logic.
  std::vector<std::uint32_t> particle_indices;
  std::vector<std::uint32_t> cell_indices;

  void clear();
};

struct ParticleActiveView {
  // Compact contiguous particle spans materialized in the transient workspace.
  std::span<const std::uint64_t> particle_id;
  std::span<const std::uint32_t> species_tag;
  std::span<const double> position_x_comoving;
  std::span<const double> position_y_comoving;
  std::span<const double> position_z_comoving;
  std::span<const double> velocity_x_peculiar;
  std::span<const double> velocity_y_peculiar;
  std::span<const double> velocity_z_peculiar;
  std::span<const double> mass_code;

  [[nodiscard]] std::size_t size() const noexcept;
};

struct CellActiveView {
  // Compact contiguous cell spans materialized in the transient workspace.
  std::span<const double> center_x_comoving;
  std::span<const double> center_y_comoving;
  std::span<const double> center_z_comoving;
  std::span<const double> mass_code;
  std::span<const std::uint32_t> patch_index;
  std::span<const double> density_code;
  std::span<const double> pressure_code;

  [[nodiscard]] std::size_t size() const noexcept;
};

struct GravityParticleKernelView {
  // Compact read/write particle hot view for gravity kernels.
  // HOT-FIELD CONTRACT (review-critical):
  //   allowed mutable physics lanes = {position_[xyz]_comoving, velocity_[xyz]_peculiar, mass_code}
  //   allowed index lane = {particle_index}
  // No IDs, species tags, owning rank, flags, or any metadata/provenance fields
  // may be added to this view.
  std::span<std::uint32_t> particle_index;
  std::span<double> position_x_comoving;
  std::span<double> position_y_comoving;
  std::span<double> position_z_comoving;
  std::span<double> velocity_x_peculiar;
  std::span<double> velocity_y_peculiar;
  std::span<double> velocity_z_peculiar;
  std::span<double> mass_code;

  [[nodiscard]] std::size_t size() const noexcept;
};

struct HydroCellKernelView {
  // Compact read/write cell hydro view for active hydrodynamics kernels.
  // HOT-FIELD CONTRACT (review-critical):
  //   allowed mutable lanes = {center_[xyz]_comoving, mass_code, density_code, pressure_code}
  //   allowed index lane = {cell_index}
  // No patch descriptors, thermodynamic cold metadata, reconstruction
  // gradients, or provenance fields may be added to this hot view.
  std::span<std::uint32_t> cell_index;
  std::span<double> center_x_comoving;
  std::span<double> center_y_comoving;
  std::span<double> center_z_comoving;
  std::span<double> mass_code;
  std::span<double> density_code;
  std::span<double> pressure_code;

  [[nodiscard]] std::size_t size() const noexcept;
};

enum class ParticleReorderMode : std::uint8_t {
  // Stable grouping by scheduler rung / hierarchical time step bin.
  kByTimeBin = 0,
  // Stable grouping by SFC key for locality-friendly traversals.
  kBySfcKey = 1,
  // Stable grouping by species for species-specialized loops/packing.
  kBySpecies = 2,
};

enum class SidecarSyncMode : std::uint8_t {
  // Sidecar rows are physically permuted with parent particle rows.
  kMoveWithParent = 0,
  // Sidecar rows stay in-place and particle_index is remapped through old->new.
  kUseParentIndirection = 1,
};

struct SidecarSyncPolicy {
  // Species sidecars are index-based by default and remap through old->new.
  SidecarSyncMode star_particles = SidecarSyncMode::kUseParentIndirection;
  SidecarSyncMode black_holes = SidecarSyncMode::kUseParentIndirection;
  SidecarSyncMode tracers = SidecarSyncMode::kUseParentIndirection;
};

struct ParticleReorderMap {
  // Explicit old/new particle index mapping used for auditable sidecar sync.
  std::vector<std::uint32_t> old_to_new_index;
  std::vector<std::uint32_t> new_to_old_index;

  [[nodiscard]] bool isConsistent(std::size_t particle_count) const noexcept;
};

class ScratchAllocator {
 public:
  virtual ~ScratchAllocator() = default;
  [[nodiscard]] virtual std::byte* allocateBytes(std::size_t bytes, std::size_t alignment) = 0;
  virtual void reset() = 0;

  template <typename T>
  [[nodiscard]] T* allocateArray(std::size_t count) {
    auto* raw = allocateBytes(sizeof(T) * count, alignof(T));
    return reinterpret_cast<T*>(raw);
  }
};

class MonotonicScratchAllocator final : public ScratchAllocator {
 public:
  // Growable monotonic byte arena for transient per-step scratch data.
  explicit MonotonicScratchAllocator(std::size_t initial_capacity_bytes = 0);

  [[nodiscard]] std::byte* allocateBytes(std::size_t bytes, std::size_t alignment) override;
  void reset() override;

  [[nodiscard]] std::size_t capacityBytes() const noexcept;

 private:
  std::vector<std::byte> m_storage;
  std::size_t m_offset_bytes = 0;
};

struct TransientStepWorkspace {
  // Compact particle active-set buffers.
  AlignedVector<std::uint64_t> particle_id;
  AlignedVector<std::uint32_t> particle_species_tag;
  AlignedVector<double> particle_position_x_comoving;
  AlignedVector<double> particle_position_y_comoving;
  AlignedVector<double> particle_position_z_comoving;
  AlignedVector<double> particle_velocity_x_peculiar;
  AlignedVector<double> particle_velocity_y_peculiar;
  AlignedVector<double> particle_velocity_z_peculiar;
  AlignedVector<double> particle_mass_code;
  AlignedVector<std::uint32_t> gravity_particle_index;

  // Compact read/write hydro kernel buffers.
  AlignedVector<std::uint32_t> hydro_cell_index;
  AlignedVector<double> hydro_cell_center_x_comoving;
  AlignedVector<double> hydro_cell_center_y_comoving;
  AlignedVector<double> hydro_cell_center_z_comoving;
  AlignedVector<double> hydro_cell_mass_code;
  AlignedVector<double> hydro_cell_density_code;
  AlignedVector<double> hydro_cell_pressure_code;

  // Compact cell active-set buffers.
  AlignedVector<double> cell_center_x_comoving;
  AlignedVector<double> cell_center_y_comoving;
  AlignedVector<double> cell_center_z_comoving;
  AlignedVector<double> cell_mass_code;
  AlignedVector<std::uint32_t> cell_patch_index;
  AlignedVector<double> cell_density_code;
  AlignedVector<double> cell_pressure_code;

  // Monotonic scratch arena reused between steps via reset().
  MonotonicScratchAllocator scratch;

  void clear();
};

[[nodiscard]] ParticleActiveView buildParticleActiveView(
    const SimulationState& state,
    std::span<const std::uint32_t> active_particle_indices,
    TransientStepWorkspace& workspace);

[[nodiscard]] CellActiveView buildCellActiveView(
    const SimulationState& state,
    std::span<const std::uint32_t> active_cell_indices,
    TransientStepWorkspace& workspace);

[[nodiscard]] GravityParticleKernelView buildGravityParticleKernelView(
    const SimulationState& state,
    std::span<const std::uint32_t> active_particle_indices,
    TransientStepWorkspace& workspace);

void scatterGravityParticleKernelView(
    const GravityParticleKernelView& view,
    SimulationState& state);

[[nodiscard]] HydroCellKernelView buildHydroCellKernelView(
    const SimulationState& state,
    std::span<const std::uint32_t> active_cell_indices,
    TransientStepWorkspace& workspace);

void scatterHydroCellKernelView(
    const HydroCellKernelView& view,
    SimulationState& state);

// Build stable old/new index maps from the selected particle ordering key.
[[nodiscard]] ParticleReorderMap buildParticleReorderMap(
    const SimulationState& state,
    ParticleReorderMode mode);

// Apply one auditable particle permutation and synchronize all sidecars.
void reorderParticles(
    SimulationState& state,
    const ParticleReorderMap& reorder_map,
    const SidecarSyncPolicy& sync_policy = {});

// Debug guard: throw on any species sidecar index no longer owned by particles.
void debugAssertNoStaleParticleIndices(const SimulationState& state);

}  // namespace cosmosim::core
