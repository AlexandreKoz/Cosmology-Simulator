#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

#include "cosmosim/core/soa_storage.hpp"

namespace cosmosim::core {

class HierarchicalTimeBinScheduler;

// Canonical species tags used in sidecar accounting and invariant checks.
enum class ParticleSpecies : std::uint8_t {
  kDarkMatter = 0,
  kGas = 1,
  kStar = 2,
  kBlackHole = 3,
  kTracer = 4,
};

struct ParticleSoa {
  // Authoritative owner for persistent gravity-hot particle truth in SimulationState.
  // Ownership lanes: position_*_comoving, velocity_*_peculiar, mass_code.
  // Deliberately no persistent acceleration lane: force accelerations are transient
  // kernel/workspace outputs and must not become a second persistent truth vector.
  AlignedVector<double> position_x_comoving;
  AlignedVector<double> position_y_comoving;
  AlignedVector<double> position_z_comoving;
  AlignedVector<double> velocity_x_peculiar;
  AlignedVector<double> velocity_y_peculiar;
  AlignedVector<double> velocity_z_peculiar;
  AlignedVector<double> mass_code;
  // Derived mirror of scheduler bin_index for diagnostics/I/O compatibility;
  // exact restart uses scheduler state.
  AlignedVector<std::uint8_t> time_bin;

  void resize(std::size_t count);
  [[nodiscard]] std::size_t size() const noexcept;
  [[nodiscard]] bool isConsistent() const noexcept;
};

struct ParticleSidecar {
  // Authoritative persistent metadata owner for particle IDs/species/rank and
  // ownership-adjacent particle metadata outside gravity-hot lanes.
  AlignedVector<std::uint64_t> particle_id;
  // Space-filling-curve key used for locality-preserving reorder/grouping.
  AlignedVector<std::uint64_t> sfc_key;
  // Authoritative species-tag ownership lane.
  AlignedVector<std::uint32_t> species_tag;
  AlignedVector<std::uint32_t> particle_flags;
  // Authoritative owning-rank lane for distributed ownership truth.
  AlignedVector<std::uint32_t> owning_rank;
  // Authoritative particle-position epoch used to predict inactive TreePM sources
  // without mutating persistent particle state during local active-bin refreshes.
  AlignedVector<double> last_drift_time_code;
  AlignedVector<double> last_drift_scale_factor;
  // Optional per-particle gravity softening value sidecar (comoving code units).
  // has_gravity_softening_override is the authoritative override mask. A populated
  // value vector without a populated mask is a materialized diagnostic/default mirror,
  // not independent per-particle override truth.
  AlignedVector<double> gravity_softening_comoving;
  AlignedVector<std::uint8_t> has_gravity_softening_override;

  void resize(std::size_t count);
  [[nodiscard]] std::size_t size() const noexcept;
  [[nodiscard]] bool isConsistent() const noexcept;
  [[nodiscard]] bool hasGravitySofteningOverride(std::size_t particle_index) const;
  [[nodiscard]] double gravitySofteningOverride(std::size_t particle_index) const;
  void setGravitySofteningOverride(std::size_t particle_index, double epsilon_comoving);
  void clearGravitySofteningOverride(std::size_t particle_index);
};

struct CellSoa {
  // Gravity-facing gas-cell skeleton kept separate from hydro thermodynamics.
  AlignedVector<double> center_x_comoving;
  AlignedVector<double> center_y_comoving;
  AlignedVector<double> center_z_comoving;
  AlignedVector<double> mass_code;
  // Derived mirror of scheduler bin_index for diagnostics/I/O compatibility;
  // not persistent authority.
  AlignedVector<std::uint8_t> time_bin;
  AlignedVector<std::uint32_t> patch_index;

  void resize(std::size_t count);
  [[nodiscard]] std::size_t size() const noexcept;
  [[nodiscard]] bool isConsistent() const noexcept;
};

struct GasCellSidecar {
  // Persistent gas-cell identity and thermodynamic state excluded from gravity-hot paths.
  // Reconstruction gradients are transient hydro/workspace scratch, not restart truth.
  // gas_cell_id is the stable nonzero identity mirror. parent_particle_id is legacy
  // lineage metadata: a value of zero mirrors a parentless GasCellIdentityRecord.
  // Cell-local velocity lanes are persistent hydro restart truth for parentless or
  // split cells; particle velocities are compatibility mirrors only when a unique
  // local parent particle is present.
  AlignedVector<std::uint64_t> gas_cell_id;
  AlignedVector<std::uint64_t> parent_particle_id;
  AlignedVector<double> velocity_x_peculiar;
  AlignedVector<double> velocity_y_peculiar;
  AlignedVector<double> velocity_z_peculiar;
  AlignedVector<double> density_code;
  AlignedVector<double> pressure_code;
  AlignedVector<double> internal_energy_code;
  AlignedVector<double> temperature_code;
  AlignedVector<double> sound_speed_code;

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

inline constexpr std::uint32_t kInvalidGasCellRow = std::numeric_limits<std::uint32_t>::max();

struct GasCellIdentityRecord {
  // Design-seam record for future AMR/moving-mesh gas ownership.
  // gas_cell_id is the stable identity and must not be derived from a particle row index.
  // A zero gas_cell_id is reserved as an uninitialized sentinel and is rejected by GasCellIdentityMap.
  std::uint64_t gas_cell_id = 0;
  std::optional<std::uint64_t> parent_particle_id;
  std::uint64_t owning_patch_id = 0;
  std::uint32_t local_cell_row = 0;
};

class GasCellIdentityMap {
 public:
  // Authoritative production row-space for gas-cell identity. Dense local rows remain transient
  // indexing; gas_cell_id is the stable key across reorder, restart, split/merge, and migration.
  // Lookup indices are rebuilt atomically on
  // assign(), and generation() changes after every successful mutation so future production views can
  // reject stale local-row mappings.
  void assign(std::vector<GasCellIdentityRecord> records);
  // Restart-authoritative import path: restore validated records together with
  // the generation value persisted at write time, so stale view checks see the
  // same identity epoch after reload. Normal runtime mutations must keep using
  // assign()/clear(), which bump the generation.
  void assignWithGeneration(std::vector<GasCellIdentityRecord> records, std::uint64_t generation);
  void clear() noexcept;
  [[nodiscard]] std::size_t size() const noexcept;
  [[nodiscard]] bool empty() const noexcept;
  [[nodiscard]] std::uint64_t generation() const noexcept;
  [[nodiscard]] bool isConsistent() const noexcept;
  [[nodiscard]] bool coversDenseLocalRows(std::size_t cell_count) const noexcept;
  void requireCoversDenseLocalRows(std::size_t cell_count, std::string_view caller) const;
  [[nodiscard]] std::span<const GasCellIdentityRecord> records() const noexcept;
  [[nodiscard]] const GasCellIdentityRecord* findByGasCellId(std::uint64_t gas_cell_id) const noexcept;
  [[nodiscard]] const GasCellIdentityRecord* findByLocalRow(std::uint32_t local_cell_row) const noexcept;
  [[nodiscard]] std::optional<std::uint32_t> rowForGasCellId(std::uint64_t gas_cell_id) const noexcept;
  [[nodiscard]] std::optional<std::uint64_t> gasCellIdForLocalRow(std::uint32_t local_cell_row) const noexcept;
  [[nodiscard]] std::optional<std::uint64_t> parentParticleIdForGasCellId(std::uint64_t gas_cell_id) const noexcept;
  [[nodiscard]] std::optional<std::uint64_t> owningPatchIdForGasCellId(std::uint64_t gas_cell_id) const noexcept;
  [[nodiscard]] std::vector<std::uint32_t> rowsForParentParticleId(std::uint64_t parent_particle_id) const;
  [[nodiscard]] std::vector<std::uint32_t> rowsForPatch(std::uint64_t owning_patch_id) const;

 private:
  bool rebuildLookupTables() noexcept;

  std::vector<GasCellIdentityRecord> m_records;
  std::unordered_map<std::uint64_t, std::size_t> m_index_by_gas_cell_id;
  std::unordered_map<std::uint32_t, std::size_t> m_index_by_local_row;
  std::uint64_t m_generation = 0;
};

// Build a dense new-row -> old-row remapping by stable gas_cell_id. New cells that are absent from
// old_map are marked with kInvalidGasCellRow so callers must deliberately initialize their hydro state.
[[nodiscard]] std::vector<std::uint32_t> buildGasCellNewToOldRowMap(
    const GasCellIdentityMap& old_map,
    const GasCellIdentityMap& new_map);

struct PatchSoa {
  // AMR patch descriptors and contiguous cell ranges [first_cell, first_cell + cell_count).
  AlignedVector<std::uint64_t> patch_id;
  AlignedVector<std::int32_t> level;
  AlignedVector<std::uint32_t> first_cell;
  AlignedVector<std::uint32_t> cell_count;
  // Restart-authoritative patch geometry. Zero extents/dimensions mark a legacy
  // non-AMR descriptor; production AMR hydro requires these lanes to be populated.
  AlignedVector<std::uint64_t> parent_patch_id;
  AlignedVector<std::uint64_t> morton_key;
  AlignedVector<double> origin_x_comoving;
  AlignedVector<double> origin_y_comoving;
  AlignedVector<double> origin_z_comoving;
  AlignedVector<double> extent_x_comoving;
  AlignedVector<double> extent_y_comoving;
  AlignedVector<double> extent_z_comoving;
  AlignedVector<std::uint16_t> cell_dim_x;
  AlignedVector<std::uint16_t> cell_dim_y;
  AlignedVector<std::uint16_t> cell_dim_z;
  // Authoritative AMR patch owner for distributed rebalancing. Patch payloads remain
  // contiguous in CellSoa; this lane is the ownership contract used by load-balance
  // planning and restart, not transient exchange scratch.
  AlignedVector<std::uint32_t> owning_rank;

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

enum class ModuleSidecarRequirementKind : std::uint8_t {
  kSparse = 0,
  kSpeciesMask = 1,
  kGasDensityAtLeast = 2,
  kBlackHoleAccretionAtLeast = 3,
  kParticleFlagMask = 4,
};

struct ModuleSidecarRequirement {
  // Structured coverage contract for particle-indexed module sidecars.
  // kSparse means optional rows.  kSpeciesMask reproduces the legacy mask.
  // Other predicates are evaluated against the authoritative particle/gas/BH
  // state during pack/commit so modules can require rows only for physically
  // active subsets instead of every particle of a species.
  ModuleSidecarRequirementKind kind = ModuleSidecarRequirementKind::kSparse;
  std::uint32_t species_mask = 0;
  std::uint32_t particle_flags_mask = 0;
  double threshold_code = 0.0;
};

struct ModuleSidecarBlock {
  // Opaque module payload with an independent schema version. Non-particle-indexed
  // blocks are run/module persistent state and are preserved as whole blocks.
  std::string module_name;
  std::uint32_t schema_version = 1;
  std::vector<std::byte> payload;

  // Optional particle-indexed row layout for module-owned persistent particle
  // sidecars. When enabled, each row is keyed by stable particle_id and the
  // corresponding byte slice must migrate with that particle.
  bool particle_indexed = false;
  std::uint32_t row_stride_bytes = 0;
  std::vector<std::uint64_t> particle_id_by_row;

  // Legacy unconditional species mask retained for compatibility. Structured
  // requirement is ORed with this mask when deciding whether a row is required.
  std::uint32_t required_species_mask = 0;
  ModuleSidecarRequirement requirement{};

  [[nodiscard]] bool isParticleIndexed() const noexcept;
  [[nodiscard]] bool requiresSpecies(ParticleSpecies species) const noexcept;
  [[nodiscard]] std::size_t rowCount() const noexcept;
  [[nodiscard]] std::span<const std::byte> rowPayload(std::size_t row) const;
};

struct ModuleParticleSidecarPayload {
  std::string module_name;
  std::uint32_t schema_version = 1;
  std::uint32_t row_stride_bytes = 0;
  std::uint32_t required_species_mask = 0;
  ModuleSidecarRequirement requirement{};
  std::vector<std::byte> payload;
};

class ModuleSidecarRegistry {
 public:
  // Insert or replace sidecar payload for a module.
  void upsert(ModuleSidecarBlock block);
  void clear() noexcept;
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
  // Derived timestep mirror copied from scheduler authority for transfer locality only.
  // Receiver-side scheduler remap/rebuild must validate or overwrite this lane
  // before exact continuation.
  AlignedVector<std::uint8_t> time_bin;
  AlignedVector<std::uint32_t> owning_rank;
  AlignedVector<double> last_drift_time_code;
  AlignedVector<double> last_drift_scale_factor;
};



struct PendingFluxRegisterRecord {
  std::uint64_t register_key = 0;
  std::uint64_t coarse_patch_id = 0;
  std::uint64_t coarse_gas_cell_id = 0;
  std::size_t coarse_cell_index = 0;
  std::uint8_t level = 0;
  std::uint8_t axis = 0;
  std::uint8_t orientation = 0;
  double expected_area_comov = 0.0;
  double coarse_area_accumulated_comov = 0.0;
  double fine_area_accumulated_comov = 0.0;
  double interval_start_code = 0.0;
  double interval_end_code = 0.0;
  double coarse_dt_code = 0.0;
  std::uint32_t expected_fine_substeps = 1;
  std::uint32_t completed_fine_substeps = 0;
  std::uint64_t fine_substep_coverage_mask = 0;
  std::uint32_t coarse_face_count = 0;
  std::uint32_t fine_face_count = 0;
  std::uint64_t gas_cell_identity_generation = 0;
  std::uint64_t patch_geometry_generation = 0;
  double coarse_mass_flux_integral_code = 0.0;
  double coarse_momentum_x_flux_integral_code = 0.0;
  double coarse_momentum_y_flux_integral_code = 0.0;
  double coarse_momentum_z_flux_integral_code = 0.0;
  double coarse_total_energy_flux_integral_code = 0.0;
  double fine_mass_flux_integral_code = 0.0;
  double fine_momentum_x_flux_integral_code = 0.0;
  double fine_momentum_y_flux_integral_code = 0.0;
  double fine_momentum_z_flux_integral_code = 0.0;
  double fine_total_energy_flux_integral_code = 0.0;

  [[nodiscard]] bool hasCoarseContribution() const noexcept { return coarse_face_count > 0U; }
  [[nodiscard]] bool hasCompleteFineSubstepCoverage() const noexcept {
    return expected_fine_substeps > 0U && completed_fine_substeps >= expected_fine_substeps;
  }
  [[nodiscard]] bool isComplete() const noexcept {
    return hasCoarseContribution() && hasCompleteFineSubstepCoverage();
  }
};

class PendingFluxRegisterStore {
 public:
  [[nodiscard]] bool empty() const noexcept { return m_records.empty(); }
  [[nodiscard]] std::size_t size() const noexcept { return m_records.size(); }
  [[nodiscard]] std::span<const PendingFluxRegisterRecord> records() const noexcept { return m_records; }
  [[nodiscard]] std::span<PendingFluxRegisterRecord> mutableRecords() noexcept { return m_records; }
  void clear() noexcept { m_records.clear(); }
  void assign(std::vector<PendingFluxRegisterRecord> records) { m_records = std::move(records); }
  [[nodiscard]] const PendingFluxRegisterRecord* findByRegisterKey(std::uint64_t register_key) const noexcept {
    for (const auto& record : m_records) {
      if (record.register_key == register_key) {
        return &record;
      }
    }
    return nullptr;
  }
  [[nodiscard]] PendingFluxRegisterRecord* findByRegisterKey(std::uint64_t register_key) noexcept {
    for (auto& record : m_records) {
      if (record.register_key == register_key) {
        return &record;
      }
    }
    return nullptr;
  }
  PendingFluxRegisterRecord& upsertByRegisterKey(const PendingFluxRegisterRecord& record) {
    if (auto* existing = findByRegisterKey(record.register_key); existing != nullptr) {
      return *existing;
    }
    m_records.push_back(record);
    return m_records.back();
  }
  void eraseCompletedByKey(std::span<const std::uint64_t> register_keys) {
    m_records.erase(
        std::remove_if(
            m_records.begin(),
            m_records.end(),
            [register_keys](const PendingFluxRegisterRecord& record) {
              return std::find(register_keys.begin(), register_keys.end(), record.register_key) != register_keys.end();
            }),
        m_records.end());
  }

 private:
  std::vector<PendingFluxRegisterRecord> m_records;
};

struct GasCellMigrationFields {
  std::uint64_t gas_cell_id = 0;
  std::uint8_t has_parent_particle = 0;
  std::uint64_t parent_particle_id = 0;
  std::uint64_t owning_patch_id = 0;
  std::uint32_t destination_local_cell_row = kInvalidGasCellRow;
  std::uint64_t gas_cell_identity_generation = 0;
  std::uint64_t ghost_hydro_epoch = 0;
  double center_x_comoving = 0.0;
  double center_y_comoving = 0.0;
  double center_z_comoving = 0.0;
  double cell_mass_code = 0.0;
  std::uint8_t cell_time_bin = 0;
  std::uint32_t patch_index = 0;
  double velocity_x_peculiar = 0.0;
  double velocity_y_peculiar = 0.0;
  double velocity_z_peculiar = 0.0;
  double density_code = 0.0;
  double pressure_code = 0.0;
  double internal_energy_code = 0.0;
  double temperature_code = 0.0;
  double sound_speed_code = 0.0;
};

struct GasCellMigrationRecord {
  // Standalone gas-cell ownership payload keyed by gas_cell_id. This is separate
  // from ParticleMigrationRecord so parentless and split/merged cells can migrate
  // without inventing a particle parent.
  std::uint32_t owning_rank = 0;
  GasCellMigrationFields fields{};
};

struct GasCellStaleGhostRecord {
  std::uint64_t gas_cell_id = 0;
  std::uint32_t local_cell_row = kInvalidGasCellRow;
  std::uint64_t gas_cell_identity_generation = 0;
  std::uint64_t ghost_hydro_epoch = 0;
};

struct GasCellMigrationCommit {
  int world_rank = 0;
  std::uint64_t expected_gas_cell_identity_generation = 0;
  std::uint64_t expected_ghost_hydro_epoch = 0;
  std::vector<std::uint32_t> outbound_local_cell_rows;
  std::vector<GasCellMigrationRecord> inbound_records;
  std::vector<GasCellStaleGhostRecord> stale_local_ghost_records;
};

struct AmrPatchMigrationFields {
  std::uint64_t patch_id = 0;
  std::uint64_t parent_patch_id = 0;
  std::uint64_t morton_key = 0;
  std::int32_t level = 0;
  std::uint32_t owning_rank = 0;
  std::uint32_t cell_count = 0;
  double origin_x_comoving = 0.0;
  double origin_y_comoving = 0.0;
  double origin_z_comoving = 0.0;
  double extent_x_comoving = 0.0;
  double extent_y_comoving = 0.0;
  double extent_z_comoving = 0.0;
  std::uint16_t cell_dim_x = 0;
  std::uint16_t cell_dim_y = 0;
  std::uint16_t cell_dim_z = 0;
};

struct AmrPatchMigrationRecord {
  // Atomic patch ownership payload. Patch descriptors and every authoritative
  // gas-cell row in the patch migrate together, keyed by stable patch_id and
  // gas_cell_id rather than by old local rows.
  AmrPatchMigrationFields patch{};
  std::vector<GasCellMigrationRecord> gas_cell_records;
};

struct AmrPatchMigrationCommit {
  int world_rank = 0;
  std::uint64_t expected_gas_cell_identity_generation = 0;
  std::uint64_t expected_ghost_hydro_epoch = 0;
  std::vector<std::uint32_t> outbound_local_patch_indices;
  std::vector<AmrPatchMigrationRecord> inbound_records;
  std::vector<GasCellStaleGhostRecord> stale_local_ghost_records;
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

struct SchedulerMigrationFields {
  // Scheduler-authoritative lanes for an element keyed by stable particle_id.
  // time_bin in ParticleSoa remains a derived mirror; these fields are the
  // migration payload needed to rebuild exact hierarchical-timestep state.
  std::uint8_t bin_index = 0;
  std::uint64_t next_activation_tick = 0;
  std::uint8_t pending_bin_index = 0xFF;
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
  // Derived timestep mirror copied with the migrating particle for diagnostics
  // and legacy transfer consumers. Exact continuation must use scheduler_fields
  // when has_scheduler_fields is true.
  std::uint8_t time_bin = 0;
  bool has_scheduler_fields = false;
  SchedulerMigrationFields scheduler_fields{};
  double last_drift_time_code = 0.0;
  double last_drift_scale_factor = 1.0;
  bool has_gravity_softening_value = false;
  bool has_gravity_softening_override = false;
  double gravity_softening_comoving = 0.0;
  bool has_gas_cell_fields = false;
  GasCellMigrationFields gas_cell_fields{};
  bool has_star_fields = false;
  StarParticleMigrationFields star_fields{};
  bool has_black_hole_fields = false;
  BlackHoleParticleMigrationFields black_hole_fields{};
  bool has_tracer_fields = false;
  TracerParticleMigrationFields tracer_fields{};
  std::vector<ModuleParticleSidecarPayload> module_sidecar_payloads;
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
  // Authoritative in-memory gas-cell identity map. During the H2 legacy transition,
  // gas_cells.gas_cell_id and gas_cells.parent_particle_id are compatibility mirrors
  // synchronized from this map for particle-bound import/restart paths.
  GasCellIdentityMap gas_cell_identity;
  PatchSoa patches;
  // Restart-authoritative AMR hydro deferred synchronization state.
  PendingFluxRegisterStore pending_flux_registers;
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
  void refreshGasCellIdentityFromParticleOrder();
  [[nodiscard]] bool gasCellIdentityMatchesParticleOrder() const;
  void refreshGasCellIdentityMapFromParticleBoundState();
  void refreshGasCellIdentityMapFromSidecarLanes();
  [[nodiscard]] bool gasCellIdentityMapMatchesParticleBoundState() const;
  [[nodiscard]] std::uint64_t gasCellIdentityGeneration() const noexcept;
  void requireGasCellIdentityMapCoversDenseRows(std::string_view caller) const;
  void requireGasCellIdentityMapFresh(std::uint64_t expected_generation, std::string_view caller) const;
  [[nodiscard]] std::optional<std::uint32_t> rowForGasCellId(std::uint64_t gas_cell_id) const noexcept;
  [[nodiscard]] std::optional<std::uint64_t> gasCellIdForLocalRow(std::uint32_t local_cell_row) const noexcept;
  [[nodiscard]] std::optional<std::uint64_t> parentParticleIdForGasCellId(std::uint64_t gas_cell_id) const noexcept;
  [[nodiscard]] std::optional<std::uint64_t> owningPatchIdForGasCellId(std::uint64_t gas_cell_id) const noexcept;
  [[nodiscard]] bool gasCellIdentityMapMatchesSidecarLanes() const;
  void requireParticleBoundGasCellContract(std::string_view caller) const;
  void legacyRequireParticleBoundGasCellContract(std::string_view caller) const;
  [[nodiscard]] std::optional<std::uint64_t> parentParticleIdForGasCellRow(std::uint32_t cell_index) const;
  [[nodiscard]] std::uint32_t gasCellRowForParticleId(std::uint64_t particle_id) const;
  [[nodiscard]] std::uint32_t gasParticleIndexForCellRow(std::uint32_t cell_index) const;

  [[nodiscard]] ParticleTransferPacket packSpeciesTransferPacket(ParticleSpecies species_tag) const;
  [[nodiscard]] std::vector<ParticleMigrationRecord> packParticleMigrationRecords(
      std::span<const std::uint32_t> local_indices,
      const HierarchicalTimeBinScheduler& scheduler) const;
  void commitParticleMigration(const ParticleMigrationCommit& commit);
  [[nodiscard]] std::vector<GasCellMigrationRecord> packGasCellMigrationRecords(
      std::span<const std::uint32_t> local_cell_rows,
      std::uint64_t ghost_hydro_epoch = 0) const;
  void commitGasCellMigration(const GasCellMigrationCommit& commit);
  [[nodiscard]] std::vector<AmrPatchMigrationRecord> packAmrPatchMigrationRecords(
      std::span<const std::uint32_t> local_patch_indices,
      std::uint64_t ghost_hydro_epoch = 0) const;
  void commitAmrPatchMigration(const AmrPatchMigrationCommit& commit);
  [[nodiscard]] std::uint64_t particleIndexGeneration() const noexcept;
  [[nodiscard]] std::uint64_t cellIndexGeneration() const noexcept;
  void bumpParticleIndexGeneration() noexcept;
  void bumpCellIndexGeneration() noexcept;

 private:
  [[nodiscard]] std::vector<ParticleMigrationRecord> packParticleMigrationRecordsCore(
      std::span<const std::uint32_t> local_indices) const;

  std::uint64_t m_particle_index_generation = 0;
  std::uint64_t m_cell_index_generation = 0;
};

template <typename T>
inline constexpr bool k_is_canonical_particle_state_owner_v = false;

template <>
inline constexpr bool k_is_canonical_particle_state_owner_v<SimulationState> = true;

struct ActiveIndexSet {
  // Per-step compact active index lists assembled by scheduler/rung logic.
  std::vector<std::uint32_t> particle_indices;
  std::vector<std::uint32_t> cell_indices;

  void clear();
};

struct ParticleActiveView {
  // Compact contiguous particle spans materialized in the transient workspace.
  // Read-only active views are transient mirrors only: they never own persistent truth
  // and consumers can compare the captured generation with SimulationState before
  // using a view across any structural transform.
  std::span<const std::uint64_t> particle_id;
  std::span<const std::uint32_t> species_tag;
  std::span<const double> position_x_comoving;
  std::span<const double> position_y_comoving;
  std::span<const double> position_z_comoving;
  std::span<const double> velocity_x_peculiar;
  std::span<const double> velocity_y_peculiar;
  std::span<const double> velocity_z_peculiar;
  std::span<const double> mass_code;
  std::uint64_t source_particle_index_generation = 0;

  [[nodiscard]] std::size_t size() const noexcept;
};

struct CellActiveView {
  // Compact contiguous cell spans materialized in the transient workspace.
  // Read-only active views are transient mirrors only: they never own persistent truth
  // and consumers can compare the captured generation with SimulationState before
  // using a view across any structural transform or gas-cell rebuild.
  std::span<const double> center_x_comoving;
  std::span<const double> center_y_comoving;
  std::span<const double> center_z_comoving;
  std::span<const double> mass_code;
  std::span<const std::uint32_t> patch_index;
  std::span<const double> density_code;
  std::span<const double> pressure_code;
  std::uint64_t source_cell_index_generation = 0;

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
  std::uint64_t source_particle_index_generation = 0;

  [[nodiscard]] std::size_t size() const noexcept;
};

struct HydroCellKernelView {
  // Compact read/write cell hydro view for active hydrodynamics kernels.
  // HOT-FIELD CONTRACT (review-critical):
  //   allowed mutable lanes = {center_[xyz]_comoving, mass_code, density_code, pressure_code}
  //   allowed identity/scatter lanes = {gas_cell_id, local_cell_row}
  // No patch descriptors, thermodynamic cold metadata, reconstruction
  // gradients, or provenance fields may be added to this hot view.
  std::span<std::uint32_t> cell_index;
  std::span<std::uint64_t> gas_cell_id;
  std::span<std::uint32_t> local_cell_row;
  std::span<double> center_x_comoving;
  std::span<double> center_y_comoving;
  std::span<double> center_z_comoving;
  std::span<double> mass_code;
  std::span<double> density_code;
  std::span<double> pressure_code;
  std::uint64_t source_cell_index_generation = 0;
  std::uint64_t source_gas_cell_identity_generation = 0;

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
  // Compact particle active-set buffers. These are reusable transient mirrors,
  // not persistent owners; clear() drops sizes but intentionally preserves
  // capacity so repeated hot-loop materialization does not churn allocations.
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

  // Compact read/write hydro kernel buffers. These carry only allowed hot lanes
  // in HydroCellKernelView and are scattered only after cell-generation checks.
  AlignedVector<std::uint32_t> hydro_cell_index;
  AlignedVector<std::uint64_t> hydro_gas_cell_id;
  AlignedVector<std::uint32_t> hydro_local_cell_row;
  AlignedVector<double> hydro_cell_center_x_comoving;
  AlignedVector<double> hydro_cell_center_y_comoving;
  AlignedVector<double> hydro_cell_center_z_comoving;
  AlignedVector<double> hydro_cell_mass_code;
  AlignedVector<double> hydro_cell_density_code;
  AlignedVector<double> hydro_cell_pressure_code;

  // Compact cell active-set buffers. These are reusable transient mirrors,
  // not persistent owners; patch/identity metadata remains outside mutable hot views.
  AlignedVector<double> cell_center_x_comoving;
  AlignedVector<double> cell_center_y_comoving;
  AlignedVector<double> cell_center_z_comoving;
  AlignedVector<double> cell_mass_code;
  AlignedVector<std::uint32_t> cell_patch_index;
  AlignedVector<double> cell_density_code;
  AlignedVector<double> cell_pressure_code;

  // Transient hydro reconstruction-gradient scratch. These lanes are derived
  // inside hydro reconstruction/source staging and are never restart truth.
  AlignedVector<double> hydro_recon_gradient_x;
  AlignedVector<double> hydro_recon_gradient_y;
  AlignedVector<double> hydro_recon_gradient_z;

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
[[nodiscard]] GravityParticleKernelView buildGravityParticleKernelViewAllParticles(
    const SimulationState& state,
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
// Debug guard: require exactly one row for each star/BH/tracer particle and no rows on other species.
void debugAssertSpeciesSidecarOwnershipInvariants(const SimulationState& state);
// Rebuild/check the temporary gas identity lanes from canonical gas-particle ordering.
void refreshGasCellIdentityFromParticleOrder(SimulationState& state);
[[nodiscard]] bool gasCellIdentityMatchesParticleOrder(const SimulationState& state);
void refreshGasCellIdentityMapFromParticleBoundState(SimulationState& state);
void refreshGasCellIdentityMapFromSidecarLanes(SimulationState& state);
[[nodiscard]] bool gasCellIdentityMapMatchesParticleBoundState(const SimulationState& state);
[[nodiscard]] bool gasCellIdentityMapMatchesSidecarLanes(const SimulationState& state);
void requireGasCellIdentityMapCoversDenseRows(const SimulationState& state, std::string_view caller);
void requireGasCellIdentityMapFresh(
    const SimulationState& state,
    std::uint64_t expected_generation,
    std::string_view caller);
// Legacy/import compatibility contract only: older particle-bound inputs require
// one local gas particle per gas-cell row before positional remaps may consume them.
void legacyRequireParticleBoundGasCellContract(const SimulationState& state, std::string_view caller);
// Backward-compatible wrapper for older tests/import callers. Production hydro
// paths should validate GasCellIdentityMap coverage and use optional parent lookup.
void requireParticleBoundGasCellContract(const SimulationState& state, std::string_view caller);
[[nodiscard]] std::optional<std::uint64_t> parentParticleIdForGasCellRow(
    const SimulationState& state,
    std::uint32_t cell_index);
[[nodiscard]] std::uint32_t gasCellRowForParticleId(const SimulationState& state, std::uint64_t particle_id);
[[nodiscard]] std::uint32_t gasParticleIndexForCellRow(const SimulationState& state, std::uint32_t cell_index);
// Backward-compatible debug alias for legacy particle-bound checks; production
// gas-cell paths should validate GasCellIdentityMap coverage instead.
void debugAssertGasCellIdentityContract(const SimulationState& state);

}  // namespace cosmosim::core
