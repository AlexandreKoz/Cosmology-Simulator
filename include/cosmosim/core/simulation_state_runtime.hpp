#pragma once

#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/core/memory_governor.hpp"

#include <memory>

namespace cosmosim::core {

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
  [[nodiscard]] bool isConsistent() const noexcept;
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
  [[nodiscard]] bool isConsistent() const noexcept;
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
  [[nodiscard]] bool isConsistent() const noexcept;
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
  [[nodiscard]] bool isConsistent() const noexcept;
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

  [[nodiscard]] bool isConsistent(std::size_t particle_count) const;
};

class ScratchAllocator {
 public:
  virtual ~ScratchAllocator() = default;
  [[nodiscard]] virtual std::byte* allocateBytes(std::size_t bytes, std::size_t alignment) = 0;
  virtual void reset() = 0;

  template <typename T>
  [[nodiscard]] T* allocateArray(std::size_t count) {
    static_assert(
        std::is_trivially_copyable_v<T> && std::is_trivially_destructible_v<T>,
        "ScratchAllocator.allocateArray supports only trivial scratch element types");
    const std::size_t bytes = checkedSizeMultiply(
        sizeof(T), count, "ScratchAllocator.allocateArray");
    auto* raw = allocateBytes(bytes, alignof(T));
    return reinterpret_cast<T*>(raw);
  }
};

class MonotonicScratchAllocator final : public ScratchAllocator {
 public:
  // Pointer-stable segmented monotonic arena for transient per-step scratch data.
  // Blocks never move after allocation; all returned pointers remain valid until
  // reset() begins the next scratch epoch or the arena is destroyed. reset()
  // reuses existing blocks instead of freeing them.
  explicit MonotonicScratchAllocator(std::size_t initial_capacity_bytes = 0);
  MonotonicScratchAllocator(
      MemoryGovernor* memory_governor,
      std::size_t initial_capacity_bytes = 0);

  MonotonicScratchAllocator(const MonotonicScratchAllocator&) = delete;
  MonotonicScratchAllocator& operator=(const MonotonicScratchAllocator&) = delete;
  MonotonicScratchAllocator(MonotonicScratchAllocator&&) noexcept = default;
  MonotonicScratchAllocator& operator=(MonotonicScratchAllocator&&) noexcept = default;

  [[nodiscard]] std::byte* allocateBytes(std::size_t bytes, std::size_t alignment) override;
  void reset() override;

  [[nodiscard]] std::size_t usedBytes() const noexcept;
  [[nodiscard]] std::size_t capacityBytes() const noexcept;
  [[nodiscard]] std::size_t logicalHighWaterBytes() const noexcept;
  [[nodiscard]] std::size_t capacityHighWaterBytes() const noexcept;
  [[nodiscard]] bool governed() const noexcept;

 private:
  struct Block {
    std::unique_ptr<std::byte[]> storage;
    std::size_t capacity_bytes = 0;
    std::size_t offset_bytes = 0;
    MemoryReservation reservation;
  };

  [[nodiscard]] Block& appendBlock(std::size_t minimum_capacity_bytes);

  MemoryGovernor* m_memory_governor = nullptr;
  std::vector<Block> m_blocks;
  std::size_t m_current_block = 0;
  std::size_t m_total_capacity_bytes = 0;
  std::size_t m_current_used_bytes = 0;
  std::size_t m_logical_high_water_bytes = 0;
  std::size_t m_capacity_high_water_bytes = 0;
  std::size_t m_next_block_capacity_bytes = 1024U;
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

  explicit TransientStepWorkspace(MemoryGovernor* memory_governor = nullptr)
      : scratch(memory_governor) {}

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
// No-copy all-active path used by the KDK orchestrator. Hot physics lanes alias
// authoritative ParticleSoa storage directly; only the uint32 scatter/index lane
// is materialized in workspace.
[[nodiscard]] GravityParticleKernelView buildGravityParticleKernelViewAllParticlesDirect(
    SimulationState& state,
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
