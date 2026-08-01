#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <new>
#include <span>
#include <string_view>
#include <vector>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/physics/star_formation.hpp"

namespace {

struct AllocationHeader {
  void* base = nullptr;
  std::size_t size = 0;
  bool tracked = false;
};

struct AllocationSnapshot {
  std::size_t count = 0;
  std::size_t bytes = 0;
  std::size_t peak_live_bytes = 0;
};

std::atomic<bool> g_track_allocations{false};
std::atomic<std::size_t> g_allocation_count{0};
std::atomic<std::size_t> g_allocation_bytes{0};
std::atomic<std::size_t> g_live_bytes{0};
std::atomic<std::size_t> g_peak_live_bytes{0};

void updatePeak(std::size_t value) noexcept {
  std::size_t previous = g_peak_live_bytes.load(std::memory_order_relaxed);
  while (value > previous &&
         !g_peak_live_bytes.compare_exchange_weak(
             previous, value, std::memory_order_relaxed, std::memory_order_relaxed)) {
  }
}

void* allocateTracked(std::size_t size, std::size_t alignment) {
  if (size == 0) {
    size = 1;
  }
  alignment = std::max(alignment, alignof(AllocationHeader));
  if ((alignment & (alignment - 1U)) != 0U) {
    throw std::bad_alloc();
  }
  if (size > std::numeric_limits<std::size_t>::max() - alignment - sizeof(AllocationHeader)) {
    throw std::bad_alloc();
  }

  const std::size_t allocation_size = size + alignment + sizeof(AllocationHeader);
  void* base = std::malloc(allocation_size);
  if (base == nullptr) {
    throw std::bad_alloc();
  }
  const auto begin = reinterpret_cast<std::uintptr_t>(base) + sizeof(AllocationHeader);
  const auto aligned = (begin + alignment - 1U) & ~(static_cast<std::uintptr_t>(alignment) - 1U);
  auto* header = reinterpret_cast<AllocationHeader*>(aligned) - 1;
  header->base = base;
  header->size = size;
  header->tracked = g_track_allocations.load(std::memory_order_relaxed);
  if (header->tracked) {
    g_allocation_count.fetch_add(1U, std::memory_order_relaxed);
    g_allocation_bytes.fetch_add(size, std::memory_order_relaxed);
    const std::size_t live = g_live_bytes.fetch_add(size, std::memory_order_relaxed) + size;
    updatePeak(live);
  }
  return reinterpret_cast<void*>(aligned);
}

void deallocateTracked(void* pointer) noexcept {
  if (pointer == nullptr) {
    return;
  }
  auto* header = reinterpret_cast<AllocationHeader*>(pointer) - 1;
  if (header->tracked) {
    g_live_bytes.fetch_sub(header->size, std::memory_order_relaxed);
  }
  std::free(header->base);
}

void beginAllocationTracking() noexcept {
  g_allocation_count.store(0, std::memory_order_relaxed);
  g_allocation_bytes.store(0, std::memory_order_relaxed);
  g_live_bytes.store(0, std::memory_order_relaxed);
  g_peak_live_bytes.store(0, std::memory_order_relaxed);
  g_track_allocations.store(true, std::memory_order_relaxed);
}

AllocationSnapshot endAllocationTracking() noexcept {
  g_track_allocations.store(false, std::memory_order_relaxed);
  return AllocationSnapshot{
      .count = g_allocation_count.load(std::memory_order_relaxed),
      .bytes = g_allocation_bytes.load(std::memory_order_relaxed),
      .peak_live_bytes = g_peak_live_bytes.load(std::memory_order_relaxed),
  };
}

struct ScenarioData {
  cosmosim::core::SimulationState state;
  std::vector<cosmosim::physics::StarFormationCellInput> inputs;
};

[[nodiscard]] cosmosim::physics::StarFormationConfig makeConfig() {
  cosmosim::physics::StarFormationConfig config;
  config.model = cosmosim::core::StarFormationModelKind::kAdaptiveBoundJeans;
  config.epsilon_ff = 1.0;
  config.bound_alpha_vir_max = 1.0;
  config.require_converging_flow = true;
  config.collapse_timescale =
      cosmosim::core::StarFormationCollapseTimescale::kMinimumFreeFallOrCompression;
  config.jeans_mass_floor_code = 0.1;
  config.star_particle_mass_policy = cosmosim::core::StarParticleMassPolicy::kFixed;
  config.target_star_particle_mass_code = 0.5;
  config.min_star_particle_mass_code = 0.1;
  config.max_star_particle_mass_code = 2.5;
  config.max_spawn_particles_per_cell_step = 8;
  config.max_fractional_mass_conversion = 0.25;
  config.min_remaining_gas_fraction = 0.01;
  config.min_remaining_gas_mass_code = 0.0;
  config.stochastic_spawning = false;
  config.random_seed = 0x51f07a11ULL;
  config.newton_g_code = 1.0;
  config.density_is_comoving = false;
  config.geometry_is_comoving = false;
  return config;
}

[[nodiscard]] ScenarioData makeScenario(
    std::size_t cell_count,
    std::size_t eligible_stride) {
  ScenarioData data;
  data.state.resizeCells(cell_count);
  data.inputs.resize(cell_count);
  std::vector<cosmosim::core::GasCellIdentityRecord> identity_records;
  identity_records.reserve(cell_count);

  for (std::size_t row = 0; row < cell_count; ++row) {
    const bool eligible = eligible_stride != 0U && row % eligible_stride == 0U;
    const std::uint64_t gas_cell_id = 1000000ULL + static_cast<std::uint64_t>(row);
    auto& input = data.inputs[row];
    input.cell_index = static_cast<std::uint32_t>(row);
    input.gas_cell_id = gas_cell_id;
    input.owning_rank = 0U;
    input.is_active = true;
    input.is_owned = true;
    input.is_leaf = true;
    input.is_ghost = false;
    input.gas_mass_code = 10.0;
    input.gas_density_code = 100.0;
    input.cell_volume_code = 0.1;
    input.gas_temperature_k = 100.0;
    input.gas_sound_speed_code = eligible ? 0.01 : 100.0;
    input.velocity_x_peculiar = 2.0;
    input.velocity_y_peculiar = -3.0;
    input.velocity_z_peculiar = 4.0;
    input.velocity_divergence_code = -20.0;
    input.velocity_gradient_frobenius_sq_code = eligible ? 0.0 : 1.0e6;
    input.gas_metal_mass_code = 0.2;
    input.center_x_comoving = static_cast<double>(row);
    input.center_y_comoving = 0.0;
    input.center_z_comoving = 0.0;

    data.state.cells.center_x_comoving[row] = input.center_x_comoving;
    data.state.cells.center_y_comoving[row] = input.center_y_comoving;
    data.state.cells.center_z_comoving[row] = input.center_z_comoving;
    data.state.cells.mass_code[row] = input.gas_mass_code;
    data.state.gas_cells.gas_cell_id[row] = gas_cell_id;
    data.state.gas_cells.density_code[row] = input.gas_density_code;
    data.state.gas_cells.pressure_code[row] = 5.0;
    data.state.gas_cells.internal_energy_code[row] = 2.0;
    data.state.gas_cells.temperature_code[row] = input.gas_temperature_k;
    data.state.gas_cells.sound_speed_code[row] = input.gas_sound_speed_code;
    data.state.gas_cells.velocity_x_peculiar[row] = input.velocity_x_peculiar;
    data.state.gas_cells.velocity_y_peculiar[row] = input.velocity_y_peculiar;
    data.state.gas_cells.velocity_z_peculiar[row] = input.velocity_z_peculiar;
    data.state.gas_cells.metal_mass_code[row] = input.gas_metal_mass_code;
    identity_records.push_back({
        .gas_cell_id = gas_cell_id,
        .parent_particle_id = std::nullopt,
        .owning_patch_id = 0U,
        .local_cell_row = static_cast<std::uint32_t>(row),
    });
  }
  data.state.replaceGasCellIdentityRecords(std::move(identity_records));
  return data;
}

void runScenario(
    std::string_view phase,
    std::size_t cell_count,
    std::size_t eligible_stride,
    std::uint64_t integration_tick) {
  auto data = makeScenario(cell_count, eligible_stride);
  const cosmosim::physics::StarFormationModel model(makeConfig());

  beginAllocationTracking();
  const auto start = std::chrono::steady_clock::now();
  const auto report = model.applyFromInputs(
      data.state,
      data.inputs,
      0.003,
      1.0,
      integration_tick);
  const auto stop = std::chrono::steady_clock::now();
  const AllocationSnapshot allocation = endAllocationTracking();

  const double elapsed_s = std::max(
      std::chrono::duration<double>(stop - start).count(), 1.0e-12);
  const double cells_per_second = static_cast<double>(cell_count) / elapsed_s;
  const double births_per_second =
      static_cast<double>(report.counters.spawned_particles) / elapsed_s;

  std::cout << "bench_star_formation_spawn"
            << " phase=" << phase
            << " build_type=" << COSMOSIM_BUILD_TYPE
            << " hardware=cpu"
            << " threads=1"
            << " model=adaptive_bound_jeans"
            << " cells=" << cell_count
            << " elapsed_s=" << elapsed_s
            << " cells_per_second=" << cells_per_second
            << " births_per_second=" << births_per_second
            << " eligible_cells=" << report.counters.eligible_cells
            << " spawned_particles=" << report.counters.spawned_particles
            << " realized_stellar_mass_code=" << report.counters.spawned_mass_code
            << " bytes_allocated_per_step=" << allocation.bytes
            << " allocation_count_per_step=" << allocation.count
            << " peak_temporary_memory_upper_bound_bytes=" << allocation.peak_live_bytes
            << '\n';
}

}  // namespace

void* operator new(std::size_t size) { return allocateTracked(size, alignof(std::max_align_t)); }
void* operator new[](std::size_t size) { return allocateTracked(size, alignof(std::max_align_t)); }
void* operator new(std::size_t size, std::align_val_t alignment) {
  return allocateTracked(size, static_cast<std::size_t>(alignment));
}
void* operator new[](std::size_t size, std::align_val_t alignment) {
  return allocateTracked(size, static_cast<std::size_t>(alignment));
}
void* operator new(std::size_t size, const std::nothrow_t&) noexcept {
  try {
    return allocateTracked(size, alignof(std::max_align_t));
  } catch (...) {
    return nullptr;
  }
}
void* operator new[](std::size_t size, const std::nothrow_t&) noexcept {
  try {
    return allocateTracked(size, alignof(std::max_align_t));
  } catch (...) {
    return nullptr;
  }
}
void operator delete(void* pointer) noexcept { deallocateTracked(pointer); }
void operator delete[](void* pointer) noexcept { deallocateTracked(pointer); }
void operator delete(void* pointer, std::size_t) noexcept { deallocateTracked(pointer); }
void operator delete[](void* pointer, std::size_t) noexcept { deallocateTracked(pointer); }
void operator delete(void* pointer, std::align_val_t) noexcept { deallocateTracked(pointer); }
void operator delete[](void* pointer, std::align_val_t) noexcept { deallocateTracked(pointer); }
void operator delete(void* pointer, std::size_t, std::align_val_t) noexcept { deallocateTracked(pointer); }
void operator delete[](void* pointer, std::size_t, std::align_val_t) noexcept { deallocateTracked(pointer); }
void operator delete(void* pointer, const std::nothrow_t&) noexcept { deallocateTracked(pointer); }
void operator delete[](void* pointer, const std::nothrow_t&) noexcept { deallocateTracked(pointer); }

int main() {
  constexpr std::size_t k_cells = 1U << 14U;
  runScenario("eligibility_no_births", k_cells, 0U, 100U);
  runScenario("eligibility_sparse_births", k_cells, 1024U, 101U);
  runScenario("dense_birth_plan_and_batch_append", k_cells, 1U, 102U);
  return 0;
}
