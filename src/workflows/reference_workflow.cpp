#include "cosmosim/workflows/reference_workflow.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <iomanip>
#include <optional>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <cstring>
#include <type_traits>
#include <utility>
#include <vector>

#include "cosmosim/analysis/diagnostics.hpp"
#include "cosmosim/core/build_config.hpp"
#include "cosmosim/core/cosmology.hpp"
#include "cosmosim/core/cuda_runtime.hpp"
#include "cosmosim/core/profiling.hpp"
#include "cosmosim/core/simulation_mode.hpp"
#include "cosmosim/core/time_integration.hpp"
#include "cosmosim/gravity/tree_pm_coupling.hpp"
#include "cosmosim/hydro/hydro_core_solver.hpp"
#include "cosmosim/hydro/hydro_reconstruction.hpp"
#include "cosmosim/hydro/hydro_riemann.hpp"
#include "cosmosim/io/ic_reader.hpp"
#include "cosmosim/io/restart_checkpoint.hpp"
#include "cosmosim/io/snapshot_hdf5.hpp"
#include "cosmosim/physics/black_hole_agn.hpp"
#include "cosmosim/physics/star_formation.hpp"

namespace cosmosim::workflows {
namespace {

constexpr double k_gamma_adiabatic = 5.0 / 3.0;
constexpr double k_pressure_floor = 1.0e-10;
constexpr double k_density_floor = 1.0e-10;
constexpr std::size_t k_default_generated_particle_axis = 6;

[[nodiscard]] gravity::PmAssignmentScheme toPmAssignmentScheme(
    core::TreePmAssignmentScheme assignment_scheme) {
  switch (assignment_scheme) {
    case core::TreePmAssignmentScheme::kCic:
      return gravity::PmAssignmentScheme::kCic;
    case core::TreePmAssignmentScheme::kTsc:
      return gravity::PmAssignmentScheme::kTsc;
  }
  throw std::runtime_error("unhandled TreePm assignment scheme enum value");
}

[[nodiscard]] std::string treePmAssignmentSchemeName(core::TreePmAssignmentScheme assignment_scheme) {
  switch (assignment_scheme) {
    case core::TreePmAssignmentScheme::kCic:
      return "cic";
    case core::TreePmAssignmentScheme::kTsc:
      return "tsc";
  }
  throw std::runtime_error("unhandled TreePm assignment scheme enum value");
}



gravity::TreePmCoordinator makeRuntimeAwareTreePmCoordinator(
    const core::SimulationConfig& config,
    const gravity::PmGridShape& pm_grid_shape) {
  parallel::MpiContext mpi_context;
  mpi_context.validateExpectedWorldSizeOrThrow(config.parallel.mpi_ranks_expected);
  const auto layout = parallel::makePmSlabLayout(
      pm_grid_shape.nx,
      pm_grid_shape.ny,
      pm_grid_shape.nz,
      mpi_context.worldSize(),
      mpi_context.worldRank());
  return gravity::TreePmCoordinator(pm_grid_shape, layout);
}

[[nodiscard]] std::string pmDecompositionModeName(core::PmDecompositionMode mode) {
  switch (mode) {
    case core::PmDecompositionMode::kSlab:
      return "slab";
    case core::PmDecompositionMode::kPencil:
      return "pencil";
  }
  throw std::runtime_error("unhandled PM decomposition mode enum value");
}


struct GasCellMigrationRecord {
  std::uint64_t particle_id = 0;
  double center_x_comoving = 0.0;
  double center_y_comoving = 0.0;
  double center_z_comoving = 0.0;
  double mass_code = 0.0;
  std::uint8_t time_bin = 0;
  std::uint32_t patch_index = 0;
  double density_code = 0.0;
  double pressure_code = 0.0;
  double internal_energy_code = 0.0;
  double temperature_code = 0.0;
  double sound_speed_code = 0.0;
  double recon_gradient_x = 0.0;
  double recon_gradient_y = 0.0;
  double recon_gradient_z = 0.0;
};

static_assert(std::is_trivially_copyable_v<GasCellMigrationRecord>);

[[nodiscard]] std::string formatIndexedRankedFileStem(
    std::string_view stem,
    std::uint64_t index,
    int world_size,
    int world_rank) {
  std::ostringstream out;
  out << stem << '_' << std::setw(3) << std::setfill('0') << index;
  if (world_size > 1) {
    out << "_rank" << std::setw(3) << std::setfill('0') << world_rank;
  }
  out << ".hdf5";
  return out.str();
}

[[nodiscard]] std::unordered_map<std::uint64_t, GasCellMigrationRecord> collectLocalGasCellRecords(
    const core::SimulationState& state,
    std::span<const std::uint32_t> gas_particle_indices) {
  std::unordered_map<std::uint64_t, GasCellMigrationRecord> records;
  const auto gas_globals = state.particle_species_index.globalIndices(core::ParticleSpecies::kGas);
  if (state.cells.size() != gas_globals.size()) {
    throw std::runtime_error("gas cell count must match local gas particle count when collecting migration records");
  }
  std::vector<std::uint8_t> keep_mask(state.particles.size(), 0U);
  for (const std::uint32_t particle_index : gas_particle_indices) {
    if (particle_index >= state.particles.size()) {
      throw std::out_of_range("gas particle index out of range while collecting gas migration records");
    }
    keep_mask[particle_index] = 1U;
  }
  for (std::size_t cell_index = 0; cell_index < gas_globals.size(); ++cell_index) {
    const std::uint32_t particle_index = gas_globals[cell_index];
    if (keep_mask[particle_index] == 0U) {
      continue;
    }
    GasCellMigrationRecord record;
    record.particle_id = state.particle_sidecar.particle_id[particle_index];
    record.center_x_comoving = state.cells.center_x_comoving[cell_index];
    record.center_y_comoving = state.cells.center_y_comoving[cell_index];
    record.center_z_comoving = state.cells.center_z_comoving[cell_index];
    record.mass_code = state.cells.mass_code[cell_index];
    record.time_bin = state.cells.time_bin[cell_index];
    record.patch_index = state.cells.patch_index[cell_index];
    record.density_code = state.gas_cells.density_code[cell_index];
    record.pressure_code = state.gas_cells.pressure_code[cell_index];
    record.internal_energy_code = state.gas_cells.internal_energy_code[cell_index];
    record.temperature_code = state.gas_cells.temperature_code[cell_index];
    record.sound_speed_code = state.gas_cells.sound_speed_code[cell_index];
    record.recon_gradient_x = state.gas_cells.recon_gradient_x[cell_index];
    record.recon_gradient_y = state.gas_cells.recon_gradient_y[cell_index];
    record.recon_gradient_z = state.gas_cells.recon_gradient_z[cell_index];
    records.emplace(record.particle_id, record);
  }
  return records;
}

[[nodiscard]] std::vector<std::uint64_t> gasParticleIdByOldCellIndex(const core::SimulationState& state) {
  const auto gas_globals = state.particle_species_index.globalIndices(core::ParticleSpecies::kGas);
  if (state.cells.size() != gas_globals.size()) {
    throw std::runtime_error("gas cell count must match local gas particle count when capturing host-cell mapping");
  }
  std::vector<std::uint64_t> cell_particle_ids(state.cells.size(), 0ULL);
  for (std::size_t cell_index = 0; cell_index < gas_globals.size(); ++cell_index) {
    cell_particle_ids[cell_index] = state.particle_sidecar.particle_id[gas_globals[cell_index]];
  }
  return cell_particle_ids;
}

void rebuildLocalGasStateFromParticleIds(
    core::SimulationState& state,
    const std::unordered_map<std::uint64_t, GasCellMigrationRecord>& gas_records_by_particle_id,
    std::span<const std::uint64_t> old_cell_particle_id) {
  const auto gas_globals = state.particle_species_index.globalIndices(core::ParticleSpecies::kGas);
  core::CellSoa rebuilt_cells;
  core::GasCellSidecar rebuilt_gas;
  rebuilt_cells.resize(gas_globals.size());
  rebuilt_gas.resize(gas_globals.size());

  std::unordered_map<std::uint64_t, std::uint32_t> new_cell_index_by_particle_id;
  new_cell_index_by_particle_id.reserve(gas_globals.size());
  for (std::size_t cell_index = 0; cell_index < gas_globals.size(); ++cell_index) {
    const std::uint64_t particle_id = state.particle_sidecar.particle_id[gas_globals[cell_index]];
    const auto found = gas_records_by_particle_id.find(particle_id);
    if (found == gas_records_by_particle_id.end()) {
      throw std::runtime_error("missing gas-cell migration record for local gas particle after ownership compaction");
    }
    const GasCellMigrationRecord& record = found->second;
    rebuilt_cells.center_x_comoving[cell_index] = record.center_x_comoving;
    rebuilt_cells.center_y_comoving[cell_index] = record.center_y_comoving;
    rebuilt_cells.center_z_comoving[cell_index] = record.center_z_comoving;
    rebuilt_cells.mass_code[cell_index] = record.mass_code;
    rebuilt_cells.time_bin[cell_index] = record.time_bin;
    rebuilt_cells.patch_index[cell_index] = record.patch_index;
    rebuilt_gas.density_code[cell_index] = record.density_code;
    rebuilt_gas.pressure_code[cell_index] = record.pressure_code;
    rebuilt_gas.internal_energy_code[cell_index] = record.internal_energy_code;
    rebuilt_gas.temperature_code[cell_index] = record.temperature_code;
    rebuilt_gas.sound_speed_code[cell_index] = record.sound_speed_code;
    rebuilt_gas.recon_gradient_x[cell_index] = record.recon_gradient_x;
    rebuilt_gas.recon_gradient_y[cell_index] = record.recon_gradient_y;
    rebuilt_gas.recon_gradient_z[cell_index] = record.recon_gradient_z;
    new_cell_index_by_particle_id.emplace(particle_id, static_cast<std::uint32_t>(cell_index));
  }

  state.cells = std::move(rebuilt_cells);
  state.gas_cells = std::move(rebuilt_gas);

  const auto remap_host_cell = [&](std::uint32_t old_cell_index) {
    if (old_cell_index >= old_cell_particle_id.size()) {
      return std::uint32_t{0};
    }
    const std::uint64_t host_particle_id = old_cell_particle_id[old_cell_index];
    const auto found = new_cell_index_by_particle_id.find(host_particle_id);
    if (found == new_cell_index_by_particle_id.end()) {
      return std::uint32_t{0};
    }
    return found->second;
  };

  for (std::size_t row = 0; row < state.black_holes.size(); ++row) {
    state.black_holes.host_cell_index[row] = remap_host_cell(state.black_holes.host_cell_index[row]);
  }
  for (std::size_t row = 0; row < state.tracers.size(); ++row) {
    state.tracers.host_cell_index[row] = remap_host_cell(state.tracers.host_cell_index[row]);
  }
}

void compactStateToCurrentOwner(
    core::SimulationState& state,
    int world_rank) {
  if (world_rank < 0) {
    throw std::invalid_argument("compactStateToCurrentOwner requires non-negative world rank");
  }
  std::vector<std::uint32_t> kept_gas_particle_indices;
  kept_gas_particle_indices.reserve(state.cells.size());
  std::vector<std::uint32_t> outbound_indices;
  outbound_indices.reserve(state.particles.size());
  for (std::size_t particle_index = 0; particle_index < state.particles.size(); ++particle_index) {
    const bool owned_here = state.particle_sidecar.owning_rank[particle_index] == static_cast<std::uint32_t>(world_rank);
    if (owned_here) {
      if (state.particle_sidecar.species_tag[particle_index] == static_cast<std::uint32_t>(core::ParticleSpecies::kGas)) {
        kept_gas_particle_indices.push_back(static_cast<std::uint32_t>(particle_index));
      }
    } else {
      outbound_indices.push_back(static_cast<std::uint32_t>(particle_index));
    }
  }
  const auto gas_records = collectLocalGasCellRecords(state, kept_gas_particle_indices);
  const auto old_cell_particle_id = gasParticleIdByOldCellIndex(state);
  core::ParticleMigrationCommit commit;
  commit.world_rank = world_rank;
  commit.outbound_local_indices = std::move(outbound_indices);
  state.commitParticleMigration(commit);
  rebuildLocalGasStateFromParticleIds(state, gas_records, old_cell_particle_id);
}

void applyInitialGravityAwareDecomposition(
    core::SimulationState& state,
    const core::SimulationConfig& config,
    int world_size,
    int world_rank) {
  if (world_size <= 1) {
    return;
  }
  constexpr std::size_t k_density_grid = 16;
  const std::size_t grid_cells = k_density_grid * k_density_grid * k_density_grid;
  std::vector<std::uint32_t> occupancy(grid_cells, 0U);
  const auto wrap = [](double x, double box) {
    if (box <= 0.0) {
      return x;
    }
    double wrapped = std::fmod(x, box);
    if (wrapped < 0.0) {
      wrapped += box;
    }
    return wrapped;
  };
  const auto density_cell_index = [&](double x, double y, double z) {
    const double box_x = config.cosmology.box_size_x_mpc_comoving;
    const double box_y = config.cosmology.box_size_y_mpc_comoving;
    const double box_z = config.cosmology.box_size_z_mpc_comoving;
    const std::size_t ix = (box_x > 0.0)
        ? std::min<std::size_t>(
              k_density_grid - 1U,
              static_cast<std::size_t>((wrap(x, box_x) / box_x) * static_cast<double>(k_density_grid)))
        : 0U;
    const std::size_t iy = (box_y > 0.0)
        ? std::min<std::size_t>(
              k_density_grid - 1U,
              static_cast<std::size_t>((wrap(y, box_y) / box_y) * static_cast<double>(k_density_grid)))
        : 0U;
    const std::size_t iz = (box_z > 0.0)
        ? std::min<std::size_t>(
              k_density_grid - 1U,
              static_cast<std::size_t>((wrap(z, box_z) / box_z) * static_cast<double>(k_density_grid)))
        : 0U;
    return (ix * k_density_grid + iy) * k_density_grid + iz;
  };
  for (std::size_t particle_index = 0; particle_index < state.particles.size(); ++particle_index) {
    const std::size_t cell = density_cell_index(
        state.particles.position_x_comoving[particle_index],
        state.particles.position_y_comoving[particle_index],
        state.particles.position_z_comoving[particle_index]);
    ++occupancy[cell];
  }
  std::vector<parallel::DecompositionItem> items;
  items.reserve(state.particles.size());
  for (std::size_t particle_index = 0; particle_index < state.particles.size(); ++particle_index) {
    parallel::DecompositionItem item;
    item.entity_id = state.particle_sidecar.particle_id[particle_index];
    item.kind = parallel::DecompositionEntityKind::kParticle;
    item.x_comov = state.particles.position_x_comoving[particle_index];
    item.y_comov = state.particles.position_y_comoving[particle_index];
    item.z_comov = state.particles.position_z_comoving[particle_index];
    const std::size_t density_cell = density_cell_index(item.x_comov, item.y_comov, item.z_comov);
    const std::uint32_t local_density = occupancy[density_cell];
    item.active_target_count_recent = local_density;
    item.remote_tree_interactions_recent = 0;
    item.work_units = 1.0 + static_cast<double>(local_density);
    item.memory_bytes = sizeof(double) * 7 + sizeof(std::uint64_t) * 2 + sizeof(std::uint32_t) * 3;
    items.push_back(item);
  }
  parallel::DecompositionConfig decomposition_config;
  decomposition_config.world_size = world_size;
  decomposition_config.domain_x_min_comov = 0.0;
  decomposition_config.domain_x_max_comov = config.cosmology.box_size_x_mpc_comoving;
  decomposition_config.domain_y_min_comov = 0.0;
  decomposition_config.domain_y_max_comov = config.cosmology.box_size_y_mpc_comoving;
  decomposition_config.domain_z_min_comov = 0.0;
  decomposition_config.domain_z_max_comov = config.cosmology.box_size_z_mpc_comoving;
  decomposition_config.owned_particle_weight = 1.0;
  decomposition_config.active_target_weight = 2.0;
  decomposition_config.remote_tree_interaction_weight = 0.0;
  decomposition_config.work_weight = 0.5;
  decomposition_config.memory_weight = 0.1;
  const auto plan = parallel::buildMortonSfcDecomposition(items, decomposition_config);
  for (std::size_t item_index = 0; item_index < items.size(); ++item_index) {
    state.particle_sidecar.owning_rank[item_index] = static_cast<std::uint32_t>(plan.owning_rank_by_item[item_index]);
  }
  compactStateToCurrentOwner(state, world_rank);
}


[[nodiscard]] double wrapPeriodicPosition(double position_comoving, double box_size_mpc_comoving) {
  if (box_size_mpc_comoving <= 0.0) {
    return position_comoving;
  }
  double wrapped = std::fmod(position_comoving, box_size_mpc_comoving);
  if (wrapped < 0.0) {
    wrapped += box_size_mpc_comoving;
  }
  if (wrapped >= box_size_mpc_comoving) {
    wrapped -= box_size_mpc_comoving;
  }
  return wrapped;
}

void seedParticleOwnershipFromPmSlabs(
    core::SimulationState& state,
    std::size_t pm_grid_nx,
    int world_size,
    double box_size_x_mpc_comoving) {
  if (pm_grid_nx == 0) {
    throw std::invalid_argument("seedParticleOwnershipFromPmSlabs requires pm_grid_nx > 0");
  }
  if (world_size <= 0) {
    throw std::invalid_argument("seedParticleOwnershipFromPmSlabs requires world_size > 0");
  }
  for (std::size_t i = 0; i < state.particles.size(); ++i) {
    const double wrapped_x = wrapPeriodicPosition(state.particles.position_x_comoving[i], box_size_x_mpc_comoving);
    std::size_t global_x = 0;
    if (box_size_x_mpc_comoving > 0.0) {
      const double scaled = (wrapped_x / box_size_x_mpc_comoving) * static_cast<double>(pm_grid_nx);
      global_x = std::min(pm_grid_nx - 1U, static_cast<std::size_t>(scaled));
    }
    state.particle_sidecar.owning_rank[i] = static_cast<std::uint32_t>(
        parallel::pmOwnerRankForGlobalX(pm_grid_nx, world_size, global_x));
  }
}

[[nodiscard]] core::ProvenanceRecord makeGravityAwareProvenanceRecord(
    const core::FrozenConfig& frozen_config,
    const core::SimulationConfig& config) {
  core::ProvenanceRecord record = core::makeProvenanceRecord(
      frozen_config.provenance.config_hash_hex, frozen_config.provenance.source_name);
  const double dx = config.cosmology.box_size_x_mpc_comoving /
      static_cast<double>(config.numerics.treepm_pm_grid_nx);
  const double dy = config.cosmology.box_size_y_mpc_comoving /
      static_cast<double>(config.numerics.treepm_pm_grid_ny);
  const double dz = config.cosmology.box_size_z_mpc_comoving /
      static_cast<double>(config.numerics.treepm_pm_grid_nz);
  const double mesh_spacing_mpc_comoving = std::cbrt(dx * dy * dz);
  const gravity::TreePmSplitPolicy split_policy = gravity::makeTreePmSplitPolicyFromMeshSpacing(
      config.numerics.treepm_asmth_cells,
      config.numerics.treepm_rcut_cells,
      mesh_spacing_mpc_comoving);
  record.gravity_treepm_pm_grid = config.numerics.treepm_pm_grid_nx;
  record.gravity_treepm_assignment_scheme =
      treePmAssignmentSchemeName(config.numerics.treepm_assignment_scheme);
  record.gravity_treepm_window_deconvolution =
      config.numerics.treepm_enable_window_deconvolution;
  record.gravity_treepm_asmth_cells = config.numerics.treepm_asmth_cells;
  record.gravity_treepm_rcut_cells = config.numerics.treepm_rcut_cells;
  record.gravity_treepm_mesh_spacing_mpc_comoving = mesh_spacing_mpc_comoving;
  record.gravity_treepm_split_scale_mpc_comoving = split_policy.split_scale_comoving;
  record.gravity_treepm_cutoff_radius_mpc_comoving = split_policy.cutoff_radius_comoving;
  record.gravity_treepm_update_cadence_steps = config.numerics.treepm_update_cadence_steps;
  record.gravity_treepm_pm_decomposition_mode =
      pmDecompositionModeName(config.numerics.treepm_pm_decomposition_mode);
  record.gravity_treepm_tree_exchange_batch_bytes =
      config.numerics.treepm_tree_exchange_batch_bytes;
  record.gravity_softening_policy = "comoving_fixed";
  record.gravity_softening_kernel = "plummer";
  record.gravity_softening_epsilon_kpc_comoving = config.numerics.gravity_softening_kpc_comoving;
  record.gravity_pm_fft_backend = gravity::PmSolver::fftBackendName();
  switch (config.mode.zoom_long_range_strategy) {
    case core::ZoomLongRangeStrategy::kDisabled:
      record.zoom_long_range_strategy = "disabled";
      break;
    case core::ZoomLongRangeStrategy::kGlobalCoarsePlusFocusedHighResCorrection:
      record.zoom_long_range_strategy = "global_coarse_plus_focused_highres_correction";
      break;
  }
  record.zoom_region_center_x_mpc_comoving = config.mode.zoom_region_center_x_mpc_comoving;
  record.zoom_region_center_y_mpc_comoving = config.mode.zoom_region_center_y_mpc_comoving;
  record.zoom_region_center_z_mpc_comoving = config.mode.zoom_region_center_z_mpc_comoving;
  record.zoom_region_radius_mpc_comoving = config.mode.zoom_region_radius_mpc_comoving;
  record.zoom_focused_pm_grid = std::to_string(config.mode.zoom_focused_pm_grid_nx) + "x" +
      std::to_string(config.mode.zoom_focused_pm_grid_ny) + "x" +
      std::to_string(config.mode.zoom_focused_pm_grid_nz);
  record.zoom_contamination_radius_mpc_comoving = config.mode.zoom_contamination_radius_mpc_comoving;
  return record;
}

[[nodiscard]] std::string pmSlabSignature(const parallel::DistributedRestartState& distributed_state) {
  std::ostringstream stream;
  for (std::size_t rank = 0; rank < distributed_state.pm_slab_begin_x_by_rank.size(); ++rank) {
    if (rank > 0) {
      stream << ';';
    }
    stream << rank << ':' << distributed_state.pm_slab_begin_x_by_rank[rank] << '-' <<
        distributed_state.pm_slab_end_x_by_rank[rank];
  }
  return stream.str();
}

[[nodiscard]] std::filesystem::path computeRunDirectory(
    const core::SimulationConfig& config,
    const std::filesystem::path* output_root_override) {
  const std::filesystem::path output_root =
      (output_root_override != nullptr && !output_root_override->empty())
      ? *output_root_override
      : std::filesystem::path(config.output.output_directory);
  return output_root / config.output.run_name;
}

[[nodiscard]] std::filesystem::path resolveConfigRelativePath(
    const core::FrozenConfig& frozen_config,
    const std::filesystem::path& candidate) {
  if (candidate.is_absolute()) {
    return candidate;
  }

  const std::filesystem::path source_path(frozen_config.provenance.source_name);
  if (!source_path.empty() && source_path.has_parent_path()) {
    return source_path.parent_path() / candidate;
  }
  return candidate;
}

[[nodiscard]] bool repeatedCanonicalOrder(const std::vector<std::string>& observed) {
  const auto canonical = core::StageScheduler::kickDriftKickOrder();
  if (observed.empty() || observed.size() % canonical.size() != 0) {
    return false;
  }

  for (std::size_t chunk = 0; chunk < observed.size(); chunk += canonical.size()) {
    for (std::size_t i = 0; i < canonical.size(); ++i) {
      if (observed[chunk + i] != core::integrationStageName(canonical[i])) {
        return false;
      }
    }
  }
  return true;
}

[[nodiscard]] std::string formatIndexedFileStem(std::string_view stem, std::uint64_t index) {
  std::ostringstream out;
  out << stem << '_' << std::setw(3) << std::setfill('0') << index << ".hdf5";
  return out.str();
}

void fnv1aMix(std::uint64_t& hash, const void* data, std::size_t size) {
  constexpr std::uint64_t k_fnv_offset_basis = 1469598103934665603ULL;
  constexpr std::uint64_t k_fnv_prime = 1099511628211ULL;
  if (hash == 0) {
    hash = k_fnv_offset_basis;
  }
  const auto* bytes = static_cast<const unsigned char*>(data);
  for (std::size_t i = 0; i < size; ++i) {
    hash ^= static_cast<std::uint64_t>(bytes[i]);
    hash *= k_fnv_prime;
  }
}

void fnv1aMix(std::uint64_t& hash, double value) { fnv1aMix(hash, &value, sizeof(value)); }

void fnv1aMix(std::uint64_t& hash, std::uint64_t value) { fnv1aMix(hash, &value, sizeof(value)); }

[[nodiscard]] std::uint64_t computeStateDigest(const core::SimulationState& state, const core::IntegratorState& integrator_state) {
  std::uint64_t hash = 0;
  fnv1aMix(hash, static_cast<std::uint64_t>(state.particles.size()));
  fnv1aMix(hash, static_cast<std::uint64_t>(state.cells.size()));
  for (const double value : state.particles.position_x_comoving) {
    fnv1aMix(hash, value);
  }
  for (const double value : state.particles.position_y_comoving) {
    fnv1aMix(hash, value);
  }
  for (const double value : state.particles.position_z_comoving) {
    fnv1aMix(hash, value);
  }
  for (const double value : state.particles.velocity_x_peculiar) {
    fnv1aMix(hash, value);
  }
  for (const double value : state.particles.velocity_y_peculiar) {
    fnv1aMix(hash, value);
  }
  for (const double value : state.particles.velocity_z_peculiar) {
    fnv1aMix(hash, value);
  }
  for (const double value : state.particles.mass_code) {
    fnv1aMix(hash, value);
  }
  for (const double value : state.gas_cells.density_code) {
    fnv1aMix(hash, value);
  }
  for (const double value : state.gas_cells.internal_energy_code) {
    fnv1aMix(hash, value);
  }
  fnv1aMix(hash, integrator_state.current_time_code);
  fnv1aMix(hash, integrator_state.current_scale_factor);
  fnv1aMix(hash, integrator_state.step_index);
  return hash;
}

[[nodiscard]] std::uint64_t computeParticleIdSum(const core::SimulationState& state) {
  std::uint64_t sum = 0;
  for (const std::uint64_t particle_id : state.particle_sidecar.particle_id) {
    sum += particle_id;
  }
  return sum;
}

[[nodiscard]] std::uint64_t computeParticleIdXor(const core::SimulationState& state) {
  std::uint64_t value = 0;
  for (const std::uint64_t particle_id : state.particle_sidecar.particle_id) {
    value ^= particle_id;
  }
  return value;
}

void ensureRunDirectory(const std::filesystem::path& run_directory) {
  std::filesystem::create_directories(run_directory);
}

void flushCommonArtifacts(
    const core::FrozenConfig& frozen_config,
    core::ProfilerSession& profiler,
    ReferenceWorkflowReport& report) {
  ensureRunDirectory(report.run_directory);

  report.normalized_config_snapshot_path = report.run_directory / "normalized_config.param.txt";
  if (!std::filesystem::exists(report.normalized_config_snapshot_path)) {
    core::writeNormalizedConfigSnapshot(frozen_config, report.run_directory);
  }
  report.normalized_config_snapshot_written = true;

  report.profiler_json_path = report.run_directory / "profile.json";
  report.profiler_csv_path = report.run_directory / "profile.csv";
  report.operational_report_json_path = report.run_directory / "operational_events.json";
  core::writeProfilerReportJson(profiler, report.profiler_json_path);
  core::writeProfilerReportCsv(profiler, report.profiler_csv_path);
  core::writeOperationalReportJson(
      profiler,
      report.operational_report_json_path,
      frozen_config.config.output.run_name,
      frozen_config.provenance.config_hash_hex);
}

[[nodiscard]] io::IcReadResult loadInitialConditions(const core::FrozenConfig& frozen_config) {
  const core::SimulationConfig& config = frozen_config.config;
  if (config.mode.ic_file == "generated") {
    return io::convertGeneratedIsolatedIcToState(config, k_default_generated_particle_axis);
  }

  const std::filesystem::path ic_path =
      resolveConfigRelativePath(frozen_config, std::filesystem::path(config.mode.ic_file));
  return io::readGadgetArepoHdf5Ic(ic_path, config);
}

void finalizeStateMetadata(const core::FrozenConfig& frozen_config, core::SimulationState& state) {
  state.metadata.run_name = frozen_config.config.output.run_name;
  state.metadata.normalized_config_hash = frozen_config.provenance.config_hash;
  state.metadata.normalized_config_hash_hex = frozen_config.provenance.config_hash_hex;
  state.metadata.snapshot_stem = frozen_config.config.output.output_stem;
  state.metadata.restart_stem = frozen_config.config.output.restart_stem;
  state.metadata.scale_factor = (state.metadata.scale_factor > 0.0)
      ? state.metadata.scale_factor
      : std::max(1.0, frozen_config.config.numerics.time_begin_code);
  state.rebuildSpeciesIndex();
}

void initializeSchedulerBins(
    const core::SimulationState& state,
    core::HierarchicalTimeBinScheduler& scheduler) {
  const std::uint32_t particle_count = static_cast<std::uint32_t>(state.particles.size());
  const std::uint32_t cell_count = static_cast<std::uint32_t>(state.cells.size());
  const std::uint32_t scheduler_count = std::max(particle_count, cell_count);
  if (scheduler_count == 0U) {
    throw std::runtime_error("reference workflow requires non-empty initial conditions");
  }

  const std::uint8_t default_bin = scheduler.maxBin() > 0 ? 1U : 0U;
  scheduler.reset(scheduler_count, default_bin, 0);

  const auto gas_globals = state.particle_species_index.globalIndices(core::ParticleSpecies::kGas);
  for (std::uint32_t cell_index = 0; cell_index < cell_count; ++cell_index) {
    scheduler.setElementBin(cell_index, 0U, scheduler.currentTick());
    if (cell_index < gas_globals.size()) {
      scheduler.setElementBin(gas_globals[cell_index], 0U, scheduler.currentTick());
    }
  }
}

void syncTimeBinsFromScheduler(
    const core::HierarchicalTimeBinScheduler& scheduler,
    core::SimulationState& state) {
  const auto persistent = scheduler.exportPersistentState();
  for (std::size_t i = 0; i < state.particles.size() && i < persistent.bin_index.size(); ++i) {
    state.particles.time_bin[i] = persistent.bin_index[i];
  }
  for (std::size_t i = 0; i < state.cells.size() && i < persistent.bin_index.size(); ++i) {
    state.cells.time_bin[i] = persistent.bin_index[i];
  }
}

class StageAuditCallback final : public core::IntegrationCallback {
 public:
  explicit StageAuditCallback(std::vector<std::string>* stage_sequence)
      : m_stage_sequence(stage_sequence) {}

  [[nodiscard]] std::string_view callbackName() const override { return "stage_audit"; }

  void onStage(core::StepContext& context) override {
    m_stage_sequence->push_back(std::string(core::integrationStageName(context.stage)));
  }

 private:
  std::vector<std::string>* m_stage_sequence = nullptr;
};

class DriftCallback final : public core::IntegrationCallback {
 public:
  [[nodiscard]] std::string_view callbackName() const override { return "drift"; }

  void onStage(core::StepContext& context) override {
    if (context.stage != core::IntegrationStage::kDrift) {
      return;
    }

    const double dt = context.integrator_state.dt_time_code;
    const double inv_a = (context.integrator_state.current_scale_factor > 0.0)
        ? (1.0 / context.integrator_state.current_scale_factor)
        : 1.0;
    parallel::MpiContext mpi_context;
    const std::uint32_t world_rank = static_cast<std::uint32_t>(mpi_context.worldRank());
    for (const std::uint32_t particle_index : context.active_set.particle_indices) {
      if (particle_index >= context.state.particles.size()) {
        throw std::out_of_range("drift callback active particle index out of range");
      }
      if (particle_index >= context.state.particle_sidecar.owning_rank.size()) {
        throw std::out_of_range("drift callback owning-rank sidecar index out of range");
      }
      if (context.state.particle_sidecar.owning_rank[particle_index] != world_rank) {
        continue;
      }
      context.state.particles.position_x_comoving[particle_index] +=
          context.state.particles.velocity_x_peculiar[particle_index] * dt * inv_a;
      context.state.particles.position_y_comoving[particle_index] +=
          context.state.particles.velocity_y_peculiar[particle_index] * dt * inv_a;
      context.state.particles.position_z_comoving[particle_index] +=
          context.state.particles.velocity_z_peculiar[particle_index] * dt * inv_a;
    }

    const auto gas_globals = context.state.particle_species_index.globalIndices(core::ParticleSpecies::kGas);
    for (std::uint32_t cell_index = 0; cell_index < context.state.cells.size() && cell_index < gas_globals.size(); ++cell_index) {
      const std::uint32_t gas_particle_index = gas_globals[cell_index];
      context.state.cells.center_x_comoving[cell_index] = context.state.particles.position_x_comoving[gas_particle_index];
      context.state.cells.center_y_comoving[cell_index] = context.state.particles.position_y_comoving[gas_particle_index];
      context.state.cells.center_z_comoving[cell_index] = context.state.particles.position_z_comoving[gas_particle_index];
    }
  }
};

class GravityStageCallback final : public core::IntegrationCallback {
 public:
  struct GravityHealthSummary {
    std::uint64_t cheap_checks_executed = 0;
    std::uint64_t heavy_checks_executed = 0;
    std::uint64_t warning_count = 0;
    std::uint64_t fatal_count = 0;
    std::uint64_t pm_field_non_finite_count = 0;
    std::uint64_t force_non_finite_count = 0;
    std::uint64_t force_abnormal_count = 0;
    std::uint64_t illegal_sync_state_count = 0;
    std::uint64_t zoom_sanity_failure_count = 0;
    std::uint64_t decomposition_sanity_failure_count = 0;
  };

  enum class PmSyncSurface : std::uint8_t {
    kKickPre = 0,
    kKickPost = 1,
  };

  struct PmLongRangeKickDecision {
    PmSyncSurface sync_surface = PmSyncSurface::kKickPre;
    std::uint64_t gravity_kick_opportunity = 0;
    bool refresh_long_range_field = false;
    std::uint64_t field_version = 0;
    std::uint64_t last_refresh_opportunity = 0;
    std::uint64_t field_built_step_index = 0;
    double field_built_scale_factor = 1.0;
  };

  GravityStageCallback(const core::SimulationConfig& config, const core::ModePolicy& mode_policy)
      : m_config(config),
        m_mode_policy(mode_policy),
        m_pm_update_cadence_steps(static_cast<std::uint64_t>(config.numerics.treepm_update_cadence_steps)),
        m_pm_grid_shape(gravity::PmGridShape{
            static_cast<std::size_t>(config.numerics.treepm_pm_grid_nx),
            static_cast<std::size_t>(config.numerics.treepm_pm_grid_ny),
            static_cast<std::size_t>(config.numerics.treepm_pm_grid_nz)}),
        m_mesh_spacing_x_mpc_comoving(
            config.cosmology.box_size_x_mpc_comoving / static_cast<double>(m_pm_grid_shape.nx)),
        m_mesh_spacing_y_mpc_comoving(
            config.cosmology.box_size_y_mpc_comoving / static_cast<double>(m_pm_grid_shape.ny)),
        m_mesh_spacing_z_mpc_comoving(
            config.cosmology.box_size_z_mpc_comoving / static_cast<double>(m_pm_grid_shape.nz)),
        m_tree_pm_coordinator(makeRuntimeAwareTreePmCoordinator(config, m_pm_grid_shape)) {
    if (!m_pm_grid_shape.isValid()) {
      throw std::runtime_error("numerics.treepm_pm_grid_n{xyz} must all be > 0");
    }
    m_tree_pm_options.pm_options.box_size_x_mpc_comoving = config.cosmology.box_size_x_mpc_comoving;
    m_tree_pm_options.pm_options.box_size_y_mpc_comoving = config.cosmology.box_size_y_mpc_comoving;
    m_tree_pm_options.pm_options.box_size_z_mpc_comoving = config.cosmology.box_size_z_mpc_comoving;
    m_tree_pm_options.pm_options.box_size_mpc_comoving = config.cosmology.box_size_x_mpc_comoving;
    m_tree_pm_options.pm_options.scale_factor = 1.0;
    m_tree_pm_options.pm_options.gravitational_constant_code = 1.0;
    m_tree_pm_options.pm_options.assignment_scheme =
        toPmAssignmentScheme(config.numerics.treepm_assignment_scheme);
    m_tree_pm_options.pm_options.enable_window_deconvolution =
        config.numerics.treepm_enable_window_deconvolution;
    m_tree_pm_options.pm_options.decomposition_mode = config.numerics.treepm_pm_decomposition_mode;
    m_tree_pm_options.pm_options.boundary_condition = mode_policy.gravity_boundary ==
            core::GravityBoundaryModel::kPeriodicPoisson
        ? gravity::PmBoundaryCondition::kPeriodic
        : gravity::PmBoundaryCondition::kIsolatedOpen;
    if (m_tree_pm_options.pm_options.boundary_condition == gravity::PmBoundaryCondition::kIsolatedOpen &&
        config.parallel.mpi_ranks_expected != 1) {
      throw std::runtime_error(
          "isolated PM gravity currently requires parallel.mpi_ranks_expected=1");
    }

    parallel::MpiContext mpi_context;
    const core::CudaRuntimeInfo cuda_runtime = core::queryCudaRuntime();
    m_runtime_topology = parallel::buildDistributedExecutionTopology(
        m_pm_grid_shape.nx,
        m_pm_grid_shape.ny,
        m_pm_grid_shape.nz,
        mpi_context,
        config.parallel.mpi_ranks_expected,
        config.parallel.gpu_devices,
        cuda_runtime.runtime_available,
        cuda_runtime.visible_device_count,
        pmDecompositionModeName(config.numerics.treepm_pm_decomposition_mode));
    if (m_runtime_topology.usesCuda()) {
      core::setCudaDeviceOrThrow(m_runtime_topology.device_assignment.assigned_device_index);
      m_tree_pm_options.pm_options.execution_policy = core::ExecutionPolicy::kCuda;
      m_tree_pm_options.pm_options.data_residency = gravity::PmDataResidencyPolicy::kPreferDevice;
    }

    m_tree_pm_options.tree_options.opening_theta = 0.7;
    m_tree_pm_options.tree_options.gravitational_constant_code = 1.0;
    m_tree_pm_options.tree_options.softening.kernel = gravity::TreeSofteningKernel::kPlummer;
    m_tree_pm_options.tree_options.softening.epsilon_comoving =
        config.numerics.gravity_softening_kpc_comoving * 1.0e-3;
    m_tree_pm_options.split_policy = gravity::makeTreePmSplitPolicyFromMeshSpacing(
        config.numerics.treepm_asmth_cells,
        config.numerics.treepm_rcut_cells,
        std::cbrt(
            m_mesh_spacing_x_mpc_comoving * m_mesh_spacing_y_mpc_comoving *
            m_mesh_spacing_z_mpc_comoving));
    m_tree_pm_options.enable_zoom_long_range_correction =
        config.mode.zoom_long_range_strategy ==
        core::ZoomLongRangeStrategy::kGlobalCoarsePlusFocusedHighResCorrection;
    m_tree_pm_options.zoom_focused_pm_shape = gravity::PmGridShape{
        static_cast<std::size_t>(std::max(config.mode.zoom_focused_pm_grid_nx, 0)),
        static_cast<std::size_t>(std::max(config.mode.zoom_focused_pm_grid_ny, 0)),
        static_cast<std::size_t>(std::max(config.mode.zoom_focused_pm_grid_nz, 0))};
    m_tree_pm_options.zoom_region_center_x_comoving = config.mode.zoom_region_center_x_mpc_comoving;
    m_tree_pm_options.zoom_region_center_y_comoving = config.mode.zoom_region_center_y_mpc_comoving;
    m_tree_pm_options.zoom_region_center_z_comoving = config.mode.zoom_region_center_z_mpc_comoving;
    m_tree_pm_options.zoom_region_radius_comoving = config.mode.zoom_region_radius_mpc_comoving;
    m_tree_pm_options.zoom_contamination_radius_comoving =
        config.mode.zoom_contamination_radius_mpc_comoving > 0.0
        ? config.mode.zoom_contamination_radius_mpc_comoving
        : config.mode.zoom_region_radius_mpc_comoving;
    m_tree_pm_options.tree_exchange_batch_bytes = config.numerics.treepm_tree_exchange_batch_bytes;
    m_pm_assignment_scheme = treePmAssignmentSchemeName(config.numerics.treepm_assignment_scheme);
    m_pm_backend = gravity::PmSolver::fftBackendName();
    if (m_runtime_topology.usesCuda()) {
      m_pm_backend += "+cuda_cic";
    }
  }

  [[nodiscard]] std::string_view callbackName() const override { return "gravity"; }
  [[nodiscard]] std::size_t pmGridSize() const noexcept { return m_pm_grid_shape.nx; }
  [[nodiscard]] const gravity::PmGridShape& pmGridShape() const noexcept { return m_pm_grid_shape; }
  [[nodiscard]] std::uint64_t longRangeRefreshCount() const noexcept { return m_long_range_refresh_count; }
  [[nodiscard]] std::uint64_t longRangeReuseCount() const noexcept { return m_long_range_reuse_count; }
  [[nodiscard]] int pmCadenceSteps() const noexcept { return static_cast<int>(m_pm_update_cadence_steps); }
  [[nodiscard]] std::span<const ReferenceWorkflowReport::TreePmCadenceRecord> cadenceRecords() const noexcept {
    return m_cadence_records;
  }
  [[nodiscard]] std::uint64_t gravityKickOpportunity() const noexcept { return m_gravity_kick_opportunity; }
  [[nodiscard]] std::uint64_t longRangeFieldVersion() const noexcept { return m_long_range_field_version; }
  [[nodiscard]] std::uint64_t lastLongRangeRefreshOpportunity() const noexcept { return m_last_long_range_refresh_opportunity; }
  [[nodiscard]] std::uint64_t lastLongRangeRefreshStepIndex() const noexcept { return m_last_long_range_refresh_step_index; }
  [[nodiscard]] double lastLongRangeRefreshScaleFactor() const noexcept { return m_last_long_range_refresh_scale_factor; }
  [[nodiscard]] const parallel::DistributedExecutionTopology& runtimeTopology() const noexcept { return m_runtime_topology; }

  [[nodiscard]] std::span<const double> cellAccelX() const noexcept { return m_cell_accel_x; }
  [[nodiscard]] std::span<const double> cellAccelY() const noexcept { return m_cell_accel_y; }
  [[nodiscard]] std::span<const double> cellAccelZ() const noexcept { return m_cell_accel_z; }

  [[nodiscard]] static PmSyncSurface toPmSyncSurface(core::IntegrationStage stage) {
    if (stage == core::IntegrationStage::kGravityKickPre) {
      return PmSyncSurface::kKickPre;
    }
    if (stage == core::IntegrationStage::kGravityKickPost) {
      return PmSyncSurface::kKickPost;
    }
    throw std::invalid_argument("TreePM sync surface is only defined on gravity kick stages");
  }

  [[nodiscard]] static std::string_view pmSyncSurfaceName(PmSyncSurface surface) {
    switch (surface) {
      case PmSyncSurface::kKickPre:
        return "kick_pre";
      case PmSyncSurface::kKickPost:
        return "kick_post";
    }
    return "unknown";
  }

  void onStage(core::StepContext& context) override {
    if (context.stage != core::IntegrationStage::kGravityKickPre &&
        context.stage != core::IntegrationStage::kGravityKickPost) {
      return;
    }

    parallel::MpiContext mpi_context;
    const std::uint64_t world_size = static_cast<std::uint64_t>(mpi_context.worldSize());
    if (m_runtime_topology.world_size != mpi_context.worldSize()) {
      throw std::runtime_error("reference workflow runtime topology world_size drifted from MPI context");
    }

    const auto requireKickConsensus = [&](std::uint64_t local_value, std::string_view name) {
      const std::uint64_t global_sum = mpi_context.allreduceSumUint64(local_value);
      if (global_sum != local_value * world_size) {
        throw std::runtime_error(
            "TreePM cadence rank-consensus failure for " + std::string(name) +
            ": local=" + std::to_string(local_value) + ", reduced_sum=" + std::to_string(global_sum) +
            ", world_size=" + std::to_string(world_size));
      }
    };

    const std::size_t particle_count = context.state.particles.size();

    // Contract: kick-opportunity cadence metadata is rank-global and must advance
    // identically on every gravity kick stage, even if a rank has no active targets.
    ++m_gravity_kick_opportunity;
    requireKickConsensus(m_gravity_kick_opportunity, "gravity_kick_opportunity");

    const bool refresh_long_range_local = !m_has_long_range_field ||
        ((m_gravity_kick_opportunity - m_last_long_range_refresh_opportunity) >= m_pm_update_cadence_steps);
    const std::uint64_t refresh_vote = refresh_long_range_local ? 1ULL : 0ULL;
    const std::uint64_t refresh_vote_sum = mpi_context.allreduceSumUint64(refresh_vote);
    if (refresh_vote_sum != 0ULL && refresh_vote_sum != world_size) {
      throw std::runtime_error(
          "TreePM long-range cadence decision diverged across ranks: refresh_vote_sum=" +
          std::to_string(refresh_vote_sum) + ", world_size=" + std::to_string(world_size));
    }
    const bool refresh_long_range = (refresh_vote_sum == world_size);
    if (refresh_long_range != refresh_long_range_local) {
      throw std::runtime_error(
          "TreePM long-range cadence local decision drifted from rank consensus");
    }

    rebuildOwnedParticleCompactView(context.state, context.active_set.particle_indices);
    m_active_accel_x.assign(m_local_active_indices.size(), 0.0);
    m_active_accel_y.assign(m_local_active_indices.size(), 0.0);
    m_active_accel_z.assign(m_local_active_indices.size(), 0.0);
    m_active_is_high_res.assign(m_local_active_indices.size(), 0U);
    m_active_slot_by_particle.assign(particle_count, -1);
    for (std::size_t i = 0; i < m_local_active_global_indices.size(); ++i) {
      m_active_slot_by_particle[m_local_active_global_indices[i]] = static_cast<int>(i);
    }
    const double box_size = m_config.cosmology.box_size_x_mpc_comoving;
    m_source_is_high_res.assign(m_local_source_x.size(), 0U);
    for (std::size_t i = 0; i < m_local_source_x.size(); ++i) {
      const double dx = m_local_source_x[i] - m_tree_pm_options.zoom_region_center_x_comoving;
      const double dy = m_local_source_y[i] - m_tree_pm_options.zoom_region_center_y_comoving;
      const double dz = m_local_source_z[i] - m_tree_pm_options.zoom_region_center_z_comoving;
      const double wrapped_dx = dx - box_size * std::nearbyint(dx / box_size);
      const double wrapped_dy = dy - box_size * std::nearbyint(dy / box_size);
      const double wrapped_dz = dz - box_size * std::nearbyint(dz / box_size);
      const double r = std::sqrt(wrapped_dx * wrapped_dx + wrapped_dy * wrapped_dy + wrapped_dz * wrapped_dz);
      m_source_is_high_res[i] = (r <= m_tree_pm_options.zoom_region_radius_comoving) ? 1U : 0U;
    }
    for (std::size_t i = 0; i < m_local_active_indices.size(); ++i) {
      m_active_is_high_res[i] = m_source_is_high_res[m_local_active_indices[i]];
    }
    m_tree_pm_options.source_is_high_res = m_source_is_high_res;
    m_tree_pm_options.active_is_high_res = m_active_is_high_res;

    gravity::TreePmForceAccumulatorView accumulator{
        .active_particle_index = m_local_active_indices,
        .accel_x_comoving = m_active_accel_x,
        .accel_y_comoving = m_active_accel_y,
        .accel_z_comoving = m_active_accel_z,
    };

    m_tree_pm_options.pm_options.scale_factor = std::max(1.0, context.integrator_state.current_scale_factor);
    if (!m_mode_policy.cosmological_comoving_frame) {
      m_tree_pm_options.pm_options.scale_factor = 1.0;
    }

    PmLongRangeKickDecision decision{
        .sync_surface = toPmSyncSurface(context.stage),
        .gravity_kick_opportunity = m_gravity_kick_opportunity,
        .refresh_long_range_field = refresh_long_range,
        .field_version = refresh_long_range ? (m_long_range_field_version + 1U) : m_long_range_field_version,
        .last_refresh_opportunity = refresh_long_range ? m_gravity_kick_opportunity : m_last_long_range_refresh_opportunity,
        .field_built_step_index = refresh_long_range
            ? context.integrator_state.step_index
            : m_last_long_range_refresh_step_index,
        .field_built_scale_factor = refresh_long_range
            ? m_tree_pm_options.pm_options.scale_factor
            : m_last_long_range_refresh_scale_factor,
    };
    requireKickConsensus(decision.field_version, "long_range_field_version");
    requireKickConsensus(decision.last_refresh_opportunity, "last_long_range_refresh_opportunity");

    m_tree_pm_coordinator.solveActiveSetWithPmCadence(
        m_local_source_x,
        m_local_source_y,
        m_local_source_z,
        m_local_source_mass,
        accumulator,
        m_tree_pm_options,
        decision.refresh_long_range_field,
        nullptr,
        &m_last_tree_pm_diagnostics,
        {});

    const bool allow_heavy_reference_checks =
        m_config.analysis.diagnostics_execution_policy ==
        core::AnalysisConfig::DiagnosticsExecutionPolicy::kAllIncludingProvisional;
    GravityHealthSummary gravity_health = validateGravityHealth(
        context,
        decision,
        allow_heavy_reference_checks);

    if (decision.refresh_long_range_field) {
      m_has_long_range_field = true;
      m_long_range_field_version = decision.field_version;
      m_last_long_range_refresh_opportunity = decision.last_refresh_opportunity;
      m_last_long_range_refresh_step_index = decision.field_built_step_index;
      m_last_long_range_refresh_scale_factor = decision.field_built_scale_factor;
      ++m_long_range_refresh_count;
    } else {
      ++m_long_range_reuse_count;
    }
    m_last_committed_field_version = m_long_range_field_version;

    const std::uint64_t inactive_particles_skipped = static_cast<std::uint64_t>(
        context.state.particles.size() - m_local_active_global_indices.size());
    m_cadence_records.push_back(ReferenceWorkflowReport::TreePmCadenceRecord{
        .step_index = context.integrator_state.step_index,
        .stage_name = std::string(core::integrationStageName(context.stage)),
        .pm_sync_surface = std::string(pmSyncSurfaceName(decision.sync_surface)),
        .gravity_kick_opportunity = decision.gravity_kick_opportunity,
        .field_version = decision.field_version,
        .last_refresh_opportunity = decision.last_refresh_opportunity,
        .field_built_step_index = decision.field_built_step_index,
        .field_built_scale_factor = decision.field_built_scale_factor,
        .field_age_in_kick_opportunities =
            decision.gravity_kick_opportunity - decision.last_refresh_opportunity,
        .active_particles_kicked = static_cast<std::uint64_t>(m_local_active_global_indices.size()),
        .inactive_particles_skipped = inactive_particles_skipped,
        .refreshed_long_range_field = decision.refresh_long_range_field,
    });
    if (context.profiler_session != nullptr) {
      context.profiler_session->recordEvent(core::RuntimeEvent{
          .event_kind = "gravity.pm_long_range_field",
          .severity = core::RuntimeEventSeverity::kInfo,
          .subsystem = "gravity.treepm",
          .step_index = context.integrator_state.step_index,
          .simulation_time_code = context.integrator_state.current_time_code,
          .scale_factor = context.integrator_state.current_scale_factor,
          .message = decision.refresh_long_range_field
              ? "PM long-range field refreshed for gravity kick"
              : "PM long-range field reused for gravity kick",
          .payload = {{"stage", std::string(core::integrationStageName(context.stage))},
                      {"pm_sync_surface", std::string(pmSyncSurfaceName(decision.sync_surface))},
                      {"gravity_kick_opportunity", std::to_string(m_gravity_kick_opportunity)},
                      {"field_version", std::to_string(m_long_range_field_version)},
                      {"field_built_step_index", std::to_string(m_last_long_range_refresh_step_index)},
                      {"field_built_scale_factor", std::to_string(m_last_long_range_refresh_scale_factor)},
                      {"pm_update_cadence_steps", std::to_string(m_pm_update_cadence_steps)},
                      {"pm_grid", std::to_string(m_pm_grid_shape.nx) + "x" + std::to_string(m_pm_grid_shape.ny) +
                              "x" + std::to_string(m_pm_grid_shape.nz)},
                      {"pm_assignment_scheme", m_pm_assignment_scheme},
                      {"pm_window_deconvolution", m_config.numerics.treepm_enable_window_deconvolution ? "true" : "false"},
                      {"asmth_cells", std::to_string(m_config.numerics.treepm_asmth_cells)},
                      {"rcut_cells", std::to_string(m_config.numerics.treepm_rcut_cells)},
                      {"mesh_spacing_x_mpc_comoving", std::to_string(m_mesh_spacing_x_mpc_comoving)},
                      {"mesh_spacing_y_mpc_comoving", std::to_string(m_mesh_spacing_y_mpc_comoving)},
                      {"mesh_spacing_z_mpc_comoving", std::to_string(m_mesh_spacing_z_mpc_comoving)},
                      {"split_scale_mpc_comoving", std::to_string(m_tree_pm_options.split_policy.split_scale_comoving)},
                      {"cutoff_radius_mpc_comoving", std::to_string(m_tree_pm_options.split_policy.cutoff_radius_comoving)},
                      {"softening_policy", "comoving_fixed"},
                      {"softening_kernel", "plummer"},
                      {"softening_epsilon_kpc_comoving", std::to_string(m_config.numerics.gravity_softening_kpc_comoving)},
                      {"pm_fft_backend", m_pm_backend},
                      {"active_particles_kicked", std::to_string(m_local_active_global_indices.size())},
                      {"inactive_particles_skipped", std::to_string(inactive_particles_skipped)},
                      {"refreshed_long_range_field", decision.refresh_long_range_field ? "true" : "false"}},
      });
      context.profiler_session->recordEvent(core::RuntimeEvent{
          .event_kind = "gravity.zoom_force_diagnostics",
          .severity = core::RuntimeEventSeverity::kInfo,
          .subsystem = "gravity.treepm",
          .step_index = context.integrator_state.step_index,
          .simulation_time_code = context.integrator_state.current_time_code,
          .scale_factor = context.integrator_state.current_scale_factor,
          .message = "zoom force decomposition and contamination diagnostics",
          .payload = {
              {"force_l2_pm_global", std::to_string(m_last_tree_pm_diagnostics.force_l2_pm_global)},
              {"force_l2_pm_zoom_correction", std::to_string(m_last_tree_pm_diagnostics.force_l2_pm_zoom_correction)},
              {"force_l2_tree_short_range", std::to_string(m_last_tree_pm_diagnostics.force_l2_tree_short_range)},
              {"force_l2_total", std::to_string(m_last_tree_pm_diagnostics.force_l2_total)},
              {"zoom_high_res_source_count", std::to_string(m_last_tree_pm_diagnostics.zoom_high_res_source_count)},
              {"zoom_low_res_source_count", std::to_string(m_last_tree_pm_diagnostics.zoom_low_res_source_count)},
              {"zoom_low_res_contamination_count", std::to_string(m_last_tree_pm_diagnostics.zoom_low_res_contamination_count)},
              {"zoom_low_res_contamination_mass_code", std::to_string(m_last_tree_pm_diagnostics.zoom_low_res_contamination_mass_code)},
          },
      });
      context.profiler_session->recordEvent(core::RuntimeEvent{
          .event_kind = "gravity.health_summary",
          .severity = (gravity_health.fatal_count > 0U)
              ? core::RuntimeEventSeverity::kFatal
              : (gravity_health.warning_count > 0U)
                    ? core::RuntimeEventSeverity::kWarning
                    : core::RuntimeEventSeverity::kInfo,
          .subsystem = "gravity.treepm",
          .step_index = context.integrator_state.step_index,
          .simulation_time_code = context.integrator_state.current_time_code,
          .scale_factor = context.integrator_state.current_scale_factor,
          .message = "gravity health summary across PM field, force, sync, zoom, and decomposition checks",
          .payload = {
              {"cheap_checks_executed", std::to_string(gravity_health.cheap_checks_executed)},
              {"heavy_checks_executed", std::to_string(gravity_health.heavy_checks_executed)},
              {"warning_count", std::to_string(gravity_health.warning_count)},
              {"fatal_count", std::to_string(gravity_health.fatal_count)},
              {"pm_field_non_finite_count", std::to_string(gravity_health.pm_field_non_finite_count)},
              {"force_non_finite_count", std::to_string(gravity_health.force_non_finite_count)},
              {"force_abnormal_count", std::to_string(gravity_health.force_abnormal_count)},
              {"illegal_sync_state_count", std::to_string(gravity_health.illegal_sync_state_count)},
              {"zoom_sanity_failure_count", std::to_string(gravity_health.zoom_sanity_failure_count)},
              {"decomposition_sanity_failure_count", std::to_string(gravity_health.decomposition_sanity_failure_count)},
              {"heavy_reference_checks_opt_in", allow_heavy_reference_checks ? "true" : "false"},
          },
      });
    }

    const double kick_factor = 0.5 * context.integrator_state.dt_time_code;
    for (std::size_t active_slot = 0; active_slot < m_local_active_global_indices.size(); ++active_slot) {
      const std::uint32_t particle_index = m_local_active_global_indices[active_slot];
      context.state.particles.velocity_x_peculiar[particle_index] += kick_factor * m_active_accel_x[active_slot];
      context.state.particles.velocity_y_peculiar[particle_index] += kick_factor * m_active_accel_y[active_slot];
      context.state.particles.velocity_z_peculiar[particle_index] += kick_factor * m_active_accel_z[active_slot];
    }

    const auto gas_globals = context.state.particle_species_index.globalIndices(core::ParticleSpecies::kGas);
    m_cell_accel_x.assign(context.state.cells.size(), 0.0);
    m_cell_accel_y.assign(context.state.cells.size(), 0.0);
    m_cell_accel_z.assign(context.state.cells.size(), 0.0);
    for (std::size_t cell_index = 0; cell_index < context.state.cells.size() && cell_index < gas_globals.size(); ++cell_index) {
      const int active_slot = m_active_slot_by_particle[gas_globals[cell_index]];
      if (active_slot < 0) {
        continue;
      }
      m_cell_accel_x[cell_index] = m_active_accel_x[static_cast<std::size_t>(active_slot)];
      m_cell_accel_y[cell_index] = m_active_accel_y[static_cast<std::size_t>(active_slot)];
      m_cell_accel_z[cell_index] = m_active_accel_z[static_cast<std::size_t>(active_slot)];
    }
  }

 private:
  void emitGravityEvent(
      core::StepContext& context,
      core::RuntimeEventSeverity severity,
      std::string message,
      std::unordered_map<std::string, std::string> payload = {}) const {
    if (context.profiler_session == nullptr) {
      return;
    }
    context.profiler_session->recordEvent(core::RuntimeEvent{
        .event_kind = "gravity.health_check",
        .severity = severity,
        .subsystem = "gravity.treepm",
        .step_index = context.integrator_state.step_index,
        .simulation_time_code = context.integrator_state.current_time_code,
        .scale_factor = context.integrator_state.current_scale_factor,
        .message = std::move(message),
        .payload = std::move(payload),
    });
  }

  [[nodiscard]] GravityHealthSummary validateGravityHealth(
      core::StepContext& context,
      const PmLongRangeKickDecision& decision,
      bool allow_heavy_reference_checks) {
    GravityHealthSummary summary;
    auto failFatal = [&](std::string message, std::unordered_map<std::string, std::string> payload = {}) {
      ++summary.fatal_count;
      emitGravityEvent(context, core::RuntimeEventSeverity::kFatal, message, payload);
      throw std::runtime_error("fatal gravity-state check failed: " + message);
    };
    auto addWarning = [&](std::string message, std::unordered_map<std::string, std::string> payload = {}) {
      ++summary.warning_count;
      emitGravityEvent(context, core::RuntimeEventSeverity::kWarning, message, payload);
    };

    ++summary.cheap_checks_executed;
    if (m_local_active_indices.size() != m_local_active_global_indices.size()) {
      ++summary.decomposition_sanity_failure_count;
      failFatal("local active index vectors diverged after compact-view rebuild");
    }
    if (m_local_active_indices.size() > m_local_source_x.size()) {
      ++summary.decomposition_sanity_failure_count;
      failFatal("active local subset exceeds local source population");
    }

    ++summary.cheap_checks_executed;
    if (!std::isfinite(m_last_tree_pm_diagnostics.force_l2_pm_global) ||
        !std::isfinite(m_last_tree_pm_diagnostics.force_l2_pm_zoom_correction) ||
        !std::isfinite(m_last_tree_pm_diagnostics.force_l2_tree_short_range) ||
        !std::isfinite(m_last_tree_pm_diagnostics.force_l2_total)) {
      ++summary.pm_field_non_finite_count;
      failFatal("PM/tree force norm diagnostics contain NaN or Inf");
    }

    ++summary.cheap_checks_executed;
    for (std::size_t i = 0; i < m_active_accel_x.size(); ++i) {
      if (!std::isfinite(m_active_accel_x[i]) ||
          !std::isfinite(m_active_accel_y[i]) ||
          !std::isfinite(m_active_accel_z[i])) {
        ++summary.force_non_finite_count;
        failFatal(
            "active-set acceleration contains NaN or Inf",
            {{"active_slot", std::to_string(i)}});
      }
    }

    ++summary.cheap_checks_executed;
    if (m_last_long_range_refresh_opportunity > m_gravity_kick_opportunity ||
        m_long_range_field_version < m_last_committed_field_version) {
      ++summary.illegal_sync_state_count;
      failFatal("gravity sync-state regressed across kick opportunities");
    }
    if (decision.refresh_long_range_field) {
      if (decision.field_version != m_long_range_field_version + 1U) {
        ++summary.illegal_sync_state_count;
        failFatal("refresh decision must increment field_version by exactly one");
      }
    } else if (decision.field_version != m_long_range_field_version) {
      ++summary.illegal_sync_state_count;
      failFatal("reuse decision must not mutate field_version");
    }

    ++summary.cheap_checks_executed;
    if (m_mode_policy.zoom_region_expected &&
        m_tree_pm_options.zoom_region_radius_comoving <= 0.0) {
      ++summary.zoom_sanity_failure_count;
      failFatal("zoom mode requires a strictly positive zoom region radius");
    }
    if (!m_mode_policy.zoom_region_expected &&
        m_tree_pm_options.enable_zoom_long_range_correction &&
        m_last_tree_pm_diagnostics.zoom_high_res_source_count == 0U &&
        m_last_tree_pm_diagnostics.zoom_low_res_source_count == 0U) {
      ++summary.zoom_sanity_failure_count;
      addWarning("zoom long-range correction enabled without detected zoom population");
    }

    if (allow_heavy_reference_checks) {
      ++summary.heavy_checks_executed;
      const double total_force = m_last_tree_pm_diagnostics.force_l2_total;
      const double pm_component = m_last_tree_pm_diagnostics.force_l2_pm_global;
      const double tree_component = m_last_tree_pm_diagnostics.force_l2_tree_short_range;
      if (total_force > 0.0 && std::isfinite(total_force)) {
        const double component_sum = pm_component + tree_component;
        const double relative_gap = std::abs(component_sum - total_force) / total_force;
        if (relative_gap > 1.0) {
          ++summary.force_abnormal_count;
          addWarning(
              "force decomposition appears abnormal under heavy reference check",
              {{"relative_gap", std::to_string(relative_gap)}});
        }
      }
      ++summary.heavy_checks_executed;
      if (m_mode_policy.zoom_region_expected &&
          m_last_tree_pm_diagnostics.zoom_low_res_contamination_count > 0U) {
        ++summary.force_abnormal_count;
        addWarning(
            "zoom contamination detected in low-resolution source particles",
            {{"zoom_low_res_contamination_count",
              std::to_string(m_last_tree_pm_diagnostics.zoom_low_res_contamination_count)}});
      }
    }

    return summary;
  }

  void rebuildOwnedParticleCompactView(const core::SimulationState& state, std::span<const std::uint32_t> active_particles) {
    const std::uint32_t world_rank = static_cast<std::uint32_t>(m_runtime_topology.world_rank);
    const std::size_t particle_count = state.particles.size();
    m_owned_local_index_by_global.assign(particle_count, -1);
    m_local_source_x.clear();
    m_local_source_y.clear();
    m_local_source_z.clear();
    m_local_source_mass.clear();
    m_local_to_global.clear();
    for (std::size_t global_index = 0; global_index < particle_count; ++global_index) {
      if (state.particle_sidecar.owning_rank[global_index] != world_rank) {
        continue;
      }
      m_owned_local_index_by_global[global_index] = static_cast<int>(m_local_to_global.size());
      m_local_to_global.push_back(static_cast<std::uint32_t>(global_index));
      m_local_source_x.push_back(state.particles.position_x_comoving[global_index]);
      m_local_source_y.push_back(state.particles.position_y_comoving[global_index]);
      m_local_source_z.push_back(state.particles.position_z_comoving[global_index]);
      m_local_source_mass.push_back(state.particles.mass_code[global_index]);
    }

    m_local_active_indices.clear();
    m_local_active_global_indices.clear();
    for (const std::uint32_t global_index : active_particles) {
      if (global_index >= particle_count) {
        throw std::out_of_range("gravity callback active particle index out of range");
      }
      const int local_index = m_owned_local_index_by_global[global_index];
      if (local_index < 0) {
        continue;
      }
      m_local_active_indices.push_back(static_cast<std::uint32_t>(local_index));
      m_local_active_global_indices.push_back(global_index);
    }
  }

  const core::SimulationConfig& m_config;
  const core::ModePolicy& m_mode_policy;
  parallel::DistributedExecutionTopology m_runtime_topology{};
  std::uint64_t m_pm_update_cadence_steps = 1;
  gravity::PmGridShape m_pm_grid_shape{};
  double m_mesh_spacing_x_mpc_comoving = 0.0;
  double m_mesh_spacing_y_mpc_comoving = 0.0;
  double m_mesh_spacing_z_mpc_comoving = 0.0;
  std::string m_pm_assignment_scheme = "unknown";
  std::string m_pm_backend = "unknown";
  gravity::TreePmCoordinator m_tree_pm_coordinator;
  gravity::TreePmOptions m_tree_pm_options;
  std::vector<double> m_active_accel_x;
  std::vector<double> m_active_accel_y;
  std::vector<double> m_active_accel_z;
  std::vector<double> m_cell_accel_x;
  std::vector<double> m_cell_accel_y;
  std::vector<double> m_cell_accel_z;
  std::vector<int> m_active_slot_by_particle;
  std::vector<int> m_owned_local_index_by_global;
  std::vector<double> m_local_source_x;
  std::vector<double> m_local_source_y;
  std::vector<double> m_local_source_z;
  std::vector<double> m_local_source_mass;
  std::vector<std::uint8_t> m_source_is_high_res;
  std::vector<std::uint8_t> m_active_is_high_res;
  std::vector<std::uint32_t> m_local_to_global;
  std::vector<std::uint32_t> m_local_active_indices;
  std::vector<std::uint32_t> m_local_active_global_indices;
  gravity::TreePmDiagnostics m_last_tree_pm_diagnostics{};
  bool m_has_long_range_field = false;
  std::uint64_t m_gravity_kick_opportunity = 0;
  std::uint64_t m_last_long_range_refresh_opportunity = 0;
  std::uint64_t m_last_long_range_refresh_step_index = 0;
  double m_last_long_range_refresh_scale_factor = 1.0;
  std::uint64_t m_long_range_field_version = 0;
  std::uint64_t m_last_committed_field_version = 0;
  std::uint64_t m_long_range_refresh_count = 0;
  std::uint64_t m_long_range_reuse_count = 0;
  std::vector<ReferenceWorkflowReport::TreePmCadenceRecord> m_cadence_records;
};

class HydroStageCallback final : public core::IntegrationCallback {
 public:
  HydroStageCallback(const core::SimulationConfig& config, const core::ModePolicy& mode_policy, const GravityStageCallback& gravity_callback)
      : m_config(config),
        m_mode_policy(mode_policy),
        m_gravity_callback(gravity_callback),
        m_solver(k_gamma_adiabatic),
        m_reconstruction(hydro::HydroReconstructionPolicy{
            .limiter = hydro::HydroSlopeLimiter::kMonotonizedCentral,
            .dt_over_dx_code = 0.0,
            .rho_floor = k_density_floor,
            .pressure_floor = k_pressure_floor,
            .enable_muscl_hancock_predictor = true,
            .adiabatic_index = k_gamma_adiabatic,
        }) {}

  [[nodiscard]] std::string_view callbackName() const override { return "hydro"; }

  void onStage(core::StepContext& context) override {
    if (context.stage != core::IntegrationStage::kHydroUpdate || context.state.cells.size() == 0) {
      return;
    }

    const auto gas_globals = context.state.particle_species_index.globalIndices(core::ParticleSpecies::kGas);
    if (gas_globals.size() < context.state.cells.size()) {
      throw std::runtime_error("hydro callback requires one gas particle per gas cell in the current state model");
    }

    rebuildGeometryIfNeeded(context.state.cells.size(), context.integrator_state.dt_time_code);
    const hydro::HydroActiveSetView active_view = buildActiveFaceView(context.active_set.cell_indices);

    m_conserved.resize(context.state.cells.size());
    for (std::size_t cell_index = 0; cell_index < context.state.cells.size(); ++cell_index) {
      const std::uint32_t gas_particle_index = gas_globals[cell_index];
      const double rho = std::max(context.state.gas_cells.density_code[cell_index], k_density_floor);
      double pressure = context.state.gas_cells.pressure_code[cell_index];
      if (pressure <= 0.0) {
        const double internal_energy = std::max(context.state.gas_cells.internal_energy_code[cell_index], k_pressure_floor);
        pressure = std::max((k_gamma_adiabatic - 1.0) * rho * internal_energy, k_pressure_floor);
      }
      const hydro::HydroPrimitiveState primitive{
          .rho_comoving = rho,
          .vel_x_peculiar = context.state.particles.velocity_x_peculiar[gas_particle_index],
          .vel_y_peculiar = context.state.particles.velocity_y_peculiar[gas_particle_index],
          .vel_z_peculiar = context.state.particles.velocity_z_peculiar[gas_particle_index],
          .pressure_comoving = pressure,
      };
      m_conserved.storeCell(cell_index, hydro::HydroCoreSolver::conservedFromPrimitive(primitive, k_gamma_adiabatic));
    }

    double hubble_rate_code = 0.0;
    if (context.cosmology_background != nullptr && context.integrator_state.current_scale_factor > 0.0) {
      hubble_rate_code =
          core::computeScaleFactorRate(*context.cosmology_background, context.integrator_state.current_scale_factor) /
          context.integrator_state.current_scale_factor;
    }

    hydro::HydroUpdateContext update{
        .dt_code = context.integrator_state.dt_time_code,
        .scale_factor = std::max(1.0, context.integrator_state.current_scale_factor),
        .hubble_rate_code = hubble_rate_code,
    };

    std::vector<double> metallicity(context.state.cells.size(), 0.0);
    std::vector<double> temperature(context.state.cells.size(), 0.0);
    std::vector<double> hydrogen_number_density(context.state.cells.size(), 0.0);
    hydro::HydroSourceContext source_context{
        .update = update,
        .gravity_accel_x_peculiar = m_gravity_callback.cellAccelX(),
        .gravity_accel_y_peculiar = m_gravity_callback.cellAccelY(),
        .gravity_accel_z_peculiar = m_gravity_callback.cellAccelZ(),
        .hydrogen_number_density_cgs = hydrogen_number_density,
        .metallicity_mass_fraction = metallicity,
        .temperature_k = temperature,
        .redshift = std::max(0.0, (update.scale_factor > 0.0 ? (1.0 / update.scale_factor - 1.0) : 0.0)),
    };

    hydro::ComovingGravityExpansionSource gravity_source;
    std::array<const hydro::HydroSourceTerm*, 1> sources{&gravity_source};
    m_solver.advancePatchActiveSetWithScratch(
        m_conserved,
        m_geometry,
        active_view,
        update,
        m_reconstruction,
        m_riemann_solver,
        sources,
        source_context,
        m_scratch,
        &m_primitive_cache,
        nullptr);

    for (std::size_t cell_index = 0; cell_index < context.state.cells.size(); ++cell_index) {
      const std::uint32_t gas_particle_index = gas_globals[cell_index];
      const hydro::HydroPrimitiveState primitive =
          hydro::HydroCoreSolver::primitiveFromConserved(m_conserved.loadCell(cell_index), k_gamma_adiabatic);
      context.state.gas_cells.density_code[cell_index] = primitive.rho_comoving;
      context.state.gas_cells.pressure_code[cell_index] = primitive.pressure_comoving;
      context.state.gas_cells.internal_energy_code[cell_index] =
          primitive.pressure_comoving / ((k_gamma_adiabatic - 1.0) * std::max(primitive.rho_comoving, k_density_floor));
      context.state.cells.mass_code[cell_index] = primitive.rho_comoving * m_geometry.cell_volume_comoving;
      context.state.particles.mass_code[gas_particle_index] = context.state.cells.mass_code[cell_index];
      context.state.particles.velocity_x_peculiar[gas_particle_index] = primitive.vel_x_peculiar;
      context.state.particles.velocity_y_peculiar[gas_particle_index] = primitive.vel_y_peculiar;
      context.state.particles.velocity_z_peculiar[gas_particle_index] = primitive.vel_z_peculiar;
    }
  }

 private:
  void rebuildGeometryIfNeeded(std::size_t cell_count, double dt_time_code) {
    if (m_cached_cell_count == cell_count && cell_count > 0) {
      return;
    }

    m_cached_cell_count = cell_count;
    m_geometry = {};
    const double box_volume_comoving =
        m_config.cosmology.box_size_x_mpc_comoving *
        m_config.cosmology.box_size_y_mpc_comoving *
        m_config.cosmology.box_size_z_mpc_comoving;
    m_geometry.cell_volume_comoving = std::max(
        1.0,
        box_volume_comoving / static_cast<double>(std::max<std::size_t>(cell_count, 1)));
    m_geometry.faces.clear();
    if (cell_count == 0) {
      return;
    }

    for (std::size_t cell_index = 0; cell_index + 1 < cell_count; ++cell_index) {
      m_geometry.faces.push_back(hydro::HydroFace{
          .owner_cell = cell_index,
          .neighbor_cell = cell_index + 1,
          .area_comoving = 1.0,
          .normal_x = 1.0,
          .normal_y = 0.0,
          .normal_z = 0.0,
      });
    }

    if (m_mode_policy.hydro_boundary == core::BoundaryCondition::kPeriodic && cell_count > 1) {
      m_geometry.faces.push_back(hydro::HydroFace{
          .owner_cell = cell_count - 1,
          .neighbor_cell = 0,
          .area_comoving = 1.0,
          .normal_x = 1.0,
          .normal_y = 0.0,
          .normal_z = 0.0,
      });
    } else {
      m_geometry.faces.push_back(hydro::HydroFace{
          .owner_cell = 0,
          .neighbor_cell = hydro::k_invalid_cell_index,
          .area_comoving = 1.0,
          .normal_x = -1.0,
          .normal_y = 0.0,
          .normal_z = 0.0,
      });
      if (cell_count > 1) {
        m_geometry.faces.push_back(hydro::HydroFace{
            .owner_cell = cell_count - 1,
            .neighbor_cell = hydro::k_invalid_cell_index,
            .area_comoving = 1.0,
            .normal_x = 1.0,
            .normal_y = 0.0,
            .normal_z = 0.0,
        });
      }
    }

    const double dx = std::max(
        1.0e-6,
        m_config.cosmology.box_size_x_mpc_comoving / static_cast<double>(std::max<std::size_t>(cell_count, 1)));
    m_reconstruction = hydro::MusclHancockReconstruction(hydro::HydroReconstructionPolicy{
        .limiter = hydro::HydroSlopeLimiter::kMonotonizedCentral,
        .dt_over_dx_code = dt_time_code / dx,
        .rho_floor = k_density_floor,
        .pressure_floor = k_pressure_floor,
        .enable_muscl_hancock_predictor = true,
        .adiabatic_index = k_gamma_adiabatic,
    });
  }

  [[nodiscard]] hydro::HydroActiveSetView buildActiveFaceView(std::span<const std::uint32_t> active_cells) {
    m_active_cells.assign(active_cells.begin(), active_cells.end());
    std::unordered_set<std::size_t> active_cell_lookup(m_active_cells.begin(), m_active_cells.end());
    m_active_faces.clear();
    for (std::size_t face_index = 0; face_index < m_geometry.faces.size(); ++face_index) {
      const hydro::HydroFace& face = m_geometry.faces[face_index];
      if (active_cell_lookup.contains(face.owner_cell) ||
          (face.neighbor_cell != hydro::k_invalid_cell_index && active_cell_lookup.contains(face.neighbor_cell))) {
        m_active_faces.push_back(face_index);
      }
    }
    return hydro::HydroActiveSetView{.active_cells = m_active_cells, .active_faces = m_active_faces};
  }

  const core::SimulationConfig& m_config;
  const core::ModePolicy& m_mode_policy;
  const GravityStageCallback& m_gravity_callback;
  hydro::HydroCoreSolver m_solver;
  hydro::MusclHancockReconstruction m_reconstruction;
  hydro::HllcRiemannSolver m_riemann_solver;
  hydro::HydroConservedStateSoa m_conserved;
  hydro::HydroScratchBuffers m_scratch;
  hydro::HydroPrimitiveCacheSoa m_primitive_cache;
  hydro::HydroPatchGeometry m_geometry;
  std::vector<std::size_t> m_active_cells;
  std::vector<std::size_t> m_active_faces;
  std::size_t m_cached_cell_count = 0;
};

void maybeWriteOutputs(
    const core::FrozenConfig& frozen_config,
    const core::SimulationConfig& config,
    const core::SimulationState& state,
    const core::IntegratorState& integrator_state,
    const core::HierarchicalTimeBinScheduler& scheduler,
    const GravityStageCallback& gravity_callback,
    ReferenceWorkflowReport& report,
    core::ProfilerSession& profiler,
    bool write_outputs_enabled) {
  if (!write_outputs_enabled) {
    return;
  }

#if !COSMOSIM_ENABLE_HDF5
  throw std::runtime_error(
      "runtime outputs requested, but this build lacks HDF5 support. Reconfigure with COSMOSIM_ENABLE_HDF5=ON.");
#else
  if (integrator_state.step_index == 0) {
    return;
  }
  if ((integrator_state.step_index % static_cast<std::uint64_t>(config.output.snapshot_interval_steps)) != 0U) {
    return;
  }

  io::SnapshotWritePayload snapshot_payload;
  snapshot_payload.state = &state;
  snapshot_payload.config = &config;
  snapshot_payload.normalized_config_text = frozen_config.normalized_text;
  snapshot_payload.provenance =
      makeGravityAwareProvenanceRecord(frozen_config, config);
  report.snapshot_path = report.run_directory / formatIndexedRankedFileStem(config.output.output_stem, integrator_state.step_index, gravity_callback.runtimeTopology().world_size, gravity_callback.runtimeTopology().world_rank);
  io::writeGadgetArepoSnapshotHdf5(report.snapshot_path, snapshot_payload);
  report.snapshot_roundtrip_executed = true;
  const io::SnapshotReadResult snapshot_read = io::readGadgetArepoSnapshotHdf5(report.snapshot_path, config);
  report.snapshot_roundtrip_ok = snapshot_read.state.particles.size() == state.particles.size();
  profiler.recordEvent(core::RuntimeEvent{
      .event_kind = "snapshot.write.complete",
      .severity = report.snapshot_roundtrip_ok ? core::RuntimeEventSeverity::kInfo : core::RuntimeEventSeverity::kWarning,
      .subsystem = "io.snapshot",
      .step_index = integrator_state.step_index,
      .simulation_time_code = integrator_state.current_time_code,
      .scale_factor = integrator_state.current_scale_factor,
      .message = "snapshot output written and verified",
      .payload = {{"path", report.snapshot_path.string()}},
  });

  if (config.output.write_restarts) {
    io::RestartWritePayload restart_payload;
    restart_payload.state = &state;
    restart_payload.integrator_state = &integrator_state;
    restart_payload.scheduler = &scheduler;
    restart_payload.provenance =
        makeGravityAwareProvenanceRecord(frozen_config, config);
    restart_payload.normalized_config_text = frozen_config.normalized_text;
    restart_payload.normalized_config_hash_hex = frozen_config.provenance.config_hash_hex;
    restart_payload.distributed_gravity_state.schema_version = 2;
    restart_payload.distributed_gravity_state.decomposition_epoch = integrator_state.step_index;
    restart_payload.distributed_gravity_state.world_size = gravity_callback.runtimeTopology().world_size;
    restart_payload.distributed_gravity_state.pm_grid_nx = gravity_callback.pmGridShape().nx;
    restart_payload.distributed_gravity_state.pm_grid_ny = gravity_callback.pmGridShape().ny;
    restart_payload.distributed_gravity_state.pm_grid_nz = gravity_callback.pmGridShape().nz;
    restart_payload.distributed_gravity_state.pm_decomposition_mode =
        pmDecompositionModeName(config.numerics.treepm_pm_decomposition_mode);
    restart_payload.distributed_gravity_state.gravity_kick_opportunity = gravity_callback.gravityKickOpportunity();
    restart_payload.distributed_gravity_state.pm_update_cadence_steps =
        static_cast<std::uint64_t>(gravity_callback.pmCadenceSteps());
    restart_payload.distributed_gravity_state.long_range_field_version =
        gravity_callback.longRangeFieldVersion();
    restart_payload.distributed_gravity_state.last_long_range_refresh_opportunity =
        gravity_callback.lastLongRangeRefreshOpportunity();
    restart_payload.distributed_gravity_state.long_range_field_built_step_index =
        gravity_callback.lastLongRangeRefreshStepIndex();
    restart_payload.distributed_gravity_state.long_range_field_built_scale_factor =
        gravity_callback.lastLongRangeRefreshScaleFactor();
    restart_payload.distributed_gravity_state.long_range_restart_policy = "deterministic_rebuild";
    restart_payload.distributed_gravity_state.owning_rank_by_item.reserve(state.particle_sidecar.owning_rank.size());
    for (const std::uint32_t owner : state.particle_sidecar.owning_rank) {
      restart_payload.distributed_gravity_state.owning_rank_by_item.push_back(static_cast<int>(owner));
    }
    const std::size_t world_size = static_cast<std::size_t>(gravity_callback.runtimeTopology().world_size);
    restart_payload.distributed_gravity_state.pm_slab_begin_x_by_rank.resize(world_size, 0);
    restart_payload.distributed_gravity_state.pm_slab_end_x_by_rank.resize(world_size, 0);
    for (std::size_t rank = 0; rank < world_size; ++rank) {
      const parallel::PmSlabRange range = parallel::pmOwnedXRangeForRank(
          gravity_callback.pmGridShape().nx,
          static_cast<int>(world_size),
          static_cast<int>(rank));
      restart_payload.distributed_gravity_state.pm_slab_begin_x_by_rank[rank] = range.begin_x;
      restart_payload.distributed_gravity_state.pm_slab_end_x_by_rank[rank] = range.end_x;
    }
    restart_payload.provenance.gravity_treepm_decomposition_epoch =
        restart_payload.distributed_gravity_state.decomposition_epoch;
    restart_payload.provenance.gravity_treepm_restart_world_size =
        restart_payload.distributed_gravity_state.world_size;
    restart_payload.provenance.gravity_treepm_restart_pm_grid =
        std::to_string(restart_payload.distributed_gravity_state.pm_grid_nx) + "x" +
        std::to_string(restart_payload.distributed_gravity_state.pm_grid_ny) + "x" +
        std::to_string(restart_payload.distributed_gravity_state.pm_grid_nz);
    restart_payload.provenance.gravity_treepm_restart_slab_signature =
        pmSlabSignature(restart_payload.distributed_gravity_state);
    restart_payload.provenance.gravity_treepm_restart_kick_opportunity =
        restart_payload.distributed_gravity_state.gravity_kick_opportunity;
    restart_payload.provenance.gravity_treepm_restart_field_version =
        restart_payload.distributed_gravity_state.long_range_field_version;
    restart_payload.provenance.gravity_treepm_long_range_restart_policy =
        restart_payload.distributed_gravity_state.long_range_restart_policy;

    report.restart_path = report.run_directory / formatIndexedRankedFileStem(config.output.restart_stem, integrator_state.step_index, gravity_callback.runtimeTopology().world_size, gravity_callback.runtimeTopology().world_rank);
    io::writeRestartCheckpointHdf5(report.restart_path, restart_payload);
    report.restart_roundtrip_executed = true;
    const io::RestartReadResult restart_read = io::readRestartCheckpointHdf5(report.restart_path);
    const auto compatibility = parallel::evaluateDistributedRestartCompatibility(
        restart_read.distributed_gravity_state,
        gravity_callback.runtimeTopology());
    const bool restart_rank_qualified_name =
        gravity_callback.runtimeTopology().world_size == 1 ||
        report.restart_path.filename().string().find("_rank") != std::string::npos;
    report.restart_roundtrip_ok = restart_read.state.particles.size() == state.particles.size() &&
        restart_read.state.cells.size() == state.cells.size() &&
        restart_read.scheduler_state.current_tick == scheduler.currentTick() &&
        restart_read.state.particle_sidecar.particle_id == state.particle_sidecar.particle_id &&
        restart_read.state.particle_sidecar.owning_rank == state.particle_sidecar.owning_rank &&
        computeParticleIdSum(restart_read.state) == computeParticleIdSum(state) &&
        computeParticleIdXor(restart_read.state) == computeParticleIdXor(state) &&
        restart_read.distributed_gravity_state.owning_rank_by_item.size() == state.particle_sidecar.owning_rank.size() &&
        restart_read.distributed_gravity_state.owning_rank_by_item == restart_payload.distributed_gravity_state.owning_rank_by_item &&
        restart_read.distributed_gravity_state.pm_slab_begin_x_by_rank == restart_payload.distributed_gravity_state.pm_slab_begin_x_by_rank &&
        restart_read.distributed_gravity_state.pm_slab_end_x_by_rank == restart_payload.distributed_gravity_state.pm_slab_end_x_by_rank &&
        restart_read.distributed_gravity_state.gravity_kick_opportunity == restart_payload.distributed_gravity_state.gravity_kick_opportunity &&
        restart_read.distributed_gravity_state.pm_update_cadence_steps == restart_payload.distributed_gravity_state.pm_update_cadence_steps &&
        restart_read.distributed_gravity_state.long_range_field_version == restart_payload.distributed_gravity_state.long_range_field_version &&
        restart_read.distributed_gravity_state.last_long_range_refresh_opportunity ==
            restart_payload.distributed_gravity_state.last_long_range_refresh_opportunity &&
        restart_read.distributed_gravity_state.long_range_field_built_step_index ==
            restart_payload.distributed_gravity_state.long_range_field_built_step_index &&
        std::abs(
            restart_read.distributed_gravity_state.long_range_field_built_scale_factor -
            restart_payload.distributed_gravity_state.long_range_field_built_scale_factor) <= 1.0e-12 &&
        restart_read.distributed_gravity_state.long_range_restart_policy ==
            restart_payload.distributed_gravity_state.long_range_restart_policy &&
        restart_rank_qualified_name &&
        compatibility.compatible();
    profiler.recordEvent(core::RuntimeEvent{
        .event_kind = "restart.write.complete",
        .severity = report.restart_roundtrip_ok ? core::RuntimeEventSeverity::kInfo : core::RuntimeEventSeverity::kWarning,
        .subsystem = "io.restart",
        .step_index = integrator_state.step_index,
        .simulation_time_code = integrator_state.current_time_code,
        .scale_factor = integrator_state.current_scale_factor,
        .message = "restart checkpoint written and verified",
        .payload = {{"path", report.restart_path.string()}},
    });
  }
#endif
}

}  // namespace

ReferenceWorkflowRunner::ReferenceWorkflowRunner(core::FrozenConfig frozen_config)
    : m_frozen_config(std::move(frozen_config)) {}

const core::FrozenConfig& ReferenceWorkflowRunner::frozenConfig() const noexcept {
  return m_frozen_config;
}

ReferenceWorkflowReport ReferenceWorkflowRunner::run(const ReferenceWorkflowOptions& options) const {
  return runImpl(nullptr, options);
}

ReferenceWorkflowReport ReferenceWorkflowRunner::run(
    const std::filesystem::path& output_root_override,
    const ReferenceWorkflowOptions& options) const {
  return runImpl(&output_root_override, options);
}

ReferenceWorkflowReport ReferenceWorkflowRunner::runImpl(
    const std::filesystem::path* output_root_override,
    const ReferenceWorkflowOptions& options) const {
  const core::SimulationConfig& config = m_frozen_config.config;

  ReferenceWorkflowReport report;
  report.run_directory = computeRunDirectory(config, output_root_override);
  report.config_compatible = true;
  report.schema_compatible =
      config.schema_version == 1 && io::gadgetArepoSchemaMap().schema_version == 2 &&
      io::isRestartSchemaCompatible(io::restartSchema().version);

  core::ProfilerSession profiler(true);
  profiler.recordEvent(core::RuntimeEvent{
      .event_kind = "config.freeze",
      .severity = core::RuntimeEventSeverity::kInfo,
      .subsystem = "core.config",
      .step_index = options.step_index,
      .simulation_time_code = config.numerics.time_begin_code,
      .scale_factor = 1.0,
      .message = "frozen configuration accepted for runtime workflow",
      .payload = {{"schema_version", std::to_string(config.schema_version)},
                  {"config_hash_hex", m_frozen_config.provenance.config_hash_hex},
                  {"source_name", m_frozen_config.provenance.source_name}},
  });

  try {
    ensureRunDirectory(report.run_directory);
    core::writeNormalizedConfigSnapshot(m_frozen_config, report.run_directory);
    report.normalized_config_snapshot_path = report.run_directory / "normalized_config.param.txt";
    report.normalized_config_snapshot_written = true;

    if (!report.schema_compatible) {
      throw std::runtime_error("runtime workflow schema compatibility validation failed");
    }

#if !COSMOSIM_ENABLE_HDF5
    if (options.write_outputs) {
      throw std::runtime_error(
          "this build lacks HDF5 support, so the config-driven runtime cannot emit snapshots/restarts; reconfigure with COSMOSIM_ENABLE_HDF5=ON or use a no-output internal test path");
    }
#endif

    const core::ModePolicy mode_policy = core::buildModePolicy(config.mode);
    core::validateModePolicy(config, mode_policy);

    io::IcReadResult ic_result = loadInitialConditions(m_frozen_config);
    core::SimulationState state = std::move(ic_result.state);
    finalizeStateMetadata(m_frozen_config, state);
    parallel::MpiContext mpi_context;
    report.world_size = mpi_context.worldSize();
    report.world_rank = mpi_context.worldRank();
    seedParticleOwnershipFromPmSlabs(
        state,
        static_cast<std::size_t>(config.numerics.treepm_pm_grid_nx),
        mpi_context.worldSize(),
        config.cosmology.box_size_x_mpc_comoving);
    applyInitialGravityAwareDecomposition(state, config, mpi_context.worldSize(), mpi_context.worldRank());

    report.local_particle_count = static_cast<std::uint64_t>(state.particles.size());
    report.local_cell_count = static_cast<std::uint64_t>(state.cells.size());
    report.local_particle_id_sum = computeParticleIdSum(state);
    report.local_particle_id_xor = computeParticleIdXor(state);
    report.global_particle_count = mpi_context.allreduceSumUint64(report.local_particle_count);
    report.global_cell_count = mpi_context.allreduceSumUint64(report.local_cell_count);
    report.global_particle_id_sum = mpi_context.allreduceSumUint64(report.local_particle_id_sum);
    report.global_particle_id_xor = mpi_context.allreduceXorUint64(report.local_particle_id_xor);
    profiler.recordEvent(core::RuntimeEvent{
        .event_kind = "parallel.ownership.partition_summary",
        .severity = core::RuntimeEventSeverity::kInfo,
        .subsystem = "parallel.ownership",
        .step_index = options.step_index,
        .simulation_time_code = config.numerics.time_begin_code,
        .scale_factor = 1.0,
        .message = "initial ownership-partition summary recorded",
        .payload = {{"world_size", std::to_string(report.world_size)},
                    {"world_rank", std::to_string(report.world_rank)},
                    {"local_particle_count", std::to_string(report.local_particle_count)},
                    {"global_particle_count", std::to_string(report.global_particle_count)},
                    {"local_cell_count", std::to_string(report.local_cell_count)},
                    {"global_cell_count", std::to_string(report.global_cell_count)},
                    {"local_particle_id_sum", std::to_string(report.local_particle_id_sum)},
                    {"global_particle_id_sum", std::to_string(report.global_particle_id_sum)},
                    {"local_particle_id_xor", std::to_string(report.local_particle_id_xor)},
                    {"global_particle_id_xor", std::to_string(report.global_particle_id_xor)}},
    });

    core::CosmologyBackgroundConfig background_config;
    background_config.hubble_param = config.cosmology.hubble_param;
    background_config.omega_matter = config.cosmology.omega_matter;
    background_config.omega_lambda = config.cosmology.omega_lambda;
    const std::optional<core::LambdaCdmBackground> background = mode_policy.cosmological_comoving_frame
        ? std::optional<core::LambdaCdmBackground>(core::LambdaCdmBackground(background_config))
        : std::nullopt;

    const std::uint32_t scheduler_count =
        static_cast<std::uint32_t>(std::max(state.particles.size(), state.cells.size()));
    core::HierarchicalTimeBinScheduler scheduler(
        static_cast<std::uint8_t>(std::max(0, std::min(config.numerics.hierarchical_max_rung, 12))));
    initializeSchedulerBins(state, scheduler);
    syncTimeBinsFromScheduler(scheduler, state);

    core::IntegratorState integrator_state;
    integrator_state.step_index = options.step_index;
    integrator_state.current_time_code = config.numerics.time_begin_code;
    integrator_state.current_scale_factor =
        (background.has_value() && config.numerics.time_begin_code > 0.0) ? config.numerics.time_begin_code : 1.0;
    integrator_state.dt_time_code = options.dt_time_code > 0.0
        ? options.dt_time_code
        : std::max(
              1.0e-6,
              (config.numerics.time_end_code - config.numerics.time_begin_code) /
                  static_cast<double>(std::max(config.numerics.max_global_steps, 1)));
    integrator_state.time_bins.hierarchical_enabled = true;
    integrator_state.time_bins.max_bin = scheduler.maxBin();

    core::StepOrchestrator orchestrator;
    StageAuditCallback stage_audit(&report.stage_sequence);
    DriftCallback drift_callback;
    GravityStageCallback gravity_callback(config, mode_policy);
    report.treepm_pm_grid = gravity_callback.pmGridSize();
    report.treepm_update_cadence_steps = gravity_callback.pmCadenceSteps();
    profiler.recordEvent(core::RuntimeEvent{
        .event_kind = "gravity.treepm_setup",
        .severity = core::RuntimeEventSeverity::kInfo,
        .subsystem = "gravity.treepm",
        .step_index = integrator_state.step_index,
        .simulation_time_code = integrator_state.current_time_code,
        .scale_factor = integrator_state.current_scale_factor,
        .message = "TreePM runtime configuration initialized",
        .payload = {
            {"pm_grid", std::to_string(config.numerics.treepm_pm_grid_nx) + "x" +
                    std::to_string(config.numerics.treepm_pm_grid_ny) + "x" +
                    std::to_string(config.numerics.treepm_pm_grid_nz)},
            {"pm_assignment_scheme", treePmAssignmentSchemeName(config.numerics.treepm_assignment_scheme)},
            {"pm_window_deconvolution", config.numerics.treepm_enable_window_deconvolution ? "true" : "false"},
            {"asmth_cells", std::to_string(config.numerics.treepm_asmth_cells)},
            {"rcut_cells", std::to_string(config.numerics.treepm_rcut_cells)},
            {"mesh_spacing_x_mpc_comoving", std::to_string(config.cosmology.box_size_x_mpc_comoving / static_cast<double>(config.numerics.treepm_pm_grid_nx))},
            {"mesh_spacing_y_mpc_comoving", std::to_string(config.cosmology.box_size_y_mpc_comoving / static_cast<double>(config.numerics.treepm_pm_grid_ny))},
            {"mesh_spacing_z_mpc_comoving", std::to_string(config.cosmology.box_size_z_mpc_comoving / static_cast<double>(config.numerics.treepm_pm_grid_nz))},
            {"split_scale_mpc_comoving", std::to_string(config.numerics.treepm_asmth_cells *
                std::cbrt((config.cosmology.box_size_x_mpc_comoving / static_cast<double>(config.numerics.treepm_pm_grid_nx)) *
                          (config.cosmology.box_size_y_mpc_comoving / static_cast<double>(config.numerics.treepm_pm_grid_ny)) *
                          (config.cosmology.box_size_z_mpc_comoving / static_cast<double>(config.numerics.treepm_pm_grid_nz))))},
            {"cutoff_radius_mpc_comoving", std::to_string(config.numerics.treepm_rcut_cells *
                std::cbrt((config.cosmology.box_size_x_mpc_comoving / static_cast<double>(config.numerics.treepm_pm_grid_nx)) *
                          (config.cosmology.box_size_y_mpc_comoving / static_cast<double>(config.numerics.treepm_pm_grid_ny)) *
                          (config.cosmology.box_size_z_mpc_comoving / static_cast<double>(config.numerics.treepm_pm_grid_nz))))},
            {"pm_update_cadence_steps", std::to_string(config.numerics.treepm_update_cadence_steps)},
            {"zoom_long_range_strategy",
                config.mode.zoom_long_range_strategy == core::ZoomLongRangeStrategy::kDisabled
                    ? "disabled"
                    : "global_coarse_plus_focused_highres_correction"},
            {"zoom_region_center_x_mpc_comoving", std::to_string(config.mode.zoom_region_center_x_mpc_comoving)},
            {"zoom_region_center_y_mpc_comoving", std::to_string(config.mode.zoom_region_center_y_mpc_comoving)},
            {"zoom_region_center_z_mpc_comoving", std::to_string(config.mode.zoom_region_center_z_mpc_comoving)},
            {"zoom_region_radius_mpc_comoving", std::to_string(config.mode.zoom_region_radius_mpc_comoving)},
            {"zoom_focused_pm_grid", std::to_string(config.mode.zoom_focused_pm_grid_nx) + "x" +
                std::to_string(config.mode.zoom_focused_pm_grid_ny) + "x" +
                std::to_string(config.mode.zoom_focused_pm_grid_nz)},
            {"zoom_contamination_radius_mpc_comoving", std::to_string(config.mode.zoom_contamination_radius_mpc_comoving)},
            {"softening_policy", "comoving_fixed"},
            {"softening_kernel", "plummer"},
            {"softening_epsilon_kpc_comoving", std::to_string(config.numerics.gravity_softening_kpc_comoving)},
            {"pm_fft_backend", gravity::PmSolver::fftBackendName()},
        },
    });
    HydroStageCallback hydro_callback(config, mode_policy, gravity_callback);
    analysis::DiagnosticsCallback diagnostics_callback(config);
    physics::StarFormationCallback star_formation_callback(
        physics::StarFormationModel(physics::makeStarFormationConfig(config.physics)));
    physics::BlackHoleAgnCallback bh_callback(
        physics::BlackHoleAgnModel(physics::makeBlackHoleAgnConfig(config.physics)));

    orchestrator.registerCallback(stage_audit);
    orchestrator.registerCallback(drift_callback);
    orchestrator.registerCallback(gravity_callback);
    orchestrator.registerCallback(hydro_callback);
    orchestrator.registerCallback(star_formation_callback);
    orchestrator.registerCallback(bh_callback);
    orchestrator.registerCallback(diagnostics_callback);

    const std::uint64_t target_step_index =
        integrator_state.step_index + static_cast<std::uint64_t>(std::max(config.numerics.max_global_steps, 0));
    while (integrator_state.step_index < target_step_index &&
           integrator_state.current_time_code < config.numerics.time_end_code) {
      const std::span<const std::uint32_t> active_scheduler_elements = scheduler.beginSubstep();
      std::vector<std::uint32_t> active_particles;
      std::vector<std::uint32_t> active_cells;
      active_particles.reserve(active_scheduler_elements.size());
      active_cells.reserve(active_scheduler_elements.size());
      for (const std::uint32_t element_index : active_scheduler_elements) {
        if (element_index < state.particles.size()) {
          active_particles.push_back(element_index);
        }
        if (element_index < state.cells.size()) {
          active_cells.push_back(element_index);
        }
      }
      if (active_particles.empty() && active_cells.empty()) {
        scheduler.endSubstep();
        continue;
      }

      core::TransientStepWorkspace workspace;
      const core::ActiveSetDescriptor active_set{
          .particle_indices = active_particles,
          .cell_indices = active_cells,
          .particles_are_subset = active_particles.size() < state.particles.size(),
          .cells_are_subset = active_cells.size() < state.cells.size(),
      };

      orchestrator.executeSingleStep(
          state,
          integrator_state,
          active_set,
          background.has_value() ? &background.value() : nullptr,
          &workspace,
          &mode_policy,
          &profiler);

      state.metadata.step_index = integrator_state.step_index;
      state.metadata.scale_factor = integrator_state.current_scale_factor;
      scheduler.endSubstep();
      syncTimeBinsFromScheduler(scheduler, state);
      maybeWriteOutputs(m_frozen_config, config, state, integrator_state, scheduler, gravity_callback, report, profiler, options.write_outputs);
    }

    report.completed_steps = integrator_state.step_index - options.step_index;
    report.final_state_digest = computeStateDigest(state, integrator_state);
    report.local_particle_count = static_cast<std::uint64_t>(state.particles.size());
    report.local_cell_count = static_cast<std::uint64_t>(state.cells.size());
    report.local_particle_id_sum = computeParticleIdSum(state);
    report.local_particle_id_xor = computeParticleIdXor(state);
    report.global_particle_count = mpi_context.allreduceSumUint64(report.local_particle_count);
    report.global_cell_count = mpi_context.allreduceSumUint64(report.local_cell_count);
    report.global_particle_id_sum = mpi_context.allreduceSumUint64(report.local_particle_id_sum);
    report.global_particle_id_xor = mpi_context.allreduceXorUint64(report.local_particle_id_xor);
    report.treepm_long_range_refresh_count = gravity_callback.longRangeRefreshCount();
    report.treepm_long_range_reuse_count = gravity_callback.longRangeReuseCount();
    report.treepm_cadence_records.assign(
        gravity_callback.cadenceRecords().begin(),
        gravity_callback.cadenceRecords().end());
    report.canonical_stage_order = repeatedCanonicalOrder(report.stage_sequence);
    profiler.recordEvent(core::RuntimeEvent{
        .event_kind = "run.complete",
        .severity = core::RuntimeEventSeverity::kInfo,
        .subsystem = "workflows.reference",
        .step_index = integrator_state.step_index,
        .simulation_time_code = integrator_state.current_time_code,
        .scale_factor = integrator_state.current_scale_factor,
        .message = "runtime workflow completed",
        .payload = {{"completed_steps", std::to_string(report.completed_steps)},
                    {"run_directory", report.run_directory.string()},
                    {"treepm_update_cadence_steps", std::to_string(report.treepm_update_cadence_steps)},
                    {"treepm_long_range_refresh_count", std::to_string(report.treepm_long_range_refresh_count)},
                    {"treepm_long_range_reuse_count", std::to_string(report.treepm_long_range_reuse_count)},
                    {"local_particle_count", std::to_string(report.local_particle_count)},
                    {"global_particle_count", std::to_string(report.global_particle_count)},
                    {"local_particle_id_sum", std::to_string(report.local_particle_id_sum)},
                    {"global_particle_id_sum", std::to_string(report.global_particle_id_sum)},
                    {"local_particle_id_xor", std::to_string(report.local_particle_id_xor)},
                    {"global_particle_id_xor", std::to_string(report.global_particle_id_xor)},
                    {"final_state_digest", std::to_string(report.final_state_digest)}},
    });
    flushCommonArtifacts(m_frozen_config, profiler, report);
    return report;
  } catch (const std::exception& ex) {
    profiler.recordEvent(core::RuntimeEvent{
        .event_kind = "run.failure",
        .severity = core::RuntimeEventSeverity::kFatal,
        .subsystem = "workflows.reference",
        .step_index = std::nullopt,
        .simulation_time_code = std::nullopt,
        .scale_factor = std::nullopt,
        .message = "runtime workflow failed",
        .payload = {{"error", ex.what()}, {"run_directory", report.run_directory.string()}},
    });
    try {
      flushCommonArtifacts(m_frozen_config, profiler, report);
    } catch (...) {
    }
    throw;
  }
}

}  // namespace cosmosim::workflows
