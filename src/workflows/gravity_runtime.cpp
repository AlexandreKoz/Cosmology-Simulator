#include "cosmosim/workflows/gravity_runtime.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <memory>
#include <numeric>
#include <optional>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "cosmosim/core/checked_arithmetic.hpp"
#include "cosmosim/core/constants.hpp"
#include "cosmosim/core/cuda_runtime.hpp"
#include "cosmosim/core/time_scheduler.hpp"
#include "cosmosim/core/units.hpp"
#include "cosmosim/gravity/gravity_memory.hpp"
#include "cosmosim/gravity/gravity_source_snapshot.hpp"
#include "cosmosim/workflows/runtime_services.hpp"
#include "workflows/internal/runtime_stage_resource_access.hpp"
#include "workflows/internal/particle_ghost_runtime.hpp"
#include "workflows/internal/gas_cell_ownership.hpp"

#if COSMOSIM_ENABLE_HDF5
#include <hdf5.h>
#endif

namespace cosmosim::workflows {
namespace {

[[nodiscard]] std::string formatRuntimeDouble(double value) {
  std::ostringstream stream;
  stream << std::scientific
         << std::setprecision(std::numeric_limits<double>::max_digits10)
         << value;
  return stream.str();
}

constexpr std::array<core::IntegrationStage, 1> k_hydro_stage = {core::IntegrationStage::kHydroUpdate};

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

[[nodiscard]] gravity::TreeOpeningCriterion toTreeOpeningCriterion(
    core::TreePmOpeningCriterion opening_criterion) {
  switch (opening_criterion) {
    case core::TreePmOpeningCriterion::kGeometric:
      return gravity::TreeOpeningCriterion::kBarnesHutGeometric;
    case core::TreePmOpeningCriterion::kComDistance:
      return gravity::TreeOpeningCriterion::kBarnesHutComDistance;
    case core::TreePmOpeningCriterion::kRelativeForceError:
      return gravity::TreeOpeningCriterion::kRelativeForceError;
  }
  throw std::runtime_error("unhandled TreePM opening criterion enum value");
}

[[nodiscard]] std::string treePmOpeningCriterionName(
    core::TreePmOpeningCriterion opening_criterion) {
  switch (opening_criterion) {
    case core::TreePmOpeningCriterion::kGeometric:
      return "geometric";
    case core::TreePmOpeningCriterion::kComDistance:
      return "com_distance";
    case core::TreePmOpeningCriterion::kRelativeForceError:
      return "relative_force_error";
  }
  throw std::runtime_error("unhandled TreePM opening criterion enum value");
}



gravity::TreePmCoordinator makeRuntimeAwareTreePmCoordinator(
    const core::SimulationConfig& config,
    const gravity::PmGridShape& pm_grid_shape,
    const RuntimeServices& services) {
  const parallel::MpiContext& mpi_context = services.mpi_context;
  mpi_context.validateExpectedWorldSizeOrThrow(config.parallel.mpi_ranks_expected);
  const auto layout = parallel::makePmSlabLayout(
      pm_grid_shape.nx,
      pm_grid_shape.ny,
      pm_grid_shape.nz,
      mpi_context.worldSize(),
      mpi_context.worldRank());
  return gravity::TreePmCoordinator(pm_grid_shape, layout, mpi_context);
}


[[nodiscard]] double newtonGCodeFromUnits(const core::UnitSystem& units) {
  return core::newtonGravitationalConstantCode(units);
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

[[nodiscard]] gravity::TreeSofteningSpeciesPolicy speciesSofteningByTag(
    const core::SimulationConfig& config) {
  gravity::TreeSofteningSpeciesPolicy values{};
  values.enabled = config.numerics.gravity_softening_gas_kpc_comoving > 0.0 ||
      config.numerics.gravity_softening_dark_matter_kpc_comoving > 0.0 ||
      config.numerics.gravity_softening_star_kpc_comoving > 0.0 ||
      config.numerics.gravity_softening_black_hole_kpc_comoving > 0.0 ||
      config.numerics.gravity_softening_tracer_kpc_comoving > 0.0;
  values.epsilon_comoving_by_species.fill(config.numerics.gravity_softening_kpc_comoving * 1.0e-3);
  values.epsilon_comoving_by_species[static_cast<std::size_t>(core::ParticleSpecies::kDarkMatter)] =
      (config.numerics.gravity_softening_dark_matter_kpc_comoving > 0.0)
      ? (config.numerics.gravity_softening_dark_matter_kpc_comoving * 1.0e-3)
      : values.epsilon_comoving_by_species[static_cast<std::size_t>(core::ParticleSpecies::kDarkMatter)];
  values.epsilon_comoving_by_species[static_cast<std::size_t>(core::ParticleSpecies::kGas)] =
      (config.numerics.gravity_softening_gas_kpc_comoving > 0.0)
      ? (config.numerics.gravity_softening_gas_kpc_comoving * 1.0e-3)
      : values.epsilon_comoving_by_species[static_cast<std::size_t>(core::ParticleSpecies::kGas)];
  values.epsilon_comoving_by_species[static_cast<std::size_t>(core::ParticleSpecies::kStar)] =
      (config.numerics.gravity_softening_star_kpc_comoving > 0.0)
      ? (config.numerics.gravity_softening_star_kpc_comoving * 1.0e-3)
      : values.epsilon_comoving_by_species[static_cast<std::size_t>(core::ParticleSpecies::kStar)];
  values.epsilon_comoving_by_species[static_cast<std::size_t>(core::ParticleSpecies::kBlackHole)] =
      (config.numerics.gravity_softening_black_hole_kpc_comoving > 0.0)
      ? (config.numerics.gravity_softening_black_hole_kpc_comoving * 1.0e-3)
      : values.epsilon_comoving_by_species[static_cast<std::size_t>(core::ParticleSpecies::kBlackHole)];
  values.epsilon_comoving_by_species[static_cast<std::size_t>(core::ParticleSpecies::kTracer)] =
      (config.numerics.gravity_softening_tracer_kpc_comoving > 0.0)
      ? (config.numerics.gravity_softening_tracer_kpc_comoving * 1.0e-3)
      : values.epsilon_comoving_by_species[static_cast<std::size_t>(core::ParticleSpecies::kTracer)];
  return values;
}

[[nodiscard]] bool hasSpeciesSpecificSoftening(const core::SimulationConfig& config) {
  return config.numerics.gravity_softening_gas_kpc_comoving > 0.0 ||
      config.numerics.gravity_softening_dark_matter_kpc_comoving > 0.0 ||
      config.numerics.gravity_softening_star_kpc_comoving > 0.0 ||
      config.numerics.gravity_softening_black_hole_kpc_comoving > 0.0 ||
      config.numerics.gravity_softening_tracer_kpc_comoving > 0.0;
}

[[nodiscard]] std::string describeSofteningPolicy(
    const core::SimulationConfig& config,
    const core::SimulationState* state = nullptr) {
  const bool per_species = hasSpeciesSpecificSoftening(config);
  const bool per_particle = state != nullptr && std::any_of(
      state->particle_sidecar.has_gravity_softening_override.begin(),
      state->particle_sidecar.has_gravity_softening_override.end(),
      [](std::uint8_t flag) { return flag != 0U; });
  if (per_particle && per_species) {
    return "comoving_species_plus_particle_override";
  }
  if (per_particle) {
    return "comoving_particle_override";
  }
  if (per_species) {
    return "comoving_species";
  }
  return "comoving_fixed";
}

void maybeInitializeParticleSofteningFromSpeciesPolicy(
    core::SimulationState& state,
    const core::SimulationConfig& config) {
  if (!state.particle_sidecar.gravity_softening_comoving.empty()) {
    return;
  }
  if (!hasSpeciesSpecificSoftening(config)) {
    return;
  }
  const auto by_species = speciesSofteningByTag(config);
  state.particle_sidecar.gravity_softening_comoving.resize(state.particles.size(), 0.0);
  for (std::size_t i = 0; i < state.particles.size(); ++i) {
    const std::size_t species_tag = static_cast<std::size_t>(state.particle_sidecar.species_tag[i]);
    if (species_tag < by_species.epsilon_comoving_by_species.size()) {
      state.particle_sidecar.gravity_softening_comoving[i] = by_species.epsilon_comoving_by_species[species_tag];
    } else {
      state.particle_sidecar.gravity_softening_comoving[i] = config.numerics.gravity_softening_kpc_comoving * 1.0e-3;
    }
  }
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

[[nodiscard]] bool hasHdf5Extension(const std::filesystem::path& path) {
  const std::string ext = path.extension().string();
  return ext == ".h5" || ext == ".hdf5";
}

[[nodiscard]] std::unordered_set<std::uint64_t> loadZoomParticleIdsFromText(const std::filesystem::path& path) {
  std::ifstream in(path);
  if (!in) {
    throw std::runtime_error("failed to open zoom region file: " + path.string());
  }
  std::unordered_set<std::uint64_t> ids;
  std::string line;
  while (std::getline(in, line)) {
    const auto comment = line.find('#');
    if (comment != std::string::npos) {
      line.erase(comment);
    }
    std::istringstream stream(line);
    std::uint64_t id = 0;
    while (stream >> id) {
      ids.insert(id);
    }
  }
  return ids;
}

#if COSMOSIM_ENABLE_HDF5
[[nodiscard]] std::unordered_set<std::uint64_t> loadZoomParticleIdsFromHdf5(const std::filesystem::path& path) {
  hid_t file = H5Fopen(path.string().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file < 0) {
    throw std::runtime_error("failed to open zoom region HDF5 file: " + path.string());
  }
  auto close_file = [&]() { H5Fclose(file); };
  const char* candidates[] = {"/Region/ParticleIDs", "/ParticleIDs", "Region/ParticleIDs", "ParticleIDs"};
  hid_t dataset = -1;
  for (const char* candidate : candidates) {
    if (H5Lexists(file, candidate, H5P_DEFAULT) > 0) {
      dataset = H5Dopen2(file, candidate, H5P_DEFAULT);
      if (dataset >= 0) {
        break;
      }
    }
  }
  if (dataset < 0) {
    close_file();
    throw std::runtime_error("zoom region HDF5 file lacks ParticleIDs dataset: " + path.string());
  }
  hid_t space = H5Dget_space(dataset);
  if (space < 0) {
    H5Dclose(dataset); close_file();
    throw std::runtime_error("failed to query zoom region HDF5 dataspace");
  }
  hsize_t dims[1] = {0};
  if (H5Sget_simple_extent_ndims(space) != 1 || H5Sget_simple_extent_dims(space, dims, nullptr) != 1) {
    H5Sclose(space); H5Dclose(dataset); close_file();
    throw std::runtime_error("zoom region ParticleIDs dataset must be 1D");
  }
  std::vector<std::uint64_t> values(static_cast<std::size_t>(dims[0]), 0ULL);
  const hid_t dtype = H5Dget_type(dataset);
  const H5T_class_t cls = H5Tget_class(dtype);
  if (cls != H5T_INTEGER) {
    H5Tclose(dtype); H5Sclose(space); H5Dclose(dataset); close_file();
    throw std::runtime_error("zoom region ParticleIDs dataset must be integer typed");
  }
  if (!values.empty() && H5Dread(dataset, H5T_NATIVE_UINT64, H5S_ALL, H5S_ALL, H5P_DEFAULT, values.data()) < 0) {
    H5Tclose(dtype); H5Sclose(space); H5Dclose(dataset); close_file();
    throw std::runtime_error("failed to read zoom region ParticleIDs dataset");
  }
  H5Tclose(dtype); H5Sclose(space); H5Dclose(dataset); close_file();
  return std::unordered_set<std::uint64_t>(values.begin(), values.end());
}
#endif


class GravityRuntimeImpl final : public GravityRuntime {
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
    kForceRefresh = 1,
    kKickPost = 2,
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

  GravityRuntimeImpl(
      const core::SimulationConfig& config,
      const core::ModePolicy& mode_policy,
      const RuntimeServices& services,
      std::filesystem::path zoom_region_path = {})
      : m_services(services),
        m_config(config),
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
        m_tree_pm_coordinator(makeRuntimeAwareTreePmCoordinator(
            config, m_pm_grid_shape, services)),
        m_zoom_region_path(std::move(zoom_region_path)) {
    if (!m_pm_grid_shape.isValid()) {
      throw std::runtime_error("numerics.treepm_pm_grid_n{xyz} must all be > 0");
    }
    m_tree_pm_options.pm_options.box_size_x_mpc_comoving = config.cosmology.box_size_x_mpc_comoving;
    m_tree_pm_options.pm_options.box_size_y_mpc_comoving = config.cosmology.box_size_y_mpc_comoving;
    m_tree_pm_options.pm_options.box_size_z_mpc_comoving = config.cosmology.box_size_z_mpc_comoving;
    m_tree_pm_options.pm_options.box_size_mpc_comoving = config.cosmology.box_size_x_mpc_comoving;
    m_tree_pm_options.pm_options.scale_factor = 1.0;
    const core::UnitSystem runtime_units = core::makeUnitSystem(
        config.units.length_unit,
        config.units.mass_unit,
        config.units.velocity_unit);
    m_gravitational_constant_code =
        core::newtonGravitationalConstantCode(runtime_units);
    m_tree_pm_options.pm_options.gravitational_constant_code =
        m_gravitational_constant_code;
    m_tree_pm_options.pm_options.assignment_scheme =
        toPmAssignmentScheme(config.numerics.treepm_assignment_scheme);
    m_tree_pm_options.pm_options.enable_window_deconvolution =
        config.numerics.treepm_enable_window_deconvolution;
    m_tree_pm_options.pm_options.decomposition_mode = config.numerics.treepm_pm_decomposition_mode;
    m_tree_pm_options.pm_options.boundary_condition = mode_policy.gravity_boundary ==
            core::GravityBoundaryModel::kPeriodicPoisson
        ? gravity::PmBoundaryCondition::kPeriodic
        : gravity::PmBoundaryCondition::kIsolatedOpen;
    m_tree_pm_options.pm_options.isolated_open_root_workspace_limit_bytes =
        config.parallel.isolated_pm_root_workspace_limit_bytes;
    const parallel::MpiContext& mpi_context = m_services.mpi_context;
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

    m_tree_pm_options.tree_options.opening_criterion =
        toTreeOpeningCriterion(config.numerics.treepm_tree_opening_criterion);
    m_tree_pm_options.tree_options.opening_theta = config.numerics.treepm_tree_opening_theta;
    m_tree_pm_options.tree_options.relative_force_tolerance =
        config.numerics.treepm_tree_relative_force_tolerance;
    m_tree_pm_options.tree_options.relative_force_acceleration_floor_code =
        config.numerics.treepm_tree_relative_force_acceleration_floor;
    m_tree_pm_options.tree_options.gravitational_constant_code =
        m_gravitational_constant_code;
    m_tree_pm_options.tree_options.softening.kernel = gravity::TreeSofteningKernel::kPlummer;
    m_tree_pm_options.tree_options.softening.epsilon_comoving =
        config.numerics.gravity_softening_kpc_comoving * 1.0e-3;
    m_tree_pm_species_softening = speciesSofteningByTag(config);
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
    m_tree_pm_options.pm_options.routing_exchange_batch_bytes =
        config.numerics.treepm_pm_exchange_batch_bytes;
    m_tree_pm_options.zoom_high_res_allgather_limit_bytes =
        config.parallel.zoom_high_res_allgather_limit_bytes;
    m_pm_assignment_scheme = treePmAssignmentSchemeName(config.numerics.treepm_assignment_scheme);
    m_pm_backend = gravity::PmSolver::fftBackendName();
    if (m_runtime_topology.usesCuda()) {
      m_pm_backend += "+cuda_cic";
    }
  }

  [[nodiscard]] std::size_t pmGridSize() const noexcept { return m_pm_grid_shape.nx; }
  [[nodiscard]] const gravity::PmGridShape& pmGridShape() const noexcept { return m_pm_grid_shape; }
  [[nodiscard]] std::uint64_t longRangeRefreshCount() const noexcept { return m_long_range_refresh_count; }
  [[nodiscard]] std::uint64_t longRangeReuseCount() const noexcept { return m_long_range_reuse_count; }
  [[nodiscard]] int pmCadenceSteps() const noexcept { return static_cast<int>(m_pm_update_cadence_steps); }
  [[nodiscard]] double gravitationalConstantCode() const noexcept {
    return m_gravitational_constant_code;
  }
  [[nodiscard]] std::span<const ReferenceWorkflowReport::TreePmCadenceRecord> cadenceRecords() const noexcept {
    return m_cadence_records;
  }
  [[nodiscard]] const parallel::DistributedExecutionTopology& runtimeTopology() const noexcept { return m_runtime_topology; }
  [[nodiscard]] core::MemoryReport memoryReport() const {
    core::MemoryReportBuilder builder;
    const auto add = [&builder](
                         core::MemorySubsystem subsystem,
                         core::MemoryLifetime lifetime,
                         std::string label,
                         const auto& container) {
      const std::uint64_t capacity_bytes = core::ownedCapacityBytesForContainer(container);
      builder.addEntry(core::MemoryEntry{
          .subsystem = subsystem,
          .lifetime = lifetime,
          .memory_class = lifetime == core::MemoryLifetime::kPersistent
                              ? core::MemoryClass::kPersistentCache
                              : core::MemoryClass::kPhaseResident,
          .label = std::move(label),
          .current_size_bytes = core::currentSizeBytesForContainer(container),
          .owned_capacity_bytes = capacity_bytes,
          .high_water_bytes = capacity_bytes,
          .estimated_next_step_bytes = 0U,
          .uncertainty_note = "retained vector capacity is the observed allocation high-water; next-step requirement is not predicted",
      });
    };
    const auto active = [&](std::string label, const auto& container) {
      add(core::MemorySubsystem::kActiveSets, core::MemoryLifetime::kTransient,
          std::move(label), container);
    };
    const auto persistent = [&](std::string label, const auto& container) {
      add(core::MemorySubsystem::kSidecars, core::MemoryLifetime::kPersistent,
          std::move(label), container);
    };

    active("gravity_runtime.active_accel_x", m_active_accel_x);
    active("gravity_runtime.active_accel_y", m_active_accel_y);
    active("gravity_runtime.active_accel_z", m_active_accel_z);
    active("gravity_runtime.active_previous_acceleration", m_active_previous_acceleration_magnitude);
    persistent("gravity_runtime.particle_accel_x", m_particle_accel_x);
    persistent("gravity_runtime.particle_accel_y", m_particle_accel_y);
    persistent("gravity_runtime.particle_accel_z", m_particle_accel_z);
    persistent("gravity_runtime.cell_accel_x", m_cell_accel_x);
    persistent("gravity_runtime.cell_accel_y", m_cell_accel_y);
    persistent("gravity_runtime.cell_accel_z", m_cell_accel_z);
    active("gravity_runtime.active_slot_by_particle", m_active_slot_by_particle);
    active("gravity_runtime.active_slot_by_cell", m_active_slot_by_cell);
    active("gravity_runtime.owned_local_index_by_particle", m_owned_local_index_by_global);
    active("gravity_runtime.owned_local_index_by_cell", m_owned_local_index_by_cell);
    active("gravity_runtime.owned_leaf_cell_mask", m_owned_leaf_cell_mask);
    active("gravity_runtime.source_x", m_local_source_x);
    active("gravity_runtime.source_y", m_local_source_y);
    active("gravity_runtime.source_z", m_local_source_z);
    active("gravity_runtime.source_mass", m_local_source_mass);
    active("gravity_runtime.source_species", m_local_source_species_tag);
    active("gravity_runtime.source_softening", m_local_source_softening_comoving);
    active("gravity_runtime.source_softening_mask", m_local_source_softening_override_mask);
    active("gravity_runtime.source_particle_row", m_local_source_particle_row);
    active("gravity_runtime.source_cell_row", m_local_source_cell_row);
    active("gravity_runtime.target_particle_row", m_active_target_particle_row);
    active("gravity_runtime.target_cell_row", m_active_target_cell_row);
    active("gravity_runtime.target_species", m_active_target_species_tag);
    active("gravity_runtime.target_softening", m_active_target_softening_comoving);
    active("gravity_runtime.target_softening_mask", m_active_target_softening_override_mask);
    active("gravity_runtime.source_high_res_mask", m_source_is_high_res);
    active("gravity_runtime.target_high_res_mask", m_active_is_high_res);
    active("gravity_runtime.local_active_source_index", m_local_active_indices);
    active("gravity_runtime.local_active_particle_row", m_local_active_global_indices);
    active("gravity_runtime.force_refresh_particle_rows", m_force_refresh_particle_indices);
    active("gravity_runtime.authoritative_top_domain_leaves", m_authoritative_top_domain_leaves);
    persistent("gravity_runtime.particle_force_cache_valid", m_particle_force_cache_valid);
    persistent("gravity_runtime.cell_force_cache_valid", m_cell_force_cache_valid);

    core::MemoryReport owned_report = std::move(builder).finish();
    core::MemoryReport coordinator_report = m_tree_pm_coordinator.memoryReport();
    const std::array<core::MemoryReport, 2> reports{
        std::move(owned_report), std::move(coordinator_report)};
    core::MemoryReport merged = core::mergeMemoryReports(reports);
    if (m_has_pre_run_memory_estimate) {
      merged.notes.push_back(
          "gravity_pre_run_known_peak_bytes=" +
          std::to_string(m_last_pre_run_memory_estimate.known_peak_bytes));
      merged.notes.push_back(
          "gravity_pre_run_external_backend_unknown_bytes=" +
          std::to_string(m_last_pre_run_memory_estimate.external_backend_unknown_bytes));
      merged.notes.push_back(
          "gravity_pre_run_budget_required_bytes=" +
          std::to_string(m_last_pre_run_memory_estimate.budget_required_bytes));
    }
    if (m_has_process_memory_estimate) {
      merged.notes.push_back(
          "dmo_process_pre_run_known_owned_peak_bytes=" +
          std::to_string(m_last_process_memory_estimate.known_owned_peak_bytes));
      merged.notes.push_back(
          "dmo_process_pre_run_external_unknown_reserve_bytes=" +
          std::to_string(m_last_process_memory_estimate.external_unknown_reserve_bytes));
      merged.notes.push_back(
          "dmo_process_pre_run_safety_margin_bytes=" +
          std::to_string(m_last_process_memory_estimate.safety_margin_bytes));
      merged.notes.push_back(
          "dmo_process_pre_run_budget_required_bytes=" +
          std::to_string(m_last_process_memory_estimate.budget_required_bytes));
      merged.notes.push_back(
          "dmo_process_pre_run_aggregate_required_bytes=" +
          std::to_string(m_last_process_memory_estimate.aggregate_required_bytes));
    }
    return merged;
  }
  [[nodiscard]] const parallel::DecompositionRuntimeMeasurements& lastRuntimeDecompositionMeasurements() const noexcept {
    return m_last_decomposition_measurements;
  }
  void shutdownMpiResources() override {
    m_tree_pm_coordinator.shutdownMpiResources();
  }
  [[nodiscard]] std::uint64_t decompositionEpoch() const noexcept { return m_decomposition_epoch.value; }

  void restoreDecompositionEpoch(std::uint64_t decomposition_epoch) {
    if (m_has_long_range_field || m_last_committed_field_version != 0U) {
      throw std::logic_error("TreePM decomposition epoch must be restored before solver use");
    }
    m_decomposition_epoch = gravity::DecompositionEpoch{decomposition_epoch};
  }

  void commitParticleDecompositionChange() {
    m_decomposition_epoch = gravity::nextGravityIdentity(
        m_decomposition_epoch, "TreePM particle-decomposition epoch overflow");
    // Acceleration lanes are dense-row mirrors and are not part of the
    // migration payload. Invalidate immediately so an output boundary after
    // migration serializes an honest cache.valid=false state instead of
    // rejecting or persisting stale row ownership.
    m_force_cache_valid = false;
    m_force_cache_particle_index_generation = 0U;
    m_force_cache_cell_index_generation = 0U;
    m_force_cache_gas_identity_generation = 0U;
    m_force_cache_source_generation = 0U;
    m_particle_force_cache_valid.clear();
    m_cell_force_cache_valid.clear();
  }

  void installAuthoritativeTopDomainLeaves(
      std::span<const parallel::TopDomainLeaf> leaves) override {
    for (const parallel::TopDomainLeaf& leaf : leaves) {
      if (leaf.owner_rank != m_services.mpi_context.worldRank()) {
        throw std::invalid_argument(
            "gravity authoritative top-domain install received a leaf owned by another rank");
      }
      if (leaf.decomposition_epoch != m_decomposition_epoch.value) {
        throw std::invalid_argument(
            "gravity authoritative top-domain install received a stale decomposition epoch");
      }
    }
    m_authoritative_top_domain_leaves.assign(leaves.begin(), leaves.end());
    m_tree_pm_options.authoritative_domain_leaves =
        std::span<const parallel::TopDomainLeaf>(
            m_authoritative_top_domain_leaves.data(),
            m_authoritative_top_domain_leaves.size());
  }

  [[nodiscard]] std::span<const double> cellAccelX() const noexcept { return m_cell_accel_x; }
  [[nodiscard]] std::span<const double> cellAccelY() const noexcept { return m_cell_accel_y; }
  [[nodiscard]] std::span<const double> cellAccelZ() const noexcept { return m_cell_accel_z; }
  [[nodiscard]] std::span<const double> particleAccelX() const noexcept { return m_particle_accel_x; }
  [[nodiscard]] std::span<const double> particleAccelY() const noexcept { return m_particle_accel_y; }
  [[nodiscard]] std::span<const double> particleAccelZ() const noexcept { return m_particle_accel_z; }

  [[nodiscard]] io::GravityForceCachePersistentState exportRestartForceCache(
      const core::SimulationState& state) const {
    io::GravityForceCachePersistentState cache;
    cache.valid = m_force_cache_valid;
    if (!cache.valid) {
      return cache;
    }
    if (m_force_cache_particle_index_generation != state.particleIndexGeneration() ||
        m_force_cache_cell_index_generation != state.cellIndexGeneration() ||
        m_force_cache_gas_identity_generation != state.gasCellIdentityGeneration() ||
        m_force_cache_source_generation != state.gravitySourceGeneration()) {
      throw std::runtime_error(
          "gravity force cache cannot be checkpointed after particle/cell ownership or gas identity changed without refresh");
    }
    if (m_particle_accel_x.size() != state.particles.size() ||
        m_particle_accel_y.size() != state.particles.size() ||
        m_particle_accel_z.size() != state.particles.size() ||
        m_cell_accel_x.size() != state.cells.size() ||
        m_cell_accel_y.size() != state.cells.size() ||
        m_cell_accel_z.size() != state.cells.size()) {
      throw std::runtime_error(
          "gravity force cache cannot be checkpointed because its dense lanes do not match SimulationState");
    }
    cache.particle_id.assign(state.particle_sidecar.particle_id.begin(), state.particle_sidecar.particle_id.end());
    cache.gas_cell_id.assign(state.gas_cells.gas_cell_id.begin(), state.gas_cells.gas_cell_id.end());
    cache.particle_accel_x_comoving = m_particle_accel_x;
    cache.particle_accel_y_comoving = m_particle_accel_y;
    cache.particle_accel_z_comoving = m_particle_accel_z;
    cache.cell_accel_x_comoving = m_cell_accel_x;
    cache.cell_accel_y_comoving = m_cell_accel_y;
    cache.cell_accel_z_comoving = m_cell_accel_z;
    return cache;
  }

  void importRestartForceCache(
      const io::GravityForceCachePersistentState& cache,
      const core::SimulationState& state) {
    if (!cache.valid) {
      m_particle_accel_x.clear();
      m_particle_accel_y.clear();
      m_particle_accel_z.clear();
      m_cell_accel_x.clear();
      m_cell_accel_y.clear();
      m_cell_accel_z.clear();
      m_particle_force_cache_valid.clear();
      m_cell_force_cache_valid.clear();
      m_force_cache_valid = false;
      m_force_cache_particle_index_generation = 0U;
      m_force_cache_cell_index_generation = 0U;
      m_force_cache_gas_identity_generation = 0U;
      m_force_cache_source_generation = 0U;
      return;
    }
    const auto requireTriplet = [](const std::vector<double>& x,
                                   const std::vector<double>& y,
                                   const std::vector<double>& z,
                                   std::size_t expected,
                                   std::string_view label) {
      if (x.size() != expected || y.size() != expected || z.size() != expected) {
        throw std::runtime_error(
            "restart gravity force cache " + std::string(label) +
            " extent does not match authoritative SimulationState");
      }
      for (const double value : x) {
        if (!std::isfinite(value)) {
          throw std::runtime_error("restart gravity force cache contains non-finite values");
        }
      }
      for (const double value : y) {
        if (!std::isfinite(value)) {
          throw std::runtime_error("restart gravity force cache contains non-finite values");
        }
      }
      for (const double value : z) {
        if (!std::isfinite(value)) {
          throw std::runtime_error("restart gravity force cache contains non-finite values");
        }
      }
    };
    const auto buildIndexById = [](std::span<const std::uint64_t> ids,
                                   std::size_t expected,
                                   std::string_view label) {
      if (ids.size() != expected) {
        throw std::runtime_error("restart gravity force cache " + std::string(label) +
                                 " identity extent does not match acceleration lanes");
      }
      std::unordered_map<std::uint64_t, std::size_t> index_by_id;
      index_by_id.reserve(ids.size());
      for (std::size_t index = 0; index < ids.size(); ++index) {
        if (ids[index] == 0U || !index_by_id.emplace(ids[index], index).second) {
          throw std::runtime_error("restart gravity force cache " + std::string(label) +
                                   " identity lane contains a zero or duplicate stable ID");
        }
      }
      return index_by_id;
    };
    requireTriplet(
        cache.particle_accel_x_comoving, cache.particle_accel_y_comoving,
        cache.particle_accel_z_comoving, state.particles.size(), "particle");
    requireTriplet(
        cache.cell_accel_x_comoving, cache.cell_accel_y_comoving,
        cache.cell_accel_z_comoving, state.cells.size(), "gas-cell");
    const auto particle_index_by_id = buildIndexById(
        cache.particle_id, state.particles.size(), "particle");
    const auto gas_index_by_id = buildIndexById(
        cache.gas_cell_id, state.cells.size(), "gas-cell");

    m_particle_accel_x.resize(state.particles.size());
    m_particle_accel_y.resize(state.particles.size());
    m_particle_accel_z.resize(state.particles.size());
    for (std::size_t row = 0; row < state.particles.size(); ++row) {
      const auto it = particle_index_by_id.find(state.particle_sidecar.particle_id[row]);
      if (it == particle_index_by_id.end()) {
        throw std::runtime_error("restart gravity force cache is missing a persistent particle ID");
      }
      const std::size_t cached_row = it->second;
      m_particle_accel_x[row] = cache.particle_accel_x_comoving[cached_row];
      m_particle_accel_y[row] = cache.particle_accel_y_comoving[cached_row];
      m_particle_accel_z[row] = cache.particle_accel_z_comoving[cached_row];
    }
    m_cell_accel_x.resize(state.cells.size());
    m_cell_accel_y.resize(state.cells.size());
    m_cell_accel_z.resize(state.cells.size());
    for (std::size_t row = 0; row < state.cells.size(); ++row) {
      const auto it = gas_index_by_id.find(state.gas_cells.gas_cell_id[row]);
      if (it == gas_index_by_id.end()) {
        throw std::runtime_error("restart gravity force cache is missing a persistent gas-cell ID");
      }
      const std::size_t cached_row = it->second;
      m_cell_accel_x[row] = cache.cell_accel_x_comoving[cached_row];
      m_cell_accel_y[row] = cache.cell_accel_y_comoving[cached_row];
      m_cell_accel_z[row] = cache.cell_accel_z_comoving[cached_row];
    }
    m_force_cache_valid = true;
    m_force_cache_particle_index_generation = state.particleIndexGeneration();
    m_force_cache_cell_index_generation = state.cellIndexGeneration();
    m_force_cache_gas_identity_generation = state.gasCellIdentityGeneration();
    m_force_cache_source_generation = state.gravitySourceGeneration();
    m_particle_force_cache_valid.assign(state.particles.size(), 1U);
    m_cell_force_cache_valid.assign(state.cells.size(), 1U);
  }

  [[nodiscard]] static PmSyncSurface toPmSyncSurface(core::IntegrationStage stage) {
    if (stage == core::IntegrationStage::kGravityKickPre) {
      return PmSyncSurface::kKickPre;
    }
    if (stage == core::IntegrationStage::kForceRefresh) {
      return PmSyncSurface::kForceRefresh;
    }
    if (stage == core::IntegrationStage::kGravityKickPost) {
      return PmSyncSurface::kKickPost;
    }
    throw std::invalid_argument("TreePM sync surface is only defined on legal gravity force-refresh/kick stages");
  }

  [[nodiscard]] static std::string_view pmSyncSurfaceName(PmSyncSurface surface) {
    switch (surface) {
      case PmSyncSurface::kKickPre:
        return "kick_pre";
      case PmSyncSurface::kForceRefresh:
        return "force_refresh";
      case PmSyncSurface::kKickPost:
        return "kick_post";
    }
    return "unknown";
  }
  [[nodiscard]] static std::string_view pmRefreshReasonName(core::PmRefreshDirective::Reason reason) {
    switch (reason) {
      case core::PmRefreshDirective::Reason::kNone:
        return "none";
      case core::PmRefreshDirective::Reason::kInitialForceBootstrap:
        return "initial_force_bootstrap";
      case core::PmRefreshDirective::Reason::kScheduledForceRefreshStage:
        return "scheduled_force_refresh_stage";
      case core::PmRefreshDirective::Reason::kSourceMutationForceRefresh:
        return "source_mutation_force_refresh";
    }
    return "unknown";
  }

  void execute(GravityStageView& view) override {
    view.requireFresh();
    core::StepContext& context = internal::RuntimeStageAccess::gravityContext(
        view,
        {{RuntimeResourceKey::kParticlePosition, RuntimeResourceAccessMode::kRead},
         {RuntimeResourceKey::kParticleVelocity, RuntimeResourceAccessMode::kReadWrite},
         {RuntimeResourceKey::kParticleGravitySource, RuntimeResourceAccessMode::kRead},
         {RuntimeResourceKey::kHydroConservedState, RuntimeResourceAccessMode::kReadWrite},
         {RuntimeResourceKey::kHydroPrimitiveState, RuntimeResourceAccessMode::kReadWrite},
         {RuntimeResourceKey::kMigrationOwnership, RuntimeResourceAccessMode::kRead},
         {RuntimeResourceKey::kSchedulerTruth, RuntimeResourceAccessMode::kRead},
         {RuntimeResourceKey::kGravityAcceleration, RuntimeResourceAccessMode::kWrite},
         {RuntimeResourceKey::kIntegratorTruth, RuntimeResourceAccessMode::kReadWrite}});
    const bool is_kick_stage = context.stage == core::IntegrationStage::kGravityKickPre ||
        context.stage == core::IntegrationStage::kGravityKickPost;
    const parallel::MpiContext& mpi_context = m_services.mpi_context;
    const std::uint64_t world_size = static_cast<std::uint64_t>(mpi_context.worldSize());
    const std::size_t particle_count = context.state.particles.size();
    const std::size_t cell_count = context.state.cells.size();
    const bool particle_cache_generation_changed =
        m_force_cache_particle_index_generation != context.state.particleIndexGeneration() ||
        m_particle_force_cache_valid.size() != particle_count;
    const bool cell_cache_generation_changed =
        m_force_cache_cell_index_generation != context.state.cellIndexGeneration() ||
        m_force_cache_gas_identity_generation != context.state.gasCellIdentityGeneration() ||
        m_cell_force_cache_valid.size() != cell_count;
    const bool source_cache_generation_changed =
        m_force_cache_source_generation != context.state.gravitySourceGeneration();

    if (particle_cache_generation_changed) {
      // Dense-row identity is not a migration payload. A row is compatible only
      // after its authoritative source was refreshed in the current generation.
      m_particle_force_cache_valid.assign(particle_count, 0U);
      m_force_cache_particle_index_generation =
          context.state.particleIndexGeneration();
      m_force_cache_valid = false;
    }
    if (cell_cache_generation_changed) {
      // AMR refinement/coarsening and gas-cell migration change both dense row
      // indices and stable gas ownership. Invalidate together so a post-AMR
      // kick cannot reuse accelerations tied to a superseded leaf population.
      m_cell_force_cache_valid.assign(cell_count, 0U);
      m_force_cache_cell_index_generation = context.state.cellIndexGeneration();
      m_force_cache_gas_identity_generation = context.state.gasCellIdentityGeneration();
      m_force_cache_valid = false;
    }
    if (source_cache_generation_changed) {
      m_particle_force_cache_valid.assign(particle_count, 0U);
      m_cell_force_cache_valid.assign(cell_count, 0U);
      m_force_cache_valid = false;
    }
    const bool local_force_cache_incompatible = !m_force_cache_valid ||
        m_force_cache_particle_index_generation != context.state.particleIndexGeneration() ||
        m_force_cache_cell_index_generation != context.state.cellIndexGeneration() ||
        m_force_cache_gas_identity_generation != context.state.gasCellIdentityGeneration() ||
        m_force_cache_source_generation != context.state.gravitySourceGeneration();
    const std::uint64_t incompatible_cache_rank_count =
        mpi_context.allreduceSumUint64(local_force_cache_incompatible ? 1ULL : 0ULL);
    const bool needs_force_cache_rebuild = incompatible_cache_rank_count > 0U;
    const bool needs_initial_force_bootstrap =
        context.stage == core::IntegrationStage::kGravityKickPre &&
        !context.integrator_state.pm_long_range_field_valid;
    if (needs_initial_force_bootstrap && !context.pm_refresh_directive.initial_cache_bootstrap_allowed) {
      throw std::runtime_error(
          "TreePM initial force bootstrap requested outside an integrator-authorized global boundary");
    }
    const bool is_force_refresh_stage = context.pm_refresh_directive.force_refresh_surface ||
        needs_initial_force_bootstrap || (is_kick_stage && needs_force_cache_rebuild);
    if (!is_kick_stage && !is_force_refresh_stage) {
      throw std::logic_error("gravity handler received an unregistered stage");
    }
    if (context.pm_refresh_directive.has_sync_event &&
        context.pm_refresh_directive.refresh_long_range_field &&
        !context.boundary.pm_refresh_allowed) {
      throw std::runtime_error("TreePM long-range PM refresh reached an illegal integration boundary");
    }
    if (context.stage == core::IntegrationStage::kForceRefresh && !context.pm_refresh_directive.force_refresh_surface) {
      throw std::runtime_error("TreePM force-refresh stage lacks an integrator-issued PM refresh directive");
    }
    if (is_kick_stage && !is_force_refresh_stage) {
      applyCachedKick(context);
      return;
    }

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
    requireKickConsensus(context.state.gravitySourceGeneration(), "gravity_source_generation");

    const internal::SolverGhostRefreshReport gravity_ghost_refresh =
        internal::refreshParticleGhostsForSolver(
        context, mpi_context, "gravity.treepm", &m_ghost_cache_lifecycle);

    std::span<const std::uint32_t> force_target_particles =
        context.active_set.particle_indices;
    SourcePredictionEpoch source_prediction_epoch =
        context.pm_refresh_directive.requires_predicted_inactive_sources
        ? SourcePredictionEpoch::kStepEnd
        : SourcePredictionEpoch::kNone;
    if (is_kick_stage && needs_force_cache_rebuild) {
      // A dense-row ownership transition invalidates acceleration rows by
      // identity, not merely by extent. Rebuild every locally owned row in one
      // coordinated solve so later local pre-kicks never mix stale cache rows
      // with current ones. Sources whose persistent coordinates lag the
      // current scheduler time are predicted to the pre-kick epoch.
      m_force_refresh_particle_indices.clear();
      m_force_refresh_particle_indices.reserve(particle_count);
      const std::uint32_t local_rank =
          static_cast<std::uint32_t>(std::max(mpi_context.worldRank(), 0));
      for (std::size_t particle_index = 0; particle_index < particle_count;
           ++particle_index) {
        if (context.state.particle_sidecar.owning_rank[particle_index] ==
            local_rank) {
          m_force_refresh_particle_indices.push_back(
              static_cast<std::uint32_t>(particle_index));
        }
      }
      force_target_particles = m_force_refresh_particle_indices;
      source_prediction_epoch = SourcePredictionEpoch::kStepBegin;
    }

    // Budget against a conservative per-rank upper bound before materializing
    // O(N) gravity staging/tree/communication buffers. Canonical state is
    // deliberately outside this gravity-only budget and is reported by core.
    const std::size_t conservative_local_source_count = core::checkedSizeAdd(
        particle_count, cell_count, "gravity pre-run local source count");
    const std::size_t conservative_local_target_count = core::checkedSizeAdd(
        force_target_particles.size(), cell_count, "gravity pre-run local target count");
    m_last_pre_run_memory_estimate = gravity::estimateGravityMemory(
        gravity::GravityMemoryEstimateInput{
            .local_source_count = static_cast<std::uint64_t>(conservative_local_source_count),
            .local_target_count = static_cast<std::uint64_t>(conservative_local_target_count),
            .local_particle_count = static_cast<std::uint64_t>(particle_count),
            .local_cell_count = static_cast<std::uint64_t>(cell_count),
            .tree_leaf_size = m_tree_pm_options.tree_options.max_leaf_size,
            .multipole_order = m_tree_pm_options.tree_options.multipole_order,
            .pm_shape = m_pm_grid_shape,
            .assignment_scheme = m_tree_pm_options.pm_options.assignment_scheme,
            .decomposition_mode = m_tree_pm_options.pm_options.decomposition_mode,
            .mpi_rank_count = static_cast<std::uint32_t>(std::max(m_runtime_topology.world_size, 1)),
            .mpi_world_rank = std::max(mpi_context.worldRank(), 0),
            .zoom_enabled = m_tree_pm_options.enable_zoom_long_range_correction,
            .zoom_pm_shape = m_tree_pm_options.zoom_focused_pm_shape,
            .periodic_tree_coordinates =
                m_tree_pm_options.pm_options.boundary_condition == gravity::PmBoundaryCondition::kPeriodic,
            .indexed_target_coordinates = true,
            .cuda_resident = m_runtime_topology.usesCuda(),
            .tree_exchange_batch_bytes = m_tree_pm_options.tree_exchange_batch_bytes,
            .pm_exchange_batch_bytes = m_tree_pm_options.pm_options.routing_exchange_batch_bytes,
            .backend_unknown_reserve_bytes =
                m_config.parallel.gravity_backend_unknown_reserve_bytes,
            .safety_margin_fraction =
                m_config.parallel.gravity_memory_safety_margin_fraction,
        });
    m_has_pre_run_memory_estimate = true;
    gravity::enforceGravityMemoryBudget(
        m_last_pre_run_memory_estimate,
        m_config.parallel.gravity_memory_budget_bytes);

    if (m_config.parallel.process_memory_budget_bytes != 0U &&
        context.particle_scheduler == nullptr) {
      throw std::runtime_error(
          "DMO process memory budget is enabled but the gravity stage has no authoritative particle scheduler for owned-memory preflight");
    }
    const std::uint64_t scheduler_current_size_bytes = context.particle_scheduler != nullptr
        ? context.particle_scheduler->logicalSizeBytes()
        : 0U;
    const std::uint64_t scheduler_owned_bytes = context.particle_scheduler != nullptr
        ? context.particle_scheduler->ownedCapacityBytes()
        : 0U;
    const std::uint64_t scheduler_high_water_bytes = context.particle_scheduler != nullptr
        ? context.particle_scheduler->ownedCapacityHighWaterBytes()
        : 0U;
    const core::MemoryReport canonical_runtime_report =
        core::collectSimulationMemoryReport(context.state, context.workspace);
    m_last_process_memory_estimate = gravity::estimateDmoProcessMemory(
        canonical_runtime_report,
        m_last_pre_run_memory_estimate,
        gravity::DmoProcessMemoryPolicy{
            .mpi_rank_count = static_cast<std::uint32_t>(
                std::max(m_runtime_topology.world_size, 1)),
            .scheduler_current_size_bytes = scheduler_current_size_bytes,
            .scheduler_owned_bytes = scheduler_owned_bytes,
            .scheduler_high_water_bytes = scheduler_high_water_bytes,
            .output_restart_overlap_bytes =
                m_config.parallel.process_output_restart_overlap_bytes,
            .mpi_external_reserve_bytes =
                m_config.parallel.process_mpi_unknown_reserve_bytes,
            .fftw_external_reserve_bytes =
                m_config.parallel.process_fftw_unknown_reserve_bytes,
            .hdf5_external_reserve_bytes =
                m_config.parallel.process_hdf5_unknown_reserve_bytes,
            .allocator_external_reserve_bytes =
                m_config.parallel.process_allocator_unknown_reserve_bytes,
            .safety_margin_fraction =
                m_config.parallel.process_memory_safety_margin_fraction,
        });

    // Replace the estimate helper's equal-rank aggregate with the exact
    // configured-topology reduction.  Particle decomposition can be uneven,
    // so production diagnostics must report the sum/max of the rank-local
    // preflights rather than local_required * world_size.  All ranks enter
    // these reductions before any rank is allowed to reject its process
    // budget, preventing a low-memory rank from throwing while peers continue
    // into later collectives.
    const std::uint64_t global_modeled_subtotal_bytes =
        mpi_context.allreduceSumUint64(
            m_last_process_memory_estimate.modeled_subtotal_bytes);
    const std::uint64_t rank_max_modeled_subtotal_bytes =
        mpi_context.allreduceMaxUint64(
            m_last_process_memory_estimate.modeled_subtotal_bytes);
    const std::uint64_t global_required_bytes =
        mpi_context.allreduceSumUint64(
            m_last_process_memory_estimate.budget_required_bytes);
    const std::uint64_t rank_max_required_bytes =
        mpi_context.allreduceMaxUint64(
            m_last_process_memory_estimate.budget_required_bytes);

    m_last_process_memory_estimate.aggregate_required_bytes =
        global_required_bytes;
    m_last_process_memory_estimate.report.notes.push_back(
        "dmo_process_runtime_exact_aggregate_required_bytes=" +
        std::to_string(global_required_bytes));
    m_last_process_memory_estimate.report.notes.push_back(
        "dmo_process_runtime_rank_max_required_bytes=" +
        std::to_string(rank_max_required_bytes));
    auto& distributed = m_last_process_memory_estimate.report.distributed;
    distributed.valid = true;
    distributed.rank_count = std::max(m_runtime_topology.world_size, 1);
    distributed.local_owned_bytes =
        m_last_process_memory_estimate.modeled_subtotal_bytes;
    distributed.global_sum_owned_bytes = global_modeled_subtotal_bytes;
    distributed.rank_max_owned_bytes = rank_max_modeled_subtotal_bytes;
    if (global_modeled_subtotal_bytes != 0U) {
      const long double mean_bytes =
          static_cast<long double>(global_modeled_subtotal_bytes) /
          static_cast<long double>(distributed.rank_count);
      distributed.max_to_mean_imbalance_ratio = static_cast<double>(
          static_cast<long double>(rank_max_modeled_subtotal_bytes) /
          mean_bytes);
    }
    m_has_process_memory_estimate = true;

    if (context.profiler_session != nullptr) {
      context.profiler_session->recordEvent(core::RuntimeEvent{
          .event_kind = "memory.dmo_process_preflight",
          .severity = core::RuntimeEventSeverity::kInfo,
          .subsystem = "core.memory",
          .step_index = context.integrator_state.step_index,
          .simulation_time_code = context.integrator_state.current_time_code,
          .scale_factor = context.integrator_state.current_scale_factor,
          .message = "authoritative DMO per-rank process memory preflight",
          .payload = {
              {"known_owned_peak_bytes",
               std::to_string(m_last_process_memory_estimate.known_owned_peak_bytes)},
              {"external_unknown_reserve_bytes",
               std::to_string(m_last_process_memory_estimate.external_unknown_reserve_bytes)},
              {"modeled_subtotal_bytes",
               std::to_string(m_last_process_memory_estimate.modeled_subtotal_bytes)},
              {"safety_margin_bytes",
               std::to_string(m_last_process_memory_estimate.safety_margin_bytes)},
              {"required_bytes",
               std::to_string(m_last_process_memory_estimate.budget_required_bytes)},
              {"configured_process_budget_bytes",
               std::to_string(m_config.parallel.process_memory_budget_bytes)},
              {"rank_max_required_bytes", std::to_string(rank_max_required_bytes)},
              {"aggregate_required_bytes", std::to_string(global_required_bytes)},
              {"rank_count", std::to_string(distributed.rank_count)},
          },
      });
    }

    const std::uint64_t process_budget_bytes =
        m_config.parallel.process_memory_budget_bytes;
    const bool local_process_budget_failure =
        process_budget_bytes != 0U &&
        m_last_process_memory_estimate.budget_required_bytes >
            process_budget_bytes;
    const std::uint64_t failed_process_budget_rank_count =
        mpi_context.allreduceSumUint64(
            local_process_budget_failure ? 1ULL : 0ULL);
    if (failed_process_budget_rank_count != 0U) {
      if (local_process_budget_failure) {
        gravity::enforceDmoProcessMemoryBudget(
            m_last_process_memory_estimate, process_budget_bytes);
      }
      throw std::runtime_error(
          "DMO process memory preflight failed on " +
          std::to_string(failed_process_budget_rank_count) +
          " rank(s); local_required_bytes=" +
          std::to_string(m_last_process_memory_estimate.budget_required_bytes) +
          ", configured_process_budget_bytes=" +
          std::to_string(process_budget_bytes) +
          ", rank_max_required_bytes=" +
          std::to_string(rank_max_required_bytes) +
          ", aggregate_required_bytes=" +
          std::to_string(global_required_bytes));
    }

    rebuildOwnedParticleCompactView(
        context,
        force_target_particles,
        source_prediction_epoch);
    m_local_kick_particle_count = 0U;
    for (const std::uint32_t particle_index :
         context.active_set.particle_indices) {
      if (particle_index < m_owned_local_index_by_global.size() &&
          m_owned_local_index_by_global[particle_index] >= 0) {
        ++m_local_kick_particle_count;
      }
    }
    m_active_accel_x.assign(m_local_active_indices.size(), 0.0);
    m_active_accel_y.assign(m_local_active_indices.size(), 0.0);
    m_active_accel_z.assign(m_local_active_indices.size(), 0.0);
    m_active_is_high_res.assign(m_local_active_indices.size(), 0U);
    m_active_slot_by_particle.assign(particle_count, -1);
    m_active_slot_by_cell.assign(cell_count, -1);
    constexpr std::uint32_t no_target_row =
        std::numeric_limits<std::uint32_t>::max();
    for (std::size_t target_slot = 0; target_slot < m_local_active_indices.size(); ++target_slot) {
      if (m_active_target_particle_row[target_slot] != no_target_row) {
        m_active_slot_by_particle[m_active_target_particle_row[target_slot]] =
            static_cast<int>(target_slot);
      }
      if (m_active_target_cell_row[target_slot] != no_target_row) {
        m_active_slot_by_cell[m_active_target_cell_row[target_slot]] =
            static_cast<int>(target_slot);
      }
    }
    ensureZoomMembershipLoaded(context.state);
    const double box_size_x = m_config.cosmology.box_size_x_mpc_comoving;
    const double box_size_y = m_config.cosmology.box_size_y_mpc_comoving;
    const double box_size_z = m_config.cosmology.box_size_z_mpc_comoving;
    m_source_is_high_res.assign(m_local_source_x.size(), 0U);
    for (std::size_t i = 0; i < m_local_source_x.size(); ++i) {
      if (!m_zoom_high_res_particle_ids.empty()) {
        const std::uint32_t particle_row = m_local_source_particle_row[i];
        if (particle_row != no_target_row) {
          const std::uint64_t particle_id =
              context.state.particle_sidecar.particle_id[particle_row];
          m_source_is_high_res[i] =
              m_zoom_high_res_particle_ids.contains(particle_id) ? 1U : 0U;
          continue;
        }
        const std::uint32_t cell_row = m_local_source_cell_row[i];
        if (cell_row != no_target_row) {
          const auto parent_id =
              context.state.parentParticleIdForGasCellRow(cell_row);
          if (parent_id.has_value()) {
            m_source_is_high_res[i] =
                m_zoom_high_res_particle_ids.contains(*parent_id) ? 1U : 0U;
            continue;
          }
        }
        // A parentless refined cell has no particle lineage authority. Fall
        // back to the explicit spatial zoom region rather than inventing one.
      }
      const double dx = m_local_source_x[i] - m_tree_pm_options.zoom_region_center_x_comoving;
      const double dy = m_local_source_y[i] - m_tree_pm_options.zoom_region_center_y_comoving;
      const double dz = m_local_source_z[i] - m_tree_pm_options.zoom_region_center_z_comoving;
      const double wrapped_dx = dx - box_size_x * std::nearbyint(dx / box_size_x);
      const double wrapped_dy = dy - box_size_y * std::nearbyint(dy / box_size_y);
      const double wrapped_dz = dz - box_size_z * std::nearbyint(dz / box_size_z);
      const double r = std::sqrt(wrapped_dx * wrapped_dx + wrapped_dy * wrapped_dy + wrapped_dz * wrapped_dz);
      m_source_is_high_res[i] = (r <= m_tree_pm_options.zoom_region_radius_comoving) ? 1U : 0U;
    }
    for (std::size_t i = 0; i < m_local_active_indices.size(); ++i) {
      m_active_is_high_res[i] = m_source_is_high_res[m_local_active_indices[i]];
    }
    m_tree_pm_options.source_is_high_res = m_source_is_high_res;
    m_tree_pm_options.active_is_high_res = m_active_is_high_res;

    m_active_previous_acceleration_magnitude.clear();
    const bool relative_mac_cache_compatible =
        m_tree_pm_options.tree_options.opening_criterion == gravity::TreeOpeningCriterion::kRelativeForceError &&
        m_force_cache_valid &&
        m_force_cache_particle_index_generation == context.state.particleIndexGeneration() &&
        m_force_cache_cell_index_generation == context.state.cellIndexGeneration() &&
        m_force_cache_gas_identity_generation == context.state.gasCellIdentityGeneration() &&
        m_particle_accel_x.size() == particle_count &&
        m_particle_accel_y.size() == particle_count &&
        m_particle_accel_z.size() == particle_count &&
        m_cell_accel_x.size() == cell_count &&
        m_cell_accel_y.size() == cell_count &&
        m_cell_accel_z.size() == cell_count;
    if (relative_mac_cache_compatible) {
      m_active_previous_acceleration_magnitude.resize(
          m_local_active_indices.size(),
          std::numeric_limits<double>::quiet_NaN());
      for (std::size_t target_slot = 0; target_slot < m_local_active_indices.size(); ++target_slot) {
        const std::uint32_t particle_row = m_active_target_particle_row[target_slot];
        if (particle_row != no_target_row &&
            particle_row < m_particle_force_cache_valid.size() &&
            m_particle_force_cache_valid[particle_row] != 0U) {
          m_active_previous_acceleration_magnitude[target_slot] = std::sqrt(
              m_particle_accel_x[particle_row] * m_particle_accel_x[particle_row] +
              m_particle_accel_y[particle_row] * m_particle_accel_y[particle_row] +
              m_particle_accel_z[particle_row] * m_particle_accel_z[particle_row]);
          continue;
        }
        const std::uint32_t cell_row = m_active_target_cell_row[target_slot];
        if (cell_row != no_target_row &&
            cell_row < m_cell_force_cache_valid.size() &&
            m_cell_force_cache_valid[cell_row] != 0U) {
          m_active_previous_acceleration_magnitude[target_slot] = std::sqrt(
              m_cell_accel_x[cell_row] * m_cell_accel_x[cell_row] +
              m_cell_accel_y[cell_row] * m_cell_accel_y[cell_row] +
              m_cell_accel_z[cell_row] * m_cell_accel_z[cell_row]);
        }
      }
    }

    gravity::TreePmForceAccumulatorView accumulator{
        .active_particle_index = m_local_active_indices,
        .accel_x_comoving = m_active_accel_x,
        .accel_y_comoving = m_active_accel_y,
        .accel_z_comoving = m_active_accel_z,
        .previous_acceleration_magnitude_code = m_active_previous_acceleration_magnitude,
    };

    // The integrator owns PM cadence state. At every rank-coordinated force
    // evaluation it issues a concrete sync event. The active particle set may
    // be local, but all ranks participate and inactive sources are predicted,
    // so production cadence one advances PM truth at that same force epoch.
    PmLongRangeKickDecision decision{
        .sync_surface = toPmSyncSurface(context.stage),
        .gravity_kick_opportunity = context.integrator_state.pm_sync_state.gravityKickOpportunity(),
        .refresh_long_range_field = false,
        .field_version = context.integrator_state.pm_sync_state.fieldVersion(),
        .last_refresh_opportunity = context.integrator_state.pm_sync_state.lastRefreshOpportunity(),
        .field_built_step_index = context.integrator_state.pm_sync_state.lastRefreshStepIndex(),
        .field_built_scale_factor = context.integrator_state.pm_sync_state.lastRefreshScaleFactor(),
    };

    if (context.pm_refresh_directive.has_sync_event) {
      if (!context.pm_refresh_directive.cadence_opportunity_allowed) {
        throw std::runtime_error(
            "TreePM received a PM sync event outside an integrator-authorized refresh opportunity");
      }
      const core::PmSyncEvent sync_event{
          .gravity_kick_opportunity = context.pm_refresh_directive.gravity_kick_opportunity,
          .refresh_long_range_field = context.pm_refresh_directive.refresh_long_range_field,
          .field_version = context.pm_refresh_directive.field_version,
          .last_refresh_opportunity = context.pm_refresh_directive.last_refresh_opportunity,
          .field_built_step_index = context.pm_refresh_directive.field_built_step_index,
          .field_built_scale_factor = context.pm_refresh_directive.field_built_scale_factor,
      };
      requireKickConsensus(sync_event.gravity_kick_opportunity, "gravity_kick_opportunity");
      const std::uint64_t refresh_vote = sync_event.refresh_long_range_field ? 1ULL : 0ULL;
      const std::uint64_t refresh_vote_sum = mpi_context.allreduceSumUint64(refresh_vote);
      if (refresh_vote_sum != 0ULL && refresh_vote_sum != world_size) {
        throw std::runtime_error(
            "TreePM long-range cadence decision diverged across ranks: refresh_vote_sum=" +
            std::to_string(refresh_vote_sum) + ", world_size=" + std::to_string(world_size));
      }
      const bool refresh_long_range = (refresh_vote_sum == world_size);
      if (refresh_long_range != sync_event.refresh_long_range_field) {
        throw std::runtime_error(
            "TreePM long-range cadence local decision drifted from rank consensus");
      }

      decision.gravity_kick_opportunity = sync_event.gravity_kick_opportunity;
      decision.refresh_long_range_field = refresh_long_range;
      decision.field_version = sync_event.field_version;
      decision.last_refresh_opportunity = sync_event.last_refresh_opportunity;
      decision.field_built_step_index = sync_event.field_built_step_index;
      decision.field_built_scale_factor = sync_event.field_built_scale_factor;
    } else {
      if (!context.integrator_state.pm_long_range_field_valid || decision.field_version == 0U) {
        throw std::runtime_error(
            "TreePM local force refresh attempted before a valid long-range PM field existed");
      }
      if (context.pm_refresh_directive.cadence_opportunity_allowed) {
        throw std::runtime_error(
            "TreePM cadence opportunity was declared without an accompanying integrator-owned PM sync event");
      }
    }

    requireKickConsensus(decision.field_version, "long_range_field_version");
    requireKickConsensus(decision.last_refresh_opportunity, "last_long_range_refresh_opportunity");
    double force_evaluation_scale_factor = decision.field_built_scale_factor;
    if (decision.refresh_long_range_field) {
      force_evaluation_scale_factor =
          context.pm_refresh_directive.force_evaluation_scale_factor;
      if (!std::isfinite(force_evaluation_scale_factor) ||
          force_evaluation_scale_factor <= 0.0 ||
          force_evaluation_scale_factor != decision.field_built_scale_factor) {
        throw std::runtime_error(
            "TreePM refresh force-evaluation scale factor disagrees with integrator field-build metadata");
      }
    }
    if (!std::isfinite(force_evaluation_scale_factor) ||
        force_evaluation_scale_factor <= 0.0) {
      throw std::runtime_error(
          "TreePM force evaluation requires a finite positive scale factor");
    }
    // At a post-drift global refresh this is timeline_step.scale_factor_end,
    // even though IntegratorState.current_scale_factor remains at the step
    // beginning until commitStep(). A reuse retains the integrator-owned field
    // build scale as validity metadata; the returned Newtonian kernel itself is
    // scale-free and the KDK integrator owns the cosmological kick factor.
    m_tree_pm_options.pm_options.scale_factor =
        m_mode_policy.cosmological_comoving_frame
        ? force_evaluation_scale_factor
        : 1.0;
    // Dense-row generations are rank-local mirrors and can legitimately differ
    // after uneven compaction. The workflow advances this shared epoch exactly
    // once, on every rank, after a committed particle-ownership transition.
    m_tree_pm_options.decomposition_epoch = m_decomposition_epoch;
    // Gravity source identity is owned by authoritative SimulationState and is
    // intentionally independent from PM cadence sequencing. PM field version
    // answers "which mesh refresh?"; source generation answers "which source
    // state?"; force epoch answers "which evaluation opportunity?".
    m_tree_pm_options.source_generation = gravity::GravitySourceGeneration{
        context.state.gravitySourceGeneration()};
    m_tree_pm_options.pm_field_version = gravity::PmFieldVersion{decision.field_version};
    m_tree_pm_options.force_epoch = gravity::ForceEvaluationEpoch{
        .sequence = decision.gravity_kick_opportunity,
        .scale_factor = force_evaluation_scale_factor,
    };
    if (!decision.refresh_long_range_field &&
        context.integrator_state.pm_source_generation !=
            context.state.gravitySourceGeneration()) {
      throw std::runtime_error(
          "TreePM PM-field reuse rejected because authoritative gravity sources changed since the committed mesh refresh");
    }

    const gravity::GravitySourceSnapshot source_snapshot{
        .pos_x_comoving = m_local_source_x,
        .pos_y_comoving = m_local_source_y,
        .pos_z_comoving = m_local_source_z,
        .mass_code = m_local_source_mass,
        .species_tag = m_local_source_species_tag,
        .generation = m_tree_pm_options.source_generation,
        .evaluation_epoch = m_tree_pm_options.force_epoch,
        .contains_predicted_coordinates = m_source_predicted_inactive_count != 0U,
    };
    source_snapshot.validate();

    if (m_particle_accel_x.size() != particle_count) {
      m_particle_accel_x.assign(particle_count, 0.0);
      m_particle_accel_y.assign(particle_count, 0.0);
      m_particle_accel_z.assign(particle_count, 0.0);
      m_particle_force_cache_valid.assign(particle_count, 0U);
    }
    if (m_cell_accel_x.size() != context.state.cells.size()) {
      m_cell_accel_x.assign(context.state.cells.size(), 0.0);
      m_cell_accel_y.assign(context.state.cells.size(), 0.0);
      m_cell_accel_z.assign(context.state.cells.size(), 0.0);
      m_cell_force_cache_valid.assign(context.state.cells.size(), 0U);
    } else if (m_cell_force_cache_valid.size() != context.state.cells.size()) {
      m_cell_force_cache_valid.assign(context.state.cells.size(), 0U);
    }
    const gravity::TreeSofteningView softening_view{
        .source_species_tag = std::span<const std::uint32_t>(m_local_source_species_tag.data(), m_local_source_species_tag.size()),
        .source_particle_epsilon_comoving = std::span<const double>(
            m_local_source_softening_comoving.empty() ? nullptr : m_local_source_softening_comoving.data(),
            m_local_source_softening_comoving.size()),
        .source_particle_epsilon_override_mask = std::span<const std::uint8_t>(
            m_local_source_softening_override_mask.empty() ? nullptr : m_local_source_softening_override_mask.data(),
            m_local_source_softening_override_mask.size()),
        .target_species_tag = std::span<const std::uint32_t>(m_active_target_species_tag.data(), m_active_target_species_tag.size()),
        .target_particle_epsilon_comoving = std::span<const double>(
            m_active_target_softening_comoving.empty() ? nullptr : m_active_target_softening_comoving.data(),
            m_active_target_softening_comoving.size()),
        .target_particle_epsilon_override_mask = std::span<const std::uint8_t>(
            m_active_target_softening_override_mask.empty() ? nullptr : m_active_target_softening_override_mask.data(),
            m_active_target_softening_override_mask.size()),
        .species_policy = m_tree_pm_species_softening,
    };
    gravity::TreePmProfileEvent tree_pm_profile;
    m_tree_pm_coordinator.solveActiveSetWithPmCadence(
        source_snapshot.pos_x_comoving,
        source_snapshot.pos_y_comoving,
        source_snapshot.pos_z_comoving,
        source_snapshot.mass_code,
        accumulator,
        m_tree_pm_options,
        decision.refresh_long_range_field,
        &tree_pm_profile,
        &m_last_tree_pm_diagnostics,
        softening_view);
    const std::uint64_t expected_pm_solve_count =
        decision.refresh_long_range_field ? 1U : 0U;
    if (m_last_tree_pm_diagnostics.pm_solve_count != expected_pm_solve_count ||
        m_last_tree_pm_diagnostics.pm_reuse_count != 1U - expected_pm_solve_count) {
      throw std::runtime_error(
          "TreePM solver PM solve/reuse outcome disagrees with the integrator-owned cadence directive");
    }

    m_last_decomposition_measurements = parallel::DecompositionRuntimeMeasurements{
        .tree_pair_evaluations_recent = m_last_tree_pm_diagnostics.residual_pair_evaluations +
            tree_pm_profile.tree_profile.particle_particle_interactions,
        .tree_remote_request_bytes_recent = m_last_tree_pm_diagnostics.residual_remote_request_bytes +
            m_last_tree_pm_diagnostics.residual_remote_response_bytes,
        .pm_mesh_cells_touched_recent = static_cast<std::uint64_t>(m_tree_pm_coordinator.slabLayout().localCellCount()),
        .pm_fft_transpose_bytes_recent = tree_pm_profile.pm_profile.fft_transpose_bytes,
        .amr_patch_cells_updated_recent = static_cast<std::uint64_t>(context.state.cells.size()),
        .hydro_face_fluxes_recent = 0,
        .ghost_exchange_bytes_recent = gravity_ghost_refresh.sent_bytes + gravity_ghost_refresh.received_bytes,
        .tree_wall_ms_recent = tree_pm_profile.tree_profile.build_ms + tree_pm_profile.tree_profile.multipole_ms +
            tree_pm_profile.tree_profile.traversal_ms + tree_pm_profile.tree_short_range_ms,
        .pm_wall_ms_recent = tree_pm_profile.pm_profile.assign_ms + tree_pm_profile.pm_profile.fft_forward_ms +
            tree_pm_profile.pm_profile.poisson_ms + tree_pm_profile.pm_profile.gradient_ms +
            tree_pm_profile.pm_profile.fft_inverse_ms + tree_pm_profile.pm_profile.fft_transpose_ms +
            tree_pm_profile.pm_profile.interpolate_ms,
        .gpu_kernel_ms_recent = tree_pm_profile.pm_profile.device_kernel_ms,
        .accelerator_occupancy_fraction_recent = (tree_pm_profile.pm_profile.device_kernel_ms > 0.0) ? 1.0 : 0.0,
        .has_measurements = true,
    };

    const bool allow_heavy_reference_checks =
        m_config.analysis.diagnostics_execution_policy ==
        core::AnalysisConfig::DiagnosticsExecutionPolicy::kAllIncludingProvisional;
    GravityHealthSummary gravity_health = validateGravityHealth(
        context,
        decision,
        allow_heavy_reference_checks);

    if (decision.refresh_long_range_field) {
      m_has_long_range_field = true;
      ++m_long_range_refresh_count;
    } else {
      ++m_long_range_reuse_count;
    }
    context.pm_refresh_directive.solver_executed = true;
    m_last_committed_field_version = decision.field_version;

    const std::uint64_t inactive_particles_skipped = static_cast<std::uint64_t>(
        context.state.particles.size() - m_local_kick_particle_count);
    m_cadence_records.push_back(ReferenceWorkflowReport::TreePmCadenceRecord{
        .step_index = context.integrator_state.step_index,
        .stage_name = std::string(core::integrationStageName(context.stage)),
        .pm_sync_surface = std::string(pmSyncSurfaceName(decision.sync_surface)),
        .pm_refresh_reason = std::string(pmRefreshReasonName(context.pm_refresh_directive.reason)),
        .gravity_kick_opportunity = decision.gravity_kick_opportunity,
        .field_version = decision.field_version,
        .last_refresh_opportunity = decision.last_refresh_opportunity,
        .field_built_step_index = decision.field_built_step_index,
        .field_built_scale_factor = decision.field_built_scale_factor,
        .field_age_in_kick_opportunities =
            decision.gravity_kick_opportunity - decision.last_refresh_opportunity,
        .active_particles_kicked = static_cast<std::uint64_t>(m_local_kick_particle_count),
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
                      {"gravity_kick_opportunity", std::to_string(decision.gravity_kick_opportunity)},
                      {"field_version", std::to_string(decision.field_version)},
                      {"field_built_step_index", std::to_string(decision.field_built_step_index)},
                      {"field_built_scale_factor", formatRuntimeDouble(decision.field_built_scale_factor)},
                      {"pm_update_cadence_steps", std::to_string(context.integrator_state.pm_sync_state.cadenceSteps())},
                      {"pm_grid", std::to_string(m_pm_grid_shape.nx) + "x" + std::to_string(m_pm_grid_shape.ny) +
                              "x" + std::to_string(m_pm_grid_shape.nz)},
                      {"pm_assignment_scheme", m_pm_assignment_scheme},
                      {"pm_window_deconvolution", m_config.numerics.treepm_enable_window_deconvolution ? "true" : "false"},
                      {"asmth_cells", formatRuntimeDouble(m_config.numerics.treepm_asmth_cells)},
                      {"rcut_cells", formatRuntimeDouble(m_config.numerics.treepm_rcut_cells)},
                      {"tree_opening_criterion", treePmOpeningCriterionName(
                          m_config.numerics.treepm_tree_opening_criterion)},
                      {"tree_opening_theta", formatRuntimeDouble(
                          m_config.numerics.treepm_tree_opening_theta)},
                      {"tree_relative_force_tolerance", formatRuntimeDouble(
                          m_config.numerics.treepm_tree_relative_force_tolerance)},
                      {"tree_relative_force_acceleration_floor", formatRuntimeDouble(
                          m_config.numerics.treepm_tree_relative_force_acceleration_floor)},
                      {"relative_mac_previous_force_available",
                          relative_mac_cache_compatible ? "true" : "false"},
                      {"mesh_spacing_x_mpc_comoving", formatRuntimeDouble(m_mesh_spacing_x_mpc_comoving)},
                      {"mesh_spacing_y_mpc_comoving", formatRuntimeDouble(m_mesh_spacing_y_mpc_comoving)},
                      {"mesh_spacing_z_mpc_comoving", formatRuntimeDouble(m_mesh_spacing_z_mpc_comoving)},
                      {"split_scale_mpc_comoving", formatRuntimeDouble(m_tree_pm_options.split_policy.split_scale_comoving)},
                      {"cutoff_radius_mpc_comoving", formatRuntimeDouble(m_tree_pm_options.split_policy.cutoff_radius_comoving)},
                      {"softening_policy", describeSofteningPolicy(m_config, &context.state)},
                      {"softening_kernel", "plummer"},
                      {"softening_epsilon_kpc_comoving", formatRuntimeDouble(m_config.numerics.gravity_softening_kpc_comoving)},
                      {"pm_fft_backend", m_pm_backend},
                      {"pm_backend_capability", std::string(gravity::pmBackendCapabilityName(gravity::pmBackendCapability()))},
                      {"gravity_acceptance_profile_id", "unverified_current_source"},
                      {"active_particles_kicked", std::to_string(m_local_active_global_indices.size())},
                      {"inactive_particles_skipped", std::to_string(inactive_particles_skipped)},
                      {"ghost_refresh_sent_bytes", std::to_string(gravity_ghost_refresh.sent_bytes)},
                      {"ghost_refresh_received_bytes", std::to_string(gravity_ghost_refresh.received_bytes)},
                      {"ghost_refresh_committed_slots", std::to_string(gravity_ghost_refresh.committed_slots)},
                      {"predicted_inactive_source_particles", std::to_string(m_source_predicted_inactive_count)},
                      {"predicted_inactive_sources_required",
                          context.pm_refresh_directive.requires_predicted_inactive_sources ? "true" : "false"},
                      {"pm_refresh_reason", std::string(pmRefreshReasonName(context.pm_refresh_directive.reason))},
                      {"pm_refresh_force_eval_scale_factor", formatRuntimeDouble(context.pm_refresh_directive.force_evaluation_scale_factor)},
                      {"refreshed_long_range_field", decision.refresh_long_range_field ? "true" : "false"}},
      });
      const auto decomposition_descriptor = gravity::describePmDecomposition(
          m_pm_grid_shape,
          m_tree_pm_coordinator.slabLayout(),
          m_tree_pm_options.pm_options.decomposition_mode);
      context.profiler_session->recordEvent(core::RuntimeEvent{
          .event_kind = "gravity.treepm_let",
          .severity = core::RuntimeEventSeverity::kInfo,
          .subsystem = "gravity.treepm",
          .step_index = context.integrator_state.step_index,
          .simulation_time_code = context.integrator_state.current_time_code,
          .scale_factor = context.integrator_state.current_scale_factor,
          .message = "locality-driven TreePM short-range communication metrics",
          .payload = {
              {"pm_decomposition_architecture", std::string(
                  gravity::pmDecompositionTopologyName(decomposition_descriptor.topology))},
              {"top_level_domain_leaf_count", std::to_string(m_last_tree_pm_diagnostics.top_level_domain_leaf_count)},
              {"authoritative_domain_leaf_count", std::to_string(m_last_tree_pm_diagnostics.authoritative_domain_leaf_count)},
              {"domain_hierarchy_node_count", std::to_string(m_last_tree_pm_diagnostics.domain_hierarchy_node_count)},
              {"domain_cache_hits", std::to_string(m_last_tree_pm_diagnostics.domain_cache_hit_count)},
              {"domain_cache_misses", std::to_string(m_last_tree_pm_diagnostics.domain_cache_miss_count)},
              {"graph_cache_hits", std::to_string(m_last_tree_pm_diagnostics.graph_cache_hit_count)},
              {"graph_cache_misses", std::to_string(m_last_tree_pm_diagnostics.graph_cache_miss_count)},
              {"let_candidate_peer_count", std::to_string(m_last_tree_pm_diagnostics.let_candidate_peer_count)},
              {"communicating_peer_count", std::to_string(m_last_tree_pm_diagnostics.communicating_peer_count)},
              {"let_exported_target_count", std::to_string(m_last_tree_pm_diagnostics.let_exported_target_count)},
              {"let_imported_target_count", std::to_string(m_last_tree_pm_diagnostics.let_imported_target_count)},
              {"let_wire_bytes_sent", std::to_string(m_last_tree_pm_diagnostics.let_wire_bytes_sent)},
              {"let_wire_bytes_received", std::to_string(m_last_tree_pm_diagnostics.let_wire_bytes_received)},
              {"let_high_water_bytes", std::to_string(m_last_tree_pm_diagnostics.let_high_water_bytes)},
              {"let_discovery_ms", formatRuntimeDouble(m_last_tree_pm_diagnostics.let_discovery_ms)},
              {"let_graph_setup_ms", formatRuntimeDouble(m_last_tree_pm_diagnostics.let_graph_setup_ms)},
              {"let_communication_ms", formatRuntimeDouble(m_last_tree_pm_diagnostics.let_communication_ms)},
              {"let_overlap_local_work_ms", formatRuntimeDouble(m_last_tree_pm_diagnostics.let_overlap_local_work_ms)},
              {"let_communication_wait_ms", formatRuntimeDouble(m_last_tree_pm_diagnostics.let_communication_wait_ms)},
              {"let_overlap_efficiency", formatRuntimeDouble(m_last_tree_pm_diagnostics.let_overlap_efficiency)},
              {"let_remote_traversal_ms", formatRuntimeDouble(m_last_tree_pm_diagnostics.let_remote_traversal_ms)},
              {"local_pair_evaluations", std::to_string(tree_pm_profile.tree_profile.particle_particle_interactions)},
              {"remote_pair_evaluations", std::to_string(m_last_tree_pm_diagnostics.residual_pair_evaluations)},
              {"remote_request_packet_imbalance_ratio", formatRuntimeDouble(
                  m_last_tree_pm_diagnostics.residual_remote_request_packet_imbalance_ratio)},
              {"pm_routed_density_records", std::to_string(tree_pm_profile.pm_profile.routed_density_records)},
              {"pm_routed_force_requests", std::to_string(tree_pm_profile.pm_profile.routed_force_requests)},
              {"pm_routed_potential_requests", std::to_string(tree_pm_profile.pm_profile.routed_potential_requests)},
              {"pm_routed_density_peer_count", std::to_string(tree_pm_profile.pm_profile.routed_density_peer_count)},
              {"pm_routed_force_peer_count", std::to_string(tree_pm_profile.pm_profile.routed_force_peer_count)},
              {"pm_routed_potential_peer_count", std::to_string(tree_pm_profile.pm_profile.routed_potential_peer_count)},
              {"pm_routed_mpi_bytes_sent", std::to_string(tree_pm_profile.pm_profile.routed_mpi_bytes_sent)},
              {"pm_routed_mpi_bytes_received", std::to_string(tree_pm_profile.pm_profile.routed_mpi_bytes_received)},
              {"pm_routed_send_buffer_high_water_bytes", std::to_string(
                  tree_pm_profile.pm_profile.routed_send_buffer_high_water_bytes)},
              {"pm_routed_receive_buffer_high_water_bytes", std::to_string(
                  tree_pm_profile.pm_profile.routed_receive_buffer_high_water_bytes)},
              {"pm_routed_combined_buffer_high_water_bytes", std::to_string(
                  tree_pm_profile.pm_profile.routed_combined_buffer_high_water_bytes)},
              {"pm_routed_workspace_high_water_bytes", std::to_string(
                  tree_pm_profile.pm_profile.routed_workspace_high_water_bytes)},
              {"pm_routed_mpi_wait_ms", formatRuntimeDouble(tree_pm_profile.pm_profile.routed_mpi_wait_ms)},
          },
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
              {"force_l2_pm_global", formatRuntimeDouble(m_last_tree_pm_diagnostics.force_l2_pm_global)},
              {"force_l2_pm_zoom_correction", formatRuntimeDouble(m_last_tree_pm_diagnostics.force_l2_pm_zoom_correction)},
              {"force_l2_tree_short_range", formatRuntimeDouble(m_last_tree_pm_diagnostics.force_l2_tree_short_range)},
              {"force_l2_total", formatRuntimeDouble(m_last_tree_pm_diagnostics.force_l2_total)},
              {"zoom_high_res_source_count", std::to_string(m_last_tree_pm_diagnostics.zoom_high_res_source_count)},
              {"zoom_low_res_source_count", std::to_string(m_last_tree_pm_diagnostics.zoom_low_res_source_count)},
              {"zoom_low_res_contamination_count", std::to_string(m_last_tree_pm_diagnostics.zoom_low_res_contamination_count)},
              {"zoom_low_res_contamination_mass_code", formatRuntimeDouble(m_last_tree_pm_diagnostics.zoom_low_res_contamination_mass_code)},
              {"zoom_membership_source", m_zoom_membership_source},
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

    if (is_kick_stage) {
      applyActiveKickFromFreshForce(context);
    }

    constexpr std::uint32_t no_scatter_row =
        std::numeric_limits<std::uint32_t>::max();
    for (std::size_t target_slot = 0; target_slot < m_local_active_indices.size(); ++target_slot) {
      const std::uint32_t particle_row = m_active_target_particle_row[target_slot];
      if (particle_row != no_scatter_row) {
        m_particle_accel_x[particle_row] = m_active_accel_x[target_slot];
        m_particle_accel_y[particle_row] = m_active_accel_y[target_slot];
        m_particle_accel_z[particle_row] = m_active_accel_z[target_slot];
        m_particle_force_cache_valid[particle_row] = 1U;
      }
      const std::uint32_t cell_row = m_active_target_cell_row[target_slot];
      if (cell_row != no_scatter_row) {
        m_cell_accel_x[cell_row] = m_active_accel_x[target_slot];
        m_cell_accel_y[cell_row] = m_active_accel_y[target_slot];
        m_cell_accel_z[cell_row] = m_active_accel_z[target_slot];
        m_cell_force_cache_valid[cell_row] = 1U;
      }
    }

    m_force_cache_particle_index_generation = context.state.particleIndexGeneration();
    m_force_cache_cell_index_generation = context.state.cellIndexGeneration();
    m_force_cache_gas_identity_generation = context.state.gasCellIdentityGeneration();
    m_force_cache_source_generation = context.state.gravitySourceGeneration();
    m_force_cache_valid = true;
    const std::uint32_t local_rank =
        static_cast<std::uint32_t>(std::max(mpi_context.worldRank(), 0));
    const std::uint32_t gas_species_tag =
        static_cast<std::uint32_t>(core::ParticleSpecies::kGas);
    for (std::size_t particle_index = 0; particle_index < particle_count;
         ++particle_index) {
      if (context.state.particle_sidecar.species_tag[particle_index] ==
          gas_species_tag) {
        continue;
      }
      if (context.state.particle_sidecar.owning_rank[particle_index] ==
              local_rank &&
          m_particle_force_cache_valid[particle_index] == 0U) {
        m_force_cache_valid = false;
        break;
      }
    }
    if (m_force_cache_valid) {
      for (std::size_t cell_index = 0; cell_index < cell_count; ++cell_index) {
        if (m_owned_leaf_cell_mask[cell_index] != 0U &&
            m_cell_force_cache_valid[cell_index] == 0U) {
          m_force_cache_valid = false;
          break;
        }
      }
    }

  }

 private:
  [[nodiscard]] static double kickFactorForStage(const core::StepContext& context) {
    if (context.stage == core::IntegrationStage::kGravityKickPre) {
      return context.timeline_step.first_kick_factor_code;
    }
    if (context.stage == core::IntegrationStage::kGravityKickPost) {
      return context.timeline_step.second_kick_factor_code;
    }
    throw std::invalid_argument("gravity kick factor requested outside KDK kick stage");
  }

  [[nodiscard]] static double hubbleDragFactorForStage(const core::StepContext& context) {
    if (context.stage == core::IntegrationStage::kGravityKickPre) {
      return context.timeline_step.first_hubble_drag_factor;
    }
    if (context.stage == core::IntegrationStage::kGravityKickPost) {
      return context.timeline_step.second_hubble_drag_factor;
    }
    throw std::invalid_argument("Hubble drag factor requested outside KDK kick stage");
  }

  static void applyPeculiarVelocityKick(
      core::SimulationState& state,
      std::uint32_t particle_index,
      double accel_x_code,
      double accel_y_code,
      double accel_z_code,
      double kick_factor_code,
      double hubble_drag_factor) {
    state.particles.velocity_x_peculiar[particle_index] =
        hubble_drag_factor * state.particles.velocity_x_peculiar[particle_index] + kick_factor_code * accel_x_code;
    state.particles.velocity_y_peculiar[particle_index] =
        hubble_drag_factor * state.particles.velocity_y_peculiar[particle_index] + kick_factor_code * accel_y_code;
    state.particles.velocity_z_peculiar[particle_index] =
        hubble_drag_factor * state.particles.velocity_z_peculiar[particle_index] + kick_factor_code * accel_z_code;
  }

  static void applyGasCellKdkKick(
      core::SimulationState& state,
      std::uint32_t cell_index,
      double accel_x_code,
      double accel_y_code,
      double accel_z_code,
      double kick_factor_code,
      double hubble_drag_factor) {
    if (cell_index >= state.gas_cells.size()) {
      throw std::out_of_range("gas KDK kick cell index is outside GasCellSidecar");
    }
    if (!std::isfinite(hubble_drag_factor) || hubble_drag_factor <= 0.0) {
      throw std::invalid_argument("gas KDK kick requires finite positive Hubble drag factor");
    }
    state.gas_cells.velocity_x_peculiar[cell_index] =
        hubble_drag_factor * state.gas_cells.velocity_x_peculiar[cell_index] +
        kick_factor_code * accel_x_code;
    state.gas_cells.velocity_y_peculiar[cell_index] =
        hubble_drag_factor * state.gas_cells.velocity_y_peculiar[cell_index] +
        kick_factor_code * accel_y_code;
    state.gas_cells.velocity_z_peculiar[cell_index] =
        hubble_drag_factor * state.gas_cells.velocity_z_peculiar[cell_index] +
        kick_factor_code * accel_z_code;

    // For gamma=5/3, homogeneous adiabatic expansion gives
    // specific internal energy u_int propto a^-2. The timeline drag factor is
    // a_begin/a_end for this half-kick, so thermal state scales as drag^2.
    // Density/mass are unchanged by a source kick; the finite-volume transport
    // stage owns those conserved quantities.
    constexpr double gamma_minus_one = 2.0 / 3.0;
    const double thermal_drag =
        std::pow(hubble_drag_factor, 3.0 * gamma_minus_one);
    state.gas_cells.internal_energy_code[cell_index] *= thermal_drag;
    state.gas_cells.pressure_code[cell_index] *= thermal_drag;
    state.gas_cells.temperature_code[cell_index] *= thermal_drag;
    state.gas_cells.sound_speed_code[cell_index] *= std::sqrt(thermal_drag);
  }

  void applyActiveKickFromFreshForce(core::StepContext& context) {
    const double kick_factor = kickFactorForStage(context);
    const double hubble_drag = hubbleDragFactorForStage(context);
    for (const std::uint32_t particle_index :
         context.active_set.particle_indices) {
      if (particle_index >= m_active_slot_by_particle.size()) {
        throw std::out_of_range(
            "gravity kick particle index is outside the fresh-force slot map");
      }
      const int active_slot = m_active_slot_by_particle[particle_index];
      if (active_slot < 0) {
        continue;
      }
      applyPeculiarVelocityKick(
          context.state,
          particle_index,
          m_active_accel_x[static_cast<std::size_t>(active_slot)],
          m_active_accel_y[static_cast<std::size_t>(active_slot)],
          m_active_accel_z[static_cast<std::size_t>(active_slot)],
          kick_factor,
          hubble_drag);
    }
    for (const std::uint32_t cell_index : context.active_set.cell_indices) {
      if (cell_index >= m_active_slot_by_cell.size()) {
        throw std::out_of_range("gas KDK fresh-force cell index is outside the active slot map");
      }
      const int active_slot = m_active_slot_by_cell[cell_index];
      if (active_slot < 0) {
        continue;
      }
      applyGasCellKdkKick(
          context.state,
          cell_index,
          m_active_accel_x[static_cast<std::size_t>(active_slot)],
          m_active_accel_y[static_cast<std::size_t>(active_slot)],
          m_active_accel_z[static_cast<std::size_t>(active_slot)],
          kick_factor,
          hubble_drag);
    }
    if (context.stage == core::IntegrationStage::kGravityKickPost &&
        !context.active_set.cell_indices.empty()) {
      internal::synchronizeParentParticleCompatibilityMirrors(
          context.state,
          static_cast<std::uint32_t>(std::max(m_services.mpi_context.worldRank(), 0)),
          "gravity post-kick gas compatibility mirror");
    }
  }

  void applyCachedKick(core::StepContext& context) {
    if (!m_force_cache_valid) {
      throw std::runtime_error("gravity kick requested before a coherent force-refresh boundary");
    }
    if (m_force_cache_particle_index_generation != context.state.particleIndexGeneration() ||
        m_force_cache_cell_index_generation != context.state.cellIndexGeneration() ||
        m_force_cache_gas_identity_generation != context.state.gasCellIdentityGeneration() ||
        m_force_cache_source_generation != context.state.gravitySourceGeneration()) {
      throw std::runtime_error(
          "gravity cached kick rejected because particle/cell source identity changed since force refresh");
    }
    const std::size_t particle_count = context.state.particles.size();
    if (m_particle_accel_x.size() != particle_count ||
        m_particle_accel_y.size() != particle_count ||
        m_particle_accel_z.size() != particle_count ||
        m_particle_force_cache_valid.size() != particle_count) {
      throw std::runtime_error(
          "gravity cached kick rejected an incomplete particle acceleration cache");
    }
    const std::uint32_t local_rank =
        static_cast<std::uint32_t>(std::max(m_services.mpi_context.worldRank(), 0));
    const std::uint32_t gas_species_tag =
        static_cast<std::uint32_t>(core::ParticleSpecies::kGas);
    const double kick_factor = kickFactorForStage(context);
    const double hubble_drag = hubbleDragFactorForStage(context);
    for (const std::uint32_t particle_index : context.active_set.particle_indices) {
      if (particle_index >= particle_count) {
        throw std::out_of_range("gravity cached kick particle index out of range");
      }
      if (context.state.particle_sidecar.owning_rank[particle_index] != local_rank ||
          context.state.particle_sidecar.species_tag[particle_index] == gas_species_tag) {
        continue;
      }
      if (m_particle_force_cache_valid[particle_index] == 0U) {
        throw std::runtime_error(
            "gravity cached kick rejected an invalid per-particle force row");
      }
      applyPeculiarVelocityKick(
          context.state,
          particle_index,
          m_particle_accel_x[particle_index],
          m_particle_accel_y[particle_index],
          m_particle_accel_z[particle_index],
          kick_factor,
          hubble_drag);
    }
    if (m_cell_accel_x.size() != context.state.cells.size() ||
        m_cell_accel_y.size() != context.state.cells.size() ||
        m_cell_accel_z.size() != context.state.cells.size() ||
        m_cell_force_cache_valid.size() != context.state.cells.size()) {
      throw std::runtime_error(
          "gravity cached kick rejected an incomplete gas-cell acceleration cache");
    }
    const auto authoritative_rows = internal::selectAuthoritativeGravitySourceRows(
        context.state,
        local_rank,
        static_cast<std::uint32_t>(std::max(m_services.mpi_context.worldSize(), 1)),
        "gas KDK cached kick");
    std::vector<std::uint8_t> is_owned_leaf(context.state.cells.size(), 0U);
    for (const std::uint32_t cell_row : authoritative_rows.gas_cell_rows) {
      is_owned_leaf[cell_row] = 1U;
    }
    for (const std::uint32_t cell_index : context.active_set.cell_indices) {
      if (cell_index >= context.state.cells.size()) {
        throw std::out_of_range("gas KDK cached-force cell index is outside CellSoa");
      }
      if (is_owned_leaf[cell_index] == 0U) {
        continue;
      }
      if (m_cell_force_cache_valid[cell_index] == 0U) {
        throw std::runtime_error(
            "gas KDK cached kick rejected an invalid per-cell force row");
      }
      applyGasCellKdkKick(
          context.state,
          cell_index,
          m_cell_accel_x[cell_index],
          m_cell_accel_y[cell_index],
          m_cell_accel_z[cell_index],
          kick_factor,
          hubble_drag);
    }
    if (context.stage == core::IntegrationStage::kGravityKickPost &&
        !context.active_set.cell_indices.empty()) {
      internal::synchronizeParentParticleCompatibilityMirrors(
          context.state, local_rank, "gravity post-kick gas compatibility mirror");
    }
  }

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
    if (m_local_active_indices.size() != m_active_target_particle_row.size() ||
        m_local_active_indices.size() != m_active_target_cell_row.size()) {
      ++summary.decomposition_sanity_failure_count;
      failFatal("local gravity target vectors diverged after compact-view rebuild");
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
    if (context.integrator_state.pm_sync_state.lastRefreshOpportunity() > context.integrator_state.pm_sync_state.gravityKickOpportunity() ||
        context.integrator_state.pm_sync_state.fieldVersion() < m_last_committed_field_version) {
      ++summary.illegal_sync_state_count;
      failFatal("gravity sync-state regressed across kick opportunities");
    }
    if (decision.refresh_long_range_field) {
      if (decision.field_version != context.integrator_state.pm_sync_state.fieldVersion() + 1U) {
        ++summary.illegal_sync_state_count;
        failFatal("refresh decision must increment field_version by exactly one");
      }
    } else if (decision.field_version != context.integrator_state.pm_sync_state.fieldVersion()) {
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

  void ensureZoomMembershipLoaded(const core::SimulationState& state) {
    if (m_zoom_membership_loaded) {
      return;
    }
    m_zoom_membership_loaded = true;
    if (!m_mode_policy.zoom_region_expected || m_config.mode.zoom_region_file.empty()) {
      m_zoom_membership_source = m_mode_policy.zoom_region_expected ? "configured_spherical_region" : "disabled";
      return;
    }
    const std::filesystem::path path = m_zoom_region_path.empty()
        ? std::filesystem::path(m_config.mode.zoom_region_file)
        : m_zoom_region_path;
    if (hasHdf5Extension(path)) {
#if COSMOSIM_ENABLE_HDF5
      m_zoom_high_res_particle_ids = loadZoomParticleIdsFromHdf5(path);
      m_zoom_membership_source = "particle_id_file_hdf5";
#else
      throw std::runtime_error("zoom_region_file requires HDF5 support in this build");
#endif
    } else {
      m_zoom_high_res_particle_ids = loadZoomParticleIdsFromText(path);
      m_zoom_membership_source = "particle_id_file_text";
    }
    for (const std::uint64_t id : state.particle_sidecar.particle_id) {
      (void)id;
    }
  }

  enum class SourcePredictionEpoch : std::uint8_t {
    kNone = 0,
    kStepBegin = 1,
    kStepEnd = 2,
  };

  void rebuildOwnedParticleCompactView(
      const core::StepContext& context,
      std::span<const std::uint32_t> active_particles,
      SourcePredictionEpoch prediction_epoch = SourcePredictionEpoch::kNone) {
    const core::SimulationState& state = context.state;
    state.requireGasCellIdentityMapCoversDenseRows(
        "gravity authoritative gas source rebuild");
    const std::uint32_t world_rank =
        static_cast<std::uint32_t>(m_runtime_topology.world_rank);
    const std::size_t particle_count = state.particles.size();
    const std::size_t cell_count = state.cells.size();
    constexpr std::uint32_t no_row = std::numeric_limits<std::uint32_t>::max();

    m_owned_local_index_by_global.assign(particle_count, -1);
    m_owned_local_index_by_cell.assign(cell_count, -1);
    m_owned_leaf_cell_mask.assign(cell_count, 0U);
    m_local_source_x.clear();
    m_local_source_y.clear();
    m_local_source_z.clear();
    m_local_source_mass.clear();
    m_local_source_species_tag.clear();
    m_local_source_softening_comoving.clear();
    m_local_source_softening_override_mask.clear();
    m_local_source_particle_row.clear();
    m_local_source_cell_row.clear();
    m_source_predicted_inactive_count = 0;

    std::vector<std::uint8_t> active_mask;
    if (prediction_epoch != SourcePredictionEpoch::kNone) {
      if (state.particle_sidecar.last_drift_time_code.size() != particle_count ||
          state.particle_sidecar.last_drift_scale_factor.size() != particle_count) {
        throw std::runtime_error(
            "PM source prediction requires per-particle drift epoch sidecars");
      }
      active_mask.assign(particle_count, 0U);
      for (const std::uint32_t global_index : active_particles) {
        if (global_index >= particle_count) {
          throw std::out_of_range(
              "gravity callback active particle index out of range");
        }
        active_mask[global_index] = 1U;
      }
    }

    const bool periodic_sources =
        m_tree_pm_options.pm_options.boundary_condition ==
        gravity::PmBoundaryCondition::kPeriodic;
    const double box_size_x = m_config.cosmology.box_size_x_mpc_comoving;
    const double box_size_y = m_config.cosmology.box_size_y_mpc_comoving;
    const double box_size_z = m_config.cosmology.box_size_z_mpc_comoving;
    const bool has_softening_values =
        !state.particle_sidecar.gravity_softening_comoving.empty();
    const bool has_softening_masks =
        !state.particle_sidecar.has_gravity_softening_override.empty();
    const std::uint32_t gas_species_tag =
        static_cast<std::uint32_t>(core::ParticleSpecies::kGas);
    const auto authoritative_source_rows =
        internal::selectAuthoritativeGravitySourceRows(
            state,
            world_rank,
            static_cast<std::uint32_t>(std::max(m_runtime_topology.world_size, 1)),
            "gravity authoritative source selection");
    const double gas_softening_comoving = m_tree_pm_species_softening.enabled
        ? m_tree_pm_species_softening.epsilon_comoving_by_species[
              static_cast<std::size_t>(core::ParticleSpecies::kGas)]
        : m_tree_pm_options.tree_options.softening.epsilon_comoving;
    constexpr double max_certified_softening_to_split_ratio = 0.20;
    const double split_scale_comoving =
        m_tree_pm_options.split_policy.split_scale_comoving;

    const auto validateCertifiedTreePmSoftening = [&] (
        std::uint32_t species_tag,
        double softening_value,
        std::uint8_t softening_override) {
      double effective_softening =
          m_tree_pm_options.tree_options.softening.epsilon_comoving;
      if (softening_override != 0U) {
        effective_softening = softening_value;
      } else if (m_tree_pm_species_softening.enabled) {
        if (!core::isValidParticleSpeciesTag(species_tag)) {
          throw std::invalid_argument(
              "gravity source species tag is outside the canonical range");
        }
        effective_softening =
            m_tree_pm_species_softening.epsilon_comoving_by_species[species_tag];
      }
      if (!std::isfinite(effective_softening) || effective_softening < 0.0) {
        throw std::invalid_argument(
            "gravity source effective softening must be finite and non-negative");
      }
      if (effective_softening >
          max_certified_softening_to_split_ratio * split_scale_comoving) {
        throw std::invalid_argument(
            "TreePM softening exceeds the certified epsilon/r_s <= 0.20 envelope; "
            "reduce gravitational softening or increase the TreePM split scale");
      }
    };

    const auto appendSource = [&](double source_x,
                                  double source_y,
                                  double source_z,
                                  double source_mass,
                                  std::uint32_t species_tag,
                                  double softening_value,
                                  std::uint8_t softening_override,
                                  std::uint32_t particle_row,
                                  std::uint32_t cell_row) {
      if (!std::isfinite(source_x) || !std::isfinite(source_y) ||
          !std::isfinite(source_z) || !std::isfinite(source_mass) ||
          source_mass < 0.0) {
        throw std::runtime_error(
            "gravity source rebuild encountered non-finite coordinates or a negative/non-finite mass");
      }
      validateCertifiedTreePmSoftening(
          species_tag, softening_value, softening_override);
      if (m_local_source_x.size() >=
          static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max())) {
        throw std::overflow_error(
            "gravity local source count exceeds TreePM 32-bit source-index contract");
      }
      if (periodic_sources) {
        source_x = wrapPeriodicPosition(source_x, box_size_x);
        source_y = wrapPeriodicPosition(source_y, box_size_y);
        source_z = wrapPeriodicPosition(source_z, box_size_z);
      }
      const auto source_index =
          static_cast<std::uint32_t>(m_local_source_x.size());
      m_local_source_x.push_back(source_x);
      m_local_source_y.push_back(source_y);
      m_local_source_z.push_back(source_z);
      m_local_source_mass.push_back(source_mass);
      m_local_source_species_tag.push_back(species_tag);
      m_local_source_particle_row.push_back(particle_row);
      m_local_source_cell_row.push_back(cell_row);
      if (has_softening_values) {
        m_local_source_softening_comoving.push_back(softening_value);
      }
      if (has_softening_masks) {
        m_local_source_softening_override_mask.push_back(softening_override);
      }
      return source_index;
    };

    // Collisionless and compact-object particle state remains particle-owned.
    // Gas-tagged generic particles are compatibility/lineage mirrors only and
    // must not duplicate the authoritative gas-cell mass in the gravity source
    // population. Remote ghost particle copies are also excluded: TreePM owns
    // its explicit remote-source exchange from each rank's authoritative
    // owned sources.
    for (const std::uint32_t global_index :
         authoritative_source_rows.particle_rows) {
      double source_x = state.particles.position_x_comoving[global_index];
      double source_y = state.particles.position_y_comoving[global_index];
      double source_z = state.particles.position_z_comoving[global_index];
      if (prediction_epoch != SourcePredictionEpoch::kNone &&
          active_mask[global_index] == 0U) {
        const double source_time =
            state.particle_sidecar.last_drift_time_code[global_index];
        const double source_scale =
            state.particle_sidecar.last_drift_scale_factor[global_index];
        const double evaluation_time =
            prediction_epoch == SourcePredictionEpoch::kStepBegin
            ? context.timeline_step.time_begin_code
            : context.timeline_step.time_end_code;
        const double evaluation_scale =
            prediction_epoch == SourcePredictionEpoch::kStepBegin
            ? context.timeline_step.scale_factor_begin
            : context.timeline_step.scale_factor_end;
        if (source_time > evaluation_time + 1.0e-12 || source_scale <= 0.0) {
          throw std::runtime_error(
              "inactive PM source has an invalid or future drift epoch");
        }
        double inactive_drift_factor = evaluation_time - source_time;
        if (context.cosmology_background != nullptr) {
          inactive_drift_factor = core::computeComovingDriftFactor(
              *context.cosmology_background,
              source_scale,
              evaluation_scale,
              64) /
              context.timeline_step.time_si_per_code;
        }
        if (!std::isfinite(inactive_drift_factor) ||
            inactive_drift_factor < -1.0e-14) {
          throw std::runtime_error(
              "inactive PM source prediction produced an invalid drift factor");
        }
        source_x += state.particles.velocity_x_peculiar[global_index] *
            inactive_drift_factor;
        source_y += state.particles.velocity_y_peculiar[global_index] *
            inactive_drift_factor;
        source_z += state.particles.velocity_z_peculiar[global_index] *
            inactive_drift_factor;
        ++m_source_predicted_inactive_count;
      }
      const double softening_value = has_softening_values
          ? state.particle_sidecar.gravity_softening_comoving[global_index]
          : 0.0;
      const std::uint8_t softening_override = has_softening_masks &&
              state.particle_sidecar.hasGravitySofteningOverride(global_index)
          ? 1U
          : 0U;
      const std::uint32_t source_index = appendSource(
          source_x,
          source_y,
          source_z,
          state.particles.mass_code[global_index],
          state.particle_sidecar.species_tag[global_index],
          softening_value,
          softening_override,
          static_cast<std::uint32_t>(global_index),
          no_row);
      m_owned_local_index_by_global[global_index] =
          static_cast<int>(source_index);
    }

    // Gas-cell state is authoritative. The shared selector excludes covered
    // coarse cells and never relies on optional parent-particle lineage.
    for (const std::uint32_t cell_index :
         authoritative_source_rows.gas_cell_rows) {
      m_owned_leaf_cell_mask[cell_index] = 1U;
      const auto gas_source = internal::authoritativeGasGravitySource(
          state, cell_index, "gravity authoritative gas source rebuild");
      const std::uint32_t source_index = appendSource(
          gas_source.x_comoving,
          gas_source.y_comoving,
          gas_source.z_comoving,
          gas_source.mass_code,
          gas_species_tag,
          gas_softening_comoving,
          0U,
          no_row,
          cell_index);
      m_owned_local_index_by_cell[cell_index] = static_cast<int>(source_index);
    }

    m_local_active_indices.clear();
    m_local_active_global_indices.clear();
    m_active_target_particle_row.clear();
    m_active_target_cell_row.clear();
    m_active_target_species_tag.clear();
    m_active_target_softening_comoving.clear();
    m_active_target_softening_override_mask.clear();

    const auto appendTarget = [&](std::uint32_t source_index,
                                  std::uint32_t particle_row,
                                  std::uint32_t cell_row) {
      m_local_active_indices.push_back(source_index);
      m_active_target_particle_row.push_back(particle_row);
      m_active_target_cell_row.push_back(cell_row);
      m_active_target_species_tag.push_back(
          m_local_source_species_tag[source_index]);
      if (has_softening_values) {
        m_active_target_softening_comoving.push_back(
            m_local_source_softening_comoving[source_index]);
      }
      if (has_softening_masks) {
        m_active_target_softening_override_mask.push_back(
            m_local_source_softening_override_mask[source_index]);
      }
    };

    for (const std::uint32_t global_index : active_particles) {
      if (global_index >= particle_count) {
        throw std::out_of_range(
            "gravity callback active particle index out of range");
      }
      if (state.particle_sidecar.species_tag[global_index] == gas_species_tag) {
        continue;
      }
      const int local_index = m_owned_local_index_by_global[global_index];
      if (local_index < 0) {
        continue;
      }
      const auto source_index = static_cast<std::uint32_t>(local_index);
      appendTarget(source_index, global_index, no_row);
      m_local_active_global_indices.push_back(global_index);
    }

    // Rung-zero production currently owns globally synchronized force
    // refreshes. Compute every owned leaf gas-cell target at those boundaries
    // so hydro and restart caches never depend on an optional lineage parent.
    // Future local gas time bins can narrow this list only after carrying a
    // per-cell force epoch through the scheduler contract.
    for (std::size_t cell_index = 0; cell_index < cell_count; ++cell_index) {
      const int local_index = m_owned_local_index_by_cell[cell_index];
      if (local_index < 0) {
        continue;
      }
      appendTarget(
          static_cast<std::uint32_t>(local_index),
          no_row,
          static_cast<std::uint32_t>(cell_index));
    }
  }

  const RuntimeServices& m_services;
  const core::SimulationConfig& m_config;
  const core::ModePolicy& m_mode_policy;
  parallel::DistributedExecutionTopology m_runtime_topology{};
  std::uint64_t m_pm_update_cadence_steps = 1;
  gravity::PmGridShape m_pm_grid_shape{};
  double m_mesh_spacing_x_mpc_comoving = 0.0;
  double m_mesh_spacing_y_mpc_comoving = 0.0;
  double m_mesh_spacing_z_mpc_comoving = 0.0;
  double m_gravitational_constant_code = 0.0;
  std::string m_pm_assignment_scheme = "unknown";
  std::string m_pm_backend = "unknown";
  gravity::TreePmCoordinator m_tree_pm_coordinator;
  gravity::TreePmOptions m_tree_pm_options;
  gravity::TreeSofteningSpeciesPolicy m_tree_pm_species_softening{};
  std::vector<double> m_active_accel_x;
  std::vector<double> m_active_accel_y;
  std::vector<double> m_active_accel_z;
  std::vector<double> m_active_previous_acceleration_magnitude;
  std::vector<double> m_particle_accel_x;
  std::vector<double> m_particle_accel_y;
  std::vector<double> m_particle_accel_z;
  std::vector<double> m_cell_accel_x;
  std::vector<double> m_cell_accel_y;
  std::vector<double> m_cell_accel_z;
  std::vector<int> m_active_slot_by_particle;
  std::vector<int> m_active_slot_by_cell;
  std::vector<int> m_owned_local_index_by_global;
  std::vector<int> m_owned_local_index_by_cell;
  std::vector<std::uint8_t> m_owned_leaf_cell_mask;
  std::vector<double> m_local_source_x;
  std::vector<double> m_local_source_y;
  std::vector<double> m_local_source_z;
  std::vector<double> m_local_source_mass;
  std::vector<std::uint32_t> m_local_source_species_tag;
  std::vector<double> m_local_source_softening_comoving;
  std::vector<std::uint8_t> m_local_source_softening_override_mask;
  std::vector<std::uint32_t> m_local_source_particle_row;
  std::vector<std::uint32_t> m_local_source_cell_row;
  std::vector<std::uint32_t> m_active_target_particle_row;
  std::vector<std::uint32_t> m_active_target_cell_row;
  std::vector<std::uint32_t> m_active_target_species_tag;
  std::vector<double> m_active_target_softening_comoving;
  std::vector<std::uint8_t> m_active_target_softening_override_mask;
  std::vector<std::uint8_t> m_source_is_high_res;
  std::vector<std::uint8_t> m_active_is_high_res;
  std::filesystem::path m_zoom_region_path;
  bool m_zoom_membership_loaded = false;
  std::unordered_set<std::uint64_t> m_zoom_high_res_particle_ids;
  std::string m_zoom_membership_source = "disabled";
  std::vector<std::uint32_t> m_local_active_indices;
  std::vector<std::uint32_t> m_local_active_global_indices;
  std::vector<std::uint32_t> m_force_refresh_particle_indices;
  std::vector<parallel::TopDomainLeaf> m_authoritative_top_domain_leaves;
  std::size_t m_local_kick_particle_count = 0U;
  gravity::TreePmDiagnostics m_last_tree_pm_diagnostics{};
  parallel::DecompositionRuntimeMeasurements m_last_decomposition_measurements{};
  bool m_has_long_range_field = false;
  std::uint64_t m_last_committed_field_version = 0;
  std::uint64_t m_long_range_refresh_count = 0;
  std::uint64_t m_long_range_reuse_count = 0;
  bool m_force_cache_valid = false;
  std::uint64_t m_force_cache_particle_index_generation = 0U;
  std::uint64_t m_force_cache_cell_index_generation = 0U;
  std::uint64_t m_force_cache_gas_identity_generation = 0U;
  std::uint64_t m_force_cache_source_generation = 0U;
  std::vector<std::uint8_t> m_particle_force_cache_valid;
  std::vector<std::uint8_t> m_cell_force_cache_valid;
  gravity::DecompositionEpoch m_decomposition_epoch{};
  gravity::GravityMemoryEstimate m_last_pre_run_memory_estimate{};
  bool m_has_pre_run_memory_estimate = false;
  gravity::DmoProcessMemoryEstimate m_last_process_memory_estimate{};
  bool m_has_process_memory_estimate = false;
  parallel::GhostCacheLifecycle m_ghost_cache_lifecycle{};
  std::uint64_t m_source_predicted_inactive_count = 0;
  std::vector<ReferenceWorkflowReport::TreePmCadenceRecord> m_cadence_records;
};

}  // namespace

void initializeGravityState(
    core::SimulationState& state,
    const core::SimulationConfig& config) {
  maybeInitializeParticleSofteningFromSpeciesPolicy(state, config);
}

std::string describeGravitySofteningPolicy(
    const core::SimulationConfig& config,
    const core::SimulationState& state) {
  return describeSofteningPolicy(config, &state);
}

std::unique_ptr<GravityRuntime> makeGravityRuntime(
    const core::SimulationConfig& config,
    const core::ModePolicy& mode_policy,
    const RuntimeServices& services,
    std::filesystem::path zoom_region_path) {
  return std::make_unique<GravityRuntimeImpl>(
      config, mode_policy, services, std::move(zoom_region_path));
}

}  // namespace cosmosim::workflows
