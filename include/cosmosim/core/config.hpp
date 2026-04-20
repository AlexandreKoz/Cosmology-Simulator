#pragma once

#include <cstdint>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <vector>

namespace cosmosim::core {

enum class SimulationMode {
  kCosmoCube,
  kZoomIn,
  kIsolatedGalaxy,
  kIsolatedCluster,
};

enum class GravitySolver {
  kTreePm,
};

enum class HydroSolver {
  kGodunovFv,
};

enum class TreePmAssignmentScheme {
  kCic,
  kTsc,
};

enum class PmDecompositionMode {
  kSlab,
  kPencil,
};

enum class CoordinateFrame {
  kComoving,
  kPhysical,
};

enum class ModeHydroBoundary {
  kAuto,
  kPeriodic,
  kOpen,
  kReflective,
};

enum class ModeGravityBoundary {
  kAuto,
  kPeriodic,
  kIsolatedMonopole,
};

enum class ZoomLongRangeStrategy {
  kDisabled,
  kGlobalCoarsePlusFocusedHighResCorrection,
};

enum class FeedbackMode {
  kThermal,
  kKinetic,
  kMomentum,
  kThermalKineticMomentum,
};

enum class FeedbackVariant {
  kNone,
  kDelayedCooling,
  kStochastic,
};

enum class UvBackgroundModel {
  kNone,
  kHm12,
  kFg20,
};

enum class SelfShieldingModel {
  kNone,
  kRahmati13Like,
};

enum class CoolingModel {
  kPrimordial,
  kPrimordialMetalLine,
};

struct CosmologyConfig {
  double omega_matter = 0.315;
  double omega_lambda = 0.685;
  double omega_baryon = 0.049;
  double hubble_param = 0.674;
  double sigma8 = 0.811;
  double scalar_index_ns = 0.965;
  double box_size_x_mpc_comoving = 50.0;
  double box_size_y_mpc_comoving = 50.0;
  double box_size_z_mpc_comoving = 50.0;
  // Legacy scalar compatibility lane; canonical runtime geometry must use axis-aware fields above.
  double box_size_mpc_comoving = 50.0;
};

struct NumericsConfig {
  double time_begin_code = 0.0;
  double time_end_code = 1.0;
  int max_global_steps = 1024;
  int hierarchical_max_rung = 12;
  int amr_max_level = 10;
  double gravity_softening_kpc_comoving = 1.0;
  GravitySolver gravity_solver = GravitySolver::kTreePm;
  HydroSolver hydro_solver = HydroSolver::kGodunovFv;
  int treepm_pm_grid_nx = 16;
  int treepm_pm_grid_ny = 16;
  int treepm_pm_grid_nz = 16;
  // Legacy scalar compatibility lane; canonical runtime geometry must use axis-aware fields above.
  int treepm_pm_grid = 16;
  double treepm_asmth_cells = 1.25;
  double treepm_rcut_cells = 4.5;
  TreePmAssignmentScheme treepm_assignment_scheme = TreePmAssignmentScheme::kCic;
  bool treepm_enable_window_deconvolution = false;
  int treepm_update_cadence_steps = 1;
  PmDecompositionMode treepm_pm_decomposition_mode = PmDecompositionMode::kSlab;
  std::uint64_t treepm_tree_exchange_batch_bytes = 4ULL * 1024ULL * 1024ULL;
};

struct PhysicsConfig {
  bool enable_cooling = true;
  bool enable_star_formation = true;
  bool enable_feedback = true;
  bool enable_stellar_evolution = true;
  std::string reionization_model = "hm12";
  UvBackgroundModel uv_background_model = UvBackgroundModel::kHm12;
  SelfShieldingModel self_shielding_model = SelfShieldingModel::kNone;
  CoolingModel cooling_model = CoolingModel::kPrimordial;
  std::string metal_line_table_path;
  double temperature_floor_k = 100.0;
  double sf_density_threshold_code = 10.0;
  double sf_temperature_threshold_k = 1.0e4;
  double sf_min_converging_flow_rate_code = 0.0;
  double sf_epsilon_ff = 0.01;
  double sf_min_star_particle_mass_code = 0.1;
  bool sf_stochastic_spawning = true;
  std::uint64_t sf_random_seed = 123456789ull;
  FeedbackMode fb_mode = FeedbackMode::kThermalKineticMomentum;
  FeedbackVariant fb_variant = FeedbackVariant::kNone;
  bool fb_use_returned_mass_budget = true;
  double fb_epsilon_thermal = 0.6;
  double fb_epsilon_kinetic = 0.3;
  double fb_epsilon_momentum = 0.1;
  double fb_sn_energy_erg_per_mass_code = 1.0e49;
  double fb_momentum_code_per_mass_code = 3.0e3;
  std::uint32_t fb_neighbor_count = 8;
  double fb_delayed_cooling_time_code = 0.0;
  double fb_stochastic_event_probability = 0.25;
  std::uint64_t fb_random_seed = 42424242ull;
  std::string stellar_evolution_table_path;
  double stellar_evolution_hubble_time_years = 1.44e10;
  bool enable_black_hole_agn = false;
  double bh_seed_halo_mass_threshold_code = 1.0e3;
  double bh_seed_mass_code = 1.0;
  std::uint32_t bh_seed_max_per_cell = 1;
  double bh_alpha_bondi = 1.0;
  bool bh_use_eddington_cap = true;
  double bh_epsilon_r = 0.1;
  double bh_epsilon_f = 0.05;
  double bh_feedback_coupling_efficiency = 1.0;
  double bh_duty_cycle_active_edd_ratio_threshold = 0.01;
  double bh_proton_mass_si = 1.67262192369e-27;
  double bh_thomson_cross_section_si = 6.6524587321e-29;
  double bh_newton_g_si = 6.67430e-11;
  double bh_speed_of_light_si = 2.99792458e8;
  bool enable_tracers = false;
  bool tracer_track_mass = true;
  double tracer_min_host_mass_code = 0.0;
};

struct OutputConfig {
  std::string run_name = "cosmosim_run";
  std::string output_directory = "outputs";
  std::string output_stem = "snapshot";
  std::string restart_stem = "restart";
  int snapshot_interval_steps = 64;
  bool write_restarts = true;
};


struct AnalysisConfig {
  enum class DiagnosticsExecutionPolicy : std::uint8_t {
    kRunHealthOnly = 0,
    kRunHealthAndLightScience = 1,
    kAllIncludingProvisional = 2,
  };

  bool enable_diagnostics = true;
  bool enable_halo_workflow = false;
  bool halo_on_the_fly = false;
  DiagnosticsExecutionPolicy diagnostics_execution_policy =
      DiagnosticsExecutionPolicy::kRunHealthAndLightScience;
  int run_health_interval_steps = 1;
  int science_light_interval_steps = 8;
  int science_heavy_interval_steps = 64;
  int retention_bundle_count = 8;
  int power_spectrum_mesh_n = 16;
  int power_spectrum_bin_count = 12;
  int sf_history_bin_count = 16;
  int quicklook_grid_n = 32;
  std::string diagnostics_stem = "diagnostics";
  std::string halo_catalog_stem = "halo_catalog";
  std::string merger_tree_stem = "merger_tree_plan";
  double halo_fof_linking_length_factor = 0.2;
  int halo_fof_min_group_size = 16;
  bool halo_include_gas = true;
  bool halo_include_stars = true;
  bool halo_include_black_holes = true;
};

struct ParallelConfig {
  int mpi_ranks_expected = 1;
  int omp_threads = 1;
  int gpu_devices = 0;
  bool deterministic_reduction = true;
};

struct UnitsConfig {
  std::string length_unit = "mpc";
  std::string mass_unit = "msun";
  std::string velocity_unit = "km_s";
  CoordinateFrame coordinate_frame = CoordinateFrame::kComoving;
};

struct ModeConfig {
  SimulationMode mode = SimulationMode::kZoomIn;
  std::string ic_file = "ics.hdf5";
  bool zoom_high_res_region = false;
  std::string zoom_region_file;
  ZoomLongRangeStrategy zoom_long_range_strategy =
      ZoomLongRangeStrategy::kDisabled;
  double zoom_region_center_x_mpc_comoving = 0.5;
  double zoom_region_center_y_mpc_comoving = 0.5;
  double zoom_region_center_z_mpc_comoving = 0.5;
  double zoom_region_radius_mpc_comoving = 0.0;
  int zoom_focused_pm_grid_nx = 0;
  int zoom_focused_pm_grid_ny = 0;
  int zoom_focused_pm_grid_nz = 0;
  double zoom_contamination_radius_mpc_comoving = 0.0;
  ModeHydroBoundary hydro_boundary = ModeHydroBoundary::kAuto;
  ModeGravityBoundary gravity_boundary = ModeGravityBoundary::kAuto;
};

struct CompatibilityConfig {
  bool allow_unknown_keys = false;
};

struct SimulationConfig {
  int schema_version = 1;
  CosmologyConfig cosmology;
  NumericsConfig numerics;
  PhysicsConfig physics;
  OutputConfig output;
  ParallelConfig parallel;
  AnalysisConfig analysis;
  UnitsConfig units;
  ModeConfig mode;
  CompatibilityConfig compatibility;
};

struct ProvenanceMetadata {
  std::string source_name;
  std::uint64_t config_hash = 0;
  std::string config_hash_hex;
  std::vector<std::string> deprecation_warnings;
};

struct FrozenConfig {
  SimulationConfig config;
  std::string normalized_text;
  ProvenanceMetadata provenance;
};

struct ParseOptions {
  bool allow_unknown_keys = false;
};

class ConfigError : public std::runtime_error {
 public:
  explicit ConfigError(const std::string& message);
};

[[nodiscard]] FrozenConfig loadFrozenConfigFromFile(
    const std::filesystem::path& path,
    const ParseOptions& options = {});

[[nodiscard]] FrozenConfig loadFrozenConfigFromString(
    const std::string& config_text,
    const std::string& source_name,
    const ParseOptions& options = {});

void writeNormalizedConfigSnapshot(
    const FrozenConfig& frozen_config,
    const std::filesystem::path& run_directory);

[[nodiscard]] std::string modeToString(SimulationMode mode);
[[nodiscard]] std::string gravitySolverToString(GravitySolver solver);
[[nodiscard]] std::string hydroSolverToString(HydroSolver solver);
[[nodiscard]] std::string coordinateFrameToString(CoordinateFrame frame);
[[nodiscard]] std::string modeHydroBoundaryToString(ModeHydroBoundary boundary);
[[nodiscard]] std::string modeGravityBoundaryToString(ModeGravityBoundary boundary);
[[nodiscard]] std::string feedbackModeToString(FeedbackMode mode);
[[nodiscard]] std::string feedbackVariantToString(FeedbackVariant variant);
[[nodiscard]] std::string uvBackgroundModelToString(UvBackgroundModel model);
[[nodiscard]] std::string selfShieldingModelToString(SelfShieldingModel model);
[[nodiscard]] std::string coolingModelToString(CoolingModel model);

}  // namespace cosmosim::core
