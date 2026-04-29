#include <cassert>
#include <cmath>
#include <string>

#include "cosmosim/core/config.hpp"
#include "cosmosim/core/provenance.hpp"

namespace {

void testCommentsWhitespaceSectionsAndUnits() {
  const std::string config_text = R"(
# top-level comment
[units]
length_unit = mpc
coordinate_frame = comoving

[mode]
mode = zoom_in
ic_file = zoom_ics.hdf5
zoom_high_res_region = true
zoom_region_file = region_zoom.hdf5
zoom_region_radius = 1.0

[cosmology]
box_size = 50000 kpc
omega_matter = 0.31
omega_lambda = 0.69

[numerics]
time_begin_code = 0.0
time_end_code = 1.0
gravity_softening = 2.5 kpc

[output]
output_stem = snapshot
restart_stem = restart
)";

  const auto frozen = cosmosim::core::loadFrozenConfigFromString(config_text, "unit_case");
  assert(frozen.config.cosmology.box_size_mpc_comoving == 50.0);
  assert(frozen.config.numerics.gravity_softening_kpc_comoving == 2.5);
  assert(frozen.config.mode.zoom_high_res_region);
  assert(frozen.config.mode.zoom_region_file == "region_zoom.hdf5");
}

void testDuplicateKeyFails() {
  const std::string config_text = "[output]\noutput_stem = snapshot\noutput_stem = other\n";
  bool threw = false;
  try {
    (void)cosmosim::core::loadFrozenConfigFromString(config_text, "dup");
  } catch (const cosmosim::core::ConfigError&) {
    threw = true;
  }
  assert(threw);
}

void testMissingValueFails() {
  const std::string config_text = "mode.mode = zoom_in\noutput.output_stem =\n";
  bool threw = false;
  try {
    (void)cosmosim::core::loadFrozenConfigFromString(config_text, "missing_value");
  } catch (const cosmosim::core::ConfigError&) {
    threw = true;
  }
  assert(threw);
}

void testUnknownKeysFailUnlessCompatibilityEnabled() {
  const std::string config_text = "mode.mode = zoom_in\nfoo.bar = 1\n";
  bool threw = false;
  try {
    (void)cosmosim::core::loadFrozenConfigFromString(config_text, "unknown");
  } catch (const cosmosim::core::ConfigError&) {
    threw = true;
  }
  assert(threw);

  const std::string compat_text =
      "mode.mode = zoom_in\ncompatibility.allow_unknown_keys = true\nfoo.bar = 1\n";
  const auto frozen = cosmosim::core::loadFrozenConfigFromString(compat_text, "compat_unknown");
  assert(frozen.config.compatibility.allow_unknown_keys);
}

void testInvalidEnumFails() {
  const std::string config_text = "mode.mode = giant_universe\n";
  bool threw = false;
  try {
    (void)cosmosim::core::loadFrozenConfigFromString(config_text, "bad_enum");
  } catch (const cosmosim::core::ConfigError&) {
    threw = true;
  }
  assert(threw);
}

void testInvalidTypedPolicyEnumsFail() {
  const std::string bad_solver = "[mode]\nmode = zoom_in\n[numerics]\ngravity_solver = spectral_magic\n";
  bool threw = false;
  try {
    (void)cosmosim::core::loadFrozenConfigFromString(bad_solver, "bad_solver");
  } catch (const cosmosim::core::ConfigError&) {
    threw = true;
  }
  assert(threw);

  const std::string bad_frame = "[mode]\nmode = zoom_in\n[units]\ncoordinate_frame = weird\n";
  threw = false;
  try {
    (void)cosmosim::core::loadFrozenConfigFromString(bad_frame, "bad_frame");
  } catch (const cosmosim::core::ConfigError&) {
    threw = true;
  }
  assert(threw);
}

void testBoundaryModeValidation() {
  const std::string bad_config = "[mode]\nmode = cosmo_cube\nhydro_boundary = open\n";
  bool threw = false;
  try {
    (void)cosmosim::core::loadFrozenConfigFromString(bad_config, "bad_boundary");
  } catch (const cosmosim::core::ConfigError&) {
    threw = true;
  }
  assert(threw);

  const std::string isolated_ok =
      "[mode]\nmode = isolated_cluster\nhydro_boundary = reflective\ngravity_boundary = isolated_monopole\n";
  const auto frozen = cosmosim::core::loadFrozenConfigFromString(isolated_ok, "isolated_ok");
  assert(frozen.config.mode.hydro_boundary == cosmosim::core::ModeHydroBoundary::kReflective);
  assert(frozen.config.mode.gravity_boundary == cosmosim::core::ModeGravityBoundary::kIsolatedMonopole);
}

void testDefaultsCanonicalizationAndDeterminism() {
  const std::string first =
      "[mode]\nmode = cosmo_cube\n[output]\nrun_name = a\noutput_stem = snapshot\nrestart_stem = "
      "restart\n";
  const std::string second =
      "[output]\nrestart_stem = restart\noutput_stem = snapshot\nrun_name = a\n[mode]\nmode = "
      "cosmo_cube\n";

  const auto frozen_first = cosmosim::core::loadFrozenConfigFromString(first, "first");
  const auto frozen_second = cosmosim::core::loadFrozenConfigFromString(second, "second");

  assert(frozen_first.config.physics.enable_cooling);
  assert(frozen_first.config.physics.fb_mode == cosmosim::core::FeedbackMode::kThermalKineticMomentum);
  assert(frozen_first.config.parallel.deterministic_reduction);
  assert(frozen_first.normalized_text == frozen_second.normalized_text);
  assert(frozen_first.provenance.config_hash_hex == frozen_second.provenance.config_hash_hex);
  assert(frozen_first.provenance.config_hash_hex == cosmosim::core::stableConfigHashHex(frozen_first.normalized_text));
  const auto reparsed =
      cosmosim::core::loadFrozenConfigFromString(frozen_first.normalized_text, "normalized_roundtrip");
  assert(reparsed.normalized_text == frozen_first.normalized_text);
  assert(reparsed.provenance.config_hash_hex == frozen_first.provenance.config_hash_hex);
}

void testConfigNormalizationHashDeterminism() {
  const std::string first =
      "[mode]\nmode = zoom_in\n[output]\nrun_name = hash_a\n";
  const std::string second =
      "[output]\nrun_name = hash_a\n[mode]\nmode = zoom_in\n";
  const auto frozen_first = cosmosim::core::loadFrozenConfigFromString(first, "hash_first");
  const auto frozen_second = cosmosim::core::loadFrozenConfigFromString(second, "hash_second");
  assert(frozen_first.normalized_text == frozen_second.normalized_text);
  assert(frozen_first.provenance.config_hash_hex == frozen_second.provenance.config_hash_hex);
  assert(frozen_first.provenance.config_hash_hex == cosmosim::core::stableConfigHashHex(frozen_first.normalized_text));
}

void testConfigRuntimeOwnership() {
  const std::string text = R"(
[mode]
mode = zoom_in
[cosmology]
box_size = 64 mpc
[numerics]
treepm_pm_grid = 32
treepm_asmth_cells = 1.25
treepm_rcut_cells = 4.5
[output]
run_name = ownership
)";
  const auto frozen = cosmosim::core::loadFrozenConfigFromString(text, "runtime_ownership");
  assert(frozen.config.numerics.treepm_pm_grid_nx == 32);
  assert(frozen.config.numerics.treepm_pm_grid_ny == 32);
  assert(frozen.config.numerics.treepm_pm_grid_nz == 32);
  assert(frozen.normalized_text.find("treepm_pm_grid =") == std::string::npos);
  assert(frozen.normalized_text.find("treepm_pm_grid_nx = 32") != std::string::npos);
  const auto reparsed = cosmosim::core::loadFrozenConfigFromString(frozen.normalized_text, "runtime_ownership_reparse");
  assert(reparsed.config.numerics.treepm_pm_grid_nx == frozen.config.numerics.treepm_pm_grid_nx);
  assert(reparsed.config.numerics.treepm_asmth_cells == frozen.config.numerics.treepm_asmth_cells);
  assert(reparsed.config.numerics.treepm_rcut_cells == frozen.config.numerics.treepm_rcut_cells);
  assert(reparsed.provenance.config_hash_hex == frozen.provenance.config_hash_hex);
}

void testDerivedRuntimeConfigOwnership() {
  const std::string text = R"(
[mode]
mode = zoom_in
[cosmology]
box_size_x = 70 mpc
box_size_y = 50 mpc
box_size_z = 30 mpc
[numerics]
time_begin_code = 0.125
time_end_code = 0.875
gravity_softening = 1.0 kpc
gravity_softening_gas = 0.5 kpc
gravity_softening_black_hole = 0.05 kpc
treepm_pm_grid_nx = 28
treepm_pm_grid_ny = 20
treepm_pm_grid_nz = 12
)";
  const auto frozen = cosmosim::core::loadFrozenConfigFromString(text, "derived_runtime");
  const auto derived = cosmosim::core::deriveRuntimeConfig(frozen);
  assert(derived.time_begin_code == 0.125);
  assert(derived.time_end_code == 0.875);
  assert(derived.box_size_mpc_comoving[0] == 70.0);
  assert(derived.box_size_mpc_comoving[1] == 50.0);
  assert(derived.box_size_mpc_comoving[2] == 30.0);
  assert(derived.treepm_pm_grid_shape[0] == 28);
  assert(derived.treepm_pm_grid_shape[1] == 20);
  assert(derived.treepm_pm_grid_shape[2] == 12);
  assert(derived.gravity_softening_kpc_comoving_by_species[0] == 1.0);
  assert(derived.gravity_softening_kpc_comoving_by_species[1] == 0.5);
  assert(derived.gravity_softening_kpc_comoving_by_species[2] == 1.0);
  assert(derived.gravity_softening_kpc_comoving_by_species[3] == 0.05);
  assert(derived.gravity_softening_kpc_comoving_by_species[4] == 1.0);
  assert(derived.normalized_config_hash == frozen.provenance.config_hash);
  assert(derived.normalized_config_hash_hex == frozen.provenance.config_hash_hex);

  const std::string bad_time = "[mode]\nmode = zoom_in\n[numerics]\ntime_begin_code = 1.0\ntime_end_code = 0.5\n";
  bool threw = false;
  try {
    const auto bad = cosmosim::core::loadFrozenConfigFromString(bad_time, "derived_runtime_bad_time");
    (void)cosmosim::core::deriveRuntimeConfig(bad);
  } catch (const cosmosim::core::ConfigError&) {
    threw = true;
  }
  assert(threw);
}

void testUrlsWindowsPathsAndQuotedHashesAreNotTruncated() {
  const std::string config_text = R"(
[mode]
mode = zoom_in
[physics]
metal_line_table_path = https://example.com/tables/metal_line_rates.h5 // fetch remotely later
stellar_evolution_table_path = "C://stellar_tables//ssp#grid.h5" # local override
)";

  const auto frozen = cosmosim::core::loadFrozenConfigFromString(config_text, "path_preservation");
  assert(frozen.config.physics.metal_line_table_path == "https://example.com/tables/metal_line_rates.h5");
  assert(frozen.config.physics.stellar_evolution_table_path == "\"C://stellar_tables//ssp#grid.h5\"");
}

void testDeprecatedAliasesAndCanonicalCollision() {
  const std::string alias_only = "mode = isolated_galaxy\n";
  const auto frozen = cosmosim::core::loadFrozenConfigFromString(alias_only, "alias_only");
  assert(frozen.config.mode.mode == cosmosim::core::SimulationMode::kIsolatedGalaxy);
  assert(!frozen.provenance.deprecation_warnings.empty());

  const std::string conflicting = "mode = zoom_in\nmode.mode = cosmo_cube\n";
  bool threw = false;
  try {
    (void)cosmosim::core::loadFrozenConfigFromString(conflicting, "alias_conflict");
  } catch (const cosmosim::core::ConfigError&) {
    threw = true;
  }
  assert(threw);
}

void testFeedbackConfigKeysAndValidation() {
  const std::string good_text = R"(
[mode]
mode = zoom_in
[physics]
fb_mode = momentum
fb_variant = stochastic
fb_stochastic_event_probability = 1.0
fb_neighbor_count = 4
)";
  const auto frozen = cosmosim::core::loadFrozenConfigFromString(good_text, "feedback_good");
  assert(frozen.config.physics.fb_mode == cosmosim::core::FeedbackMode::kMomentum);
  assert(frozen.config.physics.fb_variant == cosmosim::core::FeedbackVariant::kStochastic);
  assert(frozen.config.physics.fb_neighbor_count == 4);

  const std::string bad_text = R"(
[mode]
mode = zoom_in
[physics]
fb_mode = hidden_magic
)";
  bool threw = false;
  try {
    (void)cosmosim::core::loadFrozenConfigFromString(bad_text, "feedback_bad");
  } catch (const cosmosim::core::ConfigError&) {
    threw = true;
  }
  assert(threw);
}

void testCoolingPolicyEnumsAndValidation() {
  const std::string good_text = R"(
[mode]
mode = zoom_in
[physics]
uv_background_model = fg20
self_shielding_model = rahmati13_like
cooling_model = primordial_metal_line
)";
  const auto frozen = cosmosim::core::loadFrozenConfigFromString(good_text, "cooling_good");
  assert(frozen.config.physics.uv_background_model == cosmosim::core::UvBackgroundModel::kFg20);
  assert(frozen.config.physics.self_shielding_model == cosmosim::core::SelfShieldingModel::kRahmati13Like);
  assert(frozen.config.physics.cooling_model == cosmosim::core::CoolingModel::kPrimordialMetalLine);

  const std::string bad_text = R"(
[mode]
mode = zoom_in
[physics]
uv_background_model = cosmic_soup
)";
  bool threw = false;
  try {
    (void)cosmosim::core::loadFrozenConfigFromString(bad_text, "cooling_bad");
  } catch (const cosmosim::core::ConfigError&) {
    threw = true;
  }
  assert(threw);
}

void testDiagnosticsExecutionPolicyValidation() {
  const std::string good_text = R"(
[mode]
mode = zoom_in
[analysis]
diagnostics_execution_policy = all_including_provisional
)";
  const auto frozen = cosmosim::core::loadFrozenConfigFromString(good_text, "diag_policy_good");
  assert(
      frozen.config.analysis.diagnostics_execution_policy ==
      cosmosim::core::AnalysisConfig::DiagnosticsExecutionPolicy::kAllIncludingProvisional);
  assert(
      frozen.normalized_text.find("diagnostics_execution_policy = all_including_provisional") !=
      std::string::npos);
  assert(frozen.normalized_text.find("diagnostics_execution_policy = unknown") == std::string::npos);

  const std::string bad_text = R"(
[mode]
mode = zoom_in
[analysis]
diagnostics_execution_policy = unsupported_policy
)";
  bool threw = false;
  try {
    (void)cosmosim::core::loadFrozenConfigFromString(bad_text, "diag_policy_bad");
  } catch (const cosmosim::core::ConfigError&) {
    threw = true;
  }
  assert(threw);
}

void testEnumSerializationIsFailFastWithoutUnknownFallback() {
  bool threw = false;
  try {
    (void)cosmosim::core::modeToString(static_cast<cosmosim::core::SimulationMode>(255));
  } catch (const cosmosim::core::ConfigError&) {
    threw = true;
  }
  assert(threw);
}

void testSpeciesSofteningOverridesRoundtrip() {
  const std::string text = R"(
[mode]
mode = zoom_in
[numerics]
gravity_softening = 1.0 kpc
gravity_softening_gas = 0.5 kpc
gravity_softening_dark_matter = 2.0 kpc
gravity_softening_star = 0.2 kpc
gravity_softening_black_hole = 0.05 kpc
gravity_softening_tracer = 1.5 kpc
)";
  const auto frozen = cosmosim::core::loadFrozenConfigFromString(text, "species_softening");
  assert(std::abs(frozen.config.numerics.gravity_softening_gas_kpc_comoving - 0.5) < 1.0e-12);
  assert(std::abs(frozen.config.numerics.gravity_softening_dark_matter_kpc_comoving - 2.0) < 1.0e-12);
  assert(std::abs(frozen.config.numerics.gravity_softening_star_kpc_comoving - 0.2) < 1.0e-12);
  assert(std::abs(frozen.config.numerics.gravity_softening_black_hole_kpc_comoving - 0.05) < 1.0e-12);
  assert(std::abs(frozen.config.numerics.gravity_softening_tracer_kpc_comoving - 1.5) < 1.0e-12);
  assert(frozen.normalized_text.find("gravity_softening_gas = ") != std::string::npos);
  const auto reparsed = cosmosim::core::loadFrozenConfigFromString(frozen.normalized_text, "species_softening_roundtrip");
  assert(std::abs(reparsed.config.numerics.gravity_softening_black_hole_kpc_comoving - 0.05) < 1.0e-12);
}

void testTreePmNumericsRoundtripAndValidation() {
  const std::string good_text = R"(
[mode]
mode = zoom_in
[cosmology]
box_size = 64 mpc
[numerics]
treepm_pm_grid = 32
treepm_asmth_cells = 1.5
treepm_rcut_cells = 5.0
treepm_assignment_scheme = cic
treepm_enable_window_deconvolution = true
treepm_update_cadence_steps = 3
treepm_pm_decomposition_mode = slab
treepm_tree_exchange_batch_bytes = 8388608
)";
  const auto frozen = cosmosim::core::loadFrozenConfigFromString(good_text, "treepm_good");
  assert(frozen.config.numerics.treepm_pm_grid_nx == 32);
  assert(frozen.config.numerics.treepm_pm_grid_ny == 32);
  assert(frozen.config.numerics.treepm_pm_grid_nz == 32);
  assert(frozen.config.numerics.treepm_asmth_cells == 1.5);
  assert(frozen.config.numerics.treepm_rcut_cells == 5.0);
  assert(
      frozen.config.numerics.treepm_assignment_scheme ==
      cosmosim::core::TreePmAssignmentScheme::kCic);
  assert(frozen.config.numerics.treepm_enable_window_deconvolution);
  assert(frozen.config.numerics.treepm_update_cadence_steps == 3);
  assert(
      frozen.config.numerics.treepm_pm_decomposition_mode ==
      cosmosim::core::PmDecompositionMode::kSlab);
  assert(frozen.config.numerics.treepm_tree_exchange_batch_bytes == 8388608ULL);
  assert(frozen.normalized_text.find("treepm_pm_grid_nx = 32") != std::string::npos);
  assert(frozen.normalized_text.find("treepm_assignment_scheme = cic") != std::string::npos);
  assert(frozen.normalized_text.find("treepm_pm_decomposition_mode = slab") != std::string::npos);
  const auto reparsed =
      cosmosim::core::loadFrozenConfigFromString(frozen.normalized_text, "treepm_roundtrip");
  assert(reparsed.config.numerics.treepm_pm_grid_nx == 32);
  assert(reparsed.config.numerics.treepm_update_cadence_steps == 3);
  assert(reparsed.config.numerics.treepm_tree_exchange_batch_bytes == 8388608ULL);

  const std::string bad_grid =
      "[mode]\nmode = zoom_in\n[numerics]\ntreepm_pm_grid = 0\n";
  bool threw = false;
  try {
    (void)cosmosim::core::loadFrozenConfigFromString(bad_grid, "treepm_bad_grid");
  } catch (const cosmosim::core::ConfigError&) {
    threw = true;
  }
  assert(threw);

  const std::string bad_scheme =
      "[mode]\nmode = zoom_in\n[numerics]\ntreepm_assignment_scheme = magic\n";
  threw = false;
  try {
    (void)cosmosim::core::loadFrozenConfigFromString(bad_scheme, "treepm_bad_scheme");
  } catch (const cosmosim::core::ConfigError&) {
    threw = true;
  }
  assert(threw);

  const std::string pencil_decomposition =
      "[mode]\nmode = zoom_in\n[numerics]\ntreepm_pm_decomposition_mode = pencil\n";
  const auto parsed_pencil = cosmosim::core::loadFrozenConfigFromString(
      pencil_decomposition, "treepm_pencil_decomp");
  assert(parsed_pencil.config.numerics.treepm_pm_decomposition_mode ==
      cosmosim::core::PmDecompositionMode::kPencil);

  const std::string bad_batch_bytes =
      "[mode]\nmode = zoom_in\n[numerics]\ntreepm_tree_exchange_batch_bytes = 0\n";
  threw = false;
  try {
    (void)cosmosim::core::loadFrozenConfigFromString(bad_batch_bytes, "treepm_bad_batch");
  } catch (const cosmosim::core::ConfigError&) {
    threw = true;
  }
  assert(threw);
}

void testAxisAwareBoxAndPmGridCanonicalizationCompatibility() {
  const std::string scalar_text = R"(
[mode]
mode = zoom_in
[cosmology]
box_size = 64 mpc
[numerics]
treepm_pm_grid = 24
)";
  const auto scalar_frozen = cosmosim::core::loadFrozenConfigFromString(scalar_text, "scalar_compat");
  assert(scalar_frozen.config.cosmology.box_size_x_mpc_comoving == 64.0);
  assert(scalar_frozen.config.cosmology.box_size_y_mpc_comoving == 64.0);
  assert(scalar_frozen.config.cosmology.box_size_z_mpc_comoving == 64.0);
  assert(scalar_frozen.config.numerics.treepm_pm_grid_nx == 24);
  assert(scalar_frozen.config.numerics.treepm_pm_grid_ny == 24);
  assert(scalar_frozen.config.numerics.treepm_pm_grid_nz == 24);
  assert(scalar_frozen.normalized_text.find("box_size_x = 64 mpc") != std::string::npos);
  assert(scalar_frozen.normalized_text.find("treepm_pm_grid_nx = 24") != std::string::npos);
  assert(scalar_frozen.normalized_text.find("treepm_pm_grid =") == std::string::npos);

  const std::string axis_text = R"(
[mode]
mode = zoom_in
[cosmology]
box_size_x = 70 mpc
box_size_y = 50 mpc
box_size_z = 30 mpc
[numerics]
treepm_pm_grid_nx = 28
treepm_pm_grid_ny = 20
treepm_pm_grid_nz = 12
)";
  const auto axis_frozen = cosmosim::core::loadFrozenConfigFromString(axis_text, "axis_canonical");
  assert(axis_frozen.config.cosmology.box_size_x_mpc_comoving == 70.0);
  assert(axis_frozen.config.cosmology.box_size_y_mpc_comoving == 50.0);
  assert(axis_frozen.config.cosmology.box_size_z_mpc_comoving == 30.0);
  assert(axis_frozen.config.numerics.treepm_pm_grid_nx == 28);
  assert(axis_frozen.config.numerics.treepm_pm_grid_ny == 20);
  assert(axis_frozen.config.numerics.treepm_pm_grid_nz == 12);
}

void testBlackHoleAgnConfigKeysAndValidation() {
  const std::string good_text = R"(
[mode]
mode = zoom_in
[physics]
enable_black_hole_agn = true
bh_seed_halo_mass_threshold_code = 500
bh_seed_mass_code = 2.5
bh_seed_max_per_cell = 2
bh_alpha_bondi = 3.0
bh_use_eddington_cap = true
bh_epsilon_r = 0.15
bh_epsilon_f = 0.2
bh_feedback_coupling_efficiency = 0.8
bh_duty_cycle_active_edd_ratio_threshold = 0.03
)";
  const auto frozen = cosmosim::core::loadFrozenConfigFromString(good_text, "bh_good");
  assert(frozen.config.physics.enable_black_hole_agn);
  assert(frozen.config.physics.bh_seed_halo_mass_threshold_code == 500.0);
  assert(frozen.config.physics.bh_seed_max_per_cell == 2);
  assert(frozen.config.physics.bh_alpha_bondi == 3.0);

  const std::string bad_text =
      "[mode]\nmode = zoom_in\n[physics]\nbh_feedback_coupling_efficiency = 1.3\n";
  bool threw = false;
  try {
    (void)cosmosim::core::loadFrozenConfigFromString(bad_text, "bh_bad");
  } catch (const cosmosim::core::ConfigError&) {
    threw = true;
  }
  assert(threw);
}

void testTracerConfigKeysAndValidation() {
  const std::string good_text = R"(
[mode]
mode = zoom_in
[physics]
enable_tracers = true
tracer_track_mass = true
tracer_min_host_mass_code = 1.0e-6
)";
  const auto frozen = cosmosim::core::loadFrozenConfigFromString(good_text, "tracer_good");
  assert(frozen.config.physics.enable_tracers);
  assert(frozen.config.physics.tracer_track_mass);
  assert(frozen.config.physics.tracer_min_host_mass_code == 1.0e-6);

  const std::string bad_text = "[mode]\nmode = zoom_in\n[physics]\ntracer_min_host_mass_code = -1.0\n";
  bool threw = false;
  try {
    (void)cosmosim::core::loadFrozenConfigFromString(bad_text, "tracer_bad");
  } catch (const cosmosim::core::ConfigError&) {
    threw = true;
  }
  assert(threw);
}

void testZoomLongRangeStrategyValidationAndRoundtrip() {
  const std::string good_text = R"(
[mode]
mode = zoom_in
zoom_high_res_region = true
zoom_region_file = zoom_region.hdf5
zoom_long_range_strategy = global_coarse_plus_focused_highres_correction
zoom_region_center_x = 12.5
zoom_region_center_y = 13.5
zoom_region_center_z = 14.5
zoom_region_radius = 2.0
zoom_focused_pm_grid_nx = 32
zoom_focused_pm_grid_ny = 48
zoom_focused_pm_grid_nz = 64
zoom_contamination_radius = 3.0
)";
  const auto frozen = cosmosim::core::loadFrozenConfigFromString(good_text, "zoom_long_range_good");
  assert(frozen.config.mode.zoom_long_range_strategy ==
      cosmosim::core::ZoomLongRangeStrategy::kGlobalCoarsePlusFocusedHighResCorrection);
  assert(frozen.config.mode.zoom_region_center_x_mpc_comoving == 12.5);
  assert(frozen.config.mode.zoom_focused_pm_grid_nz == 64);
  assert(frozen.normalized_text.find("zoom_long_range_strategy = global_coarse_plus_focused_highres_correction") !=
      std::string::npos);

  const std::string bad_text = R"(
[mode]
mode = zoom_in
zoom_high_res_region = true
zoom_region_file = zoom_region.hdf5
zoom_long_range_strategy = global_coarse_plus_focused_highres_correction
zoom_region_radius = 2.0
zoom_focused_pm_grid_nx = 0
zoom_focused_pm_grid_ny = 16
zoom_focused_pm_grid_nz = 16
)";
  bool threw = false;
  try {
    (void)cosmosim::core::loadFrozenConfigFromString(bad_text, "zoom_long_range_bad");
  } catch (const cosmosim::core::ConfigError&) {
    threw = true;
  }
  assert(threw);
}

}  // namespace

int main() {
  testCommentsWhitespaceSectionsAndUnits();
  testDuplicateKeyFails();
  testMissingValueFails();
  testUnknownKeysFailUnlessCompatibilityEnabled();
  testInvalidEnumFails();
  testInvalidTypedPolicyEnumsFail();
  testBoundaryModeValidation();
  testDefaultsCanonicalizationAndDeterminism();
  testConfigNormalizationHashDeterminism();
  testConfigRuntimeOwnership();
  testDerivedRuntimeConfigOwnership();
  testUrlsWindowsPathsAndQuotedHashesAreNotTruncated();
  testDeprecatedAliasesAndCanonicalCollision();
  testFeedbackConfigKeysAndValidation();
  testCoolingPolicyEnumsAndValidation();
  testDiagnosticsExecutionPolicyValidation();
  testEnumSerializationIsFailFastWithoutUnknownFallback();
  testSpeciesSofteningOverridesRoundtrip();
  testTreePmNumericsRoundtripAndValidation();
  testAxisAwareBoxAndPmGridCanonicalizationCompatibility();
  testBlackHoleAgnConfigKeysAndValidation();
  testTracerConfigKeysAndValidation();
  testZoomLongRangeStrategyValidationAndRoundtrip();
  return 0;
}
