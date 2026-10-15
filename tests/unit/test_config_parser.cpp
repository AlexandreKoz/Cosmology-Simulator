// 20261012 regression repairs: NOT RUN — deferred to maintainer.
#include "../support/test_require.hpp"
#include <cmath>
#include <string>
#include <type_traits>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/core/config.hpp"
#include "cosmosim/core/cosmology.hpp"
#include "cosmosim/core/units.hpp"
#include "cosmosim/core/provenance.hpp"

namespace {

// Mutate one canonical numerics value after admitting the complete baseline.
// Missing mode/IC fields must never masquerade as TreePM validation coverage.
std::string withNumericsValue(
    const std::string& normalized, const std::string& key, const std::string& value) {
  auto text = normalized;
  const std::string prefix = "\n" + key + " = ";
  const auto position = text.find(prefix);
  CHUI_REQUIRE(position != std::string::npos);
  CHUI_REQUIRE(text.find(prefix, position + prefix.size()) == std::string::npos);
  const auto begin = position + prefix.size();
  const auto end = text.find('\n', begin);
  CHUI_REQUIRE(end != std::string::npos);
  text.replace(begin, end - begin, value);
  return text;
}

void requireNumericsRejected(
    const std::string& baseline, const std::string& key, const std::string& value,
    const std::string& context = {}) {
  const auto valid = cosmosim::core::loadFrozenConfigFromString(baseline, "treepm_valid_control");
  CHUI_REQUIRE(valid.normalized_text == baseline);
  bool rejected = false;
  try {
    (void)cosmosim::core::loadFrozenConfigFromString(
        withNumericsValue(baseline, key, value), "treepm_invalid_field");
  } catch (const cosmosim::core::ConfigError& error) {
    const std::string message = error.what();
    rejected = message.find(context.empty() ? "numerics." + key : context) != std::string::npos;
  }
  CHUI_REQUIRE(rejected);
}

static_assert(!std::is_default_constructible_v<cosmosim::core::SimulationConfig>,
              "SimulationConfig must not be publicly default-constructible; use FrozenConfig normalization");
static_assert(std::is_copy_constructible_v<cosmosim::core::SimulationConfig>);

void testCommentsWhitespaceSectionsAndUnits() {
  const std::string config_text = R"(
# top-level comment
[units]
length_unit = mpc
coordinate_frame = comoving

[mode]
mode = zoom_in
ic_file = zoom_ics.hdf5
ic_convention = gadget_arepo_bridge_v1
ic_bridge_source_length_unit_to_si = 3.0856775814913673e19
ic_bridge_source_mass_unit_to_si = 1.98847e30
ic_bridge_source_velocity_unit_to_si = 1000
ic_bridge_coordinate_frame = comoving
ic_bridge_velocity_convention = physical_peculiar
ic_bridge_length_hubble_exponent = 0
ic_bridge_length_scale_factor_exponent = 0
ic_bridge_mass_hubble_exponent = 0
ic_bridge_mass_scale_factor_exponent = 0
ic_bridge_velocity_hubble_exponent = 0
ic_bridge_velocity_scale_factor_exponent = 0
zoom_high_res_region = true
zoom_region_file = region_zoom.hdf5
zoom_region_radius = 1.0

[cosmology]
box_size = 50000 kpc
omega_matter = 0.31
omega_lambda = 0.69

[numerics]
t_code_begin = 0.0
t_code_end = 1.0
gravity_softening = 2.5 kpc

[output]
output_stem = snapshot
restart_stem = restart
)";

  const auto frozen = cosmosim::core::loadFrozenConfigFromString(config_text, "unit_case");
  CHUI_REQUIRE(frozen.raw_text == config_text);
  CHUI_REQUIRE(frozen.user_config.source_name == "unit_case");
  CHUI_REQUIRE(!frozen.user_config.entries.empty());
  CHUI_REQUIRE(frozen.config.cosmology.box_size_mpc_comoving == 50.0);
  CHUI_REQUIRE(frozen.config.numerics.gravity_softening_kpc_comoving == 2.5);
  CHUI_REQUIRE(frozen.config.mode.zoom_high_res_region);
  CHUI_REQUIRE(frozen.config.mode.zoom_region_file == "region_zoom.hdf5");

  std::string malformed = config_text;
  const std::string valid_value = "gravity_softening = 2.5 kpc";
  const auto offset = malformed.find(valid_value);
  CHUI_REQUIRE(offset != std::string::npos);
  malformed.replace(offset, valid_value.size(), valid_value + " trailing_garbage");
  bool trailing_threw = false;
  try {
    (void)cosmosim::core::loadFrozenConfigFromString(malformed, "trailing_tokens");
  } catch (const cosmosim::core::ConfigError&) {
    trailing_threw = true;
  }
  CHUI_REQUIRE(trailing_threw);
}

void testDuplicateKeyFails() {
  const std::string config_text = "[output]\noutput_stem = snapshot\noutput_stem = other\n";
  bool threw = false;
  try {
    (void)cosmosim::core::loadFrozenConfigFromString(config_text, "dup");
  } catch (const cosmosim::core::ConfigError&) {
    threw = true;
  }
  CHUI_REQUIRE(threw);
}

void testRequiredSectionAndFieldFailures() {
  bool threw = false;
  try {
    (void)cosmosim::core::loadFrozenConfigFromString("[numerics]\nt_code_begin=0\nt_code_end=1\n", "missing_mode");
  } catch (const cosmosim::core::ConfigError&) {
    threw = true;
  }
  CHUI_REQUIRE(threw);

  threw = false;
  try {
    (void)cosmosim::core::loadFrozenConfigFromString(
        "[mode]\nic_file=ics.hdf5\n[numerics]\nt_code_begin=0\nt_code_end=1\n[output]\noutput_stem=s\n",
        "missing_required_field");
  } catch (const cosmosim::core::ConfigError&) {
    threw = true;
  }
  CHUI_REQUIRE(threw);
}

void testMissingValueFails() {
  const std::string config_text = "mode.mode = zoom_in\noutput.output_stem =\n";
  bool threw = false;
  try {
    (void)cosmosim::core::loadFrozenConfigFromString(config_text, "missing_value");
  } catch (const cosmosim::core::ConfigError&) {
    threw = true;
  }
  CHUI_REQUIRE(threw);
}

void testUnknownKeysFailUnlessCompatibilityEnabled() {
  const std::string config_text = "mode.mode = zoom_in\nfoo.bar = 1\n";
  bool threw = false;
  try {
    (void)cosmosim::core::loadFrozenConfigFromString(config_text, "unknown");
  } catch (const cosmosim::core::ConfigError&) {
    threw = true;
  }
  CHUI_REQUIRE(threw);

  const std::string compat_text =
      "mode.mode = zoom_in\ncompatibility.allow_unknown_keys = true\nfoo.bar = 1\n";
  const auto frozen = cosmosim::core::loadFrozenConfigFromString(compat_text, "compat_unknown");
  CHUI_REQUIRE(frozen.config.compatibility.allow_unknown_keys);
}

void testAliasConflictFails() {
  const std::string config_text =
      "mode.mode = zoom_in\nnumerics.initial_redshift = 2.0\nnumerics.z_begin = 3.0\n";
  bool threw = false;
  try {
    (void)cosmosim::core::loadFrozenConfigFromString(config_text, "alias_conflict");
  } catch (const cosmosim::core::ConfigError&) {
    threw = true;
  }
  CHUI_REQUIRE(threw);
}

void testInvalidEnumFails() {
  const std::string config_text = "mode.mode = giant_universe\n";
  bool threw = false;
  try {
    (void)cosmosim::core::loadFrozenConfigFromString(config_text, "bad_enum");
  } catch (const cosmosim::core::ConfigError&) {
    threw = true;
  }
  CHUI_REQUIRE(threw);
}


void testTypedInitialConditionConfiguration() {
  const auto bridge = cosmosim::core::loadFrozenConfigFromString(
      "[mode]\n"
      "mode = zoom_in\n"
      "ic_file = input.0.hdf5\n"
      "ic_convention = gadget_arepo_bridge_v1\n"
      "ic_bridge_source_length_unit_to_si = 3.0856775814913673e19\n"
      "ic_bridge_source_mass_unit_to_si = 1.98847e30\n"
      "ic_bridge_source_velocity_unit_to_si = 1000\n"
      "ic_bridge_coordinate_frame = comoving\n"
      "ic_bridge_velocity_convention = gadget_arepo_stored_peculiar\n"
      "ic_bridge_length_hubble_exponent = 0\n"
      "ic_bridge_length_scale_factor_exponent = 0\n"
      "ic_bridge_mass_hubble_exponent = 0\n"
      "ic_bridge_mass_scale_factor_exponent = 0\n"
      "ic_bridge_velocity_hubble_exponent = 0\n"
      "ic_bridge_velocity_scale_factor_exponent = 0\n"
      "ic_chunk_particle_count = 4096\n"
      "ic_staging_particle_count = 2048\n"
      "ic_part_type2_policy = dark_matter\n"
      "ic_part_type3_policy = reject\n"
      "ic_gas_internal_energy_policy = reject\n"
      "ic_gas_density_policy = use_config_value\n"
      "ic_gas_density_value_code = 3.25\n"
      "ic_star_formation_time_policy = dialect_defined_default\n"
      "ic_star_initial_mass_policy = reject\n"
      "ic_star_metallicity_policy = reject\n"
      "ic_bh_mdot_policy = dialect_defined_default\n",
      "typed_ic_bridge");
  CHUI_REQUIRE(
      bridge.config.mode.ic_convention ==
      cosmosim::core::InitialConditionConvention::kGadgetArepoBridgeV1);
  CHUI_REQUIRE(bridge.config.mode.ic_chunk_particle_count == 4096U);
  CHUI_REQUIRE(bridge.config.mode.ic_staging_particle_count == 2048U);
  CHUI_REQUIRE(
      bridge.config.mode.ic_part_type2_policy ==
      cosmosim::core::InitialConditionSpeciesPolicy::kDarkMatter);
  CHUI_REQUIRE(
      bridge.config.mode.ic_gas_density_policy ==
      cosmosim::core::InitialConditionMissingFieldPolicy::kUseConfigValue);
  CHUI_REQUIRE(bridge.config.mode.ic_gas_density_value_code == 3.25);
  CHUI_REQUIRE(
      bridge.config.mode.ic_bh_mdot_policy ==
      cosmosim::core::InitialConditionMissingFieldPolicy::kDialectDefinedDefault);
  CHUI_REQUIRE(
      bridge.normalized_text.find(
          "ic_convention = gadget_arepo_bridge_v1") != std::string::npos);
  CHUI_REQUIRE(
      bridge.normalized_text.find("ic_chunk_particle_count = 4096") !=
      std::string::npos);
  CHUI_REQUIRE(bridge.normalized_text.find(
             "ic_bridge_coordinate_frame = comoving") != std::string::npos);
  CHUI_REQUIRE(bridge.config.mode.ic_bridge_velocity_convention ==
         cosmosim::core::InitialConditionVelocityConvention::kGadgetArepoStoredPeculiar);
  CHUI_REQUIRE(bridge.normalized_text.find(
             "ic_bridge_velocity_convention = gadget_arepo_stored_peculiar") !=
         std::string::npos);

  bool incomplete_bridge_rejected = false;
  try {
    (void)cosmosim::core::loadFrozenConfigFromString(
        "[mode]\nmode = zoom_in\nic_file = input.hdf5\n"
        "ic_convention = gadget_arepo_bridge_v1\n",
        "incomplete_bridge");
  } catch (const cosmosim::core::ConfigError&) {
    incomplete_bridge_rejected = true;
  }
  CHUI_REQUIRE(incomplete_bridge_rejected);

  bool rejected = false;
  try {
    (void)cosmosim::core::loadFrozenConfigFromString(
        "[mode]\nmode = zoom_in\nic_file = ambiguous.hdf5\n",
        "implicit_external_ic");
  } catch (const cosmosim::core::ConfigError&) {
    rejected = true;
  }
  CHUI_REQUIRE(rejected);

  rejected = false;
  try {
    (void)cosmosim::core::loadFrozenConfigFromString(
        "[mode]\nmode = zoom_in\nic_file = input.hdf5\n"
        "ic_convention = mystery_v9\n",
        "unknown_ic_convention");
  } catch (const cosmosim::core::ConfigError&) {
    rejected = true;
  }
  CHUI_REQUIRE(rejected);

  rejected = false;
  try {
    (void)cosmosim::core::loadFrozenConfigFromString(
        "[mode]\nmode = zoom_in\nic_file = input.hdf5\n"
        "ic_convention = manifest_v1\n",
        "manifest_without_path");
  } catch (const cosmosim::core::ConfigError&) {
    rejected = true;
  }
  CHUI_REQUIRE(rejected);

  const auto manifest = cosmosim::core::loadFrozenConfigFromString(
      "[mode]\nmode = zoom_in\n"
      "ic_convention = manifest_v1\n"
      "ic_manifest_file = input.audit.json\n",
      "manifest_ic");
  CHUI_REQUIRE(
      manifest.config.mode.ic_convention ==
      cosmosim::core::InitialConditionConvention::kManifestV1);
  CHUI_REQUIRE(manifest.config.mode.ic_manifest_file == "input.audit.json");
  CHUI_REQUIRE(
      manifest.normalized_text.find("ic_file = generated") !=
      std::string::npos);

  bool unsupported_reconstruct_rejected = false;
  try {
    (void)cosmosim::core::loadFrozenConfigFromString(
        "[mode]\nmode = zoom_in\nic_file = input.hdf5\n"
        "ic_convention = gadget_arepo_bridge_v1\n"
        "ic_bridge_source_length_unit_to_si = 1\n"
        "ic_bridge_source_mass_unit_to_si = 1\n"
        "ic_bridge_source_velocity_unit_to_si = 1\n"
        "ic_bridge_coordinate_frame = comoving\n"
        "ic_bridge_velocity_convention = physical_peculiar\n"
        "ic_bridge_length_hubble_exponent = 0\n"
        "ic_bridge_length_scale_factor_exponent = 0\n"
        "ic_bridge_mass_hubble_exponent = 0\n"
        "ic_bridge_mass_scale_factor_exponent = 0\n"
        "ic_bridge_velocity_hubble_exponent = 0\n"
        "ic_bridge_velocity_scale_factor_exponent = 0\n"
        "ic_gas_density_policy = reconstruct\n",
        "unsupported_ic_reconstruct");
  } catch (const cosmosim::core::ConfigError&) {
    unsupported_reconstruct_rejected = true;
  }
  CHUI_REQUIRE(unsupported_reconstruct_rejected);
}

void testInvalidTypedPolicyEnumsFail() {
  const std::string bad_solver = "[mode]\nmode = zoom_in\n[numerics]\ngravity_solver = spectral_magic\n";
  bool threw = false;
  try {
    (void)cosmosim::core::loadFrozenConfigFromString(bad_solver, "bad_solver");
  } catch (const cosmosim::core::ConfigError&) {
    threw = true;
  }
  CHUI_REQUIRE(threw);

  const std::string bad_frame = "[mode]\nmode = zoom_in\n[units]\ncoordinate_frame = weird\n";
  threw = false;
  try {
    (void)cosmosim::core::loadFrozenConfigFromString(bad_frame, "bad_frame");
  } catch (const cosmosim::core::ConfigError&) {
    threw = true;
  }
  CHUI_REQUIRE(threw);
}

void testBoundaryModeValidation() {
  const std::string bad_config = "[mode]\nmode = cosmo_cube\nhydro_boundary = open\n";
  bool threw = false;
  try {
    (void)cosmosim::core::loadFrozenConfigFromString(bad_config, "bad_boundary");
  } catch (const cosmosim::core::ConfigError&) {
    threw = true;
  }
  CHUI_REQUIRE(threw);

  const std::string isolated_ok =
      "[mode]\nmode = isolated_cluster\nhydro_boundary = reflective\ngravity_boundary = isolated_monopole\n"
      "[numerics]\ntreepm_enable_window_deconvolution = false\n";
  const auto frozen = cosmosim::core::loadFrozenConfigFromString(isolated_ok, "isolated_ok");
  CHUI_REQUIRE(frozen.config.mode.hydro_boundary == cosmosim::core::ModeHydroBoundary::kReflective);
  CHUI_REQUIRE(frozen.config.mode.gravity_boundary == cosmosim::core::ModeGravityBoundary::kIsolatedMonopole);
  CHUI_REQUIRE(!frozen.config.numerics.treepm_enable_window_deconvolution);
  CHUI_REQUIRE(frozen.normalized_text.find("treepm_enable_window_deconvolution = false") !=
      std::string::npos);

  const std::string isolated_incompatible_default =
      "[mode]\nmode = isolated_cluster\nhydro_boundary = reflective\ngravity_boundary = isolated_monopole\n";
  threw = false;
  try {
    (void)cosmosim::core::loadFrozenConfigFromString(
        isolated_incompatible_default, "isolated_incompatible_default");
  } catch (const cosmosim::core::ConfigError& ex) {
    threw = std::string(ex.what()).find(
                "treepm_enable_window_deconvolution must be false") !=
        std::string::npos;
  }
  CHUI_REQUIRE(threw);
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

  CHUI_REQUIRE(frozen_first.config.physics.enable_cooling);
  CHUI_REQUIRE(frozen_first.config.physics.fb_mode == cosmosim::core::FeedbackMode::kThermalKineticMomentum);
  CHUI_REQUIRE(frozen_first.config.parallel.deterministic_reduction);
  CHUI_REQUIRE(
      frozen_first.config.numerics.treepm_tree_opening_criterion ==
      cosmosim::core::TreePmOpeningCriterion::kComDistance);
  CHUI_REQUIRE(frozen_first.config.numerics.treepm_tree_opening_theta == 0.7);
  CHUI_REQUIRE(frozen_first.config.numerics.treepm_tree_relative_force_tolerance == 0.005);
  CHUI_REQUIRE(
      frozen_first.config.numerics.treepm_tree_relative_force_acceleration_floor == 1.0e-30);
  CHUI_REQUIRE(
      frozen_first.config.numerics.treepm_assignment_scheme ==
      cosmosim::core::TreePmAssignmentScheme::kTsc);
  CHUI_REQUIRE(frozen_first.config.numerics.treepm_enable_window_deconvolution);
  CHUI_REQUIRE(frozen_first.normalized_text.find("treepm_tree_opening_criterion = com_distance") !=
      std::string::npos);
  CHUI_REQUIRE(!frozen_first.config.parallel.decomposition_debug_exact_ownership_audit);
  CHUI_REQUIRE(frozen_first.normalized_text.find("decomposition_debug_exact_ownership_audit = false") != std::string::npos);
  CHUI_REQUIRE(frozen_first.config.parallel.isolated_pm_root_workspace_limit_bytes == 268435456ULL);
  CHUI_REQUIRE(frozen_first.config.parallel.zoom_high_res_allgather_limit_bytes == 268435456ULL);
  CHUI_REQUIRE(frozen_first.normalized_text.find("isolated_pm_root_workspace_limit_bytes = 268435456") != std::string::npos);
  CHUI_REQUIRE(frozen_first.normalized_text.find("zoom_high_res_allgather_limit_bytes = 268435456") != std::string::npos);
  CHUI_REQUIRE(frozen_first.normalized_text == frozen_second.normalized_text);
  CHUI_REQUIRE(frozen_first.provenance.config_hash_hex == frozen_second.provenance.config_hash_hex);
  CHUI_REQUIRE(frozen_first.provenance.config_hash_hex == cosmosim::core::stableConfigHashHex(frozen_first.normalized_text));
  const auto reparsed =
      cosmosim::core::loadFrozenConfigFromString(frozen_first.normalized_text, "normalized_roundtrip");
  CHUI_REQUIRE(reparsed.normalized_text == frozen_first.normalized_text);
  CHUI_REQUIRE(reparsed.provenance.config_hash_hex == frozen_first.provenance.config_hash_hex);
}

void testNormalizedConfigRoundTripPreservesDifficultDoubles() {
  const std::string text = R"(
[mode]
mode = cosmo_cube
ic_file = ../ics/chui_ic_32_10mpch_z24.hdf5
ic_convention = gadget_arepo_bridge_v1
ic_bridge_source_length_unit_to_si = 3.0856775814913673e22
ic_bridge_source_mass_unit_to_si = 1.98847e40
ic_bridge_source_velocity_unit_to_si = 1000
ic_bridge_coordinate_frame = comoving
ic_bridge_velocity_convention = gadget_arepo_stored_peculiar
ic_bridge_length_hubble_exponent = -1
ic_bridge_length_scale_factor_exponent = 0
ic_bridge_mass_hubble_exponent = -1
ic_bridge_mass_scale_factor_exponent = 0
ic_bridge_velocity_hubble_exponent = 0
ic_bridge_velocity_scale_factor_exponent = 0

[cosmology]
omega_matter = 0.3099
omega_lambda = 0.6901
omega_baryon = 0.0488911
hubble_param = 0.67742
sigma8 = 0.808992
scalar_index_ns = 0.96822
box_size_x = 14.761890702961235 mpc
box_size_y = 14.761890702961235 mpc
box_size_z = 14.761890702961235 mpc

[numerics]
a_begin = 0.04
a_end = 0.041
integrator_time_variable = scale_factor
cosmology_max_delta_ln_a = 0.007812345678901234
gravity_softening = 15 kpc
treepm_tree_opening_theta = 0.6123456789012345
treepm_tree_relative_force_tolerance = 0.0043210987654321
)";

  const auto frozen = cosmosim::core::loadFrozenConfigFromString(
      text, "normalized_difficult_doubles");
  const auto reparsed = cosmosim::core::loadFrozenConfigFromString(
      frozen.normalized_text, "normalized_difficult_doubles_roundtrip");

  CHUI_REQUIRE(reparsed.normalized_text == frozen.normalized_text);
  CHUI_REQUIRE(reparsed.provenance.config_hash_hex == frozen.provenance.config_hash_hex);
  CHUI_REQUIRE(reparsed.config.mode.ic_bridge_source_length_unit_to_si ==
         frozen.config.mode.ic_bridge_source_length_unit_to_si);
  CHUI_REQUIRE(reparsed.config.mode.ic_bridge_source_mass_unit_to_si ==
         frozen.config.mode.ic_bridge_source_mass_unit_to_si);
  CHUI_REQUIRE(reparsed.config.cosmology.omega_matter ==
         frozen.config.cosmology.omega_matter);
  CHUI_REQUIRE(reparsed.config.cosmology.omega_lambda ==
         frozen.config.cosmology.omega_lambda);
  CHUI_REQUIRE(reparsed.config.cosmology.hubble_param ==
         frozen.config.cosmology.hubble_param);
  CHUI_REQUIRE(reparsed.config.cosmology.box_size_x_mpc_comoving ==
         frozen.config.cosmology.box_size_x_mpc_comoving);
  CHUI_REQUIRE(reparsed.config.cosmology.box_size_y_mpc_comoving ==
         frozen.config.cosmology.box_size_y_mpc_comoving);
  CHUI_REQUIRE(reparsed.config.cosmology.box_size_z_mpc_comoving ==
         frozen.config.cosmology.box_size_z_mpc_comoving);
  CHUI_REQUIRE(reparsed.config.numerics.a_begin == frozen.config.numerics.a_begin);
  CHUI_REQUIRE(reparsed.config.numerics.a_end == frozen.config.numerics.a_end);
  CHUI_REQUIRE(reparsed.config.numerics.cosmology_max_delta_ln_a ==
         frozen.config.numerics.cosmology_max_delta_ln_a);
  CHUI_REQUIRE(reparsed.config.numerics.gravity_softening_kpc_comoving ==
         frozen.config.numerics.gravity_softening_kpc_comoving);
  CHUI_REQUIRE(reparsed.config.numerics.treepm_tree_opening_theta ==
         frozen.config.numerics.treepm_tree_opening_theta);
  CHUI_REQUIRE(reparsed.config.numerics.treepm_tree_relative_force_tolerance ==
         frozen.config.numerics.treepm_tree_relative_force_tolerance);
}

void testDebugExactDecompositionAuditConfig() {
  const std::string text = R"(
[mode]
mode = cosmo_cube
[parallel]
decomposition_debug_exact_ownership_audit = true
isolated_pm_root_workspace_limit_bytes = 1048576
zoom_high_res_allgather_limit_bytes = 2097152
process_memory_budget_bytes = 34359738368
process_mpi_unknown_reserve_bytes = 10485760
process_fftw_unknown_reserve_bytes = 20971520
process_hdf5_unknown_reserve_bytes = 31457280
process_allocator_unknown_reserve_bytes = 41943040
process_output_restart_overlap_bytes = 52428800
process_memory_safety_margin_fraction = 0.25
)";
  const auto frozen = cosmosim::core::loadFrozenConfigFromString(text, "debug_exact_decomposition_audit");
  CHUI_REQUIRE(frozen.config.parallel.decomposition_debug_exact_ownership_audit);
  CHUI_REQUIRE(frozen.config.parallel.isolated_pm_root_workspace_limit_bytes == 1048576ULL);
  CHUI_REQUIRE(frozen.config.parallel.zoom_high_res_allgather_limit_bytes == 2097152ULL);
  CHUI_REQUIRE(frozen.config.parallel.process_memory_budget_bytes == 34359738368ULL);
  CHUI_REQUIRE(frozen.config.parallel.process_mpi_unknown_reserve_bytes == 10485760ULL);
  CHUI_REQUIRE(frozen.config.parallel.process_fftw_unknown_reserve_bytes == 20971520ULL);
  CHUI_REQUIRE(frozen.config.parallel.process_hdf5_unknown_reserve_bytes == 31457280ULL);
  CHUI_REQUIRE(frozen.config.parallel.process_allocator_unknown_reserve_bytes == 41943040ULL);
  CHUI_REQUIRE(frozen.config.parallel.process_output_restart_overlap_bytes == 52428800ULL);
  CHUI_REQUIRE(std::abs(frozen.config.parallel.process_memory_safety_margin_fraction - 0.25) < 1.0e-15);
  CHUI_REQUIRE(frozen.normalized_text.find("process_memory_budget_bytes = 34359738368") != std::string::npos);
  CHUI_REQUIRE(frozen.normalized_text.find("decomposition_debug_exact_ownership_audit = true") != std::string::npos);
}

void testConfigNormalizationHashDeterminism() {
  const std::string first =
      "[mode]\nmode = zoom_in\n[output]\nrun_name = hash_a\n";
  const std::string second =
      "[output]\nrun_name = hash_a\n[mode]\nmode = zoom_in\n";
  const auto frozen_first = cosmosim::core::loadFrozenConfigFromString(first, "hash_first");
  const auto frozen_second = cosmosim::core::loadFrozenConfigFromString(second, "hash_second");
  CHUI_REQUIRE(frozen_first.normalized_text == frozen_second.normalized_text);
  CHUI_REQUIRE(frozen_first.provenance.config_hash_hex == frozen_second.provenance.config_hash_hex);
  CHUI_REQUIRE(frozen_first.provenance.config_hash_hex == cosmosim::core::stableConfigHashHex(frozen_first.normalized_text));
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
  CHUI_REQUIRE(frozen.config.numerics.treepm_pm_grid_nx == 32);
  CHUI_REQUIRE(frozen.config.numerics.treepm_pm_grid_ny == 32);
  CHUI_REQUIRE(frozen.config.numerics.treepm_pm_grid_nz == 32);
  CHUI_REQUIRE(frozen.normalized_text.find("treepm_pm_grid =") == std::string::npos);
  CHUI_REQUIRE(frozen.normalized_text.find("treepm_pm_grid_nx = 32") != std::string::npos);
  const auto reparsed = cosmosim::core::loadFrozenConfigFromString(frozen.normalized_text, "runtime_ownership_reparse");
  CHUI_REQUIRE(reparsed.config.numerics.treepm_pm_grid_nx == frozen.config.numerics.treepm_pm_grid_nx);
  CHUI_REQUIRE(reparsed.config.numerics.treepm_asmth_cells == frozen.config.numerics.treepm_asmth_cells);
  CHUI_REQUIRE(reparsed.config.numerics.treepm_rcut_cells == frozen.config.numerics.treepm_rcut_cells);
  CHUI_REQUIRE(reparsed.provenance.config_hash_hex == frozen.provenance.config_hash_hex);
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
t_code_begin = 0.125
t_code_end = 0.875
gravity_softening = 1.0 kpc
gravity_softening_gas = 0.5 kpc
gravity_softening_black_hole = 0.05 kpc
treepm_pm_grid_nx = 28
treepm_pm_grid_ny = 20
treepm_pm_grid_nz = 12
treepm_rcut_cells = 5.5
)";
  const auto frozen = cosmosim::core::loadFrozenConfigFromString(text, "derived_runtime");
  const auto derived = cosmosim::core::deriveRuntimeConfig(frozen);
  CHUI_REQUIRE(derived.t_code_begin == 0.125);
  CHUI_REQUIRE(derived.t_code_end == 0.875);
  CHUI_REQUIRE(derived.box_size_mpc_comoving[0] == 70.0);
  CHUI_REQUIRE(derived.box_size_mpc_comoving[1] == 50.0);
  CHUI_REQUIRE(derived.box_size_mpc_comoving[2] == 30.0);
  CHUI_REQUIRE(derived.treepm_pm_grid_shape[0] == 28);
  CHUI_REQUIRE(derived.treepm_pm_grid_shape[1] == 20);
  CHUI_REQUIRE(derived.treepm_pm_grid_shape[2] == 12);
  CHUI_REQUIRE(derived.gravity_softening_kpc_comoving_by_species[0] == 1.0);
  CHUI_REQUIRE(derived.gravity_softening_kpc_comoving_by_species[1] == 0.5);
  CHUI_REQUIRE(derived.gravity_softening_kpc_comoving_by_species[2] == 1.0);
  CHUI_REQUIRE(derived.gravity_softening_kpc_comoving_by_species[3] == 0.05);
  CHUI_REQUIRE(derived.gravity_softening_kpc_comoving_by_species[4] == 1.0);
  CHUI_REQUIRE(derived.normalized_config_hash == frozen.provenance.config_hash);
  CHUI_REQUIRE(derived.normalized_config_hash_hex == frozen.provenance.config_hash_hex);

  const std::string bad_time = "[mode]\nmode = zoom_in\n[numerics]\nt_code_begin = 1.0\nt_code_end = 0.5\n";
  bool threw = false;
  try {
    const auto bad = cosmosim::core::loadFrozenConfigFromString(bad_time, "derived_runtime_bad_time");
    (void)cosmosim::core::deriveRuntimeConfig(bad);
  } catch (const cosmosim::core::ConfigError&) {
    threw = true;
  }
  CHUI_REQUIRE(threw);
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
  CHUI_REQUIRE(frozen.config.physics.metal_line_table_path == "https://example.com/tables/metal_line_rates.h5");
  CHUI_REQUIRE(frozen.config.physics.stellar_evolution_table_path == "\"C://stellar_tables//ssp#grid.h5\"");
}

void testDeprecatedAliasesAndCanonicalCollision() {
  const std::string alias_only =
      "mode = isolated_galaxy\n"
      "numerics.treepm_enable_window_deconvolution = false\n";
  const auto frozen = cosmosim::core::loadFrozenConfigFromString(alias_only, "alias_only");
  CHUI_REQUIRE(frozen.config.mode.mode == cosmosim::core::SimulationMode::kIsolatedGalaxy);
  CHUI_REQUIRE(!frozen.provenance.deprecation_warnings.empty());
  CHUI_REQUIRE(!frozen.user_config.alias_resolution_notes.empty());

  const std::string conflicting = "mode = zoom_in\nmode.mode = cosmo_cube\n";
  bool threw = false;
  try {
    (void)cosmosim::core::loadFrozenConfigFromString(conflicting, "alias_conflict");
  } catch (const cosmosim::core::ConfigError&) {
    threw = true;
  }
  CHUI_REQUIRE(threw);
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
  CHUI_REQUIRE(frozen.config.physics.fb_mode == cosmosim::core::FeedbackMode::kMomentum);
  CHUI_REQUIRE(frozen.config.physics.fb_variant == cosmosim::core::FeedbackVariant::kStochastic);
  CHUI_REQUIRE(frozen.config.physics.fb_neighbor_count == 4);

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
  CHUI_REQUIRE(threw);
}

void testMetalDiffusionConfigurationContract() {
  const std::string good_text = R"(
[mode]
mode = zoom_in
[numerics]
hierarchical_max_rung = 0
[physics]
metal_species_mode = total_only
enable_metal_diffusion = true
metal_diffusion_model = smagorinsky
metal_diffusion_time_integrator = rkl2
metal_diffusion_coefficient = 0.075
metal_diffusion_cfl = 0.35
metal_diffusion_max_subcycles = 64
metal_diffusion_max_rkl_stages = 48
metal_diffusion_coefficient_floor_code = 0.0
metal_diffusion_coefficient_ceiling_code = 100.0
)";
  const auto frozen = cosmosim::core::loadFrozenConfigFromString(
      good_text, "metal_diffusion_good");
  CHUI_REQUIRE(frozen.config.physics.metal_species_mode ==
         cosmosim::core::MetalSpeciesMode::kTotalOnly);
  CHUI_REQUIRE(frozen.config.physics.metal_diffusion_model ==
         cosmosim::core::MetalDiffusionModel::kSmagorinsky);
  CHUI_REQUIRE(frozen.config.physics.metal_diffusion_time_integrator ==
         cosmosim::core::MetalDiffusionTimeIntegrator::kRkl2);
  CHUI_REQUIRE(frozen.normalized_text.find("metal_diffusion_coefficient = 0.075") !=
         std::string::npos);
  const auto roundtrip = cosmosim::core::loadFrozenConfigFromString(
      frozen.normalized_text, "metal_diffusion_roundtrip");
  CHUI_REQUIRE(roundtrip.provenance.config_hash == frozen.provenance.config_hash);

  for (const char* invalid : {
           "[mode]\nmode = zoom_in\n[physics]\nenable_metal_diffusion = true\nmetal_diffusion_model = none\n",
           "[mode]\nmode = zoom_in\n[physics]\nenable_metal_diffusion = true\nmetal_diffusion_model = smagorinsky\nmetal_diffusion_coefficient = 0.0\n",
           "[mode]\nmode = zoom_in\n[physics]\nmetal_species_mode = core_elements\n"}) {
    bool rejected = false;
    try {
      (void)cosmosim::core::loadFrozenConfigFromString(invalid, "metal_diffusion_bad");
    } catch (const cosmosim::core::ConfigError&) {
      rejected = true;
    }
    CHUI_REQUIRE(rejected);
  }
}

void testCoolingPolicyEnumsAndValidation() {
  const std::string good_text = R"(
[mode]
mode = zoom_in
[physics]
uv_background_model = fg20
self_shielding_model = rahmati13_like
cooling_model = primordial_metal_line
metal_line_table_path = resources/cooling/metal_line.tbl
)";
  const auto frozen = cosmosim::core::loadFrozenConfigFromString(good_text, "cooling_good");
  CHUI_REQUIRE(frozen.config.physics.uv_background_model == cosmosim::core::UvBackgroundModel::kFg20);
  CHUI_REQUIRE(frozen.config.physics.self_shielding_model == cosmosim::core::SelfShieldingModel::kRahmati13Like);
  CHUI_REQUIRE(frozen.config.physics.cooling_model == cosmosim::core::CoolingModel::kPrimordialMetalLine);

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
  CHUI_REQUIRE(threw);
}

void testDiagnosticsExecutionPolicyValidation() {
  const std::string good_text = R"(
[mode]
mode = zoom_in
[analysis]
diagnostics_execution_policy = all_including_provisional
)";
  const auto frozen = cosmosim::core::loadFrozenConfigFromString(good_text, "diag_policy_good");
  CHUI_REQUIRE(
      frozen.config.analysis.diagnostics_execution_policy ==
      cosmosim::core::AnalysisConfig::DiagnosticsExecutionPolicy::kAllIncludingProvisional);
  CHUI_REQUIRE(
      frozen.normalized_text.find("diagnostics_execution_policy = all_including_provisional") !=
      std::string::npos);
  CHUI_REQUIRE(frozen.normalized_text.find("diagnostics_execution_policy = unknown") == std::string::npos);

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
  CHUI_REQUIRE(threw);
}

void testEnumSerializationIsFailFastWithoutUnknownFallback() {
  bool threw = false;
  try {
    (void)cosmosim::core::modeToString(static_cast<cosmosim::core::SimulationMode>(255));
  } catch (const cosmosim::core::ConfigError&) {
    threw = true;
  }
  CHUI_REQUIRE(threw);
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
  CHUI_REQUIRE(std::abs(frozen.config.numerics.gravity_softening_gas_kpc_comoving - 0.5) < 1.0e-12);
  CHUI_REQUIRE(std::abs(frozen.config.numerics.gravity_softening_dark_matter_kpc_comoving - 2.0) < 1.0e-12);
  CHUI_REQUIRE(std::abs(frozen.config.numerics.gravity_softening_star_kpc_comoving - 0.2) < 1.0e-12);
  CHUI_REQUIRE(std::abs(frozen.config.numerics.gravity_softening_black_hole_kpc_comoving - 0.05) < 1.0e-12);
  CHUI_REQUIRE(std::abs(frozen.config.numerics.gravity_softening_tracer_kpc_comoving - 1.5) < 1.0e-12);
  CHUI_REQUIRE(frozen.normalized_text.find("gravity_softening_gas = ") != std::string::npos);
  const auto reparsed = cosmosim::core::loadFrozenConfigFromString(frozen.normalized_text, "species_softening_roundtrip");
  CHUI_REQUIRE(std::abs(reparsed.config.numerics.gravity_softening_black_hole_kpc_comoving - 0.05) < 1.0e-12);
}

void testHydroActiveBatchConfigurationRoundtrip() {
  const auto frozen = cosmosim::core::loadFrozenConfigFromString(
      "[mode]\nmode = zoom_in\n[numerics]\nhydro_active_batch_max_cells = 37\n",
      "hydro_active_batch_roundtrip");
  CHUI_REQUIRE(frozen.config.numerics.hydro_active_batch_max_cells == 37ULL);
  CHUI_REQUIRE(frozen.normalized_text.find("hydro_active_batch_max_cells = 37") !=
         std::string::npos);
  const auto reparsed = cosmosim::core::loadFrozenConfigFromString(
      frozen.normalized_text, "hydro_active_batch_roundtrip_reparse");
  CHUI_REQUIRE(reparsed.config.numerics.hydro_active_batch_max_cells == 37ULL);

  const auto automatic = cosmosim::core::loadFrozenConfigFromString(
      "[mode]\nmode = zoom_in\n", "hydro_active_batch_auto");
  CHUI_REQUIRE(automatic.config.numerics.hydro_active_batch_max_cells == 0ULL);
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
treepm_tree_opening_criterion = relative_force_error
treepm_tree_opening_theta = 0.65
treepm_tree_relative_force_tolerance = 0.006
treepm_tree_relative_force_acceleration_floor = 1.0e-25
treepm_assignment_scheme = cic
treepm_enable_window_deconvolution = true
treepm_update_cadence_steps = 1
treepm_pm_decomposition_mode = slab
treepm_tree_exchange_batch_bytes = 8388608
treepm_pm_exchange_batch_bytes = 12582912
)";
  const auto frozen = cosmosim::core::loadFrozenConfigFromString(good_text, "treepm_good");
  CHUI_REQUIRE(frozen.config.numerics.treepm_pm_grid_nx == 32);
  CHUI_REQUIRE(frozen.config.numerics.treepm_pm_grid_ny == 32);
  CHUI_REQUIRE(frozen.config.numerics.treepm_pm_grid_nz == 32);
  CHUI_REQUIRE(frozen.config.numerics.treepm_asmth_cells == 1.5);
  CHUI_REQUIRE(frozen.config.numerics.treepm_rcut_cells == 5.0);
  CHUI_REQUIRE(
      frozen.config.numerics.treepm_tree_opening_criterion ==
      cosmosim::core::TreePmOpeningCriterion::kRelativeForceError);
  CHUI_REQUIRE(frozen.config.numerics.treepm_tree_opening_theta == 0.65);
  CHUI_REQUIRE(frozen.config.numerics.treepm_tree_relative_force_tolerance == 0.006);
  CHUI_REQUIRE(frozen.config.numerics.treepm_tree_relative_force_acceleration_floor == 1.0e-25);
  CHUI_REQUIRE(
      frozen.config.numerics.treepm_assignment_scheme ==
      cosmosim::core::TreePmAssignmentScheme::kCic);
  CHUI_REQUIRE(frozen.config.numerics.treepm_enable_window_deconvolution);
  CHUI_REQUIRE(frozen.config.numerics.treepm_update_cadence_steps == 1);
  CHUI_REQUIRE(
      frozen.config.numerics.treepm_pm_decomposition_mode ==
      cosmosim::core::PmDecompositionMode::kSlab);
  CHUI_REQUIRE(frozen.config.numerics.treepm_tree_exchange_batch_bytes == 8388608ULL);
  CHUI_REQUIRE(frozen.config.numerics.treepm_pm_exchange_batch_bytes == 12582912ULL);
  CHUI_REQUIRE(frozen.normalized_text.find("treepm_pm_grid_nx = 32") != std::string::npos);
  CHUI_REQUIRE(frozen.normalized_text.find("treepm_tree_opening_criterion = relative_force_error") !=
      std::string::npos);
  CHUI_REQUIRE(frozen.normalized_text.find("treepm_tree_opening_theta = 0.65") != std::string::npos);
  CHUI_REQUIRE(frozen.normalized_text.find("treepm_assignment_scheme = cic") != std::string::npos);
  CHUI_REQUIRE(frozen.normalized_text.find("treepm_pm_decomposition_mode = slab") != std::string::npos);
  const auto reparsed =
      cosmosim::core::loadFrozenConfigFromString(frozen.normalized_text, "treepm_roundtrip");
  CHUI_REQUIRE(reparsed.config.numerics.treepm_pm_grid_nx == 32);
  CHUI_REQUIRE(
      reparsed.config.numerics.treepm_tree_opening_criterion ==
      cosmosim::core::TreePmOpeningCriterion::kRelativeForceError);
  CHUI_REQUIRE(reparsed.config.numerics.treepm_tree_opening_theta == 0.65);
  CHUI_REQUIRE(reparsed.config.numerics.treepm_tree_relative_force_tolerance == 0.006);
  CHUI_REQUIRE(reparsed.config.numerics.treepm_tree_relative_force_acceleration_floor == 1.0e-25);
  CHUI_REQUIRE(reparsed.config.numerics.treepm_update_cadence_steps == 1);
  CHUI_REQUIRE(reparsed.config.numerics.treepm_tree_exchange_batch_bytes == 8388608ULL);
  CHUI_REQUIRE(reparsed.config.numerics.treepm_pm_exchange_batch_bytes == 12582912ULL);

  const auto parsed_geometric = cosmosim::core::loadFrozenConfigFromString(
      "[mode]\nmode = zoom_in\n[numerics]\ntreepm_tree_opening_criterion = geometric\n",
      "treepm_geometric_opening");
  CHUI_REQUIRE(
      parsed_geometric.config.numerics.treepm_tree_opening_criterion ==
      cosmosim::core::TreePmOpeningCriterion::kGeometric);
  CHUI_REQUIRE(parsed_geometric.normalized_text.find("treepm_tree_opening_criterion = geometric") !=
      std::string::npos);

  const auto& baseline = frozen.normalized_text;
  // Every group has a typed legal control above (or below for pencil mode).
  requireNumericsRejected(baseline, "treepm_pm_grid_nx", "0", "treepm_pm_grid_n{xyz}");
  for (const auto key : {"treepm_tree_opening_theta", "treepm_tree_relative_force_tolerance",
                        "treepm_tree_relative_force_acceleration_floor", "treepm_asmth_cells", "treepm_rcut_cells"}) {
    for (const auto value : {"0", "-0.001", "nan", "inf", "-inf"})
      requireNumericsRejected(baseline, key, value);
  }
  for (const auto key : {"treepm_tree_opening_criterion", "treepm_assignment_scheme", "treepm_pm_decomposition_mode"})
    requireNumericsRejected(baseline, key, "magic");
  requireNumericsRejected(baseline, "treepm_update_cadence_steps", "3", "must be exactly 1 in production");
  for (const auto key : {"treepm_tree_exchange_batch_bytes", "treepm_pm_exchange_batch_bytes"}) {
    for (const auto value : {"0", "-1", "3.5", "nan", "inf", "18446744073709551616"})
      requireNumericsRejected(baseline, key, value);
  }
  const auto parsed_pencil = cosmosim::core::loadFrozenConfigFromString(
      withNumericsValue(baseline, "treepm_pm_decomposition_mode", "pencil"), "treepm_pencil_decomp");
  CHUI_REQUIRE(parsed_pencil.config.numerics.treepm_pm_decomposition_mode ==
      cosmosim::core::PmDecompositionMode::kPencil);

  const auto legal_cutoff = cosmosim::core::loadFrozenConfigFromString(
      "[mode]\nmode = cosmo_cube\nic_file = generated\n"
      "[cosmology]\nbox_size = 1 mpc\n"
      "[numerics]\ngravity_solver = treepm\ntreepm_pm_grid = 8\ntreepm_rcut_cells = 3.9\n",
      "treepm_cutoff_control");
  requireNumericsRejected(legal_cutoff.normalized_text, "treepm_rcut_cells", "4",
      "< half the shortest box axis");
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
  CHUI_REQUIRE(scalar_frozen.config.cosmology.box_size_x_mpc_comoving == 64.0);
  CHUI_REQUIRE(scalar_frozen.config.cosmology.box_size_y_mpc_comoving == 64.0);
  CHUI_REQUIRE(scalar_frozen.config.cosmology.box_size_z_mpc_comoving == 64.0);
  CHUI_REQUIRE(scalar_frozen.config.numerics.treepm_pm_grid_nx == 24);
  CHUI_REQUIRE(scalar_frozen.config.numerics.treepm_pm_grid_ny == 24);
  CHUI_REQUIRE(scalar_frozen.config.numerics.treepm_pm_grid_nz == 24);
  CHUI_REQUIRE(scalar_frozen.normalized_text.find("box_size_x = 64 mpc") != std::string::npos);
  CHUI_REQUIRE(scalar_frozen.normalized_text.find("treepm_pm_grid_nx = 24") != std::string::npos);
  CHUI_REQUIRE(scalar_frozen.normalized_text.find("treepm_pm_grid =") == std::string::npos);

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
treepm_rcut_cells = 5.5
)";
  const auto axis_frozen = cosmosim::core::loadFrozenConfigFromString(axis_text, "axis_canonical");
  CHUI_REQUIRE(axis_frozen.config.cosmology.box_size_x_mpc_comoving == 70.0);
  CHUI_REQUIRE(axis_frozen.config.cosmology.box_size_y_mpc_comoving == 50.0);
  CHUI_REQUIRE(axis_frozen.config.cosmology.box_size_z_mpc_comoving == 30.0);
  CHUI_REQUIRE(axis_frozen.config.numerics.treepm_pm_grid_nx == 28);
  CHUI_REQUIRE(axis_frozen.config.numerics.treepm_pm_grid_ny == 20);
  CHUI_REQUIRE(axis_frozen.config.numerics.treepm_pm_grid_nz == 12);

  const auto anisotropic_control = cosmosim::core::loadFrozenConfigFromString(
      "[mode]\nmode = zoom_in\n[cosmology]\nbox_size_x = 80 mpc\nbox_size_y = 20 mpc\nbox_size_z = 20 mpc\n"
      "[numerics]\ntreepm_pm_grid = 20\ntreepm_asmth_cells = 3\n", "anisotropic_control");
  requireNumericsRejected(anisotropic_control.normalized_text, "treepm_asmth_cells", "1.25",
      "coarsest PM cell spacing");

  const auto certified_softening = cosmosim::core::loadFrozenConfigFromString(
      R"(
[mode]
mode = zoom_in
[cosmology]
box_size = 64 mpc
[numerics]
treepm_pm_grid = 32
treepm_asmth_cells = 1.25
gravity_softening = 400 kpc
)",
      "treepm_certified_softening_envelope");
  CHUI_REQUIRE(certified_softening.config.numerics.gravity_softening_kpc_comoving == 400.0);

  requireNumericsRejected(certified_softening.normalized_text, "gravity_softening", "600 kpc",
      "epsilon/r_s <= 0.20");
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
  CHUI_REQUIRE(frozen.config.physics.enable_black_hole_agn);
  CHUI_REQUIRE(frozen.config.physics.bh_seed_halo_mass_threshold_code == 500.0);
  CHUI_REQUIRE(frozen.config.physics.bh_seed_max_per_cell == 2);
  CHUI_REQUIRE(frozen.config.physics.bh_alpha_bondi == 3.0);

  const std::string bad_text =
      "[mode]\nmode = zoom_in\n[physics]\nbh_feedback_coupling_efficiency = 1.3\n";
  bool threw = false;
  try {
    (void)cosmosim::core::loadFrozenConfigFromString(bad_text, "bh_bad");
  } catch (const cosmosim::core::ConfigError&) {
    threw = true;
  }
  CHUI_REQUIRE(threw);
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
  CHUI_REQUIRE(frozen.config.physics.enable_tracers);
  CHUI_REQUIRE(frozen.config.physics.tracer_track_mass);
  CHUI_REQUIRE(frozen.config.physics.tracer_min_host_mass_code == 1.0e-6);

  const std::string bad_text = "[mode]\nmode = zoom_in\n[physics]\ntracer_min_host_mass_code = -1.0\n";
  bool threw = false;
  try {
    (void)cosmosim::core::loadFrozenConfigFromString(bad_text, "tracer_bad");
  } catch (const cosmosim::core::ConfigError&) {
    threw = true;
  }
  CHUI_REQUIRE(threw);
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
  CHUI_REQUIRE(frozen.config.mode.zoom_long_range_strategy ==
      cosmosim::core::ZoomLongRangeStrategy::kGlobalCoarsePlusFocusedHighResCorrection);
  CHUI_REQUIRE(frozen.config.mode.zoom_region_center_x_mpc_comoving == 12.5);
  CHUI_REQUIRE(frozen.config.mode.zoom_focused_pm_grid_nz == 64);
  CHUI_REQUIRE(frozen.normalized_text.find("zoom_long_range_strategy = global_coarse_plus_focused_highres_correction") !=
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
  CHUI_REQUIRE(threw);
}


void testCosmologyScaleFactorRedshiftCanonicalizationAndValidation() {
  const auto z_only = cosmosim::core::loadFrozenConfigFromString(
      "[mode]\nmode = zoom_in\n[numerics]\nz_begin = 2.0\nz_end = 1.0\n",
      "z_only");
  CHUI_REQUIRE(std::abs(z_only.config.numerics.a_begin - (1.0 / 3.0)) < 1.0e-12);
  CHUI_REQUIRE(std::abs(z_only.config.numerics.a_end - 0.5) < 1.0e-12);
  CHUI_REQUIRE(z_only.normalized_text.find("a_begin = 0.333333") != std::string::npos);

  const auto a_only = cosmosim::core::loadFrozenConfigFromString(
      "[mode]\nmode = zoom_in\n[numerics]\na_begin = 0.25\na_end = 0.5\n",
      "a_only");
  CHUI_REQUIRE(std::abs(a_only.config.numerics.z_begin - 3.0) < 1.0e-12);
  CHUI_REQUIRE(std::abs(a_only.config.numerics.z_end - 1.0) < 1.0e-12);

  bool threw = false;
  try {
    (void)cosmosim::core::loadFrozenConfigFromString(
        "[mode]\nmode = zoom_in\n[numerics]\na_begin = 1.0\nz_begin = 2.0\n",
        "az_conflict");
  } catch (const cosmosim::core::ConfigError&) {
    threw = true;
  }
  CHUI_REQUIRE(threw);

  threw = false;
  try {
    (void)cosmosim::core::loadFrozenConfigFromString(
        "[mode]\nmode = zoom_in\n[numerics]\nz_begin = -1.0\n",
        "bad_z");
  } catch (const cosmosim::core::ConfigError&) {
    threw = true;
  }
  CHUI_REQUIRE(threw);
}

void testIntegratorTimeVariableIsTypedAndCanonical() {
  const auto frozen = cosmosim::core::loadFrozenConfigFromString(
      "[mode]\nmode = zoom_in\n[numerics]\na_begin = 0.5\na_end = 0.51\nintegrator_time_variable = ln_a\n",
      "typed_integrator");
  CHUI_REQUIRE(frozen.config.numerics.integrator_time_variable ==
      cosmosim::core::IntegratorTimeVariable::kLogScaleFactor);
  CHUI_REQUIRE(frozen.normalized_text.find("integrator_time_variable = ln_a") != std::string::npos);

  const auto alias = cosmosim::core::loadFrozenConfigFromString(
      "[mode]\nmode = zoom_in\n[numerics]\na_begin = 0.5\na_end = 0.51\nintegrator_time_variable = a\n",
      "typed_integrator_alias");
  CHUI_REQUIRE(alias.config.numerics.integrator_time_variable ==
      cosmosim::core::IntegratorTimeVariable::kScaleFactor);
  CHUI_REQUIRE(alias.normalized_text.find("integrator_time_variable = scale_factor") != std::string::npos);

  bool threw = false;
  try {
    (void)cosmosim::core::loadFrozenConfigFromString(
        "[mode]\nmode = zoom_in\n[numerics]\nintegrator_time_variable = banana\n",
        "bad_integrator_variable");
  } catch (const cosmosim::core::ConfigError&) {
    threw = true;
  }
  CHUI_REQUIRE(threw);
}

void testIntegratorEndpointAuthorityResolution() {
  const auto scale_factor = cosmosim::core::loadFrozenConfigFromString(
      "[mode]\nmode = zoom_in\n[numerics]\na_begin = 0.04\na_end = 0.041\n"
      "integrator_time_variable = scale_factor\n",
      "scale_factor_authority");
  cosmosim::core::CosmologyBackgroundConfig background_config;
  background_config.hubble_param = scale_factor.config.cosmology.hubble_param;
  background_config.omega_matter = scale_factor.config.cosmology.omega_matter;
  background_config.omega_lambda = scale_factor.config.cosmology.omega_lambda;
  const cosmosim::core::LambdaCdmBackground background(background_config);
  const auto units = cosmosim::core::makeUnitSystem(
      scale_factor.config.units.length_unit, scale_factor.config.units.mass_unit,
      scale_factor.config.units.velocity_unit);
  const double expected_t_end = scale_factor.config.numerics.t_code_begin +
      background.cosmicTimeIntervalSi(0.04, 0.041) / units.timeSiPerCode();
  CHUI_REQUIRE(std::abs(scale_factor.config.numerics.t_code_end - expected_t_end) <=
         1.0e-12 * std::max(1.0, std::abs(expected_t_end)));
  CHUI_REQUIRE(scale_factor.normalized_text.find("integrator_time_variable = scale_factor") !=
         std::string::npos);

  bool contradiction_rejected = false;
  try {
    (void)cosmosim::core::loadFrozenConfigFromString(
        "[mode]\nmode = zoom_in\n[numerics]\na_begin = 0.04\na_end = 0.041\n"
        "t_code_begin = 0\nt_code_end = 1\nintegrator_time_variable = scale_factor\n",
        "scale_factor_contradiction");
  } catch (const cosmosim::core::ConfigError&) {
    contradiction_rejected = true;
  }
  CHUI_REQUIRE(contradiction_rejected);

  const auto code_time = cosmosim::core::loadFrozenConfigFromString(
      "[mode]\nmode = zoom_in\n[numerics]\na_begin = 0.04\n"
      "t_code_begin = 0\nt_code_end = 0.001\nintegrator_time_variable = code_time\n",
      "code_time_authority");
  CHUI_REQUIRE(code_time.config.numerics.a_end > code_time.config.numerics.a_begin);
  CHUI_REQUIRE(std::abs(
      code_time.config.numerics.z_end -
      (1.0 / code_time.config.numerics.a_end - 1.0)) < 1.0e-12);

  bool physical_time_rejected = false;
  try {
    (void)cosmosim::core::loadFrozenConfigFromString(
        "[mode]\nmode = zoom_in\n[numerics]\nintegrator_time_variable = physical_time\n",
        "physical_time_unsupported");
  } catch (const cosmosim::core::ConfigError&) {
    physical_time_rejected = true;
  }
  CHUI_REQUIRE(physical_time_rejected);
}

void testAdversarialPhysicsAndCosmologyDependenciesFail() {
  const char* cases[] = {
      "[mode]\nmode = zoom_in\n[cosmology]\nomega_matter = 0.1\nomega_baryon = 0.2\n",
      "[mode]\nmode = zoom_in\n[numerics]\nt_phys_begin = 10\nt_phys_end = 1\n",
      "[mode]\nmode = zoom_in\n[physics]\nenable_star_formation = false\nenable_feedback = true\n",
      "[mode]\nmode = zoom_in\n[physics]\nenable_feedback = false\nenable_black_hole_agn = true\n",
      "[mode]\nmode = zoom_in\n[physics]\ncooling_model = primordial_metal_line\nmetal_line_table_path =\n",
  };
  for (const char* config_text : cases) {
    bool threw = false;
    try {
      (void)cosmosim::core::loadFrozenConfigFromString(config_text, "adversarial_dependency");
    } catch (const cosmosim::core::ConfigError&) {
      threw = true;
    }
    CHUI_REQUIRE(threw);
  }
}


void testFiniteNumericAndForwardCosmologyContract() {
  const char* cases[] = {
      "[mode]\nmode = zoom_in\n[cosmology]\nomega_matter = inf\n",
      "[mode]\nmode = zoom_in\n[cosmology]\nhubble_param = nan\n",
      "[mode]\nmode = zoom_in\n[cosmology]\nbox_size = inf mpc\n",
      "[mode]\nmode = zoom_in\n[numerics]\na_begin = 1.0\na_end = 0.5\n",
      "[mode]\nmode = zoom_in\n[numerics]\nz_begin = 0.0\nz_end = 1.0\n",
      "[mode]\nmode = zoom_in\n[physics]\nfb_epsilon_thermal = nan\n",
      "[mode]\nmode = zoom_in\n[analysis]\nhalo_fof_linking_length_factor = inf\n",
  };
  for (const char* config_text : cases) {
    bool threw = false;
    try {
      (void)cosmosim::core::loadFrozenConfigFromString(config_text, "finite_forward_contract");
    } catch (const cosmosim::core::ConfigError&) {
      threw = true;
    }
    CHUI_REQUIRE(threw);
  }

  const auto forward = cosmosim::core::loadFrozenConfigFromString(
      "[mode]\nmode = zoom_in\n[numerics]\na_begin = 0.25\na_end = 1.0\n",
      "forward_cosmology_ok");
  CHUI_REQUIRE(forward.config.numerics.a_end >= forward.config.numerics.a_begin);
}

void testDerivedRuntimeSerializationUsesCanonicalNames() {
  const auto frozen = cosmosim::core::loadFrozenConfigFromString(
      "[mode]\nmode = zoom_in\n[numerics]\nt_code_begin = 0.25\nt_code_end = 0.75\n",
      "derived_canonical_names");
  const std::string serialized =
      cosmosim::core::serializeDerivedRuntimeConfig(cosmosim::core::deriveRuntimeConfig(frozen));
  CHUI_REQUIRE(serialized.find("t_code_begin=0.25") != std::string::npos);
  CHUI_REQUIRE(serialized.find("t_code_end=0.75") != std::string::npos);
  CHUI_REQUIRE(serialized.find("time_begin_code") == std::string::npos);
  CHUI_REQUIRE(serialized.find("time_end_code") == std::string::npos);
}


void testStarFormationMassParametersAcceptUnits() {
  const auto frozen = cosmosim::core::loadFrozenConfigFromString(
      "[units]\n"
      "mass_unit = kg\n"
      "[mode]\n"
      "mode = zoom_in\n"
      "[physics]\n"
      "star_formation_model = adaptive_bound_jeans\n"
      "sf_jeans_mass_floor_code = 2 msun\n"
      "sf_target_star_particle_mass_code = 3 msun\n"
      "sf_min_star_particle_mass_code = 1 msun\n"
      "sf_max_star_particle_mass_code = 10 msun\n"
      "sf_min_remaining_gas_mass_code = 0.5 msun\n",
      "sf_mass_units");
  const double solar_mass_kg = 1.98847e30;
  CHUI_REQUIRE(std::abs(frozen.config.physics.sf_jeans_mass_floor_code - 2.0 * solar_mass_kg) /
             solar_mass_kg < 1.0e-12);
  CHUI_REQUIRE(std::abs(frozen.config.physics.sf_target_star_particle_mass_code - 3.0 * solar_mass_kg) /
             solar_mass_kg < 1.0e-12);
  CHUI_REQUIRE(frozen.normalized_text.find("sf_jeans_mass_floor_code = ") != std::string::npos);

  bool rejected = false;
  try {
    (void)cosmosim::core::loadFrozenConfigFromString(
        "[mode]\nmode = zoom_in\n[physics]\n"
        "sf_target_star_particle_mass_code = 4 bananas\n",
        "sf_bad_mass_unit");
  } catch (const cosmosim::core::ConfigError&) {
    rejected = true;
  }
  CHUI_REQUIRE(rejected);
}


void testEffectiveMultiphaseConfigurationContract() {
  const auto frozen = cosmosim::core::loadFrozenConfigFromString(
      "[mode]\nmode = cosmo_cube\n"
      "[physics]\n"
      "enable_star_formation = true\n"
      "star_formation_model = effective_multiphase_tng_like\n"
      "sf_effective_n_h_threshold = 0.13 cm^-3\n"
      "sf_effective_t_star_at_threshold = 2.2 Gyr\n"
      "sf_effective_q_eos = 0.3\n"
      "sf_effective_feedback_coupling = external_feedback_calibrated\n",
      "effective_multiphase_config");
  CHUI_REQUIRE(frozen.config.physics.star_formation_model ==
         cosmosim::core::StarFormationModelKind::kEffectiveMultiphaseTngLike);
  CHUI_REQUIRE(std::abs(frozen.config.physics.sf_effective_n_h_threshold_cgs - 0.13) < 1.0e-14);
  CHUI_REQUIRE(frozen.config.physics.sf_effective_t_star_at_threshold_code > 0.0);
  CHUI_REQUIRE(std::abs(frozen.config.physics.sf_effective_q_eos - 0.3) < 1.0e-14);
  CHUI_REQUIRE(frozen.normalized_text.find(
             "star_formation_model = effective_multiphase_tng_like") != std::string::npos);
  CHUI_REQUIRE(frozen.normalized_text.find("sf_effective_parameter_set = chui_sh03_tng_like_v1") !=
         std::string::npos);
  const auto reparsed = cosmosim::core::loadFrozenConfigFromString(
      frozen.normalized_text, "effective_multiphase_config_roundtrip");
  CHUI_REQUIRE(reparsed.provenance.config_hash == frozen.provenance.config_hash);

  for (const char* invalid : {
           "[mode]\nmode = cosmo_cube\n[physics]\nstar_formation_model = unknown_model\n",
           "[mode]\nmode = cosmo_cube\n[physics]\nstar_formation_model = effective_multiphase_tng_like\nsf_effective_q_eos = 1.1\n",
           "[mode]\nmode = cosmo_cube\n[physics]\nstar_formation_model = effective_multiphase_tng_like\nsf_effective_n_h_threshold = -0.1 cm^-3\n",
           "[mode]\nmode = cosmo_cube\n[physics]\nenable_feedback = true\nstar_formation_model = effective_multiphase_tng_like\nsf_effective_feedback_coupling = effective_eos_only\n"}) {
    bool rejected = false;
    try {
      (void)cosmosim::core::loadFrozenConfigFromString(invalid, "invalid_effective_config");
    } catch (const cosmosim::core::ConfigError&) {
      rejected = true;
    }
    CHUI_REQUIRE(rejected);
  }
}

void testProductionHierarchicalRungsFailClosed() {
  using cosmosim::core::loadFrozenConfigFromString;
  const auto production = loadFrozenConfigFromString(
      "[mode]\nmode = zoom_in\n", "production_single_rung_default");
  CHUI_REQUIRE(production.config.numerics.hierarchical_max_rung == 0);
  CHUI_REQUIRE(production.normalized_text.find("hierarchical_max_rung = 0") != std::string::npos);
  // Opt-in hierarchical KDK now exists, but only for source-free cosmo_cube DMO.
  requireNumericsRejected(production.normalized_text, "hierarchical_max_rung", "1",
      "requires cosmo_cube collisionless DMO with source physics disabled");
  const auto dmo = loadFrozenConfigFromString(
      "[mode]\nmode = cosmo_cube\nic_file = generated\n"
      "[physics]\nenable_cooling = false\nenable_star_formation = false\nenable_feedback = false\n"
      "enable_stellar_evolution = false\nenable_black_hole_agn = false\nenable_tracers = false\n"
      "enable_metal_diffusion = false\n", "hierarchical_dmo_control");
  CHUI_REQUIRE(dmo.config.numerics.hierarchical_max_rung == 0);
  for (const auto rung : {"1", "12"}) {
    const auto admitted = loadFrozenConfigFromString(
        withNumericsValue(dmo.normalized_text, "hierarchical_max_rung", rung), "hierarchical_dmo_opt_in");
    CHUI_REQUIRE(admitted.config.numerics.hierarchical_max_rung == std::stoi(rung));
    CHUI_REQUIRE(loadFrozenConfigFromString(admitted.normalized_text, "hierarchical_replay").normalized_text == admitted.normalized_text);
  }
  for (const auto rung : {"-1", "13", "1.5", "nan"})
    requireNumericsRejected(dmo.normalized_text, "hierarchical_max_rung", rung);
}

void testPowerSpectrumMeshMatchesCompiledFftBackend() {
  const std::string non_power_of_two =
      "[mode]\nmode = zoom_in\n[analysis]\npower_spectrum_mesh_n = 6\n";
#if COSMOSIM_ENABLE_FFTW
  const auto accepted = cosmosim::core::loadFrozenConfigFromString(
      non_power_of_two, "fftw_non_power_of_two_mesh");
  CHUI_REQUIRE(accepted.config.analysis.power_spectrum_mesh_n == 6);
#else
  bool rejected = false;
  try {
    (void)cosmosim::core::loadFrozenConfigFromString(
        non_power_of_two, "builtin_fft_non_power_of_two_mesh");
  } catch (const cosmosim::core::ConfigError& error) {
    const std::string message = error.what();
    rejected = message.find("analysis.power_spectrum_mesh_n") != std::string::npos &&
        message.find("power of two") != std::string::npos;
  }
  CHUI_REQUIRE(rejected);
#endif

  const auto valid = cosmosim::core::loadFrozenConfigFromString(
      "[mode]\nmode = zoom_in\n[analysis]\npower_spectrum_mesh_n = 8\n",
      "fft_power_of_two_mesh");
  CHUI_REQUIRE(valid.config.analysis.power_spectrum_mesh_n == 8);
}

void testFlatCosmologyClosureContract() {
  const auto fiducial = cosmosim::core::loadFrozenConfigFromString(
      "[mode]\nmode = zoom_in\n[cosmology]\n"
      "omega_matter = 0.3\n"
      "omega_lambda = 0.7\n",
      "flat_lcdm_closure");
  CHUI_REQUIRE(fiducial.config.cosmology.omega_matter == 0.3);
  CHUI_REQUIRE(fiducial.config.cosmology.omega_lambda == 0.7);
  CHUI_REQUIRE(fiducial.normalized_text.find("omega_matter = 0.3") != std::string::npos);
  CHUI_REQUIRE(fiducial.normalized_text.find("omega_lambda = 0.7") != std::string::npos);

  const auto eds = cosmosim::core::loadFrozenConfigFromString(
      "[mode]\nmode = zoom_in\n[cosmology]\n"
      "omega_matter = 1\n"
      "omega_lambda = 0\n",
      "eds_closure");
  CHUI_REQUIRE(eds.config.cosmology.omega_matter == 1.0);
  CHUI_REQUIRE(eds.config.cosmology.omega_lambda == 0.0);

  bool rejected = false;
  try {
    (void)cosmosim::core::loadFrozenConfigFromString(
        "[mode]\nmode = zoom_in\n[cosmology]\n"
        "omega_matter = 0.3\n"
        "omega_lambda = 0\n",
        "inconsistent_cosmology_closure");
  } catch (const cosmosim::core::ConfigError&) {
    rejected = true;
  }
  CHUI_REQUIRE(rejected);
}

void testOutputCodeTimeCadenceValidationAndRoundtrip() {
  const auto frozen = cosmosim::core::loadFrozenConfigFromString(
      "[mode]\nmode = zoom_in\n[output]\n"
      "snapshot_interval_steps = 0\n"
      "snapshot_interval_time_code = 0.125\n"
      "snapshot_layout = single\n"
      "snapshot_num_files = 0\n",
      "output_code_time_cadence");
  CHUI_REQUIRE(frozen.config.output.snapshot_interval_steps == 0);
  CHUI_REQUIRE(frozen.config.output.snapshot_interval_time_code == 0.125);
  CHUI_REQUIRE(frozen.config.output.snapshot_layout ==
         cosmosim::core::ScienceSnapshotLayout::kSingle);
  CHUI_REQUIRE(frozen.config.output.snapshot_num_files == 0U);
  CHUI_REQUIRE(frozen.normalized_text.find("snapshot_interval_time_code = 0.125") !=
         std::string::npos);
  CHUI_REQUIRE(frozen.normalized_text.find("snapshot_layout = single") !=
         std::string::npos);
  const auto reparsed = cosmosim::core::loadFrozenConfigFromString(
      frozen.normalized_text, "output_code_time_cadence_roundtrip");
  CHUI_REQUIRE(reparsed.normalized_text == frozen.normalized_text);

  for (const char* invalid_value : {"0.0", "-0.125", "nan"}) {
    bool rejected = false;
    try {
      (void)cosmosim::core::loadFrozenConfigFromString(
          std::string("[mode]\nmode = zoom_in\n[output]\n") +
              "snapshot_interval_steps = 0\n" +
              "snapshot_interval_time_code = " + invalid_value + "\n",
          "invalid_output_code_time_cadence");
    } catch (const cosmosim::core::ConfigError&) {
      rejected = true;
    }
    CHUI_REQUIRE(rejected);
  }

  bool invalid_single_count_rejected = false;
  try {
    (void)cosmosim::core::loadFrozenConfigFromString(
        "[mode]\nmode = zoom_in\n[output]\n"
        "snapshot_layout = single\n"
        "snapshot_num_files = 2\n",
        "invalid_single_snapshot_file_count");
  } catch (const cosmosim::core::ConfigError&) {
    invalid_single_count_rejected = true;
  }
  CHUI_REQUIRE(invalid_single_count_rejected);

  const auto aggregated = cosmosim::core::loadFrozenConfigFromString(
      "[mode]\nmode = zoom_in\n[output]\n"
      "snapshot_layout = aggregated\n"
      "snapshot_num_files = 8\n",
      "aggregated_snapshot_contract");
  CHUI_REQUIRE(aggregated.config.output.snapshot_layout ==
         cosmosim::core::ScienceSnapshotLayout::kAggregated);
  CHUI_REQUIRE(aggregated.config.output.snapshot_num_files == 8U);
}

}  // namespace

void testTreePmProductionTuningRoundTrip() {
  using cosmosim::core::loadFrozenConfigFromString;
  const auto baseline = loadFrozenConfigFromString(
      "[mode]\nmode = cosmo_cube\nic_file = generated\n"
      "[numerics]\ngravity_solver = treepm\ntreepm_pm_grid = 16\n", "tuning_baseline");
  for (const int block : {16,32,64,128,256}) {
    for (const int leaf : {1,8,1024}) {
      auto text = withNumericsValue(baseline.normalized_text, "treepm_residual_block_size", std::to_string(block));
      text = withNumericsValue(text, "treepm_tree_max_leaf_size", std::to_string(leaf));
      for (const auto key : {"treepm_full_residual_lookup_enabled", "treepm_leaf_simd_enabled",
                            "treepm_radial_coefficients_enabled"}) text = withNumericsValue(text, key, "true");
      const auto frozen = loadFrozenConfigFromString(text, "production_tuning");
      const auto replay = loadFrozenConfigFromString(frozen.normalized_text,"production_tuning_replay");
      CHUI_REQUIRE(replay.config.numerics.treepm_residual_block_size == block);
      CHUI_REQUIRE(replay.config.numerics.treepm_tree_max_leaf_size == leaf);
      CHUI_REQUIRE(replay.config.numerics.treepm_full_residual_lookup_enabled);
      CHUI_REQUIRE(replay.config.numerics.treepm_leaf_simd_enabled);
      CHUI_REQUIRE(replay.config.numerics.treepm_radial_coefficients_enabled);
      CHUI_REQUIRE(replay.normalized_text == frozen.normalized_text);
      CHUI_REQUIRE(replay.provenance.config_hash_hex == frozen.provenance.config_hash_hex);
    }
  }
  for (const auto key : {"treepm_residual_block_size", "treepm_tree_max_leaf_size"}) {
    for (const auto value : {"0", "-1", "4097", "3.5", "nan", "inf", "-inf", "2147483648"})
      requireNumericsRejected(baseline.normalized_text, key, value);
  }
  for (const auto value : {"8", "17", "512"})
    requireNumericsRejected(baseline.normalized_text, "treepm_residual_block_size", value);
  requireNumericsRejected(baseline.normalized_text, "treepm_tree_max_leaf_size", "1025");
  for (const auto value : {"0.25", "0.5", "0.75", "1"}) {
    const auto legal = loadFrozenConfigFromString(withNumericsValue(baseline.normalized_text,
        "treepm_adaptive_maximum_opening_angle", value), "adaptive_boundary_control");
    CHUI_REQUIRE(legal.config.numerics.treepm_adaptive_maximum_opening_angle == std::stod(value));
    CHUI_REQUIRE(loadFrozenConfigFromString(legal.normalized_text, "adaptive_replay").normalized_text == legal.normalized_text);
  }
  for (const auto value : {"0", "-0.1", "1.0001", "nan", "inf", "-inf"})
    requireNumericsRejected(baseline.normalized_text, "treepm_adaptive_maximum_opening_angle", value);

  // Adaptive estimator and cutoff-node policy: defaults preserve the previous
  // behavior, both alternatives round-trip, and unknown spellings fail.
  CHUI_REQUIRE(baseline.config.numerics.treepm_adaptive_error_estimator ==
      cosmosim::core::TreePmAdaptiveErrorEstimator::kResidualDerivative);
  CHUI_REQUIRE(baseline.config.numerics.treepm_cutoff_node_policy ==
      cosmosim::core::TreePmCutoffNodePolicy::kContain);
  CHUI_REQUIRE(baseline.normalized_text.find("treepm_adaptive_error_estimator = residual_derivative") != std::string::npos);
  CHUI_REQUIRE(baseline.normalized_text.find("treepm_cutoff_node_policy = contain") != std::string::npos);
  auto fast_text = withNumericsValue(baseline.normalized_text, "treepm_adaptive_error_estimator", "relative_monopole");
  fast_text = withNumericsValue(fast_text, "treepm_cutoff_node_policy", "com_distance");
  const auto fast = loadFrozenConfigFromString(fast_text, "adaptive_fast_policies");
  CHUI_REQUIRE(fast.config.numerics.treepm_adaptive_error_estimator ==
      cosmosim::core::TreePmAdaptiveErrorEstimator::kRelativeMonopole);
  CHUI_REQUIRE(fast.config.numerics.treepm_cutoff_node_policy ==
      cosmosim::core::TreePmCutoffNodePolicy::kComDistance);
  const auto fast_replay = loadFrozenConfigFromString(fast.normalized_text, "adaptive_fast_replay");
  CHUI_REQUIRE(fast_replay.normalized_text == fast.normalized_text);
  CHUI_REQUIRE(fast_replay.provenance.config_hash_hex == fast.provenance.config_hash_hex);
  CHUI_REQUIRE(fast.provenance.config_hash_hex != baseline.provenance.config_hash_hex);
  for (const auto value : {"monopole", "relative", "1"})
    requireNumericsRejected(baseline.normalized_text, "treepm_adaptive_error_estimator", value);
  for (const auto value : {"com", "containment", "true"})
    requireNumericsRejected(baseline.normalized_text, "treepm_cutoff_node_policy", value);

  // Group walk: off by default, round-trips when its prerequisites hold, and
  // is rejected (not silently ignored) without them.
  CHUI_REQUIRE(!baseline.config.numerics.treepm_group_walk_enabled);
  CHUI_REQUIRE(baseline.normalized_text.find("treepm_group_walk_enabled = false") != std::string::npos);
  auto group_text = withNumericsValue(fast.normalized_text, "treepm_group_walk_enabled", "true");
  group_text = withNumericsValue(group_text, "treepm_adaptive_acceptance_enabled", "true");
  group_text = withNumericsValue(group_text, "treepm_gaussian_pair_lookup_enabled", "true");
  const auto grouped = loadFrozenConfigFromString(group_text, "group_walk_enabled");
  CHUI_REQUIRE(grouped.config.numerics.treepm_group_walk_enabled);
  CHUI_REQUIRE(loadFrozenConfigFromString(grouped.normalized_text, "group_walk_replay").normalized_text ==
      grouped.normalized_text);
  bool group_without_prerequisites_rejected = false;
  try {
    (void)loadFrozenConfigFromString(
        withNumericsValue(baseline.normalized_text, "treepm_group_walk_enabled", "true"), "group_walk_invalid");
  } catch (const cosmosim::core::ConfigError& error) {
    group_without_prerequisites_rejected =
        std::string(error.what()).find("treepm_group_walk_enabled") != std::string::npos;
  }
  CHUI_REQUIRE(group_without_prerequisites_rejected);

  // Hierarchical coarse displacement limiter: legacy default, RMS alternative.
  CHUI_REQUIRE(baseline.config.numerics.hierarchical_coarse_displacement_criterion ==
      cosmosim::core::HierarchicalCoarseDisplacementCriterion::kMaxSpeedMesh);
  CHUI_REQUIRE(baseline.config.numerics.hierarchical_coarse_rms_displacement_factor == 0.25);
  // Encounter rung criterion: on by default, round-trips when disabled.
  CHUI_REQUIRE(baseline.config.numerics.hierarchical_encounter_timestep_enabled);
  CHUI_REQUIRE(baseline.normalized_text.find("hierarchical_encounter_timestep_enabled = true") != std::string::npos);
  {
    const auto no_encounter = loadFrozenConfigFromString(
        withNumericsValue(baseline.normalized_text, "hierarchical_encounter_timestep_enabled", "false"),
        "encounter_disabled");
    CHUI_REQUIRE(!no_encounter.config.numerics.hierarchical_encounter_timestep_enabled);
    CHUI_REQUIRE(loadFrozenConfigFromString(no_encounter.normalized_text, "encounter_replay").normalized_text ==
        no_encounter.normalized_text);
    requireNumericsRejected(baseline.normalized_text, "hierarchical_encounter_timestep_enabled", "maybe");
  }
  auto rms_text = withNumericsValue(baseline.normalized_text,
      "hierarchical_coarse_displacement_criterion", "rms_split_scale");
  rms_text = withNumericsValue(rms_text, "hierarchical_coarse_rms_displacement_factor", "0.5");
  const auto rms = loadFrozenConfigFromString(rms_text, "hierarchical_rms_limiter");
  CHUI_REQUIRE(rms.config.numerics.hierarchical_coarse_displacement_criterion ==
      cosmosim::core::HierarchicalCoarseDisplacementCriterion::kRmsSplitScale);
  CHUI_REQUIRE(rms.config.numerics.hierarchical_coarse_rms_displacement_factor == 0.5);
  CHUI_REQUIRE(loadFrozenConfigFromString(rms.normalized_text, "hierarchical_rms_replay").normalized_text ==
      rms.normalized_text);
  for (const auto value : {"rms", "max", "1"})
    requireNumericsRejected(baseline.normalized_text, "hierarchical_coarse_displacement_criterion", value);
  for (const auto value : {"0", "-0.1", "1.5", "nan", "inf"})
    requireNumericsRejected(baseline.normalized_text, "hierarchical_coarse_rms_displacement_factor", value);
}

void testTreePmPerformanceFlagsRoundTrip() {
  using cosmosim::core::loadFrozenConfigFromString;
  const auto defaults = loadFrozenConfigFromString(
      "[mode]\nmode = cosmo_cube\nic_file = generated\n"
      "[numerics]\ngravity_solver = treepm\ntreepm_pm_grid = 16\n", "flags_defaults");
  CHUI_REQUIRE(!defaults.config.numerics.treepm_gaussian_pair_lookup_enabled);
  CHUI_REQUIRE(defaults.config.numerics.treepm_full_mac_diagnostics);
  CHUI_REQUIRE(!defaults.config.numerics.treepm_full_residual_lookup_enabled);
  CHUI_REQUIRE(!defaults.config.numerics.treepm_leaf_simd_enabled);
  CHUI_REQUIRE(!defaults.config.numerics.treepm_radial_coefficients_enabled);
  CHUI_REQUIRE(!defaults.config.numerics.treepm_adaptive_acceptance_enabled);
  CHUI_REQUIRE(!defaults.config.numerics.treepm_identical_source_tree_reuse_enabled);
  CHUI_REQUIRE(!defaults.config.numerics.treepm_topology_refit_enabled);
  for (const std::string key : {"treepm_gaussian_pair_lookup_enabled", "treepm_full_mac_diagnostics",
       "treepm_full_residual_lookup_enabled", "treepm_leaf_simd_enabled", "treepm_radial_coefficients_enabled",
       "treepm_adaptive_acceptance_enabled", "treepm_identical_source_tree_reuse_enabled",
       "treepm_topology_refit_enabled", "treepm_enable_window_deconvolution"}) {
    for (const auto value : {"true", "false"}) {
      const auto frozen = loadFrozenConfigFromString(
          withNumericsValue(defaults.normalized_text, key, value), "treepm_performance_flags");
      const auto replay = loadFrozenConfigFromString(frozen.normalized_text,"flags_replay");
      CHUI_REQUIRE(frozen.normalized_text.find(key + " = " + value + "\n") != std::string::npos);
      CHUI_REQUIRE(replay.normalized_text == frozen.normalized_text);
      CHUI_REQUIRE(replay.provenance.config_hash_hex == frozen.provenance.config_hash_hex);
    }
    for (const auto value : {"maybe", "2", "nan", "inf", "yes", "on", "1", "0"})
      requireNumericsRejected(defaults.normalized_text, key, value);
  }
}

int main() {
  testTreePmProductionTuningRoundTrip();
  testTreePmPerformanceFlagsRoundTrip();
  testCommentsWhitespaceSectionsAndUnits();
  testDuplicateKeyFails();
  testRequiredSectionAndFieldFailures();
  testMissingValueFails();
  testUnknownKeysFailUnlessCompatibilityEnabled();
  testAliasConflictFails();
  testInvalidEnumFails();
  testInvalidTypedPolicyEnumsFail();
  testTypedInitialConditionConfiguration();
  testBoundaryModeValidation();
  testDefaultsCanonicalizationAndDeterminism();
  testNormalizedConfigRoundTripPreservesDifficultDoubles();
  testDebugExactDecompositionAuditConfig();
  testConfigNormalizationHashDeterminism();
  testConfigRuntimeOwnership();
  testDerivedRuntimeConfigOwnership();
  testUrlsWindowsPathsAndQuotedHashesAreNotTruncated();
  testDeprecatedAliasesAndCanonicalCollision();
  testFeedbackConfigKeysAndValidation();
  testMetalDiffusionConfigurationContract();
  testCoolingPolicyEnumsAndValidation();
  testDiagnosticsExecutionPolicyValidation();
  testEnumSerializationIsFailFastWithoutUnknownFallback();
  testSpeciesSofteningOverridesRoundtrip();
  testHydroActiveBatchConfigurationRoundtrip();
  testTreePmNumericsRoundtripAndValidation();
  testAxisAwareBoxAndPmGridCanonicalizationCompatibility();
  testBlackHoleAgnConfigKeysAndValidation();
  testTracerConfigKeysAndValidation();
  testZoomLongRangeStrategyValidationAndRoundtrip();
  testCosmologyScaleFactorRedshiftCanonicalizationAndValidation();
  testFlatCosmologyClosureContract();
  testIntegratorTimeVariableIsTypedAndCanonical();
  testIntegratorEndpointAuthorityResolution();
  testAdversarialPhysicsAndCosmologyDependenciesFail();
  testFiniteNumericAndForwardCosmologyContract();
  testDerivedRuntimeSerializationUsesCanonicalNames();
  testStarFormationMassParametersAcceptUnits();
  testEffectiveMultiphaseConfigurationContract();
  testProductionHierarchicalRungsFailClosed();
  testPowerSpectrumMeshMatchesCompiledFftBackend();
  testOutputCodeTimeCadenceValidationAndRoundtrip();
  return 0;
}
