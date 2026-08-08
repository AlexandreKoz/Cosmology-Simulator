#pragma once

#include <array>
#include <span>
#include <string_view>

namespace cosmosim::core {

constexpr std::array<std::string_view, 9> k_module_names = {
    "core", "gravity", "hydro", "amr", "physics", "io", "analysis", "parallel", "utils"};

// Compile-time composition metadata. Semicolon-delimited values name concrete
// runtime contributions; "none" is explicit and never means unknown.
struct ModuleDescriptor {
  std::string_view name;
  std::string_view config_fragment;
  std::string_view state_and_sidecar_requirements;
  std::string_view stage_and_task_contributions;
  std::string_view timestep_criteria;
  std::string_view restart_payloads;
  std::string_view migration_fields;
  std::string_view diagnostics;
};

inline constexpr std::array<ModuleDescriptor, k_module_names.size()>
    k_module_descriptors{{
        {"core", "schema;numerics;cosmology;units;mode", "particle_soa;cell_soa;scheduler_truth;integrator_state", "kdk_plan;drift;time_commit", "user_clamp", "integrator;scheduler", "particle_core;gas_cell_identity", "time_bins;memory;runtime_events"},
        {"gravity", "numerics.treepm_*;numerics.gravity_*", "particle_mass_position;force_cache;pm_sync_state", "gravity_kick_pre;force_refresh;gravity_kick_post", "gravity_acceleration", "force_cache;distributed_gravity", "force_cache_invalidated_on_migration", "treepm;force_health"},
        {"hydro", "physics.hydro_*;mode.hydro_boundary", "gas_cell_identity;hydro_primitive_conserved", "hydro_update;ghost_refresh", "hydro_cfl", "gas_cells;hydro_state", "gas_cell_fields;flux_corrections", "hydro_cfl;conservation"},
        {"amr", "physics.amr_*", "patches;gas_cells;flux_registers;temporal_ghosts", "amr_hydro;ghost_fill;reflux", "hydro_cfl;amr_level", "patches;flux_registers", "patch_payload;cell_payload", "amr_conservation;ghost_epochs"},
        {"physics", "physics.*", "star_sidecar;black_hole_sidecar;tracer_sidecar", "source_terms", "source_rate", "source_sidecars;module_sidecars", "species_sidecars", "source_conservation"},
        {"io", "output.*;mode.ic_file", "schema_maps;provenance", "initial_conditions;output_boundary;restart", "none", "all_authoritative_runtime_state", "none", "roundtrip_verification"},
        {"analysis", "analysis.*", "read_only_diagnostic_views", "analysis_hooks", "none", "diagnostic_policy", "none", "run_health;science_diagnostics"},
        {"parallel", "parallel.*", "ownership;ghost_epoch;decomposition_epoch", "collectives;ghost_exchange;migration", "none", "distributed_topology", "all_authoritative_migration_records", "communication;load_balance"},
        {"utils", "none", "none", "none", "none", "none", "none", "none"},
    }};

std::array<std::string_view, k_module_names.size()> moduleNames();
[[nodiscard]] std::span<const ModuleDescriptor> moduleDescriptors() noexcept;
[[nodiscard]] const ModuleDescriptor& requireModuleDescriptor(
    std::string_view name);

}  // namespace cosmosim::core
