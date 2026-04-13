#pragma once

#include <array>
#include <string>
#include <string_view>

#include "cosmosim/amr/amr_module.hpp"
#include "cosmosim/amr/amr_framework.hpp"
#include "cosmosim/analysis/analysis_module.hpp"
#include "cosmosim/analysis/diagnostics.hpp"
#include "cosmosim/analysis/halo_workflow.hpp"
#include "cosmosim/core/config.hpp"
#include "cosmosim/core/cosmology.hpp"
#include "cosmosim/core/provenance.hpp"
#include "cosmosim/core/profiling.hpp"
#include "cosmosim/core/soa_storage.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/core/simulation_mode.hpp"
#include "cosmosim/core/time_integration.hpp"
#include "cosmosim/core/units.hpp"
#include "cosmosim/core/module_registry.hpp"
#include "cosmosim/core/version.hpp"
#include "cosmosim/gravity/gravity_module.hpp"
#include "cosmosim/gravity/pm_solver.hpp"
#include "cosmosim/gravity/tree_gravity.hpp"
#include "cosmosim/gravity/tree_pm_coupling.hpp"
#include "cosmosim/hydro/hydro_core_solver.hpp"
#include "cosmosim/hydro/hydro_module.hpp"
#include "cosmosim/io/ic_reader.hpp"
#include "cosmosim/io/io_module.hpp"
#include "cosmosim/io/snapshot_hdf5.hpp"
#include "cosmosim/io/restart_checkpoint.hpp"
#include "cosmosim/parallel/parallel_module.hpp"
#include "cosmosim/physics/physics_module.hpp"
#include "cosmosim/physics/cooling_heating.hpp"
#include "cosmosim/physics/star_formation.hpp"
#include "cosmosim/physics/black_hole_agn.hpp"
#include "cosmosim/physics/stellar_evolution.hpp"
#include "cosmosim/physics/stellar_feedback.hpp"
#include "cosmosim/physics/tracer_support.hpp"
#include "cosmosim/utils/utils_module.hpp"
#include "cosmosim/workflows/reference_workflow.hpp"

namespace cosmosim {

constexpr std::array<std::string_view, 8> k_layer_order = {
    "core", "gravity", "hydro", "amr", "physics", "io", "analysis", "parallel"};

std::string architectureSummary();

}  // namespace cosmosim
