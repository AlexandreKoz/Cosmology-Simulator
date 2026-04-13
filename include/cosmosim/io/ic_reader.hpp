#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

#include "cosmosim/core/config.hpp"
#include "cosmosim/core/simulation_state.hpp"

namespace cosmosim::io {

enum class IcParticleKind : std::uint8_t {
  kGas = 0,
  kDarkMatter = 1,
  kDisk = 2,
  kBulge = 3,
  kStar = 4,
  kBoundary = 5,
};

// A compact per-file schema summary for transparent import auditing.
struct IcSchemaSummary {
  std::array<std::uint64_t, 6> count_by_type{};
  std::array<double, 6> mass_table{};
  double scale_factor = 1.0;
  bool velocities_are_peculiar = true;
};

struct IcImportOptions {
  bool require_velocities = true;
  bool require_particle_ids = true;
  bool allow_mass_table_fallback = true;
  std::size_t chunk_particle_count = 1u << 16;
};

struct IcImportReport {
  IcSchemaSummary schema;
  std::vector<std::string> present_aliases;
  std::vector<std::string> defaulted_fields;
  std::vector<std::string> missing_optional_fields;
  std::vector<std::string> unsupported_fields;
};

struct IcReadResult {
  core::SimulationState state;
  IcImportReport report;
};

// Read common GADGET/AREPO-style HDF5 initial conditions into SimulationState.
// This API is species-aware and writes directly into SoA buffers in chunked batches.
[[nodiscard]] IcReadResult readGadgetArepoHdf5Ic(
    const std::filesystem::path& ic_path,
    const core::SimulationConfig& config,
    const IcImportOptions& options = {});

// Internal helper path for generated isolated test problems without external files.
[[nodiscard]] IcReadResult buildGeneratedIsolatedIc(
    const core::SimulationConfig& config,
    std::size_t dark_matter_particle_count,
    std::size_t gas_particle_count,
    std::uint64_t particle_id_seed = 1u);

// A small deterministic converter for compatibility tooling and regression tests.
[[nodiscard]] IcReadResult convertGeneratedIsolatedIcToState(
    const core::SimulationConfig& config,
    std::size_t particle_count_per_axis);

}  // namespace cosmosim::io
