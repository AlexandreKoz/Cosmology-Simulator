#include "cosmosim/io/restart_checkpoint.hpp"

#include <algorithm>
#include <bit>
#include <array>
#include <cmath>
#include <cerrno>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>
#include <utility>
#include <unordered_set>
#include <span>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/io/io_contract.hpp"

#if COSMOSIM_ENABLE_HDF5
#include <fcntl.h>
#include <hdf5.h>
#include <unistd.h>
#endif

namespace cosmosim::io {
namespace {

constexpr std::uint64_t k_offset_basis = 14695981039346656037ull;
constexpr std::uint64_t k_prime = 1099511628211ull;
constexpr std::uint32_t k_restart_schema_v14 = 14;
constexpr std::uint32_t k_restart_schema_v15 = 15;
constexpr std::uint32_t k_restart_schema_v16 = 16;
constexpr std::uint32_t k_restart_schema_v17 = 17;
constexpr std::uint32_t k_restart_schema_v18 = 18;
constexpr std::uint32_t k_restart_schema_v19 = 19;
constexpr std::uint32_t k_restart_schema_v20 = 20;
constexpr std::string_view k_restart_schema_name_v14 = "cosmosim_restart_v14";
constexpr std::string_view k_restart_schema_name_v15 = "cosmosim_restart_v15";
constexpr std::string_view k_restart_schema_name_v16 = "cosmosim_restart_v16";
constexpr std::string_view k_restart_schema_name_v17 = "cosmosim_restart_v17";
constexpr std::string_view k_restart_schema_name_v18 = "cosmosim_restart_v18";
constexpr std::string_view k_restart_schema_name_v19 = "cosmosim_restart_v19";
constexpr std::string_view k_gas_identity_row_policy = "explicit_dense_local_cell_row";
constexpr std::string_view k_gas_cell_scheduler_identity_key = "gas_cell_id";

[[nodiscard]] std::string hexU64(std::uint64_t value) {
  constexpr char k_hex[] = "0123456789abcdef";
  std::string text(16, '0');
  for (int i = 15; i >= 0; --i) {
    text[static_cast<std::size_t>(i)] = k_hex[value & 0xF];
    value >>= 4;
  }
  return text;
}

[[nodiscard]] std::uint64_t fnv1aAppend(std::uint64_t seed, std::span<const std::byte> bytes) {
  std::uint64_t hash = seed;
  for (const std::byte value : bytes) {
    hash ^= static_cast<std::uint64_t>(std::to_integer<unsigned char>(value));
    hash *= k_prime;
  }
  return hash;
}

template <typename VectorLike>
[[nodiscard]] std::span<const std::byte> asBytesSpan(const VectorLike& values) {
  using ValueType = typename VectorLike::value_type;
  return {reinterpret_cast<const std::byte*>(values.data()), values.size() * sizeof(ValueType)};
}

void validateSchedulerPersistentStateForRestart(
    const core::TimeBinPersistentState& scheduler_state,
    std::size_t expected_element_count,
    std::string_view context,
    std::string_view scheduler_path) {
  if (scheduler_state.bin_index.size() != scheduler_state.next_activation_tick.size() ||
      scheduler_state.bin_index.size() != scheduler_state.active_flag.size() ||
      scheduler_state.bin_index.size() != scheduler_state.pending_bin_index.size()) {
    throw std::invalid_argument(
        std::string(context) + ": " + std::string(scheduler_path) +
        " persistent arrays must have matching sizes");
  }
  if (scheduler_state.bin_index.size() != expected_element_count) {
    throw std::invalid_argument(
        std::string(context) + ": " + std::string(scheduler_path) +
        " element count does not match the authoritative state extent");
  }
  for (std::size_t i = 0; i < scheduler_state.bin_index.size(); ++i) {
    if (scheduler_state.bin_index[i] > scheduler_state.max_bin) {
      throw std::invalid_argument(
          std::string(context) + ": " + std::string(scheduler_path) + "/bin_index exceeds max_bin");
    }
    if (scheduler_state.active_flag[i] > 1U) {
      throw std::invalid_argument(
          std::string(context) + ": " + std::string(scheduler_path) + "/active_flag must contain only 0 or 1");
    }
    const bool pending_is_unset =
        scheduler_state.pending_bin_index[i] == core::HierarchicalTimeBinScheduler::k_unset_pending_bin;
    if (!pending_is_unset && scheduler_state.pending_bin_index[i] > scheduler_state.max_bin) {
      throw std::invalid_argument(
          std::string(context) + ": " + std::string(scheduler_path) + "/pending_bin_index exceeds max_bin");
    }
  }
}

void validateGravityForceCacheForRestart(
    const GravityForceCachePersistentState& cache,
    const core::SimulationState& state,
    std::string_view context) {
  const auto validate_triplet = [&](const std::vector<double>& x,
                                    const std::vector<double>& y,
                                    const std::vector<double>& z,
                                    std::size_t expected,
                                    std::string_view label) {
    if (x.size() != y.size() || x.size() != z.size()) {
      throw std::invalid_argument(std::string(context) + ": gravity force cache " +
                                  std::string(label) + " component extents differ");
    }
    if (cache.valid && x.size() != expected) {
      throw std::invalid_argument(std::string(context) + ": gravity force cache " +
                                  std::string(label) + " extent does not match authoritative state");
    }
    if (!cache.valid && !x.empty()) {
      throw std::invalid_argument(std::string(context) +
                                  ": invalid gravity force cache must not retain stale component lanes");
    }
    for (const double value : x) {
      if (!std::isfinite(value)) {
        throw std::invalid_argument(std::string(context) + ": gravity force cache contains non-finite values");
      }
    }
    for (const double value : y) {
      if (!std::isfinite(value)) {
        throw std::invalid_argument(std::string(context) + ": gravity force cache contains non-finite values");
      }
    }
    for (const double value : z) {
      if (!std::isfinite(value)) {
        throw std::invalid_argument(std::string(context) + ": gravity force cache contains non-finite values");
      }
    }
  };
  const auto validate_identity = [&](const std::vector<std::uint64_t>& cache_ids,
                                     std::span<const std::uint64_t> state_ids,
                                     std::string_view label) {
    if (!cache.valid && !cache_ids.empty()) {
      throw std::invalid_argument(std::string(context) + ": invalid gravity force cache must not retain stale " +
                                  std::string(label) + " identity lanes");
    }
    if (!cache.valid) {
      return;
    }
    if (cache_ids.size() != state_ids.size()) {
      throw std::invalid_argument(std::string(context) + ": gravity force cache " +
                                  std::string(label) + " identity extent does not match authoritative state");
    }
    std::unordered_set<std::uint64_t> seen;
    seen.reserve(cache_ids.size());
    for (std::size_t index = 0; index < cache_ids.size(); ++index) {
      if (cache_ids[index] == 0U || !seen.insert(cache_ids[index]).second) {
        throw std::invalid_argument(std::string(context) + ": gravity force cache " +
                                    std::string(label) + " identity lane contains a zero or duplicate ID");
      }
      if (cache_ids[index] != state_ids[index]) {
        throw std::invalid_argument(std::string(context) + ": gravity force cache " +
                                    std::string(label) + " identity lane does not match authoritative dense row");
      }
    }
  };
  validate_triplet(cache.particle_accel_x_comoving, cache.particle_accel_y_comoving,
                   cache.particle_accel_z_comoving, state.particles.size(), "particles");
  validate_triplet(cache.cell_accel_x_comoving, cache.cell_accel_y_comoving,
                   cache.cell_accel_z_comoving, state.cells.size(), "gas cells");
  validate_identity(cache.particle_id, state.particle_sidecar.particle_id, "particle");
  validate_identity(cache.gas_cell_id, state.gas_cells.gas_cell_id, "gas-cell");
}

void validateRestartTimeBinMirrorsAgainstScheduler(
    const core::SimulationState& state,
    const core::TimeBinPersistentState& scheduler_state,
    std::string_view context) {
  validateSchedulerPersistentStateForRestart(
      scheduler_state, state.particles.size(), context, "/scheduler");
  if (state.particles.time_bin.size() != scheduler_state.bin_index.size()) {
    throw std::invalid_argument(
        std::string(context) + ": /state/particles/time_bin length must match /scheduler/bin_index length");
  }
  for (std::size_t i = 0; i < state.particles.time_bin.size(); ++i) {
    if (state.particles.time_bin[i] != scheduler_state.bin_index[i]) {
      throw std::invalid_argument(
          std::string(context) + ": /state/particles/time_bin is stale relative to /scheduler/bin_index");
    }
  }
}

void validateGasCellTimeBinMirrorsAgainstScheduler(
    const core::SimulationState& state,
    const core::TimeBinPersistentState& scheduler_state,
    std::span<const std::uint64_t> gas_cell_ids,
    std::string_view context) {
  validateSchedulerPersistentStateForRestart(
      scheduler_state, state.cells.size(), context, "/gas_cell_scheduler");
  state.requireGasCellIdentityMapCoversDenseRows(context);
  if (state.cells.time_bin.size() != scheduler_state.bin_index.size()) {
    throw std::invalid_argument(
        std::string(context) + ": /state/cells/time_bin length must match /gas_cell_scheduler/bin_index length");
  }
  if (gas_cell_ids.size() != state.cells.size()) {
    throw std::invalid_argument(
        std::string(context) + ": /gas_cell_scheduler/gas_cell_id length must match CellSoa extent");
  }
  for (std::uint32_t cell_index = 0; cell_index < state.cells.size(); ++cell_index) {
    const core::GasCellIdentityRecord* identity = state.gas_cell_identity.findByLocalRow(cell_index);
    if (identity == nullptr) {
      throw std::invalid_argument(std::string(context) + ": gas-cell identity map is missing a dense row");
    }
    if (gas_cell_ids[cell_index] != identity->gas_cell_id) {
      throw std::invalid_argument(
          std::string(context) + ": /gas_cell_scheduler/gas_cell_id is stale relative to GasCellIdentityMap");
    }
    if (state.cells.time_bin[cell_index] != scheduler_state.bin_index[cell_index]) {
      throw std::invalid_argument(
          std::string(context) + ": /state/cells/time_bin is stale relative to /gas_cell_scheduler/bin_index");
    }
  }
}

void rebuildRestartTimeBinMirrorsFromScheduler(
    core::SimulationState& state,
    const core::TimeBinPersistentState& scheduler_state) {
  if (state.particles.time_bin.size() != scheduler_state.bin_index.size()) {
    throw std::invalid_argument("restart reader: particle time_bin mirror length must match scheduler bin_index length");
  }
  for (std::size_t i = 0; i < state.particles.time_bin.size(); ++i) {
    state.particles.time_bin[i] = scheduler_state.bin_index[i];
  }
}

void rebuildGasCellTimeBinMirrorsFromScheduler(
    core::SimulationState& state,
    const core::TimeBinPersistentState& scheduler_state,
    std::span<const std::uint64_t> gas_cell_ids) {
  validateGasCellTimeBinMirrorsAgainstScheduler(state, scheduler_state, gas_cell_ids, "restart reader");
  for (std::size_t i = 0; i < state.cells.time_bin.size(); ++i) {
    state.cells.time_bin[i] = scheduler_state.bin_index[i];
  }
}

[[nodiscard]] std::vector<std::uint64_t> gasCellIdsByDenseLocalRow(
    const core::SimulationState& state,
    std::string_view context) {
  state.requireGasCellIdentityMapCoversDenseRows(context);
  std::vector<std::uint64_t> gas_cell_ids(state.cells.size(), 0U);
  for (std::uint32_t cell_index = 0; cell_index < state.cells.size(); ++cell_index) {
    const core::GasCellIdentityRecord* identity = state.gas_cell_identity.findByLocalRow(cell_index);
    if (identity == nullptr) {
      throw std::runtime_error(std::string(context) + ": gas-cell identity map is missing a dense row");
    }
    gas_cell_ids[cell_index] = identity->gas_cell_id;
  }
  return gas_cell_ids;
}

[[nodiscard]] core::TimeBinPersistentState legacyGasCellSchedulerStateFromMirrors(
    const core::SimulationState& state,
    const core::TimeBinPersistentState& particle_scheduler_state) {
  // Compatibility only for older direct API callers.  The reference workflow
  // always supplies a dedicated gas-cell scheduler; v19 then persists that
  // authoritative state explicitly.  This branch never maps cells through
  // parent particles and therefore remains valid for parentless cells.
  core::TimeBinPersistentState gas_scheduler_state;
  gas_scheduler_state.current_tick = particle_scheduler_state.current_tick;
  gas_scheduler_state.max_bin = particle_scheduler_state.max_bin;
  gas_scheduler_state.bin_index.assign(state.cells.time_bin.begin(), state.cells.time_bin.end());
  gas_scheduler_state.next_activation_tick.resize(state.cells.size(), particle_scheduler_state.current_tick);
  gas_scheduler_state.active_flag.assign(state.cells.size(), 0U);
  gas_scheduler_state.pending_bin_index.assign(
      state.cells.size(), core::HierarchicalTimeBinScheduler::k_unset_pending_bin);
  for (std::size_t cell_index = 0; cell_index < state.cells.size(); ++cell_index) {
    const std::uint8_t bin = gas_scheduler_state.bin_index[cell_index];
    if (bin > gas_scheduler_state.max_bin) {
      throw std::invalid_argument("restart legacy gas-cell scheduler mirror has a bin above max_bin");
    }
    const std::uint64_t interval = 1ULL << bin;
    gas_scheduler_state.next_activation_tick[cell_index] =
        ((particle_scheduler_state.current_tick / interval) + 1U) * interval;
  }
  return gas_scheduler_state;
}

void validateHydroGeometryStateForRestart(
    const core::SimulationState& state,
    std::string_view context) {
  if (state.cells.size() != state.gas_cells.size()) {
    throw std::invalid_argument(std::string(context) + ": /state/cells and /state/gas_cells sizes must match");
  }
  if (!state.cells.isConsistent()) {
    throw std::invalid_argument(std::string(context) + ": /state/cells lanes must have matching lengths");
  }
  if (!state.gas_cells.isConsistent()) {
    throw std::invalid_argument(std::string(context) + ": /state/gas_cells lanes must have matching lengths");
  }
  if (!state.patches.isConsistent()) {
    throw std::invalid_argument(std::string(context) + ": /state/patches lanes must have matching lengths");
  }
  if (state.cells.size() == 0) {
    if (state.patches.size() != 0) {
      throw std::invalid_argument(std::string(context) + ": /state/patches must be empty when no gas cells exist");
    }
    return;
  }
  // PatchSoa is mandatory for AMR-backed gas geometry, but legacy/non-AMR
  // Cartesian states may intentionally have no patch descriptors.  In that
  // case `patch_index == 0` and `owning_patch_id == 0` are sentinel metadata;
  // they are not a reason to reject a restart whose cell identity and sidecar
  // lanes are otherwise complete.  Do not manufacture a fake patch merely to
  // satisfy restart I/O.
  if (state.patches.size() == 0U) {
    for (std::size_t cell_index = 0; cell_index < state.cells.size(); ++cell_index) {
      if (state.cells.patch_index[cell_index] != 0U) {
        throw std::invalid_argument(
            std::string(context) + ": /state/cells/patch_index must be zero when PatchSoa is absent");
      }
    }
  } else {
    std::vector<std::uint8_t> cell_covered_by_patch(state.cells.size(), 0U);
    for (std::size_t patch_index = 0; patch_index < state.patches.size(); ++patch_index) {
      const std::uint64_t first_cell = state.patches.first_cell[patch_index];
      const std::uint64_t cell_count = state.patches.cell_count[patch_index];
      if (cell_count == 0U) {
        throw std::invalid_argument(std::string(context) + ": /state/patches/cell_count entries must be non-zero");
      }
      if (first_cell > state.cells.size() || cell_count > state.cells.size() - first_cell) {
        throw std::invalid_argument(std::string(context) + ": /state/patches cell range exceeds /state/cells");
      }
      for (std::uint64_t offset = 0; offset < cell_count; ++offset) {
        const std::size_t cell_index = static_cast<std::size_t>(first_cell + offset);
        if (cell_covered_by_patch[cell_index] != 0U) {
          throw std::invalid_argument(std::string(context) + ": /state/patches cell ranges overlap");
        }
        cell_covered_by_patch[cell_index] = 1U;
        if (state.cells.patch_index[cell_index] != patch_index) {
          throw std::invalid_argument(
              std::string(context) + ": /state/cells/patch_index does not match /state/patches cell ranges");
        }
      }
    }
    for (std::size_t cell_index = 0; cell_index < cell_covered_by_patch.size(); ++cell_index) {
      if (cell_covered_by_patch[cell_index] == 0U) {
        throw std::invalid_argument(std::string(context) + ": /state/patches do not cover every gas cell");
      }
    }
  }
  if (!state.gas_cell_identity.isConsistent() ||
      !state.gas_cell_identity.coversDenseLocalRows(state.cells.size())) {
    throw std::invalid_argument(
        std::string(context) + ": /state/gas_cell_identity must cover dense local cell rows");
  }
  if (!state.gasCellIdentityMapMatchesSidecarLanes()) {
    throw std::invalid_argument(
        std::string(context) + ": /state/gas_cell_identity records must match gas-cell mirror lanes and patch ownership");
  }
  for (const core::GasCellIdentityRecord& record : state.gas_cell_identity.records()) {
    if (record.parent_particle_id.has_value() && *record.parent_particle_id == 0U) {
      throw std::invalid_argument(
          std::string(context) + ": /state/gas_cell_identity/parent_particle_id must be nonzero when has_parent_particle is true");
    }
  }
}

void validateOutputCadenceStateForRestart(
    const OutputCadencePersistentState& output_state,
    const core::IntegratorState& integrator_state,
    std::string_view context) {
  if (!output_state.output_enabled) {
    return;
  }
  if (output_state.snapshot_interval_steps == 0) {
    throw std::invalid_argument(
        std::string(context) + ": /output_cadence/@snapshot_interval_steps must be > 0 when output is enabled");
  }
  if (output_state.last_completed_step_index != integrator_state.step_index) {
    throw std::invalid_argument(
        std::string(context) + ": /output_cadence/@last_completed_step_index must match /integrator/@step_index");
  }
  if (output_state.next_snapshot_step_index != 0 &&
      output_state.next_snapshot_step_index <= output_state.last_completed_step_index) {
    throw std::invalid_argument(
        std::string(context) + ": /output_cadence/@next_snapshot_step_index must be strictly after the restart step");
  }
}

void validateStochasticStateForRestart(
    const StochasticPersistentState& stochastic_state,
    const core::IntegratorState& integrator_state,
    std::string_view context) {
  std::vector<std::string> seen_modules;
  seen_modules.reserve(stochastic_state.modules.size());
  for (const StochasticModulePersistentState& module_state : stochastic_state.modules) {
    if (module_state.module_name.empty()) {
      throw std::invalid_argument(std::string(context) + ": /stochastic_state module_name must be non-empty");
    }
    if (module_state.module_name.find('/') != std::string::npos) {
      throw std::invalid_argument(
          std::string(context) + ": /stochastic_state/" + module_state.module_name +
          " module_name must not contain '/'");
    }
    if (std::find(seen_modules.begin(), seen_modules.end(), module_state.module_name) != seen_modules.end()) {
      throw std::invalid_argument(
          std::string(context) + ": /stochastic_state duplicate module_name '" + module_state.module_name + "'");
    }
    seen_modules.push_back(module_state.module_name);
    if (module_state.schema_version == 0) {
      throw std::invalid_argument(
          std::string(context) + ": /stochastic_state/" + module_state.module_name + "/@schema_version must be > 0");
    }
    if (module_state.rng_policy.empty()) {
      throw std::invalid_argument(
          std::string(context) + ": /stochastic_state/" + module_state.module_name + "/@rng_policy must be non-empty");
    }
    if (!module_state.deterministic_from_serialized_inputs) {
      throw std::invalid_argument(
          std::string(context) + ": /stochastic_state/" + module_state.module_name +
          " declares non-deterministic RNG state, but stateful RNG engine serialization is not implemented");
    }
    if (module_state.last_committed_step_index != integrator_state.step_index) {
      throw std::invalid_argument(
          std::string(context) + ": /stochastic_state/" + module_state.module_name +
          "/@last_committed_step_index must match /integrator/@step_index");
    }
  }
}

[[nodiscard]] std::vector<StochasticModulePersistentState> sortedStochasticModules(
    const StochasticPersistentState& stochastic_state) {
  std::vector<StochasticModulePersistentState> modules = stochastic_state.modules;
  std::sort(
      modules.begin(),
      modules.end(),
      [](const StochasticModulePersistentState& lhs, const StochasticModulePersistentState& rhs) {
        return lhs.module_name < rhs.module_name;
      });
  return modules;
}

[[nodiscard]] std::vector<core::GasCellIdentityRecord> gasIdentityRecordsSortedByLocalRow(
    const core::GasCellIdentityMap& identity_map) {
  std::vector<core::GasCellIdentityRecord> records(
      identity_map.records().begin(), identity_map.records().end());
  std::sort(
      records.begin(),
      records.end(),
      [](const core::GasCellIdentityRecord& lhs, const core::GasCellIdentityRecord& rhs) {
        return lhs.local_cell_row < rhs.local_cell_row;
      });
  return records;
}

#if COSMOSIM_ENABLE_HDF5

class Hdf5Handle {
 public:
  explicit Hdf5Handle(hid_t handle = -1) : m_handle(handle) {}
  Hdf5Handle(const Hdf5Handle&) = delete;
  Hdf5Handle& operator=(const Hdf5Handle&) = delete;
  Hdf5Handle(Hdf5Handle&& other) noexcept : m_handle(other.m_handle) { other.m_handle = -1; }
  Hdf5Handle& operator=(Hdf5Handle&& other) noexcept {
    if (this != &other) {
      close();
      m_handle = other.m_handle;
      other.m_handle = -1;
    }
    return *this;
  }
  ~Hdf5Handle() { close(); }

  [[nodiscard]] hid_t get() const { return m_handle; }
  [[nodiscard]] bool valid() const { return m_handle >= 0; }

 private:
  void close() {
    if (m_handle < 0) {
      return;
    }
    const H5I_type_t type = H5Iget_type(m_handle);
    if (type == H5I_FILE) {
      H5Fclose(m_handle);
    } else if (type == H5I_GROUP) {
      H5Gclose(m_handle);
    } else if (type == H5I_DATASET) {
      H5Dclose(m_handle);
    } else if (type == H5I_DATASPACE) {
      H5Sclose(m_handle);
    } else if (type == H5I_ATTR) {
      H5Aclose(m_handle);
    } else if (type == H5I_DATATYPE) {
      H5Tclose(m_handle);
    }
    m_handle = -1;
  }

  hid_t m_handle = -1;
};

void writeScalarStringAttribute(hid_t location, std::string_view key, const std::string& value) {
  Hdf5Handle scalar_space(H5Screate(H5S_SCALAR));
  Hdf5Handle string_type(H5Tcopy(H5T_C_S1));
  if (!scalar_space.valid() || !string_type.valid()) {
    throw std::runtime_error("failed to create string attribute type");
  }
  if (H5Tset_size(string_type.get(), std::max<std::size_t>(std::size_t{1}, value.size() + 1)) < 0 || H5Tset_strpad(string_type.get(), H5T_STR_NULLTERM) < 0) {
    throw std::runtime_error("failed to configure string attribute type");
  }
  Hdf5Handle attr(H5Acreate2(location, std::string(key).c_str(), string_type.get(), scalar_space.get(), H5P_DEFAULT, H5P_DEFAULT));
  if (!attr.valid() || H5Awrite(attr.get(), string_type.get(), value.c_str()) < 0) {
    throw std::runtime_error("failed to write string attribute: " + std::string(key));
  }
}

[[nodiscard]] std::string readScalarStringAttribute(hid_t location, std::string_view key) {
  Hdf5Handle attr(H5Aopen(location, std::string(key).c_str(), H5P_DEFAULT));
  if (!attr.valid()) {
    throw std::runtime_error("missing string attribute: " + std::string(key));
  }
  Hdf5Handle attr_type(H5Aget_type(attr.get()));
  const std::size_t length = H5Tget_size(attr_type.get());
  std::string value(length, '\0');
  if (H5Aread(attr.get(), attr_type.get(), value.data()) < 0) {
    throw std::runtime_error("failed reading string attribute: " + std::string(key));
  }
  while (!value.empty() && value.back() == '\0') {
    value.pop_back();
  }
  return value;
}

void writeScalarU32Attribute(hid_t location, std::string_view key, std::uint32_t value) {
  Hdf5Handle scalar_space(H5Screate(H5S_SCALAR));
  Hdf5Handle attr(H5Acreate2(location, std::string(key).c_str(), H5T_STD_U32LE, scalar_space.get(), H5P_DEFAULT, H5P_DEFAULT));
  if (!attr.valid() || H5Awrite(attr.get(), H5T_NATIVE_UINT32, &value) < 0) {
    throw std::runtime_error("failed writing u32 attribute: " + std::string(key));
  }
}

[[nodiscard]] std::uint32_t readScalarU32Attribute(hid_t location, std::string_view key) {
  Hdf5Handle attr(H5Aopen(location, std::string(key).c_str(), H5P_DEFAULT));
  std::uint32_t value = 0;
  if (!attr.valid() || H5Aread(attr.get(), H5T_NATIVE_UINT32, &value) < 0) {
    throw std::runtime_error("failed reading u32 attribute: " + std::string(key));
  }
  return value;
}

void writeScalarU64Attribute(hid_t location, std::string_view key, std::uint64_t value) {
  Hdf5Handle scalar_space(H5Screate(H5S_SCALAR));
  Hdf5Handle attr(H5Acreate2(location, std::string(key).c_str(), H5T_STD_U64LE, scalar_space.get(), H5P_DEFAULT, H5P_DEFAULT));
  if (!attr.valid() || H5Awrite(attr.get(), H5T_NATIVE_UINT64, &value) < 0) {
    throw std::runtime_error("failed writing u64 attribute: " + std::string(key));
  }
}

[[nodiscard]] std::uint64_t readScalarU64Attribute(hid_t location, std::string_view key) {
  Hdf5Handle attr(H5Aopen(location, std::string(key).c_str(), H5P_DEFAULT));
  std::uint64_t value = 0;
  if (!attr.valid() || H5Aread(attr.get(), H5T_NATIVE_UINT64, &value) < 0) {
    throw std::runtime_error("failed reading u64 attribute: " + std::string(key));
  }
  return value;
}

void writeScalarF64Attribute(hid_t location, std::string_view key, double value) {
  Hdf5Handle scalar_space(H5Screate(H5S_SCALAR));
  Hdf5Handle attr(H5Acreate2(location, std::string(key).c_str(), H5T_IEEE_F64LE, scalar_space.get(), H5P_DEFAULT, H5P_DEFAULT));
  if (!attr.valid() || H5Awrite(attr.get(), H5T_NATIVE_DOUBLE, &value) < 0) {
    throw std::runtime_error("failed writing f64 attribute: " + std::string(key));
  }
}

[[nodiscard]] double readScalarF64Attribute(hid_t location, std::string_view key) {
  Hdf5Handle attr(H5Aopen(location, std::string(key).c_str(), H5P_DEFAULT));
  double value = 0.0;
  if (!attr.valid() || H5Aread(attr.get(), H5T_NATIVE_DOUBLE, &value) < 0) {
    throw std::runtime_error("failed reading f64 attribute: " + std::string(key));
  }
  return value;
}

[[nodiscard]] bool hdf5LinkExists(hid_t root, std::string_view path) {
  return H5Lexists(root, std::string(path).c_str(), H5P_DEFAULT) > 0;
}

[[nodiscard]] bool hdf5AttributeExists(hid_t location, std::string_view key) {
  return H5Aexists(location, std::string(key).c_str()) > 0;
}

void requireHdf5Link(hid_t root, std::string_view path) {
  if (!hdf5LinkExists(root, path)) {
    throw std::runtime_error("restart schema validation missing required path: " + std::string(path));
  }
}

Hdf5Handle openRequiredGroup(hid_t root, std::string_view path) {
  requireHdf5Link(root, path);
  Hdf5Handle group(H5Gopen2(root, std::string(path).c_str(), H5P_DEFAULT));
  if (!group.valid()) {
    throw std::runtime_error("restart schema validation expected group at path: " + std::string(path));
  }
  return group;
}

void requireHdf5Attribute(hid_t location, std::string_view location_path, std::string_view key) {
  if (!hdf5AttributeExists(location, key)) {
    throw std::runtime_error(
        "restart schema validation missing required attribute: " + std::string(location_path) + "/@" +
        std::string(key));
  }
}

void requireHdf5Dataset1d(hid_t root, std::string_view path) {
  requireHdf5Link(root, path);
  Hdf5Handle dataset(H5Dopen2(root, std::string(path).c_str(), H5P_DEFAULT));
  if (!dataset.valid()) {
    throw std::runtime_error("restart schema validation expected dataset at path: " + std::string(path));
  }
  Hdf5Handle space(H5Dget_space(dataset.get()));
  if (!space.valid() || H5Sget_simple_extent_ndims(space.get()) != 1) {
    throw std::runtime_error("restart schema validation expected 1D dataset at path: " + std::string(path));
  }
}

void validateFileKindAttribute(
    hid_t file,
    std::string_view expected_kind,
    std::string_view reader_label) {
  const auto& shared_names = sharedIoContractNames();
  if (!hdf5AttributeExists(file, shared_names.file_kind_attribute)) {
    std::string hint;
    if (hdf5LinkExists(file, "/Header")) {
      hint = "; file appears to be a science snapshot, not a restart checkpoint";
    }
    throw std::runtime_error(
        std::string(reader_label) + " missing required file-kind attribute /@" +
        std::string(shared_names.file_kind_attribute) + hint);
  }
  const std::string file_kind = readScalarStringAttribute(file, shared_names.file_kind_attribute);
  if (file_kind != expected_kind) {
    throw std::runtime_error(
        std::string(reader_label) + " rejected HDF5 file-kind '" + file_kind + "'; expected '" +
        std::string(expected_kind) + "'");
  }
}

void validateRestartCheckpointSchema(hid_t file, std::uint32_t schema_version) {
  const auto& shared_names = sharedIoContractNames();
  validateFileKindAttribute(file, shared_names.restart_checkpoint_file_kind, "restart reader");

  requireHdf5Attribute(file, "/", "restart_schema_name");
  requireHdf5Attribute(file, "/", "restart_schema_version");
  requireHdf5Attribute(file, "/", "normalized_config_hash_hex");
  requireHdf5Attribute(file, "/", "payload_integrity_hash_hex");
  requireHdf5Attribute(file, "/", "payload_integrity_hash");
  requireHdf5Dataset1d(file, "/normalized_config_text");
  requireHdf5Dataset1d(file, "/provenance_record");

  requireHdf5Dataset1d(file, "/state/state_metadata");
  requireHdf5Dataset1d(file, "/state/particles/position_x_comoving");
  requireHdf5Dataset1d(file, "/state/particles/position_y_comoving");
  requireHdf5Dataset1d(file, "/state/particles/position_z_comoving");
  requireHdf5Dataset1d(file, "/state/particles/velocity_x_peculiar");
  requireHdf5Dataset1d(file, "/state/particles/velocity_y_peculiar");
  requireHdf5Dataset1d(file, "/state/particles/velocity_z_peculiar");
  requireHdf5Dataset1d(file, "/state/particles/mass_code");
  requireHdf5Dataset1d(file, "/state/particles/time_bin");
  requireHdf5Dataset1d(file, "/state/particle_sidecar/particle_id");
  requireHdf5Dataset1d(file, "/state/particle_sidecar/sfc_key");
  requireHdf5Dataset1d(file, "/state/particle_sidecar/species_tag");
  requireHdf5Dataset1d(file, "/state/particle_sidecar/particle_flags");
  requireHdf5Dataset1d(file, "/state/particle_sidecar/owning_rank");
  requireHdf5Dataset1d(file, "/state/particle_sidecar/last_drift_time_code");
  requireHdf5Dataset1d(file, "/state/particle_sidecar/last_drift_scale_factor");
  requireHdf5Dataset1d(file, "/state/cells/center_x_comoving");
  requireHdf5Dataset1d(file, "/state/cells/center_y_comoving");
  requireHdf5Dataset1d(file, "/state/cells/center_z_comoving");
  requireHdf5Dataset1d(file, "/state/cells/mass_code");
  requireHdf5Dataset1d(file, "/state/cells/time_bin");
  requireHdf5Dataset1d(file, "/state/cells/patch_index");
  requireHdf5Dataset1d(file, "/state/gas_cells/gas_cell_id");
  requireHdf5Dataset1d(file, "/state/gas_cells/parent_particle_id");
  requireHdf5Dataset1d(file, "/state/gas_cells/velocity_x_peculiar");
  requireHdf5Dataset1d(file, "/state/gas_cells/velocity_y_peculiar");
  requireHdf5Dataset1d(file, "/state/gas_cells/velocity_z_peculiar");
  requireHdf5Dataset1d(file, "/state/gas_cells/density_code");
  requireHdf5Dataset1d(file, "/state/gas_cells/pressure_code");
  requireHdf5Dataset1d(file, "/state/gas_cells/internal_energy_code");
  requireHdf5Dataset1d(file, "/state/gas_cells/temperature_code");
  requireHdf5Dataset1d(file, "/state/gas_cells/sound_speed_code");
  if (schema_version >= k_restart_schema_v18) {
    Hdf5Handle identity_group = openRequiredGroup(file, "/state/gas_cell_identity");
    requireHdf5Attribute(identity_group.get(), "/state/gas_cell_identity", "local_row_reconstruction_policy");
    requireHdf5Attribute(identity_group.get(), "/state/gas_cell_identity", "identity_generation_at_write");
    requireHdf5Dataset1d(file, "/state/gas_cell_identity/gas_cell_id");
    requireHdf5Dataset1d(file, "/state/gas_cell_identity/has_parent_particle");
    requireHdf5Dataset1d(file, "/state/gas_cell_identity/parent_particle_id");
    requireHdf5Dataset1d(file, "/state/gas_cell_identity/owning_patch_id");
    requireHdf5Dataset1d(file, "/state/gas_cell_identity/local_cell_row");
  }
  requireHdf5Dataset1d(file, "/state/patches/patch_id");
  requireHdf5Dataset1d(file, "/state/patches/level");
  requireHdf5Dataset1d(file, "/state/patches/first_cell");
  requireHdf5Dataset1d(file, "/state/patches/cell_count");
  if (schema_version >= k_restart_schema_v18) {
    requireHdf5Dataset1d(file, "/state/patches/parent_patch_id");
    requireHdf5Dataset1d(file, "/state/patches/morton_key");
    requireHdf5Dataset1d(file, "/state/patches/origin_x_comoving");
    requireHdf5Dataset1d(file, "/state/patches/origin_y_comoving");
    requireHdf5Dataset1d(file, "/state/patches/origin_z_comoving");
    requireHdf5Dataset1d(file, "/state/patches/extent_x_comoving");
    requireHdf5Dataset1d(file, "/state/patches/extent_y_comoving");
    requireHdf5Dataset1d(file, "/state/patches/extent_z_comoving");
    requireHdf5Dataset1d(file, "/state/patches/cell_dim_x");
    requireHdf5Dataset1d(file, "/state/patches/cell_dim_y");
    requireHdf5Dataset1d(file, "/state/patches/cell_dim_z");
  }
  if (schema_version >= k_restart_schema_v18) {
    Hdf5Handle pending_group = openRequiredGroup(file, "/state/amr_pending_flux_registers");
    requireHdf5Attribute(pending_group.get(), "/state/amr_pending_flux_registers", "schema_version");
    for (std::string_view dataset : {"register_key", "coarse_patch_id", "coarse_gas_cell_id",
                                     "coarse_cell_index", "level", "axis", "orientation",
                                     "expected_area_comov", "coarse_area_accumulated_comov",
                                     "fine_area_accumulated_comov", "interval_start_code",
                                     "interval_end_code", "coarse_dt_code", "expected_fine_substeps",
                                     "completed_fine_substeps", "fine_substep_coverage_mask",
                                     "coarse_face_count", "fine_face_count",
                                     "gas_cell_identity_generation", "patch_geometry_generation",
                                     "coarse_mass_flux_integral_code",
                                     "coarse_momentum_x_flux_integral_code",
                                     "coarse_momentum_y_flux_integral_code",
                                     "coarse_momentum_z_flux_integral_code",
                                     "coarse_total_energy_flux_integral_code",
                                     "fine_mass_flux_integral_code", "fine_momentum_x_flux_integral_code",
                                     "fine_momentum_y_flux_integral_code",
                                     "fine_momentum_z_flux_integral_code",
                                     "fine_total_energy_flux_integral_code"}) {
      requireHdf5Dataset1d(file, std::string("/state/amr_pending_flux_registers/") + std::string(dataset));
    }
  }
  if (schema_version >= k_restart_schema_v18) {
    Hdf5Handle temporal_group = openRequiredGroup(file, "/state/amr_temporal_boundary_history");
    requireHdf5Attribute(temporal_group.get(), "/state/amr_temporal_boundary_history", "schema_version");
    for (std::string_view dataset : {"patch_id", "patch_level", "patch_geometry_fingerprint",
                                     "gas_cell_identity_generation", "interval_start_code",
                                     "interval_end_code", "end_state_valid", "cell_offset", "cell_count",
                                     "gas_cell_id", "patch_local_cell", "start_mass_density_comoving",
                                     "start_momentum_density_x_comoving", "start_momentum_density_y_comoving",
                                     "start_momentum_density_z_comoving", "start_total_energy_density_comoving",
                                     "end_mass_density_comoving", "end_momentum_density_x_comoving",
                                     "end_momentum_density_y_comoving", "end_momentum_density_z_comoving",
                                     "end_total_energy_density_comoving"}) {
      requireHdf5Dataset1d(file, std::string("/state/amr_temporal_boundary_history/") + std::string(dataset));
    }
  }
  requireHdf5Dataset1d(file, "/state/species_count_by_species");
  requireHdf5Link(file, "/state/module_sidecars");
  requireHdf5Dataset1d(file, "/state/module_sidecar_names");

  Hdf5Handle integrator_group = openRequiredGroup(file, "/integrator");
  for (std::string_view attr : {"current_time_code", "current_scale_factor", "current_redshift",
                                "current_hubble_rate_code", "time_si_per_code", "dt_time_code",
                                "last_drift_factor_code", "last_first_kick_factor_code",
                                "last_second_kick_factor_code", "last_first_hubble_drag_factor",
                                "last_second_hubble_drag_factor", "step_index", "scheme",
                                "current_boundary_kind", "last_completed_boundary_kind", "inside_kdk_step",
                                "last_completed_restart_safe", "time_bins_hierarchical", "time_bins_active_bin",
                                "time_bins_max_bin", "pm_long_range_field_valid", "pm_cadence_steps",
                                "pm_gravity_kick_opportunity", "pm_last_refresh_opportunity",
                                "pm_field_version", "pm_last_refresh_step_index",
                                "pm_last_refresh_scale_factor", "pm_refresh_commit_pending",
                                "pm_pending_refresh_opportunity", "pm_pending_refresh_field_version"}) {
    requireHdf5Attribute(integrator_group.get(), "/integrator", attr);
  }
  if (schema_version >= k_restart_schema_v20) {
    requireHdf5Attribute(integrator_group.get(), "/integrator", "pm_refresh_enabled");
  }

  if (schema_version >= k_restart_schema_v20) {
    Hdf5Handle force_cache_group = openRequiredGroup(file, "/gravity_force_cache");
    requireHdf5Attribute(force_cache_group.get(), "/gravity_force_cache", "valid");
    for (std::string_view dataset : {"particle_id", "gas_cell_id",
                                     "particle_accel_x_comoving", "particle_accel_y_comoving",
                                     "particle_accel_z_comoving", "cell_accel_x_comoving",
                                     "cell_accel_y_comoving", "cell_accel_z_comoving"}) {
      requireHdf5Dataset1d(file, std::string("/gravity_force_cache/") + std::string(dataset));
    }
  }

  Hdf5Handle scheduler_group = openRequiredGroup(file, "/scheduler");
  requireHdf5Attribute(scheduler_group.get(), "/scheduler", "current_tick");
  requireHdf5Attribute(scheduler_group.get(), "/scheduler", "max_bin");
  requireHdf5Dataset1d(file, "/scheduler/bin_index");
  requireHdf5Dataset1d(file, "/scheduler/next_activation_tick");
  requireHdf5Dataset1d(file, "/scheduler/active_flag");
  requireHdf5Dataset1d(file, "/scheduler/pending_bin_index");
  if (schema_version >= k_restart_schema_v19) {
    Hdf5Handle gas_scheduler_group = openRequiredGroup(file, "/gas_cell_scheduler");
    requireHdf5Attribute(gas_scheduler_group.get(), "/gas_cell_scheduler", "identity_key");
    requireHdf5Attribute(gas_scheduler_group.get(), "/gas_cell_scheduler", "current_tick");
    requireHdf5Attribute(gas_scheduler_group.get(), "/gas_cell_scheduler", "max_bin");
    requireHdf5Dataset1d(file, "/gas_cell_scheduler/gas_cell_id");
    requireHdf5Dataset1d(file, "/gas_cell_scheduler/bin_index");
    requireHdf5Dataset1d(file, "/gas_cell_scheduler/next_activation_tick");
    requireHdf5Dataset1d(file, "/gas_cell_scheduler/active_flag");
    requireHdf5Dataset1d(file, "/gas_cell_scheduler/pending_bin_index");
  }

  Hdf5Handle output_group = openRequiredGroup(file, "/output_cadence");
  for (std::string_view attr : {"output_enabled", "write_restarts", "snapshot_due", "checkpoint_due",
                                "last_completed_step_index", "snapshot_interval_steps",
                                "next_snapshot_step_index", "snapshot_stem", "restart_stem"}) {
    requireHdf5Attribute(output_group.get(), "/output_cadence", attr);
  }

  Hdf5Handle stochastic_group = openRequiredGroup(file, "/stochastic_state");
  requireHdf5Attribute(stochastic_group.get(), "/stochastic_state", "module_count");
  requireHdf5Dataset1d(file, "/stochastic_state/module_names");

  Hdf5Handle diagnostics_group = openRequiredGroup(file, "/restart_diagnostics");
  for (std::string_view attr : {"restart_schema_name", "restart_schema_version", "current_boundary_kind",
                                "last_completed_boundary_kind", "restart_safe", "step_index",
                                "scheduler_current_tick", "scheduler_max_bin", "scheduler_element_count",
                                "scheduler_active_count", "scheduler_pending_transition_count",
                                "pm_cadence_steps", "pm_gravity_kick_opportunity", "pm_field_version",
                                "pm_last_refresh_opportunity", "pm_last_refresh_step_index",
                                "pm_refresh_commit_pending", "pm_long_range_field_valid", "output_enabled",
                                "output_snapshot_due", "output_checkpoint_due",
                                "output_last_completed_step_index", "output_next_snapshot_step_index",
                                "stochastic_module_count"}) {
    requireHdf5Attribute(diagnostics_group.get(), "/restart_diagnostics", attr);
  }

  if (schema_version >= k_restart_schema_v19) {
    for (std::string_view attr : {"gas_cell_scheduler_current_tick", "gas_cell_scheduler_max_bin",
                                  "gas_cell_scheduler_element_count", "gas_cell_scheduler_active_count",
                                  "gas_cell_scheduler_pending_transition_count"}) {
      requireHdf5Attribute(diagnostics_group.get(), "/restart_diagnostics", attr);
    }
  }

  requireHdf5Dataset1d(file, "/distributed_gravity/state");
}

template <typename T>
void writeDataset1d(
    hid_t group,
    std::string_view name,
    hid_t file_type,
    hid_t memory_type,
    std::span<const T> values) {
  hsize_t dims[1] = {static_cast<hsize_t>(values.size())};
  Hdf5Handle space(H5Screate_simple(1, dims, nullptr));
  Hdf5Handle dataset(H5Dcreate2(group, std::string(name).c_str(), file_type, space.get(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
  if (!dataset.valid()) {
    throw std::runtime_error("failed creating dataset: " + std::string(name));
  }
  if (!values.empty() && H5Dwrite(dataset.get(), memory_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, values.data()) < 0) {
    throw std::runtime_error("failed writing dataset: " + std::string(name));
  }
}

template <typename T, typename Allocator>
void writeDataset1d(
    hid_t group,
    std::string_view name,
    hid_t file_type,
    hid_t memory_type,
    const std::vector<T, Allocator>& values) {
  writeDataset1d<T>(group, name, file_type, memory_type, std::span<const T>(values.data(), values.size()));
}

template <typename T>
[[nodiscard]] std::vector<T> readDataset1d(hid_t group, std::string_view name, hid_t memory_type) {
  const std::string dataset_name(name);
  const htri_t exists = H5Lexists(group, dataset_name.c_str(), H5P_DEFAULT);
  if (exists < 0) {
    throw std::runtime_error("failed checking dataset presence: " + dataset_name);
  }
  if (exists == 0) {
    throw std::runtime_error("missing dataset: " + dataset_name);
  }
  Hdf5Handle dataset(H5Dopen2(group, dataset_name.c_str(), H5P_DEFAULT));
  if (!dataset.valid()) {
    throw std::runtime_error("missing dataset: " + dataset_name);
  }
  Hdf5Handle space(H5Dget_space(dataset.get()));
  hsize_t dims[1] = {0};
  if (H5Sget_simple_extent_ndims(space.get()) != 1 || H5Sget_simple_extent_dims(space.get(), dims, nullptr) != 1) {
    throw std::runtime_error("unexpected rank for dataset: " + dataset_name);
  }
  std::vector<T> values(static_cast<std::size_t>(dims[0]));
  if (!values.empty() && H5Dread(dataset.get(), memory_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, values.data()) < 0) {
    throw std::runtime_error("failed reading dataset: " + dataset_name);
  }
  return values;
}

template <typename T>
[[nodiscard]] core::AlignedVector<T> toAlignedVector(const std::vector<T>& values) {
  core::AlignedVector<T> aligned(values.size());
  std::copy(values.begin(), values.end(), aligned.begin());
  return aligned;
}

template <typename T>
[[nodiscard]] core::AlignedVector<T> readDataset1dAligned(hid_t group, std::string_view name, hid_t memory_type) {
  return toAlignedVector(readDataset1d<T>(group, name, memory_type));
}

void writeStringDataset(hid_t group, std::string_view name, const std::string& value) {
  const std::vector<std::uint8_t> bytes(value.begin(), value.end());
  writeDataset1d(group, name, H5T_STD_U8LE, H5T_NATIVE_UINT8, bytes);
}

[[nodiscard]] std::string readStringDataset(hid_t group, std::string_view name) {
  const auto bytes = readDataset1d<std::uint8_t>(group, name, H5T_NATIVE_UINT8);
  return std::string(bytes.begin(), bytes.end());
}

[[nodiscard]] hid_t openOrCreateGroup(hid_t parent, const std::string& name) {
  if (H5Lexists(parent, name.c_str(), H5P_DEFAULT) > 0) {
    return H5Gopen2(parent, name.c_str(), H5P_DEFAULT);
  }
  return H5Gcreate2(parent, name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
}

void maybeFsync(const std::filesystem::path& file_path, bool enabled) {
  if (!enabled) {
    return;
  }
#if defined(_WIN32)
  (void)file_path;
  throw std::runtime_error("restart fsync finalize is not implemented on this platform");
#else
  const int fd = ::open(file_path.c_str(), O_RDONLY);
  if (fd < 0) {
    throw std::runtime_error(
        "failed to open restart temporary file for fsync: " + file_path.string() + ": " + std::strerror(errno));
  }
  if (::fsync(fd) != 0) {
    const std::string message = std::strerror(errno);
    ::close(fd);
    throw std::runtime_error("failed to fsync restart temporary file: " + file_path.string() + ": " + message);
  }
  if (::close(fd) != 0) {
    throw std::runtime_error(
        "failed to close restart temporary file after fsync: " + file_path.string() + ": " + std::strerror(errno));
  }
#endif
}

void writeStarSidecarGroup(hid_t state_group, const core::StarParticleSidecar& stars) {
  Hdf5Handle star_group(openOrCreateGroup(state_group, "star_particles"));
  writeDataset1d(star_group.get(), "particle_index", H5T_STD_U32LE, H5T_NATIVE_UINT32, stars.particle_index);
  writeDataset1d(
      star_group.get(), "formation_scale_factor", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, stars.formation_scale_factor);
  writeDataset1d(star_group.get(), "birth_mass_code", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, stars.birth_mass_code);
  writeDataset1d(
      star_group.get(),
      "metallicity_mass_fraction",
      H5T_IEEE_F64LE,
      H5T_NATIVE_DOUBLE,
      stars.metallicity_mass_fraction);
  writeDataset1d(
      star_group.get(),
      "stellar_age_years_last",
      H5T_IEEE_F64LE,
      H5T_NATIVE_DOUBLE,
      stars.stellar_age_years_last);
  writeDataset1d(
      star_group.get(),
      "stellar_returned_mass_cumulative_code",
      H5T_IEEE_F64LE,
      H5T_NATIVE_DOUBLE,
      stars.stellar_returned_mass_cumulative_code);
  writeDataset1d(
      star_group.get(),
      "stellar_returned_metals_cumulative_code",
      H5T_IEEE_F64LE,
      H5T_NATIVE_DOUBLE,
      stars.stellar_returned_metals_cumulative_code);
  writeDataset1d(
      star_group.get(),
      "stellar_feedback_energy_cumulative_erg",
      H5T_IEEE_F64LE,
      H5T_NATIVE_DOUBLE,
      stars.stellar_feedback_energy_cumulative_erg);

  for (std::size_t channel = 0; channel < stars.stellar_returned_mass_channel_cumulative_code.size(); ++channel) {
    const std::string suffix = std::to_string(channel);
    writeDataset1d(
        star_group.get(),
        "stellar_returned_mass_channel_cumulative_code_" + suffix,
        H5T_IEEE_F64LE,
        H5T_NATIVE_DOUBLE,
        stars.stellar_returned_mass_channel_cumulative_code[channel]);
    writeDataset1d(
        star_group.get(),
        "stellar_returned_metals_channel_cumulative_code_" + suffix,
        H5T_IEEE_F64LE,
        H5T_NATIVE_DOUBLE,
        stars.stellar_returned_metals_channel_cumulative_code[channel]);
    writeDataset1d(
        star_group.get(),
        "stellar_feedback_energy_channel_cumulative_erg_" + suffix,
        H5T_IEEE_F64LE,
        H5T_NATIVE_DOUBLE,
        stars.stellar_feedback_energy_channel_cumulative_erg[channel]);
  }
}

void readStarSidecarGroup(hid_t state_group, core::StarParticleSidecar& stars) {
  Hdf5Handle star_group(H5Gopen2(state_group, "star_particles", H5P_DEFAULT));
  stars.particle_index = readDataset1dAligned<std::uint32_t>(star_group.get(), "particle_index", H5T_NATIVE_UINT32);
  stars.formation_scale_factor =
      readDataset1dAligned<double>(star_group.get(), "formation_scale_factor", H5T_NATIVE_DOUBLE);
  stars.birth_mass_code = readDataset1dAligned<double>(star_group.get(), "birth_mass_code", H5T_NATIVE_DOUBLE);
  stars.metallicity_mass_fraction =
      readDataset1dAligned<double>(star_group.get(), "metallicity_mass_fraction", H5T_NATIVE_DOUBLE);
  stars.stellar_age_years_last =
      readDataset1dAligned<double>(star_group.get(), "stellar_age_years_last", H5T_NATIVE_DOUBLE);
  stars.stellar_returned_mass_cumulative_code = readDataset1dAligned<double>(
      star_group.get(), "stellar_returned_mass_cumulative_code", H5T_NATIVE_DOUBLE);
  stars.stellar_returned_metals_cumulative_code = readDataset1dAligned<double>(
      star_group.get(), "stellar_returned_metals_cumulative_code", H5T_NATIVE_DOUBLE);
  stars.stellar_feedback_energy_cumulative_erg = readDataset1dAligned<double>(
      star_group.get(), "stellar_feedback_energy_cumulative_erg", H5T_NATIVE_DOUBLE);

  for (std::size_t channel = 0; channel < stars.stellar_returned_mass_channel_cumulative_code.size(); ++channel) {
    const std::string suffix = std::to_string(channel);
    stars.stellar_returned_mass_channel_cumulative_code[channel] = readDataset1dAligned<double>(
        star_group.get(),
        "stellar_returned_mass_channel_cumulative_code_" + suffix,
        H5T_NATIVE_DOUBLE);
    stars.stellar_returned_metals_channel_cumulative_code[channel] = readDataset1dAligned<double>(
        star_group.get(),
        "stellar_returned_metals_channel_cumulative_code_" + suffix,
        H5T_NATIVE_DOUBLE);
    stars.stellar_feedback_energy_channel_cumulative_erg[channel] = readDataset1dAligned<double>(
        star_group.get(),
        "stellar_feedback_energy_channel_cumulative_erg_" + suffix,
        H5T_NATIVE_DOUBLE);
  }
}

void writeGasCellIdentityGroup(hid_t state_group, const core::SimulationState& state) {
  Hdf5Handle identity_group(openOrCreateGroup(state_group, "gas_cell_identity"));
  writeScalarStringAttribute(identity_group.get(), "local_row_reconstruction_policy", std::string(k_gas_identity_row_policy));
  writeScalarU64Attribute(identity_group.get(), "identity_generation_at_write", state.gasCellIdentityGeneration());

  const auto records = gasIdentityRecordsSortedByLocalRow(state.gas_cell_identity);
  std::vector<std::uint64_t> gas_cell_id;
  std::vector<std::uint8_t> has_parent_particle;
  std::vector<std::uint64_t> parent_particle_id;
  std::vector<std::uint64_t> owning_patch_id;
  std::vector<std::uint32_t> local_cell_row;
  gas_cell_id.reserve(records.size());
  has_parent_particle.reserve(records.size());
  parent_particle_id.reserve(records.size());
  owning_patch_id.reserve(records.size());
  local_cell_row.reserve(records.size());
  for (const core::GasCellIdentityRecord& record : records) {
    gas_cell_id.push_back(record.gas_cell_id);
    has_parent_particle.push_back(record.parent_particle_id.has_value() ? 1U : 0U);
    parent_particle_id.push_back(record.parent_particle_id.value_or(0U));
    owning_patch_id.push_back(record.owning_patch_id);
    local_cell_row.push_back(record.local_cell_row);
  }

  writeDataset1d(identity_group.get(), "gas_cell_id", H5T_STD_U64LE, H5T_NATIVE_UINT64, gas_cell_id);
  writeDataset1d(identity_group.get(), "has_parent_particle", H5T_STD_U8LE, H5T_NATIVE_UINT8, has_parent_particle);
  writeDataset1d(identity_group.get(), "parent_particle_id", H5T_STD_U64LE, H5T_NATIVE_UINT64, parent_particle_id);
  writeDataset1d(identity_group.get(), "owning_patch_id", H5T_STD_U64LE, H5T_NATIVE_UINT64, owning_patch_id);
  writeDataset1d(identity_group.get(), "local_cell_row", H5T_STD_U32LE, H5T_NATIVE_UINT32, local_cell_row);
}

void materializeLegacyGasCellIdentityMapFromMirrors(core::SimulationState& state) {
  std::vector<core::GasCellIdentityRecord> records;
  records.reserve(state.cells.size());
  for (std::size_t cell_index = 0; cell_index < state.cells.size(); ++cell_index) {
    const std::uint64_t gas_cell_id = state.gas_cells.gas_cell_id[cell_index];
    const std::uint64_t parent_particle_id = state.gas_cells.parent_particle_id[cell_index];
    if (gas_cell_id == 0U || parent_particle_id == 0U || gas_cell_id != parent_particle_id) {
      throw std::runtime_error(
          "legacy v14 restart gas-cell identity requires nonzero gas_cell_id == parent_particle_id");
    }
    std::uint64_t owning_patch_id = 0U;
    if (!state.patches.patch_id.empty()) {
      const std::uint32_t patch_index = state.cells.patch_index[cell_index];
      if (patch_index >= state.patches.size()) {
        throw std::runtime_error("legacy v14 restart gas-cell identity has invalid patch_index");
      }
      owning_patch_id = state.patches.patch_id[patch_index];
    }
    records.push_back(core::GasCellIdentityRecord{
        .gas_cell_id = gas_cell_id,
        .parent_particle_id = parent_particle_id,
        .owning_patch_id = owning_patch_id,
        .local_cell_row = static_cast<std::uint32_t>(cell_index),
    });
  }
  state.replaceGasCellIdentityRecords(std::move(records));
}

void readGasCellIdentityGroup(hid_t state_group, core::SimulationState& state, std::uint32_t schema_version) {
  if (schema_version == k_restart_schema_v14 && H5Lexists(state_group, "gas_cell_identity", H5P_DEFAULT) <= 0) {
    materializeLegacyGasCellIdentityMapFromMirrors(state);
    return;
  }

  Hdf5Handle identity_group(H5Gopen2(state_group, "gas_cell_identity", H5P_DEFAULT));
  if (!identity_group.valid()) {
    throw std::runtime_error("restart is missing /state/gas_cell_identity");
  }
  const std::string row_policy =
      readScalarStringAttribute(identity_group.get(), "local_row_reconstruction_policy");
  if (row_policy != k_gas_identity_row_policy) {
    throw std::runtime_error("unsupported /state/gas_cell_identity local_row_reconstruction_policy");
  }
  const std::uint64_t identity_generation_at_write =
      readScalarU64Attribute(identity_group.get(), "identity_generation_at_write");

  const auto gas_cell_id =
      readDataset1d<std::uint64_t>(identity_group.get(), "gas_cell_id", H5T_NATIVE_UINT64);
  const auto has_parent_particle =
      readDataset1d<std::uint8_t>(identity_group.get(), "has_parent_particle", H5T_NATIVE_UINT8);
  const auto parent_particle_id =
      readDataset1d<std::uint64_t>(identity_group.get(), "parent_particle_id", H5T_NATIVE_UINT64);
  const auto owning_patch_id =
      readDataset1d<std::uint64_t>(identity_group.get(), "owning_patch_id", H5T_NATIVE_UINT64);
  const auto local_cell_row =
      readDataset1d<std::uint32_t>(identity_group.get(), "local_cell_row", H5T_NATIVE_UINT32);
  if (gas_cell_id.size() != has_parent_particle.size() ||
      gas_cell_id.size() != parent_particle_id.size() ||
      gas_cell_id.size() != owning_patch_id.size() ||
      gas_cell_id.size() != local_cell_row.size()) {
    throw std::runtime_error("/state/gas_cell_identity datasets must have matching lengths");
  }

  std::vector<core::GasCellIdentityRecord> records;
  records.reserve(gas_cell_id.size());
  for (std::size_t i = 0; i < gas_cell_id.size(); ++i) {
    if (has_parent_particle[i] > 1U) {
      throw std::runtime_error("/state/gas_cell_identity/has_parent_particle must contain only 0 or 1");
    }
    if (has_parent_particle[i] == 0U && parent_particle_id[i] != 0U) {
      throw std::runtime_error(
          "/state/gas_cell_identity/parent_particle_id must be 0 when has_parent_particle is false");
    }
    if (has_parent_particle[i] != 0U && parent_particle_id[i] == 0U) {
      throw std::runtime_error(
          "/state/gas_cell_identity/parent_particle_id must be nonzero when has_parent_particle is true");
    }
    records.push_back(core::GasCellIdentityRecord{
        .gas_cell_id = gas_cell_id[i],
        .parent_particle_id =
            has_parent_particle[i] != 0U ? std::optional<std::uint64_t>(parent_particle_id[i]) : std::nullopt,
        .owning_patch_id = owning_patch_id[i],
        .local_cell_row = local_cell_row[i],
    });
  }
  // Validate persisted compatibility mirrors before normalizing them from the
  // canonical map.  Otherwise a corrupt /state/cells/patch_index lane could be
  // silently overwritten during restore and evade restart integrity checks.
  if (records.size() != state.cells.size()) {
    throw std::runtime_error("/state/gas_cell_identity record count must match /state/cells extent");
  }
  for (const core::GasCellIdentityRecord& record : records) {
    if (record.local_cell_row >= state.cells.size()) {
      throw std::runtime_error("/state/gas_cell_identity/local_cell_row exceeds /state/cells extent");
    }
    const std::uint32_t row = record.local_cell_row;
    if (state.gas_cells.gas_cell_id[row] != record.gas_cell_id ||
        state.gas_cells.parent_particle_id[row] != record.parent_particle_id.value_or(0U)) {
      throw std::runtime_error("/state/gas_cell_identity records disagree with gas-cell compatibility mirrors");
    }
    if (record.owning_patch_id != 0U) {
      if (state.cells.patch_index[row] >= state.patches.size() ||
          state.patches.patch_id[state.cells.patch_index[row]] != record.owning_patch_id) {
        throw std::runtime_error("/state/cells/patch_index disagrees with /state/gas_cell_identity owning_patch_id");
      }
    }
  }
  state.restoreGasCellIdentityRecords(std::move(records), identity_generation_at_write);
}

void writePendingFluxRegisterGroup(hid_t state_group, const core::PendingFluxRegisterStore& store) {
  Hdf5Handle group(openOrCreateGroup(state_group, "amr_pending_flux_registers"));
  writeScalarU32Attribute(group.get(), "schema_version", 1U);
  const auto records = store.records();
  std::vector<std::uint64_t> register_key;
  std::vector<std::uint64_t> coarse_patch_id;
  std::vector<std::uint64_t> coarse_gas_cell_id;
  std::vector<std::uint64_t> coarse_cell_index;
  std::vector<std::uint8_t> level;
  std::vector<std::uint8_t> axis;
  std::vector<std::uint8_t> orientation;
  std::vector<double> expected_area_comov;
  std::vector<double> coarse_area_accumulated_comov;
  std::vector<double> fine_area_accumulated_comov;
  std::vector<double> interval_start_code;
  std::vector<double> interval_end_code;
  std::vector<double> coarse_dt_code;
  std::vector<std::uint32_t> expected_fine_substeps;
  std::vector<std::uint32_t> completed_fine_substeps;
  std::vector<std::uint64_t> fine_substep_coverage_mask;
  std::vector<std::uint32_t> coarse_face_count;
  std::vector<std::uint32_t> fine_face_count;
  std::vector<std::uint64_t> gas_cell_identity_generation;
  std::vector<std::uint64_t> patch_geometry_generation;
  std::vector<double> coarse_mass_flux_integral_code;
  std::vector<double> coarse_momentum_x_flux_integral_code;
  std::vector<double> coarse_momentum_y_flux_integral_code;
  std::vector<double> coarse_momentum_z_flux_integral_code;
  std::vector<double> coarse_total_energy_flux_integral_code;
  std::vector<double> fine_mass_flux_integral_code;
  std::vector<double> fine_momentum_x_flux_integral_code;
  std::vector<double> fine_momentum_y_flux_integral_code;
  std::vector<double> fine_momentum_z_flux_integral_code;
  std::vector<double> fine_total_energy_flux_integral_code;
  register_key.reserve(records.size());
  for (const core::PendingFluxRegisterRecord& record : records) {
    register_key.push_back(record.register_key);
    coarse_patch_id.push_back(record.coarse_patch_id);
    coarse_gas_cell_id.push_back(record.coarse_gas_cell_id);
    coarse_cell_index.push_back(static_cast<std::uint64_t>(record.coarse_cell_index));
    level.push_back(record.level);
    axis.push_back(record.axis);
    orientation.push_back(record.orientation);
    expected_area_comov.push_back(record.expected_area_comov);
    coarse_area_accumulated_comov.push_back(record.coarse_area_accumulated_comov);
    fine_area_accumulated_comov.push_back(record.fine_area_accumulated_comov);
    interval_start_code.push_back(record.interval_start_code);
    interval_end_code.push_back(record.interval_end_code);
    coarse_dt_code.push_back(record.coarse_dt_code);
    expected_fine_substeps.push_back(record.expected_fine_substeps);
    completed_fine_substeps.push_back(record.completed_fine_substeps);
    fine_substep_coverage_mask.push_back(record.fine_substep_coverage_mask);
    coarse_face_count.push_back(record.coarse_face_count);
    fine_face_count.push_back(record.fine_face_count);
    gas_cell_identity_generation.push_back(record.gas_cell_identity_generation);
    patch_geometry_generation.push_back(record.patch_geometry_generation);
    coarse_mass_flux_integral_code.push_back(record.coarse_mass_flux_integral_code);
    coarse_momentum_x_flux_integral_code.push_back(record.coarse_momentum_x_flux_integral_code);
    coarse_momentum_y_flux_integral_code.push_back(record.coarse_momentum_y_flux_integral_code);
    coarse_momentum_z_flux_integral_code.push_back(record.coarse_momentum_z_flux_integral_code);
    coarse_total_energy_flux_integral_code.push_back(record.coarse_total_energy_flux_integral_code);
    fine_mass_flux_integral_code.push_back(record.fine_mass_flux_integral_code);
    fine_momentum_x_flux_integral_code.push_back(record.fine_momentum_x_flux_integral_code);
    fine_momentum_y_flux_integral_code.push_back(record.fine_momentum_y_flux_integral_code);
    fine_momentum_z_flux_integral_code.push_back(record.fine_momentum_z_flux_integral_code);
    fine_total_energy_flux_integral_code.push_back(record.fine_total_energy_flux_integral_code);
  }
  writeDataset1d(group.get(), "register_key", H5T_STD_U64LE, H5T_NATIVE_UINT64, register_key);
  writeDataset1d(group.get(), "coarse_patch_id", H5T_STD_U64LE, H5T_NATIVE_UINT64, coarse_patch_id);
  writeDataset1d(group.get(), "coarse_gas_cell_id", H5T_STD_U64LE, H5T_NATIVE_UINT64, coarse_gas_cell_id);
  writeDataset1d(group.get(), "coarse_cell_index", H5T_STD_U64LE, H5T_NATIVE_UINT64, coarse_cell_index);
  writeDataset1d(group.get(), "level", H5T_STD_U8LE, H5T_NATIVE_UINT8, level);
  writeDataset1d(group.get(), "axis", H5T_STD_U8LE, H5T_NATIVE_UINT8, axis);
  writeDataset1d(group.get(), "orientation", H5T_STD_U8LE, H5T_NATIVE_UINT8, orientation);
  writeDataset1d(group.get(), "expected_area_comov", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, expected_area_comov);
  writeDataset1d(group.get(), "coarse_area_accumulated_comov", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, coarse_area_accumulated_comov);
  writeDataset1d(group.get(), "fine_area_accumulated_comov", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, fine_area_accumulated_comov);
  writeDataset1d(group.get(), "interval_start_code", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, interval_start_code);
  writeDataset1d(group.get(), "interval_end_code", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, interval_end_code);
  writeDataset1d(group.get(), "coarse_dt_code", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, coarse_dt_code);
  writeDataset1d(group.get(), "expected_fine_substeps", H5T_STD_U32LE, H5T_NATIVE_UINT32, expected_fine_substeps);
  writeDataset1d(group.get(), "completed_fine_substeps", H5T_STD_U32LE, H5T_NATIVE_UINT32, completed_fine_substeps);
  writeDataset1d(group.get(), "fine_substep_coverage_mask", H5T_STD_U64LE, H5T_NATIVE_UINT64, fine_substep_coverage_mask);
  writeDataset1d(group.get(), "coarse_face_count", H5T_STD_U32LE, H5T_NATIVE_UINT32, coarse_face_count);
  writeDataset1d(group.get(), "fine_face_count", H5T_STD_U32LE, H5T_NATIVE_UINT32, fine_face_count);
  writeDataset1d(group.get(), "gas_cell_identity_generation", H5T_STD_U64LE, H5T_NATIVE_UINT64, gas_cell_identity_generation);
  writeDataset1d(group.get(), "patch_geometry_generation", H5T_STD_U64LE, H5T_NATIVE_UINT64, patch_geometry_generation);
  writeDataset1d(group.get(), "coarse_mass_flux_integral_code", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, coarse_mass_flux_integral_code);
  writeDataset1d(group.get(), "coarse_momentum_x_flux_integral_code", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, coarse_momentum_x_flux_integral_code);
  writeDataset1d(group.get(), "coarse_momentum_y_flux_integral_code", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, coarse_momentum_y_flux_integral_code);
  writeDataset1d(group.get(), "coarse_momentum_z_flux_integral_code", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, coarse_momentum_z_flux_integral_code);
  writeDataset1d(group.get(), "coarse_total_energy_flux_integral_code", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, coarse_total_energy_flux_integral_code);
  writeDataset1d(group.get(), "fine_mass_flux_integral_code", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, fine_mass_flux_integral_code);
  writeDataset1d(group.get(), "fine_momentum_x_flux_integral_code", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, fine_momentum_x_flux_integral_code);
  writeDataset1d(group.get(), "fine_momentum_y_flux_integral_code", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, fine_momentum_y_flux_integral_code);
  writeDataset1d(group.get(), "fine_momentum_z_flux_integral_code", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, fine_momentum_z_flux_integral_code);
  writeDataset1d(group.get(), "fine_total_energy_flux_integral_code", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, fine_total_energy_flux_integral_code);
}

void readPendingFluxRegisterGroup(hid_t state_group, core::SimulationState& state, std::uint32_t schema_version) {
  if (schema_version < k_restart_schema_v19 && H5Lexists(state_group, "amr_pending_flux_registers", H5P_DEFAULT) <= 0) {
    state.pending_flux_registers.clear();
    return;
  }
  Hdf5Handle group(H5Gopen2(state_group, "amr_pending_flux_registers", H5P_DEFAULT));
  if (!group.valid()) {
    state.pending_flux_registers.clear();
    return;
  }
  const auto register_key = readDataset1d<std::uint64_t>(group.get(), "register_key", H5T_NATIVE_UINT64);
  const auto coarse_patch_id = readDataset1d<std::uint64_t>(group.get(), "coarse_patch_id", H5T_NATIVE_UINT64);
  const auto coarse_gas_cell_id = readDataset1d<std::uint64_t>(group.get(), "coarse_gas_cell_id", H5T_NATIVE_UINT64);
  const auto coarse_cell_index = readDataset1d<std::uint64_t>(group.get(), "coarse_cell_index", H5T_NATIVE_UINT64);
  const auto level = readDataset1d<std::uint8_t>(group.get(), "level", H5T_NATIVE_UINT8);
  const auto axis = readDataset1d<std::uint8_t>(group.get(), "axis", H5T_NATIVE_UINT8);
  const auto orientation = readDataset1d<std::uint8_t>(group.get(), "orientation", H5T_NATIVE_UINT8);
  const auto expected_area_comov = readDataset1d<double>(group.get(), "expected_area_comov", H5T_NATIVE_DOUBLE);
  const auto coarse_area_accumulated_comov = readDataset1d<double>(group.get(), "coarse_area_accumulated_comov", H5T_NATIVE_DOUBLE);
  const auto fine_area_accumulated_comov = readDataset1d<double>(group.get(), "fine_area_accumulated_comov", H5T_NATIVE_DOUBLE);
  const auto interval_start_code = readDataset1d<double>(group.get(), "interval_start_code", H5T_NATIVE_DOUBLE);
  const auto interval_end_code = readDataset1d<double>(group.get(), "interval_end_code", H5T_NATIVE_DOUBLE);
  const auto coarse_dt_code = readDataset1d<double>(group.get(), "coarse_dt_code", H5T_NATIVE_DOUBLE);
  const auto expected_fine_substeps = readDataset1d<std::uint32_t>(group.get(), "expected_fine_substeps", H5T_NATIVE_UINT32);
  const auto completed_fine_substeps = readDataset1d<std::uint32_t>(group.get(), "completed_fine_substeps", H5T_NATIVE_UINT32);
  const auto fine_substep_coverage_mask = readDataset1d<std::uint64_t>(group.get(), "fine_substep_coverage_mask", H5T_NATIVE_UINT64);
  const auto coarse_face_count = readDataset1d<std::uint32_t>(group.get(), "coarse_face_count", H5T_NATIVE_UINT32);
  const auto fine_face_count = readDataset1d<std::uint32_t>(group.get(), "fine_face_count", H5T_NATIVE_UINT32);
  const auto gas_cell_identity_generation = readDataset1d<std::uint64_t>(group.get(), "gas_cell_identity_generation", H5T_NATIVE_UINT64);
  const auto patch_geometry_generation = readDataset1d<std::uint64_t>(group.get(), "patch_geometry_generation", H5T_NATIVE_UINT64);
  const auto coarse_mass_flux_integral_code = readDataset1d<double>(group.get(), "coarse_mass_flux_integral_code", H5T_NATIVE_DOUBLE);
  const auto coarse_momentum_x_flux_integral_code = readDataset1d<double>(group.get(), "coarse_momentum_x_flux_integral_code", H5T_NATIVE_DOUBLE);
  const auto coarse_momentum_y_flux_integral_code = readDataset1d<double>(group.get(), "coarse_momentum_y_flux_integral_code", H5T_NATIVE_DOUBLE);
  const auto coarse_momentum_z_flux_integral_code = readDataset1d<double>(group.get(), "coarse_momentum_z_flux_integral_code", H5T_NATIVE_DOUBLE);
  const auto coarse_total_energy_flux_integral_code = readDataset1d<double>(group.get(), "coarse_total_energy_flux_integral_code", H5T_NATIVE_DOUBLE);
  const auto fine_mass_flux_integral_code = readDataset1d<double>(group.get(), "fine_mass_flux_integral_code", H5T_NATIVE_DOUBLE);
  const auto fine_momentum_x_flux_integral_code = readDataset1d<double>(group.get(), "fine_momentum_x_flux_integral_code", H5T_NATIVE_DOUBLE);
  const auto fine_momentum_y_flux_integral_code = readDataset1d<double>(group.get(), "fine_momentum_y_flux_integral_code", H5T_NATIVE_DOUBLE);
  const auto fine_momentum_z_flux_integral_code = readDataset1d<double>(group.get(), "fine_momentum_z_flux_integral_code", H5T_NATIVE_DOUBLE);
  const auto fine_total_energy_flux_integral_code = readDataset1d<double>(group.get(), "fine_total_energy_flux_integral_code", H5T_NATIVE_DOUBLE);
  const std::size_t n = register_key.size();
  const bool sizes_match = coarse_patch_id.size() == n && coarse_gas_cell_id.size() == n && coarse_cell_index.size() == n &&
      level.size() == n && axis.size() == n && orientation.size() == n && expected_area_comov.size() == n &&
      coarse_area_accumulated_comov.size() == n && fine_area_accumulated_comov.size() == n && interval_start_code.size() == n &&
      interval_end_code.size() == n && coarse_dt_code.size() == n && expected_fine_substeps.size() == n &&
      completed_fine_substeps.size() == n && fine_substep_coverage_mask.size() == n && coarse_face_count.size() == n &&
      fine_face_count.size() == n && gas_cell_identity_generation.size() == n && patch_geometry_generation.size() == n &&
      coarse_mass_flux_integral_code.size() == n && coarse_momentum_x_flux_integral_code.size() == n &&
      coarse_momentum_y_flux_integral_code.size() == n && coarse_momentum_z_flux_integral_code.size() == n &&
      coarse_total_energy_flux_integral_code.size() == n && fine_mass_flux_integral_code.size() == n &&
      fine_momentum_x_flux_integral_code.size() == n && fine_momentum_y_flux_integral_code.size() == n &&
      fine_momentum_z_flux_integral_code.size() == n && fine_total_energy_flux_integral_code.size() == n;
  if (!sizes_match) {
    throw std::runtime_error("/state/amr_pending_flux_registers datasets must have matching lengths");
  }
  std::vector<core::PendingFluxRegisterRecord> records;
  records.reserve(n);
  for (std::size_t i = 0; i < n; ++i) {
    records.push_back(core::PendingFluxRegisterRecord{
        .register_key = register_key[i],
        .coarse_patch_id = coarse_patch_id[i],
        .coarse_gas_cell_id = coarse_gas_cell_id[i],
        .coarse_cell_index = static_cast<std::size_t>(coarse_cell_index[i]),
        .level = level[i],
        .axis = axis[i],
        .orientation = orientation[i],
        .expected_area_comov = expected_area_comov[i],
        .coarse_area_accumulated_comov = coarse_area_accumulated_comov[i],
        .fine_area_accumulated_comov = fine_area_accumulated_comov[i],
        .interval_start_code = interval_start_code[i],
        .interval_end_code = interval_end_code[i],
        .coarse_dt_code = coarse_dt_code[i],
        .expected_fine_substeps = expected_fine_substeps[i],
        .completed_fine_substeps = completed_fine_substeps[i],
        .fine_substep_coverage_mask = fine_substep_coverage_mask[i],
        .coarse_face_count = coarse_face_count[i],
        .fine_face_count = fine_face_count[i],
        .gas_cell_identity_generation = gas_cell_identity_generation[i],
        .patch_geometry_generation = patch_geometry_generation[i],
        .coarse_mass_flux_integral_code = coarse_mass_flux_integral_code[i],
        .coarse_momentum_x_flux_integral_code = coarse_momentum_x_flux_integral_code[i],
        .coarse_momentum_y_flux_integral_code = coarse_momentum_y_flux_integral_code[i],
        .coarse_momentum_z_flux_integral_code = coarse_momentum_z_flux_integral_code[i],
        .coarse_total_energy_flux_integral_code = coarse_total_energy_flux_integral_code[i],
        .fine_mass_flux_integral_code = fine_mass_flux_integral_code[i],
        .fine_momentum_x_flux_integral_code = fine_momentum_x_flux_integral_code[i],
        .fine_momentum_y_flux_integral_code = fine_momentum_y_flux_integral_code[i],
        .fine_momentum_z_flux_integral_code = fine_momentum_z_flux_integral_code[i],
        .fine_total_energy_flux_integral_code = fine_total_energy_flux_integral_code[i]});
  }
  state.pending_flux_registers.assign(std::move(records));
}


void writeAmrTemporalBoundaryHistoryGroup(
    hid_t state_group,
    const core::AmrTemporalBoundaryHistoryStore& store) {
  Hdf5Handle group(openOrCreateGroup(state_group, "amr_temporal_boundary_history"));
  writeScalarU32Attribute(group.get(), "schema_version", 1U);
  std::vector<std::uint64_t> patch_id;
  std::vector<std::uint8_t> patch_level;
  std::vector<std::uint64_t> patch_geometry_fingerprint;
  std::vector<std::uint64_t> gas_cell_identity_generation;
  std::vector<double> interval_start_code;
  std::vector<double> interval_end_code;
  std::vector<std::uint8_t> end_state_valid;
  std::vector<std::uint64_t> cell_offset;
  std::vector<std::uint64_t> cell_count;
  std::vector<std::uint64_t> gas_cell_id;
  std::vector<std::uint64_t> patch_local_cell;
  std::vector<double> start_mass_density_comoving;
  std::vector<double> start_momentum_density_x_comoving;
  std::vector<double> start_momentum_density_y_comoving;
  std::vector<double> start_momentum_density_z_comoving;
  std::vector<double> start_total_energy_density_comoving;
  std::vector<double> end_mass_density_comoving;
  std::vector<double> end_momentum_density_x_comoving;
  std::vector<double> end_momentum_density_y_comoving;
  std::vector<double> end_momentum_density_z_comoving;
  std::vector<double> end_total_energy_density_comoving;
  const auto records = store.records();
  patch_id.reserve(records.size());
  std::uint64_t offset = 0U;
  for (const core::AmrTemporalBoundaryHistoryRecord& record : records) {
    patch_id.push_back(record.patch_id);
    patch_level.push_back(record.patch_level);
    patch_geometry_fingerprint.push_back(record.patch_geometry_fingerprint);
    gas_cell_identity_generation.push_back(record.gas_cell_identity_generation);
    interval_start_code.push_back(record.interval_start_code);
    interval_end_code.push_back(record.interval_end_code);
    end_state_valid.push_back(record.end_state_valid ? 1U : 0U);
    cell_offset.push_back(offset);
    cell_count.push_back(static_cast<std::uint64_t>(record.cells.size()));
    for (const core::AmrTemporalBoundaryHistoryCellRecord& cell : record.cells) {
      gas_cell_id.push_back(cell.gas_cell_id);
      patch_local_cell.push_back(static_cast<std::uint64_t>(cell.patch_local_cell));
      start_mass_density_comoving.push_back(cell.start_mass_density_comoving);
      start_momentum_density_x_comoving.push_back(cell.start_momentum_density_x_comoving);
      start_momentum_density_y_comoving.push_back(cell.start_momentum_density_y_comoving);
      start_momentum_density_z_comoving.push_back(cell.start_momentum_density_z_comoving);
      start_total_energy_density_comoving.push_back(cell.start_total_energy_density_comoving);
      end_mass_density_comoving.push_back(cell.end_mass_density_comoving);
      end_momentum_density_x_comoving.push_back(cell.end_momentum_density_x_comoving);
      end_momentum_density_y_comoving.push_back(cell.end_momentum_density_y_comoving);
      end_momentum_density_z_comoving.push_back(cell.end_momentum_density_z_comoving);
      end_total_energy_density_comoving.push_back(cell.end_total_energy_density_comoving);
      ++offset;
    }
  }
  writeDataset1d(group.get(), "patch_id", H5T_STD_U64LE, H5T_NATIVE_UINT64, patch_id);
  writeDataset1d(group.get(), "patch_level", H5T_STD_U8LE, H5T_NATIVE_UINT8, patch_level);
  writeDataset1d(group.get(), "patch_geometry_fingerprint", H5T_STD_U64LE, H5T_NATIVE_UINT64, patch_geometry_fingerprint);
  writeDataset1d(group.get(), "gas_cell_identity_generation", H5T_STD_U64LE, H5T_NATIVE_UINT64, gas_cell_identity_generation);
  writeDataset1d(group.get(), "interval_start_code", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, interval_start_code);
  writeDataset1d(group.get(), "interval_end_code", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, interval_end_code);
  writeDataset1d(group.get(), "end_state_valid", H5T_STD_U8LE, H5T_NATIVE_UINT8, end_state_valid);
  writeDataset1d(group.get(), "cell_offset", H5T_STD_U64LE, H5T_NATIVE_UINT64, cell_offset);
  writeDataset1d(group.get(), "cell_count", H5T_STD_U64LE, H5T_NATIVE_UINT64, cell_count);
  writeDataset1d(group.get(), "gas_cell_id", H5T_STD_U64LE, H5T_NATIVE_UINT64, gas_cell_id);
  writeDataset1d(group.get(), "patch_local_cell", H5T_STD_U64LE, H5T_NATIVE_UINT64, patch_local_cell);
  writeDataset1d(group.get(), "start_mass_density_comoving", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, start_mass_density_comoving);
  writeDataset1d(group.get(), "start_momentum_density_x_comoving", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, start_momentum_density_x_comoving);
  writeDataset1d(group.get(), "start_momentum_density_y_comoving", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, start_momentum_density_y_comoving);
  writeDataset1d(group.get(), "start_momentum_density_z_comoving", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, start_momentum_density_z_comoving);
  writeDataset1d(group.get(), "start_total_energy_density_comoving", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, start_total_energy_density_comoving);
  writeDataset1d(group.get(), "end_mass_density_comoving", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, end_mass_density_comoving);
  writeDataset1d(group.get(), "end_momentum_density_x_comoving", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, end_momentum_density_x_comoving);
  writeDataset1d(group.get(), "end_momentum_density_y_comoving", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, end_momentum_density_y_comoving);
  writeDataset1d(group.get(), "end_momentum_density_z_comoving", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, end_momentum_density_z_comoving);
  writeDataset1d(group.get(), "end_total_energy_density_comoving", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, end_total_energy_density_comoving);
}

void readAmrTemporalBoundaryHistoryGroup(
    hid_t state_group,
    core::SimulationState& state,
    std::uint32_t schema_version) {
  if (schema_version < k_restart_schema_v19) {
    // Older checkpoints cannot prove a mid-subcycle continuation has the
    // time-aligned coarse history required by temporal coarse-to-fine fills.
    // A legacy file is therefore safe only when no deferred AMR synchronization
    // state remains. Pending registers are the conservative observable marker
    // for an open local subcycling interval in pre-v18 checkpoints.
    if (!state.pending_flux_registers.empty()) {
      throw std::runtime_error(
          "restart reader: legacy restart with pending AMR flux registers cannot resume "
          "without v18 temporal boundary history");
    }
    state.amr_temporal_boundary_history.clear();
    return;
  }
  Hdf5Handle group(H5Gopen2(state_group, "amr_temporal_boundary_history", H5P_DEFAULT));
  if (!group.valid()) {
    throw std::runtime_error("restart reader: v18 restart is missing /state/amr_temporal_boundary_history");
  }
  if (readScalarU32Attribute(group.get(), "schema_version") != 1U) {
    throw std::runtime_error("restart reader: unsupported /state/amr_temporal_boundary_history schema_version");
  }
  const auto patch_id = readDataset1d<std::uint64_t>(group.get(), "patch_id", H5T_NATIVE_UINT64);
  const auto patch_level = readDataset1d<std::uint8_t>(group.get(), "patch_level", H5T_NATIVE_UINT8);
  const auto patch_geometry_fingerprint = readDataset1d<std::uint64_t>(group.get(), "patch_geometry_fingerprint", H5T_NATIVE_UINT64);
  const auto gas_cell_identity_generation = readDataset1d<std::uint64_t>(group.get(), "gas_cell_identity_generation", H5T_NATIVE_UINT64);
  const auto interval_start_code = readDataset1d<double>(group.get(), "interval_start_code", H5T_NATIVE_DOUBLE);
  const auto interval_end_code = readDataset1d<double>(group.get(), "interval_end_code", H5T_NATIVE_DOUBLE);
  const auto end_state_valid = readDataset1d<std::uint8_t>(group.get(), "end_state_valid", H5T_NATIVE_UINT8);
  const auto cell_offset = readDataset1d<std::uint64_t>(group.get(), "cell_offset", H5T_NATIVE_UINT64);
  const auto cell_count = readDataset1d<std::uint64_t>(group.get(), "cell_count", H5T_NATIVE_UINT64);
  const auto gas_cell_id = readDataset1d<std::uint64_t>(group.get(), "gas_cell_id", H5T_NATIVE_UINT64);
  const auto patch_local_cell = readDataset1d<std::uint64_t>(group.get(), "patch_local_cell", H5T_NATIVE_UINT64);
  const auto start_mass_density_comoving = readDataset1d<double>(group.get(), "start_mass_density_comoving", H5T_NATIVE_DOUBLE);
  const auto start_momentum_density_x_comoving = readDataset1d<double>(group.get(), "start_momentum_density_x_comoving", H5T_NATIVE_DOUBLE);
  const auto start_momentum_density_y_comoving = readDataset1d<double>(group.get(), "start_momentum_density_y_comoving", H5T_NATIVE_DOUBLE);
  const auto start_momentum_density_z_comoving = readDataset1d<double>(group.get(), "start_momentum_density_z_comoving", H5T_NATIVE_DOUBLE);
  const auto start_total_energy_density_comoving = readDataset1d<double>(group.get(), "start_total_energy_density_comoving", H5T_NATIVE_DOUBLE);
  const auto end_mass_density_comoving = readDataset1d<double>(group.get(), "end_mass_density_comoving", H5T_NATIVE_DOUBLE);
  const auto end_momentum_density_x_comoving = readDataset1d<double>(group.get(), "end_momentum_density_x_comoving", H5T_NATIVE_DOUBLE);
  const auto end_momentum_density_y_comoving = readDataset1d<double>(group.get(), "end_momentum_density_y_comoving", H5T_NATIVE_DOUBLE);
  const auto end_momentum_density_z_comoving = readDataset1d<double>(group.get(), "end_momentum_density_z_comoving", H5T_NATIVE_DOUBLE);
  const auto end_total_energy_density_comoving = readDataset1d<double>(group.get(), "end_total_energy_density_comoving", H5T_NATIVE_DOUBLE);
  const std::size_t n = patch_id.size();
  const std::size_t nc = gas_cell_id.size();
  const bool record_sizes = patch_level.size() == n && patch_geometry_fingerprint.size() == n &&
      gas_cell_identity_generation.size() == n && interval_start_code.size() == n && interval_end_code.size() == n &&
      end_state_valid.size() == n && cell_offset.size() == n && cell_count.size() == n;
  const bool cell_sizes = patch_local_cell.size() == nc && start_mass_density_comoving.size() == nc &&
      start_momentum_density_x_comoving.size() == nc && start_momentum_density_y_comoving.size() == nc &&
      start_momentum_density_z_comoving.size() == nc && start_total_energy_density_comoving.size() == nc &&
      end_mass_density_comoving.size() == nc && end_momentum_density_x_comoving.size() == nc &&
      end_momentum_density_y_comoving.size() == nc && end_momentum_density_z_comoving.size() == nc &&
      end_total_energy_density_comoving.size() == nc;
  if (!record_sizes || !cell_sizes) {
    throw std::runtime_error("/state/amr_temporal_boundary_history datasets must have matching lengths");
  }
  std::vector<core::AmrTemporalBoundaryHistoryRecord> records;
  records.reserve(n);
  for (std::size_t i = 0; i < n; ++i) {
    if (end_state_valid[i] > 1U || cell_offset[i] > nc || cell_count[i] > nc - cell_offset[i] ||
        !std::isfinite(interval_start_code[i]) || !std::isfinite(interval_end_code[i]) ||
        interval_end_code[i] <= interval_start_code[i]) {
      throw std::runtime_error("/state/amr_temporal_boundary_history has invalid interval or record bounds");
    }
    core::AmrTemporalBoundaryHistoryRecord record;
    record.patch_id = patch_id[i];
    record.patch_level = patch_level[i];
    record.patch_geometry_fingerprint = patch_geometry_fingerprint[i];
    record.gas_cell_identity_generation = gas_cell_identity_generation[i];
    record.interval_start_code = interval_start_code[i];
    record.interval_end_code = interval_end_code[i];
    record.end_state_valid = end_state_valid[i] != 0U;
    const std::size_t begin = static_cast<std::size_t>(cell_offset[i]);
    const std::size_t count = static_cast<std::size_t>(cell_count[i]);
    record.cells.reserve(count);
    for (std::size_t j = 0; j < count; ++j) {
      const std::size_t k = begin + j;
      record.cells.push_back(core::AmrTemporalBoundaryHistoryCellRecord{
          .gas_cell_id = gas_cell_id[k],
          .patch_local_cell = static_cast<std::size_t>(patch_local_cell[k]),
          .start_mass_density_comoving = start_mass_density_comoving[k],
          .start_momentum_density_x_comoving = start_momentum_density_x_comoving[k],
          .start_momentum_density_y_comoving = start_momentum_density_y_comoving[k],
          .start_momentum_density_z_comoving = start_momentum_density_z_comoving[k],
          .start_total_energy_density_comoving = start_total_energy_density_comoving[k],
          .end_mass_density_comoving = end_mass_density_comoving[k],
          .end_momentum_density_x_comoving = end_momentum_density_x_comoving[k],
          .end_momentum_density_y_comoving = end_momentum_density_y_comoving[k],
          .end_momentum_density_z_comoving = end_momentum_density_z_comoving[k],
          .end_total_energy_density_comoving = end_total_energy_density_comoving[k]});
    }
    records.push_back(std::move(record));
  }
  state.amr_temporal_boundary_history.assign(std::move(records));
}

void validateAmrTemporalBoundaryHistoryForRestart(const core::SimulationState& state) {
  for (const core::AmrTemporalBoundaryHistoryRecord& record : state.amr_temporal_boundary_history.records()) {
    if (record.patch_id == 0U || record.cells.empty() || !record.end_state_valid ||
        record.gas_cell_identity_generation != state.gasCellIdentityGeneration() ||
        !std::isfinite(record.interval_start_code) || !std::isfinite(record.interval_end_code) ||
        record.interval_end_code <= record.interval_start_code) {
      throw std::runtime_error("restart reader: active AMR temporal history is incomplete or stale");
    }
    for (std::size_t i = 0; i < record.cells.size(); ++i) {
      const core::AmrTemporalBoundaryHistoryCellRecord& cell = record.cells[i];
      const auto finite_state = [](double value) { return std::isfinite(value); };
      if (cell.gas_cell_id == 0U || !finite_state(cell.start_mass_density_comoving) ||
          !finite_state(cell.start_momentum_density_x_comoving) ||
          !finite_state(cell.start_momentum_density_y_comoving) ||
          !finite_state(cell.start_momentum_density_z_comoving) ||
          !finite_state(cell.start_total_energy_density_comoving) ||
          !finite_state(cell.end_mass_density_comoving) ||
          !finite_state(cell.end_momentum_density_x_comoving) ||
          !finite_state(cell.end_momentum_density_y_comoving) ||
          !finite_state(cell.end_momentum_density_z_comoving) ||
          !finite_state(cell.end_total_energy_density_comoving) ||
          cell.start_mass_density_comoving <= 0.0 || cell.end_mass_density_comoving <= 0.0) {
        throw std::runtime_error("restart reader: temporal history contains inadmissible conserved state");
      }
      for (std::size_t j = 0; j < i; ++j) {
        if (record.cells[j].gas_cell_id == cell.gas_cell_id ||
            record.cells[j].patch_local_cell == cell.patch_local_cell) {
          throw std::runtime_error("restart reader: temporal history has duplicate stable cell identity or patch-local index");
        }
      }
      const auto row = state.rowForGasCellId(cell.gas_cell_id);
      if (!row.has_value() || state.cells.patch_index[*row] >= state.patches.size() ||
          state.patches.patch_id[state.cells.patch_index[*row]] != record.patch_id) {
        throw std::runtime_error("restart reader: AMR temporal history references a missing or moved coarse gas cell");
      }
    }
  }
}

void writeStateGroup(hid_t root, const core::SimulationState& state) {
  Hdf5Handle state_group(openOrCreateGroup(root, "/state"));
  Hdf5Handle particles_group(openOrCreateGroup(state_group.get(), "particles"));
  writeDataset1d(particles_group.get(), "position_x_comoving", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.particles.position_x_comoving);
  writeDataset1d(particles_group.get(), "position_y_comoving", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.particles.position_y_comoving);
  writeDataset1d(particles_group.get(), "position_z_comoving", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.particles.position_z_comoving);
  writeDataset1d(particles_group.get(), "velocity_x_peculiar", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.particles.velocity_x_peculiar);
  writeDataset1d(particles_group.get(), "velocity_y_peculiar", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.particles.velocity_y_peculiar);
  writeDataset1d(particles_group.get(), "velocity_z_peculiar", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.particles.velocity_z_peculiar);
  writeDataset1d(particles_group.get(), "mass_code", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.particles.mass_code);
  writeDataset1d(particles_group.get(), "time_bin", H5T_STD_U8LE, H5T_NATIVE_UINT8, state.particles.time_bin);

  Hdf5Handle particle_sidecar_group(openOrCreateGroup(state_group.get(), "particle_sidecar"));
  writeDataset1d(particle_sidecar_group.get(), "particle_id", H5T_STD_U64LE, H5T_NATIVE_UINT64, state.particle_sidecar.particle_id);
  writeDataset1d(particle_sidecar_group.get(), "sfc_key", H5T_STD_U64LE, H5T_NATIVE_UINT64, state.particle_sidecar.sfc_key);
  writeDataset1d(particle_sidecar_group.get(), "species_tag", H5T_STD_U32LE, H5T_NATIVE_UINT32, state.particle_sidecar.species_tag);
  writeDataset1d(particle_sidecar_group.get(), "particle_flags", H5T_STD_U32LE, H5T_NATIVE_UINT32, state.particle_sidecar.particle_flags);
  writeDataset1d(particle_sidecar_group.get(), "owning_rank", H5T_STD_U32LE, H5T_NATIVE_UINT32, state.particle_sidecar.owning_rank);
  writeDataset1d(particle_sidecar_group.get(), "last_drift_time_code", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.particle_sidecar.last_drift_time_code);
  writeDataset1d(particle_sidecar_group.get(), "last_drift_scale_factor", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.particle_sidecar.last_drift_scale_factor);
  writeDataset1d(
      particle_sidecar_group.get(),
      "gravity_softening_comoving",
      H5T_IEEE_F64LE,
      H5T_NATIVE_DOUBLE,
      state.particle_sidecar.gravity_softening_comoving);
  writeDataset1d(
      particle_sidecar_group.get(),
      "has_gravity_softening_override",
      H5T_STD_U8LE,
      H5T_NATIVE_UINT8,
      state.particle_sidecar.has_gravity_softening_override);

  Hdf5Handle cells_group(openOrCreateGroup(state_group.get(), "cells"));
  writeDataset1d(cells_group.get(), "center_x_comoving", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.cells.center_x_comoving);
  writeDataset1d(cells_group.get(), "center_y_comoving", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.cells.center_y_comoving);
  writeDataset1d(cells_group.get(), "center_z_comoving", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.cells.center_z_comoving);
  writeDataset1d(cells_group.get(), "mass_code", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.cells.mass_code);
  writeDataset1d(cells_group.get(), "time_bin", H5T_STD_U8LE, H5T_NATIVE_UINT8, state.cells.time_bin);
  writeDataset1d(cells_group.get(), "patch_index", H5T_STD_U32LE, H5T_NATIVE_UINT32, state.cells.patch_index);

  Hdf5Handle gas_group(openOrCreateGroup(state_group.get(), "gas_cells"));
  writeDataset1d(gas_group.get(), "gas_cell_id", H5T_STD_U64LE, H5T_NATIVE_UINT64, state.gas_cells.gas_cell_id);
  writeDataset1d(gas_group.get(), "parent_particle_id", H5T_STD_U64LE, H5T_NATIVE_UINT64, state.gas_cells.parent_particle_id);
  writeDataset1d(gas_group.get(), "velocity_x_peculiar", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.gas_cells.velocity_x_peculiar);
  writeDataset1d(gas_group.get(), "velocity_y_peculiar", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.gas_cells.velocity_y_peculiar);
  writeDataset1d(gas_group.get(), "velocity_z_peculiar", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.gas_cells.velocity_z_peculiar);
  writeDataset1d(gas_group.get(), "density_code", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.gas_cells.density_code);
  writeDataset1d(gas_group.get(), "pressure_code", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.gas_cells.pressure_code);
  writeDataset1d(gas_group.get(), "internal_energy_code", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.gas_cells.internal_energy_code);
  writeDataset1d(gas_group.get(), "temperature_code", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.gas_cells.temperature_code);
  writeDataset1d(gas_group.get(), "sound_speed_code", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.gas_cells.sound_speed_code);

  writeGasCellIdentityGroup(state_group.get(), state);

  Hdf5Handle patches_group(openOrCreateGroup(state_group.get(), "patches"));
  writeDataset1d(patches_group.get(), "patch_id", H5T_STD_U64LE, H5T_NATIVE_UINT64, state.patches.patch_id);
  writeDataset1d(patches_group.get(), "level", H5T_STD_I32LE, H5T_NATIVE_INT32, state.patches.level);
  writeDataset1d(patches_group.get(), "first_cell", H5T_STD_U32LE, H5T_NATIVE_UINT32, state.patches.first_cell);
  writeDataset1d(patches_group.get(), "cell_count", H5T_STD_U32LE, H5T_NATIVE_UINT32, state.patches.cell_count);
  writeDataset1d(patches_group.get(), "parent_patch_id", H5T_STD_U64LE, H5T_NATIVE_UINT64, state.patches.parent_patch_id);
  writeDataset1d(patches_group.get(), "morton_key", H5T_STD_U64LE, H5T_NATIVE_UINT64, state.patches.morton_key);
  writeDataset1d(patches_group.get(), "origin_x_comoving", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.patches.origin_x_comoving);
  writeDataset1d(patches_group.get(), "origin_y_comoving", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.patches.origin_y_comoving);
  writeDataset1d(patches_group.get(), "origin_z_comoving", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.patches.origin_z_comoving);
  writeDataset1d(patches_group.get(), "extent_x_comoving", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.patches.extent_x_comoving);
  writeDataset1d(patches_group.get(), "extent_y_comoving", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.patches.extent_y_comoving);
  writeDataset1d(patches_group.get(), "extent_z_comoving", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.patches.extent_z_comoving);
  writeDataset1d(patches_group.get(), "cell_dim_x", H5T_STD_U16LE, H5T_NATIVE_UINT16, state.patches.cell_dim_x);
  writeDataset1d(patches_group.get(), "cell_dim_y", H5T_STD_U16LE, H5T_NATIVE_UINT16, state.patches.cell_dim_y);
  writeDataset1d(patches_group.get(), "cell_dim_z", H5T_STD_U16LE, H5T_NATIVE_UINT16, state.patches.cell_dim_z);
  writeDataset1d(patches_group.get(), "owning_rank", H5T_STD_U32LE, H5T_NATIVE_UINT32, state.patches.owning_rank);

  writePendingFluxRegisterGroup(state_group.get(), state.pending_flux_registers);
  writeAmrTemporalBoundaryHistoryGroup(state_group.get(), state.amr_temporal_boundary_history);

  writeDataset1d(state_group.get(), "species_count_by_species", H5T_STD_U64LE, H5T_NATIVE_UINT64, std::vector<std::uint64_t>(state.species.count_by_species.begin(), state.species.count_by_species.end()));

  writeStarSidecarGroup(state_group.get(), state.star_particles);

  Hdf5Handle bh_group(openOrCreateGroup(state_group.get(), "black_holes"));
  writeDataset1d(bh_group.get(), "particle_index", H5T_STD_U32LE, H5T_NATIVE_UINT32, state.black_holes.particle_index);
  writeDataset1d(bh_group.get(), "host_cell_index", H5T_STD_U32LE, H5T_NATIVE_UINT32, state.black_holes.host_cell_index);
  writeDataset1d(bh_group.get(), "subgrid_mass_code", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.black_holes.subgrid_mass_code);
  writeDataset1d(bh_group.get(), "accretion_rate_code", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.black_holes.accretion_rate_code);
  writeDataset1d(bh_group.get(), "feedback_energy_code", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.black_holes.feedback_energy_code);
  writeDataset1d(bh_group.get(), "eddington_ratio", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.black_holes.eddington_ratio);
  writeDataset1d(
      bh_group.get(),
      "cumulative_accreted_mass_code",
      H5T_IEEE_F64LE,
      H5T_NATIVE_DOUBLE,
      state.black_holes.cumulative_accreted_mass_code);
  writeDataset1d(
      bh_group.get(),
      "cumulative_feedback_energy_code",
      H5T_IEEE_F64LE,
      H5T_NATIVE_DOUBLE,
      state.black_holes.cumulative_feedback_energy_code);
  writeDataset1d(
      bh_group.get(),
      "duty_cycle_active_time_code",
      H5T_IEEE_F64LE,
      H5T_NATIVE_DOUBLE,
      state.black_holes.duty_cycle_active_time_code);
  writeDataset1d(
      bh_group.get(),
      "duty_cycle_total_time_code",
      H5T_IEEE_F64LE,
      H5T_NATIVE_DOUBLE,
      state.black_holes.duty_cycle_total_time_code);

  Hdf5Handle tracer_group(openOrCreateGroup(state_group.get(), "tracers"));
  writeDataset1d(tracer_group.get(), "particle_index", H5T_STD_U32LE, H5T_NATIVE_UINT32, state.tracers.particle_index);
  writeDataset1d(tracer_group.get(), "parent_particle_id", H5T_STD_U64LE, H5T_NATIVE_UINT64, state.tracers.parent_particle_id);
  writeDataset1d(tracer_group.get(), "injection_step", H5T_STD_U64LE, H5T_NATIVE_UINT64, state.tracers.injection_step);
  writeDataset1d(tracer_group.get(), "host_cell_index", H5T_STD_U32LE, H5T_NATIVE_UINT32, state.tracers.host_cell_index);
  writeDataset1d(tracer_group.get(), "mass_fraction_of_host", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.tracers.mass_fraction_of_host);
  writeDataset1d(tracer_group.get(), "last_host_mass_code", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, state.tracers.last_host_mass_code);
  writeDataset1d(
      tracer_group.get(),
      "cumulative_exchanged_mass_code",
      H5T_IEEE_F64LE,
      H5T_NATIVE_DOUBLE,
      state.tracers.cumulative_exchanged_mass_code);

  writeStringDataset(state_group.get(), "state_metadata", state.metadata.serialize());

  Hdf5Handle sidecars_group(openOrCreateGroup(state_group.get(), "module_sidecars"));
  const auto sidecar_blocks = state.sidecars.blocksSortedByName();
  std::string module_names_text;
  for (const core::ModuleSidecarBlock* block : sidecar_blocks) {
    Hdf5Handle module_group(openOrCreateGroup(sidecars_group.get(), block->module_name));
    writeScalarU32Attribute(module_group.get(), "schema_version", block->schema_version);
    writeScalarU32Attribute(module_group.get(), "particle_indexed", block->particle_indexed ? 1U : 0U);
    writeScalarU32Attribute(module_group.get(), "row_stride_bytes", block->row_stride_bytes);
    writeScalarU32Attribute(module_group.get(), "required_species_mask", block->required_species_mask);
    writeScalarU32Attribute(module_group.get(), "requirement_kind", static_cast<std::uint32_t>(block->requirement.kind));
    writeScalarU32Attribute(module_group.get(), "requirement_species_mask", block->requirement.species_mask);
    writeScalarU32Attribute(module_group.get(), "requirement_particle_flags_mask", block->requirement.particle_flags_mask);
    writeScalarF64Attribute(module_group.get(), "requirement_threshold_code", block->requirement.threshold_code);
    std::vector<std::uint8_t> payload(block->payload.size());
    for (std::size_t i = 0; i < block->payload.size(); ++i) {
      payload[i] = std::to_integer<std::uint8_t>(block->payload[i]);
    }
    writeDataset1d(module_group.get(), "payload", H5T_STD_U8LE, H5T_NATIVE_UINT8, payload);
    if (block->isParticleIndexed()) {
      writeDataset1d(
          module_group.get(),
          "particle_id_by_row",
          H5T_STD_U64LE,
          H5T_NATIVE_UINT64,
          block->particle_id_by_row);
    }
    module_names_text += block->module_name;
    module_names_text.push_back('\n');
  }
  writeStringDataset(state_group.get(), "module_sidecar_names", module_names_text);
}

void readStateGroup(hid_t root, core::SimulationState& state, std::uint32_t schema_version) {
  Hdf5Handle state_group(H5Gopen2(root, "/state", H5P_DEFAULT));
  Hdf5Handle particles_group(H5Gopen2(state_group.get(), "particles", H5P_DEFAULT));
  state.particles.position_x_comoving = readDataset1dAligned<double>(particles_group.get(), "position_x_comoving", H5T_NATIVE_DOUBLE);
  state.particles.position_y_comoving = readDataset1dAligned<double>(particles_group.get(), "position_y_comoving", H5T_NATIVE_DOUBLE);
  state.particles.position_z_comoving = readDataset1dAligned<double>(particles_group.get(), "position_z_comoving", H5T_NATIVE_DOUBLE);
  state.particles.velocity_x_peculiar = readDataset1dAligned<double>(particles_group.get(), "velocity_x_peculiar", H5T_NATIVE_DOUBLE);
  state.particles.velocity_y_peculiar = readDataset1dAligned<double>(particles_group.get(), "velocity_y_peculiar", H5T_NATIVE_DOUBLE);
  state.particles.velocity_z_peculiar = readDataset1dAligned<double>(particles_group.get(), "velocity_z_peculiar", H5T_NATIVE_DOUBLE);
  state.particles.mass_code = readDataset1dAligned<double>(particles_group.get(), "mass_code", H5T_NATIVE_DOUBLE);
  state.particles.time_bin = readDataset1dAligned<std::uint8_t>(particles_group.get(), "time_bin", H5T_NATIVE_UINT8);

  Hdf5Handle particle_sidecar_group(H5Gopen2(state_group.get(), "particle_sidecar", H5P_DEFAULT));
  state.particle_sidecar.particle_id =
      readDataset1dAligned<std::uint64_t>(particle_sidecar_group.get(), "particle_id", H5T_NATIVE_UINT64);
  state.particle_sidecar.sfc_key =
      readDataset1dAligned<std::uint64_t>(particle_sidecar_group.get(), "sfc_key", H5T_NATIVE_UINT64);
  state.particle_sidecar.species_tag =
      readDataset1dAligned<std::uint32_t>(particle_sidecar_group.get(), "species_tag", H5T_NATIVE_UINT32);
  state.particle_sidecar.particle_flags =
      readDataset1dAligned<std::uint32_t>(particle_sidecar_group.get(), "particle_flags", H5T_NATIVE_UINT32);
  state.particle_sidecar.owning_rank =
      readDataset1dAligned<std::uint32_t>(particle_sidecar_group.get(), "owning_rank", H5T_NATIVE_UINT32);
  state.particle_sidecar.last_drift_time_code =
      readDataset1dAligned<double>(particle_sidecar_group.get(), "last_drift_time_code", H5T_NATIVE_DOUBLE);
  state.particle_sidecar.last_drift_scale_factor =
      readDataset1dAligned<double>(particle_sidecar_group.get(), "last_drift_scale_factor", H5T_NATIVE_DOUBLE);
  state.particle_sidecar.gravity_softening_comoving =
      readDataset1dAligned<double>(particle_sidecar_group.get(), "gravity_softening_comoving", H5T_NATIVE_DOUBLE);
  state.particle_sidecar.has_gravity_softening_override =
      readDataset1dAligned<std::uint8_t>(
          particle_sidecar_group.get(), "has_gravity_softening_override", H5T_NATIVE_UINT8);

  Hdf5Handle cells_group(H5Gopen2(state_group.get(), "cells", H5P_DEFAULT));
  state.cells.center_x_comoving = readDataset1dAligned<double>(cells_group.get(), "center_x_comoving", H5T_NATIVE_DOUBLE);
  state.cells.center_y_comoving = readDataset1dAligned<double>(cells_group.get(), "center_y_comoving", H5T_NATIVE_DOUBLE);
  state.cells.center_z_comoving = readDataset1dAligned<double>(cells_group.get(), "center_z_comoving", H5T_NATIVE_DOUBLE);
  state.cells.mass_code = readDataset1dAligned<double>(cells_group.get(), "mass_code", H5T_NATIVE_DOUBLE);
  state.cells.time_bin = readDataset1dAligned<std::uint8_t>(cells_group.get(), "time_bin", H5T_NATIVE_UINT8);
  state.cells.patch_index = readDataset1dAligned<std::uint32_t>(cells_group.get(), "patch_index", H5T_NATIVE_UINT32);

  Hdf5Handle gas_group(H5Gopen2(state_group.get(), "gas_cells", H5P_DEFAULT));
  if (H5Lexists(gas_group.get(), "gas_cell_id", H5P_DEFAULT) > 0 &&
      H5Lexists(gas_group.get(), "parent_particle_id", H5P_DEFAULT) > 0) {
    state.gas_cells.gas_cell_id =
        readDataset1dAligned<std::uint64_t>(gas_group.get(), "gas_cell_id", H5T_NATIVE_UINT64);
    state.gas_cells.parent_particle_id =
        readDataset1dAligned<std::uint64_t>(gas_group.get(), "parent_particle_id", H5T_NATIVE_UINT64);
  } else {
    throw std::runtime_error("restart is missing gas-cell identity datasets gas_cell_id/parent_particle_id");
  }
  state.gas_cells.density_code = readDataset1dAligned<double>(gas_group.get(), "density_code", H5T_NATIVE_DOUBLE);
  state.gas_cells.velocity_x_peculiar =
      readDataset1dAligned<double>(gas_group.get(), "velocity_x_peculiar", H5T_NATIVE_DOUBLE);
  state.gas_cells.velocity_y_peculiar =
      readDataset1dAligned<double>(gas_group.get(), "velocity_y_peculiar", H5T_NATIVE_DOUBLE);
  state.gas_cells.velocity_z_peculiar =
      readDataset1dAligned<double>(gas_group.get(), "velocity_z_peculiar", H5T_NATIVE_DOUBLE);
  state.gas_cells.pressure_code = readDataset1dAligned<double>(gas_group.get(), "pressure_code", H5T_NATIVE_DOUBLE);
  state.gas_cells.internal_energy_code = readDataset1dAligned<double>(gas_group.get(), "internal_energy_code", H5T_NATIVE_DOUBLE);
  state.gas_cells.temperature_code = readDataset1dAligned<double>(gas_group.get(), "temperature_code", H5T_NATIVE_DOUBLE);
  state.gas_cells.sound_speed_code = readDataset1dAligned<double>(gas_group.get(), "sound_speed_code", H5T_NATIVE_DOUBLE);

  Hdf5Handle patches_group(H5Gopen2(state_group.get(), "patches", H5P_DEFAULT));
  state.patches.patch_id = readDataset1dAligned<std::uint64_t>(patches_group.get(), "patch_id", H5T_NATIVE_UINT64);
  state.patches.level = readDataset1dAligned<std::int32_t>(patches_group.get(), "level", H5T_NATIVE_INT32);
  state.patches.first_cell = readDataset1dAligned<std::uint32_t>(patches_group.get(), "first_cell", H5T_NATIVE_UINT32);
  state.patches.cell_count = readDataset1dAligned<std::uint32_t>(patches_group.get(), "cell_count", H5T_NATIVE_UINT32);
  if (H5Lexists(patches_group.get(), "parent_patch_id", H5P_DEFAULT) > 0) {
    state.patches.parent_patch_id = readDataset1dAligned<std::uint64_t>(patches_group.get(), "parent_patch_id", H5T_NATIVE_UINT64);
    state.patches.morton_key = readDataset1dAligned<std::uint64_t>(patches_group.get(), "morton_key", H5T_NATIVE_UINT64);
    state.patches.origin_x_comoving = readDataset1dAligned<double>(patches_group.get(), "origin_x_comoving", H5T_NATIVE_DOUBLE);
    state.patches.origin_y_comoving = readDataset1dAligned<double>(patches_group.get(), "origin_y_comoving", H5T_NATIVE_DOUBLE);
    state.patches.origin_z_comoving = readDataset1dAligned<double>(patches_group.get(), "origin_z_comoving", H5T_NATIVE_DOUBLE);
    state.patches.extent_x_comoving = readDataset1dAligned<double>(patches_group.get(), "extent_x_comoving", H5T_NATIVE_DOUBLE);
    state.patches.extent_y_comoving = readDataset1dAligned<double>(patches_group.get(), "extent_y_comoving", H5T_NATIVE_DOUBLE);
    state.patches.extent_z_comoving = readDataset1dAligned<double>(patches_group.get(), "extent_z_comoving", H5T_NATIVE_DOUBLE);
    state.patches.cell_dim_x = readDataset1dAligned<std::uint16_t>(patches_group.get(), "cell_dim_x", H5T_NATIVE_UINT16);
    state.patches.cell_dim_y = readDataset1dAligned<std::uint16_t>(patches_group.get(), "cell_dim_y", H5T_NATIVE_UINT16);
    state.patches.cell_dim_z = readDataset1dAligned<std::uint16_t>(patches_group.get(), "cell_dim_z", H5T_NATIVE_UINT16);
  } else {
    state.patches.parent_patch_id.assign(state.patches.patch_id.size(), 0U);
    state.patches.morton_key.assign(state.patches.patch_id.size(), 0U);
    state.patches.origin_x_comoving.assign(state.patches.patch_id.size(), 0.0);
    state.patches.origin_y_comoving.assign(state.patches.patch_id.size(), 0.0);
    state.patches.origin_z_comoving.assign(state.patches.patch_id.size(), 0.0);
    state.patches.extent_x_comoving.assign(state.patches.patch_id.size(), 0.0);
    state.patches.extent_y_comoving.assign(state.patches.patch_id.size(), 0.0);
    state.patches.extent_z_comoving.assign(state.patches.patch_id.size(), 0.0);
    state.patches.cell_dim_x.assign(state.patches.patch_id.size(), 0U);
    state.patches.cell_dim_y.assign(state.patches.patch_id.size(), 0U);
    state.patches.cell_dim_z.assign(state.patches.patch_id.size(), 0U);
  }
  if (H5Lexists(patches_group.get(), "owning_rank", H5P_DEFAULT) > 0) {
    state.patches.owning_rank = readDataset1dAligned<std::uint32_t>(patches_group.get(), "owning_rank", H5T_NATIVE_UINT32);
  } else {
    state.patches.owning_rank.assign(state.patches.patch_id.size(), 0U);
  }

  readGasCellIdentityGroup(state_group.get(), state, schema_version);
  readPendingFluxRegisterGroup(state_group.get(), state, schema_version);
  readAmrTemporalBoundaryHistoryGroup(state_group.get(), state, schema_version);

  const auto species_count = readDataset1d<std::uint64_t>(state_group.get(), "species_count_by_species", H5T_NATIVE_UINT64);
  if (species_count.size() != state.species.count_by_species.size()) {
    throw std::runtime_error("invalid species_count_by_species size in restart");
  }
  std::copy(species_count.begin(), species_count.end(), state.species.count_by_species.begin());

  readStarSidecarGroup(state_group.get(), state.star_particles);

  Hdf5Handle bh_group(H5Gopen2(state_group.get(), "black_holes", H5P_DEFAULT));
  state.black_holes.particle_index = readDataset1dAligned<std::uint32_t>(bh_group.get(), "particle_index", H5T_NATIVE_UINT32);
  state.black_holes.host_cell_index = readDataset1dAligned<std::uint32_t>(bh_group.get(), "host_cell_index", H5T_NATIVE_UINT32);
  state.black_holes.subgrid_mass_code = readDataset1dAligned<double>(bh_group.get(), "subgrid_mass_code", H5T_NATIVE_DOUBLE);
  state.black_holes.accretion_rate_code =
      readDataset1dAligned<double>(bh_group.get(), "accretion_rate_code", H5T_NATIVE_DOUBLE);
  state.black_holes.feedback_energy_code =
      readDataset1dAligned<double>(bh_group.get(), "feedback_energy_code", H5T_NATIVE_DOUBLE);
  state.black_holes.eddington_ratio = readDataset1dAligned<double>(bh_group.get(), "eddington_ratio", H5T_NATIVE_DOUBLE);
  state.black_holes.cumulative_accreted_mass_code =
      readDataset1dAligned<double>(bh_group.get(), "cumulative_accreted_mass_code", H5T_NATIVE_DOUBLE);
  state.black_holes.cumulative_feedback_energy_code =
      readDataset1dAligned<double>(bh_group.get(), "cumulative_feedback_energy_code", H5T_NATIVE_DOUBLE);
  state.black_holes.duty_cycle_active_time_code =
      readDataset1dAligned<double>(bh_group.get(), "duty_cycle_active_time_code", H5T_NATIVE_DOUBLE);
  state.black_holes.duty_cycle_total_time_code =
      readDataset1dAligned<double>(bh_group.get(), "duty_cycle_total_time_code", H5T_NATIVE_DOUBLE);

  Hdf5Handle tracer_group(H5Gopen2(state_group.get(), "tracers", H5P_DEFAULT));
  state.tracers.particle_index = readDataset1dAligned<std::uint32_t>(tracer_group.get(), "particle_index", H5T_NATIVE_UINT32);
  state.tracers.parent_particle_id =
      readDataset1dAligned<std::uint64_t>(tracer_group.get(), "parent_particle_id", H5T_NATIVE_UINT64);
  state.tracers.injection_step =
      readDataset1dAligned<std::uint64_t>(tracer_group.get(), "injection_step", H5T_NATIVE_UINT64);
  state.tracers.host_cell_index =
      readDataset1dAligned<std::uint32_t>(tracer_group.get(), "host_cell_index", H5T_NATIVE_UINT32);
  state.tracers.mass_fraction_of_host =
      readDataset1dAligned<double>(tracer_group.get(), "mass_fraction_of_host", H5T_NATIVE_DOUBLE);
  state.tracers.last_host_mass_code =
      readDataset1dAligned<double>(tracer_group.get(), "last_host_mass_code", H5T_NATIVE_DOUBLE);
  state.tracers.cumulative_exchanged_mass_code =
      readDataset1dAligned<double>(tracer_group.get(), "cumulative_exchanged_mass_code", H5T_NATIVE_DOUBLE);

  state.metadata = core::StateMetadata::deserialize(readStringDataset(state_group.get(), "state_metadata"));
  state.rebuildSpeciesIndex();

  Hdf5Handle sidecars_group(H5Gopen2(state_group.get(), "module_sidecars", H5P_DEFAULT));
  const std::string module_names_text = readStringDataset(state_group.get(), "module_sidecar_names");
  std::istringstream module_names_stream(module_names_text);
  std::string module_name;
  while (std::getline(module_names_stream, module_name)) {
    if (module_name.empty()) {
      continue;
    }
    Hdf5Handle module_group(H5Gopen2(sidecars_group.get(), module_name.c_str(), H5P_DEFAULT));
    core::ModuleSidecarBlock block;
    block.module_name = module_name;
    block.schema_version = readScalarU32Attribute(module_group.get(), "schema_version");
    if (hdf5AttributeExists(module_group.get(), "particle_indexed")) {
      block.particle_indexed = readScalarU32Attribute(module_group.get(), "particle_indexed") != 0U;
    }
    if (hdf5AttributeExists(module_group.get(), "row_stride_bytes")) {
      block.row_stride_bytes = readScalarU32Attribute(module_group.get(), "row_stride_bytes");
    }
    if (hdf5AttributeExists(module_group.get(), "required_species_mask")) {
      block.required_species_mask = readScalarU32Attribute(module_group.get(), "required_species_mask");
    }
    if (hdf5AttributeExists(module_group.get(), "requirement_kind")) {
      block.requirement.kind = static_cast<core::ModuleSidecarRequirementKind>(
          readScalarU32Attribute(module_group.get(), "requirement_kind"));
    }
    if (hdf5AttributeExists(module_group.get(), "requirement_species_mask")) {
      block.requirement.species_mask = readScalarU32Attribute(module_group.get(), "requirement_species_mask");
    }
    if (hdf5AttributeExists(module_group.get(), "requirement_particle_flags_mask")) {
      block.requirement.particle_flags_mask = readScalarU32Attribute(module_group.get(), "requirement_particle_flags_mask");
    }
    if (hdf5AttributeExists(module_group.get(), "requirement_threshold_code")) {
      block.requirement.threshold_code = readScalarF64Attribute(module_group.get(), "requirement_threshold_code");
    }
    const auto payload_u8 = readDataset1d<std::uint8_t>(module_group.get(), "payload", H5T_NATIVE_UINT8);
    block.payload.resize(payload_u8.size());
    for (std::size_t i = 0; i < payload_u8.size(); ++i) {
      block.payload[i] = static_cast<std::byte>(payload_u8[i]);
    }
    if (block.particle_indexed || block.row_stride_bytes != 0U ||
        hdf5LinkExists(module_group.get(), "particle_id_by_row")) {
      block.particle_indexed = true;
      block.particle_id_by_row =
          readDataset1d<std::uint64_t>(module_group.get(), "particle_id_by_row", H5T_NATIVE_UINT64);
    }
    state.sidecars.upsert(std::move(block));
  }

  validateHydroGeometryStateForRestart(state, "restart reader");
  if (!state.validateOwnershipInvariants()) {
    throw std::runtime_error("restart state failed ownership invariant validation");
  }
}

void writeDistributedGravityGroup(hid_t root, const parallel::DistributedRestartState& distributed_state) {
  Hdf5Handle group(openOrCreateGroup(root, "/distributed_gravity"));
  writeStringDataset(group.get(), "state", distributed_state.serialize());
}

void readDistributedGravityGroup(hid_t root, parallel::DistributedRestartState& distributed_state) {
  Hdf5Handle group(H5Gopen2(root, "/distributed_gravity", H5P_DEFAULT));
  distributed_state = parallel::DistributedRestartState::deserialize(readStringDataset(group.get(), "state"));
}


void writeOutputCadenceGroup(hid_t root, const OutputCadencePersistentState& output_state) {
  Hdf5Handle group(openOrCreateGroup(root, "/output_cadence"));
  writeScalarU32Attribute(group.get(), "output_enabled", output_state.output_enabled ? 1U : 0U);
  writeScalarU32Attribute(group.get(), "write_restarts", output_state.write_restarts ? 1U : 0U);
  writeScalarU32Attribute(group.get(), "snapshot_due", output_state.snapshot_due ? 1U : 0U);
  writeScalarU32Attribute(group.get(), "checkpoint_due", output_state.checkpoint_due ? 1U : 0U);
  writeScalarU64Attribute(group.get(), "last_completed_step_index", output_state.last_completed_step_index);
  writeScalarU64Attribute(group.get(), "snapshot_interval_steps", output_state.snapshot_interval_steps);
  writeScalarU64Attribute(group.get(), "next_snapshot_step_index", output_state.next_snapshot_step_index);
  writeScalarStringAttribute(group.get(), "snapshot_stem", output_state.snapshot_stem);
  writeScalarStringAttribute(group.get(), "restart_stem", output_state.restart_stem);
}

[[nodiscard]] OutputCadencePersistentState readOutputCadenceGroup(hid_t root) {
  Hdf5Handle group(H5Gopen2(root, "/output_cadence", H5P_DEFAULT));
  OutputCadencePersistentState output_state;
  output_state.output_enabled = readScalarU32Attribute(group.get(), "output_enabled") != 0U;
  output_state.write_restarts = readScalarU32Attribute(group.get(), "write_restarts") != 0U;
  output_state.snapshot_due = readScalarU32Attribute(group.get(), "snapshot_due") != 0U;
  output_state.checkpoint_due = readScalarU32Attribute(group.get(), "checkpoint_due") != 0U;
  output_state.last_completed_step_index = readScalarU64Attribute(group.get(), "last_completed_step_index");
  output_state.snapshot_interval_steps = readScalarU64Attribute(group.get(), "snapshot_interval_steps");
  output_state.next_snapshot_step_index = readScalarU64Attribute(group.get(), "next_snapshot_step_index");
  output_state.snapshot_stem = readScalarStringAttribute(group.get(), "snapshot_stem");
  output_state.restart_stem = readScalarStringAttribute(group.get(), "restart_stem");
  return output_state;
}

void writeStochasticStateGroup(hid_t root, const StochasticPersistentState& stochastic_state) {
  Hdf5Handle group(openOrCreateGroup(root, "/stochastic_state"));
  const auto modules = sortedStochasticModules(stochastic_state);
  writeScalarU32Attribute(group.get(), "module_count", static_cast<std::uint32_t>(modules.size()));
  std::string module_names;
  for (const StochasticModulePersistentState& module_state : modules) {
    if (!module_names.empty()) {
      module_names.push_back('\n');
    }
    module_names += module_state.module_name;
    Hdf5Handle module_group(openOrCreateGroup(root, "/stochastic_state/" + module_state.module_name));
    writeScalarU32Attribute(module_group.get(), "schema_version", module_state.schema_version);
    writeScalarStringAttribute(module_group.get(), "rng_policy", module_state.rng_policy);
    writeScalarU64Attribute(module_group.get(), "random_seed", module_state.random_seed);
    writeScalarU32Attribute(module_group.get(), "rank_local_seed_offset", module_state.rank_local_seed_offset);
    writeScalarU64Attribute(module_group.get(), "last_committed_step_index", module_state.last_committed_step_index);
    writeScalarU32Attribute(
        module_group.get(),
        "deterministic_from_serialized_inputs",
        module_state.deterministic_from_serialized_inputs ? 1U : 0U);
  }
  writeStringDataset(group.get(), "module_names", module_names);
}

[[nodiscard]] StochasticPersistentState readStochasticStateGroup(hid_t root) {
  Hdf5Handle group(H5Gopen2(root, "/stochastic_state", H5P_DEFAULT));
  if (!group.valid()) {
    throw std::runtime_error("restart schema validation missing required group: /stochastic_state");
  }
  const std::string module_names_text = readStringDataset(group.get(), "module_names");
  StochasticPersistentState stochastic_state;
  std::istringstream module_names_stream(module_names_text);
  std::string module_name;
  while (std::getline(module_names_stream, module_name)) {
    if (module_name.empty()) {
      continue;
    }
    Hdf5Handle module_group(H5Gopen2(group.get(), module_name.c_str(), H5P_DEFAULT));
    if (!module_group.valid()) {
      throw std::runtime_error("restart schema validation missing required group: /stochastic_state/" + module_name);
    }
    StochasticModulePersistentState module_state;
    module_state.module_name = module_name;
    module_state.schema_version = readScalarU32Attribute(module_group.get(), "schema_version");
    module_state.rng_policy = readScalarStringAttribute(module_group.get(), "rng_policy");
    module_state.random_seed = readScalarU64Attribute(module_group.get(), "random_seed");
    module_state.rank_local_seed_offset = readScalarU32Attribute(module_group.get(), "rank_local_seed_offset");
    module_state.last_committed_step_index = readScalarU64Attribute(module_group.get(), "last_committed_step_index");
    module_state.deterministic_from_serialized_inputs =
        readScalarU32Attribute(module_group.get(), "deterministic_from_serialized_inputs") != 0U;
    stochastic_state.modules.push_back(std::move(module_state));
  }
  const std::uint32_t module_count = readScalarU32Attribute(group.get(), "module_count");
  if (module_count != stochastic_state.modules.size()) {
    throw std::runtime_error("restart schema validation invalid /stochastic_state/@module_count");
  }
  return stochastic_state;
}

[[nodiscard]] RestartDiagnosticsSummary makeRestartDiagnosticsSummary(
    const core::IntegratorState& integrator_state,
    const core::TimeBinPersistentState& scheduler_state,
    const core::TimeBinPersistentState& gas_cell_scheduler_state,
    const OutputCadencePersistentState& output_state,
    const StochasticPersistentState& stochastic_state) {
  RestartDiagnosticsSummary diagnostics;
  diagnostics.restart_schema_name = restartSchema().name;
  diagnostics.restart_schema_version = restartSchema().version;
  diagnostics.current_boundary_kind = std::string(core::stepBoundaryKindName(integrator_state.current_boundary_kind));
  diagnostics.last_completed_boundary_kind =
      std::string(core::stepBoundaryKindName(integrator_state.last_completed_boundary_kind));
  diagnostics.restart_safe = !integrator_state.inside_kdk_step && integrator_state.last_completed_restart_safe &&
      core::isRestartSafeBoundary(integrator_state.last_completed_boundary_kind);
  diagnostics.step_index = integrator_state.step_index;
  diagnostics.scheduler_current_tick = scheduler_state.current_tick;
  diagnostics.scheduler_max_bin = scheduler_state.max_bin;
  diagnostics.scheduler_element_count = scheduler_state.bin_index.size();
  diagnostics.scheduler_active_count = static_cast<std::uint64_t>(
      std::count(scheduler_state.active_flag.begin(), scheduler_state.active_flag.end(), static_cast<std::uint8_t>(1U)));
  diagnostics.scheduler_pending_transition_count = static_cast<std::uint64_t>(std::count_if(
      scheduler_state.pending_bin_index.begin(),
      scheduler_state.pending_bin_index.end(),
      [](std::uint8_t pending_bin) {
        return pending_bin != core::HierarchicalTimeBinScheduler::k_unset_pending_bin;
      }));
  diagnostics.gas_cell_scheduler_current_tick = gas_cell_scheduler_state.current_tick;
  diagnostics.gas_cell_scheduler_max_bin = gas_cell_scheduler_state.max_bin;
  diagnostics.gas_cell_scheduler_element_count = gas_cell_scheduler_state.bin_index.size();
  diagnostics.gas_cell_scheduler_active_count = static_cast<std::uint64_t>(
      std::count(gas_cell_scheduler_state.active_flag.begin(), gas_cell_scheduler_state.active_flag.end(), static_cast<std::uint8_t>(1U)));
  diagnostics.gas_cell_scheduler_pending_transition_count = static_cast<std::uint64_t>(std::count_if(
      gas_cell_scheduler_state.pending_bin_index.begin(),
      gas_cell_scheduler_state.pending_bin_index.end(),
      [](std::uint8_t pending_bin) {
        return pending_bin != core::HierarchicalTimeBinScheduler::k_unset_pending_bin;
      }));

  const core::PmSynchronizationPersistentState pm_sync_state =
      integrator_state.pm_sync_state.exportPersistentState();
  diagnostics.pm_cadence_steps = pm_sync_state.cadence_steps;
  diagnostics.pm_gravity_kick_opportunity = pm_sync_state.gravity_kick_opportunity;
  diagnostics.pm_field_version = pm_sync_state.field_version;
  diagnostics.pm_last_refresh_opportunity = pm_sync_state.last_refresh_opportunity;
  diagnostics.pm_last_refresh_step_index = pm_sync_state.last_refresh_step_index;
  diagnostics.pm_refresh_commit_pending = pm_sync_state.refresh_commit_pending;
  diagnostics.pm_long_range_field_valid = integrator_state.pm_long_range_field_valid;

  diagnostics.output_enabled = output_state.output_enabled;
  diagnostics.output_snapshot_due = output_state.snapshot_due;
  diagnostics.output_checkpoint_due = output_state.checkpoint_due;
  diagnostics.output_last_completed_step_index = output_state.last_completed_step_index;
  diagnostics.output_next_snapshot_step_index = output_state.next_snapshot_step_index;
  diagnostics.stochastic_module_count = stochastic_state.modules.size();
  return diagnostics;
}

void writeRestartDiagnosticsGroup(hid_t root, const RestartDiagnosticsSummary& diagnostics) {
  Hdf5Handle group(openOrCreateGroup(root, "/restart_diagnostics"));
  writeScalarStringAttribute(group.get(), "restart_schema_name", diagnostics.restart_schema_name);
  writeScalarU32Attribute(group.get(), "restart_schema_version", diagnostics.restart_schema_version);
  writeScalarStringAttribute(group.get(), "current_boundary_kind", diagnostics.current_boundary_kind);
  writeScalarStringAttribute(group.get(), "last_completed_boundary_kind", diagnostics.last_completed_boundary_kind);
  writeScalarU32Attribute(group.get(), "restart_safe", diagnostics.restart_safe ? 1U : 0U);
  writeScalarU64Attribute(group.get(), "step_index", diagnostics.step_index);
  writeScalarU64Attribute(group.get(), "scheduler_current_tick", diagnostics.scheduler_current_tick);
  writeScalarU32Attribute(group.get(), "scheduler_max_bin", diagnostics.scheduler_max_bin);
  writeScalarU64Attribute(group.get(), "scheduler_element_count", diagnostics.scheduler_element_count);
  writeScalarU64Attribute(group.get(), "scheduler_active_count", diagnostics.scheduler_active_count);
  writeScalarU64Attribute(
      group.get(), "scheduler_pending_transition_count", diagnostics.scheduler_pending_transition_count);
  writeScalarU64Attribute(
      group.get(), "gas_cell_scheduler_current_tick", diagnostics.gas_cell_scheduler_current_tick);
  writeScalarU32Attribute(
      group.get(), "gas_cell_scheduler_max_bin", diagnostics.gas_cell_scheduler_max_bin);
  writeScalarU64Attribute(
      group.get(), "gas_cell_scheduler_element_count", diagnostics.gas_cell_scheduler_element_count);
  writeScalarU64Attribute(
      group.get(), "gas_cell_scheduler_active_count", diagnostics.gas_cell_scheduler_active_count);
  writeScalarU64Attribute(
      group.get(),
      "gas_cell_scheduler_pending_transition_count",
      diagnostics.gas_cell_scheduler_pending_transition_count);
  writeScalarU64Attribute(group.get(), "pm_cadence_steps", diagnostics.pm_cadence_steps);
  writeScalarU64Attribute(group.get(), "pm_gravity_kick_opportunity", diagnostics.pm_gravity_kick_opportunity);
  writeScalarU64Attribute(group.get(), "pm_field_version", diagnostics.pm_field_version);
  writeScalarU64Attribute(group.get(), "pm_last_refresh_opportunity", diagnostics.pm_last_refresh_opportunity);
  writeScalarU64Attribute(group.get(), "pm_last_refresh_step_index", diagnostics.pm_last_refresh_step_index);
  writeScalarU32Attribute(group.get(), "pm_refresh_commit_pending", diagnostics.pm_refresh_commit_pending ? 1U : 0U);
  writeScalarU32Attribute(group.get(), "pm_long_range_field_valid", diagnostics.pm_long_range_field_valid ? 1U : 0U);
  writeScalarU32Attribute(group.get(), "output_enabled", diagnostics.output_enabled ? 1U : 0U);
  writeScalarU32Attribute(group.get(), "output_snapshot_due", diagnostics.output_snapshot_due ? 1U : 0U);
  writeScalarU32Attribute(group.get(), "output_checkpoint_due", diagnostics.output_checkpoint_due ? 1U : 0U);
  writeScalarU64Attribute(group.get(), "output_last_completed_step_index", diagnostics.output_last_completed_step_index);
  writeScalarU64Attribute(group.get(), "output_next_snapshot_step_index", diagnostics.output_next_snapshot_step_index);
  writeScalarU64Attribute(group.get(), "stochastic_module_count", diagnostics.stochastic_module_count);
}

[[nodiscard]] RestartDiagnosticsSummary readRestartDiagnosticsGroup(
    hid_t root,
    std::uint32_t schema_version) {
  Hdf5Handle group(H5Gopen2(root, "/restart_diagnostics", H5P_DEFAULT));
  if (!group.valid()) {
    throw std::runtime_error("restart schema validation missing required group: /restart_diagnostics");
  }
  RestartDiagnosticsSummary diagnostics;
  diagnostics.restart_schema_name = readScalarStringAttribute(group.get(), "restart_schema_name");
  diagnostics.restart_schema_version = readScalarU32Attribute(group.get(), "restart_schema_version");
  diagnostics.current_boundary_kind = readScalarStringAttribute(group.get(), "current_boundary_kind");
  diagnostics.last_completed_boundary_kind = readScalarStringAttribute(group.get(), "last_completed_boundary_kind");
  diagnostics.restart_safe = readScalarU32Attribute(group.get(), "restart_safe") != 0U;
  diagnostics.step_index = readScalarU64Attribute(group.get(), "step_index");
  diagnostics.scheduler_current_tick = readScalarU64Attribute(group.get(), "scheduler_current_tick");
  diagnostics.scheduler_max_bin = readScalarU32Attribute(group.get(), "scheduler_max_bin");
  diagnostics.scheduler_element_count = readScalarU64Attribute(group.get(), "scheduler_element_count");
  diagnostics.scheduler_active_count = readScalarU64Attribute(group.get(), "scheduler_active_count");
  diagnostics.scheduler_pending_transition_count =
      readScalarU64Attribute(group.get(), "scheduler_pending_transition_count");
  if (schema_version >= k_restart_schema_v19) {
    diagnostics.gas_cell_scheduler_current_tick =
        readScalarU64Attribute(group.get(), "gas_cell_scheduler_current_tick");
    diagnostics.gas_cell_scheduler_max_bin =
        readScalarU32Attribute(group.get(), "gas_cell_scheduler_max_bin");
    diagnostics.gas_cell_scheduler_element_count =
        readScalarU64Attribute(group.get(), "gas_cell_scheduler_element_count");
    diagnostics.gas_cell_scheduler_active_count =
        readScalarU64Attribute(group.get(), "gas_cell_scheduler_active_count");
    diagnostics.gas_cell_scheduler_pending_transition_count =
        readScalarU64Attribute(group.get(), "gas_cell_scheduler_pending_transition_count");
  }
  diagnostics.pm_cadence_steps = readScalarU64Attribute(group.get(), "pm_cadence_steps");
  diagnostics.pm_gravity_kick_opportunity = readScalarU64Attribute(group.get(), "pm_gravity_kick_opportunity");
  diagnostics.pm_field_version = readScalarU64Attribute(group.get(), "pm_field_version");
  diagnostics.pm_last_refresh_opportunity = readScalarU64Attribute(group.get(), "pm_last_refresh_opportunity");
  diagnostics.pm_last_refresh_step_index = readScalarU64Attribute(group.get(), "pm_last_refresh_step_index");
  diagnostics.pm_refresh_commit_pending = readScalarU32Attribute(group.get(), "pm_refresh_commit_pending") != 0U;
  diagnostics.pm_long_range_field_valid = readScalarU32Attribute(group.get(), "pm_long_range_field_valid") != 0U;
  diagnostics.output_enabled = readScalarU32Attribute(group.get(), "output_enabled") != 0U;
  diagnostics.output_snapshot_due = readScalarU32Attribute(group.get(), "output_snapshot_due") != 0U;
  diagnostics.output_checkpoint_due = readScalarU32Attribute(group.get(), "output_checkpoint_due") != 0U;
  diagnostics.output_last_completed_step_index = readScalarU64Attribute(group.get(), "output_last_completed_step_index");
  diagnostics.output_next_snapshot_step_index = readScalarU64Attribute(group.get(), "output_next_snapshot_step_index");
  diagnostics.stochastic_module_count = readScalarU64Attribute(group.get(), "stochastic_module_count");
  return diagnostics;
}
#endif

}  // namespace

const RestartSchema& restartSchema() {
  static const RestartSchema schema{};
  return schema;
}

bool isRestartSchemaCompatible(std::uint32_t file_schema_version) {
  return file_schema_version == restartSchema().version ||
      file_schema_version == k_restart_schema_v19 ||
      file_schema_version == k_restart_schema_v18 ||
      file_schema_version == k_restart_schema_v17 || file_schema_version == k_restart_schema_v16 ||
      file_schema_version == k_restart_schema_v15 || file_schema_version == k_restart_schema_v14;
}

const std::vector<std::string_view>& exactRestartCompletenessChecklist() {
  static const std::vector<std::string_view> checklist = {
      "simulation_state_lanes_and_metadata",
      "particle_identity_softening_and_drift_epoch_lanes",
      "gas_cell_identity_lanes",
      "hydro_geometry_patch_state",
      "amr_pending_flux_register_state",
      "amr_temporal_boundary_history_state",
      "species_specific_sidecars",
      "module_sidecars_with_schema_versions",
      "integrator_state",
      "integrator_owned_pm_sync_state",
      "gravity_force_cache_at_kdk_boundary",
      "scheduler_persistent_state",
      "gas_cell_scheduler_persistent_state_keyed_by_gas_cell_id",
      "output_cadence_persistent_state",
      "stochastic_module_persistent_state",
      "restart_diagnostics_summary",
      "distributed_gravity_state",
      "normalized_config_text_and_hash",
      "provenance_record",
      "payload_integrity_hash_and_hex"};
  return checklist;
}

std::uint64_t restartPayloadIntegrityHashImpl(
    const RestartWritePayload& payload,
    bool include_gas_identity_records,
    bool include_pending_flux_registers,
    bool include_temporal_boundary_history,
    bool include_gas_cell_scheduler,
    bool include_gravity_force_cache) {
  if (payload.persistent_state.simulation_state == nullptr || payload.integrator_state == nullptr || payload.scheduler == nullptr) {
    throw std::invalid_argument("restart payload must provide state, integrator_state, and scheduler");
  }
  validateContinuationMetadata(
      payload.normalized_config_text,
      payload.normalized_config_hash_hex,
      payload.provenance,
      "restart payload");
  validateHydroGeometryStateForRestart(*payload.persistent_state.simulation_state, "restart payload");
  if (!payload.persistent_state.simulation_state->validateOwnershipInvariants()) {
    throw std::invalid_argument("restart payload state failed ownership invariant validation");
  }
  const GravityForceCachePersistentState empty_force_cache{};
  const GravityForceCachePersistentState& gravity_force_cache =
      payload.gravity_force_cache != nullptr ? *payload.gravity_force_cache : empty_force_cache;
  if (include_gravity_force_cache) {
    validateGravityForceCacheForRestart(
        gravity_force_cache, *payload.persistent_state.simulation_state, "restart payload");
  }
  const core::TimeBinPersistentState scheduler_state_for_validation = payload.scheduler->exportPersistentState();
  const core::TimeBinPersistentState gas_cell_scheduler_state_for_validation =
      payload.gas_cell_scheduler != nullptr
          ? payload.gas_cell_scheduler->exportPersistentState()
          : legacyGasCellSchedulerStateFromMirrors(
                *payload.persistent_state.simulation_state, scheduler_state_for_validation);
  const std::vector<std::uint64_t> gas_cell_scheduler_ids =
      gasCellIdsByDenseLocalRow(*payload.persistent_state.simulation_state, "restart payload");
  core::assertCanWriteCheckpointAtBoundary(*payload.integrator_state, scheduler_state_for_validation.current_tick);
  if (gas_cell_scheduler_state_for_validation.current_tick != scheduler_state_for_validation.current_tick) {
    throw std::invalid_argument("restart payload particle and gas-cell schedulers must share current_tick");
  }
  validateRestartTimeBinMirrorsAgainstScheduler(
      *payload.persistent_state.simulation_state,
      scheduler_state_for_validation,
      "restart payload");
  validateGasCellTimeBinMirrorsAgainstScheduler(
      *payload.persistent_state.simulation_state,
      gas_cell_scheduler_state_for_validation,
      gas_cell_scheduler_ids,
      "restart payload");
  if (payload.distributed_gravity_state.world_size <= 0) {
    throw std::invalid_argument("restart payload distributed_gravity_state.world_size must be positive");
  }
  if (payload.distributed_gravity_state.owning_rank_by_item.size() != payload.persistent_state.simulation_state->particles.size()) {
    throw std::invalid_argument(
        "restart payload distributed_gravity_state.owning_rank_by_item must match particle count");
  }
  if (payload.distributed_gravity_state.pm_slab_begin_x_by_rank.size() !=
          static_cast<std::size_t>(payload.distributed_gravity_state.world_size) ||
      payload.distributed_gravity_state.pm_slab_end_x_by_rank.size() !=
          static_cast<std::size_t>(payload.distributed_gravity_state.world_size)) {
    throw std::invalid_argument(
        "restart payload distributed_gravity_state PM slab ownership must match world_size");
  }
  if (payload.distributed_gravity_state.long_range_restart_policy != "deterministic_rebuild") {
    throw std::invalid_argument(
        "restart payload distributed_gravity_state.long_range_restart_policy must be deterministic_rebuild");
  }
  validateOutputCadenceStateForRestart(
      payload.output_cadence_state, *payload.integrator_state, "restart payload");
  validateStochasticStateForRestart(
      payload.stochastic_state, *payload.integrator_state, "restart payload");

  std::uint64_t hash = k_offset_basis;

  const auto append_u64 = [&hash](std::uint64_t value) {
    const std::array<std::byte, sizeof(std::uint64_t)> bytes =
        std::bit_cast<std::array<std::byte, sizeof(std::uint64_t)>>(value);
    hash = fnv1aAppend(hash, bytes);
  };

  const auto append_string = [&hash, &append_u64](const std::string& value) {
    append_u64(static_cast<std::uint64_t>(value.size()));
    hash = fnv1aAppend(hash, {reinterpret_cast<const std::byte*>(value.data()), value.size()});
  };

  append_string(payload.normalized_config_text);
  append_string(payload.normalized_config_hash_hex);
  append_string(core::serializeProvenanceRecord(payload.provenance));
  append_string(payload.persistent_state.simulation_state->metadata.serialize());

  const auto append_any_vec = [&hash, &append_u64](const auto& values) {
    const auto bytes = asBytesSpan(values);
    append_u64(static_cast<std::uint64_t>(values.size()));
    append_u64(static_cast<std::uint64_t>(bytes.size()));
    hash = fnv1aAppend(hash, bytes);
  };

  const core::SimulationState& state = *payload.persistent_state.simulation_state;
  append_any_vec(state.particles.position_x_comoving);
  append_any_vec(state.particles.position_y_comoving);
  append_any_vec(state.particles.position_z_comoving);
  append_any_vec(state.particles.velocity_x_peculiar);
  append_any_vec(state.particles.velocity_y_peculiar);
  append_any_vec(state.particles.velocity_z_peculiar);
  append_any_vec(state.particles.mass_code);
  append_any_vec(state.particles.time_bin);
  append_any_vec(state.particle_sidecar.particle_id);
  append_any_vec(state.particle_sidecar.sfc_key);
  append_any_vec(state.particle_sidecar.species_tag);
  append_any_vec(state.particle_sidecar.particle_flags);
  append_any_vec(state.particle_sidecar.owning_rank);
  append_any_vec(state.particle_sidecar.last_drift_time_code);
  append_any_vec(state.particle_sidecar.last_drift_scale_factor);
  append_any_vec(state.particle_sidecar.gravity_softening_comoving);
  append_any_vec(state.particle_sidecar.has_gravity_softening_override);

  append_any_vec(state.cells.center_x_comoving);
  append_any_vec(state.cells.center_y_comoving);
  append_any_vec(state.cells.center_z_comoving);
  append_any_vec(state.cells.mass_code);
  append_any_vec(state.cells.time_bin);
  append_any_vec(state.cells.patch_index);

  append_any_vec(state.gas_cells.gas_cell_id);
  append_any_vec(state.gas_cells.parent_particle_id);
  append_any_vec(state.gas_cells.velocity_x_peculiar);
  append_any_vec(state.gas_cells.velocity_y_peculiar);
  append_any_vec(state.gas_cells.velocity_z_peculiar);
  append_any_vec(state.gas_cells.density_code);
  append_any_vec(state.gas_cells.pressure_code);
  append_any_vec(state.gas_cells.internal_energy_code);
  append_any_vec(state.gas_cells.temperature_code);
  append_any_vec(state.gas_cells.sound_speed_code);
  if (include_gas_identity_records) {
    append_string(std::string(k_gas_identity_row_policy));
    const auto gas_identity_records = gasIdentityRecordsSortedByLocalRow(state.gas_cell_identity);
    std::vector<std::uint64_t> identity_gas_cell_id;
    std::vector<std::uint8_t> identity_has_parent_particle;
    std::vector<std::uint64_t> identity_parent_particle_id;
    std::vector<std::uint64_t> identity_owning_patch_id;
    std::vector<std::uint32_t> identity_local_cell_row;
    identity_gas_cell_id.reserve(gas_identity_records.size());
    identity_has_parent_particle.reserve(gas_identity_records.size());
    identity_parent_particle_id.reserve(gas_identity_records.size());
    identity_owning_patch_id.reserve(gas_identity_records.size());
    identity_local_cell_row.reserve(gas_identity_records.size());
    for (const core::GasCellIdentityRecord& record : gas_identity_records) {
      identity_gas_cell_id.push_back(record.gas_cell_id);
      identity_has_parent_particle.push_back(record.parent_particle_id.has_value() ? 1U : 0U);
      identity_parent_particle_id.push_back(record.parent_particle_id.value_or(0U));
      identity_owning_patch_id.push_back(record.owning_patch_id);
      identity_local_cell_row.push_back(record.local_cell_row);
    }
    append_any_vec(identity_gas_cell_id);
    append_any_vec(identity_has_parent_particle);
    append_any_vec(identity_parent_particle_id);
    append_any_vec(identity_owning_patch_id);
    append_any_vec(identity_local_cell_row);
  }

  append_any_vec(state.patches.patch_id);
  append_any_vec(state.patches.level);
  append_any_vec(state.patches.first_cell);
  append_any_vec(state.patches.cell_count);
  append_any_vec(state.patches.parent_patch_id);
  append_any_vec(state.patches.morton_key);
  append_any_vec(state.patches.origin_x_comoving);
  append_any_vec(state.patches.origin_y_comoving);
  append_any_vec(state.patches.origin_z_comoving);
  append_any_vec(state.patches.extent_x_comoving);
  append_any_vec(state.patches.extent_y_comoving);
  append_any_vec(state.patches.extent_z_comoving);
  append_any_vec(state.patches.cell_dim_x);
  append_any_vec(state.patches.cell_dim_y);
  append_any_vec(state.patches.cell_dim_z);
  append_any_vec(state.patches.owning_rank);
  if (include_pending_flux_registers) {
    std::vector<core::PendingFluxRegisterRecord> pending_flux_records(
        state.pending_flux_registers.records().begin(), state.pending_flux_registers.records().end());
    std::sort(
        pending_flux_records.begin(),
        pending_flux_records.end(),
        [](const core::PendingFluxRegisterRecord& lhs, const core::PendingFluxRegisterRecord& rhs) {
          return lhs.register_key < rhs.register_key;
        });
    append_u64(static_cast<std::uint64_t>(pending_flux_records.size()));
    for (const core::PendingFluxRegisterRecord& record : pending_flux_records) {
      append_u64(record.register_key);
      append_u64(record.coarse_patch_id);
      append_u64(record.coarse_gas_cell_id);
      append_u64(static_cast<std::uint64_t>(record.coarse_cell_index));
      append_u64(record.level);
      append_u64(record.axis);
      append_u64(record.orientation);
      append_u64(std::bit_cast<std::uint64_t>(record.expected_area_comov));
      append_u64(std::bit_cast<std::uint64_t>(record.coarse_area_accumulated_comov));
      append_u64(std::bit_cast<std::uint64_t>(record.fine_area_accumulated_comov));
      append_u64(std::bit_cast<std::uint64_t>(record.interval_start_code));
      append_u64(std::bit_cast<std::uint64_t>(record.interval_end_code));
      append_u64(std::bit_cast<std::uint64_t>(record.coarse_dt_code));
      append_u64(record.expected_fine_substeps);
      append_u64(record.completed_fine_substeps);
      append_u64(record.fine_substep_coverage_mask);
      append_u64(record.coarse_face_count);
      append_u64(record.fine_face_count);
      append_u64(record.gas_cell_identity_generation);
      append_u64(record.patch_geometry_generation);
      append_u64(std::bit_cast<std::uint64_t>(record.coarse_mass_flux_integral_code));
      append_u64(std::bit_cast<std::uint64_t>(record.coarse_momentum_x_flux_integral_code));
      append_u64(std::bit_cast<std::uint64_t>(record.coarse_momentum_y_flux_integral_code));
      append_u64(std::bit_cast<std::uint64_t>(record.coarse_momentum_z_flux_integral_code));
      append_u64(std::bit_cast<std::uint64_t>(record.coarse_total_energy_flux_integral_code));
      append_u64(std::bit_cast<std::uint64_t>(record.fine_mass_flux_integral_code));
      append_u64(std::bit_cast<std::uint64_t>(record.fine_momentum_x_flux_integral_code));
      append_u64(std::bit_cast<std::uint64_t>(record.fine_momentum_y_flux_integral_code));
      append_u64(std::bit_cast<std::uint64_t>(record.fine_momentum_z_flux_integral_code));
      append_u64(std::bit_cast<std::uint64_t>(record.fine_total_energy_flux_integral_code));
    }
  }
  if (include_temporal_boundary_history) {
    std::vector<core::AmrTemporalBoundaryHistoryRecord> temporal_records(
        state.amr_temporal_boundary_history.records().begin(),
        state.amr_temporal_boundary_history.records().end());
    std::sort(
        temporal_records.begin(), temporal_records.end(),
        [](const core::AmrTemporalBoundaryHistoryRecord& lhs,
           const core::AmrTemporalBoundaryHistoryRecord& rhs) {
          return lhs.patch_id < rhs.patch_id;
        });
    append_u64(static_cast<std::uint64_t>(temporal_records.size()));
    for (const core::AmrTemporalBoundaryHistoryRecord& record : temporal_records) {
      append_u64(record.patch_id);
      append_u64(record.patch_level);
      append_u64(record.patch_geometry_fingerprint);
      append_u64(record.gas_cell_identity_generation);
      append_u64(std::bit_cast<std::uint64_t>(record.interval_start_code));
      append_u64(std::bit_cast<std::uint64_t>(record.interval_end_code));
      append_u64(record.end_state_valid ? 1U : 0U);
      std::vector<core::AmrTemporalBoundaryHistoryCellRecord> cells = record.cells;
      std::sort(cells.begin(), cells.end(), [](const auto& lhs, const auto& rhs) {
        return lhs.patch_local_cell == rhs.patch_local_cell
            ? lhs.gas_cell_id < rhs.gas_cell_id
            : lhs.patch_local_cell < rhs.patch_local_cell;
      });
      append_u64(static_cast<std::uint64_t>(cells.size()));
      for (const core::AmrTemporalBoundaryHistoryCellRecord& cell : cells) {
        append_u64(cell.gas_cell_id);
        append_u64(static_cast<std::uint64_t>(cell.patch_local_cell));
        append_u64(std::bit_cast<std::uint64_t>(cell.start_mass_density_comoving));
        append_u64(std::bit_cast<std::uint64_t>(cell.start_momentum_density_x_comoving));
        append_u64(std::bit_cast<std::uint64_t>(cell.start_momentum_density_y_comoving));
        append_u64(std::bit_cast<std::uint64_t>(cell.start_momentum_density_z_comoving));
        append_u64(std::bit_cast<std::uint64_t>(cell.start_total_energy_density_comoving));
        append_u64(std::bit_cast<std::uint64_t>(cell.end_mass_density_comoving));
        append_u64(std::bit_cast<std::uint64_t>(cell.end_momentum_density_x_comoving));
        append_u64(std::bit_cast<std::uint64_t>(cell.end_momentum_density_y_comoving));
        append_u64(std::bit_cast<std::uint64_t>(cell.end_momentum_density_z_comoving));
        append_u64(std::bit_cast<std::uint64_t>(cell.end_total_energy_density_comoving));
      }
    }
  }
  append_any_vec(state.star_particles.particle_index);
  append_any_vec(state.star_particles.formation_scale_factor);
  append_any_vec(state.star_particles.birth_mass_code);
  append_any_vec(state.star_particles.metallicity_mass_fraction);
  append_any_vec(state.star_particles.stellar_age_years_last);
  append_any_vec(state.star_particles.stellar_returned_mass_cumulative_code);
  append_any_vec(state.star_particles.stellar_returned_metals_cumulative_code);
  append_any_vec(state.star_particles.stellar_feedback_energy_cumulative_erg);
  for (std::size_t channel = 0; channel < state.star_particles.stellar_returned_mass_channel_cumulative_code.size(); ++channel) {
    append_any_vec(state.star_particles.stellar_returned_mass_channel_cumulative_code[channel]);
    append_any_vec(state.star_particles.stellar_returned_metals_channel_cumulative_code[channel]);
    append_any_vec(state.star_particles.stellar_feedback_energy_channel_cumulative_erg[channel]);
  }
  append_any_vec(state.black_holes.particle_index);
  append_any_vec(state.black_holes.host_cell_index);
  append_any_vec(state.black_holes.subgrid_mass_code);
  append_any_vec(state.black_holes.accretion_rate_code);
  append_any_vec(state.black_holes.feedback_energy_code);
  append_any_vec(state.black_holes.eddington_ratio);
  append_any_vec(state.black_holes.cumulative_accreted_mass_code);
  append_any_vec(state.black_holes.cumulative_feedback_energy_code);
  append_any_vec(state.black_holes.duty_cycle_active_time_code);
  append_any_vec(state.black_holes.duty_cycle_total_time_code);
  append_any_vec(state.tracers.particle_index);
  append_any_vec(state.tracers.parent_particle_id);
  append_any_vec(state.tracers.injection_step);
  append_any_vec(state.tracers.host_cell_index);
  append_any_vec(state.tracers.mass_fraction_of_host);
  append_any_vec(state.tracers.last_host_mass_code);
  append_any_vec(state.tracers.cumulative_exchanged_mass_code);
  append_any_vec(state.species.count_by_species);

  const auto ordered_sidecars = state.sidecars.blocksSortedByName();
  for (const core::ModuleSidecarBlock* block : ordered_sidecars) {
    append_string(block->module_name);
    append_u64(static_cast<std::uint64_t>(block->schema_version));
    append_u64(block->particle_indexed ? 1ULL : 0ULL);
    append_u64(static_cast<std::uint64_t>(block->row_stride_bytes));
    append_u64(static_cast<std::uint64_t>(block->required_species_mask));
    append_u64(static_cast<std::uint64_t>(block->requirement.kind));
    append_u64(static_cast<std::uint64_t>(block->requirement.species_mask));
    append_u64(static_cast<std::uint64_t>(block->requirement.particle_flags_mask));
    append_u64(std::bit_cast<std::uint64_t>(block->requirement.threshold_code));
    append_u64(static_cast<std::uint64_t>(block->particle_id_by_row.size()));
    for (const std::uint64_t particle_id : block->particle_id_by_row) {
      append_u64(particle_id);
    }
    append_u64(static_cast<std::uint64_t>(block->payload.size()));
    hash = fnv1aAppend(hash, std::span<const std::byte>(block->payload.data(), block->payload.size()));
  }

  append_u64(payload.integrator_state->step_index);
  append_u64(std::bit_cast<std::uint64_t>(payload.integrator_state->current_time_code));
  append_u64(std::bit_cast<std::uint64_t>(payload.integrator_state->current_scale_factor));
  append_u64(std::bit_cast<std::uint64_t>(payload.integrator_state->current_redshift));
  append_u64(std::bit_cast<std::uint64_t>(payload.integrator_state->current_hubble_rate_code));
  append_u64(std::bit_cast<std::uint64_t>(payload.integrator_state->time_si_per_code));
  append_u64(std::bit_cast<std::uint64_t>(payload.integrator_state->dt_time_code));
  append_u64(std::bit_cast<std::uint64_t>(payload.integrator_state->last_drift_factor_code));
  append_u64(std::bit_cast<std::uint64_t>(payload.integrator_state->last_first_kick_factor_code));
  append_u64(std::bit_cast<std::uint64_t>(payload.integrator_state->last_second_kick_factor_code));
  append_u64(std::bit_cast<std::uint64_t>(payload.integrator_state->last_first_hubble_drag_factor));
  append_u64(std::bit_cast<std::uint64_t>(payload.integrator_state->last_second_hubble_drag_factor));
  append_u64(static_cast<std::uint64_t>(payload.integrator_state->scheme));
  append_u64(static_cast<std::uint64_t>(payload.integrator_state->current_boundary_kind));
  append_u64(static_cast<std::uint64_t>(payload.integrator_state->last_completed_boundary_kind));
  append_u64(payload.integrator_state->inside_kdk_step ? 1ull : 0ull);
  append_u64(payload.integrator_state->last_completed_restart_safe ? 1ull : 0ull);
  append_u64(payload.integrator_state->time_bins.hierarchical_enabled ? 1ull : 0ull);
  append_u64(static_cast<std::uint64_t>(payload.integrator_state->time_bins.active_bin));
  append_u64(static_cast<std::uint64_t>(payload.integrator_state->time_bins.max_bin));
  append_u64(payload.integrator_state->pm_long_range_field_valid ? 1ull : 0ull);
  if (include_gravity_force_cache) {
    append_u64(payload.integrator_state->pm_refresh_enabled ? 1ull : 0ull);
  }
  const core::PmSynchronizationPersistentState pm_sync_state = payload.integrator_state->pm_sync_state.exportPersistentState();
  append_u64(pm_sync_state.cadence_steps);
  append_u64(pm_sync_state.gravity_kick_opportunity);
  append_u64(pm_sync_state.last_refresh_opportunity);
  append_u64(pm_sync_state.field_version);
  append_u64(pm_sync_state.last_refresh_step_index);
  append_u64(std::bit_cast<std::uint64_t>(pm_sync_state.last_refresh_scale_factor));
  append_u64(pm_sync_state.refresh_commit_pending ? 1ull : 0ull);
  append_u64(pm_sync_state.pending_refresh_opportunity);
  append_u64(pm_sync_state.pending_refresh_field_version);
  if (include_gravity_force_cache) {
    append_u64(gravity_force_cache.valid ? 1ULL : 0ULL);
    append_any_vec(gravity_force_cache.particle_id);
    append_any_vec(gravity_force_cache.gas_cell_id);
    append_any_vec(gravity_force_cache.particle_accel_x_comoving);
    append_any_vec(gravity_force_cache.particle_accel_y_comoving);
    append_any_vec(gravity_force_cache.particle_accel_z_comoving);
    append_any_vec(gravity_force_cache.cell_accel_x_comoving);
    append_any_vec(gravity_force_cache.cell_accel_y_comoving);
    append_any_vec(gravity_force_cache.cell_accel_z_comoving);
  }

  const core::TimeBinPersistentState& scheduler_state = scheduler_state_for_validation;
  append_u64(scheduler_state.current_tick);
  append_u64(static_cast<std::uint64_t>(scheduler_state.max_bin));
  append_any_vec(scheduler_state.bin_index);
  append_any_vec(scheduler_state.next_activation_tick);
  append_any_vec(scheduler_state.active_flag);
  append_any_vec(scheduler_state.pending_bin_index);
  if (include_gas_cell_scheduler) {
    append_string(std::string(k_gas_cell_scheduler_identity_key));
    const core::TimeBinPersistentState& gas_scheduler_state = gas_cell_scheduler_state_for_validation;
    append_u64(gas_scheduler_state.current_tick);
    append_u64(static_cast<std::uint64_t>(gas_scheduler_state.max_bin));
    append_any_vec(gas_cell_scheduler_ids);
    append_any_vec(gas_scheduler_state.bin_index);
    append_any_vec(gas_scheduler_state.next_activation_tick);
    append_any_vec(gas_scheduler_state.active_flag);
    append_any_vec(gas_scheduler_state.pending_bin_index);
  }
  append_u64(payload.output_cadence_state.output_enabled ? 1ull : 0ull);
  append_u64(payload.output_cadence_state.write_restarts ? 1ull : 0ull);
  append_u64(payload.output_cadence_state.snapshot_due ? 1ull : 0ull);
  append_u64(payload.output_cadence_state.checkpoint_due ? 1ull : 0ull);
  append_u64(payload.output_cadence_state.last_completed_step_index);
  append_u64(payload.output_cadence_state.snapshot_interval_steps);
  append_u64(payload.output_cadence_state.next_snapshot_step_index);
  append_string(payload.output_cadence_state.snapshot_stem);
  append_string(payload.output_cadence_state.restart_stem);
  for (const StochasticModulePersistentState& module_state : sortedStochasticModules(payload.stochastic_state)) {
    append_string(module_state.module_name);
    append_u64(module_state.schema_version);
    append_string(module_state.rng_policy);
    append_u64(module_state.random_seed);
    append_u64(module_state.rank_local_seed_offset);
    append_u64(module_state.last_committed_step_index);
    append_u64(module_state.deterministic_from_serialized_inputs ? 1ull : 0ull);
  }
  append_string(payload.distributed_gravity_state.serialize());

  return hash;
}

std::uint64_t restartPayloadIntegrityHash(const RestartWritePayload& payload) {
  return restartPayloadIntegrityHashImpl(payload, true, true, true, true, true);
}

std::string restartPayloadIntegrityHashHex(const RestartWritePayload& payload) {
  return hexU64(restartPayloadIntegrityHash(payload));
}

void writeRestartCheckpointHdf5(
    const std::filesystem::path& output_path,
    const RestartWritePayload& payload,
    const RestartWritePolicy& policy) {
#if !COSMOSIM_ENABLE_HDF5
  (void)output_path;
  (void)payload;
  (void)policy;
  throw std::runtime_error("restart checkpoint requires COSMOSIM_ENABLE_HDF5=ON");
#else
  if (payload.persistent_state.simulation_state == nullptr || payload.integrator_state == nullptr || payload.scheduler == nullptr) {
    throw std::invalid_argument("restart write payload must include state, integrator_state, and scheduler");
  }
  validateContinuationMetadata(
      payload.normalized_config_text,
      payload.normalized_config_hash_hex,
      payload.provenance,
      "restart writer");
  if (!payload.persistent_state.simulation_state->validateOwnershipInvariants()) {
    throw std::invalid_argument("cannot checkpoint invalid simulation state");
  }
  const core::TimeBinPersistentState scheduler_state_for_validation = payload.scheduler->exportPersistentState();
  const core::TimeBinPersistentState gas_cell_scheduler_state_for_validation =
      payload.gas_cell_scheduler != nullptr
          ? payload.gas_cell_scheduler->exportPersistentState()
          : legacyGasCellSchedulerStateFromMirrors(
                *payload.persistent_state.simulation_state, scheduler_state_for_validation);
  const std::vector<std::uint64_t> gas_cell_scheduler_ids =
      gasCellIdsByDenseLocalRow(*payload.persistent_state.simulation_state, "restart writer");
  core::assertCanWriteCheckpointAtBoundary(*payload.integrator_state, scheduler_state_for_validation.current_tick);
  if (gas_cell_scheduler_state_for_validation.current_tick != scheduler_state_for_validation.current_tick) {
    throw std::invalid_argument("restart writer particle and gas-cell schedulers must share current_tick");
  }
  validateRestartTimeBinMirrorsAgainstScheduler(
      *payload.persistent_state.simulation_state,
      scheduler_state_for_validation,
      "restart writer");
  validateGasCellTimeBinMirrorsAgainstScheduler(
      *payload.persistent_state.simulation_state,
      gas_cell_scheduler_state_for_validation,
      gas_cell_scheduler_ids,
      "restart writer");
  validateOutputCadenceStateForRestart(
      payload.output_cadence_state, *payload.integrator_state, "restart writer");
  validateStochasticStateForRestart(
      payload.stochastic_state, *payload.integrator_state, "restart writer");
  const GravityForceCachePersistentState empty_force_cache{};
  const GravityForceCachePersistentState& gravity_force_cache =
      payload.gravity_force_cache != nullptr ? *payload.gravity_force_cache : empty_force_cache;
  validateGravityForceCacheForRestart(
      gravity_force_cache, *payload.persistent_state.simulation_state, "restart writer");

  std::filesystem::create_directories(output_path.parent_path());
  const std::filesystem::path temporary_path = output_path.string() + policy.temporary_suffix;

  Hdf5Handle file(H5Fcreate(temporary_path.string().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT));
  if (!file.valid()) {
    throw std::runtime_error("failed to create temporary restart file: " + temporary_path.string());
  }

  const auto& shared_names = sharedIoContractNames();
  writeScalarStringAttribute(
      file.get(), shared_names.file_kind_attribute, std::string(shared_names.restart_checkpoint_file_kind));
  writeScalarStringAttribute(file.get(), "restart_schema_name", restartSchema().name);
  writeScalarU32Attribute(file.get(), "restart_schema_version", restartSchema().version);
  writeScalarStringAttribute(file.get(), "normalized_config_hash_hex", payload.normalized_config_hash_hex);
  writeScalarStringAttribute(file.get(), "payload_integrity_hash_hex", restartPayloadIntegrityHashHex(payload));
  writeScalarU64Attribute(file.get(), "payload_integrity_hash", restartPayloadIntegrityHash(payload));
  writeStringDataset(
      file.get(),
      std::string(shared_names.normalized_config_text_dataset),
      payload.normalized_config_text);
  writeStringDataset(
      file.get(),
      std::string(shared_names.provenance_record_dataset),
      core::serializeProvenanceRecord(payload.provenance));

  writeStateGroup(file.get(), *payload.persistent_state.simulation_state);

  Hdf5Handle integrator_group(openOrCreateGroup(file.get(), "/integrator"));
  writeScalarF64Attribute(integrator_group.get(), "current_time_code", payload.integrator_state->current_time_code);
  writeScalarF64Attribute(integrator_group.get(), "current_scale_factor", payload.integrator_state->current_scale_factor);
  writeScalarF64Attribute(integrator_group.get(), "current_redshift", payload.integrator_state->current_redshift);
  writeScalarF64Attribute(integrator_group.get(), "current_hubble_rate_code", payload.integrator_state->current_hubble_rate_code);
  writeScalarF64Attribute(integrator_group.get(), "time_si_per_code", payload.integrator_state->time_si_per_code);
  writeScalarF64Attribute(integrator_group.get(), "dt_time_code", payload.integrator_state->dt_time_code);
  writeScalarF64Attribute(integrator_group.get(), "last_drift_factor_code", payload.integrator_state->last_drift_factor_code);
  writeScalarF64Attribute(integrator_group.get(), "last_first_kick_factor_code", payload.integrator_state->last_first_kick_factor_code);
  writeScalarF64Attribute(integrator_group.get(), "last_second_kick_factor_code", payload.integrator_state->last_second_kick_factor_code);
  writeScalarF64Attribute(integrator_group.get(), "last_first_hubble_drag_factor", payload.integrator_state->last_first_hubble_drag_factor);
  writeScalarF64Attribute(integrator_group.get(), "last_second_hubble_drag_factor", payload.integrator_state->last_second_hubble_drag_factor);
  writeScalarU64Attribute(integrator_group.get(), "step_index", payload.integrator_state->step_index);
  writeScalarU32Attribute(integrator_group.get(), "scheme", static_cast<std::uint32_t>(payload.integrator_state->scheme));
  writeScalarU32Attribute(integrator_group.get(), "current_boundary_kind", static_cast<std::uint32_t>(payload.integrator_state->current_boundary_kind));
  writeScalarU32Attribute(integrator_group.get(), "last_completed_boundary_kind", static_cast<std::uint32_t>(payload.integrator_state->last_completed_boundary_kind));
  writeScalarU32Attribute(integrator_group.get(), "inside_kdk_step", payload.integrator_state->inside_kdk_step ? 1U : 0U);
  writeScalarU32Attribute(integrator_group.get(), "last_completed_restart_safe", payload.integrator_state->last_completed_restart_safe ? 1U : 0U);
  writeScalarU32Attribute(integrator_group.get(), "time_bins_hierarchical", payload.integrator_state->time_bins.hierarchical_enabled ? 1U : 0U);
  writeScalarU32Attribute(integrator_group.get(), "time_bins_active_bin", payload.integrator_state->time_bins.active_bin);
  writeScalarU32Attribute(integrator_group.get(), "time_bins_max_bin", payload.integrator_state->time_bins.max_bin);
  writeScalarU32Attribute(integrator_group.get(), "pm_long_range_field_valid", payload.integrator_state->pm_long_range_field_valid ? 1U : 0U);
  writeScalarU32Attribute(integrator_group.get(), "pm_refresh_enabled", payload.integrator_state->pm_refresh_enabled ? 1U : 0U);
  const core::PmSynchronizationPersistentState pm_sync_state = payload.integrator_state->pm_sync_state.exportPersistentState();
  writeScalarU64Attribute(integrator_group.get(), "pm_cadence_steps", pm_sync_state.cadence_steps);
  writeScalarU64Attribute(integrator_group.get(), "pm_gravity_kick_opportunity", pm_sync_state.gravity_kick_opportunity);
  writeScalarU64Attribute(integrator_group.get(), "pm_last_refresh_opportunity", pm_sync_state.last_refresh_opportunity);
  writeScalarU64Attribute(integrator_group.get(), "pm_field_version", pm_sync_state.field_version);
  writeScalarU64Attribute(integrator_group.get(), "pm_last_refresh_step_index", pm_sync_state.last_refresh_step_index);
  writeScalarF64Attribute(integrator_group.get(), "pm_last_refresh_scale_factor", pm_sync_state.last_refresh_scale_factor);
  writeScalarU32Attribute(integrator_group.get(), "pm_refresh_commit_pending", pm_sync_state.refresh_commit_pending ? 1U : 0U);
  writeScalarU64Attribute(integrator_group.get(), "pm_pending_refresh_opportunity", pm_sync_state.pending_refresh_opportunity);
  writeScalarU64Attribute(integrator_group.get(), "pm_pending_refresh_field_version", pm_sync_state.pending_refresh_field_version);

  Hdf5Handle force_cache_group(openOrCreateGroup(file.get(), "/gravity_force_cache"));
  writeScalarU32Attribute(force_cache_group.get(), "valid", gravity_force_cache.valid ? 1U : 0U);
  writeDataset1d(force_cache_group.get(), "particle_id", H5T_STD_U64LE, H5T_NATIVE_UINT64, gravity_force_cache.particle_id);
  writeDataset1d(force_cache_group.get(), "gas_cell_id", H5T_STD_U64LE, H5T_NATIVE_UINT64, gravity_force_cache.gas_cell_id);
  writeDataset1d(force_cache_group.get(), "particle_accel_x_comoving", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, gravity_force_cache.particle_accel_x_comoving);
  writeDataset1d(force_cache_group.get(), "particle_accel_y_comoving", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, gravity_force_cache.particle_accel_y_comoving);
  writeDataset1d(force_cache_group.get(), "particle_accel_z_comoving", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, gravity_force_cache.particle_accel_z_comoving);
  writeDataset1d(force_cache_group.get(), "cell_accel_x_comoving", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, gravity_force_cache.cell_accel_x_comoving);
  writeDataset1d(force_cache_group.get(), "cell_accel_y_comoving", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, gravity_force_cache.cell_accel_y_comoving);
  writeDataset1d(force_cache_group.get(), "cell_accel_z_comoving", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, gravity_force_cache.cell_accel_z_comoving);

  const core::TimeBinPersistentState scheduler_state = scheduler_state_for_validation;
  Hdf5Handle scheduler_group(openOrCreateGroup(file.get(), "/scheduler"));
  writeScalarU64Attribute(scheduler_group.get(), "current_tick", scheduler_state.current_tick);
  writeScalarU32Attribute(scheduler_group.get(), "max_bin", scheduler_state.max_bin);
  writeDataset1d(scheduler_group.get(), "bin_index", H5T_STD_U8LE, H5T_NATIVE_UINT8, scheduler_state.bin_index);
  writeDataset1d(
      scheduler_group.get(),
      "next_activation_tick",
      H5T_STD_U64LE,
      H5T_NATIVE_UINT64,
      scheduler_state.next_activation_tick);
  writeDataset1d(scheduler_group.get(), "active_flag", H5T_STD_U8LE, H5T_NATIVE_UINT8, scheduler_state.active_flag);
  writeDataset1d(
      scheduler_group.get(),
      "pending_bin_index",
      H5T_STD_U8LE,
      H5T_NATIVE_UINT8,
      scheduler_state.pending_bin_index);

  Hdf5Handle gas_scheduler_group(openOrCreateGroup(file.get(), "/gas_cell_scheduler"));
  writeScalarStringAttribute(
      gas_scheduler_group.get(), "identity_key", std::string(k_gas_cell_scheduler_identity_key));
  writeScalarU64Attribute(
      gas_scheduler_group.get(), "current_tick", gas_cell_scheduler_state_for_validation.current_tick);
  writeScalarU32Attribute(
      gas_scheduler_group.get(), "max_bin", gas_cell_scheduler_state_for_validation.max_bin);
  writeDataset1d(
      gas_scheduler_group.get(), "gas_cell_id", H5T_STD_U64LE, H5T_NATIVE_UINT64, gas_cell_scheduler_ids);
  writeDataset1d(
      gas_scheduler_group.get(),
      "bin_index",
      H5T_STD_U8LE,
      H5T_NATIVE_UINT8,
      gas_cell_scheduler_state_for_validation.bin_index);
  writeDataset1d(
      gas_scheduler_group.get(),
      "next_activation_tick",
      H5T_STD_U64LE,
      H5T_NATIVE_UINT64,
      gas_cell_scheduler_state_for_validation.next_activation_tick);
  writeDataset1d(
      gas_scheduler_group.get(),
      "active_flag",
      H5T_STD_U8LE,
      H5T_NATIVE_UINT8,
      gas_cell_scheduler_state_for_validation.active_flag);
  writeDataset1d(
      gas_scheduler_group.get(),
      "pending_bin_index",
      H5T_STD_U8LE,
      H5T_NATIVE_UINT8,
      gas_cell_scheduler_state_for_validation.pending_bin_index);
  writeOutputCadenceGroup(file.get(), payload.output_cadence_state);
  writeStochasticStateGroup(file.get(), payload.stochastic_state);
  writeRestartDiagnosticsGroup(
      file.get(),
      makeRestartDiagnosticsSummary(
          *payload.integrator_state,
          scheduler_state,
          gas_cell_scheduler_state_for_validation,
          payload.output_cadence_state,
          payload.stochastic_state));
  writeDistributedGravityGroup(file.get(), payload.distributed_gravity_state);

  if (H5Fflush(file.get(), H5F_SCOPE_GLOBAL) < 0) {
    throw std::runtime_error("failed to flush temporary restart file");
  }
  file = Hdf5Handle();

  maybeFsync(temporary_path, policy.enable_fsync_finalize);

  std::error_code rename_error;
  std::filesystem::rename(temporary_path, output_path, rename_error);
  if (rename_error) {
    throw std::runtime_error(
        "failed to atomically finalize restart checkpoint from '" + temporary_path.string() + "' to '" +
        output_path.string() + "': " + rename_error.message());
  }
#endif
}

RestartReadResult readRestartCheckpointHdf5(const std::filesystem::path& input_path) {
#if !COSMOSIM_ENABLE_HDF5
  (void)input_path;
  throw std::runtime_error("restart checkpoint requires COSMOSIM_ENABLE_HDF5=ON");
#else
  Hdf5Handle file(H5Fopen(input_path.string().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
  if (!file.valid()) {
    throw std::runtime_error("failed to open restart checkpoint: " + input_path.string());
  }

  validateFileKindAttribute(
      file.get(),
      sharedIoContractNames().restart_checkpoint_file_kind,
      "restart reader");
  const std::string schema_name = readScalarStringAttribute(file.get(), "restart_schema_name");
  const std::uint32_t schema_version = readScalarU32Attribute(file.get(), "restart_schema_version");
  const bool current_schema =
      schema_name == restartSchema().name && schema_version == restartSchema().version;
  const bool legacy_v14_schema =
      schema_name == k_restart_schema_name_v14 && schema_version == k_restart_schema_v14;
  const bool legacy_v15_schema =
      schema_name == k_restart_schema_name_v15 && schema_version == k_restart_schema_v15;
  const bool legacy_v16_schema =
      schema_name == k_restart_schema_name_v16 && schema_version == k_restart_schema_v16;
  const bool legacy_v17_schema =
      schema_name == k_restart_schema_name_v17 && schema_version == k_restart_schema_v17;
  const bool legacy_v18_schema =
      schema_name == k_restart_schema_name_v18 && schema_version == k_restart_schema_v18;
  const bool legacy_v19_schema =
      schema_name == k_restart_schema_name_v19 && schema_version == k_restart_schema_v19;
  if ((!current_schema && !legacy_v14_schema && !legacy_v15_schema && !legacy_v16_schema && !legacy_v17_schema && !legacy_v18_schema && !legacy_v19_schema) ||
      !isRestartSchemaCompatible(schema_version)) {
    throw std::runtime_error(
        "restart schema is not compatible: file='" + schema_name + "' v" + std::to_string(schema_version) +
        ", expected='" + restartSchema().name + "' v" + std::to_string(restartSchema().version));
  }
  validateRestartCheckpointSchema(file.get(), schema_version);

  RestartReadResult result;
  result.normalized_config_hash_hex = readScalarStringAttribute(file.get(), "normalized_config_hash_hex");
  result.payload_hash_hex = readScalarStringAttribute(file.get(), "payload_integrity_hash_hex");
  result.payload_hash = readScalarU64Attribute(file.get(), "payload_integrity_hash");

  const auto& shared_names = sharedIoContractNames();
  result.normalized_config_text =
      readStringDataset(file.get(), std::string(shared_names.normalized_config_text_dataset));
  result.provenance = core::deserializeProvenanceRecord(
      readStringDataset(file.get(), std::string(shared_names.provenance_record_dataset)));

  readStateGroup(file.get(), result.state, schema_version);
  validateHydroGeometryStateForRestart(result.state, "restart reader");
  validateAmrTemporalBoundaryHistoryForRestart(result.state);

  Hdf5Handle integrator_group(H5Gopen2(file.get(), "/integrator", H5P_DEFAULT));
  result.integrator_state.current_time_code = readScalarF64Attribute(integrator_group.get(), "current_time_code");
  result.integrator_state.current_scale_factor = readScalarF64Attribute(integrator_group.get(), "current_scale_factor");
  result.integrator_state.current_redshift = readScalarF64Attribute(integrator_group.get(), "current_redshift");
  result.integrator_state.current_hubble_rate_code = readScalarF64Attribute(integrator_group.get(), "current_hubble_rate_code");
  result.integrator_state.time_si_per_code = readScalarF64Attribute(integrator_group.get(), "time_si_per_code");
  result.integrator_state.dt_time_code = readScalarF64Attribute(integrator_group.get(), "dt_time_code");
  result.integrator_state.last_drift_factor_code = readScalarF64Attribute(integrator_group.get(), "last_drift_factor_code");
  result.integrator_state.last_first_kick_factor_code = readScalarF64Attribute(integrator_group.get(), "last_first_kick_factor_code");
  result.integrator_state.last_second_kick_factor_code = readScalarF64Attribute(integrator_group.get(), "last_second_kick_factor_code");
  result.integrator_state.last_first_hubble_drag_factor = readScalarF64Attribute(integrator_group.get(), "last_first_hubble_drag_factor");
  result.integrator_state.last_second_hubble_drag_factor = readScalarF64Attribute(integrator_group.get(), "last_second_hubble_drag_factor");
  result.integrator_state.step_index = readScalarU64Attribute(integrator_group.get(), "step_index");
  result.integrator_state.scheme =
      static_cast<core::TimeStepScheme>(readScalarU32Attribute(integrator_group.get(), "scheme"));
  result.integrator_state.current_boundary_kind =
      static_cast<core::StepBoundaryKind>(readScalarU32Attribute(integrator_group.get(), "current_boundary_kind"));
  result.integrator_state.last_completed_boundary_kind =
      static_cast<core::StepBoundaryKind>(readScalarU32Attribute(integrator_group.get(), "last_completed_boundary_kind"));
  result.integrator_state.inside_kdk_step = readScalarU32Attribute(integrator_group.get(), "inside_kdk_step") != 0U;
  result.integrator_state.last_completed_restart_safe =
      readScalarU32Attribute(integrator_group.get(), "last_completed_restart_safe") != 0U;
  result.integrator_state.time_bins.hierarchical_enabled =
      readScalarU32Attribute(integrator_group.get(), "time_bins_hierarchical") != 0;
  result.integrator_state.time_bins.active_bin =
      static_cast<std::uint8_t>(readScalarU32Attribute(integrator_group.get(), "time_bins_active_bin"));
  result.integrator_state.time_bins.max_bin =
      static_cast<std::uint8_t>(readScalarU32Attribute(integrator_group.get(), "time_bins_max_bin"));
  result.integrator_state.pm_long_range_field_valid =
      readScalarU32Attribute(integrator_group.get(), "pm_long_range_field_valid") != 0U;
  result.integrator_state.pm_refresh_enabled =
      schema_version >= k_restart_schema_v20 &&
      readScalarU32Attribute(integrator_group.get(), "pm_refresh_enabled") != 0U;
  core::PmSynchronizationPersistentState pm_sync_state;
  pm_sync_state.cadence_steps = readScalarU64Attribute(integrator_group.get(), "pm_cadence_steps");
  pm_sync_state.gravity_kick_opportunity = readScalarU64Attribute(integrator_group.get(), "pm_gravity_kick_opportunity");
  pm_sync_state.last_refresh_opportunity = readScalarU64Attribute(integrator_group.get(), "pm_last_refresh_opportunity");
  pm_sync_state.field_version = readScalarU64Attribute(integrator_group.get(), "pm_field_version");
  pm_sync_state.last_refresh_step_index = readScalarU64Attribute(integrator_group.get(), "pm_last_refresh_step_index");
  pm_sync_state.last_refresh_scale_factor = readScalarF64Attribute(integrator_group.get(), "pm_last_refresh_scale_factor");
  pm_sync_state.refresh_commit_pending = readScalarU32Attribute(integrator_group.get(), "pm_refresh_commit_pending") != 0U;
  pm_sync_state.pending_refresh_opportunity = readScalarU64Attribute(integrator_group.get(), "pm_pending_refresh_opportunity");
  pm_sync_state.pending_refresh_field_version = readScalarU64Attribute(integrator_group.get(), "pm_pending_refresh_field_version");
  result.integrator_state.pm_sync_state.importPersistentState(pm_sync_state);
  if (schema_version >= k_restart_schema_v20) {
    Hdf5Handle force_cache_group(H5Gopen2(file.get(), "/gravity_force_cache", H5P_DEFAULT));
    result.gravity_force_cache.valid = readScalarU32Attribute(force_cache_group.get(), "valid") != 0U;
    result.gravity_force_cache.particle_id = readDataset1d<std::uint64_t>(force_cache_group.get(), "particle_id", H5T_NATIVE_UINT64);
    result.gravity_force_cache.gas_cell_id = readDataset1d<std::uint64_t>(force_cache_group.get(), "gas_cell_id", H5T_NATIVE_UINT64);
    result.gravity_force_cache.particle_accel_x_comoving = readDataset1d<double>(force_cache_group.get(), "particle_accel_x_comoving", H5T_NATIVE_DOUBLE);
    result.gravity_force_cache.particle_accel_y_comoving = readDataset1d<double>(force_cache_group.get(), "particle_accel_y_comoving", H5T_NATIVE_DOUBLE);
    result.gravity_force_cache.particle_accel_z_comoving = readDataset1d<double>(force_cache_group.get(), "particle_accel_z_comoving", H5T_NATIVE_DOUBLE);
    result.gravity_force_cache.cell_accel_x_comoving = readDataset1d<double>(force_cache_group.get(), "cell_accel_x_comoving", H5T_NATIVE_DOUBLE);
    result.gravity_force_cache.cell_accel_y_comoving = readDataset1d<double>(force_cache_group.get(), "cell_accel_y_comoving", H5T_NATIVE_DOUBLE);
    result.gravity_force_cache.cell_accel_z_comoving = readDataset1d<double>(force_cache_group.get(), "cell_accel_z_comoving", H5T_NATIVE_DOUBLE);
    validateGravityForceCacheForRestart(result.gravity_force_cache, result.state, "restart reader");
  }

  Hdf5Handle scheduler_group(H5Gopen2(file.get(), "/scheduler", H5P_DEFAULT));
  result.scheduler_state.current_tick = readScalarU64Attribute(scheduler_group.get(), "current_tick");
  result.scheduler_state.max_bin = static_cast<std::uint8_t>(readScalarU32Attribute(scheduler_group.get(), "max_bin"));
  result.scheduler_state.bin_index = readDataset1d<std::uint8_t>(scheduler_group.get(), "bin_index", H5T_NATIVE_UINT8);
  result.scheduler_state.next_activation_tick =
      readDataset1d<std::uint64_t>(scheduler_group.get(), "next_activation_tick", H5T_NATIVE_UINT64);
  result.scheduler_state.active_flag = readDataset1d<std::uint8_t>(scheduler_group.get(), "active_flag", H5T_NATIVE_UINT8);
  result.scheduler_state.pending_bin_index =
      readDataset1d<std::uint8_t>(scheduler_group.get(), "pending_bin_index", H5T_NATIVE_UINT8);
  validateRestartTimeBinMirrorsAgainstScheduler(result.state, result.scheduler_state, "restart reader");
  rebuildRestartTimeBinMirrorsFromScheduler(result.state, result.scheduler_state);
  if (schema_version >= k_restart_schema_v19) {
    Hdf5Handle gas_scheduler_group(H5Gopen2(file.get(), "/gas_cell_scheduler", H5P_DEFAULT));
    const std::string identity_key = readScalarStringAttribute(gas_scheduler_group.get(), "identity_key");
    if (identity_key != k_gas_cell_scheduler_identity_key) {
      throw std::runtime_error("restart reader: /gas_cell_scheduler identity_key must be gas_cell_id");
    }
    result.gas_cell_scheduler_state.current_tick =
        readScalarU64Attribute(gas_scheduler_group.get(), "current_tick");
    result.gas_cell_scheduler_state.max_bin =
        static_cast<std::uint8_t>(readScalarU32Attribute(gas_scheduler_group.get(), "max_bin"));
    result.gas_cell_scheduler_ids =
        readDataset1d<std::uint64_t>(gas_scheduler_group.get(), "gas_cell_id", H5T_NATIVE_UINT64);
    result.gas_cell_scheduler_state.bin_index =
        readDataset1d<std::uint8_t>(gas_scheduler_group.get(), "bin_index", H5T_NATIVE_UINT8);
    result.gas_cell_scheduler_state.next_activation_tick =
        readDataset1d<std::uint64_t>(gas_scheduler_group.get(), "next_activation_tick", H5T_NATIVE_UINT64);
    result.gas_cell_scheduler_state.active_flag =
        readDataset1d<std::uint8_t>(gas_scheduler_group.get(), "active_flag", H5T_NATIVE_UINT8);
    result.gas_cell_scheduler_state.pending_bin_index =
        readDataset1d<std::uint8_t>(gas_scheduler_group.get(), "pending_bin_index", H5T_NATIVE_UINT8);
    if (result.gas_cell_scheduler_state.current_tick != result.scheduler_state.current_tick) {
      throw std::runtime_error("restart reader: particle and gas-cell schedulers have different current_tick values");
    }
    validateGasCellTimeBinMirrorsAgainstScheduler(
        result.state,
        result.gas_cell_scheduler_state,
        result.gas_cell_scheduler_ids,
        "restart reader");
    rebuildGasCellTimeBinMirrorsFromScheduler(
        result.state, result.gas_cell_scheduler_state, result.gas_cell_scheduler_ids);
  } else {
    result.gas_cell_scheduler_state = legacyGasCellSchedulerStateFromMirrors(
        result.state, result.scheduler_state);
    result.gas_cell_scheduler_ids = gasCellIdsByDenseLocalRow(result.state, "legacy restart reader");
  }
  result.output_cadence_state = readOutputCadenceGroup(file.get());
  validateOutputCadenceStateForRestart(result.output_cadence_state, result.integrator_state, "restart reader");
  result.stochastic_state = readStochasticStateGroup(file.get());
  validateStochasticStateForRestart(result.stochastic_state, result.integrator_state, "restart reader");
  result.diagnostics = readRestartDiagnosticsGroup(file.get(), schema_version);
  if (result.diagnostics.restart_schema_name != schema_name ||
      result.diagnostics.restart_schema_version != schema_version ||
      result.diagnostics.step_index != result.integrator_state.step_index ||
      result.diagnostics.scheduler_current_tick != result.scheduler_state.current_tick ||
      result.diagnostics.scheduler_element_count != result.scheduler_state.bin_index.size() ||
      (schema_version >= k_restart_schema_v19 &&
       (result.diagnostics.gas_cell_scheduler_current_tick != result.gas_cell_scheduler_state.current_tick ||
        result.diagnostics.gas_cell_scheduler_element_count != result.gas_cell_scheduler_state.bin_index.size())) ||
      result.diagnostics.output_last_completed_step_index != result.output_cadence_state.last_completed_step_index ||
      result.diagnostics.stochastic_module_count != result.stochastic_state.modules.size()) {
    throw std::runtime_error("restart diagnostics summary is inconsistent with authoritative restart state");
  }
  readDistributedGravityGroup(file.get(), result.distributed_gravity_state);

  RestartWritePayload verify_payload;
  verify_payload.persistent_state.simulation_state = &result.state;
  verify_payload.integrator_state = &result.integrator_state;
  core::HierarchicalTimeBinScheduler verify_scheduler(result.scheduler_state.max_bin);
  verify_scheduler.importPersistentState(result.scheduler_state);
  core::HierarchicalTimeBinScheduler verify_gas_cell_scheduler(result.gas_cell_scheduler_state.max_bin);
  verify_gas_cell_scheduler.importPersistentState(result.gas_cell_scheduler_state);
  verify_payload.scheduler = &verify_scheduler;
  if (schema_version >= k_restart_schema_v20) {
    verify_payload.gravity_force_cache = &result.gravity_force_cache;
  }
  if (schema_version >= k_restart_schema_v19) {
    verify_payload.gas_cell_scheduler = &verify_gas_cell_scheduler;
  }
  verify_payload.normalized_config_hash_hex = result.normalized_config_hash_hex;
  verify_payload.normalized_config_text = result.normalized_config_text;
  verify_payload.provenance = result.provenance;
  verify_payload.distributed_gravity_state = result.distributed_gravity_state;
  verify_payload.output_cadence_state = result.output_cadence_state;
  verify_payload.stochastic_state = result.stochastic_state;

  const std::uint64_t computed_hash =
      restartPayloadIntegrityHashImpl(
          verify_payload,
          schema_version >= k_restart_schema_v15,
          schema_version >= k_restart_schema_v17,
          schema_version >= k_restart_schema_v18,
          schema_version >= k_restart_schema_v19,
          schema_version >= k_restart_schema_v20);
  if (computed_hash != result.payload_hash || hexU64(computed_hash) != result.payload_hash_hex) {
    throw std::runtime_error("restart payload integrity hash mismatch");
  }

  return result;
#endif
}

}  // namespace cosmosim::io
