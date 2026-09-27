#include "cosmosim/parallel/distributed_memory.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <iomanip>
#include <limits>
#include <numeric>
#include <optional>
#include <sstream>
#include <streambuf>
#include <unordered_map>
#include <unordered_set>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/core/memory_governor.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "parallel/internal/distributed_memory_support.hpp"

#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
#include <mpi.h>
#endif

namespace cosmosim::parallel {
namespace {

// A counting streambuf measures the exact emitted byte count without retaining
// any text. The append streambuf then writes directly into the sole result
// string. Both use checked sizes; neither materializes per-record strings.
class RestartCountingStreambuf final : public std::streambuf {
 public:
  [[nodiscard]] std::size_t count() const noexcept { return m_count; }
 protected:
  std::streamsize xsputn(const char*, std::streamsize count) override {
    if (count < 0) throw std::length_error("negative restart serialization length");
    m_count = core::checkedSizeAdd(m_count, static_cast<std::size_t>(count),
                                   "distributed restart serialization length");
    return count;
  }
  int_type overflow(int_type ch) override {
    if (!traits_type::eq_int_type(ch, traits_type::eof())) {
      m_count = core::checkedSizeAdd(m_count, 1U,
                                    "distributed restart serialization length");
    }
    return traits_type::not_eof(ch);
  }
 private:
  std::size_t m_count = 0U;
};

class RestartStringStreambuf final : public std::streambuf {
 public:
  explicit RestartStringStreambuf(std::string& output) : m_output(output) {}
 protected:
  std::streamsize xsputn(const char* bytes, std::streamsize count) override {
    if (count < 0) throw std::length_error("negative restart serialization length");
    const std::size_t size = static_cast<std::size_t>(count);
    (void)core::checkedSizeAdd(m_output.size(), size,
                               "distributed restart serialization length");
    m_output.append(bytes, size);
    return count;
  }
  int_type overflow(int_type ch) override {
    if (!traits_type::eq_int_type(ch, traits_type::eof())) {
      (void)core::checkedSizeAdd(m_output.size(), 1U,
                                 "distributed restart serialization length");
      m_output.push_back(traits_type::to_char_type(ch));
    }
    return traits_type::not_eof(ch);
  }
 private:
  std::string& m_output;
};

}  // namespace

void DistributedRestartState::serializeTo(std::ostream& stream) const {
  if (pm_slab_begin_x_by_rank.size() != pm_slab_end_x_by_rank.size()) {
    throw std::invalid_argument("distributed restart PM slab table extents differ");
  }
  stream << std::setprecision(std::numeric_limits<double>::max_digits10);
  stream << "schema_version=" << schema_version << '\n';
  stream << "decomposition_epoch=" << decomposition_epoch << '\n';
  stream << "world_size=" << world_size << '\n';
  stream << "pm_grid_nx=" << pm_grid_nx << '\n';
  stream << "pm_grid_ny=" << pm_grid_ny << '\n';
  stream << "pm_grid_nz=" << pm_grid_nz << '\n';
  stream << "pm_decomposition_mode=" << pm_decomposition_mode << '\n';
  stream << "gravity_kick_opportunity=" << gravity_kick_opportunity << '\n';
  stream << "pm_update_cadence_steps=" << pm_update_cadence_steps << '\n';
  stream << "long_range_field_version=" << long_range_field_version << '\n';
  stream << "last_long_range_refresh_opportunity=" << last_long_range_refresh_opportunity << '\n';
  stream << "long_range_field_built_step_index=" << long_range_field_built_step_index << '\n';
  stream << "long_range_field_built_scale_factor=" << long_range_field_built_scale_factor << '\n';
  stream << "long_range_restart_policy=" << long_range_restart_policy << '\n';
  stream << "item_count=" << owning_rank_by_item.size() << '\n';
  for (std::size_t i = 0; i < owning_rank_by_item.size(); ++i) {
    stream << "rank[" << i << "]=" << owning_rank_by_item[i] << '\n';
  }
  stream << "pm_slab_rank_count=" << pm_slab_begin_x_by_rank.size() << '\n';
  for (std::size_t rank = 0; rank < pm_slab_begin_x_by_rank.size(); ++rank) {
    stream << "pm_slab_begin_x[" << rank << "]=" << pm_slab_begin_x_by_rank[rank] << '\n';
    stream << "pm_slab_end_x[" << rank << "]=" << pm_slab_end_x_by_rank[rank] << '\n';
  }
 }

std::size_t DistributedRestartState::serializedSizeBytes() const {
  RestartCountingStreambuf buffer;
  std::ostream stream(&buffer);
  stream.exceptions(std::ios::badbit | std::ios::failbit);
  serializeTo(stream);
  return buffer.count();
}

std::string DistributedRestartState::serialize() const {
  const std::size_t bytes = serializedSizeBytes();
  std::string output;
  if (bytes > output.max_size()) {
    throw std::length_error("distributed restart serialization exceeds string max_size");
  }
  output.reserve(bytes);
  RestartStringStreambuf buffer(output);
  std::ostream stream(&buffer);
  stream.exceptions(std::ios::badbit | std::ios::failbit);
  serializeTo(stream);
  if (output.size() != bytes) {
    throw std::logic_error("distributed restart serialization changed size between passes");
  }
  return output;
}

DistributedRestartState DistributedRestartState::deserialize(const std::string& encoded) {
  DistributedRestartState state;
  // Parse one line at a time without retaining an O(N_item) vector of heap
  // strings. The input remains immutable so existing restart decoding and
  // integrity semantics are unchanged.
  std::size_t line_begin = 0U;
  std::size_t expected_item_count = 0;
  std::vector<bool> seen_rank_entry;
  std::size_t expected_slab_rank_count = 0;
  std::vector<bool> seen_slab_begin;
  std::vector<bool> seen_slab_end;

  while (line_begin < encoded.size()) {
    const std::size_t line_end = encoded.find('\n', line_begin);
    const std::size_t end = line_end == std::string::npos ? encoded.size() : line_end;
    const std::string line(encoded.data() + line_begin, end - line_begin);
    line_begin = line_end == std::string::npos ? encoded.size() : line_end + 1U;
    if (line.empty()) continue;

    const std::size_t eq = line.find('=');
    if (eq == std::string::npos) {
      throw std::invalid_argument("invalid restart encoding line");
    }
    const std::string key = line.substr(0, eq);
    const std::string value = line.substr(eq + 1);

    if (key == "schema_version") {
      state.schema_version = static_cast<std::uint32_t>(std::stoul(value));
    } else if (key == "decomposition_epoch") {
      state.decomposition_epoch = std::stoull(value);
    } else if (key == "world_size") {
      state.world_size = std::stoi(value);
    } else if (key == "pm_grid_nx") {
      state.pm_grid_nx = static_cast<std::size_t>(std::stoull(value));
    } else if (key == "pm_grid_ny") {
      state.pm_grid_ny = static_cast<std::size_t>(std::stoull(value));
    } else if (key == "pm_grid_nz") {
      state.pm_grid_nz = static_cast<std::size_t>(std::stoull(value));
    } else if (key == "pm_decomposition_mode") {
      state.pm_decomposition_mode = value;
    } else if (key == "gravity_kick_opportunity") {
      state.gravity_kick_opportunity = std::stoull(value);
    } else if (key == "pm_update_cadence_steps") {
      state.pm_update_cadence_steps = std::stoull(value);
    } else if (key == "long_range_field_version") {
      state.long_range_field_version = std::stoull(value);
    } else if (key == "last_long_range_refresh_opportunity") {
      state.last_long_range_refresh_opportunity = std::stoull(value);
    } else if (key == "long_range_field_built_step_index") {
      state.long_range_field_built_step_index = std::stoull(value);
    } else if (key == "long_range_field_built_scale_factor") {
      state.long_range_field_built_scale_factor = std::stod(value);
    } else if (key == "long_range_restart_policy") {
      state.long_range_restart_policy = value;
    } else if (key == "item_count") {
      expected_item_count = core::checkedIntegralNarrow<std::size_t>(
          std::stoull(value), "restart item count");
      if (expected_item_count > encoded.size() ||
          expected_item_count > state.owning_rank_by_item.max_size()) {
        throw std::length_error("restart item count exceeds encoded metadata capacity");
      }
      state.owning_rank_by_item.assign(expected_item_count, 0);
      seen_rank_entry.assign(expected_item_count, false);
    } else if (key == "pm_slab_rank_count") {
      expected_slab_rank_count = core::checkedIntegralNarrow<std::size_t>(
          std::stoull(value), "restart slab rank count");
      if (expected_slab_rank_count > encoded.size() ||
          expected_slab_rank_count > state.pm_slab_begin_x_by_rank.max_size()) {
        throw std::length_error("restart slab count exceeds encoded metadata capacity");
      }
      state.pm_slab_begin_x_by_rank.assign(expected_slab_rank_count, 0);
      state.pm_slab_end_x_by_rank.assign(expected_slab_rank_count, 0);
      seen_slab_begin.assign(expected_slab_rank_count, false);
      seen_slab_end.assign(expected_slab_rank_count, false);
    } else if (key.rfind("rank[", 0) == 0) {
      const std::size_t open = key.find('[');
      const std::size_t close = key.find(']');
      if (open == std::string::npos || close == std::string::npos || close <= open + 1) {
        throw std::invalid_argument("invalid rank entry in restart encoding");
      }
      const std::size_t index = static_cast<std::size_t>(std::stoull(key.substr(open + 1, close - open - 1)));
      if (index >= state.owning_rank_by_item.size()) {
        throw std::out_of_range("restart rank index out of bounds");
      }
      if (seen_rank_entry[index]) {
        throw std::invalid_argument("duplicate restart rank entry");
      }
      state.owning_rank_by_item[index] = std::stoi(value);
      seen_rank_entry[index] = true;
    } else if (key.rfind("pm_slab_begin_x[", 0) == 0 || key.rfind("pm_slab_end_x[", 0) == 0) {
      const bool is_begin = key.rfind("pm_slab_begin_x[", 0) == 0;
      const std::size_t open = key.find('[');
      const std::size_t close = key.find(']');
      if (open == std::string::npos || close == std::string::npos || close <= open + 1) {
        throw std::invalid_argument("invalid PM slab entry in restart encoding");
      }
      const std::size_t rank_index = static_cast<std::size_t>(std::stoull(key.substr(open + 1, close - open - 1)));
      if (rank_index >= expected_slab_rank_count) {
        throw std::out_of_range("restart PM slab rank index out of bounds");
      }
      if (is_begin) {
        if (seen_slab_begin[rank_index]) {
          throw std::invalid_argument("duplicate PM slab begin entry");
        }
        state.pm_slab_begin_x_by_rank[rank_index] = static_cast<std::size_t>(std::stoull(value));
        seen_slab_begin[rank_index] = true;
      } else {
        if (seen_slab_end[rank_index]) {
          throw std::invalid_argument("duplicate PM slab end entry");
        }
        state.pm_slab_end_x_by_rank[rank_index] = static_cast<std::size_t>(std::stoull(value));
        seen_slab_end[rank_index] = true;
      }
    }
  }

  if (state.owning_rank_by_item.size() != expected_item_count) {
    throw std::runtime_error("restart decode item count mismatch");
  }
  if (state.world_size <= 0) {
    throw std::invalid_argument("restart world_size must be positive");
  }
  if (state.pm_update_cadence_steps == 0) {
    throw std::invalid_argument("restart PM cadence must be >= 1");
  }
  if (state.last_long_range_refresh_opportunity > state.gravity_kick_opportunity) {
    throw std::invalid_argument("restart cadence state is inconsistent: last refresh opportunity exceeds kick opportunity");
  }
  if (state.long_range_field_version == 0 && state.last_long_range_refresh_opportunity != 0) {
    throw std::invalid_argument("restart cadence state is inconsistent: non-zero refresh opportunity with zero field version");
  }
  if (state.long_range_restart_policy != "deterministic_rebuild") {
    throw std::invalid_argument("restart long-range policy is unsupported: " + state.long_range_restart_policy);
  }
  if (state.schema_version >= 2) {
    if (state.pm_grid_nx == 0 || state.pm_grid_ny == 0 || state.pm_grid_nz == 0) {
      throw std::invalid_argument("restart PM grid dimensions must be > 0 for schema_version >= 2");
    }
    if (state.pm_decomposition_mode.empty()) {
      throw std::invalid_argument("restart PM decomposition mode must be non-empty");
    }
    if (expected_slab_rank_count != static_cast<std::size_t>(state.world_size)) {
      throw std::invalid_argument("restart PM slab rank count must match world_size");
    }
    for (std::size_t rank = 0; rank < expected_slab_rank_count; ++rank) {
      if (!seen_slab_begin[rank] || !seen_slab_end[rank]) {
        throw std::runtime_error("restart decode missing PM slab ownership entry");
      }
      if (state.pm_slab_end_x_by_rank[rank] < state.pm_slab_begin_x_by_rank[rank]) {
        throw std::invalid_argument("restart PM slab end_x must be >= begin_x");
      }
    }
  }
  for (bool seen : seen_rank_entry) {
    if (!seen) {
      throw std::runtime_error("restart decode missing ownership entry");
    }
  }
  for (std::size_t index = 0; index < state.owning_rank_by_item.size(); ++index) {
    const int rank = state.owning_rank_by_item[index];
    if (rank < 0 || rank >= state.world_size) {
      throw std::invalid_argument(
          "restart ownership entry rank is outside world_size at item " +
          std::to_string(index) + ": rank=" + std::to_string(rank) +
          ", world_size=" + std::to_string(state.world_size));
    }
  }
  return state;
}

DistributedRestartCompatibilityReport evaluateDistributedRestartCompatibility(
    const DistributedRestartState& restart_state,
    const DistributedExecutionTopology& runtime_topology) {
  DistributedRestartCompatibilityReport report;
  if (restart_state.schema_version != 2) {
    report.supported_schema_match = false;
    report.mismatch_messages.push_back(
        "distributed restart schema mismatch: expected=2, observed=" +
        std::to_string(restart_state.schema_version));
  }
  if (restart_state.world_size != runtime_topology.world_size) {
    report.world_size_match = false;
    report.mismatch_messages.push_back(
        "world_size mismatch: restart=" + std::to_string(restart_state.world_size) +
        ", runtime=" + std::to_string(runtime_topology.world_size));
  }
  if (restart_state.pm_slab_begin_x_by_rank.size() != restart_state.pm_slab_end_x_by_rank.size() ||
      restart_state.pm_slab_begin_x_by_rank.size() != static_cast<std::size_t>(std::max(restart_state.world_size, 0))) {
    report.pm_slab_table_shape_match = false;
    report.mismatch_messages.push_back(
        "PM slab table mismatch: begin_count=" +
        std::to_string(restart_state.pm_slab_begin_x_by_rank.size()) +
        ", end_count=" + std::to_string(restart_state.pm_slab_end_x_by_rank.size()) +
        ", restart_world_size=" + std::to_string(restart_state.world_size));
  }
  if (restart_state.pm_grid_nx != runtime_topology.pm_slab.global_nx ||
      restart_state.pm_grid_ny != runtime_topology.pm_slab.global_ny ||
      restart_state.pm_grid_nz != runtime_topology.pm_slab.global_nz) {
    report.pm_grid_shape_match = false;
    report.mismatch_messages.push_back(
        "PM grid mismatch: restart=(" + std::to_string(restart_state.pm_grid_nx) + "," +
        std::to_string(restart_state.pm_grid_ny) + "," + std::to_string(restart_state.pm_grid_nz) +
        "), runtime=(" + std::to_string(runtime_topology.pm_slab.global_nx) + "," +
        std::to_string(runtime_topology.pm_slab.global_ny) + "," +
        std::to_string(runtime_topology.pm_slab.global_nz) + ")");
  }
  if (restart_state.pm_decomposition_mode != runtime_topology.pm_decomposition_mode) {
    report.pm_decomposition_mode_match = false;
    report.mismatch_messages.push_back(
        "PM decomposition mode mismatch: restart=" + restart_state.pm_decomposition_mode +
        ", runtime=" + runtime_topology.pm_decomposition_mode);
  }
  if (restart_state.pm_update_cadence_steps == 0) {
    report.pm_cadence_steps_match = false;
    report.mismatch_messages.push_back(
        "PM cadence mismatch: restart cadence must be >= 1, got 0");
  }
  if (restart_state.last_long_range_refresh_opportunity > restart_state.gravity_kick_opportunity) {
    report.gravity_kick_state_match = false;
    report.mismatch_messages.push_back(
        "gravity kick mismatch: last refresh opportunity exceeds current kick opportunity");
  }
  if ((restart_state.long_range_field_version == 0) !=
      (restart_state.last_long_range_refresh_opportunity == 0)) {
    report.long_range_field_state_match = false;
    report.mismatch_messages.push_back(
        "long-range field mismatch: field version and refresh opportunity are inconsistent");
  }
  if (runtime_topology.world_rank < 0 ||
      runtime_topology.world_rank >= static_cast<int>(restart_state.pm_slab_begin_x_by_rank.size())) {
    report.pm_local_slab_match = false;
    report.mismatch_messages.push_back(
        "runtime world_rank is outside restart PM slab ownership table: rank=" +
        std::to_string(runtime_topology.world_rank) +
        ", table_size=" + std::to_string(restart_state.pm_slab_begin_x_by_rank.size()));
    return report;
  }
  const std::size_t rank = static_cast<std::size_t>(runtime_topology.world_rank);
  const std::size_t begin_restart = restart_state.pm_slab_begin_x_by_rank[rank];
  const std::size_t end_restart = restart_state.pm_slab_end_x_by_rank[rank];
  if (begin_restart != runtime_topology.pm_slab.owned_x.begin_x ||
      end_restart != runtime_topology.pm_slab.owned_x.end_x) {
    report.pm_local_slab_match = false;
    report.mismatch_messages.push_back(
        "PM slab mismatch for rank " + std::to_string(runtime_topology.world_rank) +
        ": restart=[" + std::to_string(begin_restart) + "," + std::to_string(end_restart) +
        "), runtime=[" + std::to_string(runtime_topology.pm_slab.owned_x.begin_x) + "," +
        std::to_string(runtime_topology.pm_slab.owned_x.end_x) + ")");
  }
  return report;
}


}  // namespace cosmosim::parallel
