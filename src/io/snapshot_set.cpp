#include "cosmosim/io/snapshot_hdf5.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <limits>
#include <regex>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/core/checked_arithmetic.hpp"
#include "cosmosim/io/io_contract.hpp"
#include "io/internal/transactional_file.hpp"

#if COSMOSIM_ENABLE_HDF5
#include <hdf5.h>
#endif

namespace cosmosim::io {
namespace {

#if COSMOSIM_ENABLE_HDF5
class H5Handle {
 public:
  explicit H5Handle(hid_t value = -1) : m_value(value) {}
  H5Handle(const H5Handle&) = delete;
  H5Handle& operator=(const H5Handle&) = delete;
  H5Handle(H5Handle&& other) noexcept : m_value(other.m_value) { other.m_value = -1; }
  H5Handle& operator=(H5Handle&& other) noexcept {
    if (this != &other) {
      close();
      m_value = other.m_value;
      other.m_value = -1;
    }
    return *this;
  }
  ~H5Handle() { close(); }
  [[nodiscard]] hid_t get() const noexcept { return m_value; }
  [[nodiscard]] bool valid() const noexcept { return m_value >= 0; }
 private:
  void close() noexcept {
    if (m_value < 0) return;
    switch (H5Iget_type(m_value)) {
      case H5I_FILE: H5Fclose(m_value); break;
      case H5I_GROUP: H5Gclose(m_value); break;
      case H5I_ATTR: H5Aclose(m_value); break;
      case H5I_DATATYPE: H5Tclose(m_value); break;
      default: break;
    }
    m_value = -1;
  }
  hid_t m_value = -1;
};

[[nodiscard]] bool readStringAttr(hid_t location, const char* key, std::string& out) {
  if (H5Aexists(location, key) <= 0) return false;
  H5Handle attr(H5Aopen(location, key, H5P_DEFAULT));
  if (!attr.valid()) return false;
  H5Handle type(H5Aget_type(attr.get()));
  if (!type.valid()) throw std::runtime_error(std::string("snapshot inspection: invalid attribute type ") + key);
  const std::size_t size = H5Tget_size(type.get());
  if (size == 0U || size > (1U << 24U)) {
    throw std::length_error(std::string("snapshot inspection: oversized string attribute ") + key);
  }
  std::string buffer(size, '\0');
  if (H5Aread(attr.get(), type.get(), buffer.data()) < 0) {
    throw std::runtime_error(std::string("snapshot inspection: failed reading attribute ") + key);
  }
  while (!buffer.empty() && buffer.back() == '\0') buffer.pop_back();
  out = std::move(buffer);
  return true;
}

[[nodiscard]] bool readU32Attr(hid_t location, const char* key, std::uint32_t& out) {
  H5Handle attr(H5Aopen(location, key, H5P_DEFAULT));
  return attr.valid() && H5Aread(attr.get(), H5T_NATIVE_UINT32, &out) >= 0;
}

[[nodiscard]] bool readDoubleAttr(hid_t location, const char* key, double& out) {
  H5Handle attr(H5Aopen(location, key, H5P_DEFAULT));
  return attr.valid() && H5Aread(attr.get(), H5T_NATIVE_DOUBLE, &out) >= 0;
}

[[nodiscard]] bool readCountArray6(
    hid_t location, const char* key, std::array<std::uint64_t, 6>& out,
    std::size_t& source_width_bytes) {
  if (H5Aexists(location, key) <= 0) return false;
  H5Handle attr(H5Aopen(location, key, H5P_DEFAULT));
  if (!attr.valid()) return false;
  H5Handle type(H5Aget_type(attr.get()));
  H5Handle space(H5Aget_space(attr.get()));
  if (!type.valid() || !space.valid() || H5Tget_class(type.get()) != H5T_INTEGER ||
      H5Sget_simple_extent_ndims(space.get()) != 1) {
    throw std::runtime_error(std::string("snapshot inspection: invalid count attribute ") + key);
  }
  hsize_t dims[1] = {0};
  if (H5Sget_simple_extent_dims(space.get(), dims, nullptr) < 0 || dims[0] != 6U) {
    throw std::runtime_error(std::string("snapshot inspection: count attribute must have six entries: ") + key);
  }
  source_width_bytes = H5Tget_size(type.get());
  if (source_width_bytes != 4U && source_width_bytes != 8U) {
    throw std::runtime_error(std::string("snapshot inspection: unsupported count width for ") + key);
  }
  if (H5Aread(attr.get(), H5T_NATIVE_UINT64, out.data()) < 0) {
    throw std::runtime_error(std::string("snapshot inspection: unreadable count attribute ") + key);
  }
  return true;
}

struct MemberHeader {
  std::filesystem::path path;
  std::array<std::uint64_t, 6> local{};
  std::array<std::uint64_t, 6> global{};
  std::uint32_t num_files = 1;
  std::uint32_t member_index = 0;
  bool has_member_index = false;
  double time = 1.0;
  double redshift = 0.0;
  double hubble = 0.0;
  std::string schema;
  std::uint32_t schema_version = 0;
  std::string generation;
  std::string file_kind;
  SnapshotDialect dialect = SnapshotDialect::kAuto;
};

[[nodiscard]] SnapshotDialect parseDialect(const std::string& value) {
  if (value == "chui_native") return SnapshotDialect::kChuiNative;
  if (value == "arepo_format3") return SnapshotDialect::kArepoFormat3;
  if (value == "gadget4_hdf5") return SnapshotDialect::kGadget4Hdf5;
  return SnapshotDialect::kAuto;
}

[[nodiscard]] MemberHeader inspectMember(const std::filesystem::path& path) {
  H5Handle file(H5Fopen(path.string().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
  if (!file.valid()) throw std::runtime_error("snapshot inspection: cannot open " + path.string());
  H5Handle header(H5Gopen2(file.get(), "/Header", H5P_DEFAULT));
  if (!header.valid()) throw std::runtime_error("snapshot inspection: missing /Header in " + path.string());

  MemberHeader result;
  result.path = path;
  const auto& names = sharedIoContractNames();
  static_cast<void>(readStringAttr(file.get(), std::string(names.file_kind_attribute).c_str(), result.file_kind));

  std::size_t local_width = 0U;
  if (!readCountArray6(header.get(), "NumPart_ThisFile", result.local, local_width)) {
    throw std::runtime_error("snapshot inspection: missing NumPart_ThisFile in " + path.string());
  }
  std::size_t global_width = 0U;
  if (!readCountArray6(header.get(), "NumPart_Total", result.global, global_width)) {
    result.global = result.local;
    global_width = local_width;
  }
  if (global_width == 4U) {
    std::array<std::uint64_t, 6> high{};
    std::size_t high_width = 0U;
    if (readCountArray6(header.get(), "NumPart_Total_HighWord", high, high_width)) {
      if (high_width != 4U) {
        throw std::runtime_error("snapshot inspection: NumPart_Total_HighWord must be 32-bit");
      }
      for (std::size_t i = 0; i < result.global.size(); ++i) {
        if (high[i] > std::numeric_limits<std::uint32_t>::max()) {
          throw std::runtime_error("snapshot inspection: invalid NumPart_Total_HighWord value");
        }
        result.global[i] |= (high[i] << 32U);
      }
    }
  }
  static_cast<void>(readU32Attr(header.get(), "NumFilesPerSnapshot", result.num_files));
  result.has_member_index =
      readU32Attr(header.get(), "CHUISnapshotMemberIndex", result.member_index);
  static_cast<void>(readDoubleAttr(header.get(), "Time", result.time));
  static_cast<void>(readDoubleAttr(header.get(), "Redshift", result.redshift));
  static_cast<void>(readDoubleAttr(header.get(), "HubbleParam", result.hubble));
  static_cast<void>(readStringAttr(header.get(), "CosmoSimSchemaName", result.schema));
  static_cast<void>(readU32Attr(header.get(), "CosmoSimSchemaVersion", result.schema_version));
  static_cast<void>(readStringAttr(header.get(), "CHUISnapshotGenerationID", result.generation));
  std::string dialect;
  if (readStringAttr(header.get(), "CHUISnapshotDialect", dialect)) result.dialect = parseDialect(dialect);
  if (result.num_files == 0U) throw std::runtime_error("snapshot inspection: NumFilesPerSnapshot is zero");
  return result;
}
#endif

[[nodiscard]] std::vector<std::filesystem::path> discoverMembers(const std::filesystem::path& input) {
  std::vector<std::filesystem::path> members;
  if (std::filesystem::is_directory(input)) {
    for (const auto& entry : std::filesystem::directory_iterator(input)) {
      if (entry.is_regular_file() && entry.path().extension() == ".hdf5") members.push_back(entry.path());
    }
  } else {
    members.push_back(input);
#if COSMOSIM_ENABLE_HDF5
    const MemberHeader first = inspectMember(input);
    if (first.num_files > 1U) {
      members.clear();
      const std::string filename = input.filename().string();
      const std::regex member_pattern(R"(^(.+)\.([0-9]+)\.hdf5$)");
      std::smatch match;
      if (!std::regex_match(filename, match, member_pattern)) {
        throw std::runtime_error(
            "snapshot inspection: multifile member name must end in .<index>.hdf5 or be opened through its snapdir");
      }
      const std::string stem = match[1].str();
      for (const auto& entry : std::filesystem::directory_iterator(input.parent_path())) {
        if (!entry.is_regular_file()) continue;
        std::smatch candidate;
        const std::string candidate_name = entry.path().filename().string();
        if (std::regex_match(candidate_name, candidate, member_pattern) && candidate[1].str() == stem) {
          members.push_back(entry.path());
        }
      }
    }
#endif
  }
  std::sort(members.begin(), members.end());
  if (members.empty()) throw std::runtime_error("snapshot inspection: no HDF5 members found");
  return members;
}

[[nodiscard]] bool completionMarkerMatches(
    const std::filesystem::path& input,
    const std::string& generation,
    const std::uint32_t num_files,
    const std::array<std::uint64_t, 6>& global_counts) {
  if (generation.empty()) return false;
  const std::filesystem::path directory =
      std::filesystem::is_directory(input) ? input : input.parent_path();
  const std::filesystem::path marker = directory / (generation + ".complete");
  std::error_code ec;
  if (!std::filesystem::is_regular_file(marker, ec) || ec) return false;
  const auto bytes = std::filesystem::file_size(marker, ec);
  if (ec || bytes == 0U || bytes > 64U * 1024U) return false;
  std::ifstream stream(marker);
  if (!stream) return false;
  std::string schema;
  std::string marker_generation;
  std::uint32_t marker_files = 0U;
  std::array<std::uint64_t, 6> marker_counts{};
  bool have_counts = false;
  std::string line;
  while (std::getline(stream, line)) {
    const auto eq = line.find('=');
    if (eq == std::string::npos) continue;
    const std::string key = line.substr(0, eq);
    const std::string value = line.substr(eq + 1U);
    try {
      if (key == "schema") schema = value;
      else if (key == "generation_id") marker_generation = value;
      else if (key == "num_files_per_snapshot") {
        const auto parsed = std::stoull(value);
        if (parsed > std::numeric_limits<std::uint32_t>::max()) return false;
        marker_files = static_cast<std::uint32_t>(parsed);
      } else if (key == "global_part_count") {
        std::size_t offset = 0U;
        for (std::size_t i = 0; i < marker_counts.size(); ++i) {
          const std::size_t comma = value.find(',', offset);
          const std::string token = value.substr(
              offset, comma == std::string::npos ? std::string::npos : comma - offset);
          if (token.empty()) return false;
          marker_counts[i] = std::stoull(token);
          if (i + 1U < marker_counts.size()) {
            if (comma == std::string::npos) return false;
            offset = comma + 1U;
          } else if (comma != std::string::npos) {
            return false;
          }
        }
        have_counts = true;
      }
    } catch (const std::exception&) {
      return false;
    }
  }
  return schema == "chui_snapshot_set_v1" && marker_generation == generation &&
      marker_files == num_files && have_counts && marker_counts == global_counts;
}

void appendState(core::SimulationState& out, const core::SimulationState& in, std::uint32_t owner_rank) {
  const std::size_t particle_offset = out.particles.size();
  const std::size_t cell_offset = out.cells.size();
  const std::size_t old_particles = out.particles.size();
  const std::size_t old_cells = out.cells.size();
  std::vector<core::GasCellIdentityRecord> identities(
      out.gas_cell_identity.records().begin(), out.gas_cell_identity.records().end());
  out.resizeParticles(core::checkedSizeAdd(old_particles, in.particles.size(), "snapshot set particles"));
  out.resizeCells(core::checkedSizeAdd(old_cells, in.cells.size(), "snapshot set cells"));

  for (std::size_t i = 0; i < in.particles.size(); ++i) {
    const std::size_t d = particle_offset + i;
    out.particles.position_x_comoving[d] = in.particles.position_x_comoving[i];
    out.particles.position_y_comoving[d] = in.particles.position_y_comoving[i];
    out.particles.position_z_comoving[d] = in.particles.position_z_comoving[i];
    out.particles.velocity_x_peculiar[d] = in.particles.velocity_x_peculiar[i];
    out.particles.velocity_y_peculiar[d] = in.particles.velocity_y_peculiar[i];
    out.particles.velocity_z_peculiar[d] = in.particles.velocity_z_peculiar[i];
    out.particles.mass_code[d] = in.particles.mass_code[i];
    out.particles.time_bin[d] = in.particles.time_bin[i];
    out.particle_sidecar.particle_id[d] = in.particle_sidecar.particle_id[i];
    out.particle_sidecar.sfc_key[d] = in.particle_sidecar.sfc_key[i];
    out.particle_sidecar.species_tag[d] = in.particle_sidecar.species_tag[i];
    out.particle_sidecar.particle_flags[d] = in.particle_sidecar.particle_flags[i];
    out.particle_sidecar.owning_rank[d] = owner_rank;
    out.particle_sidecar.last_drift_time_code[d] = in.particle_sidecar.last_drift_time_code[i];
    out.particle_sidecar.last_drift_scale_factor[d] = in.particle_sidecar.last_drift_scale_factor[i];
  }
  if (!in.particle_sidecar.gravity_softening_comoving.empty()) {
    if (out.particle_sidecar.gravity_softening_comoving.empty()) out.particle_sidecar.gravity_softening_comoving.resize(out.particles.size(), 0.0);
    else out.particle_sidecar.gravity_softening_comoving.resize(out.particles.size(), 0.0);
    for (std::size_t i = 0; i < in.particles.size(); ++i) out.particle_sidecar.gravity_softening_comoving[particle_offset+i] = in.particle_sidecar.gravity_softening_comoving[i];
  }
  if (!in.particle_sidecar.has_gravity_softening_override.empty()) {
    if (out.particle_sidecar.has_gravity_softening_override.empty()) out.particle_sidecar.has_gravity_softening_override.resize(out.particles.size(), 0U);
    else out.particle_sidecar.has_gravity_softening_override.resize(out.particles.size(), 0U);
    for (std::size_t i = 0; i < in.particles.size(); ++i) out.particle_sidecar.has_gravity_softening_override[particle_offset+i] = in.particle_sidecar.has_gravity_softening_override[i];
  }

  for (std::size_t i = 0; i < in.cells.size(); ++i) {
    const std::size_t d = cell_offset + i;
    out.cells.center_x_comoving[d] = in.cells.center_x_comoving[i];
    out.cells.center_y_comoving[d] = in.cells.center_y_comoving[i];
    out.cells.center_z_comoving[d] = in.cells.center_z_comoving[i];
    out.cells.mass_code[d] = in.cells.mass_code[i];
    out.cells.time_bin[d] = in.cells.time_bin[i];
    out.cells.patch_index[d] = in.cells.patch_index[i];
    out.gas_cells.gas_cell_id[d] = in.gas_cells.gas_cell_id[i];
    out.gas_cells.parent_particle_id[d] = in.gas_cells.parent_particle_id[i];
    out.gas_cells.velocity_x_peculiar[d] = in.gas_cells.velocity_x_peculiar[i];
    out.gas_cells.velocity_y_peculiar[d] = in.gas_cells.velocity_y_peculiar[i];
    out.gas_cells.velocity_z_peculiar[d] = in.gas_cells.velocity_z_peculiar[i];
    out.gas_cells.density_code[d] = in.gas_cells.density_code[i];
    out.gas_cells.pressure_code[d] = in.gas_cells.pressure_code[i];
    out.gas_cells.internal_energy_code[d] = in.gas_cells.internal_energy_code[i];
    out.gas_cells.metal_mass_code[d] = in.gas_cells.metal_mass_code[i];
    out.gas_cells.temperature_code[d] = in.gas_cells.temperature_code[i];
    out.gas_cells.sound_speed_code[d] = in.gas_cells.sound_speed_code[i];
  }

  for (const auto& record : in.gas_cell_identity.records()) {
    core::GasCellIdentityRecord copy = record;
    copy.local_cell_row = core::checkedLocalCellRow(cell_offset + record.local_cell_row, "snapshot set gas identity");
    identities.push_back(copy);
  }
  if (!identities.empty()) out.replaceGasCellIdentityRecords(std::move(identities));

  const std::size_t star_old = out.star_particles.size();
  out.star_particles.resize(star_old + in.star_particles.size());
  for (std::size_t i = 0; i < in.star_particles.size(); ++i) {
    const std::size_t d = star_old + i;
    out.star_particles.particle_index[d] = core::checkedLocalParticleRow(particle_offset + in.star_particles.particle_index[i], "snapshot set star index");
    out.star_particles.formation_scale_factor[d] = in.star_particles.formation_scale_factor[i];
    out.star_particles.birth_mass_code[d] = in.star_particles.birth_mass_code[i];
    out.star_particles.metallicity_mass_fraction[d] = in.star_particles.metallicity_mass_fraction[i];
    out.star_particles.birth_key[d] = in.star_particles.birth_key[i];
    out.star_particles.parent_gas_cell_id[d] = in.star_particles.parent_gas_cell_id[i];
    out.star_particles.birth_tick[d] = in.star_particles.birth_tick[i];
    out.star_particles.birth_ordinal[d] = in.star_particles.birth_ordinal[i];
    out.star_particles.stellar_age_years_last[d] = in.star_particles.stellar_age_years_last[i];
    out.star_particles.stellar_returned_mass_cumulative_code[d] = in.star_particles.stellar_returned_mass_cumulative_code[i];
    out.star_particles.stellar_returned_metals_cumulative_code[d] = in.star_particles.stellar_returned_metals_cumulative_code[i];
    out.star_particles.stellar_newly_synthesized_metals_cumulative_code[d] = in.star_particles.stellar_newly_synthesized_metals_cumulative_code[i];
    out.star_particles.stellar_feedback_energy_cumulative_erg[d] = in.star_particles.stellar_feedback_energy_cumulative_erg[i];
    out.star_particles.stellar_deposited_mass_cumulative_code[d] = in.star_particles.stellar_deposited_mass_cumulative_code[i];
    out.star_particles.stellar_deposited_metals_cumulative_code[d] = in.star_particles.stellar_deposited_metals_cumulative_code[i];
    out.star_particles.stellar_deposited_feedback_energy_cumulative_erg[d] = in.star_particles.stellar_deposited_feedback_energy_cumulative_erg[i];
  }

  const std::size_t bh_old = out.black_holes.size();
  out.black_holes.resize(bh_old + in.black_holes.size());
  for (std::size_t i = 0; i < in.black_holes.size(); ++i) {
    const std::size_t d = bh_old + i;
    out.black_holes.particle_index[d] = core::checkedLocalParticleRow(particle_offset + in.black_holes.particle_index[i], "snapshot set BH index");
    out.black_holes.host_cell_index[d] = in.black_holes.host_cell_index[i] < in.cells.size()
        ? core::checkedLocalCellRow(cell_offset + in.black_holes.host_cell_index[i], "snapshot set BH host")
        : in.black_holes.host_cell_index[i];
    out.black_holes.subgrid_mass_code[d] = in.black_holes.subgrid_mass_code[i];
    out.black_holes.accretion_rate_code[d] = in.black_holes.accretion_rate_code[i];
    out.black_holes.feedback_energy_code[d] = in.black_holes.feedback_energy_code[i];
    out.black_holes.eddington_ratio[d] = in.black_holes.eddington_ratio[i];
    out.black_holes.cumulative_accreted_mass_code[d] = in.black_holes.cumulative_accreted_mass_code[i];
    out.black_holes.cumulative_feedback_energy_code[d] = in.black_holes.cumulative_feedback_energy_code[i];
    out.black_holes.duty_cycle_active_time_code[d] = in.black_holes.duty_cycle_active_time_code[i];
    out.black_holes.duty_cycle_total_time_code[d] = in.black_holes.duty_cycle_total_time_code[i];
  }

  const std::size_t tracer_old = out.tracers.size();
  out.tracers.resize(tracer_old + in.tracers.size());
  for (std::size_t i = 0; i < in.tracers.size(); ++i) {
    const std::size_t d = tracer_old + i;
    out.tracers.particle_index[d] = core::checkedLocalParticleRow(particle_offset + in.tracers.particle_index[i], "snapshot set tracer index");
    out.tracers.parent_particle_id[d] = in.tracers.parent_particle_id[i];
    out.tracers.injection_step[d] = in.tracers.injection_step[i];
    out.tracers.host_cell_index[d] = in.tracers.host_cell_index[i] < in.cells.size()
        ? core::checkedLocalCellRow(cell_offset + in.tracers.host_cell_index[i], "snapshot set tracer host")
        : in.tracers.host_cell_index[i];
    out.tracers.mass_fraction_of_host[d] = in.tracers.mass_fraction_of_host[i];
    out.tracers.last_host_mass_code[d] = in.tracers.last_host_mass_code[i];
    out.tracers.cumulative_exchanged_mass_code[d] = in.tracers.cumulative_exchanged_mass_code[i];
  }
  out.species.count_by_species.fill(0U);
  for (const auto tag : out.particle_sidecar.species_tag) {
    if (tag >= out.species.count_by_species.size()) {
      throw std::runtime_error("snapshot set contains invalid species tag");
    }
    ++out.species.count_by_species[tag];
  }
  out.rebuildSpeciesIndex();
}

[[nodiscard]] SnapshotReadResult readSet(
    const std::filesystem::path& input_path,
    const core::SimulationConfig& config,
    SnapshotReadOptions options,
    bool require_chui) {
  const SnapshotSetInspection inspection = inspectSnapshotSet(input_path, options);
  if (options.require_complete_chui_set && require_chui && !inspection.complete) {
    throw std::runtime_error("snapshot set is incomplete or lacks its CHUI completion marker");
  }
  SnapshotReadResult merged;
  bool first = true;
  for (std::size_t member_index = 0; member_index < inspection.member_paths.size(); ++member_index) {
    SnapshotReadOptions member_options = options;
    member_options.require_complete_chui_set = false;
    SnapshotReadResult member = readGadgetArepoSnapshotHdf5(inspection.member_paths[member_index], config, member_options);
    if (require_chui && member.report.file_kind != sharedIoContractNames().science_snapshot_file_kind &&
        member.report.schema_name.rfind("gadget_arepo_v", 0) != 0) {
      throw std::runtime_error("readCosmoSimScienceSnapshotHdf5 rejected a non-CHUI snapshot member");
    }
    if (first) {
      merged.provenance = member.provenance;
      merged.normalized_config_text = member.normalized_config_text;
      merged.state.metadata = member.state.metadata;
      merged.report = member.report;
      first = false;
    } else {
      merged.report.defaulted_fields.insert(merged.report.defaulted_fields.end(), member.report.defaulted_fields.begin(), member.report.defaulted_fields.end());
      merged.report.unavailable_fields.insert(merged.report.unavailable_fields.end(), member.report.unavailable_fields.begin(), member.report.unavailable_fields.end());
    }
    appendState(merged.state, member.state, static_cast<std::uint32_t>(member_index));
  }
  merged.report.num_files_per_snapshot = inspection.num_files_per_snapshot;
  merged.report.local_part_count = inspection.global_part_count;
  merged.report.global_part_count = inspection.global_part_count;
  merged.report.generation_id = inspection.generation_id;
  merged.report.evolution_ready = merged.report.unavailable_fields.empty() && merged.state.validateUniqueParticleIds();
  if (!merged.state.validateUniqueParticleIds()) {
    throw std::runtime_error("snapshot set contains duplicate or zero particle IDs across members");
  }
  return merged;
}

}  // namespace

SnapshotSetInspection inspectSnapshotSet(
    const std::filesystem::path& input_path,
    const SnapshotReadOptions&) {
#if !COSMOSIM_ENABLE_HDF5
  static_cast<void>(input_path);
  throw std::runtime_error("COSMOSIM_ENABLE_HDF5=OFF: snapshot inspection unavailable");
#else
  const auto paths = discoverMembers(input_path);
  std::vector<MemberHeader> headers;
  headers.reserve(paths.size());
  for (const auto& path : paths) headers.push_back(inspectMember(path));
  for (std::size_t i = 0; i < headers.size(); ++i) {
    if (!headers[i].has_member_index) {
      headers[i].member_index = static_cast<std::uint32_t>(i);
    }
  }
  const MemberHeader& reference = headers.front();
  SnapshotSetInspection result;
  result.dialect = reference.dialect;
  result.schema_name = reference.schema;
  result.schema_version = reference.schema_version;
  result.num_files_per_snapshot = reference.num_files;
  result.global_part_count = reference.global;
  result.scale_factor = reference.time;
  result.redshift = reference.redshift;
  result.hubble_param = reference.hubble;
  result.generation_id = reference.generation;
  std::array<std::uint64_t, 6> summed{};
  std::set<std::uint32_t> member_indices;
  for (const auto& header : headers) {
    if (header.num_files != reference.num_files || header.global != reference.global ||
        header.schema != reference.schema || header.schema_version != reference.schema_version ||
        header.generation != reference.generation || std::abs(header.time - reference.time) > 1.0e-12 ||
        std::abs(header.redshift - reference.redshift) > 1.0e-12 ||
        std::abs(header.hubble - reference.hubble) > 1.0e-12) {
      throw std::runtime_error("snapshot set members disagree on global schema/epoch metadata");
    }
    if (!member_indices.insert(header.member_index).second) {
      throw std::runtime_error("snapshot set contains duplicate member indices");
    }
    for (std::size_t i = 0; i < 6; ++i) {
      if (header.local[i] > std::numeric_limits<std::uint64_t>::max() - summed[i]) {
        throw std::overflow_error("snapshot set local-count sum overflow");
      }
      summed[i] += header.local[i];
    }
    result.member_paths.push_back(header.path);
  }
  if (headers.size() > reference.num_files) throw std::runtime_error("snapshot set contains more members than NumFilesPerSnapshot");
  const bool counts_match = summed == reference.global;
  const bool members_complete = headers.size() == reference.num_files;
  const bool chui_multifile = reference.file_kind == sharedIoContractNames().science_snapshot_file_kind && reference.num_files > 1U;
  const bool marker_ok = !chui_multifile || completionMarkerMatches(
      input_path, reference.generation, reference.num_files, reference.global);
  result.complete = members_complete && counts_match && marker_ok;
  return result;
#endif
}

SnapshotReadResult readCosmoSimScienceSnapshotHdf5(
    const std::filesystem::path& input_path,
    const core::SimulationConfig& config,
    const SnapshotReadOptions& options) {
  return readSet(input_path, config, options, true);
}

SnapshotReadResult importExternalSnapshotHdf5(
    const std::filesystem::path& input_path,
    const core::SimulationConfig& config,
    const SnapshotReadOptions& options) {
  if (options.dialect == SnapshotDialect::kAuto) {
    throw std::invalid_argument("external snapshot import requires explicit SnapshotDialect");
  }
  SnapshotReadOptions external = options;
  external.require_complete_chui_set = false;
  return readSet(input_path, config, external, false);
}


void writeSnapshotSetCompletionMarker(
    const std::filesystem::path& snapshot_directory,
    std::string_view generation_id,
    std::uint32_t num_files_per_snapshot,
    const std::array<std::uint64_t, 6>& global_part_count,
    bool durable_publication) {
  if (generation_id.empty()) {
    throw std::invalid_argument("snapshot completion marker requires a generation id");
  }
  if (num_files_per_snapshot == 0U) {
    throw std::invalid_argument("snapshot completion marker requires at least one member");
  }
  std::string body;
  body += "schema=chui_snapshot_set_v1\n";
  body += "generation_id=" + std::string(generation_id) + "\n";
  body += "num_files_per_snapshot=" + std::to_string(num_files_per_snapshot) + "\n";
  body += "global_part_count=";
  for (std::size_t i = 0; i < global_part_count.size(); ++i) {
    if (i != 0U) body += ',';
    body += std::to_string(global_part_count[i]);
  }
  body += "\n";
  internal::writeTextFileTransactionally(
      snapshot_directory / (std::string(generation_id) + ".complete"), body,
      durable_publication ? internal::FileDurability::kDurablePublication
                          : internal::FileDurability::kAtomicVisibility);
}

}  // namespace cosmosim::io
