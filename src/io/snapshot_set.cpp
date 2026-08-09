#include "cosmosim/io/snapshot_hdf5.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <optional>
#include <regex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "core/internal/sha256.hpp"
#include "cosmosim/core/build_config.hpp"
#include "cosmosim/core/checked_arithmetic.hpp"
#include "cosmosim/core/memory_accounting.hpp"
#include "cosmosim/io/io_contract.hpp"
#include "io/internal/file_sha256.hpp"
#include "io/internal/snapshot_readiness.hpp"
#include "io/internal/snapshot_set_internal.hpp"
#include "io/internal/transactional_file.hpp"

#if COSMOSIM_ENABLE_HDF5
#include <hdf5.h>
#endif

namespace cosmosim::io {
namespace {

constexpr std::uint64_t kMaxSnapshotSetManifestBytes = 4ULL << 20U;

[[nodiscard]] std::string dialectLabel(SnapshotDialect dialect) {
  switch (dialect) {
    case SnapshotDialect::kAuto: return "auto";
    case SnapshotDialect::kChuiNative: return "chui_native";
    case SnapshotDialect::kArepoFormat3: return "arepo_format3";
    case SnapshotDialect::kGadget4Hdf5: return "gadget4_hdf5";
  }
  throw std::logic_error("unhandled SnapshotDialect");
}

[[nodiscard]] SnapshotDialect parseDialect(const std::string& value) {
  if (value == "chui_native") return SnapshotDialect::kChuiNative;
  if (value == "arepo_format3") return SnapshotDialect::kArepoFormat3;
  if (value == "gadget4_hdf5") return SnapshotDialect::kGadget4Hdf5;
  if (value.empty() || value == "auto") return SnapshotDialect::kAuto;
  throw std::runtime_error("snapshot inspection: unknown snapshot dialect '" + value + "'");
}

[[nodiscard]] bool nearlyEqual(double lhs, double rhs) {
  if (std::isnan(lhs) || std::isnan(rhs)) return std::isnan(lhs) && std::isnan(rhs);
  if (!std::isfinite(lhs) || !std::isfinite(rhs)) return lhs == rhs;
  const double scale = std::max({1.0, std::abs(lhs), std::abs(rhs)});
  return std::abs(lhs - rhs) <= 64.0 * std::numeric_limits<double>::epsilon() * scale;
}

[[nodiscard]] std::string formatDouble(double value) {
  std::ostringstream out;
  out << std::setprecision(std::numeric_limits<double>::max_digits10) << value;
  return out.str();
}

[[nodiscard]] std::string formatCountArray(const std::array<std::uint64_t, 6>& counts) {
  std::string out;
  for (std::size_t i = 0; i < counts.size(); ++i) {
    if (i != 0U) out += ',';
    out += std::to_string(counts[i]);
  }
  return out;
}

[[nodiscard]] std::array<std::uint64_t, 6> parseCountArray(std::string_view value) {
  std::array<std::uint64_t, 6> counts{};
  std::size_t offset = 0U;
  for (std::size_t i = 0; i < counts.size(); ++i) {
    const std::size_t comma = value.find(',', offset);
    const std::string token(value.substr(
        offset, comma == std::string_view::npos ? std::string_view::npos : comma - offset));
    if (token.empty()) throw std::runtime_error("snapshot set manifest contains an empty count token");
    std::size_t consumed = 0U;
    counts[i] = std::stoull(token, &consumed);
    if (consumed != token.size()) throw std::runtime_error("snapshot set manifest count is not an integer");
    if (i + 1U < counts.size()) {
      if (comma == std::string_view::npos) throw std::runtime_error("snapshot set manifest count vector is truncated");
      offset = comma + 1U;
    } else if (comma != std::string_view::npos) {
      throw std::runtime_error("snapshot set manifest count vector has extra entries");
    }
  }
  return counts;
}

void validateManifestAtom(std::string_view value, std::string_view label) {
  if (value.empty() || value.find('\n') != std::string_view::npos ||
      value.find('\r') != std::string_view::npos) {
    throw std::invalid_argument("snapshot set " + std::string(label) + " is empty or contains a line break");
  }
}

[[nodiscard]] std::map<std::string, std::string> parseKeyValueBody(std::string_view body) {
  std::map<std::string, std::string> values;
  std::istringstream stream{std::string(body)};
  std::string line;
  while (std::getline(stream, line)) {
    if (line.empty()) continue;
    const auto eq = line.find('=');
    if (eq == std::string::npos || eq == 0U) {
      throw std::runtime_error("snapshot set manifest contains a malformed line");
    }
    const std::string key = line.substr(0, eq);
    const std::string value = line.substr(eq + 1U);
    if (!values.emplace(key, value).second) {
      throw std::runtime_error("snapshot set manifest contains duplicate key: " + key);
    }
  }
  return values;
}

[[nodiscard]] std::string readBoundedTextFile(
    const std::filesystem::path& path,
    std::uint64_t max_bytes = kMaxSnapshotSetManifestBytes) {
  std::error_code ec;
  const auto bytes = std::filesystem::file_size(path, ec);
  if (ec || bytes == 0U || bytes > max_bytes) {
    throw std::runtime_error("snapshot set metadata file is missing, empty, or oversized: " + path.string());
  }
  std::ifstream input(path, std::ios::binary);
  if (!input) throw std::runtime_error("failed to open snapshot set metadata: " + path.string());
  std::string text(static_cast<std::size_t>(bytes), '\0');
  input.read(text.data(), static_cast<std::streamsize>(text.size()));
  if (!input || input.gcount() != static_cast<std::streamsize>(text.size())) {
    throw std::runtime_error("failed to read snapshot set metadata: " + path.string());
  }
  return text;
}

[[nodiscard]] std::uint32_t parseU32(const std::map<std::string, std::string>& values, const std::string& key) {
  const auto it = values.find(key);
  if (it == values.end()) throw std::runtime_error("snapshot set manifest missing key: " + key);
  std::size_t consumed = 0U;
  const auto parsed = std::stoull(it->second, &consumed);
  if (consumed != it->second.size() || parsed > std::numeric_limits<std::uint32_t>::max()) {
    throw std::runtime_error("snapshot set manifest has invalid uint32 key: " + key);
  }
  return static_cast<std::uint32_t>(parsed);
}

[[nodiscard]] std::uint64_t parseU64(const std::map<std::string, std::string>& values, const std::string& key) {
  const auto it = values.find(key);
  if (it == values.end()) throw std::runtime_error("snapshot set manifest missing key: " + key);
  std::size_t consumed = 0U;
  const auto parsed = std::stoull(it->second, &consumed);
  if (consumed != it->second.size()) throw std::runtime_error("snapshot set manifest has invalid uint64 key: " + key);
  return parsed;
}

[[nodiscard]] double parseDouble(const std::map<std::string, std::string>& values, const std::string& key) {
  const auto it = values.find(key);
  if (it == values.end()) throw std::runtime_error("snapshot set manifest missing key: " + key);
  std::size_t consumed = 0U;
  const double parsed = std::stod(it->second, &consumed);
  if (consumed != it->second.size() || !std::isfinite(parsed)) {
    throw std::runtime_error("snapshot set manifest has invalid finite double key: " + key);
  }
  return parsed;
}

[[nodiscard]] const std::string& requireString(
    const std::map<std::string, std::string>& values,
    const std::string& key) {
  const auto it = values.find(key);
  if (it == values.end()) throw std::runtime_error("snapshot set manifest missing key: " + key);
  return it->second;
}

struct MemberHeader {
  std::filesystem::path path;
  std::array<std::uint64_t, 6> local{};
  std::array<std::uint64_t, 6> global{};
  std::uint32_t num_files = 1;
  std::uint32_t member_index = 0;
  bool has_member_index = false;
  double time = std::numeric_limits<double>::quiet_NaN();
  double redshift = std::numeric_limits<double>::quiet_NaN();
  double box_x = std::numeric_limits<double>::quiet_NaN();
  double box_y = std::numeric_limits<double>::quiet_NaN();
  double box_z = std::numeric_limits<double>::quiet_NaN();
  double omega_matter = std::numeric_limits<double>::quiet_NaN();
  double omega_lambda = std::numeric_limits<double>::quiet_NaN();
  double omega_baryon = std::numeric_limits<double>::quiet_NaN();
  double hubble = std::numeric_limits<double>::quiet_NaN();
  std::string schema;
  std::uint32_t schema_version = 0;
  std::string generation;
  std::string file_kind;
  SnapshotDialect dialect = SnapshotDialect::kAuto;
  std::string unit_length;
  std::string unit_mass;
  std::string unit_velocity;
  std::string coordinate_frame;
  std::string velocity_storage_convention;
  std::string config_hash_hex;
  std::string naming_rules_version;
  std::string file_naming_rules_version;
  std::uint64_t file_size = 0;
};

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
      case H5I_DATASPACE: H5Sclose(m_value); break;
      default: break;
    }
    m_value = -1;
  }
  hid_t m_value = -1;
};

[[nodiscard]] bool readStringAttr(
    hid_t location,
    const char* key,
    std::string& out,
    std::uint64_t max_bytes) {
  if (H5Aexists(location, key) <= 0) return false;
  H5Handle attr(H5Aopen(location, key, H5P_DEFAULT));
  if (!attr.valid()) return false;
  H5Handle type(H5Aget_type(attr.get()));
  if (!type.valid() || H5Tget_class(type.get()) != H5T_STRING) {
    throw std::runtime_error(std::string("snapshot inspection: invalid string attribute type ") + key);
  }
  if (H5Tis_variable_str(type.get()) > 0) {
    char* raw = nullptr;
    if (H5Aread(attr.get(), type.get(), &raw) < 0 || raw == nullptr) {
      throw std::runtime_error(std::string("snapshot inspection: failed reading variable string attribute ") + key);
    }
    const std::size_t length = std::char_traits<char>::length(raw);
    if (length > max_bytes) {
      H5free_memory(raw);
      throw std::length_error(std::string("snapshot inspection: oversized string attribute ") + key);
    }
    out.assign(raw, length);
    H5free_memory(raw);
    return true;
  }
  const std::size_t size = H5Tget_size(type.get());
  if (size == 0U || size > max_bytes) {
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
  if (H5Aexists(location, key) <= 0) return false;
  H5Handle attr(H5Aopen(location, key, H5P_DEFAULT));
  return attr.valid() && H5Aread(attr.get(), H5T_NATIVE_UINT32, &out) >= 0;
}

[[nodiscard]] bool readDoubleAttr(hid_t location, const char* key, double& out) {
  if (H5Aexists(location, key) <= 0) return false;
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



[[nodiscard]] std::optional<std::uint32_t> filenameMemberIndex(const std::filesystem::path& path) {
  const std::regex pattern(R"(^(.+)\.([0-9]+)\.hdf5$)");
  std::smatch match;
  const std::string filename = path.filename().string();
  if (!std::regex_match(filename, match, pattern)) return std::nullopt;
  const auto parsed = std::stoull(match[2].str());
  if (parsed > std::numeric_limits<std::uint32_t>::max()) return std::nullopt;
  return static_cast<std::uint32_t>(parsed);
}

[[nodiscard]] std::string filenameMemberStem(const std::filesystem::path& path) {
  const std::regex pattern(R"(^(.+)\.([0-9]+)\.hdf5$)");
  std::smatch match;
  const std::string filename = path.filename().string();
  if (std::regex_match(filename, match, pattern)) return match[1].str();
  return filename;
}

[[nodiscard]] MemberHeader inspectMember(
    const std::filesystem::path& path,
    const SnapshotReadOptions& options) {
  H5Handle file(H5Fopen(path.string().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
  if (!file.valid()) throw std::runtime_error("snapshot inspection: cannot open " + path.string());
  H5Handle header(H5Gopen2(file.get(), "/Header", H5P_DEFAULT));
  if (!header.valid()) throw std::runtime_error("snapshot inspection: missing /Header in " + path.string());

  MemberHeader result;
  result.path = path;
  std::error_code ec;
  result.file_size = std::filesystem::file_size(path, ec);
  if (ec) throw std::runtime_error("snapshot inspection: cannot stat " + path.string());
  const auto& names = sharedIoContractNames();
  static_cast<void>(readStringAttr(
      file.get(), std::string(names.file_kind_attribute).c_str(), result.file_kind,
      options.budget.max_attribute_bytes));

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
      if (high_width != 4U) throw std::runtime_error("snapshot inspection: NumPart_Total_HighWord must be 32-bit");
      for (std::size_t i = 0; i < result.global.size(); ++i) result.global[i] |= (high[i] << 32U);
    }
  }
  static_cast<void>(readU32Attr(header.get(), "NumFilesPerSnapshot", result.num_files));
  result.has_member_index = readU32Attr(header.get(), "CHUISnapshotMemberIndex", result.member_index);
  if (!result.has_member_index) {
    if (const auto parsed_index = filenameMemberIndex(path); parsed_index.has_value()) {
      result.member_index = *parsed_index;
      result.has_member_index = true;
    }
  }
  static_cast<void>(readDoubleAttr(header.get(), "Time", result.time));
  static_cast<void>(readDoubleAttr(header.get(), "Redshift", result.redshift));
  double scalar_box = std::numeric_limits<double>::quiet_NaN();
  if (readDoubleAttr(header.get(), "BoxSize", scalar_box)) {
    result.box_x = scalar_box;
    result.box_y = scalar_box;
    result.box_z = scalar_box;
  }
  static_cast<void>(readDoubleAttr(header.get(), "CHUIBoxSizeX_MpcComoving", result.box_x));
  static_cast<void>(readDoubleAttr(header.get(), "CHUIBoxSizeY_MpcComoving", result.box_y));
  static_cast<void>(readDoubleAttr(header.get(), "CHUIBoxSizeZ_MpcComoving", result.box_z));
  static_cast<void>(readDoubleAttr(header.get(), "Omega0", result.omega_matter));
  static_cast<void>(readDoubleAttr(header.get(), "OmegaLambda", result.omega_lambda));
  static_cast<void>(readDoubleAttr(header.get(), "OmegaBaryon", result.omega_baryon));
  static_cast<void>(readDoubleAttr(header.get(), "HubbleParam", result.hubble));
  static_cast<void>(readStringAttr(header.get(), "CosmoSimSchemaName", result.schema, options.budget.max_attribute_bytes));
  static_cast<void>(readU32Attr(header.get(), "CosmoSimSchemaVersion", result.schema_version));
  static_cast<void>(readStringAttr(header.get(), "CHUISnapshotGenerationID", result.generation, options.budget.max_attribute_bytes));
  static_cast<void>(readStringAttr(header.get(), "NamingRulesVersion", result.naming_rules_version, options.budget.max_attribute_bytes));
  static_cast<void>(readStringAttr(header.get(), "FileNamingRulesVersion", result.file_naming_rules_version, options.budget.max_attribute_bytes));
  static_cast<void>(readStringAttr(header.get(), "CHUIVelocityStorageConvention", result.velocity_storage_convention, options.budget.max_attribute_bytes));
  std::string stored_dialect;
  if (readStringAttr(header.get(), "CHUISnapshotDialect", stored_dialect, options.budget.max_attribute_bytes)) {
    result.dialect = parseDialect(stored_dialect);
  } else if (options.dialect != SnapshotDialect::kAuto) {
    result.dialect = options.dialect;
  }

  H5Handle units(
      H5Lexists(file.get(), "/Units", H5P_DEFAULT) > 0
          ? H5Gopen2(file.get(), "/Units", H5P_DEFAULT)
          : -1);
  if (units.valid()) {
    static_cast<void>(readStringAttr(units.get(), "LengthUnit", result.unit_length, options.budget.max_attribute_bytes));
    static_cast<void>(readStringAttr(units.get(), "MassUnit", result.unit_mass, options.budget.max_attribute_bytes));
    static_cast<void>(readStringAttr(units.get(), "VelocityUnit", result.unit_velocity, options.budget.max_attribute_bytes));
    static_cast<void>(readStringAttr(units.get(), "CoordinateFrame", result.coordinate_frame, options.budget.max_attribute_bytes));
  } else {
    static_cast<void>(readStringAttr(header.get(), "CHUIUnitLength", result.unit_length, options.budget.max_attribute_bytes));
    static_cast<void>(readStringAttr(header.get(), "CHUIUnitMass", result.unit_mass, options.budget.max_attribute_bytes));
    static_cast<void>(readStringAttr(header.get(), "CHUIUnitVelocity", result.unit_velocity, options.budget.max_attribute_bytes));
  }
  H5Handle provenance(
      H5Lexists(file.get(), "/Provenance", H5P_DEFAULT) > 0
          ? H5Gopen2(file.get(), "/Provenance", H5P_DEFAULT)
          : -1);
  if (provenance.valid()) {
    static_cast<void>(readStringAttr(provenance.get(), "config_hash_hex", result.config_hash_hex, options.budget.max_attribute_bytes));
  }
  if (result.num_files == 0U) throw std::runtime_error("snapshot inspection: NumFilesPerSnapshot is zero");
  if (result.has_member_index && result.member_index >= result.num_files) {
    throw std::runtime_error("snapshot inspection: member index is outside NumFilesPerSnapshot");
  }
  return result;
}
#endif

[[nodiscard]] bool sameScientificIdentity(const MemberHeader& lhs, const MemberHeader& rhs) {
  return lhs.num_files == rhs.num_files && lhs.global == rhs.global &&
      lhs.schema == rhs.schema && lhs.schema_version == rhs.schema_version &&
      lhs.generation == rhs.generation && lhs.file_kind == rhs.file_kind &&
      lhs.dialect == rhs.dialect && nearlyEqual(lhs.time, rhs.time) &&
      nearlyEqual(lhs.redshift, rhs.redshift) && nearlyEqual(lhs.box_x, rhs.box_x) &&
      nearlyEqual(lhs.box_y, rhs.box_y) && nearlyEqual(lhs.box_z, rhs.box_z) &&
      nearlyEqual(lhs.omega_matter, rhs.omega_matter) &&
      nearlyEqual(lhs.omega_lambda, rhs.omega_lambda) &&
      nearlyEqual(lhs.omega_baryon, rhs.omega_baryon) &&
      nearlyEqual(lhs.hubble, rhs.hubble) && lhs.unit_length == rhs.unit_length &&
      lhs.unit_mass == rhs.unit_mass && lhs.unit_velocity == rhs.unit_velocity &&
      lhs.coordinate_frame == rhs.coordinate_frame &&
      lhs.velocity_storage_convention == rhs.velocity_storage_convention &&
      lhs.config_hash_hex == rhs.config_hash_hex &&
      lhs.naming_rules_version == rhs.naming_rules_version &&
      lhs.file_naming_rules_version == rhs.file_naming_rules_version;
}

[[nodiscard]] std::filesystem::path setDirectory(const std::filesystem::path& input) {
  if (std::filesystem::is_directory(input)) return input;
  const auto parent = input.parent_path();
  return parent.empty() ? std::filesystem::path(".") : parent;
}

[[nodiscard]] std::vector<std::filesystem::path> manifestListedMembers(
    const std::filesystem::path& directory) {
  std::vector<std::filesystem::path> manifests;
  for (const auto& entry : std::filesystem::directory_iterator(directory)) {
    if (entry.is_regular_file() && entry.path().extension() == ".complete") manifests.push_back(entry.path());
  }
  if (manifests.empty()) return {};
  if (manifests.size() != 1U) {
    throw std::runtime_error("snapshot inspection: directory contains multiple completion manifests; select a member or a specific snapdir");
  }
  const std::string text = readBoundedTextFile(manifests.front());
  const auto values = parseKeyValueBody(text);
  const auto schema_it = values.find("schema");
  if (schema_it == values.end() || schema_it->second != "chui_snapshot_set_v2") return {};
  const std::uint32_t count = parseU32(values, "member_count");
  std::vector<std::filesystem::path> paths;
  paths.reserve(count);
  for (std::uint32_t i = 0; i < count; ++i) {
    const std::string filename = requireString(values, "member." + std::to_string(i) + ".filename");
    const std::filesystem::path relative(filename);
    if (relative.is_absolute() || relative.has_parent_path() || filename == "." || filename == "..") {
      throw std::runtime_error("snapshot set manifest contains unsafe member filename");
    }
    paths.push_back(directory / relative);
  }
  return paths;
}

[[nodiscard]] std::vector<MemberHeader> discoverMemberHeaders(
    const std::filesystem::path& input,
    const SnapshotReadOptions& options,
    bool prefer_completion_manifest) {
#if !COSMOSIM_ENABLE_HDF5
  static_cast<void>(input); static_cast<void>(options); static_cast<void>(prefer_completion_manifest);
  throw std::runtime_error("COSMOSIM_ENABLE_HDF5=OFF: snapshot inspection unavailable");
#else
  std::vector<std::filesystem::path> candidates;
  if (std::filesystem::is_directory(input)) {
    if (prefer_completion_manifest) candidates = manifestListedMembers(input);
    if (candidates.empty()) {
      for (const auto& entry : std::filesystem::directory_iterator(input)) {
        if (entry.is_regular_file() && entry.path().extension() == ".hdf5") candidates.push_back(entry.path());
      }
    }
  } else {
    const MemberHeader first = inspectMember(input, options);
    if (first.num_files <= 1U) {
      candidates.push_back(input);
    } else {
      const std::string stem = filenameMemberStem(input);
      const std::filesystem::path directory = setDirectory(input);
      for (const auto& entry : std::filesystem::directory_iterator(directory)) {
        if (entry.is_regular_file() && entry.path().extension() == ".hdf5" &&
            filenameMemberStem(entry.path()) == stem) {
          candidates.push_back(entry.path());
        }
      }
    }
  }
  if (candidates.empty()) throw std::runtime_error("snapshot inspection: no HDF5 members found");
  if (candidates.size() > options.budget.max_members) throw std::length_error("snapshot inspection: member-count budget exceeded");

  std::vector<MemberHeader> headers;
  headers.reserve(candidates.size());
  for (const auto& candidate : candidates) headers.push_back(inspectMember(candidate, options));

  if (std::filesystem::is_directory(input) && !prefer_completion_manifest) {
    std::set<std::string> generations;
    bool saw_empty_generation = false;
    for (const auto& header : headers) {
      if (header.generation.empty()) saw_empty_generation = true;
      else generations.insert(header.generation);
    }
    if (!generations.empty()) {
      if (generations.size() != 1U || saw_empty_generation) {
        throw std::runtime_error("snapshot inspection: directory mixes multiple CHUI generations or CHUI/external snapshots");
      }
    } else if (headers.size() > 1U) {
      std::set<std::string> stems;
      for (const auto& header : headers) stems.insert(filenameMemberStem(header.path));
      if (stems.size() != 1U) {
        throw std::runtime_error("snapshot inspection: directory contains multiple external snapshot stems; select a member explicitly");
      }
    }
  }

  std::sort(headers.begin(), headers.end(), [](const MemberHeader& lhs, const MemberHeader& rhs) {
    if (lhs.has_member_index && rhs.has_member_index && lhs.member_index != rhs.member_index) {
      return lhs.member_index < rhs.member_index;
    }
    return lhs.path.filename().string() < rhs.path.filename().string();
  });
  return headers;
#endif
}

[[nodiscard]] std::filesystem::path memberIntegrityPath(const std::filesystem::path& member_path) {
  return std::filesystem::path(member_path.string() + ".member.integrity");
}

struct MemberIntegrity {
  std::string generation_id;
  std::string filename;
  std::uint32_t member_index = 0;
  std::uint32_t num_files = 0;
  std::uint64_t file_size = 0;
  std::string sha256;
};

[[nodiscard]] MemberIntegrity readMemberIntegrity(const std::filesystem::path& member_path) {
  const auto values = parseKeyValueBody(readBoundedTextFile(memberIntegrityPath(member_path), 64U * 1024U));
  if (requireString(values, "schema") != "chui_snapshot_member_integrity_v1") {
    throw std::runtime_error("snapshot member integrity sidecar has unsupported schema");
  }
  MemberIntegrity out;
  out.generation_id = requireString(values, "generation_id");
  out.filename = requireString(values, "filename");
  out.member_index = parseU32(values, "member_index");
  out.num_files = parseU32(values, "num_files_per_snapshot");
  out.file_size = parseU64(values, "file_size");
  out.sha256 = requireString(values, "sha256");
  if (out.sha256.size() != 64U) throw std::runtime_error("snapshot member integrity sidecar has invalid SHA-256");
  return out;
}

[[nodiscard]] std::string buildCompletionManifestBody(
    const std::vector<MemberHeader>& headers,
    const std::vector<MemberIntegrity>& integrities) {
  const MemberHeader& ref = headers.front();
  std::string body;
  body += "schema=chui_snapshot_set_v2\n";
  body += "generation_id=" + ref.generation + "\n";
  body += "num_files_per_snapshot=" + std::to_string(ref.num_files) + "\n";
  body += "schema_name=" + ref.schema + "\n";
  body += "schema_version=" + std::to_string(ref.schema_version) + "\n";
  body += "file_kind=" + ref.file_kind + "\n";
  body += "dialect=" + dialectLabel(ref.dialect) + "\n";
  body += "global_part_count=" + formatCountArray(ref.global) + "\n";
  body += "time=" + formatDouble(ref.time) + "\n";
  body += "redshift=" + formatDouble(ref.redshift) + "\n";
  body += "box_size_x=" + formatDouble(ref.box_x) + "\n";
  body += "box_size_y=" + formatDouble(ref.box_y) + "\n";
  body += "box_size_z=" + formatDouble(ref.box_z) + "\n";
  body += "omega_matter=" + formatDouble(ref.omega_matter) + "\n";
  body += "omega_lambda=" + formatDouble(ref.omega_lambda) + "\n";
  body += "omega_baryon=" + formatDouble(ref.omega_baryon) + "\n";
  body += "hubble_param=" + formatDouble(ref.hubble) + "\n";
  body += "unit_length=" + ref.unit_length + "\n";
  body += "unit_mass=" + ref.unit_mass + "\n";
  body += "unit_velocity=" + ref.unit_velocity + "\n";
  body += "coordinate_frame=" + ref.coordinate_frame + "\n";
  body += "velocity_storage_convention=" + ref.velocity_storage_convention + "\n";
  body += "config_hash_hex=" + ref.config_hash_hex + "\n";
  body += "naming_rules_version=" + ref.naming_rules_version + "\n";
  body += "file_naming_rules_version=" + ref.file_naming_rules_version + "\n";
  body += "member_count=" + std::to_string(headers.size()) + "\n";
  for (std::size_t i = 0; i < headers.size(); ++i) {
    const std::string prefix = "member." + std::to_string(i) + ".";
    body += prefix + "filename=" + integrities[i].filename + "\n";
    body += prefix + "index=" + std::to_string(headers[i].member_index) + "\n";
    body += prefix + "file_size=" + std::to_string(integrities[i].file_size) + "\n";
    body += prefix + "sha256=" + integrities[i].sha256 + "\n";
    body += prefix + "local_part_count=" + formatCountArray(headers[i].local) + "\n";
  }
  return body;
}

[[nodiscard]] bool completionManifestMatches(
    const std::filesystem::path& input,
    const std::vector<MemberHeader>& headers,
    const SnapshotReadOptions& options) {
  if (headers.empty() || headers.front().generation.empty()) return false;
  const auto marker = setDirectory(input) / (headers.front().generation + ".complete");
  std::error_code ec;
  if (!std::filesystem::is_regular_file(marker, ec) || ec) return false;
  try {
    const std::string text = readBoundedTextFile(marker);
    const std::string digest_key = "set_digest_sha256=";
    const auto digest_pos = text.rfind(digest_key);
    if (digest_pos == std::string::npos) return false;
    const std::string body = text.substr(0, digest_pos);
    std::string digest_line = text.substr(digest_pos + digest_key.size());
    while (!digest_line.empty() && (digest_line.back() == '\n' || digest_line.back() == '\r')) digest_line.pop_back();
    if (digest_line.size() != 64U || core::internal::sha256Hex(body) != digest_line) return false;
    const auto values = parseKeyValueBody(text);
    const MemberHeader& ref = headers.front();
    if (requireString(values, "schema") != "chui_snapshot_set_v2" ||
        requireString(values, "generation_id") != ref.generation ||
        parseU32(values, "num_files_per_snapshot") != ref.num_files ||
        requireString(values, "schema_name") != ref.schema ||
        parseU32(values, "schema_version") != ref.schema_version ||
        requireString(values, "file_kind") != ref.file_kind ||
        parseDialect(requireString(values, "dialect")) != ref.dialect ||
        parseCountArray(requireString(values, "global_part_count")) != ref.global ||
        !nearlyEqual(parseDouble(values, "time"), ref.time) ||
        !nearlyEqual(parseDouble(values, "redshift"), ref.redshift) ||
        !nearlyEqual(parseDouble(values, "box_size_x"), ref.box_x) ||
        !nearlyEqual(parseDouble(values, "box_size_y"), ref.box_y) ||
        !nearlyEqual(parseDouble(values, "box_size_z"), ref.box_z) ||
        !nearlyEqual(parseDouble(values, "omega_matter"), ref.omega_matter) ||
        !nearlyEqual(parseDouble(values, "omega_lambda"), ref.omega_lambda) ||
        !nearlyEqual(parseDouble(values, "omega_baryon"), ref.omega_baryon) ||
        !nearlyEqual(parseDouble(values, "hubble_param"), ref.hubble) ||
        requireString(values, "unit_length") != ref.unit_length ||
        requireString(values, "unit_mass") != ref.unit_mass ||
        requireString(values, "unit_velocity") != ref.unit_velocity ||
        requireString(values, "coordinate_frame") != ref.coordinate_frame ||
        requireString(values, "velocity_storage_convention") != ref.velocity_storage_convention ||
        requireString(values, "config_hash_hex") != ref.config_hash_hex ||
        requireString(values, "naming_rules_version") != ref.naming_rules_version ||
        requireString(values, "file_naming_rules_version") != ref.file_naming_rules_version ||
        parseU32(values, "member_count") != headers.size()) {
      return false;
    }
    for (std::size_t i = 0; i < headers.size(); ++i) {
      const std::string prefix = "member." + std::to_string(i) + ".";
      const std::string filename = requireString(values, prefix + "filename");
      if (filename != headers[i].path.filename().string() ||
          parseU32(values, prefix + "index") != headers[i].member_index ||
          parseU64(values, prefix + "file_size") != headers[i].file_size ||
          parseCountArray(requireString(values, prefix + "local_part_count")) != headers[i].local) {
        return false;
      }
      const std::string expected_sha = requireString(values, prefix + "sha256");
      if (expected_sha.size() != 64U) return false;
      if (options.verify_snapshot_set_member_hashes &&
          internal::sha256FileHex(headers[i].path) != expected_sha) {
        return false;
      }
    }
    return true;
  } catch (const std::exception&) {
    return false;
  }
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


[[nodiscard]] std::uint64_t checkedCountSum(const std::array<std::uint64_t, 6>& counts) {
  std::uint64_t total = 0U;
  for (const auto count : counts) {
    if (count > std::numeric_limits<std::uint64_t>::max() - total) {
      throw std::overflow_error("snapshot set global particle-count sum overflow");
    }
    total += count;
  }
  return total;
}

[[nodiscard]] SnapshotReadResult readSet(
    const std::filesystem::path& input_path,
    const core::SimulationConfig& config,
    SnapshotReadOptions options,
    bool require_chui) {
  const SnapshotSetInspection inspection = inspectSnapshotSet(input_path, options);
  if (options.require_complete_chui_set && require_chui && !inspection.complete) {
    throw std::runtime_error("snapshot set is incomplete or does not match its CHUI completion manifest");
  }
  if (inspection.num_files_per_snapshot > options.budget.max_members) {
    throw std::length_error("snapshot set exceeds max_members read budget");
  }
  const std::uint64_t expected_particles = checkedCountSum(inspection.global_part_count);
  if (expected_particles > options.budget.max_particles) {
    throw std::length_error("snapshot set exceeds global particle-count read budget");
  }
  if (inspection.global_part_count[0] > options.budget.max_gas_cells) {
    throw std::length_error("snapshot set exceeds global gas-cell read budget");
  }
  std::uint64_t expected_sidecars = 0U;
  for (const std::size_t type_index : {0U, 3U, 4U, 5U}) {
    if (inspection.global_part_count[type_index] >
        std::numeric_limits<std::uint64_t>::max() - expected_sidecars) {
      throw std::overflow_error("snapshot set global sidecar-count sum overflow");
    }
    expected_sidecars += inspection.global_part_count[type_index];
  }
  if (expected_sidecars > options.budget.max_sidecar_rows) {
    throw std::length_error("snapshot set exceeds global sidecar-row read budget");
  }

  SnapshotReadResult merged;
  bool first = true;
  std::uint64_t materialized_bytes = 0U;
  std::uint64_t consumed_particles = 0U;
  std::uint64_t consumed_gas = 0U;
  std::uint64_t consumed_sidecars = 0U;
  for (std::size_t member_index = 0; member_index < inspection.member_paths.size(); ++member_index) {
    SnapshotReadOptions member_options = options;
    member_options.require_complete_chui_set = false;
    member_options.budget.max_particles = options.budget.max_particles - consumed_particles;
    member_options.budget.max_gas_cells = options.budget.max_gas_cells - consumed_gas;
    member_options.budget.max_sidecar_rows = options.budget.max_sidecar_rows - consumed_sidecars;
    member_options.budget.max_materialized_bytes = options.budget.max_materialized_bytes - materialized_bytes;
    SnapshotReadResult member = readGadgetArepoSnapshotHdf5(
        inspection.member_paths[member_index], config, member_options);
    if (require_chui && member.report.file_kind != sharedIoContractNames().science_snapshot_file_kind &&
        member.report.schema_name.rfind("gadget_arepo_v", 0) != 0) {
      throw std::runtime_error("readCosmoSimScienceSnapshotHdf5 rejected a non-CHUI snapshot member");
    }
    const std::uint64_t member_particles = checkedCountSum(member.report.local_part_count);
    const std::uint64_t member_gas = member.report.local_part_count[0];
    std::uint64_t member_sidecars = 0U;
    for (const std::size_t type_index : {0U, 3U, 4U, 5U}) {
      if (member.report.local_part_count[type_index] >
          std::numeric_limits<std::uint64_t>::max() - member_sidecars) {
        throw std::overflow_error("snapshot set member sidecar-count sum overflow");
      }
      member_sidecars += member.report.local_part_count[type_index];
    }
    if (member.report.materialized_state_bytes > options.budget.max_materialized_bytes - materialized_bytes) {
      throw std::length_error("snapshot set cumulative materialized-byte budget exceeded");
    }
    materialized_bytes += member.report.materialized_state_bytes;
    consumed_particles += member_particles;
    consumed_gas += member_gas;
    consumed_sidecars += member_sidecars;

    if (first) {
      merged.provenance = member.provenance;
      merged.normalized_config_text = member.normalized_config_text;
      merged.state.metadata = member.state.metadata;
      merged.report = member.report;
      first = false;
    } else {
      merged.report.defaulted_fields.insert(
          merged.report.defaulted_fields.end(), member.report.defaulted_fields.begin(), member.report.defaulted_fields.end());
      merged.report.unavailable_fields.insert(
          merged.report.unavailable_fields.end(), member.report.unavailable_fields.begin(), member.report.unavailable_fields.end());
    }
    appendState(merged.state, member.state, static_cast<std::uint32_t>(member_index));
  }
  if (consumed_particles != expected_particles || consumed_gas != inspection.global_part_count[0]) {
    throw std::runtime_error("snapshot set materialized counts disagree with inspected global counts");
  }
  const core::MemoryReport merged_memory = core::collectSimulationMemoryReport(merged.state);
  merged.report.materialized_state_bytes = merged_memory.totals.persistent_total_bytes;
  if (merged.report.materialized_state_bytes > options.budget.max_materialized_bytes) {
    throw std::length_error("merged snapshot state exceeds max_materialized_bytes");
  }
  merged.report.num_files_per_snapshot = inspection.num_files_per_snapshot;
  merged.report.local_part_count = inspection.global_part_count;
  merged.report.global_part_count = inspection.global_part_count;
  merged.report.generation_id = inspection.generation_id;
  if (!merged.state.validatePersistentParticleIds()) {
    throw std::runtime_error("snapshot set contains duplicate or zero persistent particle IDs across members");
  }
  internal::updateSnapshotReadiness(merged.state, &merged.report);
  return merged;
}

}  // namespace

namespace internal {

void writeSnapshotMemberIntegritySidecar(
    const std::filesystem::path& member_path,
    const SnapshotSetMemberInfo& member,
    bool durable_publication) {
  if (member.num_files_per_snapshot <= 1U) return;
  validateManifestAtom(member.generation_id, "generation id");
  if (member.member_index >= member.num_files_per_snapshot) {
    throw std::invalid_argument("snapshot member integrity sidecar received an invalid member index");
  }
  std::error_code ec;
  const auto file_size = std::filesystem::file_size(member_path, ec);
  if (ec) throw std::runtime_error("failed to stat published snapshot member for integrity sidecar");
  std::string body;
  body += "schema=chui_snapshot_member_integrity_v1\n";
  body += "generation_id=" + member.generation_id + "\n";
  body += "filename=" + member_path.filename().string() + "\n";
  body += "member_index=" + std::to_string(member.member_index) + "\n";
  body += "num_files_per_snapshot=" + std::to_string(member.num_files_per_snapshot) + "\n";
  body += "file_size=" + std::to_string(file_size) + "\n";
  body += "sha256=" + sha256FileHex(member_path) + "\n";
  writeTextFileTransactionally(
      memberIntegrityPath(member_path), body,
      durable_publication ? FileDurability::kDurablePublication : FileDurability::kAtomicVisibility);
}

}  // namespace internal

SnapshotSetInspection inspectSnapshotSet(
    const std::filesystem::path& input_path,
    const SnapshotReadOptions& options) {
#if !COSMOSIM_ENABLE_HDF5
  static_cast<void>(input_path); static_cast<void>(options);
  throw std::runtime_error("COSMOSIM_ENABLE_HDF5=OFF: snapshot inspection unavailable");
#else
  std::vector<MemberHeader> headers = discoverMemberHeaders(input_path, options, true);
  if (headers.empty()) throw std::runtime_error("snapshot inspection: no members found");
  const MemberHeader& reference = headers.front();
  if (reference.num_files > options.budget.max_members) {
    throw std::length_error("snapshot inspection: NumFilesPerSnapshot exceeds max_members");
  }
  SnapshotSetInspection result;
  result.dialect = reference.dialect;
  result.file_kind = reference.file_kind;
  result.schema_name = reference.schema;
  result.schema_version = reference.schema_version;
  result.num_files_per_snapshot = reference.num_files;
  result.global_part_count = reference.global;
  result.scale_factor = reference.time;
  result.redshift = reference.redshift;
  result.box_size_x = reference.box_x;
  result.box_size_y = reference.box_y;
  result.box_size_z = reference.box_z;
  result.omega_matter = reference.omega_matter;
  result.omega_lambda = reference.omega_lambda;
  result.omega_baryon = reference.omega_baryon;
  result.hubble_param = reference.hubble;
  result.unit_length = reference.unit_length;
  result.unit_mass = reference.unit_mass;
  result.unit_velocity = reference.unit_velocity;
  result.coordinate_frame = reference.coordinate_frame;
  result.velocity_storage_convention = reference.velocity_storage_convention;
  result.config_hash_hex = reference.config_hash_hex;
  result.naming_rules_version = reference.naming_rules_version;
  result.file_naming_rules_version = reference.file_naming_rules_version;
  result.generation_id = reference.generation;

  std::array<std::uint64_t, 6> summed{};
  std::set<std::uint32_t> member_indices;
  for (const auto& header : headers) {
    if (!sameScientificIdentity(header, reference)) {
      throw std::runtime_error("snapshot set members disagree on scientific identity metadata");
    }
    if (!header.has_member_index && reference.num_files > 1U) {
      throw std::runtime_error("snapshot set multifile member lacks an explicit or filename-derived member index");
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
    if (header.file_size > std::numeric_limits<std::uint64_t>::max() - result.total_member_file_bytes) {
      throw std::overflow_error("snapshot set member file-size sum overflow");
    }
    result.total_member_file_bytes += header.file_size;
    result.member_paths.push_back(header.path);
  }
  if (headers.size() > reference.num_files) {
    throw std::runtime_error("snapshot set contains more members than NumFilesPerSnapshot");
  }
  bool contiguous = headers.size() == reference.num_files;
  if (contiguous) {
    for (std::uint32_t index = 0; index < reference.num_files; ++index) {
      if (!member_indices.contains(index)) { contiguous = false; break; }
    }
  }
  const bool counts_match = summed == reference.global;
  const bool chui_multifile =
      reference.file_kind == sharedIoContractNames().science_snapshot_file_kind &&
      reference.num_files > 1U;
  const bool manifest_ok = !chui_multifile || completionManifestMatches(input_path, headers, options);
  result.complete = contiguous && counts_match && manifest_ok;
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
  validateManifestAtom(generation_id, "generation id");
  if (num_files_per_snapshot == 0U) {
    throw std::invalid_argument("snapshot completion manifest requires at least one member");
  }
  SnapshotReadOptions options;
  options.require_complete_chui_set = false;
  options.verify_snapshot_set_member_hashes = false;
  std::vector<MemberHeader> headers = discoverMemberHeaders(snapshot_directory, options, false);
  if (headers.size() != num_files_per_snapshot) {
    throw std::runtime_error("snapshot completion manifest cannot publish before every expected member exists");
  }
  const MemberHeader& reference = headers.front();
  if (reference.generation != generation_id || reference.global != global_part_count ||
      reference.num_files != num_files_per_snapshot ||
      reference.file_kind != sharedIoContractNames().science_snapshot_file_kind) {
    throw std::runtime_error("snapshot completion manifest arguments disagree with member headers");
  }
  std::set<std::uint32_t> indices;
  std::array<std::uint64_t, 6> local_sum{};
  std::vector<MemberIntegrity> integrities;
  integrities.reserve(headers.size());
  for (const auto& header : headers) {
    if (!sameScientificIdentity(header, reference)) {
      throw std::runtime_error("snapshot completion manifest refused scientifically inconsistent members");
    }
    if (!header.has_member_index || header.member_index >= num_files_per_snapshot ||
        !indices.insert(header.member_index).second) {
      throw std::runtime_error("snapshot completion manifest refused invalid member indices");
    }
    for (std::size_t i = 0; i < 6; ++i) {
      if (header.local[i] > std::numeric_limits<std::uint64_t>::max() - local_sum[i]) {
        throw std::overflow_error("snapshot completion local-count sum overflow");
      }
      local_sum[i] += header.local[i];
    }
    MemberIntegrity integrity = readMemberIntegrity(header.path);
    if (integrity.generation_id != generation_id ||
        integrity.filename != header.path.filename().string() ||
        integrity.member_index != header.member_index ||
        integrity.num_files != num_files_per_snapshot ||
        integrity.file_size != header.file_size) {
      throw std::runtime_error("snapshot member integrity sidecar disagrees with published member");
    }
    integrities.push_back(std::move(integrity));
  }
  for (std::uint32_t index = 0; index < num_files_per_snapshot; ++index) {
    if (!indices.contains(index)) throw std::runtime_error("snapshot completion manifest found a member-index gap");
  }
  if (local_sum != global_part_count) {
    throw std::runtime_error("snapshot completion manifest local-count sum disagrees with global counts");
  }
  const std::string body = buildCompletionManifestBody(headers, integrities);
  const std::string contents = body + "set_digest_sha256=" + core::internal::sha256Hex(body) + "\n";
  internal::writeTextFileTransactionally(
      snapshot_directory / (std::string(generation_id) + ".complete"), contents,
      durable_publication ? internal::FileDurability::kDurablePublication
                          : internal::FileDurability::kAtomicVisibility);
  for (const auto& header : headers) {
    std::error_code ignored;
    std::filesystem::remove(memberIntegrityPath(header.path), ignored);
  }
}

}  // namespace cosmosim::io
