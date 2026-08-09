#include "cosmosim/io/ic_reader.hpp"
#include "io/internal/ic_file_set_common.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <cstdlib>
#include <filesystem>
#include <functional>
#include <fstream>
#include <iterator>
#include <limits>
#include <numeric>
#include <optional>
#include <set>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "cosmosim/core/build_config.hpp"
#include "io/internal/ic_byte_codec.hpp"
#include "io/internal/ic_canonical_bundle.hpp"
#include "io/internal/ic_conversion_catalog.hpp"
#include "io/internal/ic_hdf5_handle.hpp"
#include "io/internal/ic_reader_session.hpp"
#include "io/internal/ic_record_codec.hpp"
#include "io/internal/ic_stream_ingestion.hpp"
#include "cosmosim/core/units.hpp"
#include "cosmosim/parallel/distributed_memory.hpp"

#if COSMOSIM_ENABLE_HDF5
#include <hdf5.h>
#endif
#if COSMOSIM_ENABLE_MPI
#include <mpi.h>
#endif

namespace cosmosim::io::file_set_internal {

[[nodiscard]] bool nearlyEqual(double lhs, double rhs) {
  return std::abs(lhs - rhs) <= 1.0e-10 * std::max({1.0, std::abs(lhs), std::abs(rhs)});
}

[[nodiscard]] bool missingFieldContractsEqual(
    const std::vector<IcMissingFieldContract>& lhs,
    const std::vector<IcMissingFieldContract>& rhs) {
  if (lhs.size() != rhs.size()) {
    return false;
  }
  for (const auto& expected : lhs) {
    const auto observed = std::find_if(
        rhs.begin(), rhs.end(), [&](const IcMissingFieldContract& actual) {
          return actual.source_file_index == expected.source_file_index &&
              actual.field_path == expected.field_path;
        });
    if (observed == rhs.end() || observed->policy != expected.policy ||
        !nearlyEqual(
            observed->configured_value_code, expected.configured_value_code) ||
        observed->resolution != expected.resolution) {
      return false;
    }
  }
  return true;
}

[[nodiscard]] IcSpeciesPolicy mapConfiguredPolicy(
    core::InitialConditionSpeciesPolicy policy,
    std::size_t type_index) {
  switch (policy) {
    case core::InitialConditionSpeciesPolicy::kReject:
      return IcSpeciesPolicy::kReject;
    case core::InitialConditionSpeciesPolicy::kDarkMatter:
      return type_index == 2U
          ? IcSpeciesPolicy::kCollisionlessFamily2AsDarkMatter
          : IcSpeciesPolicy::kCollisionlessFamily3AsDarkMatter;
    case core::InitialConditionSpeciesPolicy::kStar:
      return IcSpeciesPolicy::kStar;
    case core::InitialConditionSpeciesPolicy::kBlackHole:
      return IcSpeciesPolicy::kBlackHole;
    case core::InitialConditionSpeciesPolicy::kTracer:
      return IcSpeciesPolicy::kTracer;
  }
  throw std::invalid_argument("unknown configured IC species policy");
}

[[nodiscard]] IcMissingFieldPolicy mapConfiguredMissingFieldPolicy(
    core::InitialConditionMissingFieldPolicy policy) {
  switch (policy) {
    case core::InitialConditionMissingFieldPolicy::kReject:
      return IcMissingFieldPolicy::kReject;
    case core::InitialConditionMissingFieldPolicy::kReconstruct:
      return IcMissingFieldPolicy::kReconstruct;
    case core::InitialConditionMissingFieldPolicy::kUseConfigValue:
      return IcMissingFieldPolicy::kUseConfigValue;
    case core::InitialConditionMissingFieldPolicy::kDialectDefinedDefault:
      return IcMissingFieldPolicy::kDialectDefinedDefault;
    case core::InitialConditionMissingFieldPolicy::kPreserveUnavailable:
      return IcMissingFieldPolicy::kPreserveUnavailable;
  }
  throw std::invalid_argument("unknown configured IC missing-field policy");
}

[[nodiscard]] std::uint32_t speciesTag(IcSpeciesPolicy policy) {
  switch (policy) {
    case IcSpeciesPolicy::kGas:
      return static_cast<std::uint32_t>(core::ParticleSpecies::kGas);
    case IcSpeciesPolicy::kDarkMatter:
    case IcSpeciesPolicy::kCollisionlessFamily2AsDarkMatter:
    case IcSpeciesPolicy::kCollisionlessFamily3AsDarkMatter:
      return static_cast<std::uint32_t>(core::ParticleSpecies::kDarkMatter);
    case IcSpeciesPolicy::kStar:
      return static_cast<std::uint32_t>(core::ParticleSpecies::kStar);
    case IcSpeciesPolicy::kBlackHole:
      return static_cast<std::uint32_t>(core::ParticleSpecies::kBlackHole);
    case IcSpeciesPolicy::kTracer:
      return static_cast<std::uint32_t>(core::ParticleSpecies::kTracer);
    case IcSpeciesPolicy::kReject:
      break;
  }
  throw std::runtime_error("attempted to materialize a rejected IC family");
}

#if COSMOSIM_ENABLE_HDF5

using internal::Hdf5Handle;
using internal::IcReaderSession;
using internal::readChunkDouble;
using internal::readChunkU64;

[[nodiscard]] bool pathExists(hid_t parent, std::string_view path) {
  return H5Lexists(parent, std::string(path).c_str(), H5P_DEFAULT) > 0;
}
[[nodiscard]] bool attributeExists(hid_t parent, std::string_view name) {
  return H5Aexists(parent, std::string(name).c_str()) > 0;
}

struct TypeDescription {
  std::string name;
  IcScalarClass scalar_class = IcScalarClass::kFloatingPoint;
  std::uint8_t width = 0;
  bool is_signed = false;
  IcByteOrder order = IcByteOrder::kNotApplicable;
};

[[nodiscard]] TypeDescription describeType(hid_t type) {
  TypeDescription description;
  const H5T_class_t scalar_class = H5Tget_class(type);
  const std::size_t width = H5Tget_size(type);
  if (width == 0U || width > 255U) {
    throw std::runtime_error("unsupported HDF5 datatype width");
  }
  description.width = static_cast<std::uint8_t>(width);
  if (scalar_class == H5T_FLOAT) {
    description.scalar_class = IcScalarClass::kFloatingPoint;
    description.name = "float" + std::to_string(width * 8U);
    description.is_signed = true;
  } else if (scalar_class == H5T_INTEGER) {
    description.scalar_class = IcScalarClass::kInteger;
    const H5T_sign_t sign = H5Tget_sign(type);
    if (sign != H5T_SGN_NONE && sign != H5T_SGN_2) {
      throw std::runtime_error("unsupported HDF5 integer sign encoding");
    }
    description.is_signed = sign == H5T_SGN_2;
    description.name =
        std::string(description.is_signed ? "int" : "uint") +
        std::to_string(width * 8U);
  } else {
    throw std::runtime_error(
        "IC scalar data must use an integer or floating-point HDF5 type");
  }
  switch (H5Tget_order(type)) {
    case H5T_ORDER_LE:
      description.order = IcByteOrder::kLittleEndian;
      break;
    case H5T_ORDER_BE:
      description.order = IcByteOrder::kBigEndian;
      break;
    case H5T_ORDER_NONE:
      description.order = IcByteOrder::kNotApplicable;
      break;
    default:
      description.order = IcByteOrder::kNative;
      break;
  }
  return description;
}

[[nodiscard]] std::vector<std::uint64_t> attributeDimensions(hid_t attribute) {
  Hdf5Handle space(H5Aget_space(attribute));
  if (!space.valid()) {
    throw std::runtime_error("failed to inspect HDF5 attribute dataspace");
  }
  const int rank = H5Sget_simple_extent_ndims(space.get());
  if (rank < 0) {
    throw std::runtime_error("failed to inspect HDF5 attribute rank");
  }
  if (rank == 0) {
    return {};
  }
  std::vector<hsize_t> raw(static_cast<std::size_t>(rank));
  if (H5Sget_simple_extent_dims(space.get(), raw.data(), nullptr) < 0) {
    throw std::runtime_error("failed to inspect HDF5 attribute dimensions");
  }
  std::vector<std::uint64_t> dimensions;
  dimensions.reserve(raw.size());
  for (const hsize_t value : raw) {
    dimensions.push_back(static_cast<std::uint64_t>(value));
  }
  return dimensions;
}

[[nodiscard]] Hdf5Handle openValidatedAttribute(
    hid_t group,
    const char* name,
    bool required,
    IcScalarClass expected_class,
    std::span<const std::uint64_t> expected_dimensions,
    bool require_unsigned_integer = false,
    std::size_t max_integer_width = 4U) {
  Hdf5Handle attribute(H5Aopen(group, name, H5P_DEFAULT));
  if (!attribute.valid()) {
    if (required) {
      throw std::runtime_error(std::string("missing Header/") + name);
    }
    return attribute;
  }
  Hdf5Handle type(H5Aget_type(attribute.get()));
  if (!type.valid()) {
    throw std::runtime_error(std::string("failed to inspect Header/") + name);
  }
  const TypeDescription description = describeType(type.get());
  if (description.scalar_class != expected_class) {
    throw std::runtime_error(
        std::string("Header/") + name + " has an invalid datatype class");
  }
  if (require_unsigned_integer && description.is_signed) {
    throw std::runtime_error(
        std::string("Header/") + name + " must be unsigned integer data");
  }
  if (expected_class == IcScalarClass::kInteger &&
      description.width > max_integer_width) {
    throw std::runtime_error(
        std::string("Header/") + name +
        " exceeds the supported integer width");
  }
  if (expected_class == IcScalarClass::kFloatingPoint &&
      description.width != 4U && description.width != 8U) {
    throw std::runtime_error(
        std::string("Header/") + name + " must use float32 or float64");
  }
  const std::vector<std::uint64_t> dimensions =
      attributeDimensions(attribute.get());
  if (!std::equal(
          dimensions.begin(), dimensions.end(), expected_dimensions.begin(),
          expected_dimensions.end())) {
    std::ostringstream message;
    message << "Header/" << name << " has shape [";
    for (std::size_t i = 0; i < dimensions.size(); ++i) {
      message << (i == 0 ? "" : ",") << dimensions[i];
    }
    message << "] but expected [";
    for (std::size_t i = 0; i < expected_dimensions.size(); ++i) {
      message << (i == 0 ? "" : ",") << expected_dimensions[i];
    }
    message << ']';
    throw std::runtime_error(message.str());
  }
  return attribute;
}

void readAttributeNonnegativeU64x6(
    hid_t group,
    const char* name,
    std::array<std::uint64_t, 6>& values) {
  static constexpr std::array<std::uint64_t, 1> kExpected{6U};
  Hdf5Handle attribute = openValidatedAttribute(
      group, name, true, IcScalarClass::kInteger, kExpected, false,
      sizeof(std::uint64_t));
  Hdf5Handle type(H5Aget_type(attribute.get()));
  if (!type.valid()) {
    throw std::runtime_error(std::string("failed to inspect Header/") + name);
  }
  const TypeDescription description = describeType(type.get());
  if (description.is_signed) {
    std::array<std::int64_t, 6> signed_values{};
    if (H5Aread(attribute.get(), H5T_NATIVE_INT64, signed_values.data()) < 0) {
      throw std::runtime_error(std::string("failed to read Header/") + name);
    }
    for (std::size_t i = 0; i < signed_values.size(); ++i) {
      if (signed_values[i] < 0) {
        throw std::runtime_error(
            std::string("Header/") + name +
            " contains a negative particle count");
      }
      values[i] = static_cast<std::uint64_t>(signed_values[i]);
    }
    return;
  }
  if (H5Aread(attribute.get(), H5T_NATIVE_UINT64, values.data()) < 0) {
    throw std::runtime_error(std::string("failed to read Header/") + name);
  }
}

void readAttributeNonnegativeU32(
    hid_t group,
    const char* name,
    std::uint32_t& value) {
  Hdf5Handle attribute = openValidatedAttribute(
      group, name, true, IcScalarClass::kInteger, {}, false,
      sizeof(std::uint64_t));
  Hdf5Handle type(H5Aget_type(attribute.get()));
  if (!type.valid()) {
    throw std::runtime_error(std::string("failed to inspect Header/") + name);
  }
  const TypeDescription description = describeType(type.get());
  std::uint64_t wide = 0U;
  if (description.is_signed) {
    std::int64_t signed_value = 0;
    if (H5Aread(attribute.get(), H5T_NATIVE_INT64, &signed_value) < 0) {
      throw std::runtime_error(std::string("failed to read Header/") + name);
    }
    if (signed_value < 0) {
      throw std::runtime_error(
          std::string("Header/") + name + " must be non-negative");
    }
    wide = static_cast<std::uint64_t>(signed_value);
  } else if (H5Aread(attribute.get(), H5T_NATIVE_UINT64, &wide) < 0) {
    throw std::runtime_error(std::string("failed to read Header/") + name);
  }
  if (wide > std::numeric_limits<std::uint32_t>::max()) {
    throw std::runtime_error(
        std::string("Header/") + name + " exceeds uint32 range");
  }
  value = static_cast<std::uint32_t>(wide);
}

void readAttributeU32x6(
    hid_t group,
    const char* name,
    std::array<std::uint32_t, 6>& values,
    bool required = true) {
  static constexpr std::array<std::uint64_t, 1> kExpected{6U};
  Hdf5Handle attribute = openValidatedAttribute(
      group, name, required, IcScalarClass::kInteger, kExpected, true);
  if (!attribute.valid()) {
    return;
  }
  if (H5Aread(attribute.get(), H5T_NATIVE_UINT32, values.data()) < 0) {
    throw std::runtime_error(std::string("failed to read Header/") + name);
  }
}

void readAttributeF64x6(
    hid_t group,
    const char* name,
    std::array<double, 6>& values) {
  static constexpr std::array<std::uint64_t, 1> kExpected{6U};
  Hdf5Handle attribute = openValidatedAttribute(
      group, name, true, IcScalarClass::kFloatingPoint, kExpected);
  if (H5Aread(attribute.get(), H5T_NATIVE_DOUBLE, values.data()) < 0) {
    throw std::runtime_error(std::string("failed to read Header/") + name);
  }
}

void readAttributeF64(hid_t group, const char* name, double& value) {
  Hdf5Handle attribute = openValidatedAttribute(
      group, name, true, IcScalarClass::kFloatingPoint, {});
  if (H5Aread(attribute.get(), H5T_NATIVE_DOUBLE, &value) < 0) {
    throw std::runtime_error(std::string("failed to read Header/") + name);
  }
}

void readAttributeU32(hid_t group, const char* name, std::uint32_t& value) {
  Hdf5Handle attribute = openValidatedAttribute(
      group, name, true, IcScalarClass::kInteger, {}, true);
  if (H5Aread(attribute.get(), H5T_NATIVE_UINT32, &value) < 0) {
    throw std::runtime_error(std::string("failed to read Header/") + name);
  }
}

[[nodiscard]] std::string readAttributeString(
    hid_t group,
    const char* name) {
  Hdf5Handle attribute(H5Aopen(group, name, H5P_DEFAULT));
  if (!attribute.valid()) {
    throw std::runtime_error(std::string("missing Header/") + name);
  }
  if (!attributeDimensions(attribute.get()).empty()) {
    throw std::runtime_error(
        std::string("Header/") + name + " must be a scalar string");
  }
  Hdf5Handle type(H5Aget_type(attribute.get()));
  if (!type.valid() || H5Tget_class(type.get()) != H5T_STRING) {
    throw std::runtime_error(std::string("Header/") + name + " must be a string");
  }
  if (H5Tis_variable_str(type.get()) > 0) {
    char* raw = nullptr;
    if (H5Aread(attribute.get(), type.get(), &raw) < 0) {
      throw std::runtime_error(std::string("failed to read Header/") + name);
    }
    std::string value = raw == nullptr ? std::string{} : std::string(raw);
    if (raw != nullptr) {
      H5free_memory(raw);
    }
    return value;
  }
  const std::size_t width = H5Tget_size(type.get());
  if (width == 0U) {
    throw std::runtime_error(
        std::string("Header/") + name + " has zero-width string type");
  }
  std::vector<char> raw(width + 1U, '\0');
  if (H5Aread(attribute.get(), type.get(), raw.data()) < 0) {
    throw std::runtime_error(std::string("failed to read Header/") + name);
  }
  const auto end = std::find(
      raw.begin(), raw.begin() + static_cast<std::ptrdiff_t>(width), '\0');
  return std::string(raw.begin(), end);
}

[[nodiscard]] IcSchemaSummary readHeader(hid_t header) {
  IcSchemaSummary summary;
  std::array<std::uint64_t, 6> local{};
  std::array<std::uint32_t, 6> low{};
  readAttributeNonnegativeU64x6(header, "NumPart_ThisFile", local);
  readAttributeU32x6(header, "NumPart_Total", low);
  readAttributeU32x6(
      header, "NumPart_Total_HighWord", summary.total_count_high_word, false);
  readAttributeF64x6(header, "MassTable", summary.mass_table);
  readAttributeF64(header, "Time", summary.scale_factor);
  readAttributeF64(header, "Redshift", summary.redshift);
  readAttributeF64(header, "BoxSize", summary.box_size);
  readAttributeF64(header, "Omega0", summary.omega_matter);
  readAttributeF64(header, "OmegaLambda", summary.omega_lambda);
  readAttributeF64(header, "HubbleParam", summary.hubble_param);
  readAttributeNonnegativeU32(
      header, "NumFilesPerSnapshot", summary.num_files_per_snapshot);
  for (std::size_t i = 0; i < 6; ++i) {
    summary.count_by_type[i] = local[i];
    summary.total_count_by_type[i] =
        static_cast<std::uint64_t>(low[i]) |
        (static_cast<std::uint64_t>(summary.total_count_high_word[i]) << 32U);
  }
  return summary;
}

[[nodiscard]] std::string headerAuditText(const IcSchemaSummary& schema) {
  std::ostringstream output;
  output.precision(std::numeric_limits<double>::max_digits10);
  output << "NumFilesPerSnapshot=" << schema.num_files_per_snapshot
         << ";NumPart_ThisFile=";
  for (std::size_t i = 0; i < schema.count_by_type.size(); ++i) {
    output << (i == 0U ? "" : ",") << schema.count_by_type[i];
  }
  output << ";NumPart_Total=";
  for (std::size_t i = 0; i < schema.total_count_by_type.size(); ++i) {
    output << (i == 0U ? "" : ",") << schema.total_count_by_type[i];
  }
  output << ";NumPart_Total_HighWord=";
  for (std::size_t i = 0; i < schema.total_count_high_word.size(); ++i) {
    output << (i == 0U ? "" : ",") << schema.total_count_high_word[i];
  }
  output << ";MassTable=";
  for (std::size_t i = 0; i < schema.mass_table.size(); ++i) {
    output << (i == 0U ? "" : ",") << schema.mass_table[i];
  }
  output << ";BoxSize=" << schema.box_size
         << ";Time=" << schema.scale_factor
         << ";Redshift=" << schema.redshift
         << ";Omega0=" << schema.omega_matter
         << ";OmegaLambda=" << schema.omega_lambda
         << ";HubbleParam=" << schema.hubble_param;
  return output.str();
}

[[nodiscard]] std::vector<std::uint64_t> dataspaceDimensions(hid_t space) {
  const int rank=H5Sget_simple_extent_ndims(space);if(rank<0)throw std::runtime_error("failed to inspect HDF5 dataspace rank");if(rank==0)return {};std::vector<hsize_t> raw(static_cast<std::size_t>(rank));if(H5Sget_simple_extent_dims(space,raw.data(),nullptr)<0)throw std::runtime_error("failed to inspect HDF5 dimensions");std::vector<std::uint64_t> dims;dims.reserve(raw.size());for(hsize_t value:raw)dims.push_back(static_cast<std::uint64_t>(value));return dims;
}

herr_t collectLinkName(
    hid_t,
    const char* name,
    const H5L_info_t*,
    void* context) {
  auto* names = static_cast<std::vector<std::string>*>(context);
  names->emplace_back(name);
  return 0;
}

[[nodiscard]] std::vector<std::string> listLinkNames(hid_t group) {
  std::vector<std::string> names;
  hsize_t index = 0U;
  if (H5Literate(
          group, H5_INDEX_NAME, H5_ITER_NATIVE, &index, collectLinkName,
          &names) < 0) {
    throw std::runtime_error("failed to enumerate HDF5 group members");
  }
  return names;
}

[[nodiscard]] H5O_type_t objectType(hid_t group, std::string_view name) {
  H5O_info_t info{};
  if (H5Oget_info_by_name(
          group, std::string(name).c_str(), &info, H5O_INFO_BASIC, H5P_DEFAULT) < 0) {
    throw std::runtime_error(
        "failed to inspect HDF5 object " + std::string(name));
  }
  return info.type;
}

[[nodiscard]] std::string selectAlias(
    hid_t group,
    const std::vector<std::string>& aliases,
    bool required,
    std::string_view canonical) {
  std::vector<std::string> present;
  for (const auto& alias : aliases) {
    if (pathExists(group, alias)) {
      present.push_back(alias);
    }
  }
  if (present.size() > 1U) {
    std::ostringstream message;
    message << "ambiguous aliases for " << canonical << ": ";
    for (std::size_t i = 0; i < present.size(); ++i) {
      message << (i == 0U ? "" : ", ") << present[i];
    }
    throw std::runtime_error(message.str());
  }
  if (!present.empty()) {
    return present.front();
  }
  if (required) {
    throw std::runtime_error(
        "required dataset missing for " + std::string(canonical));
  }
  return {};
}

[[nodiscard]] std::vector<std::filesystem::path> discoverFiles(
    const std::filesystem::path& requested,
    std::uint32_t count) {
  if (count == 0U) {
    throw std::runtime_error("NumFilesPerSnapshot must be positive");
  }
  if (count == 1U) {
    return {requested.lexically_normal()};
  }

  const std::filesystem::path directory = requested.parent_path();
  const std::string filename = requested.filename().string();
  const std::size_t extension_position = filename.rfind('.');
  if (extension_position == std::string::npos) {
    throw std::runtime_error(
        "multifile IC path requires a .hdf5 or .h5 extension");
  }
  const std::string extension = filename.substr(extension_position);
  std::string prefix = filename.substr(0U, extension_position);
  const std::size_t index_separator = prefix.rfind('.');
  if (index_separator != std::string::npos) {
    const std::string maybe_index = prefix.substr(index_separator + 1U);
    if (!maybe_index.empty() &&
        std::all_of(
            maybe_index.begin(), maybe_index.end(),
            [](char value) { return value >= '0' && value <= '9'; })) {
      prefix = prefix.substr(0U, index_separator);
    }
  }

  std::vector<std::filesystem::path> files;
  files.reserve(count);
  for (std::uint32_t file_index = 0U; file_index < count; ++file_index) {
    auto candidate =
        (directory /
         (prefix + "." + std::to_string(file_index) + extension))
            .lexically_normal();
    if (!std::filesystem::is_regular_file(candidate)) {
      throw std::runtime_error(
          "missing multifile IC member: " + candidate.string());
    }
    files.push_back(std::move(candidate));
  }
  return files;
}

using CanonicalBundleVerification =
    internal::CanonicalBundleVerification;
using internal::verifyCanonicalBundle;

using Convention = internal::IcSourceConvention;
using FieldConversionContract = internal::IcFieldConversionContract;
using internal::fieldConversionContract;
using internal::mapBridgeFrame;
using internal::mapBridgeVelocityConvention;

[[nodiscard]] Convention conventionFor(
    IcDialect dialect,
    hid_t header,
    const core::SimulationConfig& config,
    bool supplied_manifest) {
  if (dialect == IcDialect::kChuiCanonicalV1) {
    const std::string schema_name =
        readAttributeString(header, "ChuiIcSchemaName");
    std::uint32_t schema_version = 0U;
    readAttributeU32(header, "ChuiIcSchemaVersion", schema_version);
    const std::string coordinate_frame =
        readAttributeString(header, "ChuiCoordinateFrame");
    const std::string velocity_convention =
        readAttributeString(header, "ChuiVelocityConvention");
    const std::string manifest_digest =
        readAttributeString(header, "ConversionManifestSha256");
    if (schema_name != "chui_canonical_v1" || schema_version != 2U ||
        coordinate_frame != "comoving" ||
        velocity_convention != "physical_peculiar") {
      throw std::runtime_error(
          "canonical CHUÍ IC header has an unsupported schema or convention");
    }
    if (manifest_digest.size() != 64U ||
        !std::all_of(
            manifest_digest.begin(), manifest_digest.end(), [](char c) {
              return (c >= '0' && c <= '9') || (c >= 'a' && c <= 'f');
            })) {
      throw std::runtime_error(
          "canonical CHUÍ IC ConversionManifestSha256 must be 64 lowercase hexadecimal characters");
    }
    double length = 0.0;
    double mass = 0.0;
    double velocity = 0.0;
    readAttributeF64(header, "ChuiLengthUnitToSI", length);
    readAttributeF64(header, "ChuiMassUnitToSI", mass);
    readAttributeF64(header, "ChuiVelocityUnitToSI", velocity);
    if (!std::isfinite(length) || !std::isfinite(mass) ||
        !std::isfinite(velocity) || length <= 0.0 || mass <= 0.0 ||
        velocity <= 0.0) {
      throw std::runtime_error(
          "canonical CHUÍ IC unit scales must be finite and positive");
    }
    core::UnitSystem units;
    units.length_si_per_code = length;
    units.mass_si_per_code = mass;
    units.velocity_si_per_code = velocity;
    return {
        .source_units = units,
        .frame = IcCoordinateFrame::kComoving,
        .velocity = IcVelocityConvention::kPhysicalPeculiar};
  }

  if (supplied_manifest) {
    // The supplied schema-v4 manifest is the authority for every field. These
    // neutral values are used only while re-inspecting HDF5 datatype/shape.
    core::UnitSystem units;
    units.length_si_per_code = 1.0;
    units.mass_si_per_code = 1.0;
    units.velocity_si_per_code = 1.0;
    return {.source_units = units};
  }

  core::UnitSystem units;
  units.length_si_per_code =
      config.mode.ic_bridge_source_length_unit_to_si;
  units.mass_si_per_code = config.mode.ic_bridge_source_mass_unit_to_si;
  units.velocity_si_per_code =
      config.mode.ic_bridge_source_velocity_unit_to_si;
  if (!(units.length_si_per_code > 0.0) ||
      !(units.mass_si_per_code > 0.0) ||
      !(units.velocity_si_per_code > 0.0) ||
      !std::isfinite(units.length_si_per_code) ||
      !std::isfinite(units.mass_si_per_code) ||
      !std::isfinite(units.velocity_si_per_code)) {
    throw std::invalid_argument(
        "gadget_arepo_bridge_v1 requires positive explicit source unit scales");
  }
  return {
      .source_units = units,
      .frame = mapBridgeFrame(config.mode.ic_bridge_coordinate_frame),
      .velocity = mapBridgeVelocityConvention(
          config.mode.ic_bridge_velocity_convention),
      .length_hubble_exponent =
          config.mode.ic_bridge_length_hubble_exponent,
      .length_scale_factor_exponent =
          config.mode.ic_bridge_length_scale_factor_exponent,
      .mass_hubble_exponent = config.mode.ic_bridge_mass_hubble_exponent,
      .mass_scale_factor_exponent =
          config.mode.ic_bridge_mass_scale_factor_exponent,
      .velocity_hubble_exponent =
          config.mode.ic_bridge_velocity_hubble_exponent,
      .velocity_scale_factor_exponent =
          config.mode.ic_bridge_velocity_scale_factor_exponent};
}

void validateDatasetSemanticType(
    const TypeDescription& description,
    IcFieldSemantics semantics,
    std::string_view canonical_path) {
  if (semantics == IcFieldSemantics::kIdentifier) {
    if (description.scalar_class != IcScalarClass::kInteger ||
        description.is_signed || description.width > sizeof(std::uint64_t)) {
      throw std::runtime_error(
          std::string(canonical_path) +
          " must use an unsigned integer datatype no wider than uint64");
    }
    return;
  }
  if (description.scalar_class != IcScalarClass::kFloatingPoint ||
      (description.width != 4U && description.width != 8U)) {
    throw std::runtime_error(
        std::string(canonical_path) +
        " must use float32 or float64 physical data");
  }
}

[[nodiscard]] IcFieldManifest inspectDataset(
    hid_t group,
    std::uint32_t file_index,
    std::string canonical_path,
    std::string selected_alias,
    const Convention& convention,
    IcFieldSemantics semantics,
    IcVelocityConvention velocity = IcVelocityConvention::kNotVelocity) {
  Hdf5Handle dataset(
      H5Dopen2(group, selected_alias.c_str(), H5P_DEFAULT));
  if (!dataset.valid()) {
    throw std::runtime_error("failed to open dataset " + canonical_path);
  }
  Hdf5Handle type(H5Dget_type(dataset.get()));
  Hdf5Handle space(H5Dget_space(dataset.get()));
  if (!type.valid() || !space.valid()) {
    throw std::runtime_error("failed to inspect dataset " + canonical_path);
  }
  const TypeDescription type_description = describeType(type.get());
  validateDatasetSemanticType(type_description, semantics, canonical_path);
  const auto dimensions = dataspaceDimensions(space.get());
  const FieldConversionContract contract =
      fieldConversionContract(canonical_path, semantics, convention);
  return {
      .source_file_index = file_index,
      .dataset_path = std::move(canonical_path),
      .selected_alias = std::move(selected_alias),
      .scalar_type = type_description.name,
      .scalar_class = type_description.scalar_class,
      .byte_width = type_description.width,
      .is_signed = type_description.is_signed,
      .byte_order = type_description.order,
      .rank = static_cast<std::uint8_t>(dimensions.size()),
      .dimensions = dimensions,
      .record_count = dimensions.empty() ? 1U : dimensions.front(),
      .base_unit_to_si = contract.base_unit_to_si,
      .hubble_exponent = contract.hubble_exponent,
      .scale_factor_exponent = contract.scale_factor_exponent,
      .length_power = contract.length_power,
      .mass_power = contract.mass_power,
      .time_power = contract.time_power,
      .frame_scale_factor_exponent =
          contract.frame_scale_factor_exponent,
      .velocity_convention_power = contract.velocity_convention_power,
      .coordinate_frame = convention.frame,
      .velocity_convention = contract.velocity_convention_power > 0U
          ? convention.velocity
          : velocity,
      .semantics = semantics,
      .disposition = IcFieldDisposition::kConverted,
      .source_unit = contract.source_unit,
      .target_unit = contract.target_unit,
      .conversion_equation =
          "target = stored * base_unit_to_si * h^hubble_exponent * "
          "a^(scale_factor_exponent + frame_scale_factor_exponent) * "
          "velocity_convention_multiplier^velocity_convention_power / "
          "target_si_per_code(L^length_power M^mass_power T^time_power)"};
}

[[nodiscard]] IcFieldManifest inspectDroppedDataset(
    hid_t group,
    std::uint32_t file_index,
    std::string dataset_path,
    std::string selected_alias) {
  Hdf5Handle dataset(
      H5Dopen2(group, selected_alias.c_str(), H5P_DEFAULT));
  if (!dataset.valid()) {
    throw std::runtime_error("failed to open dataset " + dataset_path);
  }
  Hdf5Handle type(H5Dget_type(dataset.get()));
  Hdf5Handle space(H5Dget_space(dataset.get()));
  if (!type.valid() || !space.valid()) {
    throw std::runtime_error("failed to inspect dataset " + dataset_path);
  }
  const TypeDescription description = describeType(type.get());
  if (description.scalar_class == IcScalarClass::kFloatingPoint &&
      description.width != 4U && description.width != 8U) {
    throw std::runtime_error(
        dataset_path + " uses an unsupported floating-point width");
  }
  if (description.scalar_class == IcScalarClass::kInteger &&
      description.width > sizeof(std::uint64_t)) {
    throw std::runtime_error(dataset_path + " uses an integer wider than uint64");
  }
  const std::vector<std::uint64_t> dimensions =
      dataspaceDimensions(space.get());
  return {
      .source_file_index = file_index,
      .dataset_path = std::move(dataset_path),
      .selected_alias = std::move(selected_alias),
      .scalar_type = description.name,
      .scalar_class = description.scalar_class,
      .byte_width = description.width,
      .is_signed = description.is_signed,
      .byte_order = description.order,
      .rank = static_cast<std::uint8_t>(dimensions.size()),
      .dimensions = dimensions,
      .record_count = dimensions.empty() ? 1U : dimensions.front(),
      .base_unit_to_si = 1.0,
      .hubble_exponent = 0.0,
      .scale_factor_exponent = 0.0,
      .length_power = 0,
      .mass_power = 0,
      .time_power = 0,
      .frame_scale_factor_exponent = 0.0,
      .velocity_convention_power = 0U,
      .coordinate_frame = IcCoordinateFrame::kComoving,
      .velocity_convention = IcVelocityConvention::kNotVelocity,
      .semantics = IcFieldSemantics::kIntensive,
      .disposition = IcFieldDisposition::kDropped,
      .source_unit = "unknown",
      .target_unit = "not_imported",
      .conversion_equation = "not converted: explicitly dropped"};
}

[[nodiscard]] IcFieldManifest inspectHeaderAttribute(
    hid_t header,
    std::uint32_t file_index,
    std::string name,
    const Convention& convention,
    IcFieldSemantics semantics) {
  Hdf5Handle attribute(H5Aopen(header, name.c_str(), H5P_DEFAULT));
  if (!attribute.valid()) {
    throw std::runtime_error("failed to open Header/" + name);
  }
  Hdf5Handle type(H5Aget_type(attribute.get()));
  Hdf5Handle space(H5Aget_space(attribute.get()));
  if (!type.valid() || !space.valid()) {
    throw std::runtime_error("failed to inspect Header/" + name);
  }
  const TypeDescription type_description = describeType(type.get());
  validateDatasetSemanticType(
      type_description, semantics, "/Header/" + name);
  const auto dimensions = dataspaceDimensions(space.get());
  const FieldConversionContract contract =
      fieldConversionContract("/Header/" + name, semantics, convention);
  return {
      .source_file_index = file_index,
      .dataset_path = "/Header/" + name,
      .selected_alias = name,
      .scalar_type = type_description.name,
      .scalar_class = type_description.scalar_class,
      .byte_width = type_description.width,
      .is_signed = type_description.is_signed,
      .byte_order = type_description.order,
      .rank = static_cast<std::uint8_t>(dimensions.size()),
      .dimensions = dimensions,
      .record_count = dimensions.empty() ? 1U : dimensions.front(),
      .base_unit_to_si = contract.base_unit_to_si,
      .hubble_exponent = contract.hubble_exponent,
      .scale_factor_exponent = contract.scale_factor_exponent,
      .length_power = contract.length_power,
      .mass_power = contract.mass_power,
      .time_power = contract.time_power,
      .frame_scale_factor_exponent =
          contract.frame_scale_factor_exponent,
      .velocity_convention_power = contract.velocity_convention_power,
      .coordinate_frame = convention.frame,
      .velocity_convention = IcVelocityConvention::kNotVelocity,
      .semantics = semantics,
      .disposition = IcFieldDisposition::kConverted,
      .source_unit = contract.source_unit,
      .target_unit = contract.target_unit,
      .conversion_equation =
          "target = stored * base_unit_to_si * h^hubble_exponent * "
          "a^(scale_factor_exponent + frame_scale_factor_exponent) / "
          "target_si_per_code(L^length_power M^mass_power T^time_power)"};
}

[[nodiscard]] IcFieldManifest inspectHeaderIntegerAttribute(
    hid_t header,
    std::uint32_t file_index,
    std::string name,
    bool required = true) {
  Hdf5Handle attribute(H5Aopen(header, name.c_str(), H5P_DEFAULT));
  if (!attribute.valid()) {
    if (!required) {
      return {};
    }
    throw std::runtime_error("failed to open Header/" + name);
  }
  Hdf5Handle type(H5Aget_type(attribute.get()));
  Hdf5Handle space(H5Aget_space(attribute.get()));
  if (!type.valid() || !space.valid()) {
    throw std::runtime_error("failed to inspect Header/" + name);
  }
  const TypeDescription description = describeType(type.get());
  if (description.scalar_class != IcScalarClass::kInteger ||
      description.width > sizeof(std::uint64_t)) {
    throw std::runtime_error(
        "Header/" + name + " must use an integer type no wider than 64 bits");
  }
  const auto dimensions = dataspaceDimensions(space.get());
  return {
      .source_file_index = file_index,
      .dataset_path = "/Header/" + name,
      .selected_alias = name,
      .scalar_type = description.name,
      .scalar_class = description.scalar_class,
      .byte_width = description.width,
      .is_signed = description.is_signed,
      .byte_order = description.order,
      .rank = static_cast<std::uint8_t>(dimensions.size()),
      .dimensions = dimensions,
      .record_count = dimensions.empty() ? 1U : dimensions.front(),
      .base_unit_to_si = 1.0,
      .hubble_exponent = 0.0,
      .scale_factor_exponent = 0.0,
      .length_power = 0,
      .mass_power = 0,
      .time_power = 0,
      .frame_scale_factor_exponent = 0.0,
      .velocity_convention_power = 0U,
      .coordinate_frame = IcCoordinateFrame::kComoving,
      .velocity_convention = IcVelocityConvention::kNotVelocity,
      .semantics = IcFieldSemantics::kIdentifier,
      .disposition = IcFieldDisposition::kPreserved,
      .source_unit = "integer_count",
      .target_unit = "manifest_metadata",
      .conversion_equation = "preserved exactly after checked non-negative read"};
}


void validateCrossFileSchema(const IcManifest& manifest) {
  if (manifest.num_files_per_snapshot <= 1U) {
    return;
  }
  const auto comparable_dimensions = [](const IcFieldManifest& lhs,
                                        const IcFieldManifest& rhs) {
    if (lhs.rank != rhs.rank || lhs.dimensions.size() != rhs.dimensions.size()) {
      return false;
    }
    if (lhs.dataset_path.starts_with("/Header/")) {
      return lhs.dimensions == rhs.dimensions;
    }
    return std::equal(
        lhs.dimensions.begin() + (lhs.dimensions.empty() ? 0 : 1),
        lhs.dimensions.end(),
        rhs.dimensions.begin() + (rhs.dimensions.empty() ? 0 : 1));
  };
  const auto is_flexible_count_attribute = [](std::string_view path) {
    return path == "/Header/NumPart_ThisFile" ||
        path == "/Header/NumFilesPerSnapshot";
  };
  for (const IcFieldManifest& baseline : manifest.fields) {
    if (baseline.source_file_index != 0U) {
      continue;
    }
    for (std::uint32_t file_index = 1U;
         file_index < manifest.num_files_per_snapshot;
         ++file_index) {
      const auto candidate = std::find_if(
          manifest.fields.begin(), manifest.fields.end(),
          [&](const IcFieldManifest& field) {
            return field.source_file_index == file_index &&
                   field.dataset_path == baseline.dataset_path;
          });
      if (candidate == manifest.fields.end()) {
        throw std::runtime_error(
            "inconsistent source schema across IC files: missing " +
            baseline.dataset_path + " in file index " +
            std::to_string(file_index));
      }
      const bool flexible_count =
          is_flexible_count_attribute(baseline.dataset_path);
      const bool incompatible_type = flexible_count
          ? candidate->scalar_class != IcScalarClass::kInteger ||
              candidate->byte_width > sizeof(std::uint64_t)
          : candidate->scalar_type != baseline.scalar_type ||
              candidate->scalar_class != baseline.scalar_class ||
              candidate->byte_width != baseline.byte_width ||
              candidate->is_signed != baseline.is_signed ||
              candidate->byte_order != baseline.byte_order;
      if (candidate->selected_alias != baseline.selected_alias ||
          incompatible_type || !comparable_dimensions(baseline, *candidate)) {
        throw std::runtime_error(
            "inconsistent source schema across IC files for " +
            baseline.dataset_path);
      }
    }
  }
  for (const IcFieldManifest& field : manifest.fields) {
    if (field.source_file_index == 0U) {
      continue;
    }
    const auto baseline = std::find_if(
        manifest.fields.begin(), manifest.fields.end(),
        [&](const IcFieldManifest& candidate) {
          return candidate.source_file_index == 0U &&
                 candidate.dataset_path == field.dataset_path;
        });
    if (baseline == manifest.fields.end()) {
      throw std::runtime_error(
          "inconsistent source schema across IC files: unexpected " +
          field.dataset_path + " in file index " +
          std::to_string(field.source_file_index));
    }
  }
}

void checkedCounterAdd(
    std::uint64_t& destination,
    std::uint64_t value,
    std::string_view counter_name);

[[nodiscard]] bool isSupportedPartTypeName(
    std::string_view name,
    std::size_t& type_index) {
  static constexpr std::string_view kPrefix = "PartType";
  if (!name.starts_with(kPrefix)) {
    return false;
  }
  const std::string_view suffix = name.substr(kPrefix.size());
  if (suffix.empty() ||
      !std::all_of(suffix.begin(), suffix.end(), [](char value) {
        return value >= '0' && value <= '9';
      })) {
    return false;
  }
  std::uint64_t parsed = 0U;
  for (char value : suffix) {
    parsed = parsed * 10U + static_cast<std::uint64_t>(value - '0');
    if (parsed > std::numeric_limits<std::size_t>::max()) {
      throw std::overflow_error("PartType index overflow");
    }
  }
  type_index = static_cast<std::size_t>(parsed);
  return type_index < kParticleTypeCount;
}

[[nodiscard]] std::uint64_t logicalHeaderPayloadBytes(
    const IcSchemaSummary& schema) {
  static_cast<void>(schema);
  return 3U * 6U * sizeof(std::uint32_t) +
      6U * sizeof(double) + 6U * sizeof(double) + sizeof(std::uint32_t);
}

void recordRootSchemaDisposition(
    hid_t file,
    const IcSchemaSummary& schema,
    IcManifest& manifest) {
  for (const std::string& name : listLinkNames(file)) {
    if (name == "Header") {
      continue;
    }
    std::size_t type_index = 0U;
    if (isSupportedPartTypeName(name, type_index)) {
      if (schema.count_by_type[type_index] == 0U) {
        if (objectType(file, name) != H5O_TYPE_GROUP) {
          throw std::runtime_error("/" + name + " must be an HDF5 group");
        }
        Hdf5Handle group(H5Gopen2(file, name.c_str(), H5P_DEFAULT));
        if (!group.valid()) {
          throw std::runtime_error("failed to inspect /" + name);
        }
        if (!listLinkNames(group.get()).empty()) {
          throw std::runtime_error(
              "/" + name +
              " contains datasets although NumPart_ThisFile is zero");
        }
      }
      continue;
    }
    if (name.starts_with("PartType")) {
      throw std::runtime_error(
          "unsupported populated particle-family object /" + name);
    }
    const std::string disposition =
        "/" + name + ": auxiliary root object is not imported";
    manifest.preserved_auxiliary_fields.push_back(disposition);
    manifest.warnings.push_back(disposition);
  }
}

[[nodiscard]] SourceFileInspection inspectOneSourceFile(
    const std::filesystem::path& path,
    std::uint32_t file_index,
    IcDialect dialect,
    const std::array<IcSpeciesPolicy, kParticleTypeCount>& species_policy,
    const core::SimulationConfig& config,
    const IcImportOptions& options,
    bool has_authoritative_manifest) {
  if (!std::filesystem::is_regular_file(path)) {
    throw std::runtime_error(
        "IC source is not a regular file: " + path.string());
  }
  Hdf5Handle file(
      H5Fopen(path.string().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
  if (!file.valid()) {
    throw std::runtime_error("failed to open IC member: " + path.string());
  }
  Hdf5Handle header(H5Gopen2(file.get(), "/Header", H5P_DEFAULT));
  if (!header.valid()) {
    throw std::runtime_error("IC member missing /Header: " + path.string());
  }

  SourceFileInspection result;
  result.path = path;
  result.schema = readHeader(header.get());
  const bool source_declares_canonical =
      attributeExists(header.get(), "ChuiIcSchemaName");
  if (dialect == IcDialect::kGadgetArepoBridgeV1 &&
      source_declares_canonical) {
    throw std::runtime_error(
        "IC source declares chui_canonical_v1 but the authoritative dialect "
        "is gadget_arepo_bridge_v1");
  }
  if (dialect == IcDialect::kChuiCanonicalV1) {
    const CanonicalBundleVerification verification =
        verifyCanonicalBundle(path, file.get(), header.get());
    result.canonical_manifest_verified = verification.verified;
    result.canonical_manifest_sha256 = verification.manifest_sha256;
  }
  result.original_header_attributes = headerAuditText(result.schema);
  checkedCounterAdd(
      result.counters.logical_metadata_bytes_read,
      logicalHeaderPayloadBytes(result.schema), "logical_metadata_bytes_read");
  result.source_size_bytes = std::filesystem::file_size(path);
  result.source_sha256 = icSha256FileHex(path);
  checkedCounterAdd(
      result.counters.full_file_hash_pass_count, 1U,
      "full_file_hash_pass_count");
  checkedCounterAdd(
      result.counters.hash_bytes_read, result.source_size_bytes,
      "hash_bytes_read");

  IcManifest dispositions;
  recordRootSchemaDisposition(file.get(), result.schema, dispositions);
  result.preserved_auxiliary_fields =
      std::move(dispositions.preserved_auxiliary_fields);
  result.warnings = std::move(dispositions.warnings);

  const Convention convention = conventionFor(
      dialect, header.get(), config, has_authoritative_manifest);
  result.fields.push_back(inspectHeaderIntegerAttribute(
      header.get(), file_index, "NumPart_ThisFile"));
  result.fields.push_back(inspectHeaderIntegerAttribute(
      header.get(), file_index, "NumPart_Total"));
  if (attributeExists(header.get(), "NumPart_Total_HighWord")) {
    result.fields.push_back(inspectHeaderIntegerAttribute(
        header.get(), file_index, "NumPart_Total_HighWord"));
  }
  result.fields.push_back(inspectHeaderIntegerAttribute(
      header.get(), file_index, "NumFilesPerSnapshot"));
  result.fields.push_back(inspectHeaderAttribute(
      header.get(), file_index, "MassTable", convention,
      IcFieldSemantics::kExtensive));
  result.fields.push_back(inspectHeaderAttribute(
      header.get(), file_index, "BoxSize", convention,
      IcFieldSemantics::kCoordinate));

  for (std::size_t type = 0; type < kParticleTypeCount; ++type) {
    const std::uint64_t count = result.schema.count_by_type[type];
    if (count == 0U) {
      continue;
    }
    if (species_policy[type] == IcSpeciesPolicy::kReject) {
      throw std::invalid_argument(
          "populated PartType" + std::to_string(type) +
          " has explicit reject policy");
    }
    const std::string group_path = "/PartType" + std::to_string(type);
    Hdf5Handle group(
        H5Gopen2(file.get(), group_path.c_str(), H5P_DEFAULT));
    if (!group.valid()) {
      throw std::runtime_error(
          "missing " + group_path + " in " + path.string());
    }
    std::set<std::string> handled_aliases;
    const auto add = [&](
                         std::string canonical,
                         const std::vector<std::string>& aliases,
                         bool required,
                         IcFieldSemantics semantics,
                         IcVelocityConvention velocity =
                             IcVelocityConvention::kNotVelocity,
                         IcFieldDisposition disposition =
                             IcFieldDisposition::kConverted) {
      std::string selected;
      if (has_authoritative_manifest) {
        const std::string canonical_path = group_path + "/" + canonical;
        const auto expected = std::find_if(
            options.manifest->fields.begin(),
            options.manifest->fields.end(),
            [&](const IcFieldManifest& field) {
              return field.source_file_index == file_index &&
                  field.dataset_path == canonical_path &&
                  field.disposition != IcFieldDisposition::kDropped;
            });
        if (expected == options.manifest->fields.end()) {
          if (required) {
            throw std::runtime_error(
                "supplied IC manifest lacks required field contract " +
                canonical_path);
          }
          return false;
        }
        selected = expected->selected_alias;
        if (selected.empty() || !pathExists(group.get(), selected)) {
          throw std::runtime_error(
              "supplied IC manifest alias does not exist for " +
              canonical_path + ": " + selected);
        }
      } else {
        selected = selectAlias(
            group.get(), aliases, required, group_path + "/" + canonical);
        if (selected.empty()) {
          return false;
        }
      }
      handled_aliases.insert(selected);
      IcFieldManifest field = inspectDataset(
          group.get(), file_index, group_path + "/" + canonical, selected,
          convention, semantics, velocity);
      if (field.record_count != count) {
        throw std::runtime_error(
            "dataset record count disagrees with NumPart_ThisFile for " +
            field.dataset_path);
      }
      if ((canonical == "Coordinates" || canonical == "Velocities") &&
          (field.rank != 2U || field.dimensions.size() != 2U ||
           field.dimensions[1] != 3U)) {
        throw std::runtime_error(
            field.dataset_path + " must have dimensions [N,3]");
      }
      if (canonical != "Coordinates" && canonical != "Velocities" &&
          field.rank != 1U) {
        throw std::runtime_error(field.dataset_path + " must have rank 1");
      }
      field.disposition = disposition;
      result.fields.push_back(std::move(field));
      return true;
    };
    const auto resolve_missing = [
        &result, &options, file_index, has_authoritative_manifest](
                                     std::string field_path,
                                     core::InitialConditionMissingFieldPolicy configured_policy,
                                     double configured_value_code,
                                     std::string dialect_resolution) {
      if (has_authoritative_manifest) {
        const auto supplied = std::find_if(
            options.manifest->missing_field_contracts.begin(),
            options.manifest->missing_field_contracts.end(),
            [&](const IcMissingFieldContract& contract) {
              return contract.source_file_index == file_index &&
                  contract.field_path == field_path;
            });
        if (supplied ==
            options.manifest->missing_field_contracts.end()) {
          throw std::runtime_error(
              "supplied IC manifest lacks authoritative missing-field "
              "contract " + field_path + " for file " +
              std::to_string(file_index));
        }
        result.missing_field_contracts.push_back(*supplied);
        result.defaulted_fields.push_back(
            field_path + "=" + supplied->resolution);
        return;
      }

      const IcMissingFieldPolicy policy =
          mapConfiguredMissingFieldPolicy(configured_policy);
      if (policy == IcMissingFieldPolicy::kReject) {
        throw std::runtime_error(
            field_path +
            " is missing and its normalized missing-field policy is reject");
      }
      if (policy == IcMissingFieldPolicy::kReconstruct) {
        throw std::runtime_error(
            field_path +
            " requests reconstruction, but gadget_arepo_bridge_v1 has no "
            "validated reconstruction contract for this field");
      }
      if (policy == IcMissingFieldPolicy::kPreserveUnavailable) {
        throw std::runtime_error(
            field_path +
            " requests preserve_unavailable, but runtime sidecars do not "
            "carry availability masks");
      }
      if (policy == IcMissingFieldPolicy::kDialectDefinedDefault &&
          dialect_resolution.starts_with("no dialect-defined default")) {
        throw std::runtime_error(
            field_path +
            " requests dialect_defined_default, but this dialect defines no "
            "scientifically valid default for the field");
      }
      const std::string resolution =
          policy == IcMissingFieldPolicy::kUseConfigValue
          ? "use normalized configuration value in runtime code units"
          : std::move(dialect_resolution);
      result.missing_field_contracts.push_back(IcMissingFieldContract{
          .source_file_index = file_index,
          .field_path = field_path,
          .policy = policy,
          .configured_value_code = configured_value_code,
          .resolution = resolution});
      result.defaulted_fields.push_back(
          field_path + "=" + resolution);
    };

    add("Coordinates", {"Coordinates", "Position", "POS"}, true,
        IcFieldSemantics::kCoordinate);
    add("Velocities", {"Velocities", "Velocity", "VEL"},
        options.require_velocities, IcFieldSemantics::kVelocity,
        convention.velocity);
    add("ParticleIDs", {"ParticleIDs", "ParticleID", "ID"},
        options.require_particle_ids, IcFieldSemantics::kIdentifier);
    const bool has_masses = add(
        "Masses", {"Masses", "Mass"}, false,
        IcFieldSemantics::kExtensive);
    if (!has_masses && result.schema.mass_table[type] <= 0.0) {
      throw std::runtime_error(
          group_path + " requires Masses because MassTable is zero");
    }

    if (type == 0U) {
      if (!add("InternalEnergy", {"InternalEnergy", "U", "Internal_Energy"},
               false, IcFieldSemantics::kSpecific)) {
        resolve_missing(
            group_path + "/InternalEnergy",
            config.mode.ic_gas_internal_energy_policy,
            config.mode.ic_gas_internal_energy_value_code,
            "no dialect-defined default exists");
      }
      if (!add("Density", {"Density", "Rho"}, false,
               IcFieldSemantics::kIntensive)) {
        resolve_missing(
            group_path + "/Density", config.mode.ic_gas_density_policy,
            config.mode.ic_gas_density_value_code,
            "no dialect-defined default exists");
      }
      if (!add("Metallicity", {"Metallicity", "GFM_Metallicity"}, false,
               IcFieldSemantics::kIntensive)) {
        resolve_missing(
            group_path + "/Metallicity",
            core::InitialConditionMissingFieldPolicy::kDialectDefinedDefault,
            0.0,
            "use zero gas metallicity as the explicit dialect default");
      }
      if (add("SmoothingLength",
              {"SmoothingLength", "Hsml", "Smoothing_Length"}, false,
              IcFieldSemantics::kCoordinate,
              IcVelocityConvention::kNotVelocity,
              IcFieldDisposition::kDropped)) {
        result.dropped_fields.push_back(
            group_path + "/SmoothingLength: no gas smoothing-length lane");
      }
    }
    if (species_policy[type] == IcSpeciesPolicy::kStar) {
      if (!add("StellarFormationTime",
               {"GFM_StellarFormationTime", "StellarFormationTime",
                "BirthTime"},
               false, IcFieldSemantics::kIntensive)) {
        resolve_missing(
            group_path + "/StellarFormationTime",
            config.mode.ic_star_formation_time_policy,
            config.mode.ic_star_formation_time_value,
            "use Header/Time as the explicitly selected dialect default");
      }
      if (!add("InitialMass",
               {"GFM_InitialMass", "InitialMass", "BirthMass"}, false,
               IcFieldSemantics::kExtensive)) {
        resolve_missing(
            group_path + "/InitialMass",
            config.mode.ic_star_initial_mass_policy,
            config.mode.ic_star_initial_mass_value_code,
            "use current particle mass as the explicitly selected dialect default");
      }
      if (!add("Metallicity", {"GFM_Metallicity", "Metallicity"}, false,
               IcFieldSemantics::kIntensive)) {
        resolve_missing(
            group_path + "/Metallicity",
            config.mode.ic_star_metallicity_policy,
            config.mode.ic_star_metallicity_value,
            "use zero metallicity as the explicitly selected dialect default");
      }
    }
    if (species_policy[type] == IcSpeciesPolicy::kBlackHole) {
      add("BH_Mass", {"BH_Mass", "BlackHoleMass"}, true,
          IcFieldSemantics::kExtensive);
      if (!add("BH_Mdot", {"BH_Mdot", "BlackHoleAccretionRate"}, false,
               IcFieldSemantics::kIntensive)) {
        resolve_missing(
            group_path + "/BH_Mdot", config.mode.ic_bh_mdot_policy,
            config.mode.ic_bh_mdot_value_code,
            "use zero accretion rate as the explicitly selected dialect default");
      }
    }
    if (species_policy[type] == IcSpeciesPolicy::kTracer) {
      add("ParentParticleIDs", {"ParentParticleIDs", "TracerParentIDs"},
          true, IcFieldSemantics::kIdentifier);
      add("InjectionStep", {"InjectionStep"}, true,
          IcFieldSemantics::kIdentifier);
      if (add(
              "HostCellIndex", {"HostCellIndex"}, false,
              IcFieldSemantics::kIdentifier,
              IcVelocityConvention::kNotVelocity,
              IcFieldDisposition::kDropped)) {
        result.dropped_fields.push_back(
            group_path +
            "/HostCellIndex: source-local row is remapped from ParentParticleIDs");
      }
      add("MassFractionOfHost", {"MassFractionOfHost"}, true,
          IcFieldSemantics::kIntensive);
      add("LastHostMass", {"LastHostMass"}, true,
          IcFieldSemantics::kExtensive);
      add("CumulativeExchangedMass", {"CumulativeExchangedMass"}, true,
          IcFieldSemantics::kExtensive);
    }

    for (const std::string& dataset_name : listLinkNames(group.get())) {
      if (handled_aliases.contains(dataset_name)) {
        continue;
      }
      if (objectType(group.get(), dataset_name) != H5O_TYPE_DATASET) {
        throw std::runtime_error(
            group_path + "/" + dataset_name +
            " is not a supported scalar/vector dataset");
      }
      result.fields.push_back(inspectDroppedDataset(
          group.get(), file_index, group_path + "/" + dataset_name,
          dataset_name));
      const std::string warning =
          group_path + "/" + dataset_name +
          ": unsupported dataset explicitly dropped";
      result.dropped_fields.push_back(warning);
      result.warnings.push_back(warning);
    }
  }

  result.counters.bytes_read = result.counters.logical_metadata_bytes_read +
      result.counters.hash_bytes_read;
  return result;
}

[[nodiscard]] Inspection inspectFileSet(
    const std::filesystem::path& requested,
    const core::SimulationConfig& config,
    const IcImportOptions& options) {
  const bool has_authoritative_manifest = options.manifest != nullptr;
  const std::filesystem::path first_source =
      has_authoritative_manifest && !options.manifest->source_files.empty()
      ? options.manifest->source_files.front()
      : requested;
  Hdf5Handle first_file(
      H5Fopen(first_source.string().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
  if (!first_file.valid()) {
    throw std::runtime_error(
        "failed to open IC file: " + first_source.string());
  }
  Hdf5Handle first_header(
      H5Gopen2(first_file.get(), "/Header", H5P_DEFAULT));
  if (!first_header.valid()) {
    throw std::runtime_error("IC file missing /Header group");
  }
  const IcSchemaSummary first_schema = readHeader(first_header.get());

  std::vector<std::filesystem::path> files;
  if (has_authoritative_manifest && !options.manifest->source_files.empty()) {
    files = options.manifest->source_files;
  } else {
    files = discoverFiles(requested, first_schema.num_files_per_snapshot);
  }
  if (files.size() != first_schema.num_files_per_snapshot) {
    throw std::runtime_error(
        "manifest/discovery file count disagrees with NumFilesPerSnapshot");
  }

  Inspection inspection;
  IcManifest& manifest = inspection.manifest;
  manifest.dialect =
      config.mode.ic_convention ==
              core::InitialConditionConvention::kChuiCanonicalV1
          ? IcDialect::kChuiCanonicalV1
          : IcDialect::kGadgetArepoBridgeV1;
  manifest.dialect_version = "1";
  manifest.converter_version = "chui_runtime_inspector_v4";
  manifest.source_files = files;
  manifest.num_files_per_snapshot =
      static_cast<std::uint32_t>(files.size());
  manifest.species_policy[2] =
      mapConfiguredPolicy(config.mode.ic_part_type2_policy, 2U);
  manifest.species_policy[3] =
      mapConfiguredPolicy(config.mode.ic_part_type3_policy, 3U);
  if (options.manifest != nullptr) {
    manifest.dialect = options.manifest->dialect;
    manifest.dialect_version = options.manifest->dialect_version;
    manifest.species_policy = options.manifest->species_policy;
  }
  if (manifest.dialect_version != "1") {
    throw std::runtime_error(
        "unsupported IC dialect version: " + manifest.dialect_version);
  }
  if (config.mode.ic_convention ==
          core::InitialConditionConvention::kChuiCanonicalV1 &&
      manifest.dialect != IcDialect::kChuiCanonicalV1) {
    throw std::runtime_error(
        "chui_canonical_v1 configuration requires a canonical CHUÍ "
        "manifest/input");
  }
  if (config.mode.ic_convention ==
          core::InitialConditionConvention::kGadgetArepoBridgeV1 &&
      manifest.dialect != IcDialect::kGadgetArepoBridgeV1) {
    throw std::runtime_error(
        "gadget_arepo_bridge_v1 configuration requires a bridge "
        "manifest/input");
  }

  std::array<std::uint64_t, 6> summed{};
  for (std::size_t file_index = 0; file_index < files.size(); ++file_index) {
    SourceFileInspection source = inspectOneSourceFile(
        files[file_index], static_cast<std::uint32_t>(file_index),
        manifest.dialect, manifest.species_policy, config, options,
        has_authoritative_manifest);
    const IcSchemaSummary& schema = source.schema;
    if (schema.num_files_per_snapshot != files.size() ||
        schema.total_count_by_type != first_schema.total_count_by_type ||
        schema.total_count_high_word != first_schema.total_count_high_word ||
        schema.mass_table != first_schema.mass_table ||
        !nearlyEqual(schema.box_size, first_schema.box_size) ||
        !nearlyEqual(schema.scale_factor, first_schema.scale_factor) ||
        !nearlyEqual(schema.redshift, first_schema.redshift) ||
        !nearlyEqual(schema.omega_matter, first_schema.omega_matter) ||
        !nearlyEqual(schema.omega_lambda, first_schema.omega_lambda) ||
        !nearlyEqual(schema.hubble_param, first_schema.hubble_param)) {
      throw std::runtime_error(
          "inconsistent cosmology, box, epoch, mass table, totals, or "
          "NumFilesPerSnapshot across IC files");
    }
    for (std::size_t type = 0; type < kParticleTypeCount; ++type) {
      if (summed[type] > std::numeric_limits<std::uint64_t>::max() -
              schema.count_by_type[type]) {
        throw std::overflow_error("IC file-set particle count overflow");
      }
      summed[type] += schema.count_by_type[type];
    }
    inspection.schemas.push_back(schema);
    manifest.num_part_this_file.push_back(schema.count_by_type);
    manifest.source_file_sizes_bytes.push_back(source.source_size_bytes);
    manifest.source_sha256.push_back(source.source_sha256);
    manifest.source_provenance_ids.push_back(
        "sha256:" + source.source_sha256);
    manifest.original_header_attributes.push_back(
        std::move(source.original_header_attributes));
    manifest.fields.insert(
        manifest.fields.end(),
        std::make_move_iterator(source.fields.begin()),
        std::make_move_iterator(source.fields.end()));
    manifest.missing_field_contracts.insert(
        manifest.missing_field_contracts.end(),
        std::make_move_iterator(source.missing_field_contracts.begin()),
        std::make_move_iterator(source.missing_field_contracts.end()));
    manifest.defaulted_fields.insert(
        manifest.defaulted_fields.end(),
        std::make_move_iterator(source.defaulted_fields.begin()),
        std::make_move_iterator(source.defaulted_fields.end()));
    manifest.dropped_fields.insert(
        manifest.dropped_fields.end(),
        std::make_move_iterator(source.dropped_fields.begin()),
        std::make_move_iterator(source.dropped_fields.end()));
    manifest.preserved_auxiliary_fields.insert(
        manifest.preserved_auxiliary_fields.end(),
        std::make_move_iterator(source.preserved_auxiliary_fields.begin()),
        std::make_move_iterator(source.preserved_auxiliary_fields.end()));
    manifest.warnings.insert(
        manifest.warnings.end(),
        std::make_move_iterator(source.warnings.begin()),
        std::make_move_iterator(source.warnings.end()));
    if (source.canonical_manifest_verified) {
      if (manifest.canonical_source_manifest_verified &&
          manifest.canonical_source_manifest_sha256 !=
              source.canonical_manifest_sha256) {
        throw std::runtime_error(
            "canonical CHUÍ IC file set has inconsistent manifest bindings");
      }
      manifest.canonical_source_manifest_verified = true;
      manifest.canonical_source_manifest_sha256 =
          source.canonical_manifest_sha256;
    }
    checkedCounterAdd(
        inspection.counters.logical_metadata_bytes_read,
        source.counters.logical_metadata_bytes_read, "logical_metadata_bytes_read");
    checkedCounterAdd(
        inspection.counters.hash_bytes_read,
        source.counters.hash_bytes_read, "hash_bytes_read");
    checkedCounterAdd(
        inspection.counters.full_file_hash_pass_count,
        source.counters.full_file_hash_pass_count,
        "full_file_hash_pass_count");
  }

  validateCrossFileSchema(manifest);
  if (summed != first_schema.total_count_by_type) {
    throw std::runtime_error(
        "summed NumPart_ThisFile does not equal reconstructed 64-bit "
        "NumPart_Total");
  }
  manifest.num_part_total = first_schema.total_count_by_type;
  manifest.num_part_total_high_word = first_schema.total_count_high_word;
  manifest.mass_table = first_schema.mass_table;
  manifest.box_size = first_schema.box_size;
  manifest.scale_factor = first_schema.scale_factor;
  manifest.redshift = first_schema.redshift;
  manifest.omega_matter = first_schema.omega_matter;
  manifest.omega_lambda = first_schema.omega_lambda;
  manifest.hubble_param = first_schema.hubble_param;
  manifest.converted_fields.reserve(manifest.fields.size());
  for (const auto& field : manifest.fields) {
    if (field.disposition == IcFieldDisposition::kConverted) {
      manifest.converted_fields.push_back(field.dataset_path);
      if (std::find(
              manifest.conversion_equations.begin(),
              manifest.conversion_equations.end(), field.conversion_equation) ==
          manifest.conversion_equations.end()) {
        manifest.conversion_equations.push_back(field.conversion_equation);
      }
    }
  }

  if (has_authoritative_manifest) {
    const IcManifest& supplied = *options.manifest;
    if (supplied.num_files_per_snapshot !=
            manifest.num_files_per_snapshot ||
        supplied.source_files != manifest.source_files ||
        supplied.source_sha256 != manifest.source_sha256 ||
        supplied.source_file_sizes_bytes !=
            manifest.source_file_sizes_bytes ||
        supplied.num_part_this_file != manifest.num_part_this_file ||
        supplied.num_part_total != manifest.num_part_total ||
        supplied.num_part_total_high_word !=
            manifest.num_part_total_high_word ||
        supplied.mass_table != manifest.mass_table ||
        supplied.species_policy != manifest.species_policy ||
        !missingFieldContractsEqual(
            supplied.missing_field_contracts,
            manifest.missing_field_contracts) ||
        !nearlyEqual(supplied.box_size, manifest.box_size) ||
        !nearlyEqual(supplied.scale_factor, manifest.scale_factor) ||
        !nearlyEqual(supplied.redshift, manifest.redshift) ||
        !nearlyEqual(supplied.omega_matter, manifest.omega_matter) ||
        !nearlyEqual(supplied.omega_lambda, manifest.omega_lambda) ||
        !nearlyEqual(supplied.hubble_param, manifest.hubble_param)) {
      throw std::runtime_error(
          "supplied IC manifest provenance, scientific header, policies, or "
          "counts do not match inspected source files");
    }
    if (manifest.fields.size() != supplied.fields.size()) {
      throw std::runtime_error(
          "supplied IC manifest field count does not match actual HDF5 "
          "schema");
    }
    for (const auto& expected : supplied.fields) {
      const auto observed = std::find_if(
          manifest.fields.begin(), manifest.fields.end(),
          [&](const IcFieldManifest& actual) {
            return actual.source_file_index == expected.source_file_index &&
                actual.dataset_path == expected.dataset_path &&
                actual.selected_alias == expected.selected_alias;
          });
      if (observed == manifest.fields.end() ||
          observed->scalar_type != expected.scalar_type ||
          observed->scalar_class != expected.scalar_class ||
          observed->byte_width != expected.byte_width ||
          observed->is_signed != expected.is_signed ||
          observed->byte_order != expected.byte_order ||
          observed->rank != expected.rank ||
          observed->dimensions != expected.dimensions ||
          observed->record_count != expected.record_count ||
          observed->disposition != expected.disposition) {
        throw std::runtime_error(
            "supplied IC manifest schema does not match actual HDF5 field " +
            expected.dataset_path);
      }
    }
    manifest = supplied;
  }
  validateIcManifest(manifest);
  inspection.counters.bytes_read =
      inspection.counters.logical_metadata_bytes_read +
      inspection.counters.hash_bytes_read;
  return inspection;
}

#endif  // COSMOSIM_ENABLE_HDF5

}  // namespace cosmosim::io::file_set_internal
