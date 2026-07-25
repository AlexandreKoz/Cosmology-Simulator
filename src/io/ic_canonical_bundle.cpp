#include "io/internal/ic_canonical_bundle.hpp"

#if COSMOSIM_ENABLE_HDF5
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "cosmosim/io/ic_reader.hpp"
#include "io/internal/ic_hdf5_handle.hpp"

namespace cosmosim::io::internal {
namespace {

[[nodiscard]] std::vector<std::uint64_t> attributeDimensions(
    hid_t attribute) {
  Hdf5Handle space(H5Aget_space(attribute));
  if (!space.valid()) {
    throw std::runtime_error("failed to inspect canonical string attribute");
  }
  const int rank = H5Sget_simple_extent_ndims(space.get());
  if (rank < 0) {
    throw std::runtime_error("failed to read canonical attribute rank");
  }
  std::vector<hsize_t> raw(static_cast<std::size_t>(rank));
  if (rank > 0 && H5Sget_simple_extent_dims(space.get(), raw.data(), nullptr) < 0) {
    throw std::runtime_error("failed to read canonical attribute dimensions");
  }
  return {raw.begin(), raw.end()};
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
    throw std::runtime_error(
        std::string("Header/") + name + " must be a string");
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

[[nodiscard]] std::vector<std::uint64_t> dataspaceDimensions(hid_t space) {
  const int rank = H5Sget_simple_extent_ndims(space);
  if (rank < 0) {
    throw std::runtime_error("failed to inspect canonical dataset rank");
  }
  std::vector<hsize_t> raw(static_cast<std::size_t>(rank));
  if (rank > 0 && H5Sget_simple_extent_dims(space, raw.data(), nullptr) < 0) {
    throw std::runtime_error("failed to inspect canonical dataset dimensions");
  }
  return {raw.begin(), raw.end()};
}

[[nodiscard]] std::string readTextFileExact(
    const std::filesystem::path& path) {
  std::ifstream input(path, std::ios::binary);
  if (!input) {
    throw std::runtime_error(
        "canonical CHUÍ IC bundle is missing required file: " +
        path.string());
  }
  std::ostringstream contents;
  contents << input.rdbuf();
  if (!input.eof() && input.fail()) {
    throw std::runtime_error(
        "failed to read canonical CHUÍ IC bundle file: " + path.string());
  }
  return contents.str();
}

[[nodiscard]] std::filesystem::path resolveBundleMember(
    const std::filesystem::path& canonical_path,
    std::string_view member_name,
    std::string_view field_name) {
  const std::filesystem::path relative(member_name);
  if (member_name.empty() || relative.is_absolute() ||
      relative.has_parent_path() || relative.filename() != relative) {
    throw std::runtime_error(
        "canonical CHUÍ IC " + std::string(field_name) +
        " must name a sibling file without directory traversal");
  }
  return (canonical_path.parent_path() / relative).lexically_normal();
}

[[nodiscard]] std::string readEmbeddedManifestJson(hid_t file) {
  Hdf5Handle dataset(H5Dopen2(
      file, "/Provenance/ConversionManifestJson", H5P_DEFAULT));
  if (!dataset.valid()) {
    throw std::runtime_error(
        "canonical CHUÍ IC is missing /Provenance/ConversionManifestJson");
  }
  Hdf5Handle type(H5Dget_type(dataset.get()));
  Hdf5Handle space(H5Dget_space(dataset.get()));
  if (!type.valid() || !space.valid() ||
      H5Tget_class(type.get()) != H5T_INTEGER ||
      H5Tget_size(type.get()) != 1U ||
      H5Tget_sign(type.get()) != H5T_SGN_NONE) {
    throw std::runtime_error(
        "canonical CHUÍ IC embedded manifest must be an unsigned byte dataset");
  }
  const std::vector<std::uint64_t> dimensions =
      dataspaceDimensions(space.get());
  if (dimensions.size() != 1U || dimensions.front() == 0U ||
      dimensions.front() > std::numeric_limits<std::size_t>::max()) {
    throw std::runtime_error(
        "canonical CHUÍ IC embedded manifest must be a non-empty rank-1 byte dataset");
  }
  std::string json(static_cast<std::size_t>(dimensions.front()), '\0');
  if (H5Dread(
          dataset.get(), H5T_NATIVE_UINT8, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          json.data()) < 0) {
    throw std::runtime_error(
        "failed to read canonical CHUÍ IC embedded manifest");
  }
  return json;
}

}  // namespace

CanonicalBundleVerification verifyCanonicalBundle(
    const std::filesystem::path& canonical_path,
    hid_t file,
    hid_t header) {
  const std::string expected_digest =
      readAttributeString(header, "ConversionManifestSha256");
  const std::string sidecar_name =
      readAttributeString(header, "ConversionManifestSidecar");
  const std::string marker_name =
      readAttributeString(header, "ConversionBundleMarker");
  const std::string embedded_json = readEmbeddedManifestJson(file);
  const std::string embedded_digest = icSha256Hex(embedded_json);
  if (embedded_digest != expected_digest) {
    throw std::runtime_error(
        "canonical CHUÍ IC embedded manifest SHA-256 does not match Header/ConversionManifestSha256");
  }
  static_cast<void>(deserializeIcManifestJson(embedded_json));

  const std::filesystem::path sidecar_path = resolveBundleMember(
      canonical_path, sidecar_name, "ConversionManifestSidecar");
  const std::string sidecar_json = readTextFileExact(sidecar_path);
  if (icSha256Hex(sidecar_json) != expected_digest) {
    throw std::runtime_error(
        "canonical CHUÍ IC sidecar manifest SHA-256 does not match the embedded binding");
  }
  static_cast<void>(deserializeIcManifestJson(sidecar_json));
  if (sidecar_json != embedded_json) {
    throw std::runtime_error(
        "canonical CHUÍ IC embedded and sidecar manifests differ despite digest validation");
  }

  const std::filesystem::path marker_path = resolveBundleMember(
      canonical_path, marker_name, "ConversionBundleMarker");
  const std::string expected_marker =
      "chui_ic_bundle_v1\nsha256=" + expected_digest +
      "\ncanonical=" + canonical_path.filename().string() +
      "\nmanifest=" + sidecar_name + "\n";
  if (readTextFileExact(marker_path) != expected_marker) {
    throw std::runtime_error(
        "canonical CHUÍ IC completion marker does not match the HDF5/manifest bundle");
  }
  return {.verified = true, .manifest_sha256 = expected_digest};
}

}  // namespace cosmosim::io::internal
#endif
