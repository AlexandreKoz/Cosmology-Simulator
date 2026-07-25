#pragma once

#include "cosmosim/core/build_config.hpp"

#if COSMOSIM_ENABLE_HDF5
#include <filesystem>
#include <string>

#include <hdf5.h>

namespace cosmosim::io::internal {

struct CanonicalBundleVerification {
  bool verified = false;
  std::string manifest_sha256;
};

[[nodiscard]] CanonicalBundleVerification verifyCanonicalBundle(
    const std::filesystem::path& canonical_path,
    hid_t file,
    hid_t header);

}  // namespace cosmosim::io::internal
#endif
