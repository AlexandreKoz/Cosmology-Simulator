#include "cosmosim/core/version.hpp"

#include <sstream>

#include "cosmosim/core/build_config.hpp"

namespace cosmosim::core {

Version version() {
  return Version{COSMOSIM_VERSION_MAJOR, COSMOSIM_VERSION_MINOR, COSMOSIM_VERSION_PATCH};
}

std::string versionString() {
  const Version v = version();
  std::ostringstream stream;
  stream << v.major << "." << v.minor << "." << v.patch;
  return stream.str();
}

std::string buildProvenance() {
  std::ostringstream stream;
  stream << "project=" << projectName() << ";version=" << versionString() << ";preset="
         << COSMOSIM_BUILD_PRESET << ";build_type=" << COSMOSIM_BUILD_TYPE << ";mpi="
         << COSMOSIM_ENABLE_MPI << ";hdf5=" << COSMOSIM_ENABLE_HDF5 << ";fftw="
         << COSMOSIM_ENABLE_FFTW << ";cuda=" << COSMOSIM_ENABLE_CUDA;
  return stream.str();
}

std::string_view projectName() {
  return "cosmosim";
}

}  // namespace cosmosim::core
