#pragma once

#include <string>
#include <string_view>

namespace cosmosim::core {

struct Version {
  int major;
  int minor;
  int patch;
};

Version version();
std::string versionString();
std::string buildProvenance();
std::string_view projectName();

}  // namespace cosmosim::core
