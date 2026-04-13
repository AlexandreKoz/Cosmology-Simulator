#include "cosmosim/cosmosim.hpp"

#include <sstream>

namespace cosmosim::core {

std::array<std::string_view, k_module_names.size()> moduleNames() {
  return k_module_names;
}

}  // namespace cosmosim::core

namespace cosmosim {

std::string architectureSummary() {
  std::ostringstream stream;
  stream << "CosmoSim layers:";
  for (const std::string_view layer : k_layer_order) {
    stream << ' ' << layer;
  }
  return stream.str();
}

}  // namespace cosmosim
