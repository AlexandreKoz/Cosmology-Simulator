#include "cosmosim/cosmosim.hpp"

#include <sstream>
#include <stdexcept>
#include <string>

namespace cosmosim::core {

std::array<std::string_view, k_module_names.size()> moduleNames() {
  return k_module_names;
}

std::span<const ModuleDescriptor> moduleDescriptors() noexcept {
  return k_module_descriptors;
}

const ModuleDescriptor& requireModuleDescriptor(std::string_view name) {
  for (const ModuleDescriptor& descriptor : k_module_descriptors) {
    if (descriptor.name == name) {
      return descriptor;
    }
  }
  throw std::out_of_range("unknown runtime module descriptor: " +
                          std::string(name));
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
