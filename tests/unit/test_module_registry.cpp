#include "cosmosim/core/module_registry.hpp"

#include <cassert>
#include <stdexcept>

int main() {
  const auto names = cosmosim::core::moduleNames();
  const auto descriptors = cosmosim::core::moduleDescriptors();
  assert(descriptors.size() == names.size());
  for (std::size_t index = 0; index < descriptors.size(); ++index) {
    const auto& descriptor = descriptors[index];
    assert(descriptor.name == names[index]);
    assert(!descriptor.config_fragment.empty());
    assert(!descriptor.state_and_sidecar_requirements.empty());
    assert(!descriptor.stage_and_task_contributions.empty());
    assert(!descriptor.restart_payloads.empty());
    assert(!descriptor.migration_fields.empty());
    assert(!descriptor.diagnostics.empty());
    assert(!descriptor.capability_prerequisites.empty());
    assert(!descriptor.incompatibilities.empty());
    assert(&cosmosim::core::requireModuleDescriptor(descriptor.name) ==
           &descriptor);
  }

  bool rejected_unknown = false;
  try {
    (void)cosmosim::core::requireModuleDescriptor("not_a_module");
  } catch (const std::out_of_range&) {
    rejected_unknown = true;
  }
  assert(rejected_unknown);
  return 0;
}
