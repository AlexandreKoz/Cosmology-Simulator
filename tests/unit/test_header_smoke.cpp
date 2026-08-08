#include <cassert>
#include <string>

#include "cosmosim/cosmosim.hpp"

int main() {
  const auto modules = cosmosim::core::moduleNames();
  assert(!modules.empty());
  assert(cosmosim::core::projectName() == "cosmosim");
  assert(!cosmosim::architectureSummary().empty());
  assert(!cosmosim::core::buildProvenance().empty());
  assert(cosmosim::core::buildProvenance().find("preset=") != std::string::npos);
  return 0;
}
