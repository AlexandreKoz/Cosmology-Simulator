#include <cassert>

#include "cosmosim/cosmosim.hpp"

int main() {
  const auto version = cosmosim::core::version();
  assert(version.major >= 0);
  assert(version.minor >= 0);
  assert(version.patch >= 0);

  const auto modules = cosmosim::core::moduleNames();
  assert(modules.size() == 9);
  assert(modules.front() == "core");
  assert(modules.back() == "utils");

  return 0;
}
