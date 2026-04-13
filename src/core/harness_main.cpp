#include <iostream>

#include "cosmosim/cosmosim.hpp"

int main() {
  std::cout << cosmosim::core::projectName() << " " << cosmosim::core::versionString() << '\n';
  std::cout << cosmosim::architectureSummary() << '\n';
  return 0;
}
