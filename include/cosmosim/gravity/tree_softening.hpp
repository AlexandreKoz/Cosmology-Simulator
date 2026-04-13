#pragma once

#include <cmath>

namespace cosmosim::gravity {

// Conservative softening modes for short-range gravity force evaluation.
enum class TreeSofteningKernel {
  kPlummer,
};

struct TreeSofteningPolicy {
  TreeSofteningKernel kernel = TreeSofteningKernel::kPlummer;
  double epsilon_comoving = 0.0;
};

[[nodiscard]] inline double softenedInvR3(double squared_distance, const TreeSofteningPolicy& policy) {
  const double epsilon2 = policy.epsilon_comoving * policy.epsilon_comoving;
  const double denominator = std::pow(squared_distance + epsilon2, 1.5);
  if (denominator <= 0.0) {
    return 0.0;
  }
  return 1.0 / denominator;
}

}  // namespace cosmosim::gravity
