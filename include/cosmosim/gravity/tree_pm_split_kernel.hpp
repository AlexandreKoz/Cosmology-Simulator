#pragma once

#include <cmath>
#include <stdexcept>

namespace cosmosim::gravity {

// Explicitly documented TreePM Gaussian split metadata and utility API.
enum class TreePmSplitKernel {
  kGaussianErfc,
};

struct TreePmSplitPolicy {
  TreePmSplitKernel kernel = TreePmSplitKernel::kGaussianErfc;
  double mesh_spacing_comoving = 0.0;
  double asmth_cells = 0.0;
  double rcut_cells = 0.0;
  double split_scale_comoving = 0.0;
  double cutoff_radius_comoving = 0.0;
};

inline void validateTreePmSplitPolicy(const TreePmSplitPolicy& policy) {
  if (policy.mesh_spacing_comoving <= 0.0) {
    throw std::invalid_argument("TreePM mesh_spacing_comoving must be > 0");
  }
  if (policy.asmth_cells <= 0.0) {
    throw std::invalid_argument("TreePM asmth_cells must be > 0");
  }
  if (policy.rcut_cells <= 0.0) {
    throw std::invalid_argument("TreePM rcut_cells must be > 0");
  }
  if (policy.split_scale_comoving <= 0.0) {
    throw std::invalid_argument("TreePM split_scale_comoving must be > 0");
  }
  if (policy.cutoff_radius_comoving <= 0.0) {
    throw std::invalid_argument("TreePM cutoff_radius_comoving must be > 0");
  }
  if (policy.kernel != TreePmSplitKernel::kGaussianErfc) {
    throw std::invalid_argument("Unsupported TreePM split kernel");
  }
}

[[nodiscard]] inline TreePmSplitPolicy makeTreePmSplitPolicyFromMeshSpacing(
    double asmth_cells,
    double rcut_cells,
    double mesh_spacing_comoving,
    TreePmSplitKernel kernel = TreePmSplitKernel::kGaussianErfc) {
  TreePmSplitPolicy policy;
  policy.kernel = kernel;
  policy.mesh_spacing_comoving = mesh_spacing_comoving;
  policy.asmth_cells = asmth_cells;
  policy.rcut_cells = rcut_cells;
  policy.split_scale_comoving = asmth_cells * mesh_spacing_comoving;
  policy.cutoff_radius_comoving = rcut_cells * mesh_spacing_comoving;
  validateTreePmSplitPolicy(policy);
  return policy;
}

[[nodiscard]] inline double treePmGaussianShortRangeForceFactor(double distance_comoving, double split_scale_comoving) {
  if (distance_comoving <= 0.0) {
    return 0.0;
  }
  const double q = distance_comoving / (2.0 * split_scale_comoving);
  constexpr double inv_sqrt_pi = 0.564189583547756286948079451560772586;
  return std::erfc(q) + 2.0 * inv_sqrt_pi * q * std::exp(-q * q);
}

[[nodiscard]] inline double treePmGaussianLongRangeForceFactor(double distance_comoving, double split_scale_comoving) {
  return 1.0 - treePmGaussianShortRangeForceFactor(distance_comoving, split_scale_comoving);
}

[[nodiscard]] inline double treePmGaussianFourierLongRangeFilter(double wave_number_comoving, double split_scale_comoving) {
  const double product = wave_number_comoving * split_scale_comoving;
  return std::exp(-(product * product));
}

}  // namespace cosmosim::gravity
