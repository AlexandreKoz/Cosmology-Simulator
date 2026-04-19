#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <span>
#include <string>
#include <vector>

#include "cosmosim/gravity/pm_solver.hpp"
#include "cosmosim/gravity/tree_pm_coupling.hpp"

namespace {

struct ForceField {
  std::vector<double> ax;
  std::vector<double> ay;
  std::vector<double> az;
};

[[nodiscard]] double minimumImageDelta(double delta, double box_size_comoving) {
  return delta - box_size_comoving * std::nearbyint(delta / box_size_comoving);
}

ForceField solveTreePm(
    std::size_t pm_grid,
    std::span<const double> pos_x,
    std::span<const double> pos_y,
    std::span<const double> pos_z,
    std::span<const double> mass,
    double asmth_cells,
    double rcut_cells,
    double theta,
    std::size_t max_leaf) {
  cosmosim::gravity::TreePmOptions options;
  options.pm_options.box_size_mpc_comoving = 1.0;
  options.pm_options.scale_factor = 1.0;
  options.pm_options.gravitational_constant_code = 1.0;
  options.pm_options.enable_window_deconvolution = true;
  options.pm_options.assignment_scheme = cosmosim::gravity::PmAssignmentScheme::kTsc;
  options.tree_options.opening_theta = theta;
  options.tree_options.max_leaf_size = max_leaf;
  options.tree_options.gravitational_constant_code = 1.0;
  options.tree_options.softening.epsilon_comoving = 0.01;
  const double mesh_spacing = 1.0 / static_cast<double>(pm_grid);
  options.split_policy = cosmosim::gravity::makeTreePmSplitPolicyFromMeshSpacing(asmth_cells, rcut_cells, mesh_spacing);

  std::vector<std::uint32_t> active(pos_x.size(), 0U);
  for (std::size_t i = 0; i < active.size(); ++i) {
    active[i] = static_cast<std::uint32_t>(i);
  }

  ForceField f{
      std::vector<double>(pos_x.size(), 0.0),
      std::vector<double>(pos_x.size(), 0.0),
      std::vector<double>(pos_x.size(), 0.0)};
  cosmosim::gravity::TreePmForceAccumulatorView acc{active, f.ax, f.ay, f.az};
  cosmosim::gravity::TreePmCoordinator coordinator({pm_grid, pm_grid, pm_grid});
  coordinator.solveActiveSet(pos_x, pos_y, pos_z, mass, acc, options, nullptr, nullptr);
  return f;
}

ForceField solveFinePmLongRangeReference(
    std::size_t pm_grid,
    std::span<const double> pos_x,
    std::span<const double> pos_y,
    std::span<const double> pos_z,
    std::span<const double> mass,
    double split_scale_comoving) {
  cosmosim::gravity::PmGridStorage grid({pm_grid, pm_grid, pm_grid});
  cosmosim::gravity::PmSolver solver({pm_grid, pm_grid, pm_grid});
  cosmosim::gravity::PmSolveOptions options;
  options.box_size_mpc_comoving = 1.0;
  options.scale_factor = 1.0;
  options.gravitational_constant_code = 1.0;
  options.enable_window_deconvolution = true;
  options.assignment_scheme = cosmosim::gravity::PmAssignmentScheme::kTsc;
  options.tree_pm_split_scale_comoving = split_scale_comoving;

  ForceField field{
      std::vector<double>(pos_x.size(), 0.0),
      std::vector<double>(pos_x.size(), 0.0),
      std::vector<double>(pos_x.size(), 0.0)};
  solver.assignDensity(grid, pos_x, pos_y, pos_z, mass, options, nullptr);
  solver.solvePoissonPeriodic(grid, options, nullptr);
  solver.interpolateForces(grid, pos_x, pos_y, pos_z, field.ax, field.ay, field.az, options, nullptr);
  return field;
}

ForceField computeDirectShortRangeResidual(
    std::span<const double> pos_x,
    std::span<const double> pos_y,
    std::span<const double> pos_z,
    std::span<const double> mass,
    const cosmosim::gravity::TreePmOptions& options) {
  ForceField field{
      std::vector<double>(pos_x.size(), 0.0),
      std::vector<double>(pos_x.size(), 0.0),
      std::vector<double>(pos_x.size(), 0.0)};
  const double cutoff2 = options.split_policy.cutoff_radius_comoving * options.split_policy.cutoff_radius_comoving;
  for (std::size_t target = 0; target < pos_x.size(); ++target) {
    for (std::size_t source = 0; source < pos_x.size(); ++source) {
      if (target == source) {
        continue;
      }
      const double dx = minimumImageDelta(pos_x[source] - pos_x[target], options.pm_options.box_size_mpc_comoving);
      const double dy = minimumImageDelta(pos_y[source] - pos_y[target], options.pm_options.box_size_mpc_comoving);
      const double dz = minimumImageDelta(pos_z[source] - pos_z[target], options.pm_options.box_size_mpc_comoving);
      const double r2 = dx * dx + dy * dy + dz * dz;
      if (r2 > cutoff2) {
        continue;
      }
      const double r = std::sqrt(std::max(r2, 1.0e-30));
      const double split_factor = cosmosim::gravity::treePmGaussianShortRangeForceFactor(
          r, options.split_policy.split_scale_comoving);
      const double softened_factor = cosmosim::gravity::softenedInvR3(r2, options.tree_options.softening) *
          split_factor * options.tree_options.gravitational_constant_code;
      field.ax[target] += softened_factor * mass[source] * dx;
      field.ay[target] += softened_factor * mass[source] * dy;
      field.az[target] += softened_factor * mass[source] * dz;
    }
  }
  return field;
}

ForceField addFields(const ForceField& lhs, const ForceField& rhs) {
  ForceField out{
      std::vector<double>(lhs.ax.size(), 0.0),
      std::vector<double>(lhs.ay.size(), 0.0),
      std::vector<double>(lhs.az.size(), 0.0)};
  for (std::size_t i = 0; i < lhs.ax.size(); ++i) {
    out.ax[i] = lhs.ax[i] + rhs.ax[i];
    out.ay[i] = lhs.ay[i] + rhs.ay[i];
    out.az[i] = lhs.az[i] + rhs.az[i];
  }
  return out;
}

double relativeL2(const ForceField& lhs, const ForceField& rhs) {
  double ref2 = 0.0;
  double err2 = 0.0;
  for (std::size_t i = 0; i < lhs.ax.size(); ++i) {
    const double dx = lhs.ax[i] - rhs.ax[i];
    const double dy = lhs.ay[i] - rhs.ay[i];
    const double dz = lhs.az[i] - rhs.az[i];
    err2 += dx * dx + dy * dy + dz * dz;
    ref2 += rhs.ax[i] * rhs.ax[i] + rhs.ay[i] * rhs.ay[i] + rhs.az[i] * rhs.az[i];
  }
  return std::sqrt(err2 / std::max(ref2, 1.0e-24));
}

}  // namespace

int main() {
  constexpr std::size_t particle_count = 24;
  std::vector<double> pos_x(particle_count, 0.0);
  std::vector<double> pos_y(particle_count, 0.0);
  std::vector<double> pos_z(particle_count, 0.0);
  std::vector<double> mass(particle_count, 1.0);

  for (std::size_t i = 0; i < particle_count; ++i) {
    pos_x[i] = std::fmod((13.0 * static_cast<double>(i) + 1.0) * 0.047, 1.0);
    pos_y[i] = std::fmod((19.0 * static_cast<double>(i) + 2.0) * 0.031, 1.0);
    pos_z[i] = std::fmod((23.0 * static_cast<double>(i) + 3.0) * 0.029, 1.0);
    mass[i] = 0.9 + 0.05 * static_cast<double>(i % 5U);
  }

  const std::filesystem::path out_dir = std::filesystem::path(COSMOSIM_SOURCE_DIR) / "validation" / "artifacts";
  std::filesystem::create_directories(out_dir);
  const std::filesystem::path csv_path = out_dir / "tree_pm_force_error_map.csv";

  std::ofstream out(csv_path);
  if (!out) {
    std::cerr << "failed to open output CSV: " << csv_path << '\n';
    return 2;
  }

  out << "reference_method,pm_grid,asmth_cells,rcut_cells,relative_l2_error\n";
  out << std::setprecision(12);

  for (const std::size_t pm_grid : {16U, 24U, 32U}) {
    for (const double asmth_cells : {0.8, 1.25, 2.0}) {
      for (const double rcut_cells : {3.0, 4.5, 6.0}) {
        const ForceField candidate = solveTreePm(pm_grid, pos_x, pos_y, pos_z, mass, asmth_cells, rcut_cells, 0.55, 8);

        cosmosim::gravity::TreePmOptions ref_options;
        ref_options.pm_options.box_size_mpc_comoving = 1.0;
        ref_options.pm_options.scale_factor = 1.0;
        ref_options.pm_options.gravitational_constant_code = 1.0;
        ref_options.tree_options.gravitational_constant_code = 1.0;
        ref_options.tree_options.softening.epsilon_comoving = 0.01;
        ref_options.split_policy = cosmosim::gravity::makeTreePmSplitPolicyFromMeshSpacing(
            asmth_cells, rcut_cells, 1.0 / static_cast<double>(pm_grid));
        const ForceField long_range = solveFinePmLongRangeReference(
            96,
            pos_x,
            pos_y,
            pos_z,
            mass,
            ref_options.split_policy.split_scale_comoving);
        const ForceField short_range = computeDirectShortRangeResidual(pos_x, pos_y, pos_z, mass, ref_options);
        const ForceField reference = addFields(long_range, short_range);

        const double rel = relativeL2(candidate, reference);
        out << "periodic_spectral_direct_proxy"
            << ',' << pm_grid
            << ',' << asmth_cells
            << ',' << rcut_cells
            << ',' << rel
            << '\n';
      }
    }
  }

  out.flush();
  std::cout << "bench_tree_pm_force_error_map csv=" << csv_path.string() << "\n";
  return 0;
}
