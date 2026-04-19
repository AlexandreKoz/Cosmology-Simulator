#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <span>
#include <string>
#include <vector>

#include "cosmosim/gravity/tree_pm_coupling.hpp"

namespace {

struct ForceField {
  std::vector<double> ax;
  std::vector<double> ay;
  std::vector<double> az;
};

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
  constexpr std::size_t particle_count = 64;
  std::vector<double> pos_x(particle_count, 0.0);
  std::vector<double> pos_y(particle_count, 0.0);
  std::vector<double> pos_z(particle_count, 0.0);
  std::vector<double> mass(particle_count, 1.0);

  for (std::size_t i = 0; i < particle_count; ++i) {
    pos_x[i] = std::fmod((13.0 * static_cast<double>(i) + 1.0) * 0.023, 1.0);
    pos_y[i] = std::fmod((19.0 * static_cast<double>(i) + 2.0) * 0.019, 1.0);
    pos_z[i] = std::fmod((23.0 * static_cast<double>(i) + 3.0) * 0.017, 1.0);
    mass[i] = 0.9 + 0.05 * static_cast<double>(i % 5U);
  }

  const ForceField periodic_proxy_reference = solveTreePm(
      64,
      pos_x,
      pos_y,
      pos_z,
      mass,
      1.25,
      6.0,
      0.35,
      2);

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
        const double rel = relativeL2(candidate, periodic_proxy_reference);
        out << "periodic_proxy_treepm_highres"
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
