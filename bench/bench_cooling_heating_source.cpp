#include <chrono>
#include <cstddef>
#include <iostream>
#include <vector>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/hydro/hydro_core_solver.hpp"
#include "cosmosim/physics/cooling_heating.hpp"

int main() {
  constexpr std::size_t k_cells = 1 << 18;
  constexpr std::size_t k_iterations = 20;
  constexpr double k_gamma = 5.0 / 3.0;

  cosmosim::hydro::HydroConservedStateSoa conserved(k_cells);
  for (std::size_t i = 0; i < k_cells; ++i) {
    cosmosim::hydro::HydroPrimitiveState primitive{};
    primitive.rho_comoving = 1.0;
    primitive.pressure_comoving = 0.5;
    conserved.storeCell(i, cosmosim::hydro::HydroCoreSolver::conservedFromPrimitive(primitive, k_gamma));
  }

  cosmosim::physics::CoolingModelConfig model{};
  model.max_subcycles = 8;
  cosmosim::physics::CoolingRateProvider provider(model);
  cosmosim::physics::CoolingSourceIntegrator integrator(1.0e-8);
  cosmosim::physics::CoolingHeatingSource source(provider, integrator);

  std::vector<double> n_h(k_cells, 0.1);
  std::vector<double> z(k_cells, 0.01);
  std::vector<double> t(k_cells, 1.0e6);
  cosmosim::hydro::HydroSourceContext context{};
  context.update.dt_code = 1.0e-3;
  context.update.scale_factor = 1.0;
  context.hydrogen_number_density_cgs = n_h;
  context.metallicity_mass_fraction = z;
  context.temperature_k = t;

  const auto setup_start = std::chrono::steady_clock::now();
  volatile double sink = 0.0;
  for (std::size_t i = 0; i < k_cells; ++i) {
    const auto primitive = cosmosim::hydro::HydroCoreSolver::primitiveFromConserved(conserved.loadCell(i), k_gamma);
    sink += source.sourceForCell(i, conserved.loadCell(i), primitive, context).total_energy_density_comoving;
  }
  const auto setup_end = std::chrono::steady_clock::now();

  const auto steady_start = std::chrono::steady_clock::now();
  for (std::size_t iter = 0; iter < k_iterations; ++iter) {
    for (std::size_t i = 0; i < k_cells; ++i) {
      const auto primitive = cosmosim::hydro::HydroCoreSolver::primitiveFromConserved(conserved.loadCell(i), k_gamma);
      sink += source.sourceForCell(i, conserved.loadCell(i), primitive, context).total_energy_density_comoving;
    }
  }
  const auto steady_end = std::chrono::steady_clock::now();

  const double setup_ms = std::chrono::duration<double, std::milli>(setup_end - setup_start).count();
  const double steady_ms = std::chrono::duration<double, std::milli>(steady_end - steady_start).count();
  const double steady_s = steady_ms * 1.0e-3;
  const double cell_updates = static_cast<double>(k_cells * k_iterations);

  std::cout << "bench_cooling_heating_source"
            << " build_type=" << COSMOSIM_BUILD_TYPE
            << " hardware=cpu"
            << " threads=1"
            << " features=primordial_cooling+subcycling"
            << " setup_ms=" << setup_ms
            << " steady_ms=" << steady_ms
            << " cells=" << k_cells
            << " iterations=" << k_iterations
            << " cell_updates_per_s=" << (cell_updates / steady_s)
            << " effective_input_bandwidth_gb_s="
            << ((cell_updates * 4.0 * sizeof(double)) / steady_s * 1.0e-9)
            << " sink=" << sink
            << '\n';

  return 0;
}
