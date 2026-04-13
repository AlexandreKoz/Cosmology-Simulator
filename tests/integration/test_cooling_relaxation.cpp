#include <cassert>

#include "cosmosim/hydro/hydro_core_solver.hpp"
#include "cosmosim/core/config.hpp"
#include "cosmosim/physics/cooling_heating.hpp"

int main() {
  constexpr double k_gamma = 5.0 / 3.0;

  cosmosim::hydro::HydroPrimitiveState primitive{};
  primitive.rho_comoving = 1.0;
  primitive.pressure_comoving = 1.0;
  auto conserved = cosmosim::hydro::HydroCoreSolver::conservedFromPrimitive(primitive, k_gamma);

  cosmosim::physics::CoolingModelConfig model{};
  model.temperature_floor_k = 100.0;
  model.max_subcycles = 32;
  model.uv_background_model = cosmosim::core::UvBackgroundModel::kNone;
  cosmosim::physics::CoolingRateProvider provider(model);
  cosmosim::physics::CoolingSourceIntegrator integrator(1.0e-8);
  cosmosim::physics::CoolingHeatingSource cooling_source(provider, integrator);

  std::vector<double> nh(1, 5.0);
  std::vector<double> metal(1, 0.01);
  std::vector<double> temp(1, 1.0e6);
  cosmosim::hydro::HydroSourceContext source_context{};
  source_context.update.dt_code = 0.005;
  source_context.update.scale_factor = 1.0;
  source_context.hydrogen_number_density_cgs = nh;
  source_context.metallicity_mass_fraction = metal;
  source_context.temperature_k = temp;

  const double e_initial = conserved.total_energy_density_comoving;
  for (int step = 0; step < 16; ++step) {
    const auto source = cooling_source.sourceForCell(0, conserved, primitive, source_context);
    conserved.total_energy_density_comoving += source.total_energy_density_comoving * source_context.update.dt_code;
  }

  assert(conserved.total_energy_density_comoving <= e_initial);
  return 0;
}
