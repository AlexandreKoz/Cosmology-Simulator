#include <cassert>
#include <cmath>
#include <filesystem>
#include <fstream>

#include "cosmosim/core/config.hpp"
#include "cosmosim/physics/cooling_heating.hpp"

namespace {

void testMetalLineTableLookupInterpolation() {
  const std::filesystem::path table_path =
      std::filesystem::temp_directory_path() / "cosmosim_test_metal_table.txt";
  {
    std::ofstream out(table_path);
    out << "# log10T log10Lambda\n";
    out << "4.0 -23.0\n";
    out << "5.0 -22.0\n";
  }

  const auto table = cosmosim::physics::MetalLineCoolingTable::loadFromTextFile(table_path.string(), "unit_test");
  const double rate_mid = table.lookupCoolingRateErgCm3S(1.0e4 * std::sqrt(10.0));
  assert(std::abs(std::log10(rate_mid) + 22.5) < 1.0e-6);

  std::filesystem::remove(table_path);
}

void testCoolingIntegratorSubcyclesAndFloors() {
  cosmosim::physics::CoolingModelConfig model;
  model.temperature_floor_k = 100.0;
  model.max_fractional_energy_change_per_substep = 0.01;
  model.max_subcycles = 32;

  cosmosim::physics::CoolingRateProvider provider(model);
  cosmosim::physics::CoolingSourceIntegrator integrator(1.0e-6);

  const cosmosim::physics::CoolingRateQuery query{
      .temperature_k = 1.0e6,
      .hydrogen_number_density_cgs = 1.0e5,
      .metallicity_mass_fraction = 0.0,
      .redshift = 0.0};

  const auto result = integrator.integrateSpecificInternalEnergy(1.0, 1.0, 5.0, query, provider);
  assert(result.specific_internal_energy_code >= 1.0e-6);
  assert(result.diagnostics.subcycles_used >= 1);
  assert(!result.diagnostics.hit_subcycle_limit);
  assert(result.diagnostics.suggested_next_dt_code > 0.0);
}

void testCoolingHeatingSourceEnergySign() {
  cosmosim::physics::CoolingModelConfig model{};
  model.uv_background_model = cosmosim::core::UvBackgroundModel::kNone;
  cosmosim::physics::CoolingRateProvider provider(model);
  cosmosim::physics::CoolingSourceIntegrator integrator(1.0e-8);
  cosmosim::physics::CoolingHeatingSource source(provider, integrator);

  cosmosim::hydro::HydroPrimitiveState primitive{
      .rho_comoving = 1.0,
      .vel_x_peculiar = 0.0,
      .vel_y_peculiar = 0.0,
      .vel_z_peculiar = 0.0,
      .pressure_comoving = 0.1};
  const auto conserved = cosmosim::hydro::HydroCoreSolver::conservedFromPrimitive(primitive, 5.0 / 3.0);

  const std::array<double, 1> n_h{5.0};
  const std::array<double, 1> z{0.0};
  const std::array<double, 1> temp{1.0e6};
  cosmosim::hydro::HydroSourceContext context{};
  context.update.dt_code = 0.01;
  context.update.scale_factor = 1.0;
  context.hydrogen_number_density_cgs = n_h;
  context.metallicity_mass_fraction = z;
  context.temperature_k = temp;

  const auto source_state = source.sourceForCell(0, conserved, primitive, context);
  assert(source_state.total_energy_density_comoving <= 0.0);
}

}  // namespace

int main() {
  testMetalLineTableLookupInterpolation();
  testCoolingIntegratorSubcyclesAndFloors();
  testCoolingHeatingSourceEnergySign();
  return 0;
}
