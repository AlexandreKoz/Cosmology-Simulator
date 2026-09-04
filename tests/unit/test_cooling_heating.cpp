#include <cassert>
#include <cmath>
#include <filesystem>
#include <fstream>

#include "cosmosim/core/config.hpp"
#include "cosmosim/physics/cooling_heating.hpp"
#include "../support/test_temp_workspace.hpp"

namespace {

void testMetalLineTableLookupInterpolation() {
  const std::filesystem::path table_path =
      cosmosim::test_support::TestTempWorkspace::uniqueProcessLocalPath("cosmosim_test_metal_table.txt");
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

class StaticHydroSourcePropertyProvider final
    : public cosmosim::hydro::HydroCellSourcePropertyProvider {
 public:
  explicit StaticHydroSourcePropertyProvider(
      cosmosim::hydro::HydroCellSourceProperties properties)
      : m_properties(properties) {}

  [[nodiscard]] cosmosim::hydro::HydroCellSourceProperties sourcePropertiesForCell(
      std::size_t) const override {
    return m_properties;
  }

 private:
  cosmosim::hydro::HydroCellSourceProperties m_properties{};
};

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

  const std::array<double, 1> rho_phys_cgs{1.0};
  const std::array<double, 1> n_h{5.0};
  const std::array<double, 1> z{0.0};
  const std::array<double, 1> temp{1.0e6};
  cosmosim::hydro::HydroSourceContext context{};
  context.update.dt_code = 0.01;
  context.update.scale_factor = 1.0;
  context.mass_density_physical_cgs = rho_phys_cgs;
  context.hydrogen_number_density_cgs = n_h;
  context.metallicity_mass_fraction = z;
  context.temperature_k = temp;

  const auto source_state = source.sourceForCell(0, conserved, primitive, context);
  assert(source_state.total_energy_density_comoving <= 0.0);
}

void testCoolingHeatingSourceProviderMatchesDenseCompatibilityPath() {
  cosmosim::physics::CoolingModelConfig model{};
  model.uv_background_model = cosmosim::core::UvBackgroundModel::kNone;
  cosmosim::physics::CoolingRateProvider provider(model);
  cosmosim::physics::CoolingSourceIntegrator integrator(1.0e-8);
  cosmosim::physics::CoolingHeatingSource source(provider, integrator);

  const cosmosim::hydro::HydroPrimitiveState primitive{
      .rho_comoving = 1.0,
      .pressure_comoving = 0.1};
  const auto conserved = cosmosim::hydro::HydroCoreSolver::conservedFromPrimitive(
      primitive, 5.0 / 3.0);

  const std::array<double, 1> rho_phys_cgs{1.0};
  const std::array<double, 1> n_h{5.0};
  const std::array<double, 1> metallicity{0.02};
  const std::array<double, 1> temperature{1.0e6};
  cosmosim::hydro::HydroSourceContext dense_context{};
  dense_context.update.dt_code = 0.01;
  dense_context.update.scale_factor = 1.0;
  dense_context.mass_density_physical_cgs = rho_phys_cgs;
  dense_context.hydrogen_number_density_cgs = n_h;
  dense_context.metallicity_mass_fraction = metallicity;
  dense_context.temperature_k = temperature;

  StaticHydroSourcePropertyProvider property_provider({
      .mass_density_physical_cgs = rho_phys_cgs[0],
      .hydrogen_number_density_cgs = n_h[0],
      .metallicity_mass_fraction = metallicity[0],
      .temperature_k = temperature[0]});
  cosmosim::hydro::HydroSourceContext provider_context{};
  provider_context.update = dense_context.update;
  provider_context.source_property_provider = &property_provider;

  const auto dense = source.sourceForCell(0U, conserved, primitive, dense_context);
  const auto derived = source.sourceForCell(0U, conserved, primitive, provider_context);
  assert(std::abs(dense.total_energy_density_comoving -
                  derived.total_energy_density_comoving) < 1.0e-15);
}

}  // namespace

int main() {
  testMetalLineTableLookupInterpolation();
  testCoolingIntegratorSubcyclesAndFloors();
  testCoolingHeatingSourceEnergySign();
  testCoolingHeatingSourceProviderMatchesDenseCompatibilityPath();
  return 0;
}
