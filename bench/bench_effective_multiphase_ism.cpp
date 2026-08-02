#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <vector>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/core/config.hpp"
#include "cosmosim/core/units.hpp"
#include "cosmosim/hydro/hydro_core_solver.hpp"
#include "cosmosim/physics/effective_multiphase_ism.hpp"

namespace {

using Clock = std::chrono::steady_clock;

cosmosim::physics::EffectiveMultiphaseEosTable makeTable() {
  cosmosim::core::PhysicsConfig physics;
  physics.star_formation_model =
      cosmosim::core::StarFormationModelKind::kEffectiveMultiphaseTngLike;
  physics.sf_effective_eos_table_bins = 256U;
  physics.sf_effective_eos_max_density_ratio = 1.0e6;
  physics.sf_effective_q_eos = 0.3;
  return cosmosim::physics::EffectiveMultiphaseEosTable(
      cosmosim::physics::makeEffectiveMultiphaseEosConfig(physics),
      cosmosim::core::makeUnitSystem("kpc", "msun", "km_s"),
      cosmosim::physics::makeEffectiveIsmReferenceCoolingProvider(physics));
}

void report(std::string_view phase, std::uint64_t operations, double seconds, double checksum) {
  seconds = std::max(seconds, 1.0e-12);
  std::cout << "bench_effective_multiphase_ism"
            << " phase=" << phase
            << " build_type=" << COSMOSIM_BUILD_TYPE
            << " hardware=cpu threads=1"
            << " operations=" << operations
            << " elapsed_s=" << seconds
            << " operations_per_second=" << static_cast<double>(operations) / seconds
            << " checksum=" << checksum << '\n';
}

}  // namespace

int main() {
  constexpr std::uint64_t operations = 1000000U;
  const auto init_start = Clock::now();
  auto table = makeTable();
  const auto init_stop = Clock::now();
  report("table_initialization", table.entries().size(),
         std::chrono::duration<double>(init_stop - init_start).count(),
         static_cast<double>(table.tableHash()));

  std::vector<double> densities;
  densities.reserve(operations);
  const double threshold = table.thresholdDensityPhysCode();
  for (std::uint64_t i = 0; i < operations; ++i) {
    const double f = static_cast<double>(i % 10000U) / 9999.0;
    densities.push_back(threshold * std::exp(f * std::log(1.0e6)));
  }

  double checksum = 0.0;
  auto start = Clock::now();
  for (double density : densities) {
    const auto value = table.lookup(density);
    checksum += value.entry.pressure_phys_code + value.entry.cold_mass_fraction;
  }
  auto stop = Clock::now();
  report("table_lookup", operations, std::chrono::duration<double>(stop - start).count(), checksum);

  checksum = 0.0;
  start = Clock::now();
  for (double density : densities) {
    const auto value = table.evaluateDirect(density);
    checksum += value.pressure_phys_code + value.cold_mass_fraction;
  }
  stop = Clock::now();
  report("direct_equilibrium", operations, std::chrono::duration<double>(stop - start).count(), checksum);

  cosmosim::physics::EffectiveIsmThermodynamicClosure closure(std::move(table));
  cosmosim::hydro::HydroPrimitiveState primitive;
  primitive.rho_comoving = threshold * 10.0;
  primitive.specific_internal_energy_code = 1.0;
  primitive.pressure_comoving = (5.0 / 3.0 - 1.0) * primitive.rho_comoving;
  primitive.signal_speed_squared_code = (5.0 / 3.0) * primitive.pressure_comoving /
      primitive.rho_comoving;
  const auto conserved = cosmosim::hydro::HydroCoreSolver::conservedFromPrimitive(
      primitive, 5.0 / 3.0);
  checksum = 0.0;
  start = Clock::now();
  for (std::uint64_t i = 0; i < operations; ++i) {
    const auto value = closure.evaluate(0U, conserved, primitive, 1.0, 0.0);
    checksum += value.pressure_comoving + value.signal_speed_squared_code;
  }
  stop = Clock::now();
  report("hydro_closure", operations, std::chrono::duration<double>(stop - start).count(), checksum);
  return 0;
}
