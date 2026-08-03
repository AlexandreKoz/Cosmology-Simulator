#include <cassert>
#include <cmath>
#include <string>

#include "cosmosim/physics/stellar_evolution.hpp"

namespace {

std::string tablePath() {
  return std::string(COSMOSIM_SOURCE_DIR) +
      "/resources/stellar_evolution/test_synthetic_v2.txt";
}

void testSafeBuiltinContainsNoInventedYields() {
  const auto table = cosmosim::physics::StellarEvolutionTable::makeBuiltinReference();
  const auto state = table.evaluate(1.0e9, 0.02);
  assert(state.return_fraction_total == 0.0);
  assert(state.total_ejected_metal_fraction_total == 0.0);
  assert(!table.production_calibrated);
}

void testExactNodeAndMetallicityDependence() {
  const auto table = cosmosim::physics::StellarEvolutionTable::loadFromTextFile(tablePath());
  const auto metal_free = table.evaluate(1.0e8, 0.0);
  const auto metal_rich = table.evaluate(1.0e8, 0.03);
  assert(std::abs(metal_free.return_fraction_total - 0.153) < 1.0e-14);
  assert(std::abs(metal_free.newly_synthesized_metal_fraction_total - 0.0052) < 1.0e-14);
  assert(metal_rich.total_ejected_metal_fraction_total >
         metal_free.total_ejected_metal_fraction_total);
  assert(table.sha256 == "2ac5f1330e423f3e5f9735bcaacc8c1bd31c9832f41c4e4909fabf0482e6f513");
}

void testInterpolationAndIntervalAdditivity() {
  const auto table = cosmosim::physics::StellarEvolutionTable::loadFromTextFile(tablePath());
  const auto middle = table.evaluate(std::sqrt(1.0e7 * 1.0e8), 0.005);
  assert(middle.return_fraction_total > 0.09);
  assert(middle.return_fraction_total < 0.153);

  const auto whole = table.integrateInterval(1.0e7, 1.0e9, 2.5, 0.01);
  const auto first = table.integrateInterval(1.0e7, 1.0e8, 2.5, 0.01);
  const auto second = table.integrateInterval(1.0e8, 1.0e9, 2.5, 0.01);
  assert(std::abs(whole.returned_mass_code -
                  first.returned_mass_code - second.returned_mass_code) < 1.0e-13);
  assert(std::abs(whole.returned_metals_code -
                  first.returned_metals_code - second.returned_metals_code) < 1.0e-13);
  assert(std::abs(whole.feedback_energy_erg -
                  first.feedback_energy_erg - second.feedback_energy_erg) < 1.0e34);
  const auto zero = table.integrateInterval(1.0e8, 1.0e8, 2.5, 0.01);
  assert(zero.returned_mass_code == 0.0);
  assert(zero.returned_metals_code == 0.0);
}

void testChannelConservation() {
  const auto table = cosmosim::physics::StellarEvolutionTable::loadFromTextFile(tablePath());
  const auto interval = table.integrateInterval(0.0, 1.4e10, 1.0, 0.03);
  double returned_sum = 0.0;
  double metal_sum = 0.0;
  for (std::size_t channel = 0; channel < 3U; ++channel) {
    returned_sum += interval.returned_mass_channel_code[channel];
    metal_sum += interval.returned_metals_channel_code[channel];
  }
  assert(std::abs(interval.returned_mass_code - returned_sum) < 1.0e-13);
  assert(std::abs(interval.returned_metals_code - metal_sum) < 1.0e-13);
  assert(interval.returned_metals_code <= interval.returned_mass_code);
}

}  // namespace

int main() {
  testSafeBuiltinContainsNoInventedYields();
  testExactNodeAndMetallicityDependence();
  testInterpolationAndIntervalAdditivity();
  testChannelConservation();
  return 0;
}
