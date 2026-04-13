#include <cassert>
#include <cmath>

#include "cosmosim/physics/stellar_evolution.hpp"

namespace {

void testLifetimeInterpolationMidpoint() {
  const auto table = cosmosim::physics::StellarEvolutionTable::makeBuiltinReference();
  const auto state = table.evaluateAtAgeYears(5.0e7);
  assert(std::abs(state.return_fraction_total - 0.13555555555555557) < 1.0e-9);
  assert(state.return_fraction_channel[0] > 0.0);
}

void testIntervalConservationAndChannelConsistency() {
  const auto table = cosmosim::physics::StellarEvolutionTable::makeBuiltinReference();
  const double initial_mass = 1.0;
  const auto interval = table.integrateInterval(1.0e7, 1.0e9, initial_mass);

  const double returned_channel_sum =
      interval.returned_mass_channel_code[0] + interval.returned_mass_channel_code[1] +
      interval.returned_mass_channel_code[2];
  const double metal_channel_sum =
      interval.returned_metals_channel_code[0] + interval.returned_metals_channel_code[1] +
      interval.returned_metals_channel_code[2];

  assert(std::abs(interval.returned_mass_code - returned_channel_sum) < 1.0e-12);
  assert(std::abs(interval.returned_metals_code - metal_channel_sum) < 1.0e-12);
  assert(interval.feedback_energy_erg > 0.0);
}

}  // namespace

int main() {
  testLifetimeInterpolationMidpoint();
  testIntervalConservationAndChannelConsistency();
  return 0;
}
