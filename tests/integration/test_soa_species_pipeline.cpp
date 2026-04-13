#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <vector>
#include <utility>

#include "cosmosim/core/soa_storage.hpp"

namespace {

struct GasSpeciesContainer {
  cosmosim::core::ParticleSoaStorage particles;

  void resize(std::size_t count) { particles.resize(count); }
};

}  // namespace

int main() {
  GasSpeciesContainer gas;
  gas.resize(6);

  auto id = gas.particles.span<cosmosim::core::ParticleSoaField::kId>();
  auto pos_x = gas.particles.span<cosmosim::core::ParticleSoaField::kPosX>();
  auto pos_y = gas.particles.span<cosmosim::core::ParticleSoaField::kPosY>();
  auto pos_z = gas.particles.span<cosmosim::core::ParticleSoaField::kPosZ>();
  auto rho = gas.particles.span<cosmosim::core::ParticleSoaField::kRho>();
  auto u_int = gas.particles.span<cosmosim::core::ParticleSoaField::kUInt>();

  for (std::size_t i = 0; i < gas.particles.size(); ++i) {
    id[i] = 50U + static_cast<std::uint64_t>(i);
    pos_x[i] = static_cast<double>(i);
    pos_y[i] = static_cast<double>(i) + 0.5;
    pos_z[i] = static_cast<double>(i) + 1.0;
    rho[i] = 1.0 + static_cast<double>(i) * 0.25;
    u_int[i] = 0.2 + static_cast<double>(i) * 0.1;
  }

  const std::array<std::uint32_t, 3> active{1, 3, 5};
  std::array<double, 3> rho_active{};
  std::array<double, 3> u_int_active{};
  cosmosim::core::gatherSpan<double>(std::as_const(gas.particles).span<cosmosim::core::ParticleSoaField::kRho>(),
                                     active,
                                     rho_active);
  cosmosim::core::gatherSpan<double>(
      std::as_const(gas.particles).span<cosmosim::core::ParticleSoaField::kUInt>(),
      active,
      u_int_active);

  double cooling_checksum = 0.0;
  for (std::size_t i = 0; i < active.size(); ++i) {
    cooling_checksum += rho_active[i] * u_int_active[i];
  }

  assert(cooling_checksum > 0.0);

  const std::array<std::uint8_t, 6> keep_flags{1, 1, 0, 1, 0, 1};
  const auto compact_count = gas.particles.stableCompact(keep_flags);
  assert(compact_count == 4);

  auto id_compact = gas.particles.span<cosmosim::core::ParticleSoaField::kId>();
  assert(id_compact[0] == 50);
  assert(id_compact[1] == 51);
  assert(id_compact[2] == 53);
  assert(id_compact[3] == 55);

  return 0;
}
