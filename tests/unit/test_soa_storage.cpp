#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>
#include <utility>

#include "cosmosim/core/soa_storage.hpp"

int main() {
  cosmosim::core::ParticleSoaStorage storage;
  storage.reserve(16);
  storage.resize(8);
  assert(storage.size() == 8);
  assert(storage.capacity() >= 8);
  assert(storage.isConsistent());

  auto pos_x = storage.span<cosmosim::core::ParticleSoaField::kPosX>();
  auto pos_y = storage.span<cosmosim::core::ParticleSoaField::kPosY>();
  auto ids = storage.span<cosmosim::core::ParticleSoaField::kId>();
  assert((reinterpret_cast<std::uintptr_t>(pos_x.data()) % 64U) == 0U);
  assert((reinterpret_cast<std::uintptr_t>(pos_y.data()) % 64U) == 0U);
  assert((reinterpret_cast<std::uintptr_t>(ids.data()) % 64U) == 0U);

  for (std::size_t i = 0; i < storage.size(); ++i) {
    pos_x[i] = 100.0 + static_cast<double>(i);
    pos_y[i] = 200.0 + static_cast<double>(i);
    ids[i] = 1000U + static_cast<std::uint64_t>(i);
  }

  const auto const_view = std::as_const(storage).span<cosmosim::core::ParticleSoaField::kPosX>();
  assert(const_view[3] == 103.0);

  std::vector<std::uint32_t> gather_index{6, 2, 5};
  std::array<double, 3> gathered{};
  cosmosim::core::gatherSpan<double>(const_view, gather_index, gathered);
  assert(gathered[0] == 106.0);
  assert(gathered[1] == 102.0);
  assert(gathered[2] == 105.0);

  std::array<double, 3> delta{1.5, 2.5, 3.5};
  auto pos_x_mut = storage.span<cosmosim::core::ParticleSoaField::kPosX>();
  cosmosim::core::scatterSpan<double>(delta, gather_index, pos_x_mut);
  assert(pos_x_mut[6] == 1.5);
  assert(pos_x_mut[2] == 2.5);
  assert(pos_x_mut[5] == 3.5);

  const std::array<std::uint8_t, 8> keep_flags{1, 0, 1, 0, 1, 0, 1, 0};
  const auto kept = storage.stableCompact(keep_flags);
  assert(kept == 4);
  assert(storage.size() == 4);
  auto compact_ids = storage.span<cosmosim::core::ParticleSoaField::kId>();
  assert(compact_ids[0] == 1000);
  assert(compact_ids[1] == 1002);
  assert(compact_ids[2] == 1004);
  assert(compact_ids[3] == 1006);

  storage.swapErase(1);
  assert(storage.size() == 3);
  assert(storage.isConsistent());

  return 0;
}
