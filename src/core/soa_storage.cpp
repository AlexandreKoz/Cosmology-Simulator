#include "cosmosim/core/soa_storage.hpp"

#include <algorithm>
#include <stdexcept>

namespace cosmosim::core {

void ParticleSoaStorage::resize(std::size_t count) {
  // Maintain lock-step logical size across all field lanes.
  m_pos_x.resize(count);
  m_pos_y.resize(count);
  m_pos_z.resize(count);
  m_vel_x.resize(count);
  m_vel_y.resize(count);
  m_vel_z.resize(count);
  m_mass.resize(count);
  m_id.resize(count);
  m_rho.resize(count);
  m_u_int.resize(count);
}

void ParticleSoaStorage::reserve(std::size_t count) {
  // Reserve all lanes together to keep amortized growth behavior aligned.
  m_pos_x.reserve(count);
  m_pos_y.reserve(count);
  m_pos_z.reserve(count);
  m_vel_x.reserve(count);
  m_vel_y.reserve(count);
  m_vel_z.reserve(count);
  m_mass.reserve(count);
  m_id.reserve(count);
  m_rho.reserve(count);
  m_u_int.reserve(count);
}

std::size_t ParticleSoaStorage::size() const noexcept { return m_pos_x.size(); }

std::size_t ParticleSoaStorage::capacity() const noexcept {
  // Effective capacity is the minimum lane capacity because all fields share
  // one logical index space.
  return std::min({m_pos_x.capacity(),
                   m_pos_y.capacity(),
                   m_pos_z.capacity(),
                   m_vel_x.capacity(),
                   m_vel_y.capacity(),
                   m_vel_z.capacity(),
                   m_mass.capacity(),
                   m_id.capacity(),
                   m_rho.capacity(),
                   m_u_int.capacity()});
}

bool ParticleSoaStorage::isConsistent() const noexcept {
  // All lanes must maintain one shared row space for safe zipped iteration.
  const std::size_t expected = m_pos_x.size();
  return m_pos_y.size() == expected && m_pos_z.size() == expected && m_vel_x.size() == expected &&
         m_vel_y.size() == expected && m_vel_z.size() == expected && m_mass.size() == expected &&
         m_id.size() == expected && m_rho.size() == expected && m_u_int.size() == expected;
}

void ParticleSoaStorage::swapErase(std::size_t index) {
  // Swap-erase must be mirrored on every lane to preserve index ownership.
  m_pos_x.swapErase(index);
  m_pos_y.swapErase(index);
  m_pos_z.swapErase(index);
  m_vel_x.swapErase(index);
  m_vel_y.swapErase(index);
  m_vel_z.swapErase(index);
  m_mass.swapErase(index);
  m_id.swapErase(index);
  m_rho.swapErase(index);
  m_u_int.swapErase(index);
}

std::size_t ParticleSoaStorage::stableCompact(std::span<const std::uint8_t> keep_flags) {
  if (keep_flags.size() != size()) {
    throw std::invalid_argument("ParticleSoaStorage.stableCompact: keep_flags size mismatch");
  }

  // Compact first lane and mirror the same keep-mask on all remaining lanes.
  // This preserves row semantics without allocating temporary buffers.
  const auto kept = m_pos_x.stableCompact(keep_flags);
  (void)m_pos_y.stableCompact(keep_flags);
  (void)m_pos_z.stableCompact(keep_flags);
  (void)m_vel_x.stableCompact(keep_flags);
  (void)m_vel_y.stableCompact(keep_flags);
  (void)m_vel_z.stableCompact(keep_flags);
  (void)m_mass.stableCompact(keep_flags);
  (void)m_id.stableCompact(keep_flags);
  (void)m_rho.stableCompact(keep_flags);
  (void)m_u_int.stableCompact(keep_flags);
  return kept;
}

}  // namespace cosmosim::core
