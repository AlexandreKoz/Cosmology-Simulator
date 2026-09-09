#pragma once

#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <memory_resource>
#include <stdexcept>
#include <string_view>

#include "cosmosim/core/memory_governor.hpp"

namespace cosmosim::core {

// One physical allocation and one governor lease per phase. All descendants
// use the bounded PMR resource, with no governor locking in allocation loops.
// The caller must destroy PMR containers before the arena. The null upstream
// makes an underestimated workspace fail rather than silently grow the heap.
class GovernedScratchArena final {
 public:
  GovernedScratchArena(MemoryGovernor* governor, MemoryClass memory_class,
                       std::uint64_t bytes, std::string_view owner)
      : m_reservation(governor == nullptr ? MemoryReservation{} :
            governor->reserve(memory_class, bytes, owner)),
        m_storage(allocate(bytes)),
        m_resource(m_storage.get(), static_cast<std::size_t>(bytes),
                   std::pmr::null_memory_resource()),
        m_bytes(bytes) {
    if (m_reservation.pending()) m_reservation.commit();
  }

  GovernedScratchArena(const GovernedScratchArena&) = delete;
  GovernedScratchArena& operator=(const GovernedScratchArena&) = delete;
  GovernedScratchArena(GovernedScratchArena&&) = delete;
  GovernedScratchArena& operator=(GovernedScratchArena&&) = delete;

  [[nodiscard]] std::pmr::memory_resource* resource() noexcept { return &m_resource; }
  [[nodiscard]] std::uint64_t capacityBytes() const noexcept { return m_bytes; }

 private:
  [[nodiscard]] static std::unique_ptr<std::byte[]> allocate(std::uint64_t bytes) {
    if (bytes > std::numeric_limits<std::size_t>::max()) {
      throw std::length_error("governed scratch exceeds size_t capacity");
    }
    if (bytes == 0U) return {};
    return std::make_unique_for_overwrite<std::byte[]>(static_cast<std::size_t>(bytes));
  }

  MemoryReservation m_reservation;
  std::unique_ptr<std::byte[]> m_storage;
  std::pmr::monotonic_buffer_resource m_resource;
  std::uint64_t m_bytes = 0U;
};

}  // namespace cosmosim::core
