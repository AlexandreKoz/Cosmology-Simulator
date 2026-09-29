#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <memory_resource>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

#include "cosmosim/core/memory_governor.hpp"

namespace cosmosim::gravity {

// One rank-local physical owner for mutually exclusive TreePM communication
// phases.  The backing allocation is retained and reused; each logical lease
// resets only the monotonic allocation cursor.  The null upstream guarantees
// that an underestimated phase cannot escape to the ordinary heap.
class GravityCommunicationArena final {
 public:
  enum class Phase : std::uint8_t {
    kIdle = 0,
    kPmDensity = 1,
    kPmHalo = 2,
    kPmInterpolation = 3,
    kTreeExchange = 4,
    kCount = 5,
  };

  class Lease final {
   public:
    Lease() noexcept = default;
    ~Lease() noexcept { release(); }

    Lease(const Lease&) = delete;
    Lease& operator=(const Lease&) = delete;

    Lease(Lease&& other) noexcept
        : m_owner(other.m_owner), m_phase(other.m_phase) {
      other.m_owner = nullptr;
      other.m_phase = Phase::kIdle;
    }

    Lease& operator=(Lease&& other) noexcept {
      if (this != &other) {
        release();
        m_owner = other.m_owner;
        m_phase = other.m_phase;
        other.m_owner = nullptr;
        other.m_phase = Phase::kIdle;
      }
      return *this;
    }

    [[nodiscard]] std::pmr::memory_resource* resource() const {
      if (m_owner == nullptr) {
        throw std::logic_error("gravity communication arena lease is not active");
      }
      return m_owner->resourceFor(m_phase);
    }

    [[nodiscard]] Phase phase() const noexcept { return m_phase; }

    void recordLogicalHighWater(std::uint64_t bytes) {
      if (m_owner == nullptr) {
        throw std::logic_error("cannot record high-water on an inactive gravity communication arena lease");
      }
      m_owner->recordLogicalHighWater(m_phase, bytes);
    }

    void release() noexcept {
      if (m_owner != nullptr) {
        m_owner->endLease(m_phase);
        m_owner = nullptr;
        m_phase = Phase::kIdle;
      }
    }

   private:
    friend class GravityCommunicationArena;
    Lease(GravityCommunicationArena* owner, Phase phase) noexcept
        : m_owner(owner), m_phase(phase) {}

    GravityCommunicationArena* m_owner = nullptr;
    Phase m_phase = Phase::kIdle;
  };

  explicit GravityCommunicationArena(
      core::MemoryGovernor* governor = nullptr,
      std::string owner = "treepm.communication_arena")
      : m_governor(governor), m_owner(std::move(owner)) {}

  GravityCommunicationArena(const GravityCommunicationArena&) = delete;
  GravityCommunicationArena& operator=(const GravityCommunicationArena&) = delete;
  GravityCommunicationArena(GravityCommunicationArena&&) = delete;
  GravityCommunicationArena& operator=(GravityCommunicationArena&&) = delete;

  void configure(std::uint64_t bytes) {
    if (bytes == 0U) {
      throw std::invalid_argument("gravity communication arena capacity must be nonzero");
    }
    if (m_phase != Phase::kIdle) {
      throw std::logic_error("cannot configure gravity communication arena while a lease is active");
    }
    if (m_capacity_bytes != 0U) {
      if (bytes > m_capacity_bytes) {
        throw std::logic_error(
            "gravity communication arena cannot grow after first allocation; configure the maximum source-derived phase requirement before first use");
      }
      return;
    }
    if (bytes > static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max())) {
      throw std::length_error("gravity communication arena exceeds size_t capacity");
    }

    core::MemoryReservation reservation =
        m_governor == nullptr
            ? core::MemoryReservation{}
            : m_governor->reserve(core::MemoryClass::kCommunication, bytes, m_owner);
    auto storage = std::make_unique_for_overwrite<std::byte[]>(
        static_cast<std::size_t>(bytes));
    auto resource = std::make_unique<std::pmr::monotonic_buffer_resource>(
        storage.get(), static_cast<std::size_t>(bytes),
        std::pmr::null_memory_resource());
    if (reservation.pending()) {
      reservation.commit();
    }
    m_storage = std::move(storage);
    m_resource = std::move(resource);
    m_reservation = std::move(reservation);
    m_capacity_bytes = bytes;
  }

  [[nodiscard]] Lease begin(Phase phase) {
    if (phase == Phase::kIdle || phase == Phase::kCount) {
      throw std::invalid_argument("invalid gravity communication arena phase");
    }
    if (m_resource == nullptr || m_capacity_bytes == 0U) {
      throw std::logic_error("gravity communication arena must be configured before leasing");
    }
    if (m_phase != Phase::kIdle) {
      throw std::logic_error("gravity communication arena already has an active phase lease");
    }
    // Every prior lease has already called release(); do it again here as a
    // defensive O(1) cursor reset before publishing the new phase.
    m_resource->release();
    m_phase = phase;
    return Lease(this, phase);
  }

  [[nodiscard]] std::uint64_t capacityBytes() const noexcept {
    return m_capacity_bytes;
  }

  [[nodiscard]] bool configured() const noexcept {
    return m_capacity_bytes != 0U;
  }

  [[nodiscard]] bool governed() const noexcept {
    return m_reservation.committed();
  }

  [[nodiscard]] Phase activePhase() const noexcept { return m_phase; }

  [[nodiscard]] std::pmr::memory_resource* activeResource() {
    if (m_phase == Phase::kIdle) {
      throw std::logic_error("gravity communication arena has no active resource");
    }
    return resourceFor(m_phase);
  }

  void recordActiveLogicalHighWater(std::uint64_t bytes) {
    if (m_phase == Phase::kIdle) {
      throw std::logic_error("gravity communication arena has no active phase high-water");
    }
    recordLogicalHighWater(m_phase, bytes);
  }

  [[nodiscard]] std::uint64_t logicalHighWater(Phase phase) const noexcept {
    const std::size_t index = static_cast<std::size_t>(phase);
    return index < m_logical_high_water_bytes.size()
        ? m_logical_high_water_bytes[index]
        : 0U;
  }

 private:
  [[nodiscard]] std::pmr::memory_resource* resourceFor(Phase phase) {
    if (m_phase != phase || m_resource == nullptr) {
      throw std::logic_error("gravity communication arena resource requested outside its active lease");
    }
    return m_resource.get();
  }

  void recordLogicalHighWater(Phase phase, std::uint64_t bytes) {
    if (m_phase != phase) {
      throw std::logic_error("gravity communication arena high-water phase mismatch");
    }
    if (bytes > m_capacity_bytes) {
      throw std::logic_error("gravity communication arena logical usage exceeds physical capacity");
    }
    auto& high_water = m_logical_high_water_bytes[static_cast<std::size_t>(phase)];
    if (bytes > high_water) {
      high_water = bytes;
    }
  }

  void endLease(Phase phase) noexcept {
    if (m_phase != phase || m_resource == nullptr) {
      // Lease misuse is a programming error.  Destructors cannot throw, so put
      // the arena into a fail-closed non-active state; later begin() still
      // validates that the resource/capacity are intact.
      m_phase = Phase::kIdle;
      return;
    }
    // All PMR containers and MPI operations using the lease must be destroyed /
    // completed before Lease destruction.  release() invalidates every pointer
    // into the phase storage in O(1) without returning the backing allocation.
    m_resource->release();
    m_phase = Phase::kIdle;
  }

  core::MemoryGovernor* m_governor = nullptr;
  std::string m_owner;
  core::MemoryReservation m_reservation;
  std::unique_ptr<std::byte[]> m_storage;
  std::unique_ptr<std::pmr::monotonic_buffer_resource> m_resource;
  std::uint64_t m_capacity_bytes = 0U;
  Phase m_phase = Phase::kIdle;
  std::array<std::uint64_t, static_cast<std::size_t>(Phase::kCount)>
      m_logical_high_water_bytes{};
};

}  // namespace cosmosim::gravity
