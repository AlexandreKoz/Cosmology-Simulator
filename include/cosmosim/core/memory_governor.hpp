#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <mutex>
#include <stdexcept>
#include <string_view>

namespace cosmosim::core {

enum class MemoryClass : std::uint8_t {
  kCanonicalPersistent = 0,
  kPersistentCache = 1,
  kPhaseResident = 2,
  kScratchArena = 3,
  kCommunication = 4,
  kDiagnostic = 5,
  kExternalRuntime = 6,
  kCount = 7,
};

enum class MemoryPressure : std::uint8_t {
  kGreen = 0,
  kAmber = 1,
  kRed = 2,
  kTrip = 3,
};

[[nodiscard]] std::size_t memoryClassIndex(MemoryClass memory_class) noexcept;
[[nodiscard]] std::string_view memoryClassLabel(MemoryClass memory_class) noexcept;
[[nodiscard]] std::string_view memoryPressureLabel(MemoryPressure pressure) noexcept;

// Checked helpers are public because memory-policy construction occurs at the
// workflow composition boundary as well as inside the governor itself.
[[nodiscard]] std::uint64_t checkedMemoryBytesAdd(
    std::uint64_t lhs,
    std::uint64_t rhs,
    std::string_view context);

struct MemoryGovernorPolicy {
  // Zero preserves the repository's existing unlimited process-budget
  // semantics. Accounting and high-water telemetry remain active.
  std::uint64_t hard_limit_bytes = 0U;
  // Opaque library/runtime memory that is intentionally reserved outside
  // CHUI-owned container accounting (MPI/FFTW/HDF5/allocator/backend).
  std::uint64_t external_runtime_reserve_bytes = 0U;
  // Known owned staging that policy permits to overlap the resident runtime
  // even though it is not yet a governed physical allocation in M1B.
  std::uint64_t planned_overlap_reserve_bytes = 0U;
  // Must match parallel.process_memory_safety_margin_fraction. The governor
  // applies the same multiplicative rule as the pre-run DMO process estimate.
  double safety_margin_fraction = 0.0;
  // Internal policy only; no second user-visible configuration surface.
  std::uint16_t amber_basis_points = 8500U;
  std::uint16_t red_basis_points = 9500U;
};

struct MemoryGovernorSnapshot {
  std::uint64_t hard_limit_bytes = 0U;
  std::uint64_t baseline_owned_bytes = 0U;
  std::uint64_t external_runtime_reserve_bytes = 0U;
  std::uint64_t planned_overlap_reserve_bytes = 0U;
  std::uint64_t committed_bytes = 0U;
  std::uint64_t reserved_bytes = 0U;
  // Raw accounted demand before the configured multiplicative safety margin.
  std::uint64_t accounted_bytes = 0U;
  // Accounted demand after applying the same safety-margin rule used by
  // pre-run process memory estimation.
  std::uint64_t safety_adjusted_accounted_bytes = 0U;
  // Additional raw bytes that can be admitted while preserving the safety
  // margin. Unlimited mode reports UINT64_MAX.
  std::uint64_t headroom_bytes = 0U;
  std::uint64_t peak_committed_bytes = 0U;
  std::uint64_t peak_reserved_bytes = 0U;
  std::uint64_t peak_accounted_bytes = 0U;
  std::uint64_t rejection_count = 0U;
  MemoryPressure pressure = MemoryPressure::kGreen;
  std::array<std::uint64_t, static_cast<std::size_t>(MemoryClass::kCount)>
      committed_by_class{};
  std::array<std::uint64_t, static_cast<std::size_t>(MemoryClass::kCount)>
      reserved_by_class{};
};

class MemoryGovernor;

class MemoryReservation final {
 public:
  MemoryReservation() noexcept = default;
  ~MemoryReservation() noexcept;

  MemoryReservation(const MemoryReservation&) = delete;
  MemoryReservation& operator=(const MemoryReservation&) = delete;

  MemoryReservation(MemoryReservation&& other) noexcept;
  MemoryReservation& operator=(MemoryReservation&& other) noexcept;

  [[nodiscard]] bool valid() const noexcept;
  [[nodiscard]] bool pending() const noexcept;
  [[nodiscard]] bool committed() const noexcept;
  [[nodiscard]] std::uint64_t bytes() const noexcept;
  [[nodiscard]] MemoryClass memoryClass() const noexcept;

  // Transitions a pending reservation into physical commitment. Calling this
  // more than once is a logic error.
  void commit();
  // Atomically replaces the governor baseline with the caller-measured retained
  // capacity and releases this committed reservation under the same governor
  // lock. This transfers accounting authority without a transient double-count
  // or unaccounted gap when phase allocations become retained cache/state.
  void reconcileBaselineOwnedAndRelease(std::uint64_t baseline_owned_bytes);
  // Explicit release is idempotent. Destruction after release is harmless.
  void release() noexcept;

 private:
  friend class MemoryGovernor;
  enum class State : std::uint8_t {
    kInert = 0,
    kPending = 1,
    kCommitted = 2,
    kReleased = 3,
  };

  MemoryReservation(
      MemoryGovernor* governor,
      MemoryClass memory_class,
      std::uint64_t bytes) noexcept;

  MemoryGovernor* m_governor = nullptr;
  MemoryClass m_memory_class = MemoryClass::kCanonicalPersistent;
  std::uint64_t m_bytes = 0U;
  State m_state = State::kInert;
};

class MemoryGovernor final {
 public:
  explicit MemoryGovernor(MemoryGovernorPolicy policy = {});

  MemoryGovernor(const MemoryGovernor&) = delete;
  MemoryGovernor& operator=(const MemoryGovernor&) = delete;
  MemoryGovernor(MemoryGovernor&&) = delete;
  MemoryGovernor& operator=(MemoryGovernor&&) = delete;

  [[nodiscard]] MemoryReservation reserve(
      MemoryClass memory_class,
      std::uint64_t requested_bytes,
      std::string_view owner);

  // Reconciles physical memory already owned outside governor-managed
  // commitments. Callers must exclude MemoryEntry ranges flagged as governed
  // commitments, otherwise the same byte range would be double-counted.
  void setBaselineOwnedBytes(std::uint64_t baseline_owned_bytes);

  [[nodiscard]] MemoryGovernorSnapshot snapshot() const;
  [[nodiscard]] const MemoryGovernorPolicy& policy() const noexcept;

 private:
  friend class MemoryReservation;

  void commitReservation(MemoryClass memory_class, std::uint64_t bytes);
  void releasePendingReservation(MemoryClass memory_class, std::uint64_t bytes) noexcept;
  void releaseCommittedReservation(MemoryClass memory_class, std::uint64_t bytes) noexcept;
  void reconcileBaselineAndReleaseCommittedReservation(
      MemoryClass memory_class,
      std::uint64_t bytes,
      std::uint64_t baseline_owned_bytes);

  [[nodiscard]] std::uint64_t rawAccountedBytesLocked() const;
  [[nodiscard]] std::uint64_t safetyAdjustedBytes(std::uint64_t raw_bytes) const;
  [[nodiscard]] std::uint64_t headroomBytesLocked(std::uint64_t raw_bytes) const;
  [[nodiscard]] MemoryPressure pressureForAdjustedBytes(
      std::uint64_t adjusted_bytes) const noexcept;
  void updatePeakAccountedLocked();

  MemoryGovernorPolicy m_policy;
  mutable std::mutex m_mutex;
  std::uint64_t m_baseline_owned_bytes = 0U;
  std::uint64_t m_committed_bytes = 0U;
  std::uint64_t m_reserved_bytes = 0U;
  std::uint64_t m_peak_committed_bytes = 0U;
  std::uint64_t m_peak_reserved_bytes = 0U;
  std::uint64_t m_peak_accounted_bytes = 0U;
  std::uint64_t m_rejection_count = 0U;
  std::array<std::uint64_t, static_cast<std::size_t>(MemoryClass::kCount)>
      m_committed_by_class{};
  std::array<std::uint64_t, static_cast<std::size_t>(MemoryClass::kCount)>
      m_reserved_by_class{};
};

}  // namespace cosmosim::core
