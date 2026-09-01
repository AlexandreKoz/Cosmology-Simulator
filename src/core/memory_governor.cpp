#include "cosmosim/core/memory_governor.hpp"

#include <algorithm>
#include <cmath>
#include <exception>
#include <limits>
#include <sstream>
#include <string>
#include <utility>

namespace cosmosim::core {
namespace {

constexpr std::uint64_t k_basis_point_denominator = 10000U;

[[nodiscard]] bool validMemoryClass(MemoryClass memory_class) noexcept {
  return memoryClassIndex(memory_class) <
      static_cast<std::size_t>(MemoryClass::kCount);
}

[[nodiscard]] std::uint64_t scaledBasisPointThreshold(
    std::uint64_t hard_limit_bytes,
    std::uint16_t basis_points) noexcept {
  const std::uint64_t whole = hard_limit_bytes / k_basis_point_denominator;
  const std::uint64_t remainder = hard_limit_bytes % k_basis_point_denominator;
  return whole * static_cast<std::uint64_t>(basis_points) +
      (remainder * static_cast<std::uint64_t>(basis_points)) /
          k_basis_point_denominator;
}

[[nodiscard]] std::pair<std::uint64_t, std::uint64_t> pressureThresholds(
    std::uint64_t hard_limit_bytes,
    std::uint16_t amber_basis_points,
    std::uint16_t red_basis_points) noexcept {
  if (hard_limit_bytes <= 1U) {
    return {1U, 1U};
  }
  if (hard_limit_bytes == 2U) {
    return {1U, 2U};
  }
  const std::uint64_t amber_raw = scaledBasisPointThreshold(
      hard_limit_bytes, amber_basis_points);
  const std::uint64_t amber = std::clamp<std::uint64_t>(
      amber_raw, 1U, hard_limit_bytes - 2U);
  const std::uint64_t red_raw = scaledBasisPointThreshold(
      hard_limit_bytes, red_basis_points);
  const std::uint64_t red = std::clamp<std::uint64_t>(
      red_raw, amber + 1U, hard_limit_bytes - 1U);
  return {amber, red};
}

[[nodiscard]] std::string rejectionDiagnostic(
    MemoryClass memory_class,
    std::string_view owner,
    std::uint64_t requested_bytes,
    const MemoryGovernorSnapshot& snapshot,
    std::string_view reason) {
  std::ostringstream out;
  out << "memory reservation rejected"
      << " owner=" << owner
      << " class=" << memoryClassLabel(memory_class)
      << " requested_bytes=" << requested_bytes
      << " committed_bytes=" << snapshot.committed_bytes
      << " reserved_bytes=" << snapshot.reserved_bytes
      << " baseline_owned_bytes=" << snapshot.baseline_owned_bytes
      << " external_runtime_reserve_bytes=" << snapshot.external_runtime_reserve_bytes
      << " planned_overlap_reserve_bytes=" << snapshot.planned_overlap_reserve_bytes
      << " accounted_bytes=" << snapshot.accounted_bytes
      << " hard_limit_bytes=" << snapshot.hard_limit_bytes
      << " headroom_bytes=" << snapshot.headroom_bytes
      << " pressure=" << memoryPressureLabel(snapshot.pressure)
      << " reason=" << reason;
  return out.str();
}

}  // namespace

std::size_t memoryClassIndex(MemoryClass memory_class) noexcept {
  return static_cast<std::size_t>(memory_class);
}

std::string_view memoryClassLabel(MemoryClass memory_class) noexcept {
  switch (memory_class) {
    case MemoryClass::kCanonicalPersistent:
      return "canonical_persistent";
    case MemoryClass::kPersistentCache:
      return "persistent_cache";
    case MemoryClass::kPhaseResident:
      return "phase_resident";
    case MemoryClass::kScratchArena:
      return "scratch_arena";
    case MemoryClass::kCommunication:
      return "communication";
    case MemoryClass::kDiagnostic:
      return "diagnostic";
    case MemoryClass::kExternalRuntime:
      return "external_runtime";
    case MemoryClass::kCount:
      return "invalid";
  }
  return "invalid";
}

std::string_view memoryPressureLabel(MemoryPressure pressure) noexcept {
  switch (pressure) {
    case MemoryPressure::kGreen:
      return "green";
    case MemoryPressure::kAmber:
      return "amber";
    case MemoryPressure::kRed:
      return "red";
    case MemoryPressure::kTrip:
      return "trip";
  }
  return "trip";
}

std::uint64_t checkedMemoryBytesAdd(
    std::uint64_t lhs,
    std::uint64_t rhs,
    std::string_view context) {
  if (rhs > std::numeric_limits<std::uint64_t>::max() - lhs) {
    throw std::overflow_error(std::string(context) + ": uint64 byte addition overflow");
  }
  return lhs + rhs;
}

MemoryReservation::MemoryReservation(
    MemoryGovernor* governor,
    MemoryClass memory_class,
    std::uint64_t bytes) noexcept
    : m_governor(governor),
      m_memory_class(memory_class),
      m_bytes(bytes),
      m_state(State::kPending) {}

MemoryReservation::~MemoryReservation() noexcept {
  release();
}

MemoryReservation::MemoryReservation(MemoryReservation&& other) noexcept
    : m_governor(std::exchange(other.m_governor, nullptr)),
      m_memory_class(other.m_memory_class),
      m_bytes(std::exchange(other.m_bytes, 0U)),
      m_state(std::exchange(other.m_state, State::kInert)) {}

MemoryReservation& MemoryReservation::operator=(MemoryReservation&& other) noexcept {
  if (this == &other) {
    return *this;
  }
  release();
  m_governor = std::exchange(other.m_governor, nullptr);
  m_memory_class = other.m_memory_class;
  m_bytes = std::exchange(other.m_bytes, 0U);
  m_state = std::exchange(other.m_state, State::kInert);
  return *this;
}

bool MemoryReservation::valid() const noexcept {
  return m_state == State::kPending || m_state == State::kCommitted;
}

bool MemoryReservation::pending() const noexcept {
  return m_state == State::kPending;
}

bool MemoryReservation::committed() const noexcept {
  return m_state == State::kCommitted;
}

std::uint64_t MemoryReservation::bytes() const noexcept {
  return m_bytes;
}

MemoryClass MemoryReservation::memoryClass() const noexcept {
  return m_memory_class;
}

void MemoryReservation::commit() {
  if (m_state != State::kPending || m_governor == nullptr) {
    throw std::logic_error("MemoryReservation.commit requires one live pending reservation");
  }
  m_governor->commitReservation(m_memory_class, m_bytes);
  m_state = State::kCommitted;
}

void MemoryReservation::reconcileBaselineOwnedAndRelease(
    std::uint64_t baseline_owned_bytes) {
  if (m_state != State::kCommitted || m_governor == nullptr) {
    throw std::logic_error(
        "MemoryReservation.reconcileBaselineOwnedAndRelease requires one live committed reservation");
  }
  m_governor->reconcileBaselineAndReleaseCommittedReservation(
      m_memory_class, m_bytes, baseline_owned_bytes);
  m_state = State::kReleased;
  m_governor = nullptr;
  m_bytes = 0U;
}

void MemoryReservation::release() noexcept {
  if (m_governor == nullptr || m_state == State::kInert || m_state == State::kReleased) {
    m_state = State::kReleased;
    m_governor = nullptr;
    m_bytes = 0U;
    return;
  }
  if (m_state == State::kPending) {
    m_governor->releasePendingReservation(m_memory_class, m_bytes);
  } else if (m_state == State::kCommitted) {
    m_governor->releaseCommittedReservation(m_memory_class, m_bytes);
  }
  m_state = State::kReleased;
  m_governor = nullptr;
  m_bytes = 0U;
}

MemoryGovernor::MemoryGovernor(MemoryGovernorPolicy policy)
    : m_policy(policy) {
  if (!std::isfinite(m_policy.safety_margin_fraction) ||
      m_policy.safety_margin_fraction < 0.0 ||
      m_policy.safety_margin_fraction > 1.0) {
    throw std::invalid_argument(
        "MemoryGovernor safety margin must be finite and within [0,1]");
  }
  if (m_policy.amber_basis_points == 0U ||
      m_policy.amber_basis_points >= m_policy.red_basis_points ||
      m_policy.red_basis_points >= k_basis_point_denominator) {
    throw std::invalid_argument(
        "MemoryGovernor pressure thresholds require 0 < amber < red < 10000 basis points");
  }
  updatePeakAccountedLocked();
}

MemoryReservation MemoryGovernor::reserve(
    MemoryClass memory_class,
    std::uint64_t requested_bytes,
    std::string_view owner) {
  if (!validMemoryClass(memory_class)) {
    throw std::invalid_argument("MemoryGovernor.reserve received invalid memory class");
  }
  if (owner.empty()) {
    throw std::invalid_argument("MemoryGovernor.reserve requires a non-empty owner label");
  }

  std::lock_guard lock(m_mutex);
  const std::uint64_t current_raw = rawAccountedBytesLocked();
  std::uint64_t prospective_raw = current_raw;
  try {
    prospective_raw = checkedMemoryBytesAdd(
        current_raw, requested_bytes, "MemoryGovernor.reserve prospective demand");
  } catch (const std::overflow_error&) {
    ++m_rejection_count;
    MemoryGovernorSnapshot current;
    current.hard_limit_bytes = m_policy.hard_limit_bytes;
    current.baseline_owned_bytes = m_baseline_owned_bytes;
    current.external_runtime_reserve_bytes = m_policy.external_runtime_reserve_bytes;
    current.planned_overlap_reserve_bytes = m_policy.planned_overlap_reserve_bytes;
    current.committed_bytes = m_committed_bytes;
    current.reserved_bytes = m_reserved_bytes;
    current.accounted_bytes = current_raw;
    current.safety_adjusted_accounted_bytes = safetyAdjustedBytes(current_raw);
    current.headroom_bytes = headroomBytesLocked(current_raw);
    current.peak_committed_bytes = m_peak_committed_bytes;
    current.peak_reserved_bytes = m_peak_reserved_bytes;
    current.peak_accounted_bytes = m_peak_accounted_bytes;
    current.rejection_count = m_rejection_count;
    current.pressure = MemoryPressure::kTrip;
    current.committed_by_class = m_committed_by_class;
    current.reserved_by_class = m_reserved_by_class;
    throw std::overflow_error(rejectionDiagnostic(
        memory_class,
        owner,
        requested_bytes,
        current,
        "uint64_accounting_overflow"));
  }

  const std::uint64_t prospective_adjusted = safetyAdjustedBytes(prospective_raw);
  if (m_policy.hard_limit_bytes != 0U &&
      prospective_adjusted > m_policy.hard_limit_bytes) {
    ++m_rejection_count;
    MemoryGovernorSnapshot current;
    current.hard_limit_bytes = m_policy.hard_limit_bytes;
    current.baseline_owned_bytes = m_baseline_owned_bytes;
    current.external_runtime_reserve_bytes = m_policy.external_runtime_reserve_bytes;
    current.planned_overlap_reserve_bytes = m_policy.planned_overlap_reserve_bytes;
    current.committed_bytes = m_committed_bytes;
    current.reserved_bytes = m_reserved_bytes;
    current.accounted_bytes = current_raw;
    current.safety_adjusted_accounted_bytes = safetyAdjustedBytes(current_raw);
    current.headroom_bytes = headroomBytesLocked(current_raw);
    current.peak_committed_bytes = m_peak_committed_bytes;
    current.peak_reserved_bytes = m_peak_reserved_bytes;
    current.peak_accounted_bytes = m_peak_accounted_bytes;
    current.rejection_count = m_rejection_count;
    current.pressure = MemoryPressure::kTrip;
    current.committed_by_class = m_committed_by_class;
    current.reserved_by_class = m_reserved_by_class;
    throw std::runtime_error(rejectionDiagnostic(
        memory_class,
        owner,
        requested_bytes,
        current,
        "hard_limit_would_be_exceeded"));
  }

  const std::size_t class_index = memoryClassIndex(memory_class);
  m_reserved_bytes = checkedMemoryBytesAdd(
      m_reserved_bytes, requested_bytes, "MemoryGovernor.reserve reserved total");
  m_reserved_by_class[class_index] = checkedMemoryBytesAdd(
      m_reserved_by_class[class_index], requested_bytes,
      "MemoryGovernor.reserve class reserved total");
  m_peak_reserved_bytes = std::max(m_peak_reserved_bytes, m_reserved_bytes);
  updatePeakAccountedLocked();
  return MemoryReservation(this, memory_class, requested_bytes);
}

void MemoryGovernor::setBaselineOwnedBytes(std::uint64_t baseline_owned_bytes) {
  std::lock_guard lock(m_mutex);
  // Validate the replacement baseline before mutating state so a reconciliation
  // overflow cannot leave the governor with an unusable partially-applied value.
  std::uint64_t prospective = baseline_owned_bytes;
  prospective = checkedMemoryBytesAdd(
      prospective, m_policy.external_runtime_reserve_bytes,
      "MemoryGovernor baseline external runtime reserve");
  prospective = checkedMemoryBytesAdd(
      prospective, m_policy.planned_overlap_reserve_bytes,
      "MemoryGovernor baseline planned overlap reserve");
  prospective = checkedMemoryBytesAdd(
      prospective, m_committed_bytes, "MemoryGovernor baseline committed bytes");
  prospective = checkedMemoryBytesAdd(
      prospective, m_reserved_bytes, "MemoryGovernor baseline reserved bytes");
  m_baseline_owned_bytes = baseline_owned_bytes;
  m_peak_accounted_bytes = std::max(m_peak_accounted_bytes, prospective);
}

MemoryGovernorSnapshot MemoryGovernor::snapshot() const {
  std::lock_guard lock(m_mutex);
  const std::uint64_t raw = rawAccountedBytesLocked();
  const std::uint64_t adjusted = safetyAdjustedBytes(raw);
  MemoryGovernorSnapshot result;
  result.hard_limit_bytes = m_policy.hard_limit_bytes;
  result.baseline_owned_bytes = m_baseline_owned_bytes;
  result.external_runtime_reserve_bytes = m_policy.external_runtime_reserve_bytes;
  result.planned_overlap_reserve_bytes = m_policy.planned_overlap_reserve_bytes;
  result.committed_bytes = m_committed_bytes;
  result.reserved_bytes = m_reserved_bytes;
  result.accounted_bytes = raw;
  result.safety_adjusted_accounted_bytes = adjusted;
  result.headroom_bytes = headroomBytesLocked(raw);
  result.peak_committed_bytes = m_peak_committed_bytes;
  result.peak_reserved_bytes = m_peak_reserved_bytes;
  result.peak_accounted_bytes = m_peak_accounted_bytes;
  result.rejection_count = m_rejection_count;
  result.pressure = pressureForAdjustedBytes(adjusted);
  result.committed_by_class = m_committed_by_class;
  result.reserved_by_class = m_reserved_by_class;
  return result;
}

const MemoryGovernorPolicy& MemoryGovernor::policy() const noexcept {
  return m_policy;
}

void MemoryGovernor::commitReservation(MemoryClass memory_class, std::uint64_t bytes) {
  std::lock_guard lock(m_mutex);
  const std::size_t class_index = memoryClassIndex(memory_class);
  if (m_reserved_bytes < bytes || m_reserved_by_class[class_index] < bytes) {
    throw std::logic_error("MemoryGovernor reservation accounting underflow during commit");
  }
  m_reserved_bytes -= bytes;
  m_reserved_by_class[class_index] -= bytes;
  m_committed_bytes = checkedMemoryBytesAdd(
      m_committed_bytes, bytes, "MemoryGovernor.commit committed total");
  m_committed_by_class[class_index] = checkedMemoryBytesAdd(
      m_committed_by_class[class_index], bytes,
      "MemoryGovernor.commit class committed total");
  m_peak_committed_bytes = std::max(m_peak_committed_bytes, m_committed_bytes);
  updatePeakAccountedLocked();
}

void MemoryGovernor::releasePendingReservation(
    MemoryClass memory_class,
    std::uint64_t bytes) noexcept {
  std::lock_guard lock(m_mutex);
  const std::size_t class_index = memoryClassIndex(memory_class);
  if (m_reserved_bytes < bytes || m_reserved_by_class[class_index] < bytes) {
    std::terminate();
  }
  m_reserved_bytes -= bytes;
  m_reserved_by_class[class_index] -= bytes;
}

void MemoryGovernor::releaseCommittedReservation(
    MemoryClass memory_class,
    std::uint64_t bytes) noexcept {
  std::lock_guard lock(m_mutex);
  const std::size_t class_index = memoryClassIndex(memory_class);
  if (m_committed_bytes < bytes || m_committed_by_class[class_index] < bytes) {
    std::terminate();
  }
  m_committed_bytes -= bytes;
  m_committed_by_class[class_index] -= bytes;
}

void MemoryGovernor::reconcileBaselineAndReleaseCommittedReservation(
    MemoryClass memory_class,
    std::uint64_t bytes,
    std::uint64_t baseline_owned_bytes) {
  std::lock_guard lock(m_mutex);
  const std::size_t class_index = memoryClassIndex(memory_class);
  if (m_committed_bytes < bytes || m_committed_by_class[class_index] < bytes) {
    throw std::logic_error(
        "MemoryGovernor committed accounting underflow during baseline reconciliation");
  }

  const std::uint64_t remaining_committed = m_committed_bytes - bytes;
  std::uint64_t prospective = baseline_owned_bytes;
  prospective = checkedMemoryBytesAdd(
      prospective, m_policy.external_runtime_reserve_bytes,
      "MemoryGovernor reconciled baseline external runtime reserve");
  prospective = checkedMemoryBytesAdd(
      prospective, m_policy.planned_overlap_reserve_bytes,
      "MemoryGovernor reconciled baseline planned overlap reserve");
  prospective = checkedMemoryBytesAdd(
      prospective, remaining_committed,
      "MemoryGovernor reconciled baseline remaining committed bytes");
  prospective = checkedMemoryBytesAdd(
      prospective, m_reserved_bytes,
      "MemoryGovernor reconciled baseline reserved bytes");

  m_committed_bytes = remaining_committed;
  m_committed_by_class[class_index] -= bytes;
  m_baseline_owned_bytes = baseline_owned_bytes;
  m_peak_accounted_bytes = std::max(m_peak_accounted_bytes, prospective);
}

std::uint64_t MemoryGovernor::rawAccountedBytesLocked() const {
  std::uint64_t total = m_baseline_owned_bytes;
  total = checkedMemoryBytesAdd(
      total, m_policy.external_runtime_reserve_bytes,
      "MemoryGovernor raw external runtime reserve");
  total = checkedMemoryBytesAdd(
      total, m_policy.planned_overlap_reserve_bytes,
      "MemoryGovernor raw planned overlap reserve");
  total = checkedMemoryBytesAdd(total, m_committed_bytes, "MemoryGovernor raw committed bytes");
  total = checkedMemoryBytesAdd(total, m_reserved_bytes, "MemoryGovernor raw reserved bytes");
  return total;
}

std::uint64_t MemoryGovernor::safetyAdjustedBytes(std::uint64_t raw_bytes) const {
  if (raw_bytes == 0U || m_policy.safety_margin_fraction == 0.0) {
    return raw_bytes;
  }
  const long double scaled = static_cast<long double>(raw_bytes) *
      (1.0L + static_cast<long double>(m_policy.safety_margin_fraction));
  if (!std::isfinite(scaled) ||
      scaled > static_cast<long double>(std::numeric_limits<std::uint64_t>::max())) {
    return std::numeric_limits<std::uint64_t>::max();
  }
  return static_cast<std::uint64_t>(std::ceil(scaled));
}

std::uint64_t MemoryGovernor::headroomBytesLocked(std::uint64_t raw_bytes) const {
  if (m_policy.hard_limit_bytes == 0U) {
    return std::numeric_limits<std::uint64_t>::max();
  }
  if (safetyAdjustedBytes(raw_bytes) > m_policy.hard_limit_bytes) {
    return 0U;
  }
  if (m_policy.safety_margin_fraction == 0.0) {
    return m_policy.hard_limit_bytes - raw_bytes;
  }

  // Find the greatest raw demand whose safety-adjusted demand still fits the
  // configured process ceiling. This uses exactly the same rounding rule as
  // admission rather than relying on a lossy algebraic inversion.
  std::uint64_t low = raw_bytes;
  std::uint64_t high = m_policy.hard_limit_bytes;
  while (low < high) {
    const std::uint64_t delta = high - low;
    const std::uint64_t mid = low + delta / 2U + delta % 2U;
    if (safetyAdjustedBytes(mid) <= m_policy.hard_limit_bytes) {
      low = mid;
    } else {
      high = mid - 1U;
    }
  }
  return low - raw_bytes;
}

MemoryPressure MemoryGovernor::pressureForAdjustedBytes(
    std::uint64_t adjusted_bytes) const noexcept {
  if (m_policy.hard_limit_bytes == 0U) {
    return MemoryPressure::kGreen;
  }
  if (adjusted_bytes > m_policy.hard_limit_bytes) {
    return MemoryPressure::kTrip;
  }
  if (adjusted_bytes == 0U) {
    return MemoryPressure::kGreen;
  }
  const auto [amber, red] = pressureThresholds(
      m_policy.hard_limit_bytes,
      m_policy.amber_basis_points,
      m_policy.red_basis_points);
  if (adjusted_bytes >= red) {
    return MemoryPressure::kRed;
  }
  if (adjusted_bytes >= amber) {
    return MemoryPressure::kAmber;
  }
  return MemoryPressure::kGreen;
}

void MemoryGovernor::updatePeakAccountedLocked() {
  const std::uint64_t raw = rawAccountedBytesLocked();
  m_peak_accounted_bytes = std::max(m_peak_accounted_bytes, raw);
}

}  // namespace cosmosim::core
