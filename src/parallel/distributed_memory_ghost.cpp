#include "cosmosim/parallel/distributed_memory.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <iomanip>
#include <limits>
#include <numeric>
#include <optional>
#include <sstream>
#include <streambuf>
#include <unordered_map>
#include <unordered_set>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/core/memory_governor.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "parallel/internal/distributed_memory_support.hpp"

#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
#include <mpi.h>
#endif

namespace cosmosim::parallel {
namespace {

using internal::ghostExchangeRecordBytes;
using internal::laneIsPresentOrEmpty;
using internal::optionalLaneValue;

template <typename T>
void appendPod(std::vector<std::uint8_t>& bytes, const T& value) {
  const auto* ptr = reinterpret_cast<const std::uint8_t*>(&value);
  bytes.insert(bytes.end(), ptr, ptr + sizeof(T));
}

template <typename T>
[[nodiscard]] T readPod(const std::vector<std::uint8_t>& bytes, std::size_t* offset) {
  if (*offset > bytes.size() || sizeof(T) > bytes.size() - *offset) {
    throw std::runtime_error("distributed-memory byte payload is truncated");
  }
  T value{};
  std::memcpy(&value, bytes.data() + *offset, sizeof(T));
  *offset += sizeof(T);
  return value;
}

constexpr std::uint16_t kGhostLanePositionX = 1U << 0U;
constexpr std::uint16_t kGhostLanePositionY = 1U << 1U;
constexpr std::uint16_t kGhostLanePositionZ = 1U << 2U;
constexpr std::uint16_t kGhostLaneMass = 1U << 3U;
constexpr std::uint16_t kGhostLaneDensity = 1U << 4U;
constexpr std::uint16_t kGhostLaneVelocityX = 1U << 5U;
constexpr std::uint16_t kGhostLaneVelocityY = 1U << 6U;
constexpr std::uint16_t kGhostLaneVelocityZ = 1U << 7U;
constexpr std::uint16_t kGhostLanePressure = 1U << 8U;
constexpr std::uint16_t kGhostLaneInternalEnergy = 1U << 9U;
constexpr std::uint16_t kGhostKnownLaneMask = (1U << 10U) - 1U;

template <typename SourceView>
[[nodiscard]] std::uint16_t ghostOptionalLaneMask(const SourceView& source) noexcept {
  std::uint16_t mask = 0U;
  const auto add_if_present = [&mask](const auto& lane, std::uint16_t bit) {
    if (!lane.empty()) {
      mask = static_cast<std::uint16_t>(mask | bit);
    }
  };
  add_if_present(source.position_x_comoving, kGhostLanePositionX);
  add_if_present(source.position_y_comoving, kGhostLanePositionY);
  add_if_present(source.position_z_comoving, kGhostLanePositionZ);
  add_if_present(source.mass_code, kGhostLaneMass);
  add_if_present(source.density_code, kGhostLaneDensity);
  add_if_present(source.velocity_x_code, kGhostLaneVelocityX);
  add_if_present(source.velocity_y_code, kGhostLaneVelocityY);
  add_if_present(source.velocity_z_code, kGhostLaneVelocityZ);
  add_if_present(source.pressure_code, kGhostLanePressure);
  add_if_present(source.internal_energy_code, kGhostLaneInternalEnergy);
  return mask;
}

[[nodiscard]] double optionalLaneValue(
    std::span<const double> lane,
    std::size_t index) {
  return lane.empty() ? 0.0 : lane[index];
}

}  // namespace
std::size_t ghostRefreshPayloadRecordBytes() noexcept {
  return ghostExchangeRecordBytes();
}

void validateGhostRefreshPayloadDescriptor(const GhostTransferDescriptor& descriptor) {
  if (descriptor.intent != GhostTransferIntent::kGhostRefreshRequest &&
      descriptor.intent != GhostTransferIntent::kGhostRefreshReceiveStaging) {
    throw std::invalid_argument(
        "GhostExchangeBuffer carries ghost-refresh payloads only; ownership migration must use ParticleMigrationRecord");
  }
  if (descriptor.expected_post_transfer_residency != LocalIndexResidency::kGhost) {
    throw std::invalid_argument("ghost-refresh payload descriptor must produce ghost residency");
  }
}

int ghostExchangePairStableTag(int tag_base, int local_rank, int peer_rank) {
  if (tag_base < 0 || local_rank < 0 || peer_rank < 0 || local_rank == peer_rank) {
    throw std::invalid_argument("ghostExchangePairStableTag: invalid rank pair or tag base");
  }
  // MPI receive matching constrains source and destination ranks. Keep this tag
  // independent of local neighbor-slot order; sequence separation is provided by
  // ghostExchangeSequencedTag for overlapping or repeated phases.
  return tag_base;
}

int ghostExchangeSequencedTag(
    int tag_base,
    int local_rank,
    int peer_rank,
    std::uint64_t exchange_sequence) {
  constexpr int k_sequence_stride = 16;
  constexpr int k_sequence_window = 64;
  if (tag_base < 0) {
    throw std::invalid_argument("ghostExchangeSequencedTag: tag_base must be non-negative");
  }
  const int phased_base = tag_base +
      static_cast<int>(exchange_sequence % static_cast<std::uint64_t>(k_sequence_window)) * k_sequence_stride;
  return ghostExchangePairStableTag(phased_base, local_rank, peer_rank);
}

bool GhostExchangeBufferSoA::isConsistent() const noexcept {
  const std::size_t n = entity_id.size();
  return laneIsPresentOrEmpty(position_x_comoving.size(), n) &&
      laneIsPresentOrEmpty(position_y_comoving.size(), n) &&
      laneIsPresentOrEmpty(position_z_comoving.size(), n) &&
      laneIsPresentOrEmpty(mass_code.size(), n) &&
      laneIsPresentOrEmpty(density_code.size(), n) &&
      laneIsPresentOrEmpty(velocity_x_code.size(), n) &&
      laneIsPresentOrEmpty(velocity_y_code.size(), n) &&
      laneIsPresentOrEmpty(velocity_z_code.size(), n) &&
      laneIsPresentOrEmpty(pressure_code.size(), n) &&
      laneIsPresentOrEmpty(internal_energy_code.size(), n);
}

std::size_t GhostExchangeBufferSoA::size() const noexcept { return entity_id.size(); }

bool GhostExchangeBufferSoA::hasGravityPayload() const noexcept {
  const std::size_t n = entity_id.size();
  return position_x_comoving.size() == n && position_y_comoving.size() == n &&
      position_z_comoving.size() == n && mass_code.size() == n;
}

bool GhostExchangeBufferSoA::hasHydroPayload() const noexcept {
  const std::size_t n = entity_id.size();
  return n != 0U && density_code.size() == n && velocity_x_code.size() == n && velocity_y_code.size() == n &&
      velocity_z_code.size() == n && pressure_code.size() == n && internal_energy_code.size() == n;
}

std::size_t ReadOnlyGhostExchangeView::size() const noexcept { return entity_id.size(); }

bool ReadOnlyGhostExchangeView::isConsistent() const noexcept {
  const std::size_t n = entity_id.size();
  return laneIsPresentOrEmpty(position_x_comoving.size(), n) &&
      laneIsPresentOrEmpty(position_y_comoving.size(), n) &&
      laneIsPresentOrEmpty(position_z_comoving.size(), n) &&
      laneIsPresentOrEmpty(mass_code.size(), n) &&
      laneIsPresentOrEmpty(density_code.size(), n) &&
      laneIsPresentOrEmpty(velocity_x_code.size(), n) &&
      laneIsPresentOrEmpty(velocity_y_code.size(), n) &&
      laneIsPresentOrEmpty(velocity_z_code.size(), n) &&
      laneIsPresentOrEmpty(pressure_code.size(), n) &&
      laneIsPresentOrEmpty(internal_energy_code.size(), n);
}

bool ReadOnlyGhostExchangeView::isFresh(const GhostLayerEpoch& expected_epoch) const noexcept {
  return epoch.matches(expected_epoch);
}

ReadOnlyGhostExchangeView makeReadOnlyGhostExchangeView(const GhostExchangeBufferSoA& storage) {
  if (!storage.isConsistent()) {
    throw std::invalid_argument("ghost storage must be component-consistent before building a read-only view");
  }
  return ReadOnlyGhostExchangeView{
      .epoch = storage.epoch,
      .entity_id = storage.entity_id,
      .position_x_comoving = storage.position_x_comoving,
      .position_y_comoving = storage.position_y_comoving,
      .position_z_comoving = storage.position_z_comoving,
      .mass_code = storage.mass_code,
      .density_code = storage.density_code,
      .velocity_x_code = storage.velocity_x_code,
      .velocity_y_code = storage.velocity_y_code,
      .velocity_z_code = storage.velocity_z_code,
      .pressure_code = storage.pressure_code,
      .internal_energy_code = storage.internal_energy_code,
  };
}

void requireFreshGhostExchangeView(
    const ReadOnlyGhostExchangeView& view,
    const GhostLayerEpoch& expected_epoch) {
  if (!view.isConsistent()) {
    throw std::invalid_argument("read-only ghost view component sizes are inconsistent");
  }
  if (!view.isFresh(expected_epoch)) {
    throw std::invalid_argument("read-only ghost view is stale for the current exchange epoch");
  }
}

void GhostExchangeBuffer::clear() { m_bytes.clear(); }

std::size_t GhostExchangeBuffer::byteSize() const noexcept { return m_bytes.size(); }

std::span<const std::uint8_t> GhostExchangeBuffer::encodedBytes() const noexcept { return m_bytes; }

void GhostExchangeBuffer::replaceEncodedBytes(std::vector<std::uint8_t> bytes) {
  m_bytes = std::move(bytes);
}

void GhostExchangeBuffer::packFrom(const GhostExchangeBufferSoA& source, std::span<const std::uint32_t> local_indices) {
  packFrom(makeReadOnlyGhostExchangeView(source), local_indices);
}

void GhostExchangeBuffer::packFrom(
    const ReadOnlyGhostExchangeView& source,
    std::span<const std::uint32_t> local_indices) {
  if (!source.isConsistent()) {
    throw std::invalid_argument("ghost source view fields must have matching sizes");
  }

  m_bytes.clear();
  appendPod<std::uint64_t>(m_bytes, static_cast<std::uint64_t>(local_indices.size()));
  appendPod<std::uint16_t>(m_bytes, ghostOptionalLaneMask(source));

  for (const std::uint32_t index : local_indices) {
    if (index >= source.size()) {
      throw std::out_of_range("ghost pack local index out of range");
    }

    appendPod<std::uint64_t>(m_bytes, source.entity_id[index]);
    appendPod<double>(m_bytes, optionalLaneValue(source.position_x_comoving, index));
    appendPod<double>(m_bytes, optionalLaneValue(source.position_y_comoving, index));
    appendPod<double>(m_bytes, optionalLaneValue(source.position_z_comoving, index));
    appendPod<double>(m_bytes, optionalLaneValue(source.mass_code, index));
    appendPod<double>(m_bytes, optionalLaneValue(source.density_code, index));
    appendPod<double>(m_bytes, optionalLaneValue(source.velocity_x_code, index));
    appendPod<double>(m_bytes, optionalLaneValue(source.velocity_y_code, index));
    appendPod<double>(m_bytes, optionalLaneValue(source.velocity_z_code, index));
    appendPod<double>(m_bytes, optionalLaneValue(source.pressure_code, index));
    appendPod<double>(m_bytes, optionalLaneValue(source.internal_energy_code, index));
  }
}

void GhostExchangeBuffer::packFrom(
    const GhostTransferDescriptor& descriptor,
    const GhostExchangeBufferSoA& source,
    std::span<const std::uint32_t> local_indices) {
  packFrom(descriptor, makeReadOnlyGhostExchangeView(source), local_indices);
}

void GhostExchangeBuffer::packFrom(
    const GhostTransferDescriptor& descriptor,
    const ReadOnlyGhostExchangeView& source,
    std::span<const std::uint32_t> local_indices) {
  validateGhostRefreshPayloadDescriptor(descriptor);
  if (descriptor.local_indices.size() != local_indices.size() ||
      !std::equal(descriptor.local_indices.begin(), descriptor.local_indices.end(), local_indices.begin())) {
    throw std::invalid_argument("ghost descriptor indices must match packed local indices");
  }
  packFrom(source, local_indices);
}

void GhostExchangeBuffer::unpackAppendTo(GhostExchangeBufferSoA& destination) const {
  if (!destination.isConsistent()) {
    throw std::invalid_argument("ghost destination SoA fields must have matching sizes");
  }
  constexpr std::size_t header_bytes = sizeof(std::uint64_t) + sizeof(std::uint16_t);
  if (m_bytes.size() < header_bytes) {
    throw std::runtime_error("ghost buffer is too small");
  }

  std::size_t offset = 0;
  const std::uint64_t count = readPod<std::uint64_t>(m_bytes, &offset);
  const std::uint16_t lane_mask = readPod<std::uint16_t>(m_bytes, &offset);
  if ((lane_mask & static_cast<std::uint16_t>(~kGhostKnownLaneMask)) != 0U) {
    throw std::runtime_error("ghost buffer contains an unknown optional-lane presence bit");
  }
  const std::uint64_t expected_payload_bytes =
      static_cast<std::uint64_t>(header_bytes) + count * static_cast<std::uint64_t>(ghostExchangeRecordBytes());
  if (expected_payload_bytes != static_cast<std::uint64_t>(m_bytes.size())) {
    throw std::runtime_error("ghost buffer payload shape does not match encoded count");
  }

  const std::size_t append_count = static_cast<std::size_t>(count);
  const std::size_t existing_count = destination.entity_id.size();
  const auto require_lane_schema = [existing_count, lane_mask](
                                       const std::vector<double>& lane,
                                       std::uint16_t bit) {
    if (existing_count == 0U) {
      return;
    }
    const bool existing_present = lane.size() == existing_count;
    const bool incoming_present = (lane_mask & bit) != 0U;
    if (existing_present != incoming_present) {
      throw std::runtime_error("ghost peer payloads disagree on optional-lane presence");
    }
  };
  require_lane_schema(destination.position_x_comoving, kGhostLanePositionX);
  require_lane_schema(destination.position_y_comoving, kGhostLanePositionY);
  require_lane_schema(destination.position_z_comoving, kGhostLanePositionZ);
  require_lane_schema(destination.mass_code, kGhostLaneMass);
  require_lane_schema(destination.density_code, kGhostLaneDensity);
  require_lane_schema(destination.velocity_x_code, kGhostLaneVelocityX);
  require_lane_schema(destination.velocity_y_code, kGhostLaneVelocityY);
  require_lane_schema(destination.velocity_z_code, kGhostLaneVelocityZ);
  require_lane_schema(destination.pressure_code, kGhostLanePressure);
  require_lane_schema(destination.internal_energy_code, kGhostLaneInternalEnergy);

  destination.entity_id.reserve(destination.entity_id.size() + append_count);
  const auto reserve_if_present = [append_count, lane_mask](
                                      std::vector<double>* lane,
                                      std::uint16_t bit) {
    if ((lane_mask & bit) != 0U) {
      lane->reserve(lane->size() + append_count);
    }
  };
  reserve_if_present(&destination.position_x_comoving, kGhostLanePositionX);
  reserve_if_present(&destination.position_y_comoving, kGhostLanePositionY);
  reserve_if_present(&destination.position_z_comoving, kGhostLanePositionZ);
  reserve_if_present(&destination.mass_code, kGhostLaneMass);
  reserve_if_present(&destination.density_code, kGhostLaneDensity);
  reserve_if_present(&destination.velocity_x_code, kGhostLaneVelocityX);
  reserve_if_present(&destination.velocity_y_code, kGhostLaneVelocityY);
  reserve_if_present(&destination.velocity_z_code, kGhostLaneVelocityZ);
  reserve_if_present(&destination.pressure_code, kGhostLanePressure);
  reserve_if_present(&destination.internal_energy_code, kGhostLaneInternalEnergy);

  for (std::uint64_t i = 0; i < count; ++i) {
    destination.entity_id.push_back(readPod<std::uint64_t>(m_bytes, &offset));
    const auto append_lane = [&](std::vector<double>* lane, std::uint16_t bit) {
      const double value = readPod<double>(m_bytes, &offset);
      if ((lane_mask & bit) != 0U) {
        lane->push_back(value);
      }
    };
    append_lane(&destination.position_x_comoving, kGhostLanePositionX);
    append_lane(&destination.position_y_comoving, kGhostLanePositionY);
    append_lane(&destination.position_z_comoving, kGhostLanePositionZ);
    append_lane(&destination.mass_code, kGhostLaneMass);
    append_lane(&destination.density_code, kGhostLaneDensity);
    append_lane(&destination.velocity_x_code, kGhostLaneVelocityX);
    append_lane(&destination.velocity_y_code, kGhostLaneVelocityY);
    append_lane(&destination.velocity_z_code, kGhostLaneVelocityZ);
    append_lane(&destination.pressure_code, kGhostLanePressure);
    append_lane(&destination.internal_energy_code, kGhostLaneInternalEnergy);
  }

  if (offset != m_bytes.size()) {
    throw std::runtime_error("ghost buffer decode found trailing bytes");
  }
}

void GhostExchangeBuffer::unpackAppendTo(
    const GhostTransferDescriptor& descriptor,
    GhostExchangeBufferSoA& destination) const {
  validateGhostRefreshPayloadDescriptor(descriptor);
  if (m_bytes.size() < sizeof(std::uint64_t)) {
    throw std::runtime_error("ghost buffer is too small");
  }
  std::size_t offset = 0;
  const std::uint64_t encoded_count = readPod<std::uint64_t>(m_bytes, &offset);
  if (encoded_count != descriptor.local_indices.size()) {
    throw std::runtime_error("ghost buffer encoded count does not match receive descriptor slots");
  }
  unpackAppendTo(destination);
}


}  // namespace cosmosim::parallel
