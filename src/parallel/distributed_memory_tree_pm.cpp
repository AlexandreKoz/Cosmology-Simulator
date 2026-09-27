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

constexpr std::uint32_t k_tree_pseudo_wire_version = 1U;
constexpr std::size_t k_tree_pseudo_wire_bytes = 152U;

void appendTreeWireU32(std::vector<std::uint8_t>& bytes, std::uint32_t value) {
  for (unsigned shift = 0U; shift < 32U; shift += 8U) {
    bytes.push_back(static_cast<std::uint8_t>((value >> shift) & 0xffU));
  }
}

void appendTreeWireU64(std::vector<std::uint8_t>& bytes, std::uint64_t value) {
  for (unsigned shift = 0U; shift < 64U; shift += 8U) {
    bytes.push_back(static_cast<std::uint8_t>((value >> shift) & 0xffU));
  }
}

void appendTreeWireDouble(std::vector<std::uint8_t>& bytes, double value) {
  appendTreeWireU64(bytes, std::bit_cast<std::uint64_t>(value));
}

[[nodiscard]] std::uint32_t readTreeWireU32(std::span<const std::uint8_t> bytes, std::size_t& offset) {
  if (offset > bytes.size() || bytes.size() - offset < sizeof(std::uint32_t)) {
    throw std::runtime_error("tree pseudo-particle wire record is truncated");
  }
  std::uint32_t value = 0U;
  for (unsigned shift = 0U; shift < 32U; shift += 8U) {
    value |= static_cast<std::uint32_t>(bytes[offset++]) << shift;
  }
  return value;
}

[[nodiscard]] std::uint64_t readTreeWireU64(std::span<const std::uint8_t> bytes, std::size_t& offset) {
  if (offset > bytes.size() || bytes.size() - offset < sizeof(std::uint64_t)) {
    throw std::runtime_error("tree pseudo-particle wire record is truncated");
  }
  std::uint64_t value = 0U;
  for (unsigned shift = 0U; shift < 64U; shift += 8U) {
    value |= static_cast<std::uint64_t>(bytes[offset++]) << shift;
  }
  return value;
}

[[nodiscard]] double readTreeWireDouble(std::span<const std::uint8_t> bytes, std::size_t& offset) {
  return std::bit_cast<double>(readTreeWireU64(bytes, offset));
}

[[nodiscard, maybe_unused]] std::vector<std::uint8_t> encodeTreePseudoPackets(
    std::span<const TreePseudoParticlePacket> packets) {
  if (packets.size() > std::numeric_limits<std::size_t>::max() / k_tree_pseudo_wire_bytes) {
    throw std::overflow_error("tree pseudo-particle wire payload size overflows size_t");
  }
  std::vector<std::uint8_t> bytes;
  bytes.reserve(packets.size() * k_tree_pseudo_wire_bytes);
  for (const TreePseudoParticlePacket& packet : packets) {
    appendTreeWireU32(bytes, packet.descriptor.wire_version);
    appendTreeWireU64(bytes, packet.descriptor.pseudo_particle_id);
    appendTreeWireU32(bytes, static_cast<std::uint32_t>(packet.descriptor.source_rank));
    appendTreeWireU64(bytes, packet.descriptor.decomposition_epoch);
    appendTreeWireU32(bytes, packet.descriptor.derived_not_authoritative ? 1U : 0U);
    appendTreeWireU64(bytes, packet.descriptor.force_epoch);
    appendTreeWireU64(bytes, packet.descriptor.exchange_sequence);
    appendTreeWireU32(bytes, packet.geometry_frame);
    appendTreeWireDouble(bytes, packet.mass_code);
    appendTreeWireDouble(bytes, packet.center_x_comoving);
    appendTreeWireDouble(bytes, packet.center_y_comoving);
    appendTreeWireDouble(bytes, packet.center_z_comoving);
    appendTreeWireDouble(bytes, packet.min_x_comoving);
    appendTreeWireDouble(bytes, packet.max_x_comoving);
    appendTreeWireDouble(bytes, packet.min_y_comoving);
    appendTreeWireDouble(bytes, packet.max_y_comoving);
    appendTreeWireDouble(bytes, packet.min_z_comoving);
    appendTreeWireDouble(bytes, packet.max_z_comoving);
    appendTreeWireU64(bytes, packet.source_count);
    appendTreeWireU32(bytes, packet.hierarchy_level);
    appendTreeWireU32(bytes, packet.local_node_index);
    appendTreeWireU32(bytes, packet.child_count);
    appendTreeWireU32(bytes, packet.is_leaf);
  }
  if (bytes.size() != packets.size() * k_tree_pseudo_wire_bytes) {
    throw std::logic_error("tree pseudo-particle wire encoder size contract failed");
  }
  return bytes;
}

[[nodiscard, maybe_unused]] std::vector<TreePseudoParticlePacket> decodeTreePseudoPackets(
    std::span<const std::uint8_t> bytes) {
  if (bytes.size() % k_tree_pseudo_wire_bytes != 0U) {
    throw std::runtime_error("tree pseudo-particle wire payload is misaligned");
  }
  std::vector<TreePseudoParticlePacket> packets;
  packets.reserve(bytes.size() / k_tree_pseudo_wire_bytes);
  std::size_t offset = 0U;
  while (offset < bytes.size()) {
    TreePseudoParticlePacket packet;
    packet.descriptor.wire_version = readTreeWireU32(bytes, offset);
    packet.descriptor.pseudo_particle_id = readTreeWireU64(bytes, offset);
    packet.descriptor.source_rank = static_cast<int>(readTreeWireU32(bytes, offset));
    packet.descriptor.decomposition_epoch = readTreeWireU64(bytes, offset);
    packet.descriptor.derived_not_authoritative = readTreeWireU32(bytes, offset) != 0U;
    packet.descriptor.force_epoch = readTreeWireU64(bytes, offset);
    packet.descriptor.exchange_sequence = readTreeWireU64(bytes, offset);
    packet.geometry_frame = static_cast<std::uint8_t>(readTreeWireU32(bytes, offset));
    packet.mass_code = readTreeWireDouble(bytes, offset);
    packet.center_x_comoving = readTreeWireDouble(bytes, offset);
    packet.center_y_comoving = readTreeWireDouble(bytes, offset);
    packet.center_z_comoving = readTreeWireDouble(bytes, offset);
    packet.min_x_comoving = readTreeWireDouble(bytes, offset);
    packet.max_x_comoving = readTreeWireDouble(bytes, offset);
    packet.min_y_comoving = readTreeWireDouble(bytes, offset);
    packet.max_y_comoving = readTreeWireDouble(bytes, offset);
    packet.min_z_comoving = readTreeWireDouble(bytes, offset);
    packet.max_z_comoving = readTreeWireDouble(bytes, offset);
    packet.source_count = readTreeWireU64(bytes, offset);
    packet.hierarchy_level = readTreeWireU32(bytes, offset);
    packet.local_node_index = readTreeWireU32(bytes, offset);
    packet.child_count = static_cast<std::uint8_t>(readTreeWireU32(bytes, offset));
    packet.is_leaf = static_cast<std::uint8_t>(readTreeWireU32(bytes, offset));
    packets.push_back(packet);
  }
  return packets;
}

}  // namespace

std::vector<TreePseudoParticlePacket> executeBlockingTreePseudoParticleExchange(
    const MpiContext& mpi_context,
    const TreePseudoParticlePacket& local_packet) {
  std::exception_ptr local_validation_failure;
  try {
    validateTreePseudoParticlePacket(local_packet);
    if (local_packet.descriptor.source_rank != mpi_context.worldRank()) {
      throw std::invalid_argument("tree pseudo-particle packet source rank does not match MPI context");
    }
  } catch (...) {
    local_validation_failure = std::current_exception();
  }
  if (!mpi_context.isEnabled()) {
    if (local_validation_failure != nullptr) {
      std::rethrow_exception(local_validation_failure);
    }
    return {local_packet};
  }
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  std::uint64_t local_failure_vote = local_validation_failure != nullptr ? 1U : 0U;
  std::uint64_t global_failure_votes = 0U;
  MPI_Allreduce(
      &local_failure_vote, &global_failure_votes, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
  if (global_failure_votes != 0U) {
    throw std::runtime_error(
        "tree pseudo-particle exchange rejected invalid local input on one or more ranks");
  }
  std::vector<std::uint8_t> local_wire;
  std::exception_ptr local_encode_failure;
  try {
    local_wire = encodeTreePseudoPackets(
        std::span<const TreePseudoParticlePacket>(&local_packet, 1U));
    if (local_wire.size() != k_tree_pseudo_wire_bytes) {
      throw std::logic_error("tree pseudo-particle single-record wire size mismatch");
    }
  } catch (...) {
    local_encode_failure = std::current_exception();
  }
  local_failure_vote = local_encode_failure != nullptr ? 1U : 0U;
  global_failure_votes = 0U;
  MPI_Allreduce(
      &local_failure_vote, &global_failure_votes, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
  if (global_failure_votes != 0U) {
    throw std::runtime_error(
        "tree pseudo-particle exchange failed to encode local input on one or more ranks");
  }
  std::vector<std::uint8_t> gathered_wire(
      static_cast<std::size_t>(mpi_context.worldSize()) * k_tree_pseudo_wire_bytes,
      0U);
  MPI_Allgather(
      const_cast<std::uint8_t*>(local_wire.data()),
      static_cast<int>(k_tree_pseudo_wire_bytes),
      MPI_BYTE,
      gathered_wire.data(),
      static_cast<int>(k_tree_pseudo_wire_bytes),
      MPI_BYTE,
      MPI_COMM_WORLD);
  std::vector<TreePseudoParticlePacket> packets = decodeTreePseudoPackets(gathered_wire);
  for (int rank = 0; rank < mpi_context.worldSize(); ++rank) {
    const TreePseudoParticlePacket& packet = packets[static_cast<std::size_t>(rank)];
    validateTreePseudoParticlePacket(packet);
    if (packet.descriptor.source_rank != rank) {
      throw std::runtime_error("tree pseudo-particle exchange returned a packet with mismatched source rank");
    }
    if (packet.descriptor.exchange_sequence != local_packet.descriptor.exchange_sequence ||
        packet.descriptor.decomposition_epoch != local_packet.descriptor.decomposition_epoch ||
        packet.descriptor.force_epoch != local_packet.descriptor.force_epoch) {
      throw std::runtime_error("tree pseudo-particle exchange returned mixed protocol epochs");
    }
  }
  return packets;
#else
  throw std::runtime_error("tree pseudo-particle exchange requires MPI support when MPI context is enabled");
#endif
}

std::vector<TreePseudoParticlePacket> executeBlockingTreePseudoParticleHierarchyExchange(
    const MpiContext& mpi_context,
    std::span<const TreePseudoParticlePacket> local_packets,
    std::uint64_t exchange_sequence) {
  std::exception_ptr local_validation_failure;
  try {
    for (const TreePseudoParticlePacket& packet : local_packets) {
      validateTreePseudoParticlePacket(packet);
      if (packet.descriptor.source_rank != mpi_context.worldRank()) {
        throw std::invalid_argument("tree pseudo hierarchy packet source rank does not match MPI context");
      }
      if (packet.descriptor.exchange_sequence != exchange_sequence) {
        throw std::invalid_argument("tree pseudo hierarchy packet exchange sequence does not match exchange call");
      }
    }
  } catch (...) {
    local_validation_failure = std::current_exception();
  }
  if (!mpi_context.isEnabled()) {
    if (local_validation_failure != nullptr) {
      std::rethrow_exception(local_validation_failure);
    }
    return std::vector<TreePseudoParticlePacket>(local_packets.begin(), local_packets.end());
  }
#if defined(COSMOSIM_ENABLE_MPI) && COSMOSIM_ENABLE_MPI
  int communicator_world_size = 1;
  int communicator_world_rank = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &communicator_world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &communicator_world_rank);
  const auto coordinate_failure = [](std::exception_ptr local_failure, std::string_view phase) {
    const std::uint64_t local_failure_vote = local_failure != nullptr ? 1U : 0U;
    std::uint64_t global_failure_votes = 0U;
    MPI_Allreduce(
        &local_failure_vote, &global_failure_votes, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
    if (global_failure_votes == 0U) {
      return;
    }
    if (local_failure != nullptr) {
      std::rethrow_exception(local_failure);
    }
    throw std::runtime_error(
        "tree pseudo hierarchy exchange failed during " + std::string(phase) +
        " on one or more peer ranks");
  };

  std::exception_ptr communicator_validation_failure;
  try {
    if (mpi_context.worldSize() != communicator_world_size ||
        mpi_context.worldRank() != communicator_world_rank) {
      throw std::invalid_argument(
          "tree pseudo hierarchy MPI context does not match MPI_COMM_WORLD");
    }
  } catch (...) {
    communicator_validation_failure = std::current_exception();
  }
  coordinate_failure(communicator_validation_failure, "communicator validation");
  coordinate_failure(local_validation_failure, "local input validation");

  const int world_size = communicator_world_size;
  std::uint64_t min_exchange_sequence = 0U;
  std::uint64_t max_exchange_sequence = 0U;
  MPI_Allreduce(
      &exchange_sequence, &min_exchange_sequence, 1, MPI_UINT64_T, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(
      &exchange_sequence, &max_exchange_sequence, 1, MPI_UINT64_T, MPI_MAX, MPI_COMM_WORLD);
  if (min_exchange_sequence != max_exchange_sequence) {
    throw std::runtime_error("tree pseudo hierarchy ranks disagree on exchange sequence");
  }
  std::uint64_t local_count = 0U;
  std::vector<std::uint64_t> counts64;
  std::exception_ptr count_buffer_failure;
  try {
    if (local_packets.size() > static_cast<std::size_t>(std::numeric_limits<std::uint64_t>::max())) {
      throw std::overflow_error("tree pseudo hierarchy local packet count exceeds uint64_t");
    }
    local_count = static_cast<std::uint64_t>(local_packets.size());
    counts64.assign(static_cast<std::size_t>(world_size), 0U);
  } catch (...) {
    count_buffer_failure = std::current_exception();
  }
  coordinate_failure(count_buffer_failure, "count-buffer preparation");
  MPI_Allgather(
      &local_count,
      1,
      MPI_UINT64_T,
      counts64.data(),
      1,
      MPI_UINT64_T,
      MPI_COMM_WORLD);
  std::vector<std::uint8_t> local_wire;
  std::exception_ptr local_encode_failure;
  try {
    local_wire = encodeTreePseudoPackets(local_packets);
    const std::size_t expected_local_bytes = core::checkedSizeMultiply(
        core::checkedIntegralNarrow<std::size_t>(
            counts64[static_cast<std::size_t>(communicator_world_rank)],
            "tree pseudo hierarchy local packet count"),
        k_tree_pseudo_wire_bytes,
        "tree pseudo hierarchy encoded byte count");
    if (local_wire.size() != expected_local_bytes) {
      throw std::logic_error(
          "tree pseudo hierarchy encoded byte count disagrees with gathered packet count");
    }
  } catch (...) {
    local_encode_failure = std::current_exception();
  }
  coordinate_failure(local_encode_failure, "local wire encoding");

  std::vector<std::uint8_t> result_wire =
      mpi_context.allgatherBytesBounded(local_wire);
  std::vector<TreePseudoParticlePacket> result;
  std::exception_ptr received_wire_validation_failure;
  try {
    result = decodeTreePseudoPackets(result_wire);
    std::vector<std::uint64_t> per_rank_count(static_cast<std::size_t>(world_size), 0U);
    std::vector<std::uint32_t> per_rank_root_count(static_cast<std::size_t>(world_size), 0U);
    std::vector<std::uint8_t> per_rank_authoritative(static_cast<std::size_t>(world_size), 0U);
    std::vector<std::uint8_t> per_rank_derived(static_cast<std::size_t>(world_size), 0U);
    std::vector<std::unordered_set<std::uint64_t>> per_rank_ids(static_cast<std::size_t>(world_size));
    bool have_epoch_contract = false;
    std::uint64_t expected_decomposition_epoch = 0U;
    std::uint64_t expected_force_epoch = 0U;
    std::uint8_t expected_geometry_frame = 0U;
    for (const TreePseudoParticlePacket& packet : result) {
      validateTreePseudoParticlePacket(packet);
      if (packet.descriptor.source_rank < 0 || packet.descriptor.source_rank >= world_size) {
        throw std::runtime_error("tree pseudo hierarchy exchange returned packet with invalid source rank");
      }
      if (packet.descriptor.exchange_sequence != exchange_sequence) {
        throw std::runtime_error("tree pseudo hierarchy exchange returned a stale exchange sequence");
      }
      if (!have_epoch_contract) {
        expected_decomposition_epoch = packet.descriptor.decomposition_epoch;
        expected_force_epoch = packet.descriptor.force_epoch;
        expected_geometry_frame = packet.geometry_frame;
        have_epoch_contract = true;
      } else if (packet.descriptor.decomposition_epoch != expected_decomposition_epoch ||
                 packet.descriptor.force_epoch != expected_force_epoch ||
                 packet.geometry_frame != expected_geometry_frame) {
        throw std::runtime_error("tree pseudo hierarchy exchange returned mixed epochs or geometry frames");
      }
      const std::size_t source_rank = static_cast<std::size_t>(packet.descriptor.source_rank);
      if (packet.descriptor.derived_not_authoritative) {
        per_rank_derived[source_rank] = 1U;
      } else {
        per_rank_authoritative[source_rank] = 1U;
      }
      if (!per_rank_ids[source_rank].insert(packet.descriptor.pseudo_particle_id).second) {
        throw std::runtime_error("tree pseudo hierarchy exchange returned duplicate node identity");
      }
      ++per_rank_count[source_rank];
      if (packet.hierarchy_level == 0U) {
        ++per_rank_root_count[source_rank];
      }
    }
    for (int rank = 0; rank < world_size; ++rank) {
      if (per_rank_count[static_cast<std::size_t>(rank)] != counts64[static_cast<std::size_t>(rank)]) {
        throw std::runtime_error("tree pseudo hierarchy exchange source-rank coverage mismatch");
      }
      if (per_rank_authoritative[static_cast<std::size_t>(rank)] != 0U &&
          per_rank_derived[static_cast<std::size_t>(rank)] != 0U) {
        throw std::runtime_error(
            "tree top-domain exchange cannot mix authoritative and derived geometry for one rank");
      }
      if (per_rank_derived[static_cast<std::size_t>(rank)] != 0U &&
          per_rank_root_count[static_cast<std::size_t>(rank)] != 1U) {
        throw std::runtime_error(
            "tree pseudo hierarchy exchange requires exactly one root descriptor per participating rank");
      }
    }
  } catch (...) {
    received_wire_validation_failure = std::current_exception();
  }
  coordinate_failure(received_wire_validation_failure, "received wire decoding and validation");
  return result;
#else
  throw std::runtime_error("tree pseudo hierarchy exchange requires MPI support when MPI context is enabled");
#endif
}



}  // namespace cosmosim::parallel
