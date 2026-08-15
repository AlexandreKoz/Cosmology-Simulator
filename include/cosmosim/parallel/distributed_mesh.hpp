#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <string>
#include <vector>

#include "cosmosim/core/checked_arithmetic.hpp"

namespace cosmosim::core {
class ProfilerSession;
}

namespace cosmosim::parallel {

class MpiContext;
struct LoadBalanceMetrics;
struct PmSlabLayout;

struct PmSlabRange {
  std::size_t begin_x = 0;
  std::size_t end_x = 0;

  [[nodiscard]] std::size_t extentX() const noexcept {
    return (end_x >= begin_x) ? (end_x - begin_x) : 0;
  }
  [[nodiscard]] bool contains(std::size_t global_x) const noexcept {
    return global_x >= begin_x && global_x < end_x;
  }
};

struct PmMeshOwnershipDescriptor {
  std::string decomposition_mode = "slab";
  int owner_rank = 0;
  std::uint64_t decomposition_epoch = 0;
  std::size_t global_nx = 0;
  std::size_t global_ny = 0;
  std::size_t global_nz = 0;
  std::size_t begin_x = 0;
  std::size_t end_x = 0;
};

struct TreePseudoParticleDescriptor {
  std::uint32_t wire_version = 1U;
  std::uint64_t pseudo_particle_id = 0;
  int source_rank = 0;
  std::uint64_t decomposition_epoch = 0;
  std::uint64_t force_epoch = 0;
  std::uint64_t exchange_sequence = 0;
  // true: compatibility packet derived from a local gravity-tree node.
  // false: routing-only authoritative top-domain leaf supplied by the
  // decomposition layer. Both share the compact wire format because TreePM
  // uses only owner/epoch/bounds for peer discovery.
  bool derived_not_authoritative = true;
};

struct TreePseudoParticlePacket {
  TreePseudoParticleDescriptor descriptor{};
  double mass_code = 0.0;
  double center_x_comoving = 0.0;
  double center_y_comoving = 0.0;
  double center_z_comoving = 0.0;
  double min_x_comoving = 0.0;
  double max_x_comoving = 0.0;
  double min_y_comoving = 0.0;
  double max_y_comoving = 0.0;
  double min_z_comoving = 0.0;
  double max_z_comoving = 0.0;
  std::uint64_t source_count = 0;
  std::uint32_t hierarchy_level = 0;
  std::uint32_t local_node_index = 0;
  std::uint8_t child_count = 0;
  std::uint8_t is_leaf = 0;
  // 0: ordinary Euclidean frame; 1: contiguous periodic-unwrapped interval.
  std::uint8_t geometry_frame = 0;
};

struct HydroGhostCellDescriptor {
  std::uint64_t gas_cell_id = 0;
  int owner_rank = 0;
  int consumer_rank = 0;
  std::uint64_t hydro_sync_epoch = 0;
  std::uint64_t decomposition_epoch = 0;
  bool boundary_state_only = true;
};

struct HydroGhostCellRequest {
  HydroGhostCellDescriptor descriptor{};
  std::uint64_t face_key = 0;
  std::uint8_t axis = 0;
  std::uint8_t side = 0;
};

struct HydroGhostCellPayloadRecord {
  HydroGhostCellDescriptor descriptor{};
  std::uint64_t face_key = 0;
  double mass_density_comoving = 0.0;
  double momentum_density_x_comoving = 0.0;
  double momentum_density_y_comoving = 0.0;
  double momentum_density_z_comoving = 0.0;
  double total_energy_density_comoving = 0.0;
  double metal_mass_density_comoving = 0.0;
};

struct AmrPatchExchangeDescriptor {
  std::uint64_t patch_id = 0;
  int owner_rank = 0;
  int peer_rank = 0;
  std::uint64_t decomposition_epoch = 0;
  bool metadata_only = true;
};

struct AmrPatchPayloadRecord {
  std::uint64_t patch_id = 0;
  std::uint64_t parent_patch_id = 0;
  std::uint64_t morton_key = 0;
  int owner_rank = 0;
  std::uint32_t level = 0;
  std::uint32_t first_cell = 0;
  std::uint32_t cell_count = 0;
  double origin_x_comoving = 0.0;
  double origin_y_comoving = 0.0;
  double origin_z_comoving = 0.0;
  double extent_x_comoving = 0.0;
  double extent_y_comoving = 0.0;
  double extent_z_comoving = 0.0;
  std::uint16_t cell_dim_x = 0;
  std::uint16_t cell_dim_y = 0;
  std::uint16_t cell_dim_z = 0;
  std::uint64_t decomposition_epoch = 0;
  double cell_mass_sum_code = 0.0;
  double gas_internal_energy_sum_code = 0.0;
};


struct HydroConservativeFluxCorrectionRecord {
  std::uint64_t gas_cell_id = 0;
  std::uint64_t parent_particle_id = 0;
  int source_rank = 0;
  int owner_rank = 0;
  double delta_mass_density_comoving = 0.0;
  double delta_momentum_density_x_comoving = 0.0;
  double delta_momentum_density_y_comoving = 0.0;
  double delta_momentum_density_z_comoving = 0.0;
  double delta_total_energy_density_comoving = 0.0;
  double delta_metal_mass_density_comoving = 0.0;
};

struct AmrPatchCellPayloadRecord {
  std::uint64_t patch_id = 0;
  int owner_rank = 0;
  std::uint32_t local_cell_offset = 0;
  std::uint32_t patch_index = 0;
  double center_x_comoving = 0.0;
  double center_y_comoving = 0.0;
  double center_z_comoving = 0.0;
  double mass_code = 0.0;
  std::uint32_t time_bin = 0;
  std::uint64_t gas_cell_id = 0;
  std::uint64_t parent_particle_id = 0;
  double velocity_x_peculiar = 0.0;
  double velocity_y_peculiar = 0.0;
  double velocity_z_peculiar = 0.0;
  double density_code = 0.0;
  double pressure_code = 0.0;
  double internal_energy_code = 0.0;
  double temperature_code = 0.0;
  double sound_speed_code = 0.0;
};

struct DirectedAmrExchangeDiagnostics {
  std::uint64_t candidate_peer_count = 0;
  std::uint64_t neighbor_peer_count = 0;
  std::uint64_t directed_patch_descriptor_records_sent = 0;
  std::uint64_t directed_patch_descriptor_records_received = 0;
  std::uint64_t directed_patch_cell_records_sent = 0;
  std::uint64_t directed_patch_cell_records_received = 0;
  std::uint64_t directed_flux_records_sent = 0;
  std::uint64_t directed_flux_records_received = 0;
  std::uint64_t control_plane_bytes = 0;
  std::uint64_t patch_descriptor_bytes = 0;
  std::uint64_t patch_cell_payload_bytes = 0;
  std::uint64_t flux_payload_bytes = 0;
  std::uint64_t remote_patch_ghost_count = 0;
  std::uint64_t remote_interface_count = 0;
  std::uint64_t inbound_reflux_count = 0;
  std::uint64_t outbound_reflux_count = 0;
};

struct DirectedAmrPatchPayloadExchange {
  std::vector<AmrPatchPayloadRecord> patch_payloads_received;
  std::vector<AmrPatchCellPayloadRecord> patch_cell_payloads_received;
  DirectedAmrExchangeDiagnostics diagnostics{};
};

struct AmrFluxRegisterPayloadRecord {
  std::uint64_t register_key = 0;
  std::uint64_t coarse_patch_id = 0;
  std::uint64_t coarse_gas_cell_id = 0;
  std::uint64_t coarse_cell_index = 0;
  std::uint8_t level = 0;
  std::uint8_t axis = 0;
  std::uint8_t orientation = 0;
  int source_rank = 0;
  int owner_rank = 0;
  std::uint64_t gas_cell_identity_generation = 0;
  std::uint64_t patch_geometry_generation = 0;
  double coarse_mass_flux_code = 0.0;
  double coarse_momentum_x_flux_code = 0.0;
  double coarse_momentum_y_flux_code = 0.0;
  double coarse_momentum_z_flux_code = 0.0;
  double coarse_total_energy_flux_code = 0.0;
  double coarse_metal_mass_flux_code = 0.0;
  double fine_mass_flux_code = 0.0;
  double fine_momentum_x_flux_code = 0.0;
  double fine_momentum_y_flux_code = 0.0;
  double fine_momentum_z_flux_code = 0.0;
  double fine_total_energy_flux_code = 0.0;
  double fine_metal_mass_flux_code = 0.0;
  double face_area_comov = 0.0;
  double coarse_area_comov = 0.0;
  double fine_area_comov = 0.0;
  double dt_code = 0.0;
  std::uint32_t coarse_face_count = 0;
  std::uint32_t fine_face_count = 0;
};

void validatePmMeshOwnershipDescriptor(const PmMeshOwnershipDescriptor& descriptor);
void validateTreePseudoParticleDescriptor(const TreePseudoParticleDescriptor& descriptor);
void validateTreePseudoParticlePacket(const TreePseudoParticlePacket& packet);
void validateHydroGhostCellDescriptor(const HydroGhostCellDescriptor& descriptor);
void validateHydroGhostCellRequest(const HydroGhostCellRequest& request);
void validateHydroGhostCellPayloadRecord(const HydroGhostCellPayloadRecord& record);
void validateAmrPatchExchangeDescriptor(const AmrPatchExchangeDescriptor& descriptor);
void validateAmrPatchPayloadRecord(const AmrPatchPayloadRecord& record);
void validateAmrPatchCellPayloadRecord(const AmrPatchCellPayloadRecord& record);
void validateAmrFluxRegisterPayloadRecord(const AmrFluxRegisterPayloadRecord& record);
void validateHydroConservativeFluxCorrectionRecord(const HydroConservativeFluxCorrectionRecord& record);

[[nodiscard]] std::vector<TreePseudoParticlePacket> executeBlockingTreePseudoParticleExchange(
    const MpiContext& mpi_context,
    const TreePseudoParticlePacket& local_packet);

[[nodiscard]] std::vector<TreePseudoParticlePacket> executeBlockingTreePseudoParticleHierarchyExchange(
    const MpiContext& mpi_context,
    std::span<const TreePseudoParticlePacket> local_packets,
    std::uint64_t exchange_sequence = 0);

[[nodiscard]] std::vector<AmrPatchPayloadRecord> executeBlockingAmrPatchPayloadExchange(
    const MpiContext& mpi_context,
    std::span<const AmrPatchPayloadRecord> local_records,
    std::uint64_t exchange_sequence = 0);

[[nodiscard]] std::vector<AmrPatchCellPayloadRecord> executeBlockingAmrPatchCellPayloadExchange(
    const MpiContext& mpi_context,
    std::span<const AmrPatchCellPayloadRecord> local_records,
    std::uint64_t exchange_sequence = 0);

[[nodiscard]] DirectedAmrPatchPayloadExchange executeBlockingDirectedAmrPatchPayloadExchange(
    const MpiContext& mpi_context,
    std::span<const AmrPatchPayloadRecord> local_patch_records,
    std::span<const AmrPatchCellPayloadRecord> local_cell_records,
    std::uint64_t exchange_sequence = 0);

[[nodiscard]] std::vector<AmrFluxRegisterPayloadRecord> executeBlockingAmrFluxRegisterPayloadExchange(
    const MpiContext& mpi_context,
    std::span<const AmrFluxRegisterPayloadRecord> local_records,
    std::uint64_t exchange_sequence = 0);

[[nodiscard]] std::vector<HydroConservativeFluxCorrectionRecord> executeBlockingHydroConservativeFluxCorrectionExchange(
    const MpiContext& mpi_context,
    std::span<const HydroConservativeFluxCorrectionRecord> local_records,
    std::uint64_t exchange_sequence = 0);

[[nodiscard]] std::vector<HydroGhostCellRequest> executeBlockingHydroGhostCellRequestExchange(
    const MpiContext& mpi_context,
    std::span<const HydroGhostCellRequest> local_requests,
    std::uint64_t exchange_sequence = 0);

[[nodiscard]] std::vector<HydroGhostCellPayloadRecord> executeBlockingHydroGhostCellPayloadExchange(
    const MpiContext& mpi_context,
    std::span<const HydroGhostCellPayloadRecord> local_records,
    std::uint64_t exchange_sequence = 0);

struct PmSlabHaloExchangeResult {
  std::vector<double> left_halo;
  std::vector<double> right_halo;
  std::uint64_t sent_bytes = 0;
  std::uint64_t received_bytes = 0;
  std::size_t halo_depth_x = 0;
  int left_peer_rank = -1;
  int right_peer_rank = -1;
};

[[nodiscard]] PmSlabHaloExchangeResult executeBlockingPmSlabHaloExchange(
    const MpiContext& mpi_context,
    const PmSlabLayout& layout,
    std::span<const double> local_scalar_field,
    std::size_t halo_depth_x,
    bool periodic_x,
    std::uint64_t exchange_sequence = 0);

struct PmSlabLayout {
  std::size_t global_nx = 0;
  std::size_t global_ny = 0;
  std::size_t global_nz = 0;
  int world_size = 1;
  int world_rank = 0;
  PmSlabRange owned_x{};

  [[nodiscard]] std::size_t local_nx() const noexcept {
    return owned_x.extentX();
  }
  [[nodiscard]] std::size_t localCellCount() const {
    return core::checkedSizeProduct3(
        local_nx(), global_ny, global_nz, "PmSlabLayout::localCellCount");
  }
  [[nodiscard]] bool isValid() const noexcept;
  [[nodiscard]] bool ownsGlobalX(std::size_t global_x) const noexcept {
    return owned_x.contains(global_x);
  }
  [[nodiscard]] bool ownsGlobalCell(std::size_t global_x, std::size_t global_y, std::size_t global_z) const noexcept {
    if (global_y >= global_ny || global_z >= global_nz) {
      return false;
    }
    return ownsGlobalX(global_x);
  }
  [[nodiscard]] std::size_t localXFromGlobal(std::size_t global_x) const;
  [[nodiscard]] std::size_t globalXFromLocal(std::size_t local_x) const;
  [[nodiscard]] std::size_t localLinearIndex(std::size_t global_x, std::size_t global_y, std::size_t global_z) const;
  [[nodiscard]] bool ownsFullDomain() const noexcept {
    return owned_x.begin_x == 0 && owned_x.end_x == global_nx;
  }
  [[nodiscard]] PmMeshOwnershipDescriptor ownershipDescriptor(
      std::uint64_t decomposition_epoch = 0,
      std::string decomposition_mode = "slab") const;
};

[[nodiscard]] inline PmSlabRange pmOwnedXRangeForRank(std::size_t global_nx, int world_size, int rank) {
  if (global_nx == 0) {
    throw std::invalid_argument("global_nx must be positive");
  }
  if (world_size <= 0) {
    throw std::invalid_argument("world_size must be positive");
  }
  if (rank < 0 || rank >= world_size) {
    throw std::invalid_argument("rank must be within [0, world_size)");
  }
  const std::size_t world = static_cast<std::size_t>(world_size);
  const std::size_t rank_u = static_cast<std::size_t>(rank);
  // FFTW MPI uses a fixed ceil(global_nx / world_size) input block and
  // truncates the final block.  Keeping CHUI's authority map identical to
  // that backend partition avoids a hidden redistribution layer and, unlike
  // remainder-to-low-ranks partitioning, is correct when the remainder is
  // neither zero nor world_size - 1.
  const std::size_t block = global_nx / world + (global_nx % world != 0U ? 1U : 0U);
  // Clamp before multiplication: for extreme size_t domains the mathematical
  // fixed-block offset can lie beyond global_nx even though rank_u*block would
  // overflow the host index type.
  const std::size_t begin =
      rank_u > global_nx / block ? global_nx : std::min(rank_u * block, global_nx);
  const std::size_t extent = std::min(block, global_nx - begin);
  return {.begin_x = begin, .end_x = begin + extent};
}

[[nodiscard]] inline int pmOwnerRankForGlobalX(std::size_t global_nx, int world_size, std::size_t global_x) {
  if (global_x >= global_nx) {
    throw std::out_of_range("global_x must be within [0, global_nx)");
  }
  for (int rank = 0; rank < world_size; ++rank) {
    const PmSlabRange owned = pmOwnedXRangeForRank(global_nx, world_size, rank);
    if (owned.contains(global_x)) {
      return rank;
    }
  }
  throw std::logic_error("no owner rank found for global_x");
}

[[nodiscard]] inline int pmOwnerRankForGlobalCell(
    std::size_t global_nx,
    std::size_t global_ny,
    std::size_t global_nz,
    int world_size,
    std::size_t global_x,
    std::size_t global_y,
    std::size_t global_z) {
  if (global_y >= global_ny || global_z >= global_nz) {
    throw std::out_of_range("global cell y/z index out of range");
  }
  return pmOwnerRankForGlobalX(global_nx, world_size, global_x);
}
[[nodiscard]] inline PmSlabLayout makePmSlabLayout(
    std::size_t global_nx,
    std::size_t global_ny,
    std::size_t global_nz,
    int world_size,
    int world_rank) {
  PmSlabLayout layout{
      .global_nx = global_nx,
      .global_ny = global_ny,
      .global_nz = global_nz,
      .world_size = world_size,
      .world_rank = world_rank,
      .owned_x = pmOwnedXRangeForRank(global_nx, world_size, world_rank),
  };
  if (!layout.isValid()) {
    throw std::invalid_argument("constructed PM slab layout is invalid");
  }
  return layout;
}

inline bool PmSlabLayout::isValid() const noexcept {
  if (global_nx == 0 || global_ny == 0 || global_nz == 0) {
    return false;
  }
  if (world_size <= 0 || world_rank < 0 || world_rank >= world_size) {
    return false;
  }
  if (owned_x.begin_x > owned_x.end_x || owned_x.end_x > global_nx) {
    return false;
  }
  const std::size_t local_nx_value = owned_x.extentX();
  if (local_nx_value != 0U && global_ny > std::numeric_limits<std::size_t>::max() / local_nx_value) {
    return false;
  }
  const std::size_t local_xy = local_nx_value * global_ny;
  if (local_xy != 0U && global_nz > std::numeric_limits<std::size_t>::max() / local_xy) {
    return false;
  }
  const PmSlabRange expected = pmOwnedXRangeForRank(global_nx, world_size, world_rank);
  return expected.begin_x == owned_x.begin_x && expected.end_x == owned_x.end_x;
}

inline std::size_t PmSlabLayout::localXFromGlobal(std::size_t global_x) const {
  if (!ownsGlobalX(global_x)) {
    throw std::out_of_range("global x index is not owned by this PM slab");
  }
  return global_x - owned_x.begin_x;
}

inline std::size_t PmSlabLayout::globalXFromLocal(std::size_t local_x) const {
  if (local_x >= local_nx()) {
    throw std::out_of_range("local PM slab x index out of range");
  }
  return owned_x.begin_x + local_x;
}

inline std::size_t PmSlabLayout::localLinearIndex(std::size_t global_x, std::size_t global_y, std::size_t global_z) const {
  if (!ownsGlobalCell(global_x, global_y, global_z)) {
    throw std::out_of_range("global PM cell is not owned by this slab");
  }
  const std::size_t local_x = localXFromGlobal(global_x);
  return (local_x * global_ny + global_y) * global_nz + global_z;
}

void recordDistributedProfiling(
    core::ProfilerSession* profiler,
    const LoadBalanceMetrics& metrics,
    std::uint64_t ghost_exchange_send_bytes,
    std::uint64_t ghost_exchange_recv_bytes);

}  // namespace cosmosim::parallel
