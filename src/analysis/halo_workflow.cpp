#include "cosmosim/analysis/halo_workflow.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <numeric>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <utility>

#include "cosmosim/core/checked_arithmetic.hpp"
#include "cosmosim/parallel/distributed_memory.hpp"

namespace cosmosim::analysis {
namespace {

[[nodiscard]] double wrapPeriodicDelta(double delta_comov, double box_size_comov) {
  if (box_size_comov <= 0.0) {
    return delta_comov;
  }
  const double half_box = 0.5 * box_size_comov;
  if (delta_comov > half_box) {
    return delta_comov - box_size_comov;
  }
  if (delta_comov < -half_box) {
    return delta_comov + box_size_comov;
  }
  return delta_comov;
}

[[nodiscard]] double periodicDistanceSquared(
    const HaloParticleView& particles,
    std::uint32_t i,
    std::uint32_t j,
    double box_size_comov) {
  const double dx = wrapPeriodicDelta(
      particles.position_x_comoving[i] - particles.position_x_comoving[j],
      box_size_comov);
  const double dy = wrapPeriodicDelta(
      particles.position_y_comoving[i] - particles.position_y_comoving[j],
      box_size_comov);
  const double dz = wrapPeriodicDelta(
      particles.position_z_comoving[i] - particles.position_z_comoving[j],
      box_size_comov);
  return dx * dx + dy * dy + dz * dz;
}


void appendWireU32(std::vector<std::uint8_t>& out, std::uint32_t value) {
  for (unsigned shift = 0; shift < 32U; shift += 8U) {
    out.push_back(static_cast<std::uint8_t>((value >> shift) & 0xffU));
  }
}
void appendWireU64(std::vector<std::uint8_t>& out, std::uint64_t value) {
  for (unsigned shift = 0; shift < 64U; shift += 8U) {
    out.push_back(static_cast<std::uint8_t>((value >> shift) & 0xffU));
  }
}
void appendWireDouble(std::vector<std::uint8_t>& out, double value) {
  appendWireU64(out, std::bit_cast<std::uint64_t>(value));
}
[[nodiscard]] std::uint32_t readWireU32(std::span<const std::uint8_t> in, std::size_t& offset) {
  if (offset > in.size() || in.size() - offset < 4U) throw std::runtime_error("FOF MPI wire payload truncated (u32)");
  std::uint32_t value = 0U;
  for (unsigned shift = 0; shift < 32U; shift += 8U) value |= static_cast<std::uint32_t>(in[offset++]) << shift;
  return value;
}
[[nodiscard]] std::uint64_t readWireU64(std::span<const std::uint8_t> in, std::size_t& offset) {
  if (offset > in.size() || in.size() - offset < 8U) throw std::runtime_error("FOF MPI wire payload truncated (u64)");
  std::uint64_t value = 0U;
  for (unsigned shift = 0; shift < 64U; shift += 8U) value |= static_cast<std::uint64_t>(in[offset++]) << shift;
  return value;
}
[[nodiscard]] double readWireDouble(std::span<const std::uint8_t> in, std::size_t& offset) {
  return std::bit_cast<double>(readWireU64(in, offset));
}

[[nodiscard]] std::vector<std::uint8_t> packHaloParticleView(const HaloParticleView& p) {
  const std::size_t n = p.mass_code.size();
  if (p.position_x_comoving.size()!=n || p.position_y_comoving.size()!=n || p.position_z_comoving.size()!=n ||
      p.velocity_x_peculiar.size()!=n || p.velocity_y_peculiar.size()!=n || p.velocity_z_peculiar.size()!=n ||
      p.species_tag.size()!=n || p.particle_id.size()!=n) {
    throw std::invalid_argument("distributed FOF local view has mismatched extents");
  }
  std::vector<std::uint8_t> bytes;
  const std::size_t payload_bytes = core::checkedSizeMultiply(n, 68U, "distributed FOF local wire payload");
  bytes.reserve(core::checkedSizeAdd(16U, payload_bytes, "distributed FOF local wire payload"));
  appendWireU64(bytes, static_cast<std::uint64_t>(n));
  appendWireU64(bytes, p.normalized_config_hash);
  for (std::size_t i=0;i<n;++i) {
    appendWireDouble(bytes,p.position_x_comoving[i]); appendWireDouble(bytes,p.position_y_comoving[i]); appendWireDouble(bytes,p.position_z_comoving[i]);
    appendWireDouble(bytes,p.velocity_x_peculiar[i]); appendWireDouble(bytes,p.velocity_y_peculiar[i]); appendWireDouble(bytes,p.velocity_z_peculiar[i]);
    appendWireDouble(bytes,p.mass_code[i]); appendWireU32(bytes,p.species_tag[i]); appendWireU64(bytes,p.particle_id[i]);
  }
  return bytes;
}

class DisjointSet {
 public:
  explicit DisjointSet(std::size_t count) : m_parent(count), m_rank(count, 0) {
    std::iota(m_parent.begin(), m_parent.end(), 0U);
  }

  [[nodiscard]] std::uint32_t find(std::uint32_t x) {
    if (m_parent[x] != x) {
      m_parent[x] = find(m_parent[x]);
    }
    return m_parent[x];
  }

  bool unite(std::uint32_t a, std::uint32_t b) {
    std::uint32_t root_a = find(a);
    std::uint32_t root_b = find(b);
    if (root_a == root_b) {
      return false;
    }
    if (m_rank[root_a] < m_rank[root_b]) {
      std::swap(root_a, root_b);
    }
    m_parent[root_b] = root_a;
    if (m_rank[root_a] == m_rank[root_b]) {
      ++m_rank[root_a];
    }
    return true;
  }

 private:
  std::vector<std::uint32_t> m_parent;
  std::vector<std::uint8_t> m_rank;
};


struct SpatialCell {
  std::uint32_t x = 0;
  std::uint32_t y = 0;
  std::uint32_t z = 0;
};

[[nodiscard]] double wrapPeriodicCoordinate(double coordinate, double box_size) {
  double wrapped = std::fmod(coordinate, box_size);
  if (wrapped < 0.0) {
    wrapped += box_size;
  }
  if (wrapped >= box_size) {
    wrapped = 0.0;
  }
  return wrapped;
}

[[nodiscard]] std::uint64_t spatialCellKey(
    std::uint32_t x,
    std::uint32_t y,
    std::uint32_t z,
    std::uint32_t cells_per_axis) {
  const std::uint64_t n = cells_per_axis;
  return (static_cast<std::uint64_t>(x) * n + y) * n + z;
}

[[nodiscard]] std::uint32_t wrapCellCoordinate(
    std::int64_t coordinate,
    std::uint32_t cells_per_axis) noexcept {
  const std::int64_t n = static_cast<std::int64_t>(cells_per_axis);
  std::int64_t wrapped = coordinate % n;
  if (wrapped < 0) {
    wrapped += n;
  }
  return static_cast<std::uint32_t>(wrapped);
}

void linkAllPairsReference(
    const HaloParticleView& particles,
    std::span<const std::uint32_t> candidate_indices,
    double box_size_comov,
    double linking_length_sq,
    DisjointSet& ds,
    std::uint64_t& pair_checks,
    std::uint64_t& pair_links) {
  for (std::uint32_t a = 0; a < candidate_indices.size(); ++a) {
    for (std::uint32_t b = a + 1; b < candidate_indices.size(); ++b) {
      ++pair_checks;
      const double r2 = periodicDistanceSquared(
          particles, candidate_indices[a], candidate_indices[b], box_size_comov);
      if (r2 <= linking_length_sq && ds.unite(a, b)) {
        ++pair_links;
      }
    }
  }
}

std::uint64_t linkSpatialHash(
    const HaloParticleView& particles,
    std::span<const std::uint32_t> candidate_indices,
    double box_size_comov,
    double linking_length_comov,
    double linking_length_sq,
    DisjointSet& ds,
    std::uint64_t& pair_checks,
    std::uint64_t& pair_links) {
  if (!std::isfinite(linking_length_comov) || linking_length_comov <= 0.0 ||
      !std::isfinite(box_size_comov) || box_size_comov <= 0.0) {
    throw std::invalid_argument("FOF spatial hash requires finite positive box/linking lengths");
  }

  // The compact linear cell key is unique only while n^3 fits in uint64_t.
  // Capping below that representational limit merely makes cells wider (and
  // therefore candidate lists longer); it cannot miss a FOF neighbor because
  // the cell width remains >= linking_length_comov.
  constexpr std::uint32_t kMaxLinearKeyCellsPerAxis = 2'642'245U;
  const double raw_cells = std::floor(box_size_comov / linking_length_comov);
  const double bounded_cells = std::clamp(
      raw_cells, 1.0, static_cast<double>(kMaxLinearKeyCellsPerAxis));
  const auto cells_per_axis = static_cast<std::uint32_t>(bounded_cells);
  const double cell_width = box_size_comov / static_cast<double>(cells_per_axis);
  if (cell_width + 8.0 * std::numeric_limits<double>::epsilon() * box_size_comov <
      linking_length_comov) {
    throw std::logic_error("FOF spatial hash constructed cells narrower than linking length");
  }

  struct CellMember {
    std::uint64_t cell_key = 0;
    std::uint32_t local_index = 0;
  };
  std::vector<SpatialCell> cell_by_candidate(candidate_indices.size());
  std::vector<CellMember> cell_members;
  cell_members.reserve(candidate_indices.size());

  for (std::uint32_t local = 0; local < candidate_indices.size(); ++local) {
    const std::uint32_t particle = candidate_indices[local];
    const double x = particles.position_x_comoving[particle];
    const double y = particles.position_y_comoving[particle];
    const double z = particles.position_z_comoving[particle];
    if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z)) {
      throw std::invalid_argument("FOF requires finite particle positions");
    }
    const auto to_cell = [&](double coordinate) {
      const double wrapped = wrapPeriodicCoordinate(coordinate, box_size_comov);
      return std::min(
          cells_per_axis - 1U,
          static_cast<std::uint32_t>(wrapped / cell_width));
    };
    const SpatialCell cell{to_cell(x), to_cell(y), to_cell(z)};
    cell_by_candidate[local] = cell;
    cell_members.push_back(CellMember{
        .cell_key = spatialCellKey(cell.x, cell.y, cell.z, cells_per_axis),
        .local_index = local,
    });
  }
  std::sort(cell_members.begin(), cell_members.end(), [](const CellMember& lhs, const CellMember& rhs) {
    return lhs.cell_key < rhs.cell_key ||
        (lhs.cell_key == rhs.cell_key && lhs.local_index < rhs.local_index);
  });
  std::uint64_t occupied_cells = 0U;
  std::uint64_t previous_key = 0U;
  bool have_previous_key = false;
  for (const CellMember& member : cell_members) {
    if (!have_previous_key || member.cell_key != previous_key) {
      ++occupied_cells;
      previous_key = member.cell_key;
      have_previous_key = true;
    }
  }

  for (std::uint32_t a = 0; a < candidate_indices.size(); ++a) {
    const SpatialCell center = cell_by_candidate[a];
    std::array<std::uint64_t, 27> neighbor_keys{};
    std::size_t neighbor_count = 0;
    for (int dx = -1; dx <= 1; ++dx) {
      for (int dy = -1; dy <= 1; ++dy) {
        for (int dz = -1; dz <= 1; ++dz) {
          const std::uint32_t nx = wrapCellCoordinate(
              static_cast<std::int64_t>(center.x) + dx, cells_per_axis);
          const std::uint32_t ny = wrapCellCoordinate(
              static_cast<std::int64_t>(center.y) + dy, cells_per_axis);
          const std::uint32_t nz = wrapCellCoordinate(
              static_cast<std::int64_t>(center.z) + dz, cells_per_axis);
          neighbor_keys[neighbor_count++] = spatialCellKey(nx, ny, nz, cells_per_axis);
        }
      }
    }
    std::sort(neighbor_keys.begin(), neighbor_keys.begin() + static_cast<std::ptrdiff_t>(neighbor_count));
    const auto unique_end = std::unique(
        neighbor_keys.begin(), neighbor_keys.begin() + static_cast<std::ptrdiff_t>(neighbor_count));
    for (auto it = neighbor_keys.begin(); it != unique_end; ++it) {
      const std::uint64_t key = *it;
      const auto first = std::lower_bound(
          cell_members.begin(), cell_members.end(), key,
          [](const CellMember& member, std::uint64_t value) { return member.cell_key < value; });
      const auto last = std::upper_bound(
          first, cell_members.end(), key,
          [](std::uint64_t value, const CellMember& member) { return value < member.cell_key; });
      for (auto member = first; member != last; ++member) {
        const std::uint32_t b = member->local_index;
        if (b <= a) {
          continue;
        }
        ++pair_checks;
        const double r2 = periodicDistanceSquared(
            particles, candidate_indices[a], candidate_indices[b], box_size_comov);
        if (r2 <= linking_length_sq && ds.unite(a, b)) {
          ++pair_links;
        }
      }
    }
  }
  return occupied_cells;
}

}  // namespace

HaloParticleView buildHaloParticleView(const core::SimulationState& state) {
  return HaloParticleView{
      .position_x_comoving = state.particles.position_x_comoving,
      .position_y_comoving = state.particles.position_y_comoving,
      .position_z_comoving = state.particles.position_z_comoving,
      .velocity_x_peculiar = state.particles.velocity_x_peculiar,
      .velocity_y_peculiar = state.particles.velocity_y_peculiar,
      .velocity_z_peculiar = state.particles.velocity_z_peculiar,
      .mass_code = state.particles.mass_code,
      .species_tag = state.particle_sidecar.species_tag,
      .particle_id = state.particle_sidecar.particle_id,
      .normalized_config_hash = state.metadata.normalized_config_hash,
  };
}

std::uint32_t HaloCatalogSchema::schemaVersion() noexcept { return 2; }

std::string_view HaloCatalogSchema::catalogFormatName() noexcept { return "cosmosim_halo_catalog_v2"; }

std::string_view HaloCatalogSchema::mergerTreePlanFormatName() noexcept {
  return "cosmosim_merger_tree_plan_v2";
}

std::vector<std::string_view> HaloCatalogSchema::haloFields() {
  return {
      "halo_id",
      "snapshot_step_index",
      "snapshot_scale_factor",
      "particle_count",
      "total_mass_code",
      "center_of_mass_x_comov",
      "center_of_mass_y_comov",
      "center_of_mass_z_comov",
      "bulk_velocity_x_peculiar",
      "bulk_velocity_y_peculiar",
      "bulk_velocity_z_peculiar",
      "min_particle_id",
  };
}

std::vector<std::string_view> HaloCatalogSchema::subhaloFields() {
  return {
      "subhalo_id",
      "host_halo_id",
      "rank_in_host",
      "particle_count",
      "bound_mass_code",
  };
}

std::vector<std::string_view> HaloCatalogSchema::mergerTreeFields() {
  return {
      "tree_node_id",
      "halo_id",
      "snapshot_step_index",
      "snapshot_scale_factor",
      "descendant_tree_node_id",
      "has_descendant",
  };
}

FofHaloFinder::FofHaloFinder(FofConfig config) : m_config(std::move(config)) {}

bool FofHaloFinder::includeSpecies(core::ParticleSpecies species) const noexcept {
  switch (species) {
    case core::ParticleSpecies::kDarkMatter:
      return true;
    case core::ParticleSpecies::kGas:
      return m_config.include_gas;
    case core::ParticleSpecies::kStar:
      return m_config.include_stars;
    case core::ParticleSpecies::kBlackHole:
      return m_config.include_black_holes;
    case core::ParticleSpecies::kTracer:
      return false;
  }
  return false;
}

HaloCatalog FofHaloFinder::buildCatalogFromView(
    const HaloParticleView& particles,
    const core::SimulationConfig& config,
    std::uint64_t snapshot_step_index,
    double snapshot_scale_factor,
    FofProfilingCounters* profiling) const {
  HaloCatalog catalog;
  catalog.schema_version = HaloCatalogSchema::schemaVersion();
  catalog.snapshot_step_index = snapshot_step_index;
  catalog.snapshot_scale_factor = snapshot_scale_factor;
  catalog.run_name = config.output.run_name;
  catalog.normalized_config_hash = particles.normalized_config_hash;
  catalog.halo_finder = m_config.neighbor_search == FofNeighborSearch::kSpatialHash
      ? "fof_spatial_hash_v1"
      : "fof_all_pairs_reference_v1";
  catalog.halo_id_by_particle.assign(particles.mass_code.size(), k_unbound_halo_id);

  if (particles.position_x_comoving.size() != particles.mass_code.size() ||
      particles.position_y_comoving.size() != particles.mass_code.size() ||
      particles.position_z_comoving.size() != particles.mass_code.size() ||
      particles.velocity_x_peculiar.size() != particles.mass_code.size() ||
      particles.velocity_y_peculiar.size() != particles.mass_code.size() ||
      particles.velocity_z_peculiar.size() != particles.mass_code.size() ||
      particles.species_tag.size() != particles.mass_code.size() ||
      particles.particle_id.size() != particles.mass_code.size()) {
    throw std::invalid_argument("halo particle view has mismatched extents");
  }
  if (particles.mass_code.size() > std::numeric_limits<std::uint32_t>::max()) {
    throw std::overflow_error("FOF local particle view exceeds uint32 index representation");
  }

  std::vector<std::uint32_t> candidate_indices;
  candidate_indices.reserve(particles.mass_code.size());
  for (std::uint32_t i = 0; i < particles.mass_code.size(); ++i) {
    if (!core::isValidParticleSpeciesTag(particles.species_tag[i])) {
      throw std::invalid_argument("FOF encountered an invalid particle species tag");
    }
    const auto species = static_cast<core::ParticleSpecies>(particles.species_tag[i]);
    if (includeSpecies(species)) {
      candidate_indices.push_back(i);
    }
  }

  if (candidate_indices.empty()) {
    if (profiling != nullptr) {
      *profiling = FofProfilingCounters{};
      profiling->neighbor_search = m_config.neighbor_search;
    }
    return catalog;
  }

  const double box_size_comov = config.cosmology.box_size_mpc_comoving;
  if (!std::isfinite(box_size_comov) || box_size_comov <= 0.0) {
    throw std::invalid_argument("FOF requires a finite positive periodic box size");
  }
  if (!std::isfinite(m_config.linking_length_factor_mean_interparticle) ||
      m_config.linking_length_factor_mean_interparticle <= 0.0) {
    throw std::invalid_argument("FOF linking-length factor must be finite and positive");
  }
  const double mean_spacing_comov =
      box_size_comov / std::cbrt(static_cast<double>(candidate_indices.size()));
  const double linking_length_comov =
      m_config.linking_length_factor_mean_interparticle * mean_spacing_comov;
  const double linking_length_sq = linking_length_comov * linking_length_comov;

  DisjointSet ds(candidate_indices.size());
  std::uint64_t pair_checks = 0;
  std::uint64_t pair_links = 0;
  std::uint64_t occupied_spatial_cells = 0;
  if (m_config.neighbor_search == FofNeighborSearch::kSpatialHash) {
    occupied_spatial_cells = linkSpatialHash(
        particles,
        candidate_indices,
        box_size_comov,
        linking_length_comov,
        linking_length_sq,
        ds,
        pair_checks,
        pair_links);
  } else if (m_config.neighbor_search == FofNeighborSearch::kAllPairsReference) {
    linkAllPairsReference(
        particles,
        candidate_indices,
        box_size_comov,
        linking_length_sq,
        ds,
        pair_checks,
        pair_links);
  } else {
    throw std::invalid_argument("FOF received an unknown neighbor-search policy");
  }

  std::unordered_map<std::uint32_t, std::vector<std::uint32_t>> members_by_root;
  members_by_root.reserve(candidate_indices.size());
  for (std::uint32_t local = 0; local < candidate_indices.size(); ++local) {
    members_by_root[ds.find(local)].push_back(candidate_indices[local]);
  }

  std::vector<std::vector<std::uint32_t>> accepted_groups;
  accepted_groups.reserve(members_by_root.size());
  for (auto& [root, members] : members_by_root) {
    (void)root;
    if (members.size() >= m_config.min_group_size) {
      accepted_groups.push_back(std::move(members));
    }
  }

  std::sort(accepted_groups.begin(), accepted_groups.end(), [&particles](const auto& left, const auto& right) {
    const auto min_id = [&particles](const auto& group) {
      std::uint64_t value = std::numeric_limits<std::uint64_t>::max();
      for (const std::uint32_t particle : group) {
        value = std::min(value, particles.particle_id[particle]);
      }
      return value;
    };
    return min_id(left) < min_id(right);
  });

  catalog.halos.reserve(accepted_groups.size());
  for (std::uint64_t group_index = 0; group_index < accepted_groups.size(); ++group_index) {
    const auto& members = accepted_groups[group_index];
    const std::uint64_t halo_id = (snapshot_step_index << 32U) | (group_index + 1U);

    HaloCatalogEntry entry;
    entry.halo_id = halo_id;
    entry.snapshot_step_index = snapshot_step_index;
    entry.snapshot_scale_factor = snapshot_scale_factor;
    entry.particle_count = members.size();
    entry.min_particle_id = std::numeric_limits<std::uint64_t>::max();

    const std::uint32_t anchor_particle = members.front();
    const std::array<double, 3> anchor_position{
        particles.position_x_comoving[anchor_particle],
        particles.position_y_comoving[anchor_particle],
        particles.position_z_comoving[anchor_particle]};

    for (const std::uint32_t particle : members) {
      catalog.halo_id_by_particle[particle] = halo_id;
      const double mass = particles.mass_code[particle];
      if (!std::isfinite(mass) || mass <= 0.0) {
        throw std::invalid_argument("FOF requires finite positive particle masses");
      }
      entry.total_mass_code += mass;
      const std::array<double, 3> position{
          particles.position_x_comoving[particle],
          particles.position_y_comoving[particle],
          particles.position_z_comoving[particle]};
      for (std::size_t axis = 0; axis < 3; ++axis) {
        // Unwrap every member around one group anchor before averaging so a
        // halo crossing the periodic box face does not acquire a spurious
        // center near the middle of the box.
        const double unwrapped = anchor_position[axis] +
            wrapPeriodicDelta(position[axis] - anchor_position[axis], box_size_comov);
        entry.center_of_mass_comov[axis] += mass * unwrapped;
      }
      entry.bulk_velocity_peculiar[0] += mass * particles.velocity_x_peculiar[particle];
      entry.bulk_velocity_peculiar[1] += mass * particles.velocity_y_peculiar[particle];
      entry.bulk_velocity_peculiar[2] += mass * particles.velocity_z_peculiar[particle];
      entry.min_particle_id = std::min(entry.min_particle_id, particles.particle_id[particle]);
    }

    for (std::size_t axis = 0; axis < 3; ++axis) {
      entry.center_of_mass_comov[axis] = wrapPeriodicCoordinate(
          entry.center_of_mass_comov[axis] / entry.total_mass_code, box_size_comov);
      entry.bulk_velocity_peculiar[axis] /= entry.total_mass_code;
    }
    catalog.halos.push_back(entry);
  }

  // No physical subhalo finder exists yet. v2 deliberately emits an empty
  // subhalo collection instead of fabricating each host halo as its own subhalo.
  catalog.subhalo_candidates.clear();

  if (profiling != nullptr) {
    profiling->candidate_particle_count = candidate_indices.size();
    profiling->pair_checks = pair_checks;
    profiling->pair_links = pair_links;
    profiling->occupied_spatial_cells = occupied_spatial_cells;
    profiling->neighbor_search = m_config.neighbor_search;
    profiling->linking_length_comov = linking_length_comov;
  }
  return catalog;
}

HaloCatalog FofHaloFinder::buildCatalogReferenceAllPairsFromView(
    const HaloParticleView& particles,
    const core::SimulationConfig& config,
    std::uint64_t snapshot_step_index,
    double snapshot_scale_factor,
    FofProfilingCounters* profiling) const {
  FofConfig reference_config = m_config;
  reference_config.neighbor_search = FofNeighborSearch::kAllPairsReference;
  return FofHaloFinder(reference_config).buildCatalogFromView(
      particles, config, snapshot_step_index, snapshot_scale_factor, profiling);
}

HaloCatalog FofHaloFinder::buildCatalog(
    const core::SimulationState& state,
    const core::SimulationConfig& config,
    std::uint64_t snapshot_step_index,
    double snapshot_scale_factor,
    FofProfilingCounters* profiling) const {
  return buildCatalogFromView(
      buildHaloParticleView(state),
      config,
      snapshot_step_index,
      snapshot_scale_factor,
      profiling);
}


DistributedHaloCatalogResult FofHaloFinder::buildDistributedCatalogFromView(
    const HaloParticleView& local_particles,
    const core::SimulationConfig& config,
    const parallel::MpiContext& mpi_context,
    std::uint64_t snapshot_step_index,
    double snapshot_scale_factor,
    int root_rank,
    FofProfilingCounters* root_profiling) const {
  // Collective phases must fail coherently. A rank-local validation/allocation
  // failure is reduced before the first gather so no peer can enter a
  // collective that another rank skipped. Root-side assembly errors are
  // broadcast before the membership-label broadcast for the same reason.
  std::vector<std::uint8_t> local_wire;
  std::string local_pack_error;
  try {
    local_wire = packHaloParticleView(local_particles);
  } catch (const std::exception& error) {
    local_pack_error = error.what();
  } catch (...) {
    local_pack_error = "unknown local distributed FOF packing failure";
  }

  const std::uint64_t failed_pack_ranks = mpi_context.allreduceSumUint64(local_pack_error.empty() ? 0U : 1U);
  if (failed_pack_ranks != 0U) {
    if (!local_pack_error.empty()) {
      throw std::runtime_error("distributed FOF local preflight failed: " + local_pack_error);
    }
    throw std::runtime_error("distributed FOF local preflight failed on another MPI rank");
  }

  const std::vector<std::uint8_t> gathered = mpi_context.gatherBytesToRoot(local_wire, root_rank);

  DistributedHaloCatalogResult result;
  result.root_rank = root_rank;
  result.has_global_catalog = mpi_context.worldRank() == root_rank;
  std::vector<std::uint8_t> label_wire;
  std::vector<std::uint8_t> status_wire;

  if (result.has_global_catalog) {
    try {
      std::vector<double> px, py, pz, vx, vy, vz, mass;
      std::vector<std::uint32_t> species;
      std::vector<std::uint64_t> ids;
      std::size_t offset = 0;
      std::uint64_t expected_hash = 0U;
      for (int rank = 0; rank < mpi_context.worldSize(); ++rank) {
        const std::uint64_t n64 = readWireU64(gathered, offset);
        const std::uint64_t config_hash = readWireU64(gathered, offset);
        if (rank == 0) {
          expected_hash = config_hash;
        }
        if (config_hash != expected_hash) {
          throw std::runtime_error("distributed FOF config hash differs across ranks");
        }
        const std::size_t n = core::checkedIntegralNarrow<std::size_t>(
            n64, "distributed FOF particle count");
        const auto reserve_next = [n](auto& values, std::string_view context) {
          values.reserve(core::checkedSizeAdd(values.size(), n, context));
        };
        reserve_next(px, "distributed FOF x reserve");
        reserve_next(py, "distributed FOF y reserve");
        reserve_next(pz, "distributed FOF z reserve");
        reserve_next(vx, "distributed FOF vx reserve");
        reserve_next(vy, "distributed FOF vy reserve");
        reserve_next(vz, "distributed FOF vz reserve");
        reserve_next(mass, "distributed FOF mass reserve");
        reserve_next(species, "distributed FOF species reserve");
        reserve_next(ids, "distributed FOF ID reserve");
        for (std::size_t i = 0; i < n; ++i) {
          px.push_back(readWireDouble(gathered, offset));
          py.push_back(readWireDouble(gathered, offset));
          pz.push_back(readWireDouble(gathered, offset));
          vx.push_back(readWireDouble(gathered, offset));
          vy.push_back(readWireDouble(gathered, offset));
          vz.push_back(readWireDouble(gathered, offset));
          mass.push_back(readWireDouble(gathered, offset));
          species.push_back(readWireU32(gathered, offset));
          ids.push_back(readWireU64(gathered, offset));
        }
      }
      if (offset != gathered.size()) {
        throw std::runtime_error("distributed FOF gathered wire payload has trailing bytes");
      }

      HaloParticleView global{px, py, pz, vx, vy, vz, mass, species, ids, expected_hash};
      result.root_catalog = buildCatalogFromView(
          global, config, snapshot_step_index, snapshot_scale_factor, root_profiling);
      result.root_catalog.halo_finder = "fof_spatial_hash_mpi_root_merge_v1";

      const std::size_t label_payload = core::checkedSizeMultiply(
          ids.size(), 16U, "distributed FOF label payload");
      label_wire.reserve(core::checkedSizeAdd(8U, label_payload, "distributed FOF label payload"));
      appendWireU64(label_wire, static_cast<std::uint64_t>(ids.size()));
      for (std::size_t i = 0; i < ids.size(); ++i) {
        appendWireU64(label_wire, ids[i]);
        appendWireU64(label_wire, result.root_catalog.halo_id_by_particle[i]);
      }
      status_wire.push_back(1U);
    } catch (const std::exception& error) {
      status_wire.push_back(0U);
      const std::string message = error.what();
      status_wire.insert(status_wire.end(), message.begin(), message.end());
    } catch (...) {
      static constexpr std::string_view message = "unknown root-side distributed FOF assembly failure";
      status_wire.push_back(0U);
      status_wire.insert(status_wire.end(), message.begin(), message.end());
    }
  }

  const std::vector<std::uint8_t> status =
      mpi_context.broadcastBytesFromRoot(status_wire, root_rank);
  if (status.empty() || status.front() == 0U) {
    const std::string message = status.size() > 1U
        ? std::string(status.begin() + 1, status.end())
        : std::string("distributed FOF root assembly failed without diagnostic");
    throw std::runtime_error("distributed FOF collective assembly failed: " + message);
  }

  const std::vector<std::uint8_t> labels =
      mpi_context.broadcastBytesFromRoot(label_wire, root_rank);
  std::size_t offset = 0;
  const std::uint64_t label_count = readWireU64(labels, offset);
  const std::size_t label_count_size = core::checkedIntegralNarrow<std::size_t>(
      label_count, "distributed FOF label count");
  std::unordered_map<std::uint64_t, std::uint64_t> halo_by_id;
  halo_by_id.reserve(label_count_size);
  for (std::size_t i = 0; i < label_count_size; ++i) {
    const auto [it, inserted] = halo_by_id.emplace(
        readWireU64(labels, offset), readWireU64(labels, offset));
    if (!inserted) {
      throw std::runtime_error("distributed FOF encountered duplicate particle ID in global membership map");
    }
  }
  if (offset != labels.size()) {
    throw std::runtime_error("distributed FOF label payload has trailing bytes");
  }

  result.local_halo_id_by_particle.assign(local_particles.particle_id.size(), k_unbound_halo_id);
  for (std::size_t i = 0; i < local_particles.particle_id.size(); ++i) {
    const auto it = halo_by_id.find(local_particles.particle_id[i]);
    if (it == halo_by_id.end()) {
      throw std::runtime_error("distributed FOF global membership map omitted a local particle ID");
    }
    result.local_halo_id_by_particle[i] = it->second;
  }
  return result;
}

MergerTreePlan MergerTreePlanner::buildPlan(const HaloCatalog& current, const HaloCatalog* previous) const {
  MergerTreePlan tree_plan;
  tree_plan.schema_version = HaloCatalogSchema::schemaVersion();
  tree_plan.run_name = current.run_name;
  tree_plan.normalized_config_hash = current.normalized_config_hash;
  tree_plan.nodes.reserve(current.halos.size());

  std::unordered_map<std::uint64_t, std::uint64_t> prev_halo_to_node_id;
  if (previous != nullptr) {
    for (std::size_t i = 0; i < previous->halos.size(); ++i) {
      prev_halo_to_node_id.emplace(previous->halos[i].halo_id, (previous->halos[i].snapshot_step_index << 32U) | (i + 1U));
    }
  }

  for (std::size_t i = 0; i < current.halos.size(); ++i) {
    const auto& halo = current.halos[i];
    MergerTreeNodePlan node;
    node.tree_node_id = (halo.snapshot_step_index << 32U) | (i + 1U);
    node.halo_id = halo.halo_id;
    node.snapshot_step_index = halo.snapshot_step_index;
    node.snapshot_scale_factor = halo.snapshot_scale_factor;

    if (previous != nullptr) {
      const auto it = prev_halo_to_node_id.find(halo.halo_id - (1ULL << 32U));
      if (it != prev_halo_to_node_id.end()) {
        node.descendant_tree_node_id = it->second;
        node.has_descendant = true;
      }
    }

    tree_plan.nodes.push_back(node);
  }

  return tree_plan;
}

HaloWorkflowPlanner::HaloWorkflowPlanner(core::SimulationConfig config)
    : m_config(std::move(config)), m_fof(fofConfigFromSimulationConfig(m_config)) {}

std::string_view HaloWorkflowPlanner::ownershipBoundary() {
  return "On-the-fly planner owns only schema-governed analysis products; simulation state remains authoritative in core.";
}

HaloWorkflowMode HaloWorkflowPlanner::workflowMode() const {
  return m_config.analysis.halo_on_the_fly ? HaloWorkflowMode::kOnTheFly : HaloWorkflowMode::kPostProcess;
}

HaloWorkflowReport HaloWorkflowPlanner::runSnapshotWorkflow(
    const core::SimulationState& state,
    std::uint64_t snapshot_step_index,
    double snapshot_scale_factor,
    const HaloCatalog* previous_catalog) const {
  if (m_config.parallel.mpi_ranks_expected > 1) {
    throw std::invalid_argument(
        "HaloWorkflowPlanner serial snapshot path cannot run with mpi_ranks_expected > 1; "
        "use runSnapshotWorkflowDistributed so rank-boundary halos are merged");
  }
  HaloWorkflowReport report;
  report.catalog = m_fof.buildCatalog(
      state,
      m_config,
      snapshot_step_index,
      snapshot_scale_factor,
      &report.profiling);
  report.tree_plan = m_tree_planner.buildPlan(report.catalog, previous_catalog);
  report.halo_catalog_path = haloCatalogPath(snapshot_step_index);
  report.merger_tree_plan_path = mergerTreePlanPath(snapshot_step_index);

  writeHaloCatalog(report.catalog, report.halo_catalog_path);
  writeMergerTreePlan(report.tree_plan, report.merger_tree_plan_path);
  report.local_halo_id_by_particle = report.catalog.halo_id_by_particle;
  return report;
}

HaloWorkflowReport HaloWorkflowPlanner::runSnapshotWorkflowDistributed(
    const core::SimulationState& local_state,
    const parallel::MpiContext& mpi_context,
    std::uint64_t snapshot_step_index,
    double snapshot_scale_factor,
    const HaloCatalog* previous_root_catalog,
    int root_rank) const {
  if (!mpi_context.isEnabled() || mpi_context.worldSize() <= 1) {
    return runSnapshotWorkflow(local_state, snapshot_step_index, snapshot_scale_factor, previous_root_catalog);
  }
  if (m_config.parallel.mpi_ranks_expected > 0) {
    mpi_context.validateExpectedWorldSizeOrThrow(m_config.parallel.mpi_ranks_expected);
  }
  HaloWorkflowReport report;
  const DistributedHaloCatalogResult distributed = m_fof.buildDistributedCatalogFromView(
      buildHaloParticleView(local_state), m_config, mpi_context,
      snapshot_step_index, snapshot_scale_factor, root_rank, &report.profiling);
  report.local_halo_id_by_particle = distributed.local_halo_id_by_particle;
  report.has_global_catalog = distributed.has_global_catalog;
  if (distributed.has_global_catalog) {
    report.catalog = distributed.root_catalog;
    report.tree_plan = m_tree_planner.buildPlan(report.catalog, previous_root_catalog);
    report.halo_catalog_path = haloCatalogPath(snapshot_step_index);
    report.merger_tree_plan_path = mergerTreePlanPath(snapshot_step_index);
    writeHaloCatalog(report.catalog, report.halo_catalog_path);
    writeMergerTreePlan(report.tree_plan, report.merger_tree_plan_path);
  }
  return report;
}

void HaloWorkflowPlanner::writeHaloCatalog(const HaloCatalog& catalog, const std::filesystem::path& path) const {
  std::filesystem::create_directories(path.parent_path());
  std::ofstream out(path);
  if (!out) {
    throw std::runtime_error("failed to open halo catalog path: " + path.string());
  }

  out << "{\n";
  out << "  \"format\": \"" << HaloCatalogSchema::catalogFormatName() << "\",\n";
  out << "  \"schema_version\": " << catalog.schema_version << ",\n";
  out << "  \"run_name\": \"" << catalog.run_name << "\",\n";
  out << "  \"normalized_config_hash\": " << catalog.normalized_config_hash << ",\n";
  out << "  \"snapshot_step_index\": " << catalog.snapshot_step_index << ",\n";
  out << "  \"snapshot_scale_factor\": " << catalog.snapshot_scale_factor << ",\n";
  out << "  \"halo_finder\": \"fof\",\n";
  out << "  \"halos\": [\n";
  for (std::size_t i = 0; i < catalog.halos.size(); ++i) {
    const auto& h = catalog.halos[i];
    out << "    {\"halo_id\":" << h.halo_id << ",\"particle_count\":" << h.particle_count
        << ",\"total_mass_code\":" << h.total_mass_code << ",\"min_particle_id\":" << h.min_particle_id
        << ",\"center_of_mass_comov\":[" << h.center_of_mass_comov[0] << ',' << h.center_of_mass_comov[1] << ','
        << h.center_of_mass_comov[2] << "]"
        << ",\"bulk_velocity_peculiar\":[" << h.bulk_velocity_peculiar[0] << ',' << h.bulk_velocity_peculiar[1]
        << ',' << h.bulk_velocity_peculiar[2] << "]}";
    out << ((i + 1U < catalog.halos.size()) ? ",\n" : "\n");
  }
  out << "  ],\n";

  out << "  \"subhalo_candidates\": [\n";
  for (std::size_t i = 0; i < catalog.subhalo_candidates.size(); ++i) {
    const auto& s = catalog.subhalo_candidates[i];
    out << "    {\"subhalo_id\":" << s.subhalo_id << ",\"host_halo_id\":" << s.host_halo_id
        << ",\"rank_in_host\":" << s.rank_in_host << ",\"particle_count\":" << s.particle_count
        << ",\"bound_mass_code\":" << s.bound_mass_code << "}";
    out << ((i + 1U < catalog.subhalo_candidates.size()) ? ",\n" : "\n");
  }
  out << "  ]\n";
  out << "}\n";
}

void HaloWorkflowPlanner::writeMergerTreePlan(const MergerTreePlan& tree_plan, const std::filesystem::path& path) const {
  std::filesystem::create_directories(path.parent_path());
  std::ofstream out(path);
  if (!out) {
    throw std::runtime_error("failed to open merger-tree plan path: " + path.string());
  }

  out << "{\n";
  out << "  \"format\": \"" << HaloCatalogSchema::mergerTreePlanFormatName() << "\",\n";
  out << "  \"schema_version\": " << tree_plan.schema_version << ",\n";
  out << "  \"run_name\": \"" << tree_plan.run_name << "\",\n";
  out << "  \"normalized_config_hash\": " << tree_plan.normalized_config_hash << ",\n";
  out << "  \"nodes\": [\n";
  for (std::size_t i = 0; i < tree_plan.nodes.size(); ++i) {
    const auto& node = tree_plan.nodes[i];
    out << "    {\"tree_node_id\":" << node.tree_node_id << ",\"halo_id\":" << node.halo_id
        << ",\"snapshot_step_index\":" << node.snapshot_step_index
        << ",\"snapshot_scale_factor\":" << node.snapshot_scale_factor
        << ",\"descendant_tree_node_id\":" << node.descendant_tree_node_id
        << ",\"has_descendant\":" << (node.has_descendant ? "true" : "false") << "}";
    out << ((i + 1U < tree_plan.nodes.size()) ? ",\n" : "\n");
  }
  out << "  ]\n";
  out << "}\n";
}

std::filesystem::path HaloWorkflowPlanner::analysisDirectory() const {
  return std::filesystem::path(m_config.output.output_directory) / m_config.output.run_name / "analysis";
}

std::filesystem::path HaloWorkflowPlanner::haloCatalogPath(std::uint64_t snapshot_step_index) const {
  std::ostringstream name;
  name << m_config.analysis.halo_catalog_stem << "_step_" << snapshot_step_index << ".json";
  return analysisDirectory() / name.str();
}

std::filesystem::path HaloWorkflowPlanner::mergerTreePlanPath(std::uint64_t snapshot_step_index) const {
  std::ostringstream name;
  name << m_config.analysis.merger_tree_stem << "_step_" << snapshot_step_index << ".json";
  return analysisDirectory() / name.str();
}

FofConfig HaloWorkflowPlanner::fofConfigFromSimulationConfig(const core::SimulationConfig& config) {
  FofConfig fof;
  fof.linking_length_factor_mean_interparticle = config.analysis.halo_fof_linking_length_factor;
  fof.min_group_size = static_cast<std::uint64_t>(config.analysis.halo_fof_min_group_size);
  fof.include_gas = config.analysis.halo_include_gas;
  fof.include_stars = config.analysis.halo_include_stars;
  fof.include_black_holes = config.analysis.halo_include_black_holes;
  return fof;
}

}  // namespace cosmosim::analysis
