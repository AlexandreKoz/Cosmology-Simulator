#include "cosmosim/core/simulation_state.hpp"

#include <charconv>
#include <cstdlib>
#include <limits>
#include <new>
#include <numeric>
#include <sstream>
#include <tuple>
#include <unordered_set>

namespace cosmosim::core {
namespace {

constexpr std::size_t k_species_count = 5;

// Species tag validity helper for sidecar-to-ledger consistency checks.
[[nodiscard]] bool isValidSpeciesTag(std::uint32_t value) {
  return value <= static_cast<std::uint32_t>(ParticleSpecies::kTracer);
}

[[nodiscard]] std::string trim(std::string value) {
  const auto first = value.find_first_not_of(" \t\n\r");
  if (first == std::string::npos) {
    return {};
  }
  const auto last = value.find_last_not_of(" \t\n\r");
  return value.substr(first, last - first + 1);
}

[[nodiscard]] std::uint64_t parseUint64(std::string_view value, const std::string& key) {
  std::uint64_t parsed = 0;
  const auto* begin = value.data();
  const auto* end = value.data() + value.size();
  const auto [ptr, ec] = std::from_chars(begin, end, parsed);
  if (ec != std::errc{} || ptr != end) {
    throw std::invalid_argument("StateMetadata.deserialize: invalid integer for key '" + key + "'");
  }
  return parsed;
}

[[nodiscard]] std::uint32_t parseUint32(std::string_view value, const std::string& key) {
  const auto parsed = parseUint64(value, key);
  if (parsed > std::numeric_limits<std::uint32_t>::max()) {
    throw std::invalid_argument("StateMetadata.deserialize: integer overflow for key '" + key + "'");
  }
  return static_cast<std::uint32_t>(parsed);
}

[[nodiscard]] double parseDouble(std::string_view value, const std::string& key) {
  std::string temp(value);
  char* end = nullptr;
  const double parsed = std::strtod(temp.c_str(), &end);
  if (end != temp.c_str() + static_cast<std::ptrdiff_t>(temp.size())) {
    throw std::invalid_argument(
        "StateMetadata.deserialize: invalid floating-point value for key '" + key + "'");
  }
  return parsed;
}

}  // namespace

// Resize all particle skeleton lanes together so every index always addresses a
// full gravity-hot tuple.
void ParticleSoa::resize(std::size_t count) {
  // Keep all hot arrays in lock-step to preserve contiguous index ownership.
  position_x_comoving.resize(count);
  position_y_comoving.resize(count);
  position_z_comoving.resize(count);
  velocity_x_peculiar.resize(count);
  velocity_y_peculiar.resize(count);
  velocity_z_peculiar.resize(count);
  mass_code.resize(count);
  time_bin.resize(count);
}

// Report logical particle row count shared by all skeleton lanes.
std::size_t ParticleSoa::size() const noexcept { return position_x_comoving.size(); }

// Validate lock-step sizing so kernels can safely use a single index space.
bool ParticleSoa::isConsistent() const noexcept {
  const std::size_t expected = position_x_comoving.size();
  return position_y_comoving.size() == expected && position_z_comoving.size() == expected &&
         velocity_x_peculiar.size() == expected && velocity_y_peculiar.size() == expected &&
         velocity_z_peculiar.size() == expected && mass_code.size() == expected &&
         time_bin.size() == expected;
}

// Resize sidecar metadata lanes in lock-step with particle skeleton rows.
void ParticleSidecar::resize(std::size_t count) {
  // Sidecar arrays share the same particle index space as ParticleSoa.
  particle_id.resize(count);
  sfc_key.resize(count);
  species_tag.resize(count);
  particle_flags.resize(count);
  owning_rank.resize(count);
}

// Report logical metadata row count.
std::size_t ParticleSidecar::size() const noexcept { return particle_id.size(); }

// Validate sidecar lane consistency before ownership invariants are checked.
bool ParticleSidecar::isConsistent() const noexcept {
  const std::size_t expected = particle_id.size();
  return sfc_key.size() == expected && species_tag.size() == expected && particle_flags.size() == expected &&
         owning_rank.size() == expected;
}

// Resize gravity-facing cell skeleton lanes while keeping one shared index
// space with gas thermodynamic sidecars.
void CellSoa::resize(std::size_t count) {
  // Cell gravity skeleton remains lock-step for contiguous gather/scatter patterns.
  center_x_comoving.resize(count);
  center_y_comoving.resize(count);
  center_z_comoving.resize(count);
  mass_code.resize(count);
  time_bin.resize(count);
  patch_index.resize(count);
}

// Report logical cell row count.
std::size_t CellSoa::size() const noexcept { return center_x_comoving.size(); }

// Validate gravity-cell skeleton lock-step sizes.
bool CellSoa::isConsistent() const noexcept {
  const std::size_t expected = center_x_comoving.size();
  return center_y_comoving.size() == expected && center_z_comoving.size() == expected &&
         mass_code.size() == expected && time_bin.size() == expected && patch_index.size() == expected;
}

// Resize gas-only thermodynamic and reconstruction lanes together.
void GasCellSidecar::resize(std::size_t count) {
  // Gas-only thermodynamic fields stay separated from gravity-hot fields.
  density_code.resize(count);
  pressure_code.resize(count);
  internal_energy_code.resize(count);
  temperature_code.resize(count);
  sound_speed_code.resize(count);
  recon_gradient_x.resize(count);
  recon_gradient_y.resize(count);
  recon_gradient_z.resize(count);
}

// Report logical gas sidecar row count.
std::size_t GasCellSidecar::size() const noexcept { return density_code.size(); }

// Validate gas-sidecar lane consistency.
bool GasCellSidecar::isConsistent() const noexcept {
  const std::size_t expected = density_code.size();
  return pressure_code.size() == expected && internal_energy_code.size() == expected &&
         temperature_code.size() == expected && sound_speed_code.size() == expected &&
         recon_gradient_x.size() == expected && recon_gradient_y.size() == expected &&
         recon_gradient_z.size() == expected;
}

// Resize star metadata lanes indexed by star-local rows.
void StarParticleSidecar::resize(std::size_t count) {
  // Star sidecar rows map 1:1 to star species-local particle indices.
  particle_index.resize(count);
  formation_scale_factor.resize(count);
  birth_mass_code.resize(count);
  metallicity_mass_fraction.resize(count);
  stellar_age_years_last.resize(count);
  stellar_returned_mass_cumulative_code.resize(count);
  stellar_returned_metals_cumulative_code.resize(count);
  stellar_feedback_energy_cumulative_erg.resize(count);
  for (std::size_t channel = 0; channel < stellar_returned_mass_channel_cumulative_code.size(); ++channel) {
    stellar_returned_mass_channel_cumulative_code[channel].resize(count);
    stellar_returned_metals_channel_cumulative_code[channel].resize(count);
    stellar_feedback_energy_channel_cumulative_erg[channel].resize(count);
  }
}

// Report star sidecar row count.
std::size_t StarParticleSidecar::size() const noexcept { return particle_index.size(); }

// Validate star metadata lane consistency.
bool StarParticleSidecar::isConsistent() const noexcept {
  const std::size_t expected = particle_index.size();
  if (formation_scale_factor.size() != expected || birth_mass_code.size() != expected ||
      metallicity_mass_fraction.size() != expected || stellar_age_years_last.size() != expected ||
      stellar_returned_mass_cumulative_code.size() != expected ||
      stellar_returned_metals_cumulative_code.size() != expected ||
      stellar_feedback_energy_cumulative_erg.size() != expected) {
    return false;
  }

  for (std::size_t channel = 0; channel < stellar_returned_mass_channel_cumulative_code.size(); ++channel) {
    if (stellar_returned_mass_channel_cumulative_code[channel].size() != expected ||
        stellar_returned_metals_channel_cumulative_code[channel].size() != expected ||
        stellar_feedback_energy_channel_cumulative_erg[channel].size() != expected) {
      return false;
    }
  }
  return true;
}

// Resize black-hole metadata lanes indexed by BH-local rows.
void BlackHoleParticleSidecar::resize(std::size_t count) {
  // Black-hole sidecar rows map 1:1 to BH species-local particle indices.
  particle_index.resize(count);
  host_cell_index.resize(count);
  subgrid_mass_code.resize(count);
  accretion_rate_code.resize(count);
  feedback_energy_code.resize(count);
  eddington_ratio.resize(count);
  cumulative_accreted_mass_code.resize(count);
  cumulative_feedback_energy_code.resize(count);
  duty_cycle_active_time_code.resize(count);
  duty_cycle_total_time_code.resize(count);
}

// Report black-hole sidecar row count.
std::size_t BlackHoleParticleSidecar::size() const noexcept { return particle_index.size(); }

// Validate black-hole metadata lane consistency.
bool BlackHoleParticleSidecar::isConsistent() const noexcept {
  const std::size_t expected = particle_index.size();
  return host_cell_index.size() == expected && subgrid_mass_code.size() == expected &&
         accretion_rate_code.size() == expected && feedback_energy_code.size() == expected &&
         eddington_ratio.size() == expected && cumulative_accreted_mass_code.size() == expected &&
         cumulative_feedback_energy_code.size() == expected &&
         duty_cycle_active_time_code.size() == expected && duty_cycle_total_time_code.size() == expected;
}

// Resize tracer metadata lanes indexed by tracer-local rows.
void TracerParticleSidecar::resize(std::size_t count) {
  // Tracer sidecar rows map 1:1 to tracer species-local particle indices.
  particle_index.resize(count);
  parent_particle_id.resize(count);
  injection_step.resize(count);
  host_cell_index.resize(count);
  mass_fraction_of_host.resize(count);
  last_host_mass_code.resize(count);
  cumulative_exchanged_mass_code.resize(count);
}

// Report tracer sidecar row count.
std::size_t TracerParticleSidecar::size() const noexcept { return particle_index.size(); }

// Validate tracer metadata lane consistency.
bool TracerParticleSidecar::isConsistent() const noexcept {
  const std::size_t expected = particle_index.size();
  return parent_particle_id.size() == expected && injection_step.size() == expected &&
         host_cell_index.size() == expected && mass_fraction_of_host.size() == expected &&
         last_host_mass_code.size() == expected &&
         cumulative_exchanged_mass_code.size() == expected;
}

// Resize AMR patch descriptor lanes while preserving contiguous range contract.
void PatchSoa::resize(std::size_t count) {
  // Patch descriptors are stored in compact SoA form for traversal locality.
  patch_id.resize(count);
  level.resize(count);
  first_cell.resize(count);
  cell_count.resize(count);
}

// Report patch descriptor row count.
std::size_t PatchSoa::size() const noexcept { return patch_id.size(); }

// Validate patch descriptor lane consistency.
bool PatchSoa::isConsistent() const noexcept {
  const std::size_t expected = patch_id.size();
  return level.size() == expected && first_cell.size() == expected && cell_count.size() == expected;
}

// Sum auditable per-species ledger counts.
std::uint64_t SpeciesContainer::totalCount() const noexcept {
  std::uint64_t total = 0;
  for (const auto count : count_by_species) {
    total += count;
  }
  return total;
}

bool SpeciesContainer::isConsistentWith(const ParticleSidecar& sidecar) const noexcept {
  // Recompute measured species counts from sidecar tags and compare to ledger.
  std::array<std::uint64_t, k_species_count> measured{};
  for (const auto tag : sidecar.species_tag) {
    if (!isValidSpeciesTag(tag)) {
      return false;
    }
    ++measured.at(tag);
  }
  return measured == count_by_species;
}

// Rebuild explicit species-local/global mapping from sidecar species tags.
void ParticleSpeciesIndex::rebuild(const ParticleSidecar& sidecar) {
  for (auto& indices : global_index_by_species) {
    indices.clear();
  }

  local_index_by_global.resize(sidecar.size());
  for (std::size_t global_index = 0; global_index < sidecar.size(); ++global_index) {
    const auto tag = sidecar.species_tag[global_index];
    if (!isValidSpeciesTag(tag)) {
      throw std::invalid_argument("ParticleSpeciesIndex.rebuild: invalid species tag");
    }
    auto& species_indices = global_index_by_species[tag];
    local_index_by_global[global_index] = static_cast<std::uint32_t>(species_indices.size());
    species_indices.push_back(static_cast<std::uint32_t>(global_index));
  }
}

// Return number of particles for one species in the current mapping.
std::size_t ParticleSpeciesIndex::count(ParticleSpecies species) const noexcept {
  return global_index_by_species[static_cast<std::uint32_t>(species)].size();
}

// Return immutable global indices for one species, suitable for branch-light
// species-specific loops.
std::span<const std::uint32_t> ParticleSpeciesIndex::globalIndices(ParticleSpecies species) const noexcept {
  return global_index_by_species[static_cast<std::uint32_t>(species)];
}

// Translate global particle index -> species-local index.
std::uint32_t ParticleSpeciesIndex::localIndex(std::uint32_t global_index) const {
  if (global_index >= local_index_by_global.size()) {
    throw std::out_of_range("ParticleSpeciesIndex.localIndex: global index out of range");
  }
  return local_index_by_global[global_index];
}

// Translate species-local index -> global particle index.
std::uint32_t ParticleSpeciesIndex::globalIndex(ParticleSpecies species, std::uint32_t local_index) const {
  const auto& species_indices = global_index_by_species[static_cast<std::uint32_t>(species)];
  if (local_index >= species_indices.size()) {
    throw std::out_of_range("ParticleSpeciesIndex.globalIndex: local index out of range");
  }
  return species_indices[local_index];
}

std::string StateMetadata::serialize() const {
  // Stable key order intentionally supports deterministic snapshots/restarts.
  std::ostringstream out;
  out << "schema_version=" << schema_version << '\n';
  out << "run_name=" << run_name << '\n';
  out << "normalized_config_hash=" << normalized_config_hash << '\n';
  out << "normalized_config_hash_hex=" << normalized_config_hash_hex << '\n';
  out << "step_index=" << step_index << '\n';
  out << "scale_factor=" << scale_factor << '\n';
  out << "snapshot_stem=" << snapshot_stem << '\n';
  out << "restart_stem=" << restart_stem << '\n';
  return out.str();
}

StateMetadata StateMetadata::deserialize(std::string_view text) {
  StateMetadata metadata;
  std::istringstream in{std::string(text)};
  std::string line;

  while (std::getline(in, line)) {
    const std::string trimmed = trim(line);
    if (trimmed.empty()) {
      continue;
    }

    const auto equal_pos = trimmed.find('=');
    if (equal_pos == std::string::npos || equal_pos == 0 || equal_pos == trimmed.size() - 1) {
      throw std::invalid_argument("StateMetadata.deserialize: malformed line '" + line + "'");
    }

    const std::string key = trim(trimmed.substr(0, equal_pos));
    const std::string value = trim(trimmed.substr(equal_pos + 1));

    if (key == "schema_version") {
      metadata.schema_version = parseUint32(value, key);
    } else if (key == "run_name") {
      metadata.run_name = value;
    } else if (key == "normalized_config_hash") {
      metadata.normalized_config_hash = parseUint64(value, key);
    } else if (key == "normalized_config_hash_hex") {
      metadata.normalized_config_hash_hex = value;
    } else if (key == "step_index") {
      metadata.step_index = parseUint64(value, key);
    } else if (key == "scale_factor") {
      metadata.scale_factor = parseDouble(value, key);
    } else if (key == "snapshot_stem") {
      metadata.snapshot_stem = value;
    } else if (key == "restart_stem") {
      metadata.restart_stem = value;
    }
  }

  return metadata;
}

void ModuleSidecarRegistry::upsert(ModuleSidecarBlock block) {
  if (block.module_name.empty()) {
    throw std::invalid_argument("ModuleSidecarRegistry.upsert: module_name cannot be empty");
  }
  m_sidecars[block.module_name] = std::move(block);
}

const ModuleSidecarBlock* ModuleSidecarRegistry::find(std::string_view module_name) const {
  const auto it = m_sidecars.find(std::string(module_name));
  if (it == m_sidecars.end()) {
    return nullptr;
  }
  return &it->second;
}

std::size_t ModuleSidecarRegistry::size() const noexcept { return m_sidecars.size(); }

std::vector<const ModuleSidecarBlock*> ModuleSidecarRegistry::blocksSortedByName() const {
  std::vector<const ModuleSidecarBlock*> ordered;
  ordered.reserve(m_sidecars.size());
  for (const auto& [name, block] : m_sidecars) {
    (void)name;
    ordered.push_back(&block);
  }
  std::sort(
      ordered.begin(),
      ordered.end(),
      [](const ModuleSidecarBlock* lhs, const ModuleSidecarBlock* rhs) {
        return lhs->module_name < rhs->module_name;
      });
  return ordered;
}

// Resize shared particle skeleton and metadata sidecars together.
void SimulationState::resizeParticles(std::size_t count) {
  particles.resize(count);
  particle_sidecar.resize(count);
}

// Resize cell gravity skeleton and gas thermodynamic sidecar together.
void SimulationState::resizeCells(std::size_t count) {
  cells.resize(count);
  gas_cells.resize(count);
}

// Resize patch descriptor table.
void SimulationState::resizePatches(std::size_t count) { patches.resize(count); }

bool SimulationState::validateUniqueParticleIds() const {
  // IDs are globally unique across all species in one simulation state.
  std::unordered_set<std::uint64_t> ids;
  ids.reserve(particle_sidecar.particle_id.size());
  for (const auto id : particle_sidecar.particle_id) {
    if (!ids.insert(id).second) {
      return false;
    }
  }
  return true;
}

// Recompute species-local/global index lookup tables after sidecar updates.
void SimulationState::rebuildSpeciesIndex() { particle_species_index.rebuild(particle_sidecar); }

ParticleTransferPacket SimulationState::packSpeciesTransferPacket(ParticleSpecies species_tag) const {
  // Explicit species pack path for MPI/device transfers.
  const auto indices = particle_species_index.globalIndices(species_tag);

  ParticleTransferPacket packet;
  packet.species = species_tag;
  packet.particle_id.resize(indices.size());
  packet.position_x_comoving.resize(indices.size());
  packet.position_y_comoving.resize(indices.size());
  packet.position_z_comoving.resize(indices.size());
  packet.velocity_x_peculiar.resize(indices.size());
  packet.velocity_y_peculiar.resize(indices.size());
  packet.velocity_z_peculiar.resize(indices.size());
  packet.mass_code.resize(indices.size());
  packet.time_bin.resize(indices.size());
  packet.owning_rank.resize(indices.size());

  for (std::size_t i = 0; i < indices.size(); ++i) {
    const auto source = indices[i];
    packet.particle_id[i] = particle_sidecar.particle_id[source];
    packet.position_x_comoving[i] = particles.position_x_comoving[source];
    packet.position_y_comoving[i] = particles.position_y_comoving[source];
    packet.position_z_comoving[i] = particles.position_z_comoving[source];
    packet.velocity_x_peculiar[i] = particles.velocity_x_peculiar[source];
    packet.velocity_y_peculiar[i] = particles.velocity_y_peculiar[source];
    packet.velocity_z_peculiar[i] = particles.velocity_z_peculiar[source];
    packet.mass_code[i] = particles.mass_code[source];
    packet.time_bin[i] = particles.time_bin[source];
    packet.owning_rank[i] = particle_sidecar.owning_rank[source];
  }

  return packet;
}

// Validate structural and semantic ownership constraints across all core SoA
// blocks and sidecars.
bool SimulationState::validateOwnershipInvariants() const {
  // Structural consistency of each SoA block.
  if (!particles.isConsistent() || !particle_sidecar.isConsistent() || !cells.isConsistent() ||
      !gas_cells.isConsistent() || !patches.isConsistent() || !star_particles.isConsistent() ||
      !black_holes.isConsistent() || !tracers.isConsistent()) {
    return false;
  }

  // Particle hot/cold arrays must share one ownership cardinality.
  if (particles.size() != particle_sidecar.size()) {
    return false;
  }

  // Gas hydro sidecars are cell-owned and must match cell cardinality.
  if (cells.size() != gas_cells.size()) {
    return false;
  }

  // Species ledger must match sidecar tags exactly.
  if (!species.isConsistentWith(particle_sidecar)) {
    return false;
  }

  // Every cell must reference a valid patch index when patches exist.
  for (std::size_t i = 0; i < cells.patch_index.size(); ++i) {
    if (cells.patch_index[i] >= patches.size() && patches.size() > 0) {
      return false;
    }
  }

  // Every patch-declared cell range must fall inside the global cell index space.
  for (std::size_t patch = 0; patch < patches.size(); ++patch) {
    const std::uint64_t begin = patches.first_cell[patch];
    const std::uint64_t count = patches.cell_count[patch];
    if (begin + count > cells.size()) {
      return false;
    }
  }

  // Species sidecars must reference global particles with matching species tags.
  for (std::size_t i = 0; i < star_particles.size(); ++i) {
    const auto index = star_particles.particle_index[i];
    if (index >= particles.size()) {
      return false;
    }
    if (particle_sidecar.species_tag[index] != static_cast<std::uint32_t>(ParticleSpecies::kStar)) {
      return false;
    }
  }

  for (std::size_t i = 0; i < black_holes.size(); ++i) {
    const auto index = black_holes.particle_index[i];
    if (index >= particles.size()) {
      return false;
    }
    if (particle_sidecar.species_tag[index] != static_cast<std::uint32_t>(ParticleSpecies::kBlackHole)) {
      return false;
    }
  }

  for (std::size_t i = 0; i < tracers.size(); ++i) {
    const auto index = tracers.particle_index[i];
    if (index >= particles.size()) {
      return false;
    }
    if (tracers.host_cell_index[i] >= cells.size() && cells.size() > 0) {
      return false;
    }
    if (tracers.mass_fraction_of_host[i] < 0.0 || tracers.last_host_mass_code[i] < 0.0) {
      return false;
    }
    if (particle_sidecar.species_tag[index] != static_cast<std::uint32_t>(ParticleSpecies::kTracer)) {
      return false;
    }
  }

  return true;
}

// Clear transient scheduler-selected active lists while keeping capacity.
void ActiveIndexSet::clear() {
  particle_indices.clear();
  cell_indices.clear();
}

// Return compact particle active-view row count.
std::size_t ParticleActiveView::size() const noexcept { return particle_id.size(); }

// Return compact cell active-view row count.
std::size_t CellActiveView::size() const noexcept { return center_x_comoving.size(); }

// Return gravity kernel compact view row count.
std::size_t GravityParticleKernelView::size() const noexcept { return particle_index.size(); }

// Return hydro kernel compact view row count.
std::size_t HydroCellKernelView::size() const noexcept { return cell_index.size(); }

bool ParticleReorderMap::isConsistent(std::size_t particle_count) const noexcept {
  if (old_to_new_index.size() != particle_count || new_to_old_index.size() != particle_count) {
    return false;
  }

  // Validate bijection: each new index is hit exactly once and inverts back.
  std::vector<std::uint8_t> visited(particle_count, 0U);
  for (std::size_t i = 0; i < particle_count; ++i) {
    const auto mapped = old_to_new_index[i];
    if (mapped >= particle_count) {
      return false;
    }
    if (visited[mapped] != 0U) {
      return false;
    }
    visited[mapped] = 1U;
    if (new_to_old_index[mapped] != i) {
      return false;
    }
  }
  return true;
}

// Construct monotonic scratch allocator with optional pre-reserved capacity.
MonotonicScratchAllocator::MonotonicScratchAllocator(std::size_t initial_capacity_bytes)
    : m_storage(initial_capacity_bytes), m_offset_bytes(0) {}

std::byte* MonotonicScratchAllocator::allocateBytes(std::size_t bytes, std::size_t alignment) {
  // Require power-of-two alignment to preserve standard aligned-address arithmetic.
  if (alignment == 0 || (alignment & (alignment - 1U)) != 0) {
    throw std::invalid_argument("MonotonicScratchAllocator.allocateBytes: alignment must be power-of-two");
  }

  if (bytes == 0) {
    return m_storage.data() + m_offset_bytes;
  }

  // Round offset upward to the requested alignment boundary.
  const std::size_t aligned_offset = (m_offset_bytes + alignment - 1U) & ~(alignment - 1U);
  const std::size_t required_size = aligned_offset + bytes;

  // Geometric growth to keep amortized allocation overhead low.
  if (required_size > m_storage.size()) {
    const std::size_t grow_size = std::max(required_size, std::max<std::size_t>(1024, m_storage.size() * 2));
    m_storage.resize(grow_size);
  }

  auto* ptr = m_storage.data() + aligned_offset;
  m_offset_bytes = required_size;
  return ptr;
}

// Reset bump pointer and retain allocated capacity for next step.
void MonotonicScratchAllocator::reset() { m_offset_bytes = 0; }

// Return currently allocated scratch arena capacity in bytes.
std::size_t MonotonicScratchAllocator::capacityBytes() const noexcept { return m_storage.size(); }

void TransientStepWorkspace::clear() {
  // Preserve capacity while dropping active step data.
  particle_id.clear();
  particle_species_tag.clear();
  particle_position_x_comoving.clear();
  particle_position_y_comoving.clear();
  particle_position_z_comoving.clear();
  particle_velocity_x_peculiar.clear();
  particle_velocity_y_peculiar.clear();
  particle_velocity_z_peculiar.clear();
  particle_mass_code.clear();
  gravity_particle_index.clear();

  hydro_cell_index.clear();
  hydro_cell_center_x_comoving.clear();
  hydro_cell_center_y_comoving.clear();
  hydro_cell_center_z_comoving.clear();
  hydro_cell_mass_code.clear();
  hydro_cell_density_code.clear();
  hydro_cell_pressure_code.clear();

  cell_center_x_comoving.clear();
  cell_center_y_comoving.clear();
  cell_center_z_comoving.clear();
  cell_mass_code.clear();
  cell_patch_index.clear();
  cell_density_code.clear();
  cell_pressure_code.clear();

  scratch.reset();
}

// Gather sparse particle indices into compact contiguous SoA buffers for
// vector-friendly active kernels.
ParticleActiveView buildParticleActiveView(
    const SimulationState& state,
    std::span<const std::uint32_t> active_particle_indices,
    TransientStepWorkspace& workspace) {
  // Materialize compact particle arrays for kernel-friendly contiguous iteration.
  workspace.particle_id.resize(active_particle_indices.size());
  workspace.particle_species_tag.resize(active_particle_indices.size());
  workspace.particle_position_x_comoving.resize(active_particle_indices.size());
  workspace.particle_position_y_comoving.resize(active_particle_indices.size());
  workspace.particle_position_z_comoving.resize(active_particle_indices.size());
  workspace.particle_velocity_x_peculiar.resize(active_particle_indices.size());
  workspace.particle_velocity_y_peculiar.resize(active_particle_indices.size());
  workspace.particle_velocity_z_peculiar.resize(active_particle_indices.size());
  workspace.particle_mass_code.resize(active_particle_indices.size());

  for (const std::uint32_t source : active_particle_indices) {
    // Keep an explicit bounds pass so release builds fail fast with a clear
    // exception rather than relying on debug-only checks in gatherSpan.
    if (source >= state.particles.size()) {
      throw std::out_of_range("buildParticleActiveView: particle index out of range");
    }
  }

  // Gather hot/cold lanes into compact contiguous buffers used by step-local kernels.
  gatherSpan<std::uint64_t>(state.particle_sidecar.particle_id, active_particle_indices, workspace.particle_id);
  gatherSpan<std::uint32_t>(state.particle_sidecar.species_tag, active_particle_indices, workspace.particle_species_tag);
  gatherSpan<double>(state.particles.position_x_comoving, active_particle_indices, workspace.particle_position_x_comoving);
  gatherSpan<double>(state.particles.position_y_comoving, active_particle_indices, workspace.particle_position_y_comoving);
  gatherSpan<double>(state.particles.position_z_comoving, active_particle_indices, workspace.particle_position_z_comoving);
  gatherSpan<double>(state.particles.velocity_x_peculiar, active_particle_indices, workspace.particle_velocity_x_peculiar);
  gatherSpan<double>(state.particles.velocity_y_peculiar, active_particle_indices, workspace.particle_velocity_y_peculiar);
  gatherSpan<double>(state.particles.velocity_z_peculiar, active_particle_indices, workspace.particle_velocity_z_peculiar);
  gatherSpan<double>(state.particles.mass_code, active_particle_indices, workspace.particle_mass_code);

  return ParticleActiveView{
      .particle_id = workspace.particle_id,
      .species_tag = workspace.particle_species_tag,
      .position_x_comoving = workspace.particle_position_x_comoving,
      .position_y_comoving = workspace.particle_position_y_comoving,
      .position_z_comoving = workspace.particle_position_z_comoving,
      .velocity_x_peculiar = workspace.particle_velocity_x_peculiar,
      .velocity_y_peculiar = workspace.particle_velocity_y_peculiar,
      .velocity_z_peculiar = workspace.particle_velocity_z_peculiar,
      .mass_code = workspace.particle_mass_code,
  };
}

// Gather sparse cell indices into compact contiguous buffers, combining cell
// skeleton fields with gas thermodynamic sidecar lanes.
CellActiveView buildCellActiveView(
    const SimulationState& state,
    std::span<const std::uint32_t> active_cell_indices,
    TransientStepWorkspace& workspace) {
  // Materialize compact cell arrays for kernel-friendly contiguous iteration.
  workspace.cell_center_x_comoving.resize(active_cell_indices.size());
  workspace.cell_center_y_comoving.resize(active_cell_indices.size());
  workspace.cell_center_z_comoving.resize(active_cell_indices.size());
  workspace.cell_mass_code.resize(active_cell_indices.size());
  workspace.cell_patch_index.resize(active_cell_indices.size());
  workspace.cell_density_code.resize(active_cell_indices.size());
  workspace.cell_pressure_code.resize(active_cell_indices.size());

  for (const std::uint32_t source : active_cell_indices) {
    // Mirror particle-view policy: explicit range validation before gather.
    if (source >= state.cells.size()) {
      throw std::out_of_range("buildCellActiveView: cell index out of range");
    }
  }

  // Materialize compact cell lanes for branch-light active-set sweeps.
  gatherSpan<double>(state.cells.center_x_comoving, active_cell_indices, workspace.cell_center_x_comoving);
  gatherSpan<double>(state.cells.center_y_comoving, active_cell_indices, workspace.cell_center_y_comoving);
  gatherSpan<double>(state.cells.center_z_comoving, active_cell_indices, workspace.cell_center_z_comoving);
  gatherSpan<double>(state.cells.mass_code, active_cell_indices, workspace.cell_mass_code);
  gatherSpan<std::uint32_t>(state.cells.patch_index, active_cell_indices, workspace.cell_patch_index);
  gatherSpan<double>(state.gas_cells.density_code, active_cell_indices, workspace.cell_density_code);
  gatherSpan<double>(state.gas_cells.pressure_code, active_cell_indices, workspace.cell_pressure_code);

  return CellActiveView{
      .center_x_comoving = workspace.cell_center_x_comoving,
      .center_y_comoving = workspace.cell_center_y_comoving,
      .center_z_comoving = workspace.cell_center_z_comoving,
      .mass_code = workspace.cell_mass_code,
      .patch_index = workspace.cell_patch_index,
      .density_code = workspace.cell_density_code,
      .pressure_code = workspace.cell_pressure_code,
  };
}

GravityParticleKernelView buildGravityParticleKernelView(
    const SimulationState& state,
    std::span<const std::uint32_t> active_particle_indices,
    TransientStepWorkspace& workspace) {
  // Keep source indices alongside compact lanes so write-back is explicit.
  workspace.gravity_particle_index.resize(active_particle_indices.size());
  workspace.particle_position_x_comoving.resize(active_particle_indices.size());
  workspace.particle_position_y_comoving.resize(active_particle_indices.size());
  workspace.particle_position_z_comoving.resize(active_particle_indices.size());
  workspace.particle_velocity_x_peculiar.resize(active_particle_indices.size());
  workspace.particle_velocity_y_peculiar.resize(active_particle_indices.size());
  workspace.particle_velocity_z_peculiar.resize(active_particle_indices.size());
  workspace.particle_mass_code.resize(active_particle_indices.size());

  for (std::size_t i = 0; i < active_particle_indices.size(); ++i) {
    const auto source = active_particle_indices[i];
    if (source >= state.particles.size()) {
      throw std::out_of_range("buildGravityParticleKernelView: particle index out of range");
    }
    workspace.gravity_particle_index[i] = source;
  }

  gatherSpan<double>(
      state.particles.position_x_comoving,
      active_particle_indices,
      workspace.particle_position_x_comoving);
  gatherSpan<double>(
      state.particles.position_y_comoving,
      active_particle_indices,
      workspace.particle_position_y_comoving);
  gatherSpan<double>(
      state.particles.position_z_comoving,
      active_particle_indices,
      workspace.particle_position_z_comoving);
  gatherSpan<double>(
      state.particles.velocity_x_peculiar,
      active_particle_indices,
      workspace.particle_velocity_x_peculiar);
  gatherSpan<double>(
      state.particles.velocity_y_peculiar,
      active_particle_indices,
      workspace.particle_velocity_y_peculiar);
  gatherSpan<double>(
      state.particles.velocity_z_peculiar,
      active_particle_indices,
      workspace.particle_velocity_z_peculiar);
  gatherSpan<double>(state.particles.mass_code, active_particle_indices, workspace.particle_mass_code);

  return GravityParticleKernelView{
      .particle_index = workspace.gravity_particle_index,
      .position_x_comoving = workspace.particle_position_x_comoving,
      .position_y_comoving = workspace.particle_position_y_comoving,
      .position_z_comoving = workspace.particle_position_z_comoving,
      .velocity_x_peculiar = workspace.particle_velocity_x_peculiar,
      .velocity_y_peculiar = workspace.particle_velocity_y_peculiar,
      .velocity_z_peculiar = workspace.particle_velocity_z_peculiar,
      .mass_code = workspace.particle_mass_code,
  };
}

void scatterGravityParticleKernelView(const GravityParticleKernelView& view, SimulationState& state) {
  // Scatter changed compact lanes back to persistent hot arrays using stored indices.
  for (std::size_t i = 0; i < view.size(); ++i) {
    const auto destination = view.particle_index[i];
    if (destination >= state.particles.size()) {
      throw std::out_of_range("scatterGravityParticleKernelView: stale particle index");
    }
    state.particles.position_x_comoving[destination] = view.position_x_comoving[i];
    state.particles.position_y_comoving[destination] = view.position_y_comoving[i];
    state.particles.position_z_comoving[destination] = view.position_z_comoving[i];
    state.particles.velocity_x_peculiar[destination] = view.velocity_x_peculiar[i];
    state.particles.velocity_y_peculiar[destination] = view.velocity_y_peculiar[i];
    state.particles.velocity_z_peculiar[destination] = view.velocity_z_peculiar[i];
    state.particles.mass_code[destination] = view.mass_code[i];
  }
}

HydroCellKernelView buildHydroCellKernelView(
    const SimulationState& state,
    std::span<const std::uint32_t> active_cell_indices,
    TransientStepWorkspace& workspace) {
  // Preserve source cell indices so hydro writes remain explicit and auditable.
  workspace.hydro_cell_index.resize(active_cell_indices.size());
  workspace.hydro_cell_center_x_comoving.resize(active_cell_indices.size());
  workspace.hydro_cell_center_y_comoving.resize(active_cell_indices.size());
  workspace.hydro_cell_center_z_comoving.resize(active_cell_indices.size());
  workspace.hydro_cell_mass_code.resize(active_cell_indices.size());
  workspace.hydro_cell_density_code.resize(active_cell_indices.size());
  workspace.hydro_cell_pressure_code.resize(active_cell_indices.size());

  for (std::size_t i = 0; i < active_cell_indices.size(); ++i) {
    const auto source = active_cell_indices[i];
    if (source >= state.cells.size()) {
      throw std::out_of_range("buildHydroCellKernelView: cell index out of range");
    }
    workspace.hydro_cell_index[i] = source;
  }

  gatherSpan<double>(
      state.cells.center_x_comoving,
      active_cell_indices,
      workspace.hydro_cell_center_x_comoving);
  gatherSpan<double>(
      state.cells.center_y_comoving,
      active_cell_indices,
      workspace.hydro_cell_center_y_comoving);
  gatherSpan<double>(
      state.cells.center_z_comoving,
      active_cell_indices,
      workspace.hydro_cell_center_z_comoving);
  gatherSpan<double>(state.cells.mass_code, active_cell_indices, workspace.hydro_cell_mass_code);
  gatherSpan<double>(state.gas_cells.density_code, active_cell_indices, workspace.hydro_cell_density_code);
  gatherSpan<double>(state.gas_cells.pressure_code, active_cell_indices, workspace.hydro_cell_pressure_code);

  return HydroCellKernelView{
      .cell_index = workspace.hydro_cell_index,
      .center_x_comoving = workspace.hydro_cell_center_x_comoving,
      .center_y_comoving = workspace.hydro_cell_center_y_comoving,
      .center_z_comoving = workspace.hydro_cell_center_z_comoving,
      .mass_code = workspace.hydro_cell_mass_code,
      .density_code = workspace.hydro_cell_density_code,
      .pressure_code = workspace.hydro_cell_pressure_code,
  };
}

void scatterHydroCellKernelView(const HydroCellKernelView& view, SimulationState& state) {
  // Scatter compact hydro outputs back into persistent cell/gas storage.
  for (std::size_t i = 0; i < view.size(); ++i) {
    const auto destination = view.cell_index[i];
    if (destination >= state.cells.size()) {
      throw std::out_of_range("scatterHydroCellKernelView: stale cell index");
    }
    state.cells.center_x_comoving[destination] = view.center_x_comoving[i];
    state.cells.center_y_comoving[destination] = view.center_y_comoving[i];
    state.cells.center_z_comoving[destination] = view.center_z_comoving[i];
    state.cells.mass_code[destination] = view.mass_code[i];
    state.gas_cells.density_code[destination] = view.density_code[i];
    state.gas_cells.pressure_code[destination] = view.pressure_code[i];
  }
}

ParticleReorderMap buildParticleReorderMap(const SimulationState& state, ParticleReorderMode mode) {
  ParticleReorderMap reorder_map;
  reorder_map.new_to_old_index.resize(state.particles.size());
  std::iota(reorder_map.new_to_old_index.begin(), reorder_map.new_to_old_index.end(), 0U);

  // Comparator intentionally ties on old index to keep deterministic stability.
  const auto key_comp = [&](std::uint32_t lhs, std::uint32_t rhs) {
    if (mode == ParticleReorderMode::kByTimeBin) {
      const auto lhs_key = state.particles.time_bin[lhs];
      const auto rhs_key = state.particles.time_bin[rhs];
      return std::tuple{lhs_key, lhs} < std::tuple{rhs_key, rhs};
    }
    if (mode == ParticleReorderMode::kBySfcKey) {
      const auto lhs_key = state.particle_sidecar.sfc_key[lhs];
      const auto rhs_key = state.particle_sidecar.sfc_key[rhs];
      return std::tuple{lhs_key, lhs} < std::tuple{rhs_key, rhs};
    }
    const auto lhs_key = state.particle_sidecar.species_tag[lhs];
    const auto rhs_key = state.particle_sidecar.species_tag[rhs];
    return std::tuple{lhs_key, lhs} < std::tuple{rhs_key, rhs};
  };

  std::stable_sort(reorder_map.new_to_old_index.begin(), reorder_map.new_to_old_index.end(), key_comp);

  reorder_map.old_to_new_index.resize(state.particles.size());
  for (std::size_t new_index = 0; new_index < reorder_map.new_to_old_index.size(); ++new_index) {
    const auto old_index = reorder_map.new_to_old_index[new_index];
    reorder_map.old_to_new_index[old_index] = static_cast<std::uint32_t>(new_index);
  }

  return reorder_map;
}

template <typename T>
void reorderAlignedVector(
    AlignedVector<T>& values,
    std::span<const std::uint32_t> new_to_old_index) {
  // Rebuild destination vector by reading from old rows in new order.
  AlignedVector<T> reordered(values.size());
  for (std::size_t i = 0; i < new_to_old_index.size(); ++i) {
    reordered[i] = values[new_to_old_index[i]];
  }
  values.swap(reordered);
}

void reorderParticles(
    SimulationState& state,
    const ParticleReorderMap& reorder_map,
    const SidecarSyncPolicy& sync_policy) {
  if (!reorder_map.isConsistent(state.particles.size())) {
    throw std::invalid_argument("reorderParticles: inconsistent reorder map");
  }

  // Single permutation drives all parent-owned hot and sidecar lanes.
  const std::span<const std::uint32_t> new_to_old_index = reorder_map.new_to_old_index;

  reorderAlignedVector(state.particles.position_x_comoving, new_to_old_index);
  reorderAlignedVector(state.particles.position_y_comoving, new_to_old_index);
  reorderAlignedVector(state.particles.position_z_comoving, new_to_old_index);
  reorderAlignedVector(state.particles.velocity_x_peculiar, new_to_old_index);
  reorderAlignedVector(state.particles.velocity_y_peculiar, new_to_old_index);
  reorderAlignedVector(state.particles.velocity_z_peculiar, new_to_old_index);
  reorderAlignedVector(state.particles.mass_code, new_to_old_index);
  reorderAlignedVector(state.particles.time_bin, new_to_old_index);

  reorderAlignedVector(state.particle_sidecar.particle_id, new_to_old_index);
  reorderAlignedVector(state.particle_sidecar.sfc_key, new_to_old_index);
  reorderAlignedVector(state.particle_sidecar.species_tag, new_to_old_index);
  reorderAlignedVector(state.particle_sidecar.particle_flags, new_to_old_index);
  reorderAlignedVector(state.particle_sidecar.owning_rank, new_to_old_index);

  // Per-sidecar policy: either move rows with parents or keep rows and remap references.
  auto remap_sidecar_index = [&](AlignedVector<std::uint32_t>& particle_index, SidecarSyncMode mode) {
    if (mode == SidecarSyncMode::kMoveWithParent) {
      reorderAlignedVector(particle_index, new_to_old_index);
      return;
    }
    for (auto& index : particle_index) {
      if (index >= reorder_map.old_to_new_index.size()) {
        throw std::out_of_range("reorderParticles: sidecar particle index out of range");
      }
      index = reorder_map.old_to_new_index[index];
    }
  };

  remap_sidecar_index(state.star_particles.particle_index, sync_policy.star_particles);
  remap_sidecar_index(state.black_holes.particle_index, sync_policy.black_holes);
  remap_sidecar_index(state.tracers.particle_index, sync_policy.tracers);

  state.rebuildSpeciesIndex();
}

void debugAssertNoStaleParticleIndices(const SimulationState& state) {
  // Human-readable source labels keep invariant failures easy to triage in tests/CI.
  auto check_indices = [&](std::span<const std::uint32_t> indices, const char* name) {
    for (const auto index : indices) {
      if (index >= state.particles.size()) {
        throw std::runtime_error(std::string("debugAssertNoStaleParticleIndices: stale index in ") + name);
      }
    }
  };

  check_indices(state.star_particles.particle_index, "star_particles");
  check_indices(state.black_holes.particle_index, "black_holes");
  check_indices(state.tracers.particle_index, "tracers");
}

}  // namespace cosmosim::core
