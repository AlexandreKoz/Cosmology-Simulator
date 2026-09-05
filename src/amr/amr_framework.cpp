#include "cosmosim/amr/amr_framework.hpp"

#include "cosmosim/core/checked_arithmetic.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <utility>

namespace cosmosim::amr {

namespace {

[[nodiscard]] std::size_t cellIndex(
    std::array<std::uint16_t, 3> dims,
    std::uint16_t x,
    std::uint16_t y,
    std::uint16_t z) {
  return (static_cast<std::size_t>(z) * dims[1] + y) * dims[0] + x;
}

[[nodiscard]] CellMetrics prolongatedChildMetrics(
    const CellMetrics& parent_metrics,
    std::uint8_t fine_ordinal_within_parent) {
  CellMetrics child_metrics = parent_metrics;
  child_metrics.particle_count =
      parent_metrics.particle_count / 8U +
      (fine_ordinal_within_parent < (parent_metrics.particle_count % 8U) ? 1U : 0U);
  return child_metrics;
}

[[nodiscard]] CellMetrics restrictedParentMetrics(
    std::span<const CellMetrics> fine_metrics,
    const ConservedState& restricted_conserved,
    double parent_cell_volume_comov) {
  CellMetrics metrics;
  if (parent_cell_volume_comov > 0.0) {
    metrics.density_code = restricted_conserved.mass_code / parent_cell_volume_comov;
  }
  if (fine_metrics.empty()) {
    return metrics;
  }

  double pressure_sum = 0.0;
  double sound_speed_sum = 0.0;
  double gradient_indicator_max = 0.0;
  std::uint64_t particle_count_sum = 0;
  for (const auto& fine : fine_metrics) {
    pressure_sum += fine.pressure_code;
    sound_speed_sum += fine.sound_speed_code;
    gradient_indicator_max = std::max(gradient_indicator_max, fine.gradient_indicator);
    particle_count_sum += fine.particle_count;
  }

  const double inv_count = 1.0 / static_cast<double>(fine_metrics.size());
  metrics.pressure_code = pressure_sum * inv_count;
  metrics.sound_speed_code = sound_speed_sum * inv_count;
  metrics.gradient_indicator = gradient_indicator_max;
  metrics.particle_count = static_cast<std::uint32_t>(
      std::min<std::uint64_t>(particle_count_sum, std::numeric_limits<std::uint32_t>::max()));
  return metrics;
}

[[nodiscard]] ConservedState conservedFromHydroFlux(const hydro::HydroConservedState& flux) {
  return ConservedState{
      .mass_code = flux.mass_density_comoving,
      .momentum_x_code = flux.momentum_density_x_comoving,
      .momentum_y_code = flux.momentum_density_y_comoving,
      .momentum_z_code = flux.momentum_density_z_comoving,
      .total_energy_code = flux.total_energy_density_comoving,
      .metal_mass_code = flux.metal_mass_density_comoving};
}

void validateCompatibleRegisterRecord(
    const FluxRegisterEntry& entry,
    const hydro::HydroFluxRegisterRecord& record) {
  if (entry.coarse_patch_id != record.coarse_patch_id ||
      (entry.coarse_gas_cell_id != 0U && record.coarse_gas_cell_id != 0U &&
       entry.coarse_gas_cell_id != record.coarse_gas_cell_id) ||
      entry.coarse_cell_index != record.coarse_cell_index ||
      entry.level != static_cast<std::uint8_t>(record.level) ||
      entry.axis != record.axis ||
      entry.orientation != record.orientation) {
    throw std::invalid_argument("FluxRegisterAccumulator received incompatible metadata for register key");
  }
  if (std::abs(entry.dt_code - record.dt_code) > 1.0e-14) {
    throw std::invalid_argument("FluxRegisterAccumulator requires one dt_code per register key");
  }
}

[[nodiscard]] ConservedState averageFluxOrZero(
    const ConservedState& area_weighted_flux,
    double area_comov) {
  if (area_comov <= 0.0) {
    return {};
  }
  ConservedState flux = area_weighted_flux;
  flux *= 1.0 / area_comov;
  return flux;
}

}  // namespace

ConservedState& ConservedState::operator+=(const ConservedState& rhs) {
  mass_code += rhs.mass_code;
  momentum_x_code += rhs.momentum_x_code;
  momentum_y_code += rhs.momentum_y_code;
  momentum_z_code += rhs.momentum_z_code;
  total_energy_code += rhs.total_energy_code;
  metal_mass_code += rhs.metal_mass_code;
  return *this;
}

ConservedState& ConservedState::operator-=(const ConservedState& rhs) {
  mass_code -= rhs.mass_code;
  momentum_x_code -= rhs.momentum_x_code;
  momentum_y_code -= rhs.momentum_y_code;
  momentum_z_code -= rhs.momentum_z_code;
  total_energy_code -= rhs.total_energy_code;
  metal_mass_code -= rhs.metal_mass_code;
  return *this;
}

ConservedState& ConservedState::operator*=(double factor) {
  mass_code *= factor;
  momentum_x_code *= factor;
  momentum_y_code *= factor;
  momentum_z_code *= factor;
  total_energy_code *= factor;
  metal_mass_code *= factor;
  return *this;
}

ConservedState operator+(ConservedState lhs, const ConservedState& rhs) {
  lhs += rhs;
  return lhs;
}

ConservedState operator-(ConservedState lhs, const ConservedState& rhs) {
  lhs -= rhs;
  return lhs;
}

ConservedState operator*(ConservedState lhs, double factor) {
  lhs *= factor;
  return lhs;
}

void FluxRegisterAccumulator::recordFaceFlux(const hydro::HydroFluxRegisterRecord& record) {
  if (record.role == hydro::HydroFluxRegisterFaceRole::kNone) {
    return;
  }
  if (record.register_key == 0U ||
      record.coarse_patch_id == 0U ||
      record.coarse_cell_index == hydro::k_invalid_cell_index ||
      record.face_area_comoving <= 0.0 ||
      record.dt_code <= 0.0 ||
      record.level < 0 ||
      record.level > static_cast<int>(std::numeric_limits<std::uint8_t>::max())) {
    throw std::invalid_argument("FluxRegisterAccumulator received invalid hydro flux-register record");
  }

  std::size_t slot = 0;
  const auto found = m_slot_by_key.find(record.register_key);
  if (found == m_slot_by_key.end()) {
    slot = m_entries.size();
    m_slot_by_key.emplace(record.register_key, slot);
    AccumulatedEntry accumulated;
    accumulated.entry.register_key = record.register_key;
    accumulated.entry.coarse_patch_id = record.coarse_patch_id;
    accumulated.entry.coarse_gas_cell_id = record.coarse_gas_cell_id;
    accumulated.entry.coarse_cell_index = record.coarse_cell_index;
    accumulated.entry.level = static_cast<std::uint8_t>(record.level);
    accumulated.entry.axis = record.axis;
    accumulated.entry.orientation = record.orientation;
    accumulated.entry.face_area_comov = record.face_area_comoving;
    accumulated.entry.dt_code = record.dt_code;
    m_entries.push_back(accumulated);
  } else {
    slot = found->second;
    validateCompatibleRegisterRecord(m_entries[slot].entry, record);
  }

  AccumulatedEntry& accumulated = m_entries[slot];
  if (accumulated.entry.coarse_gas_cell_id == 0U && record.coarse_gas_cell_id != 0U) {
    accumulated.entry.coarse_gas_cell_id = record.coarse_gas_cell_id;
  }
  const ConservedState area_weighted_flux =
      conservedFromHydroFlux(record.flux_code) * record.face_area_comoving;
  if (record.role == hydro::HydroFluxRegisterFaceRole::kCoarse) {
    accumulated.coarse_area_weighted_flux += area_weighted_flux;
    accumulated.coarse_area_comov += record.face_area_comoving;
    accumulated.entry.coarse_face_count += 1U;
    accumulated.entry.face_area_comov = accumulated.coarse_area_comov;
  } else if (record.role == hydro::HydroFluxRegisterFaceRole::kFine) {
    accumulated.fine_area_weighted_flux += area_weighted_flux;
    accumulated.fine_area_comov += record.face_area_comoving;
    accumulated.entry.fine_face_count += 1U;
    if (accumulated.coarse_area_comov <= 0.0) {
      accumulated.entry.face_area_comov = accumulated.fine_area_comov;
    }
  }
}

std::vector<FluxRegisterEntry> FluxRegisterAccumulator::entries() const {
  std::vector<FluxRegisterEntry> result;
  result.reserve(m_entries.size());
  for (const AccumulatedEntry& accumulated : m_entries) {
    FluxRegisterEntry entry = accumulated.entry;
    const double register_area = accumulated.coarse_area_comov > 0.0
        ? accumulated.coarse_area_comov
        : accumulated.fine_area_comov;
    entry.face_area_comov = register_area;
    entry.coarse_area_comov = accumulated.coarse_area_comov;
    entry.fine_area_comov = accumulated.fine_area_comov;
    entry.coarse_face_flux_code =
        averageFluxOrZero(accumulated.coarse_area_weighted_flux, register_area);
    entry.fine_face_flux_code =
        averageFluxOrZero(accumulated.fine_area_weighted_flux, register_area);
    result.push_back(entry);
  }
  std::sort(result.begin(), result.end(), [](const FluxRegisterEntry& lhs, const FluxRegisterEntry& rhs) {
    return lhs.register_key < rhs.register_key;
  });
  return result;
}

void FluxRegisterAccumulator::clear() {
  m_slot_by_key.clear();
  m_entries.clear();
}

RefinementDecision RefinementEvaluator::evaluateCell(
    const CellMetrics& metrics,
    double cell_width_comov,
    const RefinementCriteria& criteria) {
  bool trigger_refine = false;
  bool allow_derefine = true;

  if (criteria.use_mass_threshold) {
    const double cell_mass_code = metrics.density_code * std::pow(cell_width_comov, 3);
    trigger_refine = trigger_refine || (cell_mass_code > criteria.mass_threshold_code);
    allow_derefine = allow_derefine &&
                     (cell_mass_code <= criteria.mass_threshold_code * criteria.derefine_hysteresis);
  }

  if (criteria.use_gradient_indicator) {
    trigger_refine = trigger_refine || (metrics.gradient_indicator > criteria.gradient_threshold);
    allow_derefine = allow_derefine &&
                     (metrics.gradient_indicator <= criteria.gradient_threshold * criteria.derefine_hysteresis);
  }

  if (criteria.use_particle_count) {
    trigger_refine = trigger_refine || (metrics.particle_count > criteria.particle_threshold);
    const double relaxed_threshold =
        static_cast<double>(criteria.particle_threshold) * criteria.derefine_hysteresis;
    allow_derefine = allow_derefine &&
                     (static_cast<double>(metrics.particle_count) <= relaxed_threshold);
  }

  if (criteria.use_jeans_condition &&
      metrics.sound_speed_code > 0.0 &&
      criteria.gravitational_constant_code > 0.0 &&
      metrics.density_code > 0.0) {
    const double jeans_length_comov =
        metrics.sound_speed_code *
        std::sqrt(std::numbers::pi_v<double> / (criteria.gravitational_constant_code * metrics.density_code));
    const double jeans_resolution = jeans_length_comov / cell_width_comov;
    trigger_refine = trigger_refine || (jeans_resolution < criteria.jeans_resolution_cells);
    allow_derefine = allow_derefine &&
                     (jeans_resolution >= criteria.jeans_resolution_cells / criteria.derefine_hysteresis);
  }

  if (trigger_refine) {
    return RefinementDecision::kRefine;
  }
  if (allow_derefine) {
    return RefinementDecision::kDerefine;
  }
  return RefinementDecision::kKeep;
}

AmrPatch::AmrPatch(PatchDescriptor descriptor)
    : m_descriptor(std::move(descriptor)) {
  if (m_descriptor.cell_dims[0] == 0U || m_descriptor.cell_dims[1] == 0U ||
      m_descriptor.cell_dims[2] == 0U) {
    throw std::invalid_argument("AMR patch cell dimensions must be nonzero");
  }
  const std::size_t cell_count = core::checkedSizeProduct3(
      static_cast<std::size_t>(m_descriptor.cell_dims[0]),
      static_cast<std::size_t>(m_descriptor.cell_dims[1]),
      static_cast<std::size_t>(m_descriptor.cell_dims[2]),
      "AMR patch cell count");
  m_conserved.resize(cell_count);
  m_metrics.resize(cell_count);
  m_gas_cell_ids.resize(cell_count);
}

const PatchDescriptor& AmrPatch::descriptor() const {
  return m_descriptor;
}

std::size_t AmrPatch::cellCount() const {
  return m_conserved.size();
}

double AmrPatch::cellVolumeComov() const {
  return (m_descriptor.extent_comov[0] * m_descriptor.extent_comov[1] * m_descriptor.extent_comov[2]) /
         static_cast<double>(cellCount());
}

double AmrPatch::cellWidthComov() const {
  return m_descriptor.extent_comov[0] / static_cast<double>(m_descriptor.cell_dims[0]);
}

std::span<ConservedState> AmrPatch::conservedView() {
  return std::span<ConservedState>(m_conserved.data(), m_conserved.size());
}

std::span<const ConservedState> AmrPatch::conservedView() const {
  return std::span<const ConservedState>(m_conserved.data(), m_conserved.size());
}

std::span<CellMetrics> AmrPatch::metricsView() {
  return std::span<CellMetrics>(m_metrics.data(), m_metrics.size());
}

std::span<const CellMetrics> AmrPatch::metricsView() const {
  return std::span<const CellMetrics>(m_metrics.data(), m_metrics.size());
}

std::span<std::uint64_t> AmrPatch::gasCellIdView() {
  return std::span<std::uint64_t>(m_gas_cell_ids.data(), m_gas_cell_ids.size());
}

std::span<const std::uint64_t> AmrPatch::gasCellIdView() const {
  return std::span<const std::uint64_t>(m_gas_cell_ids.data(), m_gas_cell_ids.size());
}

ConservedState AmrPatch::totalConserved() const {
  ConservedState total;
  for (const auto& cell : m_conserved) {
    total += cell;
  }
  return total;
}

std::size_t AmrPatch::ownedStorageCapacityBytes() const {
  std::size_t bytes = core::checkedSizeMultiply(
      m_conserved.capacity(), sizeof(ConservedState), "AMR conserved-state capacity bytes");
  bytes = core::checkedSizeAdd(
      bytes,
      core::checkedSizeMultiply(m_metrics.capacity(), sizeof(CellMetrics), "AMR cell-metrics capacity bytes"),
      "AMR patch storage capacity bytes");
  bytes = core::checkedSizeAdd(
      bytes,
      core::checkedSizeMultiply(
          m_gas_cell_ids.capacity(), sizeof(std::uint64_t), "AMR gas-cell-id capacity bytes"),
      "AMR patch storage capacity bytes");
  return bytes;
}

bool AmrPatch::isLeaf() const {
  return m_is_leaf;
}

void AmrPatch::setLeaf(bool is_leaf) {
  m_is_leaf = is_leaf;
}

std::uint64_t PatchHierarchy::createRootPatch(const PatchDescriptor& root) {
  if (!m_levels.empty()) {
    throw std::runtime_error("PatchHierarchy already initialized with a root patch.");
  }
  if (m_next_patch_id == std::numeric_limits<std::uint64_t>::max()) {
    throw std::overflow_error("AMR patch id space exhausted");
  }

  PatchDescriptor root_copy = root;
  root_copy.patch_id = m_next_patch_id;
  root_copy.parent_patch_id = 0;
  root_copy.level = 0;
  AmrPatch prepared_root(root_copy);
  std::uint64_t candidate_next_gas_cell_id = m_next_gas_cell_id;
  assignStableGasCellIds(prepared_root, candidate_next_gas_cell_id);

  std::unordered_map<std::uint64_t, std::pair<std::size_t, std::size_t>> candidate_index;
  candidate_index.reserve(1U);
  candidate_index.emplace(root_copy.patch_id, std::pair<std::size_t, std::size_t>{0U, 0U});
  std::vector<AmrPatch> root_level;
  root_level.reserve(1U);
  root_level.emplace_back(std::move(prepared_root));
  m_levels.reserve(1U);

  m_levels.emplace_back(std::move(root_level));
  m_patch_index.swap(candidate_index);
  m_next_patch_id = root_copy.patch_id + 1U;
  m_next_gas_cell_id = candidate_next_gas_cell_id;
  return root_copy.patch_id;
}

std::array<std::uint64_t, 8> PatchHierarchy::refinePatch(std::uint64_t patch_id) {
  const AmrPatch* parent_patch = findPatch(patch_id);
  if (parent_patch == nullptr) {
    throw std::runtime_error("Cannot refine missing patch id.");
  }
  if (!parent_patch->isLeaf()) {
    throw std::runtime_error("Cannot refine non-leaf patch id.");
  }

  const PatchDescriptor parent_descriptor = parent_patch->descriptor();
  const std::span<const ConservedState> parent_conserved = parent_patch->conservedView();
  const std::span<const CellMetrics> parent_metrics = parent_patch->metricsView();
  const std::size_t child_level = static_cast<std::size_t>(parent_descriptor.level) + 1U;
  if (child_level > static_cast<std::size_t>(std::numeric_limits<std::uint8_t>::max())) {
    throw std::overflow_error("AMR refinement level exceeds uint8_t range");
  }
  if (m_next_patch_id > std::numeric_limits<std::uint64_t>::max() - 8U) {
    throw std::overflow_error("AMR patch id space cannot allocate a child octet");
  }

  const std::array<double, 3> child_extent = {
      parent_descriptor.extent_comov[0] * 0.5,
      parent_descriptor.extent_comov[1] * 0.5,
      parent_descriptor.extent_comov[2] * 0.5,
  };

  std::array<std::uint64_t, 8> child_ids{};
  std::vector<AmrPatch> prepared_children;
  prepared_children.reserve(8U);
  std::uint64_t candidate_next_patch_id = m_next_patch_id;
  std::uint64_t candidate_next_gas_cell_id = m_next_gas_cell_id;
  for (std::uint8_t octant = 0; octant < 8; ++octant) {
    PatchDescriptor child;
    child.patch_id = candidate_next_patch_id++;
    child.parent_patch_id = parent_descriptor.patch_id;
    child.level = static_cast<std::uint8_t>(child_level);
    child.morton_key = (parent_descriptor.morton_key << 3U) | octant;
    child.cell_dims = parent_descriptor.cell_dims;
    child.extent_comov = child_extent;
    child.origin_comov = parent_descriptor.origin_comov;
    child.origin_comov[0] += ((octant & 1U) != 0U) ? child_extent[0] : 0.0;
    child.origin_comov[1] += ((octant & 2U) != 0U) ? child_extent[1] : 0.0;
    child.origin_comov[2] += ((octant & 4U) != 0U) ? child_extent[2] : 0.0;

    AmrPatch child_patch(child);
    assignStableGasCellIds(child_patch, candidate_next_gas_cell_id);
    auto child_conserved = child_patch.conservedView();
    auto child_metrics = child_patch.metricsView();
    for (std::uint16_t z = 0; z < child.cell_dims[2]; ++z) {
      for (std::uint16_t y = 0; y < child.cell_dims[1]; ++y) {
        for (std::uint16_t x = 0; x < child.cell_dims[0]; ++x) {
          const std::uint16_t fine_x =
              static_cast<std::uint16_t>(x + (((octant & 1U) != 0U) ? child.cell_dims[0] : 0U));
          const std::uint16_t fine_y =
              static_cast<std::uint16_t>(y + (((octant & 2U) != 0U) ? child.cell_dims[1] : 0U));
          const std::uint16_t fine_z =
              static_cast<std::uint16_t>(z + (((octant & 4U) != 0U) ? child.cell_dims[2] : 0U));
          const std::uint16_t parent_x = static_cast<std::uint16_t>(fine_x / 2U);
          const std::uint16_t parent_y = static_cast<std::uint16_t>(fine_y / 2U);
          const std::uint16_t parent_z = static_cast<std::uint16_t>(fine_z / 2U);
          const std::size_t parent_index =
              cellIndex(parent_descriptor.cell_dims, parent_x, parent_y, parent_z);
          const std::size_t child_index = cellIndex(child.cell_dims, x, y, z);
          const std::uint8_t fine_ordinal =
              static_cast<std::uint8_t>((fine_x & 1U) | ((fine_y & 1U) << 1U) | ((fine_z & 1U) << 2U));

          child_conserved[child_index] = parent_conserved[parent_index];
          child_conserved[child_index] *= 0.125;
          child_metrics[child_index] = prolongatedChildMetrics(parent_metrics[parent_index], fine_ordinal);
        }
      }
    }
    child_ids[octant] = child.patch_id;
    prepared_children.emplace_back(std::move(child_patch));
  }

  // Pre-allocate all hierarchy/index growth before changing leaf authority.
  m_levels.reserve(core::checkedSizeAdd(child_level, 1U, "AMR hierarchy level reserve"));
  if (child_level < m_levels.size()) {
    m_levels[child_level].reserve(core::checkedSizeAdd(
        m_levels[child_level].size(), prepared_children.size(), "AMR child-level patch reserve"));
  }

  std::vector<std::pair<std::uint64_t, std::uint64_t>> child_level_order;
  const std::size_t existing_child_level_size = child_level < m_levels.size() ? m_levels[child_level].size() : 0U;
  child_level_order.reserve(core::checkedSizeAdd(
      existing_child_level_size, prepared_children.size(), "AMR child-level order reserve"));
  if (child_level < m_levels.size()) {
    for (const AmrPatch& patch : m_levels[child_level]) {
      child_level_order.emplace_back(patch.descriptor().morton_key, patch.descriptor().patch_id);
    }
  }
  for (const AmrPatch& patch : prepared_children) {
    child_level_order.emplace_back(patch.descriptor().morton_key, patch.descriptor().patch_id);
  }
  std::sort(child_level_order.begin(), child_level_order.end());

  std::unordered_map<std::uint64_t, std::pair<std::size_t, std::size_t>> candidate_index;
  candidate_index.reserve(core::checkedSizeAdd(patchCount(), prepared_children.size(), "AMR patch-index reserve"));
  for (std::size_t level = 0; level < m_levels.size(); ++level) {
    if (level == child_level) {
      continue;
    }
    for (std::size_t offset = 0; offset < m_levels[level].size(); ++offset) {
      candidate_index.emplace(m_levels[level][offset].descriptor().patch_id, std::pair{level, offset});
    }
  }
  for (std::size_t offset = 0; offset < child_level_order.size(); ++offset) {
    candidate_index.emplace(child_level_order[offset].second, std::pair{child_level, offset});
  }

  if (child_level == m_levels.size()) {
    std::sort(prepared_children.begin(), prepared_children.end(), [](const AmrPatch& lhs, const AmrPatch& rhs) {
      return lhs.descriptor().morton_key < rhs.descriptor().morton_key;
    });
    m_levels.emplace_back(std::move(prepared_children));
  } else {
    auto& target_level = m_levels[child_level];
    for (AmrPatch& child : prepared_children) {
      target_level.emplace_back(std::move(child));
    }
    std::sort(target_level.begin(), target_level.end(), [](const AmrPatch& lhs, const AmrPatch& rhs) {
      return lhs.descriptor().morton_key < rhs.descriptor().morton_key;
    });
  }
  m_patch_index.swap(candidate_index);
  AmrPatch* committed_parent = findPatch(patch_id);
  if (committed_parent == nullptr) {
    throw std::logic_error("AMR refinement lost parent patch during transaction commit");
  }
  committed_parent->setLeaf(false);
  m_next_patch_id = candidate_next_patch_id;
  m_next_gas_cell_id = candidate_next_gas_cell_id;
  return child_ids;
}

bool PatchHierarchy::derefinePatch(std::uint64_t parent_patch_id) {
  AmrPatch* parent_patch = findPatch(parent_patch_id);
  if (parent_patch == nullptr) {
    return false;
  }
  const PatchDescriptor parent_descriptor = parent_patch->descriptor();
  const std::size_t child_level = static_cast<std::size_t>(parent_descriptor.level) + 1U;
  if (child_level >= m_levels.size()) {
    return false;
  }

  auto& level_patches = m_levels[child_level];
  std::array<const AmrPatch*, 8> child_by_octant{};
  std::size_t child_count = 0U;
  for (const AmrPatch& patch : level_patches) {
    if (patch.descriptor().parent_patch_id != parent_patch_id) {
      continue;
    }
    if (patch.descriptor().cell_dims != parent_descriptor.cell_dims) {
      throw std::runtime_error("Cannot conservatively derefine child patches with mismatched cell dimensions.");
    }
    const std::uint8_t octant = static_cast<std::uint8_t>(patch.descriptor().morton_key & 7U);
    if (child_by_octant[octant] != nullptr) {
      throw std::runtime_error("Cannot conservatively derefine duplicate child octant.");
    }
    child_by_octant[octant] = &patch;
    ++child_count;
  }
  if (child_count == 0U) {
    return false;
  }
  if (child_count != 8U) {
    throw std::runtime_error("Cannot conservatively derefine an incomplete octet of child patches.");
  }

  const std::size_t parent_cell_count = parent_patch->cellCount();
  std::vector<ConservedState> prepared_parent_conserved(parent_cell_count);
  std::vector<CellMetrics> prepared_parent_metrics(parent_cell_count);
  std::vector<std::uint64_t> prepared_retired_ids;
  prepared_retired_ids.reserve(core::checkedSizeMultiply(
      parent_cell_count, 8U, "AMR derefine retired gas-cell reserve"));
  const auto dims = parent_descriptor.cell_dims;
  const double parent_cell_volume_comov = parent_patch->cellVolumeComov();
  for (std::uint16_t parent_z = 0; parent_z < dims[2]; ++parent_z) {
    for (std::uint16_t parent_y = 0; parent_y < dims[1]; ++parent_y) {
      for (std::uint16_t parent_x = 0; parent_x < dims[0]; ++parent_x) {
        const std::size_t parent_index = cellIndex(dims, parent_x, parent_y, parent_z);
        std::array<ConservedState, 8> fine_conserved{};
        std::array<CellMetrics, 8> fine_metrics{};
        for (std::uint8_t fine_ordinal = 0; fine_ordinal < 8U; ++fine_ordinal) {
          const std::uint16_t fine_x_global = static_cast<std::uint16_t>(2U * parent_x + (fine_ordinal & 1U));
          const std::uint16_t fine_y_global = static_cast<std::uint16_t>(2U * parent_y + ((fine_ordinal >> 1U) & 1U));
          const std::uint16_t fine_z_global = static_cast<std::uint16_t>(2U * parent_z + ((fine_ordinal >> 2U) & 1U));
          const std::uint8_t octant = static_cast<std::uint8_t>(
              (fine_x_global >= dims[0] ? 1U : 0U) |
              (fine_y_global >= dims[1] ? 2U : 0U) |
              (fine_z_global >= dims[2] ? 4U : 0U));
          const AmrPatch* child_patch = child_by_octant[octant];
          if (child_patch == nullptr) {
            throw std::runtime_error("Cannot conservatively derefine: child octant coverage is incomplete.");
          }
          const std::uint16_t child_x = static_cast<std::uint16_t>(fine_x_global - ((octant & 1U) != 0U ? dims[0] : 0U));
          const std::uint16_t child_y = static_cast<std::uint16_t>(fine_y_global - ((octant & 2U) != 0U ? dims[1] : 0U));
          const std::uint16_t child_z = static_cast<std::uint16_t>(fine_z_global - ((octant & 4U) != 0U ? dims[2] : 0U));
          const std::size_t child_index = cellIndex(dims, child_x, child_y, child_z);
          fine_conserved[fine_ordinal] = child_patch->conservedView()[child_index];
          fine_metrics[fine_ordinal] = child_patch->metricsView()[child_index];
          prepared_retired_ids.push_back(child_patch->gasCellIdView()[child_index]);
        }
        prepared_parent_conserved[parent_index] = ConservativeTransfer::restrictToCoarse(fine_conserved);
        prepared_parent_metrics[parent_index] = restrictedParentMetrics(
            fine_metrics, prepared_parent_conserved[parent_index], parent_cell_volume_comov);
      }
    }
  }

  std::unordered_map<std::uint64_t, std::pair<std::size_t, std::size_t>> candidate_index;
  candidate_index.reserve(patchCount() - child_count);
  for (std::size_t level = 0; level < m_levels.size(); ++level) {
    std::size_t final_offset = 0U;
    for (const AmrPatch& patch : m_levels[level]) {
      if (level == child_level && patch.descriptor().parent_patch_id == parent_patch_id) {
        continue;
      }
      candidate_index.emplace(patch.descriptor().patch_id, std::pair{level, final_offset++});
    }
  }

  auto parent_conserved = parent_patch->conservedView();
  auto parent_metrics = parent_patch->metricsView();
  std::copy(prepared_parent_conserved.begin(), prepared_parent_conserved.end(), parent_conserved.begin());
  std::copy(prepared_parent_metrics.begin(), prepared_parent_metrics.end(), parent_metrics.begin());
  // Retired gas-cell identities are transaction-local invalidation evidence, not
  // an append-only history. Keeping every prior derefinement would make this
  // scaffold hierarchy grow monotonically across otherwise capacity-stable
  // regrid cycles. Swap in only the identities retired by this commit; the
  // previous transaction's storage is released when prepared_retired_ids dies.
  m_retired_gas_cell_ids.swap(prepared_retired_ids);
  level_patches.erase(
      std::remove_if(
          level_patches.begin(), level_patches.end(),
          [parent_patch_id](const AmrPatch& patch) {
            return patch.descriptor().parent_patch_id == parent_patch_id;
          }),
      level_patches.end());
  parent_patch->setLeaf(true);
  m_patch_index.swap(candidate_index);
  while (m_levels.size() > 1U && m_levels.back().empty()) {
    m_levels.pop_back();
  }
  return true;
}

AmrPatch* PatchHierarchy::findPatch(std::uint64_t patch_id) {
  const auto it = m_patch_index.find(patch_id);
  if (it == m_patch_index.end()) {
    return nullptr;
  }
  const auto [level_index, offset] = it->second;
  return &m_levels[level_index][offset];
}

const AmrPatch* PatchHierarchy::findPatch(std::uint64_t patch_id) const {
  const auto it = m_patch_index.find(patch_id);
  if (it == m_patch_index.end()) {
    return nullptr;
  }
  const auto [level_index, offset] = it->second;
  return &m_levels[level_index][offset];
}

std::span<AmrPatch> PatchHierarchy::levelView(std::size_t level) {
  if (level >= m_levels.size()) {
    return {};
  }
  return std::span<AmrPatch>(m_levels[level].data(), m_levels[level].size());
}

std::span<const AmrPatch> PatchHierarchy::levelView(std::size_t level) const {
  if (level >= m_levels.size()) {
    return {};
  }
  return std::span<const AmrPatch>(m_levels[level].data(), m_levels[level].size());
}

std::size_t PatchHierarchy::levelCount() const {
  return m_levels.size();
}

std::size_t PatchHierarchy::patchCount() const {
  std::size_t count = 0U;
  for (const auto& level : m_levels) {
    count = core::checkedSizeAdd(count, level.size(), "AMR hierarchy patch count");
  }
  return count;
}

std::size_t PatchHierarchy::allocatedCellCount() const {
  std::size_t count = 0U;
  for (const auto& level : m_levels) {
    for (const AmrPatch& patch : level) {
      count = core::checkedSizeAdd(count, patch.cellCount(), "AMR hierarchy allocated cell count");
    }
  }
  return count;
}

std::size_t PatchHierarchy::patchStorageCapacityBytes() const {
  std::size_t bytes = 0U;
  for (const auto& level : m_levels) {
    for (const AmrPatch& patch : level) {
      bytes = core::checkedSizeAdd(bytes, patch.ownedStorageCapacityBytes(), "AMR patch storage capacity total");
    }
  }
  return bytes;
}

std::span<const std::uint64_t> PatchHierarchy::retiredGasCellIds() const {
  return std::span<const std::uint64_t>(m_retired_gas_cell_ids.data(), m_retired_gas_cell_ids.size());
}

void PatchHierarchy::assignStableGasCellIds(AmrPatch& patch, std::uint64_t& next_gas_cell_id) {
  for (auto& gas_cell_id : patch.gasCellIdView()) {
    if (next_gas_cell_id == std::numeric_limits<std::uint64_t>::max()) {
      throw std::overflow_error("AMR gas-cell id space exhausted");
    }
    gas_cell_id = next_gas_cell_id++;
  }
}

void PatchHierarchy::rebuildPatchIndex() {
  m_patch_index.clear();
  for (std::size_t level = 0; level < m_levels.size(); ++level) {
    auto& patches = m_levels[level];
    std::sort(patches.begin(), patches.end(), [](const AmrPatch& lhs, const AmrPatch& rhs) {
      return lhs.descriptor().morton_key < rhs.descriptor().morton_key;
    });
    for (std::size_t i = 0; i < patches.size(); ++i) {
      m_patch_index[patches[i].descriptor().patch_id] = {level, i};
    }
  }
}

std::vector<ConservedState> ConservativeTransfer::prolongateFromCoarse(
    const ConservedState& coarse,
    std::size_t fine_cell_count) {
  if (fine_cell_count == 0) {
    return {};
  }
  std::vector<ConservedState> fine_values(fine_cell_count, coarse);
  const double scale = 1.0 / static_cast<double>(fine_cell_count);
  for (auto& fine : fine_values) {
    fine *= scale;
  }
  return fine_values;
}

ConservedState ConservativeTransfer::restrictToCoarse(std::span<const ConservedState> fine_cells) {
  ConservedState coarse;
  for (const auto& fine : fine_cells) {
    coarse += fine;
  }
  return coarse;
}

RefluxDiagnostics RefluxSynchronizer::apply(
    PatchHierarchy& hierarchy,
    std::span<const FluxRegisterEntry> entries) {
  RefluxDiagnostics diagnostics;

  for (const auto& entry : entries) {
    AmrPatch* coarse_patch = hierarchy.findPatch(entry.coarse_patch_id);
    if (coarse_patch == nullptr) {
      continue;
    }

    auto conserved = coarse_patch->conservedView();
    if (entry.coarse_cell_index >= conserved.size()) {
      continue;
    }

    const double inv_volume = 1.0 / coarse_patch->cellVolumeComov();
    ConservedState delta_flux = (entry.fine_face_flux_code - entry.coarse_face_flux_code);
    delta_flux *= (entry.face_area_comov * entry.dt_code * inv_volume);

    conserved[entry.coarse_cell_index] -= delta_flux;

    diagnostics.corrected_cells += 1;
    diagnostics.corrected_mass_code += std::abs(delta_flux.mass_code);
    diagnostics.corrected_momentum_x_code += std::abs(delta_flux.momentum_x_code);
    diagnostics.corrected_momentum_y_code += std::abs(delta_flux.momentum_y_code);
    diagnostics.corrected_momentum_z_code += std::abs(delta_flux.momentum_z_code);
    diagnostics.corrected_total_energy_code += std::abs(delta_flux.total_energy_code);
    diagnostics.corrected_metal_mass_code += std::abs(delta_flux.metal_mass_code);
    diagnostics.corrected_energy_code += std::abs(delta_flux.total_energy_code);
    diagnostics.corrected_internal_energy_code += std::abs(delta_flux.total_energy_code);
  }

  return diagnostics;
}

RefineDerefineManager::RefineDerefineManager(RefinementCriteria criteria)
    : m_criteria(std::move(criteria)) {}

RegridDiagnostics RefineDerefineManager::regrid(PatchHierarchy& hierarchy) const {
  const auto start_time = std::chrono::steady_clock::now();
  RegridDiagnostics diagnostics;

  std::vector<std::uint64_t> refine_ids;
  std::vector<std::uint64_t> derefine_parent_ids;

  for (std::size_t level = 0; level < hierarchy.levelCount(); ++level) {
    for (const auto& patch : hierarchy.levelView(level)) {
      if (!patch.isLeaf()) {
        continue;
      }
      diagnostics.touched_leaf_patch_count += 1;

      const RefinementDecision decision = evaluatePatchDecision(patch);
      if (decision == RefinementDecision::kRefine) {
        refine_ids.push_back(patch.descriptor().patch_id);
      } else if (decision == RefinementDecision::kDerefine && patch.descriptor().parent_patch_id != 0) {
        derefine_parent_ids.push_back(patch.descriptor().parent_patch_id);
      }
    }
  }

  std::sort(refine_ids.begin(), refine_ids.end());
  refine_ids.erase(std::unique(refine_ids.begin(), refine_ids.end()), refine_ids.end());
  for (const auto patch_id : refine_ids) {
    [[maybe_unused]] const auto child_ids = hierarchy.refinePatch(patch_id);
    diagnostics.refined_patch_count += 1;
  }

  std::sort(derefine_parent_ids.begin(), derefine_parent_ids.end());
  derefine_parent_ids.erase(
      std::unique(derefine_parent_ids.begin(), derefine_parent_ids.end()),
      derefine_parent_ids.end());
  for (const auto parent_id : derefine_parent_ids) {
    if (hierarchy.derefinePatch(parent_id)) {
      diagnostics.derefined_patch_count += 1;
    }
  }

  const auto end_time = std::chrono::steady_clock::now();
  diagnostics.elapsed_microseconds =
      static_cast<std::uint64_t>(
          std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count());

  return diagnostics;
}

RefinementDecision RefineDerefineManager::evaluatePatchDecision(const AmrPatch& patch) const {
  bool refine_any = false;
  bool derefine_all = true;

  for (const auto& metrics : patch.metricsView()) {
    const auto cell_decision =
        RefinementEvaluator::evaluateCell(metrics, patch.cellWidthComov(), m_criteria);
    if (cell_decision == RefinementDecision::kRefine) {
      refine_any = true;
      derefine_all = false;
      break;
    }
    if (cell_decision == RefinementDecision::kKeep) {
      derefine_all = false;
    }
  }

  if (refine_any) {
    return RefinementDecision::kRefine;
  }
  if (derefine_all) {
    return RefinementDecision::kDerefine;
  }
  return RefinementDecision::kKeep;
}

}  // namespace cosmosim::amr
