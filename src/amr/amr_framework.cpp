#include "cosmosim/amr/amr_framework.hpp"

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
      .total_energy_code = flux.total_energy_density_comoving};
}

void validateCompatibleRegisterRecord(
    const FluxRegisterEntry& entry,
    const hydro::HydroFluxRegisterRecord& record) {
  if (entry.coarse_patch_id != record.coarse_patch_id ||
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
  return *this;
}

ConservedState& ConservedState::operator-=(const ConservedState& rhs) {
  mass_code -= rhs.mass_code;
  momentum_x_code -= rhs.momentum_x_code;
  momentum_y_code -= rhs.momentum_y_code;
  momentum_z_code -= rhs.momentum_z_code;
  total_energy_code -= rhs.total_energy_code;
  return *this;
}

ConservedState& ConservedState::operator*=(double factor) {
  mass_code *= factor;
  momentum_x_code *= factor;
  momentum_y_code *= factor;
  momentum_z_code *= factor;
  total_energy_code *= factor;
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
    : m_descriptor(std::move(descriptor)),
      m_conserved(static_cast<std::size_t>(m_descriptor.cell_dims[0]) *
                  static_cast<std::size_t>(m_descriptor.cell_dims[1]) *
                  static_cast<std::size_t>(m_descriptor.cell_dims[2])),
      m_metrics(m_conserved.size()),
      m_gas_cell_ids(m_conserved.size()) {}

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

  PatchDescriptor root_copy = root;
  root_copy.patch_id = m_next_patch_id++;
  root_copy.parent_patch_id = 0;
  root_copy.level = 0;

  m_levels.emplace_back();
  m_levels[0].emplace_back(root_copy);
  assignStableGasCellIds(m_levels[0].back());
  rebuildPatchIndex();
  return root_copy.patch_id;
}

std::array<std::uint64_t, 8> PatchHierarchy::refinePatch(std::uint64_t patch_id) {
  AmrPatch* parent_patch = findPatch(patch_id);
  if (parent_patch == nullptr) {
    throw std::runtime_error("Cannot refine missing patch id.");
  }
  if (!parent_patch->isLeaf()) {
    throw std::runtime_error("Cannot refine non-leaf patch id.");
  }

  const PatchDescriptor& parent = parent_patch->descriptor();
  const PatchDescriptor parent_descriptor = parent;
  const std::vector<ConservedState> parent_conserved(
      parent_patch->conservedView().begin(),
      parent_patch->conservedView().end());
  const std::vector<CellMetrics> parent_metrics(
      parent_patch->metricsView().begin(),
      parent_patch->metricsView().end());
  std::vector<std::vector<ConservedState>> parent_prolongated;
  parent_prolongated.reserve(parent_conserved.size());
  for (const auto& parent_cell : parent_conserved) {
    parent_prolongated.push_back(ConservativeTransfer::prolongateFromCoarse(parent_cell, 8));
  }

  const std::size_t child_level = static_cast<std::size_t>(parent.level) + 1;
  if (child_level >= m_levels.size()) {
    m_levels.resize(child_level + 1);
  }

  const std::array<double, 3> child_extent = {
      parent_descriptor.extent_comov[0] * 0.5,
      parent_descriptor.extent_comov[1] * 0.5,
      parent_descriptor.extent_comov[2] * 0.5,
  };

  std::array<std::uint64_t, 8> child_ids{};
  for (std::uint8_t octant = 0; octant < 8; ++octant) {
    PatchDescriptor child;
    child.patch_id = m_next_patch_id++;
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
    assignStableGasCellIds(child_patch);
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

          child_conserved[child_index] = parent_prolongated[parent_index][fine_ordinal];
          child_metrics[child_index] = prolongatedChildMetrics(parent_metrics[parent_index], fine_ordinal);
        }
      }
    }

    m_levels[child_level].emplace_back(std::move(child_patch));
    child_ids[octant] = child.patch_id;
  }

  parent_patch->setLeaf(false);
  rebuildPatchIndex();
  return child_ids;
}

bool PatchHierarchy::derefinePatch(std::uint64_t parent_patch_id) {
  AmrPatch* parent_patch = findPatch(parent_patch_id);
  if (parent_patch == nullptr) {
    return false;
  }

  const PatchDescriptor parent_descriptor = parent_patch->descriptor();
  const std::uint8_t child_level = parent_patch->descriptor().level + 1;
  if (static_cast<std::size_t>(child_level) >= m_levels.size()) {
    return false;
  }

  auto& level_patches = m_levels[child_level];
  std::vector<const AmrPatch*> child_patches;
  child_patches.reserve(8);
  for (const auto& patch : level_patches) {
    if (patch.descriptor().parent_patch_id == parent_patch_id) {
      child_patches.push_back(&patch);
    }
  }
  if (child_patches.empty()) {
    return false;
  }
  if (child_patches.size() != 8) {
    throw std::runtime_error("Cannot conservatively derefine an incomplete octet of child patches.");
  }

  std::vector<std::vector<ConservedState>> fine_conserved_by_parent(parent_patch->cellCount());
  std::vector<std::vector<CellMetrics>> fine_metrics_by_parent(parent_patch->cellCount());
  for (auto& fine_cells : fine_conserved_by_parent) {
    fine_cells.reserve(8);
  }
  for (auto& fine_metrics : fine_metrics_by_parent) {
    fine_metrics.reserve(8);
  }

  for (const auto* child_patch : child_patches) {
    const PatchDescriptor& child = child_patch->descriptor();
    if (child.cell_dims != parent_descriptor.cell_dims) {
      throw std::runtime_error("Cannot conservatively derefine child patches with mismatched cell dimensions.");
    }

    const std::uint8_t octant = static_cast<std::uint8_t>(child.morton_key & 7U);
    const auto child_conserved = child_patch->conservedView();
    const auto child_metrics = child_patch->metricsView();
    const auto child_gas_cell_ids = child_patch->gasCellIdView();
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

          fine_conserved_by_parent[parent_index].push_back(child_conserved[child_index]);
          fine_metrics_by_parent[parent_index].push_back(child_metrics[child_index]);
          m_retired_gas_cell_ids.push_back(child_gas_cell_ids[child_index]);
        }
      }
    }
  }

  auto parent_conserved = parent_patch->conservedView();
  auto parent_metrics = parent_patch->metricsView();
  for (std::size_t parent_index = 0; parent_index < parent_conserved.size(); ++parent_index) {
    if (fine_conserved_by_parent[parent_index].size() != 8) {
      throw std::runtime_error("Cannot conservatively derefine: fine child coverage is incomplete.");
    }

    parent_conserved[parent_index] =
        ConservativeTransfer::restrictToCoarse(fine_conserved_by_parent[parent_index]);
    parent_metrics[parent_index] = restrictedParentMetrics(
        fine_metrics_by_parent[parent_index],
        parent_conserved[parent_index],
        parent_patch->cellVolumeComov());
  }

  const auto old_size = level_patches.size();
  level_patches.erase(
      std::remove_if(
          level_patches.begin(),
          level_patches.end(),
          [parent_patch_id](const AmrPatch& patch) {
            return patch.descriptor().parent_patch_id == parent_patch_id;
          }),
      level_patches.end());

  const bool removed_children = level_patches.size() != old_size;
  if (removed_children) {
    parent_patch->setLeaf(true);
    rebuildPatchIndex();
  }
  return removed_children;
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
  std::size_t count = 0;
  for (const auto& level : m_levels) {
    count += level.size();
  }
  return count;
}

std::span<const std::uint64_t> PatchHierarchy::retiredGasCellIds() const {
  return std::span<const std::uint64_t>(m_retired_gas_cell_ids.data(), m_retired_gas_cell_ids.size());
}

void PatchHierarchy::assignStableGasCellIds(AmrPatch& patch) {
  for (auto& gas_cell_id : patch.gasCellIdView()) {
    gas_cell_id = m_next_gas_cell_id++;
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
