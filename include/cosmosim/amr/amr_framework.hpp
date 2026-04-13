#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <unordered_map>
#include <vector>

namespace cosmosim::amr {

struct ConservedState {
  double mass_code = 0.0;
  double momentum_x_code = 0.0;
  double momentum_y_code = 0.0;
  double momentum_z_code = 0.0;
  double total_energy_code = 0.0;

  ConservedState& operator+=(const ConservedState& rhs);
  ConservedState& operator-=(const ConservedState& rhs);
  ConservedState& operator*=(double factor);
};

[[nodiscard]] ConservedState operator+(ConservedState lhs, const ConservedState& rhs);
[[nodiscard]] ConservedState operator-(ConservedState lhs, const ConservedState& rhs);
[[nodiscard]] ConservedState operator*(ConservedState lhs, double factor);

struct CellMetrics {
  double density_code = 0.0;
  double pressure_code = 0.0;
  double sound_speed_code = 0.0;
  double gradient_indicator = 0.0;
  std::uint32_t particle_count = 0;
};

enum class RefinementDecision : std::uint8_t {
  kKeep,
  kRefine,
  kDerefine,
};

struct RefinementCriteria {
  bool use_mass_threshold = true;
  bool use_gradient_indicator = true;
  bool use_particle_count = true;
  bool use_jeans_condition = true;

  double mass_threshold_code = 1.0;
  double gradient_threshold = 0.5;
  std::uint32_t particle_threshold = 8;
  double jeans_resolution_cells = 4.0;
  double gravitational_constant_code = 1.0;
  double derefine_hysteresis = 0.5;
};

class RefinementEvaluator {
 public:
  [[nodiscard]] static RefinementDecision evaluateCell(
      const CellMetrics& metrics,
      double cell_width_comov,
      const RefinementCriteria& criteria);
};

struct PatchDescriptor {
  std::uint64_t patch_id = 0;
  std::uint64_t parent_patch_id = 0;
  std::uint8_t level = 0;
  std::uint64_t morton_key = 0;
  std::array<double, 3> origin_comov = {0.0, 0.0, 0.0};
  std::array<double, 3> extent_comov = {1.0, 1.0, 1.0};
  std::array<std::uint16_t, 3> cell_dims = {4, 4, 4};
};

class AmrPatch {
 public:
  explicit AmrPatch(PatchDescriptor descriptor);

  [[nodiscard]] const PatchDescriptor& descriptor() const;
  [[nodiscard]] std::size_t cellCount() const;
  [[nodiscard]] double cellVolumeComov() const;
  [[nodiscard]] double cellWidthComov() const;

  [[nodiscard]] std::span<ConservedState> conservedView();
  [[nodiscard]] std::span<const ConservedState> conservedView() const;
  [[nodiscard]] std::span<CellMetrics> metricsView();
  [[nodiscard]] std::span<const CellMetrics> metricsView() const;

  [[nodiscard]] ConservedState totalConserved() const;
  [[nodiscard]] bool isLeaf() const;
  void setLeaf(bool is_leaf);

 private:
  PatchDescriptor m_descriptor;
  std::vector<ConservedState> m_conserved;
  std::vector<CellMetrics> m_metrics;
  bool m_is_leaf = true;
};

class PatchHierarchy {
 public:
  [[nodiscard]] std::uint64_t createRootPatch(const PatchDescriptor& root);
  [[nodiscard]] std::array<std::uint64_t, 8> refinePatch(std::uint64_t patch_id);
  bool derefinePatch(std::uint64_t parent_patch_id);

  [[nodiscard]] AmrPatch* findPatch(std::uint64_t patch_id);
  [[nodiscard]] const AmrPatch* findPatch(std::uint64_t patch_id) const;

  [[nodiscard]] std::span<AmrPatch> levelView(std::size_t level);
  [[nodiscard]] std::span<const AmrPatch> levelView(std::size_t level) const;

  [[nodiscard]] std::size_t levelCount() const;
  [[nodiscard]] std::size_t patchCount() const;

 private:
  std::vector<std::vector<AmrPatch>> m_levels;
  std::unordered_map<std::uint64_t, std::pair<std::size_t, std::size_t>> m_patch_index;
  std::uint64_t m_next_patch_id = 1;

  void rebuildPatchIndex();
};

class ConservativeTransfer {
 public:
  [[nodiscard]] static std::vector<ConservedState> prolongateFromCoarse(
      const ConservedState& coarse,
      std::size_t fine_cell_count);

  [[nodiscard]] static ConservedState restrictToCoarse(std::span<const ConservedState> fine_cells);
};

struct FluxRegisterEntry {
  std::uint64_t coarse_patch_id = 0;
  std::size_t coarse_cell_index = 0;
  ConservedState coarse_face_flux_code;
  ConservedState fine_face_flux_code;
  double face_area_comov = 1.0;
  double dt_code = 0.0;
};

struct RefluxDiagnostics {
  std::size_t corrected_cells = 0;
  double corrected_mass_code = 0.0;
  double corrected_energy_code = 0.0;
};

class RefluxSynchronizer {
 public:
  [[nodiscard]] static RefluxDiagnostics apply(
      PatchHierarchy& hierarchy,
      std::span<const FluxRegisterEntry> entries);
};

struct RegridDiagnostics {
  std::size_t refined_patch_count = 0;
  std::size_t derefined_patch_count = 0;
  std::size_t touched_leaf_patch_count = 0;
  std::uint64_t elapsed_microseconds = 0;
};

class RefineDerefineManager {
 public:
  explicit RefineDerefineManager(RefinementCriteria criteria);

  [[nodiscard]] RegridDiagnostics regrid(PatchHierarchy& hierarchy) const;

 private:
  RefinementCriteria m_criteria;

  [[nodiscard]] RefinementDecision evaluatePatchDecision(const AmrPatch& patch) const;
};

}  // namespace cosmosim::amr
