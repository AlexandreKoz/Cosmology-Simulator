#include "cosmosim/analysis/halo_workflow.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <utility>

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
    const core::SimulationState& state,
    std::uint32_t i,
    std::uint32_t j,
    double box_size_comov) {
  const double dx = wrapPeriodicDelta(
      state.particles.position_x_comoving[i] - state.particles.position_x_comoving[j],
      box_size_comov);
  const double dy = wrapPeriodicDelta(
      state.particles.position_y_comoving[i] - state.particles.position_y_comoving[j],
      box_size_comov);
  const double dz = wrapPeriodicDelta(
      state.particles.position_z_comoving[i] - state.particles.position_z_comoving[j],
      box_size_comov);
  return dx * dx + dy * dy + dz * dz;
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

}  // namespace

std::uint32_t HaloCatalogSchema::schemaVersion() noexcept { return 1; }

std::string_view HaloCatalogSchema::catalogFormatName() noexcept { return "cosmosim_halo_catalog_v1"; }

std::string_view HaloCatalogSchema::mergerTreePlanFormatName() noexcept {
  return "cosmosim_merger_tree_plan_v1";
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

HaloCatalog FofHaloFinder::buildCatalog(
    const core::SimulationState& state,
    const core::SimulationConfig& config,
    std::uint64_t snapshot_step_index,
    double snapshot_scale_factor,
    FofProfilingCounters* profiling) const {
  HaloCatalog catalog;
  catalog.schema_version = HaloCatalogSchema::schemaVersion();
  catalog.snapshot_step_index = snapshot_step_index;
  catalog.snapshot_scale_factor = snapshot_scale_factor;
  catalog.run_name = config.output.run_name;
  catalog.normalized_config_hash = state.metadata.normalized_config_hash;
  catalog.halo_id_by_particle.assign(state.particles.size(), k_unbound_halo_id);

  std::vector<std::uint32_t> candidate_indices;
  candidate_indices.reserve(state.particles.size());
  for (std::uint32_t i = 0; i < state.particles.size(); ++i) {
    const auto species = static_cast<core::ParticleSpecies>(state.particle_sidecar.species_tag[i]);
    if (includeSpecies(species)) {
      candidate_indices.push_back(i);
    }
  }

  if (candidate_indices.empty()) {
    if (profiling != nullptr) {
      profiling->candidate_particle_count = 0;
      profiling->pair_checks = 0;
      profiling->pair_links = 0;
      profiling->linking_length_comov = 0.0;
    }
    return catalog;
  }

  const double box_size_comov = config.cosmology.box_size_mpc_comoving;
  const double mean_spacing_comov = box_size_comov / std::cbrt(static_cast<double>(candidate_indices.size()));
  const double linking_length_comov = m_config.linking_length_factor_mean_interparticle * mean_spacing_comov;
  const double linking_length_sq = linking_length_comov * linking_length_comov;

  DisjointSet ds(candidate_indices.size());
  std::uint64_t pair_checks = 0;
  std::uint64_t pair_links = 0;

  for (std::uint32_t a = 0; a < candidate_indices.size(); ++a) {
    for (std::uint32_t b = a + 1; b < candidate_indices.size(); ++b) {
      ++pair_checks;
      const double r2 = periodicDistanceSquared(state, candidate_indices[a], candidate_indices[b], box_size_comov);
      if (r2 <= linking_length_sq && ds.unite(a, b)) {
        ++pair_links;
      }
    }
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

  std::sort(accepted_groups.begin(), accepted_groups.end(), [&state](const auto& left, const auto& right) {
    const std::uint64_t left_min = *std::min_element(
        left.begin(),
        left.end(),
        [&state](std::uint32_t a, std::uint32_t b) { return state.particle_sidecar.particle_id[a] < state.particle_sidecar.particle_id[b]; });
    const std::uint64_t right_min = *std::min_element(
        right.begin(),
        right.end(),
        [&state](std::uint32_t a, std::uint32_t b) { return state.particle_sidecar.particle_id[a] < state.particle_sidecar.particle_id[b]; });
    return state.particle_sidecar.particle_id[left_min] < state.particle_sidecar.particle_id[right_min];
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

    for (const std::uint32_t p : members) {
      catalog.halo_id_by_particle[p] = halo_id;
      const double mass = state.particles.mass_code[p];
      entry.total_mass_code += mass;
      entry.center_of_mass_comov[0] += mass * state.particles.position_x_comoving[p];
      entry.center_of_mass_comov[1] += mass * state.particles.position_y_comoving[p];
      entry.center_of_mass_comov[2] += mass * state.particles.position_z_comoving[p];
      entry.bulk_velocity_peculiar[0] += mass * state.particles.velocity_x_peculiar[p];
      entry.bulk_velocity_peculiar[1] += mass * state.particles.velocity_y_peculiar[p];
      entry.bulk_velocity_peculiar[2] += mass * state.particles.velocity_z_peculiar[p];
      entry.min_particle_id = std::min(entry.min_particle_id, state.particle_sidecar.particle_id[p]);
    }

    if (entry.total_mass_code > 0.0) {
      for (std::size_t axis = 0; axis < 3; ++axis) {
        entry.center_of_mass_comov[axis] /= entry.total_mass_code;
        entry.bulk_velocity_peculiar[axis] /= entry.total_mass_code;
      }
    }

    catalog.halos.push_back(entry);

    // Subhalo pipeline is not solved in v1; keep an explicit placeholder rank for future bound finders.
    catalog.subhalo_candidates.push_back(SubhaloCandidateEntry{
        .subhalo_id = halo_id,
        .host_halo_id = halo_id,
        .rank_in_host = 0,
        .particle_count = entry.particle_count,
        .bound_mass_code = entry.total_mass_code,
    });
  }

  if (profiling != nullptr) {
    profiling->candidate_particle_count = candidate_indices.size();
    profiling->pair_checks = pair_checks;
    profiling->pair_links = pair_links;
    profiling->linking_length_comov = linking_length_comov;
  }

  return catalog;
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
