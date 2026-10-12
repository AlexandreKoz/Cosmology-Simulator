#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <string_view>

#include "cosmosim/core/config.hpp"

namespace {
// Transform only the C++ loader's canonical serialization, then re-admit
// through that same typed loader. This is not a parameter-file parser.
void replaceNormalized(std::string& text, std::string_view section,
                       std::string_view key, const std::string& value) {
  const auto begin = text.find("["+std::string(section)+"]\n");
  if (begin == std::string::npos) throw std::runtime_error("normalized section missing");
  const auto end = text.find("\n[",begin);
  const auto position = text.find("\n"+std::string(key)+" = ",begin);
  if (position == std::string::npos || (end != std::string::npos && position >= end))
    throw std::runtime_error("normalized key missing: "+std::string(key));
  const auto value_begin = position+key.size()+4U;
  const auto value_end = text.find('\n',value_begin);
  text.replace(value_begin,value_end-value_begin,value);
}
}

int main(int argc, char** argv) {
  try {
    if (argc != 3) throw std::invalid_argument("usage: cosmosim_treepm_experiment_configs BASE.param.txt NEW_OUTPUT_DIRECTORY");
    const auto input = std::filesystem::absolute(argv[1]);
    const auto frozen = cosmosim::core::loadFrozenConfigFromFile(input);
    if (frozen.config.numerics.gravity_solver != cosmosim::core::GravitySolver::kTreePm)
      throw std::invalid_argument("base must select TreePM");
    if (frozen.config.compatibility.allow_unknown_keys)
      throw std::invalid_argument("experiment generation requires strict typed configuration");
    auto base = frozen.normalized_text;
    // Relocation must not silently point an experiment at a different input.
    const auto absolutePath = [&](std::string_view key, const std::string& value) {
      if (!value.empty() && value != "generated")
        replaceNormalized(base,"mode",key,(input.parent_path()/value).lexically_normal().string());
    };
    absolutePath("ic_file",frozen.config.mode.ic_file);
    absolutePath("ic_manifest_file",frozen.config.mode.ic_manifest_file);
    absolutePath("zoom_region_file",frozen.config.mode.zoom_region_file);
    replaceNormalized(base,"numerics","hierarchical_max_rung","0");
    replaceNormalized(base,"numerics","treepm_residual_block_size","64");
    replaceNormalized(base,"numerics","treepm_adaptive_maximum_opening_angle","0.25");
    for (const auto key : {"treepm_adaptive_acceptance_enabled", "treepm_gaussian_pair_lookup_enabled",
         "treepm_full_residual_lookup_enabled", "treepm_leaf_simd_enabled", "treepm_radial_coefficients_enabled",
         "treepm_identical_source_tree_reuse_enabled", "treepm_topology_refit_enabled"})
      replaceNormalized(base,"numerics",key,"false");
    replaceNormalized(base,"numerics","treepm_full_mac_diagnostics","true");
    replaceNormalized(base,"numerics","treepm_tree_relative_force_tolerance","0.005");
    const std::map<std::string,std::map<std::string,std::string>> experiments{
      {"reference",{}},
      {"optimized_matched",{{"treepm_full_mac_diagnostics","false"}}},
      {"matched_gaussian",{{"treepm_gaussian_pair_lookup_enabled","true"}}},
      {"matched_full_residual",{{"treepm_full_residual_lookup_enabled","true"}}},
      {"matched_simd",{{"treepm_leaf_simd_enabled","true"}}},
      {"matched_radial",{{"treepm_radial_coefficients_enabled","true"}}},
      {"matched_reuse",{{"treepm_identical_source_tree_reuse_enabled","true"}}},
      {"matched_refit",{{"treepm_topology_refit_enabled","true"}}},
      {"adaptive_025",{{"treepm_adaptive_acceptance_enabled","true"}}},
      {"adaptive_050_experimental",{{"treepm_adaptive_acceptance_enabled","true"},
          {"treepm_adaptive_maximum_opening_angle","0.5"}}},
    };
    // Validate every candidate before creating any output.
    std::map<std::string,cosmosim::core::FrozenConfig> candidates;
    for (const auto& [name,changes] : experiments) {
      auto text=base;
      for (const auto& [key,value] : changes) replaceNormalized(text,"numerics",key,value);
      replaceNormalized(text,"output","run_name","treepm_"+name);
      candidates.emplace(name,cosmosim::core::loadFrozenConfigFromString(text,"treepm_"+name));
    }
    const std::filesystem::path output(argv[2]);
    if (!std::filesystem::create_directory(output)) throw std::runtime_error("output directory must be new");
    std::ofstream manifest(output/"experiments.txt");
    manifest << "CHUI TreePM experiments: SOURCE IMPLEMENTED / QUALIFICATION PENDING\n"
             << "base=" << input.string() << "\nbase_config_hash=" << frozen.provenance.config_hash_hex << '\n';
    for (const auto& [name,candidate] : candidates) {
      std::ofstream file(output/(name+".param.txt"));
      file << candidate.normalized_text;
      file.close();
      if (!file) throw std::runtime_error("failed writing experiment "+name);
      manifest << name << " config_hash=" << candidate.provenance.config_hash_hex << '\n';
    }
    manifest.close();
    if (!manifest) throw std::runtime_error("failed writing experiment manifest");
    std::cout << "Wrote " << candidates.size() << " typed, normalized experiments to " << output << '\n';
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n'; return 1;
  }
}
