#include <chrono>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <string_view>
#include <stdexcept>
#include <vector>

namespace {

[[nodiscard]] std::string readTextFile(const std::filesystem::path& path) {
  std::ifstream input(path);
  if (!input) {
    throw std::runtime_error("unable to open file: " + path.string());
  }

  return std::string(
      std::istreambuf_iterator<char>(input),
      std::istreambuf_iterator<char>());
}

[[nodiscard]] std::size_t countMatches(
    const std::string& text,
    const std::vector<std::string_view>& patterns) {
  std::size_t count = 0;
  for (const auto pattern : patterns) {
    if (text.find(pattern) != std::string::npos) {
      ++count;
    }
  }
  return count;
}

}  // namespace

int main() {
  const std::filesystem::path repo_root = std::filesystem::current_path();
  const std::vector<std::filesystem::path> docs = {
      repo_root / "README.md",
      repo_root / "RELEASE_NOTES.md",
      repo_root / "CONTRIBUTING.md",
      repo_root / "docs/build_instructions.md",
      repo_root / "docs/configuration.md",
      repo_root / "docs/output_schema.md",
      repo_root / "docs/validation_plan.md",
      repo_root / "docs/profiling.md",
      repo_root / "docs/releases/initial_release_scope.md",
      repo_root / "docs/releases/initial_release_validation_summary.md",
      repo_root / "docs/releases/initial_release_benchmark_summary.md",
      repo_root / "docs/releases/known_issues.md",
      repo_root / "docs/releases/release_checklist.md",
      repo_root / "docs/architecture/overview.md",
      repo_root / "docs/architecture/developer_workflow_contract.md"};

  const std::vector<std::string_view> required_markers = {
      "docs/build_instructions.md",
      "docs/configuration.md",
      "docs/output_schema.md",
      "docs/validation_plan.md",
      "docs/profiling.md",
      "RELEASE_NOTES.md",
      "docs/releases/initial_release_scope.md",
      "docs/releases/initial_release_validation_summary.md",
      "docs/releases/initial_release_benchmark_summary.md",
      "docs/releases/known_issues.md",
      "docs/releases/release_checklist.md",
      "docs/architecture/overview.md",
      "docs/architecture/developer_workflow_contract.md",
      "gadget_arepo_v1",
      "cosmosim_restart_v2",
      "provenance_v1",
      "release_manifest_v1",
      "mode.mode"};

  std::size_t total_bytes = 0;
  std::size_t total_matches = 0;

  const auto start = std::chrono::steady_clock::now();
  for (const auto& path : docs) {
    const std::string content = readTextFile(path);
    total_bytes += content.size();
    total_matches += countMatches(content, required_markers);
  }
  const auto end = std::chrono::steady_clock::now();
  const auto elapsed_us =
      std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

  std::cout << "bench_docs_reference_scan"
            << " files=" << docs.size()
            << " bytes=" << total_bytes
            << " marker_hits=" << total_matches
            << " elapsed_us=" << elapsed_us << '\n';

  return 0;
}
