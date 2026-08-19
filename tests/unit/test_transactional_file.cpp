#include <cassert>
#include <filesystem>
#include <fstream>
#include <string>

#include "io/internal/transactional_file.hpp"
#include "../support/test_temp_workspace.hpp"

namespace {

std::string readText(const std::filesystem::path& path) {
  std::ifstream input(path, std::ios::binary);
  assert(input);
  return std::string(std::istreambuf_iterator<char>(input), std::istreambuf_iterator<char>());
}

}  // namespace

int main() {
  const auto root = cosmosim::test_support::TestTempWorkspace::uniqueProcessLocalPath("cosmosim_transactional_file_test");
  std::filesystem::remove_all(root);
  std::filesystem::create_directories(root);

  const auto output = root / "checkpoint.txt";
  {
    std::ofstream initial(output, std::ios::binary | std::ios::trunc);
    initial << "old";
  }
  cosmosim::io::internal::writeTextFileTransactionally(
      output, "new", cosmosim::io::internal::FileDurability::kDurablePublication);
  assert(readText(output) == "new");

  const auto previous = std::filesystem::current_path();
  std::filesystem::current_path(root);
  cosmosim::io::internal::writeTextFileTransactionally(
      "bare.txt", "bare", cosmosim::io::internal::FileDurability::kDurablePublication);
  assert(readText("bare.txt") == "bare");
  std::filesystem::current_path(previous);

  std::size_t regular_files = 0U;
  for (const auto& entry : std::filesystem::directory_iterator(root)) {
    if (entry.is_regular_file()) ++regular_files;
  }
  assert(regular_files == 2U);
  std::filesystem::remove_all(root);
  return 0;
}
