#pragma once

#include <filesystem>
#include <string_view>

namespace cosmosim::io::internal {

enum class FileDurability {
  kAtomicVisibility,
  kDurablePublication,
};

void ensureParentDirectory(const std::filesystem::path& output_path);
void validateTemporarySuffix(std::string_view suffix);
[[nodiscard]] std::filesystem::path makeUniqueTemporarySibling(
    const std::filesystem::path& output_path,
    std::string_view suffix = ".part");
void syncFileForDurability(const std::filesystem::path& path);
void syncDirectoryForDurability(const std::filesystem::path& directory);
void atomicReplaceFile(
    const std::filesystem::path& temporary_path,
    const std::filesystem::path& output_path,
    FileDurability durability);

class TransactionalFileTarget {
 public:
  TransactionalFileTarget(
      std::filesystem::path output_path,
      std::string_view temporary_suffix,
      FileDurability durability);
  TransactionalFileTarget(const TransactionalFileTarget&) = delete;
  TransactionalFileTarget& operator=(const TransactionalFileTarget&) = delete;
  TransactionalFileTarget(TransactionalFileTarget&&) = delete;
  TransactionalFileTarget& operator=(TransactionalFileTarget&&) = delete;
  ~TransactionalFileTarget();

  [[nodiscard]] const std::filesystem::path& outputPath() const noexcept;
  [[nodiscard]] const std::filesystem::path& temporaryPath() const noexcept;
  [[nodiscard]] FileDurability durability() const noexcept;
  void publish();

 private:
  std::filesystem::path m_output_path;
  std::filesystem::path m_temporary_path;
  FileDurability m_durability = FileDurability::kAtomicVisibility;
  bool m_published = false;
};

void writeTextFileTransactionally(
    const std::filesystem::path& output_path,
    std::string_view contents,
    FileDurability durability = FileDurability::kAtomicVisibility);

}  // namespace cosmosim::io::internal
