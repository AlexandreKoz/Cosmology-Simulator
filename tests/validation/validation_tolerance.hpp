#pragma once

#include <string>
#include <string_view>
#include <unordered_map>

namespace cosmosim::validation {

class ValidationToleranceTable {
 public:
  [[nodiscard]] static ValidationToleranceTable loadFromFile(const std::string& path);

  [[nodiscard]] bool has(std::string_view key) const;
  [[nodiscard]] double require(std::string_view key) const;

 private:
  std::unordered_map<std::string, double> m_values;
};

}  // namespace cosmosim::validation
