#include "validation_tolerance.hpp"

#include <cctype>
#include <charconv>
#include <fstream>
#include <stdexcept>

namespace cosmosim::validation {
namespace {

[[nodiscard]] std::string trim(std::string input) {
  while (!input.empty() && std::isspace(static_cast<unsigned char>(input.front())) != 0) {
    input.erase(input.begin());
  }
  while (!input.empty() && std::isspace(static_cast<unsigned char>(input.back())) != 0) {
    input.pop_back();
  }
  return input;
}

[[nodiscard]] double parseDouble(const std::string& value, const std::string& key) {
  if (value.empty()) {
    throw std::runtime_error("Validation tolerance key '" + key + "' has empty value");
  }
  double parsed = 0.0;
  const char* begin = value.data();
  const char* end = value.data() + value.size();
  auto result = std::from_chars(begin, end, parsed);
  if (result.ec != std::errc() || result.ptr != end) {
    throw std::runtime_error("Validation tolerance key '" + key + "' has non-numeric value '" + value + "'");
  }
  return parsed;
}

}  // namespace

ValidationToleranceTable ValidationToleranceTable::loadFromFile(const std::string& path) {
  std::ifstream input(path);
  if (!input.is_open()) {
    throw std::runtime_error("Failed to open validation tolerance file: " + path);
  }

  ValidationToleranceTable table;
  std::string line;
  std::size_t line_number = 0;
  while (std::getline(input, line)) {
    ++line_number;
    const std::string stripped = trim(line);
    if (stripped.empty() || stripped.front() == '#') {
      continue;
    }

    const std::size_t equals = stripped.find('=');
    if (equals == std::string::npos) {
      throw std::runtime_error(
          "Invalid validation tolerance entry at line " + std::to_string(line_number) + ": '" + stripped + "'");
    }

    const std::string key = trim(stripped.substr(0, equals));
    const std::string value = trim(stripped.substr(equals + 1));
    if (key.empty()) {
      throw std::runtime_error("Validation tolerance entry has empty key at line " + std::to_string(line_number));
    }

    table.m_values[key] = parseDouble(value, key);
  }

  return table;
}

bool ValidationToleranceTable::has(std::string_view key) const {
  return m_values.find(std::string(key)) != m_values.end();
}

double ValidationToleranceTable::require(std::string_view key) const {
  const auto it = m_values.find(std::string(key));
  if (it == m_values.end()) {
    throw std::runtime_error("Missing required validation tolerance key: " + std::string(key));
  }
  return it->second;
}

}  // namespace cosmosim::validation
