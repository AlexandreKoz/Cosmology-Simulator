#pragma once

#include <source_location>
#include <stdexcept>
#include <string>

namespace cosmosim::test_support {

inline void requireOrThrow(
    bool condition, const char* expression,
    const std::source_location location = std::source_location::current()) {
  if (!condition) {
    throw std::runtime_error(std::string(location.file_name()) + ':' +
        std::to_string(location.line()) + ": requirement failed: " + expression);
  }
}

}  // namespace cosmosim::test_support

// Always evaluated, including manual builds with NDEBUG. Each expression is
// evaluated exactly once; scientific tolerances remain at their call sites.
#define CHUI_REQUIRE(...) \
  ::cosmosim::test_support::requireOrThrow(static_cast<bool>((__VA_ARGS__)), #__VA_ARGS__)
