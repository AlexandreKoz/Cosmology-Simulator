#pragma once

#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>

namespace cosmosim::core {

[[nodiscard]] inline std::size_t checkedSizeAdd(
    std::size_t lhs,
    std::size_t rhs,
    std::string_view context) {
  if (rhs > std::numeric_limits<std::size_t>::max() - lhs) {
    throw std::overflow_error(std::string(context) + ": size addition overflows size_t");
  }
  return lhs + rhs;
}

[[nodiscard]] inline std::size_t checkedSizeMultiply(
    std::size_t lhs,
    std::size_t rhs,
    std::string_view context) {
  if (lhs != 0U && rhs > std::numeric_limits<std::size_t>::max() / lhs) {
    throw std::overflow_error(std::string(context) + ": size product overflows size_t");
  }
  return lhs * rhs;
}

[[nodiscard]] inline std::size_t checkedSizeProduct3(
    std::size_t a,
    std::size_t b,
    std::size_t c,
    std::string_view context) {
  return checkedSizeMultiply(checkedSizeMultiply(a, b, context), c, context);
}

[[nodiscard]] inline std::size_t checkedAlignUpSize(
    std::size_t value,
    std::size_t alignment,
    std::string_view context) {
  if (alignment == 0U || (alignment & (alignment - 1U)) != 0U) {
    throw std::invalid_argument(std::string(context) + ": alignment must be a nonzero power of two");
  }
  const std::size_t padded = checkedSizeAdd(value, alignment - 1U, context);
  return padded & ~(alignment - 1U);
}

template <typename TDestination, typename TSource>
[[nodiscard]] constexpr TDestination checkedIntegralNarrow(
    TSource value,
    std::string_view context) {
  static_assert(std::is_integral_v<TDestination> && std::is_integral_v<TSource>);
  if (!std::in_range<TDestination>(value)) {
    throw std::overflow_error(std::string(context) + ": integer value is outside destination range");
  }
  return static_cast<TDestination>(value);
}

}  // namespace cosmosim::core
