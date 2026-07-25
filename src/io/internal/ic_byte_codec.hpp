#pragma once

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace cosmosim::io::internal {

void appendLe32(std::vector<std::uint8_t>& output, std::uint32_t value);
void appendLe64(std::vector<std::uint8_t>& output, std::uint64_t value);
void appendDouble(std::vector<std::uint8_t>& output, double value);

[[nodiscard]] std::uint32_t readLe32(
    std::span<const std::uint8_t> bytes,
    std::size_t& offset);
[[nodiscard]] std::uint64_t readLe64(
    std::span<const std::uint8_t> bytes,
    std::size_t& offset);
[[nodiscard]] double readDouble(
    std::span<const std::uint8_t> bytes,
    std::size_t& offset);

}  // namespace cosmosim::io::internal
