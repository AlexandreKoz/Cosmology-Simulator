#include "io/internal/ic_byte_codec.hpp"

#include <bit>
#include <stdexcept>

namespace cosmosim::io::internal {

void appendLe32(std::vector<std::uint8_t>& output, std::uint32_t value) {
  for (unsigned shift = 0U; shift < 32U; shift += 8U) {
    output.push_back(
        static_cast<std::uint8_t>((value >> shift) & 0xffU));
  }
}

void appendLe64(std::vector<std::uint8_t>& output, std::uint64_t value) {
  for (unsigned shift = 0U; shift < 64U; shift += 8U) {
    output.push_back(
        static_cast<std::uint8_t>((value >> shift) & 0xffU));
  }
}

void appendDouble(std::vector<std::uint8_t>& output, double value) {
  appendLe64(output, std::bit_cast<std::uint64_t>(value));
}

std::uint32_t readLe32(
    std::span<const std::uint8_t> bytes,
    std::size_t& offset) {
  if (offset + 4U > bytes.size()) {
    throw std::runtime_error("truncated IC wire uint32");
  }
  std::uint32_t value = 0U;
  for (unsigned shift = 0U; shift < 32U; shift += 8U) {
    value |= static_cast<std::uint32_t>(bytes[offset++]) << shift;
  }
  return value;
}

std::uint64_t readLe64(
    std::span<const std::uint8_t> bytes,
    std::size_t& offset) {
  if (offset + 8U > bytes.size()) {
    throw std::runtime_error("truncated IC wire uint64");
  }
  std::uint64_t value = 0U;
  for (unsigned shift = 0U; shift < 64U; shift += 8U) {
    value |= static_cast<std::uint64_t>(bytes[offset++]) << shift;
  }
  return value;
}

double readDouble(
    std::span<const std::uint8_t> bytes,
    std::size_t& offset) {
  return std::bit_cast<double>(readLe64(bytes, offset));
}

}  // namespace cosmosim::io::internal
