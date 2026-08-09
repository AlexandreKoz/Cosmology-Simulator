#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>

#include "core/internal/sha256.hpp"

namespace cosmosim::io::internal {

[[nodiscard]] inline std::string sha256FileHex(const std::filesystem::path& path) {
  std::ifstream input(path, std::ios::binary);
  if (!input) {
    throw std::runtime_error("failed to open file for SHA-256: " + path.string());
  }
  core::internal::Sha256 hash;
  std::array<std::uint8_t, 1U << 16U> buffer{};
  while (input) {
    input.read(
        reinterpret_cast<char*>(buffer.data()),
        static_cast<std::streamsize>(buffer.size()));
    const std::streamsize count = input.gcount();
    if (count > 0) {
      hash.update(buffer.data(), static_cast<std::size_t>(count));
    }
  }
  if (!input.eof()) {
    throw std::runtime_error("failed while hashing file: " + path.string());
  }
  static constexpr char kHex[] = "0123456789abcdef";
  const auto digest = hash.finish();
  std::string out(64U, '0');
  for (std::size_t i = 0; i < digest.size(); ++i) {
    out[i * 2U] = kHex[digest[i] >> 4U];
    out[i * 2U + 1U] = kHex[digest[i] & 0x0fU];
  }
  return out;
}

}  // namespace cosmosim::io::internal
