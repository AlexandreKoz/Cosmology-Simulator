#include "cosmosim/io/ic_reader.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <string_view>

#include "core/internal/sha256.hpp"

namespace cosmosim::io {

std::string icSha256Hex(std::string_view value) {
  return core::internal::sha256Hex(value);
}

std::string icSha256FileHex(const std::filesystem::path& input_path) {
  std::ifstream input(input_path, std::ios::binary);
  if (!input) {
    throw std::runtime_error(
        "failed to open IC source for SHA-256: " + input_path.string());
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
    throw std::runtime_error(
        "failed while hashing IC source: " + input_path.string());
  }

  static constexpr char k_hex[] = "0123456789abcdef";
  const auto digest = hash.finish();
  std::string out(64U, '0');
  for (std::size_t i = 0; i < digest.size(); ++i) {
    out[i * 2U] = k_hex[digest[i] >> 4U];
    out[i * 2U + 1U] = k_hex[digest[i] & 0x0fU];
  }
  return out;
}

}  // namespace cosmosim::io
