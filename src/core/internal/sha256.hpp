#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>

namespace cosmosim::core::internal {

class Sha256 {
 public:
  Sha256() = default;

  void update(const std::uint8_t* data, std::size_t size) {
    for (std::size_t i = 0; i < size; ++i) {
      m_block[m_block_size++] = data[i];
      m_bit_count += 8U;
      if (m_block_size == m_block.size()) {
        transform();
        m_block_size = 0U;
      }
    }
  }

  void update(std::string_view text) {
    update(reinterpret_cast<const std::uint8_t*>(text.data()), text.size());
  }

  [[nodiscard]] std::array<std::uint8_t, 32> finish() {
    const std::uint64_t original_bits = m_bit_count;
    m_block[m_block_size++] = 0x80U;
    if (m_block_size > 56U) {
      while (m_block_size < m_block.size()) {
        m_block[m_block_size++] = 0U;
      }
      transform();
      m_block_size = 0U;
    }
    while (m_block_size < 56U) {
      m_block[m_block_size++] = 0U;
    }
    for (int shift = 56; shift >= 0; shift -= 8) {
      m_block[m_block_size++] = static_cast<std::uint8_t>((original_bits >> shift) & 0xffU);
    }
    transform();

    std::array<std::uint8_t, 32> digest{};
    for (std::size_t i = 0; i < m_state.size(); ++i) {
      digest[i * 4U + 0U] = static_cast<std::uint8_t>(m_state[i] >> 24U);
      digest[i * 4U + 1U] = static_cast<std::uint8_t>(m_state[i] >> 16U);
      digest[i * 4U + 2U] = static_cast<std::uint8_t>(m_state[i] >> 8U);
      digest[i * 4U + 3U] = static_cast<std::uint8_t>(m_state[i]);
    }
    return digest;
  }

 private:
  static constexpr std::array<std::uint32_t, 64> k_round{
      0x428a2f98U, 0x71374491U, 0xb5c0fbcfU, 0xe9b5dba5U, 0x3956c25bU,
      0x59f111f1U, 0x923f82a4U, 0xab1c5ed5U, 0xd807aa98U, 0x12835b01U,
      0x243185beU, 0x550c7dc3U, 0x72be5d74U, 0x80deb1feU, 0x9bdc06a7U,
      0xc19bf174U, 0xe49b69c1U, 0xefbe4786U, 0x0fc19dc6U, 0x240ca1ccU,
      0x2de92c6fU, 0x4a7484aaU, 0x5cb0a9dcU, 0x76f988daU, 0x983e5152U,
      0xa831c66dU, 0xb00327c8U, 0xbf597fc7U, 0xc6e00bf3U, 0xd5a79147U,
      0x06ca6351U, 0x14292967U, 0x27b70a85U, 0x2e1b2138U, 0x4d2c6dfcU,
      0x53380d13U, 0x650a7354U, 0x766a0abbU, 0x81c2c92eU, 0x92722c85U,
      0xa2bfe8a1U, 0xa81a664bU, 0xc24b8b70U, 0xc76c51a3U, 0xd192e819U,
      0xd6990624U, 0xf40e3585U, 0x106aa070U, 0x19a4c116U, 0x1e376c08U,
      0x2748774cU, 0x34b0bcb5U, 0x391c0cb3U, 0x4ed8aa4aU, 0x5b9cca4fU,
      0x682e6ff3U, 0x748f82eeU, 0x78a5636fU, 0x84c87814U, 0x8cc70208U,
      0x90befffaU, 0xa4506cebU, 0xbef9a3f7U, 0xc67178f2U};

  [[nodiscard]] static constexpr std::uint32_t rotate(std::uint32_t value, unsigned bits) noexcept {
    return (value >> bits) | (value << (32U - bits));
  }

  void transform() {
    std::array<std::uint32_t, 64> words{};
    for (std::size_t i = 0; i < 16U; ++i) {
      words[i] = (static_cast<std::uint32_t>(m_block[i * 4U]) << 24U) |
          (static_cast<std::uint32_t>(m_block[i * 4U + 1U]) << 16U) |
          (static_cast<std::uint32_t>(m_block[i * 4U + 2U]) << 8U) |
          static_cast<std::uint32_t>(m_block[i * 4U + 3U]);
    }
    for (std::size_t i = 16U; i < words.size(); ++i) {
      const std::uint32_t s0 = rotate(words[i - 15U], 7U) ^
          rotate(words[i - 15U], 18U) ^ (words[i - 15U] >> 3U);
      const std::uint32_t s1 = rotate(words[i - 2U], 17U) ^
          rotate(words[i - 2U], 19U) ^ (words[i - 2U] >> 10U);
      words[i] = words[i - 16U] + s0 + words[i - 7U] + s1;
    }

    std::uint32_t a = m_state[0];
    std::uint32_t b = m_state[1];
    std::uint32_t c = m_state[2];
    std::uint32_t d = m_state[3];
    std::uint32_t e = m_state[4];
    std::uint32_t f = m_state[5];
    std::uint32_t g = m_state[6];
    std::uint32_t h = m_state[7];
    for (std::size_t i = 0; i < words.size(); ++i) {
      const std::uint32_t s1 = rotate(e, 6U) ^ rotate(e, 11U) ^ rotate(e, 25U);
      const std::uint32_t choice = (e & f) ^ ((~e) & g);
      const std::uint32_t temp1 = h + s1 + choice + k_round[i] + words[i];
      const std::uint32_t s0 = rotate(a, 2U) ^ rotate(a, 13U) ^ rotate(a, 22U);
      const std::uint32_t majority = (a & b) ^ (a & c) ^ (b & c);
      const std::uint32_t temp2 = s0 + majority;
      h = g;
      g = f;
      f = e;
      e = d + temp1;
      d = c;
      c = b;
      b = a;
      a = temp1 + temp2;
    }
    m_state[0] += a;
    m_state[1] += b;
    m_state[2] += c;
    m_state[3] += d;
    m_state[4] += e;
    m_state[5] += f;
    m_state[6] += g;
    m_state[7] += h;
  }

  std::array<std::uint32_t, 8> m_state{
      0x6a09e667U, 0xbb67ae85U, 0x3c6ef372U, 0xa54ff53aU,
      0x510e527fU, 0x9b05688cU, 0x1f83d9abU, 0x5be0cd19U};
  std::array<std::uint8_t, 64> m_block{};
  std::size_t m_block_size = 0U;
  std::uint64_t m_bit_count = 0U;
};

[[nodiscard]] inline std::string sha256Hex(std::string_view text) {
  static constexpr char k_hex[] = "0123456789abcdef";
  Sha256 hash;
  hash.update(text);
  const auto digest = hash.finish();
  std::string out(64U, '0');
  for (std::size_t i = 0; i < digest.size(); ++i) {
    out[i * 2U] = k_hex[digest[i] >> 4U];
    out[i * 2U + 1U] = k_hex[digest[i] & 0x0fU];
  }
  return out;
}

}  // namespace cosmosim::core::internal
