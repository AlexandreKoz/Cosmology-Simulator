#pragma once

#include <cctype>
#include <stdexcept>
#include <string>
#include <string_view>

namespace cosmosim::core::internal {

[[nodiscard]] inline char hexDigit(unsigned value) noexcept {
  return value < 10U ? static_cast<char>('0' + value)
                     : static_cast<char>('A' + (value - 10U));
}

[[nodiscard]] inline unsigned hexValue(char value) {
  const unsigned char c = static_cast<unsigned char>(value);
  if (c >= static_cast<unsigned char>('0') && c <= static_cast<unsigned char>('9')) {
    return static_cast<unsigned>(c - static_cast<unsigned char>('0'));
  }
  if (c >= static_cast<unsigned char>('a') && c <= static_cast<unsigned char>('f')) {
    return 10U + static_cast<unsigned>(c - static_cast<unsigned char>('a'));
  }
  if (c >= static_cast<unsigned char>('A') && c <= static_cast<unsigned char>('F')) {
    return 10U + static_cast<unsigned>(c - static_cast<unsigned char>('A'));
  }
  throw std::invalid_argument("text codec: invalid hexadecimal escape digit");
}

// Escapes one logical key=value line. UTF-8 bytes >= 0x20 (except DEL) are
// preserved byte-for-byte; line-breaking/control bytes and grammar delimiters
// use an explicit reversible escape grammar.
[[nodiscard]] inline std::string escapeTextLine(std::string_view text) {
  std::string out;
  out.reserve(text.size());
  for (const char raw_c : text) {
    const auto c = static_cast<unsigned char>(raw_c);
    switch (c) {
      case static_cast<unsigned char>('\\'):
        out += "\\\\";
        break;
      case static_cast<unsigned char>('\n'):
        out += "\\n";
        break;
      case static_cast<unsigned char>('\r'):
        out += "\\r";
        break;
      case static_cast<unsigned char>('\t'):
        out += "\\t";
        break;
      case static_cast<unsigned char>('='):
        out += "\\=";
        break;
      default:
        if (c < 0x20U || c == 0x7fU) {
          out += "\\x";
          out.push_back(hexDigit(c >> 4U));
          out.push_back(hexDigit(c & 0x0fU));
        } else {
          out.push_back(static_cast<char>(c));
        }
        break;
    }
  }
  return out;
}

[[nodiscard]] inline std::string unescapeTextLineStrict(
    std::string_view text,
    std::string_view context = "text codec") {
  std::string out;
  out.reserve(text.size());
  for (std::size_t i = 0; i < text.size(); ++i) {
    const char c = text[i];
    if (c != '\\') {
      out.push_back(c);
      continue;
    }
    if (i + 1U >= text.size()) {
      throw std::invalid_argument(std::string(context) + ": dangling escape delimiter");
    }
    const char escaped = text[++i];
    switch (escaped) {
      case '\\':
        out.push_back('\\');
        break;
      case 'n':
        out.push_back('\n');
        break;
      case 'r':
        out.push_back('\r');
        break;
      case 't':
        out.push_back('\t');
        break;
      case '=':
        out.push_back('=');
        break;
      case 'x': {
        if (i + 2U >= text.size()) {
          throw std::invalid_argument(std::string(context) + ": truncated hexadecimal escape");
        }
        const unsigned high = hexValue(text[i + 1U]);
        const unsigned low = hexValue(text[i + 2U]);
        const unsigned value = (high << 4U) | low;
        out.push_back(static_cast<char>(static_cast<unsigned char>(value)));
        i += 2U;
        break;
      }
      default:
        throw std::invalid_argument(
            std::string(context) + ": unknown escape sequence \\" + escaped + "'");
    }
  }
  return out;
}

}  // namespace cosmosim::core::internal
