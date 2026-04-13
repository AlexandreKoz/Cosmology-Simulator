#pragma once

#include <string_view>

namespace cosmosim::core {

enum class ExecutionPolicy {
  kHostSerial,
  kCuda,
};

[[nodiscard]] constexpr std::string_view toString(ExecutionPolicy policy) {
  switch (policy) {
    case ExecutionPolicy::kHostSerial:
      return "host_serial";
    case ExecutionPolicy::kCuda:
      return "cuda";
  }
  return "unknown";
}

}  // namespace cosmosim::core
