#pragma once

#include <array>
#include <string_view>

namespace cosmosim::core {

constexpr std::array<std::string_view, 9> k_module_names = {
    "core", "gravity", "hydro", "amr", "physics", "io", "analysis", "parallel", "utils"};

std::array<std::string_view, k_module_names.size()> moduleNames();

}  // namespace cosmosim::core
