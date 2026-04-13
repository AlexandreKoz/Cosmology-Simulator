#pragma once

#include <string>
#include <string_view>

#include "cosmosim/core/provenance.hpp"

namespace cosmosim::io {

struct SharedIoContractNames {
  std::string_view normalized_config_text_dataset = "normalized_config_text";
  std::string_view normalized_config_hash_hex_attribute = "normalized_config_hash_hex";
  std::string_view provenance_record_dataset = "provenance_record";
};

[[nodiscard]] const SharedIoContractNames& sharedIoContractNames();

void validateContinuationMetadata(
    const std::string& normalized_config_text,
    const std::string& normalized_config_hash_hex,
    const core::ProvenanceRecord& provenance,
    std::string_view context_label);

}  // namespace cosmosim::io
