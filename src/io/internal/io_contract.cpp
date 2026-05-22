#include "cosmosim/io/io_contract.hpp"

#include <stdexcept>

namespace cosmosim::io {

const SharedIoContractNames& sharedIoContractNames() {
  static const SharedIoContractNames names{};
  return names;
}

void validateContinuationMetadata(
    const std::string& normalized_config_text,
    const std::string& normalized_config_hash_hex,
    const core::ProvenanceRecord& provenance,
    std::string_view context_label) {
  if (normalized_config_text.empty()) {
    throw std::invalid_argument(std::string(context_label) + " requires non-empty normalized_config_text");
  }
  if (normalized_config_hash_hex.empty()) {
    throw std::invalid_argument(std::string(context_label) + " requires non-empty normalized_config_hash_hex");
  }
  if (provenance.normalized_config_hash_hex.empty() && provenance.config_hash_hex.empty()) {
    throw std::invalid_argument(
        std::string(context_label) + " requires provenance.normalized_config_hash_hex");
  }

  const std::string computed_hash_hex = core::stableConfigHashHex(normalized_config_text);
  if (normalized_config_hash_hex != computed_hash_hex) {
    throw std::invalid_argument(
        std::string(context_label) + " normalized_config_hash_hex does not match normalized_config_text");
  }
  const std::string& provenance_hash = provenance.normalized_config_hash_hex.empty()
      ? provenance.config_hash_hex
      : provenance.normalized_config_hash_hex;
  if (provenance_hash != normalized_config_hash_hex) {
    throw std::invalid_argument(
        std::string(context_label) +
        " provenance.normalized_config_hash_hex does not match normalized_config_hash_hex");
  }
}

}  // namespace cosmosim::io
