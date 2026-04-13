#include <cassert>
#include <filesystem>
#include <string>

#include "cosmosim/core/provenance.hpp"

int main() {
  const auto run_directory = std::filesystem::temp_directory_path() / "cosmosim_provenance_roundtrip";
  std::filesystem::remove_all(run_directory);

  cosmosim::core::ProvenanceRecord in;
  in.git_sha = "abc123def";
  in.compiler_id = "test_compiler";
  in.compiler_version = "1.2.3";
  in.build_preset = "test_preset";
  in.enabled_features = "mpi=0,hdf5=0,fftw=0,cuda=0,python=0";
  in.config_hash_hex = "0011223344556677";
  in.timestamp_utc = "2026-04-05T00:00:00Z";
  in.hardware_summary = "logical_threads=8";
  in.author_rank = 0;

  cosmosim::core::writeProvenanceRecord(in, run_directory);

  const auto out = cosmosim::core::readProvenanceRecord(run_directory);
  assert(out.schema_version == in.schema_version);
  assert(out.git_sha == in.git_sha);
  assert(out.compiler_id == in.compiler_id);
  assert(out.compiler_version == in.compiler_version);
  assert(out.build_preset == in.build_preset);
  assert(out.enabled_features == in.enabled_features);
  assert(out.config_hash_hex == in.config_hash_hex);
  assert(out.timestamp_utc == in.timestamp_utc);
  assert(out.hardware_summary == in.hardware_summary);
  assert(out.author_rank == in.author_rank);

  std::filesystem::remove_all(run_directory);
  return 0;
}
