# Compatibility profiles are real files for Windows portability. Keep each
# entry point byte-identical to its canonical scientific profile.
set(cosmosim_config_alias_pairs
  "configs/adaptive_bound_jeans_isolated_galaxy.param.txt|configs/isolated_galaxy/disk_adaptive_bound_jeans_v01.param.txt"
  "configs/zoom_in.param.txt|configs/zoom_in/zoom_adaptive_bound_jeans_v01.param.txt"
  "configs/cosmo_cube.param.txt|configs/cosmo_cube/cube_effective_multiphase_tng_like_v01.param.txt")

foreach(cosmosim_pair IN LISTS cosmosim_config_alias_pairs)
  string(REPLACE "|" ";" cosmosim_pair_paths "${cosmosim_pair}")
  list(GET cosmosim_pair_paths 0 cosmosim_alias_relative)
  list(GET cosmosim_pair_paths 1 cosmosim_canonical_relative)
  set(cosmosim_alias "${COSMOSIM_SOURCE_DIR}/${cosmosim_alias_relative}")
  set(cosmosim_canonical "${COSMOSIM_SOURCE_DIR}/${cosmosim_canonical_relative}")

  if(NOT EXISTS "${cosmosim_alias}" OR NOT EXISTS "${cosmosim_canonical}")
    message(FATAL_ERROR
      "configuration alias contract is incomplete: ${cosmosim_alias_relative} -> ${cosmosim_canonical_relative}")
  endif()

  execute_process(
    COMMAND "${CMAKE_COMMAND}" -E compare_files
            "${cosmosim_alias}" "${cosmosim_canonical}"
    RESULT_VARIABLE cosmosim_compare_result)
  if(NOT cosmosim_compare_result EQUAL 0)
    message(FATAL_ERROR
      "configuration alias drift detected: ${cosmosim_alias_relative} must match ${cosmosim_canonical_relative}")
  endif()
endforeach()

message(STATUS "CosmoSim configuration compatibility aliases match canonical profiles")
