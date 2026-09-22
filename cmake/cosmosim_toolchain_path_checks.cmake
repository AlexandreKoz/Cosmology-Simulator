# Small path classifiers used by configure-time dependency provenance checks.
# Normalize separators before matching so warning-only logic cannot fail on
# Windows-style backslashes or require fragile escaped regular expressions.
function(cosmosim_normalize_dependency_path input_path output_var)
  file(TO_CMAKE_PATH "${input_path}" _cosmosim_normalized_path)
  string(TOLOWER "${_cosmosim_normalized_path}" _cosmosim_normalized_path)
  set(${output_var} "${_cosmosim_normalized_path}" PARENT_SCOPE)
endfunction()

function(cosmosim_path_looks_conda input_path output_var)
  cosmosim_normalize_dependency_path("${input_path}" _cosmosim_normalized_path)
  set(_cosmosim_is_conda FALSE)
  foreach(_cosmosim_marker IN ITEMS "/anaconda" "/miniconda" "/conda/")
    string(FIND "${_cosmosim_normalized_path}" "${_cosmosim_marker}" _cosmosim_marker_pos)
    if(NOT _cosmosim_marker_pos EQUAL -1)
      set(_cosmosim_is_conda TRUE)
      break()
    endif()
  endforeach()
  set(${output_var} "${_cosmosim_is_conda}" PARENT_SCOPE)
endfunction()
