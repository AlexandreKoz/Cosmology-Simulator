# Source identity is separate from compiler/config identity. Always hash the
# source manifest so modified checkouts and source archives remain distinguishable.
set(COSMOSIM_SOURCE_GIT "unavailable")
set(COSMOSIM_SOURCE_GIT_STATE "archive")
find_package(Git QUIET)
if(Git_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
  execute_process(COMMAND "${GIT_EXECUTABLE}" rev-parse HEAD
    WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}" OUTPUT_VARIABLE COSMOSIM_SOURCE_GIT
    OUTPUT_STRIP_TRAILING_WHITESPACE RESULT_VARIABLE _source_git_result ERROR_QUIET)
  if(NOT _source_git_result EQUAL 0)
    set(COSMOSIM_SOURCE_GIT "unavailable")
  else()
    execute_process(COMMAND "${GIT_EXECUTABLE}" status --porcelain
      WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}" OUTPUT_VARIABLE _source_git_status ERROR_QUIET)
    if(_source_git_status STREQUAL "")
      set(COSMOSIM_SOURCE_GIT_STATE "clean")
    else()
      set(COSMOSIM_SOURCE_GIT_STATE "modified")
    endif()
  endif()
endif()
file(GLOB_RECURSE _source_manifest_files RELATIVE "${PROJECT_SOURCE_DIR}" CONFIGURE_DEPENDS
  "${PROJECT_SOURCE_DIR}/include/*.hpp" "${PROJECT_SOURCE_DIR}/src/*.cpp"
  "${PROJECT_SOURCE_DIR}/src/*.hpp" "${PROJECT_SOURCE_DIR}/src/*.cu"
  "${PROJECT_SOURCE_DIR}/cmake/*.cmake" "${PROJECT_SOURCE_DIR}/cmake/*.in"
  "${PROJECT_SOURCE_DIR}/configs/*.param.txt" "${PROJECT_SOURCE_DIR}/docs/*.md"
  "${PROJECT_SOURCE_DIR}/tests/*.cpp" "${PROJECT_SOURCE_DIR}/tests/*.hpp"
  "${PROJECT_SOURCE_DIR}/tests/*.py" "${PROJECT_SOURCE_DIR}/tests/*.in"
  "${PROJECT_SOURCE_DIR}/bench/*.cpp" "${PROJECT_SOURCE_DIR}/bench/*.hpp"
  "${PROJECT_SOURCE_DIR}/tools/*.cpp" "${PROJECT_SOURCE_DIR}/tools/*.hpp"
  "${PROJECT_SOURCE_DIR}/tools/*.py" "${PROJECT_SOURCE_DIR}/scripts/*.sh")
list(APPEND _source_manifest_files CMakeLists.txt CMakePresets.json AGENTS.md README.md CONTRIBUTING.md CURRENT_STATUS.md chui)
list(SORT _source_manifest_files)
set(_source_manifest "chui_source_manifest_v1\n")
foreach(_source_file IN LISTS _source_manifest_files)
  set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${PROJECT_SOURCE_DIR}/${_source_file}")
  file(SHA256 "${PROJECT_SOURCE_DIR}/${_source_file}" _source_hash)
  string(APPEND _source_manifest "${_source_hash}  ${_source_file}\n")
endforeach()
string(SHA256 COSMOSIM_SOURCE_MANIFEST_SHA256 "${_source_manifest}")
file(WRITE "${CMAKE_BINARY_DIR}/source_manifest.sha256" "${_source_manifest}")
if(COSMOSIM_SOURCE_GIT STREQUAL "unavailable")
  set(COSMOSIM_SOURCE_IDENTITY "source-archive-manifest-sha256:${COSMOSIM_SOURCE_MANIFEST_SHA256}")
else()
  set(COSMOSIM_SOURCE_IDENTITY "git:${COSMOSIM_SOURCE_GIT}:${COSMOSIM_SOURCE_GIT_STATE}")
endif()
message(STATUS "CHUI source: ${COSMOSIM_SOURCE_IDENTITY} manifest_sha256=${COSMOSIM_SOURCE_MANIFEST_SHA256}")
