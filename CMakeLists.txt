cmake_minimum_required(VERSION 3.24)

project(
  cosmosim
  VERSION 0.1.0
  DESCRIPTION "CosmoSim: desktop-first cosmological simulation framework skeleton"
  LANGUAGES CXX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(cosmosim_feature_summary)
include(cosmosim_find_fftw)
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

function(cosmosim_require_cxx20 target_name)
  # Prefer compile features when CMake knows the compiler feature table.
  # Some CI/container toolchains report a valid compiler but leave
  # CMAKE_CXX_COMPILE_FEATURES empty; falling back to target properties keeps
  # the C++20 contract without making configure fail before tests can run.
  if("cxx_std_20" IN_LIST CMAKE_CXX_COMPILE_FEATURES)
    target_compile_features(${target_name} PUBLIC cxx_std_20)
  else()
    set_target_properties(${target_name} PROPERTIES
      CXX_STANDARD 20
      CXX_STANDARD_REQUIRED ON
      CXX_EXTENSIONS OFF)
  endif()
endfunction()

function(cosmosim_require_private_cxx20 target_name)
  if("cxx_std_20" IN_LIST CMAKE_CXX_COMPILE_FEATURES)
    target_compile_features(${target_name} PRIVATE cxx_std_20)
  else()
    set_target_properties(${target_name} PROPERTIES
      CXX_STANDARD 20
      CXX_STANDARD_REQUIRED ON
      CXX_EXTENSIONS OFF)
  endif()
endfunction()

function(cosmosim_enable_test_assertions target_name)
  if(MSVC)
    target_compile_options(${target_name} PRIVATE /UNDEBUG)
  else()
    target_compile_options(${target_name} PRIVATE -UNDEBUG)
  endif()
endfunction()

function(cosmosim_apply_test_target_defaults)
  get_property(cosmosim_build_targets DIRECTORY PROPERTY BUILDSYSTEM_TARGETS)
  foreach(cosmosim_target_name IN LISTS cosmosim_build_targets)
    get_target_property(cosmosim_target_type ${cosmosim_target_name} TYPE)
    if(cosmosim_target_name MATCHES "^test_" AND
       cosmosim_target_type STREQUAL "EXECUTABLE")
      cosmosim_enable_test_assertions(${cosmosim_target_name})
    endif()
  endforeach()
endfunction()

option(COSMOSIM_ENABLE_MPI "Enable MPI backend" OFF)
option(COSMOSIM_ENABLE_OPENMP "Enable OpenMP shared-memory backend when available" ON)
option(COSMOSIM_ENABLE_HDF5 "Enable HDF5 I/O backend" OFF)
option(COSMOSIM_ENABLE_FFTW "Enable FFTW backend" OFF)
option(COSMOSIM_ENABLE_CUDA "Enable CUDA/cuFFT backend" OFF)
option(COSMOSIM_ENABLE_PYTHON "Enable Python tooling integration" OFF)
option(COSMOSIM_ENABLE_TRACERS "Enable tracer update module" OFF)
option(COSMOSIM_ENABLE_TESTS "Build unit and integration tests" ON)
option(COSMOSIM_ENABLE_BENCHMARKS "Build benchmark/profiling hooks" ON)
option(COSMOSIM_ENABLE_PROFILING "Enable built-in profiling instrumentation" ON)
option(COSMOSIM_ENABLE_LTO "Enable interprocedural optimization (LTO)" OFF)
option(COSMOSIM_STRICT_IO_WARNINGS "Treat first-party cosmosim_io warnings as errors" OFF)

set(COSMOSIM_FFTW_PLAN_RIGOR
    "ESTIMATE"
    CACHE STRING "FFTW planner rigor for PM plans: ESTIMATE, MEASURE, PATIENT, or EXHAUSTIVE")
set_property(CACHE COSMOSIM_FFTW_PLAN_RIGOR PROPERTY STRINGS ESTIMATE MEASURE PATIENT EXHAUSTIVE)
set(_cosmosim_valid_fftw_plan_rigors ESTIMATE MEASURE PATIENT EXHAUSTIVE)
if(NOT COSMOSIM_FFTW_PLAN_RIGOR IN_LIST _cosmosim_valid_fftw_plan_rigors)
  message(FATAL_ERROR
    "COSMOSIM_FFTW_PLAN_RIGOR='${COSMOSIM_FFTW_PLAN_RIGOR}' is invalid; "
    "expected one of: ${_cosmosim_valid_fftw_plan_rigors}")
endif()
set(COSMOSIM_FFTW_PLANNER_FLAGS "FFTW_${COSMOSIM_FFTW_PLAN_RIGOR}")

set(COSMOSIM_CTEST_DEFAULT_TIMEOUT_SECONDS
    "120"
    CACHE STRING "Default timeout in seconds for ordinary CosmoSim CTest tests")
set(COSMOSIM_CTEST_LONG_TIMEOUT_SECONDS
    "300"
    CACHE STRING "Timeout in seconds for long validation, MPI, I/O, and TreePM CTest tests")

function(cosmosim_append_test_labels test_name)
  get_property(cosmosim_existing_labels TEST ${test_name} PROPERTY LABELS)
  if(NOT cosmosim_existing_labels)
    set(cosmosim_existing_labels "")
  endif()
  list(APPEND cosmosim_existing_labels ${ARGN})
  list(REMOVE_DUPLICATES cosmosim_existing_labels)
  set_tests_properties(${test_name} PROPERTIES LABELS "${cosmosim_existing_labels}")
endfunction()

function(cosmosim_mark_runtime_truth_test test_name)
  cosmosim_append_test_labels(${test_name} runtime_truth p0 ${ARGN})
endfunction()

function(cosmosim_apply_test_defaults)
  get_property(cosmosim_registered_tests DIRECTORY PROPERTY TESTS)
  foreach(cosmosim_test_name IN LISTS cosmosim_registered_tests)
    get_property(cosmosim_existing_labels TEST ${cosmosim_test_name} PROPERTY LABELS)
    if(NOT cosmosim_existing_labels)
      set(cosmosim_existing_labels "")
    endif()

    if(cosmosim_test_name MATCHES "^unit_")
      list(APPEND cosmosim_existing_labels unit)
    elseif(cosmosim_test_name MATCHES "^integration_")
      list(APPEND cosmosim_existing_labels integration)
    elseif(cosmosim_test_name MATCHES "^validation_")
      list(APPEND cosmosim_existing_labels validation)
    elseif(cosmosim_test_name MATCHES "^python_")
      list(APPEND cosmosim_existing_labels python)
    endif()

    if(COSMOSIM_ENABLE_MPI AND cosmosim_test_name MATCHES "two_rank|three_rank|four_rank|eight_rank|mpi_np[0-9]+")
      list(APPEND cosmosim_existing_labels mpi distributed)
    endif()
    if(cosmosim_test_name MATCHES "runtime_app_smoke")
      list(APPEND cosmosim_existing_labels runtime smoke)
    endif()

    set(cosmosim_timeout "${COSMOSIM_CTEST_DEFAULT_TIMEOUT_SECONDS}")
    if(cosmosim_test_name MATCHES "mpi|three_rank|two_rank|validation_|reference_workflow|tree_pm|pm_periodic|snapshot|restart|runtime_app_smoke")
      set(cosmosim_timeout "${COSMOSIM_CTEST_LONG_TIMEOUT_SECONDS}")
    endif()

    set_tests_properties(${cosmosim_test_name}
      PROPERTIES
        WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
        TIMEOUT "${cosmosim_timeout}")

    if(cosmosim_test_name MATCHES "four_rank")
      set_tests_properties(${cosmosim_test_name} PROPERTIES PROCESSORS 4)
    elseif(cosmosim_test_name MATCHES "three_rank")
      set_tests_properties(${cosmosim_test_name} PROPERTIES PROCESSORS 3)
    elseif(cosmosim_test_name MATCHES "two_rank")
      set_tests_properties(${cosmosim_test_name} PROPERTIES PROCESSORS 2)
    endif()

    if(cosmosim_existing_labels)
      list(REMOVE_DUPLICATES cosmosim_existing_labels)
      set_tests_properties(${cosmosim_test_name} PROPERTIES LABELS "${cosmosim_existing_labels}")
    endif()
  endforeach()
endfunction()

# Backward-compatibility aliases for earlier build toggles.
if(DEFINED COSMOSIM_BUILD_TESTS)
  set(COSMOSIM_ENABLE_TESTS
      "${COSMOSIM_BUILD_TESTS}"
      CACHE BOOL "Build unit and integration tests" FORCE)
endif()

if(DEFINED COSMOSIM_BUILD_BENCHMARKS)
  set(COSMOSIM_ENABLE_BENCHMARKS
      "${COSMOSIM_BUILD_BENCHMARKS}"
      CACHE BOOL "Build benchmark/profiling hooks" FORCE)
endif()

if(COSMOSIM_ENABLE_CUDA)
  enable_language(CUDA)
endif()

if(COSMOSIM_ENABLE_LTO)
  include(CheckIPOSupported)
  check_ipo_supported(RESULT cosmosim_lto_supported OUTPUT cosmosim_lto_error)
  if(NOT cosmosim_lto_supported)
    message(FATAL_ERROR "COSMOSIM_ENABLE_LTO=ON requested but IPO/LTO is not supported: ${cosmosim_lto_error}")
  endif()
endif()

if(CMAKE_BUILD_TYPE STREQUAL "")
  message(STATUS "No CMAKE_BUILD_TYPE selected; use presets for explicit Debug/Release configuration.")
endif()

# Dependencies gated by explicit feature options.
set(COSMOSIM_HAVE_OPENMP OFF)
if(COSMOSIM_ENABLE_OPENMP)
  find_package(OpenMP QUIET COMPONENTS CXX)
  if(OpenMP_CXX_FOUND)
    set(COSMOSIM_HAVE_OPENMP ON)
    message(STATUS "CosmoSim OpenMP backend enabled")
  else()
    message(STATUS "CosmoSim OpenMP backend requested but unavailable; building serial fallback")
  endif()
endif()

if(COSMOSIM_ENABLE_MPI)
  find_package(MPI REQUIRED COMPONENTS CXX)
endif()

if(COSMOSIM_ENABLE_HDF5)
  find_package(HDF5 QUIET COMPONENTS C CXX)
  if(NOT HDF5_FOUND)
    message(FATAL_ERROR
      "COSMOSIM_ENABLE_HDF5=ON but HDF5 was not found.\n"
      "Install libhdf5 development packages and reconfigure with preset 'hdf5-debug' "
      "or provide -DHDF5_ROOT=/path/to/hdf5.")
  endif()
  if(HDF5_VERSION VERSION_LESS "1.10")
    message(FATAL_ERROR
      "CosmoSim requires HDF5 >= 1.10 for the supported scientific-I/O contract; found ${HDF5_VERSION}.")
  endif()
  if(HDF5_VERSION VERSION_GREATER_EQUAL "1.15")
    message(FATAL_ERROR
      "CosmoSim's currently qualified HDF5 API range is 1.10.x through 1.14.x; found ${HDF5_VERSION}. "
      "Newer HDF5 families must be qualified explicitly before enabling production I/O.")
  endif()
endif()

if(COSMOSIM_ENABLE_FFTW)
  if(COSMOSIM_ENABLE_MPI)
    cosmosim_find_fftw(COSMOSIM_FFTW_TARGET REQUIRE_MPI)
  else()
    cosmosim_find_fftw(COSMOSIM_FFTW_TARGET)
  endif()
endif()

if(COSMOSIM_ENABLE_CUDA)
  find_package(CUDAToolkit REQUIRED)
endif()

if(COSMOSIM_ENABLE_PYTHON)
  find_package(Python3 REQUIRED COMPONENTS Interpreter Development.Module)
  find_package(pybind11 CONFIG QUIET)
  if(NOT pybind11_FOUND)
    message(FATAL_ERROR
      "COSMOSIM_ENABLE_PYTHON=ON but pybind11 was not found.\n"
      "Install pybind11 (pip/conda/system package) or provide -Dpybind11_DIR=/path/to/pybind11/share/cmake/pybind11.")
  endif()
endif()

# Build metadata and generated configuration header.
set(COSMOSIM_SELECTED_PRESET "${COSMOSIM_PRESET_NAME}")
if(COSMOSIM_SELECTED_PRESET STREQUAL "")
  if(DEFINED ENV{CMAKE_PRESET_NAME})
    set(COSMOSIM_SELECTED_PRESET "$ENV{CMAKE_PRESET_NAME}")
  else()
    set(COSMOSIM_SELECTED_PRESET "manual")
  endif()
endif()

string(TOUPPER "${CMAKE_BUILD_TYPE}" COSMOSIM_BUILD_TYPE_UPPER)
set(COSMOSIM_EFFECTIVE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
if(NOT COSMOSIM_BUILD_TYPE_UPPER STREQUAL "" AND DEFINED CMAKE_CXX_FLAGS_${COSMOSIM_BUILD_TYPE_UPPER})
  string(APPEND COSMOSIM_EFFECTIVE_CXX_FLAGS " ${CMAKE_CXX_FLAGS_${COSMOSIM_BUILD_TYPE_UPPER}}")
endif()
string(STRIP "${COSMOSIM_EFFECTIVE_CXX_FLAGS}" COSMOSIM_EFFECTIVE_CXX_FLAGS)
string(REPLACE "\\" "\\\\" COSMOSIM_EFFECTIVE_CXX_FLAGS_ESCAPED "${COSMOSIM_EFFECTIVE_CXX_FLAGS}")
string(REPLACE "\"" "\\\"" COSMOSIM_EFFECTIVE_CXX_FLAGS_ESCAPED "${COSMOSIM_EFFECTIVE_CXX_FLAGS_ESCAPED}")

set(COSMOSIM_GENERATED_INCLUDE_DIR "${CMAKE_BINARY_DIR}/generated/include")
file(MAKE_DIRECTORY "${COSMOSIM_GENERATED_INCLUDE_DIR}/cosmosim/core")

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/templates/build_config.hpp.in"
  "${COSMOSIM_GENERATED_INCLUDE_DIR}/cosmosim/core/build_config.hpp"
  @ONLY)

add_library(cosmosim_core STATIC
  src/core/version.cpp
  src/core/module_registry.cpp
  src/core/config.cpp
  src/core/provenance.cpp
  src/core/cuda_runtime.cpp
  src/core/cosmology.cpp
  src/core/units.cpp
  src/core/simulation_state.cpp
  src/core/simulation_state_structures.cpp
  src/core/simulation_state_ownership.cpp
  src/core/simulation_state_species.cpp
  src/core/simulation_state_metadata.cpp
  src/core/simulation_state_active_views.cpp
  src/core/soa_storage.cpp
  src/core/simulation_mode.cpp
  src/core/time_integration.cpp
  src/core/profiling.cpp
  src/core/memory_governor.cpp
  src/core/memory_accounting.cpp
  src/core/process_memory.cpp
  src/core/openmp_runtime.cpp)

target_include_directories(cosmosim_core
  PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${COSMOSIM_GENERATED_INCLUDE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  PRIVATE
    ${PROJECT_SOURCE_DIR}/src)

cosmosim_require_cxx20(cosmosim_core)

target_compile_definitions(cosmosim_core
  PUBLIC
    COSMOSIM_VERSION_MAJOR=${PROJECT_VERSION_MAJOR}
    COSMOSIM_VERSION_MINOR=${PROJECT_VERSION_MINOR}
    COSMOSIM_VERSION_PATCH=${PROJECT_VERSION_PATCH}
    COSMOSIM_ENABLE_PROFILING=$<BOOL:${COSMOSIM_ENABLE_PROFILING}>)


if(COSMOSIM_HAVE_OPENMP)
  target_link_libraries(cosmosim_core PUBLIC OpenMP::OpenMP_CXX)
endif()

if(COSMOSIM_ENABLE_LTO)
  set_property(TARGET cosmosim_core PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

add_library(cosmosim_gravity STATIC
  src/gravity/gravity_memory.cpp
  src/gravity/pm_solver.cpp
  src/gravity/tree_ordering.cpp
  src/gravity/tree_gravity.cpp
  src/gravity/tree_pm_coupling.cpp)
if(COSMOSIM_ENABLE_CUDA)
  target_sources(cosmosim_gravity PRIVATE src/gravity/pm_cuda_kernels.cu)
  set_target_properties(cosmosim_gravity PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_RESOLVE_DEVICE_SYMBOLS ON
    CUDA_STANDARD 20
    CUDA_STANDARD_REQUIRED ON)
endif()
add_library(cosmosim_hydro STATIC
  src/hydro/hydro_boundary_conditions.cpp
  src/hydro/hydro_cartesian_patch.cpp
  src/hydro/hydro_core_solver.cpp
  src/hydro/hydro_reconstruction.cpp
  src/hydro/hydro_riemann.cpp)
add_library(cosmosim_amr STATIC
  src/amr/amr_framework.cpp
  src/amr/amr_hydro_geometry.cpp
  src/amr/amr_ghost_fill.cpp
  src/amr/amr_hydro_orchestrator.cpp)
add_library(cosmosim_physics STATIC
  src/physics/cooling_heating.cpp
  src/physics/effective_multiphase_ism.cpp
  src/physics/black_hole_agn.cpp
  src/physics/star_formation.cpp
  src/physics/stellar_evolution.cpp
  src/physics/stellar_feedback.cpp
  src/physics/metal_diffusion.cpp
  src/physics/tracer_support.cpp)
add_library(cosmosim_io STATIC
  src/io/internal/io_contract.cpp
  src/io/transactional_file.cpp
  src/io/snapshot_conversion.cpp
  src/io/snapshot_set.cpp
  src/io/snapshot_validator.cpp
  src/io/ic_manifest.cpp
  src/io/ic_record_codec.cpp
  src/io/ic_reader.cpp
  src/io/ic_file_set_manifest.cpp
  src/io/ic_stream_ingestion.cpp
  src/io/ic_serial_ingestion.cpp
  src/io/ic_mpi_collectives.cpp
  src/io/ic_failure_protocol.cpp
  src/io/ic_distributed_audit.cpp
  src/io/ic_distributed_ingestion.cpp
  src/io/ic_canonical_bundle.cpp
  src/io/ic_conversion_catalog.cpp
  src/io/ic_reader_session.cpp
  src/io/ic_byte_codec.cpp
  src/io/ic_sha256.cpp
  src/io/snapshot_hdf5.cpp
  src/io/restart_checkpoint.cpp)
add_library(cosmosim_analysis STATIC
  src/analysis/diagnostics.cpp
  src/analysis/halo_workflow.cpp)
add_library(cosmosim_parallel STATIC
  src/parallel/distributed_memory.cpp)
target_link_libraries(cosmosim_gravity PUBLIC cosmosim_parallel)
target_link_libraries(cosmosim_analysis PUBLIC cosmosim_parallel)

add_library(cosmosim_workflows STATIC
  src/workflows/analysis_runtime.cpp
  src/workflows/gravity_runtime.cpp
  src/workflows/hydro_amr_runtime.cpp
  src/workflows/amr_migration_payload.cpp
  src/workflows/initial_condition_runtime.cpp
  src/workflows/internal/cartesian_gas_cell_layout.cpp
  src/workflows/internal/metal_diffusion_topology.cpp
  src/workflows/internal/gas_cell_ownership.cpp
  src/workflows/internal/migration_wire.cpp
  src/workflows/migration_balance_runtime.cpp
  src/workflows/output_restart_runtime.cpp
  src/workflows/particle_ghost_runtime.cpp
  src/workflows/reference_runtime_composition.cpp
  src/workflows/restart_verification.cpp
  src/workflows/runtime_capabilities.cpp
  src/workflows/runtime_module_registry.cpp
  src/workflows/runtime_resources.cpp
  src/workflows/runtime_services.cpp
  src/workflows/source_runtime.cpp
  src/workflows/time_coordinator.cpp
  src/workflows/snapshot_verification.cpp
  src/workflows/reference_workflow.cpp)

set(COSMOSIM_MODULE_LIBS
  cosmosim_gravity
  cosmosim_hydro
  cosmosim_amr
  cosmosim_physics
  cosmosim_io
  cosmosim_analysis
  cosmosim_parallel
  cosmosim_workflows)

foreach(module_lib IN LISTS COSMOSIM_MODULE_LIBS)
  target_include_directories(${module_lib}
    PUBLIC
      $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
      $<BUILD_INTERFACE:${COSMOSIM_GENERATED_INCLUDE_DIR}>
      $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    PRIVATE
      ${PROJECT_SOURCE_DIR}/src)
  target_link_libraries(${module_lib} PUBLIC cosmosim_core)
  cosmosim_require_cxx20(${module_lib})

  if(COSMOSIM_ENABLE_LTO)
    set_property(TARGET ${module_lib} PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
  endif()
endforeach()

if(COSMOSIM_STRICT_IO_WARNINGS)
  if(MSVC)
    target_compile_options(cosmosim_io PRIVATE /W4 /WX /permissive-)
  else()
    target_compile_options(cosmosim_io PRIVATE
      -Wall -Wextra -Wpedantic -Wconversion -Wsign-conversion -Wshadow
      -Wundef -Wformat=2 -Wnull-dereference -Werror)
  endif()
endif()

target_link_libraries(cosmosim_workflows
  PUBLIC
    cosmosim_analysis
    cosmosim_gravity
    cosmosim_hydro
    cosmosim_amr
    cosmosim_io
    cosmosim_parallel
    cosmosim_physics)

target_link_libraries(cosmosim_io PUBLIC cosmosim_parallel cosmosim_physics)
target_link_libraries(cosmosim_amr PUBLIC cosmosim_hydro)

if(COSMOSIM_ENABLE_MPI)
  target_link_libraries(cosmosim_parallel PUBLIC MPI::MPI_CXX)
endif()

if(COSMOSIM_ENABLE_HDF5)
  if(TARGET HDF5::HDF5)
    target_link_libraries(cosmosim_io PUBLIC HDF5::HDF5)
  else()
    target_include_directories(cosmosim_io PUBLIC ${HDF5_INCLUDE_DIRS})
    target_link_libraries(cosmosim_io PUBLIC ${HDF5_LIBRARIES})
  endif()
endif()

if(COSMOSIM_ENABLE_FFTW)
  target_compile_definitions(cosmosim_gravity PRIVATE COSMOSIM_FFTW_PLANNER_FLAGS=${COSMOSIM_FFTW_PLANNER_FLAGS})
  target_link_libraries(cosmosim_gravity PUBLIC ${COSMOSIM_FFTW_TARGET})
  target_link_libraries(cosmosim_analysis PUBLIC ${COSMOSIM_FFTW_TARGET})
  if(COSMOSIM_ENABLE_MPI)
    # pm_solver.cpp includes mpi.h/fftw3-mpi.h and calls MPI collectives in
    # distributed slab paths. Link MPI directly here instead of relying on an
    # FFTW imported target to propagate MPI's include and link interfaces;
    # distro/pkg-config FFTW MPI targets are not consistent about doing so.
    target_link_libraries(cosmosim_gravity PUBLIC MPI::MPI_CXX ${COSMOSIM_FFTW_MPI_TARGET})
  endif()
endif()

if(COSMOSIM_ENABLE_CUDA)
  target_link_libraries(cosmosim_gravity PUBLIC CUDA::cufft CUDA::cudart)
endif()

add_library(cosmosim_all INTERFACE)
target_link_libraries(cosmosim_all
  INTERFACE
    cosmosim_core
    ${COSMOSIM_MODULE_LIBS})

set_target_properties(cosmosim_core PROPERTIES EXPORT_NAME core)
set_target_properties(cosmosim_gravity PROPERTIES EXPORT_NAME gravity)
set_target_properties(cosmosim_hydro PROPERTIES EXPORT_NAME hydro)
set_target_properties(cosmosim_amr PROPERTIES EXPORT_NAME amr)
set_target_properties(cosmosim_physics PROPERTIES EXPORT_NAME physics)
set_target_properties(cosmosim_io PROPERTIES EXPORT_NAME io)
set_target_properties(cosmosim_analysis PROPERTIES EXPORT_NAME analysis)
set_target_properties(cosmosim_parallel PROPERTIES EXPORT_NAME parallel)
set_target_properties(cosmosim_workflows PROPERTIES EXPORT_NAME workflows)
set_target_properties(cosmosim_all PROPERTIES EXPORT_NAME all)

add_library(CosmoSim::core ALIAS cosmosim_core)
add_library(CosmoSim::gravity ALIAS cosmosim_gravity)
add_library(CosmoSim::hydro ALIAS cosmosim_hydro)
add_library(CosmoSim::amr ALIAS cosmosim_amr)
add_library(CosmoSim::physics ALIAS cosmosim_physics)
add_library(CosmoSim::io ALIAS cosmosim_io)
add_library(CosmoSim::analysis ALIAS cosmosim_analysis)
add_library(CosmoSim::parallel ALIAS cosmosim_parallel)
add_library(CosmoSim::workflows ALIAS cosmosim_workflows)
add_library(CosmoSim::all ALIAS cosmosim_all)

add_executable(cosmosim_harness src/core/harness_main.cpp)
target_link_libraries(cosmosim_harness PRIVATE cosmosim_all)


set(COSMOSIM_INSTALL_CMAKEDIR "${CMAKE_INSTALL_LIBDIR}/cmake/CosmoSim")
set(COSMOSIM_INSTALL_TARGETS
  cosmosim_core
  cosmosim_gravity
  cosmosim_hydro
  cosmosim_amr
  cosmosim_physics
  cosmosim_io
  cosmosim_analysis
  cosmosim_parallel
  cosmosim_workflows
  cosmosim_all)

install(TARGETS ${COSMOSIM_INSTALL_TARGETS}
  EXPORT CosmoSimTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(DIRECTORY "${PROJECT_SOURCE_DIR}/include/cosmosim"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
  FILES_MATCHING PATTERN "*.hpp")
install(FILES
  "${COSMOSIM_GENERATED_INCLUDE_DIR}/cosmosim/core/build_config.hpp"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/cosmosim/core")

configure_package_config_file(
  "${PROJECT_SOURCE_DIR}/cmake/templates/CosmoSimConfig.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/CosmoSimConfig.cmake"
  INSTALL_DESTINATION "${COSMOSIM_INSTALL_CMAKEDIR}")
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/CosmoSimConfigVersion.cmake"
  VERSION "${PROJECT_VERSION}"
  COMPATIBILITY SameMajorVersion)

install(EXPORT CosmoSimTargets
  FILE CosmoSimTargets.cmake
  NAMESPACE CosmoSim::
  DESTINATION "${COSMOSIM_INSTALL_CMAKEDIR}")
install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/CosmoSimConfig.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/CosmoSimConfigVersion.cmake"
  "${PROJECT_SOURCE_DIR}/cmake/cosmosim_find_fftw.cmake"
  "${PROJECT_SOURCE_DIR}/cmake/FindFFTW3.cmake"
  DESTINATION "${COSMOSIM_INSTALL_CMAKEDIR}")


if(COSMOSIM_ENABLE_PYTHON)
  pybind11_add_module(_cosmosim MODULE src/python_bindings.cpp)
  target_link_libraries(_cosmosim PRIVATE cosmosim_all)
  cosmosim_require_private_cxx20(_cosmosim)

  set_target_properties(_cosmosim PROPERTIES
    OUTPUT_NAME "_cosmosim"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/python/cosmosim")

  add_custom_target(cosmosim_python_package ALL
    COMMAND ${CMAKE_COMMAND} -E make_directory "${CMAKE_BINARY_DIR}/python/cosmosim"
    COMMAND ${CMAKE_COMMAND} -E copy "${PROJECT_SOURCE_DIR}/python/cosmosim/__init__.py" "${CMAKE_BINARY_DIR}/python/cosmosim/__init__.py"
    COMMAND ${CMAKE_COMMAND} -E copy "${PROJECT_SOURCE_DIR}/python/cosmosim/analysis.py" "${CMAKE_BINARY_DIR}/python/cosmosim/analysis.py"
    DEPENDS _cosmosim)
endif()

if(COSMOSIM_ENABLE_TESTS)
  enable_testing()

  add_executable(
    test_unit_assertions_enabled
    tests/unit/test_assertions_enabled.cpp)
  add_test(NAME unit_assertions_enabled COMMAND test_unit_assertions_enabled)

  add_executable(test_unit_header_smoke tests/unit/test_header_smoke.cpp)
  target_link_libraries(test_unit_header_smoke PRIVATE cosmosim_all)
  add_test(NAME unit_header_smoke COMMAND test_unit_header_smoke)

  add_executable(test_unit_benchmark_harness_smoke tests/unit/test_benchmark_harness_smoke.cpp)
  target_link_libraries(test_unit_benchmark_harness_smoke PRIVATE cosmosim_core)
  target_include_directories(test_unit_benchmark_harness_smoke PRIVATE ${PROJECT_SOURCE_DIR})
  add_test(NAME unit_benchmark_harness_smoke COMMAND test_unit_benchmark_harness_smoke)

  add_executable(test_integration_harness tests/integration/test_integration_harness.cpp)
  target_link_libraries(test_integration_harness PRIVATE cosmosim_all)
  add_test(NAME integration_harness COMMAND test_integration_harness)

  add_executable(test_integration_reference_workflow tests/integration/test_reference_workflow_end_to_end.cpp)
  target_link_libraries(test_integration_reference_workflow PRIVATE cosmosim_all)
  add_test(NAME integration_reference_workflow COMMAND test_integration_reference_workflow)

  add_executable(
    test_integration_runtime_descriptor_probe
    tests/integration/test_runtime_descriptor_probe.cpp)
  target_link_libraries(test_integration_runtime_descriptor_probe PRIVATE cosmosim_all)
  add_test(
    NAME integration_runtime_descriptor_probe
    COMMAND test_integration_runtime_descriptor_probe)
  cosmosim_mark_runtime_truth_test(
    integration_runtime_descriptor_probe workflow integration runtime_registry)

  if(COSMOSIM_ENABLE_HDF5)
  add_executable(
    test_integration_reference_workflow_hydro_row_order
    tests/integration/test_reference_workflow_hydro_row_order.cpp)
  target_link_libraries(test_integration_reference_workflow_hydro_row_order PRIVATE cosmosim_all)
  add_test(
    NAME integration_reference_workflow_hydro_row_order
    COMMAND test_integration_reference_workflow_hydro_row_order)
  cosmosim_mark_runtime_truth_test(
    integration_reference_workflow_hydro_row_order
    restart scheduler hydro gas workflow geometry integration)

  add_executable(
    test_integration_restart_equivalence_reference_workflow_hydro
    tests/integration/test_restart_equivalence_reference_workflow_hydro.cpp)
  target_link_libraries(test_integration_restart_equivalence_reference_workflow_hydro PRIVATE cosmosim_all)
  add_test(
    NAME integration_restart_equivalence_reference_workflow_hydro
    COMMAND test_integration_restart_equivalence_reference_workflow_hydro)
  cosmosim_mark_runtime_truth_test(
    integration_restart_equivalence_reference_workflow_hydro
    restart scheduler hydro gas workflow equivalence integration)
  endif()
  if(COSMOSIM_ENABLE_MPI)
    add_executable(
      test_integration_reference_workflow_distributed_treepm_mpi
      tests/integration/test_reference_workflow_distributed_treepm_mpi.cpp)
    target_link_libraries(test_integration_reference_workflow_distributed_treepm_mpi PRIVATE cosmosim_all)
    add_test(
      NAME integration_reference_workflow_distributed_treepm_mpi_two_rank
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
              $<TARGET_FILE:test_integration_reference_workflow_distributed_treepm_mpi> ${MPIEXEC_POSTFLAGS})
    cosmosim_mark_runtime_truth_test(
      integration_reference_workflow_distributed_treepm_mpi_two_rank
      mpi distributed treepm pm restart workflow integration)

    add_executable(
      test_integration_tree_pseudo_hierarchy_exchange_mpi
      tests/integration/test_tree_pseudo_hierarchy_exchange_mpi.cpp)
    target_link_libraries(test_integration_tree_pseudo_hierarchy_exchange_mpi PRIVATE cosmosim_parallel)
    add_test(
      NAME integration_tree_pseudo_hierarchy_exchange_mpi_two_rank
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
              $<TARGET_FILE:test_integration_tree_pseudo_hierarchy_exchange_mpi> ${MPIEXEC_POSTFLAGS})
    cosmosim_mark_runtime_truth_test(
      integration_tree_pseudo_hierarchy_exchange_mpi_two_rank
      mpi distributed treepm hierarchy protocol guard integration)
  endif()

  if(COSMOSIM_ENABLE_MPI)
    add_executable(
      test_integration_reference_workflow_distributed_isolated_pm_mpi
      tests/integration/test_reference_workflow_distributed_isolated_pm_mpi.cpp)
    target_link_libraries(test_integration_reference_workflow_distributed_isolated_pm_mpi PRIVATE cosmosim_all)
    add_test(
      NAME integration_reference_workflow_distributed_isolated_pm_mpi_two_rank
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
              $<TARGET_FILE:test_integration_reference_workflow_distributed_isolated_pm_mpi> ${MPIEXEC_POSTFLAGS})
    cosmosim_mark_runtime_truth_test(
      integration_reference_workflow_distributed_isolated_pm_mpi_two_rank
      mpi distributed pm isolated_pm guard workflow integration)
  endif()

  if(COSMOSIM_ENABLE_MPI AND COSMOSIM_ENABLE_HDF5)
    add_executable(
      test_integration_reference_workflow_distributed_hydro_mpi
      tests/integration/test_reference_workflow_distributed_hydro_mpi.cpp)
    target_link_libraries(test_integration_reference_workflow_distributed_hydro_mpi PRIVATE cosmosim_all)
    add_test(
      NAME integration_reference_workflow_distributed_hydro_mpi_two_rank
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
              $<TARGET_FILE:test_integration_reference_workflow_distributed_hydro_mpi> ${MPIEXEC_POSTFLAGS})
    cosmosim_mark_runtime_truth_test(
      integration_reference_workflow_distributed_hydro_mpi_two_rank
      mpi distributed restart hydro gas workflow conservation integration)
  endif()

  add_executable(test_unit_simulation_state tests/unit/test_simulation_state.cpp)
  target_link_libraries(test_unit_simulation_state PRIVATE cosmosim_all)
  add_test(NAME unit_simulation_state COMMAND test_unit_simulation_state)
  cosmosim_mark_runtime_truth_test(unit_simulation_state active_views scheduler gas sidecar migration)
  add_executable(test_unit_stage6_layout_contracts tests/unit/test_stage6_layout_contracts.cpp)
  target_link_libraries(test_unit_stage6_layout_contracts PRIVATE cosmosim_core)
  add_test(NAME unit_stage6_layout_contracts COMMAND test_unit_stage6_layout_contracts)
  cosmosim_mark_runtime_truth_test(unit_stage6_layout_contracts active_views sidecar softening)
  add_executable(test_unit_gas_cell_identity_invariants tests/unit/test_gas_cell_identity_invariants.cpp)
  target_link_libraries(test_unit_gas_cell_identity_invariants PRIVATE cosmosim_core)
  add_test(NAME unit_gas_cell_identity_invariants COMMAND test_unit_gas_cell_identity_invariants)
  cosmosim_mark_runtime_truth_test(unit_gas_cell_identity_invariants gas sidecar)

  add_executable(test_unit_config_parser tests/unit/test_config_parser.cpp)
  target_link_libraries(test_unit_config_parser PRIVATE cosmosim_all)
  add_test(NAME unit_config_parser COMMAND test_unit_config_parser)
  add_executable(test_unit_units_cosmology_provenance tests/unit/test_units_cosmology_provenance.cpp)
  target_link_libraries(test_unit_units_cosmology_provenance PRIVATE cosmosim_core)
  add_test(NAME unit_units_cosmology_provenance COMMAND test_unit_units_cosmology_provenance)

  add_executable(test_unit_soa_storage tests/unit/test_soa_storage.cpp)
  target_link_libraries(test_unit_soa_storage PRIVATE cosmosim_core)
  add_test(NAME unit_soa_storage COMMAND test_unit_soa_storage)



  add_executable(test_unit_profiling tests/unit/test_profiling.cpp)
  target_link_libraries(test_unit_profiling PRIVATE cosmosim_core)
  add_test(NAME unit_profiling COMMAND test_unit_profiling)

  add_executable(test_unit_memory_accounting tests/unit/test_memory_accounting.cpp)
  target_link_libraries(test_unit_memory_accounting PRIVATE cosmosim_core)
  add_test(NAME unit_memory_accounting COMMAND test_unit_memory_accounting)

  add_executable(test_unit_process_memory tests/unit/test_process_memory.cpp)
  target_link_libraries(test_unit_process_memory PRIVATE cosmosim_core)
  add_test(NAME unit_process_memory COMMAND test_unit_process_memory)

  add_executable(test_unit_memory_governor tests/unit/test_memory_governor.cpp)
  target_link_libraries(test_unit_memory_governor PRIVATE cosmosim_core)
  add_test(NAME unit_memory_governor COMMAND test_unit_memory_governor)

  add_executable(test_unit_scratch_memory_governor tests/unit/test_scratch_memory_governor.cpp)
  target_link_libraries(test_unit_scratch_memory_governor PRIVATE cosmosim_core)
  add_test(NAME unit_scratch_memory_governor COMMAND test_unit_scratch_memory_governor)

  add_executable(test_unit_stage6_source_diagnostics_views tests/unit/test_stage6_source_diagnostics_views.cpp)
  target_link_libraries(test_unit_stage6_source_diagnostics_views PRIVATE cosmosim_all)
  target_compile_definitions(test_unit_stage6_source_diagnostics_views PRIVATE COSMOSIM_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
  add_test(NAME unit_stage6_source_diagnostics_views COMMAND test_unit_stage6_source_diagnostics_views)
  cosmosim_mark_runtime_truth_test(unit_stage6_source_diagnostics_views active_views sidecar)

  add_executable(test_unit_stage6_final_acceptance tests/unit/test_stage6_final_acceptance.cpp)
  target_link_libraries(test_unit_stage6_final_acceptance PRIVATE cosmosim_all)
  add_test(NAME unit_stage6_final_acceptance COMMAND test_unit_stage6_final_acceptance)
  cosmosim_mark_runtime_truth_test(unit_stage6_final_acceptance active_views sidecar restart memory diagnostics)

  add_executable(test_unit_time_integration tests/unit/test_time_integration.cpp)
  target_link_libraries(test_unit_time_integration PRIVATE cosmosim_core)
  add_test(NAME unit_time_integration COMMAND test_unit_time_integration)
  cosmosim_mark_runtime_truth_test(unit_time_integration scheduler restart treepm)
  # Compatibility shim for established build invocations.
  add_custom_target(test_time_integration DEPENDS test_unit_time_integration)

  add_executable(test_unit_runtime_capabilities tests/unit/test_runtime_capabilities.cpp)
  target_link_libraries(test_unit_runtime_capabilities PRIVATE cosmosim_workflows)
  add_test(NAME unit_runtime_capabilities COMMAND test_unit_runtime_capabilities)

  add_executable(test_unit_runtime_services tests/unit/test_runtime_services.cpp)
  target_link_libraries(test_unit_runtime_services PRIVATE cosmosim_workflows)
  add_test(NAME unit_runtime_services COMMAND test_unit_runtime_services)

  add_executable(test_unit_runtime_resources tests/unit/test_runtime_resources.cpp)
  target_link_libraries(test_unit_runtime_resources PRIVATE cosmosim_workflows)
  add_test(NAME unit_runtime_resources COMMAND test_unit_runtime_resources)

  add_executable(test_unit_runtime_module_registry tests/unit/test_runtime_module_registry.cpp)
  target_link_libraries(test_unit_runtime_module_registry PRIVATE cosmosim_workflows)
  add_test(NAME unit_runtime_module_registry COMMAND test_unit_runtime_module_registry)

  add_executable(test_unit_module_registry tests/unit/test_module_registry.cpp)
  target_link_libraries(test_unit_module_registry PRIVATE cosmosim_core)
  add_test(NAME unit_module_registry COMMAND test_unit_module_registry)

  add_executable(test_unit_simulation_mode tests/unit/test_simulation_mode.cpp)
  target_link_libraries(test_unit_simulation_mode PRIVATE cosmosim_core)
  add_test(NAME unit_simulation_mode COMMAND test_unit_simulation_mode)
  add_executable(test_unit_species_state_organization tests/unit/test_species_state_organization.cpp)
  target_link_libraries(test_unit_species_state_organization PRIVATE cosmosim_core)
  add_test(NAME unit_species_state_organization COMMAND test_unit_species_state_organization)

  add_executable(test_unit_hot_cold_sidecar_layout tests/unit/test_hot_cold_sidecar_layout.cpp)
  target_link_libraries(test_unit_hot_cold_sidecar_layout PRIVATE cosmosim_core)
  add_test(NAME unit_hot_cold_sidecar_layout COMMAND test_unit_hot_cold_sidecar_layout)
  cosmosim_mark_runtime_truth_test(unit_hot_cold_sidecar_layout sidecar)


  add_executable(test_unit_parallel_distributed_memory tests/unit/test_parallel_distributed_memory.cpp)
  target_link_libraries(test_unit_parallel_distributed_memory PRIVATE cosmosim_parallel)
  add_test(NAME unit_parallel_distributed_memory COMMAND test_unit_parallel_distributed_memory)

  add_executable(test_unit_hydro_reconstruction tests/unit/test_hydro_reconstruction.cpp)
  target_link_libraries(test_unit_hydro_reconstruction PRIVATE cosmosim_all)
  add_test(NAME unit_hydro_reconstruction COMMAND test_unit_hydro_reconstruction)

  add_executable(test_unit_hydro_boundary_conditions tests/unit/test_hydro_boundary_conditions.cpp)
  target_link_libraries(test_unit_hydro_boundary_conditions PRIVATE cosmosim_hydro)
  add_test(NAME unit_hydro_boundary_conditions COMMAND test_unit_hydro_boundary_conditions)

  add_executable(test_unit_hydro_riemann tests/unit/test_hydro_riemann.cpp)
  target_link_libraries(test_unit_hydro_riemann PRIVATE cosmosim_all)
  add_test(NAME unit_hydro_riemann COMMAND test_unit_hydro_riemann)


  add_executable(test_unit_cooling_heating tests/unit/test_cooling_heating.cpp)
  target_link_libraries(test_unit_cooling_heating PRIVATE cosmosim_all)
  add_test(NAME unit_cooling_heating COMMAND test_unit_cooling_heating)


  add_executable(test_unit_star_formation tests/unit/test_star_formation.cpp)
  target_link_libraries(test_unit_star_formation PRIVATE cosmosim_all)
  add_test(NAME unit_star_formation COMMAND test_unit_star_formation)
  add_executable(test_unit_effective_multiphase_ism tests/unit/test_effective_multiphase_ism.cpp)
  target_link_libraries(test_unit_effective_multiphase_ism PRIVATE cosmosim_all)
  add_test(NAME unit_effective_multiphase_ism COMMAND test_unit_effective_multiphase_ism)
  add_executable(test_unit_black_hole_agn tests/unit/test_black_hole_agn.cpp)
  target_link_libraries(test_unit_black_hole_agn PRIVATE cosmosim_all)
  add_test(NAME unit_black_hole_agn COMMAND test_unit_black_hole_agn)

  add_executable(test_unit_stellar_evolution tests/unit/test_stellar_evolution.cpp)
  target_link_libraries(test_unit_stellar_evolution PRIVATE cosmosim_all)
  target_compile_definitions(test_unit_stellar_evolution PRIVATE COSMOSIM_SOURCE_DIR="${PROJECT_SOURCE_DIR}")
  add_test(NAME unit_stellar_evolution COMMAND test_unit_stellar_evolution)
  add_executable(test_unit_metal_diffusion tests/unit/test_metal_diffusion.cpp)
  target_link_libraries(test_unit_metal_diffusion PRIVATE cosmosim_all)
  add_test(NAME unit_metal_diffusion COMMAND test_unit_metal_diffusion)
  add_executable(test_unit_workflow_metal_diffusion_topology
    tests/unit/test_workflow_metal_diffusion_topology.cpp)
  target_include_directories(test_unit_workflow_metal_diffusion_topology PRIVATE ${PROJECT_SOURCE_DIR}/src)
  target_link_libraries(test_unit_workflow_metal_diffusion_topology PRIVATE cosmosim_all)
  add_test(NAME unit_workflow_metal_diffusion_topology COMMAND test_unit_workflow_metal_diffusion_topology)
  add_executable(test_unit_stellar_feedback tests/unit/test_stellar_feedback.cpp)
  target_link_libraries(test_unit_stellar_feedback PRIVATE cosmosim_all)
  add_test(NAME unit_stellar_feedback COMMAND test_unit_stellar_feedback)
  add_executable(test_unit_tracer_support tests/unit/test_tracer_support.cpp)
  target_link_libraries(test_unit_tracer_support PRIVATE cosmosim_all)
  add_test(NAME unit_tracer_support COMMAND test_unit_tracer_support)

  add_executable(test_integration_cooling_relaxation tests/integration/test_cooling_relaxation.cpp)
  target_link_libraries(test_integration_cooling_relaxation PRIVATE cosmosim_all)
  add_test(NAME integration_cooling_relaxation COMMAND test_integration_cooling_relaxation)

  add_executable(test_integration_stellar_evolution_ssp tests/integration/test_stellar_evolution_ssp.cpp)
  target_link_libraries(test_integration_stellar_evolution_ssp PRIVATE cosmosim_all)
  target_compile_definitions(test_integration_stellar_evolution_ssp PRIVATE COSMOSIM_SOURCE_DIR="${PROJECT_SOURCE_DIR}")
  add_test(NAME integration_stellar_evolution_ssp COMMAND test_integration_stellar_evolution_ssp)

  add_executable(test_integration_star_formation_box tests/integration/test_star_formation_box.cpp)
  target_link_libraries(test_integration_star_formation_box PRIVATE cosmosim_all)
  add_test(NAME integration_star_formation_box COMMAND test_integration_star_formation_box)

  add_executable(test_integration_star_formation_source_runtime tests/integration/test_star_formation_source_runtime.cpp)
  target_link_libraries(test_integration_star_formation_source_runtime PRIVATE cosmosim_all)
  add_test(NAME integration_star_formation_source_runtime COMMAND test_integration_star_formation_source_runtime)
  set_tests_properties(integration_star_formation_source_runtime PROPERTIES LABELS "integration;runtime_truth;star_formation;p0" RUN_SERIAL TRUE)
  foreach(star_formation_amr_test_name IN ITEMS
      integration_star_formation_amr_covered_coarse
      integration_star_formation_amr_refine_derefine
      integration_star_formation_amr_level_equivalence
      integration_star_formation_amr_reflux_ordering
      integration_star_formation_amr_patch_reorder
      integration_effective_ism_amr_threshold_invariance
      integration_effective_ism_amr_eos_restriction)
    add_test(NAME ${star_formation_amr_test_name} COMMAND test_integration_star_formation_source_runtime)
    set_tests_properties(${star_formation_amr_test_name} PROPERTIES LABELS "integration;runtime_truth;star_formation;amr" RUN_SERIAL TRUE)
  endforeach()

  if(COSMOSIM_ENABLE_MPI AND COSMOSIM_ENABLE_HDF5)
    add_executable(test_integration_star_formation_source_runtime_mpi tests/integration/test_star_formation_source_runtime_mpi.cpp)
    target_link_libraries(test_integration_star_formation_source_runtime_mpi PRIVATE cosmosim_all)

    add_executable(test_integration_star_formation_mpi tests/integration/test_star_formation_mpi.cpp)
    target_link_libraries(test_integration_star_formation_mpi PRIVATE cosmosim_all)
    foreach(star_formation_mpi_mode IN ITEMS base repartition restart)
      if(star_formation_mpi_mode STREQUAL "base")
        set(star_formation_mpi_test_name integration_star_formation_mpi_two_rank)
      elseif(star_formation_mpi_mode STREQUAL "repartition")
        set(star_formation_mpi_test_name integration_star_formation_mpi_repartition_two_rank)
      else()
        set(star_formation_mpi_test_name integration_star_formation_restart_mpi_two_rank)
      endif()
      add_test(
        NAME ${star_formation_mpi_test_name}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
                $<TARGET_FILE:test_integration_star_formation_mpi>
                ${star_formation_mpi_mode} ${MPIEXEC_POSTFLAGS})
      cosmosim_mark_runtime_truth_test(
        ${star_formation_mpi_test_name}
        mpi distributed star_formation restart ownership reproducibility integration)
    endforeach()
    foreach(gravity_gas_np IN ITEMS 2 3)
      add_test(
        NAME integration_gravity_gas_ownership_mpi_np${gravity_gas_np}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${gravity_gas_np} ${MPIEXEC_PREFLAGS}
                $<TARGET_FILE:test_integration_star_formation_source_runtime_mpi>
                gravity_gas_ownership ${MPIEXEC_POSTFLAGS})
      cosmosim_mark_runtime_truth_test(
        integration_gravity_gas_ownership_mpi_np${gravity_gas_np}
        mpi distributed gravity gas ownership integration)
    endforeach()
    if(COSMOSIM_ENABLE_FFTW)
      foreach(gravity_gas_np IN ITEMS 2 3)
        add_test(
          NAME integration_gravity_gas_force_mpi_np${gravity_gas_np}
          COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${gravity_gas_np} ${MPIEXEC_PREFLAGS}
                  $<TARGET_FILE:test_integration_star_formation_source_runtime_mpi>
                  gravity_gas_force ${MPIEXEC_POSTFLAGS})
        cosmosim_mark_runtime_truth_test(
          integration_gravity_gas_force_mpi_np${gravity_gas_np}
          mpi fftw distributed gravity gas numerical_force source_only target_only integration)
      endforeach()
    endif()
    foreach(star_formation_runtime_mpi_mode IN ITEMS source_runtime empty_rank id_registry restart one_vs_two_rank_equivalence eos_hash)
      if(star_formation_runtime_mpi_mode STREQUAL "source_runtime")
        set(star_formation_runtime_mpi_test_name integration_star_formation_source_runtime_mpi_two_rank)
      elseif(star_formation_runtime_mpi_mode STREQUAL "empty_rank")
        set(star_formation_runtime_mpi_test_name integration_star_formation_empty_rank_mpi_two_rank)
      elseif(star_formation_runtime_mpi_mode STREQUAL "id_registry")
        set(star_formation_runtime_mpi_test_name integration_star_formation_id_registry_mpi_two_rank)
      elseif(star_formation_runtime_mpi_mode STREQUAL "restart")
        set(star_formation_runtime_mpi_test_name integration_star_formation_restart_source_runtime_mpi_two_rank)
      elseif(star_formation_runtime_mpi_mode STREQUAL "one_vs_two_rank_equivalence")
        set(star_formation_runtime_mpi_test_name integration_star_formation_one_vs_two_rank_equivalence)
      else()
        set(star_formation_runtime_mpi_test_name integration_effective_ism_eos_hash_mpi_two_rank)
      endif()
      add_test(
        NAME ${star_formation_runtime_mpi_test_name}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
                $<TARGET_FILE:test_integration_star_formation_source_runtime_mpi>
                ${star_formation_runtime_mpi_mode} ${MPIEXEC_POSTFLAGS})
      cosmosim_mark_runtime_truth_test(
        ${star_formation_runtime_mpi_test_name}
        mpi distributed star_formation source_runtime amr restart ownership reproducibility integration)
    endforeach()
  endif()
  add_executable(test_integration_black_hole_agn_toy tests/integration/test_black_hole_agn_toy.cpp)
  target_link_libraries(test_integration_black_hole_agn_toy PRIVATE cosmosim_all)
  add_test(NAME integration_black_hole_agn_toy COMMAND test_integration_black_hole_agn_toy)
  add_executable(test_integration_stellar_feedback_box tests/integration/test_stellar_feedback_box.cpp)
  target_link_libraries(test_integration_stellar_feedback_box PRIVATE cosmosim_all)
  add_test(NAME integration_stellar_feedback_box COMMAND test_integration_stellar_feedback_box)

  add_executable(test_integration_simulation_state_pipeline tests/integration/test_simulation_state_pipeline.cpp)
  target_link_libraries(test_integration_simulation_state_pipeline PRIVATE cosmosim_all)
  add_test(NAME integration_simulation_state_pipeline COMMAND test_integration_simulation_state_pipeline)
  add_executable(test_integration_stage6_active_views tests/integration/test_stage6_active_views.cpp)
  target_link_libraries(test_integration_stage6_active_views PRIVATE cosmosim_core)
  add_test(NAME integration_stage6_active_views COMMAND test_integration_stage6_active_views)
  cosmosim_mark_runtime_truth_test(integration_stage6_active_views active_views)

  add_executable(test_integration_analysis_bundle tests/integration/test_analysis_bundle_pipeline.cpp)
  target_link_libraries(test_integration_analysis_bundle PRIVATE cosmosim_all)
  add_test(NAME integration_analysis_bundle COMMAND test_integration_analysis_bundle)
  add_executable(test_integration_halo_catalog_pipeline tests/integration/test_halo_catalog_pipeline.cpp)
  target_link_libraries(test_integration_halo_catalog_pipeline PRIVATE cosmosim_all)
  add_test(NAME integration_halo_catalog_pipeline COMMAND test_integration_halo_catalog_pipeline)


  add_executable(test_integration_soa_species_pipeline tests/integration/test_soa_species_pipeline.cpp)
  target_link_libraries(test_integration_soa_species_pipeline PRIVATE cosmosim_core)
  add_test(NAME integration_soa_species_pipeline COMMAND test_integration_soa_species_pipeline)

  add_executable(test_integration_species_mixed_state tests/integration/test_species_mixed_state.cpp)
  target_link_libraries(test_integration_species_mixed_state PRIVATE cosmosim_core)
  add_test(NAME integration_species_mixed_state COMMAND test_integration_species_mixed_state)

  add_executable(test_integration_reorder_compaction_sidecars tests/integration/test_reorder_compaction_sidecars.cpp)
  target_link_libraries(test_integration_reorder_compaction_sidecars PRIVATE cosmosim_core)
  add_test(NAME integration_reorder_compaction_sidecars COMMAND test_integration_reorder_compaction_sidecars)
  cosmosim_mark_runtime_truth_test(integration_reorder_compaction_sidecars sidecar active_views)

  add_executable(test_integration_species_migration_invariants tests/integration/test_species_migration_invariants.cpp)
  target_link_libraries(test_integration_species_migration_invariants PRIVATE cosmosim_core cosmosim_gravity)
  add_test(NAME integration_species_migration_invariants COMMAND test_integration_species_migration_invariants)
  cosmosim_mark_runtime_truth_test(integration_species_migration_invariants migration sidecar)

  add_executable(test_unit_migration_wire tests/unit/test_migration_wire.cpp)
  target_include_directories(test_unit_migration_wire PRIVATE ${PROJECT_SOURCE_DIR}/src)
  target_link_libraries(test_unit_migration_wire PRIVATE cosmosim_workflows)
  add_test(NAME unit_migration_wire COMMAND test_unit_migration_wire)
  cosmosim_mark_runtime_truth_test(unit_migration_wire migration memory)

  add_executable(test_integration_gas_cell_identity_migration tests/integration/test_gas_cell_identity_migration.cpp)
  target_link_libraries(test_integration_gas_cell_identity_migration PRIVATE cosmosim_core)
  add_test(NAME integration_gas_cell_identity_migration COMMAND test_integration_gas_cell_identity_migration)
  cosmosim_mark_runtime_truth_test(integration_gas_cell_identity_migration migration gas sidecar)

  add_executable(test_integration_amr_patch_migration tests/integration/test_amr_patch_migration.cpp)
  target_link_libraries(test_integration_amr_patch_migration PRIVATE cosmosim_core cosmosim_amr cosmosim_hydro)
  add_test(NAME integration_amr_patch_migration COMMAND test_integration_amr_patch_migration)
  cosmosim_mark_runtime_truth_test(integration_amr_patch_migration amr hydro gas migration)

  add_executable(test_integration_gas_cell_split_merge_remap tests/integration/test_gas_cell_split_merge_remap.cpp)
  target_link_libraries(test_integration_gas_cell_split_merge_remap PRIVATE cosmosim_core cosmosim_hydro)
  add_test(NAME integration_gas_cell_split_merge_remap COMMAND test_integration_gas_cell_split_merge_remap)
  cosmosim_mark_runtime_truth_test(integration_gas_cell_split_merge_remap gas hydro sidecar)

  add_executable(test_integration_transform_fuzz_invariants tests/integration/test_transform_fuzz_invariants.cpp)
  target_link_libraries(test_integration_transform_fuzz_invariants PRIVATE cosmosim_all)
  add_test(NAME integration_transform_fuzz_invariants COMMAND test_integration_transform_fuzz_invariants)
  cosmosim_mark_runtime_truth_test(integration_transform_fuzz_invariants migration sidecar gas softening active_views)
  add_test(NAME validation_transform_fuzz_invariants_long COMMAND test_integration_transform_fuzz_invariants --long)
  set_tests_properties(
    integration_transform_fuzz_invariants
    validation_transform_fuzz_invariants_long
    PROPERTIES RESOURCE_LOCK transform_fuzz_hdf5)
  set_tests_properties(validation_transform_fuzz_invariants_long PROPERTIES LABELS "validation;integration;fuzz")

  add_executable(test_integration_softening_ownership_invariants tests/integration/test_softening_ownership_invariants.cpp)
  target_link_libraries(test_integration_softening_ownership_invariants PRIVATE cosmosim_all)
  add_test(NAME integration_softening_ownership_invariants COMMAND test_integration_softening_ownership_invariants)
  cosmosim_mark_runtime_truth_test(integration_softening_ownership_invariants softening sidecar)


  add_executable(test_integration_profiling_mini_run tests/integration/test_profiling_mini_run.cpp)
  target_link_libraries(test_integration_profiling_mini_run PRIVATE cosmosim_core)
  add_test(NAME integration_profiling_mini_run COMMAND test_integration_profiling_mini_run)

  add_executable(test_integration_time_integration_loop tests/integration/test_time_integration_loop.cpp)
  target_link_libraries(test_integration_time_integration_loop PRIVATE cosmosim_all)
  add_test(NAME integration_time_integration_loop COMMAND test_integration_time_integration_loop)

  add_executable(test_integration_hierarchical_time_bins tests/integration/test_hierarchical_time_bins.cpp)
  target_link_libraries(test_integration_hierarchical_time_bins PRIVATE cosmosim_core)
  add_test(NAME integration_hierarchical_time_bins COMMAND test_integration_hierarchical_time_bins)
  cosmosim_mark_runtime_truth_test(integration_hierarchical_time_bins scheduler)

  add_executable(test_integration_hierarchical_timestep_regression tests/integration/test_hierarchical_timestep_regression.cpp)
  target_link_libraries(test_integration_hierarchical_timestep_regression PRIVATE cosmosim_core)
  add_test(NAME integration_hierarchical_timestep_regression COMMAND test_integration_hierarchical_timestep_regression)
  cosmosim_mark_runtime_truth_test(integration_hierarchical_timestep_regression scheduler)


  add_executable(test_unit_analysis_diagnostics tests/unit/test_analysis_diagnostics.cpp)
  target_link_libraries(test_unit_analysis_diagnostics PRIVATE cosmosim_all)
  add_test(NAME unit_analysis_diagnostics COMMAND test_unit_analysis_diagnostics)
  add_executable(test_unit_halo_workflow tests/unit/test_halo_workflow.cpp)
  target_link_libraries(test_unit_halo_workflow PRIVATE cosmosim_all)
  add_test(NAME unit_halo_workflow COMMAND test_unit_halo_workflow)

  add_executable(test_unit_pm_solver tests/unit/test_pm_solver.cpp)
  target_link_libraries(test_unit_pm_solver PRIVATE cosmosim_gravity)
  add_test(NAME unit_pm_solver COMMAND test_unit_pm_solver)
  set_tests_properties(unit_pm_solver PROPERTIES WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")

  add_executable(test_unit_tree_gravity tests/unit/test_tree_gravity.cpp)
  target_link_libraries(test_unit_tree_gravity PRIVATE cosmosim_gravity)
  target_include_directories(test_unit_tree_gravity PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
  add_test(NAME unit_tree_gravity COMMAND test_unit_tree_gravity)


  add_executable(test_unit_amr_refinement tests/unit/test_amr_refinement.cpp)
  target_link_libraries(test_unit_amr_refinement PRIVATE cosmosim_amr)
  add_test(NAME unit_amr_refinement COMMAND test_unit_amr_refinement)
  add_executable(test_unit_amr_hydro_geometry tests/unit/test_amr_hydro_geometry.cpp)
  target_link_libraries(test_unit_amr_hydro_geometry PRIVATE cosmosim_amr cosmosim_hydro)
  add_test(NAME unit_amr_hydro_geometry COMMAND test_unit_amr_hydro_geometry)
  add_executable(test_unit_amr_ghost_fill tests/unit/test_amr_ghost_fill.cpp)
  target_link_libraries(test_unit_amr_ghost_fill PRIVATE cosmosim_amr cosmosim_hydro)
  add_test(NAME unit_amr_ghost_fill COMMAND test_unit_amr_ghost_fill)
  add_executable(test_unit_amr_temporal_ghost_fill tests/unit/test_amr_temporal_ghost_fill.cpp)
  target_link_libraries(test_unit_amr_temporal_ghost_fill PRIVATE cosmosim_amr cosmosim_hydro cosmosim_core)
  add_custom_target(unit_amr_temporal_ghost_fill DEPENDS test_unit_amr_temporal_ghost_fill)
  add_test(NAME unit_amr_temporal_ghost_fill COMMAND test_unit_amr_temporal_ghost_fill)
  cosmosim_mark_runtime_truth_test(unit_amr_temporal_ghost_fill amr hydro gas unit temporal_ghost)
  add_executable(test_unit_amr_flux_register_generation tests/unit/test_amr_flux_register_generation.cpp)
  target_link_libraries(test_unit_amr_flux_register_generation PRIVATE cosmosim_amr cosmosim_hydro)
  add_test(NAME unit_amr_flux_register_generation COMMAND test_unit_amr_flux_register_generation)
  add_executable(test_unit_amr_pending_flux_registers tests/unit/test_amr_pending_flux_registers.cpp)
  target_link_libraries(test_unit_amr_pending_flux_registers PRIVATE cosmosim_amr cosmosim_hydro cosmosim_core)
  add_custom_target(unit_amr_pending_flux_registers DEPENDS test_unit_amr_pending_flux_registers)
  add_test(NAME unit_amr_pending_flux_registers COMMAND test_unit_amr_pending_flux_registers)
  cosmosim_mark_runtime_truth_test(unit_amr_pending_flux_registers amr hydro gas unit pending_flux)

  add_executable(test_unit_hydro_core_solver tests/unit/test_hydro_core_solver.cpp)
  target_link_libraries(test_unit_hydro_core_solver PRIVATE cosmosim_hydro)
  add_test(NAME unit_hydro_core_solver COMMAND test_unit_hydro_core_solver)


  add_executable(test_integration_parallel_two_rank_restart tests/integration/test_parallel_two_rank_restart.cpp)
  target_link_libraries(test_integration_parallel_two_rank_restart PRIVATE cosmosim_parallel cosmosim_core)
  add_test(NAME integration_parallel_two_rank_restart COMMAND test_integration_parallel_two_rank_restart)
  cosmosim_mark_runtime_truth_test(integration_parallel_two_rank_restart restart scheduler)

  add_executable(test_integration_hydro_sod_like tests/integration/test_hydro_sod_like.cpp)
  target_link_libraries(test_integration_hydro_sod_like PRIVATE cosmosim_hydro)
  add_test(NAME integration_hydro_sod_like COMMAND test_integration_hydro_sod_like)
  add_executable(test_integration_hydro_decoupled_gas_cells tests/integration/test_hydro_decoupled_gas_cells.cpp)
  target_link_libraries(test_integration_hydro_decoupled_gas_cells PRIVATE cosmosim_core cosmosim_hydro)
  add_test(NAME integration_hydro_decoupled_gas_cells COMMAND test_integration_hydro_decoupled_gas_cells)
  cosmosim_mark_runtime_truth_test(integration_hydro_decoupled_gas_cells hydro gas integration)

  add_executable(test_integration_hydro_conservation_periodic tests/integration/test_hydro_conservation_periodic.cpp)
  target_link_libraries(test_integration_hydro_conservation_periodic PRIVATE cosmosim_hydro)
  add_test(NAME integration_hydro_conservation_periodic COMMAND test_integration_hydro_conservation_periodic)

  add_executable(test_integration_hydro_axis_symmetry tests/integration/test_hydro_axis_symmetry.cpp)
  target_link_libraries(test_integration_hydro_axis_symmetry PRIVATE cosmosim_hydro)
  add_test(NAME integration_hydro_axis_symmetry COMMAND test_integration_hydro_axis_symmetry)



  add_executable(test_integration_amr_static_refinement_sync tests/integration/test_amr_static_refinement_sync.cpp)
  target_link_libraries(test_integration_amr_static_refinement_sync PRIVATE cosmosim_amr)
  add_test(NAME integration_amr_static_refinement_sync COMMAND test_integration_amr_static_refinement_sync)
  add_executable(test_integration_amr_patch_boundary_consistency tests/integration/test_amr_patch_boundary_consistency.cpp)
  target_link_libraries(test_integration_amr_patch_boundary_consistency PRIVATE cosmosim_amr cosmosim_hydro)
  add_test(NAME integration_amr_patch_boundary_consistency COMMAND test_integration_amr_patch_boundary_consistency)
  add_executable(test_integration_amr_conservative_refine tests/integration/test_amr_conservative_refine.cpp)
  target_link_libraries(test_integration_amr_conservative_refine PRIVATE cosmosim_amr)
  add_test(NAME integration_amr_conservative_refine COMMAND test_integration_amr_conservative_refine)
  add_executable(test_integration_amr_conservative_derefine tests/integration/test_amr_conservative_derefine.cpp)
  target_link_libraries(test_integration_amr_conservative_derefine PRIVATE cosmosim_amr)
  add_test(NAME integration_amr_conservative_derefine COMMAND test_integration_amr_conservative_derefine)
  add_executable(test_integration_amr_reflux_conservation tests/integration/test_amr_reflux_conservation.cpp)
  target_link_libraries(test_integration_amr_reflux_conservation PRIVATE cosmosim_amr cosmosim_hydro cosmosim_core)
  add_test(NAME integration_amr_reflux_conservation COMMAND test_integration_amr_reflux_conservation)
  cosmosim_mark_runtime_truth_test(integration_amr_reflux_conservation amr hydro gas integration)
  add_executable(test_integration_amr_production_hydro_integration tests/integration/test_amr_production_hydro_integration.cpp)
  target_link_libraries(test_integration_amr_production_hydro_integration PRIVATE cosmosim_amr cosmosim_hydro cosmosim_core)
  add_test(NAME integration_amr_production_hydro_integration COMMAND test_integration_amr_production_hydro_integration)
  cosmosim_mark_runtime_truth_test(integration_amr_production_hydro_integration amr hydro gas integration)
  add_executable(test_integration_amr_hydro_subcycling tests/integration/test_amr_hydro_subcycling.cpp)
  target_link_libraries(test_integration_amr_hydro_subcycling PRIVATE cosmosim_amr cosmosim_hydro cosmosim_core)
  add_custom_target(integration_amr_hydro_subcycling DEPENDS test_integration_amr_hydro_subcycling)
  add_test(NAME integration_amr_hydro_subcycling COMMAND test_integration_amr_hydro_subcycling)
  cosmosim_mark_runtime_truth_test(integration_amr_hydro_subcycling amr hydro gas integration pending_flux)
  add_executable(test_integration_amr_distributed_remote_patch_boundary tests/integration/test_amr_distributed_remote_patch_boundary.cpp)
  target_link_libraries(test_integration_amr_distributed_remote_patch_boundary PRIVATE cosmosim_amr cosmosim_hydro cosmosim_core)
  add_custom_target(integration_amr_distributed_remote_patch_boundary DEPENDS test_integration_amr_distributed_remote_patch_boundary)
  add_test(NAME integration_amr_distributed_remote_patch_boundary COMMAND test_integration_amr_distributed_remote_patch_boundary)
  cosmosim_mark_runtime_truth_test(integration_amr_distributed_remote_patch_boundary amr hydro gas integration reflux serial)
  add_executable(test_integration_amr_hydro_shock_tube tests/integration/test_amr_hydro_shock_tube.cpp)
  target_link_libraries(test_integration_amr_hydro_shock_tube PRIVATE cosmosim_amr cosmosim_hydro cosmosim_core)
  add_custom_target(integration_amr_hydro_shock_tube DEPENDS test_integration_amr_hydro_shock_tube)
  add_test(NAME integration_amr_hydro_shock_tube COMMAND test_integration_amr_hydro_shock_tube)
  cosmosim_mark_runtime_truth_test(integration_amr_hydro_shock_tube amr hydro gas integration regression smoke)
  add_executable(test_integration_amr_hydro_sedov tests/integration/test_amr_hydro_sedov.cpp)
  target_link_libraries(test_integration_amr_hydro_sedov PRIVATE cosmosim_amr cosmosim_hydro cosmosim_core)
  add_custom_target(integration_amr_hydro_sedov DEPENDS test_integration_amr_hydro_sedov)
  add_test(NAME integration_amr_hydro_sedov COMMAND test_integration_amr_hydro_sedov)
  cosmosim_mark_runtime_truth_test(integration_amr_hydro_sedov amr hydro gas integration regression smoke)
  add_executable(test_integration_amr_synchronization_stress tests/integration/test_amr_synchronization_stress.cpp)
  target_link_libraries(test_integration_amr_synchronization_stress PRIVATE cosmosim_amr cosmosim_hydro cosmosim_core)
  add_custom_target(integration_amr_synchronization_stress DEPENDS test_integration_amr_synchronization_stress)
  add_test(NAME integration_amr_synchronization_stress COMMAND test_integration_amr_synchronization_stress)
  cosmosim_mark_runtime_truth_test(integration_amr_synchronization_stress amr hydro gas integration regression smoke)

  add_executable(test_integration_simulation_mode_toy_runs tests/integration/test_simulation_mode_toy_runs.cpp)
  target_link_libraries(test_integration_simulation_mode_toy_runs PRIVATE cosmosim_core)
  add_test(NAME integration_simulation_mode_toy_runs COMMAND test_integration_simulation_mode_toy_runs)

  add_executable(test_integration_pm_periodic_mode tests/integration/test_pm_periodic_mode.cpp)
  target_link_libraries(test_integration_pm_periodic_mode PRIVATE cosmosim_gravity)
  add_test(NAME integration_pm_periodic_mode COMMAND test_integration_pm_periodic_mode)

  add_executable(test_integration_pm_scale_factor_regression tests/integration/test_pm_scale_factor_regression.cpp)
  target_link_libraries(test_integration_pm_scale_factor_regression PRIVATE cosmosim_gravity)
  add_test(NAME integration_pm_scale_factor_regression COMMAND test_integration_pm_scale_factor_regression)
  if(COSMOSIM_ENABLE_MPI)
    add_test(
      NAME integration_pm_periodic_mode_mpi_two_rank
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
              $<TARGET_FILE:test_integration_pm_periodic_mode> ${MPIEXEC_POSTFLAGS})
    cosmosim_append_test_labels(integration_pm_periodic_mode_mpi_two_rank mpi distributed pm)
    if(COSMOSIM_ENABLE_FFTW)
      cosmosim_append_test_labels(integration_pm_periodic_mode_mpi_two_rank fftw)
    endif()
    set_tests_properties(integration_pm_periodic_mode_mpi_two_rank PROPERTIES PROCESSORS 2)
    add_test(
      NAME integration_pm_zero_width_slab_mpi_four_rank
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
              $<TARGET_FILE:test_integration_pm_periodic_mode> ${MPIEXEC_POSTFLAGS})
    cosmosim_append_test_labels(integration_pm_zero_width_slab_mpi_four_rank mpi distributed pm)
    if(COSMOSIM_ENABLE_FFTW)
      cosmosim_append_test_labels(integration_pm_zero_width_slab_mpi_four_rank fftw)
    endif()
    set_tests_properties(integration_pm_zero_width_slab_mpi_four_rank PROPERTIES PROCESSORS 4)
  endif()

  if(COSMOSIM_ENABLE_MPI)
    add_executable(test_integration_mpi_collective_safety
      tests/integration/test_mpi_collective_safety.cpp)
    target_link_libraries(test_integration_mpi_collective_safety PRIVATE cosmosim_parallel)
    add_test(
      NAME integration_mpi_collective_safety_two_rank
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
              $<TARGET_FILE:test_integration_mpi_collective_safety> ${MPIEXEC_POSTFLAGS})
    add_test(
      NAME integration_mpi_collective_safety_three_rank
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
              $<TARGET_FILE:test_integration_mpi_collective_safety> ${MPIEXEC_POSTFLAGS})
    set_tests_properties(
      integration_mpi_collective_safety_two_rank
      integration_mpi_collective_safety_three_rank
      PROPERTIES LABELS "mpi;distributed;correctness;collective_safety")
    set_tests_properties(integration_mpi_collective_safety_two_rank PROPERTIES PROCESSORS 2)
    set_tests_properties(integration_mpi_collective_safety_three_rank PROPERTIES PROCESSORS 3)

    add_executable(test_integration_pm_slab_halo_exchange_mpi tests/integration/test_pm_slab_halo_exchange_mpi.cpp)
    target_link_libraries(test_integration_pm_slab_halo_exchange_mpi PRIVATE cosmosim_parallel)
    add_test(
      NAME integration_pm_slab_halo_exchange_mpi_two_rank
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
              $<TARGET_FILE:test_integration_pm_slab_halo_exchange_mpi> ${MPIEXEC_POSTFLAGS})
    add_test(
      NAME integration_pm_slab_halo_exchange_mpi_three_rank
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
              $<TARGET_FILE:test_integration_pm_slab_halo_exchange_mpi> ${MPIEXEC_POSTFLAGS})
    add_test(
      NAME integration_pm_slab_halo_exchange_mpi_four_rank
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
              $<TARGET_FILE:test_integration_pm_slab_halo_exchange_mpi> ${MPIEXEC_POSTFLAGS})
    add_test(
      NAME integration_pm_slab_halo_exchange_mpi_eight_rank
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 8 ${MPIEXEC_PREFLAGS}
              $<TARGET_FILE:test_integration_pm_slab_halo_exchange_mpi> ${MPIEXEC_POSTFLAGS})
    cosmosim_append_test_labels(integration_pm_slab_halo_exchange_mpi_two_rank mpi distributed pm)
    cosmosim_append_test_labels(integration_pm_slab_halo_exchange_mpi_three_rank mpi distributed pm)
    cosmosim_append_test_labels(integration_pm_slab_halo_exchange_mpi_four_rank mpi distributed pm)
    cosmosim_append_test_labels(integration_pm_slab_halo_exchange_mpi_eight_rank mpi distributed pm)
    set_tests_properties(integration_pm_slab_halo_exchange_mpi_two_rank PROPERTIES PROCESSORS 2)
    set_tests_properties(integration_pm_slab_halo_exchange_mpi_three_rank PROPERTIES PROCESSORS 3)
    set_tests_properties(integration_pm_slab_halo_exchange_mpi_four_rank PROPERTIES PROCESSORS 4)
    set_tests_properties(integration_pm_slab_halo_exchange_mpi_eight_rank PROPERTIES PROCESSORS 8)
    add_executable(test_integration_distributed_sfc_rebalance_mpi
      tests/integration/test_distributed_sfc_rebalance_mpi.cpp)
    target_link_libraries(test_integration_distributed_sfc_rebalance_mpi PRIVATE cosmosim_parallel)
    add_test(
      NAME integration_distributed_sfc_rebalance_mpi_two_rank
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
              $<TARGET_FILE:test_integration_distributed_sfc_rebalance_mpi> ${MPIEXEC_POSTFLAGS})
    add_test(
      NAME integration_distributed_sfc_rebalance_mpi_three_rank
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
              $<TARGET_FILE:test_integration_distributed_sfc_rebalance_mpi> ${MPIEXEC_POSTFLAGS})
    add_test(
      NAME integration_distributed_sfc_rebalance_mpi_four_rank
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
              $<TARGET_FILE:test_integration_distributed_sfc_rebalance_mpi> ${MPIEXEC_POSTFLAGS})
    set_tests_properties(
      integration_distributed_sfc_rebalance_mpi_two_rank
      integration_distributed_sfc_rebalance_mpi_three_rank
      integration_distributed_sfc_rebalance_mpi_four_rank
      PROPERTIES LABELS "mpi;distributed;decomposition;runtime_truth")
    set_tests_properties(integration_distributed_sfc_rebalance_mpi_four_rank PROPERTIES PROCESSORS 4)
  endif()

  if(COSMOSIM_ENABLE_MPI AND COSMOSIM_ENABLE_HDF5)
    add_executable(test_integration_snapshot_set_mpi
      tests/integration/test_snapshot_set_mpi.cpp)
    target_link_libraries(test_integration_snapshot_set_mpi PRIVATE cosmosim_all)
    foreach(snapshot_rank_count IN ITEMS 2 3 4 8)
      add_test(
        NAME integration_snapshot_set_mpi_${snapshot_rank_count}_rank
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${snapshot_rank_count} ${MPIEXEC_PREFLAGS}
                $<TARGET_FILE:test_integration_snapshot_set_mpi> ${MPIEXEC_POSTFLAGS})
      set_tests_properties(
        integration_snapshot_set_mpi_${snapshot_rank_count}_rank
        PROPERTIES PROCESSORS ${snapshot_rank_count} TIMEOUT 90
                   LABELS "mpi;hdf5;distributed;snapshot;archive;runtime_truth")
    endforeach()

    add_executable(test_integration_distributed_ic_reader_mpi
      tests/integration/test_distributed_ic_reader_mpi.cpp)
    target_link_libraries(test_integration_distributed_ic_reader_mpi PRIVATE cosmosim_all)
    foreach(rank_count IN ITEMS 1 2 3 4 8)
      add_test(
        NAME integration_distributed_ic_reader_mpi_${rank_count}_rank
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${rank_count} ${MPIEXEC_PREFLAGS}
                $<TARGET_FILE:test_integration_distributed_ic_reader_mpi> normal ${MPIEXEC_POSTFLAGS})
      set_tests_properties(integration_distributed_ic_reader_mpi_${rank_count}_rank
        PROPERTIES PROCESSORS ${rank_count} TIMEOUT 60
                   LABELS "mpi;hdf5;distributed;ic;mixed_species;runtime_truth")
    endforeach()
    foreach(acceptance_mode IN ITEMS
        canonical
        manifest
        canonical_manifest
        alternate
        type_policy)
      foreach(rank_count IN ITEMS 1 2 3 4 8)
        add_test(
          NAME integration_distributed_ic_reader_${acceptance_mode}_mpi_${rank_count}_rank
          COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${rank_count} ${MPIEXEC_PREFLAGS}
                  $<TARGET_FILE:test_integration_distributed_ic_reader_mpi>
                  ${acceptance_mode} ${MPIEXEC_POSTFLAGS})
        set_tests_properties(
          integration_distributed_ic_reader_${acceptance_mode}_mpi_${rank_count}_rank
          PROPERTIES PROCESSORS ${rank_count} TIMEOUT 90
                     LABELS "mpi;hdf5;distributed;ic;equivalence;runtime_truth")
      endforeach()
    endforeach()
    foreach(manifest_failure_mode IN ITEMS
        canonical_manifest_tampered_sidecar
        canonical_manifest_missing_marker
        canonical_manifest_bad_digest
        canonical_manifest_source_hash_mismatch
        canonical_manifest_source_path_mismatch
        canonical_manifest_bridge_dialect
        bridge_manifest_canonical_dialect
        invalid_manifest)
      add_test(
        NAME integration_distributed_ic_reader_${manifest_failure_mode}_mpi_two_rank
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
                $<TARGET_FILE:test_integration_distributed_ic_reader_mpi>
                ${manifest_failure_mode} ${MPIEXEC_POSTFLAGS})
      set_tests_properties(
        integration_distributed_ic_reader_${manifest_failure_mode}_mpi_two_rank
        PROPERTIES PROCESSORS 2 TIMEOUT 30
                   LABELS "mpi;hdf5;distributed;ic;manifest;malformed;runtime_truth")
    endforeach()

    foreach(rank_count IN ITEMS 2 4 8)
      add_test(
        NAME integration_distributed_ic_reader_scaling_mpi_${rank_count}_rank
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${rank_count} ${MPIEXEC_PREFLAGS}
                $<TARGET_FILE:test_integration_distributed_ic_reader_mpi> scaling ${MPIEXEC_POSTFLAGS})
      set_tests_properties(
        integration_distributed_ic_reader_scaling_mpi_${rank_count}_rank
        PROPERTIES PROCESSORS ${rank_count} TIMEOUT 90
                   LABELS "mpi;hdf5;distributed;ic;scaling;memory;runtime_truth")
    endforeach()
    add_test(
      NAME integration_distributed_ic_reader_duplicate_mpi_two_rank
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
              $<TARGET_FILE:test_integration_distributed_ic_reader_mpi> duplicate ${MPIEXEC_POSTFLAGS})
    set_tests_properties(integration_distributed_ic_reader_duplicate_mpi_two_rank
      PROPERTIES PROCESSORS 2 TIMEOUT 60
                 LABELS "mpi;hdf5;distributed;ic;malformed;runtime_truth")
    foreach(route_case IN ITEMS route_loss route_duplicate)
      add_test(
        NAME integration_distributed_ic_reader_${route_case}_mpi_two_rank
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
                $<TARGET_FILE:test_integration_distributed_ic_reader_mpi>
                ${route_case} ${MPIEXEC_POSTFLAGS})
      set_tests_properties(
        integration_distributed_ic_reader_${route_case}_mpi_two_rank
        PROPERTIES PROCESSORS 2 TIMEOUT 30
                   LABELS "mpi;hdf5;distributed;ic;reconciliation;malformed;runtime_truth")
    endforeach()
    foreach(fault_phase IN ITEMS
        manifest_digest
        reader_session_creation
        owner_serialization
        serialization_accounting
        send_layout
        post_exchange_accounting
        payload_validation
        deserialization
        source_final_reconciliation
        source_final_accounting
        global_id_audit_accounting
        exact_audit_capacity_accounting
        duplicate_audit_capacity_accounting
        audit_path_storage_accounting
        audit_peak_staging_accounting
        sidecar_append
        source_identity_revalidation
        final_state)
      add_test(
        NAME integration_distributed_ic_reader_fault_${fault_phase}_mpi_two_rank
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
                $<TARGET_FILE:test_integration_distributed_ic_reader_mpi>
                fault_${fault_phase} ${MPIEXEC_POSTFLAGS})
      set_tests_properties(
        integration_distributed_ic_reader_fault_${fault_phase}_mpi_two_rank
        PROPERTIES PROCESSORS 2 TIMEOUT 30
                   LABELS "mpi;hdf5;distributed;ic;fault_injection;malformed;runtime_truth")
    endforeach()
  endif()

  if(COSMOSIM_ENABLE_MPI AND COSMOSIM_ENABLE_HDF5)
    add_executable(test_integration_reference_workflow_distributed_amr_mpi
      tests/integration/test_reference_workflow_distributed_amr_mpi.cpp)
    target_link_libraries(test_integration_reference_workflow_distributed_amr_mpi PRIVATE cosmosim_all)
    add_test(
      NAME integration_reference_workflow_distributed_amr_mpi_two_rank
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
              $<TARGET_FILE:test_integration_reference_workflow_distributed_amr_mpi> ${MPIEXEC_POSTFLAGS})
    cosmosim_mark_runtime_truth_test(
      integration_reference_workflow_distributed_amr_mpi_two_rank
      mpi hdf5 distributed amr hydro gas workflow directed_exchange reflux integration)
  endif()

  if(COSMOSIM_ENABLE_MPI AND COSMOSIM_ENABLE_HDF5)
    add_executable(test_integration_distributed_gas_cell_migration_mpi
      tests/integration/test_distributed_gas_cell_migration_mpi.cpp)
    target_link_libraries(test_integration_distributed_gas_cell_migration_mpi PRIVATE cosmosim_all)
    add_test(
      NAME integration_distributed_gas_cell_migration_mpi_two_rank
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
              $<TARGET_FILE:test_integration_distributed_gas_cell_migration_mpi> ${MPIEXEC_POSTFLAGS})
    cosmosim_mark_runtime_truth_test(integration_distributed_gas_cell_migration_mpi_two_rank mpi hdf5 gas migration restart integration)
  endif()

  if(COSMOSIM_ENABLE_MPI AND COSMOSIM_ENABLE_HDF5)
    add_executable(test_integration_distributed_amr_patch_migration_mpi
      tests/integration/test_distributed_amr_patch_migration_mpi.cpp)
    target_include_directories(
      test_integration_distributed_amr_patch_migration_mpi PRIVATE ${PROJECT_SOURCE_DIR}/src)
    target_link_libraries(test_integration_distributed_amr_patch_migration_mpi PRIVATE cosmosim_all)
    add_test(
      NAME integration_distributed_amr_patch_migration_mpi_two_rank
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
              $<TARGET_FILE:test_integration_distributed_amr_patch_migration_mpi> ${MPIEXEC_POSTFLAGS})
    set_tests_properties(
      integration_distributed_amr_patch_migration_mpi_two_rank
      PROPERTIES PROCESSORS 2 TIMEOUT 90)
    cosmosim_mark_runtime_truth_test(
      integration_distributed_amr_patch_migration_mpi_two_rank
      mpi hdf5 distributed amr hydro gas migration restart reflux temporal_history integration)
  endif()

  if(COSMOSIM_ENABLE_CUDA)
    add_executable(test_integration_pm_cuda_smoke tests/integration/test_pm_cuda_smoke.cpp)
    target_link_libraries(test_integration_pm_cuda_smoke PRIVATE cosmosim_gravity)
    add_test(NAME integration_pm_cuda_smoke COMMAND test_integration_pm_cuda_smoke)
  endif()


  if(COSMOSIM_ENABLE_HDF5)
    add_executable(cosmosim_convert_ic tools/convert_ic.cpp)
    target_include_directories(
      cosmosim_convert_ic PRIVATE ${PROJECT_SOURCE_DIR}/src)
    target_link_libraries(cosmosim_convert_ic PRIVATE cosmosim_io)
    find_package(Python3 QUIET COMPONENTS Interpreter)
    if(Python3_Interpreter_FOUND)
      add_test(
        NAME integration_convert_ic
        COMMAND ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/tests/integration/io/test_convert_ic.py
                $<TARGET_FILE:cosmosim_convert_ic>)
      set_tests_properties(integration_convert_ic PROPERTIES
        SKIP_RETURN_CODE 77 LABELS "integration;hdf5;ic;converter;runtime_truth")
    endif()
  endif()

  if(COSMOSIM_ENABLE_HDF5)
    add_executable(test_integration_initial_condition_manifest_runtime
      tests/integration/test_initial_condition_manifest_runtime.cpp)
    target_include_directories(
      test_integration_initial_condition_manifest_runtime
      PRIVATE ${PROJECT_SOURCE_DIR}/src)
    target_link_libraries(
      test_integration_initial_condition_manifest_runtime PRIVATE cosmosim_all)
    add_test(
      NAME integration_initial_condition_manifest_runtime
      COMMAND test_integration_initial_condition_manifest_runtime)
    cosmosim_mark_runtime_truth_test(
      integration_initial_condition_manifest_runtime
      integration hdf5 ic manifest runtime_truth)
  endif()

  add_executable(test_unit_transactional_file tests/unit/test_transactional_file.cpp)
  target_include_directories(test_unit_transactional_file PRIVATE ${PROJECT_SOURCE_DIR}/src)
  target_link_libraries(test_unit_transactional_file PRIVATE cosmosim_io)
  add_test(NAME unit_transactional_file COMMAND test_unit_transactional_file)
  cosmosim_mark_runtime_truth_test(unit_transactional_file io transaction durability)

  add_executable(test_unit_ic_reader tests/unit/test_ic_reader.cpp)
  target_include_directories(
    test_unit_ic_reader PRIVATE ${PROJECT_SOURCE_DIR}/src)
  target_link_libraries(test_unit_ic_reader PRIVATE cosmosim_io)
  add_test(NAME unit_ic_reader COMMAND test_unit_ic_reader)

  add_executable(test_unit_ic_record_codec tests/unit/test_ic_record_codec.cpp)
  target_include_directories(
    test_unit_ic_record_codec PRIVATE ${PROJECT_SOURCE_DIR}/src)
  target_link_libraries(test_unit_ic_record_codec PRIVATE cosmosim_io)
  add_test(NAME unit_ic_record_codec COMMAND test_unit_ic_record_codec)

  add_executable(test_unit_snapshot_hdf5_schema tests/unit/test_snapshot_hdf5_schema.cpp)
  target_link_libraries(test_unit_snapshot_hdf5_schema PRIVATE cosmosim_io)
  add_test(NAME unit_snapshot_hdf5_schema COMMAND test_unit_snapshot_hdf5_schema)
  cosmosim_mark_runtime_truth_test(unit_snapshot_hdf5_schema snapshot)


  add_executable(test_unit_restart_checkpoint_schema tests/unit/test_restart_checkpoint_schema.cpp)
  target_link_libraries(test_unit_restart_checkpoint_schema PRIVATE cosmosim_io)
  add_test(NAME unit_restart_checkpoint_schema COMMAND test_unit_restart_checkpoint_schema)
  cosmosim_mark_runtime_truth_test(unit_restart_checkpoint_schema restart provenance)

  add_executable(test_unit_tree_pm_split_kernel tests/unit/test_tree_pm_split_kernel.cpp)
  target_include_directories(test_unit_tree_pm_split_kernel PRIVATE ${PROJECT_SOURCE_DIR}/src)
  target_link_libraries(test_unit_tree_pm_split_kernel PRIVATE cosmosim_gravity)
  add_test(NAME unit_tree_pm_split_kernel COMMAND test_unit_tree_pm_split_kernel)

  add_executable(test_integration_tree_pm_coupling_periodic tests/integration/test_tree_pm_coupling_periodic.cpp)
  target_link_libraries(test_integration_tree_pm_coupling_periodic PRIVATE cosmosim_gravity)
  add_test(NAME integration_tree_pm_coupling_periodic COMMAND test_integration_tree_pm_coupling_periodic)
  if(COSMOSIM_ENABLE_MPI)
    add_test(
      NAME integration_tree_pm_coupling_periodic_mpi_two_rank
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
              $<TARGET_FILE:test_integration_tree_pm_coupling_periodic> ${MPIEXEC_POSTFLAGS})
    add_test(
      NAME integration_tree_pm_coupling_periodic_mpi_three_rank
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
              $<TARGET_FILE:test_integration_tree_pm_coupling_periodic> ${MPIEXEC_POSTFLAGS})
    add_test(
      NAME integration_tree_pm_coupling_periodic_mpi_four_rank
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
              $<TARGET_FILE:test_integration_tree_pm_coupling_periodic> ${MPIEXEC_POSTFLAGS})
    add_test(
      NAME integration_tree_pm_coupling_periodic_mpi_eight_rank
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 8 ${MPIEXEC_PREFLAGS}
              $<TARGET_FILE:test_integration_tree_pm_coupling_periodic> ${MPIEXEC_POSTFLAGS})
    set_tests_properties(
      integration_tree_pm_coupling_periodic_mpi_two_rank
      integration_tree_pm_coupling_periodic_mpi_three_rank
      integration_tree_pm_coupling_periodic_mpi_four_rank
      integration_tree_pm_coupling_periodic_mpi_eight_rank
      PROPERTIES LABELS "mpi;distributed;treepm;pm")
    if(COSMOSIM_ENABLE_FFTW)
      foreach(cosmosim_treepm_mpi_test IN ITEMS
          integration_tree_pm_coupling_periodic_mpi_two_rank
          integration_tree_pm_coupling_periodic_mpi_three_rank
          integration_tree_pm_coupling_periodic_mpi_four_rank
          integration_tree_pm_coupling_periodic_mpi_eight_rank)
        cosmosim_append_test_labels(${cosmosim_treepm_mpi_test} fftw)
      endforeach()
    endif()
    set_tests_properties(integration_tree_pm_coupling_periodic_mpi_two_rank PROPERTIES PROCESSORS 2)
    set_tests_properties(integration_tree_pm_coupling_periodic_mpi_three_rank PROPERTIES PROCESSORS 3)
    set_tests_properties(integration_tree_pm_coupling_periodic_mpi_four_rank PROPERTIES PROCESSORS 4)
    set_tests_properties(integration_tree_pm_coupling_periodic_mpi_eight_rank PROPERTIES PROCESSORS 8)
  endif()


  add_executable(test_integration_ic_generated_pipeline tests/integration/test_ic_generated_pipeline.cpp)
  target_link_libraries(test_integration_ic_generated_pipeline PRIVATE cosmosim_all)
  add_test(NAME integration_ic_generated_pipeline COMMAND test_integration_ic_generated_pipeline)

  if(COSMOSIM_ENABLE_HDF5)
  add_executable(test_integration_snapshot_hdf5_roundtrip tests/integration/test_snapshot_hdf5_roundtrip.cpp)
  target_link_libraries(test_integration_snapshot_hdf5_roundtrip PRIVATE cosmosim_all)
  add_test(NAME integration_snapshot_hdf5_roundtrip COMMAND test_integration_snapshot_hdf5_roundtrip)
  cosmosim_mark_runtime_truth_test(integration_snapshot_hdf5_roundtrip snapshot sidecar)
  endif()
  add_executable(test_integration_tracer_advection tests/integration/test_tracer_advection.cpp)
  target_link_libraries(test_integration_tracer_advection PRIVATE cosmosim_all)
  add_test(NAME integration_tracer_advection COMMAND test_integration_tracer_advection)


  if(COSMOSIM_ENABLE_HDF5)
  add_executable(test_integration_restart_checkpoint_roundtrip tests/integration/test_restart_checkpoint_roundtrip.cpp)
  target_link_libraries(test_integration_restart_checkpoint_roundtrip PRIVATE cosmosim_all)
  add_test(NAME integration_restart_checkpoint_roundtrip COMMAND test_integration_restart_checkpoint_roundtrip)
  cosmosim_mark_runtime_truth_test(integration_restart_checkpoint_roundtrip restart provenance gas sidecar scheduler)

  add_executable(test_integration_restart_equivalence_harness tests/integration/test_restart_equivalence_harness.cpp)
  target_link_libraries(test_integration_restart_equivalence_harness PRIVATE cosmosim_all)
  add_test(NAME integration_restart_equivalence_harness COMMAND test_integration_restart_equivalence_harness)
  cosmosim_mark_runtime_truth_test(integration_restart_equivalence_harness restart scheduler output stochastic equivalence integration)

  foreach(stage8_case IN ITEMS dm_only treepm hydro_toy multirate_bins output_enabled stochastic_sources amr_hydro amr_flux_registers)
    add_executable(
      test_integration_restart_equivalence_${stage8_case}
      tests/integration/test_restart_equivalence_${stage8_case}.cpp)
    target_link_libraries(test_integration_restart_equivalence_${stage8_case} PRIVATE cosmosim_all)
    add_test(
      NAME integration_restart_equivalence_${stage8_case}
      COMMAND test_integration_restart_equivalence_${stage8_case})
  endforeach()
  add_custom_target(integration_restart_equivalence_amr_hydro DEPENDS test_integration_restart_equivalence_amr_hydro)
  add_custom_target(integration_restart_equivalence_amr_flux_registers DEPENDS test_integration_restart_equivalence_amr_flux_registers)
  add_executable(
    test_integration_restart_equivalence_amr_temporal_ghosts
    tests/integration/test_restart_equivalence_amr_temporal_ghosts.cpp)
  target_link_libraries(test_integration_restart_equivalence_amr_temporal_ghosts PRIVATE cosmosim_all)
  add_custom_target(
    integration_restart_equivalence_amr_temporal_ghosts
    DEPENDS test_integration_restart_equivalence_amr_temporal_ghosts)
  add_test(
    NAME integration_restart_equivalence_amr_temporal_ghosts
    COMMAND test_integration_restart_equivalence_amr_temporal_ghosts)
  cosmosim_mark_runtime_truth_test(integration_restart_equivalence_dm_only restart scheduler equivalence integration)
  cosmosim_mark_runtime_truth_test(integration_restart_equivalence_treepm restart scheduler treepm pm equivalence integration)
  cosmosim_mark_runtime_truth_test(integration_restart_equivalence_hydro_toy restart scheduler hydro gas equivalence integration)
  cosmosim_mark_runtime_truth_test(integration_restart_equivalence_multirate_bins restart scheduler multirate equivalence integration)
  cosmosim_mark_runtime_truth_test(integration_restart_equivalence_output_enabled restart scheduler output equivalence integration)
  cosmosim_mark_runtime_truth_test(integration_restart_equivalence_stochastic_sources restart scheduler stochastic source_terms equivalence integration)
  cosmosim_mark_runtime_truth_test(integration_restart_equivalence_amr_hydro restart scheduler hydro gas amr equivalence integration)
  cosmosim_mark_runtime_truth_test(integration_restart_equivalence_amr_flux_registers restart scheduler hydro gas amr equivalence integration pending_flux)
  cosmosim_mark_runtime_truth_test(integration_restart_equivalence_amr_temporal_ghosts restart scheduler hydro gas amr equivalence integration pending_flux temporal_ghost)

  endif()
  add_executable(test_integration_tree_gravity_vs_direct tests/integration/test_tree_gravity_vs_direct.cpp)
  target_link_libraries(test_integration_tree_gravity_vs_direct PRIVATE cosmosim_gravity)
  target_include_directories(test_integration_tree_gravity_vs_direct PRIVATE tests/support)
  add_test(NAME integration_tree_gravity_vs_direct COMMAND test_integration_tree_gravity_vs_direct)
  add_executable(test_integration_config_examples tests/integration/test_config_examples.cpp)
  target_link_libraries(test_integration_config_examples PRIVATE cosmosim_all)
  target_compile_definitions(test_integration_config_examples PRIVATE COSMOSIM_SOURCE_DIR="${PROJECT_SOURCE_DIR}")
  add_test(NAME integration_config_examples COMMAND test_integration_config_examples)
  add_test(
    NAME integration_config_alias_consistency
    COMMAND ${CMAKE_COMMAND}
      -DCOSMOSIM_SOURCE_DIR=${PROJECT_SOURCE_DIR}
      -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/check_config_aliases.cmake)
  add_executable(test_integration_provenance_roundtrip tests/integration/test_provenance_roundtrip.cpp)
  target_link_libraries(test_integration_provenance_roundtrip PRIVATE cosmosim_core)
  add_test(NAME integration_provenance_roundtrip COMMAND test_integration_provenance_roundtrip)
  cosmosim_mark_runtime_truth_test(integration_provenance_roundtrip provenance)

  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/tests/integration/test_feature_summary.cmake.in"
    "${CMAKE_BINARY_DIR}/test_feature_summary.cmake"
    @ONLY)
  add_test(NAME integration_feature_summary
           COMMAND ${CMAKE_COMMAND} -P "${CMAKE_BINARY_DIR}/test_feature_summary.cmake")

  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/tests/integration/test_configure_presets.cmake.in"
    "${CMAKE_BINARY_DIR}/test_configure_presets.cmake"
    @ONLY)
  add_test(NAME integration_configure_presets
           COMMAND ${CMAKE_COMMAND} -P "${CMAKE_BINARY_DIR}/test_configure_presets.cmake")

  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/tests/integration/test_docs_scaffold.cmake.in"
    "${CMAKE_BINARY_DIR}/test_docs_scaffold.cmake"
    @ONLY)
  add_test(NAME integration_docs_scaffold
           COMMAND ${CMAKE_COMMAND} -P "${CMAKE_BINARY_DIR}/test_docs_scaffold.cmake")

  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/tests/integration/test_core_dependency_direction.cmake.in"
    "${CMAKE_BINARY_DIR}/test_core_dependency_direction.cmake"
    @ONLY)
  add_test(NAME integration_core_dependency_direction
           COMMAND ${CMAKE_COMMAND} -P "${CMAKE_BINARY_DIR}/test_core_dependency_direction.cmake")

  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/tests/integration/test_ic_decomposition_structure.cmake.in"
    "${CMAKE_BINARY_DIR}/test_ic_decomposition_structure.cmake"
    @ONLY)
  add_test(NAME integration_ic_decomposition_structure
           COMMAND ${CMAKE_COMMAND} -P "${CMAKE_BINARY_DIR}/test_ic_decomposition_structure.cmake")
  set_tests_properties(integration_ic_decomposition_structure
    PROPERTIES LABELS "architecture;contract;ic")

  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/tests/integration/test_runtime_decomposition_structure.cmake.in"
    "${CMAKE_BINARY_DIR}/test_runtime_decomposition_structure.cmake"
    @ONLY)
  add_test(NAME integration_runtime_decomposition_structure
           COMMAND ${CMAKE_COMMAND} -P "${CMAKE_BINARY_DIR}/test_runtime_decomposition_structure.cmake")
  set_tests_properties(
    integration_runtime_decomposition_structure
    PROPERTIES LABELS "architecture;contract;workflow;runtime_registry")

  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/tests/integration/test_runtime_resource_capability_compile.cmake.in"
    "${CMAKE_BINARY_DIR}/test_runtime_resource_capability_compile.cmake"
    @ONLY)
  add_test(NAME integration_runtime_resource_capability_compile
           COMMAND ${CMAKE_COMMAND} -P "${CMAKE_BINARY_DIR}/test_runtime_resource_capability_compile.cmake")
  set_tests_properties(
    integration_runtime_resource_capability_compile
    PROPERTIES LABELS "architecture;contract;compile_only;workflow;runtime_registry")

  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/tests/integration/test_release_readiness_artifacts.cmake.in"
    "${CMAKE_BINARY_DIR}/test_release_readiness_artifacts.cmake"
    @ONLY)
  add_test(NAME integration_release_readiness_artifacts
           COMMAND ${CMAKE_COMMAND} -P "${CMAKE_BINARY_DIR}/test_release_readiness_artifacts.cmake")

  if(UNIX)
    add_test(NAME integration_source_package_completeness
             COMMAND bash "${CMAKE_CURRENT_SOURCE_DIR}/scripts/ci/test_source_package_completeness.sh")
    set_tests_properties(integration_source_package_completeness PROPERTIES
      LABELS "integration;release;packaging;source_distribution"
      RUN_SERIAL TRUE
      TIMEOUT 1800)
  endif()

  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/tests/integration/test_runtime_truth_ctest_labels.cmake.in"
    "${CMAKE_BINARY_DIR}/test_runtime_truth_ctest_labels.cmake"
    @ONLY)
  add_test(NAME integration_runtime_truth_ctest_labels
           COMMAND ${CMAKE_COMMAND} -P "${CMAKE_BINARY_DIR}/test_runtime_truth_ctest_labels.cmake")
  cosmosim_mark_runtime_truth_test(integration_runtime_truth_ctest_labels ci labels)

  set(cosmosim_stage1_runtime_truth_targets
      test_unit_simulation_state
      test_unit_gas_cell_identity_invariants
      test_unit_hot_cold_sidecar_layout
      test_integration_reorder_compaction_sidecars
      test_integration_species_migration_invariants
      test_integration_gas_cell_identity_migration
      test_integration_transform_fuzz_invariants
      test_integration_softening_ownership_invariants
      test_integration_hierarchical_time_bins
      test_integration_hierarchical_timestep_regression
      test_integration_parallel_two_rank_restart
      test_unit_snapshot_hdf5_schema
      test_unit_restart_checkpoint_schema
      test_integration_provenance_roundtrip)
  if(COSMOSIM_ENABLE_HDF5)
    list(APPEND cosmosim_stage1_runtime_truth_targets
      test_integration_snapshot_hdf5_roundtrip
      test_integration_restart_checkpoint_roundtrip
      test_integration_restart_equivalence_harness
      test_integration_restart_equivalence_dm_only
      test_integration_restart_equivalence_treepm
      test_integration_restart_equivalence_hydro_toy
      test_integration_restart_equivalence_multirate_bins
      test_integration_restart_equivalence_output_enabled
      test_integration_restart_equivalence_stochastic_sources
      test_integration_restart_equivalence_amr_hydro
      test_integration_restart_equivalence_amr_flux_registers
      test_integration_restart_equivalence_amr_temporal_ghosts
      test_integration_reference_workflow_hydro_row_order
      test_integration_restart_equivalence_reference_workflow_hydro)
  endif()
  add_custom_target(stage1_runtime_truth_targets DEPENDS ${cosmosim_stage1_runtime_truth_targets})


  # This smoke exercises the production TreePM runtime path.  Keep it out of
  # HDF5-only builds, where the runtime correctly refuses to substitute the
  # diagnostic naive-DFT backend for FFTW.
  if(COSMOSIM_ENABLE_HDF5 AND COSMOSIM_ENABLE_FFTW)
    configure_file(
      "${CMAKE_CURRENT_SOURCE_DIR}/tests/integration/test_runtime_app_smoke.cmake.in"
      "${CMAKE_BINARY_DIR}/test_runtime_app_smoke.cmake"
      @ONLY)
    add_test(NAME integration_runtime_app_smoke
             COMMAND ${CMAKE_COMMAND} -P "${CMAKE_BINARY_DIR}/test_runtime_app_smoke.cmake")
  endif()

  if(COSMOSIM_ENABLE_MPI AND COSMOSIM_ENABLE_HDF5 AND COSMOSIM_ENABLE_FFTW)
    configure_file(
      "${CMAKE_CURRENT_SOURCE_DIR}/tests/integration/test_runtime_app_mpi_treepm_smoke.cmake.in"
      "${CMAKE_BINARY_DIR}/test_runtime_app_mpi_treepm_smoke.cmake"
      @ONLY)
    add_test(NAME integration_runtime_app_mpi_treepm_smoke_two_rank
             COMMAND ${CMAKE_COMMAND} -P "${CMAKE_BINARY_DIR}/test_runtime_app_mpi_treepm_smoke.cmake")
    set_tests_properties(integration_runtime_app_mpi_treepm_smoke_two_rank PROPERTIES PROCESSORS 2)
  endif()

  add_executable(
    test_validation_unit
    tests/validation/test_validation_unit.cpp
    tests/validation/validation_tolerance.cpp)
  target_link_libraries(test_validation_unit PRIVATE cosmosim_all)
  target_include_directories(test_validation_unit PRIVATE tests/validation tests/support)
  target_compile_definitions(test_validation_unit PRIVATE COSMOSIM_SOURCE_DIR="${PROJECT_SOURCE_DIR}")
  add_test(NAME validation_unit COMMAND test_validation_unit)
  set_tests_properties(validation_unit PROPERTIES LABELS "validation;unit")

  add_executable(
    test_validation_integration
    tests/validation/test_validation_integration.cpp
    tests/validation/validation_tolerance.cpp)
  target_link_libraries(test_validation_integration PRIVATE cosmosim_all)
  target_include_directories(test_validation_integration PRIVATE tests/validation tests/support)
  target_compile_definitions(test_validation_integration PRIVATE COSMOSIM_SOURCE_DIR="${PROJECT_SOURCE_DIR}")
  add_test(NAME validation_integration COMMAND test_validation_integration)
  set_tests_properties(validation_integration PROPERTIES LABELS "validation;integration")

  add_executable(
    test_validation_convergence
    tests/validation/test_validation_convergence.cpp
    tests/validation/validation_tolerance.cpp)
  target_link_libraries(test_validation_convergence PRIVATE cosmosim_all)
  target_include_directories(test_validation_convergence PRIVATE tests/validation tests/support)
  target_compile_definitions(
    test_validation_convergence PRIVATE
    COSMOSIM_SOURCE_DIR="${PROJECT_SOURCE_DIR}"
    COSMOSIM_VALIDATION_OUTPUT_DIR="${CMAKE_BINARY_DIR}/validation_outputs")
  add_test(NAME validation_convergence COMMAND test_validation_convergence)
  set_tests_properties(validation_convergence PROPERTIES LABELS "validation;numerical_convergence;regression")
  add_custom_target(
    generate_validation_research_evidence
    COMMAND $<TARGET_FILE:test_validation_convergence> --promote-research-evidence
    DEPENDS test_validation_convergence
    WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
    COMMENT "Explicitly regenerate curated research-grade convergence evidence")

  add_executable(
    test_validation_hydro_classics
    tests/validation/test_validation_hydro_classics.cpp)
  target_link_libraries(test_validation_hydro_classics PRIVATE cosmosim_hydro)
  target_include_directories(test_validation_hydro_classics PRIVATE tests tests/validation)
  add_test(NAME validation_hydro_classics COMMAND test_validation_hydro_classics)
  set_tests_properties(validation_hydro_classics PROPERTIES LABELS "validation;integration;hydro")

  add_executable(
    test_unit_amr_hydro_validation_metric_contract
    tests/unit/test_amr_hydro_validation_metric_contract.cpp)
  target_link_libraries(test_unit_amr_hydro_validation_metric_contract PRIVATE cosmosim_amr cosmosim_hydro cosmosim_core)
  add_test(
    NAME unit_amr_hydro_validation_metric_contract
    COMMAND test_unit_amr_hydro_validation_metric_contract)
  set_tests_properties(
    unit_amr_hydro_validation_metric_contract
    PROPERTIES LABELS "unit;contract;validation_metrics;amr;hydro")

  add_executable(
    test_unit_amr_cross_code_metadata_contract
    tests/unit/test_amr_cross_code_metadata_contract.cpp)
  target_link_libraries(test_unit_amr_cross_code_metadata_contract PRIVATE cosmosim_core)
  add_test(
    NAME unit_amr_cross_code_metadata_contract
    COMMAND test_unit_amr_cross_code_metadata_contract)
  set_tests_properties(
    unit_amr_cross_code_metadata_contract
    PROPERTIES LABELS "unit;contract;validation_metadata;amr")

  add_executable(test_validation_star_formation_models tests/validation/test_validation_star_formation_models.cpp)
  target_link_libraries(test_validation_star_formation_models PRIVATE cosmosim_all)
  add_test(NAME validation_star_formation_models COMMAND test_validation_star_formation_models)
  set_tests_properties(validation_star_formation_models PROPERTIES LABELS "validation;star_formation;effective_ism;adaptive")

  add_executable(
    test_validation_regression
    tests/validation/test_validation_regression.cpp
    tests/validation/validation_tolerance.cpp)
  target_link_libraries(test_validation_regression PRIVATE cosmosim_all)
  target_include_directories(test_validation_regression PRIVATE tests/validation)
  target_compile_definitions(test_validation_regression PRIVATE COSMOSIM_SOURCE_DIR="${PROJECT_SOURCE_DIR}")
  add_test(NAME validation_regression COMMAND test_validation_regression)
  set_tests_properties(validation_regression PROPERTIES LABELS "validation;regression")

  add_executable(
    test_validation_power_spectrum_mesh_consistency
    tests/validation/test_validation_power_spectrum_mesh_consistency.cpp
    tests/validation/validation_tolerance.cpp)
  target_link_libraries(
    test_validation_power_spectrum_mesh_consistency PRIVATE cosmosim_analysis cosmosim_core)
  target_include_directories(
    test_validation_power_spectrum_mesh_consistency PRIVATE tests/validation)
  target_compile_definitions(
    test_validation_power_spectrum_mesh_consistency
    PRIVATE COSMOSIM_SOURCE_DIR="${PROJECT_SOURCE_DIR}")
  add_test(
    NAME validation_power_spectrum_mesh_consistency
    COMMAND test_validation_power_spectrum_mesh_consistency)
  set_tests_properties(
    validation_power_spectrum_mesh_consistency
    PROPERTIES LABELS "validation;analysis;power_spectrum;regression;mesh_consistency")

  # Phase-2 production TreePM validation is an FFT-backed acceptance gate.
  # Diagnostic naive DFT remains available to focused solver/unit fixtures, but
  # must not be registered as successful production-gravity validation.
  if(COSMOSIM_ENABLE_FFTW)
    add_executable(test_validation_phase2_mpi_gravity tests/validation/test_validation_phase2_mpi_gravity.cpp)
    target_link_libraries(test_validation_phase2_mpi_gravity PRIVATE cosmosim_all)
    add_test(NAME validation_phase2_mpi_gravity_single_rank COMMAND test_validation_phase2_mpi_gravity)
    set_tests_properties(
      validation_phase2_mpi_gravity_single_rank
      PROPERTIES LABELS "validation;gravity;fftw;distributed_contract")
    if(COSMOSIM_ENABLE_MPI)
      add_test(
        NAME validation_phase2_mpi_gravity_two_rank
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
                $<TARGET_FILE:test_validation_phase2_mpi_gravity> ${MPIEXEC_POSTFLAGS})
      add_test(
        NAME validation_phase2_mpi_gravity_three_rank
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
                $<TARGET_FILE:test_validation_phase2_mpi_gravity> ${MPIEXEC_POSTFLAGS})
      set_tests_properties(
        validation_phase2_mpi_gravity_two_rank
        validation_phase2_mpi_gravity_three_rank
        PROPERTIES LABELS "validation;gravity;fftw;mpi;distributed")
      set_tests_properties(validation_phase2_mpi_gravity_two_rank PROPERTIES PROCESSORS 2)
      set_tests_properties(validation_phase2_mpi_gravity_three_rank PROPERTIES PROCESSORS 3)
      add_test(
        NAME validation_phase2_mpi_gravity_four_rank
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
                $<TARGET_FILE:test_validation_phase2_mpi_gravity> ${MPIEXEC_POSTFLAGS})
      set_tests_properties(
        validation_phase2_mpi_gravity_four_rank
        PROPERTIES LABELS "validation;gravity;fftw;mpi;distributed" PROCESSORS 4)
      add_test(
        NAME validation_phase2_mpi_gravity_eight_rank
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 8 ${MPIEXEC_PREFLAGS}
                $<TARGET_FILE:test_validation_phase2_mpi_gravity> ${MPIEXEC_POSTFLAGS})
      set_tests_properties(
        validation_phase2_mpi_gravity_eight_rank
        PROPERTIES LABELS "validation;gravity;fftw;mpi;distributed" PROCESSORS 8)
    endif()
  endif()

  add_executable(
    test_validation_periodic_ewald_reference
    tests/validation/test_periodic_ewald_reference.cpp
    tests/support/periodic_ewald_reference.cpp)
  target_include_directories(test_validation_periodic_ewald_reference PRIVATE tests/support)
  add_test(NAME validation_periodic_ewald_reference COMMAND test_validation_periodic_ewald_reference)
  set_tests_properties(
    validation_periodic_ewald_reference
    PROPERTIES LABELS "validation;gravity;periodic;reference")

  add_executable(
    test_validation_tree_pm_ewald_accuracy
    tests/validation/test_tree_pm_ewald_accuracy.cpp
    tests/support/periodic_ewald_reference.cpp)
  target_include_directories(test_validation_tree_pm_ewald_accuracy PRIVATE tests/support)
  target_link_libraries(test_validation_tree_pm_ewald_accuracy PRIVATE cosmosim_gravity)
  if(COSMOSIM_ENABLE_FFTW)
    add_test(NAME validation_tree_pm_ewald_accuracy COMMAND test_validation_tree_pm_ewald_accuracy)
    set_tests_properties(
      validation_tree_pm_ewald_accuracy
      PROPERTIES LABELS "validation;gravity;treepm;periodic;accuracy;production_fft")
    add_test(
      NAME validation_tree_pm_ewald_accuracy_first_light_softening
      COMMAND test_validation_tree_pm_ewald_accuracy)
    set_tests_properties(
      validation_tree_pm_ewald_accuracy_first_light_softening
      PROPERTIES
        ENVIRONMENT "COSMOSIM_EWALD_SOFTENING_RATIO=0.016637952"
        LABELS "validation;gravity;treepm;periodic;accuracy;production_fft;finite_softening;first_light")
  endif()

  add_executable(
    test_validation_dmo_zeldovich_workflow
    tests/validation/test_dmo_zeldovich_workflow.cpp)
  target_link_libraries(test_validation_dmo_zeldovich_workflow PRIVATE cosmosim_all)
  target_compile_definitions(
    test_validation_dmo_zeldovich_workflow PRIVATE
    COSMOSIM_DMO_ARTIFACT_ROOT="${CMAKE_BINARY_DIR}/validation_artifacts/dmo_zeldovich")
  add_test(
    NAME integration_dmo_empty_rank_smoke_single_rank
    COMMAND test_validation_dmo_zeldovich_workflow --empty-rank-smoke)
  set_tests_properties(
    integration_dmo_empty_rank_smoke_single_rank
    PROPERTIES LABELS "integration;gravity;dmo;workflow;empty_rank")
  if(COSMOSIM_ENABLE_HDF5 AND COSMOSIM_ENABLE_FFTW)
    add_test(NAME validation_dmo_zeldovich_workflow_single_rank COMMAND test_validation_dmo_zeldovich_workflow)
    set_tests_properties(
      validation_dmo_zeldovich_workflow_single_rank
      PROPERTIES FIXTURES_SETUP dmo_zeldovich_rank_state_artifacts
                 LABELS "validation;gravity;dmo;workflow;restart;power_spectrum")
  endif()
  if(COSMOSIM_ENABLE_MPI)
    foreach(dmo_empty_np IN ITEMS 4 8)
      add_test(
        NAME integration_dmo_empty_rank_smoke_mpi_np${dmo_empty_np}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${dmo_empty_np} ${MPIEXEC_PREFLAGS}
                $<TARGET_FILE:test_validation_dmo_zeldovich_workflow> --empty-rank-smoke ${MPIEXEC_POSTFLAGS})
      set_tests_properties(
        integration_dmo_empty_rank_smoke_mpi_np${dmo_empty_np}
        PROPERTIES PROCESSORS ${dmo_empty_np}
                   LABELS "integration;gravity;dmo;workflow;empty_rank;mpi;distributed")
    endforeach()
  endif()
  if(COSMOSIM_ENABLE_MPI AND COSMOSIM_ENABLE_HDF5 AND COSMOSIM_ENABLE_FFTW)
    foreach(dmo_np IN ITEMS 2 3 4 8)
      add_test(
        NAME validation_dmo_zeldovich_workflow_mpi_np${dmo_np}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${dmo_np} ${MPIEXEC_PREFLAGS}
                $<TARGET_FILE:test_validation_dmo_zeldovich_workflow> ${MPIEXEC_POSTFLAGS})
      set_tests_properties(
        validation_dmo_zeldovich_workflow_mpi_np${dmo_np}
        PROPERTIES PROCESSORS ${dmo_np}
                   FIXTURES_SETUP dmo_zeldovich_rank_state_artifacts
                   LABELS "validation;gravity;dmo;workflow;restart;power_spectrum;mpi;distributed")
    endforeach()
    add_test(
      NAME validation_dmo_zeldovich_workflow_rank_equivalence
      COMMAND $<TARGET_FILE:test_validation_dmo_zeldovich_workflow> --compare-rank-artifacts)
    set_tests_properties(
      validation_dmo_zeldovich_workflow_rank_equivalence
      PROPERTIES
        DEPENDS "validation_dmo_zeldovich_workflow_single_rank;validation_dmo_zeldovich_workflow_mpi_np2;validation_dmo_zeldovich_workflow_mpi_np3;validation_dmo_zeldovich_workflow_mpi_np4;validation_dmo_zeldovich_workflow_mpi_np8"
        FIXTURES_REQUIRED dmo_zeldovich_rank_state_artifacts
        LABELS "validation;gravity;dmo;workflow;mpi;distributed;rank_equivalence")
  endif()

  if(COSMOSIM_ENABLE_PYTHON)
    add_test(
      NAME python_bindings_smoke
      COMMAND ${CMAKE_COMMAND} -E env PYTHONPATH=${CMAKE_BINARY_DIR}/python ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/tests/python/test_bindings_smoke.py)

    if(COSMOSIM_ENABLE_HDF5)
      add_test(
        NAME python_snapshot_diagnostic_integration
        COMMAND ${CMAKE_COMMAND} -E env PYTHONPATH=${CMAKE_BINARY_DIR}/python ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/tests/python/test_snapshot_diagnostic_integration.py)
    endif()
  endif()

  cosmosim_apply_test_target_defaults()
  cosmosim_apply_test_defaults()
endif()

if(COSMOSIM_ENABLE_BENCHMARKS)
  add_executable(bench_gravity_kernels bench/gravity/bench_gravity_kernels.cpp)
  target_link_libraries(bench_gravity_kernels PRIVATE cosmosim_gravity)
  target_include_directories(bench_gravity_kernels PRIVATE ${PROJECT_SOURCE_DIR})

  add_executable(bench_hydro_kernels bench/hydro/bench_hydro_kernels.cpp)
  target_link_libraries(bench_hydro_kernels PRIVATE cosmosim_hydro)
  target_include_directories(bench_hydro_kernels PRIVATE ${PROJECT_SOURCE_DIR})

  add_executable(bench_amr_regrid_kernel bench/amr/bench_amr_regrid_kernel.cpp)
  target_link_libraries(bench_amr_regrid_kernel PRIVATE cosmosim_amr)
  target_include_directories(bench_amr_regrid_kernel PRIVATE ${PROJECT_SOURCE_DIR})

  add_executable(bench_io_restart_kernel bench/io/bench_io_restart_kernel.cpp)
  target_link_libraries(bench_io_restart_kernel PRIVATE cosmosim_io)
  target_include_directories(bench_io_restart_kernel PRIVATE ${PROJECT_SOURCE_DIR})

  add_executable(bench_mini_run_pipeline bench/mini/bench_mini_run_pipeline.cpp)
  target_link_libraries(bench_mini_run_pipeline PRIVATE cosmosim_all)
  target_include_directories(bench_mini_run_pipeline PRIVATE ${PROJECT_SOURCE_DIR})

  add_executable(bench_reference_workflow bench/mini/bench_reference_workflow.cpp)
  target_link_libraries(bench_reference_workflow PRIVATE cosmosim_all)
  target_include_directories(bench_reference_workflow PRIVATE ${PROJECT_SOURCE_DIR})

  add_executable(bench_layout_smoke bench/bench_layout_smoke.cpp)
  target_link_libraries(bench_layout_smoke PRIVATE cosmosim_core)
  target_include_directories(bench_layout_smoke PRIVATE ${PROJECT_SOURCE_DIR}/src)
  target_compile_definitions(bench_layout_smoke PRIVATE COSMOSIM_BENCHMARK_BUILD=1)


  add_executable(bench_parallel_decomposition_exchange bench/bench_parallel_decomposition_exchange.cpp)
  target_link_libraries(bench_parallel_decomposition_exchange PRIVATE cosmosim_parallel)
  target_include_directories(bench_parallel_decomposition_exchange PRIVATE ${PROJECT_SOURCE_DIR})

  add_executable(bench_config_parser bench/bench_config_parser.cpp)
  target_link_libraries(bench_config_parser PRIVATE cosmosim_core)

  add_executable(bench_docs_reference_scan bench/bench_docs_reference_scan.cpp)
  target_include_directories(bench_docs_reference_scan PRIVATE ${PROJECT_SOURCE_DIR})

  add_executable(bench_simulation_state_resize bench/bench_simulation_state_resize.cpp)
  target_link_libraries(bench_simulation_state_resize PRIVATE cosmosim_core)

  add_executable(bench_soa_storage bench/bench_soa_storage.cpp)
  target_link_libraries(bench_soa_storage PRIVATE cosmosim_core)

  add_executable(bench_species_hot_cold_access bench/bench_species_hot_cold_access.cpp)
  target_link_libraries(bench_species_hot_cold_access PRIVATE cosmosim_core)



  add_executable(bench_profiling_overhead bench/bench_profiling_overhead.cpp)
  target_link_libraries(bench_profiling_overhead PRIVATE cosmosim_core)

  add_executable(bench_memory_accounting_overhead bench/bench_memory_accounting_overhead.cpp)
  target_link_libraries(bench_memory_accounting_overhead PRIVATE cosmosim_core)

  add_executable(bench_core_migration_memory bench/bench_core_migration_memory.cpp)
  target_link_libraries(bench_core_migration_memory PRIVATE cosmosim_core)

  add_executable(bench_scheduler_overhead bench/bench_scheduler_overhead.cpp)
  add_executable(bench_analysis_diagnostics_overhead bench/bench_analysis_diagnostics_overhead.cpp)
  add_executable(bench_fof_grouping bench/bench_fof_grouping.cpp)
  target_link_libraries(bench_scheduler_overhead PRIVATE cosmosim_core)
  target_link_libraries(bench_analysis_diagnostics_overhead PRIVATE cosmosim_all)
  target_link_libraries(bench_fof_grouping PRIVATE cosmosim_all)

  add_executable(bench_mode_dispatch bench/bench_mode_dispatch.cpp)
  target_link_libraries(bench_mode_dispatch PRIVATE cosmosim_core)
  add_executable(bench_active_kernel_compact_vs_full bench/bench_active_kernel_compact_vs_full.cpp)
  target_link_libraries(bench_active_kernel_compact_vs_full PRIVATE cosmosim_core)
  add_executable(bench_time_bin_scheduler bench/bench_time_bin_scheduler.cpp)
  target_link_libraries(bench_time_bin_scheduler PRIVATE cosmosim_core)


  add_executable(bench_ic_reader_streaming bench/bench_ic_reader_streaming.cpp)
  target_link_libraries(bench_ic_reader_streaming PRIVATE cosmosim_io)

  add_executable(bench_snapshot_hdf5_io bench/bench_snapshot_hdf5_io.cpp)
  target_link_libraries(bench_snapshot_hdf5_io PRIVATE cosmosim_io)

  add_executable(bench_restart_checkpoint_latency bench/bench_restart_checkpoint_latency.cpp)
  target_link_libraries(bench_restart_checkpoint_latency PRIVATE cosmosim_io)

  add_executable(bench_pm_solver bench/bench_pm_solver.cpp)
  target_link_libraries(bench_pm_solver PRIVATE cosmosim_gravity)

  add_executable(bench_pm_gpu_transfer bench/bench_pm_gpu_transfer.cpp)
  target_link_libraries(bench_pm_gpu_transfer PRIVATE cosmosim_gravity)

  add_executable(bench_tree_gravity bench/bench_tree_gravity.cpp)
  target_link_libraries(bench_tree_gravity PRIVATE cosmosim_gravity)

  add_executable(bench_tree_pm_coupling bench/bench_tree_pm_coupling.cpp)
  target_link_libraries(bench_tree_pm_coupling PRIVATE cosmosim_gravity)
  add_executable(bench_tree_pm_force_error_map bench/bench_tree_pm_force_error_map.cpp)
  target_link_libraries(bench_tree_pm_force_error_map PRIVATE cosmosim_gravity)
  target_compile_definitions(bench_tree_pm_force_error_map PRIVATE COSMOSIM_SOURCE_DIR="${PROJECT_SOURCE_DIR}")
  add_custom_target(generate_tree_pm_force_error_map
    COMMAND bench_tree_pm_force_error_map
    DEPENDS bench_tree_pm_force_error_map
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

  add_executable(bench_pm_only_scaling_mpi bench/bench_pm_only_scaling_mpi.cpp)
  target_link_libraries(bench_pm_only_scaling_mpi PRIVATE cosmosim_gravity cosmosim_parallel)

  add_executable(bench_tree_only_scaling_mpi bench/bench_tree_only_scaling_mpi.cpp)
  target_link_libraries(bench_tree_only_scaling_mpi PRIVATE cosmosim_gravity)

  if(COSMOSIM_ENABLE_MPI)
    add_custom_target(generate_mpi_gravity_scaling_artifacts
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS} $<TARGET_FILE:bench_pm_only_scaling_mpi> ${MPIEXEC_POSTFLAGS}
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS} $<TARGET_FILE:bench_pm_only_scaling_mpi> ${MPIEXEC_POSTFLAGS}
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS} $<TARGET_FILE:bench_tree_only_scaling_mpi> ${MPIEXEC_POSTFLAGS}
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS} $<TARGET_FILE:bench_tree_only_scaling_mpi> ${MPIEXEC_POSTFLAGS}
      DEPENDS bench_pm_only_scaling_mpi bench_tree_only_scaling_mpi
      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
  endif()

  add_executable(bench_amr_regrid_sync bench/bench_amr_regrid_sync.cpp)
  target_link_libraries(bench_amr_regrid_sync PRIVATE cosmosim_amr)

  add_executable(bench_hydro_face_sweep bench/bench_hydro_face_sweep.cpp)
  target_link_libraries(bench_hydro_face_sweep PRIVATE cosmosim_hydro)

  add_executable(bench_cooling_heating_source bench/bench_cooling_heating_source.cpp)
  target_link_libraries(bench_cooling_heating_source PRIVATE cosmosim_all)

  add_executable(bench_stellar_evolution_bookkeeping bench/bench_stellar_evolution_bookkeeping.cpp)
  target_link_libraries(bench_stellar_evolution_bookkeeping PRIVATE cosmosim_all)
  add_executable(bench_metal_enrichment bench/bench_metal_enrichment.cpp)
  target_link_libraries(bench_metal_enrichment PRIVATE cosmosim_all)
  target_compile_definitions(bench_metal_enrichment PRIVATE COSMOSIM_SOURCE_DIR="${PROJECT_SOURCE_DIR}")
  add_executable(bench_metal_diffusion bench/bench_metal_diffusion.cpp)
  target_link_libraries(bench_metal_diffusion PRIVATE cosmosim_all)

  add_executable(bench_star_formation_spawn bench/bench_star_formation_spawn.cpp)
  target_link_libraries(bench_star_formation_spawn PRIVATE cosmosim_all)
  add_executable(bench_effective_multiphase_ism bench/bench_effective_multiphase_ism.cpp)
  target_link_libraries(bench_effective_multiphase_ism PRIVATE cosmosim_all)
  add_executable(bench_black_hole_agn_sparse_update bench/bench_black_hole_agn_sparse_update.cpp)
  target_link_libraries(bench_black_hole_agn_sparse_update PRIVATE cosmosim_all)
  add_executable(bench_stellar_feedback_deposition bench/bench_stellar_feedback_deposition.cpp)
  target_link_libraries(bench_stellar_feedback_deposition PRIVATE cosmosim_all)
  add_executable(bench_tracer_update bench/bench_tracer_update.cpp)
  target_link_libraries(bench_tracer_update PRIVATE cosmosim_all)

  add_executable(bench_validation_ladder bench/bench_validation_ladder.cpp)
  target_link_libraries(bench_validation_ladder PRIVATE cosmosim_all)

  add_custom_target(stage6_benchmarks
    DEPENDS
      bench_active_kernel_compact_vs_full
      bench_gravity_kernels
      bench_hydro_kernels
      bench_memory_accounting_overhead
      bench_profiling_overhead)

  if(COSMOSIM_ENABLE_PYTHON AND COSMOSIM_ENABLE_HDF5)
    add_custom_target(bench_python_snapshot_access
      COMMAND ${CMAKE_COMMAND} -E env PYTHONPATH=${CMAKE_BINARY_DIR}/python ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/bench/python/bench_snapshot_access.py --particles=200000 --iterations=8
      DEPENDS cosmosim_python_package
      USES_TERMINAL)
  endif()

endif()

cosmosim_generate_feature_summary(
  OUTPUT_TEXT "${CMAKE_BINARY_DIR}/cosmosim_feature_summary.txt"
  OUTPUT_JSON "${CMAKE_BINARY_DIR}/cosmosim_build_metadata.json"
  TARGETS
    cosmosim_core
    ${COSMOSIM_MODULE_LIBS}
    cosmosim_harness)

message(STATUS "\nCosmoSim feature summary\n-----------------------")
file(STRINGS "${CMAKE_BINARY_DIR}/cosmosim_feature_summary.txt" cosmosim_summary_lines)
foreach(summary_line IN LISTS cosmosim_summary_lines)
  message(STATUS "${summary_line}")
endforeach()
message(STATUS "Build metadata JSON: ${CMAKE_BINARY_DIR}/cosmosim_build_metadata.json")
