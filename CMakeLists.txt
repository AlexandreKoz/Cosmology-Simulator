cmake_minimum_required(VERSION 3.24)

project(
  cosmosim
  VERSION 0.1.0
  DESCRIPTION "CosmoSim: desktop-first cosmological simulation framework skeleton"
  LANGUAGES CXX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(cosmosim_feature_summary)
include(cosmosim_find_fftw)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

option(COSMOSIM_ENABLE_MPI "Enable MPI backend" OFF)
option(COSMOSIM_ENABLE_HDF5 "Enable HDF5 I/O backend" OFF)
option(COSMOSIM_ENABLE_FFTW "Enable FFTW backend" OFF)
option(COSMOSIM_ENABLE_CUDA "Enable CUDA/cuFFT backend" OFF)
option(COSMOSIM_ENABLE_PYTHON "Enable Python tooling integration" OFF)
option(COSMOSIM_ENABLE_TRACERS "Enable tracer update module" OFF)
option(COSMOSIM_ENABLE_TESTS "Build unit and integration tests" ON)
option(COSMOSIM_ENABLE_BENCHMARKS "Build benchmark/profiling hooks" ON)
option(COSMOSIM_ENABLE_PROFILING "Enable built-in profiling instrumentation" ON)
option(COSMOSIM_ENABLE_LTO "Enable interprocedural optimization (LTO)" OFF)

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
endif()

if(COSMOSIM_ENABLE_FFTW)
  cosmosim_find_fftw(COSMOSIM_FFTW_TARGET)
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
  src/core/profiling.cpp)

target_include_directories(cosmosim_core
  PUBLIC
    ${PROJECT_SOURCE_DIR}/include
    ${COSMOSIM_GENERATED_INCLUDE_DIR}
  PRIVATE
    ${PROJECT_SOURCE_DIR}/src)

target_compile_features(cosmosim_core PUBLIC cxx_std_20)

target_compile_definitions(cosmosim_core
  PUBLIC
    COSMOSIM_VERSION_MAJOR=${PROJECT_VERSION_MAJOR}
    COSMOSIM_VERSION_MINOR=${PROJECT_VERSION_MINOR}
    COSMOSIM_VERSION_PATCH=${PROJECT_VERSION_PATCH}
    COSMOSIM_ENABLE_PROFILING=$<BOOL:${COSMOSIM_ENABLE_PROFILING}>)

if(COSMOSIM_ENABLE_LTO)
  set_property(TARGET cosmosim_core PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

add_library(cosmosim_gravity STATIC
  src/gravity/gravity_module.cpp
  src/gravity/pm_solver.cpp
  src/gravity/tree_ordering.cpp
  src/gravity/tree_gravity.cpp
  src/gravity/tree_pm_coupling.cpp)
if(COSMOSIM_ENABLE_CUDA)
  target_sources(cosmosim_gravity PRIVATE src/gravity/pm_cuda_kernels.cu)
  set_target_properties(cosmosim_gravity PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
endif()
add_library(cosmosim_hydro STATIC
  src/hydro/hydro_module.cpp
  src/hydro/hydro_core_solver.cpp
  src/hydro/hydro_reconstruction.cpp
  src/hydro/hydro_riemann.cpp)
add_library(cosmosim_amr STATIC
  src/amr/amr_module.cpp
  src/amr/amr_framework.cpp)
add_library(cosmosim_physics STATIC
  src/physics/physics_module.cpp
  src/physics/cooling_heating.cpp
  src/physics/black_hole_agn.cpp
  src/physics/star_formation.cpp
  src/physics/stellar_evolution.cpp
  src/physics/stellar_feedback.cpp
  src/physics/tracer_support.cpp)
add_library(cosmosim_io STATIC
  src/io/io_module.cpp
  src/io/ic_reader.cpp
  src/io/snapshot_hdf5.cpp
  src/io/restart_checkpoint.cpp)
add_library(cosmosim_analysis STATIC
  src/analysis/analysis_module.cpp
  src/analysis/diagnostics.cpp
  src/analysis/halo_workflow.cpp)
add_library(cosmosim_parallel STATIC
  src/parallel/parallel_module.cpp
  src/parallel/distributed_memory.cpp)
add_library(cosmosim_utils STATIC
  src/utils/utils_module.cpp)
add_library(cosmosim_workflows STATIC
  src/workflows/reference_workflow.cpp)

set(COSMOSIM_MODULE_LIBS
  cosmosim_gravity
  cosmosim_hydro
  cosmosim_amr
  cosmosim_physics
  cosmosim_io
  cosmosim_analysis
  cosmosim_parallel
  cosmosim_utils
  cosmosim_workflows)

foreach(module_lib IN LISTS COSMOSIM_MODULE_LIBS)
  target_include_directories(${module_lib}
    PUBLIC
      ${PROJECT_SOURCE_DIR}/include
      ${COSMOSIM_GENERATED_INCLUDE_DIR}
    PRIVATE
      ${PROJECT_SOURCE_DIR}/src)
  target_link_libraries(${module_lib} PUBLIC cosmosim_core)
  target_compile_features(${module_lib} PUBLIC cxx_std_20)

  if(COSMOSIM_ENABLE_LTO)
    set_property(TARGET ${module_lib} PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
  endif()
endforeach()

target_link_libraries(cosmosim_workflows
  PUBLIC
    cosmosim_analysis
    cosmosim_io
    cosmosim_physics)

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
  target_link_libraries(cosmosim_gravity PUBLIC ${COSMOSIM_FFTW_TARGET})
endif()

if(COSMOSIM_ENABLE_CUDA)
  target_link_libraries(cosmosim_gravity PUBLIC CUDA::cufft CUDA::cudart)
endif()

add_library(cosmosim_all INTERFACE)
target_link_libraries(cosmosim_all
  INTERFACE
    cosmosim_core
    ${COSMOSIM_MODULE_LIBS})

add_executable(cosmosim_harness src/core/harness_main.cpp)
target_link_libraries(cosmosim_harness PRIVATE cosmosim_all)


if(COSMOSIM_ENABLE_PYTHON)
  pybind11_add_module(_cosmosim MODULE src/python_bindings.cpp)
  target_link_libraries(_cosmosim PRIVATE cosmosim_all)
  target_compile_features(_cosmosim PRIVATE cxx_std_20)

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

  add_executable(test_unit_simulation_state tests/unit/test_simulation_state.cpp)
  target_link_libraries(test_unit_simulation_state PRIVATE cosmosim_all)
  add_test(NAME unit_simulation_state COMMAND test_unit_simulation_state)

  add_executable(test_unit_config_parser tests/unit/test_config_parser.cpp)
  target_link_libraries(test_unit_config_parser PRIVATE cosmosim_all)
  add_test(NAME unit_config_parser COMMAND test_unit_config_parser)

  add_executable(test_unit_soa_storage tests/unit/test_soa_storage.cpp)
  target_link_libraries(test_unit_soa_storage PRIVATE cosmosim_core)
  add_test(NAME unit_soa_storage COMMAND test_unit_soa_storage)



  add_executable(test_unit_profiling tests/unit/test_profiling.cpp)
  target_link_libraries(test_unit_profiling PRIVATE cosmosim_core)
  add_test(NAME unit_profiling COMMAND test_unit_profiling)

  add_executable(test_unit_time_integration tests/unit/test_time_integration.cpp)
  target_link_libraries(test_unit_time_integration PRIVATE cosmosim_core)
  add_test(NAME unit_time_integration COMMAND test_unit_time_integration)

  add_executable(test_unit_simulation_mode tests/unit/test_simulation_mode.cpp)
  target_link_libraries(test_unit_simulation_mode PRIVATE cosmosim_core)
  add_test(NAME unit_simulation_mode COMMAND test_unit_simulation_mode)
  add_executable(test_unit_species_state_organization tests/unit/test_species_state_organization.cpp)
  target_link_libraries(test_unit_species_state_organization PRIVATE cosmosim_core)
  add_test(NAME unit_species_state_organization COMMAND test_unit_species_state_organization)

  add_executable(test_unit_hot_cold_sidecar_layout tests/unit/test_hot_cold_sidecar_layout.cpp)
  target_link_libraries(test_unit_hot_cold_sidecar_layout PRIVATE cosmosim_core)
  add_test(NAME unit_hot_cold_sidecar_layout COMMAND test_unit_hot_cold_sidecar_layout)


  add_executable(test_unit_parallel_distributed_memory tests/unit/test_parallel_distributed_memory.cpp)
  target_link_libraries(test_unit_parallel_distributed_memory PRIVATE cosmosim_parallel)
  add_test(NAME unit_parallel_distributed_memory COMMAND test_unit_parallel_distributed_memory)

  add_executable(test_unit_hydro_reconstruction tests/unit/test_hydro_reconstruction.cpp)
  target_link_libraries(test_unit_hydro_reconstruction PRIVATE cosmosim_all)
  add_test(NAME unit_hydro_reconstruction COMMAND test_unit_hydro_reconstruction)

  add_executable(test_unit_hydro_riemann tests/unit/test_hydro_riemann.cpp)
  target_link_libraries(test_unit_hydro_riemann PRIVATE cosmosim_all)
  add_test(NAME unit_hydro_riemann COMMAND test_unit_hydro_riemann)


  add_executable(test_unit_cooling_heating tests/unit/test_cooling_heating.cpp)
  target_link_libraries(test_unit_cooling_heating PRIVATE cosmosim_all)
  add_test(NAME unit_cooling_heating COMMAND test_unit_cooling_heating)


  add_executable(test_unit_star_formation tests/unit/test_star_formation.cpp)
  target_link_libraries(test_unit_star_formation PRIVATE cosmosim_all)
  add_test(NAME unit_star_formation COMMAND test_unit_star_formation)
  add_executable(test_unit_black_hole_agn tests/unit/test_black_hole_agn.cpp)
  target_link_libraries(test_unit_black_hole_agn PRIVATE cosmosim_all)
  add_test(NAME unit_black_hole_agn COMMAND test_unit_black_hole_agn)

  add_executable(test_unit_stellar_evolution tests/unit/test_stellar_evolution.cpp)
  target_link_libraries(test_unit_stellar_evolution PRIVATE cosmosim_all)
  add_test(NAME unit_stellar_evolution COMMAND test_unit_stellar_evolution)
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
  add_test(NAME integration_stellar_evolution_ssp COMMAND test_integration_stellar_evolution_ssp)

  add_executable(test_integration_star_formation_box tests/integration/test_star_formation_box.cpp)
  target_link_libraries(test_integration_star_formation_box PRIVATE cosmosim_all)
  add_test(NAME integration_star_formation_box COMMAND test_integration_star_formation_box)
  add_executable(test_integration_black_hole_agn_toy tests/integration/test_black_hole_agn_toy.cpp)
  target_link_libraries(test_integration_black_hole_agn_toy PRIVATE cosmosim_all)
  add_test(NAME integration_black_hole_agn_toy COMMAND test_integration_black_hole_agn_toy)
  add_executable(test_integration_stellar_feedback_box tests/integration/test_stellar_feedback_box.cpp)
  target_link_libraries(test_integration_stellar_feedback_box PRIVATE cosmosim_all)
  add_test(NAME integration_stellar_feedback_box COMMAND test_integration_stellar_feedback_box)

  add_executable(test_integration_simulation_state_pipeline tests/integration/test_simulation_state_pipeline.cpp)
  target_link_libraries(test_integration_simulation_state_pipeline PRIVATE cosmosim_all)
  add_test(NAME integration_simulation_state_pipeline COMMAND test_integration_simulation_state_pipeline)

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



  add_executable(test_integration_profiling_mini_run tests/integration/test_profiling_mini_run.cpp)
  target_link_libraries(test_integration_profiling_mini_run PRIVATE cosmosim_core)
  add_test(NAME integration_profiling_mini_run COMMAND test_integration_profiling_mini_run)

  add_executable(test_integration_time_integration_loop tests/integration/test_time_integration_loop.cpp)
  target_link_libraries(test_integration_time_integration_loop PRIVATE cosmosim_all)
  add_test(NAME integration_time_integration_loop COMMAND test_integration_time_integration_loop)

  add_executable(test_integration_hierarchical_time_bins tests/integration/test_hierarchical_time_bins.cpp)
  target_link_libraries(test_integration_hierarchical_time_bins PRIVATE cosmosim_core)
  add_test(NAME integration_hierarchical_time_bins COMMAND test_integration_hierarchical_time_bins)

  add_executable(test_integration_hierarchical_timestep_regression tests/integration/test_hierarchical_timestep_regression.cpp)
  target_link_libraries(test_integration_hierarchical_timestep_regression PRIVATE cosmosim_core)
  add_test(NAME integration_hierarchical_timestep_regression COMMAND test_integration_hierarchical_timestep_regression)


  add_executable(test_unit_analysis_diagnostics tests/unit/test_analysis_diagnostics.cpp)
  target_link_libraries(test_unit_analysis_diagnostics PRIVATE cosmosim_all)
  add_test(NAME unit_analysis_diagnostics COMMAND test_unit_analysis_diagnostics)
  add_executable(test_unit_halo_workflow tests/unit/test_halo_workflow.cpp)
  target_link_libraries(test_unit_halo_workflow PRIVATE cosmosim_all)
  add_test(NAME unit_halo_workflow COMMAND test_unit_halo_workflow)

  add_executable(test_unit_pm_solver tests/unit/test_pm_solver.cpp)
  target_link_libraries(test_unit_pm_solver PRIVATE cosmosim_gravity)
  add_test(NAME unit_pm_solver COMMAND test_unit_pm_solver)

  add_executable(test_unit_tree_gravity tests/unit/test_tree_gravity.cpp)
  target_link_libraries(test_unit_tree_gravity PRIVATE cosmosim_gravity)
  add_test(NAME unit_tree_gravity COMMAND test_unit_tree_gravity)


  add_executable(test_unit_amr_refinement tests/unit/test_amr_refinement.cpp)
  target_link_libraries(test_unit_amr_refinement PRIVATE cosmosim_amr)
  add_test(NAME unit_amr_refinement COMMAND test_unit_amr_refinement)

  add_executable(test_unit_hydro_core_solver tests/unit/test_hydro_core_solver.cpp)
  target_link_libraries(test_unit_hydro_core_solver PRIVATE cosmosim_hydro)
  add_test(NAME unit_hydro_core_solver COMMAND test_unit_hydro_core_solver)


  add_executable(test_integration_parallel_two_rank_restart tests/integration/test_parallel_two_rank_restart.cpp)
  target_link_libraries(test_integration_parallel_two_rank_restart PRIVATE cosmosim_parallel cosmosim_core)
  add_test(NAME integration_parallel_two_rank_restart COMMAND test_integration_parallel_two_rank_restart)

  add_executable(test_integration_hydro_sod_like tests/integration/test_hydro_sod_like.cpp)
  target_link_libraries(test_integration_hydro_sod_like PRIVATE cosmosim_hydro)
  add_test(NAME integration_hydro_sod_like COMMAND test_integration_hydro_sod_like)



  add_executable(test_integration_amr_static_refinement_sync tests/integration/test_amr_static_refinement_sync.cpp)
  target_link_libraries(test_integration_amr_static_refinement_sync PRIVATE cosmosim_amr)
  add_test(NAME integration_amr_static_refinement_sync COMMAND test_integration_amr_static_refinement_sync)

  add_executable(test_integration_simulation_mode_toy_runs tests/integration/test_simulation_mode_toy_runs.cpp)
  target_link_libraries(test_integration_simulation_mode_toy_runs PRIVATE cosmosim_core)
  add_test(NAME integration_simulation_mode_toy_runs COMMAND test_integration_simulation_mode_toy_runs)

  add_executable(test_integration_pm_periodic_mode tests/integration/test_pm_periodic_mode.cpp)
  target_link_libraries(test_integration_pm_periodic_mode PRIVATE cosmosim_gravity)
  add_test(NAME integration_pm_periodic_mode COMMAND test_integration_pm_periodic_mode)

  if(COSMOSIM_ENABLE_CUDA)
    add_executable(test_integration_pm_cuda_smoke tests/integration/test_pm_cuda_smoke.cpp)
    target_link_libraries(test_integration_pm_cuda_smoke PRIVATE cosmosim_gravity)
    add_test(NAME integration_pm_cuda_smoke COMMAND test_integration_pm_cuda_smoke)
  endif()


  add_executable(test_unit_ic_reader tests/unit/test_ic_reader.cpp)
  target_link_libraries(test_unit_ic_reader PRIVATE cosmosim_io)
  add_test(NAME unit_ic_reader COMMAND test_unit_ic_reader)

  add_executable(test_unit_snapshot_hdf5_schema tests/unit/test_snapshot_hdf5_schema.cpp)
  target_link_libraries(test_unit_snapshot_hdf5_schema PRIVATE cosmosim_io)
  add_test(NAME unit_snapshot_hdf5_schema COMMAND test_unit_snapshot_hdf5_schema)


  add_executable(test_unit_restart_checkpoint_schema tests/unit/test_restart_checkpoint_schema.cpp)
  target_link_libraries(test_unit_restart_checkpoint_schema PRIVATE cosmosim_io)
  add_test(NAME unit_restart_checkpoint_schema COMMAND test_unit_restart_checkpoint_schema)

  add_executable(test_unit_tree_pm_split_kernel tests/unit/test_tree_pm_split_kernel.cpp)
  target_link_libraries(test_unit_tree_pm_split_kernel PRIVATE cosmosim_gravity)
  add_test(NAME unit_tree_pm_split_kernel COMMAND test_unit_tree_pm_split_kernel)

  add_executable(test_integration_tree_pm_coupling_periodic tests/integration/test_tree_pm_coupling_periodic.cpp)
  target_link_libraries(test_integration_tree_pm_coupling_periodic PRIVATE cosmosim_gravity)
  add_test(NAME integration_tree_pm_coupling_periodic COMMAND test_integration_tree_pm_coupling_periodic)


  add_executable(test_integration_ic_generated_pipeline tests/integration/test_ic_generated_pipeline.cpp)
  target_link_libraries(test_integration_ic_generated_pipeline PRIVATE cosmosim_all)
  add_test(NAME integration_ic_generated_pipeline COMMAND test_integration_ic_generated_pipeline)

  add_executable(test_integration_snapshot_hdf5_roundtrip tests/integration/test_snapshot_hdf5_roundtrip.cpp)
  target_link_libraries(test_integration_snapshot_hdf5_roundtrip PRIVATE cosmosim_all)
  add_test(NAME integration_snapshot_hdf5_roundtrip COMMAND test_integration_snapshot_hdf5_roundtrip)
  add_executable(test_integration_tracer_advection tests/integration/test_tracer_advection.cpp)
  target_link_libraries(test_integration_tracer_advection PRIVATE cosmosim_all)
  add_test(NAME integration_tracer_advection COMMAND test_integration_tracer_advection)


  add_executable(test_integration_restart_checkpoint_roundtrip tests/integration/test_restart_checkpoint_roundtrip.cpp)
  target_link_libraries(test_integration_restart_checkpoint_roundtrip PRIVATE cosmosim_all)
  add_test(NAME integration_restart_checkpoint_roundtrip COMMAND test_integration_restart_checkpoint_roundtrip)

  add_executable(test_integration_tree_gravity_vs_direct tests/integration/test_tree_gravity_vs_direct.cpp)
  target_link_libraries(test_integration_tree_gravity_vs_direct PRIVATE cosmosim_gravity)
  add_test(NAME integration_tree_gravity_vs_direct COMMAND test_integration_tree_gravity_vs_direct)
  add_executable(test_integration_config_examples tests/integration/test_config_examples.cpp)
  target_link_libraries(test_integration_config_examples PRIVATE cosmosim_all)
  target_compile_definitions(test_integration_config_examples PRIVATE COSMOSIM_SOURCE_DIR="${PROJECT_SOURCE_DIR}")
  add_test(NAME integration_config_examples COMMAND test_integration_config_examples)

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
    "${CMAKE_CURRENT_SOURCE_DIR}/tests/integration/test_release_readiness_artifacts.cmake.in"
    "${CMAKE_BINARY_DIR}/test_release_readiness_artifacts.cmake"
    @ONLY)
  add_test(NAME integration_release_readiness_artifacts
           COMMAND ${CMAKE_COMMAND} -P "${CMAKE_BINARY_DIR}/test_release_readiness_artifacts.cmake")

  add_executable(
    test_validation_unit
    tests/validation/test_validation_unit.cpp
    tests/validation/validation_tolerance.cpp)
  target_link_libraries(test_validation_unit PRIVATE cosmosim_all)
  target_include_directories(test_validation_unit PRIVATE tests/validation)
  target_compile_definitions(test_validation_unit PRIVATE COSMOSIM_SOURCE_DIR="${PROJECT_SOURCE_DIR}")
  add_test(NAME validation_unit COMMAND test_validation_unit)
  set_tests_properties(validation_unit PROPERTIES LABELS "validation;unit")

  add_executable(
    test_validation_integration
    tests/validation/test_validation_integration.cpp
    tests/validation/validation_tolerance.cpp)
  target_link_libraries(test_validation_integration PRIVATE cosmosim_all)
  target_include_directories(test_validation_integration PRIVATE tests/validation)
  target_compile_definitions(test_validation_integration PRIVATE COSMOSIM_SOURCE_DIR="${PROJECT_SOURCE_DIR}")
  add_test(NAME validation_integration COMMAND test_validation_integration)
  set_tests_properties(validation_integration PROPERTIES LABELS "validation;integration")

  add_executable(
    test_validation_convergence
    tests/validation/test_validation_convergence.cpp
    tests/validation/validation_tolerance.cpp)
  target_link_libraries(test_validation_convergence PRIVATE cosmosim_all)
  target_include_directories(test_validation_convergence PRIVATE tests/validation)
  target_compile_definitions(test_validation_convergence PRIVATE COSMOSIM_SOURCE_DIR="${PROJECT_SOURCE_DIR}")
  add_test(NAME validation_convergence COMMAND test_validation_convergence)
  set_tests_properties(validation_convergence PROPERTIES LABELS "validation;convergence")

  add_executable(
    test_validation_regression
    tests/validation/test_validation_regression.cpp
    tests/validation/validation_tolerance.cpp)
  target_link_libraries(test_validation_regression PRIVATE cosmosim_all)
  target_include_directories(test_validation_regression PRIVATE tests/validation)
  target_compile_definitions(test_validation_regression PRIVATE COSMOSIM_SOURCE_DIR="${PROJECT_SOURCE_DIR}")
  add_test(NAME validation_regression COMMAND test_validation_regression)
  set_tests_properties(validation_regression PROPERTIES LABELS "validation;regression")

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


  add_executable(bench_amr_regrid_sync bench/bench_amr_regrid_sync.cpp)
  target_link_libraries(bench_amr_regrid_sync PRIVATE cosmosim_amr)

  add_executable(bench_hydro_face_sweep bench/bench_hydro_face_sweep.cpp)
  target_link_libraries(bench_hydro_face_sweep PRIVATE cosmosim_hydro)

  add_executable(bench_cooling_heating_source bench/bench_cooling_heating_source.cpp)
  target_link_libraries(bench_cooling_heating_source PRIVATE cosmosim_all)

  add_executable(bench_stellar_evolution_bookkeeping bench/bench_stellar_evolution_bookkeeping.cpp)
  target_link_libraries(bench_stellar_evolution_bookkeeping PRIVATE cosmosim_all)

  add_executable(bench_star_formation_spawn bench/bench_star_formation_spawn.cpp)
  target_link_libraries(bench_star_formation_spawn PRIVATE cosmosim_all)
  add_executable(bench_black_hole_agn_sparse_update bench/bench_black_hole_agn_sparse_update.cpp)
  target_link_libraries(bench_black_hole_agn_sparse_update PRIVATE cosmosim_all)
  add_executable(bench_stellar_feedback_deposition bench/bench_stellar_feedback_deposition.cpp)
  target_link_libraries(bench_stellar_feedback_deposition PRIVATE cosmosim_all)
  add_executable(bench_tracer_update bench/bench_tracer_update.cpp)
  target_link_libraries(bench_tracer_update PRIVATE cosmosim_all)

  add_executable(bench_validation_ladder bench/bench_validation_ladder.cpp)
  target_link_libraries(bench_validation_ladder PRIVATE cosmosim_all)

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
