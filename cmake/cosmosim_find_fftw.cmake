include_guard(GLOBAL)

# Robust FFTW3 discovery for development and CI feature paths.
#
# Supported discovery routes, in order:
#   1. CMake package config targets from system/conda/vcpkg installations.
#   2. CosmoSim's FindFFTW3.cmake module.
#   3. pkg-config imported targets.
#
# The exported variables intentionally stay simple for call sites:
#   - <out_target> receives the serial double-precision FFTW target.
#   - COSMOSIM_FFTW_MPI_TARGET receives the MPI FFTW target when available.
#   - FFTW_VERSION is propagated for build metadata when a version is known.
function(cosmosim_find_fftw out_target)
  set(options REQUIRE_MPI)
  cmake_parse_arguments(COSMOSIM_FIND_FFTW "${options}" "" "" ${ARGN})

  set(COSMOSIM_FFTW_TARGET "")
  set(COSMOSIM_FFTW_MPI_TARGET "")

  function(_cosmosim_select_first_existing_target out_var)
    foreach(candidate IN LISTS ARGN)
      if(TARGET "${candidate}")
        set(${out_var} "${candidate}" PARENT_SCOPE)
        return()
      endif()
    endforeach()
    set(${out_var} "" PARENT_SCOPE)
  endfunction()

  function(_cosmosim_capture_fftw_version)
    foreach(candidate_version_var IN ITEMS FFTW3_VERSION FFTW_VERSION FFTW_PKG_VERSION FFTW_MPI_PKG_VERSION)
      if(DEFINED ${candidate_version_var} AND NOT "${${candidate_version_var}}" STREQUAL "")
        set(FFTW_VERSION "${${candidate_version_var}}" PARENT_SCOPE)
        return()
      endif()
    endforeach()
  endfunction()

  # Prefer package config when an installation provides imported targets.  Target
  # names are not fully standardized across FFTW package providers, so accept the
  # common spellings seen in distro, conda, vcpkg, and hand-written configs.
  find_package(FFTW3 CONFIG QUIET)
  if(FFTW3_FOUND)
    _cosmosim_select_first_existing_target(
      COSMOSIM_FFTW_TARGET
      FFTW3::fftw3
      FFTW3::fftw
      FFTW3::FFTW3
      FFTW::FFTW
      fftw3)
    _cosmosim_select_first_existing_target(
      COSMOSIM_FFTW_MPI_TARGET
      FFTW3::fftw3_mpi
      FFTW3::fftw_mpi
      FFTW3::FFTW3_MPI
      FFTW::MPI)
    _cosmosim_capture_fftw_version()
  endif()

  # Fall back to CosmoSim's module finder.  This covers Debian/Ubuntu packages,
  # where libfftw3-dev normally does not ship a CMake config package.
  if("${COSMOSIM_FFTW_TARGET}" STREQUAL "")
    find_package(FFTW3 QUIET MODULE)
    if(FFTW3_FOUND)
      _cosmosim_select_first_existing_target(
        COSMOSIM_FFTW_TARGET
        FFTW3::fftw3
        FFTW3::fftw
        FFTW3::FFTW3
        FFTW::FFTW)
      _cosmosim_select_first_existing_target(
        COSMOSIM_FFTW_MPI_TARGET
        FFTW3::fftw3_mpi
        FFTW3::fftw_mpi
        FFTW3::FFTW3_MPI
        FFTW::MPI)
      _cosmosim_capture_fftw_version()
    endif()
  endif()

  # Last resort: pkg-config.  Probe serial and MPI components independently so
  # serial PM builds do not require fftw3_mpi, while MPI+FFTW builds fail with a
  # precise missing-component diagnostic.
  find_package(PkgConfig QUIET)
  if(PkgConfig_FOUND)
    if("${COSMOSIM_FFTW_TARGET}" STREQUAL "")
      pkg_check_modules(FFTW_PKG QUIET IMPORTED_TARGET fftw3)
      if(FFTW_PKG_FOUND)
        set(COSMOSIM_FFTW_TARGET "PkgConfig::FFTW_PKG")
        _cosmosim_capture_fftw_version()
      endif()
    endif()

    if("${COSMOSIM_FFTW_MPI_TARGET}" STREQUAL "")
      pkg_check_modules(FFTW_MPI_PKG QUIET IMPORTED_TARGET fftw3_mpi)
      if(FFTW_MPI_PKG_FOUND)
        set(COSMOSIM_FFTW_MPI_TARGET "PkgConfig::FFTW_MPI_PKG")
        _cosmosim_capture_fftw_version()
      endif()
    endif()
  endif()

  if("${COSMOSIM_FFTW_TARGET}" STREQUAL "")
    message(FATAL_ERROR
      "COSMOSIM_ENABLE_FFTW=ON but FFTW3 serial double-precision library was not found.\n"
      "Tried CMake package config, CosmoSim's FindFFTW3 module, and pkg-config fallback.\n"
      "Install FFTW3 development packages and either provide FFTW3_DIR/FFTW3_ROOT or expose fftw3.pc via PKG_CONFIG_PATH.\n"
      "Debian/Ubuntu package: libfftw3-dev. Recommended validation preset: pm-hdf5-fftw-debug.")
  endif()

  if(COSMOSIM_FIND_FFTW_REQUIRE_MPI AND "${COSMOSIM_FFTW_MPI_TARGET}" STREQUAL "")
    message(FATAL_ERROR
      "COSMOSIM_ENABLE_MPI=ON with COSMOSIM_ENABLE_FFTW=ON requires FFTW MPI support, but fftw3_mpi was not found.\n"
      "Install the FFTW MPI development package and reconfigure.\n"
      "Debian/Ubuntu package: libfftw3-mpi-dev. If using a custom FFTW build, provide FFTW3_ROOT or expose fftw3_mpi.pc via PKG_CONFIG_PATH.\n"
      "Recommended validation preset: mpi-hdf5-fftw-debug.")
  endif()

  set(${out_target} "${COSMOSIM_FFTW_TARGET}" PARENT_SCOPE)
  set(COSMOSIM_FFTW_MPI_TARGET "${COSMOSIM_FFTW_MPI_TARGET}" PARENT_SCOPE)
  if(DEFINED FFTW_VERSION AND NOT "${FFTW_VERSION}" STREQUAL "")
    set(FFTW_VERSION "${FFTW_VERSION}" PARENT_SCOPE)
  endif()
endfunction()
