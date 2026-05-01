include(FindPackageHandleStandardArgs)

find_package(PkgConfig QUIET)
if(PkgConfig_FOUND)
  pkg_check_modules(FFTW3_PC QUIET fftw3)
  pkg_check_modules(FFTW3_MPI_PC QUIET fftw3_mpi)
endif()

find_path(FFTW3_INCLUDE_DIR
  NAMES fftw3.h
  HINTS
    ${FFTW3_ROOT}
    $ENV{FFTW3_ROOT}
    ${FFTW3_PC_INCLUDE_DIRS}
  PATH_SUFFIXES include)

find_library(FFTW3_LIBRARY
  NAMES fftw3 libfftw3-3
  HINTS
    ${FFTW3_ROOT}
    $ENV{FFTW3_ROOT}
    ${FFTW3_PC_LIBRARY_DIRS}
  PATH_SUFFIXES lib lib64)

find_path(FFTW3_MPI_INCLUDE_DIR
  NAMES fftw3-mpi.h
  HINTS
    ${FFTW3_ROOT}
    $ENV{FFTW3_ROOT}
    ${FFTW3_MPI_PC_INCLUDE_DIRS}
    ${FFTW3_INCLUDE_DIR}
  PATH_SUFFIXES include)

find_library(FFTW3_MPI_LIBRARY
  NAMES fftw3_mpi libfftw3_mpi
  HINTS
    ${FFTW3_ROOT}
    $ENV{FFTW3_ROOT}
    ${FFTW3_MPI_PC_LIBRARY_DIRS}
  PATH_SUFFIXES lib lib64)

set(FFTW3_LIBRARIES "${FFTW3_LIBRARY}")
set(FFTW3_INCLUDE_DIRS "${FFTW3_INCLUDE_DIR}")
set(FFTW3_MPI_LIBRARIES "${FFTW3_MPI_LIBRARY}")
set(FFTW3_MPI_INCLUDE_DIRS "${FFTW3_MPI_INCLUDE_DIR}")

if(FFTW3_PC_VERSION)
  set(FFTW3_VERSION "${FFTW3_PC_VERSION}")
elseif(FFTW3_MPI_PC_VERSION)
  set(FFTW3_VERSION "${FFTW3_MPI_PC_VERSION}")
endif()

find_package_handle_standard_args(FFTW3
  REQUIRED_VARS FFTW3_LIBRARY FFTW3_INCLUDE_DIR
  VERSION_VAR FFTW3_VERSION)

if(FFTW3_FOUND AND NOT TARGET FFTW3::fftw3)
  add_library(FFTW3::fftw3 UNKNOWN IMPORTED)
  set_target_properties(FFTW3::fftw3 PROPERTIES
    IMPORTED_LOCATION "${FFTW3_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${FFTW3_INCLUDE_DIR}")
endif()

# fftw3_mpi is intentionally optional for serial FFTW builds.  When present,
# expose it as a separate target and attach the same include directory used for
# fftw3-mpi.h.  Callers that require distributed PM should request REQUIRE_MPI
# through cosmosim_find_fftw(), which turns absence into a clear configure error.
if(FFTW3_FOUND AND FFTW3_MPI_LIBRARY AND FFTW3_MPI_INCLUDE_DIR AND NOT TARGET FFTW3::fftw3_mpi)
  add_library(FFTW3::fftw3_mpi UNKNOWN IMPORTED)
  set_target_properties(FFTW3::fftw3_mpi PROPERTIES
    IMPORTED_LOCATION "${FFTW3_MPI_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${FFTW3_MPI_INCLUDE_DIR}"
    INTERFACE_LINK_LIBRARIES FFTW3::fftw3)
endif()

mark_as_advanced(
  FFTW3_INCLUDE_DIR
  FFTW3_LIBRARY
  FFTW3_MPI_INCLUDE_DIR
  FFTW3_MPI_LIBRARY)
