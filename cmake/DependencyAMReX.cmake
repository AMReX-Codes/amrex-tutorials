# Make AMReX available, either by building it as part of this project
# (superbuild, the default) or by using an installed AMReX when AMReX_ROOT
# or AMReX_DIR is given.

set(AMReXTutorials_AMReX_REPO "https://github.com/AMReX-Codes/amrex.git"
    CACHE STRING "AMReX git repository")
set(AMReXTutorials_AMReX_TAG "development"
    CACHE STRING "AMReX git branch, tag or commit")

# Superbuild defaults. These only apply when AMReX is built here, and
# can be overridden on the command line, e.g., -DAMReX_FFT=OFF.
set(AMReX_EB             ON  CACHE BOOL "Build AMReX EB support")
set(AMReX_FFT            ON  CACHE BOOL "Build AMReX FFT support")
set(AMReX_PARTICLES      ON  CACHE BOOL "Build AMReX particle support")
set(AMReX_LINEAR_SOLVERS ON  CACHE BOOL "Build AMReX linear solvers")
set(AMReX_TINY_PROFILE   OFF CACHE BOOL "Build AMReX with TinyProfiler")
set(AMReX_INSTALL        OFF CACHE BOOL "Generate AMReX install targets")
set(AMReX_GPU_BACKEND   NONE CACHE STRING "GPU backend: NONE, CUDA, HIP or SYCL")

# The GPU language must be enabled in this scope, before AMReX is added.
if(AMReX_GPU_BACKEND STREQUAL "CUDA")
    enable_language(CUDA)
elseif(AMReX_GPU_BACKEND STREQUAL "HIP")
    enable_language(HIP)
endif()

# The Fortran tutorials need the AMReX Fortran interfaces. Fortran must be
# enabled here, in the top-level directory, for the tutorials to use it.
if(AMReX_FORTRAN_INTERFACES)
    enable_language(Fortran)
    set(AMReX_FORTRAN ON CACHE BOOL "Enable Fortran in AMReX")
endif()

set(_amrex_find_args)
if(DEFINED AMReX_ROOT OR DEFINED AMReX_DIR OR DEFINED ENV{AMReX_ROOT})
    set(_amrex_find_args FIND_PACKAGE_ARGS CONFIG REQUIRED)
endif()

include(FetchContent)

if(AMReX_SUNDIALS AND NOT _amrex_find_args)
    include(DependencySUNDIALS)
endif()

FetchContent_Declare(AMReX
    GIT_REPOSITORY ${AMReXTutorials_AMReX_REPO}
    GIT_TAG        ${AMReXTutorials_AMReX_TAG}
    GIT_SHALLOW    ON
    SYSTEM         # warnings in AMReX headers are not the tutorials' concern
    ${_amrex_find_args})
FetchContent_MakeAvailable(AMReX)

if(AMReX_FOUND)
    message(STATUS "Using installed AMReX from ${AMReX_DIR}")
    include(AMReXTargetHelpers)
    # The installed package names these AMReX_LSOLVERS and AMReX_FINTERFACES.
    set(AMReX_LINEAR_SOLVERS ${AMReX_LSOLVERS})
    if(AMReX_FORTRAN_INTERFACES AND NOT AMReX_FINTERFACES)
        message(FATAL_ERROR "The installed AMReX was built without Fortran interfaces")
    endif()
else()
    message(STATUS "Building AMReX from ${amrex_SOURCE_DIR}")
endif()

# An AMReX built for several dimensions (e.g., AMReX_SPACEDIM="2;3") can be
# used for any of them. The last one is the default.
list(GET AMReX_SPACEDIM -1 _amrex_last_dim)
set(AMReXTutorials_DIM ${_amrex_last_dim} CACHE STRING "Dimension of the tutorials")
if(NOT AMReXTutorials_DIM IN_LIST AMReX_SPACEDIM)
    message(FATAL_ERROR "AMReXTutorials_DIM=${AMReXTutorials_DIM} is not in "
                        "AMReX_SPACEDIM=${AMReX_SPACEDIM}")
endif()
message(STATUS "Tutorials are built for ${AMReXTutorials_DIM}D")
