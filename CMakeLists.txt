# AMReX tutorials: superbuild of all tutorials.
#
#   cmake -S . -B build [-DAMReX_GPU_BACKEND=CUDA] [-DAMReXTutorials_TOPICS="Base;FFT"]
#   cmake --build build -j 8
#   ctest --test-dir build
#
# AMReX is fetched from GitHub by default. To use something else:
#   local source:   -DFETCHCONTENT_SOURCE_DIR_AMREX=/path/to/amrex
#   installed:      -DAMReX_ROOT=/path/to/install
#
# Optional dependencies: -DAMReX_SUNDIALS=ON (SUNDIALS is fetched unless
# SUNDIALS_ROOT is given), -DAMReX_SIMD=ON (vir-simd is fetched unless
# installed), -DAMReX_FORTRAN_INTERFACES=ON (needs a Fortran compiler),
# -DAMReXTutorials_PYTHON=ON, -DAMReXTutorials_TORCH=ON.
cmake_minimum_required(VERSION 3.25)

project(AMReXTutorials
    DESCRIPTION "Tutorials for the AMReX block-structured AMR framework"
    HOMEPAGE_URL "https://amrex-codes.github.io/amrex/tutorials_html/"
    LANGUAGES C CXX)

if(CMAKE_BINARY_DIR STREQUAL CMAKE_SOURCE_DIR)
    message(FATAL_ERROR "In-source builds are not allowed.")
endif()

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

set(AMReXTutorials_ALL_TOPICS
    Base AMR LinearSolvers EB Particles FFT TimeIntegration
    Coupling Fortran Python ML Applications)
set(AMReXTutorials_TOPICS "${AMReXTutorials_ALL_TOPICS}"
    CACHE STRING "Semicolon-separated list of tutorial topics to build")

# Tutorials with extra dependencies. SUNDIALS is enabled with -DAMReX_SUNDIALS=ON
# and Fortran with -DAMReX_FORTRAN_INTERFACES=ON.
option(AMReXTutorials_PYTHON "Test the pyAMReX tutorials (needs the amrex Python module)" OFF)
option(AMReXTutorials_TORCH  "Build the PyTorch tutorials (needs LibTorch)" OFF)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
include(DependencyAMReX)
include(AMReXTutorials)

enable_testing()

foreach(_topic IN LISTS AMReXTutorials_TOPICS)
    if(NOT _topic IN_LIST AMReXTutorials_ALL_TOPICS)
        message(FATAL_ERROR "Unknown tutorial topic '${_topic}'. "
                            "Valid topics: ${AMReXTutorials_ALL_TOPICS}")
    endif()
    if(EXISTS ${PROJECT_SOURCE_DIR}/${_topic}/CMakeLists.txt)
        add_subdirectory(${_topic})
    endif()
endforeach()
