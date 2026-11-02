# IncFlow builds one executable per problem in Exec/: IncFlow_TaylorGreen
# and IncFlow_LidDrivenCavity.
#
# Inside the tutorials superbuild, amrex_tutorial_add creates the targets
# Applications_IncFlow_TaylorGreen and Applications_IncFlow_LidDrivenCavity.
# When this directory is copied out and built on its own, AMReX is fetched
# with FetchContent:
#
#   cmake -S . -B build [-DAMReX_SPACEDIM=2] [-DAMReX_GPU_BACKEND=CUDA]
#         [-DFETCHCONTENT_SOURCE_DIR_AMREX=/path/to/amrex]
#   cmake --build build -j 8
if(NOT COMMAND amrex_tutorial_add)
    cmake_minimum_required(VERSION 3.25)
    project(IncFlow LANGUAGES C CXX)

    if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
        set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    endif()

    # AMReX options must be set before AMReX is added.
    set(AMReX_SPACEDIM       3    CACHE STRING "Dimension of AMReX")
    set(AMReX_MPI            ON   CACHE BOOL   "Enable MPI")
    set(AMReX_OMP            OFF  CACHE BOOL   "Enable OpenMP")
    set(AMReX_GPU_BACKEND    NONE CACHE STRING "GPU backend: NONE, CUDA, HIP or SYCL")
    set(AMReX_LINEAR_SOLVERS ON   CACHE BOOL   "Build the linear solvers")
    set(AMReX_TINY_PROFILE   OFF  CACHE BOOL   "Enable TinyProfiler")
    set(AMReX_INSTALL        OFF  CACHE BOOL   "Generate AMReX install targets")

    if(AMReX_GPU_BACKEND STREQUAL "CUDA")
        enable_language(CUDA)
    elseif(AMReX_GPU_BACKEND STREQUAL "HIP")
        enable_language(HIP)
    endif()

    include(FetchContent)
    FetchContent_Declare(amrex
        GIT_REPOSITORY https://github.com/AMReX-Codes/amrex.git
        GIT_TAG        development
        GIT_SHALLOW    ON)
    FetchContent_MakeAvailable(amrex)

    set(_incflow_standalone ON)
endif()

set(_incflow_sources
    Source/main.cpp
    Source/IncFlow.cpp
    Source/Advance.cpp
    Source/Projection.cpp
    Source/IO.cpp)

foreach(_prob IN ITEMS TaylorGreen DoubleShearLayer LidDrivenCavity)
    set(_sources ${_incflow_sources} Exec/${_prob}/prob.cpp)
    if(_incflow_standalone)
        set(_target IncFlow_${_prob})
        add_executable(${_target} ${_sources})
        target_link_libraries(${_target} PRIVATE AMReX::amrex)
        if(AMReX_GPU_BACKEND STREQUAL "CUDA")
            setup_target_for_cuda_compilation(${_target})
        endif()
    else()
        set(_target Applications_IncFlow_${_prob})
        # The inputs are read from the source tree because both problems
        # have a file named inputs.
        amrex_tutorial_add(IncFlow_${_prob}
            SOURCES ${_sources}
            REQUIRES LSOLVERS
            DIMS 2 3
            TEST_ARGS ${CMAKE_CURRENT_SOURCE_DIR}/Exec/${_prob}/inputs
                      amr.plot_file=${_prob}_plt amr.chk_file=${_prob}_chk)
    endif()
    if(TARGET ${_target})
        target_include_directories(${_target} PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/Source
            ${CMAKE_CURRENT_SOURCE_DIR}/Exec/${_prob})
    endif()
endforeach()
