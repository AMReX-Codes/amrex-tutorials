# Make vir-simd, the header-only library behind AMReX_SIMD=ON, available for
# AMReX's find_package(vir-simd), fetching it unless it is installed.

set(AMReXTutorials_VIR_SIMD_TAG "v0.4.4" CACHE STRING "vir-simd git tag")

find_package(vir-simd QUIET)
if(vir-simd_FOUND)
    return()
endif()

FetchContent_Declare(vir-simd
    GIT_REPOSITORY https://github.com/mattkretz/vir-simd.git
    GIT_TAG        ${AMReXTutorials_VIR_SIMD_TAG}
    GIT_SHALLOW    ON
    SYSTEM
    OVERRIDE_FIND_PACKAGE)
FetchContent_MakeAvailable(vir-simd)

# The installed package exports vir-simd::vir-simd; the source tree only
# defines vir-simd. AMReX links it and exports its own targets from the
# build tree, which requires vir-simd to be exported too.
if(NOT TARGET vir-simd::vir-simd)
    add_library(vir-simd::vir-simd ALIAS vir-simd)
endif()
export(TARGETS vir-simd NAMESPACE vir-simd::
       FILE ${CMAKE_CURRENT_BINARY_DIR}/vir-simd-targets.cmake)
