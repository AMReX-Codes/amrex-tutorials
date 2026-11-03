# Build SUNDIALS as part of this project for AMReX_SUNDIALS=ON, unless an
# installed SUNDIALS is given with SUNDIALS_ROOT or SUNDIALS_DIR.

set(AMReXTutorials_SUNDIALS_TAG "v7.9.0" CACHE STRING "SUNDIALS git tag")

if(DEFINED SUNDIALS_ROOT OR DEFINED SUNDIALS_DIR)
    return()  # AMReX will call find_package(SUNDIALS)
endif()

set(SUNDIALS_ENABLE_ARKODE          ON  CACHE INTERNAL "")
set(SUNDIALS_ENABLE_CVODE           ON  CACHE INTERNAL "")
set(SUNDIALS_ENABLE_CVODES          OFF CACHE INTERNAL "")
set(SUNDIALS_ENABLE_IDA             OFF CACHE INTERNAL "")
set(SUNDIALS_ENABLE_IDAS            OFF CACHE INTERNAL "")
set(SUNDIALS_ENABLE_KINSOL          OFF CACHE INTERNAL "")
set(SUNDIALS_ENABLE_C_EXAMPLES      OFF CACHE INTERNAL "")
set(SUNDIALS_ENABLE_CXX_EXAMPLES    OFF CACHE INTERNAL "")
set(SUNDIALS_ENABLE_EXAMPLES_INSTALL OFF CACHE INTERNAL "")
set(BUILD_SHARED_LIBS               OFF CACHE BOOL "")
# AMReX passes its MPI communicator to SUNDIALS, so their MPI settings must match.
if(NOT DEFINED AMReX_MPI OR AMReX_MPI)
    set(SUNDIALS_ENABLE_MPI ON  CACHE INTERNAL "")
else()
    set(SUNDIALS_ENABLE_MPI OFF CACHE INTERNAL "")
endif()
if(AMReX_GPU_BACKEND STREQUAL "CUDA")
    set(SUNDIALS_ENABLE_CUDA ON CACHE INTERNAL "")
elseif(AMReX_GPU_BACKEND STREQUAL "HIP")
    set(SUNDIALS_ENABLE_HIP ON CACHE INTERNAL "")
elseif(AMReX_GPU_BACKEND STREQUAL "SYCL")
    set(SUNDIALS_ENABLE_SYCL ON CACHE INTERNAL "")
endif()

FetchContent_Declare(sundials
    GIT_REPOSITORY https://github.com/LLNL/sundials.git
    GIT_TAG        ${AMReXTutorials_SUNDIALS_TAG}
    GIT_SHALLOW    ON
    SYSTEM)
FetchContent_MakeAvailable(sundials)

# Tell AMReX that SUNDIALS (with its SUNDIALS:: targets) is already available,
# and which of its libraries to link.
set(SUNDIALS_FOUND TRUE)
set(SUNDIALS_COMPONENTS arkode cvode sunlinsolspgmr sunlinsolspfgmr
    nvecserial nvecmanyvector sunnonlinsolfixedpoint)
string(REGEX REPLACE "^v" "" SUNDIALS_VERSION ${AMReXTutorials_SUNDIALS_TAG})
