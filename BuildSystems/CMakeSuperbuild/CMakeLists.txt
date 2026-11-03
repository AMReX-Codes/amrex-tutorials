# A standalone project that builds AMReX as part of itself (superbuild).
#
#   cmake -S . -B build [-DAMReX_GPU_BACKEND=CUDA -DCMAKE_CUDA_ARCHITECTURES=native]
#   cmake --build build -j 8
#
# To use a local AMReX source tree instead of downloading it:
#   cmake -S . -B build -DFETCHCONTENT_SOURCE_DIR_AMREX=/path/to/amrex
cmake_minimum_required(VERSION 3.25)
project(MyApp LANGUAGES C CXX)

# AMReX options must be set before AMReX is added. Users can still
# override them on the command line, e.g., -DAMReX_SPACEDIM=2.
set(AMReX_SPACEDIM    3    CACHE STRING "Dimension(s) of AMReX")
set(AMReX_MPI         ON   CACHE BOOL   "Enable MPI")
set(AMReX_OMP         OFF  CACHE BOOL   "Enable OpenMP")
set(AMReX_GPU_BACKEND NONE CACHE STRING "GPU backend: NONE, CUDA, HIP or SYCL")
set(AMReX_INSTALL     OFF  CACHE BOOL   "Generate AMReX install targets")

# The GPU language must be enabled before AMReX is added.
if(AMReX_GPU_BACKEND STREQUAL "CUDA")
    enable_language(CUDA)
elseif(AMReX_GPU_BACKEND STREQUAL "HIP")
    enable_language(HIP)
endif()

include(FetchContent)
FetchContent_Declare(amrex
    GIT_REPOSITORY https://github.com/AMReX-Codes/amrex.git
    GIT_TAG        development   # or a release tag such as 26.10
    GIT_SHALLOW    ON)
FetchContent_MakeAvailable(amrex)

add_executable(MyApp main.cpp)
target_link_libraries(MyApp PRIVATE AMReX::amrex)

# With CUDA, the .cpp files must be compiled as CUDA sources.
if(AMReX_GPU_BACKEND STREQUAL "CUDA")
    setup_target_for_cuda_compilation(MyApp)
endif()
