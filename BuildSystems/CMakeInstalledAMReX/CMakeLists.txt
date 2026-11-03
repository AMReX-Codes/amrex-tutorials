# A standalone project that uses an installed AMReX.
#
# First build and install AMReX, e.g.,
#   cmake -S amrex -B amrex/build -DAMReX_SPACEDIM=3 -DCMAKE_INSTALL_PREFIX=/path/to/install
#   cmake --build amrex/build -j 8 --target install
# then
#   cmake -S . -B build -DAMReX_ROOT=/path/to/install
#   cmake --build build -j 8
cmake_minimum_required(VERSION 3.25)
project(MyApp LANGUAGES C CXX)

# Request the dimension and any optional components the code needs, e.g.,
# find_package(AMReX REQUIRED 3D EB PARTICLES).
find_package(AMReX REQUIRED 3D)

# An AMReX built with a GPU backend needs that language enabled here.
if(AMReX_GPU_BACKEND STREQUAL "CUDA")
    enable_language(CUDA)
    include(AMReXTargetHelpers)
elseif(AMReX_GPU_BACKEND STREQUAL "HIP")
    enable_language(HIP)
endif()

add_executable(MyApp main.cpp)
target_link_libraries(MyApp PRIVATE AMReX::amrex_3d)

if(AMReX_GPU_BACKEND STREQUAL "CUDA")
    setup_target_for_cuda_compilation(MyApp)
endif()
