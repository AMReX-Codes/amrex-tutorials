# CMake file for build MultiFab tutorial.
#
# Step 1: Create and enter a build directory:
#
#   mkdir build
#   cd build
#
# Step 2: Configure with CMake from inside the build directory:
#
#   cmake <compile options> ..                   #if AMReX not installed
#
#         or:
#
#   cmake -DAMReX_ROOT=/path/to/installdir <compile options> ..    #if AMReX already installed
#
# Step 3: Build the configuration:
#
#   cmake --build . -j4
#
# Done! You should now see an executable called "MultiFab". To run:
#
#   mpiexec -n 4 ./MultiFab
#
# For additional CMake compile options see
# https://amrex-codes.github.io/amrex/docs_html/BuildingAMReX.html#building-with-cmake



cmake_minimum_required(VERSION 3.16)

# Project name and source file language
project(MultiFab
        LANGUAGES C CXX)

add_executable(MultiFab main.cpp)



# To use a pre-installed AMReX build, run:
#    cmake -DAMReX_ROOT=/path/to/installdir
# Otherwise cmake will download AMReX from GitHub

if(NOT DEFINED AMReX_ROOT)
  message("-- Download and configure AMReX from GitHub")

  #Download AMReX from GitHub
  include(FetchContent)
  set(FETCHCONTENT_QUIET OFF)

  FetchContent_Declare(
    amrex_code
    GIT_REPOSITORY https://github.com/AMReX-Codes/amrex.git/
    GIT_TAG        origin/development
    )

  FetchContent_MakeAvailable(amrex_code)

else()

  # Add AMReX
  message("-- Searching for AMReX install directory at ${AMReX_ROOT}")
  find_package(AMReX PATHS ${AMReX_ROOT}/lib/cmake/AMReX/AMReXConfig.cmake)

endif()

target_link_libraries(MultiFab PRIVATE AMReX::amrex)
