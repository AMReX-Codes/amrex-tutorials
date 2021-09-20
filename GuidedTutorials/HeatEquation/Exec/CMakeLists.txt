# CMake file for build HeatEquation tutorial. 
#
# Step 1: Create and enter a build directory:
# 
#   mkdir build 
#   cd build
#
# Step 2 (AMReX Not Installed): Configure with CMake from inside the build directory:
#
#   cmake -DAMReX_SPACEDIM=2 ..       # Run with 2-dimensions
#
#         or:
#
#   cmake -DAMReX_SPACEDIM=3 ..       # Run with 3-dimensions
#
# Step 2 (AMReX Already Installed): Configure with CMake from inside the build directory:
#
#   cmake -DAMReX_ROOT=/path/to/installdir <compile options> ..    #AMReX installation determines AMReX_SPACEDIM value
#
# Step 3: Build the configuration:
#
#   cmake --build . -j4
#
# Done! You should now see an executable called "HelloWorld". To run:
#
#   mpiexec -n 4 ./HeatEquation inputs
#
# For additional CMake compile options see 
# https://amrex-codes.github.io/amrex/docs_html/BuildingAMReX.html#building-with-cmake

cmake_minimum_required(VERSION 3.16)

## Project name and source file languages
project(HeatEquation_EX0 
        LANGUAGES C CXX)

#add files from ../Source
add_executable(HeatEquation 
        ../Source/main.cpp
        ../Source/myfunc.H
        )

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
  
  FetchContent_Populate(amrex_code)
  
  # CMake will read the files in these directories and configure, build
  # and install AMReX.
  add_subdirectory(${amrex_code_SOURCE_DIR} ${amrex_code_BINARY_DIR})

else()

  # Add AMReX
  message("-- Searching for AMReX install directory at ${AMReX_ROOT}")
  find_package(AMReX REQUIRED PATHS ${AMReX_ROOT}/lib/cmake/AMReX/AMReXConfig.cmake)
  
  # Add mpi
  find_package(MPI REQUIRED)
  target_link_libraries(AMReX::amrex INTERFACE ${MPI_LIBRARIES}) 

endif()

target_link_libraries(HeatEquation PRIVATE AMReX::amrex) 

## Copy input files
file( COPY ${CMAKE_SOURCE_DIR}/inputs DESTINATION . )

