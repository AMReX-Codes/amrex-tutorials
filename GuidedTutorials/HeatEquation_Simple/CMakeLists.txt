cmake_minimum_required(VERSION 3.16)

# Project name and source file languages
project(HeatEquation_Simple
        LANGUAGES C CXX)

add_executable(HeatEquation_Simple main.cpp)

# Add AMReX
find_package(AMReX REQUIRED)
target_link_libraries(HeatEquation_Simple PRIVATE AMReX::amrex)

# Add mpi
find_package(MPI REQUIRED)
target_link_libraries(AMReX::amrex INTERFACE ${MPI_LIBRARIES})

## Copy input files to build directory
file( COPY ${CMAKE_SOURCE_DIR}/inputs DESTINATION . )
