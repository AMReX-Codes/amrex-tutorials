cmake_minimum_required(VERSION 3.16)

## Project name and source file languages
project(HeatEquation_EX0 
        LANGUAGES C CXX)

#add files from ../Source
add_executable(HeatEquation 
        ../Source/main.cpp
        ../Source/myfunc.H
        )

## Add AMReX
find_package(AMReX REQUIRED)
target_link_libraries(HeatEquation PRIVATE AMReX::amrex) 

## Add mpi
find_package(MPI REQUIRED)
target_link_libraries(AMReX::amrex INTERFACE ${MPI_LIBRARIES}) 

## Copy input files
file( COPY ${CMAKE_SOURCE_DIR}/inputs DESTINATION . )

