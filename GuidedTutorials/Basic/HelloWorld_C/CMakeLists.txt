cmake_minimum_required(VERSION 3.16)

# Project name and source file language
project(HelloWorld 
        LANGUAGES C CXX)

add_executable(HelloWorld main.cpp)

# Add AMReX
find_package(AMReX REQUIRED)
target_link_libraries(HelloWorld PRIVATE AMReX::amrex) 

# Add mpi
find_package(MPI REQUIRED)
target_link_libraries(AMReX::amrex INTERFACE ${MPI_LIBRARIES}) 
