cmake_minimum_required(VERSION 3.16)

# Project name and source file language
project(HelloWorld 
        LANGUAGES C CXX)

add_executable(HelloWorld main.cpp)







include(FetchContent)
set(FETCHCONTENT_QUIET OFF)


FetchContent_Declare(
  amrex_code
  GIT_REPOSITORY https://github.com/AMReX-Codes/amrex.git/
  GIT_TAG        origin/development
  CONFIGURE_COMMAND "-DAMREX_HOME=${amrex_code_SOURCE_DIR}"
  #BUILD_COMMAND "cmake --build . -j4"
  #INSTALL_COMMAND "make -j4 install"
  )

FetchContent_Populate(amrex_code)

add_subdirectory(${amrex_code_SOURCE_DIR} ${amrex_code_BINARY_DIR})




message("amrex_SOURCE_DIR = ${amrex_code_SOURCE_DIR}")
message("amrex_BINARY_DIR = ${amrex_code_BINARY_DIR}")






# Add AMReX
#nd_package(AMReX PATHS ${amrex_BINARY_DIR}/lib/cmake/AMReX/AMReXConfig.cmake)
target_link_libraries(HelloWorld PRIVATE AMReX::amrex) 

# Add mpi
find_package(MPI REQUIRED)
target_link_libraries(AMReX::amrex INTERFACE ${MPI_LIBRARIES}) 
