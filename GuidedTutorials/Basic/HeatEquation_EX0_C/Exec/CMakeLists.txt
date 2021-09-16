cmake_minimum_required(VERSION 3.16)

project(HeatEquation 
        LANGUAGES C CXX)


#add_executable(HeatEquation_EX0 main.cpp)

#target_link_libraries(HeatEquation_EX0 PRIVATE myfunc.H)

#target_include_directories(HeatEquation_EX0 PUBLIC ../Source)

#message( HeatEquation_EX0) 

add_subdirectory( /home/erik/AMReX/amrex-tutorials/GuidedTutorials/Basic/HeatEquation_EX0_C/Source/ )


## Project name and source file language
#project(HelloWorld 
#        LANGUAGES C CXX)
#
#add_executable(HelloWorld main.cpp)
#
## Add AMReX
#find_package(AMReX REQUIRED)
#target_link_libraries(HelloWorld PRIVATE AMReX::amrex) 
#
## Add mpi
#find_package(MPI REQUIRED)
#target_link_libraries(AMReX::amrex INTERFACE ${MPI_LIBRARIES}) 
#
## Copy input files
file( COPY ${CMAKE_CURRENT_SOURCE_DIR}/inputs DESTINATION . )
