cmake_minimum_required(VERSION 3.14)

project( AMReX-Tutorials
   DESCRIPTION "Tutorials for the AMReX adaptive mesh refinement framework"
   #VERSION  ${AMREX_PKG_VERSION}
   # Check that the line below points to the correct repo
   HOMEPAGE_URL "https://amrex-codes.github.io/amrex/tutorials_html/index.html"
   )

set( CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake )

#
# enable building against an existing AMReX library.
# Default to downloading and rebuilding the library here
#
option( AMReX_USE_LOCAL_BUILD "Use an existing build of the AMReX library" OFF )

if( NOT ${AMReX_USE_LOCAL_BUILD} )
   #
   # Fetch amrex repo
   #
   set(AMReX_GIT_BRANCH "development" CACHE STRING "The AMReX branch to checkout")
   set(AMReX_INSTALL  "NO" CACHE INTERNAL "Disable install target for amrex")

   include(FetchContent)
   set(FETCHCONTENT_QUIET OFF)  # Verbose ON

   FetchContent_Declare( amrex
       GIT_REPOSITORY https://github.com/AMReX-Codes/amrex.git/
       GIT_TAG        ${AMReX_GIT_BRANCH}
       )

   if(NOT ${amrex}_POPULATED)
       FetchContent_Populate(amrex)

       list(APPEND CMAKE_MODULE_PATH ${amrex_SOURCE_DIR}/Tools/CMake)

       # Load amrex options here so that they are
       # available to the entire project
       include(AMReXOptions)

       if (AMReX_FORTRAN)
          enable_language(Fortran)
       endif ()

       if (AMReX_GPU_BACKEND STREQUAL "CUDA")
          enable_language(CUDA)
          include(AMReX_SetupCUDA)
       endif ()

       # Bring the populated content into the build
       add_subdirectory(${amrex_SOURCE_DIR} ${amrex_BINARY_DIR})
   endif()
else()
   message(STATUS "Using existing AMReX library")

   enable_language(Fortran)
   find_package(AMReX REQUIRED)

   if (AMReX_GPU_BACKEND STREQUAL "CUDA")
      enable_language(CUDA)
      include(AMReX_SetupCUDA)
   endif()
endif()

#
# List of subdirectories to search for CMakeLists.
# For now, we do not include MUI, SDC, SWFFT
#
set( AMREX_TUTORIALS_SUBDIRS Amr Basic ForkJoin)

if (AMReX_SENSEI)
   list(APPEND AMREX_TUTORIALS_SUBDIRS SENSEI)
endif()

if (AMReX_ASCENT)
   list(APPEND AMREX_TUTORIALS_SUBDIRS Blueprint)
endif ()

if (AMReX_EB)
   list(APPEND AMREX_TUTORIALS_SUBDIRS EB)
endif ()

if (AMReX_LINEAR_SOLVERS)
   list(APPEND AMREX_TUTORIALS_SUBDIRS LinearSolvers)
endif ()

if (AMReX_PARTICLES)
   list(APPEND AMREX_TUTORIALS_SUBDIRS Particles)
endif ()

if (AMReX_FORTRAN_INTERFACES)
   list(APPEND AMREX_TUTORIALS_SUBDIRS FortranInterface)
endif ()

if (NOT (AMReX_GPU_BACKEND STREQUAL "NONE") )
   list(APPEND AMREX_TUTORIALS_SUBDIRS GPU)
endif ()

list(TRANSFORM AMREX_TUTORIALS_SUBDIRS PREPEND "${CMAKE_CURRENT_LIST_DIR}/")

#
# Search for CMakelists.txt in the subdirectories defined in the list above
# and include those that contain one into the build
#
include(SetupTutorials)
foreach (_subdir IN LISTS AMREX_TUTORIALS_SUBDIRS)
   file( GLOB_RECURSE _tests "${_subdir}/*CMakeLists.txt" )
   foreach ( _item  IN LISTS _tests)
      get_filename_component(_dir ${_item} DIRECTORY )
      add_subdirectory(${_dir})
   endforeach ()
endforeach ()
