# Fortran tutorials use the AMReX Fortran interfaces and run on CPUs only.
# They are built with -DAMReX_FORTRAN_INTERFACES=ON.
if(NOT AMReX_FORTRAN_INTERFACES)
    message(STATUS "Skipping Fortran tutorials: AMReX_FORTRAN_INTERFACES is OFF")
    return()
endif()

add_subdirectory(HeatEquation)
