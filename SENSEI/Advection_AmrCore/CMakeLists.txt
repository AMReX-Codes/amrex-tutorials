#
# Run the Single Vortex tutorial problem instrumented with 
# SENSEI either using manual instrumentation of an AmrMesh
#

#
# Does not work if amrex is built in shared mode
#
if (BUILD_SHARED_LIBS)
  return()
endif ()

set ( SRC_DIR Source ) 
set ( EXE_DIR Exec/SingleVortex ) 

#
# collect source files
#
set( _sources
  ${SRC_DIR}/AmrCoreAdv.cpp
  ${SRC_DIR}/bc_fill_nd.F90
  ${SRC_DIR}/main.cpp
  ${SRC_DIR}/AmrCoreAdv_F.H
  ${SRC_DIR}/AmrCoreAdv.H

  ${SRC_DIR}/Src_${AMReX_SPACEDIM}d/Adv_${AMReX_SPACEDIM}d.f90
  ${SRC_DIR}/Src_${AMReX_SPACEDIM}d/compute_flux_${AMReX_SPACEDIM}d.f90
  ${SRC_DIR}/Src_${AMReX_SPACEDIM}d/slope_${AMReX_SPACEDIM}d.f90
  ${EXE_DIR}/Prob.f90 
  ${EXE_DIR}/face_velocity_${AMReX_SPACEDIM}d.f90

  ${SRC_DIR}/Src_nd/Tagging_nd.f90
)

#
# Find input files
#
set(_inputs_files
  sensei
  inputs
)
list(TRANSFORM _inputs_files PREPEND "${EXE_DIR}/")

#
# Add to the "tutorial" target
#
setup_tutorial(
  _sources _inputs_files
  HAS_FORTRAN_MODULES
)
unset( _sources )
unset( _inputs_files )