#
# Run the Single Vortex tutorial problem instrumented with 
# SENSEI either using the auto-instrumentation in the Amr class
# or using manual instrumentation which allows for particle data adaptors
#

set ( SRC_DIR Source )
set ( EXE_DIR Exec/SingleVortex )

#
# collect source files
#
set(_sources
  ${SRC_DIR}/Adv_F.H
  ${SRC_DIR}/AmrLevelAdv.cpp
  ${SRC_DIR}/AmrLevelAdv.H
  ${SRC_DIR}/main.cpp
  ${SRC_DIR}/LevelBldAdv.cpp

  ${SRC_DIR}/Src_${AMReX_SPACEDIM}d/Adv_${AMReX_SPACEDIM}d.f90
  ${SRC_DIR}/Src_${AMReX_SPACEDIM}d/compute_flux_${AMReX_SPACEDIM}d.f90
  ${SRC_DIR}/Src_${AMReX_SPACEDIM}d/slope_${AMReX_SPACEDIM}d.f90

  ${EXE_DIR}/Prob.f90 
  ${EXE_DIR}/face_velocity_${AMReX_SPACEDIM}d.f90

  ${SRC_DIR}/Src_nd/Adv_nd.f90
  ${SRC_DIR}/Src_nd/Tagging_nd.f90
  ${SRC_DIR}/Src_nd/tagging_params.f90
)

#
# Find input files
#
set(_inputs_files
  sensei
  probin
  inputs
  inputs.tracers
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
