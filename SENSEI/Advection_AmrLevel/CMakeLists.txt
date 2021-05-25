#
# Single Vortex Tutorial Instrumented with SENSEI
#
# | Mini-app                | Description |
# | ----                    | ----        |
# | ExplicitAmrAndParticles | Process both mesh and particle data |
# | ImplicitAmr             | Processes mesh based data using built in instrumentation |
# | ExplicitAmr             | Processes mesh based data using explicit intrumentation |
#

set ( SRC_DIR Source )
set ( EXE_DIR Exec/SingleVortex )

#
# Tutorial executables
#
set ( _mini_apps )

if ( AMReX_PARTICLES AND AMReX_NO_SENSEI_AMR_INST )
    message ( STATUS "Tutorial SENSEI ExplicitAmrAndParticles ... enabled")
    list ( APPEND _mini_apps ExplicitAmrAndParticles )
else ()
    message ( STATUS "Tutorial SENSEI ExplicitAmrAndParticles ... disabled")
endif ()

if ( NOT AMReX_NO_SENSEI_AMR_INST )
    message ( STATUS "Tutorial SENSEI ImplcitAmr ................ enabled")
    message ( STATUS "Tutorial SENSEI ExplicitAmr ............... disabled")
    list ( APPEND _mini_apps ImplicitAmr )
else ()
    message ( STATUS "Tutorial SENSEI ImplcitAmr ................ disabled")
    message ( STATUS "Tutorial SENSEI ExplicitAmr ............... enabled")
    list ( APPEND _mini_apps ExplicitAmr )
endif ()

#
# Find input files
#
set ( _inputs_files sensei probin inputs inputs.tracers )

list(TRANSFORM _inputs_files PREPEND "${EXE_DIR}/")

foreach ( _app ${_mini_apps} )

    #
    # collect source files
    #
    set(_sources
      ${SRC_DIR}/Adv_F.H
      ${SRC_DIR}/AmrLevelAdv.cpp
      ${SRC_DIR}/AmrLevelAdv.H
      ${SRC_DIR}/LevelBldAdv.cpp
      ${SRC_DIR}/main${_app}.cpp

      ${SRC_DIR}/Src_${AMReX_SPACEDIM}d/Adv_${AMReX_SPACEDIM}d.f90
      ${SRC_DIR}/Src_${AMReX_SPACEDIM}d/compute_flux_${AMReX_SPACEDIM}d.f90
      ${SRC_DIR}/Src_${AMReX_SPACEDIM}d/slope_${AMReX_SPACEDIM}d.f90

      ${EXE_DIR}/Prob.f90
      ${EXE_DIR}/face_velocity_${AMReX_SPACEDIM}d.f90

      ${SRC_DIR}/Src_nd/Adv_nd.f90
      ${SRC_DIR}/Src_nd/Tagging_nd.f90
      ${SRC_DIR}/Src_nd/tagging_params.f90 )

    #
    # Add to the "tutorial" target
    #
    setup_tutorial( _sources _inputs_files
        BASE_NAME Advection_AmrLevel_${_app} HAS_FORTRAN_MODULES )

    unset( _sources )
    unset( _inputs_files )

endforeach()
