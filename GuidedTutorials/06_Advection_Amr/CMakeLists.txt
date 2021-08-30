if (AMReX_SPACEDIM EQUAL 1)
   return()
endif ()

# List of source files
set(_sources AdvancePhiAllLevels.cpp  AdvancePhiAtLevel.cpp  AmrCoreAdv.cpp  AmrCoreAdv.H bc_fill.H)
list(APPEND _sources DefineVelocity.cpp  face_velocity.H  Kernels.H  main.cpp  Tagging.H)
list(APPEND _sources Adv_K.H  compute_flux_${AMReX_SPACEDIM}D_K.H slope_K.H)
list(APPEND _sources Prob.H)

# List of input files
set(_input_files inputs Visualization.ipynb)

setup_tutorial(_sources _input_files)

unset( _sources )
unset( _input_files )
