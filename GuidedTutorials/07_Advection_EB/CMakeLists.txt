if (AMReX_SPACEDIM EQUAL 1)
   return()
endif ()

# List of source files
set(_sources EB_Cylinder.cpp  FluidParticleContainer.cpp  FluidParticleContainer.H)
list(APPEND _sources DefineVelocity.cpp  face_velocity.H  main.cpp  Indexing.H)
list(APPEND _sources mac_project_velocity.cpp)

# List of input files
set(_input_files inputs Visualization.ipynb)

setup_tutorial(_sources _input_files)

unset( _sources )
unset( _input_files )
