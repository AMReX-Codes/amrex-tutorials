if (AMReX_SPACEDIM EQUAL 1)
   return()
endif ()

# List of source files
set(_sources main.cpp  MyParticleContainer.cpp  MyParticleContainer.H)

# List of input files
set(_input_files inputs_3d initial_particles_3d paraview_pachinko.py Visualization.ipynb)

setup_tutorial(_sources _input_files)

unset( _sources )
unset( _input_files )
