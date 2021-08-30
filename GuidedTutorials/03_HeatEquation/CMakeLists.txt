if (AMReX_SPACEDIM EQUAL 1)
   return()
endif ()

# List of source files
set(_sources main.cpp myfunc.cpp myfunc.H mykernel.H)
set(_input_files inputs Visualization.ipynb)

setup_tutorial(_sources _input_files)

unset( _sources )
unset( _input_files   )
