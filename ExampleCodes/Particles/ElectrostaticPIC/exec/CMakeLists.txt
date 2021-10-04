if (AMReX_SPACEDIM EQUAL 1 OR NOT CMAKE_Fortran_COMPILER_LOADED)
   return()
endif ()

set(_sources ElectrostaticParticleContainer.cpp ElectrostaticParticleContainer.H Electrostatic_PIC_K.H Electrostatic_PIC_Util.H Electrostatic_PIC_Util.cpp
    main.cpp PhysConst.H Particle_Pusher.H)

set(_input_files inputs)

setup_tutorial(_sources _input_files)

unset(_sources)
unset(_input_files)
