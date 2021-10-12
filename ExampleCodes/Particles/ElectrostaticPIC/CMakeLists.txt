if (NOT AMReX_SPACEDIM EQUAL 2)
  message("-----------------------------------------------------------------------------")
  message([WARNING] "ElectrostaticPIC Tutorial only avaiable in 2D, set AMReX_SPACEDIM=2")
  message("-----------------------------------------------------------------------------")
  return()
endif ()

set(_sources main.cpp)
list(APPEND _sources PhysConst.H)

list(APPEND _sources particles/ElectrostaticParticleContainer.cpp)
list(APPEND _sources particles/ElectrostaticParticleContainer.H)
list(APPEND _sources particles/deposition/ChargeDeposition_K.H)
list(APPEND _sources particles/field_gather/FieldGather_K.H)
list(APPEND _sources particles/pusher/ParticlePusher_K.H)

list(APPEND _sources field_solver/FieldSolver.H)
list(APPEND _sources field_solver/FieldSolver.cpp)
list(APPEND _sources field_solver/FieldSolver_K.H)

list(APPEND _sources diagnostics/FieldIO.H)
list(APPEND _sources diagnostics/FieldIO.cpp)
list(APPEND _sources diagnostics/ParticleIO.cpp)

list(TRANSFORM _sources PREPEND src/)

set(_input_files example/inputs)

setup_tutorial(_sources _input_files)

unset(_sources)
unset(_input_files)
