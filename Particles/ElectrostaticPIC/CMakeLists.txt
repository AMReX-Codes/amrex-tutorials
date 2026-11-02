amrex_tutorial_add(ElectrostaticPIC
    SOURCES main.cpp
            FieldSolver.cpp FieldSolver.H
            PICParticleContainer.cpp PICParticleContainer.H
    INPUTS inputs
    REQUIRES PARTICLES LSOLVERS
    DIMS 2 3)
