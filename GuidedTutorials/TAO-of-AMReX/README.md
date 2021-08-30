# Boundary Control with 2-D Laplace Equation

## Introduction

ATPESC 2019 Numerical Software Track

[__Hands-on Lesson Page__](https://xsdk-project.github.io/MathPackagesTraining/lessons/boundary_control_tao/)

__Developed by:__

* Alp Dener, Postdoctoral Researcher, Argonne National Laboratory

* Donald E. Willcox, Postdoctoral Researcher. Lawrence-Livermore National Laboratory

## Dependencies

1. [PETSc/TAO](https://www.mcs.anl.gov/petsc/)
    - Must be configured with MPI
    - `PETSC_DIR` and `PETSC_ARCH` environment variables must be set

2. [AMReX](https://amrex-codes.github.io/amrex/)
    - Must be compiled as a library with `DIM=2` (2-D configuration)
    - Modify problem makefile to provide the correct `ATPESC_HOME` directory

3. [MPICH](https://www.mpich.org) or [OpenMPI](https://www.open-mpi.org)
    - PETSc can provide an MPICH installation if configured with `--download-mpich`

