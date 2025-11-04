#!/usr/bin/env bash
#
# Copyright 2020 The AMReX Community
#
# License: BSD-3-Clause-LBNL
# Authors: Axel Huebl

set -eu -o pipefail

brew update
brew install gfortran || true

# verify installation
gfortran-14 --version
otool -L $(which gfortran-14)

# make sure to install Open MPI with the correct Fortran compiler
export FC=$(which gfortran-14)
export F77=$FC
export F90=$FC

brew install libomp || true
brew install open-mpi || true
brew install ccache || true
