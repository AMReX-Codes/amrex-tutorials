#!/usr/bin/env bash
#
# Build a sample of tutorials with GNU Make and run them.
# Expects AMReX to be cloned next to this repository.

set -eu -o pipefail
cd "$(git rev-parse --show-toplevel)"

DIRS="
Base/HelloWorld
Base/HeatEquation
FFT/Poisson
LinearSolvers/CellPoisson
AMR/AdvectionAmrCore
EB/Poisson
Particles/Basics
Fortran/HeatEquation
Applications/CNS/Exec/Sod
Applications/IncFlow/Exec/LidDrivenCavity
BuildSystems/GNUMake/Exec
"

for dir in ${DIRS}; do
    echo "::group::${dir}"
    (
        cd "${dir}"
        make -j 4 USE_OMP=TRUE
        exe=$(ls ./*.ex | head -n 1)
        if [ -f inputs ]; then
            OMP_NUM_THREADS=2 mpiexec --oversubscribe -n 2 "${exe}" inputs
        else
            OMP_NUM_THREADS=2 mpiexec --oversubscribe -n 2 "${exe}"
        fi
        make realclean
    )
    echo "::endgroup::"
done
