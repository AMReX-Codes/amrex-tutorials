# GNUMake

A standalone project built with AMReX's GNU Make system. The application
sources are in `Source/` and the build happens in `Exec/`, a layout that
scales to codes with several problem setups (`Exec/ProblemA`, ...).

## Build and run

```
cd Exec
make -j8 AMREX_HOME=/path/to/amrex
mpiexec -n 2 ./main3d.gnu.MPI.ex
```

## Options

Options are set in `Exec/GNUmakefile` or on the command line:

| Variable   | Values                    |
|------------|---------------------------|
| `DIM`      | 1, 2, 3                   |
| `COMP`     | gnu, llvm, intel-llvm     |
| `DEBUG`    | TRUE, FALSE               |
| `USE_MPI`, `USE_OMP` | TRUE, FALSE     |
| `USE_CUDA`, `USE_HIP`, `USE_SYCL` | TRUE, FALSE (e.g., `USE_CUDA=TRUE CUDA_ARCH=80`) |

`Source/Make.package` lists the application sources. `AMREX_DIRS` in the
`GNUmakefile` lists the AMReX source directories to compile.
