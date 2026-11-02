# CMakeInstalledAMReX

A standalone CMake project that links to an AMReX installed with CMake,
using `find_package(AMReX)`.

## Build and run

Install AMReX once:

```
cmake -S /path/to/amrex -B amrex-build -DAMReX_SPACEDIM=3 -DCMAKE_INSTALL_PREFIX=$HOME/amrex-install
cmake --build amrex-build -j 8 --target install
```

Then build this project:

```
cmake -S . -B build -DAMReX_ROOT=$HOME/amrex-install
cmake --build build -j 8
mpiexec -n 2 ./build/MyApp
```

## Notes

- `find_package(AMReX REQUIRED 3D EB ...)` fails if the installed AMReX
  lacks a requested component. The available targets are
  `AMReX::amrex_1d`, `AMReX::amrex_2d` and `AMReX::amrex_3d` for each
  dimension AMReX was built for.
- The project inherits AMReX's compiler, MPI, OpenMP and GPU settings. For
  a CUDA build of AMReX, the project enables CUDA and calls
  `setup_target_for_cuda_compilation`.
