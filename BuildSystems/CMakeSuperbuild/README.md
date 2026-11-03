# CMakeSuperbuild

A standalone CMake project that downloads and builds AMReX together with
the application, using `FetchContent`. Most AMReX applications are built
this way. Copy this directory to start a new project.

## Build and run

```
cmake -S . -B build
cmake --build build -j 8
mpiexec -n 2 ./build/MyApp
```

## Options

- AMReX options are ordinary CMake cache variables, e.g.,
  `-DAMReX_SPACEDIM=2 -DAMReX_OMP=ON -DAMReX_EB=ON`.
- GPU: `-DAMReX_GPU_BACKEND=CUDA -DCMAKE_CUDA_ARCHITECTURES=native`, or `HIP`
  with `-DAMReX_AMD_ARCH=gfx90a`, or `SYCL` with Intel's `icpx` compiler.
- Use a local AMReX source tree instead of downloading:
  `-DFETCHCONTENT_SOURCE_DIR_AMREX=/path/to/amrex`.
- Pin a version by changing `GIT_TAG` to a release such as `26.10`.

Other approaches: [CMakeInstalledAMReX](../CMakeInstalledAMReX) links to an
installed AMReX; [GNUMake](../GNUMake) uses AMReX's GNU Make system.
