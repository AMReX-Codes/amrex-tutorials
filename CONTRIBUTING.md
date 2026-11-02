# Contributing a tutorial

## Where it goes

Tutorials are organized by topic. Each tutorial lives at exactly
`<Topic>/<Name>/`, e.g., `LinearSolvers/CellPoisson/`.

| Topic             | Contents                                                        |
|-------------------|-----------------------------------------------------------------|
| `Base`            | Core data structures, kernels, reductions, I/O, ParmParse       |
| `AMR`             | Adaptive mesh refinement with `AmrCore` and `AmrLevel`          |
| `LinearSolvers`   | MLMG, GMRES and related solvers                                 |
| `EB`              | Embedded boundaries                                             |
| `Particles`       | Particle containers and particle-mesh operations                |
| `FFT`             | `amrex::FFT` transforms and FFT-based solvers                   |
| `TimeIntegration` | `amrex::TimeIntegrator`, SUNDIALS                               |
| `Coupling`        | MPMD                                                            |
| `Fortran`         | Fortran applications using the AMReX Fortran interfaces         |
| `Python`          | pyAMReX                                                         |
| `ML`              | Machine-learning inference                                      |
| `Applications`    | Complete solvers combining several topics; starting points for new codes |
| `BuildSystems`    | Standalone projects showing how to build against AMReX          |

Applications may use `Source/` and `Exec/<problem>/` below `<Topic>/<Name>/`.

Each tutorial gets a level: ★ (intro), ★★ (intermediate), ★★★ (advanced), or
Blueprint (a complete code meant to be copied). Add a row for it to the
catalog in the top-level `README.md`.

## What a tutorial directory contains

```
README.md       # what it shows, key APIs, how to build and run, expected output
CMakeLists.txt  # one amrex_tutorial_add(...) call
GNUmakefile     # copy from Base/HelloWorld; adjust the Src/*/Make.package list
Make.package    # CEXE_sources / CEXE_headers
inputs          # small enough to run in CI in a few seconds
*.cpp, *.H
```

The topic's `CMakeLists.txt` lists its tutorials with `add_subdirectory`.
`amrex_tutorial_add` is documented in `cmake/AMReXTutorials.cmake`.

## Code

- C++20, no Fortran outside the `Fortran` topic. Use `amrex::Real` and
  precision-agnostic literals (`0.5_rt` or `Real(0.5)`) so that the code
  works in single and double precision.
- Write kernels with `amrex::ParallelFor`, either over a `MultiFab`
  (`ParallelFor(mf, ...)` with `mf.arrays()`) or over boxes from an `MFIter`
  loop with `TilingIfNotGPU()`. The same code must run on CPU, GPU, with MPI
  and with OpenMP.
- Do not use `Gpu::LaunchRegion`, `AMREX_LAUNCH_*`, `AMREX_FOR_*`,
  `AMREX_PARALLEL_FOR_*`, `AMREX_HOST_DEVICE_FOR_*`, old `BL_*` macros
  (`BL_PROFILE` is fine), `FabArray::copy` or Fortran files. CI rejects them.
- Write dimension-agnostic code (`AMREX_D_DECL`, `AMREX_SPACEDIM`) when it
  comes naturally, and restrict the dimensions with `DIMS` otherwise.
- Read runtime parameters with `ParmParse`, with sensible defaults.
- When there is an exact or reference answer, compute the error, print it,
  and return a non-zero exit code if it exceeds a tolerance. CTest runs
  every tutorial.
- Comments explain AMReX concepts briefly. No history or migration notes.

## Testing locally

```
cmake -S . -B build -DFETCHCONTENT_SOURCE_DIR_AMREX=/path/to/amrex \
      -DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_FLAGS="-Wall -Wextra -Wshadow -Werror"
cmake --build build -j
ctest --test-dir build --output-on-failure
```

Also check `make` in the tutorial directory with `AMREX_HOME` pointing to
AMReX, and, if possible, a GPU build (`-DAMReX_GPU_BACKEND=CUDA`).
