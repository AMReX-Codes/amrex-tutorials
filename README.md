# AMReX Tutorials

Tutorials for [AMReX](https://github.com/AMReX-Codes/amrex), a software
framework for block-structured adaptive mesh refinement. The tutorials are
organized by topic and range from introductory examples to complete codes
meant as starting points for new applications. Every tutorial runs on CPUs
and GPUs, with MPI and OpenMP, and builds with CMake or GNU Make. The
exceptions are the Fortran tutorials, which run on CPUs only, and
PyTorchInference, which supports CUDA but not HIP or SYCL.

Documentation: https://amrex-codes.github.io/amrex/tutorials_html/

## Building

With CMake, from the top of this repository (AMReX is downloaded and built
as part of the tutorials):

```
cmake -S . -B build
cmake --build build -j 8
ctest --test-dir build
```

Useful options:

| Option | Meaning |
|---|---|
| `-DAMReXTutorials_TOPICS="Base;FFT"` | Build only some topics |
| `-DAMReX_SPACEDIM=2` | Build for 1, 2 or 3 dimensions (default 3) |
| `-DAMReX_GPU_BACKEND=CUDA -DCMAKE_CUDA_ARCHITECTURES=native` | GPU build (also `HIP`, `SYCL`) |
| `-DAMReX_OMP=ON`, `-DAMReX_MPI=OFF` | OpenMP on, MPI off |
| `-DFETCHCONTENT_SOURCE_DIR_AMREX=/path/to/amrex` | Use a local AMReX source tree |
| `-DAMReX_ROOT=/path/to/install` | Use an installed AMReX |
| `-DAMReX_TESTING=ON` | Trap floating-point exceptions in the tests |
| `-DAMReX_TINY_PROFILE=ON` | Build with TinyProfiler (see [Base/HeatEquation](Base/HeatEquation)) |
| `-DAMReX_SUNDIALS=ON` | Build the SUNDIALS tutorial (SUNDIALS is downloaded) |
| `-DAMReX_FORTRAN_INTERFACES=ON` | Build the Fortran tutorials (needs a Fortran compiler) |
| `-DAMReXTutorials_PYTHON=ON` | Test the pyAMReX tutorials |
| `-DAMReXTutorials_TORCH=ON` | Build the PyTorch tutorial (needs LibTorch) |

A single tutorial is built with its target, `<Topic>_<Name>`, e.g.,
`cmake --build build --target Base_HeatEquation`.

On a GPU, AMReX preallocates 3/4 of the device memory at startup. If the
test ranks share one GPU, they may run out of memory; set a smaller initial
size through the `AMREX_DEFAULT_INIT` environment variable:

```
AMREX_DEFAULT_INIT="amrex.the_arena_init_size=268435456" ctest --test-dir build
```

With GNU Make, clone AMReX next to this repository (or set `AMREX_HOME`)
and run `make -j8` in a tutorial directory. See each tutorial's README.

The `inputs` files use small settings so that each tutorial runs in a few
seconds. Production runs usually use a much larger `max_grid_size` and
`blocking_factor` for performance, especially on GPUs.

To set up your own project, see [BuildSystems](BuildSystems).

## Catalog

Levels: ★ introductory, ★★ intermediate, ★★★ advanced, **Blueprint** a
complete code to copy as the starting point of an application.

### Base: core data structures, kernels and I/O

| Tutorial | Level | What it shows |
|---|---|---|
| [HelloWorld](Base/HelloWorld) | ★ | Initialize AMReX and print from one or all ranks |
| [ParmParse](Base/ParmParse) | ★ | Runtime parameters: arrays, expressions, enums, tables, sections |
| [MultiFab](Base/MultiFab) | ★ | `BoxArray`, `DistributionMapping`, `Geometry`, `MultiFab`, `ParallelFor`, `MFIter`, plotfiles |
| [HeatEquation](Base/HeatEquation) | ★ | Explicit heat equation on a periodic domain |
| [Communication](Base/Communication) | ★★ | `FillBoundary`, `ParallelCopy`, `ParallelAdd` and `SumBoundary` |
| [BoundaryConditions](Base/BoundaryConditions) | ★★ | Physical boundary conditions with `BCRec` and `PhysBCFunct` |
| [CheckpointRestart](Base/CheckpointRestart) | ★★ | Checkpoint, restart, and reading plotfiles back |
| [Reductions](Base/Reductions) | ★★ | `ReduceOps`, `ParReduce`, `Reducer`, MultiFab norms, MPI reductions |
| [GpuContainers](Base/GpuContainers) | ★★ | Device/pinned vectors, async copies, prefix sums, random numbers, arenas |
| [Parser](Base/Parser) | ★★ | Functions defined at runtime and evaluated in GPU kernels |
| [CompileTimeOptions](Base/CompileTimeOptions) | ★★★ | Specialized kernels from runtime options |
| [ExternalMPI](Base/ExternalMPI) | ★★ | AMReX on a communicator owned by the application |

### AMR: adaptive mesh refinement with `AmrCore` and `AmrLevel`

| Tutorial | Level | What it shows |
|---|---|---|
| [StaticHierarchy](AMR/StaticHierarchy) | ★★ | FillPatch, average down and refluxing on a fixed hierarchy |
| [AdvectionAmrCore](AMR/AdvectionAmrCore) | ★★★ / Blueprint | Complete AMR code: tagging, regridding, subcycling, checkpoints |
| [WaveAmrLevel](AMR/WaveAmrLevel) | ★★★ | Wave equation on `Amr`/`AmrLevel` with `AmrLevel::RK` (RK3, RK4) |

### LinearSolvers: multigrid and Krylov solvers

| Tutorial | Level | What it shows |
|---|---|---|
| [CellPoisson](LinearSolvers/CellPoisson) | ★ | `MLPoisson` and `MLMG` basics |
| [ABecLaplacian](LinearSolvers/ABecLaplacian) | ★★ | Variable coefficients, multi-level composite solve |
| [ImplicitHeat](LinearSolvers/ImplicitHeat) | ★★ | Backward Euler heat equation, reusing a solver |
| [NodalProjection](LinearSolvers/NodalProjection) | ★★ | `MLNodeLaplacian` velocity projection |
| [GMRES](LinearSolvers/GMRES) | ★★ | GMRES with a multigrid preconditioner |
| [AlgMG](LinearSolvers/AlgMG) | ★★★ | Algebraic multigrid on an `SpMatrix`, and as MLMG's multigrid type or bottom solver |
| [CurlCurl](LinearSolvers/CurlCurl) | ★★★ | `MLCurlCurl` for edge-centered fields |
| [CustomLinOp](LinearSolvers/CustomLinOp) | ★★★ | Writing your own multigrid operator |

### EB: embedded boundaries

| Tutorial | Level | What it shows |
|---|---|---|
| [Geometry](EB/Geometry) | ★ | Implicit functions, CSG, parser and STL geometries |
| [Poisson](EB/Poisson) | ★★ | `MLEBABecLap` with Dirichlet conditions on the EB |
| [Advection](EB/Advection) | ★★★ | Cut-cell advection with redistribution |

### Particles

| Tutorial | Level | What it shows |
|---|---|---|
| [Basics](Particles/Basics) | ★ | Pure SoA container, push, redistribute, reductions, I/O |
| [ParticleMesh](Particles/ParticleMesh) | ★★ | Deposition and gather |
| [EBTracers](Particles/EBTracers) | ★★ | Tracer particles around an EB, reflected with the signed distance |
| [Sorting](Particles/Sorting) | ★★ | Binning, sorting, partitioning, filtering |
| [NeighborList](Particles/NeighborList) | ★★★ | Molecular dynamics with neighbor lists |
| [ElectrostaticPIC](Particles/ElectrostaticPIC) | ★★★ / Blueprint | Particle-in-cell code |

### FFT: `amrex::FFT` transforms and solvers

| Tutorial | Level | What it shows |
|---|---|---|
| [Basic](FFT/Basic) | ★ | Real-to-complex transforms, spectral derivative and filter |
| [Poisson](FFT/Poisson) | ★★ | `FFT::Poisson` (periodic, Dirichlet, Neumann) and `FFT::PoissonHybrid` |
| [OpenBC](FFT/OpenBC) | ★★ | Free-space Poisson solver |
| [EnergySpectrum](FFT/EnergySpectrum) | ★★ | Kinetic energy spectrum with batched FFTs |

### TimeIntegration

| Tutorial | Level | What it shows |
|---|---|---|
| [RungeKutta](TimeIntegration/RungeKutta) | ★★ | `TimeIntegrator` with Runge-Kutta methods |
| [SUNDIALS](TimeIntegration/SUNDIALS) | ★★★ | Explicit, implicit, IMEX and multirate methods from SUNDIALS |

### Coupling

| Tutorial | Level | What it shows |
|---|---|---|
| [MPMD](Coupling/MPMD) | ★★★ | Coupling two programs (C++ or Python) in one MPI job |

### Fortran: AMReX Fortran interfaces

| Tutorial | Level | What it shows |
|---|---|---|
| [HeatEquation](Fortran/HeatEquation) | ★ | Heat equation written in Fortran: `amrex_multifab`, `amrex_mfiter`, `fill_boundary`, plotfiles |

### Python (pyAMReX)

| Tutorial | Level | What it shows |
|---|---|---|
| [MultiFab](Python/MultiFab) | ★ | MultiFabs with NumPy/CuPy views |
| [HeatEquation](Python/HeatEquation) | ★ | Heat equation written in Python |

### ML

| Tutorial | Level | What it shows |
|---|---|---|
| [PyTorchInference](ML/PyTorchInference) | ★★★ | Evaluating a PyTorch model on MultiFab data |

### Applications: complete codes

| Tutorial | Level | What it shows |
|---|---|---|
| [CNS](Applications/CNS) | Blueprint | Compressible flow on `Amr`/`AmrLevel` with GPU support |
| [IncFlow](Applications/IncFlow) | Blueprint | Incompressible flow on `AmrCore` with MAC and nodal projections: Taylor-Green, double shear layer, lid-driven cavity |

### BuildSystems: your own project

| Example | What it shows |
|---|---|
| [CMakeSuperbuild](BuildSystems/CMakeSuperbuild) | AMReX built with the project via `FetchContent` |
| [CMakeInstalledAMReX](BuildSystems/CMakeInstalledAMReX) | Linking to an installed AMReX with `find_package` |
| [GNUMake](BuildSystems/GNUMake) | A GNU Make project with `Source/` and `Exec/` |

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md).

## License

NOTICE.  This Software was developed under funding from the
U.S. Department of Energy and the U.S. Government consequently retains
certain rights. As such, the U.S. Government has been granted for
itself and others acting on its behalf a paid-up, nonexclusive,
irrevocable, worldwide license in the Software to reproduce,
distribute copies to the public, prepare derivative works, and perform
publicly and display publicly, and to permit other to do so.

The license for these tutorials can be found at [LICENSE](LICENSE).
