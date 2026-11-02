# ExternalMPI ★★

Use AMReX inside an application that manages MPI itself. The application
calls `MPI_Init`, splits `MPI_COMM_WORLD` into two groups (even and odd
ranks), and initializes AMReX on each group's communicator. It also handles
its own command line and sets AMReX parameters in a callback. Requires MPI.

**Key APIs:** `amrex::Initialize(argc, argv, build_parm_parse, comm, func)`,
`ParmParse::add`, `ParmParse::addfile`, `ParallelDescriptor::NProcs`

## Notes

- If MPI is already initialized, AMReX uses a duplicate of the given
  communicator and does not finalize MPI; the application calls
  `MPI_Finalize` after `amrex::Finalize`.
- With `build_parm_parse = false`, AMReX ignores `argc`/`argv`. Here the
  application reads the inputs file itself with `ParmParse::addfile`, so
  `name=value` arguments have no effect.
- The callback runs before AMReX reads its own parameters (e.g.
  `amrex.fpe_trap_invalid`), so it is the place to change them in code.
- Each AMReX instance only knows the ranks in its group. On GPUs, ranks of
  both groups may share a device, so the callback disables the large
  up-front memory reservation (`amrex.the_arena_init_size`).

## Build and run

CMake, from the repository root (requires `AMReX_MPI=ON`):

```
cmake -S . -B build -DAMReXTutorials_TOPICS=Base
cmake --build build --target Base_ExternalMPI -j 8
cd build/Base/ExternalMPI && mpiexec -n 4 ./ExternalMPI inputs
```

GNU Make, with AMReX cloned next to this repository:

```
make -j8
mpiexec -n 4 ./main3d.gnu.MPI.ex inputs
```

## Expected output

The two groups print independently, so the lines may interleave.

```
Group 0: AMReX runs on 2 of 4 MPI ranks
Group 0: sum = 4096 (exact 4096)
Group 0: PASSED
Group 1: AMReX runs on 2 of 4 MPI ranks
Group 1: sum = 8192 (exact 8192)
Group 1: PASSED
```
