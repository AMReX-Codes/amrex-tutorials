# MultiFab ★

Build the basic AMReX data structures and fill a `MultiFab` with data: a
`Box` for the domain, a `BoxArray` that chops it into grids, a
`DistributionMapping` that assigns grids to MPI ranks, and a `Geometry` for
the physical, periodic domain. The same function is written into two
components in two ways, ghost cells are filled with `FillBoundary`, and the
result is checked and written to a plotfile. For `FillBoundary` and the
other communication functions, see [Communication](../Communication).

**Key APIs:** `Box`, `BoxArray::maxSize`, `DistributionMapping`, `Geometry`,
`RealBox`, `MultiFab`, `ParallelFor(mf, ...)`, `MFIter`, `TilingIfNotGPU`,
`Array4`, `FillBoundary`, `MultiFab::sum/norm0/norm1/norm2`,
`WriteSingleLevelPlotfile`

## Two ways to loop over a MultiFab

- `ParallelFor(mf, ...)` with `mf.arrays()` runs one kernel over all local
  boxes. It is the simplest, and on GPUs it is fast for many small boxes.
- An `MFIter` loop gives one `Array4` per box (or tile on CPU, with
  `TilingIfNotGPU()` and OpenMP). Use it when work is organized per box,
  e.g., temporary `FArrayBox`es or several kernels on the same box.

## Build and run

CMake, from the repository root:

```
cmake -S . -B build -DAMReXTutorials_TOPICS=Base
cmake --build build --target Base_MultiFab -j 8
cd build/Base/MultiFab && mpiexec -n 2 ./MultiFab inputs
```

GNU Make, with AMReX cloned next to this repository:

```
make -j8
mpiexec -n 2 ./main3d.gnu.MPI.ex inputs
```

The plotfile `plt_multifab` can be viewed with, e.g., ParaView, VisIt or yt.

## Expected output

```
Domain ((0,0,0) (31,31,31) (0,0,0)) split into 8 boxes on 2 MPI rank(s)
integral = 1 (exact: 1)
...
max error including ghost cells = 6.661338148e-16
PASSED
```
