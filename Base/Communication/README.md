# Communication ★★

The three basic ways to exchange `MultiFab` data between boxes and MPI
ranks:

- `FillBoundary` copies valid data into the ghost cells of neighboring
  boxes, including across periodic boundaries.
- `ParallelCopy` and `ParallelAdd` copy or add data from one `MultiFab` to
  another with a different `BoxArray` and `DistributionMapping`, e.g.,
  after regridding.
- `SumBoundary` is the reverse of `FillBoundary`: it adds values in ghost
  cells to the valid cells they overlap, as needed for deposition.

Each valid cell holds a unique number, so every result is compared exactly
with a reference `MultiFab`; the program checks the full `FillBoundary` and
`SumBoundary` results.

**Key APIs:** `FillBoundary`, `FillBoundary_nowait`, `FillBoundary_finish`,
`ParallelCopy`, `ParallelAdd`, `SumBoundary`, `Geometry::periodicity`,
`MultiFab::setVal`, `MultiFab::Copy`, `MultiFab::Subtract`, `MultiFab::norm0`

## What it does

1. **FillBoundary.** A `MultiFab` with two ghost cells first fills only the
   first ghost layer (`FillBoundary(IntVect(1), period)`), leaving the
   second layer untouched. Then all ghost cells are filled with
   the non-blocking `FillBoundary_nowait` and `FillBoundary_finish`; work
   that does not need the ghost cells can be done in between, overlapping
   it with the communication.
2. **ParallelCopy and ParallelAdd.** The destination `BoxArray` is the
   domain shifted by `n_cell/4` and chopped into boxes of at most
   `dst_max_grid_size` cells, so it partly lies outside the domain. The
   copy fills its valid and ghost cells, using periodicity for the cells
   outside the domain. `ParallelAdd` then adds the source again.
3. **SumBoundary.** Every valid cell deposits 1 into the $3^d$ cells
   around it, itself included, some of them ghost cells. After
   `SumBoundary`, every valid cell holds $3^d$.

## Build and run

CMake, from the repository root:

```
cmake -S . -B build -DAMReXTutorials_TOPICS=Base
cmake --build build --target Base_Communication -j 8
cd build/Base/Communication && mpiexec -n 2 ./Communication inputs
```

GNU Make, with AMReX cloned next to this repository:

```
make -j8
mpiexec -n 2 ./main3d.gnu.MPI.ex inputs
```

Try different numbers of MPI ranks and box sizes, e.g.,
`max_grid_size=8 dst_max_grid_size=20`.

## Expected output

```
FillBoundary, first ghost layer: max error = 0
FillBoundary, both ghost layers: max error = 0
ParallelCopy: max error = 0
ParallelAdd: max error = 0
SumBoundary: max error = 0
PASSED
```
