# StaggeredGrids ★★

Data on cell centers, faces and nodes of the same grids. A `MultiFab`'s
`BoxArray` carries the index type: `convert` turns the cell-centered
`BoxArray` into one of face- or node-centered boxes. The program sets a
velocity on the faces (as a MAC projection or a finite-volume flux needs
it), takes its divergence, averages it to cell centers, averages a
cell-centered scalar to the faces, averages nodal data to cell centers, and
writes everything to a plotfile. The discrete divergence of the shear flow
is exactly zero, and the averages reproduce linear functions to round-off.
The tutorial runs in 2D and 3D.

**Key APIs:** `convert`, `IntVect::TheDimensionVector`,
`IntVect::TheNodeVector`, `IndexType`, `MultiFab::ixType`,
`MFIter::nodaltilebox`, `ParallelFor(box1, box2, box3, f1, f2, f3)`,
`OverrideSync`, `computeDivergence`, `average_face_to_cellcenter`,
`average_cellcenter_to_face`, `average_node_to_cellcenter`,
`GetArrOfPtrs`, `GetArrOfConstPtrs`, `MultiFab::Subtract`, `MultiFab::norm0`

## Notes

- Face data for direction `d` lives on
  `convert(ba, IntVect::TheDimensionVector(d))`; it has one more point
  than cells in direction `d`. Nodal data lives on
  `convert(ba, IntVect::TheNodeVector())`. `amrex::convert` and
  `amrex::surroundingNodes` do the same for a single `Box`;
  `amrex::enclosedCells` goes back to cells.
- A node at index `i` sits at `i*dx`, a cell center at `(i+1/2)*dx`, so
  the position of a point depends on the index type (`ixType()`), as in
  the `position` function of the example.
- In an `MFIter` loop over cell-centered data, `mfi.nodaltilebox(d)` is the
  box of the faces of the tile in direction `d`. The tiles' face boxes do not
  overlap: within a box, each face belongs to one tile, so face data can be
  written in the same loop as the cell data, without race conditions with
  OpenMP. `ParallelFor` takes one box and one function per direction and
  runs them in a single launch.
- A face on the boundary between two boxes is stored in both boxes.
  `OverrideSync` copies the value of one box to the others so that the
  copies agree. The example does not need it, because the value depends
  only on the index, but it is needed when the boxes compute different
  values, for instance from round-off in a sum.
- Face data can have ghost faces, filled by `FillBoundary` like ghost
  cells; `FillBoundaryAndSync` does that and `OverrideSync` in one call.
  The divergence and averaging functions read valid faces only, so the
  example has no ghost faces.
- The face-based functions take an `Array<MultiFab*,AMREX_SPACEDIM>` or
  the `const` version, made with `GetArrOfPtrs`/`GetArrOfConstPtrs`.
- Plotfile viewers usually expect cell-centered data; write face and nodal
  data after averaging them to cell centers.

## Build and run

CMake, from the repository root:

```
cmake -S . -B build -DAMReXTutorials_TOPICS=Base
cmake --build build --target Base_StaggeredGrids -j 8
cd build/Base/StaggeredGrids && mpiexec -n 2 ./StaggeredGrids inputs
```

GNU Make, with AMReX cloned next to this repository:

```
make -j8
mpiexec -n 2 ./main3d.gnu.MPI.ex inputs
```

## Expected output (3D)

```
Index types: phi (C,C,C), umac[0] (N,C,C), psi (N,N,N)
Max |div u| = 0, max error of the cell-averaged velocity = 0
Max error of the face-averaged scalar = 0
Max error of the cell-averaged nodal data = 0
PASSED
```

The plotfile `plt_staggered` holds the divergence, the cell-averaged
velocity and the cell-averaged nodal data.
