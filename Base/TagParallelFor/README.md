# TagParallelFor ★★★

Run many small kernels as one launch. Filling the ghost cells outside the
domain boundary touches a thin box per domain face of every box at the
boundary: many small boxes of different sizes, each too small to fill a GPU.
The program fills them with a reflecting boundary condition in three ways,
checks that the results agree, and times them with `BL_PROFILE` (see
[HeatEquation](../HeatEquation#profiling)):

1. one `ParallelFor` per box and face, inside an `MFIter` loop: the
   baseline, not recommended on GPU for many small boxes since no launch
   fills the device;
2. one `ParallelFor` over a `Vector` of tags, which packs all the boxes into
   a single launch;
3. the same with a `TagVector` built once and reused.

**Key APIs:** `ParallelFor(Vector<Tag>, ncomp, f)`, `ParallelFor(TagVector<Tag>, ncomp, f)`,
`TagVector`, `Array4BoxOrientationTag`, `Orientation`, `OrientationIter`,
`adjCell`, `BL_PROFILE`

## Notes

- A tag is a small struct with a `box()` member (or `size()` for 1D work);
  the kernel receives `(i, j, k, n, tag)`, so it finds the `Array4`, the box
  and any other data it needs in the tag. AMReX provides `Array4BoxTag`,
  `Array4BoxValTag`, `Array4BoxOrientationTag`, `Array4PairTag`,
  `Array4CopyTag` and `VectorTag`; any struct of your own with `box()` or
  `size()` works too.
- On GPU, the launch is one kernel over all the tags; a warp never spans two
  tags, and the work is balanced at warp granularity. On CPU, it is a loop
  over the tags that opens its own OpenMP parallel region, so do not call
  it inside one.
- `ParallelFor(Vector<Tag>, ...)` builds a `TagVector` (copies the tags to
  the device) and synchronizes the stream on every call. When the same tags
  are used repeatedly, build the `TagVector` once and rebuild it only when
  the grids or the `MultiFab` change.
- Use `ParallelFor(mf, ...)` (see [MultiFab](../MultiFab)) when the work is
  the same on every box. Tags are for irregular sets of boxes: physical
  boundary conditions, corrections at coarse/fine faces, operations on a
  subset of boxes. AMReX itself uses them in `FillBoundary`, `FluxRegister`
  and the linear solvers.

## Build and run

CMake, from the repository root:

```
cmake -S . -B build -DAMReXTutorials_TOPICS=Base
cmake --build build --target Base_TagParallelFor -j 8
cd build/Base/TagParallelFor && mpiexec -n 2 ./TagParallelFor inputs
```

GNU Make, with AMReX cloned next to this repository:

```
make -j8
mpiexec -n 2 ./main3d.gnu.MPI.ex inputs
```

For timings, build with TinyProfiler on (`-DAMReX_TINY_PROFILE=ON` or
`make TINY_PROFILE=TRUE`) and an optimized GPU build, and use e.g.
`n_cell=256 nrepeat=20 tiny_profiler.print_threshold=0` so that even the
fastest method is listed. The difference between the methods grows with the
number of boxes, i.e. with `n_cell/max_grid_size`; on a GPU, the fused
launch is typically one to two orders of magnitude faster than thousands of
tiny launches.

## Expected output

```
512 boxes, 384 boundary faces
Max difference from the expected result: per-face launches 0, tags rebuilt 0, tags cached 0
PASSED
```

With TinyProfiler on, the profiling tables list the timers
`per_face_launches`, `tags_rebuilt` and `tags_cached`. Timers below
`tiny_profiler.print_threshold` (1% by default) are folded into `Other`.
