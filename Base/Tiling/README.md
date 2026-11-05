# Tiling ★★

The boxes an `MFIter` loop works on, and the bugs that come from picking
the wrong one. On CPU, an `MFIter` constructed with tiling splits each box
of a rank into tiles, and the loop body runs once per tile, in parallel
with OpenMP. `tilebox()` is the current tile, `validbox()` the whole box
the tile belongs to, and `fabbox()` the box with its ghost cells. The
program prints these boxes for one box, then counts how many times each
cell is written for several choices of the box of an iteration, which shows
which choices are right and which are not, and finally times a stencil with
and without tiling.

**Key APIs:** `MFIter(mf, tile_size)`, `TilingIfNotGPU`, `MFIter::tilebox`,
`MFIter::validbox`, `MFIter::fabbox`, `MFIter::growntilebox`,
`MFIter::nodaltilebox`, `MFIter::index`, `MFIter::LocalTileIndex`,
`HostDevice::Atomic::Add`, `fabarray.mfiter_tile_size`

## Which box

- **`tilebox()`** for work on the cells of the tile. The tiles of a box do
  not overlap and together cover the valid box, so every cell is done once.
  This is the box for almost all kernels.
- **`validbox()`** only in a loop *without* tiling, for work that concerns
  the whole box: per-box bookkeeping such as costs for load balancing (see
  [LoadBalance](../LoadBalance)), or one launch per box on GPUs.
- **`growntilebox()`** for the tile plus the ghost cells on the sides of the
  tile that lie on the boundary of the box. The grown tiles do not overlap
  either, and together they cover the valid and ghost cells.
- **`nodaltilebox(d)`** for the faces of the tile in direction `d`, and
  `tilebox(nodal_flag)` for any index type. Each face belongs to one tile,
  including the faces between two tiles (see
  [StaggeredGrids](../StaggeredGrids)).
- **`fabbox()`** is the box of the `FArrayBox`, i.e., the valid box with its
  ghost cells; it is the same for every tile of a box.

## The bugs

- **`validbox()` in a tiled loop.** Every tile does the whole box: the work
  is multiplied by the number of tiles, a kernel that accumulates gives a
  wrong answer, and with OpenMP several threads write the same cells at
  the same time. `FArrayBox` operations such as `fab.setVal` or `fab.plus`
  inside a tiled loop are the same bug, since they act on the whole box by
  default; pass them the tile box.
- **`tilebox()` where the whole box is needed.** Per-box bookkeeping in a
  tiled loop is done once per tile instead of once per box, e.g. a value
  assigned per box is overwritten by each tile. Use a loop without tiling
  and `validbox()`, as [LoadBalance](../LoadBalance) does for its costs.
- **`grow(tilebox(), ng)`.** A tile box grown by hand overlaps the
  neighboring tiles, so the cells near tile boundaries are written by
  several tiles, from several threads. Use `growntilebox()`.
- **`surroundingNodes(tilebox(), d)`.** Converting the tile box by hand
  gives the face between two tiles to both of them. Use `nodaltilebox(d)`.

The program makes these mistakes on purpose with atomic adds, so that the
counts are exact and the program does not crash; without the atomics, the
overlaps are data races.

## Notes

- `MFIter mfi(mf, TilingIfNotGPU())` is the usual form. `TilingIfNotGPU()`
  is true on CPU and false on GPU, where a box is already split among the
  GPU threads and tiling would only make the kernels smaller. The tile
  size is the runtime parameter `fabarray.mfiter_tile_size`, by default
  `1024000 8 8` in 3D on CPU builds (tiles are not split in `x`, to keep
  the inner loop long for vectorization) and no splitting in 2D or on
  GPU builds, and `MFIter mfi(mf, IntVect(...))` sets it for one loop.
  The demonstrations use the explicit form so that they show tiles in
  2D and on GPU builds too.
- With OpenMP, the `MFIter` loop is inside a `#pragma omp parallel` region
  and the tiles of a rank are shared among the threads, so a tiled loop
  must not have an `omp parallel` or an `omp for` inside. The region is
  skipped on GPU with `if (Gpu::notInLaunchRegion())`.
- `mfi.index()` is the global index of the box and `mfi.LocalIndex()` the
  index among the boxes of this rank; `mfi.LocalTileIndex()` is the number
  of the tile within its box and `mfi.tileIndex()` numbers the iterations.
  A `LayoutData` is indexed with `mfi` or `mfi.index()`.
- Tiling helps a stencil on CPU when a box does not fit in cache; the
  benefit depends on the kernel, the box size and the tile size. The tile
  size is for tuning on the command line, e.g.
  `fabarray.mfiter_tile_size=1024000 16 16`.

## Build and run

CMake, from the repository root:

```
cmake -S . -B build -DAMReXTutorials_TOPICS=Base
cmake --build build --target Base_Tiling -j 8
cd build/Base/Tiling && mpiexec -n 2 ./Tiling inputs
```

GNU Make, with AMReX cloned next to this repository:

```
make -j8 USE_OMP=TRUE
OMP_NUM_THREADS=4 mpiexec -n 2 ./main3d.gnu.MPI.OMP.ex inputs
```

For timings, build with TinyProfiler on (`-DAMReX_TINY_PROFILE=ON` or
`make TINY_PROFILE=TRUE`) and use a larger problem, e.g.
`n_cell=256 max_grid_size=128 nrepeat=20`; the timers are
`laplacian_tiled` and `laplacian_untiled`. The difference between them
depends on the box size, the number of OpenMP threads and the cache.

## Expected output

```
Without tiling, tilebox ((0,0,0) (15,15,15) (0,0,0)) = validbox ((0,0,0) (15,15,15) (0,0,0))
  (1 iteration)
With tile size (1024000,8,8):
  tile 0 of box 0: tilebox ((0,0,0) (15,7,7) (0,0,0))
      growntilebox()  ((-1,-1,-1) (16,7,7) (0,0,0))
      nodaltilebox(2) ((0,0,0) (15,7,7) (0,0,1))
  tile 1 of box 0: tilebox ((0,8,0) (15,15,7) (0,0,0))
      growntilebox()  ((-1,8,-1) (16,16,7) (0,0,0))
      nodaltilebox(2) ((0,8,0) (15,15,7) (0,0,1))
  tile 2 of box 0: tilebox ((0,0,8) (15,7,15) (0,0,0))
      growntilebox()  ((-1,-1,8) (16,7,16) (0,0,0))
      nodaltilebox(2) ((0,0,8) (15,7,16) (0,0,1))
  tile 3 of box 0: tilebox ((0,8,8) (15,15,15) (0,0,0))
      growntilebox()  ((-1,8,8) (16,16,16) (0,0,0))
      nodaltilebox(2) ((0,8,8) (15,15,16) (0,0,1))
  validbox ((0,0,0) (15,15,15) (0,0,0)) and fabbox ((-1,-1,-1) (16,16,16) (0,0,0)) are the same for every tile

Domain ((0,0,0) (63,63,63) (0,0,0)) in 8 boxes, tile size (1024000,8,8)
  tilebox          : visits per valid cell 1 to 1, ghost cells included 0 to 1
  validbox         : visits per valid cell 16 to 16, ghost cells included 0 to 16
  grow(tilebox, 1) : visits per valid cell 1 to 4, ghost cells included 1 to 4
  growntilebox()   : visits per valid cell 1 to 1, ghost cells included 1 to 1
  nodaltilebox(2)  : visits per face 1 to 1
  surroundingNodes(tilebox, 2): visits per face 1 to 2
  Laplacian with and without tiling: max difference 0
PASSED
```
