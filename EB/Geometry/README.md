# Geometry ★

Build embedded boundary (EB) geometries and query the cut-cell data. The
program builds each geometry listed in `geometries` on `[-1,1]^D` (default
`sphere csg parser`; the `inputs` file adds `stl`):

- `sphere`: a sphere (a circle in 2D) with the fluid inside.
- `csg`: constructive solid geometry. A box with a cylindrical hole, rotated
  and translated, plus a hemisphere made by intersecting a ball with a plane.
- `parser`: an implicit function given as a runtime expression,
  `parser.function`.
- `stl`: a triangulated surface from an STL file (`octahedron.stl`), read
  through the `eb2.*` runtime parameters.

For each geometry, the program counts regular, cut and covered cells, sums
the fluid volume, the EB surface area, and the fluid area of the plane
`x = 0` from the x-face area fractions. It compares them with the exact
values and fails if the relative error is larger than 2% (at `n_cell = 32`).
It writes a plotfile `plt_<geometry>` with the cell type and the volume
fraction, and in 3D the EB surface of `surface_geometry` to `eb.pvtp` for
ParaView.

An implicit function is negative in the fluid and positive in the body.
The CSG operations (`makeUnion`, `makeIntersection`, `makeDifference`) act on
the bodies.

**Key APIs:** `EB2::SphereIF`, `EB2::BoxIF`, `EB2::CylinderIF`, `EB2::PlaneIF`,
`EB2::ParserIF`, `EB2::makeUnion`, `EB2::makeIntersection`,
`EB2::makeDifference`, `EB2::rotate`, `EB2::translate`, `EB2::makeShop`,
`EB2::Build`, `makeEBFabFactory`, `EBFArrayBoxFactory::getEBDataArrays`,
`EBData_t`, `EB_WriteSingleLevelPlotfile`, `WriteEBSurface`

## Build and run

CMake, from the repository root:

```
cmake -S . -B build -DAMReXTutorials_TOPICS=EB
cmake --build build --target EB_Geometry -j 8
cd build/EB/Geometry && mpiexec -n 2 ./Geometry inputs
```

GNU Make, with AMReX cloned next to this repository:

```
make -j8
mpiexec -n 2 ./main3d.gnu.MPI.ex inputs
```

Try `geometries=parser parser.function="x^2+y^2+z^2-0.25"` or a
finer grid, `n_cell = 128`.

## Expected output

```
Geometry: sphere
  regular, cut, covered cells: 2920, 1760, 28088
  fluid volume:        0.90020653
  EB surface area:     4.51266773
  fluid area at x = 0: 1.129279637
  expected volume:     0.9047786842, relative error 0.005053339927
  ...
PASSED
```
