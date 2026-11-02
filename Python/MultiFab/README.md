# MultiFab (Python) ★

The Python version of [Base/MultiFab](../../Base/MultiFab) with pyAMReX:
build a `MultiFab`, fill it through zero-copy NumPy (CPU) or CuPy (GPU)
views, fill ghost cells, compute reductions, and write a plotfile.

**Key APIs:** `amr.initialize`, `amr.Box`, `amr.BoxArray`,
`amr.DistributionMapping`, `amr.Geometry`, `amr.MultiFab`,
`mf.array(mfi).to_xp()`, `mf.to_xp()`, `mf.fill_boundary`, `mf.min`,
`mf.max`, `mf.sum`, `amr.write_single_level_plotfile`, `amr.ParmParse`

## Zero-copy views

`for mfi in mf:` loops over the boxes owned by this MPI rank, and
`mf.array(mfi).to_xp(copy=False)` returns a NumPy or CuPy array that shares
memory with the box's data. It includes the ghost cells, is indexed
`[x, y, z, component]`, and its indices start at 0, so the valid cells are
`[ng:-ng, ng:-ng, ng:-ng]`. `mf.to_xp()` returns the views of all local
boxes as a list. Writing to a view changes the MultiFab. Pass `order="C"`
for `[component, z, y, x]` indexing, which is usually faster with CuPy.

## Requirements

pyAMReX (`import amrex`), and CuPy for GPU builds of pyAMReX. See the
[pyAMReX documentation](https://pyamrex.readthedocs.io) for installation.

## Run

```
python3 main.py inputs
mpiexec -n 2 python3 main.py inputs n_cell=64
```

With CMake, `-DAMReXTutorials_PYTHON=ON` registers this script as the CTest
test `Python_MultiFab`; set `-DPython3_EXECUTABLE` to a Python that has
pyAMReX installed.

## Expected output

```
min including ghost cells before FillBoundary: 0.0
min including ghost cells after  FillBoundary: 1.0
max(phi)      = 1.9293757195853387
integral(phi) = 1.0055683279968117 (exact: 1.0055683279968317)
max(2 phi)    = 3.8587514391706774
Wrote plt00000
```
