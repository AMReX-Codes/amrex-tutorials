#!/usr/bin/env python3
"""MultiFab basics with pyAMReX: build a MultiFab, fill it through zero-copy
NumPy/CuPy views, fill ghost cells, compute reductions, and write a plotfile.

    python3 main.py inputs
    mpiexec -n 2 python3 main.py inputs n_cell=64
"""

import math
import sys

import amrex.space3d as amr


def query(pp, name, default):
    """The value of name in the inputs, or default. The result has the type of default."""
    params = pp.to_dict()
    return type(default)(params[name]) if name in params else default


def main():
    # NumPy on CPU, CuPy on GPU. Array views of AMReX data use the same
    # module (see to_xp below).
    if amr.Config.have_gpu:
        import cupy as xp
    else:
        import numpy as xp

    pp = amr.ParmParse("")
    n_cell = query(pp, "n_cell", 32)
    max_grid_size = query(pp, "max_grid_size", 16)

    # Index space: the domain, chopped into boxes, distributed over MPI ranks
    domain = amr.Box(amr.IntVect(0), amr.IntVect(n_cell - 1))
    ba = amr.BoxArray(domain)
    ba.max_size(max_grid_size)
    dm = amr.DistributionMapping(ba)

    # Physical space: the unit cube, periodic in all directions
    real_box = amr.RealBox([0.0, 0.0, 0.0], [1.0, 1.0, 1.0])
    geom = amr.Geometry(domain, real_box, 0, [1, 1, 1])  # 0: Cartesian
    dx = geom.data().CellSize()

    # Two components and one ghost cell
    ncomp = 2
    ng = 1
    mf = amr.MultiFab(ba, dm, ncomp, ng)
    mf.set_val(0.0)

    # Loop over the boxes owned by this rank. mf.array(mfi).to_xp() is a
    # zero-copy view of the box's data including ghost cells, indexed
    # [x, y, z, component] with local indices starting at 0.
    for mfi in mf:
        bx = mfi.validbox()
        phi = mf.array(mfi).to_xp(copy=False)
        x = (xp.arange(bx.small_end[0], bx.big_end[0] + 1) + 0.5) * dx[0]
        y = (xp.arange(bx.small_end[1], bx.big_end[1] + 1) + 0.5) * dx[1]
        z = (xp.arange(bx.small_end[2], bx.big_end[2] + 1) + 0.5) * dx[2]
        r2 = (
            (x[:, None, None] - 0.5) ** 2
            + (y[None, :, None] - 0.5) ** 2
            + (z[None, None, :] - 0.5) ** 2
        )
        # Valid cells only: skip ng ghost cells on each side
        phi[ng:-ng, ng:-ng, ng:-ng, 0] = 1.0 + xp.exp(-r2 / 0.01)

    # All local boxes at once: a list of views, one per box
    for phi in mf.to_xp():
        phi[..., 1] = 2.0 * phi[..., 0]

    # Ghost cells are still zero. FillBoundary copies them from neighboring
    # boxes and periodic images.
    amr.Print(f"min including ghost cells before FillBoundary: {mf.min(0, ng)}")
    mf.fill_boundary(geom.periodicity())
    amr.Print(f"min including ghost cells after  FillBoundary: {mf.min(0, ng)}")

    # Reductions are done in C++ (on the GPU, if any) and over all ranks
    dv = dx[0] * dx[1] * dx[2]
    integral = mf.sum(0) * dv
    exact = 1.0 + (0.01 * math.pi) ** 1.5
    amr.Print(f"max(phi)      = {mf.max(0)}")
    amr.Print(f"integral(phi) = {integral} (exact: {exact})")
    amr.Print(f"max(2 phi)    = {mf.max(1)}")

    plotfile = amr.concatenate("plt", 0, 5)
    varnames = amr.Vector_string(["phi", "two_phi"])
    amr.write_single_level_plotfile(plotfile, mf, varnames, geom, 0.0, 0)
    amr.Print(f"Wrote {plotfile}")

    ok = abs(integral - exact) < 1.0e-3 * exact
    ok = ok and abs(mf.max(1) - 2.0 * mf.max(0)) < 1.0e-12
    return 0 if ok else 1


if __name__ == "__main__":
    # Command line arguments, e.g., an inputs file, are passed to AMReX.
    amr.initialize(sys.argv[1:])
    status = main()  # AMReX objects must be gone before finalize
    if status != 0:
        amr.Print("FAILED")
    amr.finalize()
    sys.exit(status)
