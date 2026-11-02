#!/usr/bin/env python3
"""Explicit heat equation d(phi)/dt = D Laplacian(phi) on a periodic unit
cube with pyAMReX, checked against the exact solution. This is the Python
version of Base/HeatEquation.

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
    # NumPy on CPU, CuPy on GPU, matching the views returned by to_xp
    if amr.Config.have_gpu:
        import cupy as xp
    else:
        import numpy as xp

    pp = amr.ParmParse("")
    n_cell = query(pp, "n_cell", 32)
    max_grid_size = query(pp, "max_grid_size", 16)
    diffusion_coef = query(pp, "diffusion_coef", 1.0)
    cfl = query(pp, "cfl", 0.9)
    stop_time = query(pp, "stop_time", 0.01)
    plot_int = query(pp, "plot_int", -1)

    # Periodic unit domain
    domain = amr.Box(amr.IntVect(0), amr.IntVect(n_cell - 1))
    ba = amr.BoxArray(domain)
    ba.max_size(max_grid_size)
    dm = amr.DistributionMapping(ba)
    real_box = amr.RealBox([0.0, 0.0, 0.0], [1.0, 1.0, 1.0])
    geom = amr.Geometry(domain, real_box, 0, [1, 1, 1])  # 0: Cartesian
    dx = geom.data().CellSize()

    def sine_mode(bx):
        """prod_d sin(2 pi x_d) at the cell centers of bx, indexed [x, y, z]"""
        s = []
        for d in range(3):
            x = (xp.arange(bx.small_end[d], bx.big_end[d] + 1) + 0.5) * dx[d]
            s.append(xp.sin(2.0 * math.pi * x))
        return s[0][:, None, None] * s[1][None, :, None] * s[2][None, None, :]

    # One ghost cell for the 7-point stencil
    phi_old = amr.MultiFab(ba, dm, 1, 1)
    phi_new = amr.MultiFab(ba, dm, 1, 1)

    # Initial condition: phi = 1 + prod_d sin(2 pi x_d). Each view below
    # includes the ghost cells, so [1:-1] selects the valid cells.
    v = (slice(1, -1),) * 3
    for mfi in phi_new:
        phi = phi_new.array(mfi).to_xp(copy=False)[..., 0]
        phi[v] = 1.0 + sine_mode(mfi.validbox())

    # Explicit time step limit: dt <= dx^2 / (2 * dim * D)
    dt = cfl * dx[0] ** 2 / (2.0 * 3 * diffusion_coef)

    time = 0.0
    step = 0
    varnames = amr.Vector_string(["phi"])
    if plot_int > 0:
        amr.write_single_level_plotfile(
            amr.concatenate("plt", step, 5), phi_new, varnames, geom, time, step
        )

    while time < stop_time * (1.0 - 1.0e-6):
        phi_old, phi_new = phi_new, phi_old
        dt_step = min(dt, stop_time - time)

        # Fill ghost cells, here all from periodic images and neighbors.
        phi_old.fill_boundary(geom.periodicity())

        # phi_new = phi_old + dt * D * Laplacian(phi_old), with array slicing
        fac = dt_step * diffusion_coef
        for mfi in phi_new:
            po = phi_old.array(mfi).to_xp(copy=False)[..., 0]
            pn = phi_new.array(mfi).to_xp(copy=False)[..., 0]
            c = po[v]
            pn[v] = c + fac * (
                (po[2:, 1:-1, 1:-1] - 2.0 * c + po[:-2, 1:-1, 1:-1]) / dx[0] ** 2
                + (po[1:-1, 2:, 1:-1] - 2.0 * c + po[1:-1, :-2, 1:-1]) / dx[1] ** 2
                + (po[1:-1, 1:-1, 2:] - 2.0 * c + po[1:-1, 1:-1, :-2]) / dx[2] ** 2
            )

        time += dt_step
        step += 1
        if step % 10 == 0:
            amr.Print(f"Step {step}, time = {time}")
        if plot_int > 0 and step % plot_int == 0:
            amr.write_single_level_plotfile(
                amr.concatenate("plt", step, 5), phi_new, varnames, geom, time, step
            )

    # Exact solution: the sine mode decays as exp(-D |k|^2 t), |k|^2 = 3 (2 pi)^2.
    decay = math.exp(-diffusion_coef * 3 * 4.0 * math.pi**2 * time)
    error = amr.MultiFab(ba, dm, 1, 0)
    for mfi in error:
        phi = phi_new.array(mfi).to_xp(copy=False)[..., 0]
        err = error.array(mfi).to_xp(copy=False)[..., 0]
        err[...] = xp.abs(phi[v] - (1.0 + decay * sine_mode(mfi.validbox())))
    rel_err = error.max(0) / decay  # max over all ranks, relative to the amplitude
    amr.Print(
        f"Finished {step} steps at time {time}, mode amplitude = {decay}, "
        f"relative max error = {rel_err}"
    )

    if rel_err > 0.02:
        amr.Print("FAILED")
        return 1
    amr.Print("PASSED")
    return 0


if __name__ == "__main__":
    # Command line arguments, e.g., an inputs file, are passed to AMReX.
    amr.initialize(sys.argv[1:])
    status = main()  # AMReX objects must be gone before finalize
    amr.finalize()
    sys.exit(status)
