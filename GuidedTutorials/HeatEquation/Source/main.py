#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Copyright 2023 The AMReX Community
#
# This file is part of AMReX.
#
# License: BSD-3-Clause-LBNL
# Authors: Revathi Jambunathan, Edoardo Zoni, Olga Shapoval, David Grote, Axel Huebl

import amrex.space3d as amr


def load_cupy():
    if amr.Config.have_gpu:
        try:
            import cupy as cp
            xp = cp
            amr.Print("Note: found and will use cupy")
        except ImportError:
            amr.Print("Warning: GPU found but cupy not available! Trying managed memory in numpy...")
            import numpy as np
            xp = np
        if amr.Config.gpu_backend == "SYCL":
            amr.Print("Warning: SYCL GPU backend not yet implemented for Python")
            import numpy as np
            xp = np

    else:
        import numpy as np
        xp = np
        amr.Print("Note: found and will use numpy")
    return xp


def main(n_cell, max_grid_size, nsteps, plot_int, dt):
    """
    The main function, automatically called below if called as a script.
    """
    # CPU/GPU logic
    xp = load_cupy()

    # AMREX_D_DECL means "do the first X of these, where X is the dimensionality of the simulation"
    dom_lo = amr.IntVect(*amr.d_decl(       0,        0,        0))
    dom_hi = amr.IntVect(*amr.d_decl(n_cell-1, n_cell-1, n_cell-1))

    # Make a single box that is the entire domain
    domain = amr.Box(dom_lo, dom_hi)

    # Make BoxArray and Geometry:
    # ba contains a list of boxes that cover the domain,
    # geom contains information such as the physical domain size,
    # number of points in the domain, and periodicity

    # Initialize the boxarray "ba" from the single box "domain"
    ba = amr.BoxArray(domain)
    # Break up boxarray "ba" into chunks no larger than "max_grid_size" along a direction
    ba.max_size(max_grid_size)

    # This defines the physical box, [0,1] in each direction.
    real_box = amr.RealBox([*amr.d_decl( 0., 0., 0.)], [*amr.d_decl( 1., 1., 1.)])

    # This defines a Geometry object
    # periodic in all direction
    coord = 0 # Cartesian
    is_per = [*amr.d_decl(1,1,1)] # periodicity
    geom = amr.Geometry(domain, real_box, coord, is_per);

    # Extract dx from the geometry object
    dx = geom.data().CellSize()

    # Nghost = number of ghost cells for each array
    Nghost = 1

    # Ncomp = number of components for each array
    Ncomp = 1

    # How Boxes are distrubuted among MPI processes
    dm = amr.DistributionMapping(ba)

    # Allocate two phi multifabs: one will store the old state, the other the new.
    phi_old = amr.MultiFab(ba, dm, Ncomp, Nghost)
    phi_new = amr.MultiFab(ba, dm, Ncomp, Nghost)
    phi_old.set_val(0.)
    phi_new.set_val(0.)

    # time = starting time in the simulation
    time = 0.

    # Ghost cells
    ng = phi_old.nGrowVect
    ngx = ng[0]
    ngy = ng[1]
    ngz = ng[2]

    # Loop over boxes
    for mfi in phi_old:
        bx = mfi.validbox()
        # phiOld is indexed in reversed order (z,y,x) and indices are local
        phiOld = xp.array(phi_old.array(mfi), copy=False)
        # set phi = 1 + e^(-(r-0.5)^2)
        x = (xp.arange(bx.small_end[0],bx.big_end[0]+1,1) + 0.5) * dx[0]
        y = (xp.arange(bx.small_end[1],bx.big_end[1]+1,1) + 0.5) * dx[1]
        z = (xp.arange(bx.small_end[2],bx.big_end[2]+1,1) + 0.5) * dx[2]
        rsquared = ((z[:         , xp.newaxis, xp.newaxis] - 0.5)**2
                  + (y[xp.newaxis, :         , xp.newaxis] - 0.5)**2
                  + (x[xp.newaxis, xp.newaxis, :         ] - 0.5)**2) / 0.01
        phiOld[:, ngz:-ngz, ngy:-ngy, ngx:-ngx] = 1. + xp.exp(-rsquared)

    # Write a plotfile of the initial data if plot_int > 0
    if plot_int > 0:
        step = 0
        pltfile = amr.concatenate("plt", step, 5)
        varnames = amr.Vector_string(['phi'])
        amr.write_single_level_plotfile(pltfile, phi_old, varnames, geom, time, 0)

    for step in range(1, nsteps+1):
        # Fill periodic ghost cells
        phi_old.fill_boundary(geom.periodicity())

        # new_phi = old_phi + dt * Laplacian(old_phi)
        # Loop over boxes
        for mfi in phi_old:
            phiOld = xp.array(phi_old.array(mfi), copy=False)
            phiNew = xp.array(phi_new.array(mfi), copy=False)
            hix = phiOld.shape[3]
            hiy = phiOld.shape[2]
            hiz = phiOld.shape[1]
            # Advance the data by dt
            phiNew[:, ngz:-ngz,ngy:-ngy,ngx:-ngx] = (
                phiOld[:, ngz:-ngz,ngy:-ngy,ngx:-ngx]
                + dt*((   phiOld[:, ngz  :-ngz     , ngy  :-ngy     , ngx+1:hix-ngx+1]
                       -2*phiOld[:, ngz  :-ngz     , ngy  :-ngy     , ngx  :-ngx     ]
                         +phiOld[:, ngz  :-ngz     , ngy  :-ngy     , ngx-1:hix-ngx-1]) / dx[0]**2
                     +(   phiOld[:, ngz  :-ngz     , ngy+1:hiy-ngy+1, ngx  :-ngx     ]
                       -2*phiOld[:, ngz  :-ngz     , ngy  :-ngy     , ngx  :-ngx     ]
                         +phiOld[:, ngz  :-ngz     , ngy-1:hiy-ngy-1, ngx  :-ngx     ]) / dx[1]**2
                     +(   phiOld[:, ngz+1:hiz-ngz+1, ngy  :-ngy     , ngx  :-ngx     ]
                       -2*phiOld[:, ngz  :-ngz     , ngy  :-ngy     , ngx  :-ngx     ]
                         +phiOld[:, ngz-1:hiz-ngz-1, ngy  :-ngy     , ngx  :-ngx     ]) / dx[2]**2))

        # Update time
        time = time + dt

        # Copy new solution into old solution
        phi_old.copy(dst=phi_old, src=phi_new, srccomp=0, dstcomp=0, numcomp=1, nghost=0)

        # Tell the I/O Processor to write out which step we're doing
        amr.Print(f'Advanced step {step}\n')

        # Write a plotfile of the current data (plot_int was defined in the inputs file)
        if plot_int > 0 and step%plot_int == 0:
            pltfile = amr.concatenate("plt", step, 5)
            varnames = amr.Vector_string(['phi'])
            amr.write_single_level_plotfile(pltfile, phi_new, varnames, geom, time, step)


if __name__ == '__main__':
    # Initialize AMReX
    amr.initialize([])

    # TODO Implement parser
    # Simulation parameters
    # number of cells on each side of the domain
    n_cell = 32
    # size of each box (or grid)
    max_grid_size = 16
    # total steps in simulation
    nsteps = 1000
    # how often to write a plotfile
    plot_int = 100
    # time step
    dt = 1e-5

    main(n_cell, max_grid_size, nsteps, plot_int, dt)

    # Finalize AMReX
    amr.finalize()
