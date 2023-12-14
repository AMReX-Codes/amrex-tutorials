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


# Initialize AMReX
amr.initialize([])

# CPU/GPU logic
xp = load_cupy()

amr.Print(f"Hello world from pyAMReX version {amr.__version__}\n")

# Goals:
# * Define a MultiFab
# * Fill a MultiFab with data
# * Plot it

# Parameters

# Number of data components at each grid point in the MultiFab
ncomp = 1
# How many grid cells in each direction over the problem domain
n_cell = 32
# How many grid cells are allowed in each direction over each box
max_grid_size = 16

# BoxArray: abstract domain setup

# Integer vector indicating the lower coordinate bounds
dom_lo = amr.IntVect(0,0,0)
# Integer vector indicating the upper coordinate bounds
dom_hi = amr.IntVect(n_cell-1, n_cell-1, n_cell-1)
# Box containing the coordinates of this domain
domain = amr.Box(dom_lo, dom_hi)

# Will contain a list of boxes describing the problem domain
ba = amr.BoxArray(domain)

# Chop the single grid into many small boxes
ba.max_size(max_grid_size)

# Distribution Mapping
dm = amr.DistributionMapping(ba)

# Define MultiFab
mf = amr.MultiFab(ba, dm, ncomp, 0)
mf.set_val(0.)

# Geometry: physical properties for data on our domain
real_box = amr.RealBox([0., 0., 0.], [1. , 1., 1.])

coord = 0  # Cartesian
is_per = [0, 0, 0] # periodicity
geom = amr.Geometry(domain, real_box, coord, is_per)

# Calculate cell sizes
dx = geom.data().CellSize() # dx[0]=dx dx[1]=dy dx[2]=dz

# Fill a MultiFab with data
for mfi in mf:
    bx = mfi.validbox()
    # Preferred way to fill array using fast ranged operations:
    # - xp.array is indexed in reversed order (n,z,y,x),
    #   .T creates a view into the AMReX (x,y,z,n) order
    # - indices are local (range from 0 to box size)
    mf_array = xp.array(mf.array(mfi), copy=False).T
    x = (xp.arange(bx.small_end[0], bx.big_end[0]+1)+0.5)*dx[0]
    y = (xp.arange(bx.small_end[1], bx.big_end[1]+1)+0.5)*dx[1]
    z = (xp.arange(bx.small_end[2], bx.big_end[2]+1)+0.5)*dx[2]
    v = (x[xp.newaxis,xp.newaxis,:]
       + y[xp.newaxis,:,xp.newaxis]*0.1
       + z[:,xp.newaxis,xp.newaxis]*0.01)
    rsquared = ((z[xp.newaxis, xp.newaxis,          :] - 0.5)**2
              + (y[xp.newaxis,          :, xp.newaxis] - 0.5)**2
              + (x[         :, xp.newaxis, xp.newaxis] - 0.5)**2) / 0.01
    mf_array[:, :, :, 0] = 1. + xp.exp(-rsquared)

# Plot MultiFab data
plotfile = amr.concatenate(root="plt", num=1, mindigits=3)
varnames = amr.Vector_string(["comp0"])
amr.write_single_level_plotfile(plotfile, mf, varnames, geom, time=0., level_step=0)

# Finalize AMReX
amr.finalize()
