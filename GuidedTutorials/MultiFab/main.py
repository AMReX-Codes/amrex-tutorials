import numpy as np
import amrex.space3d as amr

# Initialize AMReX
amr.initialize([])

# TODO How do we print pyamrex's version?
if amr.ParallelDescriptor.IOProcessor():
    print("Hello world from pyamrex\n")

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

# TODO Should this be &real_box (reference) instead?
coord = 0  # Cartesian
is_per = [0, 0, 0] # periodicity
geom = amr.Geometry(domain, real_box, coord, is_per)

# Calculate cell sizes
dx = geom.data().CellSize() # dx[0]=dx dx[1]=dy dx[2]=dz

# Fill a MultiFab with data
mf_val = 0.
for mfi in mf:
    bx = mfi.validbox()
    # Preferred way to fill array using NumPy operations:
    # - mf_array is indexed in reversed order (n,z,y,x)
    # - indices are local (range from 0 to box size)
    mf_array = np.array(mf.array(mfi), copy=False)
    x = (np.arange(bx.small_end[0], bx.big_end[0]+1)+0.5)*dx[0]
    y = (np.arange(bx.small_end[1], bx.big_end[1]+1)+0.5)*dx[1]
    z = (np.arange(bx.small_end[2], bx.big_end[2]+1)+0.5)*dx[2]
    v = (x[np.newaxis,np.newaxis,:]
       + y[np.newaxis,:,np.newaxis]*0.1
       + z[:,np.newaxis,np.newaxis]*0.01)
    mf_array = 1. + np.exp(-v)
    # Alternative way to fill array using standard for loop over box indices:
    # - mf_array_ref is indexed in standard order (x,y,z,n)
    # - indices are global
    mf_array_ref = mf.array(mfi)
    for i, j, k in bx:
        x = (i+0.5)*dx[0]
        y = (j+0.5)*dx[1]
        z = (k+0.5)*dx[2]
        v = x + y*0.1 + z*0.01
        mf_array_ref[i,j,k,0] = 1. + np.exp(-v)
    # Check that mf_array has been filled correctly, compare it with
    # mf_array_ref (filled as Array4 and transformed to NumPy array)
    mf_array_ref_np = np.array(mf_array_ref, copy=False)
    assert np.allclose(mf_array, mf_array_ref_np, rtol=1e-16, atol=1e-16)

# Plot MultiFab data
plotfile = amr.concatenate(root="plt", num=1, mindigits=3)
varnames = amr.Vector_string(["comp0"])
amr.write_single_level_plotfile(plotfile, mf, varnames, geom, time=0., level_step=0)

# Finalize AMReX
amr.finalize()
