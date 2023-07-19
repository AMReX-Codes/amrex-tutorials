import numpy as np
import amrex.space3d as amr

# TODO How do we print pyamrex's version?
print("Hello world from pyamrex\n")

# Initialize AMReX
amr.initialize([])

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
    # mf_array is indexed in reversed order (z,y,x) and indices are local,
    # that is, indices range from 0 to the box size
    mf_array = np.array(mf.array(mfi), copy=False)
    x = (np.arange(bx.small_end[0], bx.big_end[0])+0.5)*dx[0]
    y = (np.arange(bx.small_end[1], bx.big_end[1])+0.5)*dx[1]
    z = (np.arange(bx.small_end[2], bx.big_end[2])+0.5)*dx[2]
    v = (x[np.newaxis,np.newaxis,:]*10
       + y[np.newaxis,:,np.newaxis]*100
       + z[:,np.newaxis,np.newaxis]*1000)
    mf_array = 1. + np.exp(-v)
    # Check that mf_array has been filled correctly for a given set of indices
    ii = bx.small_end[0]
    jj = bx.small_end[1] 
    kk = bx.small_end[2]
    xx = (ii+0.5)*dx[0]
    yy = (jj+0.5)*dx[1]
    zz = (kk+0.5)*dx[2]
    vv = xx*10 + yy*100 + zz*1000
    mf_refval = 1.0 + np.exp(-vv)
    assert(mf_array[0,0,0] == mf_refval)

# TODO How do we write plotfiles?
# Plot MultiFab data
# WriteSingleLevelPlotfile("plt001", mf, {"comp0"}, geom, 0., 0);

# Finalize AMReX
amr.finalize()
