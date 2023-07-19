import numpy as np
import amrex.space3d as amr

def main_main ():
    """
    TODO Add documentation
    """
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

    # time = starting time in the simulation
    time = 0.

    # Loop over boxes
    for mfi in phi_old:
        bx = mfi.validbox()
        # phiOld is indexed in reversed order (z,y,x) and indices are local
        phiOld = np.array(phi_old.array(mfi), copy=False)
        # set phi = 1 + e^(-(r-0.5)^2)
        x = (np.arange(bx.small_end[0],bx.big_end[0],1) + 0.5) * dx[0]
        y = (np.arange(bx.small_end[1],bx.big_end[1],1) + 0.5) * dx[1]
        z = (np.arange(bx.small_end[2],bx.big_end[2],1) + 0.5) * dx[2]
        rsquared = ((z[:,np.newaxis,np.newaxis]-0.5)**2
                  + (y[np.newaxis,:,np.newaxis]-0.5)**2
                  + (x[np.newaxis,np.newaxis,:]-0.5)**2) / 0.01
        phiOld = 1. + np.exp(-rsquared)

    # Write a plotfile of the initial data if plot_int > 0
    #if (plot_int > 0)
    #{
    #    int step = 0;
    #    const std::string& pltfile = amrex::Concatenate("plt",step,5);
    #    WriteSingleLevelPlotfile(pltfile, phi_old, {"phi"}, geom, time, 0);
    #}

    for step in range(nsteps):
        # Fill periodic ghost cells
        phi_old.fill_boundary(geom.periodicity())

    # new_phi = old_phi + dt * Laplacian(old_phi)
    # Loop over boxes
    for mfi in phi_old:
        #bx = mfi.validbox()
        phiOld = np.array(phi_old.array(mfi), copy=False)
        phiNew = np.array(phi_new.array(mfi), copy=False)
        # Advance the data by dt
        phiNew[ngz:-ngz,ngy:-ngy,ngx:-ngx] = (
            phiOld[ngz:-ngz,ngy:-ngy,ngx:-ngx]
            + dt*((   phiOld[ngz  :-ngz  , ngy  :-ngy  , ngx+1:-ngx+1]
                   -2*phiOld[ngz  :-ngz  , ngy  :-ngy  , ngx  :-ngx  ]
                     +phiOld[ngz  :-ngz  , ngy  :-ngy  , ngx-1:-ngx-1])/dx[0]**2
                 +(   phiOld[ngz  :-ngz  , ngy+1:-ngy+1, ngx  :-ngx  ]
                   -2*phiOld[ngz  :-ngz  , ngy  :-ngy  , ngx  :-ngx  ]
                     +phiOld[ngz  :-ngz  , ngy-1:-ngy-1, ngx  :-ngx  ])/dx[1]**2
                 +(   phiOld[ngz+1:-ngz+1, ngy  :-ngy  , ngx  :-ngx  ]
                   -2*phiOld[ngz  :-ngz  , ngy  :-ngy  , ngx  :-ngx  ]
                     +phiOld[ngz-1:-ngz-1, ngy  :-ngy  , ngx  :-ngx  ])/dx[2]**2))
        
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

main_main()

# Finalize AMReX
amr.finalize()

##include <AMReX_PlotFileUtil.H>
##include <AMReX_ParmParse.H>
#
##include "myfunc.H"
#
#void main_main ()
#{
#
#    // inputs parameters
#    {
#        // ParmParse is way of reading inputs from the inputs file
#        // pp.get means we require the inputs file to have it
#        // pp.query means we optionally need the inputs file to have it - but we must supply a default here
#        ParmParse pp;
#
#        // We need to get n_cell from the inputs file - this is the number of cells on each side of
#        //   a square (or cubic) domain.
#        pp.get("n_cell",n_cell);
#
#        // The domain is broken into boxes of size max_grid_size
#        pp.get("max_grid_size",max_grid_size);
#
#        // Default nsteps to 10, allow us to set it to something else in the inputs file
#        nsteps = 10;
#        pp.query("nsteps",nsteps);
#
#        // Default plot_int to -1, allow us to set it to something else in the inputs file
#        //  If plot_int < 0 then no plot files will be written
#        plot_int = -1;
#        pp.query("plot_int",plot_int);
#
#        // time step
#        pp.get("dt",dt);
#    }
#
#    // **********************************
#    // SIMULATION SETUP
#
#
#    // **********************************
#    // INITIALIZE DATA
#
#    for (int step = 1; step <= nsteps; ++step)
#    {
#        // new_phi = old_phi + dt * Laplacian(old_phi)
#        // loop over boxes
#        for ( MFIter mfi(phi_old); mfi.isValid(); ++mfi )
#        {
#            const Box& bx = mfi.validbox();
#
#            const Array4<Real>& phiOld = phi_old.array(mfi);
#            const Array4<Real>& phiNew = phi_new.array(mfi);
#
#            // advance the data by dt
#            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
#            {
#                phiNew(i,j,k) = phiOld(i,j,k) + dt *
#                    ( (phiOld(i+1,j,k) - 2.*phiOld(i,j,k) + phiOld(i-1,j,k)) / (dx[0]*dx[0])
#                     +(phiOld(i,j+1,k) - 2.*phiOld(i,j,k) + phiOld(i,j-1,k)) / (dx[1]*dx[1])
##if (AMREX_SPACEDIM == 3)
#                     +(phiOld(i,j,k+1) - 2.*phiOld(i,j,k) + phiOld(i,j,k-1)) / (dx[2]*dx[2])
##endif
#                        );
#            });
#        }
#
#        // update time
#        time = time + dt;
#
#        // copy new solution into old solution
#        MultiFab::Copy(phi_old, phi_new, 0, 0, 1, 0);
#
#        // Tell the I/O Processor to write out which step we're doing
#        amrex::Print() << "Advanced step " << step << "\n";
#
#        // Write a plotfile of the current data (plot_int was defined in the inputs file)
#        if (plot_int > 0 && step%plot_int == 0)
#        {
#            const std::string& pltfile = amrex::Concatenate("plt",step,5);
#            WriteSingleLevelPlotfile(pltfile, phi_new, {"phi"}, geom, time, step);
#        }
#    }
#}
