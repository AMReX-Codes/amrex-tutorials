#include <AMReX.H>

#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>

/*
#ifdef AMREX_USE_CUDA
#include <cufft.h>
#else
#include <fftw3.h>
#include <fftw3-mpi.h>
#endif
*/

int main (int argc, char* argv[])
{
    
    amrex::Initialize(argc,argv);

    // store the current time so we can later compute total run time.
    amrex::Real start_time = amrex::ParallelDescriptor::second();

    // **********************************
    // DECLARE SIMULATION PARAMETERS
    // **********************************

    // number of cells on each side of the domain
    int n_cell;

    // size of each box (or grid)
    int max_grid_size;

    // **********************************
    // READ PARAMETER VALUES FROM INPUT DATA
    // **********************************
    // inputs parameters
    {
        // ParmParse is way of reading inputs from the inputs file
        // pp.get means we require the inputs file to have it
        // pp.query means we optionally need the inputs file to have it - but we must supply a default here
        amrex::ParmParse pp;

        // We need to get n_cell from the inputs file - this is the number of cells on each side of
        //   a square (or cubic) domain.
        pp.get("n_cell",n_cell);

        // The domain is broken into boxes of size max_grid_size
        pp.get("max_grid_size",max_grid_size);
    }

    // **********************************
    // DEFINE SIMULATION SETUP AND GEOMETRY
    // **********************************

    // make BoxArray and Geometry
    // ba will contain a list of boxes that cover the domain
    // geom contains information such as the physical domain size,
    // number of points in the domain, and periodicity
    amrex::BoxArray ba;
    amrex::Geometry geom;

    // define lower and upper indices
    amrex::IntVect dom_lo(AMREX_D_DECL(       0,        0,        0));
    amrex::IntVect dom_hi(AMREX_D_DECL(n_cell-1, n_cell-1, n_cell-1));

    // Make a single box that is the entire domain
    amrex::Box domain(dom_lo, dom_hi);

    // Initialize the boxarray "ba" from the single box "domain"
    ba.define(domain);

    // Break up boxarray "ba" into chunks no larger than "max_grid_size" along a direction
    ba.maxSize(max_grid_size);

    // This defines the physical box, [0,1] in each direction.
    amrex::RealBox real_box({ AMREX_D_DECL(0., 0., 0.)},
                            { AMREX_D_DECL(1., 1., 1.)} );

    // periodic in all direction
    amrex::Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(1,1,1)};

    // This defines a Geometry object
    geom.define(domain, real_box, amrex::CoordSys::cartesian, is_periodic);

    // extract dx from the geometry object
    amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> dx = geom.CellSizeArray();

    // Nghost = number of ghost cells for each array
    int Nghost = 0;

    // Ncomp = number of components for each array
    int Ncomp = 1;

    // How Boxes are distrubuted among MPI processes
    amrex::DistributionMapping dm(ba);

    // we allocate two phi multifabs; one will store the old state, the other the new.
    amrex::MultiFab phi(ba, dm, Ncomp, Nghost);

    // **********************************
    // INITIALIZE DATA
    // **********************************

    // loop over boxes
    for (amrex::MFIter mfi(phi); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.validbox();

        const amrex::Array4<amrex::Real>& phi_ptr = phi.array(mfi);

        // set phi = 1 + e^(-(r-0.5)^2)
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            // **********************************
            // SET VALUES FOR EACH CELL
            // **********************************

            amrex::Real x = (i+0.5) * dx[0];
            amrex::Real y = (j+0.5) * dx[1];
            amrex::Real z = (k+0.5) * dx[2];
            amrex::Real rsquared = ((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)+(z-0.5)*(z-0.5))/0.01;
            phi_ptr(i,j,k) = 1. + std::exp(-rsquared);
        });
    }

    // **********************************
    // WRITE INITIAL DATA TO PLOT FILE
    // **********************************

    // Write a plotfile of the initial data 
    int step = 0;
    amrex::Real time = 0.;
    // arguments
    // 1: name of plotfile
    // 2: MultiFab containing data to plot
    // 3: variables names
    // 4: geometry object
    // 5: "time" of plotfile; not relevant in this example
    // 6: "time step" of plotfile; not relevant in this example
    WriteSingleLevelPlotfile("phi", phi, {"phi"}, geom, time, step);

    // **********************************
    // COMPUTE FFT
    // **********************************

    // number of points in the domain
    long npts = domain.numPts();
    amrex::Real sqrtnpts = std::sqrt(npts);


    // **********************************
    // WRITE FFT TO PLOT FILE
    // **********************************




    

    // Call the timer again and compute the maximum difference between the start time 
    // and stop time over all processors
    amrex::Real stop_time = amrex::ParallelDescriptor::second() - start_time;
    amrex::ParallelDescriptor::ReduceRealMax(stop_time);
    amrex::Print() << "Run time = " << stop_time << std::endl;
    
    amrex::Finalize();
}


