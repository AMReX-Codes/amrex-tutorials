/*
 * A simplified single file version of the HeatEquation_EX0_C exmaple.
 * This code is designed to be used with Demo_Tutorial.rst.
 *
 */


#include <AMReX.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_GMRES_Poisson.H>

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    {

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
    amrex::IntVect dom_lo(       0,        0,        0);
    amrex::IntVect dom_hi(n_cell-1, n_cell-1, n_cell-1);

    // Make a single box that is the entire domain
    amrex::Box domain(dom_lo, dom_hi);

    // Initialize the boxarray "ba" from the single box "domain"
    ba.define(domain);

    // Break up boxarray "ba" into chunks no larger than "max_grid_size" along a direction
    ba.maxSize(max_grid_size);

    // This defines the physical box, [0,1] in each direction.
    amrex::RealBox real_box({ 0., 0., 0.},
                            { 1., 1., 1.});

    // periodic in all direction
    amrex::Array<int,3> is_periodic{1,1,1};

    // This defines a Geometry object
    geom.define(domain, real_box, amrex::CoordSys::cartesian, is_periodic);

    // extract dx from the geometry object
    amrex::GpuArray<amrex::Real,3> dx = geom.CellSizeArray();

    // How Boxes are distrubuted among MPI processes
    amrex::DistributionMapping dm(ba);

    // we allocate two phi multifabs; one will store the old state, the other the new.
    amrex::MultiFab rhs(ba, dm, 1, 0);
    amrex::MultiFab phi(ba, dm, 1, 1);

    // **********************************
    // INITIALIZE DATA LOOP
    // **********************************

    // loop over boxes
    for (amrex::MFIter mfi(rhs); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.validbox();

        const amrex::Array4<amrex::Real>& rhs_p = rhs.array(mfi);

        // set rhs = 1 + e^(-(r-0.5)^2)
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            // **********************************
            // SET VALUES FOR EACH CELL
            // **********************************
            amrex::Real x = (i+0.5) * dx[0];
            amrex::Real y = (j+0.5) * dx[1];
            amrex::Real z = (k+0.5) * dx[2];
            amrex::Real rsquared = ((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)+(z-0.5)*(z-0.5))/0.01;
            rhs_p(i,j,k) = 1. + std::exp(-rsquared);
        });
    }

    WriteSingleLevelPlotfile("rhs", rhs, {"rhs"}, geom, 0., 0);

    amrex::GMRESPOISSON gmres_poisson(ba,dm,geom);

    // initial guess
    phi.setVal(0.);

    gmres_poisson.setVerbose(2);
    gmres_poisson.solve(phi, rhs, 1.e-12, 0.);

    WriteSingleLevelPlotfile("phi", phi, {"phi"}, geom, 0., 0);
    
    }
    amrex::Finalize();
    return 0;
}


