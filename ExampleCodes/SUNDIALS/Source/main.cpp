#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_TimeIntegrator.H>

#include "myfunc.H"

using namespace amrex;

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    main_main();

    amrex::Finalize();
    return 0;
}

void main_main ()
{

    // **********************************
    // SIMULATION PARAMETERS

    // number of cells on each side of the domain
    int n_cell;

    // size of each box (or grid)
    int max_grid_size;

    // total steps in simulation
    int nsteps;

    // how often to write a plotfile
    int plot_int;

    // time step
    Real dt_out;

    // final time
    Real t_final;

    // inputs parameters
    {
        // ParmParse is way of reading inputs from the inputs file
        // pp.get means we require the inputs file to have it
        // pp.query means we optionally need the inputs file to have it - but we must supply a default here
        ParmParse pp;

        // We need to get n_cell from the inputs file - this is the number of cells on each side of
        //   a square (or cubic) domain.
        pp.get("n_cell",n_cell);

        // The domain is broken into boxes of size max_grid_size
        pp.get("max_grid_size",max_grid_size);

        // Default nsteps to 10, allow us to set it to something else in the inputs file
        nsteps = 10;
        pp.query("nsteps",nsteps);

        // Default plot_int to -1, allow us to set it to something else in the inputs file
        //  If plot_int < 0 then no plot files will be written
        plot_int = -1;
        pp.query("plot_int",plot_int);

        // time step
        pp.get("dt_out",dt_out);

        // time step
        pp.get("t_final",t_final);
    }

    // **********************************
    // SIMULATION SETUP

    // make BoxArray and Geometry
    // ba will contain a list of boxes that cover the domain
    // geom contains information such as the physical domain size,
    //               number of points in the domain, and periodicity
    BoxArray ba;
    Geometry geom;

    // AMREX_D_DECL means "do the first X of these, where X is the dimensionality of the simulation"
    IntVect dom_lo(AMREX_D_DECL(       0,        0,        0));
    IntVect dom_hi(AMREX_D_DECL(n_cell-1, n_cell-1, n_cell-1));

    // Make a single box that is the entire domain
    Box domain(dom_lo, dom_hi);

    // Initialize the boxarray "ba" from the single box "domain"
    ba.define(domain);

    // Break up boxarray "ba" into chunks no larger than "max_grid_size" along a direction
    ba.maxSize(max_grid_size);

    // This defines the physical box, [0,1] in each direction.
    RealBox real_box({AMREX_D_DECL( 0., 0., 0.)},
                     {AMREX_D_DECL( 1., 1., 1.)});

    // periodic in all direction
    Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(1,1,1)};

    // This defines a Geometry object
    geom.define(domain, real_box, CoordSys::cartesian, is_periodic);

    // extract dx from the geometry object
    GpuArray<Real,AMREX_SPACEDIM> dx = geom.CellSizeArray();

    // Nghost = number of ghost cells for each array
    int Nghost = 1;

    // Ncomp = number of components for each array
    int Ncomp = 1;

    // How Boxes are distrubuted among MPI processes
    DistributionMapping dm(ba);

    // we allocate two phi multifabs; one will store the old state, the other the new.
    Vector<MultiFab> state;
    state.push_back(MultiFab(ba, dm, Ncomp, Nghost));
    auto& phi = state[0];

    // time = starting time in the simulation
    Real time = 0.0;

    // **********************************
    // INITIALIZE DATA

    // loop over boxes
    for (MFIter mfi(phi); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

        const Array4<Real>& phiOld = phi.array(mfi);

        // set phi = 1 + e^(-(r-0.5)^2)
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            Real x = (i+0.5) * dx[0];
            Real y = (j+0.5) * dx[1];
#if (AMREX_SPACEDIM == 2)
            Real rsquared = ((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5))/0.01;
#elif (AMREX_SPACEDIM == 3)
            Real z= (k+0.5) * dx[2];
            Real rsquared = ((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)+(z-0.5)*(z-0.5))/0.01;
#endif
            phiOld(i,j,k) = 1. + std::exp(-rsquared);
        });
    }

    // Write a plotfile of the initial data if plot_int > 0
    if (plot_int > 0)
    {
        int step = 0;
        const std::string& pltfile = amrex::Concatenate("plt",step,5);
        WriteSingleLevelPlotfile(pltfile, phi, {"phi"}, geom, time, 0);
    }

    auto pre_rhs_function = [&](Vector<MultiFab>& S_data, const Real /* time */) {
        // fill periodic ghost cells
        S_data[0].FillBoundary(geom.periodicity());
    };

    auto rhs_function = [&](Vector<MultiFab>& S_rhs,
                            const Vector<MultiFab>& S_data, const Real /* time */) {

        // loop over boxes
        auto& phi_data = S_data[0];
        auto& phi_rhs  = S_rhs[0];
        for ( MFIter mfi(phi_data); mfi.isValid(); ++mfi )
        {
            const Box& bx = mfi.validbox();

            const Array4<const Real>& phi_array = phi_data.array(mfi);
            const Array4<Real>& phi_rhs_array = phi_rhs.array(mfi);

            // fill the right-hand-side for phi
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                phi_rhs_array(i,j,k) = ( (phi_array(i+1,j,k) - 2.*phi_array(i,j,k) + phi_array(i-1,j,k)) / (dx[0]*dx[0])
                                         +(phi_array(i,j+1,k) - 2.*phi_array(i,j,k) + phi_array(i,j-1,k)) / (dx[1]*dx[1])
#if (AMREX_SPACEDIM == 3)
                                         +(phi_array(i,j,k+1) - 2.*phi_array(i,j,k) + phi_array(i,j,k-1)) / (dx[2]*dx[2])
#endif
                                         );
            });
        }
    };

    TimeIntegrator<Vector<MultiFab>> integrator(state, time);
    integrator.set_pre_rhs_update(pre_rhs_function);
    integrator.set_rhs(rhs_function);

    for (int step = 1; step <= nsteps; ++step)
    {
        // Set output time
        time += dt_out;

        // Advance to output time
        integrator.evolve(state, time);

        // Tell the I/O Processor to write out which step we're doing
        amrex::Print() << "Advanced step " << step << "\n";

        // Write a plotfile of the current data (plot_int was defined in the inputs file)
        if (plot_int > 0 && step%plot_int == 0)
        {
            const std::string& pltfile = amrex::Concatenate("plt",step,5);
            WriteSingleLevelPlotfile(pltfile, phi, {"phi"}, geom, time, step);
        }

        if (time > t_final) break;
    }
}
