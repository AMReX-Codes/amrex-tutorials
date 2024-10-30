#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>

#include <AMReX_FFT_Poisson.H>

using namespace amrex;

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv); {

    BL_PROFILE("main");

    // **********************************
    // DECLARE SIMULATION PARAMETERS
    // **********************************

    // number of cells on each side of the domain
    int n_cell_x;
    int n_cell_y;
    int n_cell_z;

    // maximum grid size in each direction
    int max_grid_size_x;
    int max_grid_size_y;
    int max_grid_size_z;

    // physical dimensions of the domain
    Real prob_lo_x = 0.;
    Real prob_lo_y = 0.;
    Real prob_lo_z = 0.;
    Real prob_hi_x = 1.;
    Real prob_hi_y = 1.;
    Real prob_hi_z = 1.;

    // **********************************
    // READ PARAMETER VALUES FROM INPUTS FILE
    // **********************************
    {
        // ParmParse is way of reading inputs from the inputs file
        // pp.get means we require the inputs file to have it
        // pp.query means we optionally need the inputs file to have it - but you should supply a default value above

        ParmParse pp;

        pp.get("n_cell_x",n_cell_x);
        pp.get("n_cell_y",n_cell_y);
        pp.get("n_cell_z",n_cell_z);

        pp.get("max_grid_size_x",max_grid_size_x);
        pp.get("max_grid_size_y",max_grid_size_y);
        pp.get("max_grid_size_z",max_grid_size_z);

        pp.query("prob_lo_x",prob_lo_x);
        pp.query("prob_lo_y",prob_lo_y);
        pp.query("prob_lo_z",prob_lo_z);

        pp.query("prob_hi_x",prob_hi_x);
        pp.query("prob_hi_y",prob_hi_y);
        pp.query("prob_hi_z",prob_hi_z);
    }

    // Determine the domain length in each direction
    Real L_x = std::abs(prob_hi_x - prob_lo_x);
    Real L_y = std::abs(prob_hi_y - prob_lo_y);
    Real L_z = std::abs(prob_hi_z - prob_lo_z);

    // define lower and upper indices of domain
    IntVect dom_lo(AMREX_D_DECL(         0,          0,          0));
    IntVect dom_hi(AMREX_D_DECL(n_cell_x-1, n_cell_y-1, n_cell_z-1));

    // Make a single box that is the entire domain
    Box domain(dom_lo, dom_hi);

    // number of points in the domain
    long npts = domain.numPts();

    // Initialize the boxarray "ba" from the single box "domain"
    BoxArray ba(domain);

    // create IntVect of max_grid_size
    IntVect max_grid_size(AMREX_D_DECL(max_grid_size_x,max_grid_size_y,max_grid_size_z));

    // Break up boxarray "ba" into chunks no larger than "max_grid_size" along a direction
    ba.maxSize(max_grid_size);

    // How Boxes are distrubuted among MPI processes
    DistributionMapping dm(ba);

    // This defines the physical box size in each direction
    RealBox real_box({ AMREX_D_DECL(prob_lo_x, prob_lo_y, prob_lo_z)},
                     { AMREX_D_DECL(prob_hi_x, prob_hi_y, prob_hi_z)} );

    // periodic in all direction
    Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(1,1,1)};

    // geometry object for real data
    Geometry geom(domain, real_box, CoordSys::cartesian, is_periodic);

    // extract dx from the geometry object
    GpuArray<Real,AMREX_SPACEDIM> dx = geom.CellSizeArray();

    amrex::FFT::Poisson<MultiFab> my_poisson(geom);

    MultiFab rhs(ba,dm,1,0);
    MultiFab soln(ba,dm,1,0);

    for (MFIter mfi(rhs); mfi.isValid(); ++mfi) {

        Array4<Real> const& rhs_ptr = rhs.array(mfi);

        const Box& bx = mfi.fabbox();

        amrex::ParallelForRNG(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k, RandomEngine const& engine) noexcept
        {

            // **********************************
            // SET VALUES FOR EACH CELL
            // **********************************

            Real x = (i+0.5) * dx[0];
            Real y = (AMREX_SPACEDIM>=2) ? (j+0.5) * dx[1] : 0.;
            Real z = (AMREX_SPACEDIM==3) ? (k+0.5) * dx[2] : 0.;

            rhs_ptr(i,j,k) = std::exp(-10.*((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)+(z-0.5)*(z-0.5)));

        });
    }

    my_poisson.solve(soln,rhs);

    // storage for variables to write to plotfile
    MultiFab plotfile(ba, dm, 2, 0);

    // copy rhs and soln into plotfile
    MultiFab::Copy(plotfile, rhs , 0, 0, 1, 0);
    MultiFab::Copy(plotfile, soln, 0, 1, 1, 0);

    // time and step are dummy variables required to WriteSingleLevelPlotfile
    Real time = 0.;
    int step = 0;

    // arguments
    // 1: name of plotfile
    // 2: MultiFab containing data to plot
    // 3: variables names
    // 4: geometry object
    // 5: "time" of plotfile; not relevant in this example
    // 6: "time step" of plotfile; not relevant in this example
    WriteSingleLevelPlotfile("plt", plotfile, {"rhs", "soln"}, geom, time, step);

    } amrex::Finalize();
}
