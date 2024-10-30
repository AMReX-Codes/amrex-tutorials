#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>

#include <AMReX_FFT.H>

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

    amrex::FFT::R2C my_fft(domain);

    // storage for input phi and phi after forward-inverse transformation
    MultiFab phi(ba,dm,1,0);
    MultiFab phi_after(ba,dm,1,0);

    // initialize phi
    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {

        Array4<Real> const& phi_ptr = phi.array(mfi);

        const Box& bx = mfi.fabbox();

        amrex::ParallelForRNG(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k, RandomEngine const& engine) noexcept
        {

            // **********************************
            // SET VALUES FOR EACH CELL
            // **********************************

            Real x = (i+0.5) * dx[0];
            Real y = (AMREX_SPACEDIM>=2) ? (j+0.5) * dx[1] : 0.;
            Real z = (AMREX_SPACEDIM==3) ? (k+0.5) * dx[2] : 0.;

            phi_ptr(i,j,k) = std::exp(-10.*((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)+(z-0.5)*(z-0.5)));

        });
    }

    // create storage for the FFT
    Box cdomain = geom.Domain();
    cdomain.setBig(0,cdomain.length(0)/2);
    Geometry cgeom(cdomain, real_box, CoordSys::cartesian, is_periodic);
    auto cba = amrex::decompose(cdomain, ParallelContext::NProcsSub(),
                                {AMREX_D_DECL(true,true,false)});
    DistributionMapping cdm = amrex::FFT::detail::make_iota_distromap(cba.size());
    FabArray<BaseFab<GpuComplex<amrex::Real> > > phi_fft(cba, cdm, 1, 0);
    
    // we will copy the real and imaginary parts of the FFT to this MultiFab
    MultiFab phi_fft_realimag(cba,cdm,2,0);

    // compute the FFT
    my_fft.forward(phi,phi_fft);

    // copy FFT into a MultiFab for plotfile visualization
    for (MFIter mfi(phi_fft); mfi.isValid(); ++mfi) {

        Array4<GpuComplex<Real>> const& phi_fft_ptr = phi_fft.array(mfi);
        Array4<Real> phi_fft_realimag_ptr = phi_fft_realimag.array(mfi);

        const Box& bx = mfi.fabbox();

        amrex::ParallelForRNG(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k, RandomEngine const& engine) noexcept
        {
            phi_fft_realimag_ptr(i,j,k,0) = phi_fft_ptr(i,j,k).real();
            phi_fft_realimag_ptr(i,j,k,1) = phi_fft_ptr(i,j,k).imag();

        });
    }

    // compute the inverse FFT and store result in phi_after
    my_fft.backward(phi_after);

    // scale phi_after by 1/n_cells so it matches the original phi
    long n_cells = n_cell_x;
    if (AMREX_SPACEDIM >= 2) n_cells *= n_cell_y;
    if (AMREX_SPACEDIM >= 3) n_cells *= n_cell_z;
    phi_after.mult(1./n_cells);

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
    WriteSingleLevelPlotfile("plt", phi, {"phi"}, geom, time, step);
    WriteSingleLevelPlotfile("plt_after", phi_after, {"phi_after"}, geom, time, step);
    WriteSingleLevelPlotfile("plt_fft", phi_fft_realimag, {"phi_fft_real", "phi_fft_imag"}, cgeom, time, step);

    } amrex::Finalize();
}
