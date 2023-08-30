#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_GpuComplex.H>
#include <AMReX_PlotFileUtil.H>

#include <heffte.h>

using namespace amrex;
//using namespace HEFFTE;

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

    int max_grid_size_x;
    int max_grid_size_y;
    int max_grid_size_z;

    // dimensions of the domain
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
        

    // define lower and upper indices
    IntVect dom_lo(AMREX_D_DECL(         0,          0,          0));
    IntVect dom_hi(AMREX_D_DECL(n_cell_x-1, n_cell_y-1, n_cell_z-1));

    // Make a single box that is the entire domain
    Box domain(dom_lo, dom_hi);

    // number of points in the domain
    long npts = domain.numPts();
    
    IntVect max_grid_size(AMREX_D_DECL(max_grid_size_x,max_grid_size_y,max_grid_size_z));
    
    // This defines the physical box size in each direction
    RealBox real_box({ AMREX_D_DECL(prob_lo_x, prob_lo_y, prob_lo_z)},
                     { AMREX_D_DECL(prob_hi_x, prob_hi_y, prob_hi_z)} );
    
    // periodic in all direction
    Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(1,1,1)};

    Geometry geom;
    geom.define(domain, real_box, CoordSys::cartesian, is_periodic);

    // Initialize the boxarray "ba" from the single box "domain"
    BoxArray ba;
    ba.define(domain);

    // Break up boxarray "ba" into chunks no larger than "max_grid_size" along a direction
    ba.maxSize(max_grid_size);

    // How Boxes are distrubuted among MPI processes
    DistributionMapping dm(ba);

    Box my_domain;
    int my_boxid;
    {
        for (int i = 0; i < ba.size(); ++i) {
            Box b = ba[i];
            // each MPI rank has its own my_domain Box and my_boxid ID
            if (ParallelDescriptor::MyProc() == dm[i]) {
                my_domain = b;
                my_boxid = i;
            }
        }
    }
    
    MultiFab real_field(ba,dm,1,0,MFInfo().SetArena(The_Device_Arena()));

    // check to make sure each MPI rank has exactly 1 box
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(real_field.local_size() == 1, "Must have one Box per process");
    
    for (MFIter mfi(real_field); mfi.isValid(); ++mfi) {
        Array4<Real> const& fab = real_field.array(mfi);
        amrex::ParallelForRNG(mfi.fabbox(),
        [=] AMREX_GPU_DEVICE (int i, int j, int k, RandomEngine const& engine) noexcept
        {
            fab(i,j,k) = amrex::Random(engine);
        });
    }

    Box r_local_box = my_domain;

    Print() << "r_local_box " << r_local_box << std::endl;
    Print() << "r_local_box.bigEnd(0) " << r_local_box.bigEnd(0) << std::endl;
    
    Box c_local_box = amrex::coarsen(r_local_box, IntVect(2,1,1));

    Print() << "c_local_box " << c_local_box << std::endl;
    Print() << "c_local_box.bigEnd(0) " << c_local_box.bigEnd(0) << std::endl;
    
    if (c_local_box.bigEnd(0) * 2 == r_local_box.bigEnd(0)) {
        // this avoids overlap for the cases when one or more r_local_box's
        // have an even cell index in the hi-x cell
        c_local_box.setBig(0,c_local_box.bigEnd(0)-1);
    }
    if (my_domain.bigEnd(0) == geom.Domain().bigEnd(0)) {
        // increase the size of boxes touching the hi-x domain by 1 in x
        // this is an (Nx x Ny x Nz) -> (Nx/2+1 x Ny x Nz) real-to-complex sizing 
        c_local_box.growHi(0,1);
    }

    BaseFab<GpuComplex<Real> > spectral_field(c_local_box, 1, The_Device_Arena());

#if (AMREX_SPACEDIM==2)

#ifdef AMREX_USE_CUDA
    heffte::fft2d_r2c<heffte::backend::cufft> fft
#elif AMREX_USE_HIP
    heffte::fft2d_r2c<heffte::backend::rocfft> fft
#else
    heffte::fft2d_r2c<heffte::backend::fftw> fft
#endif

#elif (AMREX_SPACEDIM==3)

#ifdef AMREX_USE_CUDA
    heffte::fft3d_r2c<heffte::backend::cufft> fft
#elif AMREX_USE_HIP
    heffte::fft3d_r2c<heffte::backend::rocfft> fft
#else
    heffte::fft3d_r2c<heffte::backend::fftw> fft
#endif

#endif
        ({{r_local_box.smallEnd(0),r_local_box.smallEnd(1),r_local_box.smallEnd(2)},
          {r_local_box.bigEnd(0)  ,r_local_box.bigEnd(1)  ,r_local_box.bigEnd(2)}},
         {{c_local_box.smallEnd(0),c_local_box.smallEnd(1),c_local_box.smallEnd(2)},
          {c_local_box.bigEnd(0)  ,c_local_box.bigEnd(1)  ,c_local_box.bigEnd(2)}},
         0, ParallelDescriptor::Communicator());

    using heffte_complex = typename heffte::fft_output<Real>::type;
    heffte_complex* spectral_data = (heffte_complex*) spectral_field.dataPtr();

    Real time = 0.;
    int step = 0;
    
    WriteSingleLevelPlotfile("plt_in", real_field, {"phi"}, geom, time, step);
    
    { BL_PROFILE("HEFFTE-total");
    {
        BL_PROFILE("ForwardTransform");
        fft.forward(real_field[my_boxid].dataPtr(), spectral_data);
    }

    {
        BL_PROFILE("BackwardTransform");
        fft.backward(spectral_data, real_field[my_boxid].dataPtr());
    }
    }

    // scale by 1/npts
    real_field.mult(1./npts);

    WriteSingleLevelPlotfile("plt_out", real_field, {"phi"}, geom, time, step);

    } amrex::Finalize();
}
