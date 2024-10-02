#include <heffte.h>

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_GpuComplex.H>
#include <AMReX_PlotFileUtil.H>

using namespace amrex;

AMREX_ENUM(LaplacianType,
    Unknown, Exact, Discrete
);

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

    LaplacianType laplacian_type = LaplacianType::Unknown;

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

        pp.query_enum_case_insensitive("laplacian_type",laplacian_type);
    }

    if (laplacian_type == LaplacianType::Unknown) {
        amrex::Abort("Must specify exact or discrete Laplacian");
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

    MultiFab rhs(ba,dm,1,0);
    MultiFab soln(ba,dm,1,0);

    // check to make sure each MPI rank has exactly 1 box
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(rhs.local_size() == 1, "Must have one Box per MPI process");

    Real omega = M_PI/2.0;

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

    // since there is 1 MPI rank per box, here each MPI rank obtains its local box and the associated boxid
    Box local_box;
    int local_boxid;
    {
        for (int i = 0; i < ba.size(); ++i) {
            Box b = ba[i];
            // each MPI rank has its own local_box Box and local_boxid ID
            if (ParallelDescriptor::MyProc() == dm[i]) {
                local_box = b;
                local_boxid = i;
            }
        }
    }

    // now each MPI rank works on its own box
    // for real->complex fft's, the fft is stored in an (nx/2+1) x ny x nz dataset

    // start by coarsening each box by 2 in the x-direction
    Box c_local_box = amrex::coarsen(local_box, IntVect(AMREX_D_DECL(2,1,1)));

    // if the coarsened box's high-x index is even, we shrink the size in 1 in x
    // this avoids overlap between coarsened boxes
    if (c_local_box.bigEnd(0) * 2 == local_box.bigEnd(0)) {
        c_local_box.setBig(0,c_local_box.bigEnd(0)-1);
    }
    // for any boxes that touch the hi-x domain we
    // increase the size of boxes by 1 in x
    // this makes the overall fft dataset have size (Nx/2+1 x Ny x Nz)
    if (local_box.bigEnd(0) == geom.Domain().bigEnd(0)) {
        c_local_box.growHi(0,1);
    }

    // each MPI rank gets storage for its piece of the fft
    BaseFab<GpuComplex<Real> > spectral_field(c_local_box, 1, The_Device_Arena());

    // create real->complex fft objects with the appropriate backend and data about
    // the domain size and its local box size
#if (AMREX_SPACEDIM==2)
#ifdef AMREX_USE_CUDA
    heffte::fft2d_r2c<heffte::backend::cufft> fft
#elif AMREX_USE_HIP
    heffte::fft2d_r2c<heffte::backend::rocfft> fft
#else
    heffte::fft2d_r2c<heffte::backend::fftw> fft
#endif
        ({{local_box.smallEnd(0),local_box.smallEnd(1),0},
          {local_box.bigEnd(0)  ,local_box.bigEnd(1)  ,0}},
         {{c_local_box.smallEnd(0),c_local_box.smallEnd(1),0},
          {c_local_box.bigEnd(0)  ,c_local_box.bigEnd(1)  ,0}},
         0, ParallelDescriptor::Communicator());

#elif (AMREX_SPACEDIM==3)

#ifdef AMREX_USE_CUDA
    heffte::fft3d_r2c<heffte::backend::cufft> fft
#elif AMREX_USE_HIP
    heffte::fft3d_r2c<heffte::backend::rocfft> fft
#else
    heffte::fft3d_r2c<heffte::backend::fftw> fft
#endif
        ({{local_box.smallEnd(0),local_box.smallEnd(1),local_box.smallEnd(2)},
          {local_box.bigEnd(0)  ,local_box.bigEnd(1)  ,local_box.bigEnd(2)}},
         {{c_local_box.smallEnd(0),c_local_box.smallEnd(1),c_local_box.smallEnd(2)},
          {c_local_box.bigEnd(0)  ,c_local_box.bigEnd(1)  ,c_local_box.bigEnd(2)}},
         0, ParallelDescriptor::Communicator());

#endif

    Real start_step = static_cast<Real>(ParallelDescriptor::second());
    using heffte_complex = typename heffte::fft_output<Real>::type;
    heffte_complex* spectral_data = (heffte_complex*) spectral_field.dataPtr();

    { BL_PROFILE("HEFFTE-total");
        {
            BL_PROFILE("ForwardTransform");
            fft.forward(rhs[local_boxid].dataPtr(), spectral_data);
        }

        // Now we take the standard FFT and scale it by 1/k^2
        Array4< GpuComplex<Real> > spectral = spectral_field.array();

        int use_discrete = 0;
        if (laplacian_type == LaplacianType::Discrete) {
            use_discrete = 1;
        }

        ParallelFor(c_local_box, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            Real a = 2.*M_PI*i / L_x;
            Real b = 2.*M_PI*j / L_y;
            Real c = 2.*M_PI*k / L_z;

            // the values in the upper-half of the spectral array in y and z are here interpreted as negative wavenumbers
            if (j >= n_cell_y/2) b = 2.*M_PI*(n_cell_y-j) / L_y;
            if (k >= n_cell_z/2) c = 2.*M_PI*(n_cell_z-k) / L_z;

            Real k2;

            if (use_discrete == 0) {
#if (AMREX_SPACEDIM == 1)
                k2 = -a*a;
#elif (AMREX_SPACEDIM == 2)
                k2 = -(a*a + b*b);
#elif (AMREX_SPACEDIM == 3)
                k2 = -(a*a + b*b + c*c);
#endif
            } else {
#if (AMREX_SPACEDIM == 1)
                k2 = 2*(std::cos(a*dx[0])-1.)/(dx[0]*dx[0]);
#elif (AMREX_SPACEDIM == 2)
                k2 =  2*(std::cos(a*dx[0])-1.)/(dx[0]*dx[0]) + 2*(std::cos(b*dx[1])-1.)/(dx[1]*dx[1]);
#elif (AMREX_SPACEDIM == 3)
                k2 =  2*(std::cos(a*dx[0])-1.)/(dx[0]*dx[0]) + 2*(std::cos(b*dx[1])-1.)/(dx[1]*dx[1]) + 2*(std::cos(c*dx[2])-1.)/(dx[2]*dx[2]);
#endif
            }

            if (k2 != 0.) {
                spectral(i,j,k) /= k2;
            } else {
                spectral(i,j,k) *= 0.; // interpretation here is that the average value of the solution is zero
            }

        });


        {
            BL_PROFILE("BackwardTransform");
            fft.backward(spectral_data, soln[local_boxid].dataPtr());
        }
    }

    // scale by 1/npts (both forward and inverse need sqrt(npts) scaling so I am doing it all here)
    soln.mult(1./npts);

    Real end_step = static_cast<Real>(ParallelDescriptor::second());
    // amrex::Print() << "TIME IN SOLVE " << end_step - start_step << std::endl;

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
