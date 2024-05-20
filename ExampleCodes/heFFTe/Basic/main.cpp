#include <heffte.h>

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_GpuComplex.H>
#include <AMReX_PlotFileUtil.H>
#include <cmath>

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
    Real prob_lo_x = -10.;
    Real prob_lo_y = -10.;
    Real prob_lo_z = -10.;
    Real prob_hi_x = 10.;
    Real prob_hi_y = 10.;
    Real prob_hi_z = 10.;

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


    // define lower and upper indices of domain
    IntVect dom_lo(AMREX_D_DECL(         0,          0,          0));
    IntVect dom_hi(AMREX_D_DECL(n_cell_x-1, n_cell_y-1, n_cell_z-1));

    // Make a single box that is the entire domain
    Box domain(dom_lo, dom_hi, amrex::IntVect::TheNodeVector());

    // number of points in the domain
    long npts = domain.numPts();
    Real sqrtnpts = std::sqrt(npts);

    // Initialize the boxarray "ba" from the single box "domain"
    BoxArray ba(domain);
Print() << "ba" " " << ba << " "<< std::endl;
    // create IntVect of max_grid_size
    IntVect max_grid_size(AMREX_D_DECL(max_grid_size_x,max_grid_size_y,max_grid_size_z));

    // Break up boxarray "ba" into chunks no larger than "max_grid_size" along a direction
    ba.maxSize(max_grid_size);
Print() << "ba.maxSize" " " << max_grid_size << " "<< std::endl;
Print() << "ba_after_maxSize" " " << ba << " "<< std::endl;
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

    MultiFab phi(ba,dm,1,0);

    // check to make sure each MPI rank has exactly 1 box
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(phi.local_size() == 1, "Must have one Box per MPI process");

    Real omega = M_PI/2.0;

    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {

        Array4<Real> const& fab = phi.array(mfi);

        const Box& bx = mfi.fabbox();

        amrex::ParallelForRNG(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k, RandomEngine const& engine) noexcept
        {

            // **********************************
            // SET VALUES FOR EACH CELL
            // **********************************

            Real x = prob_lo_x + (i+0.5) * dx[0];
            Real y = (AMREX_SPACEDIM>=2) ? prob_lo_y + (j+0.5) * dx[1] : 0.;
            Real z = (AMREX_SPACEDIM==3) ? prob_lo_z + (k+0.5) * dx[2] : 0.;
            fab(i,j,k) = std::exp(-x*x*0.5);
            if (AMREX_SPACEDIM >= 2) {
                fab(i,j,k) *= std::exp(-y*y*0.2);
            }
            if (AMREX_SPACEDIM == 3) {
                fab(i,j,k) *= std::sin(8*M_PI*z/prob_hi_z + omega);
            }

            // fab(i,j,k) = std::sqrt(x)*std::sqrt(y);

            // fab(i,j,k) = amrex::Random(engine);

        });
    }

    Real time = 0.;
    int step = 0;

    // write out phi to plotfile
    WriteSingleLevelPlotfile("phi", phi, {"phi"}, geom, time, step);

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

AllPrint() << "local_z" " " << local_box.smallEnd(2) << " " <<local_box.bigEnd(2) << "  " << c_local_box.smallEnd(2) << " " << c_local_box.bigEnd(2)<< " " << std::endl;
AllPrint() << "local_y" " " << local_box.smallEnd(1) << " " <<local_box.bigEnd(1) << "  " << c_local_box.smallEnd(1) << " " << c_local_box.bigEnd(1)<< " " << std::endl;
AllPrint() << "local_x" " " << local_box.smallEnd(0) << " " <<local_box.bigEnd(0) << "  " << c_local_box.smallEnd(0) << " " << c_local_box.bigEnd(0)<< " " << std::endl;

#endif

    using heffte_complex = typename heffte::fft_output<Real>::type;
    heffte_complex* spectral_data = (heffte_complex*) spectral_field.dataPtr();

    { BL_PROFILE("HEFFTE-total");
        {
            BL_PROFILE("ForwardTransform");
            fft.forward(phi[local_boxid].dataPtr(), spectral_data);
        }
        {
            BL_PROFILE("BackwardTransform");
            fft.backward(spectral_data, phi[local_boxid].dataPtr());
        }
    }

    // scale by 1/npts (both forward and inverse need 1/sqrtnpts scaling so I am doing it all here)
    phi.mult(1./npts);

    // **********************************
    // diagnostics
    // **********************************

    // create a BoxArray containing the fft boxes
    // by construction, these boxes correlate to the associated spectral_data
    // this we can copy the spectral data into this multifab since we know they are owned by the same MPI rank
    BoxArray fft_ba;
    {
        BoxList bl;
        bl.reserve(ba.size());

        for (int i = 0; i < ba.size(); ++i) {
            Box b = ba[i];

            Box r_box = b;
            Box c_box = amrex::coarsen(r_box, IntVect(AMREX_D_DECL(2,1,1)));

            // this avoids overlap for the cases when one or more r_box's
            // have an even cell index in the hi-x cell
            if (c_box.bigEnd(0) * 2 == r_box.bigEnd(0)) {
                c_box.setBig(0,c_box.bigEnd(0)-1);
            }

            // increase the size of boxes touching the hi-x domain by 1 in x
            // this is an (Nx x Ny x Nz) i-> (Nx/2+1 x Ny x Nz) real-to-complex sizing
            if (b.bigEnd(0) == geom.Domain().bigEnd(0)) {
                c_box.growHi(0,1);
            }
            bl.push_back(c_box);

        }
        fft_ba.define(std::move(bl));
    }

    // storage for real, imaginary, magnitude, and phase
    MultiFab fft_data(fft_ba,dm,4,0);

    // this copies the spectral data into a distributed MultiFab
    for (MFIter mfi(fft_data); mfi.isValid(); ++mfi) {

        Array4<Real> const& data = fft_data.array(mfi);
        Array4< GpuComplex<Real> > spectral = spectral_field.array();

        const Box& bx = mfi.fabbox();

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real re = spectral(i,j,k).real() / sqrtnpts;
            Real im = spectral(i,j,k).imag() / sqrtnpts;

            data(i,j,k,0) = re;
            data(i,j,k,1) = im;
            data(i,j,k,2) = std::sqrt(re*re + im*im);

            // Here we want to store the values of the phase angle
            // Avoid division by zero
            if (re == 0.0) {
                if (im == 0.0){
                    data(i,j,k,3) = 0.0;
                } else if (im > 0.0) {
                    data(i,j,k,3) = M_PI/2.0;
                } else {
                    data(i,j,k,3) = -M_PI/2.0;
                }
            } else {
                data(i,j,k,3) = std::atan(im/re);
            }
        });
    }

    // domain for fft data used to contruct a geometry object
    Box domain_fft = amrex::coarsen(domain, IntVect(AMREX_D_DECL(2,1,1)));
    // shrink by 1 in x in case there are an odd number of cells in the x-direction in domain
    if (domain_fft.bigEnd(0) * 2 == domain.bigEnd(0)) {
        domain_fft.setBig(0,domain_fft.bigEnd(0)-1);
    }
    // grow by 1 in the x-direction to match the size of the FFT
    domain_fft.growHi(0,1);

    Geometry geom_fft(domain_fft, real_box, CoordSys::cartesian, is_periodic);

    WriteSingleLevelPlotfile("fft_data", fft_data, {"real", "imag", "magitude", "phase"}, geom_fft, time, step);

    // **********************************
    // unpack data onto a 'full'-sized MultiFab
    // shift data so k=0 mode is at the center
    // **********************************

    BoxArray ba_onegrid(domain);
    DistributionMapping dm_onegrid(ba_onegrid);

    // real, imaginary, magnitude, phase
    MultiFab fft_data_onegrid(ba_onegrid, dm_onegrid, 4, 0);
    MultiFab fft_data_onegrid_shifted(ba_onegrid, dm_onegrid, 4, 0);

    // copy in real and imaginary parts
    fft_data_onegrid.ParallelCopy(fft_data, 0, 0, 2);

    // unpack data into a 'full'-sized MultiFab
    // this involves copying the complex conjugate from the half-sized field
    // into the appropriate place in the full MultiFab
    for (MFIter mfi(fft_data_onegrid); mfi.isValid(); ++mfi) {

        Array4<Real> const& data = fft_data_onegrid.array(mfi);

        const Box& bx = mfi.fabbox();

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
        /*
         Unpacking rules:

         For domains from (0,0,0) to (Nx-1,Ny-1,Nz-1)

         For any cells with i index > Nx/2, these values are complex conjugates of the corresponding
         entry where (Nx-i,Ny-j,Nz-k) UNLESS that index is zero, in which case you use 0.

         e.g. for an 8^3 domain, any cell with i index

         Cell (6,2,3) is complex conjugate of (2,6,5)

         Cell (4,1,0) is complex conjugate of (4,7,0)  (note that the FFT is computed for 0 <= i <= Nx/2)
        */
            if (i <= bx.length(0)/2) {
                // do nothing, data is already here
            } else {
                // copy complex conjugate
                int iloc = bx.length(0)-i;
                int jloc = 0;
                int kloc = 0;
#if (AMREX_SPACEDIM >= 2)
                jloc = (j == 0) ? 0 : bx.length(1)-j;
#endif
#if (AMREX_SPACEDIM == 3)
                kloc = (k == 0) ? 0 : bx.length(2)-k;
#endif

                data(i,j,k,0) =  data(iloc,jloc,kloc,0);
                data(i,j,k,1) = -data(iloc,jloc,kloc,1);
            }

            Real re = data(i,j,k,0);
            Real im = data(i,j,k,1);

            data(i,j,k,2) = std::sqrt(re*re + im*im);

            // Here we want to store the values of the phase angle
            // Avoid division by zero
            if (re == 0.0) {
                if (im == 0.0){
                    data(i,j,k,3) = 0.0;
                } else if (im > 0.0) {
                    data(i,j,k,3) = M_PI/2.0;
                } else {
                    data(i,j,k,3) = -M_PI/2.0;
                }
            } else {
                data(i,j,k,3) = std::atan(im/re);
            }

        });
    }

    // zero_avg=0 means set the k=0 value to zero,
    // otherwise it sets the k=0 value to the average value of the signal in real space
    //int zero_avg = 0;

    // zero k=0 mode
    /*
    if (zero_avg == 1) {
        for (MFIter mfi(fft_data_onegrid); mfi.isValid(); ++mfi) {

            const Box& bx = mfi.tilebox();

            Array4<Real> const& data = fft_data_onegrid.array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (i == 0 && j == 0 && k == 0) {
                    for (int comp=0; comp<4; ++comp) {
                        data(i,j,k,comp) = 0.;
                    }
                }
            });
        }
    }
     */
    // SHIFT DATA
    for (MFIter mfi(fft_data_onegrid); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.tilebox();

        const Array4<Real>& dft = fft_data_onegrid_shifted.array(mfi);
        const Array4<Real>& dft_temp = fft_data_onegrid.array(mfi);

        int nx = bx.length(0);
        int nxh = nx/2;
#if (AMREX_SPACEDIM >= 2)
        int ny = bx.length(1);
        int nyh = ny/2;
#endif
#if (AMREX_SPACEDIM == 3)
        int nz = bx.length(2);
        int nzh = nz/2;
#endif

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            int ip = (i+nxh)%nx;
            int jp = 0;
            int kp = 0;
#if (AMREX_SPACEDIM >= 2)
            jp = (j+nyh)%ny;
#endif
#if (AMREX_SPACEDIM == 3)
            kp = (k+nzh)%nz;
#endif
            for (int comp=0; comp<4; ++comp) {
                dft(ip,jp,kp,comp) = dft_temp(i,j,k,comp);
            }
        });

    }

    WriteSingleLevelPlotfile("plt", fft_data_onegrid_shifted,
                             {"phi_dft_real", "phi_dft_imag", "phi_dft_magitude", "phi_dft_phase"},
                             geom, time, step);

    } amrex::Finalize();
}
