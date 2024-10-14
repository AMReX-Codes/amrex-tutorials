#include <heffte.h>

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_GpuComplex.H>
#include <AMReX_PlotFileUtil.H>

using namespace amrex;

static_assert(AMREX_SPACEDIM == 3);

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

    // define lower and upper indices of domain
    IntVect dom_lo(AMREX_D_DECL(         0,          0,          0));
    IntVect dom_hi(AMREX_D_DECL(n_cell_x-1, n_cell_y-1, n_cell_z-1));

    // Make a single box that is the entire domain
    Box domain(dom_lo, dom_hi);

    // Initialize the boxarray "ba" from the single box "domain" There are
    // exactly nprocs boxes. The domain decomposition is done in the x- and
    // y-directions, but not the z-direction.
    BoxArray ba = amrex::decompose(domain, ParallelDescriptor::NProcs(),
                                   {AMREX_D_DECL(true,true,false)});

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

    for (MFIter mfi(rhs); mfi.isValid(); ++mfi) {

        Array4<Real> const& rhs_ptr = rhs.array(mfi);

        const Box& bx = mfi.fabbox();

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
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

    // Shift rhs so that its sum is zero.
    auto rhosum = rhs.sum(0);
    rhs.plus(-rhosum/geom.Domain().d_numPts(), 0, 1);

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
    using fft_r2c_t =
#ifdef AMREX_USE_CUDA
        heffte::fft2d_r2c<heffte::backend::cufft>;
#elif AMREX_USE_HIP
        heffte::fft2d_r2c<heffte::backend::rocfft>;
#else
        heffte::fft2d_r2c<heffte::backend::fftw>;
#endif

    auto lo = amrex::lbound(local_box);
    auto hi = amrex::ubound(local_box);
    auto len = amrex::length(local_box);
    auto clo = amrex::lbound(c_local_box);
    auto chi = amrex::ubound(c_local_box);
    auto clen = amrex::length(c_local_box);

    auto fft = std::make_unique<fft_r2c_t>
        (heffte::box3d({ lo.x,  lo.y, 0},
                       { hi.x,  hi.y, 0}),
         heffte::box3d({clo.x, clo.y, 0},
                       {chi.x, chi.y, 0}),
         0, ParallelDescriptor::Communicator());

    Real start_step = static_cast<Real>(ParallelDescriptor::second());
    using heffte_complex = typename heffte::fft_output<Real>::type;
    heffte_complex* spectral_data = (heffte_complex*) spectral_field.dataPtr();

    int batch_size = n_cell_z;
    Gpu::DeviceVector<heffte_complex> workspace(fft->size_workspace()*batch_size);

    { BL_PROFILE("HEFFTE-total");
        {
            BL_PROFILE("ForwardTransform");
            fft->forward(batch_size, rhs[local_boxid].dataPtr(), spectral_data,
                         workspace.data());
        }

        // Now we take the standard FFT and scale it by 1/k^2
        Array4< GpuComplex<Real> > spectral = spectral_field.array();

        FArrayBox tridiag_workspace(c_local_box,4);
        auto const& ald = tridiag_workspace.array(0);
        auto const& bd = tridiag_workspace.array(1);
        auto const& cud = tridiag_workspace.array(2);
        auto const& scratch = tridiag_workspace.array(3);

        Gpu::DeviceVector<Real> delzv(n_cell_z, dx[2]);
        auto const* delz = delzv.data();

        auto xybox = amrex::makeSlab(c_local_box, 2, 0);
        ParallelFor(xybox, [=] AMREX_GPU_DEVICE(int i, int j, int)
        {
            Real a = 2.*M_PI*i / L_x;
            Real b = 2.*M_PI*j / L_y;

            // the values in the upper-half of the spectral array in y and z are here interpreted as negative wavenumbers
            if (j >= n_cell_y/2) b = 2.*M_PI*(n_cell_y-j) / L_y;

            Real k2 = 2*(std::cos(a*dx[0])-1.)/(dx[0]*dx[0]) + 2*(std::cos(b*dx[1])-1.)/(dx[1]*dx[1]);

            // Tridiagonal solve with homogeneous Neumann
            for( int k=0; k<n_cell_z; k++) {
                if(k==0) {
                    ald(i,j,k) = 0.;
                    cud(i,j,k) = 2.0 /(delz[k]*(delz[k]+delz[k+1]));
                    bd(i,j,k) = k2 -ald(i,j,k)-cud(i,j,k);
                } else if (k == n_cell_z-1) {
                    ald(i,j,k) = 2.0 /(delz[k]*(delz[k]+delz[k-1]));
                    cud(i,j,k) = 0.;
                    bd(i,j,k) = k2 -ald(i,j,k)-cud(i,j,k);
                    if (i == 0 && j == 0) {
                        bd(i,j,k) *= 2.0;
                    }
               } else {
                    ald(i,j,k) = 2.0 /(delz[k]*(delz[k]+delz[k-1]));
                    cud(i,j,k) = 2.0 /(delz[k]*(delz[k]+delz[k+1]));
                    bd(i,j,k) = k2 -ald(i,j,k)-cud(i,j,k);
                }
            }

            scratch(i,j,0) = cud(i,j,0)/bd(i,j,0);
            spectral(i,j,0) = spectral(i,j,0)/bd(i,j,0);

            for (int k = 1; k < n_cell_z; k++) {
                if (k < n_cell_z-1){
                    scratch(i,j,k) = cud(i,j,k) / (bd(i,j,k) - ald(i,j,k) * scratch(i,j,k-1));
                }
                spectral(i,j,k) = (spectral(i,j,k) - ald(i,j,k) * spectral(i,j,k - 1)) / (bd(i,j,k) - ald(i,j,k) * scratch(i,j,k-1));
            }

            for (int k = n_cell_z - 2; k >= 0; k--) {
                spectral(i,j,k) -= scratch(i,j,k) * spectral(i,j,k + 1);
            }
        });
        Gpu::streamSynchronize();

        {
            BL_PROFILE("BackwardTransform");
            fft->backward(batch_size, spectral_data, soln[local_boxid].dataPtr(),
                          heffte::scale::full);
        }
    }

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

    {
        MultiFab phi(soln.boxArray(), soln.DistributionMap(), 1, 1);
        MultiFab res(soln.boxArray(), soln.DistributionMap(), 1, 0);
        MultiFab::Copy(phi, soln, 0, 0, 1, 0);
        phi.FillBoundary(geom.periodicity());
        auto const& res_ma = res.arrays();
        auto const& phi_ma = phi.const_arrays();
        auto const& rhs_ma = rhs.const_arrays();
        ParallelFor(res, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
        {
            auto const& phia = phi_ma[b];
            auto lap = (phia(i-1,j,k)-2.*phia(i,j,k)+phia(i+1,j,k)) / (dx[0]*dx[0])
                +      (phia(i,j-1,k)-2.*phia(i,j,k)+phia(i,j+1,k)) / (dx[1]*dx[1]);
            if (k == 0) {
                lap += (-phia(i,j,k)+phia(i,j,k+1)) / (dx[2]*dx[2]);
            } else if (k == n_cell_z-1) {
                lap += (phia(i,j,k-1)-phia(i,j,k)) / (dx[2]*dx[2]);
            } else {
                lap += (phia(i,j,k-1)-2.*phia(i,j,k)+phia(i,j,k+1)) / (dx[2]*dx[2]);
            }
            res_ma[b](i,j,k) = rhs_ma[b](i,j,k) - lap;
        });
        amrex::Print() << " rhs.min & max: " <<  rhs.min(0) << " " <<  rhs.max(0) << "\n"
                       << " res.min & max: " <<  res.min(0) << " " <<  res.max(0) << "\n";
    }

    } amrex::Finalize();
}
