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

    Geometry geom;
    BoxArray ba;
    DistributionMapping dm;
    {
        ParmParse pp;
        IntVect n_cell;
        IntVect max_grid_size;
        pp.get("n_cell", n_cell);
        pp.get("max_grid_size", max_grid_size);

        Box domain(IntVect(0),n_cell-IntVect(1));
        RealBox rb({0.,0.,0.},{1.,1.,1.});
        Array<int,3> is_periodic{1,1,1};
        geom.define(domain, rb, CoordSys::cartesian, is_periodic);

        ba.define(domain);
        ba.maxSize(max_grid_size);

        dm.define(ba);
    }

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
    
    Box c_local_box = amrex::coarsen(r_local_box, IntVect(2,1,1));

    Print() << "c_local_box " << c_local_box << std::endl;
    
    if (c_local_box.bigEnd(0) * 2 == r_local_box.bigEnd(0)) {
        c_local_box.setBig(0,c_local_box.bigEnd(0)-1);// to avoid overlap
    }
    if (my_domain.bigEnd(0) == geom.Domain().bigEnd(0)) {
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
    
    real_field.mult(1./(64.*64.*64.));

    WriteSingleLevelPlotfile("plt_out", real_field, {"phi"}, geom, time, step);

    } amrex::Finalize();
}
