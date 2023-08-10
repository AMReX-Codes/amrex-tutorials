#include <AMReX.H>

#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>

#ifdef AMREX_USE_CUDA
#include <cufft.h>
#else
#include <fftw3.h>
#include <fftw3-mpi.h>
#endif

#include "myfunc.H"

using namespace amrex;

int main (int argc, char* argv[])
{

    Initialize(argc,argv);
    {

    // store the current time so we can later compute total run time.
    Real start_time = ParallelDescriptor::second();

    // **********************************
    // DECLARE SIMULATION PARAMETERS
    // **********************************

    // number of cells on each side of the domain
    int n_cell_x;
    int n_cell_y;
    int n_cell_z;

    // dimensions of the domain
    Real prob_lo_x;
    Real prob_lo_y;
    Real prob_lo_z;

    Real prob_hi_x;
    Real prob_hi_y;
    Real prob_hi_z;

    // This is the largest size a grid can be
    int max_grid_size;

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

        pp.get("prob_lo_x",prob_lo_x);
        pp.get("prob_lo_y",prob_lo_y);
        pp.get("prob_lo_z",prob_lo_z);

        pp.get("prob_hi_x",prob_hi_x);
        pp.get("prob_hi_y",prob_hi_y);
        pp.get("prob_hi_z",prob_hi_z);

        pp.get("max_grid_size",max_grid_size);
    }

    // Determine the domain length in each direction
    Real L_x = std::abs(prob_hi_x - prob_lo_x);
    Real L_y = std::abs(prob_hi_y - prob_lo_y);
    Real L_z = std::abs(prob_hi_z - prob_lo_z);

    // **********************************
    // DEFINE SIMULATION SETUP AND GEOMETRY
    // **********************************
    // make BoxArray and Geometry
    // ba will contain a list of boxes that cover the domain
    // geom contains information such as the physical domain size,
    // number of points in the domain, and periodicity
    BoxArray ba;
    Geometry geom;

    // define lower and upper indices
    IntVect dom_lo(AMREX_D_DECL(         0,          0,          0));
    IntVect dom_hi(AMREX_D_DECL(n_cell_x-1, n_cell_y-1, n_cell_z-1));

    // Make a single box that is the entire domain
    Box domain(dom_lo, dom_hi);

    // Initialize the boxarray "ba" from the single box "domain"
    ba.define(domain);

    // Break up boxarray "ba" into chunks no larger than "max_grid_size" along a direction
    ba.maxSize(max_grid_size);

    // How Boxes are distrubuted among MPI processes
    DistributionMapping dm(ba);

    // This defines the physical box size in each direction
    RealBox real_box({ AMREX_D_DECL(prob_lo_x, prob_lo_y, prob_lo_z)},
                     { AMREX_D_DECL(prob_hi_x, prob_hi_y, prob_hi_z)} );

    // periodic in all direction
    Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(1,1,1)};

    // This defines a Geometry object
    geom.define(domain, real_box, CoordSys::cartesian, is_periodic);

    // extract dx from the geometry object
    GpuArray<Real,AMREX_SPACEDIM> dx = geom.CellSizeArray();

    // MultiFab storage for rhs and the solution to the Poisson equation, lap(soln) = rhs
    MultiFab rhs(ba, dm, 1, 0);
    MultiFab soln(ba, dm, 1, 0);

    // **********************************
    // INITIALIZE RHS
    // **********************************

    // loop over boxes
    for (MFIter mfi(rhs); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

        const Array4<Real>& rhs_ptr = rhs.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            // **********************************
            // SET VALUES FOR EACH CELL
            // **********************************

            Real x = (i+0.5) * dx[0];
            Real y = (AMREX_SPACEDIM==2) ? (j+0.5) * dx[1] : 0.;
            Real z = (AMREX_SPACEDIM==3) ? (k+0.5) * dx[2] : 0.;

            rhs_ptr(i,j,k) = std::exp(-10.*((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)+(z-0.5)*(z-0.5)));

        });
    }

    // **********************************
    // COPY RHS INTO A MULTIFAB WITH ONE BOX
    // **********************************

    // create a new BoxArray and DistributionMapping for a MultiFab with 1 box
    BoxArray ba_onegrid(geom.Domain());
    DistributionMapping dm_onegrid(ba_onegrid);

    // storage for rhs and soln
    MultiFab rhs_onegrid(ba_onegrid, dm_onegrid, 1, 0);
    MultiFab soln_onegrid(ba_onegrid, dm_onegrid, 1, 0);

    // copy rhs into rhs_onegrid
    rhs_onegrid.ParallelCopy(rhs, 0, 0, 1);

    // **********************************
    // COMPUTE FFT
    // **********************************

#ifdef AMREX_USE_CUDA
    using FFTplan = cufftHandle;
    using FFTcomplex = cuDoubleComplex;
#else
    using FFTplan = fftw_plan;
    using FFTcomplex = fftw_complex;
#endif

    // number of points in the domain
    long npts = domain.numPts();
    Real sqrtnpts = std::sqrt(npts);

    // contain to store FFT - note it is shrunk by approximately a half in x
    Vector<std::unique_ptr<BaseFab<GpuComplex<Real> > > > spectral_field;

    Vector<FFTplan> forward_plan;

    for (MFIter mfi(rhs_onegrid); mfi.isValid(); ++mfi) {

      // grab a single box
      Box realspace_bx = mfi.fabbox();

      // size of box
      IntVect fft_size = realspace_bx.length(); // This will be different for FFTs of complex data

      // this is the size of the box, except the 0th component is 'halved plus 1'
      IntVect spectral_bx_size = fft_size;
      spectral_bx_size[0] = fft_size[0]/2 + 1;

      // spectral box
      Box spectral_bx = Box(IntVect(0), spectral_bx_size - IntVect(1));

      spectral_field.emplace_back(new BaseFab<GpuComplex<Real> >(spectral_bx,1,
                                 The_Device_Arena()));
      spectral_field.back()->setVal<RunOn::Device>(0.0); // touch the memory

      FFTplan fplan;

#ifdef AMREX_USE_CUDA

#if (AMREX_SPACEDIM == 1)
      cufftResult result = cufftPlan1d(&fplan, fft_size[0], CUFFT_D2Z, 1);
      if (result != CUFFT_SUCCESS) {
          AllPrint() << " cufftplan1d forward failed! Error: "
                     << cufftErrorToString(result) << "\n";
      }
#elif (AMREX_SPACEDIM == 2)
      cufftResult result = cufftPlan2d(&fplan, fft_size[1], fft_size[0], CUFFT_D2Z);
      if (result != CUFFT_SUCCESS) {
          AllPrint() << " cufftplan2d forward failed! Error: "
                     << cufftErrorToString(result) << "\n";
      }
#elif (AMREX_SPACEDIM == 3)
      cufftResult result = cufftPlan3d(&fplan, fft_size[2], fft_size[1], fft_size[0], CUFFT_D2Z);
      if (result != CUFFT_SUCCESS) {
          AllPrint() << " cufftplan3d forward failed! Error: "
                     << cufftErrorToString(result) << "\n";
      }
#endif

#else // host

#if (AMREX_SPACEDIM == 1)
      fplan = fftw_plan_dft_r2c_1d(fft_size[0],
                   rhs_onegrid[mfi].dataPtr(),
                   reinterpret_cast<FFTcomplex*>
                   (spectral_field.back()->dataPtr()),
                   FFTW_ESTIMATE);
#elif (AMREX_SPACEDIM == 2)
      fplan = fftw_plan_dft_r2c_2d(fft_size[1], fft_size[0],
                   rhs_onegrid[mfi].dataPtr(),
                   reinterpret_cast<FFTcomplex*>
                   (spectral_field.back()->dataPtr()),
                   FFTW_ESTIMATE);
#elif (AMREX_SPACEDIM == 3)
      fplan = fftw_plan_dft_r2c_3d(fft_size[2], fft_size[1], fft_size[0],
                   rhs_onegrid[mfi].dataPtr(),
                   reinterpret_cast<FFTcomplex*>
                   (spectral_field.back()->dataPtr()),
                   FFTW_ESTIMATE);
#endif

#endif

      forward_plan.push_back(fplan);
    }

    ParallelDescriptor::Barrier();

    // ForwardTransform
    for (MFIter mfi(rhs_onegrid); mfi.isValid(); ++mfi) {
      int i = mfi.LocalIndex();
#ifdef AMREX_USE_CUDA
      cufftSetStream(forward_plan[i], Gpu::gpuStream());
      cufftResult result = cufftExecD2Z(forward_plan[i],
                    rhs_onegrid[mfi].dataPtr(),
                    reinterpret_cast<FFTcomplex*>
                    (spectral_field[i]->dataPtr()));
      if (result != CUFFT_SUCCESS) {
    AllPrint() << " forward transform using cufftExec failed! Error: "
           << cufftErrorToString(result) << "\n";
      }
#else
      fftw_execute(forward_plan[i]);
#endif
    }

    // Now we take the standard FFT and scale it by 1/k^2
    for (MFIter mfi(rhs_onegrid); mfi.isValid(); ++mfi)
    {
        Array4< GpuComplex<Real> > spectral = (*spectral_field[0]).array();

        const Box& bx = mfi.fabbox();

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            // the "spectral" data only exists for i <= nx/2
            if (i <= bx.length(0)/2) {

                Real a = 2.*M_PI*i / L_x;
                Real b = 2.*M_PI*j / L_y;
                Real c = 2.*M_PI*k / L_z;

                // the values in the upper-half of the spectral array in y and z are here interpreted as negative wavenumbers
                if (j >= n_cell_y/2) b = 2.*M_PI*(n_cell_y-j) / L_y;
                if (k >= n_cell_z/2) c = 2.*M_PI*(n_cell_z-k) / L_z;

#if (AMREX_SPACEDIM == 1)
                Real k2 = -a*a;
#elif (AMREX_SPACEDIM == 2)
                Real k2 = -(a*a + b*b);
#elif (AMREX_SPACEDIM == 3)
                Real k2 = -(a*a + b*b + c*c);
#endif

                if (k2 != 0.) {
                    spectral(i,j,k) /= k2;
                } else {
                    spectral(i,j,k) *= 0.; // interpretation here is that the average value of the solution is zero
                }
            }

        });
     }

    // Now we have completed the fft and scaled each value by 1/k^2
    // The scaled fft is inside spectral_field
    // Take inverse fft of spectral_field and put it in soln_onegrid
    Vector<FFTplan> backward_plan;

    for (MFIter mfi(soln_onegrid); mfi.isValid(); ++mfi) {

       // grab a single box including ghost cell range
       Box realspace_bx = mfi.fabbox();

       // size of box including ghost cell range
       IntVect fft_size = realspace_bx.length(); // This will be different for hybrid FFT

       FFTplan bplan;

#ifdef AMREX_USE_CUDA

#if (AMREX_SPACEDIM == 1)
       cufftResult result = cufftPlan1d(&bplan, fft_size[0], CUFFT_Z2D, 1);
      if (result != CUFFT_SUCCESS) {
          AllPrint() << " cufftplan1d forward failed! Error: "
                     << cufftErrorToString(result) << "\n";
      }
#elif (AMREX_SPACEDIM == 2)
      cufftResult result = cufftPlan2d(&bplan, fft_size[1], fft_size[0], CUFFT_Z2D);
      if (result != CUFFT_SUCCESS) {
          AllPrint() << " cufftplan2d forward failed! Error: "
                     << cufftErrorToString(result) << "\n";
      }
#elif (AMREX_SPACEDIM == 3)
      cufftResult result = cufftPlan3d(&bplan, fft_size[2], fft_size[1], fft_size[0], CUFFT_Z2D);
      if (result != CUFFT_SUCCESS) {
          AllPrint() << " cufftplan3d forward failed! Error: "
                     << cufftErrorToString(result) << "\n";
      }
#endif

#else // host

#if (AMREX_SPACEDIM == 1)
      bplan = fftw_plan_dft_c2r_1d(fft_size[0],
                   reinterpret_cast<FFTcomplex*>
                   (spectral_field.back()->dataPtr()),
                   soln_onegrid[mfi].dataPtr(),
                   FFTW_ESTIMATE);
#elif (AMREX_SPACEDIM == 2)
      bplan = fftw_plan_dft_c2r_2d(fft_size[1], fft_size[0],
                   reinterpret_cast<FFTcomplex*>
                   (spectral_field.back()->dataPtr()),
                   soln_onegrid[mfi].dataPtr(),
                   FFTW_ESTIMATE);
#elif (AMREX_SPACEDIM == 3)
      bplan = fftw_plan_dft_c2r_3d(fft_size[2], fft_size[1], fft_size[0],
                   reinterpret_cast<FFTcomplex*>
                   (spectral_field.back()->dataPtr()),
                   soln_onegrid[mfi].dataPtr(),
                   FFTW_ESTIMATE);
#endif

#endif

    backward_plan.push_back(bplan);// This adds an instance of bplan to the end of backward_plan
    }

    for (MFIter mfi(soln_onegrid); mfi.isValid(); ++mfi) {
      int i = mfi.LocalIndex();
#ifdef AMREX_USE_CUDA
      cufftSetStream(backward_plan[i], Gpu::gpuStream());
      cufftResult result = cufftExecZ2D(backward_plan[i],
                           reinterpret_cast<FFTcomplex*>
                           (spectral_field[i]->dataPtr()),
                           soln_onegrid[mfi].dataPtr());
       if (result != CUFFT_SUCCESS) {
           AllPrint() << " inverse transform using cufftExec failed! Error: "
                      << cufftErrorToString(result) << "\n";
       }
#else
      fftw_execute(backward_plan[i]);
#endif

    }

    // copy contents of soln_onegrid into soln
    soln.ParallelCopy(soln_onegrid, 0, 0, 1);

    // Must divide each point by the total number of points in the domain for properly scaled inverse FFT
    Real volinv = (AMREX_SPACEDIM == 2) ? 1. / (n_cell_x*n_cell_y) : 1. / (n_cell_x*n_cell_y*n_cell_z);
    soln.mult(volinv);

    // destroy fft plan
    for (int i = 0; i < forward_plan.size(); ++i) {
#ifdef AMREX_USE_CUDA
        cufftDestroy(forward_plan[i]);
#else
        fftw_destroy_plan(forward_plan[i]);
#endif
     }

    // destroy ifft plan
    for (int i = 0; i < backward_plan.size(); ++i) {
#ifdef AMREX_USE_CUDA
        cufftDestroy(backward_plan[i]);
#else
        fftw_destroy_plan(backward_plan[i]);
#endif

    }

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

    // Call the timer again and compute the maximum difference between the start time
    // and stop time over all processors
    Real stop_time = ParallelDescriptor::second() - start_time;
    ParallelDescriptor::ReduceRealMax(stop_time);
    Print() << "Run time = " << stop_time << std::endl;

    }
    Finalize();
}


