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

        pp.get("prob_hi_x",prob_hi_x);
        pp.get("prob_hi_y",prob_hi_y);
        pp.get("prob_hi_z",prob_hi_z);

        pp.get("max_grid_size",max_grid_size);
    }

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
    RealBox real_box({ AMREX_D_DECL(       0.,        0.,        0.)},
                     { AMREX_D_DECL(prob_hi_x, prob_hi_y, prob_hi_z)} );

    // periodic in all direction
    Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(1,1,1)};

    // This defines a Geometry object
    geom.define(domain, real_box, CoordSys::cartesian, is_periodic);

    // extract dx from the geometry object
    GpuArray<Real,AMREX_SPACEDIM> dx = geom.CellSizeArray();

    // MultiFab storage for phi, and the real and imaginary parts of the dft
    MultiFab phi         (ba, dm, 1, 0);
    MultiFab phi_dft_real(ba, dm, 1, 0);
    MultiFab phi_dft_imag(ba, dm, 1, 0);

    // **********************************
    // INITIALIZE PHI
    // **********************************

    double omega = M_PI/2.0;

    // loop over boxes
    for (MFIter mfi(phi); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

        const Array4<Real>& phi_ptr = phi.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            // **********************************
            // SET VALUES FOR EACH CELL
            // **********************************

            Real x = (i+0.5) * dx[0];
            Real y = (AMREX_SPACEDIM>=2) ? (j+0.5) * dx[1] : 0.;
            Real z = (AMREX_SPACEDIM==3) ? (k+0.5) * dx[2] : 0.;
            phi_ptr(i,j,k) = std::sin(4*M_PI*x/prob_hi_x + omega);
            if (AMREX_SPACEDIM >= 2) {
                phi_ptr(i,j,k) *= std::sin(124*M_PI*y/prob_hi_y + omega);
            }
            if (AMREX_SPACEDIM == 3) {
                phi_ptr(i,j,k) *= std::sin(2*M_PI*z/prob_hi_z + omega);
            }
        });
    }

    // **********************************
    // COPY PHI INTO A MULTIFAB WITH ONE BOX
    // **********************************

    // create a new BoxArray and DistributionMapping for a MultiFab with 1 box
    BoxArray ba_onegrid(geom.Domain());
    DistributionMapping dm_onegrid(ba_onegrid);

    // storage for phi and the dft
    MultiFab phi_onegrid         (ba_onegrid, dm_onegrid, 1, 0);
    MultiFab phi_dft_real_onegrid(ba_onegrid, dm_onegrid, 1, 0);
    MultiFab phi_dft_imag_onegrid(ba_onegrid, dm_onegrid, 1, 0);

    // copy phi into phi_onegrid
    phi_onegrid.ParallelCopy(phi, 0, 0, 1);

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

    for (MFIter mfi(phi_onegrid); mfi.isValid(); ++mfi) {

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
                   phi_onegrid[mfi].dataPtr(),
                   reinterpret_cast<FFTcomplex*>
                   (spectral_field.back()->dataPtr()),
                   FFTW_ESTIMATE);
#elif (AMREX_SPACEDIM == 2)
      fplan = fftw_plan_dft_r2c_2d(fft_size[1], fft_size[0],
                   phi_onegrid[mfi].dataPtr(),
                   reinterpret_cast<FFTcomplex*>
                   (spectral_field.back()->dataPtr()),
                   FFTW_ESTIMATE);
#elif (AMREX_SPACEDIM == 3)
      fplan = fftw_plan_dft_r2c_3d(fft_size[2], fft_size[1], fft_size[0],
                   phi_onegrid[mfi].dataPtr(),
                   reinterpret_cast<FFTcomplex*>
                   (spectral_field.back()->dataPtr()),
                   FFTW_ESTIMATE);
#endif

#endif

      forward_plan.push_back(fplan);
    }

    ParallelDescriptor::Barrier();

    // ForwardTransform
    for (MFIter mfi(phi_onegrid); mfi.isValid(); ++mfi) {
      int i = mfi.LocalIndex();
#ifdef AMREX_USE_CUDA
      cufftSetStream(forward_plan[i], Gpu::gpuStream());
      cufftResult result = cufftExecD2Z(forward_plan[i],
                    phi_onegrid[mfi].dataPtr(),
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

    // copy data to a full-sized MultiFab
    // this involves copying the complex conjugate from the half-sized field
    // into the appropriate place in the full MultiFab
    for (MFIter mfi(phi_dft_real_onegrid); mfi.isValid(); ++mfi) {

      Array4< GpuComplex<Real> > spectral = (*spectral_field[0]).array();

      Array4<Real> const& realpart = phi_dft_real_onegrid.array(mfi);
      Array4<Real> const& imagpart = phi_dft_imag_onegrid.array(mfi);

      Box bx = mfi.fabbox();

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
              // copy value
              realpart(i,j,k) = spectral(i,j,k).real();
              imagpart(i,j,k) = spectral(i,j,k).imag();
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

              realpart(i,j,k) =  spectral(iloc,jloc,kloc).real();
              imagpart(i,j,k) = -spectral(iloc,jloc,kloc).imag();
          }

          realpart(i,j,k) /= sqrtnpts;
          imagpart(i,j,k) /= sqrtnpts;
      });
    }

    // destroy fft plan
    for (int i = 0; i < forward_plan.size(); ++i) {
#ifdef AMREX_USE_CUDA
        cufftDestroy(forward_plan[i]);
#else
        fftw_destroy_plan(forward_plan[i]);
#endif
    }

    // **********************************
    // SHIFT DATA
    // **********************************

    // zero_avg=0 means set the k=0 value to zero,
    // otherwise it sets the k=0 value to the average value of the signal in real space
    int zero_avg = 1;

    // shift data
    ShiftFFT(phi_dft_real_onegrid,geom,zero_avg);
    ShiftFFT(phi_dft_imag_onegrid,geom,zero_avg);

    // **********************************
    // COPY DFT INTO THE DISTRIBUTED MULTIFAB
    // **********************************

    phi_dft_real.ParallelCopy(phi_dft_real_onegrid, 0, 0, 1);
    phi_dft_imag.ParallelCopy(phi_dft_imag_onegrid, 0, 0, 1);

    // **********************************
    // WRITE DATA AND FFT TO PLOT FILE
    // **********************************

    // storage for magnitude and phase angle
    MultiFab phi_dft_magn(ba, dm, 1, 0);
    MultiFab phi_dft_phase(ba, dm, 1, 0);

    for (MFIter mfi(phi_dft_real); mfi.isValid(); ++mfi)
    {
        // Pointers to the magnitude, phase, real, and imaginary data
        const Array4<Real>& phi_dft_magn_ptr = phi_dft_magn.array(mfi);
        const Array4<Real>& phi_dft_phase_ptr = phi_dft_phase.array(mfi);
        const Array4<Real>& phi_dft_real_ptr = phi_dft_real.array(mfi);
        const Array4<Real>& phi_dft_imag_ptr = phi_dft_imag.array(mfi);

        const Box& bx = mfi.validbox();

        // Set the value of the magnitude and phase angle using the real and imaginary parts of the dft
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            double re = phi_dft_real_ptr(i,j,k);
            double im = phi_dft_imag_ptr(i,j,k);
            phi_dft_magn_ptr(i,j,k) = std::sqrt(re*re + im*im); // Here we want to store the values of the magnitude

            // Avoid division by zero
            if (re == 0.0) {
                if (im == 0.0){
                    phi_dft_phase_ptr(i,j,k) = 0.0;
                } else if (im > 0.0) {
                    phi_dft_phase_ptr(i,j,k) = M_PI/2.0;
                } else {
                     phi_dft_phase_ptr(i,j,k) = -M_PI/2.0;
                }
            } else {
                phi_dft_phase_ptr(i,j,k) = std::atan(im/re); // Here we want to store the values of the phase angle
            }
        });
     }

     // storage for variables to write to plotfile
     MultiFab plotfile(ba, dm, 5, 0);

     // copy phi, phi_dft_real, and phi_dft_imag, phi_dft_magn, and phi_dft_phase into plotfile
     MultiFab::Copy(plotfile, phi         , 0, 0, 1, 0);
     MultiFab::Copy(plotfile, phi_dft_real, 0, 1, 1, 0);
     MultiFab::Copy(plotfile, phi_dft_imag, 0, 2, 1, 0);
     MultiFab::Copy(plotfile, phi_dft_magn, 0, 3, 1, 0);
     MultiFab::Copy(plotfile, phi_dft_phase, 0, 4, 1, 0);

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
     WriteSingleLevelPlotfile("plt", plotfile, {"phi", "phi_dft_real", "phi_dft_imag","phi_dft_magn","phi_dft_phase"}, geom, time, step);

     // Call the timer again and compute the maximum difference between the start time
     // and stop time over all processors
     Real stop_time = ParallelDescriptor::second() - start_time;
     ParallelDescriptor::ReduceRealMax(stop_time);
     Print() << "Run time = " << stop_time << std::endl;

     }
     Finalize();
}


