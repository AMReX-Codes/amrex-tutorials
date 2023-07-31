
#include "myfunc.H"

void ShiftFFT(MultiFab& dft_onegrid, const Geometry& geom, const int& zero_avg) {

   /*
    Shifting rules:

    For each direction d, shift the data in the +d direction by n_cell[d]/2 and then modulo by n_cell[d].

    e.g. for an 8^3 domain

    Cell (7,2,3) is shifted to ( (7+4)%8, (2+4)%8, (3+4)%8 ) =  (3,6,7)

  */
  MultiFab dft_onegrid_temp;
  dft_onegrid_temp.define(dft_onegrid.boxArray(), dft_onegrid.DistributionMap(), 1, 0);

  MultiFab::Copy(dft_onegrid_temp,dft_onegrid,0,0,1,0);

  // Shift DFT by N/2+1 (pi)
  for (MFIter mfi(dft_onegrid); mfi.isValid(); ++mfi) {

    const Box& bx = mfi.tilebox();

    const Array4<Real>& dft = dft_onegrid.array(mfi);
    const Array4<Real>& dft_temp = dft_onegrid_temp.array(mfi);

    if (zero_avg == 1) {
#if (AMREX_SPACEDIM == 2)
      dft_temp(bx.smallEnd(0),bx.smallEnd(1),0) = 0.;
#elif (AMREX_SPACEDIM == 3)
      dft_temp(bx.smallEnd(0),bx.smallEnd(1),bx.smallEnd(2)) = 0.;
#endif
    }

    int nx = bx.length(0);    
    int nxh = nx/2;
    int ny = bx.length(1);
    int nyh = ny/2;
#if (AMREX_SPACEDIM == 3)
    int nz = bx.length(2);
    int nzh = nz/2;
#endif

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      int ip = (i+nxh)%nx;
      int jp = (j+nyh)%ny;
      int kp = 0;
#if (AMREX_SPACEDIM == 3)
      kp = (k+nzh)%nz;
#endif
      dft(ip,jp,kp) = dft_temp(i,j,k);
    });

  }

}


#ifdef AMREX_USE_CUDA
std::string cufftErrorToString (const cufftResult& err)
{
    switch (err) {
    case CUFFT_SUCCESS:  return "CUFFT_SUCCESS";
    case CUFFT_INVALID_PLAN: return "CUFFT_INVALID_PLAN";
    case CUFFT_ALLOC_FAILED: return "CUFFT_ALLOC_FAILED";
    case CUFFT_INVALID_TYPE: return "CUFFT_INVALID_TYPE";
    case CUFFT_INVALID_VALUE: return "CUFFT_INVALID_VALUE";
    case CUFFT_INTERNAL_ERROR: return "CUFFT_INTERNAL_ERROR";
    case CUFFT_EXEC_FAILED: return "CUFFT_EXEC_FAILED";
    case CUFFT_SETUP_FAILED: return "CUFFT_SETUP_FAILED";
    case CUFFT_INVALID_SIZE: return "CUFFT_INVALID_SIZE";
    case CUFFT_UNALIGNED_DATA: return "CUFFT_UNALIGNED_DATA";
    default: return std::to_string(err) + " (unknown error code)";
    }
}
#endif
