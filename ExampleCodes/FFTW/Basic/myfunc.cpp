
#include "myfunc.H"

void ShiftFFT(MultiFab& dft_onegrid, const Geometry& geom, const int& zero_avg) {

  /*
    Shifting rules:

    For domains from (0,0,0) to (Nx-1,Ny-1,Nz-1)

    For any cells with i index >= Nx/2, these values are complex conjugates of the corresponding
    entry where (Nx-i,Ny-j,Nz-k) UNLESS that index is zero, in which case you use 0.

    e.g. for an 8^3 domain, any cell with i index

    Cell (6,2,3) is complex conjugate of (2,6,5)

    Cell (4,1,0) is complex conjugate of (4,7,0)  (note that the FFT is computed for 0 <= i <= Nx/2)
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
    int nxh = (nx+1)/2;
    int ny = bx.length(1);
    int nyh = (ny+1)/2;
#if (AMREX_SPACEDIM == 3)
    int nz = bx.length(2);
    int nzh = (nz+1)/2;
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
