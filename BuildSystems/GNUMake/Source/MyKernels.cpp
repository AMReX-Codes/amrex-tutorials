#include "MyKernels.H"

void fill (amrex::MultiFab& mf)
{
    auto const& ma = mf.arrays();
    amrex::ParallelFor(mf, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
    {
        ma[b](i,j,k) = amrex::Real(i+j+k);
    });
}
