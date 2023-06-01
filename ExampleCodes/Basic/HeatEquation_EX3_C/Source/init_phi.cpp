#include "myfunc.H"
#include "mykernel.H"

#include <AMReX_BCRec.H>
#include <AMReX_BCUtil.H>

using namespace amrex;

void init_phi(amrex::MultiFab& phi_new, amrex::Geometry const& geom)
{
    GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();
    GpuArray<Real, AMREX_SPACEDIM> prob_lo = geom.ProbLoArray();

    for (MFIter mfi(phi_new); mfi.isValid(); ++mfi)
    {
        const Box& vbx = mfi.validbox();
        auto const& phiNew = phi_new.array(mfi);
        amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    init_phi(i, j, k, phiNew, dx, prob_lo);
                });
    }
}
