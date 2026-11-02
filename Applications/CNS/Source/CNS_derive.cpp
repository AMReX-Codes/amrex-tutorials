#include "CNS_derive.H"
#include "CNS.H"

using namespace amrex;

// Pressure from rho*e
void cns_derpres (MultiFab& der, int dcomp, int /*ncomp*/,
                  const MultiFab& data, const Geometry& /*geom*/,
                  Real /*time*/, const int* /*bcrec*/, int /*level*/)
{
    auto const& p = der.arrays();
    auto const& rhoe = data.const_arrays();
    const Real gamma = CNS::parm.eos_gamma;
    ParallelFor(der, der.nGrowVect(), [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) noexcept
    {
        p[b](i,j,k,dcomp) = (gamma-Real(1.0))*rhoe[b](i,j,k);
    });
}

// Velocity from density and momentum
void cns_dervel (MultiFab& der, int dcomp, int /*ncomp*/,
                 const MultiFab& data, const Geometry& /*geom*/,
                 Real /*time*/, const int* /*bcrec*/, int /*level*/)
{
    auto const& vel = der.arrays();
    auto const& dat = data.const_arrays();
    ParallelFor(der, der.nGrowVect(), AMREX_SPACEDIM,
    [=] AMREX_GPU_DEVICE (int b, int i, int j, int k, int n) noexcept
    {
        vel[b](i,j,k,dcomp+n) = dat[b](i,j,k,1+n) / dat[b](i,j,k,0);
    });
}
