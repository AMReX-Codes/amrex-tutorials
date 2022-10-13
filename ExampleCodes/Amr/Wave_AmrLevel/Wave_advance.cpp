#include "AmrLevelWave.H"

using namespace amrex;

Real
AmrLevelWave::advance (Real time, Real dt, int iteration, int ncycle)
{
    // At the beginning of step, we make the new data from previous step the
    // old data of this step.
    for (int k = 0; k < NUM_STATE_TYPE; ++k) {
        state[k].allocOldData();
        state[k].swapTimeLevels(dt);
    }

    AmrLevel::RK(rk_order, State_Type, time, dt, iteration, ncycle,
    [&] (int /*stage*/, MultiFab& dSdt, MultiFab const& S, Real /*dtsub*/)
    {
        computeRHS(dSdt, S);
    });

    return dt;
}

void
AmrLevelWave::computeRHS (MultiFab& dSdt, MultiFab const& S)
{
    const auto dxinv = Geom().InvCellSizeArray();
    AMREX_D_TERM(Real dx2inv = dxinv[0]*dxinv[0];,
                 Real dy2inv = dxinv[1]*dxinv[1];,
                 Real dz2inv = dxinv[2]*dxinv[2]);
    auto const& sa = S.const_arrays();
    auto const& sdot = dSdt.arrays();
    amrex::ParallelFor(S,
    [=] AMREX_GPU_DEVICE (int bi, int i, int j, int k) noexcept
    {
        auto const& s = sa[bi];
        auto const& f = sdot[bi];
        f(i,j,k,0) = s(i,j,k,1);
        AMREX_D_TERM(Real lapx = dx2inv*(-2.5*s(i,j,k,0) + (4./3.)*(s(i-1,j,k,0)+s(i+1,j,k,0))
                                         -                (1./12.)*(s(i-2,j,k,0)+s(i+2,j,k,0)));,
                     Real lapy = dy2inv*(-2.5*s(i,j,k,0) + (4./3.)*(s(i,j-1,k,0)+s(i,j+1,k,0))
                                         -                (1./12.)*(s(i,j-2,k,0)+s(i,j+2,k,0)));,
                     Real lapz = dz2inv*(-2.5*s(i,j,k,0) + (4./3.)*(s(i,j,k-1,0)+s(i,j,k+1,0))
                                         -                (1./12.)*(s(i,j,k-2,0)+s(i,j,k+2,0))));
        f(i,j,k,1) = AMREX_D_TERM(lapx, +lapy, +lapz);
    });
    Gpu::streamSynchronize();
}
