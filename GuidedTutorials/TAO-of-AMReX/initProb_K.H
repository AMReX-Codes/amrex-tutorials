#ifndef INIT_PROB_K_H_
#define INIT_PROB_K_H_

#include <AMReX_FArrayBox.H>

AMREX_GPU_DEVICE AMREX_INLINE
void actual_init_poisson (int i, int j, int k,
                          amrex::Array4<amrex::Real> const& rhs,
                          amrex::Array4<amrex::Real> const& exact,
                          amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& prob_lo,
                          amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx)
{
    constexpr amrex::Real tpi = 2.*3.1415926535897932;
    constexpr amrex::Real fpi = 4.*3.1415926535897932;
    constexpr amrex::Real fac = tpi*tpi*3.;
    amrex::Real x = prob_lo[0] + dx[0] * (i + 0.5);
    amrex::Real y = prob_lo[1] + dx[1] * (j + 0.5);
    amrex::Real z = prob_lo[2] + dx[2] * (k + 0.5);

    exact(i,j,k) = (std::sin(tpi*x) * std::sin(tpi*y) * std::sin(tpi*z))
           + .25 * (std::sin(fpi*x) * std::sin(fpi*y) * std::sin(fpi*z));
                
    rhs(i,j,k) = -fac * (std::sin(tpi*x) * std::sin(tpi*y) * std::sin(tpi*z))
                 -fac * (std::sin(fpi*x) * std::sin(fpi*y) * std::sin(fpi*z));
}

#endif
