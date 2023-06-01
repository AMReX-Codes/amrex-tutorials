#ifndef MY_KERNEL_H_
#define MY_KERNEL_H_

#include <AMReX_FArrayBox.H>

using namespace amrex;

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void init_phi (int i, int j, int k,
               amrex::Array4<amrex::Real> const& phi,
               GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx,
               GpuArray<amrex::Real,AMREX_SPACEDIM> const& prob_lo)
{
    using amrex::Real;;

    Real x = prob_lo[0] + (i+Real(0.5)) * dx[0];
    Real y = prob_lo[1] + (j+Real(0.5)) * dx[1];
#if (AMREX_SPACEDIM > 2)
    Real z = prob_lo[2] + (k+Real(0.5)) * dx[2];
    Real r2 = ((x-Real(0.25))*(x-Real(0.25))+(y-Real(0.25))*(y-Real(0.25))+(z-Real(0.25))*(z-Real(0.25)))/Real(0.01);
#else
    Real z = Real(0.);
    Real r2 = ((x-Real(0.25))*(x-Real(0.25))+(y-Real(0.25))*(y-Real(0.25)))/Real(0.01);
#endif
    phi(i,j,k) = Real(1.) + std::exp(-r2);
}


#endif
