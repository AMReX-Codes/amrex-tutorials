#ifndef MY_KERNEL_H_
#define MY_KERNEL_H_

#include <AMReX_FArrayBox.H>

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


AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void compute_flux_x (int i, int j, int k,
                     amrex::Array4<amrex::Real> const& fluxx,
                     amrex::Array4<amrex::Real const> const& phi, amrex::Real dxinv,
                     int lo, int hi, int dom_lo, int dom_hi, int bc_lo, int bc_hi)
{
    if (lo == dom_lo &&
        (bc_lo == BCType::foextrap ||
         bc_lo == BCType::ext_dir))
    {
        if(i == lo)
        {
            fluxx(i,j,k) = (phi(i,j,k)-phi(i-1,j,k)) * dxinv / 0.5;
        }
        else
        {
            fluxx(i,j,k) = (phi(i,j,k)-phi(i-1,j,k)) * dxinv;
        }
    }
    else if (hi == dom_hi &&
             (bc_hi == BCType::foextrap ||
              bc_hi == BCType::ext_dir))
    {
        if(i == hi+1)
        {
            fluxx(i,j,k) = (phi(i,j,k)-phi(i-1,j,k)) * dxinv / 0.5;
        }
        else
        {
            fluxx(i,j,k) = (phi(i,j,k)-phi(i-1,j,k)) * dxinv;
        }
    }
    else
    {
        fluxx(i,j,k) = (phi(i,j,k)-phi(i-1,j,k)) * dxinv;
    }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void compute_flux_y (int i, int j, int k,
                     amrex::Array4<amrex::Real> const& fluxy,
                     amrex::Array4<amrex::Real const> const& phi, amrex::Real dyinv,
                     int lo, int hi, int dom_lo, int dom_hi, int bc_lo, int bc_hi)
{
    if (lo == dom_lo &&
        (bc_lo == BCType::foextrap ||
         bc_lo == BCType::ext_dir))
    {
        if(j == lo)
        {
            fluxy(i,j,k) = (phi(i,j,k)-phi(i,j-1,k)) * dyinv / 0.5;
        }
        else
        {
            fluxy(i,j,k) = (phi(i,j,k)-phi(i,j-1,k)) * dyinv;
        }
    }
    else if (hi == dom_hi &&
             (bc_hi == BCType::foextrap ||
              bc_hi == BCType::ext_dir))
    {
        if(j == hi+1)
        {
            fluxy(i,j,k) = (phi(i,j,k)-phi(i,j-1,k)) * dyinv / 0.5;
        }
        else
        {
            fluxy(i,j,k) = (phi(i,j,k)-phi(i,j-1,k)) * dyinv;
        }
    }
    else
    {
        fluxy(i,j,k) = (phi(i,j,k)-phi(i,j-1,k)) * dyinv;
    }
}


#if (AMREX_SPACEDIM > 2)
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void compute_flux_z (int i, int j, int k,
                     amrex::Array4<amrex::Real> const& fluxz,
                     amrex::Array4<amrex::Real const> const& phi, amrex::Real dzinv,
                     int lo, int hi, int dom_lo, int dom_hi, int bc_lo, int bc_hi)
{
    if (lo == dom_lo &&
        (bc_lo == BCType::foextrap ||
         bc_lo == BCType::ext_dir))
    {
        if(k == lo)
        {
            fluxz(i,j,k) = (phi(i,j,k)-phi(i,j,k-1)) * dzinv / 0.5;
        }
        else
        {
            fluxz(i,j,k) = (phi(i,j,k)-phi(i,j,k-1)) * dzinv;
        }
    }
    else if (hi == dom_hi &&
             (bc_hi == BCType::foextrap ||
              bc_hi == BCType::ext_dir))
    {
        if(k == hi+1)
        {
            fluxz(i,j,k) = (phi(i,j,k)-phi(i,j,k-1)) * dzinv / 0.5;
        }
        else
        {
            fluxz(i,j,k) = (phi(i,j,k)-phi(i,j,k-1)) * dzinv;
        }
    }
    else
    {
        fluxz(i,j,k) = (phi(i,j,k)-phi(i,j,k-1)) * dzinv;
    }
}
#endif

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void update_phi (int i, int j, int k,
                 amrex::Array4<amrex::Real const> const& phiold,
                 amrex::Array4<amrex::Real      > const& phinew,
                 AMREX_D_DECL(amrex::Array4<amrex::Real const> const& fluxx,
                              amrex::Array4<amrex::Real const> const& fluxy,
                              amrex::Array4<amrex::Real const> const& fluxz),
                 amrex::Real dt,
                 AMREX_D_DECL(amrex::Real dxinv,
                              amrex::Real dyinv,
                              amrex::Real dzinv))
{
    phinew(i,j,k) = phiold(i,j,k)
        + dt * dxinv * (fluxx(i+1,j  ,k  ) - fluxx(i,j,k))
        + dt * dyinv * (fluxy(i  ,j+1,k  ) - fluxy(i,j,k))
#if (AMREX_SPACEDIM > 2)
        + dt * dzinv * (fluxz(i  ,j  ,k+1) - fluxz(i,j,k));
#else
        ;
#endif
}

#endif
