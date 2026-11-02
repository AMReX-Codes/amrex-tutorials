#ifndef CNS_K_H_
#define CNS_K_H_

#include "CNS_index.H"
#include "CNS_parm.H"

#include <AMReX_Array.H>
#include <AMReX_Array4.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Math.H>

#include <cmath>
#include <limits>

// Largest stable time step (for CFL = 1) in cell (i,j,k)
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real
cns_estdt (int i, int j, int k, amrex::Array4<amrex::Real const> const& u,
           amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx,
           Parm const& parm) noexcept
{
    using amrex::Real;

    Real rhoinv = Real(1.0)/amrex::max(u(i,j,k,URHO),parm.smallr);
    Real p = amrex::max((parm.eos_gamma-Real(1.0))*u(i,j,k,UEINT), parm.smallp);
    Real cs = std::sqrt(parm.eos_gamma*p*rhoinv);
    Real dt = std::numeric_limits<Real>::max();
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        Real vel = u(i,j,k,UMX+dir)*rhoinv;
        dt = amrex::min(dt, dx[dir]/(std::abs(vel)+cs));
    }
    return dt;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void
cns_compute_temperature (int i, int j, int k, amrex::Array4<amrex::Real> const& u,
                         Parm const& parm) noexcept
{
    using amrex::Real;

    Real rhoinv = Real(1.0)/u(i,j,k,URHO);
    Real mx = u(i,j,k,UMX);
    Real my = u(i,j,k,UMY);
    Real mz = u(i,j,k,UMZ);
    u(i,j,k,UEINT) = u(i,j,k,UEDEN) - Real(0.5) * rhoinv * (mx*mx+my*my+mz*mz);
    u(i,j,k,UTEMP) = rhoinv * u(i,j,k,UEINT) / parm.cv;
}

#endif
