#ifndef CHARGE_DEPOSITION_K_H_
#define CHARGE_DEPOSITION_K_H_

#include <AMReX_Array.H>
#include <AMReX_IntVect.H>
#include <AMReX_REAL.H>

/* This routine computes the charge density due to the particles using cloud-in-cell
 * deposition. The Fab rho is assumed to be node-centered.
 *
 * Arguments:
 *     particles : a pointer to the particle array-of-structs
 *     ns        : the stride length of particle struct (the size of the struct in number of reals)
 *     np        : the number of particles
 *     weights   : the particle weights
 *     charge    : the charge of this particle species
 *     rho       : a Fab that will contain the charge density on exit
 *     lo        : a pointer to the lo corner of this valid box for rho, in index space
 *     hi        : a pointer to the hi corner of this valid box for rho, in index space
 *     plo       : the real position of the left-hand corner of the problem domain
 *     dx        : the mesh spacing
 *     ng        : the number of ghost cells in rho
*/
template <typename P>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void deposit_cic (P const& p, amrex::ParticleReal wp, amrex::Real charge,
                  amrex::Array4<amrex::Real> const& rho,
                  amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                  amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi)
{
    amrex::Real inv_vol = AMREX_D_TERM(dxi[0], *dxi[1], *dxi[2]);
    amrex::Real qp = wp * charge * inv_vol;

    amrex::Real lx = (p.pos(0) - plo[0]) * dxi[0];
    amrex::Real ly = (p.pos(1) - plo[1]) * dxi[1];

    int i = static_cast<int>(amrex::Math::floor(lx));
    int j = static_cast<int>(amrex::Math::floor(ly));

    amrex::Real wx_hi = lx - i;
    amrex::Real wy_hi = ly - j;

    amrex::Real wx_lo = amrex::Real(1.0) - wx_hi;
    amrex::Real wy_lo = amrex::Real(1.0) - wy_hi;

    amrex::Gpu::Atomic::AddNoRet(&rho(i,   j,   0, 0), wx_lo*wy_lo*qp);
    amrex::Gpu::Atomic::AddNoRet(&rho(i,   j+1, 0, 0), wx_lo*wy_hi*qp);
    amrex::Gpu::Atomic::AddNoRet(&rho(i+1, j,   0, 0), wx_hi*wy_lo*qp);
    amrex::Gpu::Atomic::AddNoRet(&rho(i+1, j+1, 0, 0), wx_hi*wy_hi*qp);
}

#endif
