#ifndef ELECTROSTATIC_PIC_2D
#define ELECTROSTATIC_PIC_2D

#include <AMReX_Array.H>
#include <AMReX_IntVect.H>
#include <AMReX_REAL.H>

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void sum_fine_to_crse_nodal (int i, int j, int,
                             amrex::Array4<amrex::Real> const& crse,
                             amrex::Array4<amrex::Real const> const& fine,
                             amrex::IntVect const& ratio) noexcept
{
    const int facx = ratio[0];
    const int facy = ratio[1];
    const int ii = i*facx;
    const int jj = j*facy;

    crse(i,j,0,0) = fine(ii,jj,0,0)                                 +
// These four fine nodes are shared by two coarse nodes...
        amrex::Real(0.5)*(fine(ii-1,jj,0,0) + fine(ii+1,jj,0,0)            +
                          fine(ii,jj-1,0,0) + fine(ii,jj+1,0,0))           +
// ... and these four are shared by four...
        amrex::Real(0.25)*(fine(ii-1,jj-1,0,0) + fine(ii-1,jj+1,0,0)       +
                           fine(ii-1,jj+1,0,0) + fine(ii+1,jj+1,0,0));
// ... note that we have 9 nodes in total...
    crse(i,j,0,0) = crse(i,j,0,0) / amrex::Real(4.0);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void zero_out_bndry (int i, int j, int,
                     amrex::Array4<amrex::Real> const& input,
                     amrex::Array4<amrex::Real> const& bndry,
                     amrex::Array4<int const> mask) noexcept
{
    if (mask(i,j,0,0) == 1) {
        bndry(i,j,0,0) = input(i,j,0,0);
        input(i,j,0,0) = amrex::Real(0.0);
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void build_mask (int i, int j, int,
                 amrex::Array4<int const> tmp_mask,
                 amrex::Array4<int> mask, int ncells) noexcept
{
    int total = 0;
    for (int jj = j-ncells; jj <= j+ncells; ++jj) {
        for (int ii = i-ncells; ii <= i+ncells; ++ii) {
            total += tmp_mask(ii, jj, 0, 0);
        }
    }

    mask(i, j, 0, 0) = (total > 0) ? 1 : 0;
}

/**
 * \brief This routine computes the node-centered electric field given a node-centered phi.
 * The gradient is computed using 2nd-order centered differences. It assumes the
 * boundary conditions have already been set and that you have one row of ghost cells.
 * Note that this routine includes the minus sign in E = - grad phi.
 *
 * Arguments:
 *     lo, hi:     The corners of the valid box over which the gradient is taken
 *     Ex, Ey:     The electric field in the x and y directions.
 *     dx:         The cell spacing
 *
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void compute_E_nodal (int i, int j, int,
                      amrex::Array4<amrex::Real> const& Ex,
                      amrex::Array4<amrex::Real> const& Ey,
                      amrex::Array4<amrex::Real const> const& phi,
                      amrex::GpuArray<amrex::Real, 2> dx) noexcept
{
    amrex::GpuArray<amrex::Real, 2> fac;
    fac[0] = amrex::Real(0.5) / dx[0];
    fac[1] = amrex::Real(0.5) / dx[1];

    Ex(i,j,0,0) = fac[0] * (phi(i-1,j,0,0) - phi(i+1,j,0,0));
    Ey(i,j,0,0) = fac[1] * (phi(i,j-1,0,0) - phi(i,j+1,0,0));
};

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

/* This routine interpolates the electric field to the particle positions
* using cloud-in-cell interpolation. The electric fields are assumed to be
* node-centered.
*
* Arguments:
*     particles : a pointer to the particle array-of-structs
*     ns        : the stride length of particle struct (the size of the struct in number of reals)
*     np        : the number of particles
*     Ex_p      : the electric field in the x-direction at the particle positions (output)
*     Ey_p      : the electric field in the y-direction at the particle positions (output)
*     Ez_p      : the electric field in the z-direction at the particle positions (output)
*     Ex, Ey, Ez: Fabs conting the electric field on the mesh
*     lo        : a pointer to the lo corner of this valid box, in index space
*     hi        : a pointer to the hi corner of this valid box, in index space
*     plo       : the real position of the left-hand corner of the problem domain
*     dx        : the mesh spacing
*     ng        : the number of ghost cells for the E-field
*/
template <typename P>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void interpolate_cic (P const& p, amrex::Real& Ex_p, amrex::Real& Ey_p,
                      amrex::Array4<amrex::Real const> const& Ex,
                      amrex::Array4<amrex::Real const> const& Ey,
                      amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                      amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi)
{
    amrex::Real lx = (p.pos(0) - plo[0]) * dxi[0];
    amrex::Real ly = (p.pos(1) - plo[1]) * dxi[1];

    int i = static_cast<int>(amrex::Math::floor(lx));
    int j = static_cast<int>(amrex::Math::floor(ly));

    amrex::Real wx_hi = lx - i;
    amrex::Real wy_hi = ly - j;

    amrex::Real wx_lo = amrex::Real(1.0) - wx_hi;
    amrex::Real wy_lo = amrex::Real(1.0) - wy_hi;

    Ex_p = wx_lo*wy_lo*Ex(i,   j  ,0,0) +
           wx_lo*wy_hi*Ex(i,   j+1,0,0) +
           wx_hi*wy_lo*Ex(i+1, j  ,0,0) +
           wx_hi*wy_hi*Ex(i+1, j+1,0,0);

    Ey_p = wx_lo*wy_lo*Ey(i,   j  ,0,0) +
           wx_lo*wy_hi*Ey(i,   j+1,0,0) +
           wx_hi*wy_lo*Ey(i+1, j  ,0,0) +
           wx_hi*wy_hi*Ey(i+1, j+1,0,0);
}

#endif
