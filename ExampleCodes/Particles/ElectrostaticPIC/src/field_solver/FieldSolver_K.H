#ifndef FIELD_SOLVER_K_H_
#define FIELD_SOLVER_K_H_

#include <AMReX.H>
#include <AMReX_Array.H>
#include <AMReX_Array4.H>
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

/**
   tmp_mask stores 1 if a ghost cell is either uncovered or outside the
   physical domain.

   On exit, mask will be 1 for cells with ncells of either of those conditions,
   zero otherwise.
 */
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

#endif
