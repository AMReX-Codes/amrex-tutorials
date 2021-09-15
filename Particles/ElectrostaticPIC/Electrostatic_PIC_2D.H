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
#endif
