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

#endif
