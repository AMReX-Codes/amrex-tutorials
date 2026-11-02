#ifndef INCFLOW_ADVECTION_K_H_
#define INCFLOW_ADVECTION_K_H_

#include <AMReX_Array4.H>
#include <AMReX_IntVect.H>
#include <AMReX_Slopes_K.H>

#include <cmath>

// Second-order upwind (MOL) face values of the cell-centered velocity.
namespace incflow {

// Walls on the low and high side of one direction, and the first and last
// cell of the domain in that direction
struct DirInfo
{
    bool wall_lo;
    bool wall_hi;
    int dlo;
    int dhi;
};

// MC-limited slope of component n in direction d at cell iv, one-sided
// next to walls
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real slope (amrex::IntVect const& iv, int n, int d,
                   amrex::Array4<amrex::Real const> const& q, DirInfo const& di)
{
    auto const c = iv.dim3();
    if (d == 0) {
        return amrex::amrex_calc_xslope_extdir(c.x, c.y, c.z, n, 2, q,
                                               di.wall_lo, di.wall_hi, di.dlo, di.dhi);
    }
#if (AMREX_SPACEDIM == 3)
    if (d == 2) {
        return amrex::amrex_calc_zslope_extdir(c.x, c.y, c.z, n, 2, q,
                                               di.wall_lo, di.wall_hi, di.dlo, di.dhi);
    }
#endif
    return amrex::amrex_calc_yslope_extdir(c.x, c.y, c.z, n, 2, q,
                                           di.wall_lo, di.wall_hi, di.dlo, di.dhi);
}

// Values of component n on the left (ql) and right (qr) of the face between
// cells iv-e_d and iv. On a wall, both are the wall value in the ghost cell.
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void face_values (amrex::IntVect const& iv, int n, int d,
                  amrex::Array4<amrex::Real const> const& q, DirInfo const& di,
                  amrex::Real& ql, amrex::Real& qr)
{
    using namespace amrex::literals;
    amrex::IntVect const im = iv - amrex::IntVect::TheDimensionVector(d);
    if (di.wall_lo && iv[d] == di.dlo) {
        ql = qr = q(im,n);
    } else if (di.wall_hi && iv[d] == di.dhi+1) {
        ql = qr = q(iv,n);
    } else {
        ql = q(im,n) + 0.5_rt * slope(im, n, d, q, di);
        qr = q(iv,n) - 0.5_rt * slope(iv, n, d, q, di);
    }
}

// Normal velocity on a face: the upwind value, or zero where the flow
// stagnates
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real normal_velocity (amrex::Real ul, amrex::Real ur)
{
    using namespace amrex::literals;
    constexpr amrex::Real small = 1.e-10_rt;
    if ((ul >= 0.0_rt && ur <= 0.0_rt) || std::abs(ul+ur) < small) {
        return 0.0_rt;
    }
    return (ul+ur > 0.0_rt) ? ul : ur;
}

// Value of a quantity on a face, upwinded with the face velocity u
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real upwind (amrex::Real ql, amrex::Real qr, amrex::Real u)
{
    using namespace amrex::literals;
    constexpr amrex::Real small = 1.e-10_rt;
    if (u > small) { return ql; }
    if (u < -small) { return qr; }
    return 0.5_rt * (ql + qr);
}

}

#endif
