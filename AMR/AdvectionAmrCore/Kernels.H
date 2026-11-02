#ifndef ADV_KERNELS_H_
#define ADV_KERNELS_H_

// Problem-specific device functions: initial data, velocity field, and the
// numerical flux. Change these to solve a different problem.

#include <AMReX_Array4.H>
#include <AMReX_Math.H>
#include <AMReX_REAL.H>

#include <cmath>

namespace adv {

using amrex::Real;
using namespace amrex::literals;

// Gaussian centered at (0.5, 0.75), independent of z.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real initial_phi (Real x, Real y) noexcept
{
    Real const r2 = (x-0.5_rt)*(x-0.5_rt) + (y-0.75_rt)*(y-0.75_rt);
    return std::exp(-r2 / 0.01_rt);
}

// Stream function of a swirling flow in the x-y plane. The flow reverses
// at t = 1, so phi returns to its initial state at t = 2.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real streamfunction (Real x, Real y, Real t) noexcept
{
    constexpr Real pi = amrex::Math::pi<Real>();
    Real const sx = std::sin(pi*x);
    Real const sy = std::sin(pi*y);
    return sx*sx * sy*sy * std::cos(0.5_rt*pi*t) / pi;
}

// Normal velocity on face (i,j) of direction dir. It is computed from the
// stream function at the nodes of the face, so its discrete divergence is zero.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real face_velocity (int dir, int i, int j,
                    amrex::GpuArray<Real,AMREX_SPACEDIM> const& problo,
                    amrex::GpuArray<Real,AMREX_SPACEDIM> const& dx, Real t) noexcept
{
    Real const x = problo[0] + i*dx[0];
    Real const y = problo[1] + j*dx[1];
    if (dir == 0) {
        return -(streamfunction(x, y+dx[1], t) - streamfunction(x, y, t)) / dx[1];
    } else if (dir == 1) {
        return  (streamfunction(x+dx[0], y, t) - streamfunction(x, y, t)) / dx[0];
    } else {
        return 0.0_rt;
    }
}

// Monotonized central (MC) limited slope.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real mc_slope (Real qm, Real q0, Real qp) noexcept
{
    Real const dl = q0 - qm;
    Real const dr = qp - q0;
    if (dl*dr <= 0.0_rt) { return 0.0_rt; }
    Real const dc = 0.5_rt*(qp - qm);
    Real const s = amrex::min(2.0_rt*std::abs(dl), 2.0_rt*std::abs(dr), std::abs(dc));
    return (dc > 0.0_rt) ? s : -s;
}

// Upwind flux u*phi on the face between cells (i,j,k)-e and (i,j,k), with
// phi reconstructed linearly from the upwind cell.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real upwind_flux (amrex::Array4<Real const> const& q, int i, int j, int k,
                  amrex::Dim3 const& e, Real u) noexcept
{
    if (u > 0.0_rt) {
        int const il = i-e.x, jl = j-e.y, kl = k-e.z;
        Real const s = mc_slope(q(il-e.x,jl-e.y,kl-e.z), q(il,jl,kl), q(i,j,k));
        return u * (q(il,jl,kl) + 0.5_rt*s);
    } else {
        Real const s = mc_slope(q(i-e.x,j-e.y,k-e.z), q(i,j,k), q(i+e.x,j+e.y,k+e.z));
        return u * (q(i,j,k) - 0.5_rt*s);
    }
}

}

#endif
