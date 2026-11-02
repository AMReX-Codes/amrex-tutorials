#ifndef CNS_PROB_H_
#define CNS_PROB_H_

#include "CNS_index.H"
#include "CNS_parm.H"

#include <AMReX_Array4.H>
#include <AMReX_Geometry.H>
#include <AMReX_Math.H>

#include <cmath>

// Rayleigh-Taylor instability: heavy fluid on top of light fluid, with
// gravity in the last direction (z in 3D, y in 2D).

// Problem parameters (prob.*).  Kernels capture a copy by value.
struct ProbParm
{
    amrex::Real rho_1 = amrex::Real(0.5);    // density below the interface
    amrex::Real rho_2 = amrex::Real(2.0);    // density above the interface
    amrex::Real p0_base = amrex::Real(5.0);  // pressure at the bottom
};

// Initial state in cell (i,j,k).  Only rho, momentum and rho*E need to be
// set; rho*e and temperature are computed from them.
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void
cns_initdata (int i, int j, int k, amrex::Array4<amrex::Real> const& state,
              amrex::GeometryData const& geomdata, Parm const& parm,
              ProbParm const& prob_parm) noexcept
{
    using amrex::Real;

    constexpr int vdir = AMREX_SPACEDIM-1;  // vertical direction
    constexpr Real pi = amrex::Math::pi<Real>();

    const auto* plo = geomdata.ProbLo();
    const auto* phi = geomdata.ProbHi();
    const auto* dx  = geomdata.CellSize();
    const amrex::IntVect iv(AMREX_D_DECL(i,j,k));

    Real pos[AMREX_SPACEDIM];
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        pos[d] = plo[d] + (Real(iv[d])+Real(0.5))*dx[d];
    }
    const Real h = pos[vdir];
    const Real hmid = Real(0.5)*(plo[vdir]+phi[vdir]);

    // Hydrostatic pressure (gravity is negative)
    Real p;
    if (h < hmid) {
        p = prob_parm.p0_base + prob_parm.rho_1*parm.gravity*(h-plo[vdir]);
    } else {
        p = prob_parm.p0_base + prob_parm.rho_1*parm.gravity*(hmid-plo[vdir])
            + prob_parm.rho_2*parm.gravity*(h-hmid);
    }

    // Distance from the vertical center line of the domain
    const Real L_x = phi[0] - plo[0];
#if (AMREX_SPACEDIM == 3)
    const Real r = amrex::min(std::hypot(pos[0]-Real(0.5)*(plo[0]+phi[0]),
                                         pos[1]-Real(0.5)*(plo[1]+phi[1])),
                              Real(0.5)*L_x);
#else
    const Real r = std::abs(pos[0]-Real(0.5)*(plo[0]+phi[0]));
#endif

    // Smooth interface with a small cosine perturbation
    const Real pertheight = hmid - Real(0.01)*std::cos(Real(2.0)*pi*r/L_x);
    const Real rho = prob_parm.rho_1 + Real(0.5)*(prob_parm.rho_2-prob_parm.rho_1)
        * (Real(1.0)+std::tanh((h-pertheight)/Real(0.005)));

    state(i,j,k,URHO ) = rho;
    state(i,j,k,UMX  ) = Real(0.0);
    state(i,j,k,UMY  ) = Real(0.0);
    state(i,j,k,UMZ  ) = Real(0.0);
    state(i,j,k,UEDEN) = p/(parm.eos_gamma-Real(1.0));
}

#endif
