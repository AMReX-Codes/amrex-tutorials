#ifndef CNS_PROB_H_
#define CNS_PROB_H_

#include "CNS_index.H"
#include "CNS_parm.H"

#include <AMReX_Array4.H>
#include <AMReX_Geometry.H>

// Sod shock tube in the x-direction

// Problem parameters (prob.*).  Kernels capture a copy by value.
struct ProbParm
{
    amrex::Real p_l = amrex::Real(1.0);
    amrex::Real p_r = amrex::Real(0.1);
    amrex::Real rho_l = amrex::Real(1.0);
    amrex::Real rho_r = amrex::Real(0.125);
    amrex::Real u_l = amrex::Real(0.0);
    amrex::Real u_r = amrex::Real(0.0);
    amrex::Real x0 = amrex::Real(0.5);  // initial position of the diaphragm
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

    const Real x = geomdata.ProbLo(0) + (Real(i)+Real(0.5))*geomdata.CellSize(0);

    const bool left = x < prob_parm.x0;
    const Real p   = left ? prob_parm.p_l   : prob_parm.p_r;
    const Real rho = left ? prob_parm.rho_l : prob_parm.rho_r;
    const Real u   = left ? prob_parm.u_l   : prob_parm.u_r;

    state(i,j,k,URHO ) = rho;
    state(i,j,k,UMX  ) = rho*u;
    state(i,j,k,UMY  ) = Real(0.0);
    state(i,j,k,UMZ  ) = Real(0.0);
    state(i,j,k,UEDEN) = p/(parm.eos_gamma-Real(1.0)) + Real(0.5)*rho*u*u;
}

#endif
