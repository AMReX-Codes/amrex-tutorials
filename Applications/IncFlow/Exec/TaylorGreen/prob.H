#ifndef INCFLOW_PROB_H_
#define INCFLOW_PROB_H_

#include <AMReX_Array4.H>
#include <AMReX_Geometry.H>

#include <cmath>

// Decaying Taylor-Green vortices on the periodic unit square, extruded in z
// in 3D. This is an exact solution of the Navier-Stokes equations:
//
//     u =  sin(2 pi x) cos(2 pi y) exp(-8 pi^2 nu t),
//     v = -cos(2 pi x) sin(2 pi y) exp(-8 pi^2 nu t),   w = 0.

class IncFlow;

// Defaults for the domain, grids and boundaries of this problem
void prob_add_defaults ();

// Returns 0 if the result is correct.
int prob_check (IncFlow const& incflow);

// Velocity of cell (i,j,k) at t = 0
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void prob_initdata (int i, int j, int k, amrex::Array4<amrex::Real> const& vel,
                    amrex::GeometryData const& geomdata)
{
    using namespace amrex::literals;
    constexpr amrex::Real twopi = 2.0_rt * amrex::Math::pi<amrex::Real>();
    amrex::Real const x = geomdata.ProbLo(0) + (i+0.5_rt)*geomdata.CellSize(0);
    amrex::Real const y = geomdata.ProbLo(1) + (j+0.5_rt)*geomdata.CellSize(1);
    vel(i,j,k,0) =  std::sin(twopi*x) * std::cos(twopi*y);
    vel(i,j,k,1) = -std::cos(twopi*x) * std::sin(twopi*y);
#if (AMREX_SPACEDIM == 3)
    vel(i,j,k,2) = 0.0_rt;
#endif
}

#endif
