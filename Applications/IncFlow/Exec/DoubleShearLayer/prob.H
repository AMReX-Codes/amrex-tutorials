#ifndef INCFLOW_PROB_H_
#define INCFLOW_PROB_H_

#include <AMReX_Array4.H>
#include <AMReX_Geometry.H>

#include <cmath>

// Double shear layer (Kelvin-Helmholtz) on the periodic unit square,
// extruded in z in 3D. Two thin shear layers at y = 1/4 and y = 3/4 are
// perturbed in v and roll up into vortices, which the fine level follows.
//
//     u = tanh((y - 1/4)/delta)  for y <= 1/2,   tanh((3/4 - y)/delta)  otherwise,
//     v = 0.05 sin(2 pi x),   w = 0,   delta = 1/30.

class IncFlow;

// Defaults for the domain, grids and boundaries of this problem
void prob_add_defaults ();

// Returns 0 if the result is sane.
int prob_check (IncFlow const& incflow);

// Velocity of cell (i,j,k) at t = 0
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void prob_initdata (int i, int j, int k, amrex::Array4<amrex::Real> const& vel,
                    amrex::GeometryData const& geomdata)
{
    using namespace amrex::literals;
    constexpr amrex::Real delta = 1.0_rt/30.0_rt;
    amrex::Real const x = geomdata.ProbLo(0) + (i+0.5_rt)*geomdata.CellSize(0);
    amrex::Real const y = geomdata.ProbLo(1) + (j+0.5_rt)*geomdata.CellSize(1);
    vel(i,j,k,0) = (y <= 0.5_rt) ? std::tanh((y-0.25_rt)/delta) : std::tanh((0.75_rt-y)/delta);
    vel(i,j,k,1) = 0.05_rt * std::sin(2.0_rt*amrex::Math::pi<amrex::Real>()*x);
#if (AMREX_SPACEDIM == 3)
    vel(i,j,k,2) = 0.0_rt;
#endif
}

#endif
