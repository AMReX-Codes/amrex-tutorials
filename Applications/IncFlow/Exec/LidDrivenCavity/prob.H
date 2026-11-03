#ifndef INCFLOW_PROB_H_
#define INCFLOW_PROB_H_

#include <AMReX_Array4.H>
#include <AMReX_Geometry.H>

// Lid-driven cavity: fluid at rest in the unit square (cube in 3D) with
// no-slip walls. The lid at the top (y = 1) moves in x with velocity 1.

class IncFlow;

// Defaults for the domain, grids and boundaries of this problem
void prob_add_defaults ();

// Returns 0 if the result is sane.
int prob_check (IncFlow const& incflow);

// Velocity of cell (i,j,k) at t = 0
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void prob_initdata (int i, int j, int k, amrex::Array4<amrex::Real> const& vel,
                    amrex::GeometryData const& /*geomdata*/)
{
    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
        vel(i,j,k,n) = amrex::Real(0.0);
    }
}

#endif
