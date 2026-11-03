#ifndef CNS_PARM_H_
#define CNS_PARM_H_

#include <AMReX_REAL.H>

// Physical parameters.  The struct is small and trivially copyable, so GPU
// kernels capture a copy of it by value.
struct Parm
{
    amrex::Real eos_gamma = amrex::Real(1.4);
    amrex::Real eos_mu = amrex::Real(28.97);  // mean molecular weight

    amrex::Real cv = 0;
    amrex::Real cp = 0;

    // Gravitational acceleration in the last coordinate direction
    amrex::Real gravity = 0;

    amrex::Real smallr = amrex::Real(1.e-19);
    amrex::Real smallp = amrex::Real(1.e-10);

    void Initialize ();
};

#endif
